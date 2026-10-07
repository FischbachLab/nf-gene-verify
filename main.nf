#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
 * nf-gene-verify -- annotation-free gene presence/absence verification.
 *
 * Nextflow port of gene_verification_docker/scripts/run_pipeline.sh.
 *
 *   BUILD_PROTEIN_DB   NCBI protein reference set for the gene list   (once)
 *   BLASTX             genome translated in 6 frames vs that set   (per sample)
 *   CALL_GENES         presence/absence calling                    (per sample)
 *
 *   --verify_nr adds:
 *   EXTRACT_LOCI       pull each called locus with flanks          (per sample)
 *   NR_VERIFICATION    blastx those loci against nr                (per sample)
 *   FINALIZE_CALLS     merge round-1 calls with nr evidence        (per sample)
 *
 * The shell version rebuilds the protein DB on every invocation. Here it is a
 * single task whose output fans out to every sample -- for a 20-genome run
 * that is 20 NCBI builds replaced by one, and every sample is then scored
 * against an identical reference set.
 *
 * Usage:
 *   nextflow run main.nf --seedfile seedfile.csv --outdir s3://.../gene-verify
 *   nextflow run main.nf --fasta LM088.fasta --project LM088 --outdir results/
 *   nextflow run main.nf --seedfile seedfile.csv --outdir ... --verify_nr true
 */

def helpMessage() {
    log.info """
    =========================================================
     nf-gene-verify : annotation-free gene verification
    =========================================================

    Input -- one of:
      --seedfile  path    CSV with header 'sampleName,fasta'
      --fasta     path    a single assembly (use with --project)

    Required:
      --outdir    path    output directory (local or s3://)

    Options:
      --project   value   label for a single --fasta run   [basename of fasta]
      --genes     path    gene list, one name per line     [bundled genes.txt]
      --verify_nr <true|false>
                          second round: extract loci and verify against nr
                                                           [${params.verify_nr}]
      --nr_db     path    prebuilt nr BLAST db prefix. Empty = NCBI remote
                          BLAST, which is queue-dependent and can take hours
                                                           [${params.nr_db}]
      --flank     int     bp of flanking sequence per locus [${params.flank}]
      --evalue    value   blastx evalue for the screen      [${params.evalue}]
      --retmax    int     NCBI records fetched per gene     [${params.retmax}]

    Seedfile (CSV, one header row):
      sampleName,fasta
      LM088-RCB001,s3://bucket/LM088-RCB001.fasta
      LM089-RCB002,s3://bucket/LM089-RCB002.fasta

      Rows whose sampleName is blank or starts with '#' are skipped.
    """.stripIndent()
}

if (params.help) {
    helpMessage()
    exit 0
}

// ── Input validation ────────────────────────────────────────────────────────
if (!params.outdir) {
    error "ERROR: --outdir is required"
}
if (!params.seedfile && !params.fasta) {
    error "ERROR: provide either --seedfile or --fasta"
}
if (params.seedfile && params.fasta) {
    error "ERROR: --seedfile and --fasta are mutually exclusive"
}

verify_nr = (params.verify_nr?.toString()?.toLowerCase() in ['true', 'yes', '1'])

// The gene list is staged as a file so the exact list used is published with
// the results, rather than being an invisible argument.
genes_f = file(params.genes, checkIfExists: true)

log.info """
─────────────────────────────────────────────
 nf-gene-verify
   input     : ${params.seedfile ?: params.fasta}
   outdir    : ${params.outdir}
   genes     : ${genes_f}
   verify_nr : ${verify_nr}
   nr_db     : ${verify_nr ? (params.nr_db ?: 'NCBI remote BLAST (slow)') : 'n/a'}
─────────────────────────────────────────────
""".stripIndent()

// ============================================================================
//  Step 1: NCBI protein reference set -- ONE task for the whole run
// ============================================================================
process BUILD_PROTEIN_DB {

    tag "${genes_f.name}"

    container params.container

    cpus   2
    memory '4 GB'
    time   '4 h'

    // Not cache = false: this is the one step worth reusing across -resume
    // runs, and it is also the only step that touches the network.
    errorStrategy { task.attempt <= 2 ? 'retry' : 'terminate' }
    maxRetries 2

    publishDir "${params.outdir}/protein_db", mode: 'copy'

    input:
      path genes

    output:
      path "protein_db/combined.fasta", emit: fasta
      path "protein_db/manifest.csv",   emit: manifest
      path "protein_db/*"

    script:
    """
    set -euo pipefail

    # One name per line; blank lines and '#' comments dropped
    mapfile -t GENES < <(grep -v '^[[:space:]]*\$' ${genes} | grep -v '^[[:space:]]*#' | tr -d '\\r')
    echo "genes: \${GENES[*]}"

    python3 /opt/pipeline/build_protein_db.py \\
        --genes "\${GENES[@]}" \\
        --retmax ${params.retmax} \\
        --out protein_db

    echo "reference sequences: \$(grep -c '^>' protein_db/combined.fasta)"
    """
}

// ============================================================================
//  Step 2: blastx screen -- genome translated in 6 frames vs the protein set
// ============================================================================
process BLASTX {

    tag "$sample"

    container params.container

    cpus   { 8 * task.attempt }
    memory { 16.GB * task.attempt }
    time   '8 h'

    errorStrategy { task.exitStatus in [104,134,137,139,143,247] ? 'retry' : 'finish' }
    maxRetries 1

    publishDir "${params.outdir}/${sample}", mode: 'copy'

    input:
      tuple val(sample), path(genome)
      path protein_fasta

    output:
      tuple val(sample), path("blastx_raw.tsv"), emit: hits

    script:
    // -max_target_seqs is deliberately absent: it silently truncates results
    // and can drop whole genes in subject-mode searches. Carried over from the
    // shell pipeline, where that was observed.
    """
    set -euo pipefail

    blastx -query ${genome} \\
           -subject ${protein_fasta} \\
           -evalue ${params.evalue} \\
           -num_threads ${task.cpus} \\
           -outfmt "6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore" \\
           -out blastx_raw.tsv

    echo "blastx hits: \$(wc -l < blastx_raw.tsv)"
    """
}

// ============================================================================
//  Step 3: presence/absence calling
// ============================================================================
process CALL_GENES {

    tag "$sample"

    container params.container

    cpus   2
    memory '8 GB'
    time   '2 h'

    publishDir "${params.outdir}/${sample}", mode: 'copy'

    input:
      tuple val(sample), path(blastx)
      path manifest
      path genes

    output:
      tuple val(sample), path("gene_calls.tsv"), emit: calls

    script:
    """
    set -euo pipefail

    mapfile -t GENES < <(grep -v '^[[:space:]]*\$' ${genes} | grep -v '^[[:space:]]*#' | tr -d '\\r')

    python3 /opt/pipeline/call_genes.py \\
        --blastx ${blastx} \\
        --manifest ${manifest} \\
        --genes "\${GENES[@]}" \\
        --identity-min ${params.identity_min} \\
        --coverage-min ${params.coverage_min} \\
        --gap ${params.gap} \\
        --out gene_calls.tsv
    """
}

// ============================================================================
//  Step 4 (optional): extract every called locus with flanks
// ============================================================================
process EXTRACT_LOCI {

    tag "$sample"

    container params.container

    cpus   2
    memory '8 GB'
    time   '2 h'

    publishDir "${params.outdir}/${sample}", mode: 'copy'

    input:
      tuple val(sample), path(calls), path(genome)

    output:
      tuple val(sample), path("nr_verification"), emit: loci
      // the individual loci, so the nr search can be batched or parallelised
      tuple val(sample), path("nr_verification/*_locus.fasta"), emit: loci_files

    script:
    """
    set -euo pipefail

    python3 /opt/pipeline/extract_loci.py \\
        --calls ${calls} \\
        --genome ${genome} \\
        --flank ${params.flank} \\
        --out nr_verification \\
        --all

    echo "loci extracted: \$(ls -1 nr_verification/*_locus.fasta 2>/dev/null | wc -l)"
    """
}

// ============================================================================
//  Step 5a (local nr): blastx a BATCH of loci against nr
//
//  The shell pipeline runs one blastx per locus, so eight loci scan all ~83 nr
//  volumes eight times. BLAST streams the whole database regardless of how
//  many query sequences it is given, so putting several loci in one query file
//  scans nr ONCE for all of them -- a far bigger saving than running the eight
//  scans concurrently, which only shifts the bottleneck onto EFS throughput.
//
//    params.nr_batch_size = 0  -> every locus in one task  (one nr scan)
//                           1  -> one task per locus       (fully parallel)
//                           N  -> batches of N
//
//  The per-gene raw TSVs this emits are byte-compatible with what
//  run_nr_verification.py writes itself, so NR_COLLECT reuses that script's
//  own parser and its skip-if-raw-exists logic rather than reimplementing it.
// ============================================================================
process NR_BLASTX {

    tag "${sample}:${genes.join(',')}"

    container params.container

    cpus   { 16 * task.attempt }
    memory { 64.GB * task.attempt }
    time   '24 h'

    // Concurrent readers of nr on EFS contend for throughput; cap them.
    maxForks params.nr_max_forks

    errorStrategy { task.attempt <= 2 ? 'retry' : 'ignore' }
    maxRetries 2

    input:
      tuple val(sample), val(genes), path(loci)

    output:
      tuple val(sample), path("raw_out/*_nr_local_raw.tsv"), emit: raw

    script:
    """
    set -euo pipefail
    mkdir -p raw_out

    # BLAST_DB_USE_MMAP=0: memory-mapping a BLAST DB on a network filesystem
    # (EFS/NFS) fails. run_nr_verification.py sets this too.
    export BLAST_DB_USE_MMAP=0

    # Read OUTFMT from the pipeline module instead of duplicating it here, so
    # the columns cannot drift from what parse_local_tsv() expects.
    OUTFMT=\$(python3 -c "import sys; sys.path.insert(0,'/opt/pipeline'); import run_nr_verification as m; print(m.OUTFMT)")
    echo "outfmt: \${OUTFMT}"

    # Consumed by the splitter below; a Groovy list joined here rather than
    # re-derived from filenames, so the gene set is exactly what was batched.
    export BATCH_GENES="${genes.join(',')}"

    cat ${loci} > batch_query.fasta
    echo "query sequences: \$(grep -c '^>' batch_query.fasta)"

    blastx -query batch_query.fasta \\
           -db ${params.nr_db} \\
           -evalue ${params.nr_evalue} \\
           -outfmt "\${OUTFMT}" \\
           -max_target_seqs ${params.nr_max_target_seqs} \\
           -max_hsps ${params.nr_max_hsps} \\
           -num_threads ${task.cpus} \\
           -out batch_raw.tsv

    # Split by qseqid into one raw TSV per gene. extract_loci.py writes headers
    # as '>{gene}_locus_{contig}_{start}-{end}', so the gene is everything
    # before '_locus_'. Genes with no hits still get a header-only file, which
    # is what run_local() produces in that case.
    python3 - "\$OUTFMT" <<'PYEOF'
import os, sys
outfmt = sys.argv[1]
header = "\\t".join(outfmt.split()[1:])

genes = [g.strip() for g in os.environ["BATCH_GENES"].split(",") if g.strip()]
fh = {g: open(f"raw_out/{g}_nr_local_raw.tsv", "w") for g in genes}
for g in genes:
    fh[g].write(header + "\\n")

n = {g: 0 for g in genes}
unmatched = 0
with open("batch_raw.tsv") as src:
    for line in src:
        if not line.strip():
            continue
        q = line.split("\\t", 1)[0]
        g = q.split("_locus_")[0]
        if g in fh:
            fh[g].write(line)
            n[g] += 1
        else:
            unmatched += 1

for g in genes:
    fh[g].close()
    print(f"  {g}: {n[g]} hit row(s)")
if unmatched:
    print(f"  WARNING: {unmatched} row(s) had a qseqid matching no gene in "
          f"this batch", file=sys.stderr)
PYEOF
    """
}

// ============================================================================
//  Step 5b (local nr): parse the per-gene raw TSVs into nr_top_hits.tsv
//
//  run_nr_verification.py skips any locus whose raw TSV already exists, so
//  staging all of them here makes it do zero BLAST work and only parse --
//  reusing its parser rather than duplicating the HSP-grouping logic.
// ============================================================================
process NR_COLLECT {

    tag "$sample"

    container params.container

    cpus   2
    memory '8 GB'
    time   '2 h'

    publishDir "${params.outdir}/${sample}/nr_verification", mode: 'copy'

    input:
      tuple val(sample), path(raws), path(loci_dir)

    output:
      tuple val(sample), path("nr_out/nr_top_hits.tsv"), emit: top_hits
      path "nr_out/*"

    script:
    """
    set -euo pipefail
    mkdir -p nr_out
    cp -f *_nr_local_raw.tsv nr_out/

    echo "raw TSVs staged: \$(ls -1 nr_out/*_nr_local_raw.tsv | wc -l)"

    # Every raw file is already present, so this only parses.
    python3 /opt/pipeline/run_nr_verification.py \\
        --loci-dir ${loci_dir} \\
        --out nr_out \\
        --db ${params.nr_db}

    [ -f nr_out/nr_top_hits.tsv ] || : > nr_out/nr_top_hits.tsv
    echo "nr_top_hits rows: \$(( \$(wc -l < nr_out/nr_top_hits.tsv) - 1 ))"
    """
}

// ============================================================================
//  Step 5 (remote nr): unchanged single task
//
//  Not parallelised: remote mode submits to NCBI's queue and sleeps 10 s
//  between submissions on purpose. The wall clock is NCBI's queue, not our
//  concurrency, and firing eight submissions at once is impolite.
// ============================================================================
process NR_VERIFICATION {

    tag "$sample"

    container params.container

    // Local nr is I/O bound streaming ~80 volumes off EFS; remote mode is
    // mostly waiting on NCBI's queue. Both want a long wall clock.
    cpus   { params.nr_db ? 16 : 2 }
    memory { params.nr_db ? 64.GB : 8.GB }
    time   '48 h'

    errorStrategy { task.attempt <= 2 ? 'retry' : 'ignore' }
    maxRetries 2

    publishDir "${params.outdir}/${sample}/nr_verification", mode: 'copy'

    input:
      tuple val(sample), path(loci_dir)

    output:
      tuple val(sample), path("nr_out/nr_top_hits.tsv"), emit: top_hits
      path "nr_out/*"

    script:
    def db_args = params.nr_db ? "--db ${params.nr_db} --threads ${task.cpus}"
                              : "--hitlist ${params.hitlist}"
    """
    set -euo pipefail
    mkdir -p nr_out

    ${params.nr_db ? "echo 'LOCAL nr: ${params.nr_db}'" : "echo 'REMOTE NCBI BLAST -- queue-dependent'"}

    python3 /opt/pipeline/run_nr_verification.py \\
        --loci-dir ${loci_dir} \\
        --out nr_out \\
        ${db_args}

    # The declared output must exist even when nothing was verified, otherwise
    # Nextflow skips output binding for the whole task and FINALIZE_CALLS hangs.
    [ -f nr_out/nr_top_hits.tsv ] || : > nr_out/nr_top_hits.tsv
    """
}

// ============================================================================
//  Step 6 (optional): merge round-1 calls with nr evidence
// ============================================================================
process FINALIZE_CALLS {

    tag "$sample"

    container params.container

    cpus   2
    memory '8 GB'
    time   '2 h'

    publishDir "${params.outdir}/${sample}", mode: 'copy'

    input:
      tuple val(sample), path(calls), path(nr_hits), path(loci_dir)
      path patterns          // optional title-pattern overrides, may be empty

    output:
      tuple val(sample), path("gene_presence_table_v2.csv"),  emit: final_table
      tuple val(sample), path("gene_presence_matrix.tsv"),    emit: matrix

    script:
    """
    set -euo pipefail

    # bin/finalize_calls.py, NOT /opt/pipeline/: Nextflow prepends the
    # pipeline's bin/ to PATH, so the fixed title-matching logic takes effect
    # without rebuilding fischbachlab/gene-verify. Keep the two in sync --
    # bin/ is a copy of gene_verification_docker/scripts/finalize_calls.py.
    finalize_calls.py \\
        --calls ${calls} \\
        --nr ${nr_hits} \\
        --loci-dir ${loci_dir} \\
        --sample ${sample} \\
        --min-bitscore-frac ${params.nr_min_bitscore_frac} \\
        --expect-genus '${params.expect_genus}' \\
        ${params.gene_patterns ? "--patterns ${file(params.gene_patterns).name}" : ''} \\
        --out gene_presence_table_v2.csv
    """
}

// ============================================================================
//  Workflow
// ============================================================================
workflow {

    // ── Inputs ──────────────────────────────────────────────────────────────
    if (params.seedfile) {
        ch_genomes = Channel
            .fromPath(params.seedfile, checkIfExists: true)
            .splitCsv(header: true)
            .filter { row -> row.sampleName && !row.sampleName.trim().startsWith('#') }
            .map { row ->
                if (!row.fasta?.trim()) {
                    error "ERROR: no fasta for sample '${row.sampleName}' in ${params.seedfile}"
                }
                tuple(row.sampleName.trim(), file(row.fasta.trim(), checkIfExists: true))
            }
    }
    else {
        def f = file(params.fasta, checkIfExists: true)
        def name = params.project ?: f.simpleName
        ch_genomes = Channel.of(tuple(name, f))
    }

    ch_genes = Channel.value(genes_f)

    // ── Step 1: one protein DB for the whole run ────────────────────────────
    BUILD_PROTEIN_DB(ch_genes)

    // first(): turn the single-item output into a value channel so it can be
    // consumed once per sample rather than being used up by the first one
    ch_db_fasta   = BUILD_PROTEIN_DB.out.fasta.first()
    ch_db_manifest = BUILD_PROTEIN_DB.out.manifest.first()

    // ── Steps 2-3 ───────────────────────────────────────────────────────────
    BLASTX(ch_genomes, ch_db_fasta)
    CALL_GENES(BLASTX.out.hits, ch_db_manifest, ch_genes)

    CALL_GENES.out.calls.view { s, f -> "calls: ${s} -> ${f.name}" }

    // ── Steps 4-6: optional nr verification ─────────────────────────────────
    if (verify_nr) {

        // join on sampleName, not positional pairing: samples finish BLASTX in
        // whatever order they finish
        EXTRACT_LOCI(CALL_GENES.out.calls.join(ch_genomes))

        if (params.nr_db) {

            /*
              Local nr. Batch the loci so nr is streamed once per batch rather
              than once per locus:
                nr_batch_size = 0 -> one batch (one scan of nr)
                              = 1 -> one task per locus
                              = N -> batches of N
            */
            def bsize = (params.nr_batch_size as int) > 0 ? (params.nr_batch_size as int)
                                                          : 10000

            ch_batches = EXTRACT_LOCI.out.loci_files
                // a single locus arrives as a bare path, several as a List
                .flatMap { sample, files ->
                    (files instanceof List ? files : [files]).collect { f -> tuple(sample, f) }
                }
                .map { sample, f ->
                    tuple(sample, f.name.replaceAll(/_locus\.fasta$/, ''), f)
                }
                .groupTuple(by: 0)                      // (sample, [genes], [files])
                .flatMap { sample, genes, files ->
                    // split into batches of bsize, keeping genes and files aligned
                    def out = []
                    [genes, files].transpose()
                        .collate(bsize)
                        .each { chunk ->
                            out << tuple(sample,
                                         chunk.collect { it[0] },
                                         chunk.collect { it[1] })
                        }
                    out
                }

            NR_BLASTX(ch_batches)

            // one task may emit several raw TSVs; flatten before grouping
            ch_raw = NR_BLASTX.out.raw
                .groupTuple(by: 0)
                .map { sample, nested -> tuple(sample, nested.flatten()) }

            NR_COLLECT(ch_raw.join(EXTRACT_LOCI.out.loci))
            ch_nr_hits = NR_COLLECT.out.top_hits
        }
        else {
            // Remote NCBI BLAST: one task, as before (see the note on NR_VERIFICATION)
            NR_VERIFICATION(EXTRACT_LOCI.out.loci)
            ch_nr_hits = NR_VERIFICATION.out.top_hits
        }

        // an empty list stages nothing but still satisfies the input
        ch_patterns = params.gene_patterns
            ? Channel.value(file(params.gene_patterns, checkIfExists: true))
            : Channel.value([])

        FINALIZE_CALLS(
            CALL_GENES.out.calls
                .join(ch_nr_hits)
                .join(EXTRACT_LOCI.out.loci),
            ch_patterns
        )

        FINALIZE_CALLS.out.final_table
            .view { s, f -> "FINAL: ${s} -> ${f.name}" }

        /*
          One presence/absence table across all samples. The per-sample files
          already carry a 'sample' column (finalize_calls.py --sample), so
          this is a plain concatenation with one header kept.

          collectFile is an operator, not a process: merging a few hundred
          bytes of TSV does not justify another Batch task.
        */
        FINALIZE_CALLS.out.matrix
            .map { sample, f -> f }
            .collectFile(name: 'gene_presence_matrix.tsv',
                         storeDir: "${params.outdir}/summary",
                         keepHeader: true, skip: 1, sort: true)
            .view { "presence matrix: ${it}" }
    }
}

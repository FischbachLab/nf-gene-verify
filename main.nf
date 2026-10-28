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

// An empty nr_db silently selects REMOTE NCBI BLAST. That is not a slower
// version of the same thing: it submits an RID to a public queue and polls it,
// so a handful of loci can run for days with nothing in the task directory but
// rids.json, and it searches FULL nr rather than the local Staphylococcus
// subset -- different evidence, not just slower evidence. One line in the
// banner was not enough to stop a re-run going out this way, so say it louder.
if (verify_nr && !params.nr_db) {
    log.warn """
    ╔════════════════════════════════════════════════════════════════════╗
    ║  --verify_nr true with NO --nr_db: using REMOTE NCBI BLAST.        ║
    ║  Expect days, not hours, and nr hits from FULL nr rather than the  ║
    ║  local Staphylococcus subset. Pass --nr_db <prefix> for the local  ║
    ║  database, or --verify_nr false to stop after round 1.             ║
    ╚════════════════════════════════════════════════════════════════════╝
    """.stripIndent()
}

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
    // runs, and it is also the only step that touches the network. The cache
    // is keyed on the gene list, the predefined dir AND the content of
    // build_protein_db.py (staged as an input below), so editing the panel
    // logic forces a rebuild while an unchanged rerun still reuses the result.
    errorStrategy { task.attempt <= 2 ? 'retry' : 'terminate' }
    maxRetries 2

    publishDir "${params.outdir}", mode: 'copy'

    input:
      path genes
      path data_dir        // predefined DBs (data/sccmec_proteins.faa); may be empty
      // Frozen per-gene panels; may be empty. Staged (not read from a host
      // path) so the panel is hashed into the task and -resume notices when
      // it changes -- the same reason build_protein_db.py itself is an input.
      path panel_dir
      // Staged as an INPUT rather than picked up from bin/ on PATH so that its
      // CONTENT is part of the task hash. This process keeps its cache across
      // -resume (it is the only networked step), and bin/ changes do not
      // reliably invalidate that cache -- so edits to the panel logic were
      // silently ignored and a whole run reproduced the previous output
      // byte-for-byte, down to the combined.fasta checksum. Hashing the script
      // means a changed filter or alias forces a rebuild, which is the
      // behaviour the comment above assumes.
      path db_script

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

    # Targets that gene-name search cannot retrieve: SCCmec is a cassette,
    # IS256 an insertion sequence. Both come from predefined FASTAs built once
    # by fetch_predefined_db.py.
    PREDEF_ARGS=""
    # args: label, filename, build_protein_db flag
    find_predef() {
        for d in . data; do
            if [ -s "\$d/\$2" ]; then
                PREDEF_ARGS="\$PREDEF_ARGS \$3 \$d/\$2"
                echo "\$1: using predefined \$d/\$2"
                return 0
            fi
        done
        if printf '%s\\n' "\${GENES[@]}" | grep -qix "\$1"; then
            echo "ERROR: '\$1' is in the gene list but \$2 was not staged." >&2
            echo "       Build it once:" >&2
            echo "         python3 bin/fetch_predefined_db.py --db \$(echo \$1 | tr A-Z a-z) --out data/\$2" >&2
            echo "       then point --predefined_dir at that directory." >&2
            return 1
        fi
        return 0
    }
    find_predef SCCmec sccmec_proteins.faa --sccmec-fasta || exit 1
    find_predef IS256  is256_proteins.faa  --is256-fasta  || exit 1

    # The STAGED bin/build_protein_db.py, not /opt/pipeline/: carries the
    # allele-panel query, product filters and diversity selection, so it takes
    # effect without rebuilding fischbachlab/gene-verify. Invoked by its staged
    # path rather than via PATH because it is also a hashed input -- see the
    # input block for why that matters to -resume.
    python3 ${db_script} \\
        --genes "\${GENES[@]}" \\
        --retmax ${params.retmax} \\
        --max-alleles ${params.max_alleles} \\
        --min-len-frac ${params.min_len_frac} \\
        --prefer-organism '${params.prefer_organism}' \\
        ${params.refseq_only ? '--refseq-only' : ''} \\
        ${params.predefined_all ? '--predefined-all' : ''} \\
        ${params.predefined_max ? "--predefined-max ${params.predefined_max}" : ''} \\
        --predefined-min-len ${params.predefined_min_len} \\
        ${params.panel_dir ? "--panel-dir ${file(params.panel_dir).name}" : ''} \\
        \${PREDEF_ARGS} \\
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

    // cpus is NOT multiplied by task.attempt. Doubling it on retry asks for
    // more cores than the node has, so Nextflow refuses to schedule the retry
    // and reports a resource error that hides the ORIGINAL failure -- which is
    // what happened to BLASTX (NIH051668): attempt 2 requested 16 cpus / 32 GB
    // on an 8-core / 28 GB node and surfaced as 'exceeds available memory'.
    //
    // Memory may escalate, since that is the usual reason for a retry, but both
    // steps must fit: 8 -> 16 GB on a 28 GB node. blastx with -subject against
    // a panel under 1 MB needs well under a gigabyte in practice.
    cpus   8
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
    // -max_target_seqs MUST be set explicitly and larger than the panel.
    //
    // Omitting it does NOT disable the limit -- blastx then uses its DEFAULT of
    // 500, and the cap applies per QUERY SEQUENCE. The query here is a whole
    // contig, so a 2.46 Mb chromosome gets at most 500 reference proteins
    // reported no matter how many the panel holds, and the ones dropped are the
    // lowest-scoring.
    //
    // That silently coupled gene detection to PANEL SIZE. Going from 45 genes
    // (743 refs) to 55 genes (908 refs) evicted the weakest subjects and made
    // agrD -- a 46 aa propeptide whose HSPs score 76-90 bits -- vanish from
    // LM088 entirely, despite having matched at 100% identity, E=7e-16, in the
    // smaller panel. ermA's background hits went the same way. Both read as
    // 'absent' with n_hits=0, indistinguishable from a real negative.
    //
    // Sized from the panel itself so it cannot be outgrown.
    """
    set -euo pipefail

    PANEL_N=\$(grep -c '^>' ${protein_fasta})
    echo "panel references: \${PANEL_N}"

    blastx -query ${genome} \\
           -subject ${protein_fasta} \\
           -seg ${params.blastx_seg} \\
           -evalue ${params.evalue} \\
           -max_target_seqs \$(( PANEL_N + 100 )) \\
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
      path calls_script     // hashed, so edits to the thresholds rebuild

    output:
      tuple val(sample), path("gene_calls.tsv"), emit: calls

    script:
    """
    set -euo pipefail

    mapfile -t GENES < <(grep -v '^[[:space:]]*\$' ${genes} | grep -v '^[[:space:]]*#' | tr -d '\\r')

    # The STAGED bin/call_genes.py, not /opt/pipeline/: carries the
    # homology_only threshold. Also a hashed input, so editing it invalidates
    # this task's cache -- see BUILD_PROTEIN_DB for why that matters.
    python3 ${calls_script} \\
        --blastx ${blastx} \\
        --manifest ${manifest} \\
        --genes "\${GENES[@]}" \\
        --identity-min ${params.identity_min} \\
        --coverage-min ${params.coverage_min} \\
        --homology-max ${params.homology_max} \\
        --frame-check-min-bp ${params.frame_check_min_bp} \\
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
      path loci_script      // hashed, so edits to --flank handling rebuild

    output:
      tuple val(sample), path("nr_verification"), emit: loci
      // the individual loci, so the nr search can be batched or parallelised
      // OPTIONAL. A sample can legitimately yield no loci at all: a gene that
      // is absent in round 1 has no coordinates, so there is nothing to
      // extract even with --all. That never happens across 55 genes, because
      // something always has a locus -- but a SINGLE-GENE re-run hits it the
      // moment the gene is absent in one genome, and a non-optional output
      // turned that into 'Missing output file(s) nr_verification/*_locus.fasta'
      // and killed the whole run (blaZ, absent in 6 of 29).
      tuple val(sample), path("nr_verification/*_locus.fasta"), emit: loci_files, optional: true
      // sample-prefixed copies, so loci from DIFFERENT samples can be staged
      // into one task without filename collisions
      path "batch/*_locus.fasta", emit: batch_files, optional: true

    script:
    // --all extracts a locus for EVERY gene, including those already called
    // 'present', so round 2 can check them for false positives -- that is how
    // mecA (round 1: present at 99.6%) and msrA (98.3%) were caught. Without
    // it, round 2 only ever rescues divergent/partial calls and never
    // contradicts a present one.
    def all_flag = params.nr_verify_present ? '--all' : ''
    """
    set -euo pipefail

    python3 ${loci_script} \\
        --calls ${calls} \\
        --genome ${genome} \\
        --flank ${params.flank} \\
        --out nr_verification \\
        ${all_flag}

    echo "loci extracted: \$(ls -1 nr_verification/*_locus.fasta 2>/dev/null | wc -l)"

    # Sample-prefixed copies for cross-sample batching. Both the FILENAME and
    # the FASTA header get the prefix: the filename so Nextflow can stage 30
    # samples' loci side by side, the header so the splitter can route each hit
    # back to the sample it came from. nr_verification/ keeps the ORIGINAL
    # headers, which is what finalize_calls.py reads for window coordinates.
    mkdir -p batch
    for f in nr_verification/*_locus.fasta; do
        [ -e "\$f" ] || continue
        sed "s|^>|>${sample}__|" "\$f" > "batch/${sample}__\$(basename \$f)"
    done
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

    tag "${labels.size()} loci"

    container params.container

    // 4 threads, NOT 32. blastx parallel efficiency against staph_nr, measured
    // on embp (30,987 bp): 100% at 1 thread, 93% at 4, 79% at 8, 44% at 16.
    // Beyond ~8 threads the extra cores do almost nothing, so one wide task
    // wastes the machine -- four concurrent 4-thread tasks deliver 102 bp/s
    // against 48 bp/s for a single 16-thread task, i.e. 8.4 h versus 17.8 h on
    // the same 3.08 Mb query. 32 threads would be worse still.
    //
    // cpus x maxForks should equal the node's core count; raise both together.
    // cpus NOT multiplied by task.attempt -- see BLASTX. 4 cpus x
    // nr_max_forks 2 = the node's 8 cores, at 93% parallel efficiency.
    cpus   8
    memory { 16.GB * task.attempt }

    // Concurrent readers contend for throughput; cap them.
    maxForks params.nr_max_forks

    // 'ignore' made sense when this was one task PER SAMPLE: losing one sample's
    // nr search still left the other 29 usable. With nr_batch_samples it is a
    // SINGLE task covering every sample, so ignoring a failure silently yields
    // no nr results at all -- finalize_calls.py would only print
    // '[WARN] nr hits file not found' and report every gene in every genome as
    // flagged/absent. That is a 30-genome run quietly producing wrong answers,
    // which is worse than a visible crash.
    // Always terminate. Loci are pooled across samples, so ANY lost chunk
    // silently removes nr evidence from many genomes at once -- and
    // finalize_calls.py would only print '[WARN] nr hits file not found'
    // before reporting those genes as flagged/absent.
    errorStrategy { task.attempt <= 2 ? 'retry' : 'terminate' }
    maxRetries 2

    input:
      tuple val(labels), path(loci)
      path nr_script        // hashed, so edits to OUTFMT/scoring rebuild

    output:
      // '<sample>__<gene>_nr_local_raw.tsv'; the workflow routes them by sample
      path "raw_out/*_nr_local_raw.tsv", emit: raw

    script:
    """
    set -euo pipefail
    mkdir -p raw_out

    # BLAST_DB_USE_MMAP=0: memory-mapping a BLAST DB on a network filesystem
    # (EFS/NFS) fails. run_nr_verification.py sets this too.
    export BLAST_DB_USE_MMAP=0

    # Read OUTFMT from the pipeline module instead of duplicating it here, so
    # the columns cannot drift from what parse_local_tsv() expects.
    # The staged copy, not /opt/pipeline/: the script is a hashed input, so it
    # is imported from the task directory by its own basename.
    OUTFMT=\$(python3 -c "import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath('${nr_script}')) or '.'); import ${nr_script.baseName} as m; print(m.OUTFMT)")
    echo "outfmt: \${OUTFMT}"

    # Consumed by the splitter below; a Groovy list joined here rather than
    # re-derived from filenames, so the gene set is exactly what was batched.
    export BATCH_LABELS="${labels.join(',')}"

    cat ${loci} > batch_query.fasta
    echo "query sequences: \$(grep -c '^>' batch_query.fasta)"

    blastx -query batch_query.fasta \\
           -db ${params.nr_db} \\
           -seg ${params.nr_seg} \\
           -evalue ${params.nr_evalue} \\
           -outfmt "\${OUTFMT}" \\
           -max_target_seqs ${params.nr_max_target_seqs} \\
           -max_hsps ${params.nr_max_hsps} \\
           -num_threads ${task.cpus} \\
           -out batch_raw.tsv

    # Split by qseqid into one raw TSV per SAMPLE and gene. The staged headers
    # are '>{sample}__{gene}_locus_{contig}_{start}-{end}'. The sample prefix is
    # STRIPPED from qseqid on the way out, so each raw TSV is byte-identical in
    # form to a single-sample run and run_nr_verification.py -- which takes the
    # gene from the FILENAME and ignores qseqid -- needs no change.
    # Labels with no hits still get a header-only file, as run_local() produces.
    python3 - "\$OUTFMT" <<'PYEOF'
import os, sys
outfmt = sys.argv[1]
header = "\\t".join(outfmt.split()[1:])

labels = [x.strip() for x in os.environ["BATCH_LABELS"].split(",") if x.strip()]
fh = {x: open(f"raw_out/{x}_nr_local_raw.tsv", "w") for x in labels}
for x in labels:
    fh[x].write(header + "\\n")

n = {x: 0 for x in labels}
unmatched = 0
with open("batch_raw.tsv") as src:
    for line in src:
        if not line.strip():
            continue
        f = line.split("\\t")
        q = f[0]
        # '<sample>__<gene>_locus_<contig>_<start>-<end>'
        label = q.split("_locus_")[0]
        if label in fh:
            # strip the sample prefix from qseqid so the raw TSV is identical
            # in form to a single-sample run
            f[0] = q.split("__", 1)[1] if "__" in q else q
            fh[label].write("\\t".join(f))
            n[label] += 1
        else:
            unmatched += 1

for x in labels:
    fh[x].close()
by_sample = {}
for x in labels:
    by_sample.setdefault(x.split("__", 1)[0], 0)
    by_sample[x.split("__", 1)[0]] += n[x]
for smp, tot in sorted(by_sample.items()):
    print(f"  {smp}: {tot} hit row(s) across "
          f"{sum(1 for x in labels if x.startswith(smp + '__'))} loci")
if unmatched:
    print(f"  WARNING: {unmatched} row(s) had a qseqid matching no label in "
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
      path nr_script        // hashed, so edits to the parser rebuild

    output:
      tuple val(sample), path("nr_out/nr_top_hits.tsv"), emit: top_hits
      path "nr_out/*"

    script:
    """
    set -euo pipefail
    mkdir -p nr_out
    # Strip the '<sample>__' prefix added for cross-sample batching:
    # run_nr_verification.py takes the gene from the FILENAME, so it must see
    # '<gene>_nr_local_raw.tsv' exactly as a single-sample run produces.
    for f in *_nr_local_raw.tsv; do
        [ -e "\$f" ] || continue
        cp -f "\$f" "nr_out/\${f#*__}"
    done

    # '|| true': zero raw files is a valid state (the gene was absent in this
    # genome, so no locus was searched), and a bare ls would write to stderr.
    echo "raw TSVs staged: \$(ls -1 nr_out/*_nr_local_raw.tsv 2>/dev/null | wc -l || true)"

    # Every raw file is already present, so this only parses.
    python3 ${nr_script} \\
        --loci-dir ${loci_dir} \\
        --out nr_out \\
        --db ${params.nr_db}

    [ -f nr_out/nr_top_hits.tsv ] || : > nr_out/nr_top_hits.tsv
    # max 0: a header-only or zero-byte file would otherwise report '-1 rows'.
    ROWS=\$(wc -l < nr_out/nr_top_hits.tsv)
    echo "nr_top_hits rows: \$(( ROWS > 0 ? ROWS - 1 : 0 ))"
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
    // Remote-NCBI path only; sized to fit a modest node either way.
    cpus   { params.nr_db ? 8 : 2 }
    memory { params.nr_db ? 16.GB : 8.GB }
    time   '48 h'

    errorStrategy { task.attempt <= 2 ? 'retry' : 'ignore' }
    maxRetries 2

    publishDir "${params.outdir}/${sample}/nr_verification", mode: 'copy'

    input:
      tuple val(sample), path(loci_dir)
      path nr_script        // hashed, so edits to the search rebuild

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

    python3 ${nr_script} \\
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
        --min-contained-frac ${params.nr_min_contained_frac} \\
        --min-overlap-frac ${params.nr_min_overlap_frac} \\
        --short-gene-aa ${params.short_gene_aa} \\
        ${params.nr_verify_present ? '' : '--present-unverified'} \\
        --expect-genus '${params.expect_genus}' \\
        ${params.gene_patterns ? "--patterns ${file(params.gene_patterns).name}" : ''} \\
        --out gene_presence_table_v2.csv
    """
}

// ============================================================================
//  Workflow
// ============================================================================
//  Step 8: wide presence matrix + prevalence summaries
// ============================================================================
process SUMMARIZE_MATRIX {

    tag "wide matrix"

    container params.container

    cpus   1
    memory '2 GB'
    time   '30 m'

    publishDir "${params.outdir}/summary", mode: 'copy'

    input:
      path long_matrix
      path genes
      path summary_script   // hashed, so edits rebuild

    output:
      path "gene_presence_matrix_wide.tsv", emit: wide
      path "gene_prevalence.tsv",           emit: prevalence
      path "genes_per_genome.tsv",          emit: per_genome

    script:
    """
    set -euo pipefail

    # --genes supplies COLUMN ORDER, so operons stay grouped (arc, opp3, kdp,
    # then resistance, then adhesins) rather than being split alphabetically.
    # Reading down a column block is how the ica operon's internal consistency
    # becomes visible: icaA/B/C/D/R all at the same count is a check that the
    # panels agree, and a single column out of step is what a bug looks like.
    python3 ${summary_script} \\
        --matrix ${long_matrix} \\
        --genes ${genes} \\
        --out-prefix .
    """
}

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
    // an empty list stages nothing but still satisfies the input
    ch_predefined = params.predefined_dir
        ? Channel.value(file(params.predefined_dir, checkIfExists: true))
        : Channel.value([])

    // Hashed as an input so edits to the panel logic invalidate the cache
    ch_db_script = Channel.value(
        file("$baseDir/bin/build_protein_db.py", checkIfExists: true))

    // A frozen panel makes the run reproducible: every gene found here is
    // built from file and nothing is fetched, so a re-run cannot resample
    // NCBI. See tools/freeze_panel.sh.
    ch_panel = params.panel_dir
        ? Channel.value(file(params.panel_dir, checkIfExists: true))
        : Channel.value([])

    BUILD_PROTEIN_DB(ch_genes, ch_predefined, ch_panel, ch_db_script)

    // No .first() needed. Every input to BUILD_PROTEIN_DB is a value channel
    // (ch_genes, ch_predefined, ch_db_script), so the process runs once and
    // Nextflow emits its outputs as VALUE channels, which can already be
    // consumed once per sample instead of being used up by the first one.
    // Applying .first() to a value channel is a no-op and Nextflow warns about
    // it. Restore it only if one of those inputs ever becomes a queue channel
    // -- a per-sample gene list, say -- which would make these queue channels
    // and reintroduce the original hazard.
    ch_db_fasta    = BUILD_PROTEIN_DB.out.fasta
    ch_db_manifest = BUILD_PROTEIN_DB.out.manifest

    // ── Steps 2-3 ───────────────────────────────────────────────────────────
    BLASTX(ch_genomes, ch_db_fasta)
    ch_calls_script = Channel.value(
        file("$baseDir/bin/call_genes.py", checkIfExists: true))

    CALL_GENES(BLASTX.out.hits, ch_db_manifest, ch_genes, ch_calls_script)

    CALL_GENES.out.calls.view { s, f -> "calls: ${s} -> ${f.name}" }

    // ── Steps 4-6: optional nr verification ─────────────────────────────────
    if (verify_nr) {

        // Staged from bin/ rather than baked into the image: the Dockerfile no
        // longer ships /opt/pipeline, so these two are hashed inputs like
        // build_protein_db.py and call_genes.py. Editing one invalidates the
        // tasks that use it, which a baked-in copy never did.
        ch_loci_script = Channel.value(
            file("$baseDir/bin/extract_loci.py", checkIfExists: true))
        ch_nr_script = Channel.value(
            file("$baseDir/bin/run_nr_verification.py", checkIfExists: true))

        // join on sampleName, not positional pairing: samples finish BLASTX in
        // whatever order they finish
        EXTRACT_LOCI(CALL_GENES.out.calls.join(ch_genomes), ch_loci_script)

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

            /*
              Loci from ALL samples go into one channel. blastx searches every
              query sequence in a file against the database in a single pass, so
              30 genomes cost ONE scan of nr rather than 30.
              
              This is the dominant cost: the query is ~100 kb per sample against
              a ~100 GB database, so essentially all of the runtime is streaming
              nr. Thresholds (--nr_evalue, --nr_max_target_seqs) do not change
              that -- they are applied while scoring and reporting hits, not to
              decide what is read off disk.
              
              Filenames and FASTA headers are sample-prefixed by EXTRACT_LOCI so
              the loci can be staged together and the hits routed back.
            */
            ch_all_loci = EXTRACT_LOCI.out.batch_files
                .flatMap { files -> (files instanceof List ? files : [files]) }
                .map { f -> tuple(f.name.replaceAll(/_locus\.fasta$/, ''), f) }

            ch_batches = (params.nr_batch_samples
                    ? ch_all_loci.toList()              // one task for everything
                    : ch_all_loci.collate(bsize).toList())
                .flatMap { rows ->
                    // rows: [[label, file], ...] (or [[[label,file],...], ...])
                    def chunks = params.nr_batch_samples ? [rows]
                                                         : rows
                    chunks.findAll { it }
                          .collect { chunk ->
                              tuple(chunk.collect { it[0] }, chunk.collect { it[1] })
                          }
                }

            NR_BLASTX(ch_batches, ch_nr_script)

            // route each '<sample>__<gene>_nr_local_raw.tsv' back to its sample
            ch_raw = NR_BLASTX.out.raw
                .flatMap { files -> (files instanceof List ? files : [files]) }
                .map { f -> tuple(f.name.tokenize('__')[0], f) }
                .groupTuple(by: 0)

            // join on the LOCI side with remainder, so a sample that produced
            // no loci still reaches NR_COLLECT. A plain inner join drops it,
            // and because FINALIZE_CALLS joins on nr hits too, that sample
            // would silently vanish from the results -- 29 genomes in, 23 out,
            // with no error anywhere. It must instead come through with an
            // empty nr_top_hits.tsv and be reported absent.
            ch_collect = EXTRACT_LOCI.out.loci
                .join(ch_raw, remainder: true)
                .map { sample, loci_dir, raws -> tuple(sample, raws ?: [], loci_dir) }

            NR_COLLECT(ch_collect, ch_nr_script)
            ch_nr_hits = NR_COLLECT.out.top_hits
        }
        else {
            // Remote NCBI BLAST: one task, as before (see the note on NR_VERIFICATION)
            NR_VERIFICATION(EXTRACT_LOCI.out.loci, ch_nr_script)
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
        ch_long_matrix = FINALIZE_CALLS.out.matrix
            .map { sample, f -> f }
            .collectFile(name: 'gene_presence_matrix.tsv',
                         storeDir: "${params.outdir}/summary",
                         keepHeader: true, skip: 1, sort: true)

        ch_long_matrix.view { "presence matrix: ${it}" }

        // Derive the wide matrix and the two prevalence summaries from the
        // merged long matrix. Takes the file as an INPUT rather than globbing
        // the summary directory: a glob there would include the file being
        // written, which is how a rebuild loop reads its own output and grows
        // without bound.
        SUMMARIZE_MATRIX(
            ch_long_matrix,
            ch_genes,
            Channel.value(file("$baseDir/bin/summarize_matrix.py",
                               checkIfExists: true))
        )
        SUMMARIZE_MATRIX.out.wide.view { "wide matrix: ${it}" }
    }
}

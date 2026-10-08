# nf-gene-verify

An Annotation-free pipeline for gene presence/absence verification.

## Workflow
```
genes.txt ──► BUILD_PROTEIN_DB ──► combined.fasta + manifest.csv   (ONCE per run)
                                        │
assembly ───────────────────────────────┴──► BLASTX ──► CALL_GENES ──► gene_calls.tsv
                                                                           │
                              --verify_nr true                             ▼
                                   EXTRACT_LOCI ──► NR_BLASTX ──► NR_COLLECT ──► FINALIZE_CALLS
                                                                   gene_presence_table_v2.csv & gene_presence_matrix.tsv
```

## Implementation notes

**The protein DB is built once, not once per genome.** The shell pipeline
rebuilds it on every invocation, so a 20-genome run meant 20 NCBI builds — and
20 reference sets that could differ from each other. Here `BUILD_PROTEIN_DB`
is a single task whose output fans out to every sample, so every genome is
scored against an identical reference and `manifest.csv` describes all of them.

**Samples run in parallel and fail independently.** One genome failing BLASTX
no longer stops the rest.

**Steps are joined on `sampleName`, not on arrival order.** `EXTRACT_LOCI` and
`FINALIZE_CALLS` use `.join()`, so samples finishing out of order can't be
paired with another sample's results.

## Usage

Single genome:

```bash
nextflow run main.nf -profile local \
  --fasta LM088-RCB001.fasta \
  --project LM088-RCB001 \
  --outdir LM088-results \
  --verify_nr true
```

Many genomes on AWS Batch:

```bash
aws batch submit-job \
  --job-name nf-gene-verify \
  --job-queue priority-maf-pipelines \
  --job-definition nextflow-production \
  --container-overrides '{"command":[
      "s3://nextflow-pipelines/nf-gene-verify",
      "--seedfile","s3://genomics-workflow-core/Results/gene-verify/seedfiles/LM088.seedfile.csv",
      "--outdir","s3://enomics-workflow-core/Results/gene-verify/20261005",
      "--genes","s3://genomics-workflow-core/Results/gene-verify/gene_files/genes.txt",
      "--verify_nr","true"]}'
```

### Seedfile format

```csv
sampleName,fasta
LM088,s3://genomics-workflow-core/Results/gene-verify/genomes/LM088-RCB001.fasta
#LM087,s3://genomics-workflow-core/Results/gene-verify/genomes/LM087.fasta
```

Rows with a blank `sampleName` or a leading `#` are skipped.

## Parameters

| Parameter | Default | Description |
|---|---|---|
| `--seedfile` | — | CSV `sampleName,fasta` (or use `--fasta`) |
| `--fasta` | — | a single assembly |
| `--project` | fasta basename | label for a single-`--fasta` run |
| `--outdir` | *required* | local or `s3://` |
| `--genes` | `genes.txt` | gene list, one name per line |
| `--retmax` | `30` | NCBI records fetched per gene |
| `--evalue` | `1e-5` | blastx evalue for the screen |
| `--identity_min` | `90.0` | passed to `call_genes.py` |
| `--coverage_min` | `80.0` | passed to `call_genes.py` |
| `--gap` | `5000` | genomic clustering gap (bp) |
| `--verify_nr` | `false` | run the nr round |
| `--nr_db` | `/mnt/efs/databases/Blast/nr/db/nr` | nr db **prefix**; empty = NCBI remote |
| `--flank` | `300` | bp of flank per extracted locus |
| `--nr_batch_size` | `0` | loci per blastx task; `0` = all in one, `1` = per locus |
| `--nr_max_forks` | `4` | concurrent nr tasks (EFS throughput cap) |
| `--nr_evalue` | `1e-10` | nr round evalue |
| `--nr_max_target_seqs` | `100` | nr round |
| `--nr_max_hsps` | `1` | nr round |

## The nr round is the slow part

The shell pipeline runs **one blastx per locus**, so eight genes scan all ~83
nr volumes **eight times**. BLAST streams the entire database regardless of
how many query sequences it is handed, so the fix is not primarily
parallelism — it is giving one blastx all eight loci at once:

```
nr_batch_size = 0   8 loci -> 1 task  -> nr scanned ONCE     (default)
nr_batch_size = 1   8 loci -> 8 tasks -> nr scanned 8x, concurrently
nr_batch_size = 3   8 loci -> 3 tasks -> nr scanned 3x
```

On EFS the database is read over the network, so the eight concurrent scans of
`nr_batch_size = 1` contend for the same throughput and typically do **not**
run 8× faster — while `nr_batch_size = 0` does ~1/8 the I/O. Start at `0`.

Use a larger batch count only if you have evidence EFS serves concurrent
readers well, and keep `--nr_max_forks` as the brake. On Batch each task also
asks for 16 cores, so `nr_batch_size 1` means 8 tasks × 16 cores across
several instances, all reading nr independently.

`NR_BLASTX` splits the combined output back per gene by `qseqid` — safe
because `extract_loci.py` writes headers as
`>{gene}_locus_{contig}_{start}-{end}`. Genes with no hits still get a
header-only raw TSV, matching what `run_local()` produces. `NR_COLLECT` then
stages those raw files and runs `run_nr_verification.py`, whose
skip-if-raw-exists logic means it does zero BLAST work and only parses — so
the HSP-grouping and coverage maths stay in one place rather than being
reimplemented here.

**Remote mode is not parallelised.** With `--nr_db ''` the original
single-task `NR_VERIFICATION` runs, because the wall clock there is NCBI's
queue, the script sleeps 10 s between submissions deliberately, and firing
eight concurrent submissions is impolite.

## Output

```
<outdir>/
  protein_db/                        combined.fasta, manifest.csv  (shared)
  summary/
    gene_presence_matrix.tsv         sample x gene, long format    (--verify_nr)
  <sample>/
    blastx_raw.tsv                   raw BLAST evidence
    gene_calls.tsv                   round-1 calls
    nr_verification/                 loci + nr hits          (--verify_nr)
    gene_presence_table_v2.csv       full table, all evidence      (--verify_nr)
    gene_presence_matrix.tsv         minimal sample/gene/present   (--verify_nr)
```

### `gene_presence_matrix.tsv`

The minimal result, written per sample and concatenated into
`summary/gene_presence_matrix.tsv` across the run. Long format:

```
sample	gene	present
LM088	agrA	1
LM088	agrB	1
...
```

Only `final_call == present` scores **1**. A `flagged` gene is **0**, because
flagged means *not verified* — scoring an unverified gene as present is how a
false positive reaches a figure. The reason behind any 0 is in
`gene_presence_table_v2.csv`.

Long rather than wide so files from separate genomes concatenate directly;
pivot when a wide table is wanted:

```r
read_tsv("gene_presence_matrix.tsv") |> pivot_wider(names_from = gene, values_from = present)
```

## Profiles

| Profile | Use |
|---|---|
| *(default)* | AWS Batch, queue `default-pipelines` |
| `-profile local` | local executor + docker on the EC2 box, `/mnt/efs` bind-mounted |
| `-profile test` | round 1 only, `--retmax 5` — fast smoke test |

## Notes

- **`--nr_db` is a prefix, not a directory**: `/mnt/efs/databases/Blast/nr/db/nr`,
  not `.../db/`. Leaving it empty switches to NCBI remote BLAST, which is
  queue-dependent and can take hours — hence the 48 h time limit on that step.
- **Local nr is I/O bound, not CPU bound.** It streams ~80 volumes off EFS, so
  the limit is EFS throughput rather than `task.cpus`. Running many samples
  concurrently against it will not help.
- **`BUILD_PROTEIN_DB` keeps its cache** (`cache = true`) while everything else
  sets `cache = false`, matching the other pipelines here. It is the only
  networked step, so `-resume` should not re-hit NCBI.
- **`NR_VERIFICATION` always creates `nr_top_hits.tsv`**, empty if nothing was
  verified. A declared output that matches nothing makes Nextflow skip output
  binding for the entire task, which would leave `FINALIZE_CALLS` waiting on a
  channel that never emits.
- `NCBI_API_KEY` in the environment raises the build step's rate limit from 3
  to 10 requests/second.

## Validation from round 1 

 gene_calls is decided by 3 parameters `--identity-min`, `--coverage-min`, `--gap` in `call_genes.py`.

- **present** — ≥90% protein identity over ≥80% of the reference length
- **divergent** — full coverage but <90% identity (likely genuine ortholog; check the locus)
- **partial** — <80% coverage (fragment, truncation, or cross-homology)
- **absent** — no hit at e-value ≤ 1e-5

## Validation from round 2

There are no identity or coverage thresholds in round 2 — that's the key thing. Those columns are reported but never gate the call.

final_call is decided by three tests in finalize_calls.py:

|#|	test|	threshold
|---|---|---|
|1	|nr hit's query span overlaps the round-1 called region	|any overlap ≥ 1 bp — g2 >= gs and g1 <= ge
|2	|hit title matches the gene's accept pattern and not reject	|pattern match, not numeric
|3	|that hit's bitscore ≥ fraction of the top overlapping hit's	|--min-bitscore-frac, default 0.80 (relative, not absolute)

- **present** - Pass all three 
- **flagged**  - Otherwise


## Notes in nr_top_hits.tsv
|column	|meaning|	computed from
|---|---|---|
|hit_len	|length of the nr subject protein in aa (slen)	|the database entry, not your genome
|bit_score	|bit score of the single best HSP	|max(hs, key=bits)
|evalue	|E-value of that same best HSP	|—
|identity_pct	|length-weighted mean % identity across all HSPs for that subject	|Σ(pident × alen) / Σ(alen)
|align_len	|aligned length in aa of the best HSP only	|—
|query_cov_pct	|always empty in local mode — deliberately unset	|see below
|hit_cov_pct	|% of the subject protein covered, merging overlapping HSPs	|union(subject spans) / slen × 100

### Example output in nr_top_hits.tsv

```bash
gene	rank	accession	title	organism	hit_len	bit_score	evalue	identity_pct	align_len	query_cov_pct	hit_cov_pct	query_start	query_end
sdrH	4	ref|WP_458554421.1|	MSCRAMM-like protein SdrH [Staphylococcus epidermidis]	Staphylococcus epidermidis	436	503.0	7.19E-170	100.0	270		61.9	865	1674


hit_len       436          the SdrH protein in nr is 436 aa
evalue        7.19E-170
bit_score     503.0        best HSP
identity_pct  100.0        perfect identity over what aligned
align_len     270          best HSP covers 270 aa
query_cov_pct              (empty)
hit_cov_pct   61.9         270/436 = 62% of the reference protein
```

### Raw BLAST output for the above sdrH 
```bash
qseqid	sseqid	stitle	staxids	sskingdoms	pident	length	qlen	slen	qstart	qend	sstart	send	evalue	bitscore
sdrH_locus_1_862515-864488	ref|WP_458554421.1|	MSCRAMM-like protein SdrH [Staphylococcus epidermidis]	1282	N/A	100.000	270	1974	436	865	1674	167	436	7.19e-170	503
```


### additional notes

- **The aggregation is mixed. bit_score, evalue and align_len come from the single best HSP; identity_pct and hit_cov_pct are aggregated across all HSPs for that subject. So align_len can be smaller than hit_cov_pct implies when a hit is split into several HSPs — they aren't describing the same thing.

- **query_cov_pct is empty by design, not a bug. The query here is a genomic locus plus 300 bp flanks — far longer than any protein — so query coverage would always be a tiny, meaningless number. Subject coverage answers the real question: how much of the reference protein did we find? That's why finalize_calls.py reads hit_cov_pct or query_cov_pct. The remote-BLAST path does populate query_cov_pct, hence the column existing at all.

- **Units differ across the row. This is blastx, so the query is nucleotide and the subject is protein: hit_len and align_len are in amino acids, while query_start/query_end (865–1674 above) are nucleotide offsets into the extracted locus. finalize_calls.py uses those two to map hits back to genome coordinates.

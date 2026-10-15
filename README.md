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
| `--retmax` | `200` | NCBI records fetched per gene before selection |
| `--max_alleles` | `8` | references kept per gene after diversity selection |
| `--refseq_only` | `false` | restrict round 1 to RefSeq (collapses alleles) |
| `--predefined_dir` | `$baseDir/data` | holds `sccmec_proteins.faa` |
| `--predefined_all` | `false` | use all predefined proteins (`core=0` too), not just mec/ccr core |
| `--predefined_max` | `0` | cap on proteins kept from a predefined set (0 = no cap) |
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

### `gene_presence_table_v2.csv` — column reference

15 columns. Columns 2–6 come from **round 1** (blastx vs the allele panel);
columns 7–13 come from **round 2** (blastx vs nr). Nothing is mixed between
rounds, so the two can be compared directly — that comparison is the point of
the table.

| # | Column | Round | Meaning |
|---|---|---|---|
| 1 | `gene` | — | gene symbol as given in `genes.txt` |
| 2 | `round1_call` | 1 | `present` / `partial` / `divergent` / `no_reference` |
| 3 | `round1_identity_pct` | 1 | % identity to the best panel reference |
| 4 | `round1_coverage_pct` | 1 | % of that reference covered |
| 5 | `contig` | 1 | contig carrying the locus |
| 6 | `coordinates` | 1 | locus start–end on that contig (bp) |
| 7 | `nr_confirming_accession` | 2 | the nr hit the call **rests on** |
| 8 | `nr_confirming_title` | 2 | its defline, truncated to 80 chars |
| 9 | `nr_confirming_organism` | 2 | its source organism |
| 10 | `nr_identity_pct` | 2 | % identity to that hit |
| 11 | `nr_hit_coverage_pct` | 2 | coverage — read with column 12 |
| 12 | `nr_coverage_basis` | 2 | `subject` or `query_locus+flank` |
| 13 | `nr_hit_rank` | 2 | that hit's rank in the nr output |
| 14 | `nr_verification` | 2 | human-readable verdict + any notes |
| 15 | `final_call` | 2 | `present` / `flagged` |

Three things that are easy to misread:

**Columns 7–9 are not the top hit.** They hold the *confirming* hit — the
highest-scoring hit whose title `GENE_PATTERNS` accepts as the gene in
question. When the top hits carry ambiguous or paralogous titles they are
skipped, so the reported hit may be rank 2, 3 or 4. The columns were named
`nr_top_hit_*` before and that was wrong; anything parsing the old names needs
updating. Reporting a top hit the pipeline had already rejected would be
actively misleading, which is why the confirming hit is what appears.

**`nr_hit_rank` is the tell.** `1` means the confirming hit was also the
top-scoring one. Anything above `1` means lower-ranked hits were accepted after
higher ones were rejected, and `nr_verification` spells it out
(`[hit rank 3, not the top hit]`). Rank > 1 is normal, not a warning — it is
the `--nr_min_bitscore_frac` fallback working.

**`nr_coverage_basis` disambiguates column 11.** `subject` = fraction of the
*nr reference protein* aligned. `query_locus+flank` = fraction of the *query*
aligned — and since the query carries `--flank` bp on each side, that figure
caps well below 100% even for a flawless full-length match (1374/1974 = 70% for
sdrH on LM088-RCB001). The two are not comparable, so check this column before
comparing coverage across genes. In practice `subject` is what you will see;
the fallback fires only when `hit_cov_pct` is absent.

Worked example — sdrH on LM088-RCB001:

```
gene                      sdrH
round1_call               present
round1_identity_pct       95.4        <- full-length identity to the panel
round1_coverage_pct       100.0
coordinates               862815-864188
nr_confirming_accession   ref|WP_002494639.1|
nr_confirming_title       MSCRAMM-like protein SdrH [Staphylococcus epidermidis]
nr_identity_pct           100.0       <- exact match
nr_hit_coverage_pct       100.0
nr_coverage_basis         subject     <- 100% of WP_002494639.1, not of the locus
nr_hit_rank               3           <- ranks 1-2 had ambiguous titles
final_call                present
```

Round 1 reads 95.4% and round 2 reads 100.0% for the same gene because the two
rounds search different databases: 95.4% is the distance to the closest allele
*in the panel*, 100.0% the distance to the closest protein *in nr*. A gap
between them means the panel lacks this strain's allele — not that the gene is
questionable. Both calls are `present`, and that agreement is the result to
read.

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

<!-->
## Round 1 uses an allele PANEL

A single reference per gene measures similarity to one strain's allele, so a
gene that varies between isolates reads as `divergent` or `partial` while
plainly present. On LM088-RCB001 that happened to half the gene set:

| gene | round1 | r1 identity | r1 coverage | nr identity | Δ |
|---|---|---|---|---|---|
| agrB | divergent | 49.2% | 96.8% | 99.5% | **+50.3** |
| fosB | divergent | 64.3% | 99.3% | 98.6% | **+34.3** |
| sdrF | partial | 42.9% | 66.3% | 99.9% | **+57.0** |
| sdrH | partial | 93.9% | **79.0%** | 100.0% | +6.1 |

The genome is identical in both rounds, so a +57-point identity jump can only
come from the reference changing — not from the genome, and not from the tool.
Causes: *agrB*/*agrD* differ by **agr type** (I–IV); SdrF varies
**780–1750 aa** by SD-repeat copy number, so a long reference yields low
coverage against a short allele; and `sdrH` missed the 80% coverage cutoff by
one point.

`build_protein_db.py` therefore builds a panel, not a single reference:

- **Product-name queries as well as `[Gene Name]`.** Most NCBI records for
  these genes carry a product name and no symbol — `accessory gene regulator B`
  contains no `agrB` token — so a symbol-only query retrieves an
  unrepresentative slice.
- **No RefSeq restriction by default** (`--refseq-only` restores it). RefSeq
  `WP_` records are *non-redundant representatives*: they deliberately collapse
  the allelic diversity a panel needs. RefSeq is still ranked first.
- **`--retmax 200`, then diversity selection down to `--max-alleles 8`.**
  Greedy farthest-point over combined **sequence and length** dissimilarity.
  The length term is essential, not cosmetic: k-mer similarity barely moves
  when a protein is extended by a low-complexity repeat, so an SD-expanded
  SdrF looks similar to a short allele — yet length is exactly what coverage
  is measured against.

`call_genes.py` needs no change: `score_cluster()` already picks the reference
with the best `(coverage, identity)` per gene, so the best-matching allele
wins automatically.

**Expected effect:** `agrB`, `fosB`, `sdrF` and `sdrH` should return `present`
in round 1, making round 2 a confirmation rather than a correction — so the
two columns stop disagreeing.

## SCCmec: a predefined cassette database

SCCmec is a mobile cassette, not a gene, so `SCCmec[Gene Name]` retrieves
nothing useful. Its references are the CDS translations of seven
whole-cassette GenBank records, prebuilt once and shipped as a file:

```
AB033763  D86934  AB037671  AB063172  AB121219  AF411935  FJ390057
```

Build it once, then commit it:

```bash
python3 bin/fetch_sccmec_db.py --out data/sccmec_proteins.faa
```

Add `SCCmec` to `genes.txt` and the builder uses that file instead of querying
NCBI. If `SCCmec` is listed but the file is missing, `BUILD_PROTEIN_DB` fails
immediately with the command to build it — rather than silently reporting the
cassette absent.

### Only the diagnostic core is used

Every protein is tagged `[core=0|1]`; **only `core=1` is used by default**.

A full cassette carries dozens of CDS — transposases, recombinase fragments,
hypothetical proteins — that also occur elsewhere in staphylococcal genomes and
on unrelated mobile elements. Scoring those under one `SCCmec` label would
report the cassette present in any genome carrying a common transposase. The
core is what SCCmec typing actually rests on:

| complex | proteins |
|---|---|
| *mec* | MecA, MecC, MecR1, MecI |
| *ccr* | CcrA, CcrB, CcrC |

**To include everything (`core=0` and `core=1`): `--predefined_all true`.**
The fetcher always writes every protein; the filter lives in
`build_protein_db.py`, not in `fetch_predefined_db.py`. `--max-alleles` does
not trim predefined sets — use `--predefined_max N` for an explicit cap.

### What `SCCmec present` means

All cassette proteins share one gene label, and `call_genes.py` scores the
best-covered single reference per gene. So `present` means **at least one core
protein matched at ≥90% identity over ≥80% of its length** — usually MecA,
which is highly conserved.

It does **not** establish which SCCmec type is present, nor that the cassette
is intact. For typing, look at which *mec* and *ccr* proteins hit in
`blastx_raw.tsv`, or use a purpose-built tool such as SCCmecFinder.
-->
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

<!--
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
-->
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

### Reproducible runs: freeze the reference panel

By default each run rebuilds the panel from a genus-wide NCBI query, and that
query drifts. Between two runs three weeks apart, **16 of 55 panels drew
different references** (`fdh` swapped 7 of 20, `aap` 5, `kdpB` 4) and **13
calls moved for that reason alone**, with no code change behind any of them.
Re-running is therefore not a way to confirm a result — without a frozen panel
it produces a new one that needs its own audit.

To make a run reproducible, freeze the panel once and point later runs at it:

```bash
tools/freeze_panel.sh results-29genomes_v5 panels/v5     # once
nextflow run main.nf -profile local --panel_dir panels/v5 ...
```

Every gene with a `<gene>.fasta` in that directory is built from the file and
**nothing is fetched from NCBI**; the rebuild is byte-identical to the run it
was frozen from, manifest included. Genes *not* in the directory fall through
to the normal query, so adding a gene to `genes.txt` without adding it to the
panel silently reintroduces drift for that gene — the build log prints
`FROZEN panel` for each gene served from file, so check the count matches.

`panels/v5/` is the panel behind the 29-genome v5 results (55 genes, 903
references). Leave `--panel_dir` unset when you *want* to pick up new alleles,
and re-freeze afterwards.

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
| `--gap` | `2000` | genomic clustering gap (bp) — see below |
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
| 15 | `frame_status` | 1 | empty, or `frameshift:Nframes` — the called region is not one ORF |
| 16 | `frameshift_at` | 1 | genomic coordinate(s) of the junction(s) |
| 17 | `final_call` | 2 | `present` / `present_frameshift` / `absent` / `unverifiable_short` / `flagged` / `no_reference` |

### `present_frameshift`: sequence verified, reading frame not

Round 1 measures sequence similarity; round 2 asks what protein that sequence
resembles. **Neither looks at the reading frame**, so both can be fully
satisfied by a pseudogene. `blastx` translates the query in six frames and
reports HSPs from any of them, and `call_genes.py` clusters by genomic distance
without checking frame.

On LM088-RCB001, `ebh`/`embp` scores 99.6% identity to an *S. epidermidis*
Ebh-like protein across 30,469 bp — and is interrupted by two frameshifts:

```
frame +0   1,365,187-1,382,946   17,760 bp
frame +2   1,382,916-1,392,902    9,987 bp
frame +1   1,392,833-1,395,655    2,823 bp
```

Two signals are required before reporting this, and the second is what matters:

1. the frames occupy **separate stretches** of the locus, so it is not one ORF;
2. their junctions **overlap** rather than gap.

At a frameshift, references align across the break in either frame, so HSPs
from both frames cover the same bases. **Adjacent genes give disjoint blocks**
instead — `SCCmec`'s two blocks sit 1 bp apart (they are blaR1 and blaI, which
is why that locus was never SCCmec) and `opp3D`'s 86 bp apart (two neighbouring
ABC transporters). Only overlapping junctions are reported. Interleaved frames
are ignored entirely: repeat proteins routinely produce alternative alignments
of the same region in several frames, which says nothing about the gene's frame.

The call scores **0** in the matrix. The matrix answers "does this genome make
this protein", and a locus interrupted by frameshifts does not. The
sequence-level finding is not lost — the row keeps
`final_call = present_frameshift` and the junction coordinates in
`frameshift_at`, so a sequence-level question is still answerable from
`gene_presence_table_v2.csv`. Scoring it 1 would place a probable pseudogene in
a presence figure as an ordinary positive, which is the harder error to catch
downstream.

**What it cannot tell you:** whether the frameshifts are a real pseudogene or
residual assembly indels. Long-read assemblies polished with short reads
(`polypolish`) retain exactly this error class, concentrated in long repeat
arrays. Tempting but invalid reasoning: "every other repeat protein in the
genome is single-frame, so this must be real." At one frameshift per ~15 kb,
the next-longest target (`aap`, 4,143 bp) would be expected to carry 0.27 — so
zero there is entirely consistent with the same length-proportional error rate.
Only read-level evidence at the two junction coordinates settles it, and the
check reduces that from a 30 kb region to two positions.

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

Only `final_call == present` scores **1**; everything else is **0**. The reason
behind any 0 is in `gene_presence_table_v2.csv`, and the four values are not
interchangeable:

| `final_call` | Meaning | Confidence in the 0 |
|---|---|---|
| `present` | verified against nr, **or** too short for nr to assess with a round-1 call of present | — (scores 1) |
| `present_frameshift` | sequence verified; locus spans multiple reading frames | **scores 0** — gene present, functional protein not supported |
| `absent` | searched against a valid panel, nothing qualified, **or** round 1 was `homology_only` and round 2 agreed | **high** — a completed negative |
| `unverifiable_short` | too short for nr to assess, and round 1 did **not** say present | **none** — no positive evidence either way |
| `flagged` | a locus was called but could not be verified | **none** — needs review |
| `no_reference` | no panel was built; untested | **none** — neither present nor absent |

`flagged` is 0 because *not verified* must never score as present — that is how
a false positive reaches a figure. But do not read it as evidence of absence:
`mecA` on LM088-RCB001 was `flagged` because its locus turned out to be a
native PBP, whereas `icaD` was `absent` because the gene is not in the genome.
Only the second is a result you can publish; the first is a to-do item.

`absent` is strongest on a closed assembly, where there is no gap for a gene to
hide in. On a draft assembly every `absent` silently carries "or it fell in a
gap between contigs".

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

<!--
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

**`--gap` affects coordinates, not calls.** It was lowered from 5000 to 2000.
At 5000 the `aap` cluster on LM088-RCB001 absorbed a detached 516 bp block of
G5/E repeats lying 3,847 bp away, reporting the locus as 33856-42361 (8,506 bp)
for a gene that is 38219-42361 (4,143 bp) — a 4.4 kb overstatement, with round-1
coverage inflated to 85.7% from 81.2%. At 2000 the coordinates match nr's
confirming hit exactly, which is independent corroboration rather than a tuned
guess.

Validated across two panels and 56 gene-calls: `aap`'s coordinates are the only
thing that changes, and **no call changes anywhere**. The 30,529 bp `ebh`/`embp`
gene — the obvious casualty of a tighter gap — is unaffected at every value
tested down to 500.

Do not go below ~1400. At 1000 the `aap` gene splits at its own internal 1,341
bp gap and loses a block. The safe window is wider than the largest *internal*
unaligned stretch of a gene and narrower than the distance to the nearest
unrelated element. Repeat proteins set both bounds, because their references
match in patches.

- **present** — ≥90% protein identity over ≥80% of the reference length
- **divergent** — full coverage but <90% identity (likely genuine ortholog; check the locus)
- **partial** — <80% coverage (fragment, truncation, or cross-homology)
- **homology_only** — <`--homology-max` identity (default 40%). Checked *first*,
  so it takes precedence over divergent/partial.
- **absent** — no hit at e-value ≤ 1e-5

### Why `homology_only` exists

`divergent` was conflating two unrelated situations. On LM088-RCB001 it gave the
same label to `arcB` at 83.1% — a real cross-species orthologue — and to `tetK`
at 24.9%, which is fold-level similarity to a native MFS efflux pump. Read
without the identity column, that table suggests eight divergent resistance
genes where there are none.

Below the threshold the hit is to a different member of the same **superfamily**,
not an allele of the gene:

| gene | id% | what it actually matched |
|---|---|---|
| `ermA`/`ermB`/`ermC` | 24–27 | one native rRNA methyltransferase, same locus for all three |
| `msrA` | 26.4 | a native ABC-F protein (the resistance Msr(A) is 488 aa; the housekeeping MsrA is ~170 aa) |
| `tetK` | 24.9 | a native MFS pump |
| `opp3D`/`opp3E` | 28–30 | unrelated ABC transporters, 435 kb and 161 kb from the opp3 cluster |
| `kdpB`/`kdpE` | 26–34 | the P-type ATPase and OmpR response-regulator superfamilies |
| `dfrG` | 39.5 | the native trimethoprim-**sensitive** DfrC, at the same coordinates |

For scale, a genuine copy matches these same panels at ~100%: `opp3B` scored
100.0% over 308/308 aa.

**40% is a judgement call, not a natural boundary.** `embp` at 42.2% is a real
*S. aureus*/*S. epidermidis* orthologue (nr-confirmed at 99.6%) and `dfrG` at
39.5% is background — 2.7 points apart. Raise the threshold and real orthologues
get mislabelled; lower it and background leaks back in. Round 2 remains the
arbiter: it independently rejected `mecA` and `msrA`. Set `--homology_max 0` to
disable.

When round 1 says `homology_only` **and** round 2 finds no overlapping nr hit
naming the gene, `final_call` is `absent`, not `flagged` — two independent lines
of evidence agreeing. `flagged` stays reserved for the genuinely unresolved.

## Validation from round 2

There are no identity or coverage thresholds in round 2 — that's the key thing. Those columns are reported but never gate the call.

final_call is decided by three tests in finalize_calls.py:

|#|	test|	threshold
|---|---|---|
|1	|nr hit's query span overlaps the round-1 called region	|any overlap ≥ 1 bp — g2 >= gs and g1 <= ge
|2	|hit title matches the gene's accept pattern and not reject	|pattern match, not numeric
|3	|that hit's bitscore ≥ fraction of the top overlapping hit's	|--min-bitscore-frac, default 0.80 (relative, not absolute)

- **present** - Pass all three
- **flagged**  - A locus WAS called but fails one of the three
- **absent** - Either (a) round 1 found no locus, so there was nothing for round
  2 to verify — round 2 only asks whether an nr hit overlaps the called region,
  and with no region nothing can overlap, so reporting that as `flagged` would
  imply a verification was attempted; or (b) round 1 called `homology_only` and
  round 2 found no overlapping nr hit naming the gene, which is two independent
  lines of evidence agreeing on absence. Note that (b) can be *overturned*: a
  confirming nr hit still wins, so `homology_only` in round 1 does not
  pre-empt round 2.
- **no_reference** - No panel was built for the gene (NCBI returned nothing, or
  `PRODUCT_FILTERS` rejected every record). Untested, so neither present nor
  absent.

Note the asymmetry: round 2 can only ever *confirm* presence, never establish
absence. Every `absent` call is round 1's verdict, resting on the panel holding
the right protein at the right length — which is why `PRODUCT_FILTERS` and the
length floors matter as much as the thresholds do. A `mecA` panel full of the
~240 aa ClpC adaptor instead of the ~668 aa PBP2a produced a confident
`present` on a methicillin-sensitive genome, and no amount of round-2 logic
would have caught it if the title had been less obviously wrong.


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
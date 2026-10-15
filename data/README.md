# Predefined reference databases

Targets that `'<gene>[Gene Name]'` cannot retrieve, because they are not genes.
Build each once with `fetch_predefined_db.py` and commit the result;
`build_protein_db.py` then uses the file whenever the target appears in the
gene list, with no network call.

```bash
python3 ../scripts/fetch_predefined_db.py --db all --outdir .
```

## `sccmec_proteins.faa`

SCCmec protein references — the CDS translations of seven whole-cassette
GenBank records:

```
AB033763  D86934  AB037671  AB063172  AB121219  AF411935  FJ390057
```

SCCmec is a mobile cassette, not a gene, so `SCCmec[Gene Name]` retrieves
nothing useful. Build this file once, commit it, and `build_protein_db.py`
will use it whenever `SCCmec` appears in the gene list — no network call.

```bash
python3 ../scripts/fetch_predefined_db.py --db sccmec --out sccmec_proteins.faa
```

Also writes `sccmec_proteins.manifest.json` recording every protein with its
source accession, gene, product and core flag.

### The diagnostic core

Each protein is tagged `[core=0|1]`. **Only `core=1` is used by default.**

A full cassette carries dozens of CDS — transposases, recombinase fragments,
hypothetical proteins — that also occur elsewhere in staphylococcal genomes
and on unrelated mobile elements. Scoring those under the single `SCCmec`
label would report the cassette present in any genome carrying a common
transposase. The core is the basis SCCmec typing actually rests on:

| complex | proteins |
|---|---|
| *mec* | MecA, MecC, MecR1, MecI |
| *ccr* | CcrA, CcrB, CcrC |

**To include everything (`core=0` and `core=1`):**

```bash
python3 ../scripts/build_protein_db.py --genes SCCmec --predefined-all --out protein_db
# Nextflow:  --predefined_all true
```

The fetcher always writes **every** protein — the filter lives in
`build_protein_db.py`, not in `fetch_predefined_db.py`. Expect more false
positives with `--predefined-all`, for the reason above.

`--max-alleles` does **not** trim predefined sets: it exists to reduce hundreds
of near-identical NCBI records to a diverse handful, and a predefined database
is already curated. Use `--predefined-max N` if you need a cap.

### Interpreting the call

All cassette proteins share one gene label, and `call_genes.py` scores the
best-covered single reference per gene. So `SCCmec present` means **at least
one core protein was found at ≥90% identity over ≥80% of its length** — most
often MecA, which is highly conserved. It does **not** establish which SCCmec
type is present, nor that the cassette is complete. For typing, inspect which
`ccr` and `mec` complex proteins hit in `blastx_raw.tsv`, or use a dedicated
tool such as SCCmecFinder.


## `is256_proteins.faa`

IS256 insertion-sequence references. IS256 is a mobile element, and its only
verifiable product is the transposase, so the references are IS256-family
transposase proteins from *Staphylococcus epidermidis*:

```bash
python3 ../scripts/fetch_predefined_db.py --db is256 --out is256_proteins.faa
```

Records are filtered to **250–500 aa** — IS256-family transposases run ~390 aa,
so anything outside that band is a fragment or a mis-annotated neighbour (an
IS1272 transposase, for instance). All are tagged `core=1`; there is no
core/non-core split, since the element has one protein.

### Interpreting the call

**IS256 is multi-copy.** `call_genes.py` clusters hits by genomic distance and
reports one call per gene, so `IS256 present` means *at least one* intact copy
was found — it says nothing about **copy number**, which is usually the
biologically interesting quantity (IS256 copy number correlates with biofilm
phase variation and is strain-variable). To count copies, cluster the `IS256`
rows of `blastx_raw.tsv` by contig and coordinate.

Also note IS256-family transposases are widespread beyond *S. epidermidis*;
restricting the reference set to *S. epidermidis* keeps the references
appropriate but does not stop a related family member in the query genome from
matching.

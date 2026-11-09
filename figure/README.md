# Fig. 1-style bacterial genome pipeline for EC2

For the included **29-genome, 44-marker dataset**, start with `RUN_29.md` and
`run_29.sh`. The sections below document the general pipeline.
The optional `Dockerfile`, `DOCKER.md` and `run_29_docker.sh` provide the same
analysis in a container without a host Conda installation.
To improve an already completed figure without rebuilding its tree, see
`REDRAW.md` and run `bash redraw_29.sh`.
For the 2026-11-04 memo grouping, NIHLM088 highlight, nosocomial emphasis and
branch distances on the completed **31-genome run**, run `bash redraw_updated.sh`;
this reads `results_31/` and writes `results_31_updated/`.
NIHLM088 is bold blue; Biofilm and Virulence panels are highlighted pale yellow.
The current redraw uses larger text/symbols and compact spacing. Seven node
rotations place NIHLM088 in row 9 and all nine nosocomial isolates within
the bottom ten rows. Clades,
genetic distances and gene calls remain unchanged; see `REDRAW.md`.
For a Conlan et al. Figure 6-style PCA of genome-wide gene-cluster content,
see `PCA.md` and run `bash run_pca.sh --metadata metadata2.tsv`. This reuses
the completed Panaroo matrix and does not repeat annotation or clustering.

This package processes a folder of assembled bacterial genomes and draws a
phylogenetic tree next to sequence types, optional sample metadata and your
reported gene-presence calls. It accepts 30 genomes without hard-coded sample
names. Genome assembly from raw sequencing reads is outside this workflow.

The analysis route is **Prokka -> Panaroo/MAFFT -> IQ-TREE -> plotting**.
The paper uses BV-BRC PGFam/codon trees with RAxML, so this is a layout adaptation
with a different tree-building method. An existing BV-BRC Newick tree can also
be supplied to reproduce that analysis route's tree panel.

## 1. Install on EC2

The bundled installer and dependency lock target **Linux x86_64 (Intel/AMD)**
with glibc **2.28 or newer**. They do not target macOS or ARM/Graviton.
The explicit environment was resolved for Linux x86_64 with glibc 2.28 and
Linux 5.15 virtual packages. Dependency resolution is a check of installability,
not evidence that the complete biological pipeline was executed on Linux.

Upload `fig1_ec2_ready.zip` to your EC2 machine, then:

```bash
unzip fig1_ec2_ready.zip
cd fig1_ec2
bash setup_linux.sh
.env/bin/python smoke_test.py
```

The installer downloads micromamba and creates `.env` inside this folder.
It requires `curl`, `tar`, `bzip2` and network access to the package repositories;
it uses no `sudo`. The environment downloads approximately **1.5 GB** of packages.
Allow additional disk space for extraction, package caching and genome analysis
outputs. Genomes are processed on the EC2 machine; the scripts do not upload them
to another analysis service or call AWS APIs.

The smoke test checks the two-genome plot and a clearly labeled synthetic
30-sample layout. It checks every plotted marker call against the corresponding
sample's input row. It does not annotate genomes or infer a phylogeny.

The installed main versions are Prokka 1.15.6, Panaroo 1.8.0, MAFFT 7.525,
IQ-TREE 2.4.0, mlst 2.35.0, Python 3.11.9, NumPy 1.26.4 and Biopython 1.85.
The exact package builds and checksums are in `ec2-linux-64.explicit.txt`.
`environment.yml` is provided for re-solving the environment if necessary;
the installer uses the explicit lock for reproducibility.

For a machine that already has Conda, an alternative installation is:

```bash
conda create --prefix "$PWD/.env" --file ec2-linux-64.explicit.txt
conda activate "$PWD/.env"
python run_pipeline.py --doctor
```

After a manual Conda installation, use `python run_pipeline.py ...` while the
environment is activated. The `run_ec2.sh` wrapper expects the micromamba installed
by `setup_linux.sh`.

## 2. Prepare your 30 genomes and tables

```text
data/
  genomes/
    NIHLM087.fna
    NIHLM095.fna
    ... remaining genomes ...
  gene_presence.tsv
  metadata.tsv
```

The genome folder should contain exactly the genomes for this run. The pipeline
reads `.fna`, `.fasta`, `.fa` and `.fas`, optionally gzip-compressed, directly in
that folder. It does not recursively scan unrelated subfolders.

Filenames define sample identifiers. The recognized FASTA suffix and optional
`.gz` suffix are removed, so `NIHLM087.fna` and `NIHLM087` match the same sample.
Identifiers must contain only letters, digits, dots, hyphens or underscores;
they must start with a letter or digit and be at most 64 characters long.
Two files or two table entries that normalize to the same key are rejected.

Use actual **tab-separated** tables, not space-aligned text.

`data/gene_presence.tsv` uses the same long format as your original table:

```tsv
sample	gene	present
NIHLM087.fna	mecA	1
NIHLM087.fna	IS256	0
NIHLM095.fna	mecA	1
NIHLM095.fna	IS256	0
```

Include one row per sample/marker combination for all 30 samples. `1` means
reported present; `0` means reported absent. Explicit `NA`, `?`, blank, or `.`
values mean unknown. Missing rows are rejected unless `--allow-missing` is set;
with that flag, missing rows appear as crosses rather than zeroes. Extra samples,
duplicate rows and nonbinary calls are rejected.

`data/metadata.tsv` can use either `sample` or `genome` as its ID column:

```tsv
genome	ST
NIHLM087	72
NIHLM095	72
```

Add the remaining 28 samples and their actual STs. `ST72`, `72` and `72.0` are
accepted and normalized to ST72. Use `NA` for an unknown ST. Optional columns
`group`, `outcome` and `biofilm` are displayed when supplied. The script preserves
their text values; it does not infer clinical outcomes or biofilm phenotypes.
Other columns can be selected with `--metadata-columns ST,group,custom_column`.

Metadata is optional. Without it, ST values display as NA. To infer missing STs
from the genomes, add `--infer-st`; for these S. epidermidis genomes the default
scheme is `sepidermidis`. Supplied STs are retained, and a conflicting known
inferred ST stops the run. Unknown/novel profiles stay unknown. The bundled mlst
software's allele database is a historical PubMLST snapshot, not a guarantee of
the newest allele definitions. Prefer your supplied ST table when it is complete.

The target gene table is independent of Panaroo's pangenome clustering table.
This package plots your reported marker calls; it does not independently re-call
AMR genes, virulence genes, IS256 or SCCmec from the genomes.

## 3. Validate and run

```bash
bash run_ec2.sh \
  --genomes-dir data/genomes \
  --genes data/gene_presence.tsv \
  --metadata data/metadata.tsv \
  --check
```

This checks the genomes and tables and prints assembly statistics. It does not
write analysis outputs or run the external bioinformatics tools. These basic
assembly checks do not confirm species identity or exclude contamination.
Use consistently prepared genomes of the same bacterial species for this
core-genome workflow; verify species and assembly quality before interpretation.

To see the planned analysis without executing it, replace `--check` with
`--dry-run`. To check the installed commands and record their versions, run
`bash run_ec2.sh --doctor`.

Run the full analysis:

```bash
bash run_ec2.sh \
  --genomes-dir data/genomes \
  --genes data/gene_presence.tsv \
  --metadata data/metadata.tsv \
  --gene-groups gene_groups.tsv \
  --threads 8 \
  --out results_30
```

Set `--threads` to suit your instance's available CPUs. Prokka runs one genome
at a time, using that thread count per genome, so annotation does not launch 30
independent jobs that compete for memory. Subsequent tools use the same thread
count. The pipeline does not impose a global memory limit.

For a long run, use an existing `tmux` session so disconnecting SSH does not
terminate it, or run the same command with `nohup` and redirect its console output
to a log. Individual tool logs are written inside the result folder.

Resume an interrupted run by repeating the same command with `--resume`:

```bash
bash run_ec2.sh \
  --genomes-dir data/genomes \
  --genes data/gene_presence.tsv \
  --metadata data/metadata.tsv \
  --gene-groups gene_groups.tsv \
  --threads 8 \
  --out results_30 \
  --resume
```

Successful stages are reused only when their input signatures and output hashes
match. Failed attempts and logs are preserved in `.failed-N` stage folders.
Changes to genome contents, tree parameters or tool versions require a fresh
output folder. Changing only gene calls, metadata, plot labels or rooting can
reuse the existing tree with `--resume`; the plot and its table join are rebuilt.

## 4. Outputs and tree interpretation

```text
results_30/
  figure1.pdf                 vector figure
  figure1.svg                 editable vector figure
  figure1.png                 preview
  figure1.matrix.tsv          marker matrix in displayed tree-tip order
  figure1.display_tree.nwk    tree with the chosen display rooting
  figure1.plot.json           sample/gene order and plot details
  tree.nwk                   normalized inference/input tree
  metadata_used.tsv          metadata used for plotting
  assembly_stats.tsv         contigs, length, N50, GC, ambiguity, checksums
  alignment_summary.json     core alignment dimensions and variation
  input_manifest.json        input paths and checksums
  software_versions.json     executed tool versions
  run_summary.json           method and run summary
  work/                      annotations, alignments, trees and tool logs
```

The pipeline requires at least four genomes for supported tree inference.
It preserves identical alignment samples with IQ-TREE's `-keep-ident`; it does
not interpret identical samples as resolved evolutionary relationships.
No-variation, empty, unequal-length or mismatched core alignments stop the run.

Defaults: bacterial Prokka annotation, genetic code 11, minimum annotated contig
length 200 bp; Panaroo strict mode, core threshold 0.95, MAFFT alignments;
Panaroo's filtered core alignment; IQ-TREE ModelFinder, 1,000 ultrafast bootstrap
replicates, 1,000 SH-aLRT replicates, seed 42. Full and filtered alignments are
retained in the Panaroo output. Use `--alignment unfiltered` in a fresh result
folder if you intentionally choose the unfiltered alignment. This workflow
does not include recombination masking; that choice should be considered for
the scientific interpretation of your bacterial population.

The initial tree is displayed in its supplied orientation. An arbitrary origin
in an unrooted tree does not establish ancestry. `--root midpoint` changes only
the display root and explicitly labels it as a midpoint choice. Alternatively,
`--outgroup SAMPLE_ID` uses one of your existing input samples as an outgroup;
it does not download or add reference genomes. Biological justification for an
outgroup remains necessary. Changing the root does not rebuild the inference.

The style configuration `gene_groups.tsv` assigns group labels, hex colors and
shapes (`square`, `circle`, `triangle`). Unlisted markers are retained in the
Other markers group. The generic `aacA` label is not equated with a specific
aminoglycoside resistance variant. SCCmec is a cassette; its binary call is
plotted exactly as provided. Presence of icaD/icaR does not establish a complete
icaADBC operon. Clinical or phenotypic information is not inferred from these
genetic markers.

## 5. Use an existing tree or replot without rebuilding

If you already have a Newick tree, supply it to the batch runner:

```bash
bash run_ec2.sh \
  --genomes-dir data/genomes \
  --genes data/gene_presence.tsv \
  --metadata data/metadata.tsv \
  --tree your_tree.nwk \
  --out results_existing_tree
```

This skips Prokka, Panaroo and IQ-TREE. Tree-tip IDs must match the genome/table
sample IDs after the documented filename normalization.

For plotting only, the genome files and bioinformatics tools are unnecessary:

```bash
.env/bin/python plot_fig1.py \
  --tree results_30/tree.nwk \
  --genes data/gene_presence.tsv \
  --metadata data/metadata.tsv \
  --gene-groups gene_groups.tsv \
  --branch-units substitutions/site \
  --show-support \
  --output revised_figure
```

For an external tree, supply the correct `--branch-units`; the standalone default
is "as supplied". Trees with absent or all-zero branch lengths are explicitly
drawn as cladograms. The standalone plotter can be installed in a small Python
environment using `pip install -r requirements-plot.txt`.

## 6. Validation supplied with this package

- The Linux x86_64 dependency environment was resolved successfully and locked
  with exact package builds/checksums; it was not installed/executed on EC2 here.
- Input validation was exercised on your original NIHLM087/NIHLM095 assemblies.
- The supplied-tree batch path was run on those two genomes.
- Two-tip and synthetic 30-tip plotting smoke tests checked every marker/sample
  join and the exported figures were visually reviewed.
- Unit tests check missing calls versus zeroes, filename collisions, duplicate
  calls, extra samples, ST normalization, invalid/no-variation alignments, stage
  integrity on resume, and recovery from failed attempts.
- External-tool orchestration was tested with explicit test doubles to check
  command arguments, pipeline outputs and stage reuse. Those doubles do not
  constitute a biological analysis or validation of Prokka/Panaroo/IQ-TREE.
- The actual 30-genome annotation/core-alignment/tree-inference run remains to be
  performed on your EC2 machine with your real dataset.

Run the unit checks after installation:

```bash
.env/bin/python -m unittest discover -s tests -v
```

## Sources

- Santos et al. (2025): https://doi.org/10.1038/s41598-025-09061-4
- Prokka: https://github.com/tseemann/prokka
- Panaroo quick start: https://github.com/gtonkinhill/panaroo/blob/master/docs/gettingstarted/quickstart.md
- Panaroo parameters: https://github.com/gtonkinhill/panaroo/blob/master/docs/gettingstarted/params.md
- IQ-TREE tutorial: https://iqtree.github.io/doc/Tutorial
- IQ-TREE options: https://iqtree.github.io/doc/Command-Reference
- mlst: https://github.com/tseemann/mlst
- BV-BRC tree service: https://www.bv-brc.org/docs/quick_references/services/phylogenetic_tree_building_service.html

Prepared 2026-10-09. Package version 1.0.0.

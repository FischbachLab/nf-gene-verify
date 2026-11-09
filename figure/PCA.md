# PCA of ortholog-cluster content, adapted from Conlan Figure 6

Added 2026-11-03. This analysis uses your completed Panaroo gene-cluster
matrix, not the 44 targeted markers plotted in Figure 1. It does not rerun
Prokka, Panaroo, MAFFT or IQ-TREE.

The reference is **Conlan et al. (2012), Figure 6**, in *Staphylococcus
epidermidis pan-genome sequence analysis reveals diversity of skin commensal
and hospital infection-associated isolates*. The Santos 2025 article used
earlier for Figure 1 has no Figure 6. Conlan used 35 genomes, OrthoMCL clusters
and R's `prcomp`. This package adapts the analysis to your 29-genome Panaroo
results; the coordinates and percentages are not expected to reproduce the
published numbers exactly.

## Run on the completed EC2 results

Upload the updated scripts and your `metadata2.tsv` to the existing package
folder, retaining its `.env` and `results_29/` directories. Then:

```bash
cd /mnt/efs/scratch/Xmeng/data/Allison/gene_verification_docker/fig1_ec2_29genomes
bash run_pca.sh --metadata metadata2.tsv
```

Outputs are saved beside your completed run, under `results_29_pca/`.
The original tree and Figure 1 files remain unchanged. No additional Python
dependencies are required: the calculation uses NumPy and Matplotlib already
present in the environment. R is not required at runtime.

The default matrix is:

```text
results_29/work/panaroo/pangenome/gene_presence_absence.Rtab
```

The metadata may use `sample`, `genome` or `strainID` as the ID column. It
must contain exactly the same samples as the matrix. Your supplied file has:

```tsv
strainID	ST	source
NIHLM095	72	commensal
NIH05001	2	nosocomial
```

All 29 rows are included in the supplied file. Metadata are joined by sample
ID, so their row order does not affect which source or ST is assigned to a point.
Sample filename suffixes are normalized as in the rest of the pipeline.
Missing source values remain "Not supplied"; clinical source is never inferred
from a sample name, ST, gene profile or PCA location.

## Method and your current results

1. Read the binary Panaroo matrix and transpose it to samples by gene clusters.
2. Exclude clusters present in every sample or absent in every sample because
   they have zero variance. By default, retain variable singletons as well as
   clusters shared by multiple samples.
3. Center each gene-cluster column across samples. Do not standardize its
   variance by default.
4. Compute PCA by singular value decomposition and label axes using the
   variance explained in your data. Axis signs are fixed using the largest
   absolute loading for reproducible drawing; their signs have no biological
   meaning.

Your matrix has **4,454 clusters**, of which **1,932** are present in all 29
genomes. The remaining **2,522 variable clusters**, including **720 singletons**,
enter the default PCA. PC1 explains **14.38%** and PC2 explains **10.34%**.

The centered, unscaled SVD calculation agrees with R 4.4.0 `stats::prcomp`
using `center=TRUE, scale.=FALSE` on the full input matrix. Removing constant
columns does not change its nonzero principal components. The paper does not
fully specify every preprocessing choice; these defaults are documented
choices, rather than a claim of an exact reconstruction of its analysis.

Panaroo Rtab encodes gene-cluster **presence/absence**, not RNA abundance,
expression, read coverage or gene copy number. Semicolon-separated IDs in
other Panaroo outputs can include merged annotations and are not treated as
copy-number measurements. Panaroo clustering also differs from OrthoMCL.

Optional sensitivity settings:

```bash
# Exclude singletons; retain clusters found in at least two genomes.
bash run_pca.sh --metadata metadata2.tsv --min-presence 2 \
  --output results_29_pca_no_singletons/figure6_pca

# Center and scale each variable feature to unit variance.
bash run_pca.sh --metadata metadata2.tsv --scale \
  --output results_29_pca_scaled/figure6_pca
```

These settings can change the PCA and must be reported if used. The default
is `--min-presence 1` without `--scale`.

## Figure styling

- Commensal: black open circles, 20 samples in the supplied metadata.
- Nosocomial, other ST: red open circles, 6 samples.
- Nosocomial, ST2: filled red circles, 3 samples.
- Unknown source: gray open circles, if present.

Labels are moved to avoid overlap and connected to their sample coordinates;
data points are not jittered. Closely related or identical gene-content
profiles can therefore overlap. The original `metadata2.tsv` has no A/B group
column, so a run with that file has no group outlines.

The 29 analyzed rows supply `PCA_group`: 22 genomes in A and 7 in B. The
current `data/metadata3.tsv` also contains two reference rows absent from the
completed 29-genome matrix. Use the matching 29-row subset to reproduce the
labeled A/B outlines:

```bash
bash run_pca.sh --metadata data/metadata_figure1_29.tsv --group-column PCA_group \
  --output results_29_pca_AB/figure6_pca_AB
```

This writes a separate figure and tables under `results_29_pca_AB/`, preserving
the earlier figure. Group names are printed beside their dashed ellipses, and
plot bounds expand to include complete outlines. Group assignments only affect
the overlay; the PCA coordinates and variance explained remain unchanged.

For other supplied groups, provide a column such as
`group` and request it explicitly:

```bash
bash run_pca.sh --metadata metadata_with_groups.tsv --group-column group
```

Dashed outlines are descriptive ellipses enclosing supplied group members,
not statistical confidence regions, hypothesis tests or inferred clusters.
At least three noncollinear points are required to draw an outline. PCA and
metadata overlays alone do not establish pathogenicity or transmission.

## Output files

- `figure6_pca.pdf`, `.svg`, `.png`: the PC1/PC2 plot.
- `.scores.tsv`: sample coordinates for all 28 nonzero components, with metadata.
- `.loadings.tsv`: gene-cluster loadings, centering values and scaling factors.
- `.variance.tsv`: eigenvalues, explained fractions and cumulative fractions.
- `.matrix.tsv`: the exact sample-by-feature matrix analyzed, before centering.
- `.features.tsv`: all input clusters, prevalence, inclusion and exclusion reasons.
- `.summary.json`: method, input checksums, feature counts and plotting categories.
- `validation_R_prcomp.json`: independent R comparison for the supplied dataset.

## Future genome sets and full-pipeline integration

For a completed 31-genome Panaroo run with matching source metadata:

```bash
bash run_pca.sh --results-dir results_31 --metadata data_31/metadata.tsv
```

This reads the 31-genome matrix and writes `results_31_pca/`. Do not append
new samples to the old Rtab manually; rerun Panaroo for the expanded dataset.

For a new full analysis, add these options to `run_ec2.sh` or `run_29.sh`:

```text
--ortholog-pca --pca-metadata metadata2.tsv
```

The full pipeline writes its PCA to `<out>/ortholog_pca/figure6_pca.*`.
Other integration options are `--pca-source-column`, `--pca-group-column`,
`--pca-min-presence`, `--pca-scale` and `--pca-no-labels`.
When `--tree` skips Panaroo, supply an existing matching `--pca-rtab` too.
The standalone command is the simplest way to add PCA to your completed run.

## Docker

Rebuild the updated image so it includes `ortholog_pca.py`. To run PCA without
annotation or tree reconstruction:

```bash
mkdir -p results_29_pca
docker run --rm --init --platform linux/amd64 \
  --user "$(id -u):$(id -g)" \
  --mount "type=bind,source=$PWD,target=/inputs,readonly" \
  --mount "type=bind,source=$PWD/results_29_pca,target=/outputs" \
  --entrypoint /usr/local/bin/_entrypoint.sh \
  fig1-29:1.0.0 python /opt/fig1/ortholog_pca.py \
  --rtab /inputs/results_29/work/panaroo/pangenome/gene_presence_absence.Rtab \
  --metadata /inputs/metadata2.tsv --output /outputs/figure6_pca
```

For PCA as part of a full Docker run, use
`bash run_29_docker.sh --threads 8 --ortholog-pca --pca-metadata /inputs/metadata2.tsv`.
Container builds still require verification on EC2; no image was built locally.

## Sources

- [Conlan et al. 2012, Figure 6 and statistical methods](https://link.springer.com/article/10.1186/gb-2012-13-7-r64).
- [Panaroo output definitions](https://github.com/gtonkinhill/panaroo/blob/master/docs/gettingstarted/output.md).
- [R `prcomp` documentation](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/prcomp.html).

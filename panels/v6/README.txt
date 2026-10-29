Frozen nf-gene-verify reference panel
  source : /sessions/eager-magical-brahmagupta/mnt/Downloads/results-29genomes_v6
  frozen : 2026-10-29T18:39:04Z
  genes  : 55

Use with:
  nextflow run main.nf -profile local --panel_dir <this directory> ...

Every gene with a <gene>.fasta here is built from that file and NOTHING is
fetched from NCBI. Genes NOT present fall through to the normal query, so
adding a gene to genes.txt without adding it here silently reintroduces drift
for that gene -- check the build log for 'FROZEN panel' on every gene you
expect to be frozen.

manifest.csv.provenance records which NCBI accessions these came from. It is
not read by the pipeline.

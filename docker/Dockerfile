FROM python:3.11-slim

RUN apt-get update \
    && apt-get install -y --no-install-recommends ncbi-blast+ \
    && rm -rf /var/lib/apt/lists/*

# NO PIPELINE CODE IS BAKED IN. This image is a runtime -- Python plus BLAST+ --
# and nothing else.
#
# It used to COPY six scripts into /opt/pipeline/. Nextflow prepends ./bin to
# PATH and stages the scripts it needs as hashed inputs, so those copies were
# shadowed for every step that mattered and were never what actually ran. The
# cost of keeping them was a second, drifting source of truth: when this image
# was last built, /opt/pipeline/build_protein_db.py still carried
# ALLOW_INCOMPLETE = {"embp","ebh","aap","bhp"} and --max-len-frac 2.0 -- the two
# guards that between them destroyed the bhp reference panel. The bug sat in the
# image for weeks while bin/ was correct, invisible because nothing ran it.
#
# Anything the pipeline executes now lives in nf-gene-verify/bin/ and is staged
# per task, so editing a script invalidates exactly the tasks that use it.
# Adding a COPY back here would reintroduce the drift.
#
# Consequence to be aware of: the standalone entry points that predate the
# Nextflow port -- run_pipeline.sh and run_no_docker.sh, i.e.
# `docker run ... bash /opt/pipeline/run_pipeline.sh` -- no longer work against
# this image. Both now live in nf-gene-verify/bin/ and can be mounted in if
# anyone still needs that path:
#   docker run --rm -v "$PWD/nf-gene-verify/bin:/opt/pipeline:ro" \
#     fischbachlab/gene-verify:latest bash /opt/pipeline/run_pipeline.sh ...
#
# BLAST+ is deliberately unpinned, and that is a known risk rather than a
# decision: `latest` currently carries 2.16.0+, while the EC2 host's own
# blastdbcmd is 2.5.0+, which is why staph_nr had to be rebuilt as a v5
# database. A rebuild today could land on a different minor version. Pin it if
# a rebuild has to reproduce a previous run exactly.

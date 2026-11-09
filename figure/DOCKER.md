# Run the 29-genome pipeline with Docker

This is an alternative to `setup_linux.sh` and `run_29.sh`. Docker installs
the locked bioinformatics dependencies inside an image; you do not need to
install or activate a host `.env` environment for this route.

It uses the same Prokka -> Panaroo/MAFFT -> IQ-TREE analysis and revised gene
groups as the existing package. It does not switch to Parsnp or whole-genome
alignment. The base image is micromamba 2.9.0 pinned by its manifest digest,
and the software packages use the existing Linux x86_64 explicit lock.

## 1. Build on EC2

Upload/extract the updated ZIP into a package folder containing the new
`Dockerfile`, `.dockerignore`, `run_29_docker.sh` and the other pipeline files.
Then run:

```bash
cd /mnt/efs/scratch/Xmeng/data/Allison/gene_verification_docker/fig1_ec2_29genomes
docker build --platform linux/amd64 -t fig1-29:1.0.0 .
```

The final dot is required: it selects this folder as the build context.
Docker must be installed and its daemon running. Use your account's usual
Docker access configuration; the wrapper calls `docker` directly.

The build downloads approximately 1.5 GB of Conda package artifacts, plus
the base image. Allow additional disk space for the extracted environment,
Docker build layers and analysis outputs. Prokka's bundled indexed databases
are included with its Bioconda package; they are not fetched at analysis time.
The build runs `prokka --listdb` and the pipeline's software-version checks,
and fails if those commands fail. These are installation checks, not a full
biological test run.

The image contains the pipeline and PCA Python modules and its software
environment. Genome FASTAs, gene tables, metadata and results are excluded
from the build context by the `.dockerignore` allowlist and are not baked
into the image. The image is built locally; these commands do not publish it.

## 2. Check and run

The wrapper already uses your EC2 genome folder:

```text
/mnt/efs/scratch/Xmeng/data/Allison/gene_verification_docker/nf-gene-verify/Segre_assemblies/
```

```bash
# Check genome FASTAs and sample IDs against the bundled tables.
bash run_29_docker.sh --check

# Run annotation, alignment, tree inference and plotting.
bash run_29_docker.sh --threads 8

# Reuse completed stages after an interruption.
bash run_29_docker.sh --threads 8 --resume
```

The wrapper mounts:

| Host location | Container location | Access |
| --- | --- | --- |
| Genome folder | `/genomes` | Read only |
| This package folder | `/inputs` | Read only |
| `results_29_docker/` inside this package | `/outputs` | Read and write |

It runs as your host numeric user/group IDs, which avoids creating root-owned
results and supports writing to EFS under your usual permissions. It also
uses Docker's `--init` to manage child processes. Results persist on EC2 when
the container exits, including `figure1.pdf`, `figure1.svg`, `figure1.png`,
`tree.nwk` and `work/` logs. The wrapper creates the output folder even when
using `--check`; the pipeline validation step itself writes no analysis files.

Keep the full output folder for resuming. Use the same image, genomes and tree
settings on a resumed run. The default Docker output folder is separate from
the host/micromamba output folder; do not mix results from different software
environments in one resumed run.

## 3. Change host paths

```bash
FIG1_DOCKER_GENOMES=/another/genome/folder \
FIG1_DOCKER_OUT=/another/result/folder \
  bash run_29_docker.sh --threads 8
```

`FIG1_DOCKER_IMAGE` selects an alternative locally built image tag.
Use `FIG1_DOCKER_OUT` rather than the pipeline's `--out` option with this
wrapper, so the selected folder is correctly mounted and results persist.
Other path options passed to the pipeline must use container paths. For
example, if you later put `actual_tree.nwk` inside this package:

```bash
FIG1_DOCKER_OUT="$PWD/results_29_docker_existing_tree" \
  bash run_29_docker.sh --tree /inputs/actual_tree.nwk
```

No tree file is required for the default run.

## 4. Inspect Prokka's databases or open a shell

The image's default entrypoint runs the Python pipeline. To run another
program, override it with micromamba's activation entrypoint:

```bash
docker run --rm --platform linux/amd64 \
  --entrypoint /usr/local/bin/_entrypoint.sh \
  fig1-29:1.0.0 prokka --listdb

docker run --rm -it --platform linux/amd64 \
  --entrypoint /usr/local/bin/_entrypoint.sh \
  fig1-29:1.0.0 bash
```

The installed environment and Prokka database directory are `/opt/conda` and
`/opt/conda/db`, respectively. No manual environment activation is required
when using the supplied entrypoints.

## Validation status

Docker is not running on the Mac where these files were prepared. Shell syntax,
Dockerfile lint, wrapper arguments/mounts and source packaging were checked
locally. Wrapper tests use an explicitly mocked Docker CLI; they do not execute
containers. The pinned base-image tag and its Linux amd64 manifest were
verified against Docker Hub. The complete Docker image build and biological
analysis still need to run on EC2. `docker_validation.json` records the checks.

References:

- [Official micromamba Docker quick start](https://micromamba-docker.readthedocs.io/en/latest/quick_start.html), including build-time activation and entrypoints.
- [Official micromamba advanced usage](https://micromamba-docker.readthedocs.io/en/latest/advanced_usage.html), including explicit locks and numeric user IDs.
- [Docker bind-mount documentation](https://docs.docker.com/engine/storage/bind-mounts/).

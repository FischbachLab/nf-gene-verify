FROM python:3.11-slim

RUN apt-get update \
    && apt-get install -y --no-install-recommends ncbi-blast+ \
    && rm -rf /var/lib/apt/lists/*

COPY scripts/build_protein_db.py scripts/call_genes.py scripts/run_pipeline.sh scripts/extract_loci.py scripts/run_nr_verification.py scripts/finalize_calls.py /opt/pipeline/
RUN chmod +x /opt/pipeline/run_pipeline.sh

#WORKDIR /work
#ENTRYPOINT ["bash", "/opt/pipeline/run_pipeline.sh"]

#!/usr/bin/env bash
# Run the two-round gene verification pipeline WITHOUT Docker.
# Requires: blastx/makeblastdb (e.g., via conda) on PATH, Python 3.8+.
#
# Usage:
#   ./run_no_docker.sh <genome.fasta> <outdir> [nr_db_path] [genes.txt]
#
#   nr_db_path  — optional; path to a pre-built nr BLAST DB (e.g.,
#                 /mnt/efs/databases/Blast/nr/db/nr). If omitted, round 2
#                 uses NCBI remote BLAST instead.
#   genes.txt   — optional; one gene name per line (defaults to the 8
#                 S. epidermidis genes).
set -euo pipefail

GENOME="$1"
OUTDIR="$2"
NR_DB="${3:-}"
GENES_FILE="${4:-}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# avoid BLAST memory-map errors on network filesystems (EFS/NFS)
export BLAST_DB_USE_MMAP=0

if [[ ! -f "$GENOME" ]]; then
    echo "ERROR: genome FASTA not found: $GENOME" >&2
    exit 1
fi
mkdir -p "$OUTDIR"

if [[ -n "$GENES_FILE" ]]; then
    mapfile -t GENES < <(grep -v '^\s*$' "$GENES_FILE" | tr -d '\r')
else
    GENES=(agrA agrB blaI blaZ fosB sarA sdrF sdrH)
fi
echo "Genes: ${GENES[*]}"

# ---- Round 1: coordinate identification ---------------------------------
python3 "$SCRIPT_DIR/build_protein_db.py" --genes "${GENES[@]}" \
    --out "$OUTDIR/protein_db"

blastx -query "$GENOME" \
       -subject "$OUTDIR/protein_db/combined.fasta" \
       -evalue 1e-5 \
       -outfmt "6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore" \
       -out "$OUTDIR/blastx_raw.tsv"
echo "round 1: $(wc -l < "$OUTDIR/blastx_raw.tsv") hits"

python3 "$SCRIPT_DIR/call_genes.py" \
    --blastx "$OUTDIR/blastx_raw.tsv" \
    --manifest "$OUTDIR/protein_db/manifest.csv" \
    --genes "${GENES[@]}" \
    --out "$OUTDIR/gene_calls.tsv"

# ---- Round 2: nr verification of all genes ------------------------------
python3 "$SCRIPT_DIR/extract_loci.py" \
    --calls "$OUTDIR/gene_calls.tsv" --genome "$GENOME" \
    --flank 300 --out "$OUTDIR/nr_verification" --all

if [[ -n "$NR_DB" ]]; then
    echo "round 2: LOCAL blastx vs $NR_DB"
    # threads: BLAST_THREADS if set, else the script auto-detects all cores
    THREAD_ARGS=()
    [[ -n "${BLAST_THREADS:-}" ]] && THREAD_ARGS=(--threads "$BLAST_THREADS")
    python3 "$SCRIPT_DIR/run_nr_verification.py" \
        --loci-dir "$OUTDIR/nr_verification" --out "$OUTDIR/nr_verification" \
        --db "$NR_DB" "${THREAD_ARGS[@]}"
else
    echo "round 2: NCBI remote BLAST (queue-dependent)"
    python3 "$SCRIPT_DIR/run_nr_verification.py" \
        --loci-dir "$OUTDIR/nr_verification" --out "$OUTDIR/nr_verification" \
        --hitlist 100
fi

# ---- Finalize ------------------------------------------------------------
python3 "$SCRIPT_DIR/finalize_calls.py" \
    --calls "$OUTDIR/gene_calls.tsv" \
    --nr "$OUTDIR/nr_verification/nr_top_hits.tsv" \
    --loci-dir "$OUTDIR/nr_verification" \
    --out "$OUTDIR/gene_presence_table_v2.csv"

echo "Done:"
echo "  $OUTDIR/gene_calls.tsv               (round 1)"
echo "  $OUTDIR/nr_verification/             (round 2 evidence)"
echo "  $OUTDIR/gene_presence_table_v2.csv   (final table)"

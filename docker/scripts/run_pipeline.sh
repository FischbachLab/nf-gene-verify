#!/usr/bin/env bash
# End-to-end annotation-free gene verification:
#   NCBI protein DB build -> blastx screen -> presence/absence calling
#
# Usage (inside the container):
#   run_pipeline.sh <genome.fasta> <outdir> [genes.txt] [--verify-nr]
#
# genes.txt is optional (one gene name per line); defaults to the 8
# S. epidermidis genes agrA agrB blaI blaZ fosB sarA sdrF sdrH.
# --verify-nr additionally extracts divergent/partial loci and verifies them
# against NCBI nr via remote BLAST (queue-dependent, minutes to hours).
set -euo pipefail

# --- environment sanitation (E2BIG defense) -------------------------------
# "Argument list too long" (E2BIG) at exec means the environment block is
# over the kernel limit (total argv+envp, or a single var > 128 KB).
# Sanitize IN-PROCESS (pure bash, no exec) before running anything external.
PATH="/usr/local/bin:/usr/bin:/bin"   # container layout; also caps a bloated PATH
for varname in $(compgen -e); do
    case "$varname" in PATH|HOME) continue ;; esac
    val="${!varname:-}"
    if (( ${#val} > 32768 )); then
        echo "[WARN] dropping oversized env var: $varname (${#val} bytes)"
        unset "$varname"
    fi
done
unset varname val

# Workaround for "Argument list too long" (E2BIG): run all external commands
# with a minimal environment as a second layer of defense.
run_clean() {
    env -i \
        PATH="${PATH}" \
        HOME="${HOME:-/root}" \
        TMPDIR="${TMPDIR:-/tmp}" \
        NCBI_API_KEY="${NCBI_API_KEY:-}" \
        NR_DB="${NR_DB:-}" \
        BLAST_THREADS="${BLAST_THREADS:-4}" \
        BLAST_DB_USE_MMAP=0 \
        "$@"
}

GENOME="$1"
shift
OUTDIR="$1"
shift
GENES_FILE=""
VERIFY_NR=0
for arg in "$@"; do
    case "$arg" in
        --verify-nr) VERIFY_NR=1 ;;
        -*) echo "[WARN] ignoring unknown flag: $arg" >&2 ;;
        *) GENES_FILE="$arg" ;;
    esac
done

if [[ ! -f "$GENOME" ]]; then
    echo "ERROR: genome FASTA not found: $GENOME" >&2
    exit 1
fi
mkdir -p "$OUTDIR"

if [[ -n "$GENES_FILE" ]]; then
    if [[ ! -f "$GENES_FILE" ]]; then
        echo "ERROR: genes file not found: $GENES_FILE" >&2
        echo "       (was the genome path passed twice, or file contents" >&2
        echo "       expanded into the arguments?)" >&2
        exit 1
    fi
    # guard: a gene list is short names, one per line — not a FASTA/sequence file
    if head -c 1 "$GENES_FILE" | grep -q '>' \
       || awk 'length($0) > 100 {found=1} END{exit !found}' "$GENES_FILE"; then
        echo "ERROR: '$GENES_FILE' does not look like a gene-name list" >&2
        echo "       (one gene name per line, e.g., 'agrA'). It looks like a" >&2
        echo "       FASTA/sequence file — check the arguments." >&2
        exit 1
    fi
    mapfile -t GENES < <(grep -v '^\s*$' "$GENES_FILE" | tr -d '\r')
else
    GENES=(agrA agrB blaI blaZ fosB sarA sdrF sdrH)
fi
echo "Genes to screen: ${GENES[*]}"

# Step 1: build the staph-restricted protein reference set from NCBI
# (set NCBI_API_KEY env var to raise the rate limit from 3 to 10 req/s)
run_clean python3 /opt/pipeline/build_protein_db.py \
    --genes "${GENES[@]}" --out "$OUTDIR/protein_db"

# Step 2: blastx screen — genome translated in 6 frames vs protein set.
# NOTE: do NOT use -max_target_seqs; it silently truncates results and can
# drop entire genes (observed with subject-mode searches).
run_clean blastx -query "$GENOME" \
       -subject "$OUTDIR/protein_db/combined.fasta" \
       -evalue 1e-5 \
       -outfmt "6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore" \
       -out "$OUTDIR/blastx_raw.tsv"
echo "blastx hits: $(wc -l < "$OUTDIR/blastx_raw.tsv")"

# Step 3: presence/absence calling
run_clean python3 /opt/pipeline/call_genes.py \
    --blastx "$OUTDIR/blastx_raw.tsv" \
    --manifest "$OUTDIR/protein_db/manifest.csv" \
    --genes "${GENES[@]}" \
    --out "$OUTDIR/gene_calls.tsv"

echo "Done. Results:"
echo "  $OUTDIR/gene_calls.tsv          (raw calls)"
echo "  $OUTDIR/blastx_raw.tsv          (raw BLAST evidence)"
echo "  $OUTDIR/protein_db/manifest.csv (reference DB provenance)"

# Optional round 2: verify ALL genes against nr (local DB or NCBI remote)
if [[ "$VERIFY_NR" -eq 1 ]]; then
    echo "Round 2: extracting ALL gene loci for nr verification..."
    run_clean python3 /opt/pipeline/extract_loci.py \
        --calls "$OUTDIR/gene_calls.tsv" --genome "$GENOME" \
        --flank 300 --out "$OUTDIR/nr_verification" --all
    if [[ -n "${NR_DB:-}" ]]; then
        # Local mode: NR_DB = path to a pre-built nr BLAST DB
        # (e.g., mount the EFS volume: -v /mnt/efs/databases/Blast/nr/db:/db
        #  and set NR_DB=/db/nr)
        echo "Running LOCAL blastx vs $NR_DB ..."
        # threads: BLAST_THREADS if set, else the script auto-detects all cores
        THREAD_ARGS=()
        [[ -n "${BLAST_THREADS:-}" ]] && THREAD_ARGS=(--threads "$BLAST_THREADS")
        run_clean python3 /opt/pipeline/run_nr_verification.py \
            --loci-dir "$OUTDIR/nr_verification" --out "$OUTDIR/nr_verification" \
            --db "$NR_DB" "${THREAD_ARGS[@]}"
    else
        # Remote mode: NCBI BLAST URL API (queue-dependent, be patient)
        echo "Submitting loci to NCBI remote BLAST (queue-dependent, be patient)..."
        run_clean python3 /opt/pipeline/run_nr_verification.py \
            --loci-dir "$OUTDIR/nr_verification" --out "$OUTDIR/nr_verification" \
            --hitlist 100
    fi
    echo "Merging round-1 calls with nr evidence..."
    run_clean python3 /opt/pipeline/finalize_calls.py \
        --calls "$OUTDIR/gene_calls.tsv" \
        --nr "$OUTDIR/nr_verification/nr_top_hits.tsv" \
        --loci-dir "$OUTDIR/nr_verification" \
        --out "$OUTDIR/gene_presence_table_v2.csv"
    echo "Final table: $OUTDIR/gene_presence_table_v2.csv"
fi

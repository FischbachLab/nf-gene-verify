#!/usr/bin/env bash
# Re-run ONE gene across a whole seedfile, then merge it into an existing
# results directory. Written for the bhp fix, but gene-agnostic.
#
#   tools/rerun_gene.sh bhp <seedfile.csv> <existing-results-dir> [outdir]
#
# Override the nr database with NR_DB=... if it lives elsewhere:
#   NR_DB=/path/to/staph_nr tools/rerun_gene.sh bhp ...
#
# Why not re-run the whole pipeline: the reference panel is a hashed input to
# BUILD_PROTEIN_DB, so changing it invalidates every downstream task and
# -resume buys nothing -- a full re-run is ~13.5 h on the 8-core node. One gene
# is a 5-sequence blastx subject and ~29 nr loci, so it is minutes.
#
# The merge is done by tools/merge_gene.py, which rewrites only that gene's
# rows and regenerates the summaries. It refuses to touch anything else.

set -euo pipefail

GENE=${1:?usage: rerun_gene.sh <gene> <seedfile.csv> <results-dir> [outdir]}
SEEDFILE=${2:?missing seedfile}
RESULTS=${3:?missing existing results dir}
OUTDIR=${4:-"$(dirname "$RESULTS")/$(basename "$RESULTS")-${GENE}-rerun"}

HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
cd "$HERE"

[[ -f "$SEEDFILE" ]] || { echo "no such seedfile: $SEEDFILE" >&2; exit 1; }
[[ -d "$RESULTS"  ]] || { echo "no such results dir: $RESULTS" >&2; exit 1; }

# ── nr database ────────────────────────────────────────────────────────────
# nextflow.config ships nr_db = '' and the pipeline reads EMPTY AS "use NCBI's
# remote BLAST service": the workflow takes the NR_VERIFICATION branch, submits
# an RID to NCBI's public queue and polls it every 10 s (the only artifact in
# the task dir is rids.json). That is days for 29 loci, and it is the WRONG
# DATABASE besides -- the rest of the table was verified against the local
# Staphylococcus subset, so remote full-nr rows cannot be spliced in beside it.
# Pass it explicitly and refuse to launch without it.
NR_DB=${NR_DB:-/mnt/efs/databases/Blast/blastdb_custom/staph_nr/staph_nr}
if ! compgen -G "${NR_DB}.p*" > /dev/null; then
  cat >&2 <<EOF
no BLAST database found at: $NR_DB
  (looked for ${NR_DB}.p*)
Set NR_DB=<prefix> to the staph_nr prefix used for the original run.
Refusing to fall back to remote NCBI BLAST -- it would take days and would
verify against a different database than every other row in the target table.
EOF
  exit 1
fi

GENELIST=$(mktemp -t "${GENE}_only.XXXXXX.txt")
echo "$GENE" > "$GENELIST"
trap 'rm -f "$GENELIST"' EXIT

echo "=== re-running '$GENE' only"
echo "    seedfile : $SEEDFILE"
echo "    outdir   : $OUTDIR"
echo "    nr_db    : $NR_DB"
echo

# --nr_verify_present true MATTERS HERE. With it false, a round-1 'present'
# call becomes final 'present' with no nr evidence at all -- and a round-1
# present call against a wrong panel is exactly the failure being fixed. The
# v4 run set it true (all 848 present calls carry a confirming accession); a
# re-run that did not would merge rows weaker than the ones they replace.
nextflow run main.nf -profile local \
  --seedfile "$SEEDFILE" \
  --genes    "$GENELIST" \
  --outdir   "$OUTDIR" \
  --verify_nr true \
  --nr_db    "$NR_DB" \
  --nr_verify_present true \
  -work-dir "work-${GENE}" \
  -ansi-log false

# ── panel sanity gate ──────────────────────────────────────────────────────
# The whole point of this re-run is that the panel was wrong. Check the panel
# BEFORE trusting the calls, and before merging anything.
python3 - "$GENE" "$OUTDIR" <<'PY'
import csv, statistics, sys, os
gene, outdir = sys.argv[1].lower(), sys.argv[2]
mf = os.path.join(outdir, "protein_db", "manifest.csv")
rows = [r for r in csv.DictReader(open(mf)) if r["gene"].lower() == gene]
if not rows:
    sys.exit(f"FAIL: no {gene} references in {mf} -- panel is empty")
lens = sorted(int(r["length"] or 0) for r in rows)
print(f"panel: n={len(rows)}  range={lens[0]}-{lens[-1]}  median={statistics.median(lens):.0f}")
for r in sorted(rows, key=lambda r: -int(r["length"] or 0)):
    print(f"  {r['length']:>6}  {r['accession']:<12} {r['organism'][:34]:<34} {r['product'][:50]}")
if gene == "bhp":
    # Bhp is ~2,500 aa. The corrupted v4 panel was 160-643 aa, every record a
    # short ', partial' fragment. Anything in that range means the fix did not
    # take -- check ALLOW_INCOMPLETE is empty and --max-len-frac is 0.
    if lens[0] < 1000 or lens[-1] < 2000:
        sys.exit(f"FAIL: bhp panel is {lens[0]}-{lens[-1]} aa; expected ~1542-3733 "
                 f"(5 refs). Short fragments are back -- the guards are not off.")
print("panel gate: OK")
PY

echo
echo "=== calls for '$GENE' (round 1 -> final, per genome)"
# csv, not awk -F',' -- the nr_verification column is a quoted free-text
# sentence that contains commas, and splitting on them shifts every field
# after it.
python3 - "$GENE" "$OUTDIR" <<'PY'
import csv, glob, os, sys
gene, outdir = sys.argv[1], sys.argv[2]
for d in sorted(glob.glob(os.path.join(outdir, "*", "gene_presence_table_v2.csv"))):
    s = os.path.basename(os.path.dirname(d))
    for r in csv.DictReader(open(d)):
        if r["gene"] != gene:
            continue
        print(f"  {s:<12} r1={r['round1_call']:<9} "
              f"{r['round1_identity_pct'] or '-':>6} id {r['round1_coverage_pct'] or '-':>6} cov"
              f"   nr={r['nr_identity_pct'] or '-':>6} cov={r['nr_hit_coverage_pct'] or '-':>6}"
              f" rank={r['nr_hit_rank'] or '-':<4} -> {r['final_call']}")
PY

cat <<EOF

=== next step (not done automatically -- inspect the calls above first)

    python3 tools/merge_gene.py --gene $GENE \\
        --from "$OUTDIR" --into "$RESULTS"

EOF

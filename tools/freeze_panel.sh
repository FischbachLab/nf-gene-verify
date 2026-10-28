#!/usr/bin/env bash
# Freeze a run's reference panel so later runs reproduce it exactly.
#
#   tools/freeze_panel.sh <results-dir> <panel-dir>
#   nextflow run main.nf ... --panel_dir <panel-dir>
#
# Why. The per-gene NCBI query is genus-wide and its result set drifts, so
# rebuilding RESAMPLES the panel rather than reproducing it: between two runs
# three weeks apart, 16 of 55 panels drew different references (fdh swapped 7
# of 20, aap 5, kdpB 4) and 13 calls moved for that reason alone, with no code
# change behind any of them. A re-run was therefore not a way to confirm a
# result -- it produced a new one needing its own audit. With the panel frozen,
# a re-run reproduces the table and any difference is attributable to the code.
#
# What gets copied is just <results>/protein_db/*.fasta: those files already
# carry the 'acc|gene|organism|product' deflines build_protein_db.py parses, so
# nothing is converted. combined.fasta is skipped -- it is a concatenation the
# build regenerates, and leaving it in would register a bogus gene named
# 'combined'.

set -euo pipefail

SRC=${1:?usage: freeze_panel.sh <results-dir> <panel-dir>}
DST=${2:?missing destination panel dir}

PDB="$SRC/protein_db"
[[ -d "$PDB" ]] || { echo "no protein_db/ under $SRC" >&2; exit 1; }

mkdir -p "$DST"
n=0
for f in "$PDB"/*.fasta; do
    b=$(basename "$f")
    [[ "$b" == "combined.fasta" ]] && continue
    cp -f "$f" "$DST/$b"
    n=$((n + 1))
done
[[ -f "$PDB/manifest.csv" ]] && cp -f "$PDB/manifest.csv" "$DST/manifest.csv.provenance"

cat > "$DST/README.txt" <<EOF
Frozen nf-gene-verify reference panel
  source : $SRC
  frozen : $(date -u +%Y-%m-%dT%H:%M:%SZ)
  genes  : $n

Use with:
  nextflow run main.nf -profile local --panel_dir <this directory> ...

Every gene with a <gene>.fasta here is built from that file and NOTHING is
fetched from NCBI. Genes NOT present fall through to the normal query, so
adding a gene to genes.txt without adding it here silently reintroduces drift
for that gene -- check the build log for 'FROZEN panel' on every gene you
expect to be frozen.

manifest.csv.provenance records which NCBI accessions these came from. It is
not read by the pipeline.
EOF

echo "froze $n gene panel(s) -> $DST"

# ── verify the frozen panel reloads identically ────────────────────────────
# A frozen panel that silently loses references is worse than no freeze, so
# check accession-for-accession against the manifest it came from.
python3 - "$PDB/manifest.csv" "$DST" <<'PY'
import csv, os, sys
from collections import defaultdict
manifest, dst = sys.argv[1], sys.argv[2]
want = defaultdict(set)
if os.path.exists(manifest):
    for r in csv.DictReader(open(manifest)):
        want[r["gene"]].add(r["accession"])
bad = 0
for gene, accs in sorted(want.items()):
    p = os.path.join(dst, f"{gene}.fasta")
    if not os.path.exists(p):
        print(f"  MISSING {gene}.fasta"); bad += 1; continue
    got = {l[1:].split("|")[0] for l in open(p) if l.startswith(">")}
    if got != accs:
        print(f"  {gene}: {len(accs)} in manifest, {len(got)} in fasta "
              f"(missing {sorted(accs - got)[:3]}, extra {sorted(got - accs)[:3]})")
        bad += 1
print(f"verified {len(want)} gene(s); {bad} mismatch(es)")
sys.exit(1 if bad else 0)
PY

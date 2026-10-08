#!/usr/bin/env python3
"""Call gene presence/absence from blastx hits (genome vs protein reference set).

Hits are clustered by genomic distance (gap <= 5 kb) so distant weak
cross-homology hits do not merge into fake loci. Within a cluster, coverage
is computed per reference accession (merged union of aligned reference
intervals); identity is the length-weighted mean.

Calling thresholds (standard): present = >=90% identity over >=80% reference
coverage; divergent = full coverage, <90% identity; partial = <80% coverage;
absent = no qualifying hit. Genes with no reference proteins in the manifest
(build_protein_db.py found none) are reported as 'no_reference' — they were
never searched, so 'absent' would be a false negative.
"""
import argparse
import csv
import re
from collections import defaultdict

ID_MIN, COV_MIN, GAP = 90.0, 80.0, 5000


def parse_hits(raw_path):
    hits = []
    with open(raw_path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 12:
                continue
            parts = f[1].split("|")
            hits.append({
                "contig": f[0], "acc": parts[0], "gene": parts[1],
                "pident": float(f[2]), "alen": int(f[3]),
                "qstart": int(f[6]), "qend": int(f[7]),
                "sstart": int(f[8]), "send": int(f[9]),
                "evalue": float(f[10]), "bits": float(f[11]),
            })
    return hits


def merge_intervals(ivs):
    ivs = sorted(ivs)
    merged = []
    for s, e in ivs:
        if merged and s <= merged[-1][1] + 1:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])
    return merged


def cluster_by_genome(hsps, gap=GAP):
    """Cluster HSPs on the same contig whose genomic spans are within `gap` bp."""
    hsps = sorted(hsps, key=lambda h: min(h["qstart"], h["qend"]))
    clusters = []
    for h in hsps:
        lo, hi = min(h["qstart"], h["qend"]), max(h["qstart"], h["qend"])
        if clusters and h["contig"] == clusters[-1]["contig"] and lo <= clusters[-1]["hi"] + gap:
            clusters[-1]["hi"] = max(clusters[-1]["hi"], hi)
            clusters[-1]["hsps"].append(h)
        else:
            clusters.append({"contig": h["contig"], "lo": lo, "hi": hi, "hsps": [h]})
    return clusters


def score_cluster(cluster, manifest):
    """Best (coverage, identity) reference accession within this cluster."""
    best = None
    by_acc = defaultdict(list)
    for h in cluster["hsps"]:
        by_acc[h["acc"]].append(h)
    for acc, hsps in by_acc.items():
        slen = int(manifest[acc]["length"])
        ivs = [[min(h["sstart"], h["send"]), max(h["sstart"], h["send"])] for h in hsps]
        covered = sum(e - s + 1 for s, e in merge_intervals(ivs))
        cov = 100.0 * covered / slen
        ident = sum(h["pident"] * h["alen"] for h in hsps) / sum(h["alen"] for h in hsps)
        if best is None or (cov, ident) > best[:2]:
            best = (cov, ident, acc)
    return best


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--blastx", required=True,
                    help="blastx outfmt 6 TSV (qseqid sseqid pident length ... bitscore)")
    ap.add_argument("--manifest", required=True, help="manifest.csv from build_protein_db.py")
    ap.add_argument("--genes", nargs="+", required=True, help="Gene names to call")
    ap.add_argument("--out", required=True, help="Output TSV path")
    ap.add_argument("--identity-min", type=float, default=ID_MIN)
    ap.add_argument("--coverage-min", type=float, default=COV_MIN)
    ap.add_argument("--gap", type=int, default=GAP, help="Genomic clustering gap (bp)")
    args = ap.parse_args()

    manifest = {r["accession"]: r for r in csv.DictReader(open(args.manifest))}
    ref_genes = {r["gene"] for r in manifest.values()}
    hits = parse_hits(args.blastx)
    by_gene = defaultdict(list)
    for h in hits:
        by_gene[h["gene"]].append(h)

    rows = []
    for gene in args.genes:
        if gene not in ref_genes:
            # No reference proteins were built for this gene, so blastx could
            # never find it. Report that explicitly instead of 'absent'.
            print(f"{gene:8s} NO REFERENCE PROTEINS — not searched. Fix with "
                  f"build_protein_db.py --alias/--extra-fasta (see README).")
            rows.append({"gene": gene, "call": "no_reference", "identity_pct": "",
                         "coverage_pct": "", "contig": "", "coords": "",
                         "best_accession": "", "organism": "", "n_hits": 0})
            continue

        ghits = by_gene.get(gene, [])
        if not ghits:
            rows.append({"gene": gene, "call": "absent", "identity_pct": "",
                         "coverage_pct": "", "contig": "", "coords": "",
                         "best_accession": "", "organism": "", "n_hits": 0})
            continue

        best = None
        for cluster in cluster_by_genome(ghits, args.gap):
            cov, ident, acc = score_cluster(cluster, manifest)
            if best is None or (cov, ident) > (best["cov"], best["ident"]):
                best = {"cov": cov, "ident": ident, "acc": acc,
                        "contig": cluster["contig"],
                        "qmin": min(min(h["qstart"], h["qend"]) for h in cluster["hsps"]),
                        "qmax": max(max(h["qstart"], h["qend"]) for h in cluster["hsps"])}

        ident, cov = best["ident"], best["cov"]
        if ident >= args.identity_min and cov >= args.coverage_min:
            call = "present"
        elif cov >= args.coverage_min:
            call = "divergent"
        else:
            call = "partial"

        org_m = re.search(r"\[([^\]]+)\]\s*$", manifest[best["acc"]]["product"])
        rows.append({
            "gene": gene, "call": call,
            "identity_pct": f"{ident:.1f}", "coverage_pct": f"{cov:.1f}",
            "contig": best["contig"],
            "coords": f"{best['qmin']}-{best['qmax']}",
            "best_accession": best["acc"],
            "organism": org_m.group(1) if org_m else "",
            "n_hits": len(ghits),
        })

    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    for r in rows:
        print(f"{r['gene']:8s} {r['call']:10s} id={r['identity_pct']:>6s}% "
              f"cov={r['coverage_pct']:>6s}% contig={r['contig']} "
              f"coords={r['coords']:>18s} ref={r['best_accession']} [{r['organism']}]")


if __name__ == "__main__":
    main()

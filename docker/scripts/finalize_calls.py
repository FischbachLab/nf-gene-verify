#!/usr/bin/env python3
"""Merge round-1 gene calls with round-2 nr verification into the final table.

Reads gene_calls.tsv (round 1: coordinates from the staph protein screen),
nr_top_hits.tsv (round 2: nr evidence per locus), and the locus FASTA
headers (which encode the extracted window as contig:start-end) to map nr
hits back to genomic coordinates.

Revision rule (mechanical, documented):
  - Take the top nr hit (by bitscore) whose query span overlaps the round-1
    called region for the gene.
  - If that hit's title contains the gene name (case-insensitive word match,
    plus a small synonym map, e.g., fosB <-> YfcC), the gene is CONFIRMED
    present — regardless of round-1 identity/coverage. This neutralizes
    artifacts from RefSeq-restricted round-1 reference sets.
  - If no nr hit overlaps the called region, or the top overlapping hit does
    not match the gene name, the gene is FLAGGED with the evidence shown;
    nothing is silently dropped.

Usage:
  finalize_calls.py --calls gene_calls.tsv --nr nr_top_hits.tsv \
                    --loci-dir nr_verification/ --out gene_presence_table_v2.csv
"""
import argparse
import csv
import os
import re

# recognized aliases: gene name -> additional title tokens that confirm it
SYNONYMS = {
    "fosb": ["fosb", "yfcc"],
    "fdh": ["fdh", "formate dehydrogenase"],
    "bhp": ["bhp", "biofilm-associated protein", "biofilm associated protein",
            "cell wall associated biofilm protein"],
}


def load_locus_windows(loci_dir):
    """Parse locus FASTA headers: >{gene}_locus_{contig}_{start}-{end}."""
    windows = {}
    for f in sorted(os.listdir(loci_dir)):
        if not f.endswith("_locus.fasta"):
            continue
        gene = f.replace("_locus.fasta", "")
        header = open(os.path.join(loci_dir, f)).readline().lstrip(">").strip()
        m = re.match(rf"{gene}_locus_(.+?)_(\d+)-(\d+)$", header)
        if not m:
            print(f"[WARN] unrecognized locus header: {header}")
            continue
        windows[gene] = {"contig": m.group(1),
                         "win_start": int(m.group(2)),
                         "win_end": int(m.group(3))}
    return windows


def title_matches(gene, title):
    tokens = SYNONYMS.get(gene.lower(), [gene.lower()])
    t = title.lower()
    return any(re.search(rf"\b{re.escape(tok)}\b", t) for tok in tokens)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--calls", required=True, help="gene_calls.tsv (round 1)")
    ap.add_argument("--nr", required=True, help="nr_top_hits.tsv (round 2)")
    ap.add_argument("--loci-dir", required=True,
                    help="dir with the *_locus.fasta files (for window coords)")
    ap.add_argument("--out", required=True,
                    help="final table path, e.g., gene_presence_table_v2.csv")
    args = ap.parse_args()

    calls = {r["gene"]: r for r in csv.DictReader(open(args.calls), delimiter="\t")}
    windows = load_locus_windows(args.loci_dir)

    # group nr hits by gene
    nr_by_gene = {}
    if os.path.exists(args.nr):
        for r in csv.DictReader(open(args.nr), delimiter="\t"):
            nr_by_gene.setdefault(r["gene"], []).append(r)
    else:
        print(f"[WARN] nr hits file not found: {args.nr}")

    rows = []
    for gene, c in calls.items():
        win = windows.get(gene)
        hits = sorted(nr_by_gene.get(gene, []),
                      key=lambda r: -float(r["bit_score"]))
        coords = c["coords"]
        gs, ge = (int(x) for x in coords.split("-")) if coords else (None, None)

        top = None
        if win and gs is not None:
            for h in hits:
                # genomic span of this hit's query alignment
                g1 = win["win_start"] + int(h["query_start"]) - 1
                g2 = win["win_start"] + int(h["query_end"]) - 1
                if g2 >= gs and g1 <= ge:
                    top = dict(h, g1=g1, g2=g2)
                    break  # hits are bitscore-sorted

        if win is None or gs is None:
            # No locus window / round-1 coordinates (absent or no_reference
            # genes): nr verification cannot apply, so keep the round-1 call
            # rather than flagging it.
            final = c["call"]
            nr_status = ("nr verification not applicable (no round-1 "
                         "coordinates); round-1 call kept")
        elif top is None:
            nr_status = "flagged: no nr hit overlaps the called region"
            final = "flagged"
        elif title_matches(gene, top["title"]):
            cov = top.get("hit_cov_pct") or top.get("query_cov_pct") or ""
            nr_status = f"confirmed by nr: {top['accession']} " \
                        f"({top['identity_pct']}% id, {cov}% cov)"
            final = "present"
        else:
            nr_status = f"flagged: top nr hit is {top['accession']} " \
                        f"({top['title'][:60]}), not {gene}"
            final = "flagged"

        cov_col = (top.get("hit_cov_pct") or top.get("query_cov_pct") or ""
                   if top else "")
        rows.append({
            "gene": gene,
            "round1_call": c["call"],
            "round1_identity_pct": c["identity_pct"],
            "round1_coverage_pct": c["coverage_pct"],
            "contig": c["contig"],
            "coordinates": coords,
            "nr_top_hit_accession": top["accession"] if top else "",
            "nr_top_hit_title": top["title"][:80] if top else "",
            "nr_organism": top["organism"] if top else "",
            "nr_identity_pct": top["identity_pct"] if top else "",
            "nr_hit_coverage_pct": cov_col,
            "nr_verification": nr_status,
            "final_call": final,
        })

    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    for r in rows:
        print(f"{r['gene']:8s} round1={r['round1_call']:10s} "
              f"final={r['final_call']:8s} {r['nr_verification']}")
    print(f"\nwritten: {args.out}")


if __name__ == "__main__":
    main()

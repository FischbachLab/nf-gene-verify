#!/usr/bin/env python3
"""Derive the wide presence matrix and two summaries from the long matrix.

Inputs the concatenated 'sample<TAB>gene<TAB>present' table and writes:

  gene_presence_matrix_wide.tsv   samples x genes, 0/1, nothing else
  gene_prevalence.tsv             gene, n_genomes, pct     (descending)
  genes_per_genome.tsv            sample, n_genes          (descending)

Three deliberate choices:

1. GENE ORDER follows genes.txt, not the alphabet, so operons stay grouped
   (arc, opp3, kdp; then resistance; then adhesins). Reading down a column
   block is how the ica operon's internal consistency becomes visible --
   icaA/B/C/D/R all at 11/29 is a check that the panels agree; 11/11/2/11/11
   is what a bug looks like.

2. The WIDE MATRIX CARRIES NO TOTALS. Row and column sums go in the two
   companion files instead, so the matrix loads straight into R or pandas
   without a totals row masquerading as a sample. Nothing to strip.

3. MISSING CELLS are written as 'NA', never as 0. A gene absent from a
   sample's table was not tested there, and recording that as absence is the
   same error as scoring an unverified call present.
"""

import argparse
import csv
import sys
from collections import Counter


def read_gene_order(path):
    """genes.txt order, comments and blanks dropped, de-duplicated."""
    order, seen = [], set()
    with open(path) as fh:
        for line in fh:
            g = line.strip()
            if not g or g.startswith("#"):
                continue
            if g.lower() in seen:
                continue
            seen.add(g.lower())
            order.append(g)
    return order


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--matrix", required=True,
                    help="long matrix: sample<TAB>gene<TAB>present")
    ap.add_argument("--genes", default="",
                    help="genes.txt, used for COLUMN ORDER (default: "
                         "alphabetical, which splits operons apart)")
    ap.add_argument("--out-prefix", default="",
                    help="directory to write into (default: alongside --matrix)")
    args = ap.parse_args()

    rows = list(csv.DictReader(open(args.matrix), delimiter="\t"))
    if not rows:
        print(f"[WARN] {args.matrix} has no data rows; nothing to summarise")
        return 0

    samples = sorted({r["sample"] for r in rows})
    present_genes = {r["gene"] for r in rows}

    order = read_gene_order(args.genes) if args.genes else []
    genes = [g for g in order if g in present_genes]
    # anything called but not in genes.txt still gets a column, at the end
    extra = sorted(present_genes - set(genes))
    if extra:
        print(f"[NOTE] {len(extra)} gene(s) not in --genes, appended: "
              f"{','.join(extra)}")
    genes += extra

    val = {(r["sample"], r["gene"]): r["present"] for r in rows}

    import os
    outdir = args.out_prefix or os.path.dirname(os.path.abspath(args.matrix))
    os.makedirs(outdir, exist_ok=True)

    # ── wide matrix ─────────────────────────────────────────────────────────
    wide = os.path.join(outdir, "gene_presence_matrix_wide.tsv")
    n_missing = 0
    with open(wide, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["sample"] + genes)
        for s in samples:
            row = []
            for g in genes:
                v = val.get((s, g))
                if v is None:
                    n_missing += 1
                    v = "NA"
                row.append(v)
            w.writerow([s] + row)

    if n_missing:
        print(f"[WARN] {n_missing} cell(s) written as NA -- those gene/sample "
              f"pairs are absent from the long matrix, meaning the gene was not "
              f"tested in that genome. NOT the same as absent.")

    # ── companions ──────────────────────────────────────────────────────────
    n = len(samples)
    prev = Counter()
    for (s, g), v in val.items():
        if v == "1":
            prev[g] += 1

    pfile = os.path.join(outdir, "gene_prevalence.tsv")
    with open(pfile, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["gene", "n_genomes", "pct"])
        for g in sorted(genes, key=lambda x: (-prev[x], x)):
            w.writerow([g, prev[g], f"{100 * prev[g] / n:.0f}"])

    gfile = os.path.join(outdir, "genes_per_genome.tsv")
    per = {s: sum(1 for g in genes if val.get((s, g)) == "1") for s in samples}
    with open(gfile, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["sample", "n_genes"])
        for s in sorted(samples, key=lambda x: (-per[x], x)):
            w.writerow([s, per[s]])

    core = [g for g in genes if prev[g] == n]
    none = [g for g in genes if prev[g] == 0]
    print(f"wide matrix: {len(samples)} samples x {len(genes)} genes -> {wide}")
    print(f"  core (all {n}): {', '.join(core) if core else 'none'}")
    print(f"  absent from all: {', '.join(none) if none else 'none'}")
    print(f"  variable: {len(genes) - len(core) - len(none)} gene(s)")
    print(f"written: {pfile}")
    print(f"written: {gfile}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

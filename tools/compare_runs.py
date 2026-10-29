#!/usr/bin/env python3
"""Diff two nf-gene-verify runs: which calls moved, and on what evidence.

    python3 tools/compare_runs.py --old <results-A> --new <results-B>
    python3 tools/compare_runs.py --old A --new B --gene bhp          # one gene
    python3 tools/compare_runs.py --old A --new B --panel             # panels too
    python3 tools/compare_runs.py --old A --new B --csv changes.csv   # for Excel

Reads <run>/<sample>/gene_presence_table_v2.csv from both and reports:

  1. SCOPE          samples and genes only in one run
  2. INTEGRITY      per run: table/matrix agreement, NA cells, present calls
                    with no nr evidence, wide-vs-long agreement
  3. PREVALENCE     per gene, n/29 before and after, biggest movers first
  4. CALL CHANGES   every changed call with round-1 and nr evidence on both
                    sides, grouped by transition
  5. PANEL          (--panel) per-gene reference count and length range

Written because a full re-run changes every number in the table and eyeballing
1,595 rows does not scale. Two habits it enforces:

  * A GENE WHOSE PREVALENCE DID NOT MOVE CAN STILL HAVE CHANGED. Two samples
    swapping present for absent nets to zero. Prevalence is reported alongside
    the per-call churn, never instead of it.

  * EVIDENCE, NOT JUST VERDICTS. 'flagged -> present' means nothing on its own;
    what matters is whether round-1 identity and the nr hit moved with it. A
    verdict that flips while the evidence stands still is a code change, and a
    verdict that flips because identity went 49.7% -> 99.9% is a panel change.
    bhp did the second; it took reading both to know which.
"""

import argparse
import csv
import glob
import os
import statistics
import sys
from collections import Counter, defaultdict

TABLE = "gene_presence_table_v2.csv"


def samples(run):
    return sorted(os.path.basename(os.path.dirname(p))
                  for p in glob.glob(os.path.join(run, "*", TABLE)))


def load(run):
    """{(sample, gene): row}"""
    out = {}
    for p in sorted(glob.glob(os.path.join(run, "*", TABLE))):
        s = os.path.basename(os.path.dirname(p))
        for r in csv.DictReader(open(p)):
            out[(s, r["gene"])] = r
    return out


def read_tsv(path):
    if not os.path.exists(path):
        return None
    with open(path, newline="") as fh:
        rows = list(csv.reader(fh, delimiter="\t"))
    return rows


def integrity(run, rows):
    """The same checks worth running before sharing either run."""
    out = []
    long_p = os.path.join(run, "summary", "gene_presence_matrix.tsv")
    wide_p = os.path.join(run, "summary", "gene_presence_matrix_wide.tsv")
    lm = read_tsv(long_p)
    if lm is None:
        out.append("summary/gene_presence_matrix.tsv missing")
    else:
        long_d = {(r[0], r[1]): r[2] for r in lm[1:]}
        bad = sum(1 for (s, g), r in rows.items()
                  if long_d.get((s, g)) != ("1" if r["final_call"] == "present" else "0"))
        out.append(f"{len(lm)-1} matrix rows; {bad} disagree with the per-sample tables")
        wm = read_tsv(wide_p)
        if wm is None:
            out.append("summary/gene_presence_matrix_wide.tsv missing")
        else:
            hdr = wm[0]
            wide_d = {(r[0], hdr[i]): r[i] for r in wm[1:] for i in range(1, len(hdr))}
            out.append(f"wide == long: {wide_d == long_d}; "
                       f"NA cells: {sum(1 for v in wide_d.values() if v == 'NA')}")
    unv = [f"{s}/{g}" for (s, g), r in rows.items()
           if r["final_call"] == "present" and not r["nr_confirming_accession"].strip()]
    out.append(f"present calls with no nr evidence: {len(unv)}"
               + (f" -- {unv[:4]}" if unv else ""))
    return out


def panel(run):
    """{gene: (n, min_aa, max_aa, median_aa)} from protein_db/manifest.csv"""
    p = os.path.join(run, "protein_db", "manifest.csv")
    if not os.path.exists(p):
        return None
    by = defaultdict(list)
    for r in csv.DictReader(open(p)):
        try:
            by[r["gene"]].append(int(r["length"] or 0))
        except ValueError:
            continue
    return {g: (len(v), min(v), max(v), statistics.median(v))
            for g, v in by.items() if v}


def ev(r):
    """One-line evidence summary for a row."""
    if r is None:
        return "-"
    return (f"r1 {r['round1_call']}/{r['round1_identity_pct'] or '-'}%"
            f"/{r['round1_coverage_pct'] or '-'}%"
            f"  nr {r['nr_identity_pct'] or '-'}%/{r['nr_hit_coverage_pct'] or '-'}%"
            f" rk{r['nr_hit_rank'] or '-'}")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--old", required=True, help="baseline run directory")
    ap.add_argument("--new", required=True, help="run to compare against it")
    ap.add_argument("--gene", default="", help="restrict to one gene")
    ap.add_argument("--panel", action="store_true",
                    help="also diff the reference panels")
    ap.add_argument("--csv", default="", help="write changed calls to this CSV")
    ap.add_argument("--quiet-evidence", action="store_true",
                    help="list transitions only, without per-call evidence")
    args = ap.parse_args()

    old, new = load(args.old), load(args.new)
    if args.gene:
        old = {k: v for k, v in old.items() if k[1] == args.gene}
        new = {k: v for k, v in new.items() if k[1] == args.gene}
    if not old or not new:
        sys.exit("no rows loaded -- check the paths (and --gene spelling)")

    print(f"old: {args.old}\nnew: {args.new}\n")

    # 1 ── scope
    os_, ns = {k[0] for k in old}, {k[0] for k in new}
    og, ng = {k[1] for k in old}, {k[1] for k in new}
    print("── scope")
    print(f"  samples: {len(os_)} -> {len(ns)}"
          + (f"   added {sorted(ns-os_)}" if ns - os_ else "")
          + (f"   dropped {sorted(os_-ns)}" if os_ - ns else ""))
    print(f"  genes:   {len(og)} -> {len(ng)}"
          + (f"   added {sorted(ng-og)}" if ng - og else "")
          + (f"   dropped {sorted(og-ng)}" if og - ng else ""))

    # 2 ── integrity, per run
    print("\n── integrity")
    for label, run, rows in (("old", args.old, old), ("new", args.new, new)):
        if args.gene:
            print(f"  {label}: skipped (--gene restricts the tables)")
            continue
        for line in integrity(run, rows):
            print(f"  {label}: {line}")

    # 3 ── distribution + prevalence
    print("\n── final_call distribution")
    co, cn = Counter(r["final_call"] for r in old.values()), \
             Counter(r["final_call"] for r in new.values())
    for k in sorted(set(co) | set(cn)):
        d = cn[k] - co[k]
        print(f"  {k:<20} {co[k]:>5} -> {cn[k]:<5} {d:+d}" if d else
              f"  {k:<20} {co[k]:>5} -> {cn[k]:<5}  =")

    shared_g = sorted(og & ng)
    print("\n── prevalence (present only; genes that moved, largest first)")
    pv = []
    for g in shared_g:
        a = sum(1 for (s, gg), r in old.items() if gg == g and r["final_call"] == "present")
        b = sum(1 for (s, gg), r in new.items() if gg == g and r["final_call"] == "present")
        if a != b:
            pv.append((abs(b - a), g, a, b))
    if pv:
        for _, g, a, b in sorted(pv, reverse=True):
            print(f"  {g:<10} {a:>3} -> {b:<3} {b-a:+d}")
    else:
        print("  no gene changed prevalence")

    # 4 ── per-call changes
    shared = sorted(set(old) & set(new))
    changed = [k for k in shared if old[k]["final_call"] != new[k]["final_call"]]
    print(f"\n── call changes: {len(changed)} of {len(shared)} shared calls")
    by_t = defaultdict(list)
    for k in changed:
        by_t[(old[k]["final_call"], new[k]["final_call"])].append(k)
    for (a, b), ks in sorted(by_t.items(), key=lambda x: -len(x[1])):
        print(f"\n  {a} -> {b}   ({len(ks)})")
        for s, g in sorted(ks):
            if args.quiet_evidence:
                print(f"    {s:<12} {g}")
            else:
                print(f"    {s:<12} {g:<10} old: {ev(old[(s,g)])}")
                print(f"    {'':<12} {'':<10} new: {ev(new[(s,g)])}")

    # A verdict that moved while its evidence did not is a CODE change; one that
    # moved with the evidence is a DATA or PANEL change. Worth separating,
    # because only the first can silently affect genes nobody was looking at.
    same_ev = [k for k in changed
               if (old[k]["round1_identity_pct"], old[k]["round1_coverage_pct"])
               == (new[k]["round1_identity_pct"], new[k]["round1_coverage_pct"])]
    if changed:
        print(f"\n  of these, {len(same_ev)} changed with IDENTICAL round-1 evidence "
              f"(a scoring/logic change, not a panel change)")
        for s, g in sorted(same_ev)[:12]:
            print(f"    {s:<12} {g:<10} {old[(s,g)]['final_call']} -> {new[(s,g)]['final_call']}")

    # 5 ── panels
    if args.panel:
        po, pn = panel(args.old), panel(args.new)
        print("\n── reference panels")
        if po is None or pn is None:
            print("  protein_db/manifest.csv missing in one run")
        else:
            rows = []
            for g in sorted(set(po) | set(pn)):
                a, b = po.get(g), pn.get(g)
                if a != b:
                    rows.append((g, a, b))
            if not rows:
                print("  identical")
            for g, a, b in rows:
                f = lambda t: "absent" if not t else f"n={t[0]} {t[1]}-{t[2]} aa (med {t[3]:.0f})"
                print(f"  {g:<10} {f(a)}  ->  {f(b)}")

    # csv
    if args.csv:
        with open(args.csv, "w", newline="") as fh:
            w = csv.writer(fh, lineterminator="\n")
            w.writerow(["sample", "gene", "old_call", "new_call",
                        "old_r1_id", "new_r1_id", "old_r1_cov", "new_r1_cov",
                        "old_nr_id", "new_nr_id", "old_nr_cov", "new_nr_cov",
                        "old_nr_rank", "new_nr_rank", "new_nr_verification"])
            for s, g in sorted(changed):
                o, n = old[(s, g)], new[(s, g)]
                w.writerow([s, g, o["final_call"], n["final_call"],
                            o["round1_identity_pct"], n["round1_identity_pct"],
                            o["round1_coverage_pct"], n["round1_coverage_pct"],
                            o["nr_identity_pct"], n["nr_identity_pct"],
                            o["nr_hit_coverage_pct"], n["nr_hit_coverage_pct"],
                            o["nr_hit_rank"], n["nr_hit_rank"],
                            n["nr_verification"]])
        print(f"\nwrote {args.csv} ({len(changed)} changed calls)")
    return 0


if __name__ == "__main__":
    sys.exit(main())

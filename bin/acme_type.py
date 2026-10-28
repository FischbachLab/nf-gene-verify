#!/usr/bin/env python3
"""Assign an ACME type from per-component gene calls.

ACME (arginine catabolic mobile element) is a composite element, not a gene.
Its types are defined by WHICH CLUSTERS are present, not by distinct sequences:

    type I    arc + opp3
    type II   arc, no opp3
    type III  opp3, no arc
    type IV   kdp, no opp3
    type V    kdp + arc + opp3

Because the types share the same genes, they cannot be separated by similarity
search -- a 'type I' panel and a 'type II' panel contain the same arc proteins,
so both would match and both would report present. Typing is a logical rule over
component presence and LOCATION, which is what this script applies.

Two traps this script exists to avoid:

1. arc and kdp are ALSO NATIVE operons in staphylococci; the ACME copy is a
   second one. On LM088-RCB001 arcA matched the native operon at 99.8%, so a
   combined 'ACME' BLAST target would have reported ACME present at 99.8%
   identity off a housekeeping operon. Presence of a cluster is therefore never
   sufficient -- it must be co-located with the element.

2. Native clusters are scattered. On LM088-RCB001 the arc genes sat at ~160 kb,
   ~277 kb and ~1291 kb -- three separate loci, 604 kb from the opp3 cluster.
   Single-linkage clustering within --window separates that from a real element.

Characterised ACMEs are 31-34 kb, which sets the default --window of 34000.

What this script CANNOT do: confirm that a cluster is carried on ACME rather
than being the native copy. ACME integrates at orfX (rlmH), and locating orfX
needs annotation this pipeline does not produce. A lone arc or kdp cluster is
therefore reported as AMBIGUOUS, never as a type.

Usage:
    acme_type.py --calls gene_calls.tsv [--out acme_type.tsv] [--window 34000]
"""

import argparse
import csv
import os
import sys

# Component clusters. argR is grouped with arc: it is the arc operon's
# regulator and sits in the element in ACME-I/II.
COMPONENTS = {
    "arc":  ["arcA", "arcB", "arcC", "arcD", "argR"],
    "opp3": ["opp3A", "opp3B", "opp3C", "opp3D", "opp3E"],
    "kdp":  ["kdpA", "kdpB", "kdpC", "kdpD", "kdpE"],
}

# Clusters whose presence is uninformative on its own because a native copy is
# conserved in the species.
#
# opp3 is in this set, which is NOT obvious. The chromosomal oligopeptide
# permease systems opp-1, opp-2, opp-3 and opp-5A are highly conserved in
# S. epidermidis, and the ACME-borne operon shares the 'opp3' label with the
# chromosomal opp-3 while being a different locus. An opp3 cluster matching
# S. epidermidis references at 100% identity is therefore at least as consistent
# with the native system as with ACME -- on LM088-RCB001 exactly that cluster
# sits at 1,894,753-1,900,441.
#
# The consequence is that ALL THREE component groups have conserved native
# counterparts, so no combination of presence calls alone can establish ACME.
# Only co-location of two groups, or the orfX position, can.
NATIVE_PRONE = {"arc", "kdp", "opp3"}

# A round-1 call counts as "the protein is there" for clustering purposes.
# 'divergent' and 'partial' are included deliberately: a real element component
# may be diverged or interrupted, and excluding them would bias towards calling
# elements absent. 'absent' and 'no_reference' carry no coordinates at all.
FOUND_CALLS = {"present", "divergent", "partial"}

# (arc, opp3, kdp) -> (type, note). Keys are booleans for "cluster is part of
# this element".
TYPE_TABLE = {
    (True,  True,  True):  ("V",   "kdp + arc + opp3"),
    (True,  True,  False): ("I",   "arc + opp3"),
    (True,  False, False): ("II",  "arc without opp3"),
    (False, True,  False): ("III", "opp3 without arc"),
    (False, False, True):  ("IV",  "kdp without opp3"),
    # Combinations with no published type
    (False, True,  True):  (None,  "opp3 + kdp without arc -- no published type"),
    (True,  False, True):  (None,  "arc + kdp without opp3 -- resembles IV/V but "
                                   "incomplete"),
}


def load_calls(path):
    """gene -> dict(call, contig, start, end) for genes with coordinates."""
    out = {}
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            coords = (r.get("coords") or "").strip()
            if r.get("call") not in FOUND_CALLS or not coords:
                out[r["gene"]] = None
                continue
            try:
                a, b = (int(x) for x in coords.split("-"))
            except ValueError:
                out[r["gene"]] = None
                continue
            out[r["gene"]] = {
                "call": r["call"],
                "contig": r.get("contig", ""),
                "start": min(a, b),
                "end": max(a, b),
                "identity": r.get("identity_pct", ""),
            }
    return out


def cluster(members, window):
    """Single-linkage grouping of located members by contig and proximity.

    Returns a list of clusters, each a dict with contig, start, end, genes.
    """
    by_contig = {}
    for gene, m in members:
        by_contig.setdefault(m["contig"], []).append((gene, m))
    clusters = []
    for contig, items in by_contig.items():
        items.sort(key=lambda t: t[1]["start"])
        cur = None
        for gene, m in items:
            if cur and m["start"] - cur["end"] <= window:
                cur["end"] = max(cur["end"], m["end"])
                cur["genes"].append(gene)
            else:
                if cur:
                    clusters.append(cur)
                cur = {"contig": contig, "start": m["start"], "end": m["end"],
                       "genes": [gene]}
        if cur:
            clusters.append(cur)
    return clusters


def gap(c1, c2):
    """Distance in bp between two clusters; None if on different contigs."""
    if c1["contig"] != c2["contig"]:
        return None
    if c1["end"] < c2["start"]:
        return c2["start"] - c1["end"]
    if c2["end"] < c1["start"]:
        return c1["start"] - c2["end"]
    return 0


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--calls", required=True, help="gene_calls.tsv from round 1")
    ap.add_argument("--out", default="",
                    help="report TSV (default: acme_type.tsv beside --calls)")
    ap.add_argument("--window", type=int, default=34000,
                    help="max bp between components of one element (default "
                         "34000; characterised ACMEs are 31-34 kb)")
    ap.add_argument("--min-cluster", type=int, default=2,
                    help="members required before a group counts as a cluster "
                         "(default 2; one lone gene is not an operon)")
    args = ap.parse_args()

    calls = load_calls(args.calls)
    missing = [g for grp in COMPONENTS.values() for g in grp if g not in calls]

    # ── locate each component cluster ───────────────────────────────────────
    found, report = {}, []
    for comp, genes in COMPONENTS.items():
        located = [(g, calls[g]) for g in genes if calls.get(g)]
        not_screened = [g for g in genes if g not in calls]
        clusters = [c for c in cluster(located, args.window)
                    if len(c["genes"]) >= args.min_cluster]
        found[comp] = clusters
        report.append({
            "component": comp,
            "genes_screened": len(genes) - len(not_screened),
            "genes_located": len(located),
            "clusters": len(clusters),
            "detail": "; ".join(
                f"{c['contig']}:{c['start']}-{c['end']} "
                f"({'+'.join(c['genes'])})" for c in clusters) or "-",
            "not_screened": ",".join(not_screened) or "-",
        })

    # ── which clusters belong to ONE element? ───────────────────────────────
    # Take every pairing of clusters from different components and keep those
    # within --window. A component is part of the element only if it is
    # co-located with another component's cluster.
    elements, notes = [], []
    comps = [c for c in COMPONENTS if found[c]]
    linked = {c: set() for c in COMPONENTS}
    for i, a in enumerate(comps):
        for b in comps[i + 1:]:
            for ca in found[a]:
                for cb in found[b]:
                    d = gap(ca, cb)
                    if d is not None and d <= args.window:
                        linked[a].add(b)
                        linked[b].add(a)
                        elements.append((a, b, ca, cb, d))

    co_located = {c for c in COMPONENTS if linked[c]}

    # ── assign ──────────────────────────────────────────────────────────────
    if not any(found.values()):
        verdict, confidence = "absent", "no component cluster detected"
    elif co_located:
        key = tuple(c in co_located for c in ("arc", "opp3", "kdp"))
        t, why = TYPE_TABLE.get(key, (None, "unrecognised combination"))
        if t:
            verdict = f"type {t}"
            confidence = f"{why}; components co-located within {args.window} bp"
        else:
            verdict, confidence = "ambiguous", why
    else:
        # Clusters exist but none are co-located: either native copies, or a
        # single-cluster element (types II, III, IV are single-cluster by
        # definition).
        present = [c for c in COMPONENTS if found[c]]
        informative = [c for c in present if c not in NATIVE_PRONE]
        if informative:
            # Reachable only if NATIVE_PRONE is narrowed; kept so a future
            # ACME-specific component yields a type rather than silence.
            verdict = "ambiguous"
            confidence = (f"cluster(s) {', '.join(informative)} detected but not "
                          f"co-located with another component")
        else:
            verdict = "ambiguous"
            confidence = (
                f"cluster(s) detected ({', '.join(present)}) but none "
                f"co-located within {args.window} bp. Every ACME component "
                f"group has a conserved native counterpart in S. epidermidis "
                f"(arc, kdp, and the chromosomal opp-3 system, which shares "
                f"the 'opp3' label with the ACME operon), so presence alone "
                f"cannot establish the element -- a scattered set of clusters "
                f"is what an ACME-NEGATIVE genome looks like")

    if verdict.startswith("type") or verdict == "ambiguous":
        notes.append("ACME integrates at orfX (rlmH); confirming that a cluster "
                     "is element-borne rather than the native copy requires the "
                     "orfX position, which this pipeline does not determine.")
    for comp in COMPONENTS:
        ns = [g for g in COMPONENTS[comp] if g not in calls]
        if ns:
            notes.append(f"{comp}: {len(ns)} gene(s) not in gene_calls.tsv "
                         f"({','.join(ns)}) -- add them to genes.txt before "
                         f"trusting a negative for this component")

    out = args.out or os.path.join(os.path.dirname(os.path.abspath(args.calls)),
                                   "acme_type.tsv")
    with open(out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["acme_type", "basis", "window_bp"])
        w.writerow([verdict, confidence, args.window])
        w.writerow([])
        w.writerow(["component", "genes_screened", "genes_located", "clusters",
                    "detail", "not_screened"])
        for r in report:
            w.writerow([r["component"], r["genes_screened"], r["genes_located"],
                        r["clusters"], r["detail"], r["not_screened"]])
        if notes:
            w.writerow([])
            for n in notes:
                w.writerow(["NOTE", n])

    print(f"ACME: {verdict}")
    print(f"  basis: {confidence}")
    for r in report:
        print(f"  {r['component']:5} {r['genes_located']}/{r['genes_screened']} "
              f"located, {r['clusters']} cluster(s): {r['detail']}")
        if r["not_screened"] != "-":
            print(f"        not screened: {r['not_screened']}")
    for a, b, ca, cb, d in elements:
        print(f"  co-located: {a} and {b} are {d} bp apart")
    for n in notes:
        print(f"  NOTE: {n}")
    print(f"written: {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

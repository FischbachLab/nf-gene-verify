#!/usr/bin/env python3
"""Flag called loci whose blastx HSPs sit in POSITIONALLY SEGREGATED frames.

blastx translates the query in six frames and reports HSPs from any of them;
call_genes.py clusters by genomic distance and never looks at frame. A
frameshifted pseudogene is therefore indistinguishable from an intact gene:
round 1 sees the right sequence over the right length, round 2 confirms the
sequence is the right protein, and both are true statements about SEQUENCE that
say nothing about READING FRAME.

An HSP's query start mod 3 is its frame. For one intact ORF every HSP shares a
frame. If instead each frame occupies its own stretch of the locus, with small
overlaps at the joins, that is a frameshift.

Mixed frames with INTERLEAVED positions are not evidence of anything -- repeat
proteins produce alternative alignments of the same region. Segregation is the
signal, so this reports the per-frame span and the overlap at each junction.
"""
import csv, sys
from collections import defaultdict

def frames_for_locus(blastx, gene, contig, gs, ge):
    hs = []
    for l in open(blastx):
        f = l.rstrip("\n").split("\t")
        if len(f) < 12 or f[1].split("|")[1] != gene or f[0] != str(contig):
            continue
        qs, qe = int(f[6]), int(f[7])
        lo, hi = min(qs, qe), max(qs, qe)
        if hi < gs or lo > ge:
            continue
        hs.append((lo, hi, "+" if qs < qe else "-", (lo - 1) % 3))
    return hs

def report(name, calls, blastx, min_bp=600):
    print(f"\n{'='*72}\n{name}\n{'='*72}")
    flagged = []
    for r in csv.DictReader(open(calls), delimiter="\t"):
        if not r.get("coords") or r["call"] in ("absent", "no_reference"):
            continue
        gs, ge = (int(x) for x in r["coords"].split("-"))
        if ge - gs + 1 < min_bp:
            continue
        hs = frames_for_locus(blastx, r["gene"], r["contig"], gs, ge)
        if not hs:
            continue
        byf = defaultdict(list)
        for lo, hi, st, fr in hs:
            byf[(st, fr)].append((lo, hi))
        spans = {k: (min(a for a, _ in v), max(b for _, b in v), len(v))
                 for k, v in byf.items()}
        if len(spans) < 2:
            print(f"  {r['gene']:7} {ge-gs+1:>7,} bp  single frame  ok")
            continue
        # segregated if sorted spans barely overlap
        ordered = sorted(spans.items(), key=lambda kv: kv[1][0])
        overlaps = [ordered[i][1][1] - ordered[i+1][1][0] + 1
                    for i in range(len(ordered) - 1)]
        # The SIGN of the junction discriminates two different situations:
        #
        #   OVERLAP (> 0)  -- HSPs from both frames span the same bases, which
        #                     is what happens at a frameshift: references align
        #                     across the break in either frame.
        #   GAP (< 0)      -- the blocks are disjoint, i.e. two ADJACENT GENES
        #                     in different frames, not one interrupted gene.
        #
        # Observed: ebh/embp junctions overlap by 31-112 bp (frameshift),
        # whereas SCCmec's blocks are 1 bp apart and opp3D's 86 bp apart --
        # both of those are gene boundaries. SCCmec's two blocks are blaR1 and
        # blaI, which is why that locus was never SCCmec to begin with.
        seg = all(o < 0.25 * (ge - gs + 1) for o in overlaps)
        if not seg:
            tag = "interleaved (repeat artefact)"
        elif all(o > 0 for o in overlaps):
            tag = "<<< FRAMESHIFT (junctions overlap)"
        else:
            tag = "adjacent genes in different frames (junctions are gaps)"
        print(f"  {r['gene']:7} {ge-gs+1:>7,} bp  {len(spans)} frames  {tag}")
        for (st, fr), (a, b, n) in ordered:
            print(f"            {st}{fr}  {a:>10,}-{b:<10,} {b-a+1:>7,} bp  n={n}")
        if len(ordered) > 1:
            print(f"            junction overlap(s): "
                  f"{', '.join(f'{o:,} bp' for o in overlaps)}")
        if seg and all(o > 0 for o in overlaps):
            flagged.append(r["gene"])
    print(f"\n  FRAMESHIFT CANDIDATES: {flagged or 'none'}")
    return flagged

if __name__ == "__main__":
    report(sys.argv[1], sys.argv[2], sys.argv[3])

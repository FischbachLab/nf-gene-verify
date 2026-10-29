#!/usr/bin/env python3
"""Call gene presence/absence from blastx hits (genome vs protein reference set).

Hits are clustered by genomic distance (gap <= 5 kb) so distant weak
cross-homology hits do not merge into fake loci. Within a cluster, coverage
is computed per reference accession (merged union of aligned reference
intervals); identity is the length-weighted mean.

Calling thresholds (standard): present = >=90% identity over >=80% reference
coverage; divergent = full coverage, <90% identity; partial = <80% coverage;
homology_only = <40% identity, i.e. a match to a different member of the same
superfamily rather than an allele of this gene (checked first, since it would
otherwise be reported as divergent/partial and read as a weak positive);
absent = no qualifying hit. Genes with no reference proteins in the manifest
(build_protein_db.py found none) are reported as 'no_reference' — they were
never searched, so 'absent' would be a false negative.
"""
import argparse
import csv
import sys
import re
from collections import defaultdict

# GAP lowered 5000 -> 2000. At 5000 the aap cluster on LM088-RCB001 absorbed a
# detached 516 bp block of G5/E repeats sitting 3,847 bp away, reporting the
# locus as 33856-42361 (8,506 bp) when the gene is 38219-42361 (4,143 bp) --
# a 4.4 kb overstatement, with coverage inflated to 85.7% from 81.2%. At 2000
# the reported coordinates match nr's confirming hit EXACTLY (38219-42361),
# which is independent corroboration rather than a tuned guess.
#
# Validated across two panels, 56 gene-calls: aap's coordinates are the only
# thing that changes and no CALL changes anywhere. The 30,529 bp ebh/embp gene
# -- the obvious casualty of a tighter gap -- is unaffected at every value
# tested down to 500.
#
# Do not go below ~1400: at 1000 the aap gene splits at its own internal
# 1,341 bp gap, losing a block. The safe window is wider than the largest
# INTERNAL unaligned stretch and narrower than the distance to the nearest
# unrelated element.
# nr/manifest fields can be enormous: a MULTISPECIES defline lists every
# redundant accession and organism in ONE field, which can exceed Python's
# default 128 KB csv cap and raise '_csv.Error: field larger than field limit'.
_lim = sys.maxsize
while True:
    try:
        csv.field_size_limit(_lim)
        break
    except OverflowError:
        _lim //= 10

ID_MIN, COV_MIN, GAP = 90.0, 80.0, 2000
# Below this identity the match is superfamily-level, not allelic. See the
# homology_only branch in main() for why 40% and why it is a judgement call.
HOMOLOGY_MAX = 40.0


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
                "slen": int(f[5]),
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


def frame_status(cluster, min_bp=600, acc=None):
    """('', '') for a clean locus, or (summary, junction_coords) for a frameshift.

    blastx translates the query in six frames and reports HSPs from any of them;
    cluster_by_genome() groups by genomic distance and never looks at frame. A
    frameshifted pseudogene is therefore indistinguishable from an intact gene:
    the right sequence is present over the right length, and nr confirms it is
    the right protein. Both are true statements about SEQUENCE and say nothing
    about READING FRAME.

    An HSP's query start mod 3 is its frame. For one intact ORF every HSP shares
    a frame. Two signals separate the cases:

      frames occupy SEPARATE stretches of the locus  -> not one ORF
      their junctions OVERLAP rather than gap        -> frameshift, not two genes

    The second is what distinguishes a frameshift from two adjacent genes. At a
    frameshift, references align across the break in either frame, so HSPs from
    both frames cover the same bases. Adjacent genes give disjoint blocks.
    Observed on LM088-RCB001: ebh/embp junctions overlap by 31-112 bp (a real
    frameshift, at identical coordinates from two independent panels), while
    SCCmec's blocks sit 1 bp apart (blaR1 + blaI) and opp3D's 86 bp apart (two
    neighbouring ABC transporters).

    Interleaved frames are NOT reported: repeat proteins routinely produce
    alternative alignments of the same region in several frames, which says
    nothing about the reading frame of the gene.
    """
    # Only the WINNING reference's HSPs. cluster_by_genome() merges HSPs from
    # every reference, so a single over-length or cross-matching reference can
    # drag in a neighbouring gene whose frame differs -- and the check then
    # reports a frameshift that does not exist. icaC did exactly this in 9
    # genomes, via a 635 aa S. pseudintermedius record that never won a call.
    # The reference that actually made the call is the only one whose frame
    # describes the gene.
    hsps = [h for h in cluster["hsps"] if acc is None or h["acc"] == acc]
    if not hsps:
        return "", ""

    # An INTACT ORF ends the question before it is asked.
    #
    # If one HSP aligns the reference END TO END -- subject residue 1 through
    # slen, in a single frame, no internal break -- then a complete, unbroken
    # copy of the protein is encoded here, whatever else the cluster contains.
    #
    # blaZ in NIH051668 is the case: the gene is DUPLICATED, one copy intact at
    # 64312-65154 (subject 1-281 of 281, 98.2%, frame +0) and one frameshifted
    # at 64979-65820 (subject 1-140 in frame +1, then 141-281 in frame +0).
    # GAP = 2000 merges the two copies into one cluster, the check sees frames
    # +0 and +1 with an overlapping junction, and a beta-lactamase with an
    # intact ORF was reported as present_frameshift -- scoring 0 for an
    # acquired resistance gene that is plainly there.
    #
    # END TO END is the whole rule; there is no coverage threshold to tune, and
    # that matters. The runners-up all fall short of a complete subject and all
    # stay flagged, correctly: atlE in NIHLM095 reaches 90.6% but its subject
    # starts at residue 137 (1-127 sits in another frame -- a real frameshift
    # at codon ~130), embp in NIHLM021 reaches 98.6% but stops at 2826 of 2867,
    # and sdrF in NIH04008 covers 1-1126 of 1784. A 90% or 98% threshold would
    # have cleared the first two wrongly.
    for h in hsps:
        s_lo, s_hi = min(h["sstart"], h["send"]), max(h["sstart"], h["send"])
        if h.get("slen") and s_lo == 1 and s_hi == h["slen"]:
            return "", ""

    lo = min(min(h["qstart"], h["qend"]) for h in hsps)
    hi = max(max(h["qstart"], h["qend"]) for h in hsps)
    if hi - lo + 1 < min_bp:
        return "", ""
    byf = {}
    bysubj = {}
    for h in hsps:
        a, b = min(h["qstart"], h["qend"]), max(h["qstart"], h["qend"])
        key = ("+" if h["qstart"] < h["qend"] else "-", (a - 1) % 3)
        s = byf.setdefault(key, [a, b])
        s[0], s[1] = min(s[0], a), max(s[1], b)
        slo, shi = min(h["sstart"], h["send"]), max(h["sstart"], h["send"])
        t = bysubj.setdefault(key, [slo, shi])
        t[0], t[1] = min(t[0], slo), max(t[1], shi)
    if len(byf) < 2:
        return "", ""

    # ── is the frame assignment interpretable at all? ───────────────────────
    #
    # Everything below this point assumes the HSPs TILE ONE ORF, so that a
    # change of frame part-way along means the reading frame changed. For a
    # repeat protein that assumption fails: blastx returns alternative
    # alignments of the same residues in several frames, and the check reads
    # the alternatives as a break. Two independent signs of that, each with a
    # wide margin on the 29-genome v5 data.
    #
    # 1. THE FRAMES ENCODE THE SAME RESIDUES. A real frameshift partitions the
    #    protein -- frame A covers reference residues 1..k, frame B covers
    #    k+1..end. If one frame's subject range is instead CONTAINED in
    #    another's, both frames encode the same part of the protein, which no
    #    single frameshift can produce. msrA in three genomes had frames
    #    spanning residues 1-479 and 34-488 of a 488 aa ABC-F protein: 98% of
    #    the shorter range inside the longer, i.e. the whole protein aligned
    #    twice. embp in NIHLM040 had 1-534 nested inside 1-6337. Every genuine
    #    case sits at 13% or below -- blaZ at 0% (1-140 | 141-281), sdrH at
    #    0.4%, atlE at 3.8%.
    #
    # 2. THE TILING IS REDUNDANT. Summing the HSP subject lengths and dividing
    #    by the residues actually covered gives ~1.0 when each part of the
    #    protein aligns once, whatever its frame. embp in NIHLM021 gave 2.13
    #    from 11 HSPs re-aligning the same region, and its whole claim to a
    #    frameshift was ONE 147-bit HSP at 38.7% identity in a second frame.
    #    NIHLM040 gave 7.06 from 18 HSPs. Every other call is at most 1.03.
    #
    # Reported as clean rather than flagged: the gene is present, and round 2
    # still has to confirm it. What is NOT claimed is that the reading frame
    # was checked -- for a 10,620 aa repeat adhesin, it cannot be, and saying
    # 'frameshift' there is a stronger statement than the data supports.
    ranges = sorted(bysubj.values())
    for i in range(len(ranges) - 1):
        for j in range(i + 1, len(ranges)):
            # NOT lo/hi -- those hold the QUERY span computed above and are
            # read again below as `span`. Shadowing them here silently turned
            # the segregation test into a comparison against a 5 bp subject
            # overlap and suppressed all 22 frameshift calls.
            o_lo = max(ranges[i][0], ranges[j][0])
            o_hi = min(ranges[i][1], ranges[j][1])
            ov = o_hi - o_lo + 1
            if ov <= 0:
                continue
            shorter = min(ranges[i][1] - ranges[i][0] + 1,
                          ranges[j][1] - ranges[j][0] + 1)
            if shorter and ov / shorter >= 0.5:
                return "", ""

    ivs = [[min(h["sstart"], h["send"]), max(h["sstart"], h["send"])] for h in hsps]
    union = sum(e - s + 1 for s, e in merge_intervals([list(i) for i in ivs]))
    total = sum(e - s + 1 for s, e in ivs)
    if union and total / union >= 1.5:
        return "", ""
    ordered = sorted(byf.items(), key=lambda kv: kv[1][0])
    overlaps = [ordered[i][1][1] - ordered[i + 1][1][0] + 1
                for i in range(len(ordered) - 1)]
    span = hi - lo + 1
    segregated = all(o < 0.25 * span for o in overlaps)
    if not (segregated and all(o > 0 for o in overlaps)):
        return "", ""
    juncs = [ordered[i + 1][1][0] for i in range(len(ordered) - 1)]
    return (f"frameshift:{len(ordered)}frames",
            ";".join(str(j) for j in juncs))


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


def orf_span(cluster, acc):
    """Genomic span of the gene, as HSPs sharing the winning reference's FRAME.

    Neither obvious choice works on its own:

      whole cluster      -- cluster_by_genome() merges HSPs from EVERY
                            reference, so one over-length or cross-matching
                            reference drags the span into a neighbouring gene.
                            icaC read 1,921 bp for a 355 aa protein because a
                            635 aa mis-annotated S. pseudintermedius record
                            aligned on through the adjacent ica gene.

      winning ref only   -- truncates any gene LONGER than its references.
                            Embp is ~10,200 aa and no single reference spans it
                            (panel 2839-10624 aa), so the winner covers a
                            fragment: v3 reported a 4,340-codon locus for a gene
                            v2 had at 10,129.

    Frame separates the two cases. An HSP's query start mod 3 is its blastx
    frame, and a neighbouring gene is in a DIFFERENT frame, while other
    references covering more of the SAME ORF are in the same one. So: take every
    HSP, from any reference, that shares the winning reference's strand and
    frame. icaC then excludes the neighbour (355 codons) and embp keeps its full
    extent.

    A genuine frameshift does truncate the span at the junction, which is
    correct -- the ORF ends there, and frame_status/frameshift_at record it.
    """
    win = [h for h in cluster["hsps"] if h["acc"] == acc]
    if not win:
        win = cluster["hsps"]
    def key(h):
        a = min(h["qstart"], h["qend"])
        return ("+" if h["qstart"] < h["qend"] else "-", (a - 1) % 3)
    frames = {key(h) for h in win}
    same = [h for h in cluster["hsps"] if key(h) in frames] or win
    lo = min(min(h["qstart"], h["qend"]) for h in same)
    hi = max(max(h["qstart"], h["qend"]) for h in same)
    return lo, hi


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
    ap.add_argument("--frame-check-min-bp", type=int, default=600,
                    help="only inspect reading frames for loci at least this "
                         "long (default 600; 0 disables the check). Short loci "
                         "give too few HSPs to tell a frameshift from noise.")
    ap.add_argument("--copy-identity-margin", type=float, default=5.0,
                    help="when a gene has several copies, an INTACT copy is "
                         "preferred over a frameshifted one if their coverage "
                         "is equal and their identity is within this many "
                         "percentage points (default 5.0; 0 disables). Guards "
                         "the case where the copies are not the same gene: "
                         "arcB's frameshifted locus (82.7%%) and the native "
                         "ornithine carbamoyltransferase (59.3%%) both reach "
                         "100%% coverage, and only the margin keeps them apart.")
    ap.add_argument("--homology-max", type=float, default=HOMOLOGY_MAX,
                    help="below this %% identity the best match is reported as "
                         "'homology_only' rather than 'divergent' or 'partial' "
                         "(default 40.0; 0 disables). It means the hit is a "
                         "different member of the same superfamily, not an "
                         "allele of this gene -- a native MFS pump matched by "
                         "tetK, a native methyltransferase matched by ermA/B/C. "
                         "Scores 0 in the presence matrix, same as absent.")
    args = ap.parse_args()

    # Normalise the gene list exactly as build_protein_db.py does. Without
    # this the two steps disagree: build_protein_db.py strips whitespace and
    # de-duplicates, so a gene listed twice (or with a trailing space) yields
    # ONE panel but TWO rows here -- two identical calls, two extracted loci and
    # two nr searches for the same gene. genes_3.txt arrived with six trailing
    # spaces and a duplicate (fbe/sdrG are the same protein), so this is not
    # hypothetical.
    seen_g, uniq = set(), []
    for g in args.genes:
        g = g.strip()
        if not g:
            continue
        if g.lower() in seen_g:
            print(f"[NOTE] gene {g!r} listed more than once; calling it once")
            continue
        seen_g.add(g.lower())
        uniq.append(g)
    args.genes = uniq

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
                         "best_accession": "", "organism": "", "n_hits": 0,
                         "frame_status": "", "frameshift_at": ""})
            continue

        ghits = by_gene.get(gene, [])
        if not ghits:
            rows.append({"gene": gene, "call": "absent", "identity_pct": "",
                         "coverage_pct": "", "contig": "", "coords": "",
                         "best_accession": "", "organism": "", "n_hits": 0,
                         "frame_status": "", "frameshift_at": ""})
            continue

        # ── which COPY of the gene to report ────────────────────────────────
        # Ranked (coverage, intact, identity), not (coverage, identity).
        #
        # A gene can be present more than once, and the copies need not be in
        # the same condition. blaZ in NIH051668 sits on two contigs: an INTACT
        # copy on NZ_AKHK01000035.1 (100% coverage, 98.6% identity, one frame,
        # a single HSP spanning subject 1-281 of 281) and a FRAMESHIFTED copy
        # on NZ_AKHK01000041.1 (100% coverage, 99.6% identity, subject 1-140 in
        # one frame then 141-281 in another). Ranking on identity alone picked
        # the broken copy over the working one for the sake of 1.0% identity,
        # and the genome was reported present_frameshift -- scoring 0 for an
        # acquired beta-lactamase that is intact and presumably functional.
        #
        # Coverage still dominates, so a short intact fragment cannot displace
        # a full-length call; 'intact' only breaks ties BETWEEN COPIES THAT
        # COVER THE REFERENCE EQUALLY WELL, where an unbroken reading frame is
        # the better evidence by any standard. Identity still breaks the
        # remaining ties.
        #
        # What this deliberately does NOT do is hide the broken copy: the
        # frameshifted locus is still in blastx_raw.tsv, and a genome carrying
        # both an intact and a disrupted copy is worth knowing about. It is the
        # presence/absence verdict that should follow the intact one.
        scored = []
        for cluster in cluster_by_genome(ghits, args.gap):
            cov, ident, acc = score_cluster(cluster, manifest)
            lo_hi = orf_span(cluster, acc)
            scored.append({"cov": cov, "ident": ident, "acc": acc,
                           "intact": not frame_status(cluster,
                                                      args.frame_check_min_bp, acc)[0],
                           "cluster": cluster,
                           "contig": cluster["contig"],
                           "qmin": lo_hi[0], "qmax": lo_hi[1]})

        best = max(scored, key=lambda c: (c["cov"], c["ident"]))

        # Prefer an INTACT copy of the same gene over a frameshifted one, but
        # only when the two are plausibly the same gene: equal coverage and
        # within --copy-identity-margin of each other.
        #
        # blaZ in NIH051668 is on two contigs -- intact on NZ_AKHK01000035.1
        # (100% cov, 98.6% id, one frame, a single HSP over subject 1-281 of
        # 281) and frameshifted on NZ_AKHK01000041.1 (100% cov, 99.6% id,
        # subject 1-140 then 141-281 in two frames). Ranking on identity alone
        # chose the broken copy for 1.0% identity and scored an intact,
        # presumably functional beta-lactamase as 0. NIH06004 is the same.
        #
        # The margin is what keeps this honest, and it was not optional: with
        # an unconditional tie-break, arcB in NIH04008 moved from a
        # frameshifted cluster at 82.7% identity to an intact one at 59.3% --
        # 23 points apart, which is a different protein (the native anabolic
        # ornithine carbamoyltransferase), not a second copy. Coverage must
        # match exactly and identity must be close, or the frameshifted call
        # stands and the gene is reported as what round 1 actually found.
        if not best["intact"]:
            rivals = [c for c in scored
                      if c["intact"] and c["cov"] == best["cov"]
                      and c["ident"] >= best["ident"] - args.copy_identity_margin]
            if rivals:
                alt = max(rivals, key=lambda c: c["ident"])
                print(f"[NOTE] {gene}: intact copy at {alt['contig']}:"
                      f"{alt['qmin']}-{alt['qmax']} ({alt['ident']:.1f}% id) preferred "
                      f"over a frameshifted copy at {best['contig']}:"
                      f"{best['qmin']}-{best['qmax']} ({best['ident']:.1f}% id)")
                best = alt

        # 'qmin'/'qmax' above are the ORF span: HSPs sharing the winning
        # reference's frame, from any reference. See orf_span().
        ident, cov = best["ident"], best["cov"]
        # ── homology_only ───────────────────────────────────────────────────
        # Checked BEFORE divergent/partial, because 'divergent' was conflating
        # two unrelated situations. On LM088-RCB001 it labelled arcB (83.1%, a
        # real cross-species orthologue) and tetK (24.9%, fold-level similarity
        # to a native MFS pump) identically. Reading that table without the
        # identity column suggests eight divergent resistance genes where there
        # are none.
        #
        # Below --homology-max the match is to a different member of the same
        # SUPERFAMILY, not to an allele of this gene: ermA/B/C landed on one
        # native rRNA methyltransferase at 24-26%, msrA on a native ABC-F
        # protein at 26%, opp3D/E on unrelated ABC transporters at ~29%. For
        # scale, a genuine copy matches these same panels at ~100% -- opp3B
        # scored 100.0% over 308/308 aa.
        #
        # The default of 40% is deliberately conservative and is a JUDGEMENT
        # CALL, not a natural boundary: embp (42.2%) is a real S. aureus/
        # S. epidermidis orthologue and dfrG (39.5%) is background, and those
        # sit 2.7 points apart. Raise it and real orthologues get mislabelled;
        # lower it and background leaks back in. Round 2 remains the arbiter --
        # it independently rejected mecA and msrA -- so this label is a
        # signpost, not a verdict.
        if ident < args.homology_max:
            call = "homology_only"
        elif ident >= args.identity_min and cov >= args.coverage_min:
            call = "present"
        elif cov >= args.coverage_min:
            call = "divergent"
        else:
            call = "partial"

        fs, fs_at = frame_status(best["cluster"], args.frame_check_min_bp,
                                 acc=best["acc"])
        if fs:
            print(f"[FRAMESHIFT] {gene}: called region {best['qmin']}-"
                  f"{best['qmax']} spans multiple reading frames "
                  f"({fs}); junction(s) near {fs_at}. The SEQUENCE is present "
                  f"and may well be confirmed by nr, but an intact full-length "
                  f"protein is not supported. Check reads at those coordinates "
                  f"before deciding pseudogene vs assembly indel.")

        org_m = re.search(r"\[([^\]]+)\]\s*$", manifest[best["acc"]]["product"])
        rows.append({
            "gene": gene, "call": call,
            "identity_pct": f"{ident:.1f}", "coverage_pct": f"{cov:.1f}",
            "contig": best["contig"],
            "coords": f"{best['qmin']}-{best['qmax']}",
            "best_accession": best["acc"],
            "organism": org_m.group(1) if org_m else "",
            "n_hits": len(ghits),
            # Empty for a clean locus. Non-empty means the called region is NOT
            # one open reading frame -- see frame_status() for how a frameshift
            # is told apart from two adjacent genes.
            "frame_status": fs,
            "frameshift_at": fs_at,
        })

    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t",
                           lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

    for r in rows:
        print(f"{r['gene']:8s} {r['call']:10s} id={r['identity_pct']:>6s}% "
              f"cov={r['coverage_pct']:>6s}% contig={r['contig']} "
              f"coords={r['coords']:>18s} ref={r['best_accession']} [{r['organism']}]")


if __name__ == "__main__":
    main()

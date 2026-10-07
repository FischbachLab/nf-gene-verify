#!/usr/bin/env python3
"""Merge round-1 gene calls with round-2 nr verification into the final table.

Reads gene_calls.tsv (round 1: coordinates from the staph protein screen),
nr_top_hits.tsv (round 2: nr evidence per locus), and the locus FASTA
headers (which encode the extracted window as contig:start-end) to map nr
hits back to genomic coordinates.

Revision rule (mechanical, documented):
  - Consider nr hits, in bitscore order, whose query span overlaps the round-1
    called region for the gene.
  - Classify each hit's title against that gene's ACCEPT and REJECT patterns.
    A hit that matches REJECT (a known paralog) can never confirm. A hit that
    matches ACCEPT confirms the gene, regardless of round-1
    identity/coverage — this neutralizes artifacts from RefSeq-restricted
    round-1 reference sets.
  - Confirmation may come from a lower-ranked hit, but only while its bitscore
    is at least --min-bitscore-frac of the top overlapping hit's. The rank used
    is recorded, so confirming on hit #2 is visible rather than silent.
  - Otherwise the gene is FLAGGED with the evidence shown; nothing is silently
    dropped.

Why patterns and not just the gene symbol: most nr deflines give a product
name, not a symbol. 'accessory gene regulator B' is agrB, 'Beta-lactamase' is
blaZ, 'serine-aspartate repeat protein F' is sdrF -- none contain the symbol
as a word. Matching the symbol alone falsely flagged 4 of 8 genes on a test
genome where every top hit was in fact correct.

Outputs:
  --out                       full table, one row per gene with all evidence
  gene_presence_matrix.tsv    minimal 'sample<TAB>gene<TAB>present' with 1/0,
                              written beside --out (override with
                              --matrix-out). Long format, so files from
                              separate genomes concatenate into one matrix.
                              Only 'present' is 1; 'flagged' is 0, because a
                              flagged gene is unverified and scoring it as
                              present is how a false positive reaches a figure.

Usage:
  finalize_calls.py --calls gene_calls.tsv --nr nr_top_hits.tsv \
                    --loci-dir nr_verification/ --out gene_presence_table_v2.csv
  finalize_calls.py ... --sample LM088 --patterns my_patterns.tsv
"""
import argparse
import csv
import os
import re
import sys

# ---------------------------------------------------------------------------
# Per-gene title patterns. Case-insensitive; applied to the nr defline.
#
#   accept -- this title identifies the gene
#   reject -- this title identifies a PARALOG and must never confirm
#
# A title matching neither is 'unknown': it does not confirm, but it also does
# not stop the scan, so an unambiguous hit a little further down can still
# confirm. sdrH is the worked example -- its top hit reads 'fibrinogen-binding
# protein', which is ambiguous (that is the usual name for SdrG/Fbe), while
# 'MSCRAMM-like protein SdrH' hits sit just below it at the same identity.
#
# Extend or override with --patterns (TSV: gene, accept, reject).
# ---------------------------------------------------------------------------
GENE_PATTERNS = {
    "agra": {
        "accept": r"\bagra\b|accessory gene regulator a\b|response regulator agra"
                  r"|quorum-?sensing response regulator",
        "reject": r"\bagr[bcd]\b|accessory gene regulator [bcd]\b",
    },
    "agrb": {
        "accept": r"\bagrb\b|accessory gene regulator b\b",
        "reject": r"\bagr[acd]\b|accessory gene regulator [acd]\b",
    },
    "blai": {
        "accept": r"\bblai\b|beta-?lactamase repressor|penicillinase repressor"
                  r"|blai/meci/copy",
        "reject": r"\bblaz\b|\bblar1?\b",
    },
    "blaz": {
        "accept": r"\bblaz\b|beta-?lactamase|penicillinase",
        "reject": r"\bblai\b|\bblar1?\b|repressor|metallo-?beta-?lactamase"
                  r"|class [cd] beta-?lactamase",
    },
    "fosb": {
        "accept": r"\bfosb\b|fosb family|metallothiol transferase"
                  r"|bacillithiol transferase|fosfomycin resistance|\byfcc\b",
        "reject": r"\bfos[acx]\b|fos[acx] family",
    },
    "sara": {
        "accept": r"\bsara\b|staphylococcal accessory regulator a\b"
                  r"|transcriptional regulator sara",
        "reject": r"\bsar[rstuvxz]\b|\brot\b",
    },
    "sdrf": {
        "accept": r"\bsdrf\b|serine-?aspartate repeat(-containing)? protein f\b",
        "reject": r"\bsdr[cdegh]\b|\bfbe\b|clumping factor|\bclf[ab]\b"
                  r"|b-like domain",
    },
    "sdrh": {
        "accept": r"\bsdrh\b|mscramm-?like protein sdrh|sdrh family",
        "reject": r"\bsdr[cdefg]\b|\bfbe\b|clumping factor|\bclf[ab]\b"
                  r"|b-like domain",
    },
}

# Kept for backward compatibility; folded into GENE_PATTERNS above.
SYNONYMS = {
    "fosb": ["fosb", "yfcc"],
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


def load_patterns(path):
    """Merge a TSV of 'gene<TAB>accept<TAB>reject' over the built-in table."""
    pats = {g: dict(v) for g, v in GENE_PATTERNS.items()}
    if not path:
        return pats
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            parts = line.split("\t")
            gene = parts[0].strip().lower()
            if not gene:
                continue
            pats[gene] = {
                "accept": parts[1].strip() if len(parts) > 1 else "",
                "reject": parts[2].strip() if len(parts) > 2 else "",
            }
    return pats


def classify_title(gene, title, pats):
    """'accept' | 'reject' | 'unknown' for this title, for this gene.

    reject is checked first: a paralog must never confirm, even when its
    defline happens to also contain an accept term (e.g. a title naming both
    SdrF and SdrG).
    """
    t = (title or "").lower()
    p = pats.get(gene.lower())

    if p is None:
        # No curated patterns: fall back to the bare symbol as a word. Weak,
        # and the caller notes it so the gap is visible rather than silent.
        return "accept" if re.search(rf"\b{re.escape(gene.lower())}\b", t) else "unknown"

    if p.get("reject") and re.search(p["reject"], t):
        return "reject"
    if p.get("accept") and re.search(p["accept"], t):
        return "accept"
    return "unknown"


def title_matches(gene, title, pats=None):
    """Backwards-compatible boolean wrapper."""
    return classify_title(gene, title, pats or GENE_PATTERNS) == "accept"


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--calls", required=True, help="gene_calls.tsv (round 1)")
    ap.add_argument("--nr", required=True, help="nr_top_hits.tsv (round 2)")
    ap.add_argument("--loci-dir", required=True,
                    help="dir with the *_locus.fasta files (for window coords)")
    ap.add_argument("--out", required=True,
                    help="final table path, e.g., gene_presence_table_v2.csv")
    ap.add_argument("--matrix-out", default="",
                    help="'sample<TAB>gene<TAB>present' TSV with 1/0 per gene "
                         "(default: gene_presence_matrix.tsv beside --out)")
    ap.add_argument("--sample", default="",
                    help="sample label for the matrix (default: the name of "
                         "the directory holding --out)")
    ap.add_argument("--patterns", default="",
                    help="TSV of 'gene<TAB>accept<TAB>reject' regexes, merged "
                         "over the built-in table")
    ap.add_argument("--min-bitscore-frac", type=float, default=0.80,
                    help="a lower-ranked hit may confirm only if its bitscore "
                         "is at least this fraction of the top overlapping "
                         "hit's (default: 0.80)")
    ap.add_argument("--expect-genus", default="Staphylococcus",
                    help="note it when the confirming hit's organism is "
                         "outside this genus ('' to disable)")
    args = ap.parse_args()

    pats = load_patterns(args.patterns)
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

        # All overlapping hits, bitscore-ordered, with their genomic spans
        overlapping = []
        if win and gs is not None:
            for rank, h in enumerate(hits, 1):
                g1 = win["win_start"] + int(h["query_start"]) - 1
                g2 = win["win_start"] + int(h["query_end"]) - 1
                if g2 >= gs and g1 <= ge:
                    overlapping.append(dict(h, g1=g1, g2=g2, _rank=rank))

        top = overlapping[0] if overlapping else None
        confirming = None
        rejected = 0
        notes = []

        if overlapping:
            top_bits = float(top["bit_score"])
            floor = top_bits * args.min_bitscore_frac
            for h in overlapping:
                verdict = classify_title(gene, h["title"], pats)
                if verdict == "reject":
                    rejected += 1
                    continue
                if verdict != "accept":
                    continue
                if float(h["bit_score"]) < floor:
                    # A matching hit exists but is much weaker than the best
                    # overlapping hit; report rather than confirm on it.
                    notes.append(
                        f"best matching hit {h['accession']} is rank "
                        f"{h['_rank']} at {h['bit_score']} bits, below "
                        f"{args.min_bitscore_frac:.0%} of top ({top['bit_score']})")
                    break
                confirming = h
                break

        if gene.lower() not in pats:
            notes.append(f"no curated title patterns for '{gene}' -- matched on "
                         f"the bare symbol only")

        if top is None:
            nr_status = "flagged: no nr hit overlaps the called region"
            final = "flagged"
        elif confirming is not None:
            cov = confirming.get("hit_cov_pct") or confirming.get("query_cov_pct") or ""
            nr_status = (f"confirmed by nr: {confirming['accession']} "
                         f"({confirming['identity_pct']}% id, {cov}% cov)")
            if confirming["_rank"] != top["_rank"]:
                nr_status += f" [hit rank {confirming['_rank']}, not the top hit]"
            genus = (confirming.get("organism") or "").split()[:1]
            if args.expect_genus and genus and genus[0] != args.expect_genus:
                notes.append(f"confirming organism is {confirming['organism']}, "
                             f"outside {args.expect_genus}")
            if rejected:
                notes.append(f"{rejected} overlapping hit(s) rejected as paralogs")
            final = "present"
        else:
            nr_status = (f"flagged: no overlapping nr hit identifies {gene}; "
                         f"top hit is {top['accession']} "
                         f"({top['title'][:60]})")
            if rejected:
                notes.append(f"{rejected} overlapping hit(s) matched a paralog "
                             f"pattern")
            notes.append(f"{len(overlapping)} overlapping hit(s) examined")
            final = "flagged"

        if notes:
            nr_status += " | " + "; ".join(notes)

        # Report the hit the decision rests on, not merely the highest-scoring
        top = confirming or top

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
            "nr_hit_rank": top.get("_rank", "") if top else "",
            "nr_verification": nr_status,
            "final_call": final,
        })

    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    # ── Minimal presence matrix ─────────────────────────────────────────────
    # sample / gene / present, long format. Long rather than wide so files
    # from separate genomes concatenate directly into one matrix; pivot when a
    # wide table is wanted.
    #
    # Only 'present' scores 1. 'flagged' is deliberately 0, because flagged
    # means NOT VERIFIED, and scoring an unverified gene as present is how a
    # false positive reaches a figure. The reason for any 0 is in --out.
    out_dir = os.path.dirname(os.path.abspath(args.out)) or "."
    sample = args.sample or os.path.basename(out_dir)
    matrix_path = args.matrix_out or os.path.join(out_dir,
                                                  "gene_presence_matrix.tsv")
    with open(matrix_path, "w", newline="") as fh:
        fh.write("sample\tgene\tpresent\n")
        for r in sorted(rows, key=lambda x: x["gene"]):
            fh.write(f"{sample}\t{r['gene']}\t"
                     f"{1 if r['final_call'] == 'present' else 0}\n")

    for r in rows:
        print(f"{r['gene']:8s} round1={r['round1_call']:10s} "
              f"final={r['final_call']:8s} {r['nr_verification']}")

    n_present = sum(1 for r in rows if r["final_call"] == "present")
    print(f"\n{n_present}/{len(rows)} gene(s) present   (sample: {sample})")
    print(f"written: {args.out}")
    print(f"written: {matrix_path}")


if __name__ == "__main__":
    main()

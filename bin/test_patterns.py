#!/usr/bin/env python3
"""Consistency checks for PRODUCT_FILTERS (round 1) and GENE_PATTERNS (round 2).

Run after editing either table:

    python3 bin/test_patterns.py [--manifest path ...]

Four checks, each guarding a bug that actually occurred:

1. KEY PARITY. bhp was fixed in PRODUCT_FILTERS and not in GENE_PATTERNS, so
   the build step and the verification step disagreed about which titles name
   the gene. Only 'sccmec' may be round-2-only: it is a predefined cassette
   with no build-time filter.

2. CROSS-TABLE AGREEMENT over every real NCBI title in the manifests. Catches
   patterns that drift apart in CONTENT rather than in key -- agrA/agrB had
   '(protein )?' in one table and not the other, so 'accessory gene regulator
   protein A' was accepted at build time and rejected in round 2.

3. NO PANEL MAY EMPTY. An over-strict accept pattern rejects every record; the
   length summary then raised IndexError and killed the whole build. bhp did
   this ('biofilm-associated protein' vs the real 'cell wall associated biofilm
   protein'), and opp3A did it again when a 'ATP[- ]?binding' deny matched the
   family name 'ATP binding cassette' in the substrate-binding subunit's title.

4. DOMAIN_ACCEPT IS READ, AND GATED. 'reject_extra' was a key I invented that
   classify_title() never read, so 31 reject lists were silently disabled; the
   same mistake in 'domain_accept' would silently confirm bhp on any Aap or
   SasG repeat. This checks the key is honoured, that the length and coverage
   floors both bite, and that 'reject' still wins over it.

5. KNOWN TRAPS. Homonyms and paralogs that produced, or nearly produced,
   reportable false positives: mecA (ClpC adaptor vs PBP2a), msrA (methionine
   sulfoxide reductase vs ABC-F macrolide resistance), dfr (native
   trimethoprim-SENSITIVE reductase), tetK vs ribosomal-protection tet genes,
   essD vs native nucleases, aae vs atlE.
"""

import argparse
import csv
import importlib.util
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))

# (gene, product title, should_be_accepted)
TRAPS = [
    ("mecA",  "adaptor protein MecA [Staphylococcus epidermidis]", False),
    ("mecA",  "PBP2a family beta-lactam-resistant peptidoglycan transpeptidase MecA", True),
    ("msrA",  "peptide-methionine (S)-S-oxide reductase MsrA", False),
    ("msrA",  "ABC-F type ribosomal protection protein Msr(A)", True),
    ("opp3A", "oligopeptide ABC superfamily ATP binding cassette transporter, binding protein", True),
    ("opp3A", "oligopeptide ABC transporter permease", False),
    ("opp3B", "oligopeptide ABC transporter permease", True),
    ("arcB",  "ornithine carbamoyltransferase, catabolic", True),
    ("arcB",  "arginine deiminase", False),
    ("dfrG",  "trimethoprim-sensitive dihydrofolate reductase DfrC", False),
    ("tetK",  "tetracycline efflux MFS transporter", True),
    ("tetK",  "tetracycline resistance protein TetM", False),
    ("essD",  "type VII secretion systems protein EssD", True),
    ("essD",  "DNA/RNA non-specific endonuclease", False),
    ("aae",   "autolysin/adhesin Aae", True),
    ("aae",   "bifunctional autolysin", False),
    ("bhp",   "cell wall associated biofilm protein", True),
    ("bhp",   "adhesin [Staphylococcus epidermidis]", False),
    # round-2 false positives caught on the first --verify_nr run
    ("opp3D", "iron ABC transporter ATP-binding protein", False),
    ("opp3E", "iron ABC transporter ATP-binding protein", False),
    ("opp3D", "nickel (Ni2+)/peptide ABC superfamily ATP binding cassette transporter, ABC protein", True),
    # round-2 false negatives from the same run: nr names these by FUNCTION
    ("atlE",  "N-acetylmuramoyl-L-alanine amidase [Staphylococcus epidermidis]", True),
    ("opp3A", "peptide ABC transporter substrate-binding protein [Staphylococcus]", True),
    # aap: nr names this family by DOMAIN only -- 0 of 100 hits said
    # 'accumulation-associated protein'. Family-level confirmation by design.
    ("aap",   "E domain-containing protein [Staphylococcus epidermidis]", True),
    ("aap",   "G5 domain-containing protein [Staphylococcus epidermidis]", True),
    ("aap",   "truncated SdrF protein", False),
    # Missing word boundary: 'E domain-containing' matched the trailing 'e' of
    # the preceding word, so three generic Gram-positive surface annotations
    # confirmed aap. In v5 that overturned five homology_only calls (round 1:
    # 16.9-26.4% identity) into 'present'.
    ("aap",   "Ig-like domain-containing protein [Staphylococcus epidermidis]", False),
    ("aap",   "KxYKxGKxW signal peptide domain-containing protein", False),
    ("aap",   "Rib/alpha-like domain-containing protein", False),
    ("aap",   "E domain-containing protein [Staphylococcus epidermidis]", True),
    ("aap",   "G5 domain-containing protein [Staphylococcus epidermidis]", True),
    # essD: 'DNA/RNA non-specific endonuclease' is deliberately NOT accepted --
    # 72 of 100 hits carry it and native nucleases do too, so accepting it would
    # confirm essD on any nuclease. essD stays flagged on LM088, correctly.
    ("essD",  "DNA/RNA non-specific endonuclease [Staphylococcus epidermidis]", False),
    ("essD",  "type VII secretion systems protein EssD", True),
    # MULTISPECIES synonym leak: classify_title() reads only the PRIMARY title
    ("opp3D", "MULTISPECIES: ABC transporter ATP-binding protein [Staphylococcus] >KTF27926.1 peptide ABC transporter ATP-binding protein", False),
    ("icaD",  "MULTISPECIES: intracellular adhesion protein IcaD [Staphylococcus] >P69519.1 Rec", True),
    # dfr family: the accept must name the gene, not the phenotype
    ("dfrG",  "MULTISPECIES: trimethoprim-resistant dihydrofolate reductase DfrS1 [Bacteria]", False),
    ("dfrG",  "trimethoprim-resistant dihydrofolate reductase DfrG [Staphylococcus aureus]", True),
    ("dfrC",  "trimethoprim-resistant dihydrofolate reductase DfrG", False),
    ("blaZ",  "Beta-lactamase [Staphylococcus aureus]", True),
    ("blaZ",  "Beta-lactamase regulator protein", False),
    # sesI: generic Gram-positive surface vocabulary must NOT confirm it.
    # These three titles overturned five homology_only calls (27.7-38.2%
    # identity) in v5. nr never names SesI, so only the identifiers accept.
    ("sesI",  "cell surface protein [Staphylococcus epidermidis]", False),
    ("sesI",  "LPXTG cell wall anchor domain-containing protein", False),
    ("sesI",  "MULTISPECIES: LPXTG cell wall anchor domain-containing protein [Staphylococcus]", False),
    ("sesI",  "hypothetical protein [Staphylococcus epidermidis]", False),
    ("sesI",  "SesI protein [Staphylococcus epidermidis]", True),
    ("sdrF",  "serine-aspartate repeat protein F", True),
    ("sdrF",  "MSCRAMM family adhesin SdrC/SdrD/SdrE-like", False),
]

ROUND2_ONLY_OK = {"sccmec"}

# Genes where the two tables are INTENTIONALLY asymmetric, and why.
#
# These are retrieved by accession/locus tag plus a seed, because their
# records carry only a general annotation. Build-time filtering must
# therefore be deny-only -- an accept gate would reject the very record the
# query exists to fetch (this emptied bhp, sesI and opp3A in turn).
#
# Round 2 must be the opposite: STRICTER. Admitting the same generic titles
# there let an unrelated ABC transporter confirm opp3D at 99.6% and an IRON
# transporter confirm opp3E at 99.2%, on loci 435 kb and 161 kb from the
# opp3 cluster, overturning two correct homology_only calls.
#
# The cost is explicit: a REAL opp3D/opp3E cannot be confirmed by title
# either, so it comes back flagged rather than present. Round-1 identity
# against the seeded reference is the evidence that counts for these
# (~100% for a genuine copy, ~29% for background).
ASYMMETRIC_OK = {"opp3d", "opp3e", "sesi", "aap"}
#
# aap is asymmetric for the opposite reason to the others. Its QUERY is
# aap[Gene Name], so every retrieved record already carries the symbol and a
# symbol-only accept is right at build time (panel: n=20, 1243-2541 aa). But nr
# annotates this family by DOMAIN -- 0 of 100 hits on the LM088 locus said
# 'accumulation-associated protein' -- so round 2 must accept 'E domain-
# containing', 'G5 domain-containing' and 'pectate lyase-like adhesive domain'.
# Widening the build filter to match would serve no purpose and would only
# loosen panel construction.




def load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--manifest", action="append", default=[],
                    help="manifest.csv with real NCBI titles (repeatable). "
                         "Checks 2 and 3 are skipped without one.")
    args = ap.parse_args()

    bp = load(os.path.join(HERE, "build_protein_db.py"), "bp")
    fc = load(os.path.join(HERE, "finalize_calls.py"), "fc")
    PF, GP = bp.PRODUCT_FILTERS, fc.GENE_PATTERNS
    fails = []

    # 1 ── key parity
    build_only = sorted(set(PF) - set(GP))
    round2_only = sorted(set(GP) - set(PF) - ROUND2_ONLY_OK)
    print(f"[1] key parity: PRODUCT_FILTERS={len(PF)} GENE_PATTERNS={len(GP)}")
    if build_only:
        fails.append(f"curated at build time but not in round 2: {build_only}")
        print(f"    FAIL build-only: {build_only}")
    if round2_only:
        fails.append(f"curated in round 2 but not at build time: {round2_only}")
        print(f"    FAIL round2-only: {round2_only}")
    if not build_only and not round2_only:
        print("    ok")

    # 2+3 ── against real titles
    titles, bygene = set(), {}
    for m in args.manifest:
        if not os.path.exists(m):
            print(f"    [skip] no such manifest: {m}")
            continue
        for r in csv.DictReader(open(m)):
            titles.add((r["gene"], r["product"]))
            bygene.setdefault(r["gene"], []).append(r["product"])

    # Seeded references carry a synthetic product ('user-supplied reference')
    # written by --extra-fasta, not an NCBI title. It cannot appear in nr, so
    # comparing the two tables on it is meaningless -- and it produces a
    # guaranteed mismatch for any gene whose build filter is deny-only
    # (opp3D/opp3E/sesI, where an accept gate would reject the very record the
    # accession query retrieves). Excluded rather than papered over by adding
    # the placeholder to the round-2 accept list, which would let any
    # unannotated record confirm those genes.
    SEED_PRODUCT = "user-supplied reference"
    titles = {(g, p) for g, p in titles if p.strip().lower() != SEED_PRODUCT}



    if titles:
        mis = [(g, p) for g, p in titles
               if g.lower() in PF and g.lower() in GP
               and g.lower() not in ASYMMETRIC_OK
               and bp.product_ok(g, p) != (fc.classify_title(g, p, GP) == "accept")]
        print(f"[2] cross-table agreement over {len(titles)} real titles: "
              f"{len(mis)} disagreement(s)")
        for g, p in mis[:10]:
            print(f"    FAIL {g}: {p[:66]}")
        if mis:
            fails.append(f"{len(mis)} title(s) classified differently by the two tables")

        zero = [g for g, ps in bygene.items()
                if g.lower() in PF and not any(bp.product_ok(g, p) for p in ps)]
        print(f"[3] panel survival: {len(zero)} gene(s) filtered to zero")
        for g in zero:
            print(f"    FAIL {g}: every record rejected -- check accept/deny "
                  f"against the real titles")
        if zero:
            fails.append(f"zero-reference genes: {zero}")
    else:
        print("[2] skipped (no manifest)")
        print("[3] skipped (no manifest)")

    # 4 ── domain_accept: honoured, and gated on both floors
    #
    # Real numbers from the 29-genome bhp re-run:
    #   NIH05003  the locus to rescue -- 99.9% id over a complete 2402 aa
    #             reference, 0 of 100 nr hits naming biofilm/Bap/Bhp
    #   NIHLM087  a shared repeat: full-length Bhp subject, 8.3% covered
    #   NIHLM018  a 'hypothetical protein, partial' subject at 100% coverage,
    #             which is why min_hit_len exists alongside min_hit_cov
    IG = "putative Ig domain-containing protein [Staphylococcus epidermidis]"
    DOMAIN_CASES = [
        # gene, title,                     hit_len, hit_cov, want, why
        ("bhp", IG,                           2402,   100.0, True,
         "NIH05003: full-length subject, fully covered"),
        ("bhp", IG,                           2491,   100.0, True,
         "NIH04003/06004/08001: the same locus by another route"),
        ("bhp", IG,                           2500,     8.3, False,
         "NIHLM087: shared repeat on a full-length subject"),
        ("bhp", IG,                           2500,    40.0, False,
         "below the 80% subject-coverage floor"),
        ("bhp", "hypothetical protein, partial [S. epidermidis]",
                                               217,   100.0, False,
         "NIHLM018: fully covered, but the subject is a fragment"),
        ("bhp", "accumulation-associated protein Aap [S. epidermidis]",
                                              2541,   100.0, False,
         "named paralog: 'reject' must win over domain_accept"),
        ("bhp", "cell wall associated biofilm protein",
                                              2402,    10.0, True,
         "named outright: plain accept, no coverage floor applies"),
        ("aap", IG,                           2541,   100.0, False,
         "aap has no domain_accept spec: must not inherit bhp's"),
    ]
    dbad = 0
    print(f"[4] domain_accept: {len(DOMAIN_CASES)} cases")
    if not any("domain_accept" in v for v in GP.values()):
        fails.append("no gene defines domain_accept -- did the key get dropped?")
        print("    FAIL: domain_accept absent from every GENE_PATTERNS entry")
    for gene, title, hit_len, hit_cov, want, why in DOMAIN_CASES:
        hit = {"title": title, "hit_len": str(hit_len), "hit_cov_pct": str(hit_cov)}
        verdict = fc.classify_title(gene, title, GP)
        # Mirror the caller: reject wins; a plain accept needs no gate.
        if verdict == "reject":
            got = False
        elif verdict == "accept":
            got = True
        else:
            got = fc.domain_accept(gene, hit, GP)
        if got != want:
            dbad += 1
            print(f"    FAIL {gene:5} got={got} want={want}  {why}")
    if dbad:
        fails.append(f"{dbad} domain_accept case(s) misclassified")
    else:
        print("    ok")

    # 5 ── known traps
    bad = 0
    print(f"[5] known traps: {len(TRAPS)} cases")
    for gene, title, want in TRAPS:
        if gene.lower() not in PF or gene.lower() not in GP:
            print(f"    SKIP {gene}: not curated in both tables")
            continue
        a = bp.product_ok(gene, title)
        b = fc.classify_title(gene, title, GP) == "accept"
        # For ASYMMETRIC_OK genes the build filter is deliberately deny-only,
        # so only the round-2 verdict is asserted here.
        a_ok = True if gene.lower() in ASYMMETRIC_OK else (a == want)
        if not a_ok or b != want:
            bad += 1
            print(f"    FAIL {gene:6} build={a} round2={b} want={want}  {title[:50]}")
    if bad:
        fails.append(f"{bad} trap case(s) misclassified")
    else:
        print("    ok")

    print()
    if fails:
        print("FAILED:")
        for f in fails:
            print(f"  - {f}")
        return 1
    print("ALL CHECKS PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())

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
    # '(protein )?' matters: NCBI writes both 'accessory gene regulator A' and
    # 'accessory gene regulator protein A'. PRODUCT_FILTERS had the optional
    # group and this table did not, so the same title was accepted at build time
    # and rejected in round 2 -- caught by the cross-table consistency test over
    # real manifest titles.
    "agra": {
        "accept": r"\bagra\b|accessory gene regulator (protein )?a\b"
                  r"|response regulator agra|quorum-?sensing response regulator",
        "reject": r"\bagr[bcd]\b|accessory gene regulator (protein )?[bcd]\b",
    },
    "agrb": {
        "accept": r"\bagrb\b|accessory gene regulator (protein )?b\b",
        "reject": r"\bagr[acd]\b|accessory gene regulator (protein )?[acd]\b",
    },
    "blai": {
        "accept": r"\bblai\b|beta-?lactamase repressor|penicillinase repressor"
                  r"|blai/meci/copy",
        "reject": r"\bblaz\b|\bblar1?\b",
    },
    "blaz": {
        "accept": r"\bblaz\b|beta-?lactamase|penicillinase",
        "reject": r"\bblai\b|\bblar1?\b|repressor|regulator|regulatory"
                  r"|sensor-?transducer|metallo-?beta-?lactamase"
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
    # ── mecA: the adaptor/PBP2a homonym ─────────────────────────────────────
    # 'adaptor protein MecA' contains the bare symbol, so symbol-only matching
    # ACCEPTS it and reports methicillin resistance in a sensitive genome. The
    # reject list below is the round-2 backstop for that; PRODUCT_FILTERS in
    # build_protein_db.py keeps it out of the panel in the first place.
    "meca": {
        "accept": r"\bmeca\b|pbp2a|pbp2'|penicillin-?binding protein 2a"
                  r"|beta-?lactam-?resistant peptidoglycan transpeptidase"
                  r"|methicillin resistance",
        "reject": r"adaptor|adapter|\bclpc\b|competence|sporulation"
                  r"|\bmec[ri]1?\b",
    },
    "icad": {
        "accept": r"\bicad\b|intercellular adhesion protein d\b"
                  r"|poly-?beta-?1,6-?n-?acetyl|pia synthesis protein icad",
        "reject": r"\bica[abcr]\b|intercellular adhesion protein [abc]\b",
    },
    "icar": {
        "accept": r"\bicar\b|ica operon|icaadbc|biofilm operon",
        "reject": r"\bica[abcd]\b|intercellular adhesion protein [abcd]\b",
    },
    "is256": {
        "accept": r"\bis256\b|is256[- ]?(like|family)",
        "reject": r"\bis(?!256)\d{1,4}\b",
    },
    "aaca": {
        # NCBI titles this enzyme both ways: spelled out ('bifunctional
        # aminoglycoside N-acetyltransferase...') and abbreviated
        # ('bifunctional AAC/APH'). Both forms must be accepted.
        "accept": r"\baaca\b|aaca-?aphd|aac\(6'\)|aph\(2''?\)|aac/aph"
                  r"|bifunctional\s+(aac|aph)"
                  r"|bifunctional.*(acetyltransferase|phosphotransferase)"
                  r"|aminoglycoside.*\(6'\).*acetyltransferase",
        "reject": r"aac\(3\)|aac\(2'\)|\baada\b|\bant\(|streptomycin"
                  r"|streptothricin",
    },
    "fdh": {
        "accept": r"\bfdh[a-c]?\b|formate dehydrogenase",
        "reject": r"\bfdh[de]\b|accessory|alcohol dehydrogenase"
                  r"|formaldehyde dehydrogenase",
    },
    "bhp": {
        # a bare 'adhesin' is NOT accepted -- see PRODUCT_FILTERS for why.
        # 'cell wall associated biofilm protein' is the title NCBI actually
        # uses for Bhp; must stay in step with PRODUCT_FILTERS['bhp'].
        "accept": r"\bbhp\b|bap[- ]?(homolog\w*|family)"
                  r"|(cell wall[- ]?associated )?biofilm[- ]?(associated )?protein",
        "reject": r"\bsdr[cdefgh]\b|\baap\b|\bembp\b|clumping factor"
                  r"|fibronectin|elastin|serine-?aspartate",
        # nr annotates this family by DOMAIN and, for most records, never by
        # name: 93 of 100 hits on NIH05003's locus say 'putative Ig
        # domain-containing protein' and 0 say biofilm/Bap/Bhp. Accepted ONLY
        # on a full-length subject -- see domain_accept() for the measurement
        # that fixes these two numbers. Deliberately NOT mirrored in
        # PRODUCT_FILTERS: the build-time query retrieves bhp BY TITLE, so
        # widening that accept would only admit Aap and SasG to the panel.
        "domain_accept": {
            "pattern": r"ig domain-containing|g5 domain-containing"
                       r"|pectate lyase-?like adhesive domain",
            "min_hit_len": 1500,   # Bhp is ~2,500 aa; excludes fragment records
            "min_hit_cov": 80.0,   # % of the SUBJECT covered
        },
    },
    "ebh": {
        "accept": r"\bebh\b|ebh-?like|hyperosmolarity resistance",
        "reject": r"\bembp\b",
    },
    # ── SCCmec: a cassette, not a gene ──────────────────────────────────────
    # MecI/MecR1 are homologues of BlaI/BlaR1 at ~60-65% identity, so a cassette
    # panel cross-matches the bla operon: on LM088-RCB001 the SCCmec 'locus'
    # landed at 903916-906025, 107 bp downstream of blaZ, with BlaR1 as the top
    # hit at 65.7% round-1 identity. Rejecting the bla operon is what stops that
    # becoming a confident 'SCCmec present'.
    "sccmec": {
        "accept": r"\bmeca\b|pbp2a|\bmec[ri]1?\b|\bccr[abc]\b"
                  r"|cassette chromosome recombinase|staphylococcal cassette",
        "reject": r"\bbla[izr]1?\b|penicillinase|beta-?lactamase repressor"
                  r"|regulatory protein blar|adaptor|adapter|\bclpc\b",
    },
    # ════════════════════════════════════════════════════════════════════════
    # Written from the product titles actually returned by NCBI on the 40-gene
    # run (TEST_40_genes/protein_db/manifest.csv), not from memory. Guessing a
    # title emptied the bhp panel and killed a build after five genes had
    # succeeded, so every accept term below appears verbatim in that manifest
    # except where marked UNVERIFIED.
    # ════════════════════════════════════════════════════════════════════════

    # ── ACME arc cluster ────────────────────────────────────────────────────
    # Each denies the other three: all four are arginine-pathway enzymes and a
    # shared pathway term would let one confirm another's locus.
    "arca": {"accept": r"\barcA\b|arginine deiminase",
             "reject":   r"\barc[BCD]\b|ornithine carbamoyltransferase"
                       r"|carbamate kinase|antiporter"},
    # 'catabolic' matters: ArcB is the CATABOLIC ornithine carbamoyltransferase,
    # and the anabolic isoenzyme ArgF is a native arginine-biosynthesis gene
    # present in every genome.
    "arcb": {"accept": r"\barcB\b|ornithine carbamoyltransferase",
             "reject":   r"\barc[ACD]\b|\bargF\b|anabolic|arginine deiminase"
                       r"|carbamate kinase"},
    "arcc": {"accept": r"\barcC\b|carbamate kinase",
             "reject":   r"\barc[ABD]\b|arginine deiminase"
                       r"|ornithine carbamoyltransferase"},
    "arcd": {"accept": r"\barcD\b|arginine[-/ ]ornithine antiporter",
             "reject":   r"\barc[ABC]\b|arginine deiminase|carbamate kinase"},
    "argr": {"accept": r"\bargR\b|\bahrC\b|AhrC/ArgR|arginine repressor",
             "reject":   r"\barc[ABCD]\b"},

    # ── ACME opp3 cluster ───────────────────────────────────────────────────
    # NOTE: NCBI gives opp3B and opp3C the IDENTICAL title 'oligopeptide ABC
    # transporter permease', so no pattern can separate them. Their panels
    # overlap by design. That is acceptable because ACME typing needs the
    # CLUSTER detected, not each subunit resolved -- acme_type.py requires
    # --min-cluster members, so a permease matching both still yields one
    # cluster rather than two spurious ones.
    "opp3a": {"accept": r"\bopp-?3A\b|\boppA\b"
                        r"|(oligo)?peptide ABC.*substrate-?binding"
                        r"|oligopeptide.*binding,? protein"
                        r"|oligopeptide.*substrate-?binding",
              "reject":   r"permease|ATP[- ]?binding protein"},
    "opp3b": {"accept": r"\bopp-?3B\b|\boppB\b|oligopeptide ABC transporter permease",
              "reject":   r"ATP[- ]?binding|substrate-?binding"},
    "opp3c": {"accept": r"\bopp-?3C\b|\boppC\b|oligopeptide ABC transporter permease",
              "reject":   r"ATP[- ]?binding|substrate-?binding"},
    # TIGHTENED after a false positive. The previous accept list admitted
    # 'ABC transporter ATP-binding' and 'ATP-binding cassette
    # domain-containing', which matches essentially every ABC transporter in a
    # genome. Round 2 then "confirmed" opp3D at 99% on a background locus
    # 435 kb from the opp3 cluster -- for opp3E the confirming hit was an
    # IRON transporter -- overturning a correct homology_only call.
    #
    # The unavoidable consequence: the TRUE opp3D is also generically titled in
    # nr, so a real copy can no longer be confirmed by title either. That is the
    # honest position -- these two subunits are UNVERIFIABLE in round 2, and a
    # real one will come back flagged rather than present. Round 1 identity
    # against the seeded ASM79718.1 reference is the evidence that counts
    # here (a genuine copy scores ~100%; background scores ~29%).
    "opp3d": {"accept": r"\bopp-?3D\b|ASM79718|WP_002490823"
                       r"|(oligo)?peptide ABC.*ATP[- ]?binding"
                       r"|nickel.*peptide ABC",
              # Other substrate specificities: an iron ABC transporter is
              # what wrongly confirmed opp3E. Merged into "reject" -- a
              # "reject_extra" key would be silently ignored by
              # classify_title(), which reads only accept/reject.
              "reject":   r"permease|substrate-?binding"
                          r"|iron|zinc|manganese|molybdate|sulfate"
                          r"|phosphate|spermidine|polyamine|sugar"},
    # TIGHTENED after a false positive. The previous accept list admitted
    # 'ABC transporter ATP-binding' and 'ATP-binding cassette
    # domain-containing', which matches essentially every ABC transporter in a
    # genome. Round 2 then "confirmed" opp3E at 99% on a background locus
    # 161 kb from the opp3 cluster -- for opp3E the confirming hit was an
    # IRON transporter -- overturning a correct homology_only call.
    #
    # The unavoidable consequence: the TRUE opp3E is also generically titled in
    # nr, so a real copy can no longer be confirmed by title either. That is the
    # honest position -- these two subunits are UNVERIFIABLE in round 2, and a
    # real one will come back flagged rather than present. Round 1 identity
    # against the seeded ASM79719.1 reference is the evidence that counts
    # here (a genuine copy scores ~100%; background scores ~29%).
    "opp3e": {"accept": r"\bopp-?3E\b|ASM79719|WP_002490830"
                       r"|(oligo)?peptide ABC.*ATP[- ]?binding"
                       r"|nickel.*peptide ABC",
              # Other substrate specificities: an iron ABC transporter is
              # what wrongly confirmed opp3E. Merged into "reject" -- a
              # "reject_extra" key would be silently ignored by
              # classify_title(), which reads only accept/reject.
              "reject":   r"permease|substrate-?binding"
                          r"|iron|zinc|manganese|molybdate|sulfate"
                          r"|phosphate|spermidine|polyamine|sugar"},

    # ── ACME kdp operon (types IV, V) ───────────────────────────────────────
    # UNVERIFIED: kdp was not in the 40-gene run, so these accept terms are the
    # standard NCBI product names but have not been checked against real output.
    # The per-gene counts in the next build log will show whether they worked.
    # Like arc, kdp is NATIVE in staphylococci -- presence is not ACME evidence,
    # which is why acme_type.py requires co-location.
    "kdpa": {"accept": r"\bkdpA\b|(potassium|K\(\+\))-?transporting ATPase.*subunit A",
             "reject":   r"\bkdp[BCDEF]\b|subunit [BCDF]\b"},
    "kdpb": {"accept": r"\bkdpB\b|(potassium|K\(\+\))-?transporting ATPase.*subunit B",
             "reject":   r"\bkdp[ACDEF]\b|subunit [ACDF]\b"},
    "kdpc": {"accept": r"\bkdpC\b|(potassium|K\(\+\))-?transporting ATPase.*subunit C",
             "reject":   r"\bkdp[ABDEF]\b|subunit [ABDF]\b"},
    "kdpd": {"accept": r"\bkdpD\b|sensor histidine kinase KdpD"
                       r"|potassium-?transporting ATPase KdpD",
             "reject":   r"\bkdp[ABCE]\b|response regulator"},
    "kdpe": {"accept": r"\bkdpE\b|response regulator KdpE|KdpE family",
             "reject":   r"\bkdp[ABCD]\b|histidine kinase"},

    # ── acquired resistance ─────────────────────────────────────────────────
    "aadd": {"accept": r"\baadD\b|ANT\(4'\)"
                       r"|kanamycin.*(nucleotidyl|adenylyl)transferase"
                       r"|aminoglycoside adenyl",
             "reject":   r"\baad[ABE]\b|ANT\(3|ANT\(6|ANT\(9|streptomycin"
                       r"|spectinomycin"},
    "mphc": {"accept": r"\bmph\(?C\)?\b|macrolide 2'-phosphotransferase",
             "reject":   r"\bmph\(?[ABDE]\)?\b"},
    # tetL is the close homologue; tetM/tetO are ribosomal protection proteins,
    # a different mechanism that should never confirm an efflux gene.
    "tetk": {"accept": r"\btetK\b|tetracycline efflux MFS"
                       r"|tetracycline resistance protein",
             "reject":   r"\btet[ABCDEGHLMOSTVWZ]\b|ribosomal protection"},
    # norB/norC/mdeA/mepA are the other staphylococcal MFS/MATE pumps
    "nora": {"accept": r"\bnorA\b|multidrug efflux MFS transporter NorA",
             "reject":   r"\bnor[BCDG]\b|\bmdeA\b|\bmepA\b|\bsepA\b|\btet[KL]\b"},
    "uhpt": {"accept": r"\buhpT\b|hexose-6-phosphate",
             "reject":   r"\buhp[ABC]\b|\bglpT\b"},

    # ── adhesins / biofilm / secretion ──────────────────────────────────────
    # Aae and AtlE are BOTH autolysins. 'autolysin/adhesin Aae' vs 'bifunctional
    # autolysin' is the only thing separating them, so each denies the other's
    # exact title.
    "aae":  {"accept": r"\baae\b|autolysin/adhesin",
             "reject":   r"\batlE?\b|bifunctional autolysin"},
    # NOTE: 'bifunctional autolysin' is also the title of S. aureus Atl, which
    # no pattern can distinguish from AtlE. Species preference in the panel is
    # what keeps this honest, not the title.
    # AtlE is a BIFUNCTIONAL autolysin: an N-acetylmuramoyl-L-alanine amidase
    # plus a glucosaminidase. nr names it by DOMAIN, not by 'autolysin', so the
    # round-1 call (99.3% over 100% of AAW53968, S. epidermidis RP62A AtlE) was
    # flagged by round 2 purely because the accept list lacked the domain names.
    "atle": {"accept": r"\batlE\b|bifunctional autolysin"
                       r"|N-acetylmuramoyl-L-alanine amidase|glucosaminidase",
             "reject":   r"\baae\b|autolysin/adhesin|\batlA\b"},
    # Aap shares B-repeats with Bhp, Embp and the Sdr family, and -seg no makes
    # that cross-matching more likely, so every one of them is denied.
    # nr does not annotate this family by symbol: of 100 hits on the LM088
    # locus, ZERO said 'accumulation-associated protein'. They are named by
    # DOMAIN -- 66 'E domain-containing protein', 22 'pectate lyase-like
    # adhesive domain', 5 'G5 domain-containing protein' -- which together are
    # precisely the Aap architecture: a pectate lyase-like A domain followed by
    # G5 repeats separated by E spacers.
    #
    # LIMITATION, stated plainly: these domain names are FAMILY-level, shared
    # with SasG, Pls and Bhp. A confirmation here means "a G5-E repeat adhesin
    # of the Aap/SasG family", not "Aap specifically". The deny list cannot help
    # because the domain titles never mention the paralogues by name. Round-1
    # identity against the Aap panel is what narrows it to Aap.
    # '\bE domain-containing', NOT 'E domain-containing'. Without the word
    # boundary the alternative matches the trailing 'e' of the PRECEDING word,
    # so it accepted 'Ig-likE domain-containing', 'signal peptidE
    # domain-containing' and 'Rib/alpha-likE domain-containing'. In the
    # 29-genome v5 run that confirmed aap in five genomes whose round 1 had
    # called it homology_only at 16.9-26.4% identity -- i.e. round 1 said this
    # is NOT aap and round 2 overturned it on a generic Gram-positive surface
    # motif. A one-letter alternative needs an anchor on both sides.
    "aap":  {"accept": r"\baap\b|accumulation[- ]?associated protein"
                       r"|\bE domain-containing|\bG5 domain-containing"
                       r"|pectate lyase-?like adhesive domain",
             "reject":   r"\bbhp\b|\bsdr[CDEFGH]\b|\bembp\b|biofilm-?associated"
                       r"|clumping factor"},
    "agrc": {"accept": r"\bagrC\b|quorum-?sensing sensor histidine kinase"
                       r"|accessory gene regulator (protein )?C\b",
             "reject":   r"\bagr[ABD]\b|accessory gene regulator (protein )?[ABD]\b"},
    # AgrD is a ~46 aa propeptide; most records carry no symbol, only the
    # product name or the AIP type.
    "agrd": {"accept": r"\bagrD\b|cyclic lactone autoinducer|autoinducing peptide"
                       r"|AIP-I{1,3}V?|accessory gene regulator (protein )?D\b",
             "reject":   r"\bagr[ABC]\b|accessory gene regulator (protein )?[ABC]\b"},
    "icaa": {"accept": r"\bicaA\b|intercellular adhesion protein A\b"
                       r"|N-acetyl-D-glucosamine synthase",
             "reject":   r"\bica[BCDR]\b|intercellular adhesion protein [BCD]\b"
                       r"|N-deacetylase|biosynthesis/export"},
    "icab": {"accept": r"\bicaB\b|polysaccharide N-deacetylase"
                       r"|adhesin biosynthesis polysaccharide N-deacetyl",
             "reject":   r"\bica[ACDR]\b|intercellular adhesion protein [ACD]\b"
                       r"|biosynthesis/export|glucosamine synthase"},
    "icac": {"accept": r"\bicaC\b|biosynthesis/export protein"
                       r"|intercellular adhesion protein C\b",
             "reject":   r"\bica[ABDR]\b|intercellular adhesion protein [ABD]\b"
                       r"|N-deacetylase|glucosamine synthase"},
    # SdrG == Fbe. Do not list both symbols in genes.txt.
    "sdrg": {"accept": r"\bsdrG\b|\bfbe\b|fibrinogen-?binding adhesin"
                       r"|MSCRAMM family adhesin SdrG",
             "reject":   r"\bsdr[CDEFH]\b|clumping factor|\bclf[AB]\b"
                       r"|B-like domain"},
    # EssD records split between the T7SS name and a bare 'DNA/RNA non-specific
    # endonuclease', which is also the title of native nucleases. Accepting only
    # the secretion-system wording filters the native ones out.
    "essd": {"accept": r"\bessD\b|ESAT-6 secretion|type VII secretion",
             "reject":   r"\bess[ABC]\b|\besx[AB]\b|\besa[ABCD]\b"},
    # UNVERIFIED: sesI returned no records under its symbol on the 40-gene run,
    # so no real titles were observed. Needs ALIASES before it can work at all.
    # Round 2 cannot verify sesI on its title: nr carries the same general
    # annotation that defeated the symbol query. The accept list therefore
    # admits generic surface-protein wording, and deliberately does NOT admit a
    # bare 'hypothetical protein' -- if that is what nr returns, sesI will come
    # back flagged, which is the honest outcome rather than a title match that
    # proves nothing. Paste the real title and this can be pinned exactly.
    # The real title, from the RP62A record: 'cell wall surface anchor family
    # protein'. My guess of 'cell wall-?anchored' did NOT match it ('anchor
    # family', not 'anchored'), so this is pinned to the observed wording.
    # ── sesI: ONLY identifiers that name this gene ──────────────────────────
    # The generic terms removed from here -- 'surface protein', 'LPXTG',
    # 'cell wall-?anchored', 'surface anchor family', 'cell wall surface
    # anchor' -- are shared with SdrC/D/E, SasG, ClfA/B, Pls and FnBP, and an
    # LPXTG anchor is carried by most Gram-positive surface proteins. In the
    # 29-genome v5 run, five of seven sesI positives were loci where ROUND 1
    # SAID THIS IS NOT THE GENE (homology_only, 27.7-38.2% identity) and round
    # 2 overturned it on exactly those words: 'cell surface protein' in three
    # genomes, 'LPXTG cell wall anchor domain-containing protein' in two.
    #
    # A coverage floor could not rescue them, which is why the terms had to go
    # rather than be gated the way bhp's domain titles are: the subject is
    # ~200 aa, so a generic hit reaches 100% subject coverage exactly as a real
    # one does. The genuine NIH04008 call and the spurious NIH05001 call were
    # BOTH 100% identity over 100% coverage.
    #
    # The consequence is the same bargain as opp3D/opp3E. Across all 9 sesI
    # loci, 0 of 900 nr hits name SesI, SERP1654 or AAW54982 -- nr does not
    # name this protein -- so a real copy now comes back FLAGGED rather than
    # present, and round-1 identity against the seeded reference is the
    # evidence that counts: ~100% for a genuine copy (NIH04008 100.0%,
    # NIH051668 95.1%) against 27.7-38.2% for background.
    "sesi": {"accept": r"\bsesI\b|SesI|\bSERP1654\b|AAW54982",
             "reject": r"\bses[ABCDEFGHJ]\b|\bsdr[CDEFGH]\b|\baap\b"
                       r"|\bembp\b|\bbhp\b|clumping factor"},


    # ── mirrored from PRODUCT_FILTERS in build_protein_db.py ────────────────
    # Round 2 must use the same definition of "this title names this gene" as
    # the build step, or the two rounds disagree. These nine were curated at
    # build time first and were missing here -- the same split that made bhp
    # behave differently in each round.

    "msra": {
        "accept": r"""\bmsr\(?A\)?\b|ABC-?F type ribosomal protection|macrolide.*(efflux|resistance)""",
        "reject": r"""peptide-?methionine|methionine sulfoxide|S-?oxide reductase|\bmsrB\b""",
    },
    "dfrc": {
        "accept": r"""\bdfrC\b|dihydrofolate reductase""",
        "reject": r"""\bdfr(S1|S2|[ABDGK])\b|\bfolA\b""",
    },
    # ACCEPT REQUIRES THE NAME. 'trimethoprim-resistant dihydrofolate
    # reductase' alone also matches DfrA, DfrK and DfrS1 -- and DfrS1 is what
    # confirmed dfrG at 97-100% in 4 genomes whose locus is the native
    # trimethoprim-SENSITIVE dfrC: identical coordinates to the dfrC call,
    # round-1 identity 38.9-40.7%, confirming organism [Bacteria], ranks 7-22.
    "dfrg": {
        "accept": r"""\bdfrG\b""",
        "reject": r"""\bdfr(S1|S2|[ABCDK])\b|trimethoprim-?sensitive|\bfolA\b""",
    },
    "erma": {
        "accept": r"""\berm\(?A\)?\b|23S rRNA .*adenine.*methyltransferase|rRNA adenine N-6-methyltransferase""",
        "reject": r"""\berm\(?[BCTFGY]\)?\b""",
    },
    "ermb": {
        "accept": r"""\berm\(?B\)?\b|23S rRNA .*adenine.*methyltransferase""",
        "reject": r"""\berm\(?[ACTFGY]\)?\b""",
    },
    "ermc": {
        "accept": r"""\berm\(?C\)?\b|23S rRNA .*adenine.*methyltransferase""",
        "reject": r"""\berm\(?[ABTFGY]\)?\b""",
    },
    "fusb": {
        "accept": r"""\bfusB\b|fusidic acid resistance""",
        "reject": r"""\bfus[CD]\b""",
    },
    "fusc": {
        "accept": r"""\bfusC\b|fusidic acid resistance""",
        "reject": r"""\bfus[BD]\b""",
    },
    "embp": {
        "accept": r"""\bembp\b|\bebh\b|Ebh-?like|hyperosmolarity resistance|extracellular matrix[- ]?binding protein""",
        "reject": r"""\bsdr[CDEFGH]\b|\baap\b|\bbhp\b|clumping factor""",
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


# Titles that name the right gene but flag the RECORD as incomplete. These
# still count as accepting -- ", partial" is extremely common in nr and often
# the best available evidence (blaI's top hit is 'BlaI/MecI/CopY family
# transcriptional regulator, partial') -- but an intact record is preferred
# when one exists. Confirming sdrF on 'sdrF protein, truncation' while
# full-length SdrF records sit in the same hit list is weaker than it looks.
# nr deflines can be enormous. A MULTISPECIES WP_ record lists every redundant
# accession and organism in ONE title field -- the IS6-family transposase
# WP_002319817.1 spans dozens of species on a single line -- and makeblastdb
# preserves that whole string, so stitle can exceed Python's default 128 KB csv
# field cap and raise '_csv.Error: field larger than field limit'. Raising the
# limit is the right fix rather than truncating the title: classify_title()
# searches the ENTIRE defline, and a gene name or paralog term can appear late
# in it (blaI confirms on 'BlaI/MecI/CopY family ... , partial'), so truncation
# would silently change calls.
_lim = sys.maxsize
while True:
    try:
        csv.field_size_limit(_lim)
        break
    except OverflowError:          # 32-bit C long on some platforms
        _lim //= 10

TRUNCATED_RE = re.compile(r",\s*partial\b|\btruncation\b|\btruncated\b"
                          r"|\bfragment\b", re.I)


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
            # Preserve keys the TSV cannot express. A plain assignment here
            # would silently drop 'domain_accept', so a --patterns file that
            # only meant to tweak bhp's wording would disable its coverage-gated
            # domain rule and turn NIH05003 back into a flag, with nothing in
            # the output saying why.
            entry = dict(pats.get(gene, {}))
            entry["accept"] = parts[1].strip() if len(parts) > 1 else ""
            entry["reject"] = parts[2].strip() if len(parts) > 2 else ""
            pats[gene] = entry
    return pats


def classify_title(gene, title, pats):
    """'accept' | 'reject' | 'unknown' for this title, for this gene.

    reject is checked first: a paralog must never confirm, even when its
    defline happens to also contain an accept term (e.g. a title naming both
    SdrF and SdrG).
    """
    # Only the PRIMARY title, not the redundant synonyms.
    #
    # An nr MULTISPECIES record aggregates every collapsed accession into one
    # defline: 'MULTISPECIES: ABC transporter ATP-binding protein
    # [Staphylococcus] >EGS79756.1 ... >KTF27926.1 peptide ABC transporter
    # ATP-binding protein ...'. Searching the whole string turns "the title
    # names this gene" into "ANY of the N synonyms names it", which is a much
    # weaker claim.
    #
    # That is how opp3D was confirmed at 100% identity in 8 genomes whose round 1
    # said homology_only at ~29%: the primary annotation is the generic 'ABC
    # transporter ATP-binding protein', and the accept pattern matched a
    # SYNONYM further down the defline. The call then flipped between runs
    # depending on which record ranked first -- 10 of 29 genomes changed between
    # v2 and v3 with identical round-1 calls.
    #
    # The first title is the representative annotation for the accession that
    # was actually hit; the rest describe other accessions that happen to share
    # the sequence. Matching only the first is the stricter and more honest test.
    t = (title or "").split(" >")[0].lower()
    p = pats.get(gene.lower())

    if p is None:
        # No curated patterns: fall back to the bare symbol as a word. Weak,
        # and the caller notes it so the gap is visible rather than silent.
        return "accept" if re.search(rf"\b{re.escape(gene.lower())}\b", t) else "unknown"

    # re.I so that these patterns can be written IDENTICALLY to the
    # PRODUCT_FILTERS table in build_protein_db.py, which uses re.I in
    # product_ok(). Two tables with different casing conventions drift: bhp was
    # fixed in one and not the other, and the two rounds then disagreed about
    # which titles name the gene. A consistency test over real manifest titles
    # guards the rest.
    if p.get("reject") and re.search(p["reject"], t, re.I):
        return "reject"
    if p.get("accept") and re.search(p["accept"], t, re.I):
        return "accept"
    return "unknown"


def title_matches(gene, title, pats=None):
    """Backwards-compatible boolean wrapper."""
    return classify_title(gene, title, pats or GENE_PATTERNS) == "accept"


def domain_accept(gene, hit, pats):
    """Second chance for a hit whose title names a DOMAIN, not the gene.

    Returns True only if the title matches the gene's 'domain_accept' pattern
    AND the hit covers a full-length subject. Callers use it only after
    classify_title() has already returned something other than 'accept', and
    only after 'reject' has been honoured.

    Why this exists. nr annotates some repeat-adhesin families by domain and
    never by name. For NIH05003's bhp locus -- a single HSP, 99.917% identity,
    aligned 1-2402 of the 2402 aa AAW53225 reference, no frameshift, as
    unambiguous as this pipeline gets -- all 100 nr hits came back titled
    'putative Ig domain-containing protein' and not one named biofilm, Bap or
    Bhp. A title-only test cannot confirm that locus at any threshold.

    Why it is gated on COVERAGE OF THE SUBJECT rather than just accepting the
    domain title. 'Ig domain-containing' is shared by Bhp, Aap, SasG, Embp and
    the Sdr family, which is precisely what the 17 flagged bhp loci are. Those
    loci hit full-length Bhp records at 8-40% subject coverage, because they
    are a shared repeat and not the gene; a genuine full-length copy covers
    ~100%. Over the 29 genomes the floor below confirms 4 loci -- the three
    already confirmed by a biofilm-titled hit, plus NIH05003 -- and none of the
    17 flagged or 7 absent. The separation is the evidence for the rule.

    min_hit_len matters as much as min_hit_cov: NIHLM018's top hit is a
    'hypothetical protein, partial' at 100% identity and 100% subject coverage
    because the SUBJECT is a tiny fragment record. Coverage alone would pass it.
    """
    p = (pats or {}).get(gene.lower()) or {}
    spec = p.get("domain_accept")
    if not spec:
        return False
    t = (hit.get("title") or "").split(" >")[0].lower()
    if not re.search(spec["pattern"], t, re.I):
        return False
    try:
        hit_len = int(hit.get("hit_len") or 0)
        hit_cov = float(hit.get("hit_cov_pct") or 0)
    except (TypeError, ValueError):
        return False
    return (hit_len >= spec.get("min_hit_len", 0)
            and hit_cov >= spec.get("min_hit_cov", 0))


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
    ap.add_argument("--present-unverified", action="store_true",
                    help="set when extract_loci.py ran WITHOUT --all, so genes "
                         "already called present have no locus and were never "
                         "nr-searched. Carries round 1's 'present' through "
                         "instead of reporting 'flagged'. Without this flag a "
                         "missing locus for a present gene is treated as a "
                         "genuine verification failure.")
    ap.add_argument("--short-gene-aa", type=int, default=80,
                    help="a gene whose called region translates to fewer than "
                         "this many residues is reported 'unverifiable_short' "
                         "when no nr hit overlaps, instead of 'flagged' "
                         "(default 80; 0 disables). A short alignment cannot "
                         "reach the nr E-value threshold -- agrD is a 46-codon "
                         "propeptide -- so round 2 can say nothing either way, "
                         "which is not the same as the locus looking doubtful.")
    ap.add_argument("--min-overlap-frac", type=float, default=0.0,
                    help="an nr hit counts only if it covers at least this "
                         "fraction of the called region (default 0.0 = off). "
                         "OFF BY DEFAULT because it was the wrong hypothesis: "
                         "it was added to fix dfrC, but every hit there covers "
                         "0.98 of the called region (the neighbouring DegV "
                         "protein included), so it discriminates nothing -- "
                         "--min-contained-frac is what fixes dfrC. At 0.5 it "
                         "changed no call on a 45-gene panel while discarding "
                         "all nr evidence for aap, whose hits are short "
                         "relative to an 8.5 kb locus, leaving a less "
                         "informative message. Kept as a lever for edge-clip "
                         "hits on SHORT genes, where it does discriminate.")
    ap.add_argument("--min-contained-frac", type=float, default=0.80,
                    help="an nr hit counts only if at least this fraction of "
                         "its own query span lies inside the called region "
                         "(default 0.80; 0 disables). Excludes a hit that runs "
                         "out through --flank into a NEIGHBOURING gene: dfrC's "
                         "top hit aligned 261 aa from the gene into the "
                         "adjacent DegV ORF (containment 0.62) and outscored "
                         "the real DfrC's 161 aa fully-contained alignment on "
                         "bitscore alone, 428 vs 333, despite 78.2% identity "
                         "against 100.0%.")
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
        # Approximate protein length from the called region. Good enough to
        # recognise a gene too short for nr to verify, and avoids making
        # finalize_calls depend on the manifest just for a length.
        ref_aa = ((ge - gs + 1) // 3) if coords else 0

        # All overlapping hits, bitscore-ordered, with their genomic spans.
        #
        # 'Overlapping' used to mean >= 1 bp, which let a NEIGHBOURING gene
        # count as a hit on this one. The locus query carries --flank bp on each
        # side, so for a short gene the flank dominates: dfrC is 492 bp with 600
        # bp of flank, 55% of the query, and an adjacent DegV family protein
        # became the top overlapping hit at 428 bits. That set the bitscore
        # ceiling, and the real dihydrofolate reductase at 333 bits (77.8%) fell
        # below the 80% floor -- so a gene called present at 99.4%/100% in round
        # 1 was reported flagged, on the strength of its neighbour.
        #
        # Two tests, catching opposite failures:
        #
        #   coverage   -- how much of the CALLED REGION the hit spans. OFF by
        #                 default; see --min-overlap-frac for why.
        #   containment-- how much of the HIT lies inside the called region.
        #                 Excludes a hit that runs out through the flank into a
        #                 neighbouring gene. This is the one that works.
        #
        # Containment is what matters for dfrC, and it is not what I first
        # assumed. Every hit there covers 0.98 of the called region,
        # DegV included, so coverage alone discriminates nothing. What separates
        # them is alignment length: the real DfrC aligns 161 aa, entirely inside
        # the gene (containment 1.00), while the DegV hit aligns 261 aa from the
        # gene's start out through 300 bp of flank into the adjacent DegV ORF
        # (containment 0.62). The longer alignment wins on BITSCORE -- 428 vs
        # 333 -- despite 78.2% identity against DfrC's 100.0%, and that inflated
        # ceiling then pushed the true hit below the 80% bitscore floor.
        #
        # So a long hit extending past the gene is NOT automatically evidence
        # about the gene; past some point it is evidence about its neighbour.
        overlapping, near_misses = [], 0
        if win and gs is not None:
            glen = max(ge - gs + 1, 1)
            for rank, h in enumerate(hits, 1):
                g1 = win["win_start"] + int(h["query_start"]) - 1
                g2 = win["win_start"] + int(h["query_end"]) - 1
                ov = min(g2, ge) - max(g1, gs) + 1
                if ov <= 0:
                    continue
                span = max(g2 - g1 + 1, 1)
                if (ov / glen < args.min_overlap_frac
                        or ov / span < args.min_contained_frac):
                    near_misses += 1
                    continue
                overlapping.append(dict(h, g1=g1, g2=g2, _rank=rank,
                                        _ov_frac=ov / glen,
                                        _contained=ov / span))

        top = overlapping[0] if overlapping else None
        confirming = None
        rejected = 0
        notes = []

        if overlapping:
            top_bits = float(top["bit_score"])
            floor = top_bits * args.min_bitscore_frac
            fallback = None          # accepting, but a truncated record
            for h in overlapping:
                verdict = classify_title(gene, h["title"], pats)
                if verdict == "reject":
                    rejected += 1
                    continue
                if verdict != "accept":
                    # A domain-named record may still confirm, but only if it
                    # covers a full-length subject -- see domain_accept(). The
                    # 'reject' check above still runs first, so a paralog named
                    # outright is never rescued this way.
                    if not domain_accept(gene, h, pats):
                        continue
                    notes.append(
                        f"confirmed on a domain-level annotation "
                        f"({h['title'].split(' >')[0][:46]}) covering "
                        f"{float(h['hit_cov_pct']):.0f}% of a {h['hit_len']} aa "
                        f"subject; nr does not name this gene")
                if float(h["bit_score"]) < floor:
                    # A matching hit exists but is much weaker than the best
                    # overlapping hit; report rather than confirm on it.
                    notes.append(
                        f"best matching hit {h['accession']} is rank "
                        f"{h['_rank']} at {h['bit_score']} bits, below "
                        f"{args.min_bitscore_frac:.0%} of top ({top['bit_score']})")
                    break
                if TRUNCATED_RE.search(h["title"]):
                    # keep looking for an intact record of the same gene
                    if fallback is None:
                        fallback = h
                    continue
                confirming = h
                break

            if confirming is None and fallback is not None:
                confirming = fallback

            # Overturning a homology_only call takes the TOP overlapping hit,
            # not a lower-ranked fallback.
            #
            # The bitscore fallback exists for an ambiguous TOP title -- sdrH
            # confirms on rank 4 because ranks 1-3 are the uninformative
            # 'fibrinogen-binding protein'. Contradicting round 1 is a different
            # matter: there, round 1 has positive evidence the gene is NOT here.
            #
            # dfrG showed why. Its locus is the native dfrC gene: round 1 scored
            # it 99.4% against the DfrC panel and 39.5% against DfrG. Round 2
            # then confirmed dfrG on a RANK 47 hit at 98.1%, organism
            # 'Bacteria'. A record titled DfrG that is 98% identical to a locus
            # only 39.5% similar to every DfrG reference is mis-annotated, not
            # evidence -- and accepting it would have reported an acquired
            # trimethoprim-resistance gene in a susceptible genome.
            if (confirming is not None
                    and c["call"] == "homology_only"
                    and confirming["_rank"] != top["_rank"]):
                notes.append(
                    f"round 1 called this homology_only at "
                    f"{c['identity_pct']}% identity, and the matching nr hit "
                    f"{confirming['accession']} is rank {confirming['_rank']}, "
                    f"not the top overlapping hit -- too weak to overturn "
                    f"round 1, so it does not confirm")
                confirming = None
                notes.append("confirming record is annotated as "
                             "truncated/partial; no intact record overlapped")

        if gene.lower() not in pats:
            notes.append(f"no curated title patterns for '{gene}' -- matched on "
                         f"the bare symbol only")

        # ── round 1 found no locus at all ───────────────────────────────────
        # Round 2 only ever asks 'does an nr hit overlap the called region?'.
        # With no called region there is nothing to overlap, so the gene used to
        # fall into the 'flagged' branch below -- the same label given to a gene
        # whose locus WAS searched and came back unverifiable. Those are
        # different claims: 'absent' is a completed negative (panel built, right
        # protein, right length, zero qualifying hits), whereas 'flagged' means
        # something aligned and a human needs to look. On a closed genome the
        # former is a publishable result and the latter is a to-do item.
        #
        # 'no_reference' is kept separate again: there, nothing was searched,
        # so the gene is neither present nor absent -- merely untested.
        if not coords and c["call"] in ("absent", "no_reference"):
            final = c["call"]
            nr_status = (
                "absent in round 1 (no locus to verify); round 2 not applicable"
                if c["call"] == "absent" else
                "no reference proteins were built for this gene; neither "
                "present nor absent -- untested"
            )
        elif top is None and ref_aa and ref_aa < args.short_gene_aa:
            # Too short for nr to verify at all. agrD is a 46-codon autoinducer
            # propeptide; a 46 aa alignment cannot reach --nr_evalue 1e-10, so
            # no hit overlaps and the gene was reported flagged despite scoring
            # 100% identity over 100% coverage in round 1.
            #
            # Reported as its own state rather than fixed by lowering the
            # e-value, which would loosen verification for all 45 genes to
            # rescue one. 'flagged' would wrongly imply the locus looks
            # questionable; it does not -- round 2 simply cannot speak to it.
            # Round 1's verdict stands when it says present: it is the only
            # evidence that exists, it is strong (agrD: 100% identity over 100%
            # coverage against 20 references), and the agr operon corroborates it
            # (agrB, agrC and agrA all present and nr-confirmed, and the operon
            # is agrBDCA). 'unverifiable_short' is kept only for a short gene
            # whose round-1 call was NOT present, where there is no positive
            # evidence to carry forward.
            #
            # Consequence to be aware of: 'present' therefore no longer means
            # "nr-verified" in every case. The nr_verification text below always
            # states when round 2 could not assess the gene, so the weaker
            # evidence stays auditable in the table.
            final = "present" if c["call"] == "present" else "unverifiable_short"
            nr_status = (f"round 1 called this present at "
                         f"{c['identity_pct']}% identity over "
                         f"{c['coverage_pct']}% coverage, but the reference is "
                         f"only {ref_aa} aa -- too short for an alignment to "
                         f"reach the nr E-value threshold, so round 2 cannot "
                         f"verify it either way. Trust round 1 here.")
        elif (top is None and args.present_unverified
              and c["call"] == "present"):
            # Round 2 was deliberately not run for this gene.
            #
            # extract_loci.py only extracts a locus for a 'present' call when
            # --all is given. Without it there is no locus, so no nr hit can
            # overlap, and this would otherwise fall through to 'flagged' --
            # turning EVERY present call into a flag. Round 1's verdict is
            # carried through instead, with the gap stated.
            #
            # What this gives up: round 2's false-positive check. mecA was
            # 'present' at 99.6% in round 1 and only round 2 showed the locus
            # was a native PBP; msrA was 'present' at 98.3% on the housekeeping
            # methionine sulfoxide reductase. Both were ultimately fixed in the
            # PANEL, but it was round 2 that exposed them. Set
            # --nr_verify_present true to restore the check.
            final = "present"
            nr_status = ("present on round 1 only; round 2 not run for genes "
                         "already called present (nr_verify_present=false). No "
                         "nr evidence either way -- the false-positive check "
                         "that caught mecA and msrA is not applied to this "
                         "call.")
        elif top is None and c["call"] == "homology_only":
            # Superfamily-level match in round 1 and nothing in nr overlapping
            # it either -- the same two-line agreement as the branch below, so
            # it gets the same verdict rather than landing on a review list.
            final = "absent"
            nr_status = (f"absent: round 1 matched at {c['identity_pct']}% "
                         f"identity (below the superfamily threshold) and no "
                         f"nr hit overlaps the locus")
        elif top is None:
            nr_status = ("flagged: a locus was called but no nr hit overlaps it")
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
        elif c["call"] == "homology_only":
            # Two independent lines of evidence now agree the gene is not here:
            # round 1 matched only at superfamily level (< --homology-max), and
            # round 2 found no overlapping nr hit that names the gene. Calling
            # that 'flagged' would put eight background genes on a review list
            # they do not belong on. 'flagged' is reserved for the genuinely
            # unresolved -- a locus that looks like the gene but could not be
            # verified.
            final = "absent"
            nr_status = (f"absent: round 1 matched at {c['identity_pct']}% "
                         f"identity (below the superfamily threshold) and no "
                         f"overlapping nr hit identifies {gene}; the locus is "
                         f"{top['accession']} ({top['title'][:50]})")
        else:
            nr_status = (f"flagged: no overlapping nr hit identifies {gene}; "
                         f"top hit is {top['accession']} "
                         f"({top['title'][:60]})")
            if rejected:
                notes.append(f"{rejected} overlapping hit(s) matched a paralog "
                             f"pattern")
            notes.append(f"{len(overlapping)} overlapping hit(s) examined")
            final = "flagged"

        # A frameshifted locus must not be reported as a plain 'present'.
        # It scores 0 in the matrix -- see the matrix section for why.
        #
        # Round 1 measures sequence similarity and round 2 asks what protein the
        # sequence resembles. NEITHER looks at reading frame, so both can be
        # fully satisfied by a pseudogene: ebh/embp on LM088-RCB001 scores 99.6%
        # identity to an S. epidermidis Ebh-like protein across 30.5 kb while
        # being interrupted by two frameshifts, at coordinates reproduced
        # exactly by two independent panels.
        #
        # 'present_frameshift' keeps that visible instead of laundering it into
        # 'present'. It still scores 1 in the matrix, because the sequence IS
        # there and scoring 0 would assert absence of something plainly present;
        # the distinction between "gene present" and "functional protein
        # present" is what the label carries.
        if c.get("frame_status") and final == "present":
            final = "present_frameshift"
            notes.append(
                f"the called region spans multiple reading frames "
                f"({c['frame_status']}, junction(s) near "
                f"{c.get('frameshift_at', '?')}), so an intact full-length "
                f"protein is NOT supported -- nr confirms the sequence, not the "
                f"reading frame. Check reads at those coordinates to tell a "
                f"pseudogene from an assembly indel")

        if notes:
            nr_status += " | " + "; ".join(notes)

        # Report the hit the decision rests on, not merely the highest-scoring.
        # The columns below are therefore named nr_confirming_*, NOT nr_top_hit_*:
        # when GENE_PATTERNS rejects ranks 1-2 as paralogs or ambiguous titles,
        # the reported hit is rank 3 (sdrH on LM088-RCB001 did exactly that).
        # nr_hit_rank == 1 is the only case where the confirming hit is also the
        # top-scoring one.
        top = confirming or top

        # NOTE: 'hit_cov_pct' is coverage of the SUBJECT protein, 'query_cov_pct'
        # coverage of the QUERY locus -- which carries --flank bp on each side,
        # capping query coverage well below 100% even for a full-length match
        # (1374/1974 = 70% for sdrH). The two are not comparable, so the basis
        # actually used is recorded in nr_coverage_basis rather than left
        # ambiguous.
        if top and top.get("hit_cov_pct"):
            cov_col, cov_basis = top["hit_cov_pct"], "subject"
        elif top and top.get("query_cov_pct"):
            cov_col, cov_basis = top["query_cov_pct"], "query_locus+flank"
        else:
            cov_col, cov_basis = "", ""

        rows.append({
            "gene": gene,
            "round1_call": c["call"],
            "round1_identity_pct": c["identity_pct"],
            "round1_coverage_pct": c["coverage_pct"],
            "contig": c["contig"],
            "coordinates": coords,
            "nr_confirming_accession": top["accession"] if top else "",
            "nr_confirming_title": top["title"][:80] if top else "",
            "nr_confirming_organism": top["organism"] if top else "",
            "nr_identity_pct": top["identity_pct"] if top else "",
            "nr_hit_coverage_pct": cov_col,
            "nr_coverage_basis": cov_basis,
            "nr_hit_rank": top.get("_rank", "") if top else "",
            "frame_status": c.get("frame_status", ""),
            "frameshift_at": c.get("frameshift_at", ""),
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
    # final_call values: present | absent | flagged | no_reference
    #   present      verified by nr, OR too short for nr to assess with a
    #                round-1 call of present (see --short-gene-aa)
    #   present_frameshift  sequence verified, but the locus spans multiple
    #                reading frames -- gene present, functional PROTEIN not
    #                supported. Scores 0.
    #
    #                Scores 0 because the matrix answers "does this genome make
    #                this protein", and a locus interrupted by frameshifts does
    #                not. The sequence-level finding is not lost: the row keeps
    #                final_call='present_frameshift' plus the junction
    #                coordinates in frameshift_at, so a sequence-level question
    #                can still be answered from --out. Scoring it 1 would put a
    #                probable pseudogene in a presence figure as an ordinary
    #                positive, which is the harder error to notice.
    #   absent       searched against a valid panel, nothing qualified
    #   flagged      a locus was called but could not be verified -- review
    #   no_reference no panel was built; untested, neither present nor absent
    # 'present' scores 1. 'flagged' is deliberately 0, because flagged means NOT
    # VERIFIED, and scoring an unverified gene as present is how a false
    # positive reaches a figure. The reason for any 0 is in --out.
    #
    # ONE EXCEPTION: 'unverifiable_short' with a round-1 call of 'present'
    # also scores 1. Round 2 is structurally incapable of assessing a gene that
    # short, so scoring it 0 would not be caution -- it would discard the only
    # evidence that exists and assert absence on no grounds. agrD is the case:
    # 100% identity over 100% coverage against 20 references, inside an
    # otherwise intact agr operon (agrB, agrC and agrA all present and
    # nr-confirmed), where agrBDCA argues strongly for agrD. The exception is
    # narrow by design: round 1 must say 'present', and the row still carries
    # 'unverifiable_short' so the weaker evidence is visible rather than
    # laundered into a plain 'present'.
    out_dir = os.path.dirname(os.path.abspath(args.out)) or "."
    sample = args.sample or os.path.basename(out_dir)
    matrix_path = args.matrix_out or os.path.join(out_dir,
                                                  "gene_presence_matrix.tsv")
    with open(matrix_path, "w", newline="") as fh:
        fh.write("sample\tgene\tpresent\n")
        for r in sorted(rows, key=lambda x: x["gene"]):
            # Only 'present' scores 1. A short gene whose round 1 said present
            # is now final_call='present' itself, so no special case is needed;
            # 'unverifiable_short' only survives where round 1 did NOT say
            # present, and correctly scores 0.
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

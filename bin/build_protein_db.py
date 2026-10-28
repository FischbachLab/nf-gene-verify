#!/usr/bin/env python3
"""Build a staph-restricted protein reference DB by gene name from NCBI.

Queries NCBI Protein with '<gene>[Gene Name] AND Staphylococcus[Organism]
AND srcdb_refseq[PROP]', prefers RefSeq records, deduplicates by sequence
hash, drops ', partial' records, and writes per-gene FASTAs, a combined
FASTA, and a manifest CSV.

Genes with no RefSeq records (common for plasmid-borne resistance genes such
as aacA, which RefSeq annotates under the fused name aacA-aphD) are retried
without the RefSeq restriction; aliases can be added with --alias and
user-supplied reference proteins with --extra-fasta. A gene for which no
reference can be built is reported as 'no_reference' by call_genes.py —
distinct from 'absent' (searched, no hit).
"""
import argparse
import csv
import hashlib
import json
import os
import re
import sys
import time
import urllib.parse
import urllib.request

BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
DEFAULT_GENES = ["agrA", "agrB", "blaI", "blaZ", "fosB", "sarA", "sdrF", "sdrH"]

# Known gene-name fusions and PRODUCT-NAME tokens. Keys are lowercase.
# Tokens containing '[' or '"' are used verbatim as field-restricted queries;
# bare tokens become '<token>[Gene Name]'.
#
# Product names matter as much as the symbol. Most NCBI protein records for
# these genes carry a product name and no gene symbol -- 'accessory gene
# regulator B' and 'serine-aspartate repeat protein F' contain no 'agrB' or
# 'sdrF' token at all -- so a bare '<gene>[Gene Name]' query retrieves a
# narrow, unrepresentative slice of the available alleles.
ALIASES = {
    "aaca": ["aacA-aphD"],  # bifunctional AAC/APH (Tn4001); plasmid-borne,
                            # essentially absent from RefSeq gene-name records
    "fdh": ['"formate dehydrogenase"'],
        # bare 'fdh[All Fields]' wrongly matches zinc-dependent ALCOHOL
        # dehydrogenases; the [Protein Name] field index is unreliable
        # (1 hit), so use the untagged quoted phrase
    "is256": ['"IS256 family transposase"[Protein Name]'],
        # IS256 is an insertion sequence; its verifiable gene is the transposase

    # ── mecA: TWO UNRELATED PROTEINS SHARE THIS SYMBOL ──────────────────────
    # 1. the ClpC adaptor protein MecA (~240 aa), a housekeeping regulator
    #    present in every staphylococcus and bacillus
    # 2. PBP2a (~668 aa), the methicillin-resistance transpeptidase
    # A bare 'mecA' query returns overwhelmingly the former: on LM088-RCB001 the
    # panel came back 17 adaptors to 3 PBP2a, the adaptor won the call at 99.6%
    # identity over 100% coverage of a 240-codon locus, and the pipeline
    # reported mecA PRESENT for a methicillin-SENSITIVE genome. Querying the
    # product name instead of the symbol is what keeps the adaptor out.
    "meca": ['"penicillin-binding protein 2a"[Title]',
             '"penicillin-binding protein PBP2a"[Title]',
             '"PBP2a family beta-lactam-resistant peptidoglycan transpeptidase"[Title]',
             '"beta-lactam-resistant peptidoglycan transpeptidase"[Title]'],

    # ── ica operon ──────────────────────────────────────────────────────────
    # IcaD is ~101 aa and IcaR ~186 aa, so neither tolerates a paralogue from
    # elsewhere in the operon: icaA is ~412 aa and would dominate on length.
    "icad": ['"intercellular adhesion protein D"[Title]',
             '"poly-beta-1,6-N-acetyl-D-glucosamine synthesis protein IcaD"[Title]',
             '"biofilm PIA synthesis protein IcaD"[Title]'],
    "icar": ['"biofilm operon icaADBC HTH-type negative transcriptional regulator IcaR"[Title]',
             '"ica operon transcriptional regulator IcaR"[Title]',
             '"biofilm operon repressor IcaR"[Title]'],

    # ── ACME components ─────────────────────────────────────────────────────
    # opp3D/opp3E returned NO records under their symbols on the 40-gene run:
    # NCBI indexes these by product. Note D and E are BOTH ATP-binding subunits,
    # so their panels overlap heavily -- acceptable here, because ACME typing
    # only needs the CLUSTER detected, not each subunit resolved individually.
    # The records spell these HYPHENATED -- 'opp-3D', not 'opp3D' -- so the
    # unhyphenated symbol query returned nothing on four consecutive runs. The
    # accessions below are from S. epidermidis I11OR1 and sit INSIDE a
    # characterised ACME type III element, which makes them ACME-borne
    # references rather than generic oligopeptide permease subunits:
    #   opp3D = ASM79718.1  (768 bp -> ~255 aa)
    #   opp3E = ASM79719.1  (636 bp -> ~211 aa)
    "opp3d": ['"opp-3D"[Gene Name]', 'ASM79718[All Fields]',
              'WP_002490823[All Fields]',
              '"oligopeptide ABC transporter ATP-binding protein"[Title]',
              '"oligopeptide transport ATP-binding protein OppD"[Title]'],
    "opp3e": ['"opp-3E"[Gene Name]', 'ASM79719[All Fields]',
              'WP_002490830[All Fields]',
              '"oligopeptide ABC transporter ATP-binding protein"[Title]',
              '"oligopeptide transport ATP-binding protein OppF"[Title]'],
    # kdp: required to detect ACME types IV and V. CAUTION: like arc, the kdp
    # operon is ALSO native in staphylococci -- the ACME copy is a second one.
    # A kdp hit is therefore not ACME evidence on its own; only co-location
    # with the element is. acme_type.py enforces that.
    "kdpa": ['"potassium-transporting ATPase subunit A"[Title]',
             '"K(+)-transporting ATPase subunit A"[Title]'],
    "kdpb": ['"potassium-transporting ATPase subunit B"[Title]',
             '"K(+)-transporting ATPase subunit B"[Title]'],
    "kdpc": ['"potassium-transporting ATPase subunit C"[Title]',
             '"K(+)-transporting ATPase subunit C"[Title]'],
    "kdpd": ['"sensor histidine kinase KdpD"[Title]',
             '"potassium-transporting ATPase KdpD"[Title]'],
    "kdpe": ['"response regulator KdpE"[Title]',
             '"KdpE family response regulator"[Title]'],

    # ── parenthesised official names ────────────────────────────────────────
    # The nomenclature-committee spellings are erm(A), mph(C), msr(A). Those
    # parentheses are Entrez GROUPING operators, so 'erm(A)[Gene Name]' is not
    # the literal string intended, and the token also lands on the filesystem as
    # <gene>.fasta / <gene>_locus.fasta where parens and spaces break the shell
    # portions of the pipeline. genes.txt therefore uses ermA/mphC/msrA, and the
    # official form is restored here as a QUOTED query term -- quoting makes the
    # parentheses literal to Entrez.
    "erma": ['"erm(A)"[Gene Name]', '"rRNA adenine N-6-methyltransferase"[Title]'],
    "ermb": ['"erm(B)"[Gene Name]'],
    "ermc": ['"erm(C)"[Gene Name]'],
    "mphc": ['"mph(C)"[Gene Name]', '"macrolide 2\'-phosphotransferase"[Title]'],
    "msra": ['"msr(A)"[Gene Name]', '"ABC-F type ribosomal protection"[Title]'],

    # ── sesI: retrieved by LOCUS TAG, not by symbol ─────────────────────────
    # 'sesI[Gene Name]' returns nothing. The S. epidermidis RP62A record
    # (CP000029.1, 1693235-1693840 complement, 606 bp = ~201 aa) carries only a
    # general protein annotation, so the symbol never appears in the title.
    # The locus tag SERP1654 and the protein accession AAW54982.1 are the
    # handles that work.
    #
    # Consequence: this query is PRECISE but NARROW -- expect a very small
    # panel, possibly n=1. A one-allele reference set is the exact problem that
    # made agrB read 49.2% in round 1 while nr said 99.5%, so if sesI matters,
    # seed it (data/sesI_seed.faa) and widen with a product-name alias once the
    # manifest shows what NCBI actually calls this protein.
    "sesi": ['SERP1654[All Fields]', 'AAW54982[All Fields]'],

    # Bhp is the Bap homologue; most records carry the product name, not 'bhp'.
    "bhp": ['"biofilm-associated protein"[Title]',
            '"Bap/Bhp family"[Title]'],
    # Ebh records are titled by product; the symbol rarely appears.
    "ebh": ['"hyperosmolarity resistance protein Ebh"[Title]',
            '"Ebh-like protein"[Title]'],

    # ── embp == ebh ─────────────────────────────────────────────────────────
    # Embp (S. epidermidis) and Ebh (S. aureus) are orthologues of the same
    # ~1 MDa ECM-binding giant (~10,200 aa / ~30.5 kb of DNA). DO NOT list both
    # symbols in genes.txt -- they build two panels and emit two rows for one
    # locus, exactly as fbe/sdrG did.
    #
    # 'embp[Gene Name]' returns NOTHING: NCBI indexes the S. epidermidis protein
    # by product, and titles it 'Ebh-like' even in that species. That empty
    # query is why embp came back no_reference on the 40-gene run while the
    # SAME gene had already been called under 'ebh' at 1365187-1395715
    # (30,529 bp = 10,176 codons), nr-confirmed at 99.6%.
    #
    # Querying by product makes the panel S. epidermidis Embp rather than
    # S. aureus Ebh, which should lift round-1 identity from 48.2% (the
    # between-species distance) to near 100% and turn 'divergent' into
    # 'present'.
    "embp": ['"extracellular matrix-binding protein"[Title]',
             '"extracellular matrix binding protein"[Title]',
             '"Embp"[Title]',
             '"Ebh-like protein"[Title]',
             '"hyperosmolarity resistance"[Title]'],

    # ── product-name tokens for the default S. epidermidis gene set ──────────
    "agra": ['"accessory gene regulator A"[Title]',
             '"response regulator AgrA"[Title]',
             '"quorum-sensing response regulator"[Title]'],
    "agrb": ['"accessory gene regulator B"[Title]'],
        # agrB/agrD define the agr type (I-IV); a single allele is not enough.
        # LM088-RCB001 scored 49.2% identity in round 1 against a one-allele
        # reference and 99.5% against nr -- a different agr type, not a
        # different gene.
    "blai": ['"beta-lactamase repressor"[Title]',
             '"penicillinase repressor"[Title]',
             '"BlaI/MecI/CopY family transcriptional regulator"[Title]'],
    "blaz": ['"beta-lactamase BlaZ"[Title]', '"penicillinase"[Title]'],
    "fosb": ['"metallothiol transferase"[Title]',
             '"bacillithiol transferase"[Title]',
             '"fosfomycin resistance"[Title]'],
    "sara": ['"staphylococcal accessory regulator A"[Title]',
             '"transcriptional regulator SarA"[Title]'],
    "sdrf": ['"serine-aspartate repeat protein F"[Title]',
             '"serine-aspartate repeat-containing protein F"[Title]'],
        # SdrF length varies 780-1750 aa between strains (SD-repeat copy
        # number). A long reference gives low COVERAGE for a short allele and
        # vice versa, so the panel needs both ends of the range.
    "sdrh": ['"MSCRAMM-like protein SdrH"[Title]', '"SdrH family protein"[Title]'],
        # most SdrH records carry no gene symbol at all
}


def _open_with_retry(url, timeout, retries=5):
    """GET with exponential backoff on HTTP 429/5xx (NCBI rate limits)."""
    delay = 5
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=timeout) as r:
                return r.read()
        except urllib.error.HTTPError as e:
            if e.code == 429 or e.code >= 500:
                if attempt == retries - 1:
                    raise
                wait = delay * (2 ** attempt)
                print(f"[WARN] HTTP {e.code}, retrying in {wait}s "
                      f"({attempt + 1}/{retries})")
                time.sleep(wait)
            else:
                raise
    raise RuntimeError("unreachable")


def esearch(term, retmax, api_key):
    params = {"db": "protein", "term": term, "retmax": retmax, "retmode": "json"}
    if api_key:
        params["api_key"] = api_key
    url = f"{BASE}esearch.fcgi?{urllib.parse.urlencode(params)}"
    return json.loads(_open_with_retry(url, 30))["esearchresult"]["idlist"]


def esummary(ids, api_key):
    params = {"db": "protein", "id": ",".join(ids), "retmode": "json"}
    if api_key:
        params["api_key"] = api_key
    url = f"{BASE}esummary.fcgi?{urllib.parse.urlencode(params)}"
    return json.loads(_open_with_retry(url, 30))["result"]


def efetch_fasta(ids, api_key):
    """Return {unversioned_accession: sequence} from a batch FASTA fetch."""
    params = {"db": "protein", "id": ",".join(ids),
              "rettype": "fasta", "retmode": "text"}
    if api_key:
        params["api_key"] = api_key
    url = f"{BASE}efetch.fcgi?{urllib.parse.urlencode(params)}"
    text = _open_with_retry(url, 60).decode()
    seqs, acc, chunks = {}, None, []
    for line in text.splitlines():
        if line.startswith(">"):
            if acc is not None:
                seqs[acc] = "".join(chunks)
            acc = line[1:].split()[0].split(".")[0]  # strip version suffix
            chunks = []
        else:
            chunks.append(line.strip())
    if acc is not None:
        seqs[acc] = "".join(chunks)
    return seqs


def gene_names(gene, aliases):
    """Full set of query tokens for `gene` (gene + aliases)."""
    return [gene] + aliases.get(gene.lower(), [])


# ---------------------------------------------------------------------------
# Diversity selection
#
# Exact-sequence dedup alone leaves a reference set dominated by whichever
# allele happens to be commonest in NCBI. What round 1 needs is a set that
# SPANS the allelic range, so that whichever allele a genome carries, some
# reference matches it well. call_genes.py already picks the best
# (coverage, identity) reference per gene, so extra alleles cost nothing but
# search time.
#
# Representatives are chosen greedily by farthest-point: start from the
# RefSeq-preferred longest sequence, then repeatedly add whichever remaining
# sequence is least similar to everything already chosen. Similarity is k-mer
# Jaccard, which needs no alignment and no external tool (cd-hit is not in the
# pipeline image) and tracks protein identity closely enough for this purpose.
# ---------------------------------------------------------------------------

def kmer_set(seq, k=5):
    return {seq[i:i + k] for i in range(len(seq) - k + 1)} or {seq}


def jaccard(a, b):
    union = len(a | b)
    return (len(a & b) / union) if union else 0.0


def select_diverse(records, n_max, reserve_frac=0.34):
    """Up to n_max records spanning the sequence diversity of `records`.

    `records` must already be ordered by preference (RefSeq first, then
    longest).

    Two mechanisms, because farthest-point alone proved insufficient. It is a
    k-center heuristic -- it minimises the distance from any allele to its
    nearest representative -- but with few slots it favours OUTLIERS, and on
    LM088-RCB001 it discarded the near-exact references for agrA, sarA and
    sdrH, dropping round-1 identity from 99.2% to 84.5% for agrA. So:

      1. the first `reserve_frac` of the slots go to the top of the preference
         order (RefSeq, longest) -- the canonical alleles most likely to match
         a typical strain;
      2. the rest are filled by farthest-point, which covers the outliers.

    Raising n_max is the other half of the fix: k-center error shrinks as k
    grows, and BLAST against a few more short proteins costs nothing.
    """
    if len(records) <= n_max:
        return records

    sets = [kmer_set(r["sequence"]) for r in records]
    lens = [len(r["sequence"]) for r in records]

    def distance(i, j):
        """Combined sequence and LENGTH dissimilarity.

        Length is not a refinement here, it is essential. k-mer Jaccard barely
        moves when a protein is extended by a low-complexity repeat -- an
        SD-expanded SdrF shares nearly all its k-mers with a short allele,
        because 'SDSDSD...' contributes only a handful of distinct k-mers. Yet
        length is precisely what round-1 calling is sensitive to, since
        coverage is measured against reference length: a 1733 aa reference
        gives 66% coverage for a 783 aa allele and gets called 'partial'.

        Taking the max means a reference is kept if it differs in EITHER
        sequence or length.
        """
        seq_d = 1.0 - jaccard(sets[i], sets[j])
        len_d = 1.0 - (min(lens[i], lens[j]) / max(lens[i], lens[j], 1))
        return max(seq_d, len_d)

    n_reserved = min(max(1, int(round(n_max * reserve_frac))), n_max)
    chosen = list(range(n_reserved))      # canonical, by preference order

    while len(chosen) < n_max:
        best_i, best_dist = None, -1.0
        for i in range(len(records)):
            if i in chosen:
                continue
            # distance to the nearest already-chosen reference
            d = min(distance(i, c) for c in chosen)
            if d > best_dist:
                best_dist, best_i = d, i
        if best_i is None:
            break
        chosen.append(best_i)
    return [records[i] for i in sorted(chosen)]


def name_query(names):
    """OR-query over tokens; tokens with '[' or '"' are verbatim queries."""
    parts = [n if ("[" in n or '"' in n) else f"{n}[Gene Name]" for n in names]
    return "(" + " OR ".join(parts) + ")"


# ---------------------------------------------------------------------------
# Product filters: keep only records whose PRODUCT actually names the gene.
#
# The queries above are deliberately loose so that records carrying a product
# name and no gene symbol are retrieved. The cost is that a loose term also
# pulls in neighbours: '"penicillinase"[Title]' matches 'penicillinase
# repressor' (BlaI) and 'penicillinase regulatory protein' (BlaR1), so without
# a filter blaZ's reference set contains its own operon partners.
#
# That is not a cosmetic problem. On LM088-RCB001 those BlaI/BlaR1 references
# hit 2.7 kb away from blaZ, call_genes.py merged them into one cluster (gap
# 5 kb), and the "blaZ locus" became the entire bla operon -- 902,967-906,040
# instead of 902,967-903,773. Round 2 then saw a locus dominated by BlaR1 and
# flagged blaZ.
#
# Same accept/deny concept as finalize_calls.py, applied at build time.
# Records whose product marks them incomplete. Such a reference distorts
# coverage-first ranking far more than it helps.
INCOMPLETE_RE = re.compile(r",\s*partial\b|\btruncation\b|\btruncated\b"
                           r"|\bfragment\b", re.I)

# ---------------------------------------------------------------------------
# Per-gene organism restriction for the NCBI query.
#
# The default is genus-wide ('Staphylococcus[Organism]') because most of these
# genes are shared and a wider query gives a better allele panel.
#
# Override only where the genus-wide query returns essentially NOTHING from the
# species of interest, so that --prefer-organism has nothing to rank. Narrowing
# the query trades allelic breadth for species accuracy, so it is the wrong
# default -- it was narrowing to RefSeq that collapsed the panels to one allele
# per gene in the first place.
# ---------------------------------------------------------------------------
GENE_ORGANISM = {
    # Deliberately EMPTY. embp was tried here and it backfired: restricting to
    # S. epidermidis returned only 3 records, all of them ', partial', so after
    # INCOMPLETE_RE the panel was empty and a working call (divergent 42.2% at
    # the correct locus, nr-confirmed 99.6%) became no_reference. The fail-open
    # fallback did not help because it guards "no IDs returned", and 3 IDs were
    # returned -- they were lost later, during filtering.
    #
    # Lesson for anyone adding an entry here: a species restriction only pays
    # off when that species has ENOUGH COMPLETE records. Check first, e.g.
    #   esearch -db protein -query '<terms> AND "<species>"[Organism]' | esummary
    # and confirm the count and that the titles are not all partial. Otherwise
    # ALLOW_INCOMPLETE and --prefer-organism are the safer levers, since they
    # widen the panel rather than narrowing it.
}

# Genes exempt from INCOMPLETE_RE.
#
# INCOMPLETE_RE drops ', partial' / 'truncated' / 'fragment' records because a
# short reference wins outright under coverage-first ranking -- that is what
# moved the sdrF call to a 222 bp locus. But for proteins of several thousand
# residues, nearly every GenBank record is partial: all 3 S. epidermidis embp
# records are, and a 'partial' 5,000 aa Embp is far more useful than a
# full-length S. aureus Ebh at 42% identity.
#
# Exempting a gene is safe only because --min-len-frac still applies: fragments
# are then excluded by LENGTH relative to the gene's median, which is the
# property that actually matters, rather than by a word in the title.
# Emptied. The exemption was added for embp on the assumption that a ~10,200 aa
# protein has only partial records -- but that was true only under the
# S. epidermidis query restriction, which was itself reverted. Genus-wide, embp
# returns 20 FULL-LENGTH references (2839-10624 aa) and 0 partials, so the
# exemption does nothing for the gene it was written for.
#
# For bhp it was actively harmful. Bhp is ~2,500 aa, but its 'partial' records
# are SHORT FRAGMENTS (160-297 aa), not long partials. Admitting them collapsed
# the panel median from ~2,500 to 296 aa -- and --max-len-frac then evicted the
# genuine full-length references as "over-length" at 2.0x of a median set by
# junk. Two guards compounding to destroy the panel.
#
# Re-add a gene here only with evidence that ITS partials are long: check
# lengths of the ', partial' records before trusting the reasoning.
ALLOW_INCOMPLETE = set()

PRODUCT_FILTERS = {
    "agra": {"accept": r"\bagrA\b|accessory gene regulator (protein )?A\b"
                       r"|response regulator AgrA|quorum-?sensing response regulator",
             "deny":   r"\bagr[BCD]\b|accessory gene regulator (protein )?[BCD]\b"},
    "agrb": {"accept": r"\bagrB\b|accessory gene regulator (protein )?B\b",
             "deny":   r"\bagr[ACD]\b|accessory gene regulator (protein )?[ACD]\b"},
    "blai": {"accept": r"\bblaI\b|beta-?lactamase repressor|penicillinase repressor"
                       r"|BlaI/MecI/CopY",
             "deny":   r"\bblaZ\b|\bblaR1?\b|regulatory|sensor-?transducer"},
    "blaz": {"accept": r"\bblaZ\b|beta-?lactamase|penicillinase",
             # the operon partners, which 'penicillinase' alone would admit
             "deny":   r"\bblaI\b|\bblaR1?\b|repressor|regulator|regulatory"
                       r"|sensor-?transducer|metallo-?beta-?lactamase"
                       r"|class [CD] beta-?lactamase"},
    "fosb": {"accept": r"\bfosB\b|metallothiol transferase|bacillithiol transferase"
                       r"|fosfomycin resistance|\byfcC\b",
             "deny":   r"\bfos[ACX]\b|fos[ACX] family"},
    "sara": {"accept": r"\bsarA\b|accessory regulator A\b|regulator SarA",
             "deny":   r"\bsar[RSTUVXZ]\b|\brot\b"},
    "sdrf": {"accept": r"\bsdrF\b|serine-?aspartate repeat(-containing)? protein F\b",
             "deny":   r"\bsdr[CDEGH]\b|\bfbe\b|clumping factor|\bclf[AB]\b"
                       r"|B-like domain"},
    "sdrh": {"accept": r"\bsdrH\b|MSCRAMM-?like protein SdrH|SdrH family",
             "deny":   r"\bsdr[CDEFG]\b|\bfbe\b|clumping factor|\bclf[AB]\b"
                       r"|B-like domain"},

    # ── mecA ────────────────────────────────────────────────────────────────
    # 'adaptor' is the single most important token in this whole table: it is
    # what separates the ~240 aa ClpC adaptor MecA from the ~668 aa PBP2a. With
    # the adaptors denied the panel median becomes ~669 aa, so --min-len-frac
    # then independently excludes any 240 aa record that slips through.
    "meca": {"accept": r"\bmecA\b|PBP2a|PBP2'|penicillin-?binding protein 2a"
                       r"|beta-?lactam-?resistant peptidoglycan transpeptidase"
                       r"|methicillin resistance",
             "deny":   r"adaptor|adapter|\bClpC\b|competence|sporulation"
                       r"|\bmec[RI]1?\b|negative regulator of genetic competence"},

    # ── ica operon ──────────────────────────────────────────────────────────
    "icad": {"accept": r"\bicaD\b|intercellular adhesion protein D\b"
                       r"|poly-?beta-?1,6-?N-?acetyl|PIA synthesis protein IcaD",
             "deny":   r"\bica[ABCR]\b|intercellular adhesion protein [ABC]\b"},
    "icar": {"accept": r"\bicaR\b|ica operon|icaADBC|biofilm operon",
             # icaADBC appears in IcaR's own product name, so the deny list
             # must not use a bare 'icaA' substring -- \bica[ABCD]\b requires a
             # word boundary after the letter, which 'icaADBC' does not provide.
             "deny":   r"\bica[ABCD]\b|intercellular adhesion protein [ABCD]\b"},

    # ── IS256 ───────────────────────────────────────────────────────────────
    # Deny every OTHER numbered IS family; a bare 'IS[0-9]+ family' deny would
    # reject IS256 itself.
    "is256": {"accept": r"\bIS256\b|IS256[- ]?(like|family)",
              "deny":   r"\bIS(?!256)\d{1,4}\b"},

    # ── aacA (aacA-aphD, Tn4001) ────────────────────────────────────────────
    # Bifunctional, so APH tokens are legitimate here; other AAC classes aren't.
    # NCBI titles this enzyme both spelled out and abbreviated ('bifunctional
    # AAC/APH'); both forms must be accepted.
    "aaca": {"accept": r"\baacA\b|aacA-?aphD|AAC\(6'\)|APH\(2''?\)|AAC/APH"
                       r"|bifunctional\s+(AAC|APH)"
                       r"|bifunctional.*(acetyltransferase|phosphotransferase)"
                       r"|aminoglycoside.*(6').*acetyltransferase",
             "deny":   r"AAC\(3\)|AAC\(2'\)|\baadA\b|\bant\(|streptomycin"
                       r"|streptothricin"},

    # ── fdh ─────────────────────────────────────────────────────────────────
    # The accessory subunits FdhD/FdhE are short and would drag the median down;
    # 'fdh' untagged also matches zinc-dependent ALCOHOL dehydrogenases.
    "fdh": {"accept": r"\bfdh[A-C]?\b|formate dehydrogenase",
            "deny":   r"\bfdh[DE]\b|accessory|alcohol dehydrogenase"
                      r"|formaldehyde dehydrogenase|\badhB?\b"},

    # ── bhp (Bap homologue protein) ─────────────────────────────────────────
    # Deliberately does NOT accept a bare 'adhesin': S. epidermidis carries Aap,
    # Embp, SdrC-G and others, all titled 'adhesin', and they share repeat
    # domains with Bhp. On LM088-RCB001 a 582 bp locus matched an 'adhesin' at
    # 99.5% identity but only 18-21% coverage -- a shared repeat, not the gene.
    # 'cell wall associated biofilm protein' is the title NCBI actually uses for
    # Bhp in S. epidermidis -- 'biofilm-associated protein' alone matched none of
    # the 6 records in the panel and emptied the gene entirely.
    "bhp": {"accept": r"\bbhp\b|Bap[- ]?(homolog\w*|family)"
                      r"|(cell wall[- ]?associated )?biofilm[- ]?(associated )?protein",
            "deny":   r"\bsdr[CDEFGH]\b|\baap\b|\bembp\b|clumping factor"
                      r"|fibronectin|elastin|serine-?aspartate|\batl[E]?\b"},

    # ── ebh / embp: the SAME protein, two names ─────────────────────────────
    # The earlier deny of '\bembp\b' here was wrong. It was written to stop
    # Embp 'leaking' into the Ebh panel, on the assumption they were distinct
    # paralogues. They are orthologues of one ~10,200 aa ECM-binding giant, and
    # in S. epidermidis the Embp records ARE the correct ones -- denying them
    # left the panel S. aureus-dominated, which is why round 1 scored 48.2%
    # against a locus nr matched at 99.6%.
    #
    # Both entries are deliberately identical so either symbol works. List only
    # ONE of them in genes.txt; embp is preferred for S. epidermidis.
    "ebh": {"accept": r"\bebh\b|\bembp\b|Ebh-?like|hyperosmolarity resistance"
                      r"|extracellular matrix[- ]?binding protein",
            "deny":   r"\bsdr[CDEFGH]\b|\baap\b|\bbhp\b|clumping factor"},
    "embp": {"accept": r"\bembp\b|\bebh\b|Ebh-?like|hyperosmolarity resistance"
                       r"|extracellular matrix[- ]?binding protein",
             "deny":   r"\bsdr[CDEFGH]\b|\baap\b|\bbhp\b|clumping factor"},

    # ── msrA: TWO UNRELATED PROTEINS SHARE THIS SYMBOL ──────────────────────
    # 1. peptide-methionine (S)-S-oxide reductase MsrA (~170 aa), a universal
    #    methionine-repair enzyme
    # 2. ABC-F type ribosomal protection protein Msr(A) (~488 aa), the
    #    macrolide/streptogramin resistance determinant
    # On the 40-gene run the panel came back 19 housekeeping records to 1
    # resistance record, median 172 aa, and the call landed on a 171-codon
    # locus at 98.3% -- reported as msrA PRESENT for a genome carrying no
    # macrolide resistance gene at all. Identical in structure to mecA/ClpC.
    "msra": {"accept": r"\bmsr\(?A\)?\b|ABC-?F type ribosomal protection"
                       r"|macrolide.*(efflux|resistance)",
             "deny":   r"peptide-?methionine|methionine sulfoxide"
                       r"|S-?oxide reductase|\bmsrB\b"},

    # ── dfrC / dfrG ─────────────────────────────────────────────────────────
    # Every genome carries a native dihydrofolate reductase, and NCBI titles
    # the chromosomal S. epidermidis enzyme 'trimethoprim-SENSITIVE
    # dihydrofolate reductase DfrC'. So a dfrC 'present' call is the wild-type,
    # susceptible state -- NOT a resistance finding. Acquired dfrG is a
    # different gene: its panel cross-matched the native dfrC locus at 73.0%
    # (same coordinates, 1401771-1402247), correctly rejected by the 90%
    # threshold.
    "dfrc": {"accept": r"\bdfrC\b|dihydrofolate reductase",
             "deny":   r"\bdfr(S1|S2|[ABDGK])\b|\bfolA\b"},
    # ACCEPT REQUIRES THE NAME. 'trimethoprim-resistant dihydrofolate
    # reductase' alone also matches DfrA, DfrK and DfrS1 -- and DfrS1 is what
    # confirmed dfrG at 97-100% in 4 genomes whose locus is the native
    # trimethoprim-SENSITIVE dfrC: identical coordinates to the dfrC call,
    # round-1 identity 38.9-40.7%, confirming organism [Bacteria], ranks 7-22.
    "dfrg": {
        "accept": r"""\bdfrG\b""",
        "deny": r"""\bdfr(S1|S2|[ABCDK])\b|trimethoprim-?sensitive|\bfolA\b""",
    },

    # ── erm family: mutually ~60-80% identical ──────────────────────────────
    # Each must deny the other two, or one erm gene's panel confirms another's
    # locus. On LM088-RCB001 all three called ONE locus (2333913-2334713) at
    # 24-26% identity -- a native rRNA methyltransferase, not an erm gene.
    "erma": {"accept": r"\berm\(?A\)?\b|23S rRNA .*adenine.*methyltransferase"
                       r"|rRNA adenine N-6-methyltransferase",
             "deny":   r"\berm\(?[BCTFGY]\)?\b"},
    "ermb": {"accept": r"\berm\(?B\)?\b|23S rRNA .*adenine.*methyltransferase",
             "deny":   r"\berm\(?[ACTFGY]\)?\b"},
    "ermc": {"accept": r"\berm\(?C\)?\b|23S rRNA .*adenine.*methyltransferase",
             "deny":   r"\berm\(?[ABTFGY]\)?\b"},

    # fusB and fusC are ~45% identical to each other
    "fusb": {"accept": r"\bfusB\b|fusidic acid resistance",
             "deny":   r"\bfus[CD]\b"},
    "fusc": {"accept": r"\bfusC\b|fusidic acid resistance",
             "deny":   r"\bfus[BD]\b"},
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
             "deny":   r"\barc[BCD]\b|ornithine carbamoyltransferase"
                       r"|carbamate kinase|antiporter"},
    # 'catabolic' matters: ArcB is the CATABOLIC ornithine carbamoyltransferase,
    # and the anabolic isoenzyme ArgF is a native arginine-biosynthesis gene
    # present in every genome.
    "arcb": {"accept": r"\barcB\b|ornithine carbamoyltransferase",
             "deny":   r"\barc[ACD]\b|\bargF\b|anabolic|arginine deiminase"
                       r"|carbamate kinase"},
    "arcc": {"accept": r"\barcC\b|carbamate kinase",
             "deny":   r"\barc[ABD]\b|arginine deiminase"
                       r"|ornithine carbamoyltransferase"},
    "arcd": {"accept": r"\barcD\b|arginine[-/ ]ornithine antiporter",
             "deny":   r"\barc[ABC]\b|arginine deiminase|carbamate kinase"},
    "argr": {"accept": r"\bargR\b|\bahrC\b|AhrC/ArgR|arginine repressor",
             "deny":   r"\barc[ABCD]\b"},

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
              "deny":   r"permease|ATP[- ]?binding protein"},
    "opp3b": {"accept": r"\bopp-?3B\b|\boppB\b|oligopeptide ABC transporter permease",
              "deny":   r"ATP[- ]?binding|substrate-?binding"},
    "opp3c": {"accept": r"\bopp-?3C\b|\boppC\b|oligopeptide ABC transporter permease",
              "deny":   r"ATP[- ]?binding|substrate-?binding"},
    # NO accept gate, for the same reason as sesI: the records carry GENERIC
    # titles -- 'ABC transporter ATP-binding protein', 'ATP-binding cassette
    # domain-containing protein', even 'hypothetical protein'. Nothing in them
    # says 'opp' or 'oligopeptide', so any product gate rejects the very record
    # the accession query and the seed exist to supply. Retrieval is anchored on
    # ASM79718.1 / WP_002490823.1 plus data/opp3D_seed.faa, so precision comes
    # from the query, not from the title.
    "opp3d": {"deny": r"permease|substrate-?binding"},
    "opp3e": {"deny": r"permease|substrate-?binding"},

    # ── ACME kdp operon (types IV, V) ───────────────────────────────────────
    # UNVERIFIED: kdp was not in the 40-gene run, so these accept terms are the
    # standard NCBI product names but have not been checked against real output.
    # The per-gene counts in the next build log will show whether they worked.
    # Like arc, kdp is NATIVE in staphylococci -- presence is not ACME evidence,
    # which is why acme_type.py requires co-location.
    "kdpa": {"accept": r"\bkdpA\b|(potassium|K\(\+\))-?transporting ATPase.*subunit A",
             "deny":   r"\bkdp[BCDEF]\b|subunit [BCDF]\b"},
    "kdpb": {"accept": r"\bkdpB\b|(potassium|K\(\+\))-?transporting ATPase.*subunit B",
             "deny":   r"\bkdp[ACDEF]\b|subunit [ACDF]\b"},
    "kdpc": {"accept": r"\bkdpC\b|(potassium|K\(\+\))-?transporting ATPase.*subunit C",
             "deny":   r"\bkdp[ABDEF]\b|subunit [ABDF]\b"},
    "kdpd": {"accept": r"\bkdpD\b|sensor histidine kinase KdpD"
                       r"|potassium-?transporting ATPase KdpD",
             "deny":   r"\bkdp[ABCE]\b|response regulator"},
    "kdpe": {"accept": r"\bkdpE\b|response regulator KdpE|KdpE family",
             "deny":   r"\bkdp[ABCD]\b|histidine kinase"},

    # ── acquired resistance ─────────────────────────────────────────────────
    "aadd": {"accept": r"\baadD\b|ANT\(4'\)"
                       r"|kanamycin.*(nucleotidyl|adenylyl)transferase"
                       r"|aminoglycoside adenyl",
             "deny":   r"\baad[ABE]\b|ANT\(3|ANT\(6|ANT\(9|streptomycin"
                       r"|spectinomycin"},
    "mphc": {"accept": r"\bmph\(?C\)?\b|macrolide 2'-phosphotransferase",
             "deny":   r"\bmph\(?[ABDE]\)?\b"},
    # tetL is the close homologue; tetM/tetO are ribosomal protection proteins,
    # a different mechanism that should never confirm an efflux gene.
    "tetk": {"accept": r"\btetK\b|tetracycline efflux MFS"
                       r"|tetracycline resistance protein",
             "deny":   r"\btet[ABCDEGHLMOSTVWZ]\b|ribosomal protection"},
    # norB/norC/mdeA/mepA are the other staphylococcal MFS/MATE pumps
    "nora": {"accept": r"\bnorA\b|multidrug efflux MFS transporter NorA",
             "deny":   r"\bnor[BCDG]\b|\bmdeA\b|\bmepA\b|\bsepA\b|\btet[KL]\b"},
    "uhpt": {"accept": r"\buhpT\b|hexose-6-phosphate",
             "deny":   r"\buhp[ABC]\b|\bglpT\b"},

    # ── adhesins / biofilm / secretion ──────────────────────────────────────
    # Aae and AtlE are BOTH autolysins. 'autolysin/adhesin Aae' vs 'bifunctional
    # autolysin' is the only thing separating them, so each denies the other's
    # exact title.
    "aae":  {"accept": r"\baae\b|autolysin/adhesin",
             "deny":   r"\batlE?\b|bifunctional autolysin"},
    # NOTE: 'bifunctional autolysin' is also the title of S. aureus Atl, which
    # no pattern can distinguish from AtlE. Species preference in the panel is
    # what keeps this honest, not the title.
    "atle": {"accept": r"\batlE\b|bifunctional autolysin"
                       r"|N-acetylmuramoyl-L-alanine amidase|glucosaminidase",
             "deny":   r"\baae\b|autolysin/adhesin|\batlA\b"},
    # Aap shares B-repeats with Bhp, Embp and the Sdr family, and -seg no makes
    # that cross-matching more likely, so every one of them is denied.
    "aap":  {"accept": r"\baap\b|accumulation[- ]?associated protein",
             "deny":   r"\bbhp\b|\bsdr[CDEFGH]\b|\bembp\b|biofilm-?associated"
                       r"|clumping factor"},
    "agrc": {"accept": r"\bagrC\b|quorum-?sensing sensor histidine kinase"
                       r"|accessory gene regulator (protein )?C\b",
             "deny":   r"\bagr[ABD]\b|accessory gene regulator (protein )?[ABD]\b"},
    # AgrD is a ~46 aa propeptide; most records carry no symbol, only the
    # product name or the AIP type.
    "agrd": {"accept": r"\bagrD\b|cyclic lactone autoinducer|autoinducing peptide"
                       r"|AIP-I{1,3}V?|accessory gene regulator (protein )?D\b",
             "deny":   r"\bagr[ABC]\b|accessory gene regulator (protein )?[ABC]\b"},
    "icaa": {"accept": r"\bicaA\b|intercellular adhesion protein A\b"
                       r"|N-acetyl-D-glucosamine synthase",
             "deny":   r"\bica[BCDR]\b|intercellular adhesion protein [BCD]\b"
                       r"|N-deacetylase|biosynthesis/export"},
    "icab": {"accept": r"\bicaB\b|polysaccharide N-deacetylase"
                       r"|adhesin biosynthesis polysaccharide N-deacetyl",
             "deny":   r"\bica[ACDR]\b|intercellular adhesion protein [ACD]\b"
                       r"|biosynthesis/export|glucosamine synthase"},
    "icac": {"accept": r"\bicaC\b|biosynthesis/export protein"
                       r"|intercellular adhesion protein C\b",
             "deny":   r"\bica[ABDR]\b|intercellular adhesion protein [ABD]\b"
                       r"|N-deacetylase|glucosamine synthase"},
    # SdrG == Fbe. Do not list both symbols in genes.txt.
    "sdrg": {"accept": r"\bsdrG\b|\bfbe\b|fibrinogen-?binding adhesin"
                       r"|MSCRAMM family adhesin SdrG",
             "deny":   r"\bsdr[CDEFH]\b|clumping factor|\bclf[AB]\b"
                       r"|B-like domain"},
    # EssD records split between the T7SS name and a bare 'DNA/RNA non-specific
    # endonuclease', which is also the title of native nucleases. Accepting only
    # the secretion-system wording filters the native ones out.
    "essd": {"accept": r"\bessD\b|ESAT-6 secretion|type VII secretion",
             "deny":   r"\bess[ABC]\b|\besx[AB]\b|\besa[ABCD]\b"},
    # UNVERIFIED: sesI returned no records under its symbol on the 40-gene run,
    # so no real titles were observed. Needs ALIASES before it can work at all.
    # sesI is fetched by LOCUS TAG (SERP1654 / AAW54982.1), because the record
    # carries only a general protein annotation and the symbol never appears in
    # the title. A product-name accept gate would therefore reject the very
    # record the query was written to retrieve -- which is how bhp and opp3A
    # each emptied their panel. Retrieval is already precise, so only a deny
    # list is applied here: it costs nothing and still blocks the repeat-sharing
    # surface proteins if the query ever widens.
    "sesi": {"deny": r"\bses[ABCDEFGHJ]\b|\bsdr[CDEFGH]\b|\baap\b"
                     r"|\bembp\b|\bbhp\b|clumping factor"},

}


def product_ok(gene, product):
    """True if `product` names `gene` and not a paralog/operon partner."""
    f = PRODUCT_FILTERS.get(gene.lower())
    if not f:
        return True              # no curated filter: keep whatever was fetched
    t = product or ""
    if f.get("deny") and re.search(f["deny"], t, re.I):
        return False
    if f.get("accept") and not re.search(f["accept"], t, re.I):
        return False
    return True


# ---------------------------------------------------------------------------
# Predefined cassette databases
#
# SCCmec is a mobile cassette, not a gene: '<gene>[Gene Name]' retrieves
# nothing useful. Its reference proteins are the CDS translations of seven
# whole-cassette GenBank records (AB033763, D86934, AB037671, AB063172,
# AB121219, AF411935, FJ390057), prebuilt by fetch_sccmec_db.py and shipped as
# a file -- so this step needs no network call for SCCmec.
#
# Only the diagnostic core (mec and ccr complexes, tagged core=1 by the
# fetcher) is used by default. A full cassette carries dozens of transposases
# and hypothetical proteins that also occur elsewhere in staphylococcal genomes
# and on unrelated mobile elements; scoring those under the 'SCCmec' label
# would report the cassette present in any genome carrying a common
# transposase.
# ---------------------------------------------------------------------------
PREDEFINED = {
    "sccmec": "sccmec_proteins.faa",
    "is256":  "is256_proteins.faa",
}


def panel_dir_path(gene, panel_dir):
    """Locate a frozen panel FASTA for `gene` in --panel-dir, or None.

    Case-insensitive on the gene name, because a directory copied between
    macOS and Linux can arrive with either casing, and a silent miss here
    means the gene goes back to querying NCBI -- which is the one thing a
    frozen panel exists to prevent.
    """
    if not panel_dir or not os.path.isdir(panel_dir):
        return None
    for ext in (".fasta", ".faa", ".fa"):
        cand = os.path.join(panel_dir, gene + ext)
        if os.path.exists(cand):
            return cand
    want = gene.lower()
    for fn in sorted(os.listdir(panel_dir)):
        stem, ext = os.path.splitext(fn)
        if ext.lower() in (".fasta", ".faa", ".fa") and stem.lower() == want:
            return os.path.join(panel_dir, fn)
    return None


def predefined_path(gene, explicit, script_dir):
    """Locate the predefined FASTA for `gene`, or None."""
    if explicit:
        return explicit
    fname = PREDEFINED.get(gene.lower())
    if not fname:
        return None
    for cand in (os.path.join(script_dir, "..", "data", fname),
                 os.path.join(script_dir, "data", fname),
                 os.path.join(os.getcwd(), "data", fname),
                 os.path.join(os.getcwd(), fname)):
        cand = os.path.normpath(cand)
        if os.path.exists(cand):
            return cand
    return None


def load_predefined(path, gene, core_only, min_len=100):
    """Read a predefined FASTA into build_protein_db's record shape.

    min_len drops short CDS fragments. Whole-cassette GenBank records contain
    partial/interrupted copies of the core genes, and because score_cluster()
    ranks COVERAGE before identity, the shortest reference in the panel wins:
    on LM088-RCB001 a 67 aa mecI fragment (AB037671_cds130) reached 100%
    coverage and dragged the SCCmec call onto the bla operon 107 bp downstream
    of blaZ -- MecI being a homologue of BlaI. An ABSOLUTE floor is used rather
    than the median-relative --min-len-frac because a cassette panel holds
    genuinely different proteins (MecI ~123 aa, CcrB ~383, PBP2a ~668), so a
    median-relative floor would discard real short components.

    INCOMPLETE_RE is deliberately NOT applied here: for a cassette, 'truncated'
    is not a defect but a diagnostic feature -- a truncated mecR1 (dmecR1) is
    what distinguishes SCCmec types III and IV.
    """
    recs, hdr, chunks = [], None, []

    def flush():
        if hdr is None:
            return
        seq = "".join(chunks)
        if not seq:
            return
        parts = hdr.split("|")
        acc = parts[0]
        product = parts[3] if len(parts) > 3 else ""
        organism = parts[2] if len(parts) > 2 else ""
        m = re.search(r"\[core=([01])\]", hdr)
        is_core = (m.group(1) == "1") if m else True
        if core_only and not is_core:
            return
        if min_len and len(seq) < min_len:
            print(f"[PREDEF] {gene}: dropped {acc} ({len(seq)} aa < {min_len} aa "
                  f"floor) -- {product[:48]}")
            return
        recs.append({
            "uid": "", "accession": acc, "organism": organism,
            "product": re.sub(r"\s*\[(core|len|src)=[^\]]*\]", "", product).strip(),
            "length": len(seq),
            "is_refseq": False,
            "sequence": seq,
        })

    for line in open(path):
        line = line.rstrip("\n")
        if line.startswith(">"):
            flush()
            hdr, chunks = line[1:], []
        elif hdr is not None:
            chunks.append(line.strip())
    flush()
    return recs


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--genes", nargs="+", default=DEFAULT_GENES,
                    help="Gene names to retrieve (default: 8 S. epidermidis genes)")
    ap.add_argument("--out", required=True, help="Output directory")
    ap.add_argument("--retmax", type=int, default=200,
                    help="Max NCBI records fetched per gene before diversity "
                         "selection (default 200). The old default of 30, "
                         "combined with a RefSeq-only query, produced "
                         "one-allele reference sets.")
    ap.add_argument("--min-len-frac", type=float, default=0.5,
                    help="drop references shorter than this fraction of the "
                         "gene's median reference length (default 0.5; 0 "
                         "disables). Guards against a fragment hijacking the "
                         "call, since call_genes.py ranks coverage first.")
    ap.add_argument("--max-len-frac", type=float, default=0.0,
                    help="drop references LONGER than this multiple of the "
                         "gene's median reference length (0 = OFF, the "
                         "default). Disabled: it never fired where it was aimed "
                         "(icaC's bad reference is 1.79x the median, under a 2.0x "
                         "ceiling) and it DID fire on bhp, where a median "
                         "corrupted by short partials made it evict the real "
                         "full-length references. orf_span() in call_genes.py "
                         "fixes the neighbour-merging problem properly. An "
                         "gene and its neighbour, inflating the reported locus "
                         "and making adjacent genes look like a frameshift.")
    ap.add_argument("--prefer-organism", default="Staphylococcus epidermidis",
                    help="rank references from this organism first (substring "
                         "match on the NCBI organism field). The gene query is "
                         "genus-wide, so without this the panel is dominated "
                         "by whichever species has most records -- usually "
                         "S. aureus. '' disables.")
    ap.add_argument("--max-alleles", type=int, default=20,
                    help="Max reference proteins KEPT per gene after dedup and "
                         "diversity selection (default 20). Raised from 8 "
                         "after 8 slots proved too few: the near-exact "
                         "reference was dropped for 3 of 8 genes, costing "
                         "agrA 14.7 identity points in round 1.")
    ap.add_argument("--sccmec-fasta", default="",
                    help="predefined SCCmec protein FASTA (default: look for "
                         "data/sccmec_proteins.faa next to this script). Build "
                         "it once with fetch_sccmec_db.py.")
    ap.add_argument("--is256-fasta", default="",
                    help="predefined IS256 transposase FASTA (default: look "
                         "for data/is256_proteins.faa next to this script). "
                         "Build it with fetch_predefined_db.py --db is256.")
    ap.add_argument("--panel-dir", default="",
                    help="directory of FROZEN per-gene reference FASTAs "
                         "(<gene>.fasta, as written to protein_db/ by a "
                         "previous run). Any gene found here is built from "
                         "that file and NOTHING is fetched from NCBI, so the "
                         "panel -- and therefore the calls -- reproduce "
                         "exactly. Genes absent from the directory fall "
                         "through to the normal query. Make one with "
                         "tools/freeze_panel.sh.")
    ap.add_argument("--predefined-max", type=int, default=0,
                    help="cap on proteins kept from a PREDEFINED database "
                         "(0 = no cap, the default). --max-alleles applies "
                         "only to genes fetched by name, not to curated "
                         "predefined sets.")
    ap.add_argument("--predefined-min-len", type=int, default=100,
                    help="drop proteins shorter than this from a PREDEFINED "
                         "database (default 100 aa; 0 disables). Whole-cassette "
                         "GenBank records carry 36-67 aa CDS fragments of the "
                         "core genes, and since call_genes.py ranks coverage "
                         "before identity the shortest reference wins: a 67 aa "
                         "mecI fragment put the SCCmec call on the bla operon.")
    ap.add_argument("--predefined-all", "--sccmec-all", action="store_true",
                    dest="predefined_all",
                    help="use EVERY protein in a predefined database, "
                         "including core=0 entries. Default is core=1 only. "
                         "For SCCmec that means all cassette proteins rather "
                         "than just the mec/ccr complexes -- which raises "
                         "false positives, since cassettes carry transposases "
                         "and hypothetical proteins found all over "
                         "staphylococcal genomes. (--sccmec-all is the old "
                         "name for this flag.)")
    ap.add_argument("--refseq-only", action="store_true",
                    help="Restrict to RefSeq (the old behaviour). RefSeq WP_ "
                         "records are non-redundant representatives, so they "
                         "COLLAPSE allelic diversity -- which is exactly what "
                         "a panel needs to keep. Off by default; RefSeq is "
                         "still ranked first.")
    ap.add_argument("--api-key", default=os.environ.get("NCBI_API_KEY", ""),
                    help="NCBI API key (or set NCBI_API_KEY env var)")
    ap.add_argument("--alias", action="append", default=[],
                    metavar="GENE=ALIAS1,ALIAS2",
                    help="Extra gene-name tokens to query for GENE "
                         "(repeatable). Built-in: aacA=aacA-aphD")
    ap.add_argument("--seed-dir", default="",
                    help="directory searched for <gene>_seed.faa files "
                         "(default: data/, data/seeds/ and the working "
                         "directory). Seeds are added on top of the fetched "
                         "panel and are never trimmed by --max-alleles.")
    ap.add_argument("--extra-fasta", action="append", default=[],
                    metavar="GENE=PATH",
                    help="Add user-supplied reference proteins for GENE from "
                         "a FASTA (headers: >acc description). Use for genes "
                         "NCBI gene-name search cannot find. Repeatable.")
    args = ap.parse_args()

    # ── normalise gene tokens ───────────────────────────────────────────────
    # Every gene name reaches the filesystem as <gene>.fasta here and as
    # <gene>_locus.fasta in extract_loci.py, and both are handled by shell code
    # in main.nf. A token carrying whitespace or shell metacharacters therefore
    # breaks the run some processes later, far from the cause. genes_3.txt
    # arrived with six trailing spaces ('sdrG ', 'erm(A) ') which mapfile
    # preserved, yielding 'sdrG .fasta'.
    cleaned, dropped = [], []
    for g in args.genes:
        g2 = g.strip()
        if not g2:
            continue
        if g2 != g:
            print(f"[NOTE] gene {g!r}: stripped surrounding whitespace -> {g2!r}")
        if re.search(r"[^\w.\-]", g2):
            print(f"[WARN] gene {g2!r} contains characters that are unsafe in a "
                  f"filename or an Entrez query. Rename it in genes.txt (e.g. "
                  f"'erm(A)' -> 'ermA') and add the official spelling to "
                  f"ALIASES as a QUOTED query term.")
            dropped.append(g2)
            continue
        cleaned.append(g2)
    if dropped:
        print(f"[WARN] skipping {len(dropped)} unusable gene name(s): "
              f"{', '.join(dropped)}. They are reported as no_reference, "
              f"never as absent.")
    # dedup, preserving order -- a repeated gene builds the panel twice and
    # emits two rows for one locus (genes_3.txt listed both fbe and sdrG, which
    # are the same protein)
    seen_g, uniq = set(), []
    for g in cleaned:
        if g.lower() in seen_g:
            print(f"[NOTE] gene {g!r} listed more than once; keeping one copy")
            continue
        seen_g.add(g.lower())
        uniq.append(g)
    args.genes = uniq

    # ── auto-discover seed FASTAs ───────────────────────────────────────────
    # data/<gene>_seed.faa is picked up without a per-gene --extra-fasta flag,
    # so pipeline runs honour seeds too: Nextflow stages the predefined
    # directory as ./data, which is one of the locations searched below.
    #
    # Seeds are APPENDED after diversity selection, so they always survive --
    # that is the point of seeding a known-matching allele.
    seed_dirs = [args.seed_dir] if args.seed_dir else [
        os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data"),
        os.path.join(os.getcwd(), "data"),
        os.path.join(os.getcwd(), "data", "seeds"),
        os.getcwd(),
    ]
    already = {spec.partition("=")[0].strip().lower() for spec in args.extra_fasta}
    for gene in args.genes:
        if gene.lower() in already:
            continue
        for d in seed_dirs:
            cand = os.path.normpath(os.path.join(d, f"{gene}_seed.faa"))
            if os.path.exists(cand):
                args.extra_fasta.append(f"{gene}={cand}")
                print(f"[SEED] {gene}: {cand}")
                break

    aliases = dict(ALIASES)
    for spec in args.alias:
        g, _, al = spec.partition("=")
        toks = [a.strip() for a in al.split(",") if a.strip()]
        if g.strip() and toks:
            aliases[g.strip().lower()] = toks

    os.makedirs(args.out, exist_ok=True)
    sleep = 0.15 if args.api_key else 0.5  # 10 req/s with key, 3/s without
    manifest_rows = []
    written = set()

    script_dir = os.path.dirname(os.path.abspath(__file__))

    # Provenance for a frozen panel. The FASTA defline is the only thing a
    # freeze preserves, and for SEEDED references it carries no organism or
    # product at all ('>ASM79718.1|opp3D|| [len=255]') -- those fields only
    # ever existed in the manifest. Reloading from the FASTA alone would
    # therefore blank them, losing the record of where a reference came from.
    # freeze_panel.sh copies the old manifest alongside for exactly this.
    frozen_meta = {}
    if args.panel_dir:
        prov = os.path.join(args.panel_dir, "manifest.csv.provenance")
        if os.path.exists(prov):
            with open(prov) as fh:
                for r in csv.DictReader(fh):
                    frozen_meta[(r["gene"].lower(), r["accession"])] = r
            print(f"[NOTE] frozen panel: read provenance for "
                  f"{len(frozen_meta)} reference(s)")

    for gene in args.genes:
        # ── frozen panel (--panel-dir) ──────────────────────────────────────
        # Takes precedence over everything, including the cassette registry:
        # the point of a frozen panel is that NOTHING is fetched, so the run
        # is reproducible.
        #
        # Why this exists. The gene query is genus-wide and NCBI's result set
        # drifts, so rebuilding resamples the panel: between two runs three
        # weeks apart, 16 of 55 panels drew different references (fdh swapped
        # 7 of 20, aap 5, kdpB 4) and 13 calls moved for that reason alone,
        # none of them from a code change. Re-running was therefore not a way
        # to confirm a result -- it produced a new one that had to be audited
        # from scratch. Freeze the panel and a re-run reproduces the table
        # exactly, so any difference is attributable to the code.
        #
        # The frozen files are just the previous run's protein_db/<gene>.fasta,
        # which already carry the 'acc|gene|organism|product' deflines
        # load_predefined() parses. See tools/freeze_panel.sh.
        frozen = panel_dir_path(gene, args.panel_dir) if args.panel_dir else None
        if frozen:
            # min_len=0, NOT --predefined-min-len. That floor is 100 aa and
            # exists to drop CDS fragments from whole-cassette GenBank records;
            # applied to a frozen allele panel it would delete all 20 agrD
            # references (44-99 aa -- it is a 46-codon propeptide) and 6 of 20
            # for icaD, silently turning two present genes into no_reference.
            kept = load_predefined(frozen, gene, core_only=False, min_len=0)
            if not kept:
                print(f"[WARN] {gene}: {frozen} yielded no proteins; "
                      f"reported as no_reference, not absent")
                continue
            # NOT sorted. The frozen file is reproduced in its original order,
            # byte for byte, so a frozen rebuild is identical to the run it
            # came from rather than merely equivalent. Sorting by length here
            # left all 55 panels with the same accessions in a different
            # order, which is the kind of "shouldn't matter" that has to be
            # reasoned about -- score_cluster() keeps the FIRST reference on a
            # (coverage, identity) tie -- instead of simply not arising.
            with open(os.path.join(args.out, f"{gene}.fasta"), "w") as fh:
                for r in kept:
                    fh.write(f">{r['accession']}|{gene}|{r['organism']}|"
                             f"{r['product']} [len={r['length']}]\n")
                    for i in range(0, len(r["sequence"]), 60):
                        fh.write(r["sequence"][i:i + 60] + "\n")
            written.add(gene)
            for r in kept:
                row = {"gene": gene,
                       **{k: v for k, v in r.items() if k != "sequence"}}
                old = frozen_meta.get((gene.lower(), r["accession"]))
                if old:
                    # Restore only what the defline could not carry; the
                    # sequence and its length come from the FASTA.
                    for k in ("uid", "organism", "product", "is_refseq"):
                        if not row.get(k) and old.get(k):
                            row[k] = old[k]
                manifest_rows.append(row)
            lens = sorted(r["length"] for r in kept)
            print(f"{gene}: FROZEN panel {os.path.basename(frozen)}, "
                  f"{len(kept)} references ({lens[0]}-{lens[-1]} aa)")
            continue

        # ── predefined cassette databases (SCCmec) ───────────────────────────
        explicit = (args.sccmec_fasta if gene.lower() == 'sccmec'
                    else args.is256_fasta if gene.lower() == 'is256' else '')
        pre = predefined_path(gene, explicit, script_dir)
        if gene.lower() in PREDEFINED:
            if not pre:
                print(f"[WARN] {gene}: no predefined FASTA found "
                      f"(expected data/{PREDEFINED.get(gene.lower(), '?')}). "
                      f"Build it with fetch_sccmec_db.py; reported as "
                      f"no_reference, not absent")
                continue
            kept = load_predefined(pre, gene, core_only=not args.predefined_all,
                                   min_len=args.predefined_min_len)
            n_total = len(load_predefined(pre, gene, core_only=False, min_len=0))
            n_dropped = n_total - len(kept)
            if n_dropped:
                print(f"[NOTE] {gene}: {n_dropped} of {n_total} proteins are "
                      f"core=0 and were excluded; --predefined-all includes them")
            if not kept:
                print(f"[WARN] {gene}: {pre} yielded no usable proteins")
                continue
            kept.sort(key=lambda x: -x["length"])
            n_unique = len(kept)
            # --max-alleles does NOT apply here by default. It exists to trim
            # hundreds of near-identical NCBI records down to a diverse set; a
            # predefined database is already a curated list, and silently
            # dropping most of a cassette would defeat --predefined-all.
            if args.predefined_max and n_unique > args.predefined_max:
                kept = select_diverse(kept, args.predefined_max)
            with open(os.path.join(args.out, f"{gene}.fasta"), "w") as fh:
                for r in kept:
                    fh.write(f">{r['accession']}|{gene}|{r['organism']}|"
                             f"{r['product']} [len={r['length']}]\n")
                    for i in range(0, len(r["sequence"]), 60):
                        fh.write(r["sequence"][i:i + 60] + "\n")
            written.add(gene)
            for r in kept:
                manifest_rows.append({"gene": gene,
                                      **{k: v for k, v in r.items()
                                         if k != "sequence"}})
            scope = ("all proteins" if args.predefined_all or gene.lower() != "sccmec"
                     else "mec/ccr core")
            lens = sorted(r["length"] for r in kept)
            print(f"{gene}: predefined {os.path.basename(pre)} ({scope}), "
                  f"{n_unique} available, {len(kept)} kept "
                  f"(lengths {lens[0]}-{lens[-1]} aa)")
            continue

        names = gene_names(gene, aliases)
        name_q = name_query(names)
        refseq_filter = " AND srcdb_refseq[PROP]" if args.refseq_only else ""
        # Per-gene organism restriction. --prefer-organism only RANKS results;
        # it cannot help when the genus-wide query returns no records from the
        # target species at all. embp is that case: querying by product name
        # returned 20 records, all S. aureus Ebh (median 6,337 aa), so round 1
        # scored the S. epidermidis locus at 42.2% -- the between-species
        # distance -- for a gene nr confirms at 99.6%.
        org = GENE_ORGANISM.get(gene.lower(), "Staphylococcus")
        term = f'{name_q} AND "{org}"[Organism]{refseq_filter}'
        if org != "Staphylococcus":
            print(f"[NOTE] {gene}: query restricted to {org} "
                  f"(see GENE_ORGANISM)")
        # Per-gene try/except: a network failure on ONE gene used to abort the
        # whole build, discarding every gene already fetched AND the seed
        # FASTAs, which are written at the end. A gene that cannot be fetched
        # is reported as no_reference by call_genes.py -- the same treatment as
        # a gene with no records -- which is never confused with 'absent'.
        try:
            ids = esearch(term, args.retmax, args.api_key)
        except Exception as e:                                   # noqa: BLE001
            print(f"[WARN] {gene}: NCBI query failed ({e}); reported as "
                  f"no_reference, not absent")
            continue
        time.sleep(sleep)

        # Organism restriction must FAIL OPEN. Narrowing embp to
        # S. epidermidis returned zero records -- the species' Embp is not
        # titled like the S. aureus Ebh records the product query was written
        # for -- turning a usable call (divergent, 42.2%, correct locus,
        # nr-confirmed at 99.6%) into no_reference. A genus-wide panel with a
        # between-species identity gap is strictly more informative than no
        # panel at all, so fall back rather than lose the gene.
        if not ids and org != "Staphylococcus":
            term = f'{name_q} AND "Staphylococcus"[Organism]{refseq_filter}'
            try:
                ids = esearch(term, args.retmax, args.api_key)
            except Exception as e:                               # noqa: BLE001
                print(f"[WARN] {gene}: genus-wide retry failed ({e}); "
                      f"reported as no_reference")
                continue
            time.sleep(sleep)
            print(f"[NOTE] {gene}: no records under {org}; fell back to "
                  f"genus-wide ({len(ids)} records). Round-1 identity will "
                  f"reflect the between-species distance, so read a low "
                  f"identity here as 'different species', not 'different gene'.")
        if not ids and args.refseq_only:
            # Resistance genes (aacA, ant, erm, ...) are plasmid/transposon
            # borne and often missing from RefSeq gene-name records; retry
            # without the RefSeq restriction before giving up.
            term = f"{name_q} AND Staphylococcus[Organism]"
            try:
                ids = esearch(term, args.retmax, args.api_key)
            except Exception as e:                               # noqa: BLE001
                print(f"[WARN] {gene}: NCBI retry failed ({e}); reported as "
                      f"no_reference")
                continue
            time.sleep(sleep)
            if ids:
                print(f"[NOTE] {gene}: no RefSeq records under "
                      f"{'/'.join(names)}; using {len(ids)} non-RefSeq "
                      f"(GenBank) records")
            else:
                print(f"[WARN] no records for {gene} (will be reported as "
                      f"no_reference, not absent)")
                continue

        # An empty search result was previously unguarded unless --refseq-only
        # was set, so esummary([]) was called, NCBI returned a payload with no
        # 'result' key, and the gene was reported as "esummary failed
        # ('result')" -- a message that reads like an API outage when the real
        # answer is simply that the query matched nothing. opp3D and opp3E hit
        # this on every run.
        if not ids:
            print(f"[WARN] {gene}: no records matched "
                  f"{name_q} AND \"{org}\"[Organism]. Reported as "
                  f"no_reference, not absent -- the gene was never searched "
                  f"for. Add a product-name or locus-tag alias (see ALIASES; "
                  f"sesI needed SERP1654/AAW54982 because its record carries "
                  f"only a general annotation).")
            continue

        try:
            summary = esummary(ids, args.api_key)
        except Exception as e:                                   # noqa: BLE001
            print(f"[WARN] {gene}: esummary failed ({e}); reported as "
                  f"no_reference")
            continue
        time.sleep(sleep)
        records = []
        for uid in ids:
            s = summary.get(uid, {})
            acc = s.get("caption", "")
            records.append({
                "uid": uid, "accession": acc,
                "organism": s.get("organism", ""),
                "product": s.get("title", ""),
                "length": s.get("slen", 0),
                "is_refseq": acc.startswith("WP_") or acc.startswith("NP_"),
            })

        try:
            seqs = efetch_fasta([r["uid"] for r in records], args.api_key)
        except Exception as e:                                   # noqa: BLE001
            print(f"[WARN] {gene}: efetch failed ({e}); reported as "
                  f"no_reference")
            continue
        time.sleep(sleep)

        # per-gene dedup by sequence hash; drop ', partial'.
        #
        # Preference order, most significant first:
        #   1. the organism the samples come from
        #   2. RefSeq over GenBank
        #   3. longest
        #
        # Organism comes FIRST because the query is genus-wide
        # ('AND Staphylococcus[Organism]'), so broadening it from RefSeq-only
        # flooded the sets with S. aureus. AgrA differs ~15% between S. aureus
        # and S. epidermidis, which is exactly the 99.2% -> 84.5% round-1 drop
        # observed on LM088-RCB001 when the panel went genus-wide.
        prefer = args.prefer_organism.lower()

        def pref_key(x):
            org = (x.get("organism") or "").lower()
            return (not (prefer and prefer in org), not x["is_refseq"], -x["length"])

        seen, seen_acc, kept = set(), set(), []
        n_offtarget = n_incomplete = n_short = 0
        for r in sorted(records, key=pref_key):
            seq = seqs.get(r["accession"])
            if not seq or r["accession"] in seen_acc:
                continue
            # ', partial' | 'truncation' | 'truncated' | 'fragment':
            # incomplete records. A truncated SdrF reference reached the panel
            # and, because score_cluster() ranks (coverage, identity) with
            # COVERAGE FIRST, a short reference trivially wins -- it moved the
            # called sdrF locus from the real gene at 52 kb to a 222 bp hit at
            # 2.43 Mb.
            if (INCOMPLETE_RE.search(r["product"])
                    and gene.lower() not in ALLOW_INCOMPLETE):
                n_incomplete += 1
                continue
            if not product_ok(gene, r["product"]):
                n_offtarget += 1
                continue
            h = hashlib.sha256(seq.encode()).hexdigest()
            if h in seen:
                continue
            seen.add(h)
            seen_acc.add(r["accession"])
            r["sequence"] = seq
            kept.append(r)

        # Relative length floor. Coverage-first ranking means a short
        # reference always outscores a long one, so a fragment that slipped
        # through the product filter would hijack the call. The floor is
        # relative to the MEDIAN, not the max, because some of these genes are
        # genuinely length-polymorphic: SdrF runs 783-1733 aa, and a
        # max-relative floor would discard the short allele that is often the
        # correct one.
        # Upper bound, symmetric with --min-len-frac. An OVER-length reference
        # is as damaging as a short one, just differently: it aligns across the
        # target gene AND its neighbour, which extends the genomic cluster and
        # inflates the reported locus. icaC's panel carried CAI6018687, a 635 aa
        # 'intercellular adhesion protein C [S. pseudintermedius]' against a
        # 354 aa median (1.80x) -- a mis-annotated gene model. It never won a
        # single call, but its HSPs stretched the icaC cluster to 1,921 bp for a
        # 355 aa protein, spanning into the adjacent ica gene in a different
        # reading frame, and the frame check then reported a frameshift in 9
        # genomes. There was no frameshift.
        if kept and args.max_len_frac > 0:
            lens = sorted(len(r["sequence"]) for r in kept)
            median = lens[len(lens) // 2]
            ceiling = median * args.max_len_frac
            before = len(kept)
            kept = [r for r in kept if len(r["sequence"]) <= ceiling]
            n_long = before - len(kept)
            if n_long:
                print(f"[NOTE] {gene}: dropped {n_long} reference(s) over "
                      f"{ceiling:.0f} aa ({args.max_len_frac:.1f}x of the "
                      f"{median} aa median) -- an over-length reference merges "
                      f"the gene with its neighbour")

        if kept and args.min_len_frac > 0:
            lens = sorted(len(r["sequence"]) for r in kept)
            median = lens[len(lens) // 2]
            floor = median * args.min_len_frac
            before = len(kept)
            kept = [r for r in kept if len(r["sequence"]) >= floor]
            n_short = before - len(kept)
            if n_short:
                print(f"[NOTE] {gene}: dropped {n_short} reference(s) under "
                      f"{floor:.0f} aa ({args.min_len_frac:.0%} of the "
                      f"{median} aa median)")

        bits = []
        if n_offtarget:  bits.append(f"{n_offtarget} off-target")
        if n_incomplete: bits.append(f"{n_incomplete} incomplete")
        if n_short:      bits.append(f"{n_short} too short")
        off = (", " + ", ".join(bits)) if bits else ""

        # An over-strict PRODUCT_FILTERS entry can reject every record, and the
        # length summary below then raised IndexError on an empty list, killing
        # the WHOLE build over one bad pattern -- the same class of failure the
        # per-gene try/except around the NCBI calls exists to prevent. bhp did
        # exactly this: NCBI titles it 'cell wall associated biofilm protein',
        # the filter expected 'biofilm-associated protein', all 6 records were
        # rejected and BUILD_PROTEIN_DB died after five genes had succeeded.
        #
        # Checked BEFORE the FASTA is written so no empty <gene>.fasta is left
        # behind and the gene is not added to `written`. With no manifest rows,
        # call_genes.py reports it as no_reference -- never as 'absent'.
        if not kept:
            print(f"[WARN] {gene}: {len(records)} fetched but 0 survived "
                  f"filtering{off}. Reported as no_reference, not absent. "
                  f"Check PRODUCT_FILTERS['{gene.lower()}'] against the real "
                  f"NCBI product names before trusting this.")
            continue

        n_unique = len(kept)
        kept = select_diverse(kept, args.max_alleles)

        with open(os.path.join(args.out, f"{gene}.fasta"), "w") as fh:
            for r in kept:
                fh.write(f">{r['accession']}|{gene}|{r['organism']}|{r['product']}"
                         f" [len={r['length']}]\n")
                for i in range(0, len(r["sequence"]), 60):
                    fh.write(r["sequence"][i:i + 60] + "\n")
        written.add(gene)

        for r in kept:
            manifest_rows.append({"gene": gene,
                                  **{k: v for k, v in r.items() if k != "sequence"}})
        lens = sorted(len(r["sequence"]) for r in kept)
        print(f"{gene}: {len(records)} fetched, {n_unique} unique{off}, "
              f"{len(kept)} kept  (lengths {lens[0]}-{lens[-1]} aa)")

    # user-supplied reference proteins (escape hatch for genes that NCBI
    # gene-name search cannot find at all)
    for spec in args.extra_fasta:
        g, _, path = spec.partition("=")
        g, path = g.strip(), path.strip()
        if not g or not path or not os.path.exists(path):
            print(f"[WARN] ignoring --extra-fasta {spec!r} (need gene=path "
                  f"with an existing file)")
            continue
        # A frozen panel is complete, and its seeds are already in it -- they
        # were appended when the run being frozen was built. Appending them
        # again would duplicate every seeded reference (agrB, sdrF, sdrH,
        # opp3D, opp3E, sesI), inflate the panel and skew the median that
        # --min-len-frac is computed from.
        if args.panel_dir and panel_dir_path(g, args.panel_dir):
            print(f"[NOTE] {g}: seed {os.path.basename(path)} not appended "
                  f"-- the frozen panel already contains it")
            continue
        seqs, acc, chunks = {}, None, []
        for line in open(path):
            if line.startswith(">"):
                if acc is not None:
                    seqs[acc] = "".join(chunks)
                head = line[1:].strip()
                acc = head.split()[0] if head else "unknown"
                chunks = []
            else:
                chunks.append(line.strip())
        if acc is not None:
            seqs[acc] = "".join(chunks)

        mode = "a" if g in written else "w"
        n_added = 0
        with open(os.path.join(args.out, f"{g}.fasta"), mode) as fh:
            for acc, seq in seqs.items():
                if not seq:
                    continue
                fh.write(f">{acc}|{g}|| [len={len(seq)}]\n")
                for i in range(0, len(seq), 60):
                    fh.write(seq[i:i + 60] + "\n")
                manifest_rows.append({"uid": "", "accession": acc, "organism": "",
                                      "product": "user-supplied reference",
                                      "length": len(seq), "is_refseq": False,
                                      "gene": g})
                n_added += 1
        written.add(g)
        print(f"{g}: +{n_added} user-supplied reference proteins from {path}")

    # combined FASTA for BLAST
    all_genes = list(args.genes) + [s.partition("=")[0].strip()
                                    for s in args.extra_fasta]
    seen_genes = set()
    with open(os.path.join(args.out, "combined.fasta"), "w") as out:
        for gene in all_genes:
            if gene in seen_genes:
                continue
            seen_genes.add(gene)
            path = os.path.join(args.out, f"{gene}.fasta")
            if os.path.exists(path):
                out.write(open(path).read())

    if manifest_rows:
        with open(os.path.join(args.out, "manifest.csv"), "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(manifest_rows[0].keys()))
            w.writeheader()
            w.writerows(manifest_rows)
    print(f"Total manifest rows: {len(manifest_rows)} -> {args.out}/manifest.csv")


if __name__ == "__main__":
    main()

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


def load_predefined(path, gene, core_only):
    """Read a predefined FASTA into build_protein_db's record shape."""
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
    ap.add_argument("--predefined-max", type=int, default=0,
                    help="cap on proteins kept from a PREDEFINED database "
                         "(0 = no cap, the default). --max-alleles applies "
                         "only to genes fetched by name, not to curated "
                         "predefined sets.")
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

    for gene in args.genes:
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
            kept = load_predefined(pre, gene, core_only=not args.predefined_all)
            n_total = len(load_predefined(pre, gene, core_only=False))
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
        term = f"{name_q} AND Staphylococcus[Organism]{refseq_filter}"
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
            if INCOMPLETE_RE.search(r["product"]):
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
        bits = []
        if n_offtarget:  bits.append(f"{n_offtarget} off-target")
        if n_incomplete: bits.append(f"{n_incomplete} incomplete")
        if n_short:      bits.append(f"{n_short} too short")
        off = (", " + ", ".join(bits)) if bits else ""
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

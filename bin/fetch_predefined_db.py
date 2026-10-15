#!/usr/bin/env python3
"""Build the predefined reference databases that gene-name search cannot.

Some targets are not genes and cannot be retrieved by '<gene>[Gene Name]':

  sccmec   a mobile CASSETTE carrying dozens of CDS. References are the CDS
           translations of seven whole-cassette GenBank records:
             AB033763 D86934 AB037671 AB063172 AB121219 AF411935 FJ390057
           Each protein is tagged core=1 (mec and ccr complexes) or core=0
           (transposases, hypothetical proteins).

THIS SCRIPT WRITES EVERY PROTEIN, core=0 and core=1 alike. The filtering
happens downstream: build_protein_db.py uses core=1 only by default, and
--predefined-all makes it use everything.

  is256    an INSERTION SEQUENCE. Its only verifiable product is the
           transposase, so references are IS256-family transposase proteins
           from Staphylococcus epidermidis. All are tagged core=1.

Run once per database and commit the output; build_protein_db.py then uses the
file whenever the target appears in the gene list, with no network call.

Usage:
    python3 fetch_predefined_db.py --db sccmec --out ../data/sccmec_proteins.faa
    python3 fetch_predefined_db.py --db is256  --out ../data/is256_proteins.faa
    python3 fetch_predefined_db.py --db all    --outdir ../data
"""

import argparse
import hashlib
import json
import os
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request

BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"

# ── sccmec ──────────────────────────────────────────────────────────────────
SCCMEC_ACCESSIONS = [
    "AB033763", "D86934", "AB037671", "AB063172",
    "AB121219", "AF411935", "FJ390057",
]

# Diagnostic core: the mec and ccr complexes, which is what SCCmec typing
# rests on. Everything else in a cassette (transposases, recombinase
# fragments, hypothetical proteins) also occurs elsewhere in staphylococcal
# genomes and on unrelated mobile elements.
SCCMEC_CORE_PATTERNS = [
    r"\bmecA\b", r"\bmecC\b", r"\bmecR1?\b", r"\bmecI\b",
    r"\bccrA\d?\b", r"\bccrB\d?\b", r"\bccrC\d?\b",
    r"penicillin-binding protein 2a", r"penicillin-binding protein 2'",
    r"cassette chromosome recombinase", r"methicillin resistance",
]

# ── is256 ───────────────────────────────────────────────────────────────────
# Searched by PRODUCT name: most records carry no 'IS256' gene symbol, and the
# [Protein Name] index alone is unreliable, so several title forms are OR'd.
IS256_TERM = (
    '("IS256 family transposase"[Protein Name] '
    'OR "IS256 family transposase"[Title] '
    'OR "transposase IS256"[Title] '
    'OR "IS256"[Title]) '
    'AND "Staphylococcus epidermidis"[Organism]'
)
# IS256-family transposases run ~390 aa; anything far off that is a fragment
# or a mis-annotated neighbour.
IS256_MIN_LEN, IS256_MAX_LEN = 250, 500


def fetch(url, timeout=120, retries=4):
    delay = 5
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=timeout) as r:
                return r.read().decode("utf-8", errors="replace")
        except urllib.error.HTTPError as e:
            if (e.code == 429 or e.code >= 500) and attempt < retries - 1:
                wait = delay * (2 ** attempt)
                print(f"  HTTP {e.code}, retrying in {wait}s", file=sys.stderr)
                time.sleep(wait); continue
            raise
        except (urllib.error.URLError, TimeoutError):
            if attempt < retries - 1:
                wait = delay * (2 ** attempt)
                print(f"  network error, retrying in {wait}s", file=sys.stderr)
                time.sleep(wait); continue
            raise
    raise RuntimeError("unreachable")


def _params(extra, api_key):
    p = dict(extra)
    if api_key:
        p["api_key"] = api_key
    return urllib.parse.urlencode(p)


# ── sccmec: nucleotide records -> CDS translations ──────────────────────────

def parse_cds_aa(text, source_acc, core_re, n_so_far=0):
    """Parse efetch rettype=fasta_cds_aa output.

    Headers look like:
      >lcl|AB033763.1_prot_BAB41255.1_2 [gene=mecA] [protein=penicillin-binding
       protein 2'] [protein_id=BAB41255.1] [location=...] [gbkey=CDS]
    """
    recs, hdr, chunks = [], None, []

    def flush():
        if hdr is None:
            return
        seq = "".join(chunks).replace("*", "").strip()
        if not seq:
            return
        gene = (re.search(r"\[gene=([^\]]+)\]", hdr) or [None, ""])[1]
        protein = (re.search(r"\[protein=([^\]]+)\]", hdr) or [None, ""])[1]
        pid = (re.search(r"\[protein_id=([^\]]+)\]", hdr) or [None, ""])[1]
        if not pid:
            m = re.search(r"_prot_([A-Z0-9_.]+)_\d+", hdr)
            pid = m.group(1) if m else f"{source_acc}_cds{n_so_far + len(recs) + 1}"
        label = " ".join(x for x in (gene, protein) if x) or "unnamed CDS"
        recs.append({"protein_id": pid, "gene": gene, "product": protein,
                     "label": label, "source": source_acc, "sequence": seq,
                     "core": bool(core_re.search(label))})

    for line in text.splitlines():
        if line.startswith(">"):
            flush(); hdr, chunks = line[1:], []
        elif hdr is not None:
            chunks.append(line.strip())
    flush()
    return recs


def build_sccmec(accessions, api_key, sleep):
    core_re = re.compile("|".join(SCCMEC_CORE_PATTERNS), re.I)
    out, seen = [], set()
    for acc in accessions:
        print(f"{acc} ...", file=sys.stderr)
        try:
            text = fetch(f"{BASE}efetch.fcgi?" + _params(
                {"db": "nuccore", "id": acc,
                 "rettype": "fasta_cds_aa", "retmode": "text"}, api_key))
        except Exception as e:                                   # noqa: BLE001
            print(f"  FAILED: {e}", file=sys.stderr); continue
        recs = parse_cds_aa(text, acc, core_re, len(out))
        kept = 0
        for r in recs:
            h = hashlib.sha256(r["sequence"].encode()).hexdigest()
            if h in seen:
                continue
            seen.add(h); out.append(r); kept += 1
        print(f"  {len(recs)} CDS ({sum(x['core'] for x in recs)} core), "
              f"{kept} new after dedup", file=sys.stderr)
        time.sleep(sleep)
    return out, {"accessions": accessions,
                 "core_patterns": SCCMEC_CORE_PATTERNS}


# ── is256: protein query ────────────────────────────────────────────────────

def build_is256(retmax, api_key, sleep):
    print(f"query: {IS256_TERM}", file=sys.stderr)
    data = json.loads(fetch(f"{BASE}esearch.fcgi?" + _params(
        {"db": "protein", "term": IS256_TERM,
         "retmax": retmax, "retmode": "json"}, api_key)))
    ids = data["esearchresult"].get("idlist", [])
    count = data["esearchresult"].get("count", "0")
    print(f"  matched {count}, retrieving {len(ids)}", file=sys.stderr)
    if not ids:
        return [], {"query": IS256_TERM}
    time.sleep(sleep)

    text = fetch(f"{BASE}efetch.fcgi?" + _params(
        {"db": "protein", "id": ",".join(ids),
         "rettype": "fasta", "retmode": "text"}, api_key))

    out, seen, skipped = [], set(), 0
    hdr, chunks = None, []

    def flush():
        nonlocal skipped
        if hdr is None:
            return
        seq = "".join(chunks).replace("*", "").strip()
        if not seq:
            return
        acc = hdr.split()[0]
        desc = hdr.split(None, 1)[1] if len(hdr.split(None, 1)) > 1 else ""
        org = ""
        if desc.rstrip().endswith("]"):
            org = desc.rstrip()[:-1].rsplit("[", 1)[-1]
            desc = desc.rstrip()[:desc.rstrip().rfind("[")].strip()
        # length sanity: IS256-family transposases are ~390 aa
        if not (IS256_MIN_LEN <= len(seq) <= IS256_MAX_LEN):
            skipped += 1
            return
        h = hashlib.sha256(seq.encode()).hexdigest()
        if h in seen:
            return
        seen.add(h)
        out.append({"protein_id": acc, "gene": "IS256", "product": desc,
                    "label": desc or "IS256 family transposase",
                    "source": org or "Staphylococcus epidermidis",
                    "sequence": seq, "core": True})

    for line in text.splitlines():
        if line.startswith(">"):
            flush(); hdr, chunks = line[1:], []
        elif hdr is not None:
            chunks.append(line.strip())
    flush()
    print(f"  {len(out)} unique after dedup "
          f"({skipped} dropped outside {IS256_MIN_LEN}-{IS256_MAX_LEN} aa)",
          file=sys.stderr)
    return out, {"query": IS256_TERM,
                 "length_filter": [IS256_MIN_LEN, IS256_MAX_LEN]}


# ── output ──────────────────────────────────────────────────────────────────

def write_db(recs, label, out_path, extra_meta):
    out_dir = os.path.dirname(os.path.abspath(out_path)) or "."
    os.makedirs(out_dir, exist_ok=True)

    # Header format is what call_genes.py parses: ACC|GENE|ORGANISM|PRODUCT
    with open(out_path, "w") as fh:
        for r in recs:
            fh.write(f">{r['protein_id']}|{label}|{r['source']}|"
                     f"{r['label']} [core={1 if r['core'] else 0}] "
                     f"[len={len(r['sequence'])}]\n")
            for i in range(0, len(r["sequence"]), 60):
                fh.write(r["sequence"][i:i + 60] + "\n")

    manifest = out_path.rsplit(".", 1)[0] + ".manifest.json"
    with open(manifest, "w") as fh:
        json.dump({"label": label, "n_proteins": len(recs),
                   "n_core": sum(r["core"] for r in recs),
                   **extra_meta,
                   "proteins": [{k: v for k, v in r.items() if k != "sequence"}
                                for r in recs]}, fh, indent=1)

    lens = sorted(len(r["sequence"]) for r in recs) or [0]
    n_core = sum(r["core"] for r in recs)
    n_noncore = len(recs) - n_core
    print(f"\nwrote {out_path}", file=sys.stderr)
    print(f"  {len(recs)} proteins written: {n_core} core=1, "
          f"{n_noncore} core=0; lengths {lens[0]}-{lens[-1]} aa",
          file=sys.stderr)
    print(f"wrote {manifest}", file=sys.stderr)
    if n_noncore:
        # This script writes EVERYTHING; the core filter lives downstream.
        print(f"\n  NOTE: all {len(recs)} proteins are in the file. "
              f"build_protein_db.py uses core=1 only by default\n"
              f"        ({n_core} of {len(recs)} here) -- pass "
              f"--predefined-all to score the {n_noncore} core=0 protein(s) "
              f"too.", file=sys.stderr)


def main(argv=None):
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--db", required=True, choices=["sccmec", "is256", "all"])
    p.add_argument("--out", default="", help="output FASTA (single --db)")
    p.add_argument("--outdir", default="", help="output directory (--db all)")
    p.add_argument("--accessions", default=",".join(SCCMEC_ACCESSIONS),
                   help="sccmec: comma-separated nucleotide accessions")
    p.add_argument("--retmax", type=int, default=200,
                   help="is256: NCBI records to retrieve (default 200)")
    p.add_argument("--api-key", default=os.environ.get("NCBI_API_KEY", ""))
    args = p.parse_args(argv)

    sleep = 0.15 if args.api_key else 0.4
    targets = ["sccmec", "is256"] if args.db == "all" else [args.db]

    if args.db == "all" and not args.outdir:
        print("ERROR: --db all needs --outdir", file=sys.stderr)
        return 2
    if args.db != "all" and not args.out:
        print("ERROR: --out is required", file=sys.stderr)
        return 2

    rc = 0
    for t in targets:
        out = (args.out if args.db != "all"
               else os.path.join(args.outdir, f"{t}_proteins.faa"))
        print(f"\n=== {t} ===", file=sys.stderr)
        if t == "sccmec":
            accs = [a.strip() for a in args.accessions.split(",") if a.strip()]
            recs, meta = build_sccmec(accs, args.api_key, sleep)
            label = "SCCmec"
        else:
            recs, meta = build_is256(args.retmax, args.api_key, sleep)
            label = "IS256"
        if not recs:
            print(f"ERROR: {t}: nothing retrieved", file=sys.stderr)
            rc = 1
            continue
        write_db(recs, label, out, meta)
    return rc


if __name__ == "__main__":
    sys.exit(main())

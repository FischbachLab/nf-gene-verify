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
import time
import urllib.parse
import urllib.request

BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
DEFAULT_GENES = ["agrA", "agrB", "blaI", "blaZ", "fosB", "sarA", "sdrF", "sdrH"]

# Known gene-name fusions: RefSeq annotates these under the fused name, so a
# bare '<gene>[Gene Name]' query returns nothing. Keys are lowercase.
# Tokens containing '[' are used verbatim as field-restricted queries.
ALIASES = {
    "aaca": ["aacA-aphD"],  # bifunctional AAC/APH (Tn4001); plasmid-borne,
                            # essentially absent from RefSeq gene-name records
    "fdh": ['"formate dehydrogenase"'],
        # bare 'fdh[All Fields]' wrongly matches zinc-dependent ALCOHOL
        # dehydrogenases; the [Protein Name] field index is unreliable
        # (1 hit), so use the untagged quoted phrase
    "is256": ['"IS256 family transposase"[Protein Name]'],
        # IS256 is an insertion sequence; its verifiable gene is the transposase
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


def name_query(names):
    """OR-query over tokens; tokens with '[' or '"' are verbatim queries."""
    parts = [n if ("[" in n or '"' in n) else f"{n}[Gene Name]" for n in names]
    return "(" + " OR ".join(parts) + ")"


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--genes", nargs="+", default=DEFAULT_GENES,
                    help="Gene names to retrieve (default: 8 S. epidermidis genes)")
    ap.add_argument("--out", required=True, help="Output directory")
    ap.add_argument("--retmax", type=int, default=30,
                    help="Max NCBI records per gene (default 30)")
    ap.add_argument("--api-key", default=os.environ.get("NCBI_API_KEY", ""),
                    help="NCBI API key (or set NCBI_API_KEY env var)")
    ap.add_argument("--alias", action="append", default=[],
                    metavar="GENE=ALIAS1,ALIAS2",
                    help="Extra gene-name tokens to query for GENE "
                         "(repeatable). Built-in: aacA=aacA-aphD")
    ap.add_argument("--extra-fasta", action="append", default=[],
                    metavar="GENE=PATH",
                    help="Add user-supplied reference proteins for GENE from "
                         "a FASTA (headers: >acc description). Use for genes "
                         "NCBI gene-name search cannot find. Repeatable.")
    args = ap.parse_args()

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

    for gene in args.genes:
        names = gene_names(gene, aliases)
        name_q = name_query(names)
        term = f"{name_q} AND Staphylococcus[Organism] AND srcdb_refseq[PROP]"
        ids = esearch(term, args.retmax, args.api_key)
        time.sleep(sleep)
        if not ids:
            # Resistance genes (aacA, ant, erm, ...) are plasmid/transposon
            # borne and often missing from RefSeq gene-name records; retry
            # without the RefSeq restriction before giving up.
            term = f"{name_q} AND Staphylococcus[Organism]"
            ids = esearch(term, args.retmax, args.api_key)
            time.sleep(sleep)
            if ids:
                print(f"[NOTE] {gene}: no RefSeq records under "
                      f"{'/'.join(names)}; using {len(ids)} non-RefSeq "
                      f"(GenBank) records")
            else:
                print(f"[WARN] no records for {gene} (will be reported as "
                      f"no_reference, not absent)")
                continue

        summary = esummary(ids, args.api_key)
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

        seqs = efetch_fasta([r["uid"] for r in records], args.api_key)
        time.sleep(sleep)

        # per-gene dedup by sequence hash; prefer RefSeq, drop ', partial'
        seen, seen_acc, kept = set(), set(), []
        for r in sorted(records, key=lambda x: (not x["is_refseq"], -x["length"])):
            seq = seqs.get(r["accession"])
            if not seq or r["accession"] in seen_acc:
                continue
            if ", partial" in r["product"].lower():
                continue
            h = hashlib.sha256(seq.encode()).hexdigest()
            if h in seen:
                continue
            seen.add(h)
            seen_acc.add(r["accession"])
            r["sequence"] = seq
            kept.append(r)

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
        print(f"{gene}: {len(records)} fetched, {len(kept)} unique kept")

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

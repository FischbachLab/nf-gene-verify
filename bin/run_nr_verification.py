#!/usr/bin/env python3
"""Verify candidate gene loci against NCBI nr.

Two modes:

1. Local (--db PATH): run blastx locally against a pre-built BLAST database
   (e.g., /mnt/efs/databases/Blast/nr/db/nr). Fast, no queue; requires the
   nr DB volumes (nr.Psq etc.) and optionally the taxonomy files
   (taxdb.* / nodes.dmp) for organism names.

2. Remote (default): submit to NCBI's BLAST URL API and poll. Queue-dependent
   (minutes to hours); no local DB needed.

Both modes write the raw result per locus plus a combined parsed TSV
(nr_top_hits.tsv) with per-hit identity, coverage, and description.

Usage:
  run_nr_verification.py --loci-dir DIR --out DIR [--db /path/to/nr]
                         [--expect 1e-10] [--hitlist 100]
                         [--max-target-seqs 100] [--max-wait-minutes 120]
"""
import argparse
import csv
import json
import os
import shutil
import subprocess
import time
import urllib.parse
import urllib.request

BLAST_URL = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
# Columns of nr_top_hits.tsv. Defined once: finalize_calls.py reads these by
# name, and the no-loci path below writes a header-only file that must match
# the populated one exactly.
TOP_HIT_COLUMNS = ["gene", "rank", "accession", "title", "organism", "hit_len",
                   "bit_score", "evalue", "identity_pct", "align_len",
                   "query_cov_pct", "hit_cov_pct", "query_start", "query_end"]

OUTFMT = ("6 qseqid sseqid stitle staxids sskingdoms pident length qlen slen "
          "qstart qend sstart send evalue bitscore")


# ---------------------------------------------------------------- local mode

def run_local(loci_dir, out_dir, db_path, expect, max_target_seqs, threads,
              max_hsps, seg="no"):
    """Run blastx per locus against a local DB; save raw TSVs; return genes done."""
    env = os.environ.copy()
    env["BLAST_DB_USE_MMAP"] = "0"  # avoid memory-map errors on network FS (EFS/NFS)
    if threads <= 0:  # auto-detect: use all available cores
        threads = os.cpu_count() or 1
    done = []
    for f in sorted(os.listdir(loci_dir)):
        if not f.endswith("_locus.fasta"):
            continue
        gene = f.replace("_locus.fasta", "")
        raw = os.path.join(out_dir, f"{gene}_nr_local_raw.tsv")
        if os.path.exists(raw):
            done.append(gene)
            continue
        tmp = raw + ".tmp"
        # -seg: blastx masks low-complexity tracts in the translated query by
        # default and EXCLUDES those residues from the alignment. For repeat
        # proteins (Sdr family: serine-aspartate repeats) that silently
        # truncates the alignment -- sdrH lost its N-terminal 166 aa, reading
        # as 'partial' for a full-length gene. Must match round 1's setting or
        # the two rounds report identity over different residue ranges.
        cmd = ["blastx", "-query", os.path.join(loci_dir, f), "-db", db_path,
               "-seg", seg,
               "-evalue", expect, "-outfmt", OUTFMT,
               "-max_target_seqs", str(max_target_seqs),
               "-max_hsps", str(max_hsps),
               "-num_threads", str(threads), "-out", tmp]
        print(f"[local] blastx {gene} vs {db_path} ...")
        subprocess.run(cmd, check=True, env=env)
        # prepend a column-header line so the raw TSV is self-describing
        with open(raw, "w") as out, open(tmp) as src:
            out.write("\t".join(OUTFMT.split()[1:]) + "\n")
            shutil.copyfileobj(src, out)
        os.remove(tmp)
        done.append(gene)
    return done


def parse_local_tsv(path):
    """Parse one locus's tabular blastx output into hit records.

    OUTFMT columns: qseqid sseqid stitle staxids sskingdoms pident length
    qlen slen qstart qend sstart send evalue bitscore

    Groups HSP rows by subject; hit coverage = union of subject spans /
    subject length; identity = length-weighted mean across HSPs. Also keeps
    the query (nt) span per hit so finalize_calls.py can map hits back to
    genomic coordinates.
    """
    hsps = {}
    titles = {}
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        if len(f) < 15 or f[0] == "qseqid":  # skip header line
            continue
        sseqid, stitle, staxids, kingdom = f[1], f[2], f[3], f[4]
        pident, alen = float(f[5]), int(f[6])
        slen = int(f[8])
        q1, q2 = sorted((int(f[9]), int(f[10])))
        s1, s2 = sorted((int(f[11]), int(f[12])))
        evalue, bits = f[13], float(f[14])
        hsps.setdefault(sseqid, []).append(
            {"q1": q1, "q2": q2, "s1": s1, "s2": s2,
             "pident": pident, "alen": alen, "evalue": evalue, "bits": bits})
        titles[sseqid] = (stitle, staxids, kingdom, slen)

    rows = []
    for sseqid, hs in hsps.items():
        stitle, staxids, kingdom, slen = titles[sseqid]
        merged = []
        for s, e in sorted((h["s1"], h["s2"]) for h in hs):
            if merged and s <= merged[-1][1] + 1:
                merged[-1][1] = max(merged[-1][1], e)
            else:
                merged.append([s, e])
        cov = 100.0 * sum(e - s + 1 for s, e in merged) / max(1, slen)
        ident = (sum(h["pident"] * h["alen"] for h in hs)
                 / sum(h["alen"] for h in hs))
        best = max(hs, key=lambda h: h["bits"])
        # organism: trailing [bracket] in the defline, else taxid/kingdom
        org = ""
        if stitle.rstrip().endswith("]"):
            org = stitle.rstrip()[:-1].rsplit("[", 1)[-1]
        rows.append({
            "accession": sseqid,
            "title": stitle,
            "organism": org or (f"taxid:{staxids}" if staxids else kingdom),
            "hit_len": slen,
            "bit_score": round(best["bits"], 1),
            "evalue": best["evalue"],
            "identity_pct": round(ident, 1),
            "align_len": best["alen"],
            "query_cov_pct": "",  # subject coverage reported instead; see hit_cov_pct
            "hit_cov_pct": round(cov, 1),
            # query (nt) span of this hit's HSPs, within the locus window —
            # used by finalize_calls.py to map hits back to genomic coordinates
            "query_start": min(h["q1"] for h in hs),
            "query_end": max(h["q2"] for h in hs),
        })
    rows.sort(key=lambda r: -r["bit_score"])
    return rows


# --------------------------------------------------------------- remote mode

def http_post(params):
    data = urllib.parse.urlencode(params).encode()
    req = urllib.request.Request(BLAST_URL, data=data)
    return urllib.request.urlopen(req, timeout=60).read().decode()


def http_get(url):
    return urllib.request.urlopen(url, timeout=60).read().decode()


def submit_locus(fasta_path, expect, hitlist, seg="no"):
    query = "".join(l for l in open(fasta_path) if not l.startswith(">"))
    rid = None
    for attempt in range(5):
        try:
            html = http_post({
                "CMD": "Put", "PROGRAM": "blastx", "DATABASE": "nr",
                "QUERY": query, "EXPECT": expect, "HITLIST_SIZE": hitlist,
                # URL API spells -seg as FILTER: 'L' enables low-complexity
                # masking, 'F' (or empty) disables it. Keep in step with the
                # local path, or remote and local results are not comparable.
                "FILTER": "L" if seg != "no" else "F",
            })
            for line in html.splitlines():
                if line.startswith("    RID ="):
                    rid = line.split("=")[1].strip()
                    break
            if rid:
                return rid
        except Exception as e:
            wait = 30 * (2 ** attempt)
            print(f"[WARN] submit {fasta_path}: {e}; retry in {wait}s")
            time.sleep(wait)
    raise RuntimeError(f"could not submit {fasta_path}")


def poll(rid, max_wait_minutes):
    """Return result text when ready, or None if timed out.

    When a remote BLAST result is ready the API returns the payload directly
    (JSON for FORMAT_TYPE=JSON2_S) with no Status marker; while queued it
    returns an HTML/JSON page containing Status=WAITING.
    """
    deadline = time.time() + max_wait_minutes * 60
    delay = 30
    while time.time() < deadline:
        time.sleep(delay)
        delay = min(delay * 1.3, 120)
        try:
            url = (f"{BLAST_URL}?CMD=Get&FORMAT_TYPE=JSON2_S&RID={rid}"
                   f"&ALIGNMENTS=50&DESCRIPTIONS=50")
            text = http_get(url)
        except Exception as e:
            print(f"[WARN] poll {rid}: {e}")
            continue
        stripped = text.lstrip()
        if stripped.startswith("{") and "BlastOutput2" in text:
            return text  # result payload delivered directly
        if "Status=UNKNOWN" in text:
            raise RuntimeError(f"RID {rid} expired/unknown")
        if "Status=FAILED" in text:
            raise RuntimeError(f"RID {rid} failed on NCBI side")
        # else still WAITING
    return None


def parse_json2(text):
    """Return list of parsed hits: dicts with description, identity, coverage."""
    data = json.loads(text)
    out = []
    for item in data["BlastOutput2"]:
        # JSON2_S layout: item['report']['results']['search']['hits'];
        # each hit's 'description' is a list (one per subject sequence)
        try:
            hits = item["report"]["results"]["search"]["hits"]
        except (KeyError, TypeError):
            continue
        for h in hits:
            descs = h.get("description", [])
            desc = descs[0] if descs else {}
            hsps = h.get("hsps", [])
            if not hsps:
                continue
            best = max(hsps, key=lambda x: x.get("bit_score", 0))
            # query coverage: union of hsp query spans (query is nucleotide; /3 for aa)
            spans = sorted((min(x["query_from"], x["query_to"]),
                            max(x["query_from"], x["query_to"])) for x in hsps)
            merged = []
            for s, e in spans:
                if merged and s <= merged[-1][1] + 1:
                    merged[-1][1] = max(merged[-1][1], e)
                else:
                    merged.append([s, e])
            cov_nt = sum(e - s + 1 for s, e in merged)
            out.append({
                "accession": desc.get("accession", ""),
                "title": desc.get("title", ""),
                "organism": desc.get("sciname", ""),
                "hit_len": h.get("len", 0),
                "bit_score": round(best.get("bit_score", 0), 1),
                "evalue": best.get("evalue", ""),
                "identity_pct": round(100.0 * best.get("identity", 0)
                                      / max(1, best.get("align_len", 1)), 1),
                "align_len": best.get("align_len", 0),
                "query_cov_pct": round(100.0 * cov_nt / 3 / max(1, h.get("len", 1)), 1),
                # query (nt) span of this hit's HSPs, within the locus window
                "query_start": min(min(x["query_from"], x["query_to"]) for x in hsps),
                "query_end": max(max(x["query_from"], x["query_to"]) for x in hsps),
            })
    return out


# --------------------------------------------------------------------- main

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--loci-dir", required=True, help="dir with *_locus.fasta files")
    ap.add_argument("--out", required=True, help="output dir for raw + parsed results")
    ap.add_argument("--db", default="",
                    help="local BLAST DB path (e.g., /mnt/efs/databases/Blast/nr/db/nr); "
                         "enables local mode instead of remote submission")
    ap.add_argument("--expect", default="1e-10")
    ap.add_argument("--hitlist", default="100",
                    help="remote mode: HITLIST_SIZE (default 100)")
    ap.add_argument("--max-target-seqs", default="100",
                    help="local mode: -max_target_seqs (default 100; matches "
                         "remote hitlist — finalize_calls.py only needs the "
                         "top overlapping hit)")
    ap.add_argument("--seg", default="no", choices=["no", "yes"],
                    help="blastx low-complexity masking of the translated "
                         "query (default no). 'yes' EXCLUDES masked residues "
                         "from the alignment, which truncates repeat proteins: "
                         "sdrH lost its N-terminal 166 aa and scored 100%% "
                         "identity over only the unmasked 270 residues while "
                         "the full-length figure was 89.7%%. Keep this equal "
                         "to round 1's setting.")
    ap.add_argument("--max-hsps", default="1",
                    help="local mode: -max_hsps (default 1; one HSP per "
                         "subject, collapses redundant per-subject rows)")
    ap.add_argument("--threads", type=int, default=0,
                    help="local mode: blastx threads (default 0 = auto-detect "
                         "all available cores)")
    ap.add_argument("--max-wait-minutes", type=int, default=120,
                    help="remote mode: polling timeout per locus")
    args = ap.parse_args()

    os.makedirs(args.out, exist_ok=True)

    # A MISSING directory is an error; an EXISTING directory with no loci is
    # not. The two used to share one exit(1), and that killed a pipeline run:
    # a gene absent in round 1 has no coordinates, so extract_loci.py writes an
    # empty nr_verification/ -- a valid state that means "nothing to verify
    # here", which is exactly what blaZ does in 6 of 29 genomes. Exiting
    # non-zero there made NR_COLLECT fail and took the whole run with it.
    #
    # Pointing at the wrong path still fails loudly, which is the case the
    # original check was written for.
    if not os.path.isdir(args.loci_dir):
        raise SystemExit(f"no such loci directory: {args.loci_dir}")
    fastas = sorted(f for f in os.listdir(args.loci_dir) if f.endswith("_locus.fasta"))
    if not fastas:
        print(f"no *_locus.fasta files in {args.loci_dir} -- nothing to verify")
        # Write the header-only output the caller expects, so downstream steps
        # read an empty result rather than a missing file.
        with open(os.path.join(args.out, "nr_top_hits.tsv"), "w") as fh:
            fh.write("\t".join(TOP_HIT_COLUMNS) + "\n")
        return 0

    if args.db:
        genes = run_local(args.loci_dir, args.out, args.db, args.expect,
                          args.max_target_seqs, args.threads, args.max_hsps,
                          args.seg)
        raw_name = "{gene}_nr_local_raw.tsv"
    else:
        # submit all, save RIDs immediately (resumable)
        rids_path = os.path.join(args.out, "rids.json")
        rids = json.load(open(rids_path)) if os.path.exists(rids_path) else {}
        for f in fastas:
            gene = f.replace("_locus.fasta", "")
            if gene in rids:
                continue
            rid = submit_locus(os.path.join(args.loci_dir, f), args.expect,
                               args.hitlist, args.seg)
            rids[gene] = rid
            json.dump(rids, open(rids_path, "w"), indent=1)
            print(f"submitted {gene}: RID={rid}")
            time.sleep(10)  # be polite to the queue

        # poll all pending
        pending = dict(rids)
        while pending:
            for gene, rid in list(pending.items()):
                raw_path = os.path.join(args.out, f"{gene}_nr_raw.json")
                if os.path.exists(raw_path):
                    del pending[gene]
                    continue
                print(f"polling {gene} ({rid})...")
                text = poll(rid, args.max_wait_minutes)
                if text is None:
                    print(f"[WARN] {gene}: timed out after {args.max_wait_minutes} min")
                    continue
                with open(raw_path, "w") as fh:
                    fh.write(text)
                print(f"{gene}: result saved ({len(text)} bytes)")
                del pending[gene]
        genes = [f.replace("_locus.fasta", "") for f in fastas]
        raw_name = "{gene}_nr_raw.json"

    # parse whatever results exist into a combined TSV
    rows = []
    for gene in genes:
        if args.db:
            raw_path = os.path.join(args.out, raw_name.format(gene=gene))
            if not os.path.exists(raw_path):
                continue
            for rank, h in enumerate(parse_local_tsv(raw_path), 1):
                rows.append({"gene": gene, "rank": rank, **h})
        else:
            raw_path = os.path.join(args.out, raw_name.format(gene=gene))
            if not os.path.exists(raw_path):
                continue
            for rank, h in enumerate(parse_json2(open(raw_path).read()), 1):
                rows.append({"gene": gene, "rank": rank, **h})

    tsv = os.path.join(args.out, "nr_top_hits.tsv")
    if rows:
        with open(tsv, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=TOP_HIT_COLUMNS, delimiter="\t",
                               extrasaction="ignore")
            w.writeheader()
            w.writerows(rows)
        print(f"parsed {len(rows)} hits -> {tsv}")


if __name__ == "__main__":
    main()

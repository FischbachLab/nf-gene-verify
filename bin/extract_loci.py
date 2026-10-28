#!/usr/bin/env python3
"""Extract divergent/partial gene loci (plus flanks) from the assembly.

Reads gene_calls.tsv, selects loci whose call is not 'present', and writes
one FASTA per locus for nr verification.

Use a SMALL flank (~300 bp): wide flanks (e.g., 5 kb) let a single neighbor
gene's redundant nr entries saturate the remote-BLAST hit list and mask the
target ORF entirely.

Usage:
  extract_loci.py --calls gene_calls.tsv --genome genome.fasta \
                  --flank 300 --out nr_verification/
"""
import argparse
import csv
import os


def load_contigs(path):
    contigs, name, chunks = {}, None, []
    for line in open(path):
        if line.startswith(">"):
            if name is not None:
                contigs[name] = "".join(chunks)
            name = line[1:].split()[0]
            chunks = []
        else:
            chunks.append(line.strip())
    if name is not None:
        contigs[name] = "".join(chunks)
    return contigs


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--calls", required=True, help="gene_calls.tsv from call_genes.py")
    ap.add_argument("--genome", required=True, help="assembly FASTA")
    ap.add_argument("--flank", type=int, default=5000)
    ap.add_argument("--out", required=True)
    ap.add_argument("--all", action="store_true",
                    help="extract loci for ALL genes, including 'present' calls "
                         "(round-2 verification of every gene)")
    args = ap.parse_args()

    contigs = load_contigs(args.genome)
    os.makedirs(args.out, exist_ok=True)

    n = 0
    for r in csv.DictReader(open(args.calls), delimiter="\t"):
        if not args.all and (r["call"] == "present" or not r["coords"]):
            continue
        if not r["coords"]:
            continue
        gene = r["gene"]
        contig = r["contig"]
        s, e = (int(x) for x in r["coords"].split("-"))
        seq = contigs[contig]
        s2, e2 = max(1, s - args.flank), min(len(seq), e + args.flank)
        sub = seq[s2 - 1:e2]
        out = os.path.join(args.out, f"{gene}_locus.fasta")
        with open(out, "w") as fh:
            fh.write(f">{gene}_locus_{contig}_{s2}-{e2}\n")
            for i in range(0, len(sub), 70):
                fh.write(sub[i:i + 70] + "\n")
        print(f"{gene}: {e2 - s2 + 1} bp ({contig}:{s2}-{e2}) -> {out}")
        n += 1
    print(f"{n} loci extracted")


if __name__ == "__main__":
    main()

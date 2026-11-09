#!/usr/bin/env python3
"""Apply supplied group overrides, preserving every unlisted marker style."""
import argparse
import csv
from pathlib import Path
from fig1lib import read_tsv,save_json,sha256


def update(base,memo,output):
    fields,rows=read_tsv(base,['gene','group','color','shape'])
    _,overrides=read_tsv(memo,['gene','group'])
    lookup={row['gene'].casefold():row for row in rows}
    if len(lookup)!=len(rows):raise ValueError('Base gene labels are ambiguous or duplicated')
    groups={row['group'].casefold():(row['group'],row['color'],row['shape']) for row in rows}
    groups.update(biofilm=('Biofilm','#00897b','circle'),colonization=('Colonization','#8c564b','circle'))
    changes=[];seen=set()
    for supplied in overrides:
        key=supplied['gene'].casefold()
        if key not in lookup or key in seen:raise ValueError(f'Unknown or duplicated memo gene: {supplied["gene"]}')
        group=supplied['group'].casefold()
        if group not in groups:raise ValueError(f'No plotting style defined for group: {supplied["group"]}')
        seen.add(key);row=lookup[key];before=row.copy()
        row.update(zip(['group','color','shape'],groups[group]))
        changes.append(dict(memo_gene=supplied['gene'],gene=row['gene'],before=before,after=row.copy()))
    output=Path(output);output.parent.mkdir(parents=True,exist_ok=True)
    with output.open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=fields,delimiter='\t',lineterminator='\n')
        writer.writeheader();writer.writerows(rows)
    save_json(str(output)+'.memo_update.json',dict(base=str(Path(base).resolve()),memo=str(Path(memo).resolve()),
        memo_sha256=sha256(memo),listed_genes=len(changes),unlisted_genes_preserved=len(rows)-len(changes),
        matching='Whitespace trimmed; gene and group names matched case-insensitively; original gene labels retained',changes=changes))
    print(f'Applied {len(changes)} memo entries; preserved {len(rows)-len(changes)} unlisted marker styles.')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--base',type=Path,required=True)
    parser.add_argument('--memo',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    update(args.base,args.memo,args.output)

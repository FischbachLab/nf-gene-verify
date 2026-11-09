#!/usr/bin/env python3
"""Test table/tree joining and 30-tip rendering; no phylogeny is inferred here."""
import argparse
import csv
import json
import subprocess
import sys
from pathlib import Path

from fig1lib import read_tsv, sample_id


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,default=Path('smoke_results'))
    args = parser.parse_args()
    package = Path(__file__).resolve().parent
    args.out.mkdir(parents=True,exist_ok=True)
    for case in ['two_genomes','synthetic_30']:
        source = package/'examples'/case
        prefix = args.out/case
        subprocess.run([sys.executable,str(package/'plot_fig1.py'),'--tree',str(source/'tree.nwk'),
                        '--genes',str(source/'genes.tsv'),'--metadata',str(source/'metadata.tsv'),
                        '--gene-groups',str(package/'gene_groups.tsv'),'--output',str(prefix),
                        '--title',('Two-genome test: NIHLM087 and NIHLM095' if case == 'two_genomes' else 'SYNTHETIC 30-sample layout test (not biological results)'),
                        '--branch-units',('mismatches / aligned base' if case == 'two_genomes' else 'synthetic units')],check=True)
        result = json.loads(Path(str(prefix)+'.plot.json').read_text())
        _, raw = read_tsv(source/'genes.tsv',['sample','gene','present'])
        expected = {(sample_id(r['sample']),r['gene']):r['present'] for r in raw}
        _, plotted = read_tsv(Path(str(prefix)+'.matrix.tsv'),['sample'])
        for row in plotted:
            for gene in result['gene_order']:
                if row[gene] != expected[row['sample'],gene]:
                    raise AssertionError(f'Tree/table joining failed: {case}, {row["sample"]}, {gene}')
        wanted = 2 if case == 'two_genomes' else 30
        if len(result['sample_order']) != wanted or len(plotted) != wanted:
            raise AssertionError('Wrong number of plotted samples')
        for suffix in ['pdf','svg','png']:
            if not Path(str(prefix)+'.'+suffix).is_file():
                raise AssertionError(f'Missing {suffix} output')
        print(f'PASS: {case}; every marker call matches its sample in tree-tip order.')
    print('Smoke test passed. This checks plotting; it does not exercise genome annotation or tree inference.')


if __name__ == '__main__':
    main()

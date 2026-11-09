#!/usr/bin/env python3
"""Batch bacterial core-genome tree and Fig. 1-style marker plot for Linux EC2.

Inputs: a folder of assembled nucleotide genomes, a long-form gene TSV, and
optional sample metadata. Run --check first. See README.md for installation.
"""
import argparse
import csv
import hashlib
import json
import os
import platform
import shlex
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from Bio import SeqIO, Phylo

from fig1lib import (assembly_stats,discover_genomes,read_genes,read_metadata,
                     read_fasta,read_tree,sample_id,save_json,sha256,write_metadata)

VERSION = '1.0.0'
OWNER = 'fig1_ec2_pipeline'


def fingerprint(data):
    return hashlib.sha256(json.dumps(data,sort_keys=True).encode()).hexdigest()


def tool_versions(infer_st=False, supplied_tree=False):
    versions, executables = {}, {}
    specs = [] if supplied_tree else [('prokka',['--version']),('panaroo',['--version']),
                                      ('mafft',['--version']),('iqtree',['--version'])]
    if infer_st:
        specs.append(('mlst',['--version']))
    for name,flags in specs:
        candidates = ['iqtree2','iqtree'] if name == 'iqtree' else [name]
        executable = next((shutil.which(n) for n in candidates if shutil.which(n)),None)
        if executable is None:
            raise ValueError(f'{name} was not found. Activate the fig1 environment; see README.md.')
        completed = subprocess.run([executable]+flags,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,
                                   text=True,timeout=30,check=False)
        if completed.returncode:
            raise ValueError(f'Cannot run {name}: {completed.stdout[-1500:]}')
        versions[name] = completed.stdout.strip()
        executables[name] = executable
    versions['python'] = sys.version
    versions['platform'] = platform.platform()
    return executables,versions


def alignment_summary(path, samples):
    with Path(path).open() as stream:
        records = list(SeqIO.parse(stream,'fasta'))
    labels = [sample_id(r.id) for r in records]
    if len(set(labels)) != len(labels) or set(labels) != set(samples):
        raise ValueError(f'Core alignment sample names do not match the input genomes: {path}')
    lengths = {len(r.seq) for r in records}
    if len(lengths) != 1 or next(iter(lengths),0) == 0:
        raise ValueError(f'Core alignment is empty or has unequal sequence lengths: {path}')
    sequences = [str(r.seq).upper() for r in records]
    allowed = set('ACGTRYSWKMBDHVN-?')
    if any(set(s)-allowed for s in sequences):
        raise ValueError('Core alignment has an unexpected nucleotide character')
    length = len(sequences[0])
    variable = 0
    for start in range(0,length,10000):
        block = np.vstack([np.frombuffer(s[start:start+10000].encode(),dtype='S1') for s in sequences])
        alleles = sum(np.any(block == b,axis=0).astype(np.int8) for b in [b'A',b'C',b'G',b'T'])
        variable += int((alleles>1).sum())
    if variable == 0:
        raise ValueError('No variable canonical sites in the chosen core alignment. A supported phylogeny cannot be inferred from this alignment.')
    return dict(samples=len(samples),columns=length,variable_canonical_columns=variable,
                identical_alignment_groups=len(sequences)-len(set(sequences)))


class Stages:
    """Resume only when inputs, parameters, versions and output hashes agree."""
    def __init__(self,output,resume):
        self.output,self.resume = output,resume

    def run(self,label,signature,command_factory,expected,validator=None):
        directory = self.output/'work'/label
        state_path = directory/'stage.json'
        state = json.loads(state_path.read_text()) if state_path.exists() else None
        if directory.exists():
            if not self.resume:
                raise ValueError(f'Stage already exists: {directory}. Use --resume or a fresh output directory.')
            if not state or state.get('owner') != OWNER:
                raise ValueError(f'Unrecognized stage directory; refusing to overwrite: {directory}')
            if state.get('signature') != signature:
                raise ValueError(f'Inputs, parameters or software changed for {label}. Choose a fresh output directory.')
            if state.get('status') == 'complete':
                valid = all((directory/p).is_file() and sha256(directory/p) == state['outputs'].get(p) for p in expected)
                if not valid:
                    raise ValueError(f'Saved outputs for {label} are missing or changed. Choose a fresh output directory.')
                if validator:
                    validator(directory)
                print(f'[resume] {label}',flush=True)
                return directory
            # A failed attempt is preserved instead of mixing partial outputs into a new run.
            index = 1
            while directory.with_name(directory.name+f'.failed-{index}').exists():
                index += 1
            directory.rename(directory.with_name(directory.name+f'.failed-{index}'))
        directory.mkdir(parents=True,exist_ok=False)
        command = [str(v) for v in command_factory(directory)]
        state = dict(owner=OWNER,status='running',signature=signature,command=command,
                     started_utc=datetime.now(timezone.utc).isoformat())
        save_json(state_path,state)
        print(f'[run] {label}\n      {shlex.join(command)}',flush=True)
        try:
            with (directory/'tool.log').open('w') as log:
                log.write(shlex.join(command)+'\n\n')
                log.flush()
                subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,check=True)
            for relative in expected:
                file = directory/relative
                if not file.is_file() or file.stat().st_size == 0:
                    raise ValueError(f'{label} did not produce its expected nonempty output: {relative}')
            if validator:
                validator(directory)
            state.update(status='complete',outputs={p:sha256(directory/p) for p in expected},
                         completed_utc=datetime.now(timezone.utc).isoformat())
            save_json(state_path,state)
        except (Exception,KeyboardInterrupt) as error:
            state.update(status='failed',error=str(error))
            save_json(state_path,state)
            if isinstance(error,KeyboardInterrupt):
                raise
            raise ValueError(f'{label} failed: {error}. See {directory}/tool.log; rerun with --resume after fixing the cause.') from error
        return directory


def prepared_genomes(genomes,directory):
    directory.mkdir(parents=True,exist_ok=True)
    prepared = {}
    for sample,path in genomes.items():
        records = read_fasta(path)
        for index,record in enumerate(records,1):
            record.id = f'contig{index:06d}'
            record.name = record.id
            record.description = ''
        destination = directory/(sample+'.fna')
        SeqIO.write(records,destination,'fasta')
        prepared[sample] = destination
    return prepared


def merge_mlst(path,metadata,samples):
    seen, evidence = set(), []
    for line in Path(path).read_text().splitlines():
        fields = line.split('\t')
        if len(fields)<3:
            raise ValueError('Malformed mlst output')
        sample = sample_id(fields[0])
        if sample not in samples or sample in seen:
            raise ValueError(f'Unexpected or duplicated sample in mlst output: {sample}')
        seen.add(sample)
        inferred = fields[2] if fields[2].isdigit() else ''
        supplied = metadata[sample].get('ST','')
        if supplied and inferred and supplied != inferred:
            raise ValueError(f'Supplied ST{supplied} conflicts with inferred ST{inferred} for {sample}. Resolve the metadata before plotting.')
        if not supplied:
            metadata[sample]['ST'] = inferred
        evidence.append(dict(sample=sample,scheme=fields[1],inferred_ST=inferred or None,
                             supplied_ST=supplied or None,alleles=fields[3:]))
    if seen != set(samples):
        raise ValueError('mlst did not return a call for every sample')
    return evidence


def plan(args,genomes):
    if args.tree:
        print(f'Use supplied Newick tree: {args.tree}')
    else:
        print(f'Annotate {len(genomes)} genomes with Prokka, sequentially, using {args.threads} threads per genome.')
        print(f'Panaroo core alignment: clean-mode {args.clean_mode}; core threshold {args.core_threshold}; MAFFT; {args.threads} threads.')
        print(f'IQ-TREE: {args.alignment} core alignment; MFP; 1000 UFBoot + 1000 SH-aLRT; seed {args.seed}; {args.threads} threads.')
    if args.infer_st:
        print(f'Infer STs with mlst --scheme {args.mlst_scheme}; compare with supplied STs.')
    print(f'Plot tree + {args.genes} + metadata; output prefix {args.out}/figure1.')
    if args.ortholog_pca:
        print(f'PCA of genome-wide Panaroo cluster presence/absence; output prefix {args.out}/ortholog_pca/figure6_pca.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--genomes-dir',type=Path)
    parser.add_argument('--genes',type=Path)
    parser.add_argument('--metadata',type=Path)
    parser.add_argument('--gene-groups',type=Path)
    parser.add_argument('--out',type=Path,default=Path('results'))
    parser.add_argument('--threads',type=int,default=8)
    parser.add_argument('--core-threshold',type=float,default=.95)
    parser.add_argument('--clean-mode',choices=['strict','moderate','sensitive'],default='strict')
    parser.add_argument('--alignment',choices=['filtered','unfiltered'],default='filtered')
    parser.add_argument('--seed',type=int,default=42)
    parser.add_argument('--infer-st',action='store_true')
    parser.add_argument('--mlst-scheme',default='sepidermidis')
    parser.add_argument('--tree',type=Path,help='Use an existing tree instead of Prokka/Panaroo/IQ-TREE')
    parser.add_argument('--root',choices=['as-is','midpoint'],default='as-is')
    parser.add_argument('--outgroup',help='An existing sample identifier; adds no genomes')
    parser.add_argument('--title',default='Genome phylogeny and marker profiles')
    parser.add_argument('--tree-layout',choices=['phylogram','cladogram'],default='phylogram')
    parser.add_argument('--tree-width',type=float,default=4.4)
    parser.add_argument('--row-height',type=float,default=.28)
    parser.add_argument('--keep-tip-order',action='store_true')
    parser.add_argument('--support-style',choices=['symbols','values','ufboot'],default='symbols')
    parser.add_argument('--ortholog-pca',action='store_true',help='Also analyze genome-wide Panaroo cluster presence/absence')
    parser.add_argument('--pca-metadata',type=Path,help='Optional PCA metadata with source/group columns')
    parser.add_argument('--pca-rtab',type=Path,help='Existing Panaroo matrix; required for PCA when --tree skips Panaroo')
    parser.add_argument('--pca-source-column',default='source')
    parser.add_argument('--pca-group-column')
    parser.add_argument('--pca-min-presence',type=int,default=1)
    parser.add_argument('--pca-scale',action='store_true')
    parser.add_argument('--pca-no-labels',action='store_true')
    parser.add_argument('--metadata-columns',help='Comma-separated columns to display')
    parser.add_argument('--sample-labels',type=Path,help='Display alias/color/bold table for Figure 1')
    parser.add_argument('--highlight-nosocomial',action='store_true')
    parser.add_argument('--show-branch-lengths',action='store_true')
    parser.add_argument('--highlight-groups',help='Comma-separated gene-group names to shade pale yellow')
    parser.add_argument('--rotate-node',action='append',help='Reverse child order at the MRCA of comma-separated input sample names; repeat for multiple nodes')
    parser.add_argument('--compact',action='store_true',help='Larger text and symbols with tighter gene/panel spacing')
    parser.add_argument('--allow-missing',action='store_true')
    parser.add_argument('--resume',action='store_true')
    parser.add_argument('--check',action='store_true',help='Validate inputs without running or writing analysis outputs')
    parser.add_argument('--dry-run',action='store_true',help='Validate inputs and show the plan without running or writing outputs')
    parser.add_argument('--doctor',action='store_true',help='Check installed external commands and versions, then exit')
    args = parser.parse_args()
    try:
        if args.threads<1 or not 0<args.core_threshold<=1 or args.seed<1:
            raise ValueError('threads and seed must be positive; core-threshold must be in (0,1]')
        if not np.isfinite(args.tree_width) or args.tree_width < 1 or not np.isfinite(args.row_height) or args.row_height < .15:
            raise ValueError('tree-width must be >=1 inch and row-height >=0.15 inch')
        if args.outgroup and args.root != 'as-is':
            raise ValueError('Use either --outgroup or --root midpoint')
        if args.doctor:
            _,versions = tool_versions(True,args.tree is not None)
            print(json.dumps(versions,indent=2))
            return
        if args.genomes_dir is None or args.genes is None:
            raise ValueError('--genomes-dir and --genes are required (except with --doctor)')
        genomes = discover_genomes(args.genomes_dir)
        samples = list(genomes)
        stats = [assembly_stats(s,p) for s,p in genomes.items()]
        calls,genes,missing = read_genes(args.genes,samples,args.allow_missing)
        metadata,columns = read_metadata(args.metadata,samples)
        if args.ortholog_pca:
            from ortholog_pca import checked_metadata,read_rtab
            if len(samples)<3 or not 1<=args.pca_min_presence<len(samples):
                raise ValueError('PCA needs at least three genomes and pca-min-presence between 1 and n-1')
            pca_columns = checked_metadata(args.pca_metadata,samples)[1] if args.pca_metadata else columns
            if args.pca_group_column and args.pca_group_column not in pca_columns:
                raise ValueError('Requested PCA group column is missing from the PCA metadata')
            if args.tree and args.pca_rtab is None:
                raise ValueError('With --tree, PCA requires an existing --pca-rtab because Panaroo is skipped')
            if args.pca_rtab and set(read_rtab(args.pca_rtab)[0])!=set(samples):
                raise ValueError('PCA matrix sample IDs do not match the input genomes')
        from plot_fig1 import gene_styles
        gene_styles(genes,args.gene_groups)
        if args.metadata_columns and set(args.metadata_columns.split(','))-set(columns):
            raise ValueError('Requested metadata-columns do not exist in the metadata table')
        if args.outgroup and sample_id(args.outgroup) not in samples:
            raise ValueError('Outgroup is not one of the input genomes')
        if args.tree:
            read_tree(args.tree,samples)
        elif len(samples)<4 and not args.check and not args.dry_run:
            raise ValueError('Tree inference with support requires at least four genomes. For the two-genome smoke test, supply --tree examples/two_genomes/tree.nwk.')
        print(f'Validated {len(samples)} genomes, {len(genes)} markers, {len(calls)} supplied calls; {missing} absent table rows.',flush=True)
        for row in stats:
            print(f'  {row["sample"]}: {row["bases"]:,} bp, {row["contigs"]} contigs, N50 {row["N50"]:,}',flush=True)
        if args.check or args.dry_run:
            if args.dry_run:
                plan(args,genomes)
            return
        executables,versions = tool_versions(args.infer_st,args.tree is not None)
        output = args.out.resolve()
        state_path = output/'run_state.json'
        base = dict(version=VERSION,genomes={r['sample']:r['sha256'] for r in stats},
                    tree_sha256=sha256(args.tree) if args.tree else None,
                    core_threshold=args.core_threshold,clean_mode=args.clean_mode,
                    alignment=args.alignment,seed=args.seed,infer_st=args.infer_st,
                    mlst_scheme=args.mlst_scheme,tool_versions=versions)
        base_signature = fingerprint(base)
        if output.exists() and any(output.iterdir()):
            if not args.resume or not state_path.is_file():
                raise ValueError('Output folder is nonempty. Use a fresh --out directory, or --resume for a previous pipeline run.')
            previous = json.loads(state_path.read_text())
            if previous.get('owner') != OWNER or previous.get('signature') != base_signature:
                raise ValueError('Genome inputs, tree settings, or tool versions changed. Use a fresh --out directory.')
        output.mkdir(parents=True,exist_ok=True)
        save_json(state_path,dict(owner=OWNER,signature=base_signature,configuration=base,status='running'))
        save_json(output/'input_manifest.json',dict(assemblies=stats,genes_sha256=sha256(args.genes),
                                                  metadata_sha256=sha256(args.metadata) if args.metadata else None,
                                                  gene_groups_sha256=sha256(args.gene_groups) if args.gene_groups else None))
        with (output/'assembly_stats.tsv').open('w',newline='') as stream:
            writer = csv.DictWriter(stream,fieldnames=list(stats[0]),delimiter='\t')
            writer.writeheader(); writer.writerows(stats)
        save_json(output/'software_versions.json',versions)
        prepared = prepared_genomes(genomes,output/'work'/'prepared_genomes')
        stages = Stages(output,args.resume)
        mlst_evidence = None
        if args.infer_st:
            signature = fingerprint(dict(genomes=base['genomes'],scheme=args.mlst_scheme,version=versions['mlst']))
            directory = stages.run('mlst',signature,
                lambda d:[executables['mlst'],'--scheme',args.mlst_scheme,'--outfile',d/'calls.tsv']+list(prepared.values()),
                ['calls.tsv'])
            mlst_evidence = merge_mlst(directory/'calls.tsv',metadata,samples)
            save_json(output/'mlst_evidence.json',mlst_evidence)
        write_metadata(output/'metadata_used.tsv',metadata,columns)
        alignment_info = None
        if args.tree:
            tree = read_tree(args.tree,samples)
            branch_units = 'as supplied'
            method = 'User-supplied Newick tree; see the source analysis for phylogenetic methodology'
        else:
            gffs = []
            for sample,genome in prepared.items():
                prefix = 'G'+hashlib.sha256(sample.encode()).hexdigest()[:10].upper()
                signature = fingerprint(dict(genome=base['genomes'][sample],version=versions['prokka'],
                                             sample=sample,minimum_contig_length=200,locustag=prefix))
                directory = stages.run('prokka/'+sample,signature,
                    lambda d,s=sample,g=genome,l=prefix:[executables['prokka'],'--outdir',d/'annotation',
                        '--prefix',s,'--locustag',l,'--kingdom','Bacteria','--gcode','11',
                        '--mincontiglen','200','--cpus',str(args.threads),g],
                    ['annotation/'+sample+'.gff'])
                gffs.append(directory/'annotation'/(sample+'.gff'))
            alignment_name = 'core_gene_alignment_filtered.aln' if args.alignment == 'filtered' else 'core_gene_alignment.aln'
            signature = fingerprint(dict(gffs={p.name:sha256(p) for p in gffs},version=versions['panaroo'],
                                         mafft=versions['mafft'],clean_mode=args.clean_mode,core=args.core_threshold))
            directory = stages.run('panaroo',signature,
                lambda d:[executables['panaroo'],'-i']+gffs+['-o',d/'pangenome','--clean-mode',args.clean_mode,
                    '-a','core','--aligner','mafft','--core_threshold',str(args.core_threshold),'-t',str(args.threads)],
                ['pangenome/'+alignment_name],
                lambda d:alignment_summary(d/'pangenome'/alignment_name,samples))
            alignment = directory/'pangenome'/alignment_name
            alignment_info = alignment_summary(alignment,samples)
            save_json(output/'alignment_summary.json',alignment_info)
            signature = fingerprint(dict(alignment_sha256=sha256(alignment),version=versions['iqtree'],
                                         seed=args.seed,model='MFP',UFBoot=1000,SH_aLRT=1000,keep_ident=True))
            directory = stages.run('iqtree',signature,
                lambda d:[executables['iqtree'],'-s',alignment,'-m','MFP','-B','1000','-alrt','1000',
                          '-keep-ident','-seed',str(args.seed),'-T',str(args.threads),'--prefix',d/'core'],
                ['core.treefile','core.iqtree'],lambda d:read_tree(d/'core.treefile',samples))
            tree = read_tree(directory/'core.treefile',samples)
            branch_units = 'substitutions/site'
            method = 'Prokka annotation; Panaroo core-gene alignment with MAFFT; IQ-TREE MFP, UFBoot 1000, SH-aLRT 1000; no recombination-masking stage'
        Phylo.write(tree,output/'tree.nwk','newick',format_branch_length='%1.10g')
        from plot_fig1 import render
        description = render(output/'tree.nwk',args.genes,output/'figure1',output/'metadata_used.tsv',
                             args.gene_groups,args.allow_missing,args.title,
                             args.metadata_columns.split(',') if args.metadata_columns else None,
                             args.root,args.outgroup,branch_units,show_support=True,
                             tree_layout=args.tree_layout,tree_width=args.tree_width,
                             row_height=args.row_height,ladderize=not args.keep_tip_order,
                             support_style=args.support_style,sample_labels_path=args.sample_labels,
                             highlight_nosocomial=args.highlight_nosocomial,show_branch_lengths=args.show_branch_lengths,
                             highlight_groups=args.highlight_groups.split(',') if args.highlight_groups else None,
                             rotate_nodes=args.rotate_node,compact=args.compact)
        pca_description = None
        if args.ortholog_pca:
            from ortholog_pca import analyze
            pca_matrix = args.pca_rtab or output/'work'/'panaroo'/'pangenome'/'gene_presence_absence.Rtab'
            pca_description = analyze(pca_matrix,output/'ortholog_pca'/'figure6_pca',
                args.pca_metadata or output/'metadata_used.tsv',args.pca_source_column,
                args.pca_group_column,args.pca_min_presence,args.pca_scale,
                not args.pca_no_labels,expected_samples=samples)
        save_json(output/'run_summary.json',dict(pipeline_version=VERSION,samples=len(samples),markers=len(genes),
                                                method=method,alignment=alignment_info,
                                                ST_source='Provided metadata; optional mlst fills missing STs' if args.infer_st else 'Provided metadata only',
                                                plot=description,ortholog_pca=pca_description))
        save_json(state_path,dict(owner=OWNER,signature=base_signature,configuration=base,status='complete'))
        print(f'Completed: {output}/figure1.pdf',flush=True)
    except (ValueError,OSError,subprocess.SubprocessError) as error:
        parser.exit(2,f'Error: {error}\n')


if __name__ == '__main__':
    main()

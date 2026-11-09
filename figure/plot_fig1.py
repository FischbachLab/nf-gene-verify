#!/usr/bin/env python3
"""Draw a sample-aligned tree, metadata and reported gene-presence panel."""
import argparse
import copy
import json
import math
import os
from pathlib import Path

# Set before importing pyplot; works without a display or writable home directory.
os.environ.setdefault('MPLCONFIGDIR', str(Path(os.environ.get('TMPDIR', '/tmp')) / 'fig1-matplotlib'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.ticker import MaxNLocator
from matplotlib.textpath import TextPath
from matplotlib.font_manager import FontProperties
from matplotlib.text import Text
from matplotlib.transforms import Bbox

from Bio import Phylo
from fig1lib import read_tree, read_genes, read_metadata, read_tsv, save_json, write_metadata

DEFAULTS = {
    'blaZ': ('Resistance', '#df3d3d', 'square'),
    'mecA': ('Resistance', '#df3d3d', 'square'),
    **{gene: ('Virulence-associated', '#172b3a', 'circle') for gene in ['agrA','bhp','ebh','icaD','icaR']},
    'IS256': ('IS256', '#c4a000', 'triangle'),
}
MARKERS = {'square': 's', 'circle': 'o', 'triangle': '>'}


def text_width_inches(text,size):
    if not text:
        return 0.
    return TextPath((0,0),str(text),size=size,prop=FontProperties(family='DejaVu Sans')).get_extents().width/72


def gene_styles(genes, config=None):
    styles = {g: DEFAULTS.get(g, ('Other markers', '#64748b', 'square')) for g in genes}
    if config:
        _, rows = read_tsv(config, ['gene','group','color','shape'])
        seen = set()
        for row in rows:
            gene = row['gene']
            if gene in seen:
                raise ValueError(f'Duplicated gene in group configuration: {gene}')
            seen.add(gene)
            import re
            if not row['group'] or not re.fullmatch(r'#[0-9A-Fa-f]{6}', row['color']) or row['shape'] not in MARKERS:
                raise ValueError(f'Invalid group/color/shape for {gene}; colors need #RRGGBB and shapes square/circle/triangle')
            if gene in genes:
                styles[gene] = (row['group'], row['color'], row['shape'])
    groups = {}
    for gene in genes:
        label, color, shape = styles[gene]
        if label in groups and (groups[label]['color'],groups[label]['shape']) != (color,shape):
            raise ValueError(f'All genes in group {label!r} must use the same color and shape')
        groups.setdefault(label, dict(color=color,shape=shape,genes=[]))['genes'].append(gene)
    # Keep Fig. 1 drug classes together even when input rows start with IS256
    # or a virulence marker. Retain the earlier generic group names as well.
    preferred = ['Beta-lactams','Quinolones','Aminoglycosides','MLSB',
                 'Tetracycline','Other ARGs','Resistance','Biofilm','Colonization','Virulence',
                 'Virulence-associated','Other markers','Additional markers']
    ordered = [name for name in preferred if name in groups]
    ordered += [name for name in groups if name not in ordered and name != 'IS256']
    if 'IS256' in groups:
        ordered.append('IS256')
    return {name:groups[name] for name in ordered}


def sample_display_styles(samples,path=None):
    styles={sample:dict(label=sample,color='#172b3a',bold=False) for sample in samples}
    if path:
        _,rows=read_tsv(path,['sample','label','color','bold'])
        seen=set()
        for row in rows:
            sample=row['sample']
            if sample not in styles or sample in seen:
                raise ValueError(f'Unknown or duplicated styled sample: {sample}')
            if not row['label'] or not matplotlib.colors.is_color_like(row['color']) or row['bold'] not in ('0','1'):
                raise ValueError(f'Invalid label/color/bold for {sample}')
            seen.add(sample)
            styles[sample]=dict(label=row['label'],color=row['color'],bold=row['bold']=='1')
    labels=[styles[sample]['label'] for sample in samples]
    if len(set(labels))!=len(labels):
        raise ValueError('Display sample labels must be unique')
    return styles


def branch_distance_labels(fig,ax,records,depths,ys,fontsize=6.2):
    """Label edge midpoints, using short leaders when small edges are crowded."""
    fig.canvas.draw()
    renderer=fig.canvas.get_renderer()
    allowed=ax.get_window_extent(renderer).padded(-3)
    points=ax.transData.transform([(depths[c],ys[c]) for c in depths])
    occupied=[Bbox.from_bounds(x-2,y-2,4,4) for x,y in points]
    for record in sorted(records,key=lambda r:(r['y'],-r['length'])):
        anchor=((record['start']+record['end'])/2,record['y'])
        text=ax.annotate(f'{record["length"]:.3g}',anchor,xytext=(0,3),
            textcoords='offset points',ha='center',va='bottom',fontsize=fontsize,
            color='#334155',zorder=7,
            bbox=dict(facecolor='white',edgecolor='none',alpha=.9,pad=.4),
            arrowprops=dict(arrowstyle='-',color='#94a3b8',lw=.4))
        candidates=[(0,3),(0,-10)]
        candidates += [(dx,dy) for radius in [8,14,22,32,45,60,80,105]
                       for dx,dy in [(radius,3),(-radius,3),(radius,-10),(-radius,-10),
                                     (0,radius),(0,-radius-6),(radius,radius),(-radius,-radius-6)]]
        for dx,dy in candidates:
            text.set_position((dx,dy));text.update_positions(renderer)
            box=Text.get_window_extent(text,renderer).padded(1.2)
            if not (allowed.contains(box.x0,box.y0) and allowed.contains(box.x1,box.y1)):
                continue
            if any(box.overlaps(other) for other in occupied):continue
            occupied.append(box)
            record['label_offset_points']=[dx,dy]
            break
        else:
            raise ValueError('Branch-distance labels are crowded; increase --tree-width or --row-height')


def support_label(clade):
    if clade.confidence is not None:
        return f'{clade.confidence:g}'
    if clade.name:
        # IQ-TREE writes paired SH-aLRT/UFBoot supports as an internal-node label.
        try:
            values = [float(v) for v in clade.name.split('/')]
            return clade.name if all(math.isfinite(v) for v in values) else ''
        except ValueError:
            pass
    return ''


def support_values(clade):
    label = support_label(clade)
    if not label:
        return []
    return [float(value) for value in label.split('/')]


def tree_display_depths(tree, layout):
    """Return display coordinates without modifying branch lengths."""
    if layout not in ('phylogram','cladogram'):
        raise ValueError('tree-layout must be phylogram or cladogram')
    incomplete = any(c.branch_length is None for c in tree.find_clades() if c is not tree.root)
    if layout == 'phylogram' and not incomplete:
        depths = tree.depths()
        origin = depths[tree.root]
        depths = {clade:depth-origin for clade,depth in depths.items()}
        if max(depths.values()) > 0:
            return depths, True
    heights = {}
    def height(clade):
        heights[clade] = 0 if clade.is_terminal() else 1+max(height(c) for c in clade.clades)
        return heights[clade]
    maximum = height(tree.root)
    return {clade:maximum-value for clade,value in heights.items()}, False


def apply_node_rotations(tree,requests=None):
    """Reverse child order at selected MRCAs without changing any tree edge."""
    samples={tip.name for tip in tree.get_terminals()}
    targets=[]
    for request in requests or []:
        anchors=[name.strip() for name in request.split(',')]
        if len(anchors)<2 or len(set(anchors))!=len(anchors) or not set(anchors)<=samples:
            raise ValueError(f'Rotation needs at least two distinct input sample names in the tree: {request}')
        clade=tree.common_ancestor(*anchors)
        if any(clade is target[0] for target in targets):
            raise ValueError(f'The same node was requested for rotation more than once: {request}')
        targets.append((clade,anchors))
    records=[]
    for clade,anchors in targets:
        before=[[tip.name for tip in child.get_terminals()] for child in clade.clades]
        clade.clades.reverse()
        records.append(dict(anchors=anchors,descendant_samples=sorted(t.name for t in clade.get_terminals()),
                            child_order_before=before,child_order_after=list(reversed(before))))
    return records


def render(tree_path, genes_path, output, metadata_path=None, groups_path=None,
           allow_missing=False, title='Genome phylogeny and marker profiles',
           metadata_columns=None, root='as-is', outgroup=None,
           branch_units='as supplied', show_support=False,
           tree_layout='phylogram', tree_width=4.4, row_height=.28,
           ladderize=True, support_style='symbols',sample_labels_path=None,
           highlight_nosocomial=False,show_branch_lengths=False,highlight_groups=None,
           rotate_nodes=None,compact=False):
    if not math.isfinite(tree_width) or tree_width < 1 or not math.isfinite(row_height) or row_height < .15:
        raise ValueError('tree-width must be >=1 inch and row-height >=0.15 inch')
    if support_style not in ('symbols','values','ufboot'):
        raise ValueError('support-style must be symbols, values or ufboot')
    sample_font=11 if compact else 9
    metadata_font=10.5 if compact else 8.5
    gene_font=11 if compact else 9
    group_font=10 if compact else 8.5
    title_font=17 if compact else 15
    note_font=9 if compact else 8
    branch_font=7.2 if compact else 6.2
    gene_spacing=.23 if compact else .25
    group_gap=.12 if compact else .22
    tree = copy.deepcopy(read_tree(tree_path))
    if outgroup:
        from fig1lib import sample_id
        outgroup = sample_id(outgroup)
        if outgroup not in {t.name for t in tree.get_terminals()}:
            raise ValueError(f'Outgroup {outgroup!r} is not a tree tip')
        if root != 'as-is':
            raise ValueError('Choose either --outgroup or --root midpoint')
        tree.root_with_outgroup(outgroup)
        root_description = f'Outgroup: {outgroup}'
    elif root == 'midpoint':
        if any(c.branch_length is None for c in tree.find_clades() if c is not tree.root) or max(tree.depths().values()) <= 0:
            raise ValueError('Midpoint rooting requires complete, positive branch-length distances')
        tree.root_at_midpoint()
        root_description = 'Midpoint display root; ancestral root not inferred'
    else:
        root_description = 'Tree orientation: as supplied'
    if ladderize:
        tree.ladderize()
    rotation_records=apply_node_rotations(tree,rotate_nodes)
    if rotation_records:
        root_description+=f'; {len(rotation_records)} display node rotations'
    samples = [tip.name for tip in tree.get_terminals()]
    calls, genes, missing_rows = read_genes(genes_path, samples, allow_missing)
    metadata, available = read_metadata(metadata_path, samples)
    display_styles=sample_display_styles(samples,sample_labels_path)
    display_names={s:display_styles[s]['label'] for s in samples}
    if metadata_columns:
        columns = metadata_columns
        if set(columns)-set(available):
            raise ValueError(f'Unknown metadata columns: {sorted(set(columns)-set(available))}')
    else:
        columns = [c for c in ['ST','source','group','outcome','biofilm'] if c in available]
    groups = gene_styles(genes, groups_path)
    highlight_groups=list(highlight_groups or [])
    if set(highlight_groups)-set(groups):
        raise ValueError(f'Unknown groups requested for highlighting: {sorted(set(highlight_groups)-set(groups))}')
    positions, blocks, cursor = {}, [], 0.
    for label, group in groups.items():
        block_width = max(gene_spacing*len(group['genes']), text_width_inches(label,group_font)+.1)
        for index, gene in enumerate(group['genes']):
            positions[gene] = cursor + (index+.5)*block_width/len(group['genes'])
        blocks.append((cursor,cursor+block_width,label,group))
        cursor += block_width+group_gap
    gene_width = max(cursor-group_gap, .3)
    label_width = max(1.15, max(text_width_inches(display_names[s],sample_font)*1.13 for s in samples)+.12)
    metadata_widths = []
    for column in columns:
        values = [column] + [metadata[s].get(column,'') or 'NA' for s in samples]
        if column == 'ST':
            values += ['ST'+metadata[s]['ST'] for s in samples if metadata[s].get('ST')]
        metadata_widths.append(max(.55,max(text_width_inches(v,metadata_font)*(1.10 if compact else 1) for v in values)+(.12 if compact else .18)))
    margin, gap = (.30,.08) if compact else (.35,.12)
    label_start = margin+tree_width+gap
    meta_start = label_start+label_width+gap
    gene_start = meta_start+sum(metadata_widths)+((.12 if compact else .22) if columns else 0)
    width = max(gene_start+gene_width+margin,text_width_inches(title,title_font)+2*margin)
    body_height = max(row_height*len(samples), .6)
    header_height = max(.85, max(text_width_inches(g,gene_font) for g in genes)+.45)
    footer_height, title_height = (1.55 if highlight_nosocomial or sample_labels_path or show_branch_lengths else 1.15), .72
    height = footer_height+body_height+header_height+title_height
    fig = plt.figure(figsize=(width,height),facecolor='white')
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':metadata_font,
                         'pdf.fonttype':42,'svg.fonttype':'none'})
    def panel(left, panel_width, header=False):
        bottom = footer_height+body_height if header else footer_height
        size = header_height if header else body_height
        ax = fig.add_axes([left/width,bottom/height,panel_width/width,size/height])
        ax.set_xlim(0,panel_width)
        if header:
            ax.set_ylim(0,header_height)
        else:
            ax.set_ylim(len(samples)-.5,-.5)
        ax.axis('off')
        return ax
    tree_ax = panel(margin,tree_width)
    labels_ax = panel(label_start,label_width)
    gene_ax = panel(gene_start,gene_width)
    gene_header = panel(gene_start,gene_width,True)
    meta_axes, meta_headers = [], []
    x = meta_start
    for column,column_width in zip(columns,metadata_widths):
        meta_axes.append(panel(x,column_width))
        meta_headers.append(panel(x,column_width,True))
        meta_headers[-1].text(column_width/2,.12,column,ha='center',va='bottom',weight='bold',fontsize=metadata_font)
        x += column_width
    for index,sample in enumerate(samples):
        nosocomial=highlight_nosocomial and metadata[sample].get('source','').lower()=='nosocomial'
        style=display_styles[sample]
        custom=style['label']!=sample or style['bold'] or style['color']!='#172b3a'
        for ax,panel_width in [(labels_ax,label_width),(gene_ax,gene_width)] + list(zip(meta_axes,metadata_widths)):
            if nosocomial or custom or index % 2 == 0:
                tint=tuple(.08*channel+.92 for channel in matplotlib.colors.to_rgb(style['color']))
                shade=matplotlib.colors.to_hex(tint) if custom else '#fff0f0' if nosocomial else '#f3f6f9'
                ax.add_patch(Rectangle((0,index-.5),panel_width,1,facecolor=shade,edgecolor='none',zorder=-5))
        labels_ax.text(0,index,display_names[sample],va='center',fontsize=sample_font,
                       color=style['color'] if custom else '#c62828' if nosocomial else style['color'],
                       weight='bold' if style['bold'] or nosocomial else 'normal')
        for ax,column,column_width in zip(meta_axes,columns,metadata_widths):
            value = metadata[sample].get(column,'')
            text = 'ST'+value if column == 'ST' and value else value or 'NA'
            ax.text(column_width/2,index,text,ha='center',va='center',fontsize=metadata_font,
                    color='#c62828' if nosocomial and column=='source' else '#172b3a' if value else '#7a8792',
                    weight='bold' if nosocomial and column=='source' else 'normal')
    for left,right,label,group in blocks:
        if label in highlight_groups:
            gene_ax.add_patch(Rectangle((left,-.5),right-left,len(samples),
                facecolor='#fff2b3',edgecolor='none',zorder=-4,gid='highlight-group-'+label))
            gene_header.add_patch(Rectangle((left,0),right-left,header_height,
                facecolor='#fff2b3',edgecolor='none',zorder=-4,gid='highlight-header-'+label))
        line_y = header_height-.23
        gene_header.plot([left,right],[line_y,line_y],color=group['color'],lw=1.8)
        gene_header.text((left+right)/2,line_y+.07,label,ha='center',va='bottom',fontsize=group_font,color=group['color'])
        for gene in group['genes']:
            gx = positions[gene]
            gene_header.text(gx,.08,gene,rotation=90,ha='center',va='bottom',fontsize=gene_font)
            for index,sample in enumerate(samples):
                value = calls.get((sample,gene))
                if value is None:
                    gene_ax.scatter(gx,index,marker='x',s=30 if compact else 19,color='#8995a0',linewidths=.8)
                else:
                    area=(64 if group['shape']=='square' else 54) if compact else 36
                    gene_ax.scatter(gx,index,marker=MARKERS[group['shape']],s=area,
                                    facecolors=group['color'] if value else 'white',
                                    edgecolors=group['color'],linewidths=.8)
    depths, proportional = tree_display_depths(tree,tree_layout)
    maximum = max(depths.values())
    tip_y = {tip:float(index) for index,tip in enumerate(tree.get_terminals())}
    ys = dict(tip_y)
    support_counts = dict(high=0,lower=0)
    branch_records=[]
    def draw(clade):
        if clade.clades:
            for child in clade.clades:
                draw(child)
            ys[clade] = (ys[clade.clades[0]]+ys[clade.clades[-1]])/2
            tree_ax.plot([depths[clade],depths[clade]],[ys[clade.clades[0]],ys[clade.clades[-1]]],color='#253446',lw=.8)
            for child in clade.clades:
                tree_ax.plot([depths[clade],depths[child]],[ys[child],ys[child]],color='#253446',lw=.8)
                if show_branch_lengths and child.branch_length is not None:
                    branch_records.append(dict(descendants=[display_names[t.name] for t in child.get_terminals()],
                        length=float(child.branch_length),start=depths[clade],end=depths[child],y=ys[child]))
            support = support_label(clade)
            if show_support and support and clade is not tree.root:
                values = support_values(clade)
                high = values[-1] >= 95 and (len(values) == 1 or values[0] >= 80)
                support_counts['high' if high else 'lower'] += 1
                if support_style == 'symbols':
                    tree_ax.scatter(depths[clade],ys[clade],s=16 if compact else 12,
                                    facecolors='#253446' if high else 'white',
                                    edgecolors='#253446',linewidths=.65,zorder=5)
                else:
                    label = f'{values[-1]:g}' if support_style == 'ufboot' else support
                    tree_ax.annotate(label,(depths[clade],ys[clade]),xytext=(2,3),
                                     textcoords='offset points',fontsize=6.5 if compact else 5.5)
    draw(tree.root)
    if proportional:
        for tip in tree.get_terminals():
            tree_ax.plot([depths[tip],maximum*1.035],[ys[tip],ys[tip]],
                         color='#cbd2da',lw=.55,ls=(0,(1.5,2.5)),zorder=-2)
    tree_ax.set_xlim(-maximum*.02,maximum*1.04)
    if show_branch_lengths:
        branch_distance_labels(fig,tree_ax,branch_records,depths,ys,branch_font)
    if proportional:
        scale = maximum/4
        tree_ax.plot([0,scale],[0,0],transform=tree_ax.get_xaxis_transform(),color='#253446',lw=1,clip_on=False)
        tree_ax.annotate(f'{scale:.2g} {branch_units}',(0,0),xytext=(0,-8 if show_branch_lengths else -20),
                         xycoords=tree_ax.get_xaxis_transform(),textcoords='offset points',ha='left',fontsize=8 if compact else 7)
        scale_note = branch_units
    else:
        scale_note = ('Cladogram: branch lengths not displayed' if tree_layout == 'cladogram'
                      else 'Cladogram: branch lengths unavailable or zero')
        fig.text(margin/width,.86/height,scale_note,fontsize=note_font,color='#586372')
    if show_support and support_style == 'symbols' and sum(support_counts.values()):
        fig.text(margin/width,.64/height,
                 'Internal nodes: filled = UFBoot >=95 and SH-aLRT >=80 when paired; open = below these thresholds.',
                 fontsize=note_font,color='#586372')
    if proportional:
        fig.text((margin+tree_width+.3)/width,.86/height,
                 'Dotted lines connect tips to labels; they do not represent branch length.',fontsize=note_font,color='#586372')
    fig.text(margin/width,(height-.27)/height,title,fontsize=title_font,weight='bold',va='top')
    fig.text(margin/width,(height-.54)/height,f'{len(samples)} genomes | {len(genes)} markers | {root_description}',fontsize=10.5 if compact else 9,color='#586372',va='top')
    fig.text(margin/width,.36/height,'Filled = reported present    Open = reported absent    Cross = missing / unknown',fontsize=9.5 if compact else 8.5)
    fig.text(margin/width,.16/height,'Gene calls and sample metadata are plotted as supplied; phenotypes are not inferred from genomes.',fontsize=note_font,color='#586372')
    if highlight_nosocomial or sample_labels_path:
        fig.text(margin/width,1.08/height,'Red bold names/source and pale red rows: nosocomial isolates.  Colored bold name: highlighted isolate.',
                 fontsize=note_font,color='#586372')
    if show_branch_lengths:
        fig.text(margin/width,1.30/height,f'Numbers on branches: {branch_units} (3 significant figures); leaders connect displaced labels to branch midpoints.',
                 fontsize=note_font,color='#586372')
    output = Path(output)
    output.parent.mkdir(parents=True,exist_ok=True)
    for suffix in ['pdf','svg','png']:
        fig.savefig(str(output)+'.'+suffix,dpi=180)
    plt.close(fig)
    ordered_genes = [g for _,_,_,group in blocks for g in group['genes']]
    import csv
    with Path(str(output)+'.matrix.tsv').open('w',newline='') as stream:
        writer = csv.writer(stream,delimiter='\t')
        writer.writerow(['sample']+ordered_genes)
        for sample in samples:
            writer.writerow([display_names[sample]]+['NA' if calls.get((sample,g)) is None else calls[sample,g] for g in ordered_genes])
    display_tree=copy.deepcopy(tree)
    for tip in display_tree.get_terminals():tip.name=display_names[tip.name]
    Phylo.write(display_tree,str(output)+'.display_tree.nwk','newick',format_branch_length='%1.10g')
    write_metadata(str(output)+'.metadata.tsv',{display_names[s]:metadata[s] for s in samples},columns)
    if show_branch_lengths:save_json(str(output)+'.branch_lengths.json',branch_records)
    description = dict(sample_order=[display_names[s] for s in samples],input_sample_order=samples,
                       node_rotations=rotation_records,
                       compact_layout=compact,
                       font_points=dict(sample=sample_font,metadata=metadata_font,gene=gene_font,group=group_font,branch_distance=branch_font),
                       marker_area_points2=dict(square=64,circle=54,triangle=54) if compact else dict(square=36,circle=36,triangle=36),
                       gene_spacing_inches=gene_spacing,group_gap_inches=group_gap,
                       sample_display_styles=display_styles,highlight_nosocomial=highlight_nosocomial,
                       highlighted_gene_groups=highlight_groups,gene_group_highlight_color='#fff2b3' if highlight_groups else None,
                       labeled_branch_distances=len(branch_records),gene_order=ordered_genes,metadata_columns=columns,
                       metadata_source=str(Path(metadata_path).resolve()) if metadata_path else None,
                       missing_rows=missing_rows,unknown_calls=sum(calls.get((s,g)) is None for s in samples for g in genes),
                       rooting=root_description,branch_units=scale_note,
                       figure_inches=[width,height],source_tree=str(Path(tree_path).resolve()),
                       tree_layout='phylogram' if proportional else 'cladogram',
                       source_branch_units=branch_units,tree_width_inches=tree_width,
                       row_height_inches=row_height,ladderized=ladderize,
                       support_style=support_style if show_support else 'hidden',support_counts=support_counts)
    save_json(str(output)+'.plot.json',description)
    print(f'Wrote {output}.pdf, .svg, .png; {len(samples)} samples, {len(genes)} markers.')
    return description


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--tree',type=Path,required=True)
    parser.add_argument('--genes',type=Path,required=True)
    parser.add_argument('--metadata',type=Path)
    parser.add_argument('--gene-groups',type=Path)
    parser.add_argument('--output',type=Path,required=True,help='Output prefix, without an extension')
    parser.add_argument('--allow-missing',action='store_true')
    parser.add_argument('--metadata-columns',help='Comma-separated column names; default ST,source,group,outcome,biofilm when present')
    parser.add_argument('--title',default='Genome phylogeny and marker profiles')
    parser.add_argument('--root',choices=['as-is','midpoint'],default='as-is')
    parser.add_argument('--outgroup',help='An existing sample in the tree; does not add genomes')
    parser.add_argument('--branch-units',default='as supplied')
    parser.add_argument('--show-support',action='store_true')
    parser.add_argument('--tree-layout',choices=['phylogram','cladogram'],default='phylogram')
    parser.add_argument('--tree-width',type=float,default=4.4,help='Width of the tree panel in inches')
    parser.add_argument('--row-height',type=float,default=.28,help='Row spacing in inches')
    parser.add_argument('--keep-tip-order',action='store_true',help='Keep supplied tip order instead of ladderizing')
    parser.add_argument('--support-style',choices=['symbols','values','ufboot'],default='symbols')
    parser.add_argument('--sample-labels',type=Path,help='TSV with sample,label,color,bold; input identifiers remain unchanged')
    parser.add_argument('--highlight-nosocomial',action='store_true')
    parser.add_argument('--show-branch-lengths',action='store_true',help='Label original edge distances in either tree layout')
    parser.add_argument('--highlight-groups',help='Comma-separated gene-group names to shade pale yellow')
    parser.add_argument('--rotate-node',action='append',help='Reverse child order at the MRCA of comma-separated input sample names; repeat for multiple nodes')
    parser.add_argument('--compact',action='store_true',help='Larger text and symbols with tighter gene/panel spacing')
    args = parser.parse_args()
    try:
        render(args.tree,args.genes,args.output,args.metadata,args.gene_groups,args.allow_missing,
               args.title,args.metadata_columns.split(',') if args.metadata_columns else None,
               args.root,args.outgroup,args.branch_units,args.show_support,
               args.tree_layout,args.tree_width,args.row_height,not args.keep_tip_order,args.support_style,
               args.sample_labels,args.highlight_nosocomial,args.show_branch_lengths,
               args.highlight_groups.split(',') if args.highlight_groups else None,args.rotate_node,args.compact)
    except (ValueError,OSError) as error:
        parser.exit(2,f'Error: {error}\n')


if __name__ == '__main__':
    main()

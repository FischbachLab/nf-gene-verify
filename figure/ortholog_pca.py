#!/usr/bin/env python3
"""PCA of genome-wide Panaroo gene-cluster presence/absence, with sample metadata."""
import argparse
import csv
import math
import os
from pathlib import Path

os.environ.setdefault('MPLCONFIGDIR',str(Path(os.environ.get('TMPDIR','/tmp'))/'fig1-matplotlib'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Ellipse
from matplotlib.text import Text
from matplotlib.transforms import Bbox
import numpy as np

from fig1lib import read_metadata,read_tsv,sample_id,save_json,sha256,MISSING


def read_rtab(path):
    with Path(path).open(newline='',encoding='utf-8-sig') as stream:
        reader=csv.reader(stream,delimiter='\t')
        header=next(reader,[])
        if len(header)<4 or header[0].strip()!='Gene':
            raise ValueError('Expected a Panaroo Rtab: Gene followed by at least three sample columns')
        samples=[sample_id(value) for value in header[1:]]
        if len(set(samples))!=len(samples):
            raise ValueError('Duplicated sample identifier in ortholog matrix')
        genes,rows,seen=[],[],set()
        for line,row in enumerate(reader,2):
            if len(row)!=len(header) or not row[0].strip():
                raise ValueError(f'Malformed ortholog matrix row {line}')
            gene=row[0].strip()
            if gene in seen:
                raise ValueError(f'Duplicated gene cluster: {gene}')
            values=[v.strip() for v in row[1:]]
            if any(v not in ('0','1') for v in values):
                raise ValueError(f'Ortholog matrix must be binary with no missing values; row {line}')
            seen.add(gene);genes.append(gene);rows.append([int(v) for v in values])
    if not genes:
        raise ValueError('Ortholog matrix contains no gene clusters')
    return samples,genes,np.asarray(rows,dtype=np.float64).T


def fit_pca(matrix,min_presence=1,scale=False):
    n,p=matrix.shape
    if n<3 or not 1<=min_presence<n:
        raise ValueError('PCA requires at least three samples and min-presence between 1 and n-1')
    prevalence=matrix.sum(axis=0).astype(int)
    selected=(prevalence>=min_presence)&(prevalence<n)
    if selected.sum()<2:
        raise ValueError('Fewer than two variable clusters remain after filtering')
    x=matrix[:,selected]
    center=x.mean(axis=0)
    divisor=x.std(axis=0,ddof=1) if scale else np.ones(x.shape[1])
    x=(x-center)/divisor
    u,s,vt=np.linalg.svd(x,full_matrices=False)
    tolerance=max(x.shape)*np.finfo(float).eps*s[0]
    rank=int(np.sum(s>tolerance))
    if rank<2:
        raise ValueError('Gene-content matrix has fewer than two independent principal components')
    scores=u[:,:rank]*s[:rank]
    loadings=vt[:rank].T.copy()
    for index in range(rank):
        # PCA axis signs are arbitrary. Fix each using its largest absolute loading.
        if loadings[np.argmax(np.abs(loadings[:,index])),index]<0:
            scores[:,index]*=-1;loadings[:,index]*=-1
    eigenvalues=s[:rank]**2/(n-1)
    ratios=s[:rank]**2/np.sum(s**2)
    reconstruction=float(np.linalg.norm(x-scores@loadings.T)/np.linalg.norm(x))
    if reconstruction>1e-10:
        raise ValueError('PCA numerical reconstruction check failed')
    return dict(selected=selected,prevalence=prevalence,center=center,divisor=divisor,
                scores=scores,loadings=loadings,eigenvalues=eigenvalues,ratios=ratios,
                rank=rank,reconstruction_relative_error=reconstruction)


def checked_metadata(path,samples):
    metadata,columns=read_metadata(path,samples)
    if path is not None:
        fields,rows=read_tsv(path,[])
        key=next(c for c in ('sample','genome','strainID') if c in fields)
        if {sample_id(row[key]) for row in rows}!=set(samples):
            raise ValueError('PCA metadata must contain exactly one row per input sample; use explicit NA for unknown values')
    return metadata,columns


def source_value(row,column):
    value=row.get(column,'').strip()
    if value.upper() in MISSING:
        return 'Not supplied'
    if value.lower() in ('commensal','nosocomial'):
        return value.lower()
    return value


def point_style(row,column,other_colors):
    source=source_value(row,column)
    if source=='commensal':
        return 'Commensal','#17212b','white'
    if source=='nosocomial':
        if not row.get('ST'):
            return 'Nosocomial, ST unknown','#df3232','white'
        return ('Nosocomial, ST2','#df3232','#df3232') if row.get('ST')=='2' else ('Nosocomial, other ST','#df3232','white')
    if source=='Not supplied':
        return source,'#7a8792','white'
    return source,other_colors[source],'white'


def place_sample_labels(fig,ax,coordinates,samples,reserved=()):
    """Offset labels, never data points, and check label/point collisions."""
    fig.canvas.draw()
    renderer=fig.canvas.get_renderer()
    xy=ax.transData.transform(coordinates)
    allowed=ax.get_window_extent(renderer).padded(-5)
    point_boxes=[Bbox.from_bounds(x-5,y-5,10,10) for x,y in xy]
    occupancy=[artist.get_window_extent(renderer).padded(3) for artist in reserved]
    crowding=np.sum(np.linalg.norm(xy[:,None,:]-xy[None,:,:],axis=2)<50,axis=1)
    order=sorted(range(len(samples)),key=lambda i:(-crowding[i],samples[i]))
    angles=[45,-45,135,-135,0,90,180,-90,25,65,115,155,205,245,295,335]
    for index in order:
        annotation=ax.annotate(samples[index],coordinates[index],xytext=(8,8),
            textcoords='offset points',fontsize=8.5,color='#596b7b',
            arrowprops=dict(arrowstyle='-',color='#b7c1ca',lw=.5),zorder=4)
        found=False
        for radius in [8,12,18,25,34,45,60,78,100,125,155,190]:
            for angle in angles:
                radians=math.radians(angle)
                dx,dy=radius*math.cos(radians),radius*math.sin(radians)
                annotation.set_position((dx,dy))
                annotation.set_ha('left' if dx>=0 else 'right')
                annotation.set_va('bottom' if dy>=0 else 'top')
                annotation.update_positions(renderer)
                bounds=Text.get_window_extent(annotation,renderer).padded(2)
                if not (allowed.contains(bounds.x0,bounds.y0) and allowed.contains(bounds.x1,bounds.y1)):
                    continue
                if any(bounds.overlaps(other) for other in occupancy+point_boxes):
                    continue
                occupancy.append(bounds);found=True;break
            if found:break
        if not found:
            raise ValueError(f'Cannot place sample label {samples[index]} without overlap; use --no-labels for a crowded dataset')


def descriptive_outlines(ax,coordinates,samples,metadata,column):
    drawn=[]
    groups={}
    for index,sample in enumerate(samples):
        value=metadata[sample].get(column,'')
        if value.upper() not in MISSING:
            groups.setdefault(value,[]).append(index)
    for name,indices in groups.items():
        if len(indices)<3:continue
        points=coordinates[indices]
        covariance=np.cov(points,rowvar=False)
        eigenvalues,eigenvectors=np.linalg.eigh(covariance)
        if eigenvalues[0]<=np.finfo(float).eps*max(1,eigenvalues[-1]):continue
        mean=points.mean(axis=0)
        standardized=(points-mean)@eigenvectors/np.sqrt(eigenvalues)
        radius=np.sqrt(np.max(np.sum(standardized**2,axis=1)))*1.10
        angle=np.degrees(np.arctan2(eigenvectors[1,1],eigenvectors[0,1]))
        ax.add_patch(Ellipse(mean,2*radius*np.sqrt(eigenvalues[1]),2*radius*np.sqrt(eigenvalues[0]),
            angle=angle,facecolor='#64748b',alpha=.06,edgecolor='none',zorder=0))
        ax.add_patch(Ellipse(mean,2*radius*np.sqrt(eigenvalues[1]),2*radius*np.sqrt(eigenvalues[0]),
            angle=angle,facecolor='none',edgecolor='#64748b',lw=1.1,ls='--',zorder=0))
        # Exact axis-aligned bounds of the rotated ellipse, including its padding.
        extent=radius*np.sqrt((eigenvectors**2)@eigenvalues)
        top=mean+radius*covariance[:,1]/np.sqrt(covariance[1,1])
        ax.annotate(name,top,xytext=(0,5),textcoords='offset points',
                    ha='center',va='bottom',fontsize=13,weight='bold',color='#475569',
                    bbox=dict(facecolor='white',edgecolor='none',alpha=.85,pad=1.5),
                    zorder=5,gid='pca-group-label')
        drawn.append(dict(group=name,samples=len(indices),center=mean.tolist(),
                          bounds=[(mean-extent).tolist(),(mean+extent).tolist()],
                          description='Data-enclosing descriptive ellipse; not a confidence region'))
    return drawn


def analyze(matrix_path,output,metadata_path=None,source_column='source',group_column=None,
            min_presence=1,scale=False,labels=True,expected_samples=None):
    samples,genes,matrix=read_rtab(matrix_path)
    if expected_samples is not None and set(samples)!=set(expected_samples):
        raise ValueError('PCA matrix sample IDs do not match the genome dataset')
    metadata,columns=checked_metadata(metadata_path,samples)
    if group_column and group_column not in columns:
        raise ValueError(f'Metadata lacks requested group column: {group_column}')
    fit=fit_pca(matrix,min_presence,scale)
    selected_genes=[g for g,keep in zip(genes,fit['selected']) if keep]
    coordinates=fit['scores'][:,:2]
    output=Path(output);output.parent.mkdir(parents=True,exist_ok=True)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'pdf.fonttype':42,'svg.fonttype':'none'})
    fig=plt.figure(figsize=(10.4,8.8),facecolor='white')
    ax=fig.add_axes([.10,.17,.72,.70])
    span=np.ptp(coordinates,axis=0)
    for axis in range(2):
        low,high=coordinates[:,axis].min(),coordinates[:,axis].max()
        limits=(low-.30*span[axis],high+.30*span[axis])
        (ax.set_xlim if axis==0 else ax.set_ylim)(limits)
    ax.set_aspect('equal',adjustable='box')
    for spine in ax.spines.values():spine.set_color('#aeb9c3')
    ax.axhline(0,color='#e2e7ec',lw=.7,zorder=-3)
    ax.axvline(0,color='#e2e7ec',lw=.7,zorder=-3)
    outline_records=descriptive_outlines(ax,coordinates,samples,metadata,group_column) if group_column else []
    if outline_records:
        bounds=np.asarray([record['bounds'] for record in outline_records])
        low=np.minimum(coordinates.min(axis=0),bounds[:,0,:].min(axis=0))
        high=np.maximum(coordinates.max(axis=0),bounds[:,1,:].max(axis=0))
        padding=.10*(high-low)
        ax.set_xlim(low[0]-padding[0],high[0]+padding[0])
        ax.set_ylim(low[1]-padding[1],high[1]+padding[1])
    extra=sorted({source_value(row,source_column) for row in metadata.values()}-{'commensal','nosocomial','Not supplied'})
    other_colors={value:matplotlib.colors.to_hex(plt.get_cmap('tab10')(i%10)) for i,value in enumerate(extra)}
    legend={}
    category_counts={}
    for index,sample in enumerate(samples):
        label,edge,face=point_style(metadata[sample],source_column,other_colors)
        ax.scatter(*coordinates[index],s=51,marker='o',edgecolors=edge,facecolors=face,linewidths=1.1,zorder=3)
        legend[label]=Line2D([],[],marker='o',ls='none',markeredgecolor=edge,markerfacecolor=face,markersize=6)
        category_counts[label]=category_counts.get(label,0)+1
    ordered=[v for v in ['Commensal','Nosocomial, other ST','Nosocomial, ST2','Not supplied'] if v in legend]
    ordered += [v for v in legend if v not in ordered]
    source_legend=ax.legend([legend[v] for v in ordered],[f'{v} (n={category_counts[v]})' for v in ordered],
              frameon=False,fontsize=8.5,loc='upper left',borderaxespad=.8)
    ax.set_xlabel(f'PC1 ({100*fit["ratios"][0]:.2f}% variance)',labelpad=9)
    ax.set_ylabel(f'PC2 ({100*fit["ratios"][1]:.2f}% variance)',labelpad=9)
    if labels:
        reserved=[source_legend]+[text for text in ax.texts if text.get_gid()=='pca-group-label']
        place_sample_labels(fig,ax,coordinates,samples,reserved)
    fig.text(.10,.965,'PCA of ortholog-cluster content',fontsize=16,weight='bold',va='top')
    fig.text(.10,.925,f'{len(samples)} genomes | {len(genes):,} Panaroo clusters | {len(selected_genes):,} variable clusters analyzed',
             fontsize=10,color='#596b7b',va='top')
    fig.text(.10,.095,'Binary gene-cluster presence/absence; columns centered; '+('unit-variance scaled.' if scale else 'no variance scaling.'),fontsize=8.5,color='#596b7b')
    fig.text(.10,.067,'Labels are offset for readability; sample coordinates are not jittered.',fontsize=8.5,color='#596b7b')
    if group_column:
        fig.text(.10,.039,'Dashed group outlines are descriptive, not confidence regions.',fontsize=8.5,color='#596b7b')
    for extension in ['pdf','svg','png']:
        fig.savefig(str(output)+'.'+extension,dpi=200)
    plt.close(fig)

    def table(suffix,header,rows):
        with Path(str(output)+suffix).open('w',newline='') as stream:
            writer=csv.writer(stream,delimiter='\t',lineterminator='\n');writer.writerow(header);writer.writerows(rows)
    pcs=['PC'+str(i+1) for i in range(fit['rank'])]
    table('.scores.tsv',['sample']+columns+pcs,
          ([s]+[metadata[s].get(c,'') or 'NA' for c in columns]+list(fit['scores'][i]) for i,s in enumerate(samples)))
    table('.loadings.tsv',['gene','center','scale']+pcs,
          ([g,fit['center'][i],fit['divisor'][i]]+list(fit['loadings'][i]) for i,g in enumerate(selected_genes)))
    table('.variance.tsv',['component','eigenvalue','explained_variance_ratio','cumulative_ratio'],
          ([pc,fit['eigenvalues'][i],fit['ratios'][i],np.cumsum(fit['ratios'])[i]] for i,pc in enumerate(pcs)))
    table('.matrix.tsv',['sample']+selected_genes,
          ([s]+list(matrix[i,fit['selected']].astype(int)) for i,s in enumerate(samples)))
    table('.features.tsv',['gene','present_in_samples','included','reason'],
          ([g,int(fit['prevalence'][i]),int(fit['selected'][i]),
            'included' if fit['selected'][i] else 'constant_present' if fit['prevalence'][i]==len(samples)
            else 'constant_absent' if fit['prevalence'][i]==0 else 'below_min_presence'] for i,g in enumerate(genes)))
    summary=dict(samples=len(samples),sample_order=samples,total_clusters=len(genes),
        analyzed_clusters=len(selected_genes),constant_clusters=int(np.sum((fit['prevalence']==0)|(fit['prevalence']==len(samples)))),
        singleton_clusters=int(np.sum(fit['prevalence']==1)),min_presence=min_presence,
        centered=True,unit_variance_scaled=scale,rank=fit['rank'],
        PC1_variance_percent=float(100*fit['ratios'][0]),PC2_variance_percent=float(100*fit['ratios'][1]),
        category_counts=category_counts,source_column=source_column,source_column_present=source_column in columns,
        group_column=group_column,descriptive_outlines=outline_records,
        matrix_source=str(Path(matrix_path).resolve()),matrix_sha256=sha256(matrix_path),
        metadata_source=str(Path(metadata_path).resolve()) if metadata_path else None,
        metadata_sha256=sha256(metadata_path) if metadata_path else None,
        method='SVD PCA of centered genome-by-Panaroo-cluster binary presence/absence matrix',
        reconstruction_relative_error=fit['reconstruction_relative_error'],
        reference='Conlan et al. 2012 Figure 6; adaptation using Panaroo rather than OrthoMCL; not an exact reproduction',
        source_and_group_labels_used_for_plotting_only=True)
    save_json(str(output)+'.summary.json',summary)
    print(f'PCA: {len(samples)} genomes, {len(selected_genes)} variable clusters; PC1 {summary["PC1_variance_percent"]:.2f}%, PC2 {summary["PC2_variance_percent"]:.2f}%.')
    print(f'Wrote {output}.pdf, .svg and .png plus scores, loadings, variance, matrix and feature tables.')
    return summary


def main():
    package=Path(__file__).resolve().parent
    default_metadata=package/'metadata2.tsv'
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--results-dir',type=Path,default=package/'results_29')
    parser.add_argument('--rtab',type=Path,help='Override the results folder Panaroo Rtab path')
    parser.add_argument('--metadata',type=Path,default=default_metadata if default_metadata.is_file() else None)
    parser.add_argument('--output',type=Path,help='Output prefix; default sibling results_FOLDER_pca/figure6_pca')
    parser.add_argument('--source-column',default='source')
    parser.add_argument('--group-column',help='Optional supplied group column for descriptive outlines; groups are never inferred')
    parser.add_argument('--min-presence',type=int,default=1,help='Minimum genomes containing a cluster; 2 excludes singletons')
    parser.add_argument('--scale',action='store_true',help='Scale variable clusters to unit variance; default is centered only')
    parser.add_argument('--no-labels',action='store_true')
    args=parser.parse_args()
    rtab=args.rtab or args.results_dir/'work/panaroo/pangenome/gene_presence_absence.Rtab'
    output=args.output or args.results_dir.parent/(args.results_dir.name+'_pca')/'figure6_pca'
    try:
        analyze(rtab,output,args.metadata,args.source_column,args.group_column,args.min_presence,args.scale,not args.no_labels)
    except (OSError,ValueError) as error:
        parser.exit(2,f'Error: {error}\n')


if __name__=='__main__':
    main()

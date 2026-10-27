"""Standalone command line entry point; outputs an offline interactive plot."""
from pathlib import Path
import argparse, platform, time, importlib.metadata
import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist,squareform
from scipy.cluster.hierarchy import dendrogram
from tfcds_core import *

def interactive_plot(x,pops,z,target,title):
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots
    dd=dendrogram(z,no_plot=True)
    order=dd['leaves'];fig=make_subplots(rows=1,cols=2,column_widths=[.2,.8],horizontal_spacing=.015)
    for yy,xx in zip(dd['icoord'],dd['dcoord']):
        fig.add_trace(go.Scatter(x=xx,y=[(a-5)/10 for a in yy],mode='lines',line=dict(color='#395770',width=1.5),hoverinfo='skip',showlegend=False),row=1,col=1)
    fig.add_trace(go.Heatmap(x=list(range(x.shape[1])),y=list(range(len(pops))),z=x[order],colorscale='Magma',
                  colorbar=dict(title='Profile value'),hovertemplate='Population index %{y}<br>Input start position %{x}<br>Value %{z}<extra></extra>'),row=1,col=2)
    fig.update_xaxes(autorange='reversed',title_text='Complete-linkage distance',row=1,col=1)
    fig.update_yaxes(visible=False,row=1,col=1)
    fig.update_yaxes(tickvals=list(range(len(pops))),ticktext=[pops[i] for i in order],side='right',row=1,col=2)
    fig.update_xaxes(title_text='Zero-based start position in each input sequence',row=1,col=2)
    fig.update_layout(title=title,height=max(450,len(pops)*28),width=1400,template='plotly_white',margin=dict(r=160))
    fig.write_html(target,include_plotlyjs=True,full_html=True)

def main():
    ap=argparse.ArgumentParser(description='Sequence k-mer profiling and hierarchical clustering.')
    ap.add_argument('--version',action='version',version=f'Genetical Mapping {VERSION}')
    ap.add_argument('--fasta',required=True,type=Path);ap.add_argument('--out',required=True,type=Path)
    ap.add_argument('--k',type=int,default=20);ap.add_argument('--lower',type=int,default=1);ap.add_argument('--upper',type=int)
    ap.add_argument('--clusters',type=int,default=2)
    ap.add_argument('--representation',choices=['position_raw','position_l1','kmer_frequency'],default='position_raw')
    ap.add_argument('--variant',choices=['nonempty','exact_unique','gene_representative','complete_only','complete_unique'],default='nonempty')
    ap.add_argument('--legacy-end-drop',action='store_true');ap.add_argument('--canonical',action='store_true')
    args=ap.parse_args()
    if args.out.exists() and any(args.out.iterdir()):raise SystemExit('Output directory must be new or empty; existing files are never silently replaced.')
    if args.lower<1 or (args.upper is not None and args.upper<args.lower):raise SystemExit('Invalid occurrence bounds')
    rr0=read_fasta(args.fasta);rr=choose_records(rr0,args.variant)
    if set(r.population for r in rr)!=set(r.population for r in rr0):raise SystemExit('Filtering removes an entire population; inspect inputs instead of silent exclusion.')
    t=time.perf_counter();data=ProfileData(rr,args.k,args.canonical)
    if not 2<=args.clusters<data.n:raise SystemExit('clusters must be >=2 and < number of populations')
    if args.representation=='kmer_frequency':
        if args.legacy_end_drop:raise SystemExit('legacy-end-drop applies only to positional profiles')
        matrix,ret=data.frequency(lower=args.lower,upper=args.upper)
        if np.any(np.asarray(matrix.sum(axis=1)).ravel()==0):raise SystemExit('Zero-profile population after filtering')
        d=frequency_distance(matrix)
    else:
        matrix,ret=data.positional(lower=args.lower,upper=args.upper,legacy_end_drop=args.legacy_end_drop)
        if np.any(matrix.sum(axis=1)==0):raise SystemExit('Zero-profile population after filtering')
        if args.representation=='position_l1':matrix=row_l1(matrix)
        d=squareform(pdist(matrix))
    info,labels,z=analyze_distance(d,args.clusters)
    args.out.mkdir(parents=True,exist_ok=True)
    pd.DataFrame({'population':data.populations,'cluster':labels}).to_csv(args.out/'cluster_assignments.tsv',sep='\t',index=False)
    pd.DataFrame(d,index=data.populations,columns=data.populations).to_csv(args.out/'distance_matrix.tsv',sep='\t')
    np.save(args.out/'linkage.npy',z)
    if sparse.issparse(matrix):sparse.save_npz(args.out/'feature_matrix.npz',matrix)
    else:
        pd.DataFrame(matrix,index=data.populations).to_csv(args.out/'feature_matrix.tsv',sep='\t')
        interactive_plot(matrix,data.populations,z,args.out/'interactive_heatmap.html',f'Sequence profiling {VERSION}; {args.representation}')
    pd.DataFrame({'kmer':np.array(data.words)[ret],'pooled_count':data.full_pool[ret]}).to_csv(args.out/'retained_kmers.tsv',sep='\t',index=False)
    pd.DataFrame(data.support(),index=data.populations).to_csv(args.out/'position_support.tsv',sep='\t')
    manifest={k:str(v) if isinstance(v,Path) else v for k,v in vars(args).items()}
    manifest.update(info,input_sha256=sha256(args.fasta),version=VERSION,python=platform.python_version(),
                    n_input_headers=len(rr0),n_selected_records=len(rr),n_retained_kmers=int(ret.sum()),feature_dimensions=list(matrix.shape),
                    seconds=time.perf_counter()-t,metric='sqrt-frequency Euclidean' if args.representation=='kmer_frequency' else 'Euclidean',linkage_method='complete',
                    populations=data.populations,dependencies={n:importlib.metadata.version(n) for n in ['numpy','scipy','pandas','scikit-learn','plotly']})
    json_write(args.out/'run_manifest.json',manifest)
    print('SUCCESS',args.out)
if __name__=='__main__':main()

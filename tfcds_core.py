"""Sequence k-mer profiling and complete-linkage clustering.

Positional profiles sum pooled occurrence weights of distinct k-mers present
in each population/start-position cell. The legacy_end_drop option sets each
population's last available start position to zero for compatibility.
"""
from __future__ import annotations
from dataclasses import dataclass, replace
from pathlib import Path
import hashlib, json, re
import numpy as np
from scipy import sparse
from scipy.cluster.hierarchy import linkage, cut_tree
from scipy.spatial.distance import pdist, squareform
from sklearn.metrics import adjusted_rand_score, silhouette_score

VERSION = '1.1'

@dataclass(frozen=True)
class Record:
    id: str
    population: str
    sequence: str
    header: str
    gene: str
    orf_type: str
    header_line: int

def sha256(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()

def read_fasta(path):
    records=[]; header=None; seq=[]; line=0; seen=set()
    def emit():
        if header is None:return
        rid=header.split()[0]
        if rid in seen:raise ValueError('Duplicate FASTA identifier: '+rid)
        seen.add(rid)
        sequence=''.join(seq).upper()
        bad=set(sequence)-set('ACGT')
        if bad:raise ValueError(f'{rid}: unsupported symbols {sorted(bad)}')
        gene=re.sub(r'_i\d+(?:\.p\d+)?$', '',rid)
        typ=re.search(r'type:(\S+)',header)
        records.append(Record(rid,rid.split('_',1)[0],sequence,header,gene,typ.group(1) if typ else 'unknown',line))
    for no,s in enumerate(Path(path).read_text(encoding='utf-8-sig').splitlines(),1):
        s=s.strip()
        if s.startswith('>'):
            emit();header=s[1:];seq=[];line=no
        elif s:
            if header is None:raise ValueError('Sequence before first FASTA header')
            seq.append(s)
    emit()
    return records

def write_fasta(records,path):
    with open(path,'w',encoding='utf-8',newline='\n') as f:
        for r in records:
            f.write('>'+r.header+'\n')
            for j in range(0,len(r.sequence),80):f.write(r.sequence[j:j+80]+'\n')

def choose_records(records,variant='nonempty'):
    rr=[r for r in records if r.sequence]
    if variant.startswith('complete'):rr=[r for r in rr if r.orf_type=='complete']
    if variant in ['exact_unique','complete_unique']:
        seen=set();keep=[]
        for r in rr:
            key=(r.population,r.sequence)
            if key not in seen:seen.add(key);keep.append(r)
        rr=keep
    if variant=='gene_representative':
        groups={}
        for r in rr:groups.setdefault((r.population,r.gene),[]).append(r)
        rr=[sorted(v,key=lambda r:(r.orf_type!='complete',-len(r.sequence),r.id))[0] for v in groups.values()]
    if variant not in ['nonempty','complete_only','complete_unique','exact_unique','gene_representative']:
        raise ValueError('Unknown variant '+variant)
    return rr

class ProfileData:
    """Sparse sufficient statistics permit full-pipeline record subsampling."""
    def __init__(self,records,k=20,canonical=False,populations=None):
        if not 1<=k<=64:raise ValueError('k must lie between 1 and 64')
        self.records=records;self.k=k;self.canonical=canonical
        self.populations=populations or list(dict.fromkeys(r.population for r in records))
        self.n=len(self.populations);self.nr=len(records)
        popmap={p:i for i,p in enumerate(self.populations)}
        self.pop_index=np.array([popmap[r.population] for r in records],dtype=np.int32)
        self.lengths=np.array([max(0,len(r.sequence)-k+1) for r in records],dtype=np.int32)
        self.width=int(self.lengths.max(initial=0))
        if not self.width:raise ValueError('No valid kmer windows')
        words={};tokens=[];positions=[];record_ids=[]
        trans=str.maketrans('ACGT','TGCA')
        for i,r in enumerate(records):
            for j in range(self.lengths[i]):
                t=r.sequence[j:j+k]
                if canonical:t=min(t,t.translate(trans)[::-1])
                idx=words.get(t)
                if idx is None:idx=len(words);words[t]=idx
                tokens.append(idx);positions.append(j);record_ids.append(i)
        self.words=list(words);v=len(words);self.v=v
        tok=np.asarray(tokens,dtype=np.int32);pos=np.asarray(positions,dtype=np.int32);rid=np.asarray(record_ids,dtype=np.int32)
        self.record_counts=sparse.csr_matrix((np.ones(len(tok),dtype=np.int32),(rid,tok)),shape=(self.nr,v))
        cells=self.pop_index[rid].astype(np.int64)*self.width+pos
        code=cells*v+tok
        triple,inv=np.unique(code,return_inverse=True)
        self.triple_token=(triple%v).astype(np.int32)
        self.triple_cell=(triple//v).astype(np.int32)
        self.triple_records=sparse.csr_matrix((np.ones(len(tok),dtype=np.int8),(inv,rid)),shape=(len(triple),self.nr))
        self.triple_records.data[:]=1
        self.full_pool=np.asarray(self.record_counts.sum(axis=0)).ravel().astype(np.float64)
        self.full_selection=np.ones(self.nr,dtype=np.int32)

    def counts_by_population(self,selection=None):
        selection=self.full_selection if selection is None else np.asarray(selection,dtype=np.int32)
        group=sparse.csr_matrix((selection,(self.pop_index,np.arange(self.nr))),shape=(self.n,self.nr))
        return (group@self.record_counts).tocsr()

    def pooled(self,selection=None):
        if selection is None:return self.full_pool
        return np.asarray(self.record_counts.T@selection).ravel().astype(np.float64)

    def support(self,selection=None):
        sel=self.full_selection if selection is None else selection
        out=np.zeros((self.n,self.width),dtype=np.int32)
        for i,p in enumerate(self.pop_index):
            if sel[i]:out[p,:self.lengths[i]]+=int(sel[i])
        return out

    def positional(self,selection=None,lower=1,upper=None,legacy_end_drop=False,fixed_pool=False):
        sel=self.full_selection if selection is None else np.asarray(selection,dtype=np.int32)
        pool=self.full_pool if fixed_pool else self.pooled(selection)
        retained=(pool>=lower)&(pool<=(np.inf if upper is None else upper))&(pool>0)
        active=np.asarray(self.triple_records@sel).ravel()>0
        active &= retained[self.triple_token]
        x=np.bincount(self.triple_cell[active],weights=pool[self.triple_token[active]],minlength=self.n*self.width).reshape(self.n,self.width)
        if legacy_end_drop:
            for p in range(self.n):
                available=self.lengths[(self.pop_index==p)&(sel>0)]
                end=int(available.max(initial=0))-1
                if end>=0:x[p,end]=0
        return x,retained

    def frequency(self,selection=None,lower=1,upper=None,equal_cds=False):
        sel=self.full_selection if selection is None else np.asarray(selection,dtype=np.int32)
        pool=self.pooled(selection)
        retained=(pool>=lower)&(pool<=(np.inf if upper is None else upper))&(pool>0)
        # Denominator is ALL valid windows before feature filtering.
        if equal_cds:
            weights=np.divide(sel,self.lengths,out=np.zeros(self.nr,dtype=float),where=self.lengths>0)
            g=sparse.csr_matrix((weights,(self.pop_index,np.arange(self.nr))),shape=(self.n,self.nr))
            counts=(g@self.record_counts).tocsr()
            denominator=np.bincount(self.pop_index,weights=(sel*(self.lengths>0)),minlength=self.n)
        else:
            counts=self.counts_by_population(sel)
            denominator=np.asarray(counts.sum(axis=1)).ravel()
        if np.any(denominator==0):raise ValueError('Population without eligible windows')
        freq=sparse.diags(1/denominator)@counts[:,retained]
        return freq.tocsr(),retained

def row_l1(x):
    d=np.sum(x,axis=1)
    if np.any(d<=0):raise ValueError('Zero-feature population; cannot normalize')
    return x/d[:,None]

def frequency_distance(freq):
    r=freq.copy();r.data=np.sqrt(r.data)
    gram=(r@r.T).toarray()
    norms=np.asarray(r.multiply(r).sum(axis=1)).ravel()
    squared=norms[:,None]+norms[None,:]-2*gram
    # Gram subtraction suffers cancellation for identical compositions. Re-evaluate
    # only near-zero pairs by sparse vector subtraction, not a biological cutoff.
    near=np.argwhere(np.triu(squared<1e-10,1))
    for i,j in near:
        delta=r.getrow(int(i))-r.getrow(int(j))
        exact=float(delta.multiply(delta).sum())
        squared[i,j]=squared[j,i]=exact
    d=np.sqrt(np.maximum(squared,0))
    np.fill_diagonal(d,0)
    return (d+d.T)/2

def analyze_distance(d,K,historical=None):
    d=np.asarray(d,dtype=float)
    if not np.all(np.isfinite(d)):raise ValueError('Non-finite distances')
    if len(d)<3 or not np.any(d>0):raise ValueError('Insufficient distinct profiles')
    if not isinstance(K,(int,np.integer)) or not 2<=K<len(d):
        raise ValueError('K must be an integer >=2 and < number of populations')
    z=linkage(squareform(d,checks=True),method='complete')
    labels=cut_tree(z,n_clusters=K).ravel()+1
    scores={}
    for kk in range(2,min(5,len(d)-1)+1):
        ll=cut_tree(z,n_clusters=kk).ravel()
        scores[kk]=float(silhouette_score(d,ll,metric='precomputed'))
    best=max(scores,key=lambda kk:(scores[kk],-kk))
    # The fixed-K request is independent of the limited diagnostic scan.
    requested_silhouette=scores[K] if K in scores else float(silhouette_score(d,labels,metric='precomputed'))
    # A strict scalar-height interval may not exist when merge heights tie.
    low=float(z[len(d)-K-1,2]) if len(d)>K else 0.0
    high=float(z[len(d)-K,2]) if K>1 else None
    info={'n':len(d),'K':K,'silhouette':requested_silhouette,'best_K_2to5':best,'best_silhouette':scores[best],
          'cut_height_lower_inclusive':low,'cut_height_upper_exclusive':high,'unique_height_cut_possible':bool(high is not None and high>low),
          'cut_rule':'scipy.cut_tree fixed number of clusters; distance interval reported separately',
          'ari_to_historical':float(adjusted_rand_score(historical,labels)) if historical is not None else None,
          'silhouette_scan':scores}
    return info,labels,z

def threshold_bounds(name,pool):
    if name=='all':return 1,None
    if name=='count2plus':return 2,None
    if name=='count5plus':return 5,None
    if name=='count10plus':return 10,None
    if name=='central_2_q99':return 2,max(2,int(np.ceil(np.quantile(pool[pool>0],.99))))
    raise ValueError(name)

def json_write(path,data):
    Path(path).write_text(json.dumps(data,ensure_ascii=False,indent=2,allow_nan=False),encoding='utf-8')

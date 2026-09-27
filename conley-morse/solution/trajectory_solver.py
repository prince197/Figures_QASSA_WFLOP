#!/usr/bin/env python3
import json, os, sys
from collections import defaultdict, Counter
from functools import lru_cache
from fractions import Fraction

sys.setrecursionlimit(200000)

class GF2Basis:
    def __init__(self):
        self.piv={}; self.rank=0
    def add(self, vec, label=0):
        v=int(vec); lab=int(label)
        while v:
            p=v.bit_length()-1
            if p in self.piv:
                vv,ll=self.piv[p]; v^=vv; lab^=ll
            else:
                self.piv[p]=(v,lab); self.rank+=1
                return True,(v,lab)
        return False,(0,lab)
    def reduce(self, vec, label=0):
        v=int(vec); lab=int(label)
        while v:
            p=v.bit_length()-1
            if p not in self.piv: break
            vv,ll=self.piv[p]; v^=vv; lab^=ll
        return v,lab

def independent_vectors(vecs):
    b=GF2Basis(); out=[]
    for v in vecs:
        ok,_=b.add(v)
        if ok: out.append(v)
    return out

def kernel_basis_from_columns(cols):
    b=GF2Basis(); deps=[]
    for j,col in enumerate(cols):
        v=int(col); coeff=1<<j
        while v:
            p=v.bit_length()-1
            if p in b.piv:
                vv,ll=b.piv[p]; v^=vv; coeff^=ll
            else:
                b.piv[p]=(v,coeff); b.rank+=1; break
        if v==0: deps.append(coeff)
    return deps



def _ceil_fraction(x):
    return -((-x.numerator)//x.denominator)

def _floor_fraction(x):
    return x.numerator//x.denominator

def _logistic_interval(num,scale,k,n):
    def endpoint(t):
        return Fraction(num*t*(n-t),scale*n*n)
    vals=[endpoint(k),endpoint(k+1)]
    if 2*k <= n <= 2*(k+1):
        vals.append(Fraction(num,4*scale))
    return min(vals),max(vals)

def _target_range(lo,hi,n):
    # Closed source and target boxes intersect at shared boundaries.
    a=max(0,_ceil_fraction(lo*n-1))
    b=min(n-1,_floor_fraction(hi*n))
    return a,b

def build_frames(data):
    nx,ny=data['phase_grid']['nx'],data['phase_grid']['ny']
    pg=data['parameter_grid']; model=data['model']
    scale=model['parameter_scale']; cnum=model['coupling_num']; cden=model['coupling_den']
    rnums=pg['r_num']; snums=pg['s_num']
    xr={(ri,i):_logistic_interval(r,scale,i,nx) for ri,r in enumerate(rnums) for i in range(nx)}
    yr={(sj,j):_logistic_interval(s,scale,j,ny) for sj,s in enumerate(snums) for j in range(ny)}
    frames=[]
    for sj,s in enumerate(snums):
        for ri,r in enumerate(rnums):
            offsets=[0]; targets=[]
            for j in range(ny):
                ly,hy=yr[(sj,j)]
                for i in range(nx):
                    lx,hx=xr[(ri,i)]
                    xlo=Fraction((cden-cnum),cden)*lx + Fraction(cnum,cden)*ly
                    xhi=Fraction((cden-cnum),cden)*hx + Fraction(cnum,cden)*hy
                    ylo=Fraction(cnum,cden)*lx + Fraction((cden-cnum),cden)*ly
                    yhi=Fraction(cnum,cden)*hx + Fraction((cden-cnum),cden)*hy
                    i0,i1=_target_range(xlo,xhi,nx); j0,j1=_target_range(ylo,yhi,ny)
                    for jj in range(j0,j1+1):
                        for ii in range(i0,i1+1):
                            targets.append(ii+nx*jj)
                    offsets.append(len(targets))
            frames.append({'offsets':offsets,'targets':targets})
    return frames

def tarjan_scc(n, offsets, targets):
    index=0; stack=[]; on=[False]*n; indices=[-1]*n; low=[0]*n; comps=[]
    def dfs(v):
        nonlocal index
        indices[v]=low[v]=index; index+=1; stack.append(v); on[v]=True
        for kk in range(offsets[v],offsets[v+1]):
            w=targets[kk]
            if indices[w]<0:
                dfs(w); low[v]=min(low[v],low[w])
            elif on[w]:
                low[v]=min(low[v],indices[w])
        if low[v]==indices[v]:
            c=[]
            while True:
                w=stack.pop(); on[w]=False; c.append(w)
                if w==v: break
            comps.append(c)
    for v in range(n):
        if indices[v]<0: dfs(v)
    return comps

def cells_of_boxes(boxes,nx):
    verts=set(); edges=set(); squares=set()
    for q in boxes:
        i,j=q%nx,q//nx
        squares.add((i,j))
        verts.update(((i,j),(i+1,j),(i,j+1),(i+1,j+1)))
        edges.update(((0,i,j),(0,i,j+1),(1,i,j),(1,i+1,j)))
    return verts,edges,squares

def homology_data_boxes(boxes,nx):
    KN=cells_of_boxes(boxes,nx)
    bases=[sorted(KN[d]) for d in range(3)]
    idx=[{c:i for i,c in enumerate(bases[d])} for d in range(3)]
    cols1=[]
    for axis,i,j in bases[1]:
        faces=((i,j),(i+1,j)) if axis==0 else ((i,j),(i,j+1))
        v=0
        for f in faces: v ^= 1<<idx[0][f]
        cols1.append(v)
    cols2=[]
    for i,j in bases[2]:
        faces=((0,i,j),(0,i,j+1),(1,i,j),(1,i+1,j))
        v=0
        for f in faces: v ^= 1<<idx[1][f]
        cols2.append(v)
    boundaries=[None,cols1,cols2]
    out=[]
    for d in range(3):
        Z=[1<<j for j in range(len(bases[0]))] if d==0 else kernel_basis_from_columns(boundaries[d])
        B=independent_vectors(boundaries[d+1] if d<2 else [])
        span=GF2Basis(); Bind=[]
        for v in B:
            ok,_=span.add(v)
            if ok: Bind.append(v)
        H=[]
        for z in Z:
            rem,_=span.reduce(z)
            if rem:
                span.add(z); H.append(z)
        coord=GF2Basis()
        for i,v in enumerate(Bind): coord.add(v,1<<i)
        for j,v in enumerate(H): coord.add(v,1<<(len(Bind)+j))
        out.append({'basis':bases[d],'B':Bind,'H':H,'coord':coord})
    return out

def relative_homology_data_boxes(P,E,nx):
    KP=cells_of_boxes(P,nx); KE=cells_of_boxes(E,nx)
    bases=[sorted(set(KP[d])-set(KE[d])) for d in range(3)]
    idx=[{c:i for i,c in enumerate(bases[d])} for d in range(3)]
    cols1=[]
    for axis,i,j in bases[1]:
        faces=((i,j),(i+1,j)) if axis==0 else ((i,j),(i,j+1))
        v=0
        for f in faces:
            q=idx[0].get(f)
            if q is not None: v ^= 1<<q
        cols1.append(v)
    cols2=[]
    for i,j in bases[2]:
        faces=((0,i,j),(0,i,j+1),(1,i,j),(1,i+1,j))
        v=0
        for f in faces:
            q=idx[1].get(f)
            if q is not None: v ^= 1<<q
        cols2.append(v)
    boundaries=[None,cols1,cols2]
    out=[]
    for d in range(3):
        Z=[1<<j for j in range(len(bases[0]))] if d==0 else kernel_basis_from_columns(boundaries[d])
        B=independent_vectors(boundaries[d+1] if d<2 else [])
        span=GF2Basis(); Bind=[]
        for v in B:
            ok,_=span.add(v)
            if ok: Bind.append(v)
        H=[]
        for z in Z:
            rem,_=span.reduce(z)
            if rem:
                span.add(z); H.append(z)
        coord=GF2Basis()
        for i,v in enumerate(Bind): coord.add(v,1<<i)
        for j,v in enumerate(H): coord.add(v,1<<(len(Bind)+j))
        out.append({'basis':bases[d],'B':Bind,'H':H,'coord':coord})
    return out

def relative_homology_betti_boxes(P,E,nx):
    KP=cells_of_boxes(P,nx); KE=cells_of_boxes(E,nx)
    bases=[sorted(set(KP[d])-set(KE[d])) for d in range(3)]
    idx=[{c:i for i,c in enumerate(bases[d])} for d in range(3)]
    cols1=[]
    for axis,i,j in bases[1]:
        faces=((i,j),(i+1,j)) if axis==0 else ((i,j),(i,j+1))
        v=0
        for f in faces:
            q=idx[0].get(f)
            if q is not None: v ^= 1<<q
        cols1.append(v)
    cols2=[]
    for i,j in bases[2]:
        faces=((0,i,j),(0,i,j+1),(1,i,j),(1,i+1,j))
        v=0
        for f in faces:
            q=idx[1].get(f)
            if q is not None: v ^= 1<<q
        cols2.append(v)
    r1=GF2Basis(); r2=GF2Basis()
    for v in cols1: r1.add(v)
    for v in cols2: r2.add(v)
    return [len(bases[0])-r1.rank, len(cols1)-r1.rank-r2.rank, len(cols2)-r2.rank]

@lru_cache(maxsize=4096)
def cached_homology_data_boxes(nx,boxes_tuple):
    return homology_data_boxes(set(boxes_tuple),nx)

def get_homology_data_boxes(boxes,nx):
    return cached_homology_data_boxes(nx,tuple(sorted(boxes)))

@lru_cache(maxsize=16384)
def cached_relative_homology_data_boxes(nx,P_tuple,E_tuple):
    return relative_homology_data_boxes(set(P_tuple),set(E_tuple),nx)

def get_relative_homology_data_boxes(P,E,nx):
    return cached_relative_homology_data_boxes(nx,tuple(sorted(P)),tuple(sorted(E)))

def apply_box_chain_map(vec,src_basis,tgt_index):
    out=0; x=vec
    while x:
        bit=x&-x; j=bit.bit_length()-1; x^=bit
        q=tgt_index.get(src_basis[j])
        if q is not None: out ^= 1<<q
    return out

def induced_box_map(srcH,tgtH,d):
    tgt_index={c:i for i,c in enumerate(tgtH[d]['basis'])}; bdim=len(tgtH[d]['B']); cols=[]
    for h in srcH[d]['H']:
        y=apply_box_chain_map(h,srcH[d]['basis'],tgt_index)
        rem,lab=tgtH[d]['coord'].reduce(y)
        if rem: raise RuntimeError('inclusion did not map cycles to represented cycles')
        cols.append(lab>>bdim)
    return cols

def map_rank(cols):
    b=GF2Basis()
    for c in cols: b.add(c)
    return b.rank

def common_image_rank(cols_a,cols_b):
    a=GF2Basis()
    for c in cols_a: a.add(c)
    ra=a.rank
    b=GF2Basis()
    for c in cols_b: b.add(c)
    rb=b.rank
    both=GF2Basis()
    for c in cols_a: both.add(c)
    for c in cols_b: both.add(c)
    return ra+rb-both.rank

def recurrent_structure(n,offsets,targets):
    comps=tarjan_scc(n,offsets,targets)
    recurrent=[]
    for c in comps:
        if len(c)>1: recurrent.append(sorted(c))
        else:
            v=c[0]
            if any(targets[k]==v for k in range(offsets[v],offsets[v+1])): recurrent.append([v])
    recurrent.sort(key=lambda c:c[0])
    cid={v:i for i,c in enumerate(comps) for v in c}
    cadj=[set() for _ in comps]
    for v in range(n):
        cv=cid[v]
        for kk in range(offsets[v],offsets[v+1]):
            cw=cid[targets[kk]]
            if cw!=cv: cadj[cv].add(cw)
    rcids=[cid[c[0]] for c in recurrent]
    desc=[set() for _ in recurrent]
    for i,rc in enumerate(rcids):
        seen={rc}; st=[rc]
        while st:
            x=st.pop()
            for y in cadj[x]:
                if y not in seen: seen.add(y); st.append(y)
        for j,rr in enumerate(rcids):
            if i!=j and rr in seen: desc[i].add(j)
    hasse=[]
    for i in range(len(recurrent)):
        for j in sorted(desc[i]):
            if not any(k!=i and k!=j and k in desc[i] and j in desc[k] for k in range(len(recurrent))):
                hasse.append([i,j])
    terminals=[i for i,d in enumerate(desc) if not d]
    return recurrent,desc,hasse,terminals,cid

def analyze_frame(frame,nx,ny):
    n=nx*ny; offsets=frame['offsets']; targets=frame['targets']
    recurrent,desc,edges,terminals,cid=recurrent_structure(n,offsets,targets)
    nodes=[]; pairs=[]
    for i,c in enumerate(recurrent):
        H=get_homology_data_boxes(set(c),nx)
        P=set(c); st=list(c)
        while st:
            v=st.pop()
            for kk in range(offsets[v],offsets[v+1]):
                w=targets[kk]
                if w not in P:
                    P.add(w); st.append(w)
        E=P-set(c)
        ib=relative_homology_betti_boxes(P,E,nx)
        nodes.append({'morse_id':i,'min_box_id':c[0],'size':len(c),
                      'block_betti_F2':[len(H[d]['H']) for d in range(3)],
                      'index_pair_betti_F2':ib,
                      'attractor':i in terminals})
        pairs.append({'P':P,'E':E})
    return {'recurrent':recurrent,'desc':desc,'nodes':nodes,'pairs':pairs,'edges':edges,'terminals':terminals}

def node_label(node):
    return (tuple(node['block_betti_F2']),tuple(node['index_pair_betti_F2']),bool(node['attractor']))

def refined_signatures(analyses):
    labels=[[node_label(n) for n in a['nodes']] for a in analyses]
    palette={x:i for i,x in enumerate(sorted({x for row in labels for x in row}))}
    colors=[[palette[x] for x in row] for row in labels]
    while True:
        keys=[]; node_keys=[]
        for a,crow in zip(analyses,colors):
            pred=[[] for _ in crow]; succ=[[] for _ in crow]
            for u,v in a['edges']:
                succ[u].append(crow[v]); pred[v].append(crow[u])
            row=[]
            for i,c in enumerate(crow):
                k=(c,tuple(sorted(pred[i])),tuple(sorted(succ[i])))
                row.append(k); keys.append(k)
            node_keys.append(row)
        pal={x:i for i,x in enumerate(sorted(set(keys)))}
        newcolors=[[pal[k] for k in row] for row in node_keys]
        if len(pal)==len({c for row in colors for c in row}):
            colors=newcolors; break
        colors=newcolors
    sigs=[]
    for a,crow in zip(analyses,colors):
        nc=Counter(crow); ec=Counter((crow[u],crow[v]) for u,v in a['edges'])
        sigs.append((tuple(sorted(nc.items())),tuple(sorted((u,v,n) for (u,v),n in ec.items()))))
    return sigs

def union_edge_structure(frame_a,ana_a,frame_b,ana_b,n):
    oa,ta=frame_a['offsets'],frame_a['targets']; ob,tb=frame_b['offsets'],frame_b['targets']
    offsets=[0]; targets=[]
    for q in range(n):
        s=set(ta[oa[q]:oa[q+1]]); s.update(tb[ob[q]:ob[q+1]])
        targets.extend(sorted(s)); offsets.append(len(targets))
    recurrent,desc,hasse,terminals,cid=recurrent_structure(n,offsets,targets)
    urid={cid[c[0]]:i for i,c in enumerate(recurrent)}
    amap=[]
    for A in ana_a['recurrent']:
        cc=cid[A[0]]; amap.append(urid[cc])
        if any(cid[v]!=cc for v in A): raise RuntimeError('source SCC split under union')
    bmap=[]
    for B in ana_b['recurrent']:
        cc=cid[B[0]]; bmap.append(urid[cc])
        if any(cid[v]!=cc for v in B): raise RuntimeError('target SCC split under union')
    rel=sorted([[i,j] for i,u in enumerate(amap) for j,v in enumerate(bmap) if u==v])
    return rel,{'recurrent':recurrent,'desc':desc,'source_map':amap,'target_map':bmap}

def aggregate_transport_for_edge(ana_a,ana_b,rel,nx):
    src=sorted({i for i,j in rel}); tgt=sorted({j for i,j in rel})
    def make_pair(ana,ids):
        B=set(); P=set()
        for mid in ids:
            B.update(ana['recurrent'][mid]); P.update(ana['pairs'][mid]['P'])
        return P,P-B
    Pa,Ea=make_pair(ana_a,src); Pb,Eb=make_pair(ana_b,tgt)
    Pu=Pa|Pb; Eu=Ea|Eb
    Ha=get_relative_homology_data_boxes(Pa,Ea,nx)
    Hb=get_relative_homology_data_boxes(Pb,Eb,nx)
    Hu=get_relative_homology_data_boxes(Pu,Eu,nx)
    ar=[]; br=[]; cr=[]
    for d in range(3):
        fa=induced_box_map(Ha,Hu,d); fb=induced_box_map(Hb,Hu,d)
        ar.append(map_rank(fa)); br.append(map_rank(fb)); cr.append(common_image_rank(fa,fb))
    return {'source_morse_ids':src,'target_morse_ids':tgt,
            'source_betti_F2':[len(Ha[d]['H']) for d in range(3)],
            'target_betti_F2':[len(Hb[d]['H']) for d in range(3)],
            'bridge_betti_F2':[len(Hu[d]['H']) for d in range(3)],
            'source_to_bridge_rank_F2':ar,'target_to_bridge_rank_F2':br,
            'common_image_rank_F2':cr}

def classify_event(ana_a,ana_b,rel):
    da=[0]*len(ana_a['nodes']); db=[0]*len(ana_b['nodes'])
    for i,j in rel: da[i]+=1; db[j]+=1
    birth=any(x==0 for x in db); death=any(x==0 for x in da)
    split=any(x>1 for x in da); merge=any(x>1 for x in db)
    if split and merge: return 'split-merge'
    if split: return 'split'
    if merge: return 'merge'
    if birth and death: return 'birth-death'
    if birth: return 'birth'
    if death: return 'death'
    mapping={i:j for i,j in rel}
    same_nodes=len(rel)==len(ana_a['nodes'])==len(ana_b['nodes']) and all(node_label(ana_a['nodes'][i])==node_label(ana_b['nodes'][j]) for i,j in rel)
    mapped_edges=sorted((mapping[a],mapping[b]) for a,b in ana_a['edges']) if len(mapping)==len(ana_a['nodes']) else []
    return 'continuation' if same_nodes and mapped_edges==sorted(map(tuple,ana_b['edges'])) else 'rewiring'

def relation_for_oriented_edge(p,q,edge_rel):
    if (p,q) in edge_rel: return edge_rel[(p,q)]
    return [[j,i] for i,j in edge_rel[(q,p)]]

def relation_space(dims,maps,a,b):
    offs=[0]
    for i in range(a,b+1): offs.append(offs[-1]+dims[i])
    rels=[]
    for k in range(a,b):
        direction,cols=maps[k]
        src,tgt=(k,k+1) if direction=='forward' else (k+1,k)
        so,to=offs[src-a],offs[tgt-a]
        for j,col in enumerate(cols):
            v=1<<(so+j); x=col
            while x:
                bit=x&-x; r=bit.bit_length()-1; x^=bit; v^=1<<(to+r)
            rels.append(v)
    return offs,independent_vectors(rels)

def limit_space(dims,maps,a,b):
    offs=[0]
    for i in range(a,b+1): offs.append(offs[-1]+dims[i])
    cbase=[]; ncon=0
    for k in range(a,b):
        direction,_=maps[k]; td=dims[k+1] if direction=='forward' else dims[k]
        cbase.append(ncon); ncon+=td
    cols=[0]*offs[-1]
    for kk,k in enumerate(range(a,b)):
        direction,mcols=maps[k]; co=cbase[kk]
        src,tgt=(k,k+1) if direction=='forward' else (k+1,k)
        so,to=offs[src-a],offs[tgt-a]
        for j,col in enumerate(mcols):
            x=col
            while x:
                bit=x&-x; r=bit.bit_length()-1; x^=bit; cols[so+j]^=1<<(co+r)
        for r in range(dims[tgt]): cols[to+r]^=1<<(co+r)
    return offs,kernel_basis_from_columns(cols)

def interval_rank(dims,maps,a,b):
    if a==b: return dims[a]
    offs,rels=relation_space(dims,maps,a,b); _,lim=limit_space(dims,maps,a,b)
    basis=GF2Basis()
    for v in rels: basis.add(v)
    r0=basis.rank; mask=(1<<dims[a])-1 if dims[a] else 0
    for x in lim: basis.add(x&mask)
    return basis.rank-r0

def _rank_bit_columns(cols):
    basis=GF2Basis()
    for x in cols:
        basis.add(x)
    return basis.rank

def _edge_linear_relation(du,dv,direction,cols):
    """Basis for the edge relation inside F2^du direct_sum F2^dv."""
    out=[]
    if direction=='forward':
        for i,col in enumerate(cols):
            out.append((1<<i)|(col<<du))
    else:
        for j,col in enumerate(cols):
            out.append(col|((1<<j)<<du))
    return independent_vectors(out)

def _compose_linear_relations(R,du,dv,S,dw):
    """Compose R subset U x V with S subset V x W."""
    nr=len(R)
    vm=(1<<dv)-1 if dv else 0
    constraints=[(x>>du)&vm for x in R]
    constraints.extend(x&vm for x in S)
    coeff_kernel=kernel_basis_from_columns(constraints)
    out=[]
    um=(1<<du)-1 if du else 0
    for coeff in coeff_kernel:
        u=w=0
        x=coeff
        while x:
            bit=x&-x
            j=bit.bit_length()-1
            x^=bit
            if j<nr:
                u^=R[j]&um
            else:
                w^=S[j-nr]>>dv
        out.append(u|(w<<du))
    return independent_vectors(out)

def _relation_regular_rank(R,du,dv):
    if not R:
        return 0
    um=(1<<du)-1 if du else 0
    left_rank=_rank_bit_columns([x&um for x in R])
    right_rank=_rank_bit_columns([x>>du for x in R])
    return left_rank+right_rank-len(R)

def barcode_for_zigzag(dims,maps):
    """Compute every generalized interval rank by incremental linear relations.

    For a path representation, the generalized rank on [a,b] is the regular
    rank of the composed endpoint relation.  This avoids rebuilding a full
    limit and colimit for every interval while producing the same rank table.
    """
    n=len(dims)
    rho=[[0]*(n-a) for a in range(n)]
    edge_rel=[_edge_linear_relation(dims[i],dims[i+1],maps[i][0],maps[i][1]) for i in range(n-1)]
    for a in range(n):
        rho[a][0]=dims[a]
        if dims[a]==0:
            continue
        # Identity relation from V_a to itself.  Generalized interval rank is
        # monotone under interval extension, so after it reaches zero the rest
        # of this row is forced to zero.
        R=[(1<<i)|((1<<i)<<dims[a]) for i in range(dims[a])]
        right_dim=dims[a]
        for b in range(a+1,n):
            R=_compose_linear_relations(R,dims[a],right_dim,edge_rel[b-1],dims[b])
            right_dim=dims[b]
            rr=_relation_regular_rank(R,dims[a],dims[b])
            rho[a][b-a]=rr
            if rr==0:
                break
    def RANK(a,b):
        if a<0 or b>=n or a>b:
            return 0
        return rho[a][b-a]
    bars=[]
    for a in range(n):
        for b in range(a,n):
            m=RANK(a,b)-RANK(a-1,b)-RANK(a,b+1)+RANK(a-1,b+1)
            if m>0:
                bars.append({'start':a,'end':b,'multiplicity':m})
            elif m<0:
                raise RuntimeError('negative zigzag multiplicity')
    return rho,bars

def selection_pair(pid,sel,frames,analyses):
    if not sel: return set(),set(),set()
    B=set(); P=set()
    for mid in sel:
        B.update(analyses[pid]['recurrent'][mid]); P.update(analyses[pid]['pairs'][mid]['P'])
    E=P-B
    return B,P,E

def path_analysis(pathrec,frames,analyses,edge_rel,nx):
    ids=pathrec['parameter_ids']
    seed=pathrec['seed_box_id']
    start=next(i for i,c in enumerate(analyses[ids[0]]['recurrent']) if seed in c)
    selected=[[start]]
    for p,q in zip(ids,ids[1:]):
        rel=relation_for_oriented_edge(p,q,edge_rel)
        prev=set(selected[-1])
        selected.append(sorted({j for i,j in rel if i in prev}))

    # Construct the zigzag maps sequentially so long trajectories do not retain
    # hundreds of full cubical chain complexes at once.
    dims_by_d=[[],[],[]]
    maps_by_d=[[],[],[]]
    Bcur,Pcur,Ecur=selection_pair(ids[0],selected[0],frames,analyses)
    Hcur=get_relative_homology_data_boxes(Pcur,Ecur,nx)
    for d in range(3):
        dims_by_d[d].append(len(Hcur[d]['H']))
    for k in range(len(ids)-1):
        Bnext,Pnext,Enext=selection_pair(ids[k+1],selected[k+1],frames,analyses)
        Hnext=get_relative_homology_data_boxes(Pnext,Enext,nx)
        Pb=Pcur|Pnext
        Eb=Ecur|Enext
        Hbridge=get_relative_homology_data_boxes(Pb,Eb,nx)
        for d in range(3):
            dims_by_d[d].append(len(Hbridge[d]['H']))
            dims_by_d[d].append(len(Hnext[d]['H']))
            maps_by_d[d].append(('forward',induced_box_map(Hcur,Hbridge,d)))
            maps_by_d[d].append(('backward',induced_box_map(Hnext,Hbridge,d)))
        Pcur,Ecur,Hcur=Pnext,Enext,Hnext

    dimsout=[]
    for d in range(3):
        rho,bars=barcode_for_zigzag(dims_by_d[d],maps_by_d[d])
        dimsout.append({'dimension':d,'node_dimensions_F2':dims_by_d[d],
                        'rank_rows_F2':rho,'interval_multiplicities_F2':bars})
    return {'path_id':pathrec['path_id'],'sweep_id':pathrec['sweep_id'],'parameter_ids':ids,'seed_box_id':seed,
            'selected_morse_ids_by_frame':selected,'dimensions':dimsout}




def _relation_from_span(left_cols,right_cols):
    """Kernel relation {(x,y): f(x)=g(y)} in left direct_sum right."""
    return independent_vectors(kernel_basis_from_columns(list(left_cols)+list(right_cols)))

def _subspace_intersection_basis(A,B):
    if not A or not B:
        return []
    na=len(A)
    deps=kernel_basis_from_columns(list(A)+list(B))
    out=[]
    for dep in deps:
        v=0; x=dep
        while x:
            bit=x&-x; j=bit.bit_length()-1; x^=bit
            if j<na:
                v^=A[j]
        if v:
            out.append(v)
    return independent_vectors(out)

def _relation_projection_ranks(R,du):
    if not R:
        return 0,0
    mask=(1<<du)-1 if du else 0
    return _rank_bit_columns([x&mask for x in R]), _rank_bit_columns([x>>du for x in R])

def _quiver_rank_invariant(node_dims, arrows):
    """Return (limit dim, colimit dim, rank(lim->colim)) for a finite quiver.

    node_dims is a list.  Each arrow is (source_node,target_node,columns), where
    columns are target-coordinate bit vectors for the source basis.
    """
    offs=[0]
    for d in node_dims:
        offs.append(offs[-1]+d)
    total=offs[-1]

    # Inverse-limit constraints f(x_source)+x_target=0, one target block per arrow.
    coff=[]; ncon=0
    for _,t,_ in arrows:
        coff.append(ncon); ncon+=node_dims[t]
    ccols=[0]*total
    for k,(src,tgt,cols) in enumerate(arrows):
        co=coff[k]; so=offs[src]; to=offs[tgt]
        if len(cols)!=node_dims[src]:
            raise RuntimeError('quiver map column count mismatch')
        for j,col in enumerate(cols):
            x=col
            while x:
                bit=x&-x; r=bit.bit_length()-1; x^=bit
                ccols[so+j]^=1<<(co+r)
        for r in range(node_dims[tgt]):
            ccols[to+r]^=1<<(co+r)
    lim=kernel_basis_from_columns(ccols)

    # Colimit relations identify x_source with f(x_source) at the target.
    rel=[]
    for src,tgt,cols in arrows:
        so=offs[src]; to=offs[tgt]
        for j,col in enumerate(cols):
            v=1<<(so+j); x=col
            while x:
                bit=x&-x; r=bit.bit_length()-1; x^=bit
                v^=1<<(to+r)
            rel.append(v)
    rel=independent_vectors(rel)
    span=GF2Basis()
    for v in rel: span.add(v)
    r0=span.rank
    for v in lim: span.add(v)
    return len(lim), total-r0, span.rank-r0

def _patch_vertices_edges(rec,na):
    i0,j0,i1,j1=rec['rectangle']
    verts=[i+na*j for j in range(j0,j1+1) for i in range(i0,i1+1)]
    vset=set(verts); edges=[]
    for p in verts:
        i=p%na; j=p//na
        if i<i1 and p+1 in vset: edges.append((p,p+1))
        if j<j1 and p+na in vset: edges.append((p,p+na))
    return verts,edges

def parameter_patch_analysis(rec,frames,analyses,edge_rel,nx,na):
    verts,edges=_patch_vertices_edges(rec,na)
    if rec['seed_parameter_id'] not in set(verts):
        raise RuntimeError('patch seed parameter outside rectangle')
    seedp=rec['seed_parameter_id']; seedbox=rec['seed_box_id']
    try:
        seedmid=next(i for i,c in enumerate(analyses[seedp]['recurrent']) if seedbox in c)
    except StopIteration:
        raise RuntimeError('patch seed box is not recurrent')

    selected={p:set() for p in verts}; selected[seedp].add(seedmid)
    changed=True
    while changed:
        changed=False
        for p,q in edges:
            rel=relation_for_oriented_edge(p,q,edge_rel)
            addq={j for i,j in rel if i in selected[p]}
            addp={i for i,j in rel if j in selected[q]}
            if not addq.issubset(selected[q]): selected[q].update(addq); changed=True
            if not addp.issubset(selected[p]): selected[p].update(addp); changed=True

    pair={}; H={}
    for p in verts:
        _,P,E=selection_pair(p,sorted(selected[p]),frames,analyses)
        pair[p]=(P,E)
        H[p]=get_relative_homology_data_boxes(P,E,nx)

    edgeH={}
    for p,q in edges:
        Pp,Ep=pair[p]; Pq,Eq=pair[q]
        edgeH[(p,q)]=get_relative_homology_data_boxes(Pp|Pq,Ep|Eq,nx)

    degree_records=[]
    for d in range(3):
        # Endpoint nodes first, then one bridge node per edge in the stated edge order.
        node_dims=[len(H[p][d]['H']) for p in verts]
        node_dims.extend(len(edgeH[e][d]['H']) for e in edges)
        index={p:i for i,p in enumerate(verts)}
        arrows=[]
        for ei,(p,q) in enumerate(edges):
            bnode=len(verts)+ei; Hb=edgeH[(p,q)]
            arrows.append((index[p],bnode,induced_box_map(H[p],Hb,d)))
            arrows.append((index[q],bnode,induced_box_map(H[q],Hb,d)))
        ld,cd,gr=_quiver_rank_invariant(node_dims,arrows)
        degree_records.append({'dimension':d,'limit_dimension_F2':ld,'colimit_dimension_F2':cd,'generalized_rank_F2':gr})

    return {'patch_id':rec['patch_id'],'rectangle':rec['rectangle'],'seed_parameter_id':seedp,'seed_box_id':seedbox,
            'parameter_ids':verts,'selected_morse_ids_by_parameter':[sorted(selected[p]) for p in verts],
            'degrees':degree_records}

def continuation_square_coherence(analyses,edge_struct,nx,na,nb):
    # Use the aggregate endpoint pairs participating in each union SCC relation.
    endpoint_pair={}; endpoint_H={}
    for p,ana in enumerate(analyses):
        B=set(); P=set()
        for mid in range(len(ana['recurrent'])):
            B.update(ana['recurrent'][mid]); P.update(ana['pairs'][mid]['P'])
        E=P-B; endpoint_pair[p]=(P,E)
        endpoint_H[p]=get_relative_homology_data_boxes(P,E,nx)

    edge_R={}
    for (p,q),_U in edge_struct.items():
        Pp,Ep=endpoint_pair[p]; Pq,Eq=endpoint_pair[q]
        Hb=get_relative_homology_data_boxes(Pp|Pq,Ep|Eq,nx)
        rec=[]
        for d in range(3):
            lp=induced_box_map(endpoint_H[p],Hb,d)
            rq=induced_box_map(endpoint_H[q],Hb,d)
            rec.append(_relation_from_span(lp,rq))
        edge_R[(p,q)]=rec

    out=[]
    for j in range(nb-1):
        for i in range(na-1):
            p=i+na*j; q=p+1; r=p+na; s=r+1
            deg=[]
            for d in range(3):
                du=len(endpoint_H[p][d]['H']); dq=len(endpoint_H[q][d]['H']); dr=len(endpoint_H[r][d]['H']); dw=len(endpoint_H[s][d]['H'])
                Rhv=_compose_linear_relations(edge_R[(p,q)][d],du,dq,edge_R[(q,s)][d],dw)
                Rvh=_compose_linear_relations(edge_R[(p,r)][d],du,dr,edge_R[(r,s)][d],dw)
                I=_subspace_intersection_basis(Rhv,Rvh)
                lh,rh=_relation_projection_ranks(I,du)
                deg.append({'dimension':d,
                            'horizontal_then_vertical_relation_dimension_F2':len(Rhv),
                            'vertical_then_horizontal_relation_dimension_F2':len(Rvh),
                            'common_relation_dimension_F2':len(I),
                            'horizontal_then_vertical_regular_rank_F2':_relation_regular_rank(Rhv,du,dw),
                            'vertical_then_horizontal_regular_rank_F2':_relation_regular_rank(Rvh,du,dw),
                            'common_regular_rank_F2':_relation_regular_rank(I,du,dw),
                            'common_left_projection_rank_F2':lh,
                            'common_right_projection_rank_F2':rh})
            out.append({'lower_left_parameter_id':p,'lower_right_parameter_id':q,
                        'upper_left_parameter_id':r,'upper_right_parameter_id':s,'degrees':deg})
    return out

def sink_first_morse_order(ana):
    m=len(ana['nodes']); chosen=[]; used=set()
    while len(chosen)<m:
        ready=[i for i in range(m) if i not in used and ana['desc'][i].issubset(used)]
        if not ready: raise RuntimeError('Morse order is not acyclic')
        x=min(ready); chosen.append(x); used.add(x)
    return chosen

def attractor_filtration_sets(ana):
    order=sink_first_morse_order(ana); out=[]; cur=set()
    for mid in order:
        cur.update(ana['pairs'][mid]['P'])
        out.append(set(cur))
    return order,out

def cell_boundary_F2(cell,d):
    if d==1:
        axis,i,j=cell
        return ((i,j),(i+1,j)) if axis==0 else ((i,j),(i,j+1))
    if d==2:
        i,j=cell
        return ((0,i,j),(0,i,j+1),(1,i,j),(1,i+1,j))
    return ()

def connecting_map_rank(srcH,tgtH,d):
    if d<=0: return 0
    tgt_index={c:i for i,c in enumerate(tgtH[d-1]['basis'])}; bdim=len(tgtH[d-1]['B']); cols=[]
    for h in srcH[d]['H']:
        y=0; x=h
        while x:
            bit=x&-x; j=bit.bit_length()-1; x^=bit
            for face in cell_boundary_F2(srcH[d]['basis'][j],d):
                q=tgt_index.get(face)
                if q is not None: y ^= 1<<q
        rem,lab=tgtH[d-1]['coord'].reduce(y)
        if rem: raise RuntimeError('connecting boundary is not represented in target homology')
        cols.append(lab>>bdim)
    return map_rank(cols)

def connection_braid_analysis(pid,ana,nx,pairs):
    order,A=attractor_filtration_sets(ana); graded=[]
    for k,X in enumerate(A):
        Y=set() if k==0 else A[k-1]
        H=relative_homology_data_boxes(X,Y,nx)
        graded.append([len(H[d]['H']) for d in range(3)])
    con=[]
    target_cache={}
    for a,b in pairs:
        if not (0<=a<b<len(A)): raise RuntimeError('invalid braid query pair')
        if a not in target_cache:
            Z=set() if a==0 else A[a-1]; Y=A[a]
            target_cache[a]=relative_homology_data_boxes(Y,Z,nx)
        tgt=target_cache[a]; Y=A[a]
        src=relative_homology_data_boxes(A[b],Y,nx)
        sb=[len(src[d]['H']) for d in range(3)]; tb=[len(tgt[d]['H']) for d in range(3)]
        ranks=[0,connecting_map_rank(src,tgt,1),connecting_map_rank(src,tgt,2)]
        con.append({'a':a,'b':b,'source_betti_F2':sb,'target_betti_F2':tb,'connecting_rank_F2':ranks})
    return {'parameter_id':pid,'morse_order':order,'filtration_box_counts':[len(x) for x in A],
            'graded_betti_F2':graded,'connecting_maps':con}

def graph_homology_data(nverts,edges):
    edges=sorted({tuple(sorted(e)) for e in edges if e[0]!=e[1]})
    b0=list(range(nverts)); b1=edges
    cols=[]
    for u,v in b1: cols.append((1<<u)^(1<<v))
    out=[]
    # H0
    B0=independent_vectors(cols); span=GF2Basis(); Bind=[]
    for v in B0:
        ok,_=span.add(v)
        if ok: Bind.append(v)
    H0=[]
    for j in range(nverts):
        z=1<<j; rem,_=span.reduce(z)
        if rem: span.add(z); H0.append(z)
    coord0=GF2Basis()
    for i,v in enumerate(Bind): coord0.add(v,1<<i)
    for j,v in enumerate(H0): coord0.add(v,1<<(len(Bind)+j))
    out.append({'basis':b0,'B':Bind,'H':H0,'coord':coord0})
    # H1, graph has no 2-cells
    H1=kernel_basis_from_columns(cols); coord1=GF2Basis()
    for j,v in enumerate(H1): coord1.add(v,1<<j)
    out.append({'basis':b1,'B':[],'H':H1,'coord':coord1})
    return out

def induced_graph_map(srcH,tgtH,d,vertex_map):
    if d==0:
        tgt_index={v:i for i,v in enumerate(tgtH[0]['basis'])}; bdim=len(tgtH[0]['B']); cols=[]
        for h in srcH[0]['H']:
            y=0; x=h
            while x:
                bit=x&-x; j=bit.bit_length()-1; x^=bit; y^=1<<tgt_index[vertex_map[j]]
            rem,lab=tgtH[0]['coord'].reduce(y)
            if rem: raise RuntimeError('graph H0 map failure')
            cols.append(lab>>bdim)
        return cols
    tgt_index={e:i for i,e in enumerate(tgtH[1]['basis'])}; cols=[]
    for h in srcH[1]['H']:
        y=0; x=h
        while x:
            bit=x&-x; j=bit.bit_length()-1; x^=bit
            u,v=srcH[1]['basis'][j]; a,b=vertex_map[u],vertex_map[v]
            if a!=b:
                e=tuple(sorted((a,b)))
                if e not in tgt_index: raise RuntimeError('graph map does not preserve edges')
                y^=1<<tgt_index[e]
        rem,lab=tgtH[1]['coord'].reduce(y)
        if rem: raise RuntimeError('graph H1 map failure')
        cols.append(lab)
    return cols

def order_graph_edges(desc):
    return [(i,j) for i,d in enumerate(desc) for j in sorted(d) if i!=j]

def graph_path_analysis(sweep,analyses,edge_struct):
    ids=sweep['parameter_ids']; Hnodes=[]; maps_by_d=[[],[]]
    endpointH={}
    for p in ids:
        endpointH[p]=graph_homology_data(len(analyses[p]['recurrent']),order_graph_edges(analyses[p]['desc']))
    unionHs=[]
    unions=[]
    for p,q in zip(ids,ids[1:]):
        key=(p,q) if (p,q) in edge_struct else (q,p)
        U=edge_struct[key]
        unionHs.append(graph_homology_data(len(U['recurrent']),order_graph_edges(U['desc'])))
        unions.append((key,U))
    for i,p in enumerate(ids):
        Hnodes.append(endpointH[p])
        if i+1<len(ids): Hnodes.append(unionHs[i])
    for d in range(2):
        maps=[]
        for i,(p,q) in enumerate(zip(ids,ids[1:])):
            key,U=unions[i]
            if key==(p,q): pmap,qmap=U['source_map'],U['target_map']
            else: pmap,qmap=U['target_map'],U['source_map']
            maps.append(('forward',induced_graph_map(endpointH[p],unionHs[i],d,pmap)))
            maps.append(('backward',induced_graph_map(endpointH[q],unionHs[i],d,qmap)))
        maps_by_d[d]=maps
    dimsout=[]
    for d in range(2):
        dims=[len(H[d]['H']) for H in Hnodes]; rho,bars=barcode_for_zigzag(dims,maps_by_d[d])
        dimsout.append({'dimension':d,'node_dimensions_F2':dims,'rank_rows_F2':rho,'interval_multiplicities_F2':bars})
    return {'sweep_id':sweep['sweep_id'],'parameter_ids':ids,'dimensions':dimsout}



def morse_graph_probe_analysis(rec,analyses,edge_struct):
    """Zigzag persistence of H0,H1 for Morse reachability graphs."""
    tmp={'sweep_id':rec['probe_id'],'parameter_ids':rec['parameter_ids']}
    ans=graph_path_analysis(tmp,analyses,edge_struct)
    return {'probe_id':rec['probe_id'],'parameter_ids':rec['parameter_ids'],'dimensions':ans['dimensions']}


def morse_graph_square_coherence(analyses,edge_struct,na,nb):
    """Two-route coherence on H0,H1 of the Morse reachability graphs."""
    endpoint_H={}
    for p,a in enumerate(analyses):
        endpoint_H[p]=graph_homology_data(len(a['recurrent']),order_graph_edges(a['desc']))
    edge_R={}
    for (p,q),U in edge_struct.items():
        Hu=graph_homology_data(len(U['recurrent']),order_graph_edges(U['desc']))
        rec=[]
        for d in range(2):
            fp=induced_graph_map(endpoint_H[p],Hu,d,U['source_map'])
            fq=induced_graph_map(endpoint_H[q],Hu,d,U['target_map'])
            rec.append(_relation_from_span(fp,fq))
        edge_R[(p,q)]=rec
    out=[]
    for j in range(nb-1):
        for i in range(na-1):
            p=i+na*j; q=p+1; r=p+na; s=r+1
            deg=[]
            for d in range(2):
                du=len(endpoint_H[p][d]['H']); dq=len(endpoint_H[q][d]['H'])
                dr=len(endpoint_H[r][d]['H']); dw=len(endpoint_H[s][d]['H'])
                Rhv=_compose_linear_relations(edge_R[(p,q)][d],du,dq,edge_R[(q,s)][d],dw)
                Rvh=_compose_linear_relations(edge_R[(p,r)][d],du,dr,edge_R[(r,s)][d],dw)
                I=_subspace_intersection_basis(Rhv,Rvh)
                lh,rh=_relation_projection_ranks(I,du)
                deg.append({'dimension':d,
                            'horizontal_then_vertical_relation_dimension_F2':len(Rhv),
                            'vertical_then_horizontal_relation_dimension_F2':len(Rvh),
                            'common_relation_dimension_F2':len(I),
                            'horizontal_then_vertical_regular_rank_F2':_relation_regular_rank(Rhv,du,dw),
                            'vertical_then_horizontal_regular_rank_F2':_relation_regular_rank(Rvh,du,dw),
                            'common_regular_rank_F2':_relation_regular_rank(I,du,dw),
                            'common_left_projection_rank_F2':lh,
                            'common_right_projection_rank_F2':rh})
            out.append({'lower_left_parameter_id':p,'lower_right_parameter_id':q,
                        'upper_left_parameter_id':r,'upper_right_parameter_id':s,'degrees':deg})
    return out

def solve(data):
    nx,ny=data['phase_grid']['nx'],data['phase_grid']['ny']; frames=build_frames(data); n=nx*ny
    analyses=[analyze_frame(f,nx,ny) for f in frames]
    sigtuples=refined_signatures(analyses); unique=sorted(set(sigtuples)); sig_id={s:i for i,s in enumerate(unique)}
    signatures=[]
    for sid,s in enumerate(unique):
        nodes,edges=s
        signatures.append({'signature_id':sid,
                           'node_color_counts':[[c,n] for c,n in nodes],
                           'edge_color_counts':[[u,v,n] for u,v,n in edges]})
    frameout=[]
    for pid,a in enumerate(analyses):
        frameout.append({'parameter_id':pid,'morse_sets':a['nodes'],'morse_edges':a['edges'],
                         'terminal_morse_ids':a['terminals'],'signature_id':sig_id[sigtuples[pid]]})
    edge_rel={}; edge_struct={}; edge_events=[]
    for p,q in data['parameter_edges']:
        rel,U=union_edge_structure(frames[p],analyses[p],frames[q],analyses[q],n)
        edge_rel[(p,q)]=rel; edge_struct[(p,q)]=U
        edge_events.append({'source_parameter_id':p,'target_parameter_id':q,'event_type':classify_event(analyses[p],analyses[q],rel),
                            'continuation_relation':rel})
    transport_queries=[]
    for qid,eidx in enumerate(data.get('transport_edge_indices',[])):
        p,q=data['parameter_edges'][eidx]
        rec=aggregate_transport_for_edge(analyses[p],analyses[q],edge_rel[(p,q)],nx)
        rec={'query_id':qid,'parameter_edge_index':eidx,'source_parameter_id':p,'target_parameter_id':q,**rec}
        transport_queries.append(rec)
    padj=[[] for _ in frames]
    for p,q in data['parameter_edges']: padj[p].append(q); padj[q].append(p)
    seen=set(); regions=[]; param_region=[-1]*len(frames)
    for p in range(len(frames)):
        if p in seen: continue
        sid=sig_id[sigtuples[p]]; comp=[]; st=[p]; seen.add(p)
        while st:
            x=st.pop(); comp.append(x)
            for y in padj[x]:
                if y not in seen and sig_id[sigtuples[y]]==sid:
                    seen.add(y); st.append(y)
        comp.sort(); rid=len(regions)
        for x in comp: param_region[x]=rid
        regions.append({'region_id':rid,'signature_id':sid,'representative_parameter_id':comp[0],'parameter_ids':comp})
    event_names=['continuation','birth','death','split','merge','split-merge','birth-death','rewiring']
    agg={}
    for idx,ev in enumerate(edge_events):
        p,q=ev['source_parameter_id'],ev['target_parameter_id']; a,b=param_region[p],param_region[q]
        if a==b: continue
        u,v=sorted((a,b)); rec=agg.setdefault((u,v),{'region_u':u,'region_v':v,'parameter_edge_count':0,'event_type_counts':Counter(),'witness_parameter_edge':None,'witness_index':10**9})
        rec['parameter_edge_count']+=1; rec['event_type_counts'][ev['event_type']]+=1
        if idx<rec['witness_index']:
            rec['witness_index']=idx; rec['witness_parameter_edge']=[p,q]
    regadj=[]
    for key in sorted(agg):
        r=agg[key]
        regadj.append({'region_u':r['region_u'],'region_v':r['region_v'],'parameter_edge_count':r['parameter_edge_count'],
                       'event_type_counts':{k:r['event_type_counts'].get(k,0) for k in event_names},
                       'witness_parameter_edge':r['witness_parameter_edge']})
    conley_paths=[]
    for pr in data['paths']:
        conley_paths.append(path_analysis(pr,frames,analyses,edge_rel,nx))
        cached_relative_homology_data_boxes.cache_clear()
        cached_homology_data_boxes.cache_clear()
    morse_graph_paths=[morse_graph_probe_analysis(rec,analyses,edge_struct) for rec in data.get('morse_graph_probes',[])]
    morse_graph_coherence=morse_graph_square_coherence(analyses,edge_struct,data['parameter_grid']['na'],data['parameter_grid']['nb'])
    parameter_patches=[parameter_patch_analysis(rec,frames,analyses,edge_rel,nx,data['parameter_grid']['na']) for rec in data.get('patches',[])]
    cached_relative_homology_data_boxes.cache_clear(); cached_homology_data_boxes.cache_clear()
    continuation_coherence=continuation_square_coherence(analyses,edge_struct,nx,data['parameter_grid']['na'],data['parameter_grid']['nb'])
    cached_relative_homology_data_boxes.cache_clear(); cached_homology_data_boxes.cache_clear()
    connection_braids=[connection_braid_analysis(rec['parameter_id'],analyses[rec['parameter_id']],nx,rec['pairs']) for rec in data.get('braid_frames',[])]
    bistable=[]; cycle=[]
    for r in regions:
        a=analyses[r['representative_parameter_id']]; attrs=[x for x in a['nodes'] if x['attractor']]
        if len(attrs)>=2: bistable.append(r['region_id'])
        if any(x['block_betti_F2'][1]>0 for x in attrs): cycle.append(r['region_id'])
    counts=Counter(ev['event_type'] for ev in edge_events)
    summary={'num_signatures':len(signatures),'num_regions':len(regions),'bistable_region_ids':bistable,
             'cycle_rich_attractor_region_ids':cycle,'parameter_edge_event_counts':{k:counts.get(k,0) for k in event_names}}
    return {'phase_grid':[nx,ny],'parameter_grid':[data['parameter_grid']['na'],data['parameter_grid']['nb']],
            'frames':frameout,'signatures':signatures,'regions':regions,'parameter_edge_events':edge_events,
            'continuation_transport':transport_queries,'region_adjacencies':regadj,'conley_index_barcodes':conley_paths,
            'morse_graph_barcodes':morse_graph_paths,'morse_graph_coherence':morse_graph_coherence,
            'parameter_patch_invariants':parameter_patches,'continuation_coherence':continuation_coherence,'connection_braids':connection_braids,
            'scientific_summary':summary}



def trajectory_core(query,frames,analyses,edge_rel,nx):
    ids=query['parameter_ids']; seed=query['seed_box_id']
    try:
        start=next(i for i,c in enumerate(analyses[ids[0]]['recurrent']) if seed in c)
    except StopIteration:
        raise ValueError('seed box is not recurrent at the first parameter')
    selected=[[start]]
    for p,q in zip(ids,ids[1:]):
        rel=relation_for_oriented_edge(p,q,edge_rel)
        prev=set(selected[-1])
        nxt=sorted({j for i,j in rel if i in prev})
        if not nxt: raise ValueError('trajectory selection dies')
        selected.append(nxt)
    dims_by_d=[[],[],[]]; maps_by_d=[[],[],[]]
    Bcur,Pcur,Ecur=selection_pair(ids[0],selected[0],frames,analyses)
    Hcur=get_relative_homology_data_boxes(Pcur,Ecur,nx)
    for d in range(3): dims_by_d[d].append(len(Hcur[d]['H']))
    for k in range(len(ids)-1):
        Bnext,Pnext,Enext=selection_pair(ids[k+1],selected[k+1],frames,analyses)
        Hnext=get_relative_homology_data_boxes(Pnext,Enext,nx)
        Pb=Pcur|Pnext; Eb=Ecur|Enext
        Hbridge=get_relative_homology_data_boxes(Pb,Eb,nx)
        for d in range(3):
            dims_by_d[d].append(len(Hbridge[d]['H']))
            dims_by_d[d].append(len(Hnext[d]['H']))
            maps_by_d[d].append(('forward',induced_box_map(Hcur,Hbridge,d)))
            maps_by_d[d].append(('backward',induced_box_map(Hnext,Hbridge,d)))
        Pcur,Ecur,Hcur=Pnext,Enext,Hnext
    return selected,dims_by_d,maps_by_d


def _relation_between_nodes(dims,maps,left,right):
    if left==right:
        return [(1<<i)|((1<<i)<<dims[left]) for i in range(dims[left])]
    R=[(1<<i)|((1<<i)<<dims[left]) for i in range(dims[left])]
    rd=dims[left]
    for k in range(left,right):
        E=_edge_linear_relation(dims[k],dims[k+1],maps[k][0],maps[k][1])
        # No early exit: the zero relation is not absorbing under relation
        # composition ({0} o S = {0} x S(0) can be nonzero).
        R=_compose_linear_relations(R,dims[left],rd,E,dims[k+1]); rd=dims[k+1]
    return R


def _diag_intersection_dim(R,d):
    diag=[(1<<i)|((1<<i)<<d) for i in range(d)]
    return len(_subspace_intersection_basis(R,diag))


def ranks_for_windows(dims,maps,windows):
    # Efficiently answer sparse generalized-rank queries by reusing prefixes
    # with the same left endpoint.
    by_start={}
    for idx,(a,b) in enumerate(windows): by_start.setdefault(a,[]).append((b,idx))
    ans=[0]*len(windows)
    for a,items in by_start.items():
        items=sorted(items); targets={b for b,_ in items}
        for b,idx in items:
            if b==a: ans[idx]=dims[a]
        if all(b==a for b,_ in items): continue
        R=[(1<<i)|((1<<i)<<dims[a]) for i in range(dims[a])]
        rd=dims[a]
        idxs={}
        for b,idx in items: idxs.setdefault(b,[]).append(idx)
        maxb=max(b for b,_ in items)
        for b in range(a+1,maxb+1):
            E=_edge_linear_relation(dims[b-1],dims[b],maps[b-1][0],maps[b-1][1])
            R=_compose_linear_relations(R,dims[a],rd,E,dims[b]); rd=dims[b]
            if b in targets:
                rr=_relation_regular_rank(R,dims[a],dims[b])
                for idx in idxs[b]: ans[idx]=rr
            if not R:
                # Only the *regular rank* is needed here.  Once the running
                # relation is zero its left projection is zero, and every later
                # composite {(0,w)} still has zero left projection, so all later
                # regular ranks are 0 (generalized ranks are monotone).  Loop and
                # holonomy relations never stop early: there the relation itself
                # is reported and {0} composed with S can be nonzero.
                break
    return ans



def relation_loop_signature(window,dims_by_d,maps_by_d,with_powers=True):
    a,b=window; left=2*a; right=2*b; out=[]
    for d in range(len(dims_by_d)):
        dims=dims_by_d[d]; maps=maps_by_d[d]
        if dims[left]!=dims[right]: raise ValueError('loop endpoint homology dimensions differ')
        R=_relation_between_nodes(dims,maps,left,right); n=dims[left]
        lp,rp=_relation_projection_ranks(R,n); rr=_relation_regular_rank(R,n,n)
        fixed=_diag_intersection_dim(R,n)
        rec={'dimension':d,'relation_dimension_F2':len(R),
             'left_projection_rank_F2':lp,'right_projection_rank_F2':rp,
             'regular_rank_F2':rr,'fixed_subspace_dimension_F2':fixed,
             'kronecker_invariants_F2':kronecker_type(R,n)}
        if with_powers:
            powers=[rr]; P=R
            for _ in range(3):
                P=_compose_linear_relations(P,n,n,P,n)
                powers.append(_relation_regular_rank(P,n,n))
            rec['power_regular_ranks_F2']=powers
        out.append(rec)
    return out



# Upper-case letters are the two closed-loop endorelations; lower-case letters
# are their converse relations.  Words are composed from left to right.
HOLONOMY_WORDS=['A','B','AB','BA','ABA','BAB','AABB','BBAA','ABAB','BABA','Ab','aB','ABab','abAB']

def holonomy_word_signature(probe,dims_by_d,maps_by_d,with_powers=True):
    """Basis-free invariants of the semigroup generated by two closed continuation loops.

    probe is [a,b,c] with the same endpoint state at all three occurrences.
    A is the relation carried by a..b and B by b..c.  The reported words expose
    information that cannot be recovered from the ranks of A and B separately.
    """
    a,b,c=probe
    words=HOLONOMY_WORDS
    out=[]
    for d in range(len(dims_by_d)):
        dims=dims_by_d[d]; maps=maps_by_d[d]
        ia,ib,ic=2*a,2*b,2*c
        n=dims[ia]
        if dims[ib]!=n or dims[ic]!=n:
            raise ValueError('holonomy probe endpoint homology dimensions differ')
        A=_relation_between_nodes(dims,maps,ia,ib)
        B=_relation_between_nodes(dims,maps,ib,ic)
        rels={'A':A,'B':B,'a':kr_converse(A,n),'b':kr_converse(B,n)}
        sig=[]
        for word in words:
            R=[(1<<i)|((1<<i)<<n) for i in range(n)]
            for ch in word:
                R=_compose_linear_relations(R,n,n,rels[ch],n)
            lp,rp=_relation_projection_ranks(R,n)
            sig.append({'word':word,'relation_dimension_F2':len(R),
                        'left_projection_rank_F2':lp,'right_projection_rank_F2':rp,
                        'regular_rank_F2':_relation_regular_rank(R,n,n),
                        'fixed_subspace_dimension_F2':_diag_intersection_dim(R,n),
                        'kronecker_invariants_F2':kronecker_type(R,n)})
        out.append({'dimension':d,'space_dimension_F2':n,'word_signatures':sig})
    return out



# ---------------------------------------------------------------------------
# Kronecker invariants of an endorelation R subset U (+) U over F2.
# R is a list of bit vectors: low n bits = left coordinate, next n bits = right.
# The invariants are computed from four functorial chains of subspaces
#   M_j = R^j(0), K_j = R^{-j}(0), R^k(U), R^{-k}(U)
# whose dimensions are additive over the Kronecker decomposition, and from the
# elementary divisors of the automorphism induced on the regular core.
# ---------------------------------------------------------------------------

def _k_ech(vecs):
    piv = {}
    for v in vecs:
        while v:
            p = v.bit_length() - 1
            if p in piv:
                v ^= piv[p]
            else:
                piv[p] = v
                break
    return piv

def _k_basis(vecs):
    return list(_k_ech(vecs).values())

def _k_dim(vecs):
    return len(_k_ech(vecs))

def _k_reduce(piv, v):
    while v:
        p = v.bit_length() - 1
        if p in piv:
            v ^= piv[p]
        else:
            break
    return v

def _k_intersect(A, B):
    """Basis of span(A) cap span(B)."""
    A = _k_basis(A); B = _k_basis(B)
    if not A or not B:
        return []
    # Track combinations of A that land in span(B).
    piv = {}
    out = []
    for v in B:
        lab = 0
        while v:
            p = v.bit_length() - 1
            if p in piv:
                vv, ll = piv[p]; v ^= vv; lab ^= ll
            else:
                piv[p] = (v, lab); break
    for a in A:
        v = a; lab = a
        while v:
            p = v.bit_length() - 1
            if p in piv:
                vv, ll = piv[p]; v ^= vv; lab ^= ll
            else:
                piv[p] = (v, lab); break
        if v == 0:
            out.append(lab)
    return _k_basis(out)

def _k_image(R, n, S):
    """R(S) = {y : (x,y) in R, x in S}."""
    m = (1 << n) - 1
    piv = {}
    for s in S:
        v = s; lab = 0
        while v:
            p = v.bit_length() - 1
            if p in piv:
                vv, ll = piv[p]; v ^= vv; lab ^= ll
            else:
                piv[p] = (v, lab); break
    out = []
    for r in R:
        v = r & m; lab = r >> n
        while v:
            p = v.bit_length() - 1
            if p in piv:
                vv, ll = piv[p]; v ^= vv; lab ^= ll
            else:
                piv[p] = (v, lab); break
        if v == 0 and lab:
            out.append(lab)
    return _k_basis(out)

def kr_converse(R, n):
    m = (1 << n) - 1
    return [(r >> n) | ((r & m) << n) for r in R]

def _k_chain(R, n, start):
    seq = [_k_basis(start)]
    while True:
        nxt = _k_image(R, n, seq[-1])
        if len(nxt) == len(seq[-1]):
            # Chains used here are monotone, so equal dimension means stable.
            return seq
        seq.append(nxt)

def _k_counts_from_cumulative(f):
    """f[j] = sum_i min(j, t_i).  Return sorted list of t_i (t_i >= 1)."""
    ge = [f[j] - f[j - 1] for j in range(1, len(f))]  # ge[j-1] = #{t >= j}
    out = []
    for j in range(1, len(ge) + 1):
        c = ge[j - 1] - (ge[j] if j < len(ge) else 0)
        out.extend([j] * c)
    return out

# ---------- F2[x] polynomial arithmetic (ints, bit i = coeff of x^i) ----------

def kr_pdeg(a):
    return a.bit_length() - 1

def kr_pmul(a, b):
    r = 0
    while b:
        if b & 1:
            r ^= a
        a <<= 1; b >>= 1
    return r

def kr_pdivmod(a, b):
    q = 0; db = kr_pdeg(b)
    while a and kr_pdeg(a) >= db:
        s = kr_pdeg(a) - db
        q ^= 1 << s; a ^= b << s
    return q, a

def kr_pmod(a, b):
    db = kr_pdeg(b)
    while a and kr_pdeg(a) >= db:
        a ^= b << (kr_pdeg(a) - db)
    return a

def kr_pgcd(a, b):
    while b:
        a, b = b, kr_pmod(a, b)
    return a

def kr_pmulmod(a, b, f):
    return kr_pmod(kr_pmul(a, b), f)

def kr_pderiv(a):
    r = 0; k = 1
    a >>= 1
    while a:
        if a & 1 and k & 1:
            r |= 1 << (k - 1)
        a >>= 1; k += 1
    return r

def kr_psqrt(a):
    r = 0; k = 0
    while a:
        if a & 1:
            r |= 1 << (k // 2)
        a >>= 1; k += 1
    return r

def kr_squarefree_decomposition(f):
    """Return list of (g, e) with f = prod g^e, g squarefree (not nec. irreducible)."""
    out = []
    if kr_pdeg(f) <= 0:
        return out
    d = kr_pderiv(f)
    if d == 0:
        for g, e in kr_squarefree_decomposition(kr_psqrt(f)):
            out.append((g, 2 * e))
        return out
    c = kr_pgcd(f, d)
    w = kr_pdivmod(f, c)[0]
    i = 1
    while kr_pdeg(w) > 0:
        y = kr_pgcd(w, c)
        z = kr_pdivmod(w, y)[0]
        if kr_pdeg(z) > 0:
            out.append((z, i))
        i += 1
        w = y
        c = kr_pdivmod(c, y)[0]
    if kr_pdeg(c) > 0:
        for g, e in kr_squarefree_decomposition(kr_psqrt(c)):
            out.append((g, 2 * e))
    return out

def kr_ddf(f):
    out = []
    h = 2  # x
    d = 0
    while kr_pdeg(f) >= 2 * (d + 1):
        d += 1
        h = kr_pmulmod(h, h, f)
        g = kr_pgcd(h ^ 2, f)
        if kr_pdeg(g) > 0:
            out.append((g, d))
            f = kr_pdivmod(f, g)[0]
            h = kr_pmod(h, f)
    if kr_pdeg(f) > 0:
        out.append((f, kr_pdeg(f)))
    return out

def kr_edf(f, d, seed=[0x9E3779B97F4A7C15]):
    if kr_pdeg(f) == d:
        return [f]
    n = kr_pdeg(f)
    while True:
        seed[0] = (seed[0] * 6364136223846793005 + 1442695040888963407) & ((1 << 64) - 1)
        a = (seed[0] ^ (seed[0] << 64) ^ (seed[0] << 128) ^ (seed[0] << 192)) & ((1 << n) - 1)
        if kr_pdeg(a) <= 0:
            continue
        t = a; s = a
        for _ in range(d - 1):
            s = kr_pmulmod(s, s, f); t ^= s
        g = kr_pgcd(t, f)
        if 0 < kr_pdeg(g) < n:
            return kr_edf(g, d) + kr_edf(kr_pdivmod(f, g)[0], d)

def kr_factor(f):
    """Irreducible factorisation of f over F2: dict irreducible -> multiplicity."""
    out = {}
    for g, e in kr_squarefree_decomposition(f):
        for h, d in kr_ddf(g):
            for p in kr_edf(h, d):
                out[p] = out.get(p, 0) + e
    return out

# ---------- matrices over F2 as column lists ----------

def kr_mat_apply(cols, v):
    r = 0; j = 0
    while v:
        if v & 1:
            r ^= cols[j]
        v >>= 1; j += 1
    return r

def kr_mat_mul(A, B):
    """(A B) as columns; A, B column lists of square matrices."""
    return [kr_mat_apply(A, c) for c in B]

def kr_charpoly(cols, m):
    """Characteristic polynomial via Danilevsky-free Krylov/Frobenius blocks."""
    # Build a cyclic decomposition V = sum Krylov subspaces (not nec. canonical);
    # char poly is product of the local min polys of the quotient chain.
    piv = {}
    total = 1
    remaining = list(range(m))
    for s in range(m):
        e = 1 << s
        if _k_reduce({k: v for k, (v, _) in piv.items()}, e) == 0:
            continue
        # Krylov sequence of e modulo the span so far
        seq_piv = dict(piv)
        v = e; k = 0
        # labels track polynomial coefficient of the Krylov vector
        while True:
            vv = v; lab = 1 << k
            while vv:
                p = vv.bit_length() - 1
                if p in seq_piv:
                    a, l = seq_piv[p]; vv ^= a; lab ^= l
                else:
                    break
            if vv == 0:
                # lab has x^k + lower (relative to base span, whose labels are 0)
                total = kr_pmul(total, lab)
                break
            seq_piv[vv.bit_length() - 1] = (vv, lab)
            v = kr_mat_apply(cols, v); k += 1
        # the new base span: keep vectors but labels zero
        piv = {p: (a, 0) for p, (a, _) in seq_piv.items()}
    return total

def kr_poly_of_matrix(p, cols, m):
    """Columns of p(T)."""
    res = [0] * m
    ident = [1 << i for i in range(m)]
    for k in range(kr_pdeg(p), -1, -1):
        res = kr_mat_mul(cols, res) if any(res) else [0] * m
        if (p >> k) & 1:
            res = [res[i] ^ ident[i] for i in range(m)]
    return res

def kr_elementary_divisors(cols, m):
    """Sorted list of [p_code, exponent] for the linear map with given columns."""
    if m == 0:
        return []
    chi = kr_charpoly(cols, m)
    out = []
    for p, e in sorted(kr_factor(chi).items()):
        P = kr_poly_of_matrix(p, cols, m)
        d = kr_pdeg(p)
        kers = [0]
        Q = [1 << i for i in range(m)]
        while True:
            Q = kr_mat_mul(P, Q)
            k = m - _k_dim(Q)
            if k == kers[-1]:
                break
            kers.append(k)
        # kers[j] = dim ker p(T)^j = d * sum_t min(j,t)
        f = [k // d for k in kers]
        for t in _k_counts_from_cumulative(f):
            out.append([p, t])
        assert sum(t for _, t in out if _ == p) == e, (p, e, kers)
    out.sort()
    return out

def kronecker_type(R, n):
    """Complete conjugacy invariant of the endorelation R on F2^n."""
    R = _k_basis(R)
    Rc = kr_converse(R, n)
    full = [1 << i for i in range(n)]
    M = _k_chain(R, n, [])        # M_j = R^j(0)
    K = _k_chain(Rc, n, [])       # K_j = R^{-j}(0)
    Ran = _k_chain(R, n, full)    # R^k(U)
    D = _k_chain(Rc, n, full)     # R^{-k}(U)
    Minf, Kinf, Raninf, Dinf = M[-1], K[-1], Ran[-1], D[-1]
    # right minimal indices
    f = [len(_k_intersect(Mj, Kinf)) for Mj in M]
    right = _k_counts_from_cumulative(f)
    # left minimal indices: h(k) = n - dim(R^k(U) + D_inf); #{eta >= k-1} = h(k)-h(k-1)
    h = [n - _k_dim(Rk + Dinf) for Rk in Ran]
    left = []
    ge = [h[k] - h[k - 1] for k in range(1, len(h))]
    for k in range(1, len(ge) + 1):
        c = ge[k - 1] - (ge[k] if k < len(ge) else 0)
        left.extend([k - 1] * c)
    # infinite elementary divisors
    finf = [len(Mj) - len(_k_intersect(Mj, Kinf)) for Mj in M]
    inf = _k_counts_from_cumulative(finf)
    # nilpotent (p = x) part
    f0 = [len(Kj) - len(_k_intersect(Kj, Minf)) for Kj in K]
    nil = _k_counts_from_cumulative(f0)
    # invertible finite part: T on X/Y, X = Dinf cap Raninf, Y = (Minf+Kinf) cap X
    X = _k_intersect(Dinf, Raninf)
    Y = _k_intersect(_k_basis(Minf + Kinf), X)
    Ypiv = _k_ech(Y)
    comp = []
    cp = dict(Ypiv)
    for x in X:
        v = _k_reduce(cp, x)
        if v:
            cp[v.bit_length() - 1] = v; comp.append(v)
    m = len(comp)
    coord = {}
    for y in Y:
        v = y; lab = 0
        while v:
            p = v.bit_length() - 1
            if p in coord:
                a, l = coord[p]; v ^= a; lab ^= l
            else:
                coord[p] = (v, lab); break
    for i, c in enumerate(comp):
        v = c; lab = 1 << i
        while v:
            p = v.bit_length() - 1
            if p in coord:
                a, l = coord[p]; v ^= a; lab ^= l
            else:
                coord[p] = (v, lab); break
    def co(v):
        lab = 0
        while v:
            p = v.bit_length() - 1
            a, l = coord[p]; v ^= a; lab ^= l
        return lab
    XX = [x for x in X] + [x << n for x in X]
    RX = _k_intersect(R, XX)
    mask = (1 << n) - 1
    piv = {}
    for r in RX:
        a = co(r & mask); b = co(r >> n)
        while a:
            p = a.bit_length() - 1
            if p in piv:
                aa, bb = piv[p]; a ^= aa; b ^= bb
            else:
                piv[p] = (a, b); break
        if a == 0:
            assert b == 0, 'core relation is not a graph'
    assert len(piv) == m, 'core relation not total'
    # reduce to columns: solve for e_i
    cols = []
    for i in range(m):
        a = 1 << i; b = 0
        while a:
            p = a.bit_length() - 1
            aa, bb = piv[p]; a ^= aa; b ^= bb
        cols.append(b)
    fin = kr_elementary_divisors(cols, m)
    fin = sorted(fin + [[2, t] for t in nil])
    # consistency
    dimR = len(R)
    assert n == sum(right) + sum(e + 1 for e in left) + sum(kr_pdeg(p) * t for p, t in fin) + sum(inf), 'U dim'
    assert dimR == sum(e + 1 for e in right) + sum(left) + sum(kr_pdeg(p) * t for p, t in fin) + sum(inf), 'R dim'
    return {'right_minimal_indices': sorted(right), 'left_minimal_indices': sorted(left),
            'infinite_elementary_divisors': sorted(inf), 'finite_elementary_divisors': fin}


# ---------------------------------------------------------------------------
# Full zigzag barcode.
# ---------------------------------------------------------------------------
# Full zigzag barcode of a path zigzag by a right sweep over the two flags
#   E_b = ran(R_{b,k})   and   N_b = mul(R_{b,k})   inside V_k,
# where R_{b,k} is the composed linear relation from node b to node k.
# rho(b,k) = dim E_b - dim N_b is the generalized rank on [b,k].

def _z_ech_insert(piv, v):
    while v:
        p = v.bit_length() - 1
        if p in piv:
            v ^= piv[p]
        else:
            piv[p] = v
            return True
    return False

def _z_basis(vecs):
    piv = {}
    for v in vecs:
        _z_ech_insert(piv, v)
    return tuple(sorted(piv.values()))

def _z_canon(vecs):
    """Reduced row echelon form: a canonical, hashable representative."""
    piv = {}
    for v in vecs:
        _z_ech_insert(piv, v)
    keys = sorted(piv)
    for i, p in enumerate(keys):
        for q in keys[i + 1:]:
            if (piv[q] >> p) & 1:
                piv[q] ^= piv[p]
    return tuple(piv[p] for p in keys)

def _z_apply(cols, v):
    r = 0; j = 0
    while v:
        if v & 1:
            r ^= cols[j]
        v >>= 1; j += 1
    return r

def _z_forward_image(cols, S):
    return _z_canon([_z_apply(cols, v) for v in S])

def _z_preimage(cols, dsrc, S):
    """{v in F2^dsrc : g(v) in span S}, g given by columns."""
    piv = {}
    for s in S:
        v = s
        while v:
            p = v.bit_length() - 1
            if p in piv:
                v ^= piv[p][0]
            else:
                piv[p] = (v, 0); break
    out = []
    for j in range(dsrc):
        v = cols[j]; lab = 1 << j
        while v:
            p = v.bit_length() - 1
            if p in piv:
                vv, ll = piv[p]; v ^= vv; lab ^= ll
            else:
                piv[p] = (v, lab); break
        if v == 0:
            out.append(lab)
    return _z_canon(out)

def zigzag_barcode(dims, maps):
    """dims: node dimensions; maps[k] = ('forward', cols V_k->V_{k+1}) or
    ('backward', cols V_{k+1}->V_k).  Returns sorted [[start,end,mult],...]."""
    n = len(dims)
    if n == 0:
        return []
    # segs: list of [b_lo, E, N] for b in [b_lo, next b_lo), increasing b.
    full = lambda d: _z_canon([1 << i for i in range(d)])
    segs = [[0, full(dims[0]), ()]]
    bars = {}

    def rho_profile(sg):
        return [(b, len(E) - len(N)) for b, E, N in sg]

    for k in range(n):
        if k + 1 < n:
            direction, cols = maps[k]
            new = []
            for b, E, N in segs:
                if direction == 'forward':
                    E2 = _z_forward_image(cols, E); N2 = _z_forward_image(cols, N)
                else:
                    E2 = _z_preimage(cols, dims[k + 1], E); N2 = _z_preimage(cols, dims[k + 1], N)
                if new and new[-1][1] == E2 and new[-1][2] == N2:
                    continue
                new.append([b, E2, N2])
            F = full(dims[k + 1])
            if new and new[-1][1] == F and new[-1][2] == ():
                pass
            else:
                new.append([k + 1, F, ()])
        else:
            new = []
        # m(b,k) = [rho_k(b)-rho_k(b-1)] - [rho_{k+1}(b)-rho_{k+1}(b-1)], b<=k
        def jumps(sg, upto):
            out = {}; prev = 0
            for b, E, N in sg:
                if b > upto:
                    break
                r = len(E) - len(N)
                if r != prev:
                    out[b] = out.get(b, 0) + r - prev
                prev = r
            return out
        jk = jumps(segs, k); jn = jumps(new, k)
        for b in set(jk) | set(jn):
            m = jk.get(b, 0) - jn.get(b, 0)
            if m < 0:
                raise RuntimeError('negative zigzag multiplicity')
            if m:
                bars[(b, k)] = m
        segs = new
    return [[a, b, m] for (a, b), m in sorted(bars.items())]


def graph_trajectory_core(ids,analyses,edge_struct):
    endpointH={}
    for p in set(ids):
        endpointH[p]=graph_homology_data(len(analyses[p]['recurrent']),order_graph_edges(analyses[p]['desc']))
    union_cache={}
    for p,q in zip(ids,ids[1:]):
        a,b=sorted((p,q))
        if (a,b) in union_cache: continue
        key=(a,b) if (a,b) in edge_struct else (b,a)
        U=edge_struct[key]
        union_cache[(a,b)]=(key,U,graph_homology_data(len(U['recurrent']),order_graph_edges(U['desc'])))
    Hnodes=[]
    for i,p in enumerate(ids):
        Hnodes.append(endpointH[p])
        if i+1<len(ids):
            a,b=sorted((p,ids[i+1])); Hnodes.append(union_cache[(a,b)][2])
    maps_by_d=[[],[]]
    for d in range(2):
        for p,q in zip(ids,ids[1:]):
            a,b=sorted((p,q)); key,U,HU=union_cache[(a,b)]
            if key==(p,q): pmap,qmap=U['source_map'],U['target_map']
            else: pmap,qmap=U['target_map'],U['source_map']
            maps_by_d[d].append(('forward',induced_graph_map(endpointH[p],HU,d,pmap)))
            maps_by_d[d].append(('backward',induced_graph_map(endpointH[q],HU,d,qmap)))
    dims_by_d=[[len(H[d]['H']) for H in Hnodes] for d in range(2)]
    return dims_by_d,maps_by_d



def solve_case(data):
    nx=data['phase_grid']['nx']; ny=data['phase_grid']['ny']
    frames=build_frames(data); n=nx*ny
    analyses=[analyze_frame(f,nx,ny) for f in frames]
    edge_rel={}; edge_struct={}
    for p,q in data['parameter_edges']:
        rel,U=union_edge_structure(frames[p],analyses[p],frames[q],analyses[q],n)
        edge_rel[(p,q)]=rel; edge_struct[(p,q)]=U
    queries=[]
    for q in data['queries']:
        selected,dims,maps=trajectory_core(q,frames,analyses,edge_rel,nx)
        z=2*len(q['parameter_ids'])-1
        windows=q.get('rank_windows',[])
        for a,b in windows:
            if not (0<=a<=b<z): raise ValueError('rank window bounds')
        rank_by_d=[ranks_for_windows(dims[d],maps[d],windows) for d in range(3)]
        rank_records=[]
        for k,(a,b) in enumerate(windows):
            rank_records.append({'window_id':k,'start_node':a,'end_node':b,
                                 'rank_F2':[rank_by_d[d][k] for d in range(3)]})
        loops=[]
        for lid,w in enumerate(q.get('loop_windows',[])):
            a,b=w
            if not (0<=a<b<len(q['parameter_ids'])): raise ValueError('loop window bounds')
            if q['parameter_ids'][a]!=q['parameter_ids'][b]: raise ValueError('loop is not closed in parameter space')
            if selected[a]!=selected[b]: raise ValueError('loop endpoints do not carry the same Morse selection')
            loops.append({'loop_id':lid,'start_occurrence':a,'end_occurrence':b,
                          'degrees':relation_loop_signature(w,dims,maps)})
        barcodes=[zigzag_barcode(dims[d],maps[d]) for d in range(3)]
        holonomy=[]
        for hid,probe in enumerate(q.get('holonomy_probes',[])):
            a,b,c=probe
            if not (0<=a<b<c<len(q['parameter_ids'])): raise ValueError('holonomy probe bounds')
            if not (q['parameter_ids'][a]==q['parameter_ids'][b]==q['parameter_ids'][c]):
                raise ValueError('holonomy probe is not closed at one parameter')
            if not (selected[a]==selected[b]==selected[c]):
                raise ValueError('holonomy probe endpoints do not carry the same Morse selection')
            holonomy.append({'probe_id':hid,'occurrences':[a,b,c],
                             'degrees':holonomy_word_signature(probe,dims,maps)})

        gdims,gmaps=graph_trajectory_core(q['parameter_ids'],analyses,edge_struct)
        gwindows=q.get('graph_rank_windows',[])
        for a,b in gwindows:
            if not (0<=a<=b<z): raise ValueError('graph rank window bounds')
        grank_by_d=[ranks_for_windows(gdims[d],gmaps[d],gwindows) for d in range(2)]
        granks=[]
        for k,(a,b) in enumerate(gwindows):
            granks.append({'window_id':k,'start_node':a,'end_node':b,
                           'rank_F2':[grank_by_d[d][k] for d in range(2)]})
        gbarcodes=[zigzag_barcode(gdims[d],gmaps[d]) for d in range(2)]
        gloops=[]
        for lid,w in enumerate(q.get('graph_loop_windows',q.get('loop_windows',[])[:4])):
            a,b=w
            if q['parameter_ids'][a]!=q['parameter_ids'][b]: raise ValueError('graph loop is not closed')
            gloops.append({'loop_id':lid,'start_occurrence':a,'end_occurrence':b,
                           'degrees':relation_loop_signature(w,gdims,gmaps,False)})
        gholonomy=[]
        for hid,probe in enumerate(q.get('graph_holonomy_probes',[])):
            a,b,c=probe
            if not (0<=a<b<c<len(q['parameter_ids'])): raise ValueError('graph holonomy probe bounds')
            if not (q['parameter_ids'][a]==q['parameter_ids'][b]==q['parameter_ids'][c]):
                raise ValueError('graph holonomy probe is not closed at one parameter')
            gholonomy.append({'probe_id':hid,'occurrences':[a,b,c],
                              'degrees':holonomy_word_signature(probe,gdims,gmaps,False)})

        queries.append({'query_id':q['query_id'],'parameter_ids':q['parameter_ids'],
                        'seed_box_id':q['seed_box_id'],'selected_morse_ids_by_frame':selected,
                        'conley_node_dimensions_F2':[dims[d] for d in range(3)],
                        'conley_generalized_rank_queries':rank_records,
                        'conley_zigzag_barcodes_F2':barcodes,'conley_loop_signatures':loops,
                        'conley_holonomy_word_signatures':holonomy,
                        'morse_graph_node_dimensions_F2':[gdims[d] for d in range(2)],
                        'morse_graph_generalized_rank_queries':granks,
                        'morse_graph_zigzag_barcodes_F2':gbarcodes,'morse_graph_loop_signatures':gloops,
                        'morse_graph_holonomy_word_signatures':gholonomy})
        cached_relative_homology_data_boxes.cache_clear(); cached_homology_data_boxes.cache_clear()
    return {'case_id':data['case_id'],'queries':queries}


def main():
    import argparse
    ap=argparse.ArgumentParser(); ap.add_argument('input'); ap.add_argument('output')
    args=ap.parse_args()
    with open(args.input) as f: data=json.load(f)
    ans=solve_case(data)
    with open(args.output,'w') as f: json.dump(ans,f,separators=(',',':'),sort_keys=False)

if __name__=='__main__': main()

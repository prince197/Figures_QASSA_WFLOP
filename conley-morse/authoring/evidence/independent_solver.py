#!/usr/bin/env python3
"""Independent verification solver for the Conley-Morse continuation task.

Independent of solution/trajectory_solver.py in every layer that the graded
records depend on beyond cubical homology:
  * box maps use pure integer arithmetic (no Fraction objects);
  * strongly connected components use an iterative Kosaraju route;
  * every relation between zigzag nodes is the projection of the section space
    (inverse limit) of the restricted zigzag, computed chunk-wise and glued by
    pullbacks, with no early exit of any kind;
  * Kronecker invariants use kernel dimensions of the pencil beta - s*alpha over
    polynomial quotient rings and polynomial-kernel (Toeplitz) ranks, not the
    subspace chains used by the reference;
  * zigzag barcodes use per-start generalized-rank sweeps and inclusion-exclusion,
    not the reference's single right-to-left flag sweep.
The relative cubical homology and the chain-level inclusion maps descend from the
original author's implementation; they are cross-checked by the generalized-rank,
dimension and barcode gates, which every external solver has reproduced.
"""
import json, sys
from functools import lru_cache
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

def tarjan_scc(n, offsets, targets):
    # Independent SCC route: iterative Kosaraju rather than Tarjan.
    adj=[targets[offsets[v]:offsets[v+1]] for v in range(n)]
    radj=[[] for _ in range(n)]
    for v in range(n):
        for w in adj[v]: radj[w].append(v)
    seen=[False]*n; order=[]
    for s in range(n):
        if seen[s]: continue
        seen[s]=True; st=[(s,0)]
        while st:
            v,k=st[-1]
            if k<len(adj[v]):
                w=adj[v][k]; st[-1]=(v,k+1)
                if not seen[w]: seen[w]=True; st.append((w,0))
            else:
                order.append(v); st.pop()
    seen=[False]*n; comps=[]
    for s in reversed(order):
        if seen[s]: continue
        seen[s]=True; comp=[]; st=[s]
        while st:
            v=st.pop(); comp.append(v)
            for w in radj[v]:
                if not seen[w]: seen[w]=True; st.append(w)
        comps.append(comp)
    return comps

def cells_of_boxes(boxes,nx):
    verts=set(); edges=set(); squares=set()
    for q in boxes:
        i,j=q%nx,q//nx
        squares.add((i,j))
        verts.update(((i,j),(i+1,j),(i,j+1),(i+1,j+1)))
        edges.update(((0,i,j),(0,i,j+1),(1,i,j),(1,i+1,j)))
    return verts,edges,squares

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

def relation_for_oriented_edge(p,q,edge_rel):
    if (p,q) in edge_rel: return edge_rel[(p,q)]
    return [[j,i] for i,j in edge_rel[(q,p)]]

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

# ===========================================================================
# Independent routes.  Nothing below is shared with solution/trajectory_solver.py.
# ===========================================================================

# ---------------------------------------------------------------------------
# Exact box map with integer arithmetic only (no Fraction objects).
# x-factor on source box i: g(t) = r_num * t * (nx - t) / (scale * nx^2), t in [i, i+1].
# Image coordinate X = ((cd-cn)*gx + cn*gy)/cd.  Target box ii intersects [lo,hi]
# (closed) iff ii <= hi*n and ii+1 >= lo*n.
# ---------------------------------------------------------------------------

def _factor_range_int(num, n, k):
    """Return (lo, hi) of t*(n-t) over t in [k, k+1] as integers scaled by 4
    (so that the critical value n^2/4 stays integral)."""
    a = 4 * k * (n - k); b = 4 * (k + 1) * (n - k - 1)
    vals = [a, b]
    if 2 * k <= n <= 2 * k + 2:
        vals.append(n * n)
    return min(vals) * num, max(vals) * num


def build_frames_int(data):
    nx = data['phase_grid']['nx']; ny = data['phase_grid']['ny']
    pg = data['parameter_grid']; model = data['model']
    scale = model['parameter_scale']; cn = model['coupling_num']; cd = model['coupling_den']
    # value of factor = F / (4 * scale * n^2) with F from _factor_range_int.
    # X * n = ((cd-cn)*Fx/(4 scale nx^2) + cn*Fy/(4 scale ny^2)) * n / cd
    frames = []
    for sj, s in enumerate(pg['s_num']):
        yr = [_factor_range_int(s, ny, j) for j in range(ny)]
        for ri, r in enumerate(pg['r_num']):
            xr = [_factor_range_int(r, nx, i) for i in range(nx)]
            offsets = [0]; targets = []
            for j in range(ny):
                ly, hy = yr[j]
                for i in range(nx):
                    lx, hx = xr[i]
                    rng = []
                    for (w1, w2, n) in (((cd - cn), cn, nx), (cn, (cd - cn), ny)):
                        # common denominator: 4*scale*nx^2*ny^2*cd ; value*n compared to integers
                        den = 4 * scale * nx * nx * ny * ny * cd
                        lo_num = (w1 * lx * ny * ny + w2 * ly * nx * nx) * n
                        hi_num = (w1 * hx * ny * ny + w2 * hy * nx * nx) * n
                        # smallest box index ii with ii+1 >= lo*n  -> ii >= lo*n - 1
                        a = -((-(lo_num - den)) // den)
                        b = hi_num // den
                        rng.append((max(0, a), min(n - 1, b)))
                    (i0, i1), (j0, j1) = rng
                    for jj in range(j0, j1 + 1):
                        for ii in range(i0, i1 + 1):
                            targets.append(ii + nx * jj)
                    offsets.append(len(targets))
            frames.append({'offsets': offsets, 'targets': targets})
    return frames


def analyze_frame_ind(frame, nx, ny):
    n = nx * ny; off = frame['offsets']; tg = frame['targets']
    recurrent, desc, _, _, _ = recurrent_structure(n, off, tg)
    pairs = []
    for c in recurrent:
        P = set(c); st = list(c)
        while st:
            v = st.pop()
            for k in range(off[v], off[v + 1]):
                w = tg[k]
                if w not in P:
                    P.add(w); st.append(w)
        pairs.append({'P': P, 'E': P - set(c)})
    return {'recurrent': recurrent, 'desc': desc, 'pairs': pairs}


def propagate_selection(ids, seed, analyses, edge_rel):
    start = [i for i, c in enumerate(analyses[ids[0]]['recurrent']) if seed in c]
    if len(start) != 1:
        raise ValueError('seed box is not recurrent at the first parameter')
    cur = set(start); out = [sorted(cur)]
    for p, q in zip(ids, ids[1:]):
        rel = relation_for_oriented_edge(p, q, edge_rel)
        cur = {j for i, j in rel if i in cur}
        if not cur:
            raise ValueError('selection dies')
        out.append(sorted(cur))
    return out


def conley_zigzag_ind(ids, selected, analyses, nx):
    def pair(pid, sel):
        B = set(); P = set()
        for m in sel:
            B |= set(analyses[pid]['recurrent'][m]); P |= analyses[pid]['pairs'][m]['P']
        return P, P - B
    dims = [[], [], []]; maps = [[], [], []]
    P0, E0 = pair(ids[0], selected[0]); H0 = get_relative_homology_data_boxes(P0, E0, nx)
    for d in range(3):
        dims[d].append(len(H0[d]['H']))
    for k in range(len(ids) - 1):
        P1, E1 = pair(ids[k + 1], selected[k + 1]); H1 = get_relative_homology_data_boxes(P1, E1, nx)
        Hb = get_relative_homology_data_boxes(P0 | P1, E0 | E1, nx)
        for d in range(3):
            dims[d].append(len(Hb[d]['H'])); dims[d].append(len(H1[d]['H']))
            maps[d].append(('forward', induced_box_map(H0, Hb, d)))
            maps[d].append(('backward', induced_box_map(H1, Hb, d)))
        P0, E0, H0 = P1, E1, H1
    return dims, maps


def graph_zigzag_ind(ids, analyses, edge_struct):
    H = {p: graph_homology_data(len(analyses[p]['recurrent']), order_graph_edges(analyses[p]['desc'])) for p in set(ids)}
    dims = [[], []]; maps = [[], []]
    for d in range(2):
        dims[d].append(len(H[ids[0]][d]['H']))
    for p, q in zip(ids, ids[1:]):
        if (p, q) in edge_struct:
            U = edge_struct[(p, q)]; pm, qm = U['source_map'], U['target_map']
        else:
            U = edge_struct[(q, p)]; pm, qm = U['target_map'], U['source_map']
        HU = graph_homology_data(len(U['recurrent']), order_graph_edges(U['desc']))
        for d in range(2):
            dims[d].append(len(HU[d]['H'])); dims[d].append(len(H[q][d]['H']))
            maps[d].append(('forward', induced_graph_map(H[p], HU, d, pm)))
            maps[d].append(('backward', induced_graph_map(H[q], HU, d, qm)))
    return dims, maps


# ---------------------------------------------------------------------------
# Linear relations as explicit section spaces.
# A relation between nodes a and b is the projection to (x_a, x_b) of the space
# of compatible families (x_a, ..., x_b) on the restricted zigzag.  Long windows
# are split into chunks; chunk relations are glued by a pullback over the shared
# node.  There is no early exit anywhere: every chunk is always glued.
# ---------------------------------------------------------------------------

class Rel:
    __slots__ = ('du', 'dv', 'vecs')
    def __init__(self, du, dv, vecs):
        self.du = du; self.dv = dv
        piv = {}
        for v in vecs:
            while v:
                p = v.bit_length() - 1
                if p in piv:
                    v ^= piv[p]
                else:
                    piv[p] = v; break
        self.vecs = [piv[p] for p in sorted(piv)]
    def left(self, v):
        return v & ((1 << self.du) - 1)
    def right(self, v):
        return v >> self.du


def _nullspace(cols):
    """Basis (as coefficient masks) of {c : sum c_j cols[j] = 0}."""
    piv = {}; out = []
    for j, v in enumerate(cols):
        lab = 1 << j
        while v:
            p = v.bit_length() - 1
            if p in piv:
                a, l = piv[p]; v ^= a; lab ^= l
            else:
                piv[p] = (v, lab); break
        if v == 0:
            out.append(lab)
    return out


def sections_relation(dims, maps, a, b):
    """Endpoint relation of the zigzag restricted to nodes a..b via its section space."""
    offs = [0]
    for k in range(a, b + 1):
        offs.append(offs[-1] + dims[k])
    # one constraint block per arrow, living in the arrow target
    cons_off = []; nc = 0
    for k in range(a, b):
        tgt = k + 1 if maps[k][0] == 'forward' else k
        cons_off.append(nc); nc += dims[tgt]
    cols = [0] * offs[-1]
    for t, k in enumerate(range(a, b)):
        direction, mc = maps[k]
        src, tgt = (k, k + 1) if direction == 'forward' else (k + 1, k)
        so = offs[src - a]; to = offs[tgt - a]; co = cons_off[t]
        for j, c in enumerate(mc):
            cols[so + j] ^= c << co
        for r in range(dims[tgt]):
            cols[to + r] ^= 1 << (co + r)
    du = dims[a]; dv = dims[b]; ro = offs[-2]
    out = []
    for sec in _nullspace(cols):
        # sec is a coefficient mask over the concatenated node coordinates
        x = sec & ((1 << du) - 1); y = (sec >> ro) & ((1 << dv) - 1)
        out.append(x | (y << du))
    return Rel(du, dv, out)


def glue(R, S):
    """Pullback composition {(x,z): exists y, (x,y) in R, (y,z) in S}."""
    assert R.dv == S.du
    dy = R.dv
    ym = (1 << dy) - 1
    cols = [R.right(v) for v in R.vecs] + [S.left(v) & ym for v in S.vecs]
    nr = len(R.vecs)
    out = []
    for c in _nullspace(cols):
        x = 0; z = 0; j = 0
        while c:
            if c & 1:
                if j < nr:
                    x ^= R.left(R.vecs[j])
                else:
                    z ^= S.right(S.vecs[j - nr])
            c >>= 1; j += 1
        out.append(x | (z << R.du))
    return Rel(R.du, S.dv, out)


def window_relation(dims, maps, a, b, budget=400):
    if a == b:
        return Rel(dims[a], dims[a], [(1 << i) | (1 << (i + dims[a])) for i in range(dims[a])])
    cuts = [a]; acc = 0
    for k in range(a, b):
        acc += dims[k]
        if acc >= budget or k + 1 == b:
            cuts.append(k + 1); acc = 0
    R = None
    for u, v in zip(cuts, cuts[1:]):
        C = sections_relation(dims, maps, u, v)
        R = C if R is None else glue(R, C)
    return R


def rel_converse(R):
    return Rel(R.dv, R.du, [R.right(v) | (R.left(v) << R.dv) for v in R.vecs])


def _rank(vecs):
    piv = set(); red = {}
    for v in vecs:
        while v:
            p = v.bit_length() - 1
            if p in red:
                v ^= red[p]
            else:
                red[p] = v; break
    return len(red)


def rel_stats(R):
    lp = _rank([R.left(v) for v in R.vecs]); rp = _rank([R.right(v) for v in R.vecs])
    return len(R.vecs), lp, rp, lp + rp - len(R.vecs)


def rel_fixed_dim(R):
    n = R.du
    diag = [(1 << i) | (1 << (i + n)) for i in range(n)]
    return len(R.vecs) + len(diag) - _rank(R.vecs + diag)


# ---------------------------------------------------------------------------
# Kronecker invariants via kernels of the pencil  P(s) = beta - s*alpha
# (alpha = left projection, beta = right projection of a basis of R),
# acting on (F2[s]/(q))^r  ->  (F2[s]/(q))^n  and on polynomial vectors of
# bounded degree.  This is the classical rank characterisation of the Kronecker
# canonical form of a pencil; it never uses subspace chains of the relation.
# ---------------------------------------------------------------------------

def _pmul(a, b):
    r = 0
    while b:
        if b & 1:
            r ^= a
        a <<= 1; b >>= 1
    return r


def _pmod(a, m):
    dm = m.bit_length() - 1
    while a and a.bit_length() - 1 >= dm:
        a ^= m << (a.bit_length() - 1 - dm)
    return a


def _spread(v, stride, k):
    out = 0; u = 0
    while v:
        if v & 1:
            out |= 1 << (u * stride + k)
        v >>= 1; u += 1
    return out


def _kernel_dim_mod(alpha, beta, n, q):
    dq = q.bit_length() - 1
    cols = []
    for x, y in zip(alpha, beta):
        for l in range(dq):
            c = _spread(y, dq, l)
            t = _pmod(1 << (l + 1), q)
            k = 0
            while t:
                if t & 1:
                    c ^= _spread(x, dq, k)
                t >>= 1; k += 1
            cols.append(c)
    return len(cols) - _rank(cols)


def _poly_kernel_dim(alpha, beta, n, k):
    """dim{(x_0..x_k): beta x_0 = 0, beta x_i = alpha x_{i-1}, alpha x_k = 0}."""
    cols = []
    for i in range(k + 1):
        for x, y in zip(alpha, beta):
            # contributes beta*x at block i and alpha*x at block i+1
            cols.append((y << (i * n)) ^ (x << ((i + 1) * n)))
    return len(cols) - _rank(cols)


def _minimal_indices(alpha, beta, n, mass_bound):
    """Multiset of minimal indices from polynomial-kernel growth.

    N_k = sum_{e <= k} (k - e + 1) m_e, so N_k - N_{k-1} = #{e <= k}.
    A block of index e contributes e to rank(alpha) and to rank(beta), so the
    total index mass of blocks not yet found is at most mass_bound minus the
    mass already found; once that is below k+1 no block of index > k remains."""
    Nk = []; out = []
    k = 0
    while True:
        Nk.append(_poly_kernel_dim(alpha, beta, n, k))
        g = Nk[k] - (Nk[k - 1] if k else 0)          # #{e <= k}
        gprev = (Nk[k - 1] - (Nk[k - 2] if k >= 2 else 0)) if k >= 1 else 0
        out.extend([k] * (g - gprev))
        if mass_bound - sum(out) < k + 1:
            return out
        k += 1


def _irreducibles(max_deg):
    for d in range(1, max_deg + 1):
        for low in range(1 << d):
            p = (1 << d) | low
            if not (p & 1) and p != 2:
                continue
            ok = True
            for e in range(1, d // 2 + 1):
                for low2 in range(1 << e):
                    f = (1 << e) | low2
                    if _pmod(p, f) == 0:
                        ok = False; break
                if not ok:
                    break
            if ok:
                yield p


def kronecker_invariants_pencil(R):
    n = R.du; assert R.dv == n
    alpha = [R.left(v) for v in R.vecs]; beta = [R.right(v) for v in R.vecs]
    r = len(R.vecs)
    mass = min(_rank(alpha), _rank(beta))
    right = _minimal_indices(alpha, beta, n, mass)
    # left indices: transposed pencil  beta^T - s alpha^T : F2^n -> F2^r
    def transpose(cols, rows):
        out = []
        for u in range(rows):
            v = 0
            for j, c in enumerate(cols):
                if (c >> u) & 1:
                    v |= 1 << j
            out.append(v)
        return out
    aT = transpose(alpha, n); bT = transpose(beta, n)
    left = _minimal_indices(aT, bT, r, mass - sum(right)) if n else []
    nright = len(right)
    # infinite divisors: reversed pencil alpha - s*beta at s = 0
    inf = []
    budget = n - sum(right) - sum(e + 1 for e in left)
    prev = 0; ge = []
    j = 1
    while True:
        kj = _kernel_dim_mod(beta, alpha, n, 1 << j) - j * nright
        ge.append(kj - prev); prev = kj
        if ge[-1] == 0:
            ge.pop(); break
        j += 1
    for t in range(1, len(ge) + 1):
        inf.extend([t] * (ge[t - 1] - (ge[t] if t < len(ge) else 0)))
    budget -= sum(inf)
    fin = []
    if budget:
        for p in _irreducibles(budget):
            d = p.bit_length() - 1
            prev = 0; ge = []; q = 1
            while True:
                q = _pmul(q, p)
                kj = _kernel_dim_mod(alpha, beta, n, q) - (q.bit_length() - 1) * nright
                assert kj % d == 0
                ge.append(kj // d - prev); prev = kj // d
                if ge[-1] == 0:
                    ge.pop(); break
            for t in range(1, len(ge) + 1):
                fin.extend([[p, t]] * (ge[t - 1] - (ge[t] if t < len(ge) else 0)))
            budget -= sum(d * t for pp, t in fin if pp == p)
            if budget == 0:
                break
            assert budget > 0
    assert r == sum(e + 1 for e in right) + sum(left) + sum(inf) + sum((p.bit_length() - 1) * t for p, t in fin)
    return {'right_minimal_indices': sorted(right), 'left_minimal_indices': sorted(left),
            'infinite_elementary_divisors': sorted(inf), 'finite_elementary_divisors': sorted(fin)}


_KR_CACHE = {}

def kronecker_cached(R):
    key = (R.du, tuple(R.vecs))
    if key not in _KR_CACHE:
        _KR_CACHE[key] = kronecker_invariants_pencil(R)
    return _KR_CACHE[key]


# ---------------------------------------------------------------------------
# Generalized ranks and barcode: rho(a,b) = dim ran R_{a,b} - dim mul R_{a,b},
# computed by a per-start sweep; bars by inclusion-exclusion.
# ---------------------------------------------------------------------------

def _push(direction, cols, dsrc_next, S):
    """Image of subspace S of V_k under the edge relation k -> k+1."""
    if direction == 'forward':
        return [_apply(cols, v) for v in S]
    # backward: g : V_{k+1} -> V_k ; image = g^{-1}(S)
    piv = {}
    for v in S:
        while v:
            p = v.bit_length() - 1
            if p in piv:
                v ^= piv[p][0]
            else:
                piv[p] = (v, 0); break
    out = []
    for j in range(dsrc_next):
        v = cols[j]; lab = 1 << j
        while v:
            p = v.bit_length() - 1
            if p in piv:
                a, l = piv[p]; v ^= a; lab ^= l
            else:
                piv[p] = (v, lab); break
        if v == 0:
            out.append(lab)
    return out


def _apply(cols, v):
    r = 0; j = 0
    while v:
        if v & 1:
            r ^= cols[j]
        v >>= 1; j += 1
    return r


def rank_rows(dims, maps):
    z = len(dims); rho = {}
    for a in range(z):
        E = [1 << i for i in range(dims[a])]; N = []
        r = dims[a]; b = a
        while r > 0:
            rho[(a, b)] = r
            if b + 1 >= z:
                break
            direction, cols = maps[b]
            E = _push(direction, cols, dims[b + 1], E); N = _push(direction, cols, dims[b + 1], N)
            b += 1
            r = _rank(E) - _rank(N)
    return rho


def barcode_ind(dims, maps):
    rho = rank_rows(dims, maps)
    R = lambda a, b: rho.get((a, b), 0)
    bars = {}
    for (a, b) in rho:
        for (x, y) in ((a, b), (a + 1, b), (a, b - 1), (a + 1, b - 1)):
            if x <= y and (x, y) not in bars:
                m = R(x, y) - R(x - 1, y) - R(x, y + 1) + R(x - 1, y + 1)
                bars[(x, y)] = m
    out = []
    for (x, y), m in sorted(bars.items()):
        if m < 0:
            raise RuntimeError('negative multiplicity')
        if m:
            out.append([x, y, m])
    return out


# ---------------------------------------------------------------------------
# Full barcode by a RIGHT-TO-LEFT sweep over different flags:
#   D_e = dom R_{k,e}  and  K_e = ker R_{k,e}  inside V_k,  for every end e >= k,
# so rho(k,e) = dim D_e - dim K_e.  Moving from k+1 to k pulls every subspace back
# through the edge relation S_k between V_k and V_{k+1}:  X -> S_k^{-1}(X), i.e.
# f^{-1}(X) for a forward map f and g(X) for a backward map g.
# ---------------------------------------------------------------------------

def _rref(vecs):
    piv = {}
    for v in vecs:
        while v:
            p = v.bit_length() - 1
            if p in piv:
                v ^= piv[p]
            else:
                piv[p] = v; break
    keys = sorted(piv)
    for i, p in enumerate(keys):
        for q in keys[i + 1:]:
            if (piv[q] >> p) & 1:
                piv[q] ^= piv[p]
    return tuple(piv[p] for p in keys)


def _pull_forward(cols, dsrc, X):
    """f^{-1}(X) for f: F2^dsrc -> target given by columns."""
    piv = {}
    for x in X:
        v = x
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
                a, l = piv[p]; v ^= a; lab ^= l
            else:
                piv[p] = (v, lab); break
        if v == 0:
            out.append(lab)
    return _rref(out)


def _push_backward(cols, X):
    """g(X) for g given by columns."""
    return _rref([_apply(cols, x) for x in X])


def barcode_reverse(dims, maps):
    n = len(dims)
    if n == 0:
        return []
    full = lambda d: _rref([1 << i for i in range(d)])
    # segs over e (ascending): [e_lo, D, K], valid for e in [e_lo, next e_lo)
    segs = [[n - 1, full(dims[n - 1]), ()]]
    bars = {}

    def jumps(sg, lo):
        """rho(k,e) - rho(k,e+1) for e >= lo, from a segment list over e."""
        out = {}
        vals = [(e, len(D) - len(K)) for e, D, K in sg if True]
        # value at e is that of the last segment with e_lo <= e; rho(k, n) = 0
        for i, (e0, r) in enumerate(vals):
            e_last = (vals[i + 1][0] - 1) if i + 1 < len(vals) else n - 1
            nxt = vals[i + 1][1] if i + 1 < len(vals) else 0
            if e_last >= lo and r != nxt:
                out[e_last] = out.get(e_last, 0) + r - nxt
        return out

    for k in range(n - 1, -1, -1):
        if k > 0:
            direction, cols = maps[k - 1]
            new = [[k - 1, full(dims[k - 1]), ()]]
            for e, D, K in segs:
                if direction == 'forward':
                    D2 = _pull_forward(cols, dims[k - 1], D); K2 = _pull_forward(cols, dims[k - 1], K)
                else:
                    D2 = _push_backward(cols, D); K2 = _push_backward(cols, K)
                if new[-1][1] == D2 and new[-1][2] == K2:
                    continue
                new.append([e, D2, K2])
        else:
            new = []
        # m(k,e) = [rho(k,e)-rho(k,e+1)] - [rho(k-1,e)-rho(k-1,e+1)] for e >= k
        jk = jumps(segs, k); jn = jumps(new, k) if new else {}
        for e in set(jk) | set(jn):
            m = jk.get(e, 0) - jn.get(e, 0)
            if m < 0:
                raise RuntimeError('negative multiplicity (reverse sweep)')
            if m:
                bars[(k, e)] = m
        segs = new
    return [[a, b, m] for (a, b), m in sorted(bars.items())]


# ---------------------------------------------------------------------------
# Conley index maps (independent route): a different chain selector (top-right
# vertex choice, vertical-first edge paths, top-down square fill) and elementary
# divisors from the Kronecker pencil of the graph relation {(x, Ix)}.
# ---------------------------------------------------------------------------

def _pair_for_interval(frame, recurrent, S):
    off = frame['offsets']; tg = frame['targets']
    core = set()
    for m in S:
        core |= set(recurrent[m])
    fwd = set(core); todo = list(core)
    while todo:
        v = todo.pop()
        for w in tg[off[v]:off[v + 1]]:
            if w not in fwd:
                fwd.add(w); todo.append(w)
    back = {}
    for v in fwd:
        for w in tg[off[v]:off[v + 1]]:
            back.setdefault(w, set()).add(v)
    can = set(core); todo = list(core)
    while todo:
        w = todo.pop()
        for v in back.get(w, ()):
            if v not in can:
                can.add(v); todo.append(v)
    return fwd, {b for b in fwd if b not in can}


def _rects(frame, nx):
    off = frame['offsets']; tg = frame['targets']; out = []
    for q in range(len(off) - 1):
        xs = [t % nx for t in tg[off[q]:off[q + 1]]]; ys = [t // nx for t in tg[off[q]:off[q + 1]]]
        out.append((min(xs), max(xs) + 1, min(ys), max(ys) + 1))
    return out


def _carrier(cell, dim, P, rects, nx):
    if dim == 0:
        i, j = cell; boxes = [(i - 1, j - 1), (i, j - 1), (i - 1, j), (i, j)]
    elif dim == 1:
        a, i, j = cell; boxes = [(i, j - 1), (i, j)] if a == 0 else [(i - 1, j), (i, j)]
    else:
        boxes = [cell]
    lo_x = lo_y = -1; hi_x = hi_y = 1 << 30; seen = False
    for i, j in boxes:
        if 0 <= i < nx and j >= 0 and (i + nx * j) in P:
            r = rects[i + nx * j]; seen = True
            lo_x = max(lo_x, r[0]); hi_x = min(hi_x, r[1]); lo_y = max(lo_y, r[2]); hi_y = min(hi_y, r[3])
    if not seen or lo_x > hi_x or lo_y > hi_y:
        raise RuntimeError('empty carrier')
    return lo_x, hi_x, lo_y, hi_y


def _index_map_columns(P, E, frame, nx):
    rects = _rects(frame, nx)
    Hd = get_relative_homology_data_boxes(P, E, nx)
    vmemo = {}; ememo = {}; smemo = {}
    def phi_v(v):
        if v not in vmemo:
            c = _carrier(v, 0, P, rects, nx); vmemo[v] = (c[1], c[3])
        return vmemo[v]
    def phi_e(e):
        if e not in ememo:
            a, i, j = e
            u = phi_v((i, j)); w = phi_v((i + 1, j) if a == 0 else (i, j + 1))
            ch = set()
            y0, y1 = sorted((u[1], w[1]))
            for y in range(y0, y1):
                ch ^= {(1, u[0], y)}
            x0, x1 = sorted((u[0], w[0]))
            for x in range(x0, x1):
                ch ^= {(0, x, w[1])}
            ememo[e] = ch
        return ememo[e]
    def phi_s(sq):
        if sq not in smemo:
            i, j = sq; z = set()
            for e in ((0, i, j), (0, i, j + 1), (1, i, j), (1, i + 1, j)):
                z ^= phi_e(e)
            lx, hx, ly, hy = _carrier(sq, 2, P, rects, nx)
            out = set()
            if hx > lx and hy > ly:
                for x in range(lx, hx):
                    inside = False
                    for y in range(hy - 1, ly - 1, -1):
                        if (0, x, y + 1) in z:
                            inside = not inside
                        if inside:
                            out.add((x, y))
            chk = set()
            for x, y in out:
                chk ^= {(0, x, y), (0, x, y + 1), (1, x, y), (1, x + 1, y)}
            if chk != z:
                raise RuntimeError('square fill failed')
            smemo[sq] = out
        return smemo[sq]
    img = [lambda c: {phi_v(c)}, phi_e, phi_s]
    res = []
    for d in range(3):
        H = Hd[d]; basis = H['basis']; pos = {c: t for t, c in enumerate(basis)}; nb = len(H['B'])
        cols = []
        for h in H['H']:
            y = 0; t = 0; x = h
            while x:
                if x & 1:
                    for c in img[d](basis[t]):
                        u = pos.get(c)
                        if u is not None:
                            y ^= 1 << u
                x >>= 1; t += 1
            rem, lab = H['coord'].reduce(y)
            if rem:
                raise RuntimeError('selector image not a relative cycle')
            cols.append(lab >> nb)
        res.append(cols)
    return res


def index_map_record(frame, recurrent, S, nx):
    P, E = _pair_for_interval(frame, recurrent, S)
    out = []
    for d, cols in enumerate(_index_map_columns(P, E, frame, nx)):
        n = len(cols)
        R = Rel(n, n, [(1 << i) | (c << n) for i, c in enumerate(cols)])
        fin = kronecker_cached(R)['finite_elementary_divisors'] if n else []
        out.append({'dimension': d, 'space_dimension_F2': n, 'rank_F2': _rank(cols),
                    'elementary_divisors_F2': fin})
    return out


# ---------------------------------------------------------------------------

WORDS = ['A', 'B', 'AB', 'BA', 'ABA', 'BAB', 'AABB', 'BBAA', 'ABAB', 'BABA', 'Ab', 'aB', 'ABab', 'abAB']


def loop_record(R, with_powers):
    dim, lp, rp, rr = rel_stats(R)
    rec = {'dimension': None, 'relation_dimension_F2': dim, 'left_projection_rank_F2': lp,
           'right_projection_rank_F2': rp, 'regular_rank_F2': rr,
           'fixed_subspace_dimension_F2': rel_fixed_dim(R)}
    if with_powers:
        pw = [rr]; P = R
        for _ in range(3):
            P = glue(P, P); pw.append(rel_stats(P)[3])
        rec['power_regular_ranks_F2'] = pw
    rec['kronecker_invariants_F2'] = kronecker_cached(R)
    return rec


def holonomy_record(probe, dims_by_d, maps_by_d):
    a, b, c = probe; out = []
    for d in range(len(dims_by_d)):
        dims = dims_by_d[d]; maps = maps_by_d[d]; n = dims[2 * a]
        A = window_relation(dims, maps, 2 * a, 2 * b); B = window_relation(dims, maps, 2 * b, 2 * c)
        L = {'A': A, 'B': B, 'a': rel_converse(A), 'b': rel_converse(B)}
        sig = []
        for w in WORDS:
            R = L[w[0]]
            for ch in w[1:]:
                R = glue(R, L[ch])
            dim, lp, rp, rr = rel_stats(R)
            sig.append({'word': w, 'relation_dimension_F2': dim, 'left_projection_rank_F2': lp,
                        'right_projection_rank_F2': rp, 'regular_rank_F2': rr,
                        'fixed_subspace_dimension_F2': rel_fixed_dim(R),
                        'kronecker_invariants_F2': kronecker_cached(R)})
        out.append({'dimension': d, 'space_dimension_F2': n, 'word_signatures': sig})
    return out


def solve_case_ind(data):
    nx = data['phase_grid']['nx']; ny = data['phase_grid']['ny']; n = nx * ny
    frames = build_frames_int(data)
    analyses = [analyze_frame_ind(f, nx, ny) for f in frames]
    edge_rel = {}; edge_struct = {}
    for p, q in data['parameter_edges']:
        rel, U = union_edge_structure(frames[p], analyses[p], frames[q], analyses[q], n)
        edge_rel[(p, q)] = rel; edge_struct[(p, q)] = U
    queries = []
    import os
    only = os.environ.get('IND_QUERIES')
    only = None if not only else {int(x) for x in only.split(',')}
    for q in data['queries']:
        if only is not None and q['query_id'] not in only:
            continue
        ids = q['parameter_ids']
        sel = propagate_selection(ids, q['seed_box_id'], analyses, edge_rel)
        dims, maps = conley_zigzag_ind(ids, sel, analyses, nx)
        ranks = []
        for k, (a, b) in enumerate(q['rank_windows']):
            rk = []
            for d in range(3):
                R = window_relation(dims[d], maps[d], a, b)
                rk.append(rel_stats(R)[3])
            ranks.append({'window_id': k, 'start_node': a, 'end_node': b, 'rank_F2': rk})
        bars = [barcode_ind(dims[d], maps[d]) for d in range(3)]
        loops = []
        for lid, (a, b) in enumerate(q['loop_windows']):
            degs = []
            for d in range(3):
                rec = loop_record(window_relation(dims[d], maps[d], 2 * a, 2 * b), True)
                rec['dimension'] = d; degs.append(rec)
            loops.append({'loop_id': lid, 'start_occurrence': a, 'end_occurrence': b, 'degrees': degs})
        hol = [{'probe_id': h, 'occurrences': list(pr), 'degrees': holonomy_record(pr, dims, maps)}
               for h, pr in enumerate(q['holonomy_probes'])]
        gdims, gmaps = graph_zigzag_ind(ids, analyses, edge_struct)
        gbars = [barcode_reverse(gdims[d], gmaps[d]) for d in range(2)]
        granks = []
        for k, (a, b) in enumerate(q['graph_rank_windows']):
            granks.append({'window_id': k, 'start_node': a, 'end_node': b,
                           'rank_F2': [rel_stats(window_relation(gdims[d], gmaps[d], a, b))[3] for d in range(2)]})
        gloops = []
        for lid, (a, b) in enumerate(q['graph_loop_windows']):
            degs = []
            for d in range(2):
                rec = loop_record(window_relation(gdims[d], gmaps[d], 2 * a, 2 * b), False)
                rec['dimension'] = d; degs.append(rec)
            gloops.append({'loop_id': lid, 'start_occurrence': a, 'end_occurrence': b, 'degrees': degs})
        ghol = [{'probe_id': h, 'occurrences': list(pr), 'degrees': holonomy_record(pr, gdims, gmaps)}
                for h, pr in enumerate(q['graph_holonomy_probes'])]
        imaps = []
        for h, (k, S) in enumerate(q.get('index_probes', [])):
            par = ids[k]; S = sorted(set(S))
            imaps.append({'probe_id': h, 'occurrence': k, 'morse_ids': S,
                          'degrees': index_map_record(frames[par], analyses[par]['recurrent'], S, nx)})
        queries.append({'query_id': q['query_id'], 'parameter_ids': ids, 'seed_box_id': q['seed_box_id'],
                        'selected_morse_ids_by_frame': sel, 'conley_node_dimensions_F2': dims,
                        'conley_generalized_rank_queries': ranks, 'conley_zigzag_barcodes_F2': bars,
                        'conley_loop_signatures': loops, 'conley_holonomy_word_signatures': hol,
                        'morse_graph_node_dimensions_F2': gdims, 'morse_graph_generalized_rank_queries': granks,
                        'morse_graph_zigzag_barcodes_F2': gbars,
                        'morse_graph_loop_signatures': gloops, 'morse_graph_holonomy_word_signatures': ghol,
                        'conley_index_maps': imaps})
        cached_relative_homology_data_boxes.cache_clear()
    return {'case_id': data['case_id'], 'queries': queries}


def main():
    import argparse
    ap = argparse.ArgumentParser(); ap.add_argument('input'); ap.add_argument('output')
    a = ap.parse_args()
    with open(a.input) as f:
        data = json.load(f)
    with open(a.output, 'w') as f:
        json.dump(solve_case_ind(data), f, separators=(',', ':'))


if __name__ == '__main__':
    main()

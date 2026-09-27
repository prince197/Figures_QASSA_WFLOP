#!/usr/bin/env python3
from pathlib import Path
import importlib.util, json, random, time
from collections import Counter, defaultdict

ROOT = Path(__file__).resolve().parents[2]
SOLVER = ROOT / 'solution' / 'trajectory_solver.py'
spec = importlib.util.spec_from_file_location('trajectory_ref', SOLVER)
ref = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ref)

PUBLIC = ROOT / 'environment' / 'data' / 'example.json'
PUBLIC_EXPECTED = ROOT / 'environment' / 'data' / 'example_expected.json'
HIDDEN_DIR = ROOT / 'tests' / 'heldout'
ORACLE_DIR = ROOT / 'tests' / 'oracles'

CONFIGS = [
    ('public', 33, 14, 12, 1000, 2880, 27, 2920, 26, 26, 1000, 111),
    ('heldout_0', 36, 15, 15, 1000, 2890, 25, 2910, 24, 27, 1000, 202),
    ('heldout_1', 35, 14, 16, 1000, 2875, 27, 2925, 23, 24, 1000, 303),
    ('heldout_2', 34, 16, 14, 1000, 2885, 26, 2920, 25, 27, 1000, 404),
    ('heldout_3', 35, 15, 14, 1000, 2885, 26, 2915, 25, 25, 1000, 505),
]

LENGTHS = {
    'public': (1200, 1600),
    'heldout_0': (3000, 3900, 4800),
    'heldout_1': (2800, 3700, 4500),
    'heldout_2': (2900, 3800, 4600),
    'heldout_3': (2800, 3600, 4400),
}


def lattice_edges(na, nb):
    out = []
    for j in range(nb):
        for i in range(na):
            p = i + na * j
            if i + 1 < na:
                out.append([p, p + 1])
            if j + 1 < nb:
                out.append([p, p + na])
    return out


def make_base(name, nx, na, nb, scale, r0, dr, s0, ds, cnum, cden):
    return {
        'case_id': name,
        'model': {
            'family': 'coupled_logistic',
            'formula': [
                '(1-c)*r*x*(1-x)+c*s*y*(1-y)',
                'c*r*x*(1-x)+(1-c)*s*y*(1-y)',
            ],
            'phase_bounds': [0, 1, 0, 1],
            'parameter_scale': scale,
            'coupling_num': cnum,
            'coupling_den': cden,
        },
        'phase_grid': {'nx': nx, 'ny': nx},
        'parameter_grid': {
            'na': na,
            'nb': nb,
            'r_num': [r0 + dr * i for i in range(na)],
            's_num': [s0 + ds * j for j in range(nb)],
        },
        'parameter_edges': lattice_edges(na, nb),
        'queries': [],
    }


def analyze_base(data):
    nx = data['phase_grid']['nx']
    ny = data['phase_grid']['ny']
    n = nx * ny
    frames = ref.build_frames(data)
    analyses = [ref.analyze_frame(f, nx, ny) for f in frames]
    edge_rel = {}
    edge_struct = {}
    edge_score = {}
    adj = defaultdict(list)
    for p, q in data['parameter_edges']:
        rel, U = ref.union_edge_structure(frames[p], analyses[p], frames[q], analyses[q], n)
        edge_rel[(p, q)] = rel
        edge_struct[(p, q)] = U
        ca = Counter(i for i, j in rel)
        cb = Counter(j for i, j in rel)
        branch = sum(max(0, v - 1) for v in ca.values()) + sum(max(0, v - 1) for v in cb.values())
        score = (
            len(rel)
            + 7 * branch
            + 3 * sum(v > 1 for v in ca.values())
            + 3 * sum(v > 1 for v in cb.values())
        )
        edge_score[frozenset((p, q))] = score
        adj[p].append(q)
        adj[q].append(p)
    return frames, analyses, edge_rel, edge_struct, edge_score, adj


def coverage_walk(rng, adj, edge_score, start, length, na, nb):
    path = [start]
    prev = None
    cur = start
    evis = Counter()
    dvis = Counter()
    nvis = Counter([start])
    previous_direction = None
    for t in range(length - 1):
        opts = adj[cur]
        weights = []
        for v in opts:
            sc = edge_score[frozenset((cur, v))]
            w = 1.0 + (sc + 1.0) ** 2
            if prev is not None and v == prev:
                w *= 0.012
            undirected = tuple(sorted((cur, v)))
            w *= 1.0 + 5.0 / (1 + evis[undirected])
            w *= 1.0 + 3.0 / (1 + dvis[(cur, v)])
            w *= 1.0 + 1.8 / (1 + nvis[v])
            direction = ((v % na) - (cur % na), (v // na) - (cur // na))
            if previous_direction is not None:
                w *= 1.30 if direction != previous_direction else 0.90
            if t % 37 in (0, 1, 2):
                x, y = v % na, v // na
                boundary_distance = min(x, na - 1 - x, y, nb - 1 - y)
                w *= 1.0 + 0.45 / (1 + boundary_distance)
            weights.append(w)
        nxt = rng.choices(opts, weights=weights, k=1)[0]
        undirected = tuple(sorted((cur, nxt)))
        evis[undirected] += 1
        dvis[(cur, nxt)] += 1
        nvis[nxt] += 1
        previous_direction = ((nxt % na) - (cur % na), (nxt // na) - (cur // na))
        prev, cur = cur, nxt
        path.append(cur)
    return path


def propagate(path, analyses, edge_rel, start_mid):
    selected = [[start_mid]]
    total = 1
    peak = 1
    changes = 0
    for p, q in zip(path, path[1:]):
        rel = ref.relation_for_oriented_edge(p, q, edge_rel)
        prev = set(selected[-1])
        nxt = sorted({j for i, j in rel if i in prev})
        if not nxt:
            return None
        changes += int(nxt != selected[-1])
        selected.append(nxt)
        total += len(nxt)
        peak = max(peak, len(nxt))
    return selected, total, peak, changes


def choose_diverse_pairs(scored_pairs, cap, min_sep):
    out = []
    used_span_bins = Counter()
    for score, a, b in sorted(scored_pairs, reverse=True):
        span = b - a
        span_bin = span // max(25, min_sep)
        if used_span_bins[span_bin] >= 2:
            continue
        if any(abs(a - c) < min_sep and abs(b - d) < min_sep for c, d in out):
            continue
        out.append((a, b))
        used_span_bins[span_bin] += 1
        if len(out) >= cap:
            break
    return out


def _type_key(R, n):
    kt = ref.kronecker_type(R, n)
    return (n, tuple(kt['right_minimal_indices']), tuple(kt['left_minimal_indices']),
            tuple(kt['infinite_elementary_divisors']),
            tuple(tuple(x) for x in kt['finite_elementary_divisors']))


def _identity_relation(n):
    return [(1 << i) | ((1 << i) << n) for i in range(n)]


def window_profiles(dims_by_d, maps_by_d, windows):
    """For occurrence windows (a,b) return per-degree Kronecker keys and whether the
    running composition passes through the zero relation before ending nonzero."""
    by_start = defaultdict(set)
    for a, b in windows:
        by_start[a].add(b)
    out = {}
    for a, ends in by_start.items():
        info = {b: [] for b in ends}
        for d in range(len(dims_by_d)):
            dims = dims_by_d[d]; maps = maps_by_d[d]
            left = 2 * a; n = dims[left]
            R = _identity_relation(n); rd = n; hit_zero = False
            last = 2 * max(ends)
            for k in range(left, last):
                E = ref._edge_linear_relation(dims[k], dims[k + 1], maps[k][0], maps[k][1])
                R = ref._compose_linear_relations(R, n, rd, E, dims[k + 1]); rd = dims[k + 1]
                if not R:
                    hit_zero = True
                if (k + 1) % 2 == 0 and (k + 1) // 2 in info:
                    info[(k + 1) // 2].append((_type_key(R, n), bool(hit_zero and R)))
        for b in ends:
            out[(a, b)] = info[b]
    return out


def word_profiles(dims_by_d, maps_by_d, probe):
    a, b, c = probe
    keys = []; crossing = False
    for d in range(len(dims_by_d)):
        dims = dims_by_d[d]; maps = maps_by_d[d]; n = dims[2 * a]
        A = ref._relation_between_nodes(dims, maps, 2 * a, 2 * b)
        B = ref._relation_between_nodes(dims, maps, 2 * b, 2 * c)
        rels = {'A': A, 'B': B, 'a': ref.kr_converse(A, n), 'b': ref.kr_converse(B, n)}
        for w in ref.HOLONOMY_WORDS:
            R = _identity_relation(n); z = False
            for ch in w:
                R = ref._compose_linear_relations(R, n, n, rels[ch], n)
                if not R:
                    z = True
            crossing |= bool(z and R)
            keys.append((d, w, _type_key(R, n)))
    return keys, crossing


def select_by_structure(candidates, profile_of, cap, too_close, rarity_boost=6.0):
    """Greedy selection that first covers rare relation types and zero-crossing
    compositions, then falls back to the original score."""
    freq = Counter()
    prof = {}
    for sc, *w in candidates:
        keys, crossing = profile_of(tuple(w))
        prof[tuple(w)] = (keys, crossing)
        freq.update(set(keys))
    chosen = []; seen = set(); crossings = 0
    def gain(item):
        sc, *w = item
        keys, crossing = prof[tuple(w)]
        new = [k for k in set(keys) if k not in seen]
        g = sum(rarity_boost / freq[k] for k in new)
        if crossing and crossings < max(2, cap // 4):
            g += 50.0
        return g + 0.001 * sc
    pool = list(candidates)
    while pool and len(chosen) < cap:
        pool.sort(key=gain, reverse=True)
        item = pool.pop(0)
        w = tuple(item[1:])
        if any(too_close(w, x) for x in chosen):
            continue
        keys, crossing = prof[w]
        seen.update(keys); crossings += int(crossing)
        chosen.append(w)
    return chosen


def choose_query(qid, length, rng, frames, analyses, edge_rel, edge_struct, edge_score, adj, na, nb, nx, public=False):
    hard = sorted(edge_score.items(), key=lambda kv: kv[1], reverse=True)
    start_candidates = []
    for ek, _ in hard[:max(24, len(hard) // 4)]:
        start_candidates.extend(tuple(ek))

    best = None
    attempts = 6 if public else 10
    for _ in range(attempts):
        start = rng.choice(start_candidates)
        path = coverage_walk(rng, adj, edge_score, start, length, na, nb)
        unique_edges = len(set(tuple(sorted((a, b))) for a, b in zip(path, path[1:])))
        unique_parameters = len(set(path))
        reversals = sum(path[i + 1] == path[i - 1] for i in range(1, len(path) - 1))
        for mid, comp in enumerate(analyses[start]['recurrent']):
            pr = propagate(path, analyses, edge_rel, mid)
            if pr is None:
                continue
            selected, total, peak, changes = pr
            by_state = defaultdict(list)
            for k, (pid, s) in enumerate(zip(path, selected)):
                by_state[(pid, tuple(s))].append(k)
            long_repeats = sum(
                1 for inds in by_state.values()
                if len(inds) > 1 and inds[-1] - inds[0] >= length // 3
            )
            score = (
                2.5 * unique_edges
                + 1.5 * unique_parameters
                + 0.18 * changes
                + 0.025 * total
                + 7 * peak
                + 4 * long_repeats
                - 20 * reversals
            )
            candidate = (score, path, min(comp), selected)
            if best is None or score > best[0]:
                best = candidate
    if best is None:
        raise RuntimeError('failed to build persistent coverage trajectory')

    _, path, seed, selected = best
    q = {'query_id': qid, 'parameter_ids': path, 'seed_box_id': seed}
    sel, dims, maps = ref.trajectory_core(q, frames, analyses, edge_rel, nx)
    if sel != selected:
        raise RuntimeError('internal propagation mismatch')

    by_state = defaultdict(list)
    for k, (pid, s) in enumerate(zip(path, sel)):
        by_state[(pid, tuple(s))].append(k)
    loop_candidates = []
    for inds in by_state.values():
        if len(inds) < 2:
            continue
        stride = max(1, len(inds) // 14)
        sample = sorted(set(inds[::stride] + [inds[-1]]))
        for ii, a in enumerate(sample):
            for b in sample[ii + 1:]:
                span = b - a
                if span < max(40, length // 18):
                    continue
                h1 = dims[1][2 * a]
                loop_candidates.append((2.5 * h1 + 0.004 * span + 0.2 * len(sel[a]), a, b))
    loop_cap = 8 if public else 20
    # Keep a broad shortlist, then choose windows that exercise rare continuation
    # relation types (zero, multivalued, kernel, full, graph) and compositions that
    # pass through the zero relation before becoming nonzero again.
    shortlist = choose_diverse_pairs(loop_candidates, 6 * loop_cap, max(6, length // 300))
    extra = sorted(loop_candidates, key=lambda t: rng.random())[:4 * loop_cap]
    pool = {}
    for item in [(sc, a, b) for sc, a, b in loop_candidates if (a, b) in set(shortlist)] + extra:
        pool[(item[1], item[2])] = item
    pool = list(pool.values())
    prof = window_profiles(dims, maps, [(a, b) for _, a, b in pool])
    sep = max(12, length // 100)
    loops = select_by_structure(
        pool,
        lambda w: ([(d, k) for d, (k, _) in enumerate(prof[w])], any(c for _, c in prof[w])),
        loop_cap,
        lambda w, x: abs(w[0] - x[0]) < sep and abs(w[1] - x[1]) < sep,
    )
    if len(loops) < max(4, loop_cap // 2):
        raise RuntimeError('not enough Conley loop windows')
    q['loop_windows'] = [list(x) for x in loops]

    # Continuation holonomy probes use three returns to exactly the same
    # parameter and propagated Morse selection.  The two consecutive loops
    # generate an endorelation semigroup whose word signatures retain
    # information that ordinary loop ranks do not.
    holonomy_candidates = []
    min_gap = max(55, length // 45)
    for inds in by_state.values():
        if len(inds) < 3:
            continue
        if len(inds) > 28:
            sample = [inds[round(k * (len(inds)-1) / 27)] for k in range(28)]
            sample = sorted(set(sample))
        else:
            sample = inds
        for j in range(1, len(sample)-1):
            # Pair a few separated predecessors and successors around each middle return.
            lefts = sample[max(0, j-5):j]
            rights = sample[j+1:min(len(sample), j+6)]
            for a in lefts:
                for c in rights:
                    b = sample[j]
                    if b-a < min_gap or c-b < min_gap:
                        continue
                    h1 = dims[1][2*a]
                    score = 3.0*h1 + 0.003*(c-a) + 0.4*len(sel[a]) + 0.001*abs((b-a)-(c-b))
                    holonomy_candidates.append((score,a,b,c))
    holonomy_cap = 5 if public else 10
    hpool = []
    for item in sorted(holonomy_candidates, reverse=True):
        _, a, b, c = item
        if any(abs(a-x)<min_gap//3 and abs(b-y)<min_gap//3 and abs(c-z)<min_gap//3 for _,x,y,z in hpool):
            continue
        hpool.append(item)
        if len(hpool) >= 8 * holonomy_cap:
            break
    extra = sorted(holonomy_candidates, key=lambda t: rng.random())[:4 * holonomy_cap]
    hpool = list({(t[1], t[2], t[3]): t for t in hpool + extra}.values())
    hprof = {tuple(t[1:]): word_profiles(dims, maps, tuple(t[1:])) for t in hpool}
    holonomy = select_by_structure(
        hpool, lambda w: hprof[w], holonomy_cap,
        lambda w, x: abs(w[0]-x[0])<min_gap//2 and abs(w[1]-x[1])<min_gap//2 and abs(w[2]-x[2])<min_gap//2,
    )
    if len(holonomy) < max(2, holonomy_cap//2):
        raise RuntimeError('not enough Conley holonomy probes')
    q['holonomy_probes'] = [list(x) for x in holonomy]

    z = 2 * length - 1
    conley_candidates = set()
    starts = sorted(set([0, z - 1] + [rng.randrange(z) for _ in range(54 if public else 96)]))
    for a in starts:
        for span in (0, 1, 2, 3, 7, 15, 31, 63, 127, 255, 511, 1023, 1535):
            conley_candidates.add((a, min(z - 1, a + span)))
        for frac in (0.17, 0.33, 0.51, 0.68, 0.84, 0.97):
            conley_candidates.add((a, a + int((z - 1 - a) * frac)))
    for _ in range(140 if public else 360):
        a = rng.randrange(z)
        remain = z - 1 - a
        span = rng.randint(0, remain) if remain else 0
        conley_candidates.add((a, a + span))
    conley_candidates = sorted(conley_candidates)
    vals = [ref.ranks_for_windows(dims[d], maps[d], conley_candidates) for d in range(3)]
    nonzero = []
    zero = []
    for k, w in enumerate(conley_candidates):
        rank = tuple(vals[d][k] for d in range(3))
        item = (sum(rank) + 0.0015 * (w[1] - w[0]) + 0.08 * len(set(rank)), w)
        (nonzero if any(rank) else zero).append(item)
    nonzero.sort(reverse=True)
    zero.sort(reverse=True)
    rank_cap = 160 if public else 384
    target_nz = 120 if public else 300
    target_z = rank_cap - target_nz
    chosen = [w for _, w in nonzero[:target_nz]] + [w for _, w in zero[:target_z]]
    seen = set(chosen)
    for _, w in nonzero[target_nz:] + zero[target_z:]:
        if w not in seen:
            chosen.append(w)
            seen.add(w)
        if len(chosen) >= rank_cap:
            break
    q['rank_windows'] = [list(w) for w in sorted(chosen[:rank_cap])]

    gdims, gmaps = ref.graph_trajectory_core(path, analyses, edge_struct)
    gh1 = gdims[1]
    graph_candidates = set()
    point_count = 14 if public else 32
    for i in sorted(range(z), key=lambda k: (gh1[k], -k), reverse=True)[:point_count]:
        graph_candidates.add((i, i))
    spread_n = 12 if public else 24
    spread = [int((z - 1) * k / (spread_n - 1)) for k in range(spread_n)]
    for a in spread:
        for span in (1, 2, 5, 11, 23, 47):
            graph_candidates.add((a, min(z - 1, a + span)))
        for frac in (0.12, 0.27, 0.46, 0.69, 0.91):
            graph_candidates.add((a, a + int((z - 1 - a) * frac)))
    for _ in range(70 if public else 170):
        a = rng.randrange(z)
        remain = z - 1 - a
        if remain == 0:
            graph_candidates.add((a, a))
            continue
        if rng.random() < 0.75:
            lo = min(remain, max(1, z // 12))
            span = rng.randint(lo, remain) if remain >= lo else remain
        else:
            span = rng.randint(0, min(remain, 64))
        graph_candidates.add((a, a + span))
    graph_candidates = sorted(graph_candidates)
    gvals = [ref.ranks_for_windows(gdims[d], gmaps[d], graph_candidates) for d in range(2)]
    gnonzero = []
    gzero = []
    for k, w in enumerate(graph_candidates):
        rank = tuple(gvals[d][k] for d in range(2))
        item = (sum(rank) + 0.001 * (w[1] - w[0]) + 0.03 * max(rank), w)
        (gnonzero if any(rank) else gzero).append(item)
    gnonzero.sort(reverse=True)
    gzero.sort(reverse=True)
    graph_cap = 40 if public else 128
    graph_nz = 30 if public else 96
    graph_z = graph_cap - graph_nz
    graph_windows = [w for _, w in gnonzero[:graph_nz]] + [w for _, w in gzero[:graph_z]]
    seen = set(graph_windows)
    for _, w in gnonzero[graph_nz:] + gzero[graph_z:]:
        if w not in seen:
            graph_windows.append(w)
            seen.add(w)
        if len(graph_windows) >= graph_cap:
            break
    q['graph_rank_windows'] = [list(w) for w in sorted(graph_windows[:graph_cap])]

    by_parameter = defaultdict(list)
    for i, pid in enumerate(path):
        by_parameter[pid].append(i)
    graph_loop_candidates = []
    min_span = max(120, length // 14)
    max_span = min(1800, max(min_span + 1, length // 2))
    for inds in by_parameter.values():
        for ii, a in enumerate(inds):
            for b in inds[ii + 1:]:
                span = b - a
                if not (min_span <= span <= max_span):
                    continue
                graph_loop_candidates.append((0.04 * gh1[2 * a] + 0.0035 * span, a, b))
    graph_loop_cap = 6 if public else 10
    by_start_dim = defaultdict(list)
    for score, a, b in graph_loop_candidates:
        by_start_dim[gh1[2 * a]].append((score, a, b))
    available_dims = sorted(by_start_dim)
    graph_loops = []
    if available_dims:
        if graph_loop_cap == 1:
            targets = [available_dims[-1]]
        else:
            targets = [
                available_dims[round((len(available_dims) - 1) * k / (graph_loop_cap - 1))]
                for k in range(graph_loop_cap)
            ]
        for dim in reversed(targets):
            for _, a, b in sorted(by_start_dim[dim], reverse=True):
                if all(abs(a - c) >= max(14, length // 110) or abs(b - d) >= max(14, length // 110) for c, d in graph_loops):
                    graph_loops.append((a, b))
                    break
        if len(graph_loops) < graph_loop_cap:
            leftovers = sorted(graph_loop_candidates, reverse=True)
            for _, a, b in leftovers:
                if (a, b) in graph_loops:
                    continue
                graph_loops.append((a, b))
                if len(graph_loops) >= graph_loop_cap:
                    break
    if len(graph_loops) < max(3, graph_loop_cap // 2):
        raise RuntimeError('not enough graph loop windows')
    q['graph_loop_windows'] = [list(x) for x in graph_loops[:graph_loop_cap]]

    graph_holonomy_candidates = []
    gmin_gap = max(65, length // 42)
    for inds in by_parameter.values():
        if len(inds) < 3:
            continue
        if len(inds) > 30:
            sample = [inds[round(k * (len(inds)-1) / 29)] for k in range(30)]
            sample = sorted(set(sample))
        else:
            sample = inds
        for j in range(1, len(sample)-1):
            lefts = sample[max(0,j-5):j]
            rights = sample[j+1:min(len(sample),j+6)]
            for a in lefts:
                for c in rights:
                    b = sample[j]
                    if b-a < gmin_gap or c-b < gmin_gap:
                        continue
                    score = 0.06*gh1[2*a] + 0.0025*(c-a)
                    graph_holonomy_candidates.append((score,a,b,c))
    graph_holonomy=[]
    graph_holonomy_cap = 3 if public else 6
    for _,a,b,c in sorted(graph_holonomy_candidates, reverse=True):
        if any(abs(a-x)<gmin_gap//2 and abs(b-y)<gmin_gap//2 and abs(c-z)<gmin_gap//2 for x,y,z in graph_holonomy):
            continue
        graph_holonomy.append((a,b,c))
        if len(graph_holonomy)>=graph_holonomy_cap:
            break
    if len(graph_holonomy) < max(2, graph_holonomy_cap//2):
        raise RuntimeError('not enough graph holonomy probes')
    q['graph_holonomy_probes']=[list(x) for x in graph_holonomy]
    return q


def _index_key(rec):
    return tuple((d['dimension'], d['space_dimension_F2'], tuple(map(tuple, d['elementary_divisors_F2'])))
                 for d in rec if d['space_dimension_F2'])


def choose_index_probes(q, frames, analyses, nx, seed, public=False):
    """Conley index-map probes [occurrence, morse_ids].  Public probes use single
    Morse sets only; hidden probes are mostly multi-set Morse intervals, where the
    exit set E is the part of P that can no longer reach the interval."""
    rng = random.Random(seed)
    path = q['parameter_ids']
    occ = sorted(rng.sample(range(len(path)), min(len(path), 40 if public else 70)))
    seen_par = set(); cands = []
    for k in occ:
        p = path[k]
        if p in seen_par:
            continue
        seen_par.add(p)
        ana = analyses[p]; desc = ana['desc']; n = len(ana['recurrent'])
        ints = {(m,) for m in range(n)}
        if True:
            for a in range(n):
                ints.add(tuple(sorted({a, *desc[a]})))
                ints.add(tuple(sorted({a, *[c for c in range(n) if a in desc[c]]})))
                for b in desc[a]:
                    ints.add(tuple(sorted({a, b, *[c for c in range(n) if c in desc[a] and b in desc[c]]})))
            ints.add(tuple(range(n)))
        for I in ints:
            if not ref.is_morse_interval(ana, I):
                continue
            rec = ref.conley_index_record(frames[p], ana, list(I), nx)
            cands.append((k, list(I), _index_key(rec)))
    if public:
        singles_only = [c for c in cands if len(c[1]) == 1]
        multis = [c for c in cands if len(c[1]) > 1]
        return sorted(_pick_index(singles_only, 8, 8, rng) + _pick_index(multis, 6, 0, rng))
    return _pick_index(cands, 24, 4, rng)


def _pick_index(cands, cap, max_single, rng):
    freq = Counter(c[2] for c in cands)
    chosen = []; seen = set(); singles = 0; empties = 0
    order = sorted(cands, key=lambda c: (-(len(c[2]) > 0), 1.0 / freq[c[2]] * -1, -len(c[1]), rng.random()))
    # first pass: one probe per distinct nontrivial type, rarest first
    for k, I, key in sorted(cands, key=lambda c: (freq[c[2]], -len(c[1]), c[0])):
        if len(chosen) >= cap:
            break
        if key in seen:
            continue
        if not key and empties >= 1:
            continue
        if len(I) == 1 and singles >= max_single:
            continue
        seen.add(key); chosen.append([k, I]); singles += len(I) == 1; empties += not key
    # second pass: fill with remaining multi-set intervals of nontrivial type
    for k, I, key in order:
        if len(chosen) >= cap:
            break
        if [k, I] in chosen or not key or (len(I) == 1 and singles >= max_single):
            continue
        chosen.append([k, I]); singles += len(I) == 1
    chosen.sort()
    return chosen


def build_case(cfg):
    name, nx, na, nb, scale, r0, dr, s0, ds, cnum, cden, seed = cfg
    data = make_base(name, nx, na, nb, scale, r0, dr, s0, ds, cnum, cden)
    frames, analyses, edge_rel, edge_struct, edge_score, adj = analyze_base(data)
    rng = random.Random(seed + 7000)
    data['queries'] = [
        choose_query(
            i, L, rng, frames, analyses, edge_rel, edge_struct, edge_score, adj,
            na, nb, nx, public=(name == 'public')
        )
        for i, L in enumerate(LENGTHS[name])
    ]
    add_index_probes(data, frames, analyses, seed, name == 'public')
    return data


def add_index_probes(data, frames, analyses, seed, public):
    nx = data['phase_grid']['nx']
    for q in data['queries']:
        q['index_probes'] = choose_index_probes(q, frames, analyses, nx, seed + 9100 + q['query_id'], public)


def withhold_public_interval_divisors(out):
    """The public example shows multi-set interval probes with their index-space
    dimensions and index-map ranks, but withholds the elementary divisors (null)."""
    out = json.loads(json.dumps(out))
    for q in out['queries']:
        for r in q.get('conley_index_maps', []):
            if len(r['morse_ids']) > 1:
                for d in r['degrees']:
                    d['elementary_divisors_F2'] = None
    return out


def write_json(path, obj):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, separators=(',', ':')) + '\n', encoding='utf-8')


def generate_named_case(name):
    cfg = next(c for c in CONFIGS if c[0] == name)
    t = time.time()
    data = build_case(cfg)
    out = ref.solve_case(data)
    if name == 'public':
        write_json(PUBLIC, data)
        write_json(PUBLIC_EXPECTED, withhold_public_interval_divisors(out))
    else:
        write_json(HIDDEN_DIR / f'{name}.json', data)
        write_json(ORACLE_DIR / f'{name}.json', out)
    print(
        name,
        'generated', round(time.time() - t, 2), 's',
        'lengths', [len(q['parameter_ids']) for q in data['queries']],
        'conley_windows', [len(q['rank_windows']) for q in data['queries']],
        'graph_windows', [len(q['graph_rank_windows']) for q in data['queries']],
        'graph_loops', [len(q['graph_loop_windows']) for q in data['queries']],
        'holonomy', [len(q['holonomy_probes']) for q in data['queries']],
        'graph_holonomy', [len(q['graph_holonomy_probes']) for q in data['queries']],
        flush=True,
    )


def main():
    import argparse
    import subprocess
    import sys
    ap = argparse.ArgumentParser()
    ap.add_argument('--case', choices=[c[0] for c in CONFIGS])
    args = ap.parse_args()
    HIDDEN_DIR.mkdir(parents=True, exist_ok=True)
    ORACLE_DIR.mkdir(parents=True, exist_ok=True)
    if args.case:
        generate_named_case(args.case)
        return
    for cfg in CONFIGS:
        subprocess.run([sys.executable, str(Path(__file__).resolve()), '--case', cfg[0]], check=True)


if __name__ == '__main__':
    main()

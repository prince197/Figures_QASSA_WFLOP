"""Revision 3 (item C2 / R8): direct optimization with 1-degree direction bins (experiment rev3_fine).

Usage:  python3 rev3_fine.py SHARD NSHARDS [--procs P] [--rose 15]
        python3 rev3_fine.py --validate          (checks (1) and (2) of rev3_fine_manifest.md; ~3 min, 1 process)
        python3 rev3_fine.py --timing [SEEDS]     (times runs per method at the largest N; scratch output only)
Output: rev3_fine_s<SHARD>of<NSHARDS>.csv (1-deg arm, primary) or, with --rose 15, rev3_fine15_s<SHARD>of<NSHARDS>.csv
        (optional 15-deg control arm). Rows are appended as runs finish; a restarted shard skips the keys already in
        its file (resumable). The tasks are ordered heaviest case first and split round-robin (tasks[SHARD::NSHARDS]),
        so every shard gets the same mix of cases and methods (sizes differ by at most one run).

Design (prespecified in rev3_fine_manifest.md): 8 cases (data sets I/II x (500 m, 10), (750 m, 6), (750 m, 12),
(1000 m, 15)), 5 methods (PSOBV = PSO-VNS, PSOC = PSO, GA, SLSQP = MS-SLSQP, RSDVNS = RSD-VNS), seeds 31-60,
6,030 objective evaluations, random starts; 1,200 runs per arm.

Objective: the benchmark objective with 15 sub-bins per 15-deg bin (360 bins of 1 deg; each sub-bin keeps the
parent bin's Weibull parameters and gets 1/15 of its frequency) = mpce_direction.bench_objective(xy, ds, 15), the
evaluator of the paper's 1-deg re-evaluation. FineObjective computes the same function faster:
  * per ordered turbine pair, the directions in which the in-wake test can hold are found analytically (front cone:
    |phi| < atan K + asin(c), back part of the cone between apex and turbine: |phi - pi| < asin(c) - atan K, with
    c = R / (D sqrt(1 + K^2)); all directions if c >= 1), with a 1e-6 rad safety margin, and the EXACT test of
    wflop_model.jensen_deficits (arccos(clip(num / den)) < arctan K) is evaluated on these candidates only, with the
    same floating-point expressions;
  * the squared deficits are summed over the same dense (direction, turbine, turbine) array as the reference, and
    the expected power is computed only for (direction, turbine) entries with a nonzero deficit (elsewhere it equals
    the wake-free value of the sub-bin), with the elementwise arithmetic of wflop_model.expected_power_linear (the
    survival function is evaluated once per grid speed instead of twice);
  * the farm sum is the same dense weighted sum as the reference.
Result: bit-for-bit equal to bench_objective (validation (2)) at about 1/4 of its cost (see --timing).
The penalized objective is authors_objective.make_objective with farm_objective replaced by the 1-deg objective
(same penalty expression); MS-SLSQP minimizes the unpenalized 1-deg wake loss. Methods are run through the unchanged
mpce_experiments.run_method / run_method_x and rev2_ga.run_method_ga; with --rose 15 the driver is run_grid itself
(validation (1) reproduces stored runs bit for bit).
"""
import os, sys, time, platform
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import init_hook                                            # noqa: E402
import wflop_model as W                                     # noqa: E402
from wflop_model import farm_objective, min_spacing, R     # noqa: E402
from authors_objective import make_objective, PENALTY      # noqa: E402
import mpce_experiments as mx                               # noqa: E402
import rev2_ga as ga                                        # noqa: E402
from record_io import encode_coordinates                    # noqa: E402

EXP = "rev3_fine"
CASES = [(ds, r, n) for ds in (1, 2) for r, n in ((500, 10), (750, 6), (750, 12), (1000, 15))]
METHODS = ["PSOBV", "PSOC", "GA", "SLSQP", "RSDVNS"]
SEEDS = range(31, 61)
BUDGET = 6030
SUB = 15                       # sub-bins per 15-deg bin -> 1-deg bins
NP = mx.NP
KEY = ["Algorithm", "Dataset", "Radius", "Turbines", "Seed", "Budget", "Init", "Rose"]


# ====================================================================== 1-deg objective
class FineObjective:
    """Benchmark objective with `sub` sub-bins per 15-deg bin; equal (bit for bit) to
    mpce_direction.bench_objective(xy, ds, sub) for the Jensen wake."""

    def __init__(self, ds, sub=SUB):
        k, psi, w = W.DATASETS[int(ds)]
        j = np.repeat(np.arange(24), sub)
        m = np.tile(np.arange(sub), 24)
        deg = 15.0 * j + (m + 0.5) * 15.0 / sub              # as mpce_direction.bench_bins
        self.th = np.deg2rad(deg)
        self.nt = len(self.th)
        self.kk, self.pp, self.ww = k[j], psi[j], w[j] / sub
        self.cos, self.sin = np.cos(self.th), np.sin(self.th)
        self.ep0 = W.expected_power_linear(self.kk, self.pp)
        self.step = 15.0 / sub
        self.sgrid = np.arange(W.S_CI, W.S_R + 1e-9, 0.5)
        s = self.sgrid
        self.mid = 0.5 * (s[:-1] + s[1:])
        self.ww_col = self.ww[:, None]
        self._pairs = {}
        self._ideal = {}

    def ideal(self, n):
        if n not in self._ideal:
            self._ideal[n] = float(W.BIN_WIDTH * n * np.sum(self.ww * W.expected_power_linear(self.kk, self.pp)))
        return self._ideal[n]

    def _epl(self, k, psi):
        """wflop_model.expected_power_linear for 1-D k, psi (same elementwise arithmetic; survival function evaluated
        once per grid speed: the lo / hi arrays of the reference share 20 of their 21 points)."""
        G = np.exp(-(self.sgrid / psi[:, None]) ** np.broadcast_to(k[:, None], (len(k), len(self.sgrid))))
        lin = W.LAM * np.sum(self.mid * (G[:, :-1] - G[:, 1:]), axis=-1)
        lin += W.ETA * (G[:, 0] - G[:, -1])
        return lin + W.P_RATED * G[:, -1]

    def _candidates(self, xy):
        n = len(xy)
        if n not in self._pairs:
            self._pairs[n] = np.nonzero(~np.eye(n, dtype=bool))
        ii, jj = self._pairs[n]
        dx = xy[ii, 0] - xy[jj, 0]                    # rho_i - rho_j as in jensen_deficits (i waked by j)
        dy = xy[ii, 1] - xy[jj, 1]
        dd = np.hypot(dx, dy)
        K = W.K
        with np.errstate(divide="ignore", invalid="ignore"):
            c = np.where(dd > 0, W.R / (dd * np.sqrt(1.0 + K * K)), np.inf)
        full = c >= 1.0 - 1e-9
        asn = np.arcsin(np.minimum(c, 1.0))
        atk = np.arctan(K)
        marg = 1e-6
        h0 = np.where(full, np.pi + 1.0, atk + asn + marg)       # front part of the cone (i downstream of j)
        h1 = np.where(full, -1.0, asn - atk + marg)             # back part (i between the cone apex and j)
        phi = np.rad2deg(np.arctan2(dy, dx))
        ps, ts = [], []
        for cen, h in ((phi, h0), (phi + 180.0, h1)):
            hd = np.rad2deg(h)
            lo = np.ceil((cen - hd) / self.step - 0.5).astype(np.int64)
            hi = np.floor((cen + hd) / self.step - 0.5).astype(np.int64)
            cnt = np.where(h > 0, np.clip(hi - lo + 1, 0, self.nt), 0)
            tot = int(cnt.sum())
            if tot == 0:
                continue
            ps.append(np.repeat(np.arange(len(dx)), cnt))
            off = np.arange(tot) - np.repeat(np.cumsum(cnt) - cnt, cnt)
            ts.append((np.repeat(lo, cnt) + off) % self.nt)
        if not ps:
            return ii, jj, dx, dy, np.zeros(0, np.int64), np.zeros(0, np.int64)
        return ii, jj, dx, dy, np.concatenate(ps), np.concatenate(ts)

    def deficits(self, xy):
        """Combined Jensen deficit, shape (n_theta, N) = wflop_model.jensen_deficits(xy, self.th)."""
        xy = np.asarray(xy, float)
        n = len(xy)
        ii, jj, dx, dy, p, t = self._candidates(xy)
        sq = np.zeros((self.nt, n, n))
        if len(p):
            ddx, ddy = dx[p], dy[p]
            c, s = self.cos[t], self.sin[t]
            a = R / W.K
            num = ddx * c + ddy * s + a
            den = np.sqrt((ddx + a * c) ** 2 + (ddy + a * s) ** 2)
            beta = np.arccos(np.clip(num / den, -1.0, 1.0))
            inw = beta < np.arctan(W.K)
            d = np.abs(ddx[inw] * c[inw] + ddy[inw] * s[inw])
            dij = (1 - np.sqrt(1 - W.CT)) / (1 + W.K * d / R) ** 2
            pw = p[inw]
            sq[t[inw], ii[pw], jj[pw]] = dij ** 2
        return np.sqrt(np.sum(sq, axis=2))

    def __call__(self, xy):
        """(objective, ideal) in legacy units (15 x expected farm power, kW), as farm_objective."""
        xy = np.asarray(xy, float)
        delta = self.deficits(xy)
        ep = np.broadcast_to(self.ep0[:, None], delta.shape).copy()
        wm = delta > 0
        if wm.any():
            tt = np.nonzero(wm)[0]
            ep[wm] = self._epl(self.kk[tt], self.pp[tt] * (1 - delta[wm]))
        return float(W.BIN_WIDTH * np.sum(self.ww_col * ep)), self.ideal(len(xy))


_FINE = {}


def fine_objective(ds):
    if ds not in _FINE:
        _FINE[ds] = FineObjective(ds)
    return _FINE[ds]


def make_objective_with(evalf, radius, smin=8 * R):
    """authors_objective.make_objective (penalty="authors") with farm_objective replaced by evalf(xy) ->
    (objective, ideal); same penalty expression."""
    def f(x):
        xy = x.reshape(-1, 2)
        obj, ideal = evalf(xy)
        n = len(xy)
        dd = np.sqrt(((xy[:, None] - xy[None]) ** 2).sum(-1))[np.triu_indices(n, 1)]
        gb = (xy ** 2).sum(1) - radius ** 2
        gs = smin - dd
        pen = (((1 + PENALTY * gb[gb > 0]) ** 2).sum() + ((1 + PENALTY * gs[gs > 0]) ** 2).sum())
        return ideal - obj + pen
    return f


# ====================================================================== one run
def dispatch(alg):
    return ga.run_method_ga if alg == "GA" else mx.run_method_x      # run_method_x: RSDVNS / PSOBV90, else run_method


def run_case(task):
    """One run; task = (alg, ds, r, n, seed, budget, rose) with rose in {"1deg", "15deg"}. Mirrors
    mpce_experiments.run_grid (random initialization)."""
    alg, ds, rad, n, seed, budget, rose = task
    smin = 8 * R
    init_hook.GEN = None
    fo = fine_objective(ds)
    if rose == "15deg":
        f = make_objective(ds, rad)                       # the legacy objective object of run_grid

        def evalf(xy):
            return farm_objective(xy, ds)
    elif rose == "1deg":
        evalf = fo
        f = make_objective_with(evalf, rad, smin)
    else:
        raise ValueError(rose)

    def wake(x):
        o, i = evalf(x.reshape(-1, 2))
        return i - o

    def feasible(x):
        xy = x.reshape(-1, 2)
        return (n == 1 or min_spacing(xy) >= smin - 1e-6) and np.sqrt((xy ** 2).sum(1)).max() <= rad + 1e-6
    t0 = time.perf_counter()
    pos, tr = dispatch(alg)(alg, seed, budget, f, wake, feasible, 2 * n, -rad, rad, rad, smin)
    sec = time.perf_counter() - t0
    xy = np.asarray(pos, float).reshape(-1, 2)
    obj, ideal = evalf(xy)
    o15, i15 = farm_objective(xy, ds)
    of, iff = fo(xy)
    curve = ideal - np.array(tr.curve)
    return dict(Algorithm=alg, Dataset=ds, Radius=rad, Turbines=n, Seed=seed, Budget=budget, Init="random",
                Calls=tr.calls, Objective=obj, Ideal=ideal, WakeLoss=ideal - obj, Feasible=bool(feasible(pos)),
                MinSpacing=min_spacing(xy), Seconds=sec, Coordinates=encode_coordinates(xy),
                Curve=";".join("nan" if not np.isfinite(c) else f"{c:.3f}" for c in curve),
                Rose=rose, SubBins=SUB if rose == "1deg" else 1, Obj15=o15, Ideal15=i15, ObjFine=of, IdealFine=iff,
                Experiment=EXP if rose == "1deg" else EXP + "15")


def tasks(rose="1deg"):
    """All runs of one arm, heaviest first (cases by N descending, then SLSQP / PSOBV / ... ), seeds 31-60."""
    cases = sorted(CASES, key=lambda c: (-c[2], -c[1], c[0]))
    return [(a, ds, r, n, s, BUDGET, rose) for (ds, r, n) in cases for a in METHODS for s in SEEDS]


def _key(row):
    return (str(row["Algorithm"]), int(row["Dataset"]), int(row["Radius"]), int(row["Turbines"]), int(row["Seed"]),
            int(row["Budget"]), str(row["Init"]), str(row["Rose"]))


def platform_info():
    flags = ""
    try:
        with open("/proc/cpuinfo") as fh:
            for line in fh:
                if line.startswith("model name"):
                    flags = line.split(":", 1)[1].strip(); break
    except OSError:
        pass
    import scipy
    try:
        from numpy._core._multiarray_umath import __cpu_features__ as cf
        simd = ",".join(k for k in ("AVX2", "AVX512F", "AVX512_SKX") if cf.get(k))
    except ImportError:
        simd = "?"
    return (f"python {platform.python_version()}, numpy {np.__version__}, scipy {scipy.__version__}, "
            f"pandas {pd.__version__}, cpu {flags}, numpy SIMD {simd}")


def selftest(n_layouts=20, seed=12345):
    """Quick check on this machine: FineObjective == mpce_direction.bench_objective(., ., 15), bit for bit."""
    import mpce_direction as MD
    rng = np.random.default_rng(seed)
    bad = 0
    for ds in (1, 2):
        fo = FineObjective(ds)
        for i in range(n_layouts):
            n = int(rng.integers(2, 16)); r = float(rng.choice([500.0, 750.0, 1000.0]))
            xy = rng.uniform(-r, r, (n, 2))
            if i % 5 == 0:
                xy[1] = xy[0] + rng.normal(0, 10, 2)
            bad += fo(xy) != MD.bench_objective(xy, ds, 15)
    return bad == 0


def run_shard(shard, nsh, procs, rose):
    tl = tasks(rose)[shard::nsh]
    tag = EXP if rose == "1deg" else EXP + "15"
    fn = os.path.join(HERE, f"{tag}_s{shard}of{nsh}.csv")
    done = set()
    if os.path.exists(fn) and os.path.getsize(fn) > 0:
        old = pd.read_csv(fn, usecols=KEY)
        done = {_key(r) for _, r in old.iterrows()}
    todo = [t for t in tl if (t[0], t[1], t[2], t[3], t[4], t[5], "random", t[6]) not in done]
    print(f"{tag} shard {shard}/{nsh}: {len(tl)} runs, {len(tl) - len(todo)} already in {os.path.basename(fn)}, "
          f"{len(todo)} to do; procs {procs}; {platform_info()}", flush=True)
    if not selftest():
        raise SystemExit("selftest failed: FineObjective differs from mpce_direction.bench_objective on this machine")
    t0 = time.time()
    header = not (os.path.exists(fn) and os.path.getsize(fn) > 0)

    def write(row):
        nonlocal header
        pd.DataFrame([row]).to_csv(fn, mode="a", header=header, index=False)
        header = False
    if procs <= 1:
        it = map(run_case, todo)
        pool = None
    else:
        from multiprocessing import Pool
        pool = Pool(procs)
        it = pool.imap_unordered(run_case, todo, chunksize=1)
    k = 0
    for row in it:
        write(row); k += 1
        print(f"  [{k}/{len(todo)}] {row['Algorithm']} ds{row['Dataset']} r{row['Radius']} N{row['Turbines']} "
              f"seed {row['Seed']}: loss {100 * row['WakeLoss'] / row['Ideal']:.4f}% feas {row['Feasible']} "
              f"calls {row['Calls']} {row['Seconds']:.1f}s (elapsed {time.time() - t0:.0f}s)", flush=True)
    if pool is not None:
        pool.close(); pool.join()
    print(f"{tag} shard {shard}/{nsh} done: {k} runs in {time.time() - t0:.0f} s", flush=True)


# ====================================================================== validation and timing
VAL_CASES = [(1, 750, 6), (2, 1000, 15)]
VAL_SRC = {"PSOBV": "mpce_psobv", "PSOC": "mpce_psoc", "SLSQP": "mpce_slsqp", "RSDVNS": "mpce_rsdisc", "GA": "rev2_ga"}


def _stored(alg, case, seed=1):
    import glob
    ds, r, n = case
    files = sorted(glob.glob(os.path.join(HERE, VAL_SRC[alg] + "_s*of*.csv")))
    for f in files:
        d = pd.read_csv(f, float_precision="round_trip")         # the default parser can be off by 1 ulp
        d = d[(d.Algorithm == alg) & (d.Dataset.astype(str) == str(ds)) & (d.Radius == r) & (d.Turbines == n)
              & (d.Seed == seed) & (d.Budget == 6030 if "Budget" in d else True)]
        if "Init" in d:
            d = d[d.Init == "random"]
        if len(d):
            return d.iloc[0]
    return None


def validate():
    import mpce_direction as MD
    print("platform:", platform_info(), flush=True)
    ok_all = True
    print("(1) 15-deg objective swapped back: reproduction of stored seed-1 runs", flush=True)
    for case in VAL_CASES:
        for alg in METHODS:
            st = _stored(alg, case)
            row = run_case((alg, *case, 1, BUDGET, "15deg"))
            c3 = ";".join(f"{a:.3f} {b:.3f}" for a, b in np.array(W.parse_coords(row["Coordinates"])))
            checks = dict(Objective=row["Objective"] == float(st.Objective), Ideal=row["Ideal"] == float(st.Ideal),
                          Calls=int(row["Calls"]) == int(st.Calls),
                          Feasible=bool(row["Feasible"]) == (str(st.Feasible).lower() == "true"),
                          Curve=row["Curve"] == st.Curve, Coordinates3=c3 == st.Coordinates)
            ok = all(checks.values()); ok_all &= ok
            print(f"  {alg:7s} ds{case[0]} r{case[1]} N{case[2]}: {'BIT-IDENTICAL' if ok else 'DIFFERENT'} "
                  f"objective {row['Objective']!r} vs stored {float(st.Objective)!r}; {checks}; {row['Seconds']:.1f}s",
                  flush=True)
    print("(2) 1-deg objective of stored layouts vs mpce_direction.bench_objective(xy, ds, 15)", flush=True)
    L = pd.read_csv(os.path.join(HERE, "mpce_direction_layouts.csv.gz"))
    import glob
    parts = []
    for alg, src in VAL_SRC.items():
        for f in sorted(glob.glob(os.path.join(HERE, src + "_s*of*.csv"))):
            d = pd.read_csv(f, usecols=["Algorithm", "Dataset", "Radius", "Turbines", "Seed", "Budget", "Init",
                                        "Feasible", "Coordinates"])
            parts.append(d[(d.Algorithm == alg) & d.Dataset.astype(str).isin(["1", "2"]) & (d.Budget == 6030)
                           & (d.Init == "random")])
    S = pd.concat(parts, ignore_index=True)
    S["Dataset"] = S.Dataset.astype(int)
    S = S[[c in CASES for c in zip(S.Dataset, S.Radius, S.Turbines)]]
    nbit = nall = 0; maxrel = 0.0; t_ref = t_fast = 0.0; dev_j15 = []
    Lk = L.set_index(["Algorithm", "Dataset", "Radius", "Turbines", "Seed"])
    for _, r in S.iterrows():
        xy = W.parse_coords(r.Coordinates)
        a = time.perf_counter(); ref = MD.bench_objective(xy, r.Dataset, 15); t_ref += time.perf_counter() - a
        a = time.perf_counter(); mine = fine_objective(r.Dataset)(xy); t_fast += time.perf_counter() - a
        nall += 1; nbit += (ref == mine)
        maxrel = max(maxrel, abs(ref[0] - mine[0]) / abs(ref[0]))
        k = (r.Algorithm, r.Dataset, r.Radius, r.Turbines, r.Seed)
        if str(r.Feasible).lower() == "true" and k in Lk.index:
            dev_j15.append(float(f"{mine[0]:.8g}") == float(Lk.loc[k].J15))      # stored with float_format %.8g
    ok2 = nbit == nall
    ok_all &= ok2
    print(f"  {nall} stored layouts (5 methods x 8 cases x seeds 1-30, all labels): bit-identical {nbit}/{nall}, "
          f"max rel. dev {maxrel:.2e}; vs stored J15 column of mpce_direction_layouts.csv.gz (written with %.8g; "
          f"{len(dev_j15)} feasible layouts of the 4 methods in it): equal after rounding to 8 significant digits "
          f"{sum(dev_j15)}/{len(dev_j15)}; time per evaluation: reference "
          f"{1e3 * t_ref / nall:.2f} ms, fast {1e3 * t_fast / nall:.2f} ms", flush=True)
    ok_all &= all(dev_j15)
    print("VALIDATION", "PASSED" if ok_all else "FAILED", flush=True)
    return ok_all


def timing(seeds=(31, 32)):
    """Times runs of every method at the largest N (data set II, 1000 m, N = 15) with both roses (scratch only)."""
    import mpce_direction as MD
    rng = np.random.default_rng(1)
    for ds in (1, 2):
        fo = fine_objective(ds)
        for n, r in ((15, 1000), (12, 750), (10, 500), (6, 750)):
            xs = [rng.uniform(-r, r, 2 * n) for _ in range(200)]
            fl, ff = make_objective(ds, r), make_objective_with(fo, r)
            tt = {}
            for name, fn in (("legacy15", lambda x: fl(x)), ("fine1", lambda x: ff(x)),
                             ("ref1", lambda x: MD.bench_objective(x.reshape(-1, 2), ds, 15))):
                a = time.perf_counter()
                for x in xs:
                    fn(x)
                tt[name] = 1e3 * (time.perf_counter() - a) / len(xs)
            print(f"  objective ds{ds} N{n}: penalized 15-deg {tt['legacy15']:.3f} ms, penalized 1-deg (fast) "
                  f"{tt['fine1']:.3f} ms ({tt['fine1'] / tt['legacy15']:.1f}x), reference bench_objective(15) "
                  f"{tt['ref1']:.3f} ms ({tt['ref1'] / tt['legacy15']:.1f}x)", flush=True)
    rows = []
    for alg in METHODS:
        for s in seeds:
            row = run_case((alg, 2, 1000, 15, s, BUDGET, "1deg"))
            rows.append(row)
            print(f"  run {alg:7s} ds2 r1000 N15 seed {s}: {row['Seconds']:.1f} s, loss "
                  f"{100 * row['WakeLoss'] / row['Ideal']:.4f}%, feasible {row['Feasible']}, calls {row['Calls']}",
                  flush=True)
    return pd.DataFrame(rows)


if __name__ == "__main__":
    argv = sys.argv[1:]
    if "--validate" in argv:
        sys.exit(0 if validate() else 1)
    if "--timing" in argv:
        rest = [int(a) for a in argv if a.isdigit()]
        T = timing(tuple(rest) if rest else (31, 32))
        out = os.environ.get("REV3_FINE_TIMING_OUT")
        if out:
            T.to_csv(out, index=False)
        sys.exit(0)
    procs, rose, args, i = 1, "1deg", [], 0
    while i < len(argv):
        a = argv[i]
        if a.startswith("--procs="):
            procs = int(a.split("=")[1])
        elif a == "--procs":
            procs = int(argv[i + 1]); i += 1
        elif a.startswith("--rose="):
            rose = {"1": "1deg", "15": "15deg"}[a.split("=")[1]]
        elif a == "--rose":
            rose = {"1": "1deg", "15": "15deg"}[argv[i + 1]]; i += 1
        else:
            args.append(a)
        i += 1
    if len(args) != 2:
        raise SystemExit(__doc__)
    run_shard(int(args[0]), int(args[1]), procs, rose)

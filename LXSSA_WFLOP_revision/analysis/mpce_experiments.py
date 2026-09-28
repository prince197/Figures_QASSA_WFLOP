"""Additional experiments for the MPCE resubmission (SSA-VNS study).

Usage:  python mpce_experiments.py EXP [SHARD NSHARDS] [--procs P]
Output: mpce_<EXP>_s<SHARD>of<NSHARDS>.csv (same columns as fresh_*.csv plus Budget, Init).

Method ids: SSABV = SSA-VNS (proposed), LXBV = LX-SSA-VNS, RSVNS = random-multistart VNS
(random sampling phase 1 + VNS), BVNS = basic VNS, SSA, LXSSA, PSOC = PSO with Clerc-Kennedy
constriction coefficients (w = 0.7298, c1 = c2 = 1.49618), DE, SLSQP = MS-SLSQP.

Experiments
  rsvns    RS-VNS on the 68 benchmark cases, 30 seeds, 6,030 calls
  slsqp    MS-SLSQP rerun on the 68 cases (SLSQP results are platform dependent; rerun for consistency)
  psobv    PSO-VNS (constriction PSO phase 1 + VNS) on the 68 cases and HR16, 30 seeds
  psoc     constriction PSO on the 68 cases, 30 seeds, 6,030 calls
  ssasplit SSA-VNS with 25% / 75% budget split on the 12 split cases
  hr16new  PSOC and RSVNS on the Horns Rev 1 16-turbine block, 30 seeds
  feas     feasibility-preserving initialization, 6,030 calls, 9 methods, 6 largest cases + HR16
  b30k     random initialization, 30,030 calls, 9 methods, 6 largest cases + HR16 (30 seeds)
  b120k    random initialization, 120,030 calls, 9 methods, 6 largest cases (30 seeds) + HR16 (10 seeds)
Convergence curves have 201 checkpoints (every (B-30)/200 calls).
"""
import sys, time
import numpy as np, pandas as pd
from multiprocessing import Pool
import init_hook
from feasible_init import make_generator
from authors_optimizers import SSA, LXSSA, PSO, DE
from extra_baselines import MSSLSQP
from original_vns import BVNS
from hybrid_lxssa_bvns import HybridBVNS
from rs_vns import RSVNS
from authors_objective import make_objective
from wflop_model import farm_objective, min_spacing, R
import hornsrev_model as hr

GRID = [(ds, r, n) for ds in (1, 2) for r, nmax in ((500, 10), (750, 12), (1000, 15))
        for n in range(2, nmax + 1)]
SPLITCASES = [(ds, r, n) for ds in (1, 2) for r, n in ((500, 6), (750, 8), (1000, 10), (500, 10), (750, 12), (1000, 15))]
LARGE = [(ds, r, n) for ds in (1, 2) for r, n in ((500, 10), (750, 12), (1000, 15))]
M9 = ["SSABV", "LXBV", "RSVNS", "BVNS", "SSA", "LXSSA", "PSOC", "DE", "SLSQP"]
NP = 30


class Tracker:
    def __init__(self, feasible, check):
        self.feasible = feasible; self.check = check
        self.calls = 0; self.best = np.inf; self.curve = []

    def see(self, x, wake_loss):
        self.calls += 1
        if wake_loss < self.best and self.feasible(x):
            self.best = wake_loss
        if self.calls % self.check == 0:
            self.curve.append(self.best)


def run_method(alg, seed, budget, f, wake, feasible, dim, lb, ub, radius, smin, bcons=None):
    tr = Tracker(feasible, max(1, (budget - NP) // 200))

    def fp(x):
        v = f(x); tr.see(x, v); return v

    def wk(x):
        v = wake(x); tr.see(x, v); return v
    it1 = (budget - NP) // NP                  # SSA, PSO, DE iterations (one call per member)
    if alg in ("SSABV", "SSABV25", "SSABV75", "LXBV", "LXBV25", "LXBV75"):
        split = {"25": 0.25, "75": 0.75}.get(alg[-2:], 0.5)
        pos, _, _ = HybridBVNS(NP, budget, split, 0.0, 1.0, radius, "SSA" if alg.startswith("SSA") else "LXSSA",
                               seed=seed).optimize(fp, dim, lb, ub)
    elif alg in ("PSOBV", "PSOBV25", "PSOBV75"):
        split = {"25": 0.25, "75": 0.75}.get(alg[-2:], 0.5)
        pos, _, _ = HybridBVNS(NP, budget, split, 0.0, 1.0, radius, "PSOC", seed=seed).optimize(fp, dim, lb, ub)
    elif alg == "RSVNS":
        pos, _, _ = RSVNS(NP, budget, 0.5, radius, seed=seed).optimize(fp, dim, lb, ub)
    elif alg == "BVNS":
        pos, _, _ = BVNS(NP, budget, radius, seed=seed).optimize(fp, dim, lb, ub)
    elif alg == "SLSQP":
        pos, _, _ = MSSLSQP(wk, radius, smin, NP, budget, seed=seed, boundary_cons=bcons).optimize(fp, dim, lb, ub)
    elif alg == "LXSSA":
        pos, _, _ = LXSSA(NP, (budget - NP) // (2 * NP), 0.0, 1.0, seed=seed).optimize(fp, dim, lb, ub)
    elif alg == "SSA":
        pos, _, _ = SSA(NP, it1, seed=seed).optimize(fp, dim, lb, ub)
    elif alg == "PSOC":
        pos, _, _ = PSO(NP, it1, w=0.7298, c1=1.49618, c2=1.49618, seed=seed).optimize(fp, dim, lb, ub)
    elif alg == "PSO":
        pos, _, _ = PSO(NP, it1, seed=seed).optimize(fp, dim, lb, ub)
    elif alg == "DE":
        pos, _, _ = DE(NP, it1, seed=seed).optimize(fp, dim, lb, ub)
    else:
        raise ValueError(alg)
    return pos, tr


def run_grid(task):
    alg, ds, rad, n, seed, budget, init = task
    smin = 8 * R
    init_hook.GEN = make_generator(smin, circle_r=rad) if init == "feasible" else None
    f = make_objective(ds, rad)

    def wake(x):
        o, i = farm_objective(x.reshape(-1, 2), ds)
        return i - o

    def feasible(x):
        xy = x.reshape(-1, 2)
        return (n == 1 or min_spacing(xy) >= smin - 1e-6) and np.sqrt((xy ** 2).sum(1)).max() <= rad + 1e-6
    t0 = time.perf_counter()
    pos, tr = run_method(alg, seed, budget, f, wake, feasible, 2 * n, -rad, rad, rad, smin)
    sec = time.perf_counter() - t0
    xy = pos.reshape(-1, 2)
    obj, ideal = farm_objective(xy, ds)
    curve = ideal - np.array(tr.curve)
    return dict(Algorithm=alg, Dataset=ds, Radius=rad, Turbines=n, Seed=seed, Budget=budget, Init=init,
                Calls=tr.calls, Objective=obj, Ideal=ideal, WakeLoss=ideal - obj, Feasible=bool(feasible(pos)),
                MinSpacing=min_spacing(xy), Seconds=sec,
                Coordinates=";".join(f"{a:.3f} {b:.3f}" for a, b in xy),
                Curve=";".join("nan" if not np.isfinite(c) else f"{c:.3f}" for c in curve))


def run_hr(task):
    alg, n, seed, budget, init = task
    f, ideal, poly = hr.make_objective(n)
    smin = 4 * hr.D
    init_hook.GEN = make_generator(smin, poly=poly) if init == "feasible" else None
    half = np.abs(poly).max()
    a, b = poly, np.roll(poly, -1, axis=0)
    e = b - a
    nrm = np.c_[e[:, 1], -e[:, 0]] / np.linalg.norm(e, axis=1)[:, None]

    def wake(x):
        return ideal - hr.aep_gwh(x.reshape(-1, 2))

    def bcons(p):
        return (-np.einsum("ijk,jk->ij", p[:, None, :] - a[None], nrm) / half).ravel()

    def feasible(x):
        xy = x.reshape(-1, 2)
        return min_spacing(xy) >= smin - 1e-6 and hr.outside_distance(xy, poly).max() <= 1e-6
    t0 = time.perf_counter()
    pos, tr = run_method(alg, seed, budget, f, wake, feasible, 2 * n, -half, half, half, smin, bcons)
    sec = time.perf_counter() - t0
    xy = pos.reshape(-1, 2)
    curve = ideal - np.array(tr.curve)
    return dict(Algorithm=alg, Dataset="HR", Radius=0, Turbines=n, Seed=seed, Budget=budget, Init=init,
                Calls=tr.calls, Objective=hr.aep_gwh(xy), Ideal=ideal, WakeLoss=ideal - hr.aep_gwh(xy),
                Feasible=bool(feasible(pos)), MinSpacing=min_spacing(xy), Seconds=sec,
                Coordinates=";".join(f"{u:.2f} {v:.2f}" for u, v in xy),
                Curve=";".join("nan" if not np.isfinite(c) else f"{c:.4f}" for c in curve))


def tasks(exp):
    """List of (function, task) pairs; heaviest first for better load balancing."""
    S30 = range(1, 31)
    if exp == "rsvns":
        return [(run_grid, ("RSVNS", *c, s, 6030, "random")) for c in GRID for s in S30]
    if exp == "slsqp":
        return [(run_grid, ("SLSQP", *c, s, 6030, "random")) for c in GRID for s in S30]
    if exp == "psosplit":
        return [(run_grid, (a, *c, s, 6030, "random")) for c in SPLITCASES for a in ("PSOBV25", "PSOBV75") for s in S30]
    if exp == "psobv":
        return [(run_hr, ("PSOBV", 16, s, 6030, "random")) for s in S30] + \
               [(run_grid, ("PSOBV", *c, s, 6030, "random")) for c in GRID for s in S30]
    if exp == "psoc":
        return [(run_grid, ("PSOC", *c, s, 6030, "random")) for c in GRID for s in S30]
    if exp == "ssasplit":
        return [(run_grid, (a, *c, s, 6030, "random")) for c in SPLITCASES for a in ("SSABV25", "SSABV75") for s in S30]
    if exp == "hr16new":
        return [(run_hr, (a, 16, s, 6030, "random")) for a in ("PSOC", "RSVNS") for s in S30]
    if exp == "hrfix":
        # all Horns Rev 1 16-turbine runs again with the corrected direction binning of hornsrev_model
        # (2026-09-28): 10 methods (M9 + PSO-VNS); 6,030 random / 6,030 feasible (no RS-VNS) /
        # 30,030 random (30 seeds each) and 120,030 random (10 seeds); heaviest first
        m10 = M9 + ["PSOBV"]
        tl = [(run_hr, (a, 16, s, 120030, "random")) for a in m10 for s in range(1, 11)]
        tl += [(run_hr, (a, 16, s, 30030, "random")) for a in m10 for s in S30]
        tl += [(run_hr, (a, 16, s, 6030, "random")) for a in m10 for s in S30]
        tl += [(run_hr, (a, 16, s, 6030, "feasible")) for a in m10 if a != "RSVNS" for s in S30]
        return tl
    if exp == "feasx":
        # feasible initialization without RS-VNS: with feasible initialization RS-VNS must build
        # ~3,000 feasible random layouts per run (one packing solve each), about 100x the cost of
        # the other methods; the remaining 8 methods of M9 are run here
        return [t for t in tasks("feas") if t[1][0] != "RSVNS"]
    if exp.rstrip("p") in ("feas", "b30k", "b120k"):
        # feasp / b30kp / b120kp: the same design for the PSO-VNS arm only
        methods = ["PSOBV"] if exp.endswith("p") else M9
        base = exp.rstrip("p")
        budget = {"feas": 6030, "b30k": 30030, "b120k": 120030}[base]
        init = "feasible" if base == "feas" else "random"
        hseeds = range(1, 11) if base == "b120k" else S30
        tl = [(run_hr, (a, 16, s, budget, init)) for a in methods for s in hseeds]
        tl += [(run_grid, (a, *c, s, budget, init)) for c in LARGE[::-1] for a in methods for s in S30]
        return tl
    raise ValueError(exp)


def _call(pair):
    fn, t = pair
    return fn(t)


if __name__ == "__main__":
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    procs = 4
    for a in sys.argv[1:]:
        if a.startswith("--procs="):
            procs = int(a.split("=")[1])
    exp = args[0]
    shard, nsh = (int(args[1]), int(args[2])) if len(args) >= 3 else (0, 1)
    tl = tasks(exp)[shard::nsh]
    t0 = time.time()
    with Pool(procs) as pool:
        rows = pool.map(_call, tl, chunksize=1)
    out = pd.DataFrame(rows)
    fn = f"mpce_{exp}_s{shard}of{nsh}.csv"
    out.to_csv(fn, index=False)
    print(exp, fn, len(out), "runs in", round(time.time() - t0), "s; calls", out.Calls.min(), out.Calls.max(), flush=True)

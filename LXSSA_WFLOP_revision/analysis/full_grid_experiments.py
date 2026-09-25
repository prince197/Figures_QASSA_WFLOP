"""Fresh numerical study for the revised manuscript.

grid : all 68 benchmark cases (Data Sets I and II; 500 m: N = 2..10, 750 m: N = 2..12,
       1000 m: N = 2..15), six methods, 30 seed-paired runs, 6,030 objective calls each.
hr16 / hr80 : Horns Rev 1 site case (real V80 turbine, wind climate and outline), six methods,
       6,030 calls (30 seeds for the 16-turbine block, 10 for the full 80-turbine farm).

Every run records a convergence curve: the best feasible objective found so far at every 30
calls (201 checkpoints), counted over all objective calls including SLSQP gradient calls.
"""
import sys, time
import numpy as np, pandas as pd
from multiprocessing import Pool
from authors_optimizers import SSA, LXSSA, PSO, DE
from extra_baselines import VNS, MSSLSQP
from hybrid_lxssa_vns import LXSSAVNS
from original_vns import BVNS
from authors_objective import make_objective
from wflop_model import farm_objective, min_spacing, R
import hornsrev_model as hr

BUDGET = 6030
CHECK = 30
ALGS = ["LXSSA", "SSA", "PSO", "DE", "VNS", "SLSQP"]
GRID = [(ds, r, n) for ds in (1, 2) for r, nmax in ((500, 10), (750, 12), (1000, 15))
        for n in range(2, nmax + 1)]
ITERS = {"LXSSA": 100, "SSA": 200, "PSO": 200, "DE": 200}


class Tracker:
    """Counts calls and records the best feasible wake loss found so far."""

    def __init__(self, feasible):
        self.feasible = feasible
        self.calls = 0
        self.best = np.inf
        self.curve = []

    def see(self, x, wake_loss):
        self.calls += 1
        if wake_loss < self.best and self.feasible(x):
            self.best = wake_loss
        if self.calls % CHECK == 0:
            self.curve.append(self.best)


def run_method(alg, seed, f, wake, feasible, dim, lb, ub, radius, smin, bcons=None):
    """f: penalised objective, wake: wake loss only; both are tracked."""
    tr = Tracker(feasible)

    def fp(x):
        v = f(x)
        tr.see(x, v)          # when feasible, F_p equals the wake loss
        return v

    def wk(x):
        v = wake(x)
        tr.see(x, v)
        return v
    if alg.startswith("LXVNS"):
        split = {"LXVNS": 0.5, "LXVNS25": 0.25, "LXVNS75": 0.75}[alg]
        pos, _, _ = LXSSAVNS(30, BUDGET, split, 0.0, 1.0, radius, seed=seed).optimize(fp, dim, lb, ub)
    elif alg == "BVNS":
        pos, _, _ = BVNS(30, BUDGET, radius, seed=seed).optimize(fp, dim, lb, ub)
    elif alg == "VNS":
        pos, _, _ = VNS(30, BUDGET, radius, seed=seed).optimize(fp, dim, lb, ub)
    elif alg == "SLSQP":
        opt = MSSLSQP(wk, radius, smin, 30, BUDGET, seed=seed, boundary_cons=bcons)
        pos, _, _ = opt.optimize(fp, dim, lb, ub)
    else:
        o = LXSSA(30, ITERS[alg], 0.0, 1.0, seed=seed) if alg == "LXSSA" else \
            {"SSA": SSA, "PSO": PSO, "DE": DE}[alg](30, ITERS[alg], seed=seed)
        pos, _, _ = o.optimize(fp, dim, lb, ub)
    return pos, tr


def run_grid(task):
    alg, ds, rad, n, seed = task
    smin = 8 * R
    f = make_objective(ds, rad)

    def wake(x):
        o, i = farm_objective(x.reshape(-1, 2), ds)
        return i - o

    def feasible(x):
        xy = x.reshape(-1, 2)
        return (n == 1 or min_spacing(xy) >= smin - 1e-6) and np.sqrt((xy ** 2).sum(1)).max() <= rad + 1e-6
    t0 = time.perf_counter()
    pos, tr = run_method(alg, seed, f, wake, feasible, 2 * n, -rad, rad, rad, smin)
    sec = time.perf_counter() - t0
    xy = pos.reshape(-1, 2)
    obj, ideal = farm_objective(xy, ds)
    curve = ideal - np.array(tr.curve)
    return dict(Algorithm=alg, Dataset=ds, Radius=rad, Turbines=n, Seed=seed, Calls=tr.calls,
                Objective=obj, Ideal=ideal, WakeLoss=ideal - obj, Feasible=bool(feasible(pos)),
                MinSpacing=min_spacing(xy), Seconds=sec,
                Coordinates=";".join(f"{a:.3f} {b:.3f}" for a, b in xy),
                Curve=";".join("nan" if not np.isfinite(c) else f"{c:.3f}" for c in curve))


def run_hr(task):
    alg, n, seed = task
    f, ideal, poly = hr.make_objective(n)
    smin = 4 * hr.D
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
    pos, tr = run_method(alg, seed, f, wake, feasible, 2 * n, -half, half, half, smin, bcons)
    sec = time.perf_counter() - t0
    xy = pos.reshape(-1, 2)
    curve = ideal - np.array(tr.curve)
    return dict(Algorithm=alg, Turbines=n, Seed=seed, Calls=tr.calls, AEP=hr.aep_gwh(xy),
                IdealAEP=ideal, Feasible=bool(feasible(pos)), MinSpacing=min_spacing(xy),
                Seconds=sec, Coordinates=";".join(f"{u:.2f} {v:.2f}" for u, v in xy),
                Curve=";".join("nan" if not np.isfinite(c) else f"{c:.4f}" for c in curve))


if __name__ == "__main__":
    exp = sys.argv[1]
    SPLITCASES = [(ds, r, n) for ds in (1, 2) for r, n in ((500, 6), (750, 8), (1000, 10), (500, 10), (750, 12), (1000, 15))]
    if exp == "grid":
        fn, tl = run_grid, [(a, *c, s) for c in GRID for a in ALGS for s in range(1, 31)]
    elif exp == "hgrid":
        fn, tl = run_grid, [("LXVNS", *c, s) for c in GRID for s in range(1, 31)]
    elif exp == "hsplit":
        fn, tl = run_grid, [(a, *c, s) for c in SPLITCASES for a in ("LXVNS25", "LXVNS75") for s in range(1, 31)]
    elif exp == "vgrid":
        fn, tl = run_grid, [("BVNS", *c, s) for c in GRID for s in range(1, 31)]
    elif exp == "vhr16":
        fn, tl = run_hr, [("BVNS", 16, s) for s in range(1, 31)]
    elif exp == "vhr80":
        fn, tl = run_hr, [("BVNS", 80, s) for s in range(1, 11)]
    elif exp == "hhr16":
        fn, tl = run_hr, [("LXVNS", 16, s) for s in range(1, 31)]
    elif exp == "hhr80":
        fn, tl = run_hr, [("LXVNS", 80, s) for s in range(1, 11)]
    elif exp == "hr16":
        fn, tl = run_hr, [(a, 16, s) for a in ALGS for s in range(1, 31)]
    else:
        fn, tl = run_hr, [(a, 80, s) for a in ALGS for s in range(1, 11)]
    t0 = time.time()
    with Pool(4) as pool:
        out = pd.DataFrame(pool.map(fn, tl, chunksize=1))
    out.to_csv(f"fresh_{exp}.csv", index=False)
    print(exp, len(out), "runs in", round(time.time() - t0), "s; calls",
          out.Calls.min(), out.Calls.max(), flush=True)

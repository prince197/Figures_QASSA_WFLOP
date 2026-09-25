"""R3-3 and R1-9 experiments.

  extra     : VNS and multistart SLSQP on the six representative cases (budgets 3,030 and 6,030)
              and the six crowded cases (6,030); same seeds / initial populations as
              authors_experiments.py, so results pair with authors_runs_budget/crowded.csv.
  hr16/hr80 : Horns Rev 1 with the real V80 turbine, wind climate and farm outline, all six
              methods at 6,030 calls (16-turbine block: 30 seeds; full 80-turbine farm: 10 seeds).
"""
import sys, time
import numpy as np, pandas as pd
from multiprocessing import Pool
from authors_optimizers import SSA, LXSSA, PSO, DE
from extra_baselines import VNS, MSSLSQP
from authors_objective import make_objective
from wflop_model import farm_objective, min_spacing, R
import hornsrev_model as hr

REP = [(1, 500, 4), (1, 750, 8), (1, 1000, 8), (2, 500, 4), (2, 750, 8), (2, 1000, 8)]
CROWDED = [(ds, r, n) for ds in (1, 2) for r, n in ((500, 10), (750, 12), (1000, 15))]
ITERS = {3030: {"LXSSA": 50, "SSA": 100, "PSO": 100, "DE": 100},
         6030: {"LXSSA": 100, "SSA": 200, "PSO": 200, "DE": 200}}


def optimise(alg, budget, seed, f, wake, dim, lb, ub, radius, smin, bcons=None):
    calls = [0]

    def counted(x):
        calls[0] += 1
        return f(x)
    if alg == "VNS":
        pos, _, _ = VNS(30, budget, radius, seed=seed).optimize(counted, dim, lb, ub)
    elif alg == "SLSQP":
        opt = MSSLSQP(wake, radius, smin, 30, budget, seed=seed, boundary_cons=bcons)
        pos, _, _ = opt.optimize(counted, dim, lb, ub)
        calls[0] = opt.calls
    else:
        it = ITERS[budget][alg]
        o = LXSSA(30, it, 0.0, 1.0, seed=seed) if alg == "LXSSA" else \
            {"SSA": SSA, "PSO": PSO, "DE": DE}[alg](30, it, seed=seed)
        pos, _, _ = o.optimize(counted, dim, lb, ub)
    return pos, calls[0]


def run_synthetic(task):
    exp, alg, budget, ds, rad, n, seed = task
    smin = 8 * R
    f = make_objective(ds, rad)

    def wake(x):
        o, i = farm_objective(x.reshape(-1, 2), ds)
        return i - o
    t0 = time.perf_counter()
    pos, calls = optimise(alg, budget, seed, f, wake, 2 * n, -rad, rad, rad, smin)
    sec = time.perf_counter() - t0
    xy = pos.reshape(-1, 2)
    obj, ideal = farm_objective(xy, ds)
    feas = (min_spacing(xy) >= smin - 1e-6) and (np.sqrt((xy ** 2).sum(1)).max() <= rad + 1e-6)
    return dict(Experiment=exp, Algorithm=alg, Budget=budget, Calls=calls, Spacing="4D", Dataset=ds,
                Radius=rad, Turbines=n, Seed=seed, Objective=obj, Ideal=ideal, WakeLoss=ideal - obj,
                Feasible=feas, MinSpacing=min_spacing(xy), Seconds=sec,
                Coordinates=";".join(f"{a:.4f} {b:.4f}" for a, b in xy))


def run_hornsrev(task):
    exp, alg, budget, n, seed = task
    f, ideal, poly = hr.make_objective(n)
    smin = 4 * hr.D
    half = np.abs(poly).max()

    def wake(x):
        return ideal - hr.aep_gwh(x.reshape(-1, 2))

    a, b = poly, np.roll(poly, -1, axis=0)
    e = b - a
    nrm = np.c_[e[:, 1], -e[:, 0]] / np.linalg.norm(e, axis=1)[:, None]

    def bcons(p):   # >= 0 inside the (convex) outline, scaled to O(1)
        return (-np.einsum("ijk,jk->ij", p[:, None, :] - a[None], nrm) / half).ravel()
    t0 = time.perf_counter()
    pos, calls = optimise(alg, budget, seed, f, wake, 2 * n, -half, half, half, smin, bcons)
    sec = time.perf_counter() - t0
    xy = pos.reshape(-1, 2)
    aep = hr.aep_gwh(xy)
    feas = (min_spacing(xy) >= smin - 1e-6) and (hr.outside_distance(xy, poly).max() <= 1e-6)
    return dict(Experiment=exp, Algorithm=alg, Budget=budget, Calls=calls, Turbines=n, Seed=seed,
                AEP=aep, IdealAEP=ideal, Feasible=feas, MinSpacing=min_spacing(xy), Seconds=sec,
                Coordinates=";".join(f"{u:.2f} {v:.2f}" for u, v in xy))


def tasks(exp):
    seeds = range(1, 31)
    if exp == "extra":
        t = [(exp, a, b, *c, s) for b in (3030, 6030) for c in REP for a in ("VNS", "SLSQP") for s in seeds]
        t += [(exp, a, 6030, *c, s) for c in CROWDED for a in ("VNS", "SLSQP") for s in seeds]
        return run_synthetic, t
    algs = ["LXSSA", "SSA", "PSO", "DE", "VNS", "SLSQP"]
    if exp == "hr16":
        return run_hornsrev, [(exp, a, 6030, 16, s) for a in algs for s in seeds]
    if exp == "hr80":
        return run_hornsrev, [(exp, a, 6030, 80, s) for a in algs for s in range(1, 11)]
    raise ValueError(exp)


if __name__ == "__main__":
    exp = sys.argv[1]
    fn, tl = tasks(exp)
    t0 = time.time()
    with Pool(4) as pool:
        out = pd.DataFrame(pool.map(fn, tl, chunksize=1))
    out.to_csv(f"extra_runs_{exp}.csv", index=False)
    print(exp, len(out), "runs in", round(time.time() - t0), "s; calls:",
          out.groupby("Algorithm").Calls.agg(["min", "max"]).to_dict())

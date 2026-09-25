"""Controlled re-runs with the authors' optimizer code (authors_optimizers.py) on the validated
benchmark evaluator.

Design
  * SSA, LX-SSA (phi = 0, chi = 1), PSO, DE; population 30; seeds 1..30.
  * The same seed gives every algorithm the same initial population (each optimizer seeds NumPy
    in its constructor and draws the initial population first), so runs are paired by seed.
  * Budgets are matched in objective-function calls: 3,030 (LX-SSA 50 iterations, others 100)
    and 6,030 (LX-SSA 100 iterations, others 200). Calls are counted, not assumed.
  * Each run is timed individually (4 worker processes in parallel on the same machine).

Experiments (argument 1)
  budget   : six representative cases, both datasets, both budgets
  crowded  : largest tabulated N per radius (500 m/10, 750 m/12, 1000 m/15), budget 6,030
  spacing  : SSA and LX-SSA at 5D and 6D, six representative cases, budget 6,030
  hornsrev : measured Horns Rev 1 wind climate (dataset 3), budget 6,030
"""
import sys, time
import numpy as np, pandas as pd
from multiprocessing import Pool
from authors_optimizers import SSA, LXSSA, PSO, DE
from authors_objective import make_objective
from wflop_model import farm_objective, min_spacing, R

POP = 30
REP = [(1, 500, 4), (1, 750, 8), (1, 1000, 8), (2, 500, 4), (2, 750, 8), (2, 1000, 8)]
CROWDED = [(ds, r, n) for ds in (1, 2) for r, n in ((500, 10), (750, 12), (1000, 15))]
HORNSREV = [(3, 500, 4), (3, 750, 8), (3, 1000, 8), (3, 1000, 15)]
ITERS = {3030: {"LXSSA": 50, "SSA": 100, "PSO": 100, "DE": 100},
         6030: {"LXSSA": 100, "SSA": 200, "PSO": 200, "DE": 200}}


def build(alg, iters, seed):
    if alg == "LXSSA":
        return LXSSA(POP, iters, phi=0.0, chi=1.0, seed=seed)
    return {"SSA": SSA, "PSO": PSO, "DE": DE}[alg](POP, iters, seed=seed)


def run(task):
    exp, alg, budget, spacing, ds, rad, n, seed = task
    smin = {"4D": 8 * R, "5D": 10 * R, "6D": 12 * R}[spacing]
    f = make_objective(ds, rad, smin=smin, penalty="linear")
    calls = [0]

    def counted(x):
        calls[0] += 1
        return f(x)
    t0 = time.perf_counter()
    pos, score, curve = build(alg, ITERS[budget][alg], seed).optimize(counted, 2 * n, -rad, rad)
    sec = time.perf_counter() - t0
    xy = pos.reshape(-1, 2)
    obj, ideal = farm_objective(xy, ds)
    rr = np.sqrt((xy ** 2).sum(1)).max()
    feas = (min_spacing(xy) >= smin - 1e-6) and (rr <= rad + 1e-6)
    return dict(Experiment=exp, Algorithm=alg, Budget=budget, Calls=calls[0], Spacing=spacing,
                Dataset=ds, Radius=rad, Turbines=n, Seed=seed, Objective=obj, Ideal=ideal,
                WakeLoss=ideal - obj, Feasible=feas, MinSpacing=min_spacing(xy), Seconds=sec,
                Coordinates=";".join(f"{a:.4f} {b:.4f}" for a, b in xy))


def tasks(exp):
    algs = ["LXSSA", "SSA", "PSO", "DE"]
    seeds = range(1, 31)
    if exp == "budget":
        return [(exp, a, b, "4D", *c, s) for b in (3030, 6030) for c in REP for a in algs for s in seeds]
    if exp == "crowded":
        return [(exp, a, 6030, "4D", *c, s) for c in CROWDED for a in algs for s in seeds]
    if exp == "spacing":
        return [(exp, a, 6030, sp, *c, s) for sp in ("5D", "6D") for c in REP
                for a in ("LXSSA", "SSA") for s in seeds]
    if exp == "hornsrev":
        return [(exp, a, 6030, "4D", *c, s) for c in HORNSREV for a in algs for s in seeds]
    raise ValueError(exp)


if __name__ == "__main__":
    exp = sys.argv[1]
    t0 = time.time()
    with Pool(4) as pool:
        out = pd.DataFrame(pool.map(run, tasks(exp), chunksize=2))
    out.to_csv(f"authors_runs_{exp}.csv", index=False)
    print(exp, len(out), "runs in", round(time.time() - t0), "s; calls by budget/alg:",
          out.groupby(["Budget", "Algorithm"]).Calls.unique().to_dict())

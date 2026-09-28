"""IEA Wind Task 37 WFLO benchmark runs (Case Study 1, 16- and 36-turbine scenarios; iea37_model.py).

Usage:  python3 iea37_experiments.py EXP [SHARD NSHARDS] [--procs=P]
Output: mpce_<EXP>_s<SHARD>of<NSHARDS>.csv, same columns as run_grid in mpce_experiments.py.

Experiments
  iea16   16 turbines, circle r = 1300 m, 9 methods (M9), 30 seeds, 6,030 and 30,030 calls, random init
  iea36   36 turbines, circle r = 2000 m, same design
Minimum spacing 2D = 260 m. Objective / Ideal / WakeLoss / Curve are in MWh (Objective = AEP of the
returned layout computed exactly as the official iea37-aepcalc.py; Ideal = wake-free AEP).
The minimised function is (Ideal - AEP) + the authors' quadratic penalties (iea37_model.make_objective).
"""
import sys, time
import numpy as np, pandas as pd
from multiprocessing import Pool
import init_hook
from mpce_experiments import run_method, M9
import iea37_model as iea

CASES = {"iea16": 16, "iea36": 36}
# Dataset labels requested by the lead. NB: officially both scenarios belong to IEA37 "Case Study 1"
# (optimization only; 16/36/64 turbines); the official "Case Study 2" is a 9-turbine combined study.
DATASET = {16: "IEA37cs1", 36: "IEA37cs2"}
BUDGETS = (30030, 6030)
SEEDS = range(1, 31)
# Heaviest first, from seed-1 timings of the 36-turbine case at 30,030 calls on 4 busy cores:
# SLSQP 23.8 s, DE 20.5, SSA 16.3, LXSSA 15.6, RSVNS 15.3, SSABV 15.0, PSOC 14.6, LXBV 14.2, BVNS 13.2.
COST_ORDER = ["SLSQP", "DE", "SSA", "LXSSA", "RSVNS", "SSABV", "PSOC", "LXBV", "BVNS"]


def run_iea(task):
    alg, n, seed, budget, init = task
    f, ideal, radius, smin = iea.make_objective(n)
    init_hook.GEN = None                    # random initialization (uniform in the bounding box)

    def wake(x):
        return ideal - iea.aep(x.reshape(-1, 2))

    def feasible(x):
        xy = x.reshape(-1, 2)
        return iea.min_spacing(xy) >= smin - 1e-6 and np.sqrt((xy ** 2).sum(1)).max() <= radius + 1e-6
    t0 = time.perf_counter()
    pos, tr = run_method(alg, seed, budget, f, wake, feasible, 2 * n, -radius, radius, radius, smin)
    sec = time.perf_counter() - t0
    xy = pos.reshape(-1, 2)
    obj = iea.aep(xy)
    curve = ideal - np.array(tr.curve)
    return dict(Algorithm=alg, Dataset=DATASET[n], Radius=radius, Turbines=n, Seed=seed, Budget=budget, Init=init,
                Calls=tr.calls, Objective=obj, Ideal=ideal, WakeLoss=ideal - obj, Feasible=bool(feasible(pos)),
                MinSpacing=iea.min_spacing(xy), Seconds=sec,
                Coordinates=";".join(f"{a:.3f} {b:.3f}" for a, b in xy),
                Curve=";".join("nan" if not np.isfinite(c) else f"{c:.3f}" for c in curve))


def tasks(exp):
    """Heaviest first: larger budget first, then methods in measured cost order, then seeds."""
    n = CASES[exp]
    assert sorted(COST_ORDER) == sorted(M9)
    return [(alg, n, s, b, "random") for b in BUDGETS for alg in COST_ORDER for s in SEEDS]


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
        rows = pool.map(run_iea, tl, chunksize=1)
    out = pd.DataFrame(rows)
    fn = f"mpce_{exp}_s{shard}of{nsh}.csv"
    out.to_csv(fn, index=False)
    print(exp, fn, len(out), "runs in", round(time.time() - t0), "s; calls", out.Calls.min(), out.Calls.max(), flush=True)

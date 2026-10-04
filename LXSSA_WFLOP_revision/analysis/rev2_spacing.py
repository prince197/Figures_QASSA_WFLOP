"""Review round 2, reviewer R1.2 (the 4D minimum spacing is untested): re-optimization at 5D and 6D.

Usage:  python3 rev2_spacing.py EXP [SHARD NSHARDS] [--procs P | --procs=P]
Output: rev2_<EXP>_s<SHARD>of<NSHARDS>.csv (columns of mpce_experiments.run_grid plus Spacing, SMin)

Experiments
  spacing   the 12 split cases (mpce_experiments.SPLITCASES), 30 seeds, 6,030 calls, random initialization,
            minimum spacing ell_min = 5D (= 10R = 385 m) and 6D (= 12R = 462 m), 9 methods: the 8 of the main
            comparison (PSOBV = PSO-VNS, PSOC = PSO, SSABV = SSA-VNS, SSA, LXSSA, DE, BVNS = VNS,
            SLSQP = MS-SLSQP) and the disc-sampling control RSDVNS (RSD-VNS);
            9 x 12 x 2 x 30 = 6,480 runs, heaviest first (largest N first; SLSQP first within a case)
  spacing4d verification only (not part of the design): the same tasks at 4D; with Spacing = "4D" every
            run is identical to mpce_experiments.run_grid / run_grid_x (e.g. PSOC, DS 1, 1000 m, N = 15,
            seed 1 -> WakeLoss 4380.629934684286)

How the spacing is applied (everything that uses 8R / 4D in the 68-case code path, now smin = SMIN[spacing]):
  * penalty:     authors_objective.make_objective(ds, rad, smin=smin)
  * feasibility: feasible() (Tracker best-so-far / convergence curve and the Feasible column), min_spacing >= smin - 1e-6
  * MS-SLSQP:    run_method passes smin to extra_baselines.MSSLSQP (explicit constraints d^2 >= smin^2 and its
                 internal feasible-best bookkeeping)
  * feasible initialization would use make_generator(smin, ...) (not used here: Init = "random").
Nothing else depends on the spacing: BVNS (original_vns), RS-VNS / RSD-VNS (rs_vns, mpce_experiments.RSDVNS),
the hybrids (hybrid_lxssa_bvns.HybridBVNS) and SSA / LX-SSA / PSO / DE (authors_optimizers) see the constraint
only through the penalized objective; their step sizes (rho_k = 0.1 k r, h0 = 0.05 r, h_min = 1e-3 r) scale with the
farm radius only. The method code paths are exactly mpce_experiments.run_method / run_method_x (unchanged).

Geometric capacity (packing_capacity.csv, constructive lower bounds that coincide with the proven optimal
circle packings for r = 500 m): 500 m holds at most 8 turbines at 5D and 7 at 6D, so the split cases
(DS 1 and 2, 500 m, N = 10) have NO feasible layout at 5D or at 6D (4 of the 24 case x spacing cells).
All other split cases are feasible (750 m: 19 at 5D, 13 at 6D >= 12; 1000 m: 29 / 21 >= 15). They are still
run; the ranking is feasibility-aware (infeasible runs are ranked by constraint violation / MinSpacing).
"""
from record_io import encode_coordinates

import sys, time
import numpy as np, pandas as pd
from multiprocessing import Pool
import init_hook
import mpce_experiments as E
from feasible_init import make_generator
from authors_objective import make_objective
from wflop_model import farm_objective, min_spacing, R

SMIN = {"4D": 8 * R, "5D": 10 * R, "6D": 12 * R}
METHODS = ["SLSQP", "PSOBV", "SSABV", "RSDVNS", "BVNS", "PSOC", "SSA", "LXSSA", "DE"]
CAPACITY = {(500, "4D"): 13, (500, "5D"): 8, (500, "6D"): 7, (750, "4D"): 26, (750, "5D"): 19, (750, "6D"): 13,
            (1000, "4D"): 42, (1000, "5D"): 29, (1000, "6D"): 21}      # packing_capacity.csv (MaxNConstructed)


def run_grid_sp(task):
    """mpce_experiments.run_grid with the minimum spacing as a task parameter (smin = SMIN[spacing]) and the
    Phase-6 dispatcher run_method_x (RSDVNS; every other label falls through to run_method unchanged)."""
    alg, ds, rad, n, seed, budget, init, spacing = task
    smin = SMIN[spacing]
    init_hook.GEN = make_generator(smin, circle_r=rad) if init == "feasible" else None
    f = make_objective(ds, rad, smin=smin)

    def wake(x):
        o, i = farm_objective(x.reshape(-1, 2), ds)
        return i - o

    def feasible(x):
        xy = x.reshape(-1, 2)
        return (n == 1 or min_spacing(xy) >= smin - 1e-6) and np.sqrt((xy ** 2).sum(1)).max() <= rad + 1e-6
    t0 = time.perf_counter()
    pos, tr = E.run_method_x(alg, seed, budget, f, wake, feasible, 2 * n, -rad, rad, rad, smin)
    sec = time.perf_counter() - t0
    xy = pos.reshape(-1, 2)
    obj, ideal = farm_objective(xy, ds)
    curve = ideal - np.array(tr.curve)
    return dict(Algorithm=alg, Dataset=ds, Radius=rad, Turbines=n, Seed=seed, Budget=budget, Init=init,
                Calls=tr.calls, Objective=obj, Ideal=ideal, WakeLoss=ideal - obj, Feasible=bool(feasible(pos)),
                MinSpacing=min_spacing(xy), Seconds=sec,
                Coordinates=encode_coordinates(xy),
                Curve=";".join("nan" if not np.isfinite(c) else f"{c:.3f}" for c in curve),
                Spacing=spacing, SMin=smin)


def infeasible_cells(spacings=("5D", "6D")):
    """(Dataset, Radius, N, Spacing) of the split cases whose N exceeds the constructed packing capacity."""
    return [(ds, r, n, sp) for sp in spacings for ds, r, n in E.SPLITCASES if n > CAPACITY[(r, sp)]]


def tasks(exp):
    S30 = range(1, 31)
    if exp in ("spacing", "spacing4d"):
        sps = ("6D", "5D") if exp == "spacing" else ("4D",)
        cases = sorted(E.SPLITCASES, key=lambda c: (-c[2], -c[1], c[0]))          # largest N first
        return [(run_grid_sp, (a, *c, s, 6030, "random", sp)) for c in cases for sp in sps for a in METHODS
                for s in S30]
    raise ValueError(exp)


def _call(pair):
    fn, t = pair
    return fn(t)


if __name__ == "__main__":
    argv = sys.argv[1:]
    procs, args, i = 2, [], 0
    while i < len(argv):
        a = argv[i]
        if a.startswith("--procs="):
            procs = int(a.split("=")[1])
        elif a == "--procs":
            procs = int(argv[i + 1]); i += 1
        else:
            args.append(a)
        i += 1
    exp = args[0]
    shard, nsh = (int(args[1]), int(args[2])) if len(args) >= 3 else (0, 1)
    tl = tasks(exp)[shard::nsh]
    print(exp, "shard", shard, "of", nsh, ":", len(tl), "runs; split cases without a feasible layout:",
          infeasible_cells(), flush=True)
    t0 = time.time()
    with Pool(procs) as pool:
        rows = pool.map(_call, tl, chunksize=1)
    out = pd.DataFrame(rows)
    fn = f"rev2_{exp}_s{shard}of{nsh}.csv"
    out.to_csv(fn, index=False)
    print(exp, fn, len(out), "runs in", round(time.time() - t0), "s; calls", out.Calls.min(), out.Calls.max(),
          flush=True)

"""Revision 2 (R1.9, AE.8): second measured-wind site, the Lillgrund 16-turbine block (rev2_site_model.py).

Usage:  python3 rev2_site.py EXP SHARD NSHARDS [--procs=P]
Output: rev2_<EXP>_s<SHARD>of<NSHARDS>.csv with the columns of mpce_experiments.run_hr (Dataset = "LG") plus Spacing.

Experiments (random initialization, seeds 1..30, N_p = 30, every objective call counted):
  lg16     9 methods x 30 seeds at 6,030 evaluations   (270 runs)
  lg16b    9 methods x 30 seeds at 30,030 evaluations  (270 runs)
  lg16all  lg16b followed by lg16 (540 runs, heaviest first)
Methods: PSOBV (PSO-VNS), PSOC, SSABV (SSA-VNS), SSA, LXSSA, DE, BVNS, SLSQP (MS-SLSQP) -- the eight methods of the main
comparison, run through mpce_experiments.run_method unchanged -- and RSDVNS, the random-sampling control whose Phase-1
samples are uniform in the SITE (here the block parallelogram, sampled exactly by an affine map of two uniform draws
per turbine; for the circular farms of experiment rsdisc the site is the disc). As in mpce_experiments.RSDVNS the
seeded initial population of 30 comes from init_pop, and Phase 2 is the unchanged basic VNS.

The runner run_lg is mpce_experiments.run_hr with the Lillgrund model in place of hornsrev_model: search box of
half-width max|poly| centred on the block (also the r of the VNS radii), the same boundary constraints for MS-SLSQP,
feasibility = min spacing >= l_min - 1e-6 and max outside distance <= 1e-6. The minimum spacing is
rev2_site_model.SMIN = 3 D (4 D is geometrically impossible for 16 turbines in this block; see rev2_site_model).
"""
from record_io import encode_coordinates

import sys, time
import numpy as np, pandas as pd
from multiprocessing import Pool
import init_hook
import mpce_experiments as mx
from mpce_experiments import NP, Tracker
from wflop_model import min_spacing
import rev2_site_model as lg

METHODS = ["PSOBV", "PSOC", "SSABV", "SSA", "LXSSA", "DE", "BVNS", "SLSQP", "RSDVNS"]


class RSPVNS(mx.RSDVNS):
    """mpce_experiments.RSDVNS with Phase-1 samples uniform in a parallelogram site (vertices p0, p1, p2, p3):
    x = p0 + s (p1 - p0) + t (p3 - p0), s, t ~ U(0, 1) per turbine."""

    def __init__(self, pop_size, budget, split, radius, poly, seed=None):
        super().__init__(pop_size, budget, split, radius, seed=seed)
        self.p0, self.e1, self.e2 = poly[0], poly[1] - poly[0], poly[3] - poly[0]
        assert np.allclose(poly[2], poly[0] + self.e1 + self.e2), "site is not a parallelogram"

    def optimize(self, obj_fun, dim, lb, ub):
        from original_vns import BVNS as _BVNS, BudgetExhausted as _BE
        from init_hook import init_pop as _init_pop
        calls = [0]

        def f(x):
            if calls[0] >= self.budget:
                raise _BE
            calls[0] += 1
            return obj_fun(x)

        pop = _init_pop(self.pop_size, dim, lb, ub)
        fit = np.array([f(p) for p in pop])
        best, fbest = pop[np.argmin(fit)].copy(), fit.min()
        n = dim // 2
        while calls[0] < self.n1:
            s = np.random.uniform(0.0, 1.0, n)
            t = np.random.uniform(0.0, 1.0, n)
            x = (self.p0[None] + s[:, None] * self.e1[None] + t[:, None] * self.e2[None]).ravel()
            fx = f(x)
            if fx < fbest:
                best, fbest = x.copy(), fx
        self.phase1_calls = calls[0]
        return _BVNS(self.pop_size, self.budget, self.radius).search(f, best, fbest, dim, lb, ub)


def run_method_lg(alg, seed, budget, f, wake, feasible, dim, lb, ub, radius, smin, bcons, poly):
    if alg != "RSDVNS":
        return mx.run_method(alg, seed, budget, f, wake, feasible, dim, lb, ub, radius, smin, bcons)
    tr = Tracker(feasible, max(1, (budget - NP) // 200))

    def fp(x):
        v = f(x); tr.see(x, v); return v
    pos, _, _ = RSPVNS(NP, budget, 0.5, radius, poly, seed=seed).optimize(fp, dim, lb, ub)
    return pos, tr


def run_lg(task):
    alg, n, seed, budget, init = task
    if init != "random":
        raise ValueError("only random initialization is designed for this site")
    smin = lg.SMIN
    f, ideal, poly = lg.make_objective(n, smin)
    init_hook.GEN = None
    half = np.abs(poly).max()
    a, b = poly, np.roll(poly, -1, axis=0)
    e = b - a
    nrm = np.c_[e[:, 1], -e[:, 0]] / np.linalg.norm(e, axis=1)[:, None]

    def wake(x):
        return ideal - lg.aep_gwh(x.reshape(-1, 2))

    def bcons(p):
        return (-np.einsum("ijk,jk->ij", p[:, None, :] - a[None], nrm) / half).ravel()

    def feasible(x):
        xy = x.reshape(-1, 2)
        return min_spacing(xy) >= smin - 1e-6 and lg.outside_distance(xy, poly).max() <= 1e-6
    t0 = time.perf_counter()
    pos, tr = run_method_lg(alg, seed, budget, f, wake, feasible, 2 * n, -half, half, half, smin, bcons, poly)
    sec = time.perf_counter() - t0
    xy = pos.reshape(-1, 2)
    curve = ideal - np.array(tr.curve)
    aep = lg.aep_gwh(xy)
    return dict(Algorithm=alg, Dataset="LG", Radius=0, Turbines=n, Seed=seed, Budget=budget, Init=init,
                Calls=tr.calls, Objective=aep, Ideal=ideal, WakeLoss=ideal - aep,
                Feasible=bool(feasible(pos)), MinSpacing=min_spacing(xy), Seconds=sec,
                Coordinates=encode_coordinates(xy),
                Curve=";".join("nan" if not np.isfinite(c) else f"{c:.4f}" for c in curve), Spacing=smin)


def tasks(exp):
    S30 = range(1, 31)
    if exp == "lg16":
        return [(run_lg, (a, 16, s, 6030, "random")) for a in METHODS for s in S30]
    if exp == "lg16b":
        return [(run_lg, (a, 16, s, 30030, "random")) for a in METHODS for s in S30]
    if exp == "lg16all":
        return tasks("lg16b") + tasks("lg16")
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
        elif a.startswith("--procs"):
            raise SystemExit("use --procs=P")
    exp = args[0]
    shard, nsh = (int(args[1]), int(args[2])) if len(args) >= 3 else (0, 1)
    tl = tasks(exp)[shard::nsh]
    t0 = time.time()
    with Pool(procs) as pool:
        rows = pool.map(_call, tl, chunksize=1)
    out = pd.DataFrame(rows)
    fn = f"rev2_{exp}_s{shard}of{nsh}.csv"
    out.to_csv(fn, index=False)
    print(exp, fn, len(out), "runs in", round(time.time() - t0), "s; calls", out.Calls.min(), out.Calls.max(), flush=True)

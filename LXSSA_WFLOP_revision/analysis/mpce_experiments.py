"""Additional experiments of the MPCE study (SSA-VNS study).

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
  feas     design of the feasibility-preserving-initialization study (6,030 calls, 6 largest cases + HR16): the
           9 methods of M9. NOT run as such: RS-VNS would need one packing solve per random sample. The data
           were produced by `feasx` (the same design without RS-VNS, 8 methods, 8 shards), whose shards were
           merged into mpce_feas_s0of1.csv (the file mpce_results.py reads); PSO-VNS comes from `feasp`.
  feasx    feas without RS-VNS (8 methods) -- the experiment actually run
  b30k     random initialization, 30,030 calls, 9 methods, 6 largest cases + HR16 (30 seeds)
  b120k    random initialization, 120,030 calls, 9 methods, 6 largest cases (30 seeds) + HR16 (10 seeds)
  feasp / b30kp / b120kp   the PSO-VNS arm of feas / b30k / b120k;  psosplit  PSO-VNS with 25 % / 75 % split
  hrfix    all Horns Rev 1 16-turbine runs again after the direction-binning fix of hornsrev_model (2026-09-28):
           10 methods; 6,030 random / 6,030 feasible (no RS-VNS) / 30,030 random (30 seeds), 120,030 (10 seeds)
  omega90  (additional) PSO-VNS with split 0.9 (PSOBV90) on the 12 split cases, 30 seeds, 6,030 calls
  rsdisc   (additional) RS-VNS with Phase-1 samples uniform in the farm disc (RSDVNS), 68 cases, 30 seeds, 6,030 calls
  csweep   (sensitivity study) stand-alone PSO on the 12 split cases, 30 seeds, 6,030 calls: w = 0.7 with
           c1 = c2 = c in {1.2, 1.4, 1.6, 1.7, 1.8, 1.9, 2.0} (PSOW07C12 ... PSOW07C20), the old setting with velocity
           zeroed on clip (PSOOLD_VZERO) or clamped to 0.2 (ub - lb) (PSOOLD_VMAX); 3,240 runs; see CSWEEP
Convergence curves: one checkpoint every (B-30)//200 calls, i.e. 201 checkpoints at 6,030 calls and 200 at
30,030 / 120,030 calls (the first checkpoint is at call 30 only for B = 6,030).
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
    # ---- additional controls: new labels only, run through run_grid_x (see below) ----
    if exp == "omega90":
        # PSO-VNS with split omega = 0.9 (label PSOBV90) on the 12 split cases, 30 seeds, 6,030 calls: completes
        # the split curve 0.25 / 0.5 / 0.75 / 0.9 / 1 (= PSO) of the psosplit design
        return [(run_grid_x, ("PSOBV90", *c, s, 6030, "random")) for c in SPLITCASES for s in S30]
    if exp == "rsdisc":
        # RS-VNS with Phase-1 samples uniform in the farm disc instead of the bounding square (label RSDVNS),
        # 68 cases, 30 seeds, 6,030 calls; heaviest cases first
        return [(run_grid_x, ("RSDVNS", *c, s, 6030, "random")) for c in GRID[::-1] for s in S30]
    # ---- sensitivity study: PSO coefficient sweep and bound handling, run through run_grid_cs ----
    if exp == "csweep":
        # stand-alone PSO on the 12 split cases (psosplit design), 30 seeds, 6,030 calls, random starts:
        # w = 0.7, c1 = c2 = c in CS_C (PSOW07C12 ... PSOW07C20; c = 2 = old setting), the old setting with
        # velocity zeroed on clip (PSOOLD_VZERO) and with velocity clamping |v| <= 0.2 (ub - lb) (PSOOLD_VMAX);
        # 9 labels x 12 cases x 30 seeds = 3,240 runs. The constriction reference is the stored PSOC data
        # (mpce_psoc); label PSOCREF (not in this list) reruns it through PSOBH for verification only
        return [(run_grid_cs, (a, *c, s, 6030, "random")) for c in SPLITCASES for a in CSWEEP if a != "PSOCREF"
                for s in S30]
    raise ValueError(exp)


# ---------------------------------------------------------------------------------------------------------
_RUN_METHOD = run_method          # the original function (run_grid_x swaps the module-level name)
# Additional controls. Nothing above is changed: the existing labels keep their code path (run_method /
# run_grid); the two new labels are handled by run_method_x, and run_grid_x runs the unchanged run_grid with
# run_method_x in place of run_method (every other label falls through to run_method unchanged).
class RSDVNS(RSVNS):
    """RS-VNS whose Phase-1 samples are uniform in the farm disc (polar sampling: angle U(0, 2 pi), radius
    r sqrt(U(0, 1)) per turbine) instead of the bounding square, a stronger random-sampling control.
    The seeded initial population of 30 is the common one drawn by init_pop (so runs stay seed-paired with
    all other methods); the remaining round(split B) - 30 Phase-1 samples are disc samples. Every turbine of
    a disc sample lies inside the site, so only the spacing constraint can be violated. Phase 2 is the
    unchanged basic VNS (original_vns.BVNS) from the best sample (lowest F_p)."""

    def __init__(self, pop_size=30, budget=6030, split=0.5, radius=500.0, disc_r=None, seed=None):
        super().__init__(pop_size, budget, split, radius, seed=seed)
        self.disc_r = radius if disc_r is None else disc_r

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
            ang = np.random.uniform(0.0, 2 * np.pi, n)
            rad = self.disc_r * np.sqrt(np.random.uniform(0.0, 1.0, n))
            x = np.c_[rad * np.cos(ang), rad * np.sin(ang)].ravel()
            fx = f(x)
            if fx < fbest:
                best, fbest = x.copy(), fx
        self.phase1_calls = calls[0]
        return _BVNS(self.pop_size, self.budget, self.radius).search(f, best, fbest, dim, lb, ub)


def run_method_x(alg, seed, budget, f, wake, feasible, dim, lb, ub, radius, smin, bcons=None):
    """run_method for the additional labels PSOBV90 (PSO-VNS, omega = 0.9) and RSDVNS (disc-sampling RS-VNS);
    any other label is passed to run_method unchanged."""
    if alg not in ("PSOBV90", "RSDVNS"):
        return _RUN_METHOD(alg, seed, budget, f, wake, feasible, dim, lb, ub, radius, smin, bcons)
    tr = Tracker(feasible, max(1, (budget - NP) // 200))

    def fp(x):
        v = f(x); tr.see(x, v); return v
    if alg == "PSOBV90":
        pos, _, _ = HybridBVNS(NP, budget, 0.9, 0.0, 1.0, radius, "PSOC", seed=seed).optimize(fp, dim, lb, ub)
    else:
        pos, _, _ = RSDVNS(NP, budget, 0.5, radius, seed=seed).optimize(fp, dim, lb, ub)
    return pos, tr


def run_grid_x(task):
    """run_grid (unchanged) with run_method_x in place of run_method (one task at a time per process)."""
    g = globals()
    orig = g["run_method"]
    g["run_method"] = run_method_x
    try:
        return run_grid(task)
    finally:
        g["run_method"] = orig


# ---------------------------------------------------------------------------------------------------------
# Sensitivity study: experiment csweep. Nothing above is changed (the only edit above is
# the new `if exp == "csweep"` block in tasks()). The new labels are handled by run_method_cs; run_grid_cs runs
# the unchanged run_grid with run_method_cs in place of run_method and appends the dynamics columns logged by
# PSOBH (Spread, VelMean, ClipPct, ClipLatePct, FeasEvalLatePct) after the standard columns.
CS_C = (1.2, 1.4, 1.6, 1.7, 1.8, 1.9, 2.0)
CSWEEP = {f"PSOW07C{int(round(10 * c))}": (0.7, c, "clip") for c in CS_C}      # label: (w, c1 = c2, bound)
CSWEEP["PSOOLD_VZERO"] = (0.7, 2.0, "vzero")
CSWEEP["PSOOLD_VMAX"] = (0.7, 2.0, "vmax")
CSWEEP["PSOCREF"] = (0.7298, 1.49618, "clip")
VMAX_FRAC = 0.2
_CS_LOG = {}                      # dynamics of the last PSOBH run in this process (one task at a time)


class PSOBH(PSO):
    """authors_optimizers.PSO with a choice of bound handling and passive logging.

    bound = "clip"  : position clipped to [lb, ub], velocity kept (the implemented PSO; bit-identical to PSO)
            "vzero" : position clipped, and every clipped velocity component set to zero
            "vmax"  : velocity clamped to |v| <= vmax_frac (ub - lb) per coordinate before the position update
                      (as in early PSO, Kennedy-Eberhart 1995 / Shi-Eberhart 1998), then the position clipped
    optimize() is a verbatim copy of PSO.optimize with the two bound-handling lines and [log] lines added; the
    log lines only read state and draw no random numbers, so the random stream (and the seeded initial swarm)
    is the same for every setting and bound. Logged (in self.log): final swarm spread (mean over particles of
    the RMS turbine distance to the global best, m, as in mpce_diagnostics T1), final mean |V| per coordinate,
    share of coordinates outside the box before clipping (all iterations / iterations 101-200), share of
    particle evaluations with a feasible layout in iterations 101-200 (needs `feasible`)."""

    def __init__(self, pop_size=50, max_iter=100, w=0.7, c1=2.0, c2=2.0, bound="clip", vmax_frac=VMAX_FRAC,
                 feasible=None, seed=None):
        super().__init__(pop_size, max_iter, w, c1, c2, seed=seed)
        if bound not in ("clip", "vzero", "vmax"):
            raise ValueError(bound)
        self.bound = bound; self.vmax_frac = vmax_frac; self.feasible = feasible

    def optimize(self, obj_fun, dim, lb, ub):
        from init_hook import init_pop as _init_pop
        half = self.max_iter // 2; nclip = [0, 0]; ncoord = [0, 0]; nfe = 0; nev = 0              # [log]
        vmax = self.vmax_frac * (np.asarray(ub, float) - np.asarray(lb, float))
        X = _init_pop(self.pop_size, dim, lb, ub)
        V = np.zeros((self.pop_size, dim))
        pbest = X.copy()
        pbest_score = np.array([obj_fun(x) for x in X])
        gbest_idx = np.argmin(pbest_score)
        gbest = pbest[gbest_idx].copy(); gbest_score = pbest_score[gbest_idx]
        curve = [gbest_score]
        for it in range(self.max_iter):
            late = int(it >= half)                                                                   # [log]
            for i in range(self.pop_size):
                r1 = np.random.rand(dim); r2 = np.random.rand(dim)
                V[i] = (self.w * V[i] + self.c1 * r1 * (pbest[i] - X[i])
                        + self.c2 * r2 * (gbest - X[i]))
                if self.bound == "vmax":
                    V[i] = np.clip(V[i], -vmax, vmax)
                X[i] = X[i] + V[i]
                out = (X[i] < lb) | (X[i] > ub)
                nclip[late] += int(out.sum()); ncoord[late] += dim                                   # [log]
                X[i] = np.clip(X[i], lb, ub)
                if self.bound == "vzero":
                    V[i][out] = 0.0
                score = obj_fun(X[i])
                if late and self.feasible is not None:                                               # [log]
                    nev += 1; nfe += int(bool(self.feasible(X[i].copy())))
                if score < pbest_score[i]:
                    pbest_score[i] = score; pbest[i] = X[i].copy()
                if score < gbest_score:
                    gbest_score = score; gbest = X[i].copy()
            curve.append(gbest_score)
        n = dim // 2                                                                                 # [log]
        self.log = dict(Spread=float(np.mean(np.sqrt(((X - gbest) ** 2).sum(1) / n))), VelMean=float(np.abs(V).mean()),
                        ClipPct=100.0 * sum(nclip) / max(1, sum(ncoord)), ClipLatePct=100.0 * nclip[1] / max(1, ncoord[1]),
                        FeasEvalLatePct=100.0 * nfe / nev if nev else np.nan)
        return gbest, gbest_score, np.array(curve)


def run_method_cs(alg, seed, budget, f, wake, feasible, dim, lb, ub, radius, smin, bcons=None):
    """run_method for the csweep labels (CSWEEP); any other label is passed to run_method unchanged."""
    if alg not in CSWEEP:
        return _RUN_METHOD(alg, seed, budget, f, wake, feasible, dim, lb, ub, radius, smin, bcons)
    tr = Tracker(feasible, max(1, (budget - NP) // 200))

    def fp(x):
        v = f(x); tr.see(x, v); return v
    w, c, bound = CSWEEP[alg]
    opt = PSOBH(NP, (budget - NP) // NP, w=w, c1=c, c2=c, bound=bound, feasible=feasible, seed=seed)
    pos, _, _ = opt.optimize(fp, dim, lb, ub)
    _CS_LOG.clear(); _CS_LOG.update(opt.log)
    return pos, tr


def run_grid_cs(task):
    """run_grid (unchanged) with run_method_cs in place of run_method; dynamics columns appended."""
    g = globals()
    orig = g["run_method"]
    g["run_method"] = run_method_cs
    _CS_LOG.clear()
    try:
        row = run_grid(task)
    finally:
        g["run_method"] = orig
    for k in ("Spread", "VelMean", "ClipPct", "ClipLatePct", "FeasEvalLatePct"):
        row[k] = _CS_LOG.get(k, np.nan)
    return row


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

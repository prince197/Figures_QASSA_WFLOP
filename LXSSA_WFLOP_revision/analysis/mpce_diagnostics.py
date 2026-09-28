"""Phase-6 diagnostics (workstream W3): instrumented runs that explain the mechanisms behind the component
analysis. No new optimizer; the instrumentation only observes.

Usage (from analysis/):  python3 mpce_diagnostics.py [--procs=2] [--cache=DIR]
Outputs: mpce_summary_diag.json, mpce_numbers_diag.tex (\\ND... macros), mpce_supp_diag.tex (tables/figure,
labels tab:D-*, fig:D-*), ../figures_mpce/diag_pso_dynamics.pdf (+ .png). Checks: mpce_check_diag.py.

Studies
  T1 PSO dynamics (R2 #1, #2). PSO with the old setting (w = 0.7, c1 = c2 = 2, label PSO) and the constriction
     setting (PSOC) on 6 cases x 10 seeds at 6,030 calls (200 iterations). Per iteration: swarm spread (mean
     over particles of the RMS turbine distance to the global best, m), mean |velocity| per coordinate (m),
     share of coordinates clipped to the box in that iteration (position outside [-r, r] before clipping),
     share of particle evaluations with a feasible layout / a turbine outside the circle / a spacing
     violation, and the global best. Plus, on all 2,040 stored runs of each setting: how the infeasible final
     layouts violate the constraints and how many final layouts have a coordinate on the box bound.
  T2 VNS internals (R2 #9). Phase 2 of PSO-VNS, SSA-VNS and RS-VNS on the 12 split cases x 10 seeds: shaking
     steps, accepted moves per neighbourhood k, evaluations of the local search, the improvement of the first
     descent from the switch point vs. the shake + local-search cycles, feasibility changes; also whether the
     start handed to VNS is the best feasible Phase-1 layout (question of an "rsbest" control).
  T3 Feasible starts (R2 #15). PSO (constriction) and DE with feasibility-preserving initialization on
     3 cases x 3 seeds: every candidate is logged (feasibility, spacing, comparison with the bests), the
     first update is dissected, and a counterfactual computation (separate RNG, after the run) measures how
     often the first-update PSO move between two feasible layouts is feasible with the original turbine
     labels and after relabelling one layout to match the other (label symmetry). Plus the 210 stored runs.
  T4 CPU estimates for the cloud experiments omega90 and rsdisc of mpce_experiments.py (stored run times).

Instrumentation and bit-identity. Instrumented classes are verbatim copies of the optimize()/search()
methods of authors_optimizers.PSO / DE and original_vns.BVNS with logging lines added; logging only reads
state (copies), draws no random numbers and changes no array. The instrumented classes are swapped into the
modules (mpce_experiments.PSO / .DE, hybrid_lxssa_bvns.BVNS, rs_vns.BVNS) and the runs go through the
unchanged mpce_experiments.run_method. Verification (reported in the summary, checks D01-D02): (a) every
instrumented run reproduces the stored run of the paper (same WakeLoss float, same coordinate string, same
convergence-curve string, same number of calls); (b) for seed 1 of every case and method, the same run is
repeated with the original classes in the same process and gives bit-identical final positions, curve and
call count.
"""
import os, sys, glob, json, time, pickle, argparse
from contextlib import contextmanager
import numpy as np, pandas as pd
from multiprocessing import Pool
from scipy.optimize import linear_sum_assignment

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import init_hook
from init_hook import init_pop
import authors_optimizers as AO
import original_vns as OV
from original_vns import BudgetExhausted
import hybrid_lxssa_bvns as HB
import rs_vns as RSm
import mpce_experiments as E
import hornsrev_model as hr
from feasible_init import make_generator
from authors_objective import make_objective
from wflop_model import farm_objective, min_spacing, R

SMIN = 8 * R
NPOP = 30
BUDGET = 6030
T1_CASES = [(1, 500, 6), (1, 1000, 10), (1, 750, 12), (2, 750, 8), (2, 500, 10), (2, 1000, 15)]
T1_SEEDS = range(1, 11)
T1_ALGS = ["PSO", "PSOC"]                     # old setting, constriction setting
T2_CASES = [(ds, r, n) for ds in (1, 2) for r, n in ((500, 6), (750, 8), (1000, 10), (500, 10), (750, 12), (1000, 15))]
T2_SEEDS = range(1, 11)
T2_ALGS = ["PSOBV", "SSABV", "RSVNS"]
T3_CASES = [(ds, r, n) for ds in (1, 2) for r, n in ((500, 10), (750, 12), (1000, 15))] + [("HR", 0, 16)]
T3_SEEDS = range(1, 4)
T3_ALGS = ["PSOC", "DE"]
LAB = {"PSO": "PSO (old setting)", "PSOC": "PSO (constriction)", "PSOBV": "PSO-VNS", "SSABV": "SSA-VNS",
       "RSVNS": "RS-VNS", "DE": "DE"}
MAC = {"PSOBV": "PSOVNS", "SSABV": "SSAVNS", "RSVNS": "RSVNS", "PSOC": "PSO", "DE": "DE"}
STORED = {"PSO": ("fresh_grid.csv",), "PSOC": ("mpce_psoc_s0of1.csv",), "PSOBV": ("mpce_psobv_s0of2.csv", "mpce_psobv_s1of2.csv"),
          "SSABV": ("fresh_bgrid.csv",), "RSVNS": ("mpce_rsvns_s0of1.csv",)}
STORED_FEAS = "mpce_feas_s0of1.csv"
COL_OLD, COL_NEW = "#eb6834", "#1baf7a"       # validated pair (dataviz validator: CVD dE 9.2, normal dE 27.6)
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#e4e3df"


# ------------------------------------------------------------------ problem (mirrors mpce_experiments.run_grid)
class Ctx:
    """Problem context read by the instrumented classes (set per run, one run at a time per process)."""
    def __init__(self, ds, rad, n):
        self.ds, self.rad, self.n = ds, rad, n
        self.ideal = farm_objective(np.zeros((n, 2)), ds)[1]
        self.calls = []                          # (F, feasible) of every objective call (T2 only)
        self.record_calls = False
        self.smin, self.box = SMIN, rad

    def bviol(self, x):
        return bool(np.sqrt((x.reshape(-1, 2) ** 2).sum(1)).max() > self.rad + 1e-6)

    def sviol(self, x):
        return bool(self.n > 1 and min_spacing(x.reshape(-1, 2)) < SMIN - 1e-6)

    def feasible(self, x):
        return not (self.bviol(x) or self.sviol(x))

    def pct(self, F):
        return 100.0 * F / self.ideal


class HRCtx(Ctx):
    """Horns Rev 1 block (mirrors mpce_experiments.run_hr): parallelogram site, spacing 4D."""
    def __init__(self, n):
        self.ds, self.rad, self.n = "HR", 0, n
        _, self.ideal, self.poly = hr.make_objective(n)
        self.calls = []; self.record_calls = False
        self.smin, self.box = 4 * hr.D, float(np.abs(self.poly).max())

    def bviol(self, x):
        return bool(hr.outside_distance(x.reshape(-1, 2), self.poly).max() > 1e-6)

    def sviol(self, x):
        return bool(min_spacing(x.reshape(-1, 2)) < self.smin - 1e-6)


CTX = None


def setup_hr(n, init):
    """f, wake, feasible, box half-width, smin, boundary constraints exactly as in mpce_experiments.run_hr."""
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
    return f, wake, feasible, half, smin, bcons


def setup(ds, rad, n, init):
    """f, wake, feasible exactly as in mpce_experiments.run_grid."""
    smin = 8 * R
    init_hook.GEN = make_generator(smin, circle_r=rad) if init == "feasible" else None
    f = make_objective(ds, rad)

    def wake(x):
        o, i = farm_objective(x.reshape(-1, 2), ds)
        return i - o

    def feasible(x):
        xy = x.reshape(-1, 2)
        return (n == 1 or min_spacing(xy) >= smin - 1e-6) and np.sqrt((xy ** 2).sum(1)).max() <= rad + 1e-6
    return f, wake, feasible


def rms_dist(a, b, n):
    return float(np.sqrt(((a - b) ** 2).sum() / n))


# ------------------------------------------------------------------ instrumented PSO (copy of AO.PSO.optimize)
class PSOInstr(AO.PSO):
    def optimize(self, obj_fun, dim, lb, ub):
        C = CTX; n = dim // 2; LOG = []; self.log = LOG; first = {}; PSOInstr.LAST = self   # [instr]
        X = init_pop(self.pop_size, dim, lb, ub)
        V = np.zeros((self.pop_size, dim))
        pbest = X.copy()
        pbest_score = np.array([obj_fun(x) for x in X])
        gbest_idx = np.argmin(pbest_score)
        gbest = pbest[gbest_idx].copy(); gbest_score = pbest_score[gbest_idx]
        curve = [gbest_score]
        fe0 = [C.feasible(x) for x in X]                                          # [instr]
        LOG.append(self._rec(0, X, V, gbest, gbest_score, n, dim, 0, fe0, [C.bviol(x) for x in X],
                             [C.sviol(x) for x in X], 0, 0, np.zeros(self.pop_size), 0, 0))
        self.init_X = X.copy(); self.init_gbest_idx = int(gbest_idx)             # [instr]
        for it in range(self.max_iter):
            nclip = 0; fe = []; bv = []; sv = []; npb = 0; ngb = 0; disp = []; nfb = 0; nfp = 0; msp = []  # [instr]
            for i in range(self.pop_size):
                r1 = np.random.rand(dim); r2 = np.random.rand(dim)
                V[i] = (self.w * V[i] + self.c1 * r1 * (pbest[i] - X[i])
                        + self.c2 * r2 * (gbest - X[i]))
                xprev = X[i].copy()                                               # [instr]
                X[i] = X[i] + V[i]
                nclip += int(((X[i] < lb) | (X[i] > ub)).sum())                   # [instr] clipped coordinates
                X[i] = np.clip(X[i], lb, ub)
                score = obj_fun(X[i])
                xi = X[i].copy(); fz = C.feasible(xi)                             # [instr]
                fe.append(fz); bv.append(C.bviol(xi)); sv.append(C.sviol(xi)); disp.append(rms_dist(xi, xprev, n))
                msp.append(min_spacing(xi.reshape(-1, 2)) if n > 1 else np.inf)
                if fz and score < gbest_score:                                    # [instr]
                    nfb += 1
                if fz and score < pbest_score[i]:                                 # [instr]
                    nfp += 1
                if it == 0:                                                       # [instr] first update
                    first.setdefault("v_zero", []).append(bool(np.all(V[i] == 0)))
                    first.setdefault("moved", []).append(bool(np.any(xi != xprev)))
                    first.setdefault("is_gbest_particle", []).append(bool(np.all(xprev == gbest)))
                if score < pbest_score[i]:
                    pbest_score[i] = score; pbest[i] = X[i].copy()
                    npb += 1                                                      # [instr]
                if score < gbest_score:
                    gbest_score = score; gbest = X[i].copy()
                    ngb += 1                                                      # [instr]
            curve.append(gbest_score)
            LOG.append(self._rec(it + 1, X, V, gbest, gbest_score, n, dim, nclip, fe, bv, sv, npb, ngb,  # [instr]
                                 np.array(disp), nfb, nfp, msp))
            if it == 0:                                                           # [instr]
                first["min_spacing"] = msp; first["feasible"] = fe
        self.first = first                                                        # [instr]
        return gbest, gbest_score, np.array(curve)

    @staticmethod
    def _rec(t, X, V, gbest, gscore, n, dim, nclip, fe, bv, sv, npb, ngb, disp, nfb, nfp, msp=None):
        C = CTX
        spread = float(np.mean(np.sqrt(((X - gbest) ** 2).sum(1) / n)))
        gf = C.feasible(gbest)
        return dict(t=t, spread=spread, vel=float(np.abs(V).mean()), clip=nclip / (len(X) * dim),
                    feas=float(np.mean(fe)), bviol=float(np.mean(bv)), sviol=float(np.mean(sv)),
                    npbest=npb, ngbest=ngb, disp=float(np.mean(disp)), nzero_v=int((np.abs(V).sum(1) == 0).sum()),
                    gfeas=bool(gf), gwl=C.pct(gscore) if gf else np.nan, gF=float(gscore),
                    b_only=int(np.sum(np.array(bv) & ~np.array(sv))), s_only=int(np.sum(~np.array(bv) & np.array(sv))),
                    both=int(np.sum(np.array(bv) & np.array(sv))), nfeas=int(np.sum(fe)), ncand=len(fe),
                    nfeas_better_g=nfb, nfeas_better_p=nfp,
                    msp_med=float(np.median(msp)) if msp else np.nan)


# ------------------------------------------------------------------ instrumented DE (copy of AO.DE.optimize)
class DEInstr(AO.DE):
    def optimize(self, obj_fun, dim, lb, ub):
        C = CTX; n = dim // 2; LOG = []; self.log = LOG; first = {}; DEInstr.LAST = self    # [instr]
        pop = init_pop(self.pop_size, dim, lb, ub)
        fitness = np.array([obj_fun(ind) for ind in pop])
        best_idx = np.argmin(fitness)
        best_pos = pop[best_idx].copy(); best_score = fitness[best_idx]
        curve = [best_score]
        self.init_X = pop.copy()                                                  # [instr]
        for it in range(self.max_iter):
            fe = []; acc = 0; nb = 0; msp = []; bv = []; sv = []                  # [instr]
            for i in range(self.pop_size):
                idxs = list(range(self.pop_size)); idxs.remove(i)
                r1, r2, r3 = np.random.choice(idxs, 3, replace=False)
                mutant = pop[r1] + self.F * (pop[r2] - pop[r3])
                mutant = np.clip(mutant, lb, ub)
                trial = pop[i].copy()
                j_rand = np.random.randint(dim)
                for j in range(dim):
                    if (np.random.rand() < self.CR) or (j == j_rand):
                        trial[j] = mutant[j]
                trial = np.clip(trial, lb, ub)
                trial_fit = obj_fun(trial)
                fe.append(C.feasible(trial)); msp.append(min_spacing(trial.reshape(-1, 2)))  # [instr]
                bv.append(C.bviol(trial)); sv.append(C.sviol(trial))                       # [instr]
                if trial_fit < fitness[i]:
                    acc += 1                                                      # [instr]
                    pop[i] = trial; fitness[i] = trial_fit
                    if trial_fit < best_score:
                        nb += 1                                                   # [instr]
                        best_score = trial_fit; best_pos = trial.copy()
            curve.append(best_score)
            bv_, sv_ = np.array(bv), np.array(sv)                                 # [instr]
            LOG.append(dict(t=it + 1, feas=float(np.mean(fe)), accepted=acc, nbest=nb,                # [instr]
                            msp_med=float(np.median(msp)), nfeas=int(np.sum(fe)), ncand=len(fe),
                            b_only=int(np.sum(bv_ & ~sv_)), s_only=int(np.sum(~bv_ & sv_)), both=int(np.sum(bv_ & sv_)),
                            bviol=float(np.mean(bv_)), sviol=float(np.mean(sv_))))
            if it == 0:                                                           # [instr]
                first["min_spacing"] = msp; first["feasible"] = fe
        self.first = first                                                        # [instr]
        return best_pos, best_score, np.array(curve)


# ------------------------------------------------------------------ instrumented VNS (copy of OV.BVNS.search)
class BVNSInstr(OV.BVNS):
    LAST = None

    def search(self, f0, x, fx, dim, lb, ub):
        C = CTX                                                                   # [instr]
        L = dict(k_max=self.k_max, cycles=[], stage_evals=dict(descent=0, shake=0, ls=0), descent_moves=0,  # [instr]
                 descent_done=False, stop_stage=None)
        BVNSInstr.LAST = L
        stage = ["descent"]                                                       # [instr]
        # [instr] start handed to VNS vs. the best feasible layout among the Phase-1 calls
        p1 = C.calls if C.record_calls else []
        feasF = [F for F, fz in p1 if fz]
        L.update(p1_calls=len(p1), p1_any_feasible=bool(feasF), p1_best_feasible_F=min(feasF) if feasF else None,
                 x0_feasible=C.feasible(x), x0_F=float(fx), x0_wl=C.pct(fx) if C.feasible(x) else None)

        def f(z):                                                                 # [instr] counting wrapper
            v = f0(z)
            L["stage_evals"][stage[0]] += 1
            return v
        best = [x.copy(), fx]                               # best point evaluated so far

        def shake(x, k):
            lo, hi = self.rho[k - 1], self.rho[k]
            while True:                       # uniform in the l_inf shell (rejection)
                d = np.random.uniform(-hi, hi, dim)
                if np.abs(d).max() > lo:
                    return np.clip(x + d, lb, ub)

        def best_improvement(y, fy):
            h = self.h0
            while h >= self.hmin:
                cand_x, cand_f = None, fy
                for i in range(dim):
                    for s in (h, -h):
                        z = y.copy(); z[i] = min(max(z[i] + s, lb), ub)
                        fz = f(z)
                        if fz < best[1]:
                            best[0], best[1] = z.copy(), fz
                        if fz < cand_f:
                            cand_x, cand_f = z, fz
                if cand_x is None:
                    h *= 0.5                  # no improving move: refine the step
                else:
                    y, fy = cand_x, cand_f    # best-improvement move
                    if stage[0] == "descent":                                     # [instr]
                        L["descent_moves"] += 1
            return y, fy

        try:
            x, fx = best_improvement(x, fx)   # descend from the initial solution
            L.update(descent_done=True, x1_F=float(fx), x1_feasible=C.feasible(x))                    # [instr]
            k = 1
            while True:
                stage[0] = "shake"                                                # [instr]
                y = shake(x, k)
                fy = f(y)
                cyc = dict(k=k, F_before=float(fx), F_shake=float(fy), shake_feasible=C.feasible(y),  # [instr]
                           inc_feasible_before=C.feasible(x), complete=False)
                L["cycles"].append(cyc)                                           # [instr]
                if fy < best[1]:
                    best[0], best[1] = y.copy(), fy
                stage[0] = "ls"                                                   # [instr]
                e0 = L["stage_evals"]["ls"]                                       # [instr]
                y, fy = best_improvement(y, fy)
                cyc.update(F_ls=float(fy), ls_evals=L["stage_evals"]["ls"] - e0, complete=True,        # [instr]
                           accepted=bool(fy < fx), ls_feasible=C.feasible(y))
                if fy < fx:
                    x, fx, k = y, fy, 1
                else:
                    k = k + 1 if k < self.k_max else 1
        except BudgetExhausted:
            L["stop_stage"] = stage[0]                                            # [instr]
        # return the best point evaluated (equals x unless the budget ran out mid-descent)
        L.update(inc_F=float(fx), inc_feasible=C.feasible(x), best_F=float(best[1]),                   # [instr]
                 best_feasible=C.feasible(best[0]))
        if "x1_F" not in L:                                                       # [instr] descent cut off
            L.update(x1_F=float(fx), x1_feasible=C.feasible(x))
        return best[0], best[1], None


@contextmanager
def patched(pairs):
    old = [(m, a, getattr(m, a)) for m, a, _ in pairs]
    for m, a, v in pairs:
        setattr(m, a, v)
    try:
        yield
    finally:
        for m, a, v in old:
            setattr(m, a, v)


def run_one(alg, ds, rad, n, seed, init, instrument):
    """One run through mpce_experiments.run_method (instrumented classes swapped in when `instrument`)."""
    global CTX
    if ds == "HR":
        CTX = HRCtx(n)
        f, wake, feasible, box, smin, bcons = setup_hr(n, init)
    else:
        CTX = Ctx(ds, rad, n)
        f, wake, feasible = setup(ds, rad, n, init)
        box, smin, bcons = rad, SMIN, None
    pairs = []
    if instrument:
        pairs = [(E, "PSO", PSOInstr), (E, "DE", DEInstr), (HB, "BVNS", BVNSInstr), (RSm, "BVNS", BVNSInstr)]
        if alg in ("PSOBV", "SSABV", "RSVNS"):
            CTX.record_calls = True

            def fr(x, _f=f):
                v = _f(x)
                CTX.calls.append((float(v), CTX.feasible(x)))
                return v
            fuse = fr
        else:
            fuse = f
    else:
        fuse = f
    PSOInstr.LAST = DEInstr.LAST = BVNSInstr.LAST = None
    t0 = time.perf_counter()
    with patched(pairs):
        pos, tr = E.run_method(alg, seed, BUDGET, fuse, wake, feasible, 2 * n, -box, box, box, smin, bcons)
    sec = time.perf_counter() - t0
    xy = pos.reshape(-1, 2)
    if ds == "HR":
        ideal = CTX.ideal; obj = hr.aep_gwh(xy); cf, pf = "{:.2f}", "{:.4f}"
    else:
        obj, ideal = farm_objective(xy, ds); cf, pf = "{:.3f}", "{:.3f}"
    curve = ideal - np.array(tr.curve)
    out = dict(alg=alg, ds=ds, rad=rad, n=n, seed=seed, init=init, pos=pos.copy(), calls=tr.calls,
               curve=np.array(tr.curve), WakeLoss=ideal - obj, sec=sec,
               Coordinates=";".join(f"{cf.format(a)} {cf.format(b)}" for a, b in xy),
               Curve=";".join("nan" if not np.isfinite(c) else pf.format(c) for c in curve),
               Feasible=bool(feasible(pos)))
    if instrument:
        o = PSOInstr.LAST if alg in ("PSO", "PSOC") else DEInstr.LAST if alg == "DE" else None
        if o is not None:
            out["log"] = o.log; out["first"] = o.first
            if alg in ("PSOC", "PSO") and init == "feasible":
                out["init_X"] = o.init_X
                out["init_gbest_idx"] = getattr(o, "init_gbest_idx", None)
            if alg == "DE" and init == "feasible":
                out["init_X"] = o.init_X
        if alg in ("PSOBV", "SSABV", "RSVNS"):
            out["vns"] = BVNSInstr.LAST
    init_hook.GEN = None
    return out


def job(t):
    kind, alg, ds, rad, n, seed, init, verify = t
    r = run_one(alg, ds, rad, n, seed, init, True)
    r["kind"] = kind
    if verify:
        u = run_one(alg, ds, rad, n, seed, init, False)
        r["bitident"] = bool(np.array_equal(u["pos"], r["pos"]) and u["calls"] == r["calls"]
                             and np.array_equal(u["curve"], r["curve"], equal_nan=True)
                             and u["WakeLoss"] == r["WakeLoss"])
    return r


# ------------------------------------------------------------------ run everything
def all_jobs():
    J = []
    for alg in T1_ALGS:
        for c in T1_CASES:
            for s in T1_SEEDS:
                J.append(("T1", alg, *c, s, "random", s == 1))
    for alg in T2_ALGS:
        for c in T2_CASES:
            for s in T2_SEEDS:
                J.append(("T2", alg, *c, s, "random", s == 1))
    for alg in T3_ALGS:
        for c in T3_CASES:
            for s in T3_SEEDS:
                J.append(("T3", alg, *c, s, "feasible", s == 1))
    J.sort(key=lambda t: -t[4])            # larger N first (load balance)
    return J


def stored_feasinit():
    """Stored feasible-initialization runs of PSO and DE at 6,030 calls as used in the paper: the 6 largest
    cases from mpce_feas_s0of1.csv; Horns Rev 16 from the rerun with the corrected model (mpce_hrfix)."""
    d = pd.read_csv(os.path.join(HERE, STORED_FEAS)); d = d[d.Algorithm.isin(T3_ALGS) & (d.Dataset.astype(str) != "HR")]
    h = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(os.path.join(HERE, "mpce_hrfix_s*of24.csv")))])
    h = h[h.Algorithm.isin(T3_ALGS) & (h.Init == "feasible") & (h.Budget == 6030)]
    x = pd.concat([d, h], ignore_index=True); x["Dataset"] = x.Dataset.astype(str)
    return x


def stored_table():
    rows = []
    for alg, files in STORED.items():
        for fn in files:
            d = pd.read_csv(os.path.join(HERE, fn))
            d = d[d.Algorithm == alg].copy()
            if "Init" not in d:
                d["Init"] = "random"
            if "Budget" not in d:
                d["Budget"] = 6030
            rows.append(d)
    rows.append(stored_feasinit())
    S = pd.concat(rows, ignore_index=True)
    S["Dataset"] = S.Dataset.astype(str)
    S = S[S.Budget == 6030]
    return S.set_index(["Algorithm", "Dataset", "Radius", "Turbines", "Seed", "Init"]).sort_index()


def compare_stored(res, S):
    """(a) instrumented run == stored run of the paper."""
    out = []
    for r in res:
        key = (r["alg"], str(r["ds"]), r["rad"], r["n"], r["seed"], r["init"])
        try:
            s = S.loc[key]
        except KeyError:
            out.append((key, None)); continue
        if isinstance(s, pd.DataFrame):
            s = s.iloc[0]
        # the stored CSVs hold the WakeLoss float with 16-17 significant digits (last digit may be rounded), so the
        # float is compared to a relative tolerance of 1e-13; the coordinate and curve strings are compared exactly
        ok = (abs(float(s.WakeLoss) - float(r["WakeLoss"])) <= 1e-13 * abs(float(r["WakeLoss"])) and s.Coordinates == r["Coordinates"]
              and int(s.Calls) == int(r["calls"]) and s.Curve == r["Curve"])
        out.append((key, bool(ok)))
    return out


# ------------------------------------------------------------------ T1 summary
def t1_summary(res):
    R1 = [r for r in res if r["kind"] == "T1"]
    fields = ["spread", "vel", "clip", "feas", "bviol", "sviol", "disp", "gwl"]
    arr = {a: {f: np.array([[rec[f] for rec in r["log"]] for r in R1 if r["alg"] == a]) for f in fields} for a in T1_ALGS}
    gfe = {a: np.array([[rec["gfeas"] for rec in r["log"]] for r in R1 if r["alg"] == a]) for a in T1_ALGS}
    keys = {a: [(r["ds"], r["rad"], r["n"], r["seed"]) for r in R1 if r["alg"] == a] for a in T1_ALGS}
    assert sorted(keys["PSO"]) == sorted(keys["PSOC"])
    # align runs by key
    order = {a: np.argsort([str(k) for k in keys[a]]) for a in T1_ALGS}
    for a in T1_ALGS:
        for f in fields:
            arr[a][f] = arr[a][f][order[a]]
        gfe[a] = gfe[a][order[a]]
    ks = [keys["PSO"][i] for i in order["PSO"]]
    T = arr["PSO"]["spread"].shape[1] - 1
    late = slice(T // 2 + 1, T + 1)
    allit = slice(1, T + 1)
    S = dict(cases=[list(c) for c in T1_CASES], seeds=len(T1_SEEDS), iterations=int(T), runs_per_setting=len(ks))
    for a in T1_ALGS:
        A = arr[a]
        S[a] = dict(
            spread_final_median_m=float(np.median(A["spread"][:, -1])),
            spread_final_over_r_median=float(np.median(A["spread"][:, -1] / np.array([k[1] for k in ks]))),
            vel_final_median_m=float(np.median(A["vel"][:, -1])),
            vel_max_iter_median_m=float(np.median(A["vel"].max(1))),
            clip_pct_all=100 * float(A["clip"][:, allit].mean()), clip_pct_late=100 * float(A["clip"][:, late].mean()),
            clip_pct_first10=100 * float(A["clip"][:, 1:11].mean()),
            feas_particles_pct_all=100 * float(A["feas"][:, allit].mean()),
            feas_particles_pct_late=100 * float(A["feas"][:, late].mean()),
            bviol_particles_pct_late=100 * float(A["bviol"][:, late].mean()),
            sviol_particles_pct_late=100 * float(A["sviol"][:, late].mean()),
            bviol_particles_pct_all=100 * float(A["bviol"][:, allit].mean()),
            sviol_particles_pct_all=100 * float(A["sviol"][:, allit].mean()),
            disp_late_median_m=float(np.median(A["disp"][:, late].mean(1))),
            gbest_feasible_final=int(gfe[a][:, -1].sum()),
            gbest_first_feasible_iter_median=float(np.nanmedian([np.argmax(g) if g.any() else np.nan for g in gfe[a]])),
        )
    for a in T1_ALGS:       # infeasible particle evaluations (iterations 1-200) by violated constraint
        recs = [rec for r in R1 if r["alg"] == a for rec in r["log"][1:]]
        nc = sum(rec["ncand"] for rec in recs); ninf = nc - sum(rec["nfeas"] for rec in recs)
        S[a].update(evals=int(nc), infeasible_evals_pct=100 * ninf / nc,
                    infeasible_boundary_only_pct=100 * sum(rec["b_only"] for rec in recs) / ninf,
                    infeasible_spacing_only_pct=100 * sum(rec["s_only"] for rec in recs) / ninf,
                    infeasible_both_pct=100 * sum(rec["both"] for rec in recs) / ninf)
    ratio_spread = arr["PSO"]["spread"][:, -1] / arr["PSOC"]["spread"][:, -1]
    ratio_vel = arr["PSO"]["vel"][:, -1] / arr["PSOC"]["vel"][:, -1]
    S["paired"] = dict(spread_ratio_median=float(np.median(ratio_spread)), spread_ratio_min=float(ratio_spread.min()),
                       vel_ratio_median=float(np.median(ratio_vel)), vel_ratio_min=float(ratio_vel.min()),
                       n_old_clip_higher_late=int((arr["PSO"]["clip"][:, late].mean(1) > arr["PSOC"]["clip"][:, late].mean(1)).sum()),
                       n_old_feas_lower_late=int((arr["PSO"]["feas"][:, late].mean(1) < arr["PSOC"]["feas"][:, late].mean(1)).sum()))
    # per case
    pc = []
    for c in T1_CASES:
        idx = [i for i, k in enumerate(ks) if k[:3] == c]
        row = dict(case=list(c))
        for a in T1_ALGS:
            A = arr[a]
            row[a] = dict(spread=float(np.median(A["spread"][idx, -1])), vel=float(np.median(A["vel"][idx, -1])),
                          clip=100 * float(A["clip"][idx][:, allit].mean()),
                          feas=100 * float(A["feas"][idx][:, late].mean()),
                          bviol=100 * float(A["bviol"][idx][:, late].mean()),
                          sviol=100 * float(A["sviol"][idx][:, late].mean()),
                          gfeas=int(gfe[a][idx, -1].sum()),
                          gwl=float(np.nanmedian(A["gwl"][idx, -1])) if np.isfinite(A["gwl"][idx, -1]).any() else None)
        pc.append(row)
    S["per_case"] = pc
    curves = {a: {f: dict(med=np.median(arr[a][f], 0).tolist(), q1=np.percentile(arr[a][f], 25, 0).tolist(),
                          q3=np.percentile(arr[a][f], 75, 0).tolist()) for f in ("spread", "vel", "clip", "feas")}
              for a in T1_ALGS}
    return S, curves


def stored_violation_summary():
    """Constraint violations of the final layouts of all stored runs of both PSO settings (68 cases x 30)."""
    out = {}
    for a, fn in (("PSO", "fresh_grid.csv"), ("PSOC", "mpce_psoc_s0of1.csv")):
        d = pd.read_csv(os.path.join(HERE, fn)); d = d[d.Algorithm == a]
        nb = ns = both = inf = atbox = atbox_inf = corner_out = 0
        for _, r in d.iterrows():
            xy = np.array([[float(v) for v in p.split()] for p in r.Coordinates.split(";")])
            rad = float(r.Radius); n = len(xy)
            b = np.sqrt((xy ** 2).sum(1)).max() > rad + 1e-3          # coordinates stored to 1 mm
            s = n > 1 and min_spacing(xy) < SMIN - 1e-3
            onbox = bool((np.abs(xy) >= rad - 5e-4).any())
            atbox += onbox
            if b or s:
                inf += 1; nb += b and not s; ns += s and not b; both += b and s; atbox_inf += onbox
            # turbines outside the circle that sit on the box bound (clipped coordinate)
            out_t = np.sqrt((xy ** 2).sum(1)) > rad + 1e-3
            corner_out += int((out_t & (np.abs(xy) >= rad - 5e-4).any(1)).sum())
        out[a] = dict(runs=len(d), infeasible=inf, boundary_only=nb, spacing_only=ns, both=both,
                      any_boundary=nb + both, final_on_box=atbox, infeasible_on_box=atbox_inf,
                      outside_turbines_on_box=corner_out, feasible_flag=int(d.Feasible.sum()))
    return out


# ------------------------------------------------------------------ T2 summary
def t2_summary(res):
    R2 = [r for r in res if r["kind"] == "T2"]
    S = dict(cases=[list(c) for c in T2_CASES], seeds=len(T2_SEEDS))
    for a in T2_ALGS:
        rs = [r for r in R2 if r["alg"] == a]
        per = []
        for r in rs:
            L = r["vns"]; C = Ctx(r["ds"], r["rad"], r["n"])
            cyc = L["cycles"]; comp = [c for c in cyc if c["complete"]]; accd = [c for c in comp if c["accepted"]]
            ev = L["stage_evals"]; p2 = sum(ev.values())
            d = dict(n=r["n"], ds=r["ds"], rad=r["rad"], seed=r["seed"], p1_calls=L["p1_calls"],
                     phase2_evals=p2, descent_evals=ev["descent"], shake_evals=ev["shake"], ls_evals=ev["ls"],
                     descent_done=L["descent_done"], descent_moves=L["descent_moves"], shakes=len(cyc),
                     cycles_complete=len(comp), accepted=len(accd),
                     acc_k=[sum(c["k"] == k for c in accd) for k in range(1, 6)],
                     try_k=[sum(c["k"] == k for c in comp) for k in range(1, 6)],
                     max_k=max([c["k"] for c in cyc], default=0),
                     stop_stage=L["stop_stage"], x0_feasible=L["x0_feasible"], x1_feasible=L["x1_feasible"],
                     inc_feasible=L["inc_feasible"], best_feasible=L["best_feasible"],
                     p1_any_feasible=L["p1_any_feasible"],
                     start_is_best_feasible=(L["p1_best_feasible_F"] is not None and L["x0_feasible"]
                                             and L["x0_F"] <= L["p1_best_feasible_F"]),
                     shake_better=sum(c["F_shake"] < c["F_before"] for c in accd),
                     shake_feasible_pct=100 * np.mean([c["shake_feasible"] for c in cyc]) if cyc else np.nan)
            # feasibility changes: along the incumbent chain (switch -> after descent -> after accepted cycles)
            chain = [L["x0_feasible"], L["x1_feasible"]] + [c["ls_feasible"] for c in accd]
            d["made_feasible_stage"] = (None if L["x0_feasible"] else
                                        ("descent" if L["x1_feasible"] else ("cycle" if L["inc_feasible"] else "never")))
            d["lost_feasibility"] = any(chain[i] and not chain[i + 1] for i in range(len(chain) - 1))
            if L["x0_feasible"] and L["x1_feasible"]:
                # a descent cut off by the budget has not returned; its progress is in best_F (not in x1_F / inc_F)
                x1 = L["x1_F"] if L["descent_done"] else L["best_F"]
                d["imp_total_pp"] = C.pct(L["x0_F"] - L["best_F"])
                d["imp_descent_pp"] = C.pct(L["x0_F"] - x1)
                d["imp_cycles_pp"] = C.pct(sum(c["F_before"] - c["F_ls"] for c in accd))
                d["imp_cutoff_pp"] = C.pct(L["inc_F"] - L["best_F"]) if L["descent_done"] else 0.0
            per.append(d)
        P = pd.DataFrame(per)
        fe = P[P.imp_total_pp.notna()] if "imp_total_pp" in P else P.iloc[:0]
        tot = fe.imp_total_pp.sum()
        acc_k = np.sum(np.vstack(P.acc_k.values), 0); try_k = np.sum(np.vstack(P.try_k.values), 0)
        S[a] = dict(
            runs=len(P), switch_feasible=int(P.x0_feasible.sum()),
            p1_calls=int(P.p1_calls.median()), phase2_evals_median=float(P.phase2_evals.median()),
            descent_evals_median=float(P.descent_evals.median()), descent_evals_max=int(P.descent_evals.max()),
            descent_not_finished=int((~P.descent_done).sum()),
            descent_evals_pct_of_phase2=100 * float(P.descent_evals.sum() / P.phase2_evals.sum()),
            ls_evals_pct_of_phase2=100 * float((P.descent_evals.sum() + P.ls_evals.sum()) / P.phase2_evals.sum()),
            shake_evals_pct_of_phase2=100 * float(P.shake_evals.sum() / P.phase2_evals.sum()),
            shakes_median=float(P.shakes.median()), shakes_min=int(P.shakes.min()), shakes_max=int(P.shakes.max()),
            shakes_median_by_n={int(k): float(v) for k, v in P.groupby("n").shakes.median().items()},
            ls_evals_per_cycle_median=float(P.ls_evals.sum() / max(1, P.cycles_complete.sum())),
            accepted_mean=float(P.accepted.mean()), accepted_median=float(P.accepted.median()),
            accepted_total=int(P.accepted.sum()), cycles_total=int(P.cycles_complete.sum()),
            accept_rate_pct=100 * float(P.accepted.sum() / max(1, P.cycles_complete.sum())),
            accepted_by_k=acc_k.tolist(), tried_by_k=try_k.tolist(),
            accepted_k1_pct=100 * float(acc_k[0] / max(1, acc_k.sum())),
            max_k_median=float(P.max_k.median()), runs_reaching_kmax=int((P.max_k == 5).sum()),
            runs_no_accepted_cycle=int((P.accepted == 0).sum()),
            shake_better_pct_of_accepted=100 * float(P.shake_better.sum() / max(1, P.accepted.sum())),
            stop_in_ls_pct=100 * float(P.stop_stage.isin(["ls", "descent"]).mean()),
            imp_runs=int(len(fe)), imp_total_mean_pp=float(fe.imp_total_pp.mean()) if len(fe) else None,
            imp_descent_share_pct=100 * float(fe.imp_descent_pp.sum() / tot) if tot > 0 else None,
            imp_cycles_share_pct=100 * float(fe.imp_cycles_pp.sum() / tot) if tot > 0 else None,
            imp_cutoff_share_pct=100 * float(fe.imp_cutoff_pp.sum() / tot) if tot > 0 else None,
            imp_runs_zero=int((fe.imp_total_pp <= 0).sum()) if len(fe) else None,
            runs_descent_whole_phase2=int(((~P.descent_done) & (P.shakes == 0)).sum()),
            descent_not_finished_by_n={int(k): int((~g.descent_done).sum()) for k, g in P.groupby("n")},
            made_feasible_descent=int((P.made_feasible_stage == "descent").sum()),
            made_feasible_cycle=int((P.made_feasible_stage == "cycle").sum()),
            never_feasible=int((P.made_feasible_stage == "never").sum()),
            lost_feasibility=int(P.lost_feasibility.sum()),
            final_feasible=int(P.best_feasible.sum()),
            p1_any_feasible=int(P.p1_any_feasible.sum()),
            start_not_best_feasible=int((P.p1_any_feasible & ~P.start_is_best_feasible.astype(bool)).sum()),
        )
        S[a]["per_n"] = {int(k): dict(shakes=float(g.shakes.median()), descent=float(g.descent_evals.median()),
                                      accepted=float(g.accepted.median()),
                                      descent_share=(100 * float(g.imp_descent_pp.sum() / g.imp_total_pp.sum())
                                                     if "imp_total_pp" in g and g.imp_total_pp.sum() > 0 else None))
                         for k, g in P.groupby("n")}
    return S


# ------------------------------------------------------------------ T3 summary
MIX_MODES = ("raw", "matched", "matched_convex", "matched_scalar")


def mixture_stats(P, G, ctx, c2, rng, mode, draws=200):
    """Counterfactual first-update candidates between a particle's feasible start P (= its pbest; V = 0) and the
    global best G (separate RNG; nothing here touches the optimizer). Modes:
      raw             P + c2 U(0,1)^dim (G - P): the actual first PSO move (independent weight per coordinate,
                      up to c2 = 1.496, i.e. overshooting G), original turbine labels, clipped to the box
      matched         the same after relabelling the turbines of P to the nearest turbines of G (assignment)
      matched_convex  matched, independent weights in [0, 1] (no overshoot)
      matched_scalar  matched, one weight in [0, 1] for all coordinates (a convex combination of the layouts)
    Returns (share feasible, share with a turbine outside the site, median minimum spacing (m))."""
    n = ctx.n
    p = P.reshape(n, 2).copy(); g = G.reshape(n, 2)
    if mode != "raw":
        cost = ((p[:, None] - g[None]) ** 2).sum(-1)
        ri, ci = linear_sum_assignment(cost)
        q = np.empty_like(p); q[ci] = p[ri]; p = q
    ok = bo = 0; sp = []
    for _ in range(draws):
        if mode == "matched_scalar":
            w = rng.random()
        else:
            w = rng.random((n, 2)) * (c2 if mode in ("raw", "matched") else 1.0)
        x = np.clip(p + w * (g - p), -ctx.box, ctx.box).ravel()
        b = ctx.bviol(x); sv = ctx.sviol(x)
        ok += not (b or sv); bo += b; sp.append(min_spacing(x.reshape(-1, 2)))
    return ok / draws, bo / draws, float(np.median(sp))


def matched_rms(P, G, n):
    p = P.reshape(n, 2); g = G.reshape(n, 2)
    cost = ((p[:, None] - g[None]) ** 2).sum(-1)
    ri, ci = linear_sum_assignment(cost)
    return float(np.sqrt(cost[ri, ci].mean())), float(np.sqrt(((p - g) ** 2).sum(1).mean()))


def t3_summary(res):
    R3 = [r for r in res if r["kind"] == "T3"]
    S = dict(cases=[list(c) for c in T3_CASES], seeds=len(T3_SEEDS))
    rng = np.random.default_rng(20260928)
    for a in T3_ALGS:
        rs = [r for r in R3 if r["alg"] == a]
        first_feas = []; first_msp = []; improved = 0; nbest = 0; nfb = 0; nfp = 0; npb = 0; acc = 0
        ncand = nfeas = nb_only = ns_only = nboth = 0
        frozen_runs = 0; gb_particle_still = 0; mix = {m: [] for m in MIX_MODES}; mrms = []; rawrms = []
        it_one_feas = 0; it_total = 0; disp_late = []; vel_late = []; per_case = {}
        for r in rs:
            L = r["log"]; n = r["n"]
            ctx = HRCtx(n) if r["ds"] == "HR" else Ctx(r["ds"], r["rad"], n)
            first_feas += list(r["first"]["feasible"]); first_msp += list(r["first"]["min_spacing"])
            curve = r["curve"]; imp = int(np.nanmin(curve) < curve[0]); improved += imp
            it = L[1:] if a == "PSOC" else L
            c_nc = sum(rec["ncand"] for rec in it); c_nf = sum(rec["nfeas"] for rec in it)
            ncand += c_nc; nfeas += c_nf
            nb_only += sum(rec["b_only"] for rec in it); ns_only += sum(rec["s_only"] for rec in it)
            nboth += sum(rec["both"] for rec in it)
            key = f"{r['ds']}-{r['rad']}-{n}"
            pc = per_case.setdefault(key, dict(runs=0, improved=0, cand=0, feas=0))
            pc["runs"] += 1; pc["improved"] += imp; pc["cand"] += c_nc; pc["feas"] += c_nf
            if a == "PSOC":
                nbest += sum(rec["ngbest"] for rec in it); npb += sum(rec["npbest"] for rec in it)
                nfb += sum(rec["nfeas_better_g"] for rec in it); nfp += sum(rec["nfeas_better_p"] for rec in it)
                frozen_runs += int(all(rec["nzero_v"] >= 1 for rec in it))
                # iterations in which exactly one of the 30 candidates is feasible (the G-particle re-evaluating G)
                it_one_feas += sum(rec["nfeas"] == 1 for rec in it); it_total += len(it)
                fi = r["first"]
                gb_particle_still += int(any(g and v and not m for g, v, m in zip(fi["is_gbest_particle"], fi["v_zero"], fi["moved"])))
                disp_late.append(np.mean([rec["disp"] for rec in L[-50:]])); vel_late.append(np.mean([rec["vel"] for rec in L[-50:]]))
                X0 = r["init_X"]; gi = r["init_gbest_idx"]; G = X0[gi]
                for i in range(len(X0)):
                    if i == gi:
                        continue
                    for m in MIX_MODES:
                        mix[m].append(mixture_stats(X0[i], G, ctx, 1.49618, rng, m))
                    m1, m0 = matched_rms(X0[i], G, n); mrms.append(m1); rawrms.append(m0)
            else:
                nbest += sum(rec["nbest"] for rec in it); acc += sum(rec["accepted"] for rec in it)
        ninf = ncand - nfeas
        d = dict(runs=len(rs), improved_runs=improved, gbest_updates=int(nbest), candidates=int(ncand),
                 feasible_candidates=int(nfeas), feasible_candidates_pct=100 * nfeas / ncand,
                 infeasible_boundary_only_pct=100 * nb_only / max(1, ninf), infeasible_spacing_only_pct=100 * ns_only / max(1, ninf),
                 infeasible_both_pct=100 * nboth / max(1, ninf),
                 first_update_feasible_pct=100 * float(np.mean(first_feas)),
                 first_update_min_spacing_median_m=float(np.median(first_msp)),
                 first_update_min_spacing_max_m=float(np.max(first_msp)), per_case=per_case)
        if a == "PSOC":
            d.update(feasible_better_than_gbest=int(nfb), feasible_better_than_pbest=int(nfp),
                     pbest_updates=int(npb), pbest_updates_per_particle=npb / (30 * len(rs)),
                     feasible_candidates_not_G_reevaluation=int(nfeas - it_total),
                     runs_gbest_particle_zero_velocity_throughout=frozen_runs,
                     runs_gbest_particle_not_moved_first_update=gb_particle_still,
                     iterations_exactly_one_feasible=int(it_one_feas), iterations=int(it_total),
                     mixture={m: dict(feasible_pct=100 * float(np.mean([v[0] for v in mix[m]])),
                                      outside_pct=100 * float(np.mean([v[1] for v in mix[m]])),
                                      min_spacing_median_m=float(np.median([v[2] for v in mix[m]]))) for m in MIX_MODES},
                     init_matched_rms_median_m=float(np.median(mrms)), init_matched_rms_min_m=float(np.min(mrms)),
                     init_unmatched_rms_median_m=float(np.median(rawrms)),
                     disp_last50_median_m=float(np.median(disp_late)), vel_last50_median_m=float(np.median(vel_late)))
        else:
            d.update(accepted_trials=int(acc), feasible_better_than_gbest=0 if nbest == 0 else None)
        S[a] = d
    # all stored feasible-start runs of the paper (6 largest cases, mpce_feas; Horns Rev 16, mpce_hrfix)
    x = stored_feasinit(); st = {}
    for a in T3_ALGS:
        y = x[x.Algorithm == a]
        first = y.Curve.str.split(";").str[0].astype(float).values
        last = y.Curve.str.split(";").str[-1].astype(float).values
        # the curve stores the objective (AEP for Horns Rev) of the best feasible layout so far; its first checkpoint
        # (call 30) is the best initial layout. Final objective == first checkpoint <=> no improvement.
        st[a] = dict(runs=len(y), runs_final_differs_from_initial_best=int((np.abs(y.Objective.values - first) > 1e-3).sum()),
                     curve_moves=int((first != last).sum()), all_feasible=bool(y.Feasible.all()))
    S["stored"] = st
    return S


# ------------------------------------------------------------------ extra: calls to first feasibility (stored)
def first_feasible_calls():
    """From the stored convergence curves (checkpoint every 30 calls at 6,030): calls until the first feasible
    layout, Phase 1 of the hybrids (identical to their stand-alone Phase 1 up to the switch), over 68 x 30 runs."""
    out = {}
    for a, files in (("PSOBV", STORED["PSOBV"]), ("SSABV", STORED["SSABV"]), ("RSVNS", STORED["RSVNS"])):
        d = pd.concat([pd.read_csv(os.path.join(HERE, f)) for f in files])
        d = d[(d.Algorithm == a) & (d.Dataset.astype(str) != "HR")]
        calls = []
        for cv in d.Curve:
            v = np.array([np.nan if s == "nan" else float(s) for s in cv.split(";")])
            k = np.argmax(np.isfinite(v)) if np.isfinite(v).any() else None
            calls.append(np.nan if k is None else 30 * (k + 1))
        calls = np.array(calls, float)
        out[a] = dict(runs=len(calls), median_calls=float(np.nanmedian(calls)),
                      pct_by_switch=100 * float(np.mean(calls <= 3030)), pct_by_30=100 * float(np.mean(calls <= 30)),
                      never=int(np.isnan(calls).sum()))
    return out


# ------------------------------------------------------------------ RS-VNS Phase 1 replay (geometry only)
def rs_replay(seeds=range(1, 31), n1=3015, disc_seed=20260928):
    """Replay the seeded random stream of RS-VNS Phase 1 (random initialization): np.random.seed(seed), the common
    initial population uniform(-r, r, (30, 2N)), then n1 - 30 samples uniform(-r, r, (1, 2N)) -- sequential draws
    from the legacy stream, so one uniform(-r, r, (n1, 2N)) call gives the same numbers. Counts the feasible
    samples (1e-6 m tolerance) and the samples with all turbines inside the circle, per case over 30 seeds; compares
    with (pi/4)^N. No objective call. Also a Monte Carlo estimate (separate RNG, same number of samples) of the
    feasibility rate of disc samples (the Phase-1 distribution of the planned control RSDVNS)."""
    out = {}
    rng = np.random.default_rng(disc_seed)
    for ds, rad, n in E.GRID:
        if ds != 1:                       # geometry does not depend on the wind data set: Data Set II identical
            continue
        tot_f = tot_in = 0; runs_any = 0; per_seed = []
        for sd in seeds:
            np.random.seed(sd)
            X = np.random.uniform(-rad, rad, (n1, 2 * n)).reshape(n1, n, 2)
            inside = np.sqrt((X ** 2).sum(-1)).max(-1) <= rad + 1e-6
            if n > 1:
                dd = np.sqrt(((X[:, :, None] - X[:, None]) ** 2).sum(-1))
                iu = np.triu_indices(n, 1)
                spaced = dd[:, iu[0], iu[1]].min(-1) >= SMIN - 1e-6
            else:
                spaced = np.ones(n1, bool)
            fz = inside & spaced
            tot_f += int(fz.sum()); tot_in += int(inside.sum()); runs_any += int(fz.any()); per_seed.append(int(fz.sum()))
        m = n1 * len(seeds)
        ang = rng.uniform(0, 2 * np.pi, (m, n)); rr = rad * np.sqrt(rng.uniform(0, 1, (m, n)))
        Y = np.stack([rr * np.cos(ang), rr * np.sin(ang)], -1)
        if n > 1:
            dd = np.sqrt(((Y[:, :, None] - Y[:, None]) ** 2).sum(-1)); iu = np.triu_indices(n, 1)
            disc_f = int((dd[:, iu[0], iu[1]].min(-1) >= SMIN - 1e-6).sum())
        else:
            disc_f = m
        out[f"{rad}-{n}"] = dict(rad=rad, n=n, samples=m, feasible=tot_f, inside=tot_in, runs_with_feasible=runs_any,
                                 p_feasible=tot_f / m, p_inside=tot_in / m, p_inside_theory=(np.pi / 4) ** n,
                                 disc_feasible=disc_f, p_feasible_disc=disc_f / m, feasible_per_seed=per_seed)
    return out


# ------------------------------------------------------------------ T4 CPU estimates
def cpu_estimates():
    sp = pd.read_csv(os.path.join(HERE, "mpce_psosplit_s0of1.csv"))
    rv = pd.read_csv(os.path.join(HERE, "mpce_rsvns_s0of1.csv"))
    s75 = float(sp[sp.Algorithm == "PSOBV75"].Seconds.mean())
    srs = float(rv.Seconds.mean())
    return dict(
        omega90=dict(runs=len(E.tasks("omega90")), sec_per_run=s75, cpu_hours=len(E.tasks("omega90")) * s75 / 3600,
                     basis="mean Seconds of PSOBV75 in mpce_psosplit (same 12 cases)",
                     command="python3 mpce_experiments.py omega90 <i> <n> --procs=4",
                     output="mpce_omega90_s<i>of<n>.csv", label="PSOBV90"),
        rsdisc=dict(runs=len(E.tasks("rsdisc")), sec_per_run=srs, cpu_hours=len(E.tasks("rsdisc")) * srs / 3600,
                    basis="mean Seconds of RSVNS in mpce_rsvns (same 68 cases)",
                    command="python3 mpce_experiments.py rsdisc <i> <n> --procs=4",
                    output="mpce_rsdisc_s<i>of<n>.csv", label="RSDVNS"),
    )


# ------------------------------------------------------------------ figure
def figure(curves, S1):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.family": "serif", "font.size": 8, "axes.edgecolor": MUTED, "axes.labelcolor": INK,
                         "xtick.color": MUTED, "ytick.color": MUTED, "axes.grid": True, "grid.color": GRID,
                         "grid.linewidth": 0.5, "axes.spines.top": False, "axes.spines.right": False,
                         "legend.frameon": False, "lines.linewidth": 1.4, "pdf.fonttype": 42})
    fig, axes = plt.subplots(1, 3, figsize=(7.1, 2.15))
    spec = [("spread", "swarm spread to $\\mathbf{G}$ (m)", True, 1.0),
            ("clip", "coordinates clipped (%)", False, 100.0),
            ("feas", "particles feasible (%)", False, 100.0)]
    for ax, (f, yl, logy, sc) in zip(axes, spec):
        for a, col, lsty in (("PSO", COL_OLD, "--"), ("PSOC", COL_NEW, "-")):
            c = curves[a][f]; x = np.arange(len(c["med"]))
            ax.fill_between(x, np.array(c["q1"]) * sc, np.array(c["q3"]) * sc, color=col, alpha=0.15, lw=0)
            ax.plot(x, np.array(c["med"]) * sc, color=col, ls=lsty, lw=1.5 if a == "PSOC" else 1.3,
                    label="old: $w=0.7$, $c_1=c_2=2$" if a == "PSO" else "constriction")
        if logy:
            ax.set_yscale("log")
        ax.set_xlabel("iteration"); ax.set_ylabel(yl); ax.set_xlim(0, len(x) - 1)
    h, l = axes[0].get_legend_handles_labels()
    fig.tight_layout(w_pad=1.2, rect=(0, 0, 1, 0.9))
    fig.legend(h, l, loc="upper center", ncol=2, fontsize=7.5, bbox_to_anchor=(0.5, 1.0))
    d = os.path.join(HERE, "..", "figures_mpce")
    fig.savefig(os.path.join(d, "diag_pso_dynamics.pdf"), bbox_inches="tight")
    fig.savefig(os.path.join(d, "diag_pso_dynamics.png"), dpi=160, bbox_inches="tight")
    plt.close(fig)


# ------------------------------------------------------------------ outputs
def fmt(v, d=1):
    if v is None or (isinstance(v, float) and not np.isfinite(v)):
        return "--"
    if isinstance(v, (int, np.integer)):
        s = f"{int(v):,}"
    else:
        s = f"{v:,.{d}f}"
    return s.replace(",", "{,}")


def macros(Sm):
    S1, V, S2, S3, FF, C4 = Sm["T1"], Sm["T1_stored"], Sm["T2"], Sm["T3"], Sm["first_feasible"], Sm["T4"]
    M = {}
    M["NDDynCases"] = len(S1["cases"]); M["NDDynSeeds"] = S1["seeds"]; M["NDDynRuns"] = S1["runs_per_setting"]
    for a, s in (("PSO", "Old"), ("PSOC", "New")):
        A = S1[a]
        M[f"NDSpread{s}"] = fmt(A["spread_final_median_m"], 1)
        M[f"NDVel{s}"] = fmt(A["vel_final_median_m"], 1)
        M[f"NDClip{s}"] = fmt(A["clip_pct_all"], 1)
        M[f"NDClipLate{s}"] = fmt(A["clip_pct_late"], 1)
        M[f"NDFeasPart{s}"] = fmt(A["feas_particles_pct_late"], 1)
        M[f"NDBViolPart{s}"] = fmt(A["bviol_particles_pct_late"], 1)
        M[f"NDSViolPart{s}"] = fmt(A["sviol_particles_pct_late"], 1)
        M[f"NDGbestFeas{s}"] = A["gbest_feasible_final"]
    M["NDSpreadRatio"] = fmt(S1["paired"]["spread_ratio_median"], 0 if S1["paired"]["spread_ratio_median"] >= 10 else 1)
    M["NDVelRatio"] = fmt(S1["paired"]["vel_ratio_median"], 0 if S1["paired"]["vel_ratio_median"] >= 10 else 1)
    for a, s in (("PSO", "Old"), ("PSOC", "New")):
        v = V[a]
        M[f"NDStoredInfeas{s}"] = v["infeasible"]; M[f"NDStoredRuns{s}"] = v["runs"]
        M[f"NDStoredBoundOnly{s}"] = v["boundary_only"]; M[f"NDStoredSpacingOnly{s}"] = v["spacing_only"]
        M[f"NDStoredBoth{s}"] = v["both"]; M[f"NDStoredOnBox{s}"] = v["final_on_box"]
        M[f"NDStoredInfeasOnBox{s}"] = v["infeasible_on_box"]
    M["NDVnsCases"] = len(S2["cases"]); M["NDVnsSeeds"] = S2["seeds"]
    for a in T2_ALGS:
        A = S2[a]; p = MAC[a]
        M[f"NDShakes{p}"] = fmt(A["shakes_median"], 0 if float(A["shakes_median"]).is_integer() else 1)
        M[f"NDShakesMin{p}"] = A["shakes_min"]; M[f"NDShakesMax{p}"] = A["shakes_max"]
        M[f"NDDescentEvals{p}"] = fmt(A["descent_evals_median"], 0)
        M[f"NDDescentEvalsPct{p}"] = fmt(A["descent_evals_pct_of_phase2"], 1)
        M[f"NDLSEvalsPct{p}"] = fmt(A["ls_evals_pct_of_phase2"], 2)
        M[f"NDShakeEvalsPct{p}"] = fmt(A["shake_evals_pct_of_phase2"], 2)
        M[f"NDAccepted{p}"] = fmt(A["accepted_mean"], 1)
        M[f"NDAcceptRate{p}"] = fmt(A["accept_rate_pct"], 1)
        M[f"NDAccKOne{p}"] = fmt(A["accepted_k1_pct"], 1)
        M[f"NDNoAccepted{p}"] = A["runs_no_accepted_cycle"]
        M[f"NDDescentShare{p}"] = fmt(A["imp_descent_share_pct"], 1)
        M[f"NDCycleShare{p}"] = fmt(A["imp_cycles_share_pct"], 1)
        M[f"NDShakeBetter{p}"] = fmt(A["shake_better_pct_of_accepted"], 1)
        M[f"NDSwitchInfeas{p}"] = A["runs"] - A["switch_feasible"]
        M[f"NDMadeFeas{p}"] = A["made_feasible_descent"] + A["made_feasible_cycle"]
        M[f"NDMadeFeasDescent{p}"] = A["made_feasible_descent"]
        M[f"NDLostFeas{p}"] = A["lost_feasibility"]
        M[f"NDStartNotBestFeas{p}"] = A["start_not_best_feasible"]
        M[f"NDVnsRuns{p}"] = A["runs"]
        nn = A["shakes_median_by_n"]
        M[f"NDShakesNSix{p}"] = fmt(nn[6], 0 if float(nn[6]).is_integer() else 1)
        M[f"NDShakesNFifteen{p}"] = fmt(nn[15], 0 if float(nn[15]).is_integer() else 1)
        dn = A["descent_not_finished_by_n"]
        M[f"NDDescentUnfinishedLarge{p}"] = dn[12] + dn[15]
        M[f"NDRunsLarge{p}"] = 2 * S2["seeds"] * 2           # N = 12 and N = 15, both data sets
    M["NDFsCases"] = len(S3["cases"]); M["NDFsSeeds"] = S3["seeds"]
    for a, p in (("PSOC", "PSO"), ("DE", "DE")):
        A = S3[a]
        M[f"NDFsRuns{p}"] = A["runs"]; M[f"NDFsImproved{p}"] = A["improved_runs"]
        M[f"NDFsCand{p}"] = fmt(A["candidates"])
        M[f"NDFsFeasCand{p}"] = fmt(A["feasible_candidates_pct"], 2)
        M[f"NDFsFirstFeas{p}"] = fmt(A["first_update_feasible_pct"], 1)
        M[f"NDFsFirstMinSp{p}"] = fmt(A["first_update_min_spacing_median_m"], 0)
        M[f"NDFsInfBoundOnly{p}"] = fmt(A["infeasible_boundary_only_pct"], 1)
        M[f"NDFsInfSpacingOnly{p}"] = fmt(A["infeasible_spacing_only_pct"], 1)
        M[f"NDFsInfBoth{p}"] = fmt(A["infeasible_both_pct"], 1)
        M[f"NDFsStoredRuns{p}"] = S3["stored"][a]["runs"]
        M[f"NDFsStoredMoved{p}"] = S3["stored"][a]["runs_final_differs_from_initial_best"]
    M["NDFsAcceptedDE"] = S3["DE"]["accepted_trials"]; M["NDFsFeasCandCountDE"] = S3["DE"]["feasible_candidates"]
    A = S3["PSOC"]
    for m, nm in zip(MIX_MODES, ("Raw", "Matched", "Convex", "Scalar")):
        M[f"NDFsMix{nm}"] = fmt(A["mixture"][m]["feasible_pct"], 1)
        M[f"NDFsMixOut{nm}"] = fmt(A["mixture"][m]["outside_pct"], 1)
        M[f"NDFsMixMinSp{nm}"] = fmt(A["mixture"][m]["min_spacing_median_m"], 0)
    M["NDFsOneFeasIter"] = fmt(A["iterations_exactly_one_feasible"]); M["NDFsIter"] = fmt(A["iterations"])
    M["NDFsNotG"] = A["feasible_candidates_not_G_reevaluation"]
    M["NDFsPbestUpd"] = A["pbest_updates"]
    M["NDFsInitRms"] = fmt(A["init_matched_rms_median_m"], 0)
    M["NDFsInitRmsMin"] = fmt(A["init_matched_rms_min_m"], 0)
    M["NDFsFrozen"] = A["runs_gbest_particle_zero_velocity_throughout"]
    M["NDFsFeasBetterG"] = A["feasible_better_than_gbest"]
    M["NDFsFeasBetterP"] = A["feasible_better_than_pbest"]
    M["NDFsDispLate"] = fmt(A["disp_last50_median_m"], 1)
    for a, s_ in (("PSO", "Old"), ("PSOC", "New")):
        M[f"NDInfBoundOnly{s_}"] = fmt(S1[a]["infeasible_boundary_only_pct"], 1)
        M[f"NDInfSpacingOnly{s_}"] = fmt(S1[a]["infeasible_spacing_only_pct"], 1)
        M[f"NDInfBoth{s_}"] = fmt(S1[a]["infeasible_both_pct"], 1)
        M[f"NDInfEvals{s_}"] = fmt(S1[a]["infeasible_evals_pct"], 1)
    RP = Sm["rs_replay"]
    big = RP["1000-15"]
    M["NDRsInsideNFifteen"] = fmt(100 * big["p_inside"], 2); M["NDRsInsideTheoryNFifteen"] = fmt(100 * big["p_inside_theory"], 2)
    M["NDRsFeasNFifteen"] = fmt(big["feasible"]); M["NDRsSamples"] = fmt(big["samples"])
    M["NDRsDiscFeasNFifteen"] = f"{100 * big['p_feasible_disc']:.1g}"          # percent
    M["NDRsGeomNoFeas"] = sum(v["feasible"] == 0 for v in RP.values()); M["NDRsGeom"] = len(RP)
    M["NDRsCasesNoFeasSample"] = 2 * sum(v["runs_with_feasible"] == 0 for v in RP.values())
    M["NDRsRunsNoFeasSample"] = fmt(2 * sum(30 - v["runs_with_feasible"] for v in RP.values()))
    M["NDRsReplayChecked"] = Sm["verification"]["rs_replay_checked"]; M["NDRsReplayAgree"] = Sm["verification"]["rs_replay_agree"]
    M["NDSmin"] = fmt(SMIN, 0)
    for a in ("PSOBV", "SSABV", "RSVNS"):
        M[f"NDFirstFeasCalls{MAC[a]}"] = fmt(FF[a]["median_calls"], 0)
        M[f"NDFirstFeasBySwitch{MAC[a]}"] = fmt(FF[a]["pct_by_switch"], 1)
        M[f"NDFirstFeasAtInit{MAC[a]}"] = fmt(FF[a]["pct_by_30"], 1)
    M["NDOmegaNinetyRuns"] = fmt(C4["omega90"]["runs"]); M["NDOmegaNinetyCPUh"] = fmt(C4["omega90"]["cpu_hours"], 2)
    M["NDRsDiscRuns"] = fmt(C4["rsdisc"]["runs"]); M["NDRsDiscCPUh"] = fmt(C4["rsdisc"]["cpu_hours"], 2)
    ver = Sm["verification"]
    M["NDVerRuns"] = fmt(ver["instrumented_runs"]); M["NDVerStoredMatch"] = fmt(ver["stored_match"])
    M["NDVerBitRuns"] = fmt(ver["bitident_runs"]); M["NDVerBitMatch"] = fmt(ver["bitident_match"])
    return M


def write_tex(M, path):
    with open(path, "w") as fh:
        fh.write("% generated by mpce_diagnostics.py from mpce_summary_diag.json -- do not edit by hand\n")
        for k, v in M.items():
            fh.write(f"\\newcommand{{\\{k}}}{{{v}}}\n")


def supp_tex(Sm, path):
    S1, V, S2, S3, FF = Sm["T1"], Sm["T1_stored"], Sm["T2"], Sm["T3"], Sm["first_feasible"]
    L = ["%% generated by mpce_diagnostics.py -- do not edit by hand",
         "%% Supplementary diagnostics (Phase 6, W3); requires booktabs, graphicx and mpce_numbers_diag.tex", ""]
    DSR = {1: "I", 2: "II"}
    # Figure
    L += ["\\begin{figure*}[!t]", "\\centering",
          "\\includegraphics[width=\\textwidth]{figures_mpce/diag_pso_dynamics.pdf}",
          "\\caption{PSO dynamics with the old setting ($w=0.7$, $c_1=c_2=2$) and the constriction setting "
          f"on {len(S1['cases'])} cases $\\times$ {S1['seeds']} seeds at 6,030 evaluations (200 iterations): median (line) "
          "and interquartile range (band) over the runs. Left: swarm spread, the mean over particles of the RMS distance "
          "of the turbines to the global best $\\mathbf G$ (log scale). Middle: share of the coordinates that leave the "
          "bounding square in an iteration and are clipped (velocities are kept). Right: share of particle evaluations "
          "whose layout is feasible. Per-case values: Table~\\ref{tab:D-psodyn}.}",
          "\\label{fig:D-psodyn}", "\\end{figure*}", ""]
    # Table PSO dynamics per case
    L += ["\\begin{table*}[!t]", "\\centering",
          "\\caption{PSO dynamics per case, old setting / constriction setting (10 seeds each, 6,030 evaluations). "
          "Spread: median over the runs of the swarm spread at the last iteration (m); $|V|$: median mean absolute "
          "velocity per coordinate at the last iteration (m); clipped: share of coordinates clipped per iteration "
          "(iterations 1--200); feasible, outside, spacing: share of particle evaluations in iterations 101--200 with a "
          "feasible layout, with a turbine outside the circle, with two turbines closer than $8R$; $\\mathbf G$ feas.: "
          "runs whose final global best is feasible; $L$: median wake loss of the feasible final global bests (\\%). "
          "Note below the table: violations of the infeasible evaluations and of the final layouts of all stored runs of both settings "
          "(68 cases $\\times$ 30 seeds).}",
          "\\label{tab:D-psodyn}", "\\scriptsize\\setlength{\\tabcolsep}{3pt}",
          "\\begin{tabular}{lcccccccc}", "\\toprule",
          "Case (DS, $r$, $N$) & Spread (m) & $|V|$ (m) & Clipped (\\%) & Feasible (\\%) & Outside (\\%) & Spacing (\\%) & $\\mathbf G$ feas. & $L$ (\\%) \\\\",
          "\\midrule"]
    for row in S1["per_case"]:
        ds, r, n = row["case"]
        o, c = row["PSO"], row["PSOC"]
        L.append(f"{DSR[ds]}, {r}, {n} & {fmt(o['spread'],1)} / {fmt(c['spread'],1)} & {fmt(o['vel'],1)} / {fmt(c['vel'],1)} & "
                 f"{fmt(o['clip'],1)} / {fmt(c['clip'],1)} & {fmt(o['feas'],1)} / {fmt(c['feas'],1)} & "
                 f"{fmt(o['bviol'],1)} / {fmt(c['bviol'],1)} & {fmt(o['sviol'],1)} / {fmt(c['sviol'],1)} & "
                 f"{o['gfeas']} / {c['gfeas']} & {fmt(o['gwl'],2)} / {fmt(c['gwl'],2)} \\\\")
    a, b = S1["PSO"], S1["PSOC"]
    L.append("\\midrule")
    L.append(f"All {len(S1['cases'])} cases & {fmt(a['spread_final_median_m'],1)} / {fmt(b['spread_final_median_m'],1)} & "
             f"{fmt(a['vel_final_median_m'],1)} / {fmt(b['vel_final_median_m'],1)} & {fmt(a['clip_pct_all'],1)} / {fmt(b['clip_pct_all'],1)} & "
             f"{fmt(a['feas_particles_pct_late'],1)} / {fmt(b['feas_particles_pct_late'],1)} & "
             f"{fmt(a['bviol_particles_pct_late'],1)} / {fmt(b['bviol_particles_pct_late'],1)} & "
             f"{fmt(a['sviol_particles_pct_late'],1)} / {fmt(b['sviol_particles_pct_late'],1)} & "
             f"{a['gbest_feasible_final']} / {b['gbest_feasible_final']} & \\\\")
    L += ["\\bottomrule", "\\end{tabular}", "\\\\[3pt]",
          "\\parbox{\\textwidth}{\\scriptsize Infeasible particle evaluations (iterations 1--200; old / constriction): "
          f"{fmt(a['infeasible_evals_pct'],1)} / {fmt(b['infeasible_evals_pct'],1)}\\% of all; of these, outside the circle only "
          f"{fmt(a['infeasible_boundary_only_pct'],1)} / {fmt(b['infeasible_boundary_only_pct'],1)}\\%, spacing only "
          f"{fmt(a['infeasible_spacing_only_pct'],1)} / {fmt(b['infeasible_spacing_only_pct'],1)}\\%, both "
          f"{fmt(a['infeasible_both_pct'],1)} / {fmt(b['infeasible_both_pct'],1)}\\%. "
          "Final layouts of all stored runs (old / constriction): "
          f"infeasible {V['PSO']['infeasible']} / {V['PSOC']['infeasible']} of {fmt(V['PSO']['runs'])}; of these, "
          f"outside the circle only {V['PSO']['boundary_only']} / {V['PSOC']['boundary_only']}, "
          f"spacing only {V['PSO']['spacing_only']} / {V['PSOC']['spacing_only']}, both {V['PSO']['both']} / {V['PSOC']['both']}; "
          "final layouts with a coordinate on the box bound: "
          f"{fmt(V['PSO']['final_on_box'])} / {fmt(V['PSOC']['final_on_box'])} (infeasible ones: {V['PSO']['infeasible_on_box']} / {V['PSOC']['infeasible_on_box']}).}}",
          "\\end{table*}", ""]
    # Table VNS internals
    heads = " & ".join(LAB[a] for a in T2_ALGS)
    L += ["\\begin{table}[!t]", "\\centering",
          f"\\caption{{Inside the VNS phase (Phase~2) of the hybrids at 6,030 evaluations: {len(S2['cases'])} split cases "
          f"(both data sets; $r=500, 750, 1000$~m with $N=6, 8, 10$ and $N=10, 12, 15$) $\\times$ {S2['seeds']} seeds. "
          "A cycle is one shaking step followed by the complete local search (compass search); the first descent is the "
          "local search from the switch point, before any shaking. Improvement shares refer to the reduction of the "
          "wake loss in Phase~2 in the runs whose switch point is feasible; the remainder is found in a cycle cut off by "
          "the budget.}",
          "\\label{tab:D-vns}", "\\scriptsize\\setlength{\\tabcolsep}{3pt}",
          "\\resizebox{\\ifdim\\width>\\columnwidth\\columnwidth\\else\\width\\fi}{!}{%", "\\begin{tabular}{lccc}", "\\toprule", f" & {heads} \\\\", "\\midrule"]

    def row(name, key, d=1, f=None):
        vals = [S2[a][key] if f is None else f(S2[a]) for a in T2_ALGS]
        L.append(f"{name} & " + " & ".join(fmt(v, d) if not isinstance(v, str) else v for v in vals) + " \\\\")
    row("Runs with feasible switch point", "switch_feasible", 0, lambda s: f"{s['switch_feasible']} of {s['runs']}")
    row("Phase-2 evaluations (median)", "phase2_evals_median", 0)
    row("First descent: evaluations (median)", "descent_evals_median", 0)
    row("\\quad share of Phase-2 evaluations (\\%)", "descent_evals_pct_of_phase2", 1)
    row("Shaking steps per run (median)", "shakes_median", 1)
    row("\\quad range over runs", None, 0, lambda s: f"{s['shakes_min']}--{s['shakes_max']}")
    for nn in sorted(S2[T2_ALGS[0]]["shakes_median_by_n"]):
        row(f"\\quad median for $N={nn}$", None, 1, lambda s, nn=nn: s["shakes_median_by_n"][nn])
    row("Shaking evaluations (\\% of Phase~2; rest: local search)", "shake_evals_pct_of_phase2", 2)
    row("Accepted cycles per run (mean)", "accepted_mean", 1)
    row("\\quad acceptance rate of cycles (\\%)", "accept_rate_pct", 1)
    row("\\quad share accepted with $k=1$ (\\%)", "accepted_k1_pct", 1)
    row("\\quad accepted with $k=1,\\dots,5$", None, 0, lambda s: "/".join(str(v) for v in s["accepted_by_k"]))
    row("\\quad runs without an accepted cycle", "runs_no_accepted_cycle", 0)
    row("\\quad shaken point already better (\\% of accepted)", "shake_better_pct_of_accepted", 1)
    row("Improvement from the first descent (\\%)", "imp_descent_share_pct", 1)
    row("Improvement from accepted cycles (\\%)", "imp_cycles_share_pct", 1)
    row("Infeasible at switch, feasible after descent", "made_feasible_descent", 0)
    row("Infeasible at switch, feasible after a cycle", "made_feasible_cycle", 0)
    row("Feasible incumbent made infeasible", "lost_feasibility", 0)
    row("Start is not the best feasible Phase-1 layout$^a$", "start_not_best_feasible", 0)
    L += ["\\bottomrule", "\\end{tabular}}",
          "\\\\[2pt]\\parbox{\\columnwidth}{\\scriptsize $^a$ runs in which Phase~1 evaluated at least one feasible "
          "layout but handed an infeasible or worse layout to VNS (the penalized comparison picks the best feasible layout "
          "whenever one exists).}",
          "\\end{table}", ""]
    # Table feasible starts
    P, D = S3["PSOC"], S3["DE"]
    mx = P["mixture"]
    ml = lambda m: f"{fmt(mx[m]['feasible_pct'],1)} / {fmt(mx[m]['outside_pct'],1)} / {fmt(mx[m]['min_spacing_median_m'],0)}"
    L += ["\\begin{table}[!t]", "\\centering",
          f"\\caption{{Why PSO (constriction) and DE do not move from feasible starts: the six largest cases and the Horns Rev~1 "
          f"16-turbine block $\\times$ {S3['seeds']} seeds with feasibility-preserving initialization, 6,030 evaluations; every "
          "candidate (trial) is logged; feasibility with the 1e-6~m tolerance of the paper. Mixture test: first-update PSO "
          "candidates $\\mathbf P^i+c_2\\,\\mathbf r\\odot(\\mathbf G-\\mathbf P^i)$ between a particle's feasible start and "
          "the global best (200 draws per particle, separate random stream; $c_2=1.496$, independent $\\mathbf r\\sim U(0,1)$ "
          "per coordinate, clipped to the box), with the original turbine labels and after relabelling the turbines of "
          "$\\mathbf P^i$ to their nearest turbines of $\\mathbf G$ (assignment problem); then, relabelled, without overshoot "
          "(weights in $[0,1]$) and with one common weight (a convex combination of the two layouts).}",
          "\\label{tab:D-feasstart}", "\\scriptsize\\setlength{\\tabcolsep}{3pt}",
          "\\resizebox{\\ifdim\\width>\\columnwidth\\columnwidth\\else\\width\\fi}{!}{%", "\\begin{tabular}{lcc}", "\\toprule", " & PSO & DE \\\\", "\\midrule",
          f"Runs with an improved global best & {P['improved_runs']} of {P['runs']} & {D['improved_runs']} of {D['runs']} \\\\",
          f"Candidates evaluated after the initial population & {fmt(P['candidates'])} & {fmt(D['candidates'])} \\\\",
          f"\\quad feasible (\\%) & {fmt(P['feasible_candidates_pct'],2)} & {fmt(D['feasible_candidates_pct'],2)} \\\\",
          f"\\quad feasible, other than re-evaluations of $\\mathbf G$ & {P['feasible_candidates_not_G_reevaluation']} & {D['feasible_candidates']} \\\\",
          f"Infeasible candidates: outside the site only / spacing only / both (\\%) & "
          f"{fmt(P['infeasible_boundary_only_pct'],1)} / {fmt(P['infeasible_spacing_only_pct'],1)} / {fmt(P['infeasible_both_pct'],1)} & "
          f"{fmt(D['infeasible_boundary_only_pct'],1)} / {fmt(D['infeasible_spacing_only_pct'],1)} / {fmt(D['infeasible_both_pct'],1)} \\\\",
          f"Feasible at the first update (\\%) & {fmt(P['first_update_feasible_pct'],1)} & {fmt(D['first_update_feasible_pct'],1)} \\\\",
          f"Median min.\\ spacing at the first update (m; limit 308, Horns Rev 320) & {fmt(P['first_update_min_spacing_median_m'],0)} & {fmt(D['first_update_min_spacing_median_m'],0)} \\\\",
          f"Personal-best updates / accepted DE trials & {P['pbest_updates']} & {D['accepted_trials']} \\\\",
          f"Feasible candidates better than $\\mathbf G$ & {P['feasible_better_than_gbest']} & 0 \\\\",
          f"Runs in which the $\\mathbf G$-particle keeps $\\mathbf V=\\mathbf 0$ & {P['runs_gbest_particle_zero_velocity_throughout']} of {P['runs']} & -- \\\\",
          f"Iterations with exactly one feasible candidate ($\\mathbf G$ re-evaluated) & {fmt(P['iterations_exactly_one_feasible'])} of {fmt(P['iterations'])} & -- \\\\",
          "\\midrule",
          "\\multicolumn{3}{l}{Mixture test (feasible \\% / turbine outside the site \\% / median min.\\ spacing, m):} \\\\",
          f"\\quad first PSO move, original labels & \\multicolumn{{2}}{{l}}{{{ml('raw')}}} \\\\",
          f"\\quad first PSO move, relabelled & \\multicolumn{{2}}{{l}}{{{ml('matched')}}} \\\\",
          f"\\quad relabelled, weights in $[0,1]$ & \\multicolumn{{2}}{{l}}{{{ml('matched_convex')}}} \\\\",
          f"\\quad relabelled, one weight (convex comb.) & \\multicolumn{{2}}{{l}}{{{ml('matched_scalar')}}} \\\\",
          f"Initial layouts: RMS turbine distance to $\\mathbf G$ after relabelling (m, median / min) & {fmt(P['init_matched_rms_median_m'],0)} / {fmt(P['init_matched_rms_min_m'],0)} & \\\\",
          "\\midrule",
          f"All stored runs of the paper (30 seeds): final $\\ne$ best initial layout & {S3['stored']['PSOC']['runs_final_differs_from_initial_best']} of {S3['stored']['PSOC']['runs']} & {S3['stored']['DE']['runs_final_differs_from_initial_best']} of {S3['stored']['DE']['runs']} \\\\",
          "\\bottomrule", "\\end{tabular}}", "\\end{table}", ""]
    # Table RS-VNS Phase-1 replay
    RP = Sm["rs_replay"]
    L += ["\\begin{table}[!t]", "\\centering",
          "\\caption{Phase~1 of RS-VNS replayed (geometry only, the seeded random stream of the runs): feasible samples among "
          "the 3,015 uniform samples in the bounding square, over 30 seeds (90,450 samples per row; the geometry is the same "
          "for both data sets). Inside: all $N$ turbines inside the circle, observed rate and $(\\pi/4)^N$. Disc: Monte Carlo "
          "feasibility rate of the same number of samples uniform in the disc (the Phase-1 distribution of the planned "
          "control RSDVNS). Rows with $N\\le 5$ at 750 and 1000~m omitted (all 30 seeds have feasible samples).}",
          "\\label{tab:D-rsreplay}", "\\scriptsize\\setlength{\\tabcolsep}{3pt}",
          "\\begin{tabular}{rrrrrrr}", "\\toprule",
          "$r$ (m) & $N$ & feasible & seeds with one & inside (\\%) & $(\\pi/4)^N$ (\\%) & disc feasible (\\%) \\\\", "\\midrule"]
    for k, v in RP.items():
        if v["rad"] > 500 and v["n"] <= 5:
            continue
        L.append(f"{v['rad']} & {v['n']} & {fmt(v['feasible'])} & {v['runs_with_feasible']} & {fmt(100*v['p_inside'],2)} & "
                 f"{fmt(100*v['p_inside_theory'],2)} & {100*v['p_feasible_disc']:.3g} \\\\")
    L += ["\\bottomrule", "\\end{tabular}", "\\end{table}", ""]
    with open(path, "w") as fh:
        fh.write("\n".join(L) + "\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--procs", default="2")
    ap.add_argument("--cache", default=None, help="directory for the raw instrumented results (pickle)")
    ap.add_argument("--rerun", default=None, help="with --cache: comma-separated studies to recompute (T1,T2,T3)")
    a = ap.parse_args()
    procs = int(a.procs)
    cache = os.path.join(a.cache, "diag_raw.pkl") if a.cache else None
    rerun = set(a.rerun.split(",")) if a.rerun else set()
    res = []
    if cache and os.path.exists(cache):
        res = [r for r in pickle.load(open(cache, "rb")) if r["kind"] not in rerun]
        print("loaded", len(res), "instrumented runs from", cache, flush=True)
    else:
        rerun = {"T1", "T2", "T3"}
    J = [j for j in all_jobs() if j[0] in rerun]
    if J:
        t0 = time.time()
        with Pool(procs) as pool:
            new = pool.map(job, J, chunksize=1)
        print(len(new), "instrumented runs in", round(time.time() - t0), "s", flush=True)
        res = res + new
        if cache:
            os.makedirs(a.cache, exist_ok=True); pickle.dump(res, open(cache, "wb"))
    res.sort(key=lambda r: (r["kind"], r["alg"], str(r["ds"]), r["rad"], r["n"], r["seed"]))
    S = stored_table()
    cmp = compare_stored(res, S)
    bit = [r["bitident"] for r in res if "bitident" in r]
    ver = dict(instrumented_runs=len(res), stored_found=sum(v is not None for _, v in cmp),
               stored_match=sum(bool(v) for _, v in cmp), stored_mismatch=[list(map(str, k)) for k, v in cmp if v is False],
               stored_missing=[list(map(str, k)) for k, v in cmp if v is None],
               bitident_runs=len(bit), bitident_match=int(sum(bit)),
               note="stored_match: WakeLoss (relative difference <= 1e-13, the precision of the CSV), coordinate string, curve string and call count equal to the stored run "
                    "of the paper; bitident: seed 1 of every case and method rerun uninstrumented in the same process, "
                    "final position array, tracker curve, calls and WakeLoss bit-identical")
    RP = rs_replay()
    chk = [(r["vns"]["p1_any_feasible"], RP[f"{r['rad']}-{r['n']}"]["feasible_per_seed"][r["seed"] - 1] > 0)
           for r in res if r["kind"] == "T2" and r["alg"] == "RSVNS"]
    ver.update(rs_replay_checked=len(chk), rs_replay_agree=int(sum(x == y for x, y in chk)),
               rs_replay_note="RS-VNS runs of T2: 'Phase 1 evaluated a feasible sample' (recorded during the run) equals "
                              "'the replayed stream has a feasible sample'")
    S1, curves = t1_summary(res)
    Sm = dict(verification=ver, T1=S1, T1_stored=stored_violation_summary(), T2=t2_summary(res), T3=t3_summary(res),
              first_feasible=first_feasible_calls(), T4=cpu_estimates(), rs_replay=RP,
              design=dict(T1=dict(cases=T1_CASES, seeds=list(T1_SEEDS), algs=T1_ALGS),
                          T2=dict(cases=T2_CASES, seeds=list(T2_SEEDS), algs=T2_ALGS),
                          T3=dict(cases=T3_CASES, seeds=list(T3_SEEDS), algs=T3_ALGS), budget=BUDGET))
    Sm["T1"]["curves_median"] = {a: {f: curves[a][f]["med"] for f in curves[a]} for a in curves}
    json.dump(Sm, open(os.path.join(HERE, "mpce_summary_diag.json"), "w"), indent=1,
              default=lambda o: o.item() if hasattr(o, "item") else (list(o) if isinstance(o, (tuple, range)) else str(o)))
    figure(curves, S1)
    M = macros(Sm)
    write_tex(M, os.path.join(HERE, "mpce_numbers_diag.tex"))
    supp_tex(Sm, os.path.join(HERE, "mpce_supp_diag.tex"))
    print(json.dumps(ver, indent=1, default=str)[:2000])
    for k, v in M.items():
        print(f"\\{k} = {v}")


if __name__ == "__main__":
    main()

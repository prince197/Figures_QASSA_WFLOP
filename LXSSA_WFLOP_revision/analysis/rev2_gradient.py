"""Review round 2, Reviewer 3 issue 3: gradient-based baselines with EXACT (analytic) gradients on IEA37 Case Study 1.

Usage:  python3 rev2_gradient.py EXP [SHARD NSHARDS] [--procs=P]
        python3 rev2_gradient.py check      # analytic gradient vs central finite differences + timing (c_g)
Output: rev2_<EXP>_s<SHARD>of<NSHARDS>.csv, same columns as iea37_experiments.run_iea (= run_grid of
        mpce_experiments.py) plus FunCalls, GradCalls, CG, Unused.

Experiment
  grad   IEA37 Case Study 1, 16 turbines (r = 1300 m) and 36 turbines (r = 2000 m), 2 methods, 30 seeds (1..30),
         budgets 30,030 and 6,030, random initialization: 2 x 2 x 2 x 30 = 240 runs.

Methods (new labels; every other label falls through to mpce_experiments.run_method unchanged)
  SLSQPX     multistart SLSQP identical to MSSLSQP (extra_baselines.py: same seeded initial population of 30 evaluated
             with the penalised objective, starts = that population sorted by penalised objective and then new uniform
             points, SLSQP on the wake loss with the same spacing / boundary constraints, bounds, maxiter = 200,
             ftol = 1e-9, same best-feasible bookkeeping) except that the objective gradient is the exact analytic
             gradient of the IEA37 AEP (aep_grad below) instead of forward differences, and the constraint Jacobian is
             also given analytically (it was forward-differenced with eps = 1 m; constraint calls are not counted in
             either method, as in MSSLSQP).
  PSOSLSQPX  two-phase hybrid. Phase 1 = PSO with Clerc-Kennedy constriction (authors_optimizers.PSO, w = 0.7298,
             c1 = c2 = 1.49618, N_p = 30), T1 = round((omega B - N_p) / N_p) iterations with omega = 0.5, i.e. exactly
             the Phase 1 of PSO-VNS (PSOBV) and the first 3,030 / 15,030 calls of PSOC (bitwise, same seed and stream).
             Phase 2 = SLSQP with the exact gradient (same constraints/options as SLSQPX) started from the PSO global
             best; then restarts alternating between (a) a perturbed copy of the incumbent (best feasible layout found
             so far; PSO best if none): ceil(n/4) random turbines displaced by N(0, s_min^2) per coordinate, clipped to
             the box, and (b) the next personal best of the final swarm in order of penalised objective (PSO pbest
             memories, reconstructed from the evaluation log without touching PSO); once the swarm members are used up,
             only (a); until the budget is used.

BUDGET ACCOUNTING (identical definition of "Calls" for all methods, reported in the Calls column)
  * every objective evaluation counts 1 (initial population, PSO particles, every SLSQP function evaluation, including
    line-search evaluations; scipy's cache means a repeated point is evaluated once and counted once);
  * every exact gradient evaluation is charged c_g objective evaluations, where c_g = measured wall-clock time of one
    gradient evaluation (aep_grad, which also computes the AEP) divided by the time of one objective evaluation (the
    wake-loss function called by SLSQP), rounded UP to an integer (CG below, measured by `check`; the common alternative
    charge c_g = 1 is available as label suffix "1", e.g. SLSQPX1 / PSOSLSQPX1, not part of `grad`);
  * the convergence Tracker is advanced by c_g calls per gradient (best value unchanged by a gradient), so the Curve
    checkpoints are on the same charged-call axis as for every other method;
  * a gradient requested with fewer than c_g calls left is not computed: the remaining calls are charged as unused
    (column Unused, < c_g) and the run stops, so Calls = FunCalls + c_g GradCalls + Unused = B exactly;
  * constraint and constraint-Jacobian evaluations are free (as in MSSLSQP).
  For reference, MS-SLSQP (FD) spends 2n + 1 function calls per iteration (forward differences), i.e. c_g = 2n = 32 / 72.

Non-differentiable points (handled consistently with the model's own branches in iea37_model.dir_power):
  * abeam pairs (|dx| == 0 in a direction frame) contribute no deficit and no derivative; the downwind indicator
    (dx > 0) is piecewise constant, so its derivative is taken as zero (the model is discontinuous at dx = 0, where a
    wake switches on: no gradient method can see that jump);
  * sigma = K |dx| + D/sqrt 8 uses d|dx|/ddx = sign(dx) (dx != 0 whenever the pair contributes);
  * sqrt(1 - CT D^2 / (8 sigma^2)) >= 1/3 > 0 always (sigma >= D/sqrt 8), no singularity;
  * root-sum-square loss_r = sqrt(sum d2): where loss_r == 0 (turbine waked by nobody, or all deficits underflow to 0)
    the turbine is in the rated branch (v = 9.8 m/s) and its derivative is set to 0 (all incoming d2 and their
    derivatives are 0 there anyway);
  * power curve: dP/dv = 3 P_rated (v - v_ci)^2 / (v_r - v_ci)^3 on CUT_IN <= v < RATED, 0 elsewhere (same masks).
"""
from record_io import encode_coordinates

import sys, time, math
import numpy as np, pandas as pd
from multiprocessing import Pool
from scipy.optimize import minimize
import init_hook
from init_hook import init_pop
import mpce_experiments as mpe
import iea37_experiments as ie
import iea37_model as iea
from authors_optimizers import PSO
from extra_baselines import MSSLSQP, BudgetExhausted

NP = mpe.NP
# c_g = ceil(t_grad / t_obj), measured by `python3 rev2_gradient.py check` (2 shared cores, numpy 2.4.4):
# t_grad / t_obj = 2.2-2.3 for both 16 and 36 turbines (median of 7 x 300 calls, random layouts) -> c_g = 3
CG = {16: 3, 36: 3}
OMEGA = 0.5
BUDGETS = (30030, 6030)
SEEDS = range(1, 31)
NEW = ("SLSQPX", "PSOSLSQPX", "SLSQPX1", "PSOSLSQPX1")
_B = iea.CT * iea.D ** 2 / 8.0
_CK = iea.HRS * iea.WIND_FREQ / 1.E6                        # MWh per W, per direction
_LOG = {}


# ------------------------------------------------------------------------------------------ analytic gradient
def aep_grad(xy):
    """(AEP [MWh], dAEP/dxy [(n, 2), MWh per m]) of iea37_model.aep, same equations and branches."""
    xy = np.asarray(xy, float).reshape(-1, 2)
    n = len(xy)
    if n not in iea._PAIRS:
        iu, ju = np.triu_indices(n, 1)
        iea._PAIRS[n] = (iu, ju, np.arange(len(iea.WIND_DIR))[:, None] * n)
    iu, ju, off = iea._PAIRS[n]
    C, S_ = iea._COS[:, None], iea._SIN[:, None]
    x, y = xy[:, 0], xy[:, 1]
    fx = x[None, :] * C - y[None, :] * S_
    fy = x[None, :] * S_ + y[None, :] * C
    dx = fx[:, iu] - fx[:, ju]
    dy = fy[:, iu] - fy[:, ju]
    adx = np.abs(dx)
    sigma = iea.K * adx + iea.D / np.sqrt(8.)
    sq = np.sqrt(1. - iea.CT / (8. * sigma ** 2 / iea.D ** 2))
    a = 1. - sq
    e = np.exp(-0.5 * (dy / sigma) ** 2)
    g = a * e
    d2 = g ** 2
    live = adx != 0.
    d2[~live] = 0.0
    rec = np.where(dx > 0., iu, ju) + off
    tot = np.bincount(rec.ravel(), d2.ravel(), minlength=len(iea.WIND_DIR) * n).reshape(-1, n)
    loss = np.sqrt(tot)
    v = iea.WIND_SPEED * (1. - loss)
    cub = (iea.CUT_IN <= v) & (v < iea.RATED_WS)
    rng = iea.RATED_WS - iea.CUT_IN
    p = np.where(cub, iea.RATED_PWR * ((v - iea.CUT_IN) / rng) ** 3, 0.0)
    p = np.where((iea.RATED_WS <= v) & (v < iea.CUT_OUT), iea.RATED_PWR, p)
    aep = (iea.HRS * (iea.WIND_FREQ * p.sum(-1)) / 1.E6).sum()
    # backward pass
    dpdv = np.where(cub, 3. * iea.RATED_PWR * (v - iea.CUT_IN) ** 2 / rng ** 3, 0.0)
    dA_dloss = -_CK[:, None] * dpdv * iea.WIND_SPEED                                  # (16, n)
    pos = tot > 0.
    dA_dtot = np.where(pos, dA_dloss / (2. * np.where(pos, loss, 1.)), 0.0)          # (16, n)
    W = dA_dtot.ravel()[rec]                                                          # (16, pairs)
    s3 = sigma ** 3
    dd2_dsig = 2. * g * e * (-_B / (s3 * sq) + a * dy ** 2 / s3)
    Gx = np.where(live, W * dd2_dsig * iea.K * np.sign(dx), 0.0)                      # dA/d(dx)
    Gy = np.where(live, W * (-2. * d2 * dy / sigma ** 2), 0.0)                        # dA/d(dy)
    gpx = (Gx * C + Gy * S_).sum(0)                                                   # dA/dx_i per pair (i = iu)
    gpy = (-Gx * S_ + Gy * C).sum(0)
    gx = np.bincount(iu, gpx, n) - np.bincount(ju, gpx, n)
    gy = np.bincount(iu, gpy, n) - np.bincount(ju, gpy, n)
    return aep, np.c_[gx, gy]


def wake_grad_fn(n):
    """Gradient of the wake loss (Ideal - AEP) w.r.t. the flat coordinate vector."""
    def gr(x):
        return -aep_grad(np.asarray(x, float).reshape(-1, 2))[1].ravel()
    return gr


# ------------------------------------------------------------------------------------------ optimizers
class _SLSQPCore(MSSLSQP):
    """MSSLSQP machinery (constraints, bookkeeping) with an exact objective gradient charged c_g calls."""

    def __init__(self, wake_fun, grad_fun, charge, cg, radius, smin, pop_size=30, budget=3030, seed=None,
                 boundary_cons=None):
        super().__init__(wake_fun, radius, smin, pop_size, budget, 1.0, seed, boundary_cons)
        self.grad_fun = grad_fun; self.charge = charge; self.cg = int(cg)
        self.calls = 0; self.fcalls = 0; self.gcalls = 0; self.unused = 0

    def _cons_jac(self, x):
        p = x.reshape(-1, 2)
        n = len(p)
        iu, ju = np.triu_indices(n, 1)
        m = len(iu)
        J = np.zeros((m + n, 2 * n))
        dd = 2. * (p[iu] - p[ju]) / self.smin ** 2                     # d/dp_i of (|p_i - p_j|^2 - s^2)/s^2
        r = np.arange(m)
        J[r, 2 * iu] = dd[:, 0]; J[r, 2 * iu + 1] = dd[:, 1]
        J[r, 2 * ju] = -dd[:, 0]; J[r, 2 * ju + 1] = -dd[:, 1]
        rb = m + np.arange(n)
        J[rb, 2 * np.arange(n)] = -2. * p[:, 0] / self.radius ** 2
        J[rb, 2 * np.arange(n) + 1] = -2. * p[:, 1] / self.radius ** 2
        return J

    def _setup(self, best):
        def feasible(x):
            return self._cons(x).min() >= -1e-9

        def wrapped(x):
            if self.calls >= self.budget:
                raise BudgetExhausted
            self.calls += 1; self.fcalls += 1
            w = self.wake_fun(x)
            if feasible(x):
                if best[1] > 1e9 or w < best[1]:
                    best[0], best[1] = x.copy(), w
            return w

        def grad(x):
            if self.calls + self.cg > self.budget:
                self.unused = self.budget - self.calls
                self.charge(self.unused); self.calls = self.budget
                raise BudgetExhausted
            self.calls += self.cg; self.gcalls += 1
            self.charge(self.cg)
            return self.grad_fun(x)
        cons = [{"type": "ineq", "fun": self._cons}] if self.boundary_cons is not None else \
            [{"type": "ineq", "fun": self._cons, "jac": self._cons_jac}]
        return wrapped, grad, cons

    def _local(self, wrapped, grad, cons, x0, dim, lb, ub):
        minimize(wrapped, x0, jac=grad, method="SLSQP", bounds=[(lb, ub)] * dim, constraints=cons,
                 options=dict(maxiter=200, ftol=1e-9))


class SLSQPX(_SLSQPCore):
    def optimize(self, obj_fun, dim, lb, ub):
        pop = init_pop(self.pop_size, dim, lb, ub)
        fit = np.array([obj_fun(p) for p in pop])
        self.calls = self.fcalls = self.pop_size
        order = list(np.argsort(fit))
        best = [pop[order[0]].copy(), fit[order[0]]]
        wrapped, grad, cons = self._setup(best)
        starts = [pop[i] for i in order]
        while self.calls < self.budget:
            x0 = starts.pop(0) if starts else init_pop(1, dim, lb, ub)[0]
            try:
                self._local(wrapped, grad, cons, x0, dim, lb, ub)
            except BudgetExhausted:
                break
        return best[0], best[1], None


class PSOSLSQPX(_SLSQPCore):
    def __init__(self, *a, split=OMEGA, **k):
        super().__init__(*a, **k)
        self.split = split
        self.iters1 = max(1, int(round((split * self.budget - self.pop_size) / self.pop_size)))

    def optimize(self, obj_fun, dim, lb, ub):
        log = []

        def f1(x):
            self.calls += 1; self.fcalls += 1
            v = obj_fun(x); log.append((x.copy(), v)); return v
        # ---- Phase 1: constriction PSO (continues the seeded stream, as HybridBVNS phase1 = "PSOC") ----
        h, fh, _ = PSO(self.pop_size, self.iters1, w=0.7298, c1=1.49618, c2=1.49618, seed=None).optimize(f1, dim, lb, ub)
        self.phase1_calls = self.calls
        pb = [None] * self.pop_size; pbs = np.full(self.pop_size, np.inf)
        for m, (x, v) in enumerate(log):               # PSO evaluates particle m % N_p at call m (init, then sweeps)
            i = m % self.pop_size
            if v < pbs[i]:
                pbs[i] = v; pb[i] = x
        members = [pb[i] for i in np.argsort(pbs, kind="stable") if not np.array_equal(pb[i], h)]
        best = [np.array(h, float).copy(), fh]
        wrapped, grad, cons = self._setup(best)
        n = dim // 2
        k = max(1, int(math.ceil(n / 4)))
        # ---- Phase 2: SLSQP with exact gradient from the PSO best, then restarts ----
        x0 = np.array(h, float).copy()
        turn = 0
        while self.calls < self.budget:
            try:
                self._local(wrapped, grad, cons, x0, dim, lb, ub)
            except BudgetExhausted:
                break
            if turn % 2 == 0 or not members:
                x0 = best[0].copy()
                for t in np.random.choice(n, k, replace=False):
                    x0[2 * t:2 * t + 2] += np.random.normal(0.0, self.smin, 2)
                x0 = np.clip(x0, lb, ub)
            else:
                x0 = members.pop(0).copy()
            turn += 1
        return best[0], best[1], None


def run_method_g(alg, seed, budget, f, wake, feasible, dim, lb, ub, radius, smin, bcons=None):
    """run_method for SLSQPX / PSOSLSQPX (IEA37 only: the gradient is that of iea37_model.aep); any other label is
    passed to mpce_experiments.run_method unchanged."""
    if alg not in NEW:
        return mpe._RUN_METHOD(alg, seed, budget, f, wake, feasible, dim, lb, ub, radius, smin, bcons)
    tr = mpe.Tracker(feasible, max(1, (budget - NP) // 200))

    def fp(x):
        v = f(x); tr.see(x, v); return v

    def wk(x):
        v = wake(x); tr.see(x, v); return v

    def charge(k):                                   # advance the call counter / checkpoints, best unchanged
        for _ in range(k):
            tr.calls += 1
            if tr.calls % tr.check == 0:
                tr.curve.append(tr.best)
    n = dim // 2
    cg = 1 if alg.endswith("1") else CG[n]
    cls = SLSQPX if alg.startswith("SLSQPX") else PSOSLSQPX
    opt = cls(wk, wake_grad_fn(n), charge, cg, radius, smin, NP, budget, seed=seed, boundary_cons=bcons)
    pos, _, _ = opt.optimize(fp, dim, lb, ub)
    assert tr.calls == opt.calls == opt.fcalls + cg * opt.gcalls + opt.unused == budget, (tr.calls, opt.calls)
    _LOG.clear(); _LOG.update(FunCalls=opt.fcalls, GradCalls=opt.gcalls, CG=cg, Unused=opt.unused)
    _LOG["Coordinates"] = encode_coordinates(np.asarray(pos).reshape(-1, 2))
    return pos, tr


def run_iea_g(task):
    """iea37_experiments.run_iea (unchanged) with run_method_g in place of run_method; extra columns appended."""
    g = vars(ie)
    orig = g["run_method"]
    g["run_method"] = run_method_g
    _LOG.clear()
    try:
        row = ie.run_iea(task)
    finally:
        g["run_method"] = orig
    if "Coordinates" in _LOG:
        row["Coordinates"] = _LOG["Coordinates"]
    for k in ("FunCalls", "GradCalls", "CG", "Unused"):
        row[k] = _LOG.get(k, np.nan)
    return row


def tasks(exp):
    if exp == "grad":
        # heaviest first: 36 turbines before 16, larger budget first
        return [(a, n, s, b, "random") for n in (36, 16) for b in BUDGETS for a in ("SLSQPX", "PSOSLSQPX")
                for s in SEEDS]
    raise ValueError(exp)


# ------------------------------------------------------------------------------------------ verification
def check(nlay=10, h=1e-3, seed=0):
    rs = np.random.RandomState(seed)
    for n in (16, 36):
        r = iea.RADIUS[n]
        errs, nerrs, ierrs = [], [], []
        for _ in range(nlay):
            ang = rs.uniform(0, 2 * np.pi, n); rad = r * np.sqrt(rs.uniform(0, 1, n))
            xy = np.c_[rad * np.cos(ang), rad * np.sin(ang)]
            a, gr = aep_grad(xy)
            assert abs(a - iea.aep(xy)) <= 1e-9 * a
            fd = np.zeros_like(xy)
            for i in range(n):
                for c in range(2):
                    xp = xy.copy(); xp[i, c] += h; xm = xy.copy(); xm[i, c] -= h
                    fd[i, c] = (iea.aep(xp) - iea.aep(xm)) / (2 * h)
            big = np.abs(gr) > 1e-6 * np.abs(gr).max()
            errs.append((np.abs(gr - fd)[big] / np.abs(gr)[big]).max())
            nerrs.append(np.linalg.norm(gr - fd) / np.linalg.norm(gr))
            ierrs.append(np.abs(gr - fd).max() / np.abs(gr).max())
        # timing: objective as called by SLSQP (wake loss) vs gradient
        f, ideal, _, _ = iea.make_objective(n)
        x = xy.ravel()
        wake = lambda z: ideal - iea.aep(z.reshape(-1, 2))
        gfun = wake_grad_fn(n)
        tf, tg = [], []
        for _ in range(7):
            t0 = time.perf_counter()
            for _ in range(300):
                wake(x)
            tf.append((time.perf_counter() - t0) / 300)
            t0 = time.perf_counter()
            for _ in range(300):
                gfun(x)
            tg.append((time.perf_counter() - t0) / 300)
        tf, tg = np.median(tf), np.median(tg)
        print(f"n={n}: {nlay} random layouts, central FD h={h} m: max componentwise rel.err (|g_i| > 1e-6 |g|_inf) {max(errs):.2e}, "
              f"max normwise rel.err {max(nerrs):.2e}, max |g - fd|_inf / |g|_inf {max(ierrs):.2e}; t_obj {tf * 1e6:.0f} us, t_grad {tg * 1e6:.0f} us, "
              f"ratio {tg / tf:.2f} -> c_g = {math.ceil(tg / tf)} (CG used: {CG[n]})", flush=True)


if __name__ == "__main__":
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    procs = 4
    for a in sys.argv[1:]:
        if a.startswith("--procs"):
            procs = int(a.split("=")[1]) if "=" in a else int(sys.argv[sys.argv.index(a) + 1])
    exp = args[0]
    if exp == "check":
        check()
        sys.exit()
    shard, nsh = (int(args[1]), int(args[2])) if len(args) >= 3 and args[1].isdigit() else (0, 1)
    tl = tasks(exp)[shard::nsh]
    t0 = time.time()
    with Pool(procs) as pool:
        rows = pool.map(run_iea_g, tl, chunksize=1)
    out = pd.DataFrame(rows)
    fn = f"rev2_{exp}_s{shard}of{nsh}.csv"
    out.to_csv(fn, index=False)
    print(exp, fn, len(out), "runs in", round(time.time() - t0), "s; calls", out.Calls.min(), out.Calls.max(),
          flush=True)

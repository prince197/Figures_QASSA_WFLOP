"""Review round 2 (R1.7, R2.5, AE.6): a standard real-coded genetic algorithm baseline (label GA).

Usage:  python3 rev2_ga.py EXP SHARD NSHARDS [--procs P | --procs=P]
Output: rev2_<EXP>_s<SHARD>of<NSHARDS>.csv, same columns as mpce_experiments.run_grid / run_hr and
        iea37_experiments.run_iea (Algorithm, Dataset, Radius, Turbines, Seed, Budget, Init, Calls, Objective,
        Ideal, WakeLoss, Feasible, MinSpacing, Seconds, Coordinates, Curve).

Experiments (label GA, random initialization, seeds 1..30, N_p = 30)
  ga     68 benchmark cases (mpce_experiments.GRID) x 30 seeds, 6,030 calls                    2,040 runs
  gahr   Horns Rev 1 16-turbine block (run_hr, as hrfix), 30,030 and 6,030 calls x 30 seeds         60 runs
  gaiea  IEA37 Case Study 1, 16 and 36 turbines (iea37_experiments.run_iea), 30,030 and 6,030 calls
         x 30 seeds                                                                              120 runs
The search procedure is unchanged; future coordinate output uses 17 significant digits: the GA label is handled by run_method_ga, and the unchanged
run_grid / run_hr / run_iea are run with run_method_ga in place of run_method (as run_grid_x in
mpce_experiments.py), so objective, penalty, Tracker (every call counts, feasible-best, 201 / 200 convergence
checkpoints), bounds and output columns are exactly those of the other methods.

GA (RCGA): N_p = 30, binary tournament selection, SBX (eta_c = 15, p_c = 0.9, per-variable exchange
probability 0.5), polynomial mutation (eta_m = 20, p_m = 1/n with n = 2 N_turbines variables), clipping to the
box, generational replacement with elitism. Budget: N_p calls for the initial population + N_p calls per
generation, T = (B - N_p) // N_p generations (200 at B = 6,030; 1,000 at B = 30,030), i.e. exactly B calls.
"""
from record_io import encode_coordinates

import sys, time
import numpy as np, pandas as pd
from multiprocessing import Pool
from init_hook import init_pop
import mpce_experiments as mx
import iea37_experiments as ix

_FULL_COORDINATES = {}

NP = mx.NP
S30 = range(1, 31)


class RCGA:
    """Standard real-coded GA (minimisation).

    References
      SBX: K. Deb, R.B. Agrawal, Simulated binary crossover for continuous search space, Complex Systems 9 (1995)
           115-148.
      Polynomial mutation: K. Deb, M. Goyal, A combined genetic adaptive search (GeneAS) for engineering design,
           Computer Science and Informatics 26(4) (1996) 30-45.
      Operator settings and per-variable SBX probability 0.5 as in the NSGA-II reference implementation
           (K. Deb, A. Pratap, S. Agarwal, T. Meyarivan, IEEE TEVC 6(2) (2002) 182-197; nsga2-v1.1.6 realcross).

    Generation: N_p offspring from N_p/2 parent pairs, each parent the winner of a binary tournament (two
    distinct members drawn uniformly, lower objective wins). With probability p_c a pair undergoes SBX: each
    variable independently with probability 0.5 (and |y1 - y2| > 1e-14) gets
        beta_q = (2u)^(1/(eta_c+1))            if u <= 0.5
               = (1/(2(1-u)))^(1/(eta_c+1))     otherwise,
        c1 = 0.5[(1+beta_q) y1 + (1-beta_q) y2],  c2 = 0.5[(1-beta_q) y1 + (1+beta_q) y2]
    (otherwise the parents are copied). Each child variable is then mutated with probability p_m = 1/dim by
        delta = (2u)^(1/(eta_m+1)) - 1          if u < 0.5
              = 1 - (2(1-u))^(1/(eta_m+1))      otherwise,     y <- y + delta (ub - lb),
    and the child is clipped to [lb, ub] (the bound handling of all other methods in the study; the unbounded
    textbook operators are used so that the boundary treatment is the same clipping as for SSA / PSO / DE).
    All N_p offspring are evaluated and replace the population (generational); elitism: if no offspring is
    better than the best individual found so far, the worst offspring is replaced by it. Returns the best
    individual found (always the elite).

    Seeding: np.random.seed(seed) in the constructor, as the authors' SSA / PSO / DE classes, and the initial
    population is the common init_pop draw, so the seed-s initial population is identical to theirs."""

    def __init__(self, pop_size=30, max_iter=200, pc=0.9, eta_c=15.0, eta_m=20.0, pm=None, seed=None):
        self.pop_size = pop_size; self.max_iter = max_iter; self.pc = pc
        self.eta_c = eta_c; self.eta_m = eta_m; self.pm = pm
        if seed is not None:
            np.random.seed(seed)

    def _tournament(self, fit):
        i1, i2 = np.random.choice(self.pop_size, 2, replace=False)
        return i1 if fit[i1] < fit[i2] else i2

    def _sbx(self, y1, y2):
        c1, c2 = y1.copy(), y2.copy()
        if np.random.rand() < self.pc:
            dim = y1.size
            do = (np.random.rand(dim) < 0.5) & (np.abs(y1 - y2) > 1e-14)
            u = np.random.rand(dim)
            e = 1.0 / (self.eta_c + 1.0)
            bq = np.where(u <= 0.5, (2.0 * u) ** e, (1.0 / (2.0 * (1.0 - u))) ** e)
            c1[do] = 0.5 * ((1 + bq[do]) * y1[do] + (1 - bq[do]) * y2[do])
            c2[do] = 0.5 * ((1 - bq[do]) * y1[do] + (1 + bq[do]) * y2[do])
        return c1, c2

    def _mutate(self, y, lb, ub, pm):
        dim = y.size
        m = np.random.rand(dim) < pm
        u = np.random.rand(dim)
        e = 1.0 / (self.eta_m + 1.0)
        delta = np.where(u < 0.5, (2.0 * u) ** e - 1.0, 1.0 - (2.0 * (1.0 - u)) ** e)
        y = y + np.where(m, delta * (ub - lb), 0.0)
        return np.clip(y, lb, ub)

    def optimize(self, obj_fun, dim, lb, ub):
        pop = init_pop(self.pop_size, dim, lb, ub)          # the common seeded draw (same call as SSA/PSO/DE)
        lb = np.broadcast_to(np.asarray(lb, float), (dim,)); ub = np.broadcast_to(np.asarray(ub, float), (dim,))
        pm = 1.0 / dim if self.pm is None else self.pm
        fit = np.array([obj_fun(x) for x in pop])
        b = int(np.argmin(fit)); best, fbest = pop[b].copy(), fit[b]
        curve = [fbest]
        for _ in range(self.max_iter):
            kids = []
            while len(kids) < self.pop_size:
                p1 = pop[self._tournament(fit)]; p2 = pop[self._tournament(fit)]
                c1, c2 = self._sbx(p1, p2)
                kids.append(self._mutate(c1, lb, ub, pm))
                if len(kids) < self.pop_size:
                    kids.append(self._mutate(c2, lb, ub, pm))
            pop = np.array(kids)
            fit = np.array([obj_fun(x) for x in pop])
            b = int(np.argmin(fit))
            if fit[b] < fbest:
                best, fbest = pop[b].copy(), fit[b]
            else:                                   # elitism: best-so-far replaces the worst offspring
                w = int(np.argmax(fit)); pop[w] = best; fit[w] = fbest
            curve.append(fbest)
        return best, fbest, np.array(curve)


def run_method_ga(alg, seed, budget, f, wake, feasible, dim, lb, ub, radius, smin, bcons=None):
    """run_method for the GA label; any other label is passed to mpce_experiments.run_method unchanged."""
    if alg != "GA":
        return mx._RUN_METHOD(alg, seed, budget, f, wake, feasible, dim, lb, ub, radius, smin, bcons)
    tr = mx.Tracker(feasible, max(1, (budget - NP) // 200))

    def fp(x):
        v = f(x); tr.see(x, v); return v
    pos, _, _ = RCGA(NP, (budget - NP) // NP, seed=seed).optimize(fp, dim, lb, ub)
    _FULL_COORDINATES["value"] = encode_coordinates(np.asarray(pos).reshape(-1, 2))
    return pos, tr


def _with_ga(fn, task):
    """Run the unchanged fn (run_grid / run_hr / run_iea) with run_method_ga in place of run_method."""
    _FULL_COORDINATES.clear()
    orig_m, orig_i = mx.run_method, ix.run_method
    mx.run_method = run_method_ga; ix.run_method = run_method_ga
    try:
        row = fn(task)
        if "value" in _FULL_COORDINATES:
            row["Coordinates"] = _FULL_COORDINATES["value"]
        return row
    finally:
        mx.run_method = orig_m; ix.run_method = orig_i


def tasks(exp):
    """List of (function, task) pairs; heaviest first."""
    if exp == "ga":
        cases = sorted(mx.GRID, key=lambda c: (-c[2], -c[1], c[0]))          # most turbines first
        return [(mx.run_grid, ("GA", *c, s, 6030, "random")) for c in cases for s in S30]
    if exp == "gahr":
        return [(mx.run_hr, ("GA", 16, s, b, "random")) for b in (30030, 6030) for s in S30]
    if exp == "gaiea":
        return [(ix.run_iea, ("GA", n, s, b, "random")) for b, n in ((30030, 36), (30030, 16), (6030, 36), (6030, 16))
                for s in S30]
    raise ValueError(exp)


def _call(pair):
    fn, t = pair
    return _with_ga(fn, t)


if __name__ == "__main__":
    argv = sys.argv[1:]
    procs = 2
    args = []
    i = 0
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
    t0 = time.time()
    with Pool(procs) as pool:
        rows = pool.map(_call, tl, chunksize=1)
    out = pd.DataFrame(rows)
    fn = f"rev2_{exp}_s{shard}of{nsh}.csv"
    out.to_csv(fn, index=False)
    print(exp, fn, len(out), "runs in", round(time.time() - t0), "s; calls", out.Calls.min(), out.Calls.max(), flush=True)

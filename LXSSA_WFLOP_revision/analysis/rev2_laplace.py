"""Additional experiment: role and parameterization of the Laplace step of LX-SSA (ablation).

Usage:  python3 rev2_laplace.py EXP [SHARD NSHARDS] [--procs P | --procs=P]
Output: rev2_<EXP>_s<SHARD>of<NSHARDS>.csv  (columns of mpce_experiments.run_grid plus Iters = iterations run)

Experiments
  laplace  6 new LX-SSA variants x 12 split cases (SPLITCASES) x 30 seeds, 6,030 calls, random init = 2,160 runs
  verify   re-runs 4 stored runs of fresh_grid.csv (LXSSA and SSA, case (1, 500, 6), seeds 1-2) through this
           framework and prints whether WakeLoss / Coordinates / Curve are bit-identical (not part of `laplace`)

LX-SSA as published (authors_optimizers.LXSSA), per iteration: leaders (first N_p/2) as in SSA; every follower
evaluates BOTH the SSA midpoint (X_i + X_{i-1})/2 and the Laplace candidate X_i + gamma (H - X_i),
gamma ~ Laplace(phi = 0, chi = 1), and keeps the better; then the whole population is evaluated again
-> N_p + 2 (N_p/2) = 2 N_p = 60 calls per iteration, T = (B - N_p) // (2 N_p) = 100 at B = 6,030.
SSA: N_p = 30 calls per iteration, T = 200.

New variants (labels; nothing in authors_optimizers.py is changed; random-number call order as in LXSSA:
leaders 2 rand() per coordinate, then followers 1 rand() each, in the same order)
  LXNR    LX-SSA without the redundant re-evaluation: leaders evaluated once, each follower keeps the already
          known fitness of the chosen point -> N_p/2 + 2 (N_p/2) = 1.5 N_p = 45 calls per iteration.
          (B - N_p) / 45 = 133.33, so T = ceil(6000 / 45) = 134 (used in the r1 schedule): 133 full iterations
          (5,985 calls) + a 134th iteration in which only the 15 leaders are evaluated (the budget is then
          exhausted; followers keep their position and fitness) -> exactly 6,030 calls. Leaders are evaluated
          before the followers (clipped copies; the follower midpoint uses the unclipped leader exactly as in
          LXSSA), so with the same T the trajectory is identical to LXSSA (checked by `selftest`).
  LXREP   SSA cost (N_p = 30 calls per iteration, T = 200): the follower takes the Laplace candidate
          X_i + gamma (H - X_i) instead of the midpoint (no comparison, no extra calls).
  LXU     LX-SSA as published (60 calls per iteration, T = 100, same comparison) with gamma ~ Uniform(-2, 2)
          instead of Laplace(0, 1): symmetric, E|gamma| = 1 = E|gamma| of Laplace(0, 1), light-tailed. NOTE: U(0, 2) removes the heavy
          tail AND the negative sign (Laplace(0, 1) gives gamma < 0, a move away from H, with probability 0.5).
  LXC025, LXC05, LXC2   LX-SSA as published (60 calls per iteration, T = 100) with chi = 0.25, 0.5, 2 (phi = 0).
The chi = 1 LX-SSA and SSA runs of these cases are the stored LXSSA / SSA rows of fresh_grid.csv (`verify`).
"""
from record_io import encode_coordinates

import sys, time, math
import numpy as np, pandas as pd
from multiprocessing import Pool
import mpce_experiments as M
from mpce_experiments import Tracker, NP, SPLITCASES
from authors_optimizers import LXSSA
from init_hook import init_pop

B = 6030
CHI = {"LXC025": 0.25, "LXC05": 0.5, "LXC2": 2.0}
LAPLACE = ["LXNR", "LXREP", "LXU", "LXC025", "LXC05", "LXC2"]
EVALS_PER_ITER = {"LXNR": 1.5 * NP, "LXREP": NP, "LXU": 2 * NP, "LXC025": 2 * NP, "LXC05": 2 * NP, "LXC2": 2 * NP}
_LOG = {}


class LXNR(LXSSA):
    """LX-SSA without the redundant re-evaluation (1.5 N_p calls per iteration); stops at `budget` calls."""

    def __init__(self, pop_size=50, max_iter=100, phi=0.0, chi=1.0, budget=None, seed=None):
        super().__init__(pop_size, max_iter, phi, chi, seed=seed)
        self.budget = np.inf if budget is None else budget

    def optimize(self, obj_fun, dim, lb, ub):
        calls = [0]

        def f(x):
            calls[0] += 1
            return obj_fun(x)
        pop = init_pop(self.pop_size, dim, lb, ub)
        fitness = np.array([f(ind) for ind in pop])
        idx = np.argsort(fitness); pop = pop[idx]; fitness = fitness[idx]
        food_position = pop[0].copy(); food_fitness = fitness[0]
        curve = [food_fitness]
        self.iters = 0
        for l in range(1, self.max_iter + 1):
            if calls[0] >= self.budget:
                break
            r1 = 2 * np.exp(-(4 * l / self.max_iter) ** 2)
            new_pop = pop.copy(); new_fit = fitness.copy()
            half = self.pop_size // 2
            for i in range(half):                                      # leaders (as LXSSA)
                for j in range(dim):
                    r2 = np.random.rand(); r3 = np.random.rand()
                    step = r1 * ((ub - lb) * r2 + lb)
                    if r3 >= 0.5:
                        new_pop[i, j] = food_position[j] + step
                    else:
                        new_pop[i, j] = food_position[j] - step
            lead = np.clip(new_pop[:half], lb, ub)                     # evaluated once (clipped, as LXSSA's
            for i in range(half):                                      # re-evaluation of the clipped new_pop)
                if calls[0] < self.budget:
                    new_fit[i] = f(lead[i])
                else:                                                  # budget exhausted: keep the old salp
                    lead[i] = pop[i]; new_pop[i] = pop[i]; new_fit[i] = fitness[i]
            for i in range(half, self.pop_size):                       # followers (as LXSSA)
                follower = (pop[i] + new_pop[i - 1]) / 2.0
                gamma = self.laplace_random()
                current_salp = pop[i]
                candidate = current_salp + gamma * (food_position - current_salp)
                follower = np.clip(follower, lb, ub)
                candidate = np.clip(candidate, lb, ub)
                if calls[0] + 2 > self.budget:                         # cannot afford both: keep the old salp
                    continue                                           # (never happens for B = 6,030)
                fc = f(candidate); ff = f(follower)                    # same call order as LXSSA
                if fc < ff:
                    new_pop[i] = candidate; new_fit[i] = fc
                else:
                    new_pop[i] = follower; new_fit[i] = ff
            new_pop = np.clip(new_pop, lb, ub)
            pop = new_pop; fitness = new_fit                          # no re-evaluation
            idx = np.argsort(fitness); pop = pop[idx]; fitness = fitness[idx]
            if fitness[0] < food_fitness:
                food_fitness = fitness[0]; food_position = pop[0].copy()
            curve.append(food_fitness)
            self.iters = l
        return food_position, food_fitness, np.array(curve)


class LXREP(LXSSA):
    """SSA cost: the follower takes the Laplace candidate instead of the SSA midpoint (no comparison)."""

    def optimize(self, obj_fun, dim, lb, ub):
        pop = init_pop(self.pop_size, dim, lb, ub)
        fitness = np.array([obj_fun(ind) for ind in pop])
        idx = np.argsort(fitness); pop = pop[idx]; fitness = fitness[idx]
        food_position = pop[0].copy(); food_fitness = fitness[0]
        curve = [food_fitness]
        for l in range(1, self.max_iter + 1):
            r1 = 2 * np.exp(-(4 * l / self.max_iter) ** 2)
            new_pop = pop.copy()
            half = self.pop_size // 2
            for i in range(half):                                      # leaders (as LXSSA)
                for j in range(dim):
                    r2 = np.random.rand(); r3 = np.random.rand()
                    step = r1 * ((ub - lb) * r2 + lb)
                    if r3 >= 0.5:
                        new_pop[i, j] = food_position[j] + step
                    else:
                        new_pop[i, j] = food_position[j] - step
            for i in range(half, self.pop_size):                       # followers: Laplace move only
                gamma = self.laplace_random()
                current_salp = pop[i]
                new_pop[i] = np.clip(current_salp + gamma * (food_position - current_salp), lb, ub)
            new_pop = np.clip(new_pop, lb, ub)
            pop = new_pop
            fitness = np.array([obj_fun(ind) for ind in pop])
            idx = np.argsort(fitness); pop = pop[idx]; fitness = fitness[idx]
            if fitness[0] < food_fitness:
                food_fitness = fitness[0]; food_position = pop[0].copy()
            curve.append(food_fitness)
        self.iters = self.max_iter
        return food_position, food_fitness, np.array(curve)


class LXU(LXSSA):
    """LX-SSA as published with gamma = phi + 4 U(0, 1) - 2 ~ Uniform(-2, 2) (phi = 0): symmetric like Laplace(0, 1)
    and with the same E|gamma| = 1, but light-tailed (|gamma| <= 2); one rand() per draw as before. Isolates the
    heavy tail of the Laplace step."""

    def laplace_random(self):
        return self.phi + 4.0 * np.random.rand() - 2.0


def make_opt(alg, seed, budget):
    T2 = (budget - NP) // (2 * NP)                                     # 100 at B = 6,030
    if alg == "LXNR":
        return LXNR(NP, math.ceil((budget - NP) / (1.5 * NP)), 0.0, 1.0, budget=budget, seed=seed)
    if alg == "LXREP":
        return LXREP(NP, (budget - NP) // NP, 0.0, 1.0, seed=seed)
    if alg == "LXU":
        return LXU(NP, T2, 0.0, 1.0, seed=seed)
    if alg in CHI:
        return LXSSA(NP, T2, 0.0, CHI[alg], seed=seed)
    raise ValueError(alg)


def run_method_lap(alg, seed, budget, f, wake, feasible, dim, lb, ub, radius, smin, bcons=None):
    """run_method for the laplace labels; any other label (LXSSA, SSA, ...) goes to the unchanged run_method."""
    if alg not in LAPLACE:
        _LOG["Iters"] = {"LXSSA": (budget - NP) // (2 * NP), "SSA": (budget - NP) // NP}.get(alg, np.nan)
        pos, tr = M._RUN_METHOD(alg, seed, budget, f, wake, feasible, dim, lb, ub, radius, smin, bcons)
        _LOG["Coordinates"] = encode_coordinates(np.asarray(pos).reshape(-1, 2))
        return pos, tr
    tr = Tracker(feasible, max(1, (budget - NP) // 200))

    def fp(x):
        v = f(x); tr.see(x, v); return v
    opt = make_opt(alg, seed, budget)
    pos, _, _ = opt.optimize(fp, dim, lb, ub)
    _LOG["Coordinates"] = encode_coordinates(np.asarray(pos).reshape(-1, 2))
    _LOG["Iters"] = getattr(opt, "iters", opt.max_iter)
    return pos, tr


def run_grid_lap(task):
    """mpce_experiments.run_grid (unchanged) with run_method_lap in place of run_method; Iters appended."""
    orig = M.run_method
    M.run_method = run_method_lap
    _LOG.clear()
    try:
        row = M.run_grid(task)
    finally:
        M.run_method = orig
    if "Coordinates" in _LOG:
        row["Coordinates"] = _LOG["Coordinates"]
    row["Iters"] = _LOG.get("Iters", np.nan)
    return row


def tasks(exp):
    S30 = range(1, 31)
    if exp == "laplace":
        cases = sorted(SPLITCASES, key=lambda c: -c[2])                # heaviest (most turbines) first, stable
        return [(ds, r, n, a, s) for (ds, r, n) in cases for a in LAPLACE for s in S30]
    raise ValueError(exp)


def _run(t):
    ds, r, n, a, s = t
    return run_grid_lap((a, ds, r, n, s, B, "random"))


def verify():
    ref = pd.read_csv("fresh_grid.csv")
    for a in ("LXSSA", "SSA"):
        for s in (1, 2):
            row = run_grid_lap((a, 1, 500, 6, s, B, "random"))
            q = ref[(ref.Algorithm == a) & (ref.Dataset == 1) & (ref.Radius == 500) & (ref.Turbines == 6) & (ref.Seed == s)].iloc[0]
            print(a, s, "calls", row["Calls"], "WakeLoss", row["WakeLoss"], q.WakeLoss,
                  "bitwise:", row["WakeLoss"] == q.WakeLoss, ";".join(f"{a:.3f} {b:.3f}" for a,b in np.asarray([[float(z) for z in t.split()] for t in row["Coordinates"].split(";")])) == q.Coordinates, row["Curve"] == q.Curve,
                  flush=True)


def selftest():
    """LXNR with T = 100 and no budget cap must follow LXSSA (T = 100) exactly, with 30 + 100 x 45 calls."""
    rng = np.random.RandomState(0); w = rng.rand(12)
    obj = lambda x: float(np.sum(w * np.sin(x / 50.0) ** 2) + 1e-4 * np.sum(x))
    ncall = [0]

    def cnt(x):
        ncall[0] += 1; return obj(x)
    p1, f1, c1 = LXSSA(NP, 100, 0.0, 1.0, seed=3).optimize(cnt, 12, -500, 500); n1 = ncall[0]; ncall[0] = 0
    p2, f2, c2 = LXNR(NP, 100, 0.0, 1.0, seed=3).optimize(cnt, 12, -500, 500); n2 = ncall[0]
    print("selftest LXNR==LXSSA:", np.array_equal(p1, p2), f1 == f2, np.array_equal(c1, c2), "calls", n1, n2, flush=True)


if __name__ == "__main__":
    argv = sys.argv[1:]
    procs = 4; args = []; k = 0
    while k < len(argv):
        a = argv[k]
        if a.startswith("--procs="):
            procs = int(a.split("=")[1])
        elif a == "--procs":
            procs = int(argv[k + 1]); k += 1
        else:
            args.append(a)
        k += 1
    exp = args[0]
    if exp == "verify":
        verify(); sys.exit()
    if exp == "selftest":
        selftest(); sys.exit()
    shard, nsh = (int(args[1]), int(args[2])) if len(args) >= 3 else (0, 1)
    tl = tasks(exp)[shard::nsh]
    t0 = time.time()
    with Pool(procs) as pool:
        rows = pool.map(_run, tl, chunksize=1)
    out = pd.DataFrame(rows)
    fn = f"rev2_{exp}_s{shard}of{nsh}.csv"
    out.to_csv(fn, index=False)
    print(exp, fn, len(out), "runs in", round(time.time() - t0), "s; calls", out.Calls.min(), out.Calls.max(), flush=True)

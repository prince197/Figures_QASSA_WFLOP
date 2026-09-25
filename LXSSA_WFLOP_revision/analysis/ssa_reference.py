"""Independent reference SSA (Section IV eqs. 12 and 14) on the validated benchmark evaluator.

Used for (i) calibration against the recorded SSA runs, (ii) 4D/5D/6D spacing sensitivity and
(iii) a doubled-budget (6,030-evaluation) SSA reference. LX-SSA is NOT re-implemented here
because its Laplace parameters (phi, chi) are not available in the archive.
"""
import time, platform
import numpy as np, pandas as pd, scipy
from wflop_model import (DATASETS, THETA, R, K, CT, BIN_WIDTH, expected_power_linear)

A = R / K
ALPHA = np.arctan(K)


def batch_objective(X, dataset):
    """X: (P, N, 2) -> objective (P,), ideal (scalar)."""
    k, psi, w = DATASETS[dataset]
    dx = X[:, :, None, 0] - X[:, None, :, 0]
    dy = X[:, :, None, 1] - X[:, None, :, 1]
    c = np.cos(THETA)[None, :, None, None]; s = np.sin(THETA)[None, :, None, None]
    dx, dy = dx[:, None], dy[:, None]
    num = dx * c + dy * s + A
    den = np.sqrt((dx + A * c) ** 2 + (dy + A * s) ** 2)
    beta = np.arccos(np.clip(num / den, -1, 1))
    d = np.abs(dx * c + dy * s)
    n = X.shape[1]
    inwake = (beta < ALPHA) & ~np.eye(n, dtype=bool)[None, None]
    dij = (1 - np.sqrt(1 - CT)) / (1 + K * d / R) ** 2
    delta = np.sqrt(np.sum(np.where(inwake, dij, 0.0) ** 2, axis=3))       # (P, 24, N)
    psi_i = psi[None, :, None] * (1 - delta)
    ep = expected_power_linear(np.full(psi_i.shape, 2.0), psi_i)
    obj = BIN_WIDTH * np.sum(w[None, :, None] * ep, axis=(1, 2))
    ideal = BIN_WIDTH * n * np.sum(w * expected_power_linear(k, psi))
    return obj, ideal


def violation(X, smin):
    d2 = ((X[:, :, None, :] - X[:, None, :, :]) ** 2).sum(-1)
    iu = np.triu_indices(X.shape[1], 1)
    return np.maximum(0.0, smin ** 2 - d2[:, iu[0], iu[1]]).sum(1)


def project(X, r):
    rad = np.sqrt((X ** 2).sum(-1, keepdims=True))
    return np.where(rad > r, X * (r / np.maximum(rad, 1e-12)), X)


def better(o1, v1, o2, v2):
    """Feasibility-first comparison (Section III-C): is (o1,v1) better than (o2,v2)?"""
    f1, f2 = v1 <= 0, v2 <= 0
    return np.where(f1 & ~f2, True, np.where(~f1 & f2, False,
                    np.where(f1 & f2, o1 > o2, v1 < v2)))


def rank_key(o, v):
    return np.lexsort((-o, np.where(v > 0, v, 0.0), v > 0))


def ssa(dataset, radius, n, smin, pop=30, iters=100, seed=0, greedy=False):
    rng = np.random.default_rng(seed)
    X = project(rng.uniform(-radius, radius, (pop, n, 2)), radius)
    o, ideal = batch_objective(X, dataset); v = violation(X, smin); evals = pop
    idx = rank_key(o, v); X, o, v = X[idx], o[idx], v[idx]
    H, Ho, Hv = X[0].copy(), o[0], v[0]
    lb, ub = -radius, radius
    for l in range(1, iters + 1):
        r1 = 2 * np.exp(-(4 * l / iters) ** 2)
        Y = X.copy()
        half = pop // 2
        r2 = rng.random((half, n, 2)); r3 = rng.random((half, n, 2))
        step = r1 * ((ub - lb) * r2 + lb)
        Y[:half] = np.where(r3 >= 0.5, H + step, H - step)
        for i in range(half, pop):
            Y[i] = 0.5 * (Y[i] + Y[i - 1])
        Y = project(Y, radius)
        yo, _ = batch_objective(Y, dataset); yv = violation(Y, smin); evals += pop
        if greedy:
            keep = better(yo, yv, o, v)
            X = np.where(keep[:, None, None], Y, X); o = np.where(keep, yo, o); v = np.where(keep, yv, v)
        else:
            X, o, v = Y, yo, yv
        idx = rank_key(o, v); X, o, v = X[idx], o[idx], v[idx]
        if better(o[:1], v[:1], np.array([Ho]), np.array([Hv]))[0]:
            H, Ho, Hv = X[0].copy(), o[0], v[0]
    return H, Ho, Hv, ideal, evals


CASES = [(1, 500, 4), (1, 750, 8), (1, 1000, 8), (2, 500, 4), (2, 750, 8), (2, 1000, 8)]

if __name__ == "__main__":
    configs = [("4D", 8 * R, 100, False), ("4D", 8 * R, 100, True), ("5D", 10 * R, 100, False),
               ("6D", 12 * R, 100, False), ("4D", 8 * R, 200, False)]
    rows = []
    t0 = time.time()
    for label, smin, iters, greedy in configs:
        for ci, (ds, rad, n) in enumerate(CASES):
            for run in range(1, 31):
                ts = time.perf_counter()
                H, Ho, Hv, ideal, ev = ssa(ds, rad, n, smin, iters=iters, seed=100000 * ci + run, greedy=greedy)
                rows.append(dict(Spacing=label, Greedy=greedy, Iterations=iters, Evaluations=ev,
                                 Dataset=ds, Radius=rad, Turbines=n, Run=run, Objective=Ho,
                                 WakeLoss=ideal - Ho, Violation=Hv, Feasible=Hv <= 0,
                                 Seconds=time.perf_counter() - ts,
                                 Coordinates=";".join(f"{a:.4f} {b:.4f}" for a, b in H)))
        print(label, iters, greedy, "done", round(time.time() - t0, 1), "s", flush=True)
    out = "ssa_reference_runs.csv"
    pd.DataFrame(rows).to_csv(out, index=False)
    print("env:", platform.platform(), platform.processor(), "python", platform.python_version(),
          "numpy", np.__version__, "scipy", scipy.__version__)

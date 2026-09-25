"""Largest turbine count N for which a layout satisfying the circular boundary and a minimum
spacing s could be constructed (multi-start penalty minimisation with analytic gradients).
The result is a constructive lower bound on geometric capacity."""
import numpy as np, pandas as pd
from multiprocessing import Pool
from scipy.optimize import minimize
from wflop_model import R, min_spacing


def _fg(z, n, r, s):
    p = z.reshape(n, 2)
    diff = p[:, None] - p[None]
    d2 = (diff ** 2).sum(-1)
    g = np.maximum(0, s ** 2 - d2)
    np.fill_diagonal(g, 0)
    rad2 = (p ** 2).sum(1)
    b = np.maximum(0, rad2 - r ** 2)
    f = 0.25 * (g ** 2).sum() + (b ** 2).sum()
    grad = -2 * (g[:, :, None] * diff).sum(1) + 4 * (b[:, None] * p)
    return f / s ** 4, grad.ravel() / s ** 4


def try_pack(n, r, s, restarts, seed=0):
    rng = np.random.default_rng(seed)
    for _ in range(restarts):
        ang = rng.uniform(0, 2 * np.pi, n)
        rad = r * np.sqrt(rng.uniform(0, 1, n))
        z0 = np.c_[rad * np.cos(ang), rad * np.sin(ang)].ravel()
        z = minimize(_fg, z0, args=(n, r, 1.0005 * s), jac=True, method="L-BFGS-B",
                     options=dict(maxiter=3000)).x.reshape(n, 2)
        rn = np.sqrt((z ** 2).sum(1))
        z[rn > r] *= (r / rn[rn > r])[:, None]
        if min_spacing(z) >= s - 1e-6:
            return True
    return False


def capacity(args):
    r, lab, s = args
    n = 2
    while try_pack(n + 1, r, s, restarts=200):
        n += 1
    return dict(Radius=r, Spacing=lab, MinDistance=s, MaxNConstructed=n)


if __name__ == "__main__":
    jobs = [(r, lab, s) for r in (1000, 750, 500) for lab, s in (("4D", 8 * R), ("5D", 10 * R), ("6D", 12 * R))]
    with Pool(4) as pool:
        rows = pool.map(capacity, jobs, chunksize=1)
    df = pd.DataFrame(rows).sort_values(["Radius", "MinDistance"])
    df.to_csv("packing_capacity.csv", index=False)
    print(df.to_string(index=False))

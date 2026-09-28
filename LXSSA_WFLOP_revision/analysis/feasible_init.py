"""Feasibility-preserving initialization (geometry only; no objective calls).

Each layout starts from uniform random turbine positions inside the site and is pushed apart
by minimising a smooth packing penalty (spacing overlaps + boundary violations) with L-BFGS-B.
The result is accepted if it satisfies the minimum spacing and the boundary; otherwise a new
random start is drawn. Only the global NumPy stream is used, so runs remain seed-paired.
Works for a circle of radius r or for a convex polygon (counter-clockwise vertices).
"""
import numpy as np
from scipy.optimize import minimize


def _penalty(z, n, s, circle_r, poly):
    p = z.reshape(n, 2)
    diff = p[:, None] - p[None]
    d2 = (diff ** 2).sum(-1)
    g = np.maximum(0.0, s ** 2 - d2)
    np.fill_diagonal(g, 0.0)
    f = 0.25 * (g ** 2).sum()
    grad = -2 * (g[:, :, None] * diff).sum(1)
    if poly is None:
        rad2 = (p ** 2).sum(1)
        b = np.maximum(0.0, rad2 - circle_r ** 2)
        f += (b ** 2).sum()
        grad += 4 * (b[:, None] * p)
    else:
        a, e = poly, np.roll(poly, -1, axis=0) - poly
        nrm = np.c_[e[:, 1], -e[:, 0]] / np.linalg.norm(e, axis=1)[:, None]
        sd = np.einsum("ijk,jk->ij", p[:, None, :] - a[None], nrm)          # (n, edges)
        v = np.maximum(0.0, sd)
        sc = s ** 2                                                          # same units as g
        f += (sc * v ** 2).sum()
        grad += 2 * sc * (v[:, :, None] * nrm[None]).sum(1)
    return f / s ** 4, grad.ravel() / s ** 4


def _uniform_in_site(n, circle_r, poly):
    if poly is None:
        ang = np.random.uniform(0, 2 * np.pi, n)
        rad = circle_r * np.sqrt(np.random.uniform(0, 1, n))
        return np.c_[rad * np.cos(ang), rad * np.sin(ang)]
    lo, hi = poly.min(0), poly.max(0)
    out = []
    a, e = poly, np.roll(poly, -1, axis=0) - poly
    nrm = np.c_[e[:, 1], -e[:, 0]]
    while len(out) < n:
        q = np.random.uniform(lo, hi)
        if (np.einsum("jk,jk->j", q[None] - a, nrm) <= 0).all():
            out.append(q)
    return np.array(out)


def _feasible(p, s, circle_r, poly, tol=1e-6):
    n = len(p)
    if n > 1:
        d = np.sqrt(((p[:, None] - p[None]) ** 2).sum(-1))[np.triu_indices(n, 1)]
        if d.min() < s - tol:
            return False
    if poly is None:
        return np.sqrt((p ** 2).sum(1)).max() <= circle_r + tol
    a, e = poly, np.roll(poly, -1, axis=0) - poly
    nrm = np.c_[e[:, 1], -e[:, 0]] / np.linalg.norm(e, axis=1)[:, None]
    return np.einsum("ijk,jk->ij", p[:, None, :] - a[None], nrm).max() <= tol


def feasible_layout(n, s, circle_r=None, poly=None, max_tries=200):
    """One random feasible layout of n turbines with spacing s."""
    for _ in range(max_tries):
        z0 = _uniform_in_site(n, circle_r, poly).ravel()
        if n == 1:
            return z0.reshape(1, 2)
        z = minimize(_penalty, z0, args=(n, 1.0005 * s, None if circle_r is None else circle_r * (1 - 1e-7), poly),
                     jac=True, method="L-BFGS-B", options=dict(maxiter=3000, ftol=1e-15, gtol=1e-12)).x
        p = z.reshape(n, 2)
        if poly is None:                       # pull back any residual boundary overshoot
            rn = np.sqrt((p ** 2).sum(1))
            p[rn > circle_r] *= (circle_r / rn[rn > circle_r])[:, None]
        if _feasible(p, s, circle_r, poly):
            return p
    raise RuntimeError(f"no feasible layout found for n={n}")


def make_generator(s, circle_r=None, poly=None):
    """Generator for init_hook.GEN: gen(n_pop, dim) -> (n_pop, dim) feasible layouts."""
    def gen(n_pop, dim):
        return np.array([feasible_layout(dim // 2, s, circle_r, poly).ravel() for _ in range(n_pop)])
    return gen

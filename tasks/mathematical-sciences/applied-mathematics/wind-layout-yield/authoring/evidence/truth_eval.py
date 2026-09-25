"""Verifier-side truth evaluator: mean farm power of a layout over the clean hub-height record."""
import numpy as np

D = 77.0
R0 = D / 2
CT = 0.8
K_WAKE = 0.05
A0 = 1.0 - np.sqrt(1.0 - CT)


def load_curve(path):
    d = np.loadtxt(path, delimiter=",", skiprows=1)
    return d[:, 0], d[:, 1]


def farm_power_series(xy, U, theta, curve, chunk=4000):
    v, p = curve
    xy = np.asarray(xy, float)
    n = len(xy)
    out = np.empty(len(U))
    dx = xy[None, :, 0] - xy[:, None, 0]   # [j, i] = x_i - x_j  (j upstream candidate, i target)
    dy = xy[None, :, 1] - xy[:, None, 1]
    for a in range(0, len(U), chunk):
        th = np.deg2rad(theta[a:a + chunk])[:, None, None]
        # wind blows FROM theta -> travels toward theta+180; unit downwind vector:
        ex, ey = -np.sin(th), -np.cos(th)
        down = dx * ex + dy * ey            # downstream distance of i from j
        lat = np.abs(-dx * ey + dy * ex)
        rw = R0 + K_WAKE * down
        inw = (down > 0) & (lat <= rw)
        d = np.where(inw, A0 * (R0 / np.where(inw, rw, 1.0)) ** 2, 0.0)
        tot = np.sqrt((d ** 2).sum(axis=1))  # sum over upstream j -> per target i
        ui = U[a:a + chunk, None] * (1.0 - tot)
        pw = np.interp(ui, v, p, right=0.0)
        pw[ui > v[-1]] = 0.0
        out[a:a + chunk] = pw.sum(axis=1)
    return out


def farm_power(xy, U, theta, curve):
    return float(farm_power_series(xy, U, theta, curve).mean())


def free_power(U, curve):
    v, p = curve
    pw = np.interp(U, v, p, right=0.0)
    pw[U > v[-1]] = 0.0
    return float(pw.mean())

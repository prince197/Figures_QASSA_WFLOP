"""Experimental binned-direction layout optimizer (random search with incremental evaluation)."""
import numpy as np

D = 77.0
R0 = D / 2
CT = 0.8
K_WAKE = 0.05
A0 = 1.0 - np.sqrt(1.0 - CT)
SMIN = 4 * D

S_GRID = np.linspace(0.3, 1.0, 141)


def build_bins(U, theta, curve, width, mode="center"):
    """mode='center': all records in a sector evaluated at the sector centre direction.
    mode='fine': width is small; still centre (fine resolution)."""
    v, p = curve
    nb = int(round(360 / width))
    idx = np.floor(np.mod(theta + width / 2, 360) / width).astype(int) % nb
    centers = np.arange(nb) * width
    w = np.bincount(idx, minlength=nb) / len(U)
    G = np.zeros((nb, len(S_GRID)))
    for b in range(nb):
        u = U[idx == b]
        if len(u) == 0:
            continue
        us = u[None, :] * S_GRID[:, None]
        pw = np.interp(us, v, p, right=0.0)
        pw[us > v[-1]] = 0.0
        G[b] = pw.mean(axis=1)
    keep = w > 0
    return centers[keep], w[keep], G[keep]


class Evaluator:
    def __init__(self, centers, weights, G):
        th = np.deg2rad(centers)
        self.ex, self.ey = -np.sin(th), -np.cos(th)
        self.w = weights
        self.G = G
        self.ds = S_GRID[1] - S_GRID[0]

    def pair_d2(self, xs, ys, xt, yt):
        """deficit^2 of targets t caused by sources s, shape [B, S, T]."""
        dx = xt[None, :] - xs[:, None]
        dy = yt[None, :] - ys[:, None]
        ex = self.ex[:, None, None]
        ey = self.ey[:, None, None]
        down = dx * ex + dy * ey
        lat = np.abs(-dx * ey + dy * ex)
        rw = R0 + K_WAKE * down
        inw = (down > 0) & (lat <= rw)
        return np.where(inw, (A0 * (R0 / np.where(inw, rw, 1.0)) ** 2) ** 2, 0.0)

    def gval(self, s):
        # s: [B, N] speed factor -> per-bin power via linear interp on S_GRID
        f = (np.clip(s, S_GRID[0], 1.0) - S_GRID[0]) / self.ds
        i0 = np.minimum(f.astype(int), len(S_GRID) - 2)
        fr = f - i0
        rows = np.arange(self.G.shape[0])[:, None]
        return self.G[rows, i0] * (1 - fr) + self.G[rows, i0 + 1] * fr

    def full(self, xy):
        x, y = xy[:, 0], xy[:, 1]
        self.D2 = self.pair_d2(x, y, x, y)  # [B, src, tgt]
        return self.power_from(self.D2)

    def power_from(self, D2):
        s = 1.0 - np.sqrt(D2.sum(axis=1))
        return float((self.w[:, None] * self.gval(s)).sum())

    def try_move(self, xy, m, newp):
        x, y = xy[:, 0].copy(), xy[:, 1].copy()
        x[m], y[m] = newp
        D2 = self.D2.copy()
        D2[:, m, :] = self.pair_d2(x[m:m + 1], y[m:m + 1], x, y)[:, 0, :]
        D2[:, :, m] = self.pair_d2(x, y, x[m:m + 1], y[m:m + 1])[:, :, 0]
        return self.power_from(D2), D2


def feasible_point(xy, m, p, R):
    if p[0] ** 2 + p[1] ** 2 > R ** 2:
        return False
    d2 = (xy[:, 0] - p[0]) ** 2 + (xy[:, 1] - p[1]) ** 2
    d2[m] = np.inf
    return d2.min() >= SMIN ** 2


def random_layout(N, R, rng, tries=20000):
    pts = []
    for _ in range(tries):
        r = R * np.sqrt(rng.random())
        a = rng.random() * 2 * np.pi
        p = np.array([r * np.cos(a), r * np.sin(a)])
        if all((p - q) @ (p - q) >= SMIN ** 2 for q in pts):
            pts.append(p)
            if len(pts) == N:
                return np.array(pts)
    raise RuntimeError("could not place")


def optimize(N, R, ev, seed, iters=20000, xy0=None, verbose=False):
    rng = np.random.default_rng(seed)
    xy = random_layout(N, R, rng) if xy0 is None else xy0.copy()
    f = ev.full(xy)
    for it in range(iters):
        frac = it / iters
        m = rng.integers(N)
        if rng.random() < 0.1 * (1 - frac):
            r = R * np.sqrt(rng.random()); a = rng.random() * 2 * np.pi
            newp = np.array([r * np.cos(a), r * np.sin(a)])
        else:
            step = 200 * (1 - frac) + 5
            newp = xy[m] + rng.normal(0, step, 2)
        if not feasible_point(xy, m, newp, R):
            continue
        fn, D2 = ev.try_move(xy, m, newp)
        if fn > f:
            f, xy[m], ev.D2 = fn, newp, D2
        if verbose and it % 5000 == 0:
            print(it, f)
    return xy, f

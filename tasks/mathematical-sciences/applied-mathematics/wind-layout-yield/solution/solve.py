"""Reference solution: hub-height resource from a two-level mast, then wake-aware layout optimisation.

Choices (see README): drop logger gaps and iced-sensor episodes; extrapolate the 60 m cup
to hub height with the shear exponent measured between 40 m and 60 m; keep the wind rose at
1-degree resolution so the optimiser cannot hide turbines between coarse sector centre-lines;
multi-start random search with incremental evaluation; report the yield of the final layout
evaluated record-by-record (no binning) on the cleaned hub-height series.
"""
import json
import os
import sys

import numpy as np

DATA = os.environ.get("WLY_DATA", "/app/data")
OUT = os.environ.get("WLY_OUT", "/app/output")
SEED = int(os.environ.get("WLY_SEED", "20260925"))
N_STARTS = int(os.environ.get("WLY_STARTS", "6"))
ITERS = int(os.environ.get("WLY_ITERS", "60000"))

SITE = json.load(open(os.path.join(DATA, "site.json")))
D = SITE["rotor_diameter_m"]
R0 = D / 2
HUB = SITE["hub_height_m"]
CT = SITE["thrust_coefficient"]
K = SITE["wake_decay_constant"]
N = SITE["n_turbines"]
RB = SITE["boundary_radius_m"]
SMIN = SITE["min_spacing_m"]
A0 = 1.0 - np.sqrt(1.0 - CT)
S_GRID = np.linspace(0.3, 1.0, 141)


# ---------------------------------------------------------------- resource
def load_mast():
    raw = np.genfromtxt(os.path.join(DATA, "mast_timeseries.csv"), delimiter=",", skip_header=1,
                        usecols=(1, 2, 3, 4, 5))
    ws40, ws60, wd, temp = raw[:, 0], raw[:, 1], raw[:, 2], raw[:, 3]
    ok = ~np.isnan(ws40) & ~np.isnan(ws60) & ~np.isnan(wd)
    # iced cups: near-zero speeds on both cups with a frozen vane for >= 3 h; the event
    # must touch sub-zero air somewhere (ice lingers after the air warms up)
    frozen = np.r_[False, np.diff(wd) == 0]
    cand = ok & (ws60 < 1.0) & (ws40 < 1.0)
    flag = np.zeros(len(wd), bool)
    i, n = 0, len(wd)
    while i < n:
        if not cand[i]:
            i += 1
            continue
        j = i + 1
        while j < n and cand[j] and frozen[j]:
            j += 1
        if j - i >= 3 and np.nanmin(temp[i:j]) < 1.0:
            flag[i:j] = True
        i = j
    good = ok & ~flag
    both = good & (ws40 > 3.0) & (ws60 > 3.0)
    alpha = np.log(ws60[both].mean() / ws40[both].mean()) / np.log(60.0 / 40.0)
    U = ws60[good] * (HUB / 60.0) ** alpha
    return U, wd[good], dict(alpha=float(alpha), n_good=int(good.sum()), n_iced=int(flag.sum()),
                             n_total=len(wd))


def load_curve():
    d = np.loadtxt(os.path.join(DATA, "power_curve.csv"), delimiter=",", skiprows=1)
    return d[:, 0], d[:, 1]


def pcurve(u, curve):
    v, p = curve
    out = np.interp(u, v, p)
    out[u > v[-1]] = 0.0
    return out


# ---------------------------------------------------------------- binned wake model
class Model:
    def __init__(self, U, wd, curve, width=1.0):
        nb = int(round(360 / width))
        idx = np.floor(np.mod(wd + width / 2, 360) / width).astype(int) % nb
        w = np.bincount(idx, minlength=nb) / len(U)
        keep = np.where(w > 0)[0]
        th = np.deg2rad(keep * width)
        self.ex, self.ey = -np.sin(th), -np.cos(th)
        self.w = w[keep]
        self.G = np.array([pcurve(U[idx == b][None, :] * S_GRID[:, None], curve).mean(axis=1)
                           for b in keep])
        self.ds = S_GRID[1] - S_GRID[0]
        self.rows = np.arange(len(keep))[:, None]

    def d2(self, xs, ys, xt, yt):
        dx = xt[None, :] - xs[:, None]
        dy = yt[None, :] - ys[:, None]
        ex, ey = self.ex[:, None, None], self.ey[:, None, None]
        down = dx * ex + dy * ey
        lat = np.abs(dy * ex - dx * ey)
        rw = R0 + K * down
        inw = (down > 0) & (lat <= rw)
        return np.where(inw, (A0 * (R0 / np.where(inw, rw, 1.0)) ** 2) ** 2, 0.0)

    def value(self, S):
        s = 1.0 - np.sqrt(np.maximum(S, 0.0))
        f = (np.clip(s, S_GRID[0], 1.0) - S_GRID[0]) / self.ds
        i0 = np.minimum(f.astype(int), len(S_GRID) - 2)
        fr = f - i0
        g = self.G[self.rows, i0] * (1 - fr) + self.G[self.rows, i0 + 1] * fr
        return float((self.w[:, None] * g).sum())

    def reset(self, xy):
        self.xy = xy.copy()
        self.D2 = self.d2(xy[:, 0], xy[:, 1], xy[:, 0], xy[:, 1])   # [B, src, tgt]
        self.S = self.D2.sum(axis=1)
        self.f = self.value(self.S)
        return self.f

    def propose(self, m, p):
        x, y = self.xy[:, 0].copy(), self.xy[:, 1].copy()
        x[m], y[m] = p
        row = self.d2(x[m:m + 1], y[m:m + 1], x, y)[:, 0, :]      # m -> others
        col = self.d2(x, y, x[m:m + 1], y[m:m + 1])[:, :, 0]      # others -> m
        S = self.S - self.D2[:, m, :] + row
        S[:, m] = col.sum(axis=1)
        return self.value(S), (m, p, row, col, S)

    def accept(self, f, st):
        m, p, row, col, S = st
        self.xy[m] = p
        self.D2[:, m, :] = row
        self.D2[:, :, m] = col
        self.S, self.f = S, f


def feasible(xy, m, p):
    if p @ p > RB ** 2:
        return False
    d2 = ((xy - p) ** 2).sum(axis=1)
    d2[m] = np.inf
    return d2.min() >= SMIN ** 2


def random_layout(rng):
    while True:
        pts = []
        for _ in range(20000):
            r, a = RB * np.sqrt(rng.random()), 2 * np.pi * rng.random()
            p = np.array([r * np.cos(a), r * np.sin(a)])
            if all(((p - q) ** 2).sum() >= SMIN ** 2 for q in pts):
                pts.append(p)
                if len(pts) == N:
                    return np.array(pts)


def search(model, xy, rng, iters, step0=250.0, jump=0.1):
    f = model.reset(xy)
    for it in range(iters):
        frac = it / iters
        m = rng.integers(N)
        if rng.random() < jump * (1 - frac):
            r, a = RB * np.sqrt(rng.random()), 2 * np.pi * rng.random()
            p = np.array([r * np.cos(a), r * np.sin(a)])
        else:
            p = model.xy[m] + rng.normal(0, step0 * (1 - frac) + 3.0, 2)
        if not feasible(model.xy, m, p):
            continue
        fn, st = model.propose(m, p)
        if fn > model.f:
            model.accept(fn, st)
        if it % 5000 == 4999:
            model.reset(model.xy)
    model.reset(model.xy)
    return model.xy.copy(), model.f


# ---------------------------------------------------------------- exact yield of a layout
def farm_mean_power(xy, U, wd, curve, chunk=4000):
    dx = xy[None, :, 0] - xy[:, None, 0]
    dy = xy[None, :, 1] - xy[:, None, 1]
    tot = 0.0
    for a in range(0, len(U), chunk):
        th = np.deg2rad(wd[a:a + chunk])[:, None, None]
        ex, ey = -np.sin(th), -np.cos(th)
        down = dx * ex + dy * ey
        lat = np.abs(dy * ex - dx * ey)
        rw = R0 + K * down
        inw = (down > 0) & (lat <= rw)
        d = np.where(inw, A0 * (R0 / np.where(inw, rw, 1.0)) ** 2, 0.0)
        u = U[a:a + chunk, None] * (1.0 - np.sqrt((d ** 2).sum(axis=1)))
        tot += pcurve(u, curve).sum()
    return tot / len(U)


def main():
    rng = np.random.default_rng(SEED)
    U, wd, info = load_mast()
    curve = load_curve()
    model = Model(U, wd, curve, width=1.0)
    best = None
    for s in range(N_STARTS):
        xy, f = search(model, random_layout(rng), rng, ITERS)
        if best is None or f > best[1]:
            best = (xy, f)
    xy, f = search(model, best[0], rng, ITERS, step0=60.0, jump=0.0)
    net = farm_mean_power(xy, U, wd, curve)
    gross = N * pcurve(U, curve).mean()
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, "layout.csv"), "w") as fh:
        fh.write("turbine_id,x_m,y_m\n")
        for i, (x, y) in enumerate(xy, 1):
            fh.write(f"{i},{x:.3f},{y:.3f}\n")
    with open(os.path.join(OUT, "yield.json"), "w") as fh:
        json.dump({"net_mean_power_kw": round(float(net), 2),
                   "gross_mean_power_kw": round(float(gross), 2)}, fh, indent=2)
    print(json.dumps(dict(info, model_power=f, net=net, gross=gross)), file=sys.stderr)


if __name__ == "__main__":
    main()

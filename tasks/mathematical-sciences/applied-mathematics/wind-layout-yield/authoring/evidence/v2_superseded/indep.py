"""Independent implementation (different route from solution/solve.py).

- icing QC: rolling 3-h flat-line test on the 60 m cup plus vane freeze, seeded by sub-zero air
- shear: per-record exponent from the cup pair, sector-median (30 deg) exponent applied per record
- climate model: parametric -- Weibull (MLE) per 10 deg sector for hub speed, and a
  von Mises kernel density (3 deg bandwidth) for direction, evaluated on a 1 deg grid
- optimiser: simulated annealing on the parametric model with full re-evaluation
- yield: reported from the parametric model (not the time series)
"""
import json
import os
import sys

import numpy as np
from scipy import stats
from scipy.special import gamma as G_

DATA = os.environ.get("WLY_DATA", "/app/data")
OUT = os.environ.get("WLY_OUT", "/app/output")
SEED = int(os.environ.get("WLY_SEED", "7"))
SITE = json.load(open(os.path.join(DATA, "site.json")))
D, HUB, CT, K = (SITE[k] for k in ("rotor_diameter_m", "hub_height_m", "thrust_coefficient",
                                   "wake_decay_constant"))
N, RB, SMIN = SITE["n_turbines"], SITE["boundary_radius_m"], SITE["min_spacing_m"]
R0 = D / 2
A0 = 1 - np.sqrt(1 - CT)


def read():
    import pandas as pd
    df = pd.read_csv(os.path.join(DATA, "mast_timeseries.csv"))
    pc = pd.read_csv(os.path.join(DATA, "power_curve.csv")).to_numpy()
    return df, pc


def clean(df):
    df = df.dropna().reset_index(drop=True)
    flat = df.ws_60m.rolling(3, center=True, min_periods=3).apply(lambda s: s.max() - s.min(), raw=True)
    vane = df.wd_58m.rolling(3, center=True, min_periods=3).apply(lambda s: s.max() - s.min(), raw=True)
    seed_ = (flat <= 0.02) & (vane == 0) & (df.ws_60m < 1.0)
    # grow each flagged window to its full run of identical readings; require cold at some point
    ice = np.zeros(len(df), bool)
    idx = np.where(seed_.fillna(False).to_numpy())[0]
    ws, wd, T = df.ws_60m.to_numpy(), df.wd_58m.to_numpy(), df.temp_2m.to_numpy()
    for i in idx:
        if ice[i]:
            continue
        a = i
        while a > 0 and ws[a - 1] == ws[i] and wd[a - 1] == wd[i]:
            a -= 1
        b = i
        while b < len(ws) - 1 and ws[b + 1] == ws[i] and wd[b + 1] == wd[i]:
            b += 1
        if T[a:b + 1].min() < 1.0:
            ice[a:b + 1] = True
    return df[~ice].reset_index(drop=True), int(ice.sum())


def hub_speed(df):
    ok = (df.ws_40m > 2) & (df.ws_60m > 2)
    a_rec = np.log(df.ws_60m / df.ws_40m) / np.log(1.5)
    sec = (np.floor(np.mod(df.wd_58m + 15, 360) / 30).astype(int)) % 12
    med = np.array([np.median(a_rec[ok & (sec == s)]) for s in range(12)])
    return df.ws_60m.to_numpy() * (HUB / 60.0) ** med[sec.to_numpy()], med


def climate(U, wd, nsec=36, bw=3.0):
    grid = np.arange(360.0)
    kap = 1 / np.deg2rad(bw) ** 2
    th = np.deg2rad(wd)
    f = np.zeros(360)
    for chunk in np.array_split(np.arange(len(wd)), 20):
        f += np.exp(kap * (np.cos(np.deg2rad(grid)[:, None] - th[None, chunk]) - 1)).sum(axis=1)
    f /= f.sum()
    width = 360 / nsec
    sec_rec = (np.floor(np.mod(wd + width / 2, 360) / width).astype(int)) % nsec
    kc = np.zeros((nsec, 2))
    for s in range(nsec):
        u = U[(sec_rec == s) & (U > 0)]
        k, _, c = stats.weibull_min.fit(u, floc=0)
        p0 = np.mean(U[sec_rec == s] == 0)
        kc[s] = k, c
    sec_grid = (np.floor(np.mod(grid + width / 2, 360) / width).astype(int)) % nsec
    return grid, f, kc[sec_grid]


def power_tables(kc, pc, sgrid):
    """E[P(s U)] for Weibull U, by Gauss-Legendre on the speed axis."""
    x = np.linspace(0, 40, 1601)
    out = np.zeros((len(kc), len(sgrid)))
    for i, (k, c) in enumerate(kc):
        pdf = stats.weibull_min.pdf(x, k, scale=c)
        for j, s in enumerate(sgrid):
            us = s * x
            p = np.interp(us, pc[:, 0], pc[:, 1])
            p[us > pc[-1, 0]] = 0
            out[i, j] = np.trapezoid(p * pdf, x)
    return out


class Park:
    def __init__(self, grid, f, T, sgrid):
        keep = f > 1e-6
        th = np.deg2rad(grid[keep])
        self.ex, self.ey = -np.sin(th), -np.cos(th)
        self.f, self.T, self.s = f[keep], T[keep], sgrid

    def __call__(self, xy):
        dx = xy[None, :, 0] - xy[:, None, 0]
        dy = xy[None, :, 1] - xy[:, None, 1]
        down = dx[None] * self.ex[:, None, None] + dy[None] * self.ey[:, None, None]
        lat = np.abs(dy[None] * self.ex[:, None, None] - dx[None] * self.ey[:, None, None])
        rw = R0 + K * down
        inw = (down > 0) & (lat <= rw)
        d = np.where(inw, A0 * (R0 / np.where(inw, rw, 1)) ** 2, 0)
        sp = 1 - np.sqrt((d ** 2).sum(axis=1))
        pw = np.array([np.interp(sp[b], self.s, self.T[b]) for b in range(len(self.f))])
        return float((self.f[:, None] * pw).sum())


def ok_layout(xy):
    if ((xy ** 2).sum(1) > RB ** 2).any():
        return False
    d = np.sqrt(((xy[:, None] - xy[None]) ** 2).sum(-1)) + np.eye(len(xy)) * 1e9
    return d.min() >= SMIN


def start_layout(rng):
    # concentric-ring seed: 1 + 7 + 16 = 24 turbines, rings rotated at random
    pts = [np.zeros(2)]
    for r, m in ((460.0, 7), (RB - 0.5, 16)):
        a0 = rng.random() * 2 * np.pi
        for q in range(m):
            a = a0 + 2 * np.pi * q / m
            pts.append([r * np.cos(a), r * np.sin(a)])
    xy = np.array(pts, float)
    assert ok_layout(xy)
    return xy


def anneal(model, rng, iters=24000):
    xy = start_layout(rng)
    f = model(xy)
    best = (xy.copy(), f)
    T0 = 0.002 * f
    for it in range(iters):
        T = T0 * (1 - it / iters) ** 2 + 1e-9
        m = rng.integers(N)
        cand = xy.copy()
        cand[m] += rng.normal(0, 150 * (1 - it / iters) + 5, 2)
        if not ok_layout(cand):
            continue
        fc = model(cand)
        if fc > f or rng.random() < np.exp((fc - f) / T):
            xy, f = cand, fc
            if f > best[1]:
                best = (xy.copy(), f)
    return best


def main():
    rng = np.random.default_rng(SEED)
    df, pc = read()
    df, n_ice = clean(df)
    U, med = hub_speed(df)
    grid, f, kc = climate(U, df.wd_58m.to_numpy())
    sgrid = np.linspace(0.3, 1.0, 71)
    model = Park(grid, f, power_tables(kc, pc, sgrid), sgrid)
    best = None
    for _ in range(int(os.environ.get("WLY_STARTS", "2"))):
        r = anneal(model, rng, int(os.environ.get("WLY_ITERS", "24000")))
        if best is None or r[1] > best[1]:
            best = r
    xy, net = best
    gross = N * model.T[:, -1] @ model.f
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, "layout.csv"), "w") as fh:
        fh.write("turbine_id,x_m,y_m\n")
        for i, (x, y) in enumerate(xy, 1):
            fh.write(f"{i},{x:.3f},{y:.3f}\n")
    json.dump({"net_mean_power_kw": float(net), "gross_mean_power_kw": float(gross)},
              open(os.path.join(OUT, "yield.json"), "w"))
    print(json.dumps(dict(n_ice=n_ice, shear=med.round(3).tolist(), net=net, gross=gross)), file=sys.stderr)


if __name__ == "__main__":
    main()

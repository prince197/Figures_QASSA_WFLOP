"""Independent implementation v3 (different route from solution/solve.py at every step).

- icing QC: rolling 3-h flat-line test on the 60 m A cup plus vane freeze, seeded by sub-zero air
- tower shadow: sector-wise A/B ratio decides which 60 m cup to trust; shear per record from the 40/60A pair
- reference homogeneity: CUSUM break in the monthly log-ratio reanalysis/station
- MCP: variance ratio per 45-deg sector on an equal-weight normalised reanalysis+station predictor; veer-corrected reference directions
- climate model: empirical 1-deg direction bins of the analog long-term series
- optimiser: simulated annealing on the parametric model with full re-evaluation
- yield: reported from the binned empirical model
"""
import json
import os
import sys

import numpy as np
import pandas as pd
from scipy import stats
from scipy.special import gamma as G_

DATA = os.environ.get("WLY_DATA", "/root/data")
OUT = os.environ.get("WLY_OUT", "/root/output")
SEED = int(os.environ.get("WLY_SEED", "7"))
SITE = json.load(open(os.path.join(DATA, "site.json")))
D, HUB, CT, K = (SITE[k] for k in ("rotor_diameter_m", "hub_height_m", "thrust_coefficient",
                                   "wake_decay_constant"))
N, RB, SMIN = SITE["n_turbines"], SITE["boundary_radius_m"], SITE["min_spacing_m"]
R0 = D / 2
A0 = 1 - np.sqrt(1 - CT)


def clean_mast():
    df = pd.read_csv(os.path.join(DATA, "mast_timeseries.csv"), parse_dates=["timestamp_utc"])
    df = df.dropna().reset_index(drop=True)
    flat = df.ws_60m_a.rolling(3, center=True, min_periods=3).apply(lambda s: s.max() - s.min(), raw=True)
    vane = df.wd_58m.rolling(3, center=True, min_periods=3).apply(lambda s: s.max() - s.min(), raw=True)
    seed_ = (flat <= 0.02) & (vane == 0) & (df.ws_60m_a < 1.0)
    ice = np.zeros(len(df), bool)
    ws, wd, T = df.ws_60m_a.to_numpy(), df.wd_58m.to_numpy(), df.temp_2m.to_numpy()
    for i in np.where(seed_.fillna(False).to_numpy())[0]:
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
    """Data-driven tower-shadow handling: median A/B ratio per 10-deg vane sector; where one cup
    reads >3% low, use the other, otherwise average. Per-record shear from the same-boom 40/60A
    pair, falling back to the hour-of-day median exponent at low speed."""
    a, b = df.ws_60m_a.to_numpy(), df.ws_60m_b.to_numpy()
    sec = (np.floor(np.mod(df.wd_58m.to_numpy() + 5, 360) / 10).astype(int)) % 36
    good = (a > 3) & (b > 3)
    rat = np.array([np.median(a[good & (sec == s)] / b[good & (sec == s)]) if (good & (sec == s)).sum() > 20 else 1.0
                    for s in range(36)])
    rs = rat[sec]
    u60 = np.where(rs < 0.97, b, np.where(rs > 1.03, a, 0.5 * (a + b)))
    ok = ((df.ws_40m > 2) & (df.ws_60m_a > 2)).to_numpy()
    a_rec = np.log(df.ws_60m_a.clip(lower=0.01) / df.ws_40m.clip(lower=0.01)).to_numpy() / np.log(1.5)
    hod = df.timestamp_utc.dt.hour.to_numpy()
    med = np.array([np.median(a_rec[ok & (hod == h)]) for h in range(24)])
    alpha = np.where(ok, np.clip(a_rec, -0.1, 0.7), med[hod])
    return u60 * (HUB / 60.0) ** alpha


def reference():
    r = pd.read_csv(os.path.join(DATA, "reference_timeseries.csv"), parse_dates=["timestamp_utc"])
    st = pd.read_csv(os.path.join(DATA, "station_timeseries.csv"), parse_dates=["timestamp_utc"])
    r = r.merge(st, on="timestamp_utc")
    mo = r.timestamp_utc.dt.to_period("M")
    q = np.log(r.ws_100m.groupby(mo).mean() / r.ws_10m.groupby(mo).mean()).to_numpy()
    cus = np.cumsum(q - q.mean())                     # CUSUM break location
    k = int(np.argmax(np.abs(cus))) + 1
    f = np.exp(np.median(q[k:]) - np.median(q[:k]))
    months = mo.unique()
    pre = (mo < months[k]).to_numpy()
    r.loc[pre, "ws_100m"] *= f
    return r, str(months[k]), float(f)


def vr_mcp(df, U_hub, r, rng=None):
    """Variance-ratio MCP on an equal-weight predictor: reanalysis and station speeds, each
    normalised by its concurrent-period mean, summed. Eight 45-deg sectors of the
    veer-corrected reanalysis direction; directions from the veer-corrected reanalysis."""
    c = df[["timestamp_utc", "wd_58m"]].assign(U=U_hub).merge(r, on="timestamp_utc")
    veer = np.rad2deg(np.angle(np.mean(np.exp(1j * np.deg2rad(c.wd_58m - c.wd_100m)))))
    sec = lambda x: (np.floor(np.mod(x + veer + 22.5, 360) / 45).astype(int)) % 8
    sc, sr = sec(c.wd_100m.to_numpy()), sec(r.wd_100m.to_numpy())
    mr, ms = c.ws_100m.mean(), c.ws_10m.mean()
    pc = (c.ws_100m / mr + c.ws_10m / ms).to_numpy()
    pl = (r.ws_100m / mr + r.ws_10m / ms).to_numpy()
    U_lt = np.empty(len(r))
    for s in range(8):
        mc, ml = sc == s, sr == s
        x, y = pc[mc], c.U.to_numpy()[mc]
        U_lt[ml] = y.mean() + y.std() / x.std() * (pl[ml] - x.mean())
    D_lt = np.mod(r.wd_100m.to_numpy() + veer, 360.0)
    return np.clip(U_lt, 0, None), D_lt, float(veer)


def empirical_tables(U, wd, pc, sgrid):
    """1-deg direction frequencies and E[P(s*U) | direction bin] from the long-term series."""
    b = np.mod(np.round(wd), 360).astype(int)
    f = np.bincount(b, minlength=360) / len(b)
    T = np.zeros((360, len(sgrid)))
    for d in np.where(f > 0)[0]:
        us = U[b == d][None, :] * sgrid[:, None]
        p = np.interp(us, pc[:, 0], pc[:, 1])
        p[us > pc[-1, 0]] = 0
        T[d] = p.mean(axis=1)
    return np.arange(360.0), f, T


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
    pc = pd.read_csv(os.path.join(DATA, "power_curve.csv")).to_numpy()
    df, n_ice = clean_mast()
    U_hub = hub_speed(df)
    r, brk, f = reference()
    U, wd, veer = vr_mcp(df, U_hub, r, rng=rng)
    sgrid = np.linspace(0.3, 1.0, 71)
    grid, fr, T = empirical_tables(U, wd, pc, sgrid)
    model = Park(grid, fr, T, sgrid)
    best = None
    for _ in range(int(os.environ.get("WLY_STARTS", "2"))):
        res = anneal(model, rng, int(os.environ.get("WLY_ITERS", "24000")))
        if best is None or res[1] > best[1]:
            best = res
    xy, net = best
    gross = N * model.T[:, -1] @ model.f
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, "layout.csv"), "w") as fh:
        fh.write("turbine_id,x_m,y_m\n")
        for i, (x, y) in enumerate(xy, 1):
            fh.write(f"{i},{x:.3f},{y:.3f}\n")
    json.dump({"net_mean_power_kw": float(net), "gross_mean_power_kw": float(gross)},
              open(os.path.join(OUT, "yield.json"), "w"))
    print(json.dumps(dict(n_ice=n_ice, brk=brk, f=f, veer=veer, net=net, gross=gross)), file=sys.stderr)


if __name__ == "__main__":
    main()

"""Reference solution (v3): long-term hub-height climate from a short mast campaign, then
wake-aware layout optimisation.

Choices (see README):
  * mast QC: drop logger gaps and iced-cup episodes (all cups near zero, vane frozen >= 3 h,
    event touching sub-zero air);
  * tower shadow: at 60 m use the cup whose boom faces the wind;
  * shear: power-law exponent by hour of day from the same-boom 40 m / 60 m A pair;
  * reference homogeneity: find the step in the monthly reanalysis/station speed ratio and
    rescale the reanalysis before the step;
  * MCP: veer from the concurrent mast/reanalysis directions; per 30-degree sector of the
    veer-corrected reanalysis direction, a least-squares combination of reanalysis and station
    speed is mapped to hub height by variance ratio (variance-preserving, so the windy
    concurrent period is not regressed into the long term); long-term series 2004-2023;
  * layout: 1-degree wind rose, multi-start random search with incremental evaluation, then a
    simulated-annealing polish;
  * yield: the final layout evaluated hour by hour on the long-term series.
"""
import json
import os
import sys

import numpy as np
import pandas as pd

DATA = os.environ.get("WLY_DATA", "/app/data")
OUT = os.environ.get("WLY_OUT", "/app/output")
SEED = int(os.environ.get("WLY_SEED", "20260925"))
N_STARTS = int(os.environ.get("WLY_STARTS", "6"))
ITERS = int(os.environ.get("WLY_ITERS", "60000"))
SA_ITERS = int(os.environ.get("WLY_SA_ITERS", "150000"))

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
def angdiff(a, b):
    return (np.asarray(a) - np.asarray(b) + 180.0) % 360.0 - 180.0


def mast_hub():
    m = pd.read_csv(os.path.join(DATA, "mast_timeseries.csv"), parse_dates=["timestamp_utc"])
    ws = m[["ws_40m", "ws_60m_a", "ws_60m_b"]].to_numpy()
    wd = m.wd_58m.to_numpy()
    T = m.temp_2m.to_numpy()
    ok = ~np.isnan(ws).any(axis=1) & ~np.isnan(wd)
    cand = ok & (ws < 1.0).all(axis=1)
    frozen = np.r_[False, np.diff(wd) == 0]
    ice = np.zeros(len(m), bool)
    i, n = 0, len(m)
    while i < n:
        if not cand[i]:
            i += 1
            continue
        j = i + 1
        while j < n and cand[j] and frozen[j]:
            j += 1
        if j - i >= 3 and np.nanmin(T[i:j]) < 1.0:
            ice[i:j] = True
        i = j
    d = m[ok & ~ice]
    wd = d.wd_58m.to_numpy()
    ba = SITE["mast"]["ws_60m_a"]["boom_azimuth_deg"]
    bb = SITE["mast"]["ws_60m_b"]["boom_azimuth_deg"]
    use_a = np.abs(angdiff(wd, ba)) <= np.abs(angdiff(wd, bb))
    u60 = np.where(use_a, d.ws_60m_a, d.ws_60m_b)
    pair = ((d.ws_40m > 3) & (d.ws_60m_a > 3)).to_numpy()
    hod = d.timestamp_utc.dt.hour.to_numpy()
    a_h = np.array([np.log(d.ws_60m_a.to_numpy()[pair & (hod == h)].mean()
                           / d.ws_40m.to_numpy()[pair & (hod == h)].mean()) for h in range(24)]) / np.log(1.5)
    U = u60 * (HUB / 60.0) ** a_h[hod]
    info = dict(n_mast=len(m), n_iced=int(ice.sum()), alpha_day_night=[float(a_h[14]), float(a_h[2])])
    return pd.DataFrame({"t": d.timestamp_utc.to_numpy(), "U": U, "wd": wd}), info


def homogenized_reference():
    r = pd.read_csv(os.path.join(DATA, "reference_timeseries.csv"), parse_dates=["timestamp_utc"])
    st = pd.read_csv(os.path.join(DATA, "station_timeseries.csv"), parse_dates=["timestamp_utc"])
    r = r.merge(st, on="timestamp_utc", how="left")
    mo = r.timestamp_utc.dt.to_period("M")
    q = (r.ws_100m.groupby(mo).mean() / r.ws_10m.groupby(mo).mean()).to_numpy()
    best = None
    for k in range(12, len(q) - 12):
        sse = ((q[:k] - q[:k].mean()) ** 2).sum() + ((q[k:] - q[k:].mean()) ** 2).sum()
        if best is None or sse < best[0]:
            best = (sse, k)
    k = best[1]
    months = mo.unique()
    f = q[k:].mean() / q[:k].mean()
    pre = (mo < months[k]).to_numpy()
    r.loc[pre, "ws_100m"] = r.loc[pre, "ws_100m"] * f
    return r, dict(step_month=str(months[k]), step_factor=float(f))


def long_term_series():
    h, info = mast_hub()
    r, hinfo = homogenized_reference()
    c = h.join(r.set_index("timestamp_utc"), on="t", how="inner")
    veer = float(np.rad2deg(np.angle(np.mean(np.exp(1j * np.deg2rad(c.wd - c.wd_100m))))))
    rdir = np.mod(r.wd_100m.to_numpy() + veer, 360)
    cdir = np.mod(c.wd_100m.to_numpy() + veer, 360)
    sec = lambda x: (np.floor(np.mod(x + 15, 360) / 30).astype(int)) % 12
    sr, sc = sec(rdir), sec(cdir)
    # predictor: per-sector least-squares combination of reanalysis and station speed; the
    # long-term mapping is variance-preserving (variance ratio on that predictor), not a
    # regression, so the anomalous concurrent period is not carried into the long term
    XL = np.c_[np.ones(len(r)), r.ws_100m.to_numpy(), r.ws_10m.to_numpy()]
    XC = np.c_[np.ones(len(c)), c.ws_100m.to_numpy(), c.ws_10m.to_numpy()]
    yc = c.U.to_numpy()
    U = np.empty(len(r))
    for s in range(12):
        mc, ml = sc == s, sr == s
        b, *_ = np.linalg.lstsq(XC[mc], yc[mc], rcond=None)
        p, pl = XC[mc] @ b, XL[ml] @ b
        U[ml] = yc[mc].mean() + yc[mc].std() / p.std() * (pl - p.mean())
    info.update(hinfo, veer_deg=veer, n_long_term=len(r))
    return np.clip(U, 0.0, None), rdir, info


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


def anneal(model, xy, rng, iters, T0=2.0, step0=60.0):
    """Metropolis polish on the binned model; keeps the best layout seen."""
    f = model.reset(xy)
    best = (model.xy.copy(), f)
    for it in range(iters):
        frac = it / iters
        T = T0 * (1 - frac) ** 2 + 1e-6
        m = rng.integers(N)
        p = model.xy[m] + rng.normal(0, step0 * (1 - frac) + 2.0, 2)
        if not feasible(model.xy, m, p):
            continue
        fn, st = model.propose(m, p)
        if fn > model.f or rng.random() < np.exp((fn - model.f) / T):
            model.accept(fn, st)
            if model.f > best[1]:
                best = (model.xy.copy(), model.f)
        if it % 5000 == 4999:
            model.reset(model.xy)
    model.reset(best[0])
    return model.xy.copy(), model.f


def main():
    rng = np.random.default_rng(SEED)
    U, wd, info = long_term_series()
    curve = load_curve()
    model = Model(U, wd, curve, width=1.0)
    best = None
    for s in range(N_STARTS):
        xy, f = search(model, random_layout(rng), rng, ITERS)
        if best is None or f > best[1]:
            best = (xy, f)
    xy, f = search(model, best[0], rng, ITERS, step0=60.0, jump=0.0)
    xy, f = anneal(model, xy, rng, SA_ITERS)
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

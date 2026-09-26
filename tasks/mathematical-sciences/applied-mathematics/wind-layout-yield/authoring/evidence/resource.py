"""Long-term hub-height series from mast + reference, with switchable (right and wrong) routes."""
import json
import os

import numpy as np
import pandas as pd

HUB = 80.0


def angdiff(a, b):
    return (np.asarray(a) - np.asarray(b) + 180.0) % 360.0 - 180.0


def load(data):
    m = pd.read_csv(os.path.join(data, "mast_timeseries.csv"), parse_dates=["timestamp_utc"])
    r = pd.read_csv(os.path.join(data, "reference_timeseries.csv"), parse_dates=["timestamp_utc"])
    site = json.load(open(os.path.join(data, "site.json")))
    st = pd.read_csv(os.path.join(data, "station_timeseries.csv"), parse_dates=["timestamp_utc"])
    r = r.merge(st, on="timestamp_utc", how="left")
    return m, r, site


def homogenize(r):
    """Find a single step in the monthly reanalysis/station ratio and rescale the earlier segment."""
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
    out = r.copy()
    pre = mo < months[k]
    out.loc[pre, "ws_100m"] = out.loc[pre, "ws_100m"] * f
    return out, str(months[k]), f


def icing_mask(m):
    ws = m[["ws_40m", "ws_60m_a", "ws_60m_b"]].to_numpy()
    wd = m.wd_58m.to_numpy()
    T = m.temp_2m.to_numpy()
    ok = ~np.isnan(ws).any(axis=1) & ~np.isnan(wd)
    cand = ok & (ws < 1.0).all(axis=1)
    frozen = np.r_[False, np.diff(wd) == 0]
    flag = np.zeros(len(m), bool)
    i, n = 0, len(m)
    while i < n:
        if not cand[i]:
            i += 1
            continue
        j = i + 1
        while j < n and cand[j] and frozen[j]:
            j += 1
        if j - i >= 3 and np.nanmin(T[i:j]) < 1.0:
            flag[i:j] = True
        i = j
    return ok, flag


def mast_hub(m, site, qc=True, cups="boom", shear="diurnal"):
    ok, ice = icing_mask(m)
    good = ok & ~ice if qc else ok
    d = m[good].copy()
    ba = site["mast"]["ws_60m_a"]["boom_azimuth_deg"]
    bb = site["mast"]["ws_60m_b"]["boom_azimuth_deg"]
    wd = d.wd_58m.to_numpy()
    if cups == "boom":        # use the cup whose boom faces the wind (the other sits in the tower wake)
        use_a = np.abs(angdiff(wd, ba)) <= np.abs(angdiff(wd, bb))
        u60 = np.where(use_a, d.ws_60m_a, d.ws_60m_b)
    elif cups == "avg":
        u60 = 0.5 * (d.ws_60m_a + d.ws_60m_b).to_numpy()
    elif cups == "a":
        u60 = d.ws_60m_a.to_numpy()
    elif cups == "max":
        u60 = np.maximum(d.ws_60m_a, d.ws_60m_b).to_numpy()
    # shear from the same-boom pair (40 m and 60 m A share boom A, so shadow cancels in the ratio)
    pair = (d.ws_40m > 3) & (d.ws_60m_a > 3)
    hod = d.timestamp_utc.dt.hour.to_numpy()
    if shear == "diurnal":
        a_h = np.array([np.log(d.ws_60m_a[pair & (hod == h)].mean() / d.ws_40m[pair & (hod == h)].mean())
                        for h in range(24)]) / np.log(1.5)
        alpha = a_h[hod]
    elif shear == "const":
        alpha = np.log(d.ws_60m_a[pair].mean() / d.ws_40m[pair].mean()) / np.log(1.5)
    elif shear == "seventh":
        alpha = 1 / 7
    U = u60 * (HUB / 60.0) ** alpha
    return pd.DataFrame({"t": d.timestamp_utc.to_numpy(), "U": U, "wd": wd})


def long_term(m, r, site, mcp="vr2", veer=True, homog=True, refcol="ws_100m", dircol="wd_100m", **kw):
    h = mast_hub(m, site, **kw)
    if homog and refcol == "ws_100m":
        r, _, _ = homogenize(r)
    r = r.rename(columns={refcol: "x_ref", dircol: "d_ref"})
    if mcp == "none":
        return h.U.to_numpy(), h.wd.to_numpy()
    rr = r.set_index("timestamp_utc")
    c = h.join(rr, on="t", how="inner")
    dv = np.rad2deg(np.angle(np.mean(np.exp(1j * np.deg2rad(c.wd - c.d_ref))))) if veer else 0.0
    rdir = np.mod(r.d_ref.to_numpy() + dv, 360)
    cdir = np.mod(c.d_ref.to_numpy() + dv, 360)
    sec = lambda x: (np.floor(np.mod(x + 15, 360) / 30).astype(int)) % 12
    sr, sc = sec(rdir), sec(cdir)
    x_lt = r.x_ref.to_numpy()
    U_lt = np.empty(len(r))
    for s in range(12):
        mc, ml = sc == s, sr == s
        x, y = c.x_ref.to_numpy()[mc], c.U.to_numpy()[mc]
        if mcp == "vr":
            U_lt[ml] = y.mean() + y.std() / x.std() * (x_lt[ml] - x.mean())
        elif mcp in ("vr2", "ols2_resid", "knn2", "qm2"):
            XC = np.c_[np.ones(mc.sum()), x, c.ws_10m.to_numpy()[mc]]
            XL = np.c_[np.ones(ml.sum()), x_lt[ml], r.ws_10m.to_numpy()[ml]]
            b, *_ = np.linalg.lstsq(XC, y, rcond=None)
            p, pl = XC @ b, XL @ b
            rng = np.random.default_rng(s)
            if mcp == "vr2":
                U_lt[ml] = y.mean() + y.std() / p.std() * (pl - p.mean())
            elif mcp == "qm2":
                qs = np.linspace(0, 1, 401)
                U_lt[ml] = np.interp(pl, np.quantile(p, qs), np.quantile(y, qs))
            else:
                o = np.argsort(p)
                pick = np.clip(np.searchsorted(p[o], pl) + rng.integers(-20, 20, len(pl)), 0, len(p) - 1)
                U_lt[ml] = pl + (y - p)[o][pick] if mcp == "ols2_resid" else y[o][pick]
        elif mcp == "ols":
            b, a = np.polyfit(x, y, 1)
            U_lt[ml] = a + b * x_lt[ml]
        elif mcp == "ratio":
            U_lt[ml] = x_lt[ml] * y.mean() / x.mean()
    return np.clip(U_lt, 0, None), rdir

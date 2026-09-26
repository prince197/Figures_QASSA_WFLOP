"""Seeded generator v3 for wind-layout-yield: short mast campaign + long reanalysis reference.

Truth: 20 years (2004-2023) of clean hourly hub-height (80 m) wind speed and direction at the
site, produced first, before any measurement or reference model. Nothing here reads a solver.

Agent-visible products
  mast_timeseries.csv       2021-2023 met mast: 40 m cup, two 60 m cups on opposite booms
                            (tower shadow), 58 m vane, temperature, RH; icing and outages
  reference_timeseries.csv  2004-2023 reanalysis node at 100 m: smoothed, scaled, noisy,
                            direction veered relative to the site vane
  power_curve.csv, site.json
"""
import argparse
import json
import os

import numpy as np
from scipy.special import ndtr

HUB = 80.0
LT_START = np.datetime64("2004-01-01T00:00")
LT_END = np.datetime64("2024-01-01T00:00")
MAST_START = np.datetime64("2021-01-01T00:00")
BOOM_A = 232.0      # 40 m cup and 60 m cup "A"
BOOM_B = 52.0       # 60 m cup "B"
SITE = {
    "n_turbines": 24,
    "boundary_radius_m": 1000.0,
    "min_spacing_m": 308.0,
    "rotor_diameter_m": 77.0,
    "hub_height_m": HUB,
    "thrust_coefficient": 0.8,
    "wake_decay_constant": 0.05,
    "long_term_period_utc": ["2004-01-01T00:00", "2023-12-31T23:00"],
    "mast": {
        "ws_40m": {"height_m": 40.0, "boom_azimuth_deg": BOOM_A},
        "ws_60m_a": {"height_m": 60.0, "boom_azimuth_deg": BOOM_A},
        "ws_60m_b": {"height_m": 60.0, "boom_azimuth_deg": BOOM_B},
        "wd_58m": {"height_m": 58.0},
        "tower": "triangular lattice, 0.9 m face width",
    },
    "reference": {"source": "reanalysis grid node, 28 km from site", "height_m": 100.0},
    "station": {"source": "airport synoptic station, 14 km from site", "height_m": 10.0,
                "speed_resolution_ms": 0.5, "direction_resolution_deg": 10},
}


def default_cfg():
    return dict(
        lobes=[(232.0, 12.0, 0.50), (38.0, 8.0, 0.28)], w_bg=0.22,
        c0=8.1, c1=1.2, c_dir=232.0, k0=2.15, k1=0.25,
        year_sd=0.03, mast_years_mult=1.045, w2_sd=0.04, mast_w2_shift=0.06,
        alpha0=0.19, alpha_diurnal=0.08, alpha_sd=0.04,
        phi_dir=0.97, phi_spd=0.93, vane_sd=3.0,
        shadow_depth=0.12, shadow_width=14.0,
        ref_scale=0.92, ref_mult_sd=0.12, ref_add_sd=0.9, ref_veer=14.0, ref_dir_sd=4.0,
        ref_step_date='2012-07-01T00:00', ref_step=0.94,
        st_mult_sd=0.15, st_add_sd=0.6, st_dir_sd=12.0,
        ice_rate=0.004, ice_len=(8, 72), gap_blocks=4, gap_len=(6, 96),
    )


def ar1(rng, n, phi):
    e = rng.standard_normal(n)
    z = np.empty(n)
    z[0] = e[0]
    s = np.sqrt(1 - phi ** 2)
    for t in range(1, n):
        z[t] = phi * z[t - 1] + s * e[t]
    return z


def dir_cdf(grid, lobes, w_bg):
    th = np.deg2rad(grid)
    f = np.full_like(th, w_bg / (2 * np.pi))
    for mu, kap, w in lobes:
        f += w * np.exp(kap * np.cos(th - np.deg2rad(mu))) / (2 * np.pi * np.i0(kap))
    c = np.cumsum(f)
    return c / c[-1]


def angdiff(a, b):
    return (a - b + 180.0) % 360.0 - 180.0


def runs(mask):
    d = np.diff(np.r_[0, mask.astype(int), 0])
    return list(zip(np.where(d == 1)[0], np.where(d == -1)[0]))


def generate(seed, cfg=None):
    cfg = cfg or default_cfg()
    rng = np.random.default_rng(seed)
    t = np.arange(LT_START, LT_END, np.timedelta64(1, "h"))
    n = len(t)
    year = t.astype("datetime64[Y]").astype(int) + 1970
    years = np.unique(year)
    mast_years = years >= 2021

    # --- inter-annual variability: speed multiplier and secondary-lobe weight per year
    ymult = 1.0 + cfg["year_sd"] * rng.standard_normal(len(years))
    ymult[mast_years] = cfg["mast_years_mult"] + 0.01 * rng.standard_normal(mast_years.sum())
    ymult[~mast_years] *= 1.0 / ymult[~mast_years].mean()        # long-term mean 1 outside mast years
    w2 = cfg["lobes"][1][2] + cfg["w2_sd"] * rng.standard_normal(len(years))
    w2[mast_years] = cfg["lobes"][1][2] + cfg["mast_w2_shift"] + 0.01 * rng.standard_normal(mast_years.sum())

    # --- true hub direction: latent AR(1) through a year-specific marginal
    grid = np.arange(0, 360, 0.05)
    u = ndtr(ar1(rng, n, cfg["phi_dir"]))
    theta = np.empty(n)
    (mu1, k1, w1), (mu2, k2, _) = cfg["lobes"]
    for yi, y in enumerate(years):
        m = year == y
        bg = 1.0 - w1 - w2[yi]
        cdf = dir_cdf(grid, [(mu1, k1, w1), (mu2, k2, w2[yi])], bg)
        theta[m] = np.interp(u[m], cdf, grid)
    theta = np.mod(theta + rng.uniform(0, 0.05, n), 360.0)

    # --- true hub speed: direction-dependent Weibull, year multiplier
    th = np.deg2rad(theta)
    c = (cfg["c0"] + cfg["c1"] * np.cos(th - np.deg2rad(cfg["c_dir"]))) * ymult[np.searchsorted(years, year)]
    k = cfg["k0"] + cfg["k1"] * np.cos(th - np.deg2rad(cfg["c_dir"]))
    u2 = np.clip(ndtr(ar1(rng, n, cfg["phi_spd"])), 1e-12, 1 - 1e-12)
    U = c * (-np.log1p(-u2)) ** (1.0 / k)

    # --- reference node: 3-h smoothing, scale, multiplicative + additive noise, veer
    Us = np.convolve(np.r_[U[0], U, U[-1]], np.ones(3) / 3, mode="valid")
    ref_ws = cfg["ref_scale"] * Us * (1 + cfg["ref_mult_sd"] * ar1(rng, n, 0.6)) + cfg["ref_add_sd"] * rng.standard_normal(n)
    ref_ws = np.round(np.clip(ref_ws, 0.0, None), 2)
    ref_wd = np.mod(np.round(theta + cfg["ref_veer"] + cfg["ref_dir_sd"] * rng.standard_normal(n)), 360.0)
    # reanalysis inhomogeneity: an input-stream change makes the node read low before the step date
    pre = t < np.datetime64(cfg["ref_step_date"])
    ref_ws = np.round(np.where(pre, ref_ws * cfg["ref_step"], ref_ws), 2)

    # --- airport station 14 km away: 10 m, noisier, directions reported to 10 deg, homogeneous
    st_ws = U * (10 / HUB) ** 0.26 * (1 + cfg["st_mult_sd"] * ar1(rng, n, 0.7)) + cfg["st_add_sd"] * rng.standard_normal(n)
    st_ws = np.round(np.clip(st_ws, 0.0, None) * 2) / 2          # reported to 0.5 m/s
    st_wd = np.mod(np.round((theta - 6 + cfg["st_dir_sd"] * rng.standard_normal(n)) / 10) * 10, 360.0)
    st_wd[st_ws == 0] = 0

    # --- mast campaign
    mm = t >= MAST_START
    tm = t[mm]
    nm = mm.sum()
    Um, thm = U[mm], theta[mm]
    hod = (tm - tm.astype("datetime64[D]")).astype("timedelta64[h]").astype(float)
    alpha = cfg["alpha0"] + cfg["alpha_diurnal"] * np.cos(2 * np.pi * (hod - 2) / 24.0) + cfg["alpha_sd"] * ar1(rng, nm, 0.8)
    alpha = np.clip(alpha, 0.0, 0.5)

    def shadow(boom):
        d = angdiff(thm, boom + 180.0)
        return 1.0 - cfg["shadow_depth"] * np.exp(-(d / cfg["shadow_width"]) ** 2)

    def cup(x):
        y = x * (1 + 0.01 * rng.standard_normal(nm)) + 0.05 * rng.standard_normal(nm)
        return np.round(np.where(y < 0.3, 0.0, y), 2)

    ws40 = cup(Um * (40 / HUB) ** alpha * shadow(BOOM_A))
    ws60a = cup(Um * (60 / HUB) ** alpha * shadow(BOOM_A))
    ws60b = cup(Um * (60 / HUB) ** alpha * shadow(BOOM_B))
    wd = np.mod(np.round(thm + cfg["vane_sd"] * rng.standard_normal(nm)), 360.0)

    doy = (tm - tm.astype("datetime64[Y]")).astype("timedelta64[h]").astype(float) / 24.0
    temp = np.round(5.5 - 11.0 * np.cos(2 * np.pi * (doy - 15) / 365.25)
                    - 3.0 * np.cos(2 * np.pi * (hod - 3) / 24.0) + 3.0 * ar1(rng, nm, 0.97), 1)
    rh = np.clip(np.round(80 + 12 * ar1(rng, nm, 0.9)), 30, 100)
    iced = np.zeros(nm, bool)
    i = 0
    while i < nm:
        if temp[i] < -0.5 and rh[i] > 88 and rng.random() < cfg["ice_rate"] * 6:
            L = rng.integers(*cfg["ice_len"])
            iced[i:i + L] = True
            i += L + 24
        else:
            i += 1
    for a, b in runs(iced):
        s = np.round(rng.uniform(0.0, 0.4), 2)
        ws60a[a:b] = s
        ws60b[a:b] = np.round(max(0.0, s + rng.normal(0, 0.05)), 2)
        ws40[a:b] = np.round(max(0.0, s + rng.normal(0, 0.05)), 2)
        wd[a:b] = wd[a]
    gap = np.zeros(nm, bool)
    for _ in range(cfg["gap_blocks"]):
        a = rng.integers(0, nm - 100)
        gap[a:a + rng.integers(*cfg["gap_len"])] = True

    return dict(t=t, U=U, theta=theta, ref_ws=ref_ws, ref_wd=ref_wd, st_ws=st_ws, st_wd=st_wd, tm=tm,
                ws40=ws40, ws60a=ws60a, ws60b=ws60b, wd=wd, temp=temp, rh=rh,
                iced=iced, gap=gap, ymult=ymult, years=years)


def power_curve_table():
    v = np.arange(0.0, 25.01, 0.5)
    p = 1500.0 / (1.0 + np.exp(-(v - 8.3) / 1.15))
    p[v < 3.5] = 0.0
    p[v >= 13.5] = 1500.0
    return v, np.round(p, 1)


def write(seed, env_dir, truth_path, dev_path=None):
    g = generate(seed)
    os.makedirs(env_dir, exist_ok=True)
    ts = np.datetime_as_string(g["tm"], unit="m")
    with open(os.path.join(env_dir, "mast_timeseries.csv"), "w") as f:
        f.write("timestamp_utc,ws_40m,ws_60m_a,ws_60m_b,wd_58m,temp_2m,rh_2m\n")
        for i in range(len(ts)):
            if g["gap"][i]:
                f.write(f"{ts[i]},,,,,,\n")
            else:
                f.write(f"{ts[i]},{g['ws40'][i]:.2f},{g['ws60a'][i]:.2f},{g['ws60b'][i]:.2f},"
                        f"{int(g['wd'][i])},{g['temp'][i]:.1f},{int(g['rh'][i])}\n")
    tr = np.datetime_as_string(g["t"], unit="m")
    with open(os.path.join(env_dir, "reference_timeseries.csv"), "w") as f:
        f.write("timestamp_utc,ws_100m,wd_100m\n")
        for i in range(len(tr)):
            f.write(f"{tr[i]},{g['ref_ws'][i]:.2f},{int(g['ref_wd'][i])}\n")
    with open(os.path.join(env_dir, "station_timeseries.csv"), "w") as f:
        f.write("timestamp_utc,ws_10m,wd_10m\n")
        for i in range(len(tr)):
            f.write(f"{tr[i]},{g['st_ws'][i]:.1f},{int(g['st_wd'][i])}\n")
    v, p = power_curve_table()
    with open(os.path.join(env_dir, "power_curve.csv"), "w") as f:
        f.write("wind_speed_ms,power_kw\n")
        for a, b in zip(v, p):
            f.write(f"{a:.1f},{b:.1f}\n")
    with open(os.path.join(env_dir, "site.json"), "w") as f:
        json.dump(SITE, f, indent=2)
    os.makedirs(os.path.dirname(truth_path) or ".", exist_ok=True)
    np.savez_compressed(truth_path, U=g["U"], theta=g["theta"])
    if dev_path:
        np.savez_compressed(dev_path, iced=g["iced"], gap=g["gap"], ymult=g["ymult"], years=g["years"])
    return dict(seed=int(seed), n_lt=len(g["t"]), n_mast=len(g["tm"]), iced=float(g["iced"].mean()),
                gap=float(g["gap"].mean()), ymult_mast=g["ymult"][-3:].round(3).tolist())


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, required=True)
    ap.add_argument("--env-dir", required=True)
    ap.add_argument("--truth", required=True)
    ap.add_argument("--dev-masks", default=None)
    a = ap.parse_args()
    print(json.dumps(write(a.seed, a.env_dir, a.truth, a.dev_masks)))

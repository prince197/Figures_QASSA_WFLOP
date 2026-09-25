"""Seeded generator for the wind-layout-yield task.

Produces
  <env_dir>/mast_timeseries.csv   agent-visible met-mast record (noisy, iced, gappy)
  <env_dir>/power_curve.csv       turbine power curve table
  <truth_path>                    hub-height truth (clean speed, true direction), npz

Truth = the clean hub-height wind that the mast record was derived from. It is
produced here, before any measurement model is applied, and never by a solver.
"""
import argparse
import json
import os

import numpy as np

HUB = 80.0
SITE = {
    "n_turbines": 24,
    "boundary_radius_m": 1000.0,
    "min_spacing_m": 308.0,
    "rotor_diameter_m": 77.0,
    "hub_height_m": HUB,
    "thrust_coefficient": 0.8,
    "wake_decay_constant": 0.05,
    "mast_anemometer_heights_m": [40.0, 60.0],
    "mast_vane_height_m": 58.0,
}
Z_LO, Z_HI, Z_VANE = 40.0, 60.0, 58.0
START = np.datetime64("2014-01-01T00:00")
YEARS = 10


def power_curve_table():
    v = np.arange(0.0, 25.01, 0.5)
    p = 1500.0 / (1.0 + np.exp(-(v - 8.3) / 1.15))
    p[v < 3.5] = 0.0
    p[v >= 13.5] = 1500.0
    p = np.round(p, 1)
    return v, p


def direction_marginal(theta_deg, cfg):
    """Unnormalised density of the true hub direction (deg, from-north, clockwise)."""
    th = np.deg2rad(theta_deg)
    f = np.full_like(th, cfg["w_bg"] / (2 * np.pi))
    for mu, kap, w in cfg["lobes"]:
        m = np.deg2rad(mu)
        f += w * np.exp(kap * np.cos(th - m)) / (2 * np.pi * np.i0(kap))
    return f


def weibull_params(theta_deg, cfg):
    th = np.deg2rad(theta_deg)
    c = cfg["c0"] + cfg["c1"] * np.cos(th - np.deg2rad(cfg["c_dir"]))
    k = cfg["k0"] + cfg["k1"] * np.cos(th - np.deg2rad(cfg["c_dir"]))
    return c, k


def ar1(rng, n, phi):
    e = rng.standard_normal(n)
    z = np.empty(n)
    z[0] = e[0]
    s = np.sqrt(1 - phi ** 2)
    for t in range(1, n):
        z[t] = phi * z[t - 1] + s * e[t]
    return z


def default_cfg():
    return dict(
        lobes=[(232.0, 12.0, 0.50), (38.0, 8.0, 0.28)],
        w_bg=0.22,
        c0=8.1, c1=1.2, c_dir=232.0,
        k0=2.15, k1=0.25,
        alpha0=0.20, alpha_sd=0.05,
        phi_dir=0.97, phi_spd=0.93,
        vane_sd=3.0,
        ice_rate=0.004, ice_len=(8, 72),
        gap_blocks=8, gap_len=(6, 96),
    )


def generate(seed, cfg=None):
    from scipy.special import ndtr

    cfg = cfg or default_cfg()
    rng = np.random.default_rng(seed)
    n = int(YEARS * 365.25 * 24)
    t = START + np.arange(n).astype("timedelta64[h]")

    # --- true hub-height direction: latent AR(1) mapped through the marginal CDF
    grid = np.arange(0, 360, 0.05)
    pdf = direction_marginal(grid, cfg)
    cdf = np.cumsum(pdf)
    cdf /= cdf[-1]
    u = ndtr(ar1(rng, n, cfg["phi_dir"]))
    theta = np.interp(u, cdf, grid) + rng.uniform(0, 0.05, n)
    theta = np.mod(theta, 360.0)

    # --- true hub-height speed: Weibull conditional on direction, AR(1) latent
    c, k = weibull_params(theta, cfg)
    u2 = np.clip(ndtr(ar1(rng, n, cfg["phi_spd"])), 1e-12, 1 - 1e-12)
    U = c * (-np.log1p(-u2)) ** (1.0 / k)

    # --- shear: per-hour power-law exponent
    alpha = np.clip(cfg["alpha0"] + cfg["alpha_sd"] * ar1(rng, n, 0.8), 0.02, 0.45)
    true60 = U * (Z_HI / HUB) ** alpha
    true40 = U * (Z_LO / HUB) ** alpha

    def cup(x):
        y = x * (1 + 0.01 * rng.standard_normal(n)) + 0.05 * rng.standard_normal(n)
        y = np.where(y < 0.3, 0.0, y)
        return np.round(y, 2)

    ws60, ws40 = cup(true60), cup(true40)
    wd = np.mod(np.round(theta + cfg["vane_sd"] * rng.standard_normal(n)), 360.0)

    # --- 2 m temperature: seasonal + diurnal + weather noise
    doy = (t - t.astype("datetime64[Y]")).astype("timedelta64[h]").astype(float) / 24.0
    hod = (t - t.astype("datetime64[D]")).astype("timedelta64[h]").astype(float)
    temp = (5.5 - 11.0 * np.cos(2 * np.pi * (doy - 15) / 365.25)
            - 3.0 * np.cos(2 * np.pi * (hod - 3) / 24.0) + 3.0 * ar1(rng, n, 0.97))
    temp = np.round(temp, 1)

    # --- icing: starts only in sub-zero, humid hours; unrelated to wind
    rh = np.clip(np.round(80 + 12 * ar1(rng, n, 0.9)), 30, 100)
    iced = np.zeros(n, bool)
    i = 0
    while i < n:
        if temp[i] < -0.5 and rh[i] > 88 and rng.random() < cfg["ice_rate"] * 6:
            L = rng.integers(*cfg["ice_len"])
            iced[i:i + L] = True
            i += L + 24
        else:
            i += 1
    for a, b in _runs(iced):
        stuck = np.round(rng.uniform(0.0, 0.4), 2)
        ws60[a:b] = stuck
        ws40[a:b] = np.round(max(0.0, stuck + rng.normal(0, 0.05)), 2)
        wd[a:b] = wd[a]

    # --- logger outages
    gap = np.zeros(n, bool)
    for _ in range(cfg["gap_blocks"]):
        a = rng.integers(0, n - 300)
        gap[a:a + rng.integers(*cfg["gap_len"])] = True

    rec = dict(t=t, ws40=ws40, ws60=ws60, wd=wd, temp=temp, rh=rh)
    truth = dict(U=U.astype(np.float64), theta=theta.astype(np.float64))
    meta = dict(seed=int(seed), n=n, iced_frac=float(iced.mean()), gap_frac=float(gap.mean()),
                alpha_mean=float(alpha.mean()))
    return rec, gap, truth, meta, iced


def _runs(mask):
    d = np.diff(np.r_[0, mask.astype(int), 0])
    return list(zip(np.where(d == 1)[0], np.where(d == -1)[0]))


def write(seed, env_dir, truth_path, dev_path=None):
    rec, gap, truth, meta, iced = generate(seed)
    os.makedirs(env_dir, exist_ok=True)
    ts = np.datetime_as_string(rec["t"], unit="m")
    with open(os.path.join(env_dir, "mast_timeseries.csv"), "w") as f:
        f.write("timestamp_utc,ws_40m,ws_60m,wd_58m,temp_2m,rh_2m\n")
        for i in range(len(ts)):
            if gap[i]:
                f.write(f"{ts[i]},,,,,\n")
            else:
                f.write(f"{ts[i]},{rec['ws40'][i]:.2f},{rec['ws60'][i]:.2f},"
                        f"{int(rec['wd'][i])},{rec['temp'][i]:.1f},{int(rec['rh'][i])}\n")
    v, p = power_curve_table()
    with open(os.path.join(env_dir, "power_curve.csv"), "w") as f:
        f.write("wind_speed_ms,power_kw\n")
        for a, b in zip(v, p):
            f.write(f"{a:.1f},{b:.1f}\n")
    with open(os.path.join(env_dir, "site.json"), "w") as f:
        json.dump(SITE, f, indent=2)
    os.makedirs(os.path.dirname(truth_path) or ".", exist_ok=True)
    np.savez_compressed(truth_path, U=truth["U"], theta=truth["theta"])
    if dev_path:
        np.savez_compressed(dev_path, iced=iced, gap=gap)
    return meta


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, required=True)
    ap.add_argument("--env-dir", required=True)
    ap.add_argument("--truth", required=True)
    ap.add_argument("--dev-masks", default=None, help="optional npz of the icing/gap masks (authoring only)")
    a = ap.parse_args()
    print(json.dumps(write(a.seed, a.env_dir, a.truth, a.dev_masks)))

"""Ablation routes built from the reference pipeline's pieces; each skips or swaps one decision."""
import json
import os
import sys

import numpy as np

import solve as S


def mast_variant(qc=True, shear="measured"):
    raw = np.genfromtxt(os.path.join(S.DATA, "mast_timeseries.csv"), delimiter=",", skip_header=1,
                        usecols=(1, 2, 3, 4, 5))
    if qc:
        U, wd, _ = S.load_mast()
        if shear == "measured":
            return U, wd
        # recover ws60 from hub speed to re-extrapolate
    ws40, ws60, wd_all = raw[:, 0], raw[:, 1], raw[:, 2]
    ok = ~np.isnan(ws40) & ~np.isnan(ws60) & ~np.isnan(wd_all)
    if qc:
        # reuse the reference QC mask by re-running its logic
        U_ref, wd_ref, info = S.load_mast()
        both = ok & (ws40 > 3) & (ws60 > 3)
        a = np.log(ws60[both].mean() / ws40[both].mean()) / np.log(1.5)
        ws60_good = U_ref / (S.HUB / 60.0) ** info["alpha"]
        base, wd = ws60_good, wd_ref
    else:
        both = ok & (ws40 > 3) & (ws60 > 3)
        a = np.log(ws60[both].mean() / ws40[both].mean()) / np.log(1.5)
        base, wd = ws60[ok], wd_all[ok]
    alpha = {"measured": a, "none": 0.0, "seventh": 1 / 7}[shear]
    return base * (S.HUB / 60.0) ** alpha, wd


def run(route, out):
    rng = np.random.default_rng(S.SEED)
    curve = S.load_curve()
    qc = route != "no_icing_qc"
    shear = {"no_shear": "none", "shear_one_seventh": "seventh"}.get(route, "measured")
    U, wd = mast_variant(qc, shear)
    width = {"sector30_model": 30, "sector30_series": 30, "sector10_model": 10,
             "sector10_series": 10, "sector5_series": 5}.get(route, 1)
    model = S.Model(U, wd, curve, width=float(width))
    starts, iters = S.N_STARTS, S.ITERS
    if route == "weak_optimiser":
        starts, iters = 1, 3000
    if route in ("random_layout", "ring_layout"):
        if route == "random_layout":
            xy = S.random_layout(rng)
        else:
            import indep
            xy = indep.start_layout(rng)
        f = model.reset(xy)
    else:
        best = None
        for _ in range(starts):
            xy, f = S.search(model, S.random_layout(rng), rng, iters)
            if best is None or f > best[1]:
                best = (xy, f)
        xy, f = S.search(model, best[0], rng, iters, step0=60.0, jump=0.0) if route != "weak_optimiser" else best
    if route.endswith("_model"):
        net = f
    elif route == "sector_mean_speed":
        net = mean_speed_predict(xy, U, wd, curve)
    elif route == "random_layout":
        net = S.N * S.pcurve(U, curve).mean()          # naive: no wake losses claimed
    else:
        net = S.farm_mean_power(xy, U, wd, curve)
    os.makedirs(out, exist_ok=True)
    np.savetxt(os.path.join(out, "layout.csv"), np.c_[np.arange(1, S.N + 1), xy], delimiter=",",
               header="turbine_id,x_m,y_m", comments="", fmt=["%d", "%.3f", "%.3f"])
    json.dump({"net_mean_power_kw": float(net), "gross_mean_power_kw": float(S.N * S.pcurve(U, curve).mean())},
              open(os.path.join(out, "yield.json"), "w"))


def mean_speed_predict(xy, U, wd, curve):
    """Classic shortcut: 12 sectors, power evaluated at each sector's mean hub speed."""
    sec = (np.floor(np.mod(wd + 15, 360) / 30).astype(int)) % 12
    tot = 0.0
    for s in range(12):
        m = sec == s
        if not m.any():
            continue
        Ubar = U[m].mean()
        fake_U = np.full(1, Ubar)
        tot += m.mean() * S.farm_mean_power(xy, fake_U, np.array([30.0 * s]), curve)
    return tot


if __name__ == "__main__":
    run(sys.argv[1], sys.argv[2])

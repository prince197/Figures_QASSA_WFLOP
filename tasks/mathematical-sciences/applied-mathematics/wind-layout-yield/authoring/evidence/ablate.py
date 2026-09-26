"""Ablation routes: the v3 reference pipeline with one decision skipped or swapped."""
import json, os, sys
import numpy as np
route, out = sys.argv[1], sys.argv[2]
import solve as S
import resource as R
kw = {"no_homog": dict(homog=False), "mast_only": dict(mcp="none"), "avg_cups": dict(cups="avg"),
      "cup_a_only": dict(cups="a"), "no_icing_qc": dict(qc=False), "shear_one_seventh": dict(shear="seventh"),
      "no_veer": dict(veer=False), "ols_mcp": dict(mcp="ols"), "ols2_resid_mcp": dict(mcp="ols2_resid"),
      "knn2_mcp": dict(mcp="knn2"), "qm2_mcp": dict(mcp="qm2"), "vr_single_mcp": dict(mcp="vr"), "station_mcp": dict(refcol="ws_10m", dircol="wd_10m", mcp="vr"),
      }.get(route, {})
m, r, site = R.load(S.DATA)
U, wd = R.long_term(m, r, site, **kw)
curve = S.load_curve()
rng = np.random.default_rng(S.SEED)
width = 30.0 if route == "sector30_model" else 1.0
model = S.Model(U, wd, curve, width=width)
if route == "ring_layout":
    pts = [np.zeros(2)]
    for rad, k in ((460.0, 7), (S.RB - 0.5, 16)):
        a0 = rng.random() * 2 * np.pi
        pts += [[rad * np.cos(a0 + 2 * np.pi * q / k), rad * np.sin(a0 + 2 * np.pi * q / k)] for q in range(k)]
    xy = np.array(pts)
    f = model.reset(xy)
elif route == "weak_optimiser":
    xy, f = S.search(model, S.random_layout(rng), rng, 3000)
else:
    best = None
    for _ in range(3):
        xy, f = S.search(model, S.random_layout(rng), rng, 40000)
        if best is None or f > best[1]:
            best = (xy, f)
    xy, f = S.search(model, best[0], rng, 40000, step0=60.0, jump=0.0)
    xy, f = S.anneal(model, xy, rng, 60000)
net = f if route == "sector30_model" else S.farm_mean_power(xy, U, wd, curve)
os.makedirs(out, exist_ok=True)
np.savetxt(os.path.join(out, "layout.csv"), np.c_[np.arange(1, S.N + 1), xy], delimiter=",",
           header="turbine_id,x_m,y_m", comments="", fmt=["%d", "%.3f", "%.3f"])
json.dump({"net_mean_power_kw": float(net), "gross_mean_power_kw": float(S.N * S.pcurve(U, curve).mean())},
          open(os.path.join(out, "yield.json"), "w"))

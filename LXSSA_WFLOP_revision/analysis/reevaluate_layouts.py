"""Re-evaluation of the 720 stored DS1/DS2 follow-up layouts (no re-optimisation).

Outputs: layout_reevaluation.csv (per run) and robustness_summary.csv (per case/algorithm).
"""
import numpy as np, pandas as pd
from wflop_model import farm_objective, parse_coords, min_spacing, R

df = pd.read_csv("../selected_30_run_data.csv")
rows = []
for _, r in df.iterrows():
    xy, ds = parse_coords(r.Coordinates), int(r.Dataset)
    lin, ideal = farm_objective(xy, ds)
    cub, ideal_c = farm_objective(xy, ds, curve="cubic")
    cubco, ideal_cc = farm_objective(xy, ds, curve="cubic_cutout")
    g04, _ = farm_objective(xy, ds, wake="gaussian", kstar=0.04)
    g022, _ = farm_objective(xy, ds, wake="gaussian", kstar=0.022)
    rows.append(dict(Dataset=ds, Radius=r.Radius, Turbines=r.Turbines, Algorithm=r.Algorithm,
                     Seed=r.Seed, Linear=lin, Ideal=ideal, Cubic=cub, IdealCubic=ideal_c,
                     CubicCutout=cubco, IdealCubicCutout=ideal_cc,
                     GaussLoss004=ideal - g04, GaussLoss0022=ideal - g022, JensenLoss=ideal - lin,
                     MinSpacing=min_spacing(xy)))
out = pd.DataFrame(rows)
out["Feas5D"] = out.MinSpacing >= 10 * R - 1e-6
out["Feas6D"] = out.MinSpacing >= 12 * R - 1e-6
out.to_csv("layout_reevaluation.csv", index=False)
g = out.groupby(["Dataset", "Radius", "Turbines", "Algorithm"])
summ = g.agg(Linear=("Linear", "mean"), Cubic=("Cubic", "mean"), CubicCutout=("CubicCutout", "mean"),
             JensenLoss=("JensenLoss", "mean"), GaussLoss004=("GaussLoss004", "mean"),
             GaussLoss0022=("GaussLoss0022", "mean"), Feas5D=("Feas5D", "mean"),
             Feas6D=("Feas6D", "mean"), MinSpacingMin=("MinSpacing", "min")).reset_index()
summ["RelChangeCubicPct"] = 100 * (summ.Cubic / summ.Linear - 1)
summ["RelChangeCubicCutoutPct"] = 100 * (summ.CubicCutout / summ.Linear - 1)
summ.to_csv("robustness_summary.csv", index=False)
pd.set_option("display.width", 250); pd.set_option("display.max_columns", 30)
print(summ.round(3).to_string())
# rank agreement per case
for key, sub in out.groupby(["Dataset", "Radius", "Turbines"]):
    m = sub.groupby("Algorithm")[["Linear", "Cubic", "CubicCutout", "GaussLoss004", "GaussLoss0022"]].mean()
    print(key, "lin:", list(m.Linear.sort_values(ascending=False).index),
          "cubic:", list(m.Cubic.sort_values(ascending=False).index),
          "cubCO:", list(m.CubicCutout.sort_values(ascending=False).index),
          "gauss.04 (low loss first):", list(m.GaussLoss004.sort_values().index),
          "gauss.022:", list(m.GaussLoss0022.sort_values().index))
    from scipy.stats import spearmanr
    print("   spearman lin-cubic", round(spearmanr(sub.Linear, sub.Cubic)[0], 4),
          "jensen-gauss", round(spearmanr(sub.JensenLoss, sub.GaussLoss004)[0], 3))
print("ideal per turbine kW DS1 lin/cub/cubco", out[out.Dataset == 1].iloc[0][["Ideal", "IdealCubic", "IdealCubicCutout"]].values / 4 / 15)
print("ideal per turbine kW DS2 lin/cub/cubco", out[out.Dataset == 2].iloc[0][["Ideal", "IdealCubic", "IdealCubicCutout"]].values / 4 / 15)

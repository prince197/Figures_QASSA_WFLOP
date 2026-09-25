"""Robustness of the fresh results to the benchmark model: every feasible final layout of
fresh_grid.csv is re-evaluated with a cubic power curve (with and without a 25 m/s cut-out)
and with the Bastankhah--Porte-Agel Gaussian wake model (k* = 0.04). No re-optimization."""
import numpy as np, pandas as pd
from multiprocessing import Pool
from scipy.stats import rankdata, kendalltau
from wflop_model import farm_objective

ALGS = ["LXVNS", "LXSSA", "SSA", "PSO", "DE", "VNS", "SLSQP"]
LAB = {"LXVNS": "LX-SSA-VNS", "LXSSA": "LX-SSA", "SSA": "SSA", "PSO": "PSO", "DE": "DE", "VNS": "VNS", "SLSQP": "MS-SLSQP"}


def reval(row):
    xy = np.array([[float(v) for v in p.split()] for p in row["Coordinates"].split(";")])
    ds = int(row["Dataset"])
    lin, il = farm_objective(xy, ds)
    cub, ic = farm_objective(xy, ds, curve="cubic")
    cco, icc = farm_objective(xy, ds, curve="cubic_cutout")
    gau, _ = farm_objective(xy, ds, wake="gaussian", kstar=0.04)
    return dict(Linear=lin, Cubic=cub, CubicCutout=cco, Gauss=gau, IdealCubic=ic, IdealCubicCutout=icc)


if __name__ == "__main__":
    G = pd.concat([pd.read_csv("fresh_grid.csv"), pd.read_csv("fresh_hgrid.csv")], ignore_index=True)
    G = G[G.Feasible].reset_index(drop=True)
    with Pool(4) as pool:
        R = pd.DataFrame(pool.map(reval, G.to_dict("records"), chunksize=50))
    G = pd.concat([G, R], axis=1)
    G.drop(columns=["Coordinates", "Curve"]).to_csv("hybrid_reevaluation.csv", index=False)
    out, lines = {}, []
    for col, lab in (("Linear", "Benchmark (linear curve, Jensen)"), ("Cubic", "Cubic power curve"),
                     ("CubicCutout", "Cubic curve with 25 m/s cut-out"), ("Gauss", "Gaussian wake ($k^*=0.04$)")):
        m = G.groupby(["Dataset", "Radius", "Turbines", "Algorithm"])[col].mean().unstack()[ALGS]
        cnt = G.groupby(["Dataset", "Radius", "Turbines", "Algorithm"]).size().unstack()[ALGS].fillna(0)
        m = m.where(cnt >= 5)                      # only methods with >= 5 feasible runs
        ranks = np.vstack([rankdata(-np.nan_to_num(r, nan=-np.inf)) for r in m.values])
        base = G.groupby(["Dataset", "Radius", "Turbines", "Algorithm"])["Linear"].mean().unstack()[ALGS].where(cnt >= 5)
        taus = [kendalltau(a[np.isfinite(a) & np.isfinite(b)], b[np.isfinite(a) & np.isfinite(b)])[0]
                for a, b in zip(base.values, m.values)]
        same_top = np.mean([np.nanargmax(a) == np.nanargmax(b) for a, b in zip(base.values, m.values)])
        rel = {ds: 100 * (G[G.Dataset == ds][col] / G[G.Dataset == ds]["Linear"] - 1).mean() for ds in (1, 2)}
        out[col] = dict(avg_rank=dict(zip(ALGS, ranks.mean(0))), tau=float(np.nanmean(taus)),
                        same_top=float(same_top), rel=rel)
        lines.append(f"{lab} & {rel[1]:+.1f} & {rel[2]:+.1f} & " +
                     " & ".join(f"{v:.2f}" for v in ranks.mean(0)) +
                     (" & -- & -- \\\\" if col == "Linear" else f" & {np.nanmean(taus):.2f} & {100*same_top:.0f} \\\\"))
    tex = r"""\begin{table*}[!t]
\centering
\caption{Robustness of the fresh results to the benchmark model. All feasible final layouts are re-evaluated (not re-optimized) with alternative power-curve and wake models. $\Delta$: mean relative change of the objective (expected power) with respect to the benchmark model; average rank of each method over the 68 cases (methods with at least five feasible runs); $\bar\tau$: mean Kendall rank correlation between the benchmark and alternative orderings of the methods within a case; ``Same best'': percentage of cases in which the best method is unchanged.}
\label{tab:robust-hybrid}
\scriptsize\setlength{\tabcolsep}{3pt}
\begin{tabular}{lccccccccccc}
\toprule
& \multicolumn{2}{c}{$\Delta$ (\%)} & \multicolumn{7}{c}{Average rank} & & \\
\cmidrule(lr){2-3}\cmidrule(lr){4-10}
Model & DS I & DS II & LX-SSA-VNS & LX-SSA & SSA & PSO & DE & VNS & MS-SLSQP & $\bar\tau$ & Same best (\%) \\
\midrule
""" + "\n".join(lines) + r"""
\bottomrule
\end{tabular}
\end{table*}
"""
    open("hybrid_robust_table.tex", "w").write(tex)
    import json
    json.dump(out, open("hybrid_robust_summary.json", "w"), indent=1, default=float)
    print(json.dumps(out, indent=1, default=float))

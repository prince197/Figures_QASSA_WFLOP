"""Robustness re-evaluation for the MPCE (SSA-VNS) study.

Every feasible final layout is re-evaluated (not re-optimized) with
  * the benchmark model (linear power curve, Jensen wake)              -> "Linear"
  * a cubic power curve                                                -> "Cubic"
  * a cubic power curve with a 25 m/s cut-out                          -> "CubicCutout"
  * the Bastankhah--Porte-Agel Gaussian wake (k* = 0.04), linear curve  -> "Gauss"
Same logic as robustness_final.py. Re-evaluations are cached in `mpce_reevaluation_cache.csv`
(key: md5 of dataset + coordinate string), so a rerun only evaluates new layouts; the cache is
seeded from final_reevaluation.csv (layouts from fresh_grid/vgrid/bgrid) when that file exists.
Used by mpce_results.py (section 8); can also be run on its own:  python mpce_robustness.py
"""
import os, hashlib
import numpy as np, pandas as pd
from multiprocessing import Pool
from scipy.stats import rankdata, kendalltau

MODELS = (("Objective", "Benchmark (linear curve, Jensen)"), ("Cubic", "Cubic power curve"),
          ("CubicCutout", "Cubic curve with 25 m/s cut-out"), ("Gauss", "Gaussian wake ($k^*=0.04$)"))
VALS = ["Linear", "Cubic", "CubicCutout", "Gauss", "IdealCubic", "IdealCubicCutout"]


def _key(ds, coords):
    return hashlib.md5(f"{ds}|{coords}".encode()).hexdigest()


def _reval(args):
    ds, coords = args
    from wflop_model import farm_objective
    xy = np.array([[float(v) for v in p.split()] for p in coords.split(";")])
    ds = int(ds)
    lin, _ = farm_objective(xy, ds)
    cub, ic = farm_objective(xy, ds, curve="cubic")
    cco, icc = farm_objective(xy, ds, curve="cubic_cutout")
    gau, _ = farm_objective(xy, ds, wake="gaussian", kstar=0.04)
    return dict(Linear=lin, Cubic=cub, CubicCutout=cco, Gauss=gau, IdealCubic=ic, IdealCubicCutout=icc)


def _seed_cache(here):
    """Keys + values recovered from final_reevaluation.csv (same layouts as the fresh_* files)."""
    fn = os.path.join(here, "final_reevaluation.csv")
    if not os.path.exists(fn):
        return pd.DataFrame(columns=["Key"] + VALS)
    old = pd.read_csv(fn)
    src = pd.concat([pd.read_csv(os.path.join(here, f), usecols=["Algorithm", "Dataset", "Radius", "Turbines",
                                                                "Seed", "Feasible", "Objective", "Coordinates"])
                     for f in ("fresh_grid.csv", "fresh_vgrid.csv", "fresh_bgrid.csv")], ignore_index=True)
    src = src[src.Feasible]
    k = ["Algorithm", "Dataset", "Radius", "Turbines", "Seed"]
    m = old.merge(src, on=k, suffixes=("", "_src"))
    m = m[np.abs(m.Linear - m.Objective_src) <= 1e-6 * np.abs(m.Objective_src).clip(lower=1)]
    m["Key"] = [_key(d, c) for d, c in zip(m.Dataset.astype(str), m.Coordinates)]
    return m[["Key"] + VALS].drop_duplicates("Key")


def reevaluate(G, here, procs=4, log=print):
    """G: feasible benchmark rows (Dataset in {1,2} as str/int, Coordinates). Returns G + model columns."""
    cache_fn = os.path.join(here, "mpce_reevaluation_cache.csv")
    if os.path.exists(cache_fn):
        cache = pd.read_csv(cache_fn)
    else:
        cache = _seed_cache(here)
        log(f"  robustness cache seeded from final_reevaluation.csv: {len(cache)} layouts")
    G = G.reset_index(drop=True).copy()
    G["Key"] = [_key(str(d), c) for d, c in zip(G.Dataset, G.Coordinates)]
    have = set(cache.Key)
    todo = G[~G.Key.isin(have)].drop_duplicates("Key")
    if len(todo):
        log(f"  re-evaluating {len(todo)} new layouts on {procs} processes ...")
        with Pool(procs) as pool:
            res = pool.map(_reval, list(zip(todo.Dataset.astype(str), todo.Coordinates)), chunksize=25)
        new = pd.DataFrame(res); new.insert(0, "Key", todo.Key.values)
        cache = pd.concat([cache, new], ignore_index=True)
        cache.to_csv(cache_fn, index=False)
    log(f"  robustness: {len(G)} layouts ({len(todo)} newly evaluated, cache {len(cache)})")
    return G.merge(cache, on="Key", how="left")


def robust_table(G, methods, lab, rank_fn, n_runs=30):
    """Table rows + summary. rank_fn(fmean, nfeas, nruns) implements the ranking rule of mpce_results."""
    out, lines = {}, []
    keys = ["Dataset", "Radius", "Turbines"]
    cnt = G.groupby(keys + ["Algorithm"]).size().unstack().reindex(columns=methods).fillna(0)
    # Base ordering = recorded objective of the runs (identical to the main case-level table). The stored
    # coordinates are rounded to 1 mm, and optimized layouts often sit on a wake-cone edge, so a
    # re-evaluation of the benchmark model ("Linear") can differ from the recorded objective; relative
    # changes Delta are therefore taken with respect to the re-evaluated benchmark ("Linear"), i.e. the
    # same rounded coordinates under both models.
    base = G.groupby(keys + ["Algorithm"])["Objective"].mean().unstack().reindex(columns=methods)
    dev = (G.Linear - G.Objective).abs() / G.Objective
    out["rounding_check"] = dict(median_rel_dev=float(dev.median()), max_rel_dev=float(dev.max()),
                                 n_rel_dev_above_1e_4=int((dev > 1e-4).sum()), n_layouts=int(len(G)))
    for col, name in MODELS:
        m = G.groupby(keys + ["Algorithm"])[col].mean().unstack().reindex(columns=methods)
        ranks = np.vstack([rank_fn(mv, cv, np.full(len(methods), n_runs)) for mv, cv in zip(m.values, cnt.values)])
        qual = cnt.values >= np.ceil(n_runs / 2)
        taus, same = [], []
        for a, b, q in zip(base.values, m.values, qual):
            ok = q & np.isfinite(a) & np.isfinite(b)
            if ok.sum() >= 2:
                taus.append(kendalltau(a[ok], b[ok])[0])
            if ok.sum() >= 1:
                same.append(np.argmax(np.where(ok, a, -np.inf)) == np.argmax(np.where(ok, b, -np.inf)))
        num = "Linear" if col == "Objective" else col
        rel = {ds: float(100 * (G[G.Dataset.astype(str) == ds][num] / G[G.Dataset.astype(str) == ds]["Linear"] - 1).mean())
               for ds in ("1", "2")}
        out[col] = dict(avg_rank=dict(zip(methods, map(float, ranks.mean(0)))), tau=float(np.nanmean(taus)),
                        same_best_pct=float(100 * np.mean(same)), rel_change_pct=rel)
        lines.append(f"{name} & {rel['1']:+.1f} & {rel['2']:+.1f} & " +
                     " & ".join(f"{v:.2f}" for v in ranks.mean(0)) +
                     (" & -- & -- \\\\" if col == "Objective" else f" & {np.nanmean(taus):.2f} & {100*np.mean(same):.0f} \\\\"))
    k = len(methods)
    tex = r"""\begin{table*}[!t]
\centering
\caption{Robustness of the results to the benchmark model. All feasible final layouts of the seven methods are re-evaluated (not re-optimized) with alternative power-curve and wake models. $\Delta$: mean relative change of the objective (expected power) with respect to the benchmark model (both evaluated on the stored coordinates, rounded to 1~mm); average rank of each method over the 68 cases (ranking rule of Table~\ref{tab:friedman68}: methods with fewer than 15 feasible runs are ranked last, by their number of feasible runs); $\bar\tau$: mean Kendall rank correlation between the benchmark and the alternative ordering of the methods with at least 15 feasible runs within a case; ``Same best'': percentage of cases in which the best method is unchanged.}
\label{tab:robust}
\scriptsize\setlength{\tabcolsep}{3pt}
\begin{tabular}{l""" + "c" * (k + 4) + r"""}
\toprule
& \multicolumn{2}{c}{$\Delta$ (\%)} & \multicolumn{""" + str(k) + r"""}{c}{Average rank} & & \\
\cmidrule(lr){2-3}\cmidrule(lr){4-""" + str(3 + k) + r"""}
Model & DS I & DS II & """ + " & ".join(lab[a] for a in methods) + r""" & $\bar\tau$ & Same best (\%) \\
\midrule
""" + "\n".join(lines) + r"""
\bottomrule
\end{tabular}
\end{table*}
"""
    return tex, out


if __name__ == "__main__":
    import mpce_results
    mpce_results.main(["--only-robust"])

"""Summary of the GA Lillgrund runs (rev4_galg*.csv) next to the nine archived methods (rev2_lg16*.csv).

Usage:  python3 rev4_galg_summary.py      -> prints the table and writes rev4_galg_summary.json
Per budget and method: n, mean / SD (ddof = 1) over all 30 seeds of Objective (AEP, GWh/yr) and WakeLoss
(GWh/yr), the same over the feasible runs only, feasible count, mean Seconds, Calls range. Paired test PSO-VNS
(PSOBV) vs GA: exact two-sided Wilcoxon signed-rank test on WakeLoss over the seeds where both final layouts are
feasible (WakeLoss = Ideal - Objective with a constant Ideal, so the test on Objective is identical), plus the
number of those seeds in which each method has the lower wake loss.
"""
import glob, json
import numpy as np, pandas as pd
from scipy.stats import wilcoxon

FILES = {6030: ("rev2_lg16_s*of4.csv", "rev4_galg_s0of1.csv"), 30030: ("rev2_lg16b_s*of30.csv", "rev4_galgb_s0of1.csv")}
ORDER = ["GA", "PSOBV", "PSOC", "SSABV", "SSA", "LXSSA", "DE", "BVNS", "SLSQP", "RSDVNS"]
LABEL = {"PSOBV": "PSO-VNS", "PSOC": "PSO", "SSABV": "SSA-VNS", "SLSQP": "MS-SLSQP"}


def sd(v):
    return float(np.std(v, ddof=1)) if len(v) > 1 else None


def mean(v):
    return float(np.mean(v)) if len(v) else None


def main():
    out = {"source": {str(b): {"archived": sorted(glob.glob(a)), "GA": g} for b, (a, g) in FILES.items()},
           "units": "Objective = AEP and WakeLoss in GWh/yr; Seconds per run", "budgets": {}}
    for b, (pat, gfile) in FILES.items():
        d = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(pat))] + [pd.read_csv(gfile)], ignore_index=True)
        assert (d.Budget == b).all() and (d.Calls == b).all() and (d.Dataset == "LG").all()
        rows = {}
        print(f"\nLillgrund 4x4 block, budget {b}: mean (SD) over 30 seeds; feasible-only mean (SD)")
        print(f"{'method':9s} {'n':>3s} {'feas':>5s} {'AEP':>18s} {'WakeLoss':>16s} {'AEP|feas':>18s} {'WL|feas':>16s} {'sec':>6s}")
        for a in ORDER:
            g = d[d.Algorithm == a].sort_values("Seed")
            assert list(g.Seed) == list(range(1, 31)), a
            fz = g[g.Feasible]
            r = dict(label=LABEL.get(a, a), n=len(g), feasible=int(g.Feasible.sum()),
                     AEP_mean=mean(g.Objective), AEP_sd=sd(g.Objective),
                     WakeLoss_mean=mean(g.WakeLoss), WakeLoss_sd=sd(g.WakeLoss),
                     AEP_feasible_mean=mean(fz.Objective), AEP_feasible_sd=sd(fz.Objective),
                     WakeLoss_feasible_mean=mean(fz.WakeLoss), WakeLoss_feasible_sd=sd(fz.WakeLoss),
                     Seconds_mean=mean(g.Seconds), Seconds_total=float(g.Seconds.sum()),
                     Calls_min=int(g.Calls.min()), Calls_max=int(g.Calls.max()), Ideal=float(g.Ideal.iloc[0]))
            rows[a] = r
            f4 = lambda m, s: "      -          " if m is None else f"{m:8.3f} ({s if s is not None else float('nan'):6.3f})"
            print(f"{r['label']:9s} {r['n']:3d} {r['feasible']:3d}/30 {f4(r['AEP_mean'], r['AEP_sd']):>18s} "
                  f"{f4(r['WakeLoss_mean'], r['WakeLoss_sd']):>16s} {f4(r['AEP_feasible_mean'], r['AEP_feasible_sd']):>18s} "
                  f"{f4(r['WakeLoss_feasible_mean'], r['WakeLoss_feasible_sd']):>16s} {r['Seconds_mean']:6.1f}")
        p = d[d.Algorithm == "PSOBV"].set_index("Seed").sort_index()
        q = d[d.Algorithm == "GA"].set_index("Seed").sort_index()
        both = p.Feasible & q.Feasible
        x, y = p.WakeLoss[both].to_numpy(), q.WakeLoss[both].to_numpy()
        diff = x - y
        test = dict(seeds_both_feasible=[int(s) for s in p.index[both]], n_pairs=int(both.sum()),
                    PSOBV_better=int((diff < 0).sum()), GA_better=int((diff > 0).sum()), ties=int((diff == 0).sum()),
                    median_diff_WakeLoss_PSOBV_minus_GA=float(np.median(diff)) if len(diff) else None,
                    mean_diff_WakeLoss_PSOBV_minus_GA=float(np.mean(diff)) if len(diff) else None)
        if len(diff) and (diff != 0).any():
            w = wilcoxon(x, y, alternative="two-sided", method="exact")
            test.update(W=float(w.statistic), p_exact_two_sided=float(w.pvalue), method="scipy.stats.wilcoxon exact")
        allp = p.WakeLoss - q.WakeLoss
        test["all30_PSOBV_lower_WakeLoss"] = int((allp < 0).sum())
        test["all30_GA_lower_WakeLoss"] = int((allp > 0).sum())
        test["feasible_only_PSOBV_not_GA"] = int((p.Feasible & ~q.Feasible).sum())
        test["feasible_only_GA_not_PSOBV"] = int((~p.Feasible & q.Feasible).sum())
        print(f"PSO-VNS vs GA (WakeLoss, seeds both feasible, n = {test['n_pairs']}): PSO-VNS better {test['PSOBV_better']}, "
              f"GA better {test['GA_better']}, ties {test['ties']}; exact Wilcoxon p = {test.get('p_exact_two_sided')}")
        out["budgets"][str(b)] = dict(methods=rows, PSOBV_vs_GA=test)
    with open("rev4_galg_summary.json", "w") as fh:
        json.dump(out, fh, indent=1)
    print("\nwrote rev4_galg_summary.json")


if __name__ == "__main__":
    main()

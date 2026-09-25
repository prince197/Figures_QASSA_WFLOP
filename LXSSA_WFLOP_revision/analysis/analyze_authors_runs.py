"""Statistics and LaTeX tables for the controlled re-runs (authors_runs_*.csv).

Ranking rule for paired tests: feasible runs rank ahead of infeasible ones; among feasible runs a
larger objective is better; among infeasible runs a smaller total spacing shortfall is better.
Wilcoxon signed-rank (LX-SSA vs each baseline, paired by seed, Holm within case), matched-pairs
rank-biserial correlation, and a run-level Friedman test with seeds as blocks."""
import numpy as np, pandas as pd
from scipy.stats import wilcoxon, friedmanchisquare, rankdata, mannwhitneyu

ALGS = ["LXSSA", "SSA", "PSO", "DE"]
LAB = {"LXSSA": "LX-SSA", "SSA": "SSA", "PSO": "PSO", "DE": "DE"}


def score(df):
    """Scalar 'goodness' honouring feasibility-first ordering (larger is better)."""
    short = np.maximum(0, 308.0 - df.MinSpacing)
    return np.where(df.Feasible, df.Objective, -1e9 - short)


def holm(p):
    p = np.asarray(p, float); o = np.argsort(p); m = len(p); out = np.empty(m); run = 0
    for i, k in enumerate(o):
        run = max(run, min(1.0, (m - i) * p[k])); out[k] = run
    return out


def paired_tests(sub):
    """sub: runs of one case/budget. Returns Friedman stats and per-baseline Wilcoxon rows."""
    piv = sub.assign(S=score(sub)).pivot(index="Seed", columns="Algorithm", values="S")[ALGS]
    chi, pf = friedmanchisquare(*[piv[a] for a in ALGS])
    ranks = np.vstack([rankdata(-r) for r in piv.values]).mean(0)
    rows = []
    for b in ALGS[1:]:
        d = piv["LXSSA"] - piv[b]
        nz = d[d != 0]
        if len(nz) == 0:
            p, rb = 1.0, 0.0
        else:
            p = wilcoxon(nz).pvalue
            rk = rankdata(np.abs(nz))
            rb = (rk[nz > 0].sum() - rk[nz < 0].sum()) / rk.sum()
        fm = sub[sub.Feasible].groupby("Algorithm").Objective.mean()
        rows.append(dict(Baseline=b, MeanDiff=fm.get("LXSSA", np.nan) - fm.get(b, np.nan),
                         Wins=int((d > 0).sum()), Losses=int((d < 0).sum()), P=p, RB=rb))
    ph = holm([r["P"] for r in rows])
    for r, h in zip(rows, ph):
        r["PHolm"] = h
    return chi, pf, dict(zip(ALGS, ranks)), rows


def fmt_ms(mean, sd, scale=1.0):
    if not np.isfinite(mean):
        return "--"
    if not np.isfinite(sd):
        return f"{mean/scale:.1f}"
    return f"{mean/scale:.1f} $\\pm$ {sd/scale:.1f}"


def fmt_p(p):
    if p < 1e-3:
        m, e = f"{p:.2e}".split("e")
        return f"${m}\\times10^{{{int(e)}}}$"
    return f"{p:.3f}"


def main():
    B = pd.read_csv("authors_runs_budget.csv")
    C = pd.read_csv("authors_runs_crowded.csv")
    S = pd.read_csv("authors_runs_spacing.csv")
    H = pd.read_csv("authors_runs_hornsrev.csv")
    rec = pd.read_csv("../selected_30_run_data.csv")
    tex = []
    summary = []

    # ---------- calibration: native budgets vs recorded ----------
    cal = []
    for (ds, rad, n), sub in B.groupby(["Dataset", "Radius", "Turbines"]):
        for a in ALGS:
            native = 6030 if a == "LXSSA" else 3030
            x = sub[(sub.Algorithm == a) & (sub.Budget == native)].Objective
            y = rec[(rec.Algorithm == a) & (rec.Dataset == ds) & (rec.Radius == rad) & (rec.Turbines == n)].EnergyProduction
            cal.append(dict(Dataset=ds, Radius=rad, Turbines=n, Algorithm=a, Rerun=x.mean(),
                            Recorded=y.mean(), P=mannwhitneyu(x, y).pvalue))
    cal = pd.DataFrame(cal); cal.to_csv("calibration_vs_recorded.csv", index=False)
    summary.append(f"calibration: {int((cal.P >= 0.05).sum())}/{len(cal)} comparisons p>=0.05; "
                   f"min p {cal.P.min():.3g}; max |rel diff| {100*np.max(np.abs(cal.Rerun/cal.Recorded-1)):.3f}%")

    # ---------- equal-budget table + tests ----------
    rows_t, rows_s = [], []
    allstats = []
    for budget in (3030, 6030):
        for (ds, rad, n), sub in B[B.Budget == budget].groupby(["Dataset", "Radius", "Turbines"]):
            chi, pf, ranks, rows = paired_tests(sub)
            m = sub[sub.Feasible].groupby("Algorithm").Objective.agg(["mean", "std"]).reindex(ALGS)
            feas = sub.groupby("Algorithm").Feasible.mean()
            best = m["mean"].idxmax()
            cells = []
            for a in ALGS:
                c = fmt_ms(m.loc[a, 'mean'], m.loc[a, 'std'])
                cells.append(f"\\textbf{{{c}}}" if a == best else c)
                if feas[a] < 1:
                    cells[-1] += f"$^{{{int(round(30*feas[a]))}}}$"
            rows_t.append(f"{budget:,} & {ds} & {rad} & {n} & " + " & ".join(cells) + " \\\\")
            for r in rows:
                allstats.append(dict(Budget=budget, Dataset=ds, Radius=rad, Turbines=n, Chi2=chi, PF=pf, **r))
            rows_s.append((budget, ds, rad, n, chi, pf, ranks, rows))
        rows_t.append("\\midrule")
    rows_t = rows_t[:-1]
    pd.DataFrame(allstats).to_csv("equal_budget_tests.csv", index=False)
    tex.append(r"""\begin{table*}[!t]
\centering
\caption{Controlled equal-budget re-runs with the authors' optimizer code: benchmark objective (mean $\pm$ SD of the feasible runs among 30 seed-paired runs). Budget = objective-function calls per run (LX-SSA 50 or 100 iterations; SSA, PSO and DE 100 or 200 iterations; population 30). Bold: highest mean. A superscript gives the number of feasible runs when fewer than 30; infeasible runs are excluded from the mean.}
\label{tab:equalbudget}
\footnotesize\setlength{\tabcolsep}{4pt}
\begin{tabular}{cccccccc}
\toprule
Budget & DS & Radius (m) & $N$ & LX-SSA & SSA & PSO & DE \\
\midrule
""" + "\n".join(rows_t) + r"""
\bottomrule
\end{tabular}
\end{table*}
""")
    lines = []
    for budget, ds, rad, n, chi, pf, ranks, rows in rows_s:
        first = True
        for r in rows:
            pre = (f"{budget:,} & {ds} & {rad} & {n} & {chi:.2f} & {fmt_p(pf)} & {ranks['LXSSA']:.2f}"
                   if first else " & & & & & & ")
            lines.append(pre + f" & {r['Baseline']} & {r['MeanDiff']:+.1f} & {r['Wins']}/{r['Losses']} & "
                         f"{fmt_p(r['PHolm'])} & {r['RB']:+.2f} \\\\")
            first = False
        lines.append("\\midrule")
    lines = lines[:-1]
    tex.append(r"""\begin{table*}[!t]
\centering
\caption{Paired statistics for the equal-budget re-runs (runs paired by seed, i.e., by identical initial population). $\chi^2_F$, $p_F$: run-level Friedman test over the four algorithms (30 blocks); $\bar r$: average Friedman rank of LX-SSA (1 = best). Infeasible runs rank below feasible ones (and among themselves by spacing shortfall). For each baseline: difference of feasible-run means (LX-SSA minus baseline), LX-SSA wins/losses over the 30 pairs, Holm-adjusted two-sided Wilcoxon signed-rank $p$, and matched-pairs rank-biserial correlation (positive favors LX-SSA).}
\label{tab:equalbudget-tests}
\scriptsize\setlength{\tabcolsep}{3pt}
\begin{tabular}{cccccccccccc}
\toprule
Budget & DS & Radius & $N$ & $\chi^2_F$ & $p_F$ & $\bar r$ & Baseline & $\Delta$ & W/L & $p_{\rm Holm}$ & $r_{rb}$ \\
\midrule
""" + "\n".join(lines) + r"""
\bottomrule
\end{tabular}
\end{table*}
""")

    # ---------- crowded cases ----------
    lines = []
    for (ds, rad, n), sub in C.groupby(["Dataset", "Radius", "Turbines"]):
        cells = []
        for a in ALGS:
            s = sub[sub.Algorithm == a]
            f = s[s.Feasible]
            val = f"{f.Objective.mean():.1f}" if len(f) else "--"
            cells.append(f"{len(f)}/30 & {val}")
        lines.append(f"{ds} & {rad} & {n} & " + " & ".join(cells) + " \\\\")
    tex.append(r"""\begin{table*}[!t]
\centering
\caption{Equal-budget re-runs (6,030 calls) for the largest turbine count tabulated per radius: number of feasible final layouts out of 30 and mean benchmark objective of the feasible layouts.}
\label{tab:crowded}
\footnotesize\setlength{\tabcolsep}{4pt}
\begin{tabular}{ccccccccccc}
\toprule
& & & \multicolumn{2}{c}{LX-SSA} & \multicolumn{2}{c}{SSA} & \multicolumn{2}{c}{PSO} & \multicolumn{2}{c}{DE} \\
\cmidrule(lr){4-5}\cmidrule(lr){6-7}\cmidrule(lr){8-9}\cmidrule(lr){10-11}
DS & Radius (m) & $N$ & Feas. & Objective & Feas. & Objective & Feas. & Objective & Feas. & Objective \\
\midrule
""" + "\n".join(lines) + r"""
\bottomrule
\end{tabular}
\end{table*}
""")

    # ---------- spacing with LX-SSA and SSA ----------
    S4 = B[(B.Budget == 6030) & B.Algorithm.isin(["LXSSA", "SSA"])].assign(Spacing="4D")
    SS = pd.concat([S4, S], ignore_index=True)
    SS.to_csv("authors_spacing_all.csv", index=False)
    lines = []
    for (ds, rad, n), sub in SS.groupby(["Dataset", "Radius", "Turbines"]):
        cells = []
        for sp in ("4D", "5D", "6D"):
            for a in ("LXSSA", "SSA"):
                s = sub[(sub.Spacing == sp) & (sub.Algorithm == a)]
                c = f"{s[s.Feasible].Objective.mean():.1f}"
                if s.Feasible.mean() < 1:
                    c += f"$^{{{int(s.Feasible.sum())}}}$"
                cells.append(c)
            s = sub[sub.Spacing == sp]
            piv = s.assign(Sc=score(s)).pivot(index="Seed", columns="Algorithm", values="Sc")
            d = piv["LXSSA"] - piv["SSA"]; nz = d[d != 0]
            cells.append(fmt_p(wilcoxon(nz).pvalue) if len(nz) else "1.000")
        wl = [sub[(sub.Spacing == sp) & (sub.Algorithm == "LXSSA")] for sp in ("4D", "6D")]
        pct = [100 * w.WakeLoss.mean() / w.Ideal.iloc[0] for w in wl]
        lines.append(f"{ds} & {rad} & {n} & " + " & ".join(cells) + f" & {pct[0]:.2f} $\\rightarrow$ {pct[1]:.2f} \\\\")
    tex.append(r"""\begin{table*}[!t]
\centering
\caption{Minimum-spacing sensitivity with the authors' LX-SSA and SSA code (6,030 calls, 30 seed-paired runs per cell): mean benchmark objective, Wilcoxon signed-rank $p$ for LX-SSA vs.\ SSA, and the mean LX-SSA wake loss as a percentage of the ideal objective at $4D$ and $6D$. A superscript gives the number of feasible runs when fewer than 30.}
\label{tab:spacing-authors}
\footnotesize\setlength{\tabcolsep}{3pt}
\begin{tabular}{ccccccccccccc}
\toprule
& & & \multicolumn{3}{c}{$4D$ (308 m)} & \multicolumn{3}{c}{$5D$ (385 m)} & \multicolumn{3}{c}{$6D$ (462 m)} & \\
\cmidrule(lr){4-6}\cmidrule(lr){7-9}\cmidrule(lr){10-12}
DS & Radius & $N$ & LX-SSA & SSA & $p$ & LX-SSA & SSA & $p$ & LX-SSA & SSA & $p$ & LX-SSA loss (\%) \\
\midrule
""" + "\n".join(lines) + r"""
\bottomrule
\end{tabular}
\end{table*}
""")

    # ---------- Horns Rev 1 ----------
    lines = []
    hstats = []
    for (ds, rad, n), sub in H.groupby(["Dataset", "Radius", "Turbines"]):
        chi, pf, ranks, rows = paired_tests(sub)
        m = sub[sub.Feasible].groupby("Algorithm").Objective.agg(["mean", "std"]).reindex(ALGS)
        feas = sub.groupby("Algorithm").Feasible.sum()
        best = m["mean"].idxmax()
        cells = []
        for a in ALGS:
            c = fmt_ms(m.loc[a, "mean"], m.loc[a, "std"], 15.0)
            c = f"\\textbf{{{c}}}" if a == best else c
            if feas[a] < 30:
                c += f"$^{{{int(feas[a])}}}$"
            cells.append(c)
        ideal_kw = sub.Ideal.iloc[0] / 15
        lines.append(f"{rad} & {n} & {ideal_kw:.1f} & " + " & ".join(cells) + f" & {fmt_p(pf)} & "
                     + " & ".join(fmt_p(r['PHolm']) for r in rows) + " \\\\")
        for r in rows:
            hstats.append(dict(Radius=rad, Turbines=n, Chi2=chi, PF=pf, **r))
    pd.DataFrame(hstats).to_csv("hornsrev_tests.csv", index=False)
    tex.append(r"""\begin{table*}[!t]
\centering
\caption{Measured wind climate (Horns Rev~1 wind rose, 12 sectors): expected farm power in kW (mean $\pm$ SD of the feasible runs among 30 seed-paired runs, 6,030 calls each; a superscript gives the number of feasible runs when fewer than 30; AEP in MWh $=8.76\times$ kW). ``Ideal'' is the wake-free farm power. Bold: highest mean. $p_F$: run-level Friedman test; last three columns: Holm-adjusted Wilcoxon signed-rank $p$ of LX-SSA vs.\ SSA, PSO and DE.}
\label{tab:hornsrev}
\footnotesize\setlength{\tabcolsep}{3pt}
\begin{tabular}{cccccccccccc}
\toprule
& & & & & & & & \multicolumn{3}{c}{LX-SSA vs.\ ($p_{\rm Holm}$)} \\
\cmidrule(lr){9-11}
Radius (m) & $N$ & Ideal & LX-SSA & SSA & PSO & DE & $p_F$ & SSA & PSO & DE \\
\midrule
""" + "\n".join(lines) + r"""
\bottomrule
\end{tabular}
\end{table*}
""")

    # ---------- runtime ----------
    rt = pd.concat([B[B.Budget == 6030], C, H]).groupby(["Algorithm", "Turbines"]).Seconds.mean().unstack()
    rt.to_csv("authors_runtime_6030.csv")
    summary.append("runtime s/run @6030 by N:\n" + rt.round(2).to_string())
    open("authors_tables.tex", "w").write("\n".join(tex))
    print("\n".join(summary))


if __name__ == "__main__":
    main()

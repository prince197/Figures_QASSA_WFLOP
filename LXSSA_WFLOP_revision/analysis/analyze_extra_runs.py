"""Statistics and LaTeX tables for the R3-3 baselines (VNS, multistart SLSQP) and the Horns Rev 1
site case. Runs are seed-paired with authors_runs_budget/crowded.csv."""
import numpy as np, pandas as pd
from scipy.stats import wilcoxon, friedmanchisquare, rankdata
from analyze_authors_runs import holm, fmt_p, fmt_ms

ALL = ["LXSSA", "SSA", "PSO", "DE", "VNS", "SLSQP"]
LAB = {"LXSSA": "LX-SSA", "SSA": "SSA", "PSO": "PSO", "DE": "DE", "VNS": "VNS", "SLSQP": "MS-SLSQP"}


def goodness(df, val="Objective"):
    return np.where(df.Feasible, df[val], -1e12 - np.maximum(0, 1e3 - df.MinSpacing))


def lx_tests(sub, baselines, val="Objective"):
    piv = sub.assign(S=goodness(sub, val)).pivot(index="Seed", columns="Algorithm", values="S")
    out = []
    for b in baselines:
        d = piv["LXSSA"] - piv[b]; nz = d[d != 0]
        p = wilcoxon(nz).pvalue if len(nz) else 1.0
        rk = rankdata(np.abs(nz)) if len(nz) else np.array([1.0])
        rb = (rk[nz > 0].sum() - rk[nz < 0].sum()) / rk.sum() if len(nz) else 0.0
        out.append(dict(Baseline=b, P=p, RB=rb))
    for r, h in zip(out, holm([r["P"] for r in out])):
        r["PHolm"] = h
    return out, piv


def main():
    A = pd.concat([pd.read_csv("authors_runs_budget.csv"), pd.read_csv("authors_runs_crowded.csv")])
    X = pd.read_csv("extra_runs_extra.csv")
    S = pd.concat([A, X], ignore_index=True)
    tex, notes = [], []

    # ---------- representative cases: VNS / SLSQP vs LX-SSA ----------
    lines, ranks = [], {3030: [], 6030: []}
    for budget in (3030, 6030):
        for (ds, rad, n), sub in S[(S.Budget == budget) & (S.Turbines <= 8)].groupby(["Dataset", "Radius", "Turbines"]):
            tests, piv = lx_tests(sub, ["VNS", "SLSQP"])
            ranks[budget].append(np.vstack([rankdata(-r) for r in piv[ALL].values]).mean(0))
            f = sub[sub.Feasible].groupby("Algorithm").Objective.agg(["mean", "std"]).reindex(ALL)
            cnt = sub.groupby("Algorithm").Feasible.sum()
            cells = []
            for a in ("LXSSA", "VNS", "SLSQP"):
                c = fmt_ms(f.loc[a, "mean"], f.loc[a, "std"])
                if cnt[a] < 30:
                    c += f"$^{{{int(cnt[a])}}}$"
                cells.append(c)
            lines.append(f"{budget:,} & {ds} & {rad} & {n} & " + " & ".join(cells) + " & "
                         + " & ".join(f"{fmt_p(t['PHolm'])} ({t['RB']:+.2f})" for t in tests) + " \\\\")
        lines.append("\\midrule")
    lines = lines[:-1]
    avg = {b: dict(zip(ALL, np.mean(ranks[b], 0))) for b in ranks}
    pd.DataFrame(avg).to_csv("six_method_average_ranks.csv")
    notes.append("six-method average ranks: " + str({b: {k: round(v, 2) for k, v in avg[b].items()} for b in avg}))
    rank_line = "; ".join(f"{b:,} calls: " + ", ".join(f"{LAB[a]} {avg[b][a]:.2f}" for a in sorted(ALL, key=lambda a: avg[b][a]))
                          for b in (3030, 6030))
    tex.append(r"""\begin{table*}[!t]
\centering
\caption{Additional baselines requested by Reviewer~3: variable neighbourhood search (VNS) and multistart SLSQP (MS-SLSQP, gradient-based with explicit constraints), seed-paired with LX-SSA at equal numbers of objective calls. Entries: mean $\pm$ SD of the benchmark objective over the feasible runs (superscript: number of feasible runs when fewer than 30). Last two columns: Holm-adjusted Wilcoxon signed-rank $p$ of LX-SSA vs.\ the baseline (matched-pairs rank-biserial correlation in parentheses; positive favors LX-SSA).}
\label{tab:vns-slsqp}
\footnotesize\setlength{\tabcolsep}{3pt}
\begin{tabular}{ccccccccc}
\toprule
Budget & DS & Radius (m) & $N$ & LX-SSA & VNS & MS-SLSQP & vs.\ VNS & vs.\ MS-SLSQP \\
\midrule
""" + "\n".join(lines) + r"""
\bottomrule
\end{tabular}
\end{table*}
""")

    # ---------- crowded: feasibility ----------
    lines = []
    for (ds, rad, n), sub in S[(S.Budget == 6030) & (S.Turbines >= 10)].groupby(["Dataset", "Radius", "Turbines"]):
        cells = []
        for a in ALL:
            s = sub[sub.Algorithm == a]; fz = s[s.Feasible]
            cells.append(f"{len(fz)}" + (f" ({fz.Objective.mean():.0f})" if len(fz) else ""))
        lines.append(f"{ds} & {rad} & {n} & " + " & ".join(cells) + " \\\\")
    tex.append(r"""\begin{table*}[!t]
\centering
\caption{Crowded cases (largest tabulated $N$, 6,030 calls): number of feasible layouts out of 30 for all six methods, with the mean benchmark objective of the feasible layouts in parentheses.}
\label{tab:crowded-all}
\footnotesize\setlength{\tabcolsep}{4pt}
\begin{tabular}{ccccccccc}
\toprule
DS & Radius (m) & $N$ & LX-SSA & SSA & PSO & DE & VNS & MS-SLSQP \\
\midrule
""" + "\n".join(lines) + r"""
\bottomrule
\end{tabular}
\end{table*}
""")

    # ---------- Horns Rev 1 site case ----------
    import hornsrev_model as hr
    H16 = pd.read_csv("extra_runs_hr16.csv"); H80 = pd.read_csv("extra_runs_hr80.csv")
    inst = {n: hr.aep_gwh(hr.site(n)[0]) for n in (16, 80)}
    ideal = {n: H.IdealAEP.iloc[0] for n, H in ((16, H16), (80, H80))}
    t16, _ = lx_tests(H16, ["SSA", "PSO", "DE", "VNS", "SLSQP"], val="AEP")
    t16 = {t["Baseline"]: t for t in t16}
    lines = [f"Installed layout & {inst[16]:.2f} & {100*(1-inst[16]/ideal[16]):.2f} & -- & -- & "
             f"{inst[80]:.1f} & {100*(1-inst[80]/ideal[80]):.2f} & -- \\\\", "\\midrule"]
    for a in ALL:
        c = []
        for n, H in ((16, H16), (80, H80)):
            s = H[H.Algorithm == a]; fz = s[s.Feasible]; ns = len(s)
            if len(fz):
                c += [f"{fz.AEP.mean():.2f} $\\pm$ {fz.AEP.std():.2f}" if n == 16 else f"{fz.AEP.mean():.1f} $\\pm$ {fz.AEP.std():.1f}",
                      f"{100*(1-fz.AEP.mean()/ideal[n]):.2f}", f"{len(fz)}/{ns}"]
            else:
                c += ["--", "--", f"0/{ns}"]
            if n == 16:
                c.append("--" if a == "LXSSA" else fmt_p(t16[a]["PHolm"]))
        lines.append(f"{LAB[a]} & " + " & ".join(c) + " \\\\")
    tex.append(r"""\begin{table*}[!t]
\centering
\caption{Horns Rev~1 site case with the real Vestas V80 power and thrust curves, the measured 12-sector wind climate and the installed farm outline (Jensen, $k=0.04$, minimum spacing $4D=320$~m). AEP in GWh/yr (mean $\pm$ SD over feasible runs; 6,030 calls per run; 30 seeds for the 16-turbine block, 10 for the 80-turbine farm) and mean wake loss relative to the wake-free AEP ("Loss", \%). $p_{\rm Holm}$: Holm-adjusted Wilcoxon signed-rank test of LX-SSA vs.\ each method for the 16-turbine block.}
\label{tab:hornsrev-site}
\footnotesize\setlength{\tabcolsep}{4pt}
\begin{tabular}{lccccccc}
\toprule
& \multicolumn{4}{c}{16-turbine block (wake-free %.2f GWh)} & \multicolumn{3}{c}{80-turbine farm (wake-free %.1f GWh)} \\
\cmidrule(lr){2-5}\cmidrule(lr){6-8}
Layout / method & AEP & Loss & Feas. & $p_{\rm Holm}$ & AEP & Loss & Feas. \\
\midrule
""".replace("%.2f GWh", f"{ideal[16]:.2f} GWh").replace("%.1f GWh", f"{ideal[80]:.1f} GWh") + "\n".join(lines) + r"""
\bottomrule
\end{tabular}
\end{table*}
""")
    open("extra_tables.tex", "w").write("\n".join(tex))
    for n, H in ((16, H16), (80, H80)):
        g = H[H.Feasible].groupby("Algorithm").AEP.agg(["mean", "std", "max", "count"])
        notes.append(f"HR{n} installed {inst[n]:.3f} ideal {ideal[n]:.3f}\n{g.round(3)}")
        notes.append(f"HR{n} runtime s/run: " + str(H.groupby('Algorithm').Seconds.mean().round(1).to_dict()))
    notes.append("HR16 tests: " + str({k: (round(v['PHolm'], 4), round(v['RB'], 2)) for k, v in t16.items()}))
    print("\n".join(notes))


if __name__ == "__main__":
    main()

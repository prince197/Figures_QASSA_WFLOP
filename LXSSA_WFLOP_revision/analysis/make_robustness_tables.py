"""Generate the LaTeX tables of the robustness section from the analysis CSVs."""
import pandas as pd

ALG = ["LXSSA", "SSA", "PSO", "DE"]
LAB = {"LXSSA": "LX", "SSA": "SSA", "PSO": "PSO", "DE": "DE"}
CASES = [(1, 500, 4), (1, 750, 8), (1, 1000, 8), (2, 500, 4), (2, 750, 8), (2, 1000, 8)]
rs = pd.read_csv("robustness_summary.csv").set_index(["Dataset", "Radius", "Turbines", "Algorithm"])
out = []

# power curve table
out.append(r"""\begin{table*}[!t]
\centering
\caption{Re-evaluation of the 720 stored follow-up layouts with a nonlinear (cubic) power curve. $\bar P$ is the mean expected farm power of the LX-SSA layouts in kW (objective/15). $\Delta_{\rm c}$ and $\Delta_{\rm co}$ are the relative changes of the mean expected power (range over the four algorithms) for the cubic curve without and with a 25~m/s cut-out. The last two columns list the algorithms in decreasing order of mean expected power.}
\label{tab:powercurve}
\footnotesize\setlength{\tabcolsep}{4pt}
\begin{tabular}{cccccccc}
\toprule
DS & Radius (m) & $N$ & $\bar P$ linear (kW) & $\bar P$ cubic (kW) & $\Delta_{\rm c}$ (\%) & $\Delta_{\rm co}$ (\%) & Order: linear $\rightarrow$ cubic \\
\midrule""")
for c in CASES:
    sub = rs.loc[c]
    order = lambda col: ", ".join(LAB[a] for a in sub[col].sort_values(ascending=False).index)
    out.append(f"{c[0]} & {c[1]} & {c[2]} & {sub.loc['LXSSA','Linear']/15:.1f} & {sub.loc['LXSSA','Cubic']/15:.1f} & "
               f"[{sub.RelChangeCubicPct.min():.1f}, {sub.RelChangeCubicPct.max():.1f}] & "
               f"[{sub.RelChangeCubicCutoutPct.min():.1f}, {sub.RelChangeCubicCutoutPct.max():.1f}] & "
               f"{order('Linear')} $\\rightarrow$ {order('Cubic')} \\\\")
out.append(r"""\bottomrule
\end{tabular}
\end{table*}
""")

# Gaussian table
out.append(r"""\begin{table*}[!t]
\centering
\caption{Mean wake loss (benchmark-objective units) of the stored follow-up layouts under the Jensen model used for optimization and under the Gaussian model of Bastankhah and Port\'e-Agel ($k^*=0.04$), 30 layouts per algorithm and case. The layouts were optimized with the Jensen model and are only re-evaluated here.}
\label{tab:gaussian}
\footnotesize\setlength{\tabcolsep}{4pt}
\begin{tabular}{ccccccccccc}
\toprule
& & & \multicolumn{4}{c}{Jensen wake loss} & \multicolumn{4}{c}{Gaussian wake loss ($k^*=0.04$)} \\
\cmidrule(lr){4-7}\cmidrule(lr){8-11}
DS & Radius (m) & $N$ & LX-SSA & SSA & PSO & DE & LX-SSA & SSA & PSO & DE \\
\midrule""")
for c in CASES:
    sub = rs.loc[c]
    def cells(col):
        m = sub[col]; best = m.idxmin()
        return " & ".join((r"\textbf{%.1f}" if a == best else "%.1f") % m[a] for a in ALG)
    out.append(f"{c[0]} & {c[1]} & {c[2]} & {cells('JensenLoss')} & {cells('GaussLoss004')} \\\\")
out.append(r"""\bottomrule
\multicolumn{11}{l}{\footnotesize Bold: lowest mean wake loss in the case under the respective model.}
\end{tabular}
\end{table*}
""")

# packing table
pk = pd.read_csv("packing_capacity.csv")
out.append(r"""\begin{table}[!t]
\centering
\caption{Largest turbine count for which a layout satisfying the circular boundary and the minimum spacing was constructed (multi-start packing search; a constructive lower bound on geometric capacity), compared with the largest $N$ in the historical tables.}
\label{tab:capacity}
\footnotesize
\begin{tabular}{ccccc}
\toprule
Radius (m) & $4D$ (308 m) & $5D$ (385 m) & $6D$ (462 m) & Largest $N$ tabulated \\
\midrule""")
tab_max = {500: 10, 750: 12, 1000: 15}
for r in (500, 750, 1000):
    v = pk[pk.Radius == r].set_index("Spacing").MaxNConstructed
    out.append(f"{r} & {v.get('4D','--')} & {v.get('5D','--')} & {v.get('6D','--')} & {tab_max[r]} \\\\")
out.append(r"""\bottomrule
\end{tabular}
\end{table}
""")

# reference SSA spacing table
d = pd.read_csv("ssa_reference_runs.csv"); d = d[~d.Greedy]
out.append(r"""\begin{table*}[!t]
\centering
\caption{Spacing sensitivity with the independent reference SSA (30 runs per cell, 3,030 objective evaluations per run unless stated). Entries are mean $\pm$ SD of the benchmark objective with mean wake loss as a percentage of the ideal objective in parentheses. All 540 final layouts are feasible for their spacing. The last column repeats the $4D$ case with 200 iterations (6,030 evaluations).}
\label{tab:spacing-reopt}
\footnotesize\setlength{\tabcolsep}{4pt}
\begin{tabular}{ccccccc}
\toprule
DS & Radius (m) & $N$ & $4D$ & $5D$ & $6D$ & $4D$, 6,030 evaluations \\
\midrule""")
for c in CASES:
    sub = d[(d.Dataset == c[0]) & (d.Radius == c[1]) & (d.Turbines == c[2])]
    cells = []
    for sp, it in (("4D", 100), ("5D", 100), ("6D", 100), ("4D", 200)):
        s = sub[(sub.Spacing == sp) & (sub.Iterations == it)]
        assert len(s) == 30 and s.Feasible.all()
        ideal = (s.Objective + s.WakeLoss).iloc[0]
        cells.append(f"${s.Objective.mean():.1f}\\pm{s.Objective.std():.1f}$ ({100*s.WakeLoss.mean()/ideal:.2f})")
    out.append(f"{c[0]} & {c[1]} & {c[2]} & " + " & ".join(cells) + r" \\")
out.append(r"""\bottomrule
\end{tabular}
\end{table*}
""")
open("robustness_tables.tex", "w").write("\n".join(out))
print("\n".join(out))

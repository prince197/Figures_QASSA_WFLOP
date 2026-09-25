"""Original basic VNS (BVNS) against each other method: per-case two-sided Wilcoxon signed-rank
tests on the 30 seed-paired runs, Holm-adjusted over the seven comparisons of each case.
Outputs bvns_tables.tex (tab:bvns), bvns_summary.json and bvns_case_tests.csv."""
import json
import numpy as np, pandas as pd
from hybrid_results import ALGS, LAB, DSN, RADII, goodness, holm
from scipy.stats import wilcoxon, rankdata

FOC = "BVNS"
OTH = [a for a in ALGS if a != FOC]
HDR = {"LXVNS": r"\begin{tabular}{@{}c@{}}LX-SSA-\\VNS\end{tabular}",
       "SLSQP": r"\begin{tabular}{@{}c@{}}MS-\\SLSQP\end{tabular}"}
DSS = {1: "I", 2: "II"}


def tests(sub, val="Objective"):
    piv = sub.assign(S=goodness(sub, val)).pivot(index="Seed", columns="Algorithm", values="S")
    rows = []
    for b in OTH:
        d = piv[FOC] - piv[b]; nz = d[np.abs(d) > 1e-9]
        if len(nz) == 0:
            p, rb = 1.0, 0.0
        else:
            p = wilcoxon(nz).pvalue; rk = rankdata(np.abs(nz))
            rb = (rk[nz > 0].sum() - rk[nz < 0].sum()) / rk.sum()
        rows.append(dict(Baseline=b, P=p, RB=rb))
    for r, h in zip(rows, holm([r["P"] for r in rows])):
        r["PHolm"] = h
    return rows


cols = lambda c: c not in ("Coordinates", "Curve")
G = pd.concat([pd.read_csv(f, usecols=cols) for f in ("fresh_grid.csv", "fresh_hgrid.csv", "fresh_vgrid.csv")],
              ignore_index=True)
rows = []
for (ds, r, n), sub in G.groupby(["Dataset", "Radius", "Turbines"]):
    for x in tests(sub):
        rows.append(dict(Dataset=ds, Radius=r, Turbines=n, **x))
C = pd.DataFrame(rows)
C["Outcome"] = np.where(C.PHolm < 0.05, np.where(C.RB > 0, "W", "L"), "T")
C.to_csv("bvns_case_tests.csv", index=False)

lines = []
for ds in (1, 2):
    for r in RADII:
        s = C[(C.Dataset == ds) & (C.Radius == r)]
        cells = []
        for b in OTH:
            x = s[s.Baseline == b].Outcome.value_counts()
            cells.append(f"{x.get('W',0)}/{x.get('T',0)}/{x.get('L',0)}")
        lines.append(f"{DSS[ds]} & {r} & {s.Turbines.nunique()} & " + " & ".join(cells) + " \\\\")
tot = []
for b in OTH:
    x = C[C.Baseline == b].Outcome.value_counts()
    tot.append(f"{x.get('W',0)}/{x.get('T',0)}/{x.get('L',0)}")
lines.append("\\midrule\nAll & & 68 & " + " & ".join(tot) + " \\\\")
tex = r"""\begin{table}[!t]
\centering
\caption{Pairwise outcome of the original basic VNS (BVNS) against each method: number of cases in which BVNS is significantly better / not significantly different / significantly worse (two-sided Wilcoxon signed-rank test on 30 seed-paired runs, Holm-adjusted over the seven comparisons of each case, $\alpha=0.05$).}
\label{tab:bvns}
\scriptsize\setlength{\tabcolsep}{1.6pt}
\begin{tabular}{llcccccccc}
\toprule
DS & $r$ (m) & Cases & """ + " & ".join(HDR.get(b, LAB[b]) for b in OTH) + r""" \\
\midrule
""" + "\n".join(lines) + r"""
\bottomrule
\end{tabular}
\end{table}
"""
open("bvns_tables.tex", "w").write(tex)

summ = {"wtl": {b: C[C.Baseline == b].Outcome.value_counts().to_dict() for b in OTH}}
for b in ("VNS", "LXVNS"):
    d = C[C.Baseline == b]
    summ[f"vs_{b}"] = {o: d[d.Outcome == o][["Dataset", "Radius", "Turbines"]].values.tolist() for o in "WL"}
    summ[f"vs_{b}_rb_median"] = {o: float(d[d.Outcome == o].RB.median()) for o in "WL" if (d.Outcome == o).any()}
# Horns Rev 16: BVNS vs VNS and vs hybrid
H = pd.concat([pd.read_csv(f"fresh_hr16.csv", usecols=cols), pd.read_csv("fresh_hhr16.csv", usecols=cols),
               pd.read_csv("fresh_vhr16.csv", usecols=cols)], ignore_index=True)
summ["hr16"] = {t["Baseline"]: [float(t["PHolm"]), float(t["RB"])]
                for t in tests(H.rename(columns={"AEP": "Objective"}))}
# convergence: calls until the median run is within 1% of its final wake loss
F = pd.concat([pd.read_csv(f, usecols=["Algorithm", "Feasible", "Calls"]) for f in ("fresh_vgrid.csv",)])
json.dump(summ, open("bvns_summary.json", "w"), indent=1)
print(json.dumps(summ, indent=1))

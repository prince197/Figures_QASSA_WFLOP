"""Tables, figures and statistics for the fresh numerical study (full_grid_experiments.py).

Outputs go to ../figures_fresh/ (PDF figures) and fresh_tables.tex / fresh_stats.txt here.
"""
import os, json
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker
from scipy.stats import wilcoxon, friedmanchisquare, rankdata, norm

ALGS = ["LXSSA", "SSA", "PSO", "DE", "VNS", "SLSQP"]
LAB = {"LXSSA": "LX-SSA", "SSA": "SSA", "PSO": "PSO", "DE": "DE", "VNS": "VNS", "SLSQP": "MS-SLSQP"}
# categorical palette (validated: scripts/validate_palette.js, light mode), fixed order
COL = {"LXSSA": "#2a78d6", "SSA": "#eb6834", "PSO": "#1baf7a", "DE": "#eda100",
       "VNS": "#e87ba4", "SLSQP": "#008300"}
MRK = {"LXSSA": "o", "SSA": "s", "PSO": "^", "DE": "v", "VNS": "D", "SLSQP": "P"}
LS = {"LXSSA": "-", "SSA": "--", "PSO": "-.", "DE": ":", "VNS": (0, (5, 1)), "SLSQP": (0, (3, 1, 1, 1))}
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#e4e3df"
RADII = {500: 10, 750: 12, 1000: 15}
DSN = {1: "Data Set I", 2: "Data Set II"}
FIG = "../figures_fresh"
CHECK = 30

plt.rcParams.update({"font.family": "serif", "font.size": 8, "axes.edgecolor": MUTED,
                     "axes.labelcolor": INK, "xtick.color": MUTED, "ytick.color": MUTED,
                     "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.5,
                     "axes.spines.top": False, "axes.spines.right": False,
                     "legend.frameon": False, "lines.linewidth": 1.4, "pdf.fonttype": 42})


# ------------------------------------------------------------------ helpers
def holm(p):
    p = np.asarray(p, float); o = np.argsort(p); m = len(p); out = np.empty(m); run = 0
    for i, k in enumerate(o):
        run = max(run, min(1.0, (m - i) * p[k])); out[k] = run
    return out


def fmt_p(p):
    if p < 1e-3:
        m, e = f"{p:.2e}".split("e")
        return f"${m}\\times10^{{{int(e)}}}$"
    return f"{p:.3f}"


def goodness(df, val):
    return np.where(df.Feasible, df[val], -1e12 - np.maximum(0, 1e4 - df.MinSpacing))


def paired(sub, val="Objective"):
    """Seed-paired statistics for one case: Friedman + Holm-Wilcoxon LX-SSA vs each baseline."""
    piv = sub.assign(S=goodness(sub, val)).pivot(index="Seed", columns="Algorithm", values="S")[ALGS]
    if np.all(piv.values == piv.values[:, :1]):
        chi, pf = 0.0, 1.0
    else:
        chi, pf = friedmanchisquare(*[piv[a] for a in ALGS])
    ranks = np.vstack([rankdata(-r) for r in piv.values]).mean(0)
    rows = []
    for b in ALGS[1:]:
        d = piv["LXSSA"] - piv[b]; nz = d[np.abs(d) > 1e-9]
        if len(nz) == 0:
            p, rb = 1.0, 0.0
        else:
            p = wilcoxon(nz).pvalue
            rk = rankdata(np.abs(nz))
            rb = (rk[nz > 0].sum() - rk[nz < 0].sum()) / rk.sum()
        rows.append(dict(Baseline=b, P=p, RB=rb))
    for r, h in zip(rows, holm([r["P"] for r in rows])):
        r["PHolm"] = h
    return chi, pf, dict(zip(ALGS, ranks)), rows


def curves(sub):
    """Matrix (runs x checkpoints) of best-feasible objective (NaN before first feasible)."""
    return np.array([[float(v) for v in c.split(";")] for c in sub.Curve])


def save(fig, name):
    os.makedirs(FIG, exist_ok=True)
    fig.savefig(f"{FIG}/{name}.pdf", bbox_inches="tight")
    fig.savefig(f"{FIG}/{name}.png", dpi=160, bbox_inches="tight")
    plt.close(fig)


def legend_row(fig, y=1.02):
    h = [plt.Line2D([], [], color=COL[a], ls=LS[a], marker=MRK[a], ms=4, label=LAB[a]) for a in ALGS]
    fig.legend(handles=h, loc="lower center", ncol=6, bbox_to_anchor=(0.5, y), fontsize=8)


# ------------------------------------------------------------------ main
def main():
    G = pd.read_csv("fresh_grid.csv")
    G["LossPct"] = 100 * G.WakeLoss / G.Ideal
    tex, stats = [], {}

    # ---------- per-case statistics ----------
    caserows = []
    for (ds, r, n), sub in G.groupby(["Dataset", "Radius", "Turbines"]):
        chi, pf, ranks, rows = paired(sub)
        feas = sub.groupby("Algorithm").Feasible.sum()
        fm = sub[sub.Feasible].groupby("Algorithm").Objective.mean()
        for x in rows:
            caserows.append(dict(Dataset=ds, Radius=r, Turbines=n, Chi2=chi, PF=pf, **x,
                                 MeanDiff=fm.get("LXSSA", np.nan) - fm.get(x["Baseline"], np.nan)))
        stats[(ds, r, n)] = dict(ranks=ranks, feas=feas.to_dict(), fm=fm.to_dict())
    C = pd.DataFrame(caserows)
    C["Outcome"] = np.where(C.PHolm < 0.05, np.where(C.RB > 0, "W", "L"), "T")
    C.to_csv("fresh_case_tests.csv", index=False)

    # ---------- Table: detailed results per dataset x radius ----------
    for ds in (1, 2):
        for r, nmax in RADII.items():
            sub = G[(G.Dataset == ds) & (G.Radius == r)]
            lines = []
            for n in range(2, nmax + 1):
                s = sub[sub.Turbines == n]
                ideal = s.Ideal.iloc[0]
                f = s[s.Feasible].groupby("Algorithm").Objective.agg(["mean", "std"]).reindex(ALGS)
                cnt = s.groupby("Algorithm").Feasible.sum().reindex(ALGS)
                best = f["mean"].idxmax()
                cells = []
                for a in ALGS:
                    if not np.isfinite(f.loc[a, "mean"]):
                        cells.append("--$^{0}$"); continue
                    sd = f.loc[a, "std"]
                    c = f"{f.loc[a,'mean']:.1f} ({0 if not np.isfinite(sd) else sd:.1f})"
                    if a == best:
                        c = f"\\textbf{{{c}}}"
                    if cnt[a] < 30:
                        c += f"$^{{{int(cnt[a])}}}$"
                    cells.append(c)
                lines.append(f"{n} & {ideal:.1f} & " + " & ".join(cells) + " \\\\")
            tex.append(r"""\begin{table*}[!t]
\centering
\caption{%s, %d-m farm: benchmark objective, mean (SD) over the feasible runs of 30 seed-paired runs at 6,030 objective calls. Bold: highest mean. A superscript gives the number of feasible runs when fewer than 30; ``--'' means no feasible run.}
\label{tab:res-%d-%d}
\scriptsize\setlength{\tabcolsep}{3pt}
\begin{tabular}{cccccccc}
\toprule
$N$ & Ideal & LX-SSA & SSA & PSO & DE & VNS & MS-SLSQP \\
\midrule
""" % (DSN[ds], r, ds, r) + "\n".join(lines) + r"""
\bottomrule
\end{tabular}
\end{table*}
""")

    # ---------- Table: W/T/L summary ----------
    lines = []
    for ds in (1, 2):
        for r in RADII:
            s = C[(C.Dataset == ds) & (C.Radius == r)]
            cells = []
            for b in ALGS[1:]:
                x = s[s.Baseline == b].Outcome.value_counts()
                cells.append(f"{x.get('W',0)}/{x.get('T',0)}/{x.get('L',0)}")
            lines.append(f"{DSN[ds]} & {r} & {s.Turbines.nunique()} & " + " & ".join(cells) + " \\\\")
    tot = []
    for b in ALGS[1:]:
        x = C[C.Baseline == b].Outcome.value_counts()
        tot.append(f"{x.get('W',0)}/{x.get('T',0)}/{x.get('L',0)}")
    lines.append("\\midrule\nAll & & 68 & " + " & ".join(tot) + " \\\\")
    tex.append(r"""\begin{table}[!t]
\centering
\caption{Pairwise outcome of LX-SSA against each method over all cases: number of cases in which LX-SSA is significantly better / not significantly different / significantly worse (two-sided Wilcoxon signed-rank test on 30 seed-paired runs, Holm-adjusted over the five comparisons of each case, $\alpha=0.05$).}
\label{tab:wtl}
\scriptsize\setlength{\tabcolsep}{2.5pt}
\begin{tabular}{llcccccc}
\toprule
Data set & $r$ (m) & Cases & SSA & PSO & DE & VNS & MS-SLSQP \\
\midrule
""" + "\n".join(lines) + r"""
\bottomrule
\end{tabular}
\end{table}
""")

    # ---------- case-level Friedman over all 68 cases ----------
    R = np.array([[stats[k]["ranks"][a] for a in ALGS] for k in sorted(stats)])
    means = pd.DataFrame(
        [[stats[k]["fm"].get(a, -np.inf) for a in ALGS] for k in sorted(stats)], columns=ALGS)
    cr = np.vstack([rankdata(-row) for row in means.values])      # ranks of case means
    chi, pf = friedmanchisquare(*[cr[:, j] for j in range(len(ALGS))])
    n_c, k = cr.shape
    ff = (n_c - 1) * chi / (n_c * (k - 1) - chi)
    avg = cr.mean(0)
    se = np.sqrt(k * (k + 1) / (6 * n_c))
    best_i = int(np.argmin(avg))
    zp = [(ALGS[j], (avg[j] - avg[best_i]) / se) for j in range(k) if j != best_i]
    ph = holm([2 * norm.sf(abs(z)) for _, z in zp])
    lx_i = ALGS.index("LXSSA")
    zl = [(ALGS[j], (avg[j] - avg[lx_i]) / se) for j in range(k) if j != lx_i]
    phl = holm([2 * norm.sf(abs(z)) for _, z in zl])
    rows = []
    for j, a in enumerate(ALGS):
        vb = "--" if j == best_i else fmt_p(ph[[x[0] for x in zp].index(a)])
        vl = "--" if j == lx_i else fmt_p(phl[[x[0] for x in zl].index(a)])
        wins = int((cr[:, j] == 1).sum())
        feas = G[G.Algorithm == a].Feasible.mean() * 100
        rows.append(f"{LAB[a]} & {avg[j]:.2f} & {wins} & {feas:.1f} & {vb} & {vl} \\\\")
    order = np.argsort(avg)
    tex.append(r"""\begin{table}[!t]
\centering
\caption{Case-level Friedman analysis over all 68 benchmark cases (blocks: cases; ranks of the mean feasible objective, 1 = best). Friedman $\chi^2_F=%.1f$ (5 d.f.), $p=%s$; Iman--Davenport $F_F=%.1f$. Holm-adjusted $p$ values of the average-rank $z$ test against the best-ranked method (%s) and against LX-SSA; ``Best'': cases with the highest mean; ``Feas.'': percentage of feasible runs.}
\label{tab:friedman68}
\scriptsize\setlength{\tabcolsep}{3pt}
\begin{tabular}{lccccc}
\toprule
Method & Avg.\ rank & Best & Feas.\ (\%%) & $p$ vs.\ %s & $p$ vs.\ LX-SSA \\
\midrule
""" % (chi, fmt_p(pf).strip("$") if pf >= 1e-3 else fmt_p(pf).strip("$"), ff, LAB[ALGS[best_i]], LAB[ALGS[best_i]]) +
        "\n".join(rows[i] for i in order) + r"""
\bottomrule
\end{tabular}
\end{table}
""")
    summary = dict(friedman_chi2=chi, friedman_p=pf, iman_davenport=ff,
                   avg_rank={a: float(avg[j]) for j, a in enumerate(ALGS)},
                   best_count={a: int((cr[:, j] == 1).sum()) for j, a in enumerate(ALGS)},
                   feasible_pct={a: float(G[G.Algorithm == a].Feasible.mean() * 100) for a in ALGS},
                   wtl={b: C[C.Baseline == b].Outcome.value_counts().to_dict() for b in ALGS[1:]})

    # ---------- Figure: average ranks ----------
    fig, ax = plt.subplots(figsize=(3.4, 1.9))
    for pos, j in enumerate(order[::-1]):
        a = ALGS[j]
        ax.plot([1, avg[j]], [pos, pos], color=GRID, lw=2, zorder=1)
        ax.scatter(avg[j], pos, s=36, color=COL[a], marker=MRK[a], zorder=3, edgecolor="white", lw=1)
        ax.text(avg[j] + 0.08, pos, f"{avg[j]:.2f}", va="center", fontsize=7, color=INK)
    ax.set_yticks(range(k)); ax.set_yticklabels([LAB[ALGS[j]] for j in order[::-1]])
    ax.set_xlim(1, 6.3); ax.set_xlabel("Average rank over 68 cases (1 = best)")
    ax.grid(axis="y", visible=False)
    save(fig, "avg_ranks")

    # ---------- Figure: wake loss vs N ----------
    fig, axes = plt.subplots(2, 3, figsize=(7.1, 4.0), sharey="row")
    for i, ds in enumerate((1, 2)):
        for j, (r, nmax) in enumerate(RADII.items()):
            ax = axes[i, j]
            sub = G[(G.Dataset == ds) & (G.Radius == r) & G.Feasible]
            for a in ALGS:
                m = sub[sub.Algorithm == a].groupby("Turbines").LossPct.mean()
                ax.plot(m.index, m.values, color=COL[a], ls=LS[a], marker=MRK[a], ms=3.5, label=LAB[a])
            ax.set_title(f"{DSN[ds]}, $r$ = {r} m", fontsize=8, color=INK)
            ax.set_xticks(range(2, nmax + 1, 2 if nmax > 10 else 1))
            if i == 1: ax.set_xlabel("Number of turbines $N$")
            if j == 0: ax.set_ylabel("Mean wake loss (% of ideal)")
    legend_row(fig, 1.0)
    fig.tight_layout()
    save(fig, "wakeloss_vs_n")

    # ---------- Figure: feasibility vs N ----------
    fig, axes = plt.subplots(2, 3, figsize=(7.1, 3.6), sharey=True)
    for i, ds in enumerate((1, 2)):
        for j, (r, nmax) in enumerate(RADII.items()):
            ax = axes[i, j]
            sub = G[(G.Dataset == ds) & (G.Radius == r)]
            for a in ALGS:
                m = sub[sub.Algorithm == a].groupby("Turbines").Feasible.mean() * 100
                ax.plot(m.index, m.values, color=COL[a], ls=LS[a], marker=MRK[a], ms=3.5)
            ax.set_ylim(-5, 105)
            ax.set_title(f"{DSN[ds]}, $r$ = {r} m", fontsize=8, color=INK)
            ax.set_xticks(range(2, nmax + 1, 2 if nmax > 10 else 1))
            if i == 1: ax.set_xlabel("Number of turbines $N$")
            if j == 0: ax.set_ylabel("Feasible runs (%)")
    legend_row(fig, 1.0)
    fig.tight_layout()
    save(fig, "feasibility_vs_n")

    # ---------- Figure: convergence (two N levels) ----------
    levels = {"mid": {500: 6, 750: 8, 1000: 10}, "max": RADII}
    for tag, lv in levels.items():
        fig, axes = plt.subplots(2, 3, figsize=(7.1, 4.0))
        for i, ds in enumerate((1, 2)):
            for j, r in enumerate(RADII):
                ax = axes[i, j]; n = lv[r]
                sub = G[(G.Dataset == ds) & (G.Radius == r) & (G.Turbines == n)]
                ideal = sub.Ideal.iloc[0]
                x = np.arange(1, 202) * CHECK
                for a in ALGS:
                    M = curves(sub[sub.Algorithm == a])
                    L = 100 * (ideal - M) / ideal                       # wake loss %, NaN if none
                    frac = np.mean(np.isfinite(L), 0)
                    med = np.where(frac >= 0.5, np.nanmedian(np.where(np.isfinite(L), L, np.inf), 0), np.nan)
                    q1 = np.where(frac >= 0.5, np.nanpercentile(np.where(np.isfinite(L), L, np.inf), 25, 0), np.nan)
                    q3 = np.where(frac >= 0.75, np.nanpercentile(np.where(np.isfinite(L), L, np.inf), 75, 0), np.nan)
                    ax.plot(x, med, color=COL[a], ls=LS[a], lw=1.2)
                    ax.fill_between(x, q1, q3, color=COL[a], alpha=0.12, lw=0)
                    k_last = np.where(np.isfinite(med))[0]
                    if len(k_last):
                        ax.plot(x[k_last[-1]], med[k_last[-1]], marker=MRK[a], color=COL[a], ms=4)
                ax.set_yscale("log")
                fmt = matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:g}")
                ax.yaxis.set_major_formatter(fmt); ax.yaxis.set_minor_formatter(fmt)
                ax.tick_params(axis="y", which="minor", labelsize=6)
                ax.set_title(f"{DSN[ds]}, $r$ = {r} m, $N$ = {n}", fontsize=8, color=INK)
                if i == 1: ax.set_xlabel("Objective-function calls")
                if j == 0: ax.set_ylabel("Best wake loss (% of ideal)")
        legend_row(fig, 1.0)
        fig.tight_layout()
        save(fig, f"convergence_{tag}")

    # ---------- Figure: box plots (moderate N) ----------
    fig, axes = plt.subplots(2, 3, figsize=(7.1, 3.8))
    for i, ds in enumerate((1, 2)):
        for j, r in enumerate(RADII):
            ax = axes[i, j]; n = levels["mid"][r]
            sub = G[(G.Dataset == ds) & (G.Radius == r) & (G.Turbines == n) & G.Feasible]
            data = [sub[sub.Algorithm == a].LossPct.values for a in ALGS]
            bp = ax.boxplot(data, widths=0.6, patch_artist=True, showfliers=True,
                            flierprops=dict(marker=".", ms=3, mec=MUTED),
                            medianprops=dict(color=INK, lw=1), whiskerprops=dict(color=MUTED),
                            capprops=dict(color=MUTED))
            for patch, a in zip(bp["boxes"], ALGS):
                patch.set_facecolor(COL[a]); patch.set_alpha(0.55); patch.set_edgecolor(COL[a])
            ax.set_xticks(range(1, 7)); ax.set_xticklabels([LAB[a] for a in ALGS], rotation=35, fontsize=6.5)
            ax.set_title(f"{DSN[ds]}, $r$ = {r} m, $N$ = {n}", fontsize=8, color=INK)
            if j == 0: ax.set_ylabel("Wake loss (% of ideal)")
            ax.grid(axis="x", visible=False)
    fig.tight_layout()
    save(fig, "boxplots_mid")

    # ---------- Figure: best layouts (largest N) ----------
    import wflop_model as wm
    fig, axes = plt.subplots(2, 3, figsize=(7.1, 4.8))
    best_rows = []
    for i, ds in enumerate((1, 2)):
        for j, (r, n) in enumerate(RADII.items()):
            ax = axes[i, j]
            sub = G[(G.Dataset == ds) & (G.Radius == r) & (G.Turbines == n) & G.Feasible]
            lx = sub[sub.Algorithm == "LXSSA"].sort_values("Objective").iloc[-1]
            top = sub.sort_values("Objective").iloc[-1]
            for row, mk, colr, lab in ((top, "s", COL[top.Algorithm], f"best overall ({LAB[top.Algorithm]})"),
                                       (lx, "o", COL["LXSSA"], "best LX-SSA")):
                xy = np.array([[float(v) for v in p.split()] for p in row.Coordinates.split(";")])
                ax.scatter(xy[:, 0], xy[:, 1], marker=mk, s=22 if mk == "s" else 12,
                           facecolor="none" if mk == "s" else colr, edgecolor=colr, lw=1, label=lab, zorder=3)
            for row in (lx, top):
                best_rows.append(dict(Dataset=ds, Radius=r, Turbines=n, Algorithm=row.Algorithm, Seed=row.Seed,
                                      Objective=row.Objective, WakeLoss=row.WakeLoss, Coordinates=row.Coordinates))
            t = np.linspace(0, 2 * np.pi, 200)
            ax.plot(r * np.cos(t), r * np.sin(t), color=MUTED, lw=0.8)
            ax.set_aspect("equal"); ax.set_xlim(-1.08 * r, 1.08 * r); ax.set_ylim(-1.08 * r, 1.08 * r)
            ax.set_title(f"{DSN[ds]}, $r$ = {r} m, $N$ = {n}", fontsize=8, color=INK)
            ax.tick_params(labelsize=6)
            ax.legend(fontsize=6, loc="lower left", bbox_to_anchor=(0, -0.02), handletextpad=0.2)
    fig.tight_layout()
    save(fig, "layouts_max")
    pd.DataFrame(best_rows).drop_duplicates().to_csv("fresh_best_layouts_maxN.csv", index=False)

    # ---------- runtime ----------
    rt = (G.Seconds / G.Calls * 1000).groupby(G.Algorithm).mean().reindex(ALGS)
    rn = G.groupby(["Algorithm", "Turbines"]).Seconds.mean().unstack().reindex(ALGS)
    summary["ms_per_call"] = rt.round(3).to_dict()
    summary["sec_per_run_range"] = {a: [float(rn.loc[a].min()), float(rn.loc[a].max())] for a in ALGS}

    # ---------- Horns Rev ----------
    import hornsrev_model as hr
    H = {n: pd.read_csv(f"fresh_hr{n}.csv") for n in (16, 80)}
    inst = {n: hr.aep_gwh(hr.site(n)[0]) for n in (16, 80)}
    ideal = {n: H[n].IdealAEP.iloc[0] for n in (16, 80)}
    _, pf16, rk16, t16 = paired(H[16].rename(columns={"AEP": "Objective"}), "Objective")
    t16 = {t["Baseline"]: t for t in t16}
    lines = [f"Installed layout & {inst[16]:.2f} & {100*(1-inst[16]/ideal[16]):.2f} & -- & -- & "
             f"{inst[80]:.1f} & {100*(1-inst[80]/ideal[80]):.2f} & -- \\\\", "\\midrule"]
    for a in ALGS:
        c = []
        for n in (16, 80):
            s = H[n][H[n].Algorithm == a]; fz = s[s.Feasible]
            if len(fz):
                sd = fz.AEP.std() if len(fz) > 1 else 0.0
                c += [f"{fz.AEP.mean():.2f} ({sd:.2f})" if n == 16 else f"{fz.AEP.mean():.1f} ({sd:.1f})",
                      f"{100*(1-fz.AEP.mean()/ideal[n]):.2f}", f"{len(fz)}/{len(s)}"]
            else:
                c += ["--", "--", f"0/{len(s)}"]
            if n == 16:
                c.append("--" if a == "LXSSA" else fmt_p(t16[a]["PHolm"]))
        lines.append(f"{LAB[a]} & " + " & ".join(c) + " \\\\")
    tex.append((r"""\begin{table*}[!t]
\centering
\caption{Horns Rev~1 site case (real Vestas V80 power and thrust curves, measured 12-sector wind climate, installed farm outline; Jensen wake with $k=0.04$; minimum spacing $4D=320$~m): AEP in GWh/yr, mean (SD) over the feasible runs at 6,030 objective calls (30 seeds for the 16-turbine block, 10 for the 80-turbine farm), mean wake loss relative to the wake-free AEP (``Loss'', \%%), feasible runs, and Holm-adjusted Wilcoxon signed-rank $p$ of LX-SSA vs.\ each method (16-turbine block; run-level Friedman $p=%s$).}
\label{tab:hr-site}
\footnotesize\setlength{\tabcolsep}{4pt}
\begin{tabular}{lccccccc}
\toprule
& \multicolumn{4}{c}{16-turbine block (wake-free %.2f GWh/yr)} & \multicolumn{3}{c}{80-turbine farm (wake-free %.1f GWh/yr)} \\
\cmidrule(lr){2-5}\cmidrule(lr){6-8}
Layout / method & AEP & Loss & Feas. & $p_{\rm Holm}$ & AEP & Loss & Feas. \\
\midrule
""" % (fmt_p(pf16).strip("$"), ideal[16], ideal[80])) + "\n".join(lines) + r"""
\bottomrule
\end{tabular}
\end{table*}
""")
    summary["hr"] = {n: dict(installed=inst[n], ideal=ideal[n],
                             mean={a: float(H[n][(H[n].Algorithm == a) & H[n].Feasible].AEP.mean()) for a in ALGS},
                             feasible={a: int(H[n][H[n].Algorithm == a].Feasible.sum()) for a in ALGS})
                     for n in (16, 80)}
    summary["hr16_tests"] = {k: (float(v["PHolm"]), float(v["RB"])) for k, v in t16.items()}

    # Horns Rev convergence
    fig, axes = plt.subplots(1, 2, figsize=(7.1, 2.4))
    for ax, n in zip(axes, (16, 80)):
        x = np.arange(1, 202) * CHECK
        for a in ALGS:
            M = curves(H[n][H[n].Algorithm == a])
            frac = np.mean(np.isfinite(M), 0)
            med = np.where(frac >= 0.5, np.nanmedian(np.where(np.isfinite(M), M, -np.inf), 0), np.nan)
            ax.plot(x, med, color=COL[a], ls=LS[a], lw=1.2, label=LAB[a])
            k_last = np.where(np.isfinite(med))[0]
            if len(k_last):
                ax.plot(x[k_last[-1]], med[k_last[-1]], marker=MRK[a], color=COL[a], ms=4)
        ax.axhline(inst[n], color=INK, lw=0.9, ls=(0, (1, 1)))
        ax.text(x[5], inst[n], "installed layout", va="bottom", fontsize=6.5, color=INK)
        ax.set_title(f"Horns Rev 1, {n} turbines", fontsize=8, color=INK)
        ax.set_xlabel("Objective-function calls"); ax.set_ylabel("Best feasible AEP (GWh/yr)")
    legend_row(fig, 1.0)
    fig.tight_layout()
    save(fig, "hr_convergence")

    # Horns Rev layouts
    fig, axes = plt.subplots(1, 2, figsize=(7.1, 2.8))
    for ax, n in zip(axes, (16, 80)):
        xy0, poly = hr.site(n)
        pp = np.vstack([poly, poly[:1]])
        ax.plot(pp[:, 0], pp[:, 1], color=MUTED, lw=0.8)
        ax.scatter(xy0[:, 0], xy0[:, 1], marker="x", s=14, color=INK, lw=0.8, label="installed", zorder=3)
        fz = H[n][H[n].Feasible]
        if len(fz):
            top = fz.sort_values("AEP").iloc[-1]
            xy = np.array([[float(v) for v in p.split()] for p in top.Coordinates.split(";")])
            ax.scatter(xy[:, 0], xy[:, 1], marker=MRK[top.Algorithm], s=16, facecolor="none",
                       edgecolor=COL[top.Algorithm], lw=1, label=f"best optimized ({LAB[top.Algorithm]})", zorder=3)
        ax.set_aspect("equal"); ax.tick_params(labelsize=6)
        ax.set_title(f"Horns Rev 1, {n} turbines (m)", fontsize=8, color=INK)
        ax.legend(fontsize=6, loc="upper left", bbox_to_anchor=(0, -0.08), ncol=2)
    fig.tight_layout()
    save(fig, "hr_layouts")

    open("fresh_tables.tex", "w").write("\n".join(tex))
    json.dump(summary, open("fresh_summary.json", "w"), indent=1, default=float)
    print(json.dumps(summary, indent=1, default=float)[:4000])


if __name__ == "__main__":
    main()

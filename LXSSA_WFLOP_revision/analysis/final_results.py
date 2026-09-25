"""Final tables, figures and statistics: hybrid LX-SSA-VNS (LX-SSA + original basic VNS, method id "LXBV")
against LX-SSA, SSA, PSO, DE, the original basic VNS ("BVNS", shown as "VNS") and MS-SLSQP, plus the
component ablation (SSA, LX-SSA, VNS, SSA-VNS, LX-SSA-VNS).

Outputs go to ../figures_fresh/ (PDF figures) and fresh_tables.tex / fresh_stats.txt here.
"""
import os, json
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker
from scipy.stats import wilcoxon, friedmanchisquare, rankdata, norm

ALGS = ["LXBV", "LXSSA", "SSA", "PSO", "DE", "BVNS", "SLSQP"]
FOCUS = "LXBV"
LAB = {"LXBV": "LX-SSA-VNS", "SSABV": "SSA-VNS", "LXVNS": "LX-SSA-VNS (mod.)", "LXSSA": "LX-SSA", "SSA": "SSA", "PSO": "PSO", "DE": "DE", "VNS": "VNS",
       "SLSQP": "MS-SLSQP", "BVNS": "VNS"}
# categorical palette (validated: scripts/validate_palette.js, light mode, legend order)
COL = {"LXBV": "#2a78d6", "SSABV": "#e34948", "LXVNS": "#2a78d6", "LXSSA": "#4a3aa7", "SSA": "#eb6834", "PSO": "#1baf7a", "DE": "#eda100",
       "VNS": "#e87ba4", "SLSQP": "#008300", "BVNS": "#e87ba4"}
MRK = {"LXBV": "*", "SSABV": "X", "LXVNS": "*", "LXSSA": "o", "SSA": "s", "PSO": "^", "DE": "v", "VNS": "D", "SLSQP": "P", "BVNS": "D"}
LS = {"LXBV": "-", "SSABV": (0, (4, 2, 1, 2)), "LXVNS": "-", "LXSSA": (0, (1, 1)), "SSA": "--", "PSO": "-.", "DE": ":", "VNS": (0, (5, 1)),
      "SLSQP": (0, (3, 1, 1, 1)), "BVNS": (0, (5, 1))}
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#e4e3df"
RADII = {500: 10, 750: 12, 1000: 15}
DSN = {1: "Data Set I", 2: "Data Set II"}
FIG = "../figures_final"
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
        d = piv[FOCUS] - piv[b]; nz = d[np.abs(d) > 1e-9]
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


def legend_row(fig, y=1.02, algs=None):
    algs = algs or ALGS
    h = [plt.Line2D([], [], color=COL[a], ls=LS[a], marker=MRK[a], ms=7 if MRK[a] == "*" else 4, label=LAB[a])
         for a in algs]
    fig.legend(handles=h, loc="lower center", ncol=len(algs), bbox_to_anchor=(0.5, y), fontsize=7.5)


# ------------------------------------------------------------------ main
def main():
    GA = pd.concat([pd.read_csv("fresh_grid.csv"), pd.read_csv("fresh_vgrid.csv"), pd.read_csv("fresh_bgrid.csv")],
                   ignore_index=True)
    GA["LossPct"] = 100 * GA.WakeLoss / GA.Ideal
    G = GA[GA.Algorithm.isin(ALGS)].reset_index(drop=True)
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
                                 MeanDiff=fm.get(FOCUS, np.nan) - fm.get(x["Baseline"], np.nan)))
        stats[(ds, r, n)] = dict(ranks=ranks, feas=feas.to_dict(), fm=fm.to_dict())
    C = pd.DataFrame(caserows)
    C["Outcome"] = np.where(C.PHolm < 0.05, np.where(C.RB > 0, "W", "L"), "T")
    C.to_csv("final_case_tests.csv", index=False)

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
\scriptsize\setlength{\tabcolsep}{2pt}
\resizebox{\textwidth}{!}{%%
\begin{tabular}{ccccccccc}
\toprule
$N$ & Ideal & LX-SSA-VNS & LX-SSA & SSA & PSO & DE & VNS & MS-SLSQP \\
\midrule
""" % (DSN[ds], r, ds, r) + "\n".join(lines) + r"""
\bottomrule
\end{tabular}}
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
            lines.append(f"{'I' if ds == 1 else 'II'} & {r} & {s.Turbines.nunique()} & " + " & ".join(cells) + " \\\\")
    tot = []
    for b in ALGS[1:]:
        x = C[C.Baseline == b].Outcome.value_counts()
        tot.append(f"{x.get('W',0)}/{x.get('T',0)}/{x.get('L',0)}")
    lines.append("\\midrule\nAll & & 68 & " + " & ".join(tot) + " \\\\")
    tex.append(r"""\begin{table}[!t]
\centering
\caption{Pairwise outcome of the hybrid LX-SSA-VNS against each method: number of cases in which LX-SSA-VNS is significantly better / not significantly different / significantly worse (two-sided Wilcoxon signed-rank test on 30 seed-paired runs, Holm-adjusted over the six comparisons of each case, $\alpha=0.05$).}
\label{tab:wtl}
\scriptsize\setlength{\tabcolsep}{2.5pt}
\begin{tabular}{llccccccc}
\toprule
DS & $r$ (m) & Cases & LX-SSA & SSA & PSO & DE & VNS & \begin{tabular}{@{}c@{}}MS-\\SLSQP\end{tabular} \\
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
    lx_i = ALGS.index(FOCUS)
    zl = [(ALGS[j], (avg[j] - avg[lx_i]) / se) for j in range(k) if j != lx_i]
    phl = holm([2 * norm.sf(abs(z)) for _, z in zl])
    rows = []
    for j, a in enumerate(ALGS):
        vb = "--" if j == best_i else fmt_p(ph[[x[0] for x in zp].index(a)])
        vl = "--" if j == lx_i else fmt_p(phl[[x[0] for x in zl].index(a)])
        wins = int((cr[:, j] == 1).sum())
        feas = G[G.Algorithm == a].Feasible.mean() * 100
        rows.append(f"{LAB[a]} & {avg[j]:.2f} & {wins} & {feas:.1f} & {vl} \\\\")
    order = np.argsort(avg)
    tex.append(r"""\begin{table}[!t]
\centering
\caption{Case-level Friedman analysis over all 68 benchmark cases (blocks: cases; ranks of the mean feasible objective, 1 = best). Friedman $\chi^2_F=%.1f$ (%d d.f.), $p=%s$; Iman--Davenport $F_F=%.1f$. Last column: Holm-adjusted $p$ value of the average-rank $z$ test against LX-SSA-VNS (best-ranked method: %s); ``Best'': cases in which the method alone has the highest mean (remaining cases are ties); ``Feas.'': percentage of feasible runs.}
\label{tab:friedman68}
\scriptsize\setlength{\tabcolsep}{3pt}
\begin{tabular}{lcccc}
\toprule
Method & Avg.\ rank & Best & Feas.\ (\%%) & $p$ vs.\ LX-SSA-VNS \\
\midrule
""" % (chi, k - 1, fmt_p(pf).strip("$"), ff, LAB[ALGS[best_i]]) +
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
    fig, ax = plt.subplots(figsize=(3.4, 2.1))
    for pos, j in enumerate(order[::-1]):
        a = ALGS[j]
        ax.plot([1, avg[j]], [pos, pos], color=GRID, lw=2, zorder=1)
        ax.scatter(avg[j], pos, s=90 if MRK[a] == "*" else 36, color=COL[a], marker=MRK[a], zorder=3,
                   edgecolor="white", lw=0.6)
        ax.text(avg[j] + 0.08, pos, f"{avg[j]:.2f}", va="center", fontsize=7, color=INK)
    ax.set_yticks(range(k)); ax.set_yticklabels([LAB[ALGS[j]] for j in order[::-1]])
    ax.set_xlim(1, len(ALGS) + 0.3); ax.set_xlabel("Average rank over 68 cases (1 = best)")
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
                ax.plot(m.index, m.values, color=COL[a], ls=LS[a], marker=MRK[a],
                        ms=6 if MRK[a] == "*" else 3.5, label=LAB[a])
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
                ax.plot(m.index, m.values, color=COL[a], ls=LS[a], marker=MRK[a], ms=6 if MRK[a] == "*" else 3.5)
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
                        ax.plot(x[k_last[-1]], med[k_last[-1]], marker=MRK[a], color=COL[a],
                                ms=7 if MRK[a] == "*" else 4)
                ax.set_yscale("log")
                fmt = matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:g}")
                ax.yaxis.set_major_locator(matplotlib.ticker.LogLocator(subs=(1.0, 2.0, 3.0, 5.0)))
                ax.yaxis.set_major_formatter(fmt)
                ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
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
            ax.set_xticks(range(1, len(ALGS) + 1)); ax.set_xticklabels([LAB[a] for a in ALGS], rotation=40, fontsize=6)
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
            lx = sub[sub.Algorithm == FOCUS].sort_values("Objective").iloc[-1]
            top = sub.sort_values("Objective").iloc[-1]
            for row, mk, colr, lab in ((top, "s", COL[top.Algorithm], f"best overall ({LAB[top.Algorithm]})"),
                                       (lx, "o", COL[FOCUS], "best LX-SSA-VNS")):
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
    pd.DataFrame(best_rows).drop_duplicates().to_csv("final_best_layouts_maxN.csv", index=False)

    # ---------- runtime ----------
    rt = (G.Seconds / G.Calls * 1000).groupby(G.Algorithm).mean().reindex(ALGS)
    rn = G.groupby(["Algorithm", "Turbines"]).Seconds.mean().unstack().reindex(ALGS)
    summary["ms_per_call"] = rt.round(3).to_dict()
    summary["sec_per_run_range"] = {a: [float(rn.loc[a].min()), float(rn.loc[a].max())] for a in ALGS}

    # ---------- Horns Rev ----------
    import hornsrev_model as hr
    H = {n: pd.concat([pd.read_csv(f"fresh_hr{n}.csv"), pd.read_csv(f"fresh_vhr{n}.csv"),
                    pd.read_csv(f"fresh_bhr{n}.csv")], ignore_index=True).query("Algorithm in @ALGS").reset_index(drop=True)
         for n in (16, 80)}
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
                c.append("--" if a == FOCUS else fmt_p(t16[a]["PHolm"]))
        lines.append(f"{LAB[a]} & " + " & ".join(c) + " \\\\")
    tex.append((r"""\begin{table*}[!t]
\centering
\caption{Horns Rev~1 site case (real Vestas V80 power and thrust curves, measured 12-sector wind climate, installed farm outline; Jensen wake with $k=0.04$; minimum spacing $4D=320$~m): AEP in GWh/yr, mean (SD) over the feasible runs at 6,030 objective calls (30 seeds for the 16-turbine block, 10 for the 80-turbine farm), mean wake loss relative to the wake-free AEP (``Loss'', \%%), feasible runs, and Holm-adjusted Wilcoxon signed-rank $p$ of LX-SSA-VNS vs.\ each method (16-turbine block; run-level Friedman $p=%s$).}
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

    # ---------- ablation: budget split and phase contribution ----------
    Sp = pd.concat([pd.read_csv("fresh_bsplit.csv"), G[G.Algorithm == FOCUS]], ignore_index=True)
    lines, split_rows = [], []
    for (ds, r, n), s in Sp.groupby(["Dataset", "Radius", "Turbines"]):
        if s.Algorithm.nunique() < 3:
            continue
        piv = s.assign(S=goodness(s, "Objective")).pivot(index="Seed", columns="Algorithm", values="S")
        fm = s[s.Feasible].groupby("Algorithm").Objective.mean()
        fc = s.groupby("Algorithm").Feasible.sum()
        ps = []
        for b in ("LXBV25", "LXBV75"):
            d = piv[FOCUS] - piv[b]; nz = d[np.abs(d) > 1e-9]
            ps.append(wilcoxon(nz).pvalue if len(nz) else 1.0)
        ph = holm(ps)
        cells = []
        for a in ("LXBV25", FOCUS, "LXBV75"):
            v = fm.get(a, np.nan)
            c = "--" if not np.isfinite(v) else f"{v:.1f}"
            if fc.get(a, 0) < 30:
                c += f"$^{{{int(fc.get(a, 0))}}}$"
            cells.append(c)
        best = max(("LXBV25", FOCUS, "LXBV75"), key=lambda a: fm.get(a, -np.inf))
        split_rows.append(dict(Dataset=ds, Radius=r, Turbines=n, best=best, p25=ph[0], p75=ph[1]))
        lines.append(f"{ds} & {r} & {n} & " + " & ".join(cells) + f" & {fmt_p(ph[0])} & {fmt_p(ph[1])} \\\\")
    tex.append(r"""\begin{table}[!t]
\centering
\caption{Sensitivity of LX-SSA-VNS to the budget split between the LX-SSA and VNS phases (25\%%, 50\%% and 75\%% of the 6,030 calls for LX-SSA): mean benchmark objective of the feasible runs (superscript: feasible runs when fewer than 30) and Holm-adjusted Wilcoxon signed-rank $p$ of the 50\%% split against the 25\%% and 75\%% splits.}
\label{tab:split}
\scriptsize\setlength{\tabcolsep}{2.5pt}
\begin{tabular}{cccccccc}
\toprule
DS & $r$ & $N$ & 25\%% & 50\%% & 75\%% & $p$ (25) & $p$ (75) \\\\
\midrule
""".replace("%%", "%") + "\n".join(lines) + r"""
\bottomrule
\end{tabular}
\end{table}
""")
    summary["split"] = pd.DataFrame(split_rows).to_dict("records")
    # phase contribution: wake loss at the end of the LX-SSA phase (call 3,030) vs final
    Hy = G[(G.Algorithm == FOCUS) & G.Feasible]
    gains, still_infeasible = [], 0
    for _, row in Hy.iterrows():
        c = np.array([float(v) for v in row.Curve.split(";")])
        at = c[3030 // CHECK - 1]
        if not np.isfinite(at):
            still_infeasible += 1; continue
        l1, l2 = row.Ideal - at, row.Ideal - row.Objective
        if l1 > 1e-9:
            gains.append(100 * (l1 - l2) / l1)
    summary["phase2_loss_reduction_pct"] = dict(mean=float(np.mean(gains)), median=float(np.median(gains)),
                                                n=len(gains), infeasible_at_switch=still_infeasible)
    # ---------- component ablation: SSA, LX-SSA, VNS, SSA-VNS, LX-SSA-VNS ----------
    ABL = ["SSABV", "LXBV", "BVNS", "LXSSA", "SSA"]          # legend order (palette validated)
    CONTR = [("LXBV", "LXSSA", "VNS phase (after LX-SSA)"),
             ("LXBV", "BVNS", "LX-SSA start vs.\\ best initial point"),
             ("LXBV", "SSABV", "Laplace step inside the hybrid"),
             ("SSABV", "SSA", "VNS phase (after SSA)"),
             ("SSABV", "BVNS", "SSA start vs.\\ best initial point"),
             ("LXSSA", "SSA", "Laplace step alone")]
    GB = GA[GA.Algorithm.isin(ABL)]
    arows = []
    for (ds, r, n), sub in GB.groupby(["Dataset", "Radius", "Turbines"]):
        piv = sub.assign(S=goodness(sub, "Objective")).pivot(index="Seed", columns="Algorithm", values="S")
        lm = sub[sub.Feasible].groupby("Algorithm").LossPct.mean()
        rr = []
        for a, b, _ in CONTR:
            d = piv[a] - piv[b]; nz = d[np.abs(d) > 1e-9]
            if len(nz) == 0:
                p, rb = 1.0, 0.0
            else:
                p = wilcoxon(nz).pvalue; rk = rankdata(np.abs(nz))
                rb = (rk[nz > 0].sum() - rk[nz < 0].sum()) / rk.sum()
            rr.append(dict(Dataset=ds, Radius=r, Turbines=n, A=a, B=b, P=p, RB=rb,
                           DLoss=lm.get(a, np.nan) - lm.get(b, np.nan)))
        for x, h in zip(rr, holm([x["P"] for x in rr])):
            x["PHolm"] = h
        arows += rr
    A = pd.DataFrame(arows)
    A["Outcome"] = np.where(A.PHolm < 0.05, np.where(A.RB > 0, "W", "L"), "T")
    A.to_csv("final_ablation_tests.csv", index=False)
    am = pd.DataFrame([[GB[(GB.Dataset == k[0]) & (GB.Radius == k[1]) & (GB.Turbines == k[2]) & GB.Feasible
                           & (GB.Algorithm == a)].Objective.mean() for a in ABL] for k in sorted(stats)], columns=ABL)
    acr = np.vstack([rankdata(-np.nan_to_num(row, nan=-np.inf)) for row in am.values])
    achi, apf = friedmanchisquare(*[acr[:, j] for j in range(len(ABL))])
    aavg = acr.mean(0)
    comp = {"SSA": ("SSA", "--"), "LXSSA": ("LX-SSA", "--"), "BVNS": ("--", "VNS"),
            "SSABV": ("SSA", "VNS"), "LXBV": ("LX-SSA", "VNS")}
    rows1 = []
    for j in np.argsort(aavg):
        a = ABL[j]
        rows1.append(f"{LAB[a]} & {comp[a][0]} & {comp[a][1]} & {aavg[j]:.2f} & "
                     f"{GB[GB.Algorithm == a].Feasible.mean() * 100:.1f} \\\\")
    rows2 = []
    for a, b, what in CONTR:
        x = A[(A.A == a) & (A.B == b)]
        oc = x.Outcome.value_counts()
        rows2.append(f"{LAB[a]} vs.\\ {LAB[b]} & {what} & {oc.get('W', 0)}/{oc.get('T', 0)}/{oc.get('L', 0)} & "
                     f"${x.DLoss.mean():+.2f}$ \\\\")
    tex.append(r"""\begin{table}[!t]
\centering
\caption{Component ablation over the 68 benchmark cases (30 seed-paired runs, 6,030 calls). Top: phase-1 and phase-2 components, average rank of the mean feasible objective among the five variants (Friedman $\chi^2_F=%.1f$, 4 d.f., $p=%s$) and percentage of feasible runs. Bottom: planned contrasts, number of cases in which the first variant is significantly better / not different / worse (Wilcoxon signed-rank, Holm-adjusted over the six contrasts of each case, $\alpha=0.05$), and mean difference of the wake loss $\overline{\Delta L}$ (percentage points; negative = first variant better).}
\label{tab:ablation}
\scriptsize\setlength{\tabcolsep}{2.5pt}
\begin{tabular}{lcccc}
\toprule
Variant & Phase 1 & Phase 2 & Avg.\ rank & Feas.\ (\%%) \\
\midrule
""" % (achi, fmt_p(apf).strip("$")) + "\n".join(rows1) + r"""
\bottomrule
\end{tabular}

\smallskip
\begin{tabular}{l>{\raggedright\arraybackslash}p{2.6cm}cc}
\toprule
Contrast & Isolates & W/T/L & $\overline{\Delta L}$ \\
\midrule
""" + "\n".join(rows2) + r"""
\bottomrule
\end{tabular}
\end{table}
""")
    summary["ablation"] = dict(chi2=achi, p=apf, avg_rank={a: float(aavg[j]) for j, a in enumerate(ABL)},
                               feasible_pct={a: float(GB[GB.Algorithm == a].Feasible.mean() * 100) for a in ABL},
                               contrasts={f"{a}-{b}": dict(A[(A.A == a) & (A.B == b)].Outcome.value_counts().to_dict(),
                                                            dloss=float(A[(A.A == a) & (A.B == b)].DLoss.mean()))
                                          for a, b, _ in CONTR})
    # phase-2 gain of SSA-VNS as well
    for key, alg in (("phase2_ssabv", "SSABV"),):
        Hy2 = GA[(GA.Algorithm == alg) & GA.Feasible]; g2, inf2 = [], 0
        for _, row in Hy2.iterrows():
            c = np.array([float(v) for v in row.Curve.split(";")]); at = c[3030 // CHECK - 1]
            if not np.isfinite(at):
                inf2 += 1; continue
            l1, l2 = row.Ideal - at, row.Ideal - row.Objective
            if l1 > 1e-9:
                g2.append(100 * (l1 - l2) / l1)
        summary[key] = dict(mean=float(np.mean(g2)), median=float(np.median(g2)), infeasible_at_switch=inf2)
    # ablation convergence figure (largest N)
    fig, axes = plt.subplots(2, 3, figsize=(7.1, 4.0))
    for i, ds in enumerate((1, 2)):
        for j, (r, n) in enumerate(RADII.items()):
            ax = axes[i, j]
            sub = GB[(GB.Dataset == ds) & (GB.Radius == r) & (GB.Turbines == n)]
            ideal = sub.Ideal.iloc[0]; x = np.arange(1, 202) * CHECK
            for a in ABL:
                M = curves(sub[sub.Algorithm == a]); L = 100 * (ideal - M) / ideal
                frac = np.mean(np.isfinite(L), 0)
                med = np.where(frac >= 0.5, np.nanmedian(np.where(np.isfinite(L), L, np.inf), 0), np.nan)
                ax.plot(x, med, color=COL[a], ls=LS[a], lw=1.2)
                k_last = np.where(np.isfinite(med))[0]
                if len(k_last):
                    ax.plot(x[k_last[-1]], med[k_last[-1]], marker=MRK[a], color=COL[a], ms=7 if MRK[a] == "*" else 4)
            ax.axvline(3030, color=MUTED, lw=0.7, ls=(0, (1, 2)))
            ax.set_yscale("log")
            ax.yaxis.set_major_locator(matplotlib.ticker.LogLocator(subs=(1.0, 2.0, 3.0, 5.0)))
            ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:g}"))
            ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
            ax.set_title(f"{DSN[ds]}, $r$ = {r} m, $N$ = {n}", fontsize=8, color=INK)
            if i == 1: ax.set_xlabel("Objective-function calls")
            if j == 0: ax.set_ylabel("Median best wake loss (%)")
    legend_row(fig, 1.0, ABL)
    fig.tight_layout()
    save(fig, "ablation_convergence")
    open("final_tables.tex", "w").write("\n".join(tex))
    json.dump(summary, open("final_summary.json", "w"), indent=1, default=float)
    print(json.dumps(summary, indent=1, default=float)[:4000])


if __name__ == "__main__":
    main()

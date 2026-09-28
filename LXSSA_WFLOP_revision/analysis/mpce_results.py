"""Results and statistics pipeline for the MPCE resubmission: proposed method SSA-VNS (id "SSABV").

Usage:  python mpce_results.py [--data-dir D] [--out-dir O] [--fig-dir F] [--procs 4]
                               [--partial] [--skip-robust] [--only-robust]

Inputs (all runs at 6,030 calls with random initialization unless stated)
  fresh_grid.csv  (LXSSA, SSA, PSO [old settings; development fallback for PSOC], DE,
                   SLSQP [old platform; development fallback for the SLSQP rerun]; VNS ignored)
  fresh_vgrid.csv (BVNS = basic VNS, "VNS"),  fresh_bgrid.csv (LXBV, SSABV)
  fresh_hr16.csv, fresh_vhr16.csv, fresh_bhr16.csv (Horns Rev 1, 16 turbines; AEP/IdealAEP)
  mpce_<exp>_s<i>of<k>.csv from mpce_experiments.py (read from --data-dir; all shards of an
  experiment are merged; an experiment is used only when all k shards are present, unless
  --partial is given): rsvns, psoc, slsqp, ssasplit, hr16new, feas, b30k, b120k, iea16, iea36.
  iea37_published_results.csv (optional; columns case, participant/algorithm, AEP).

Outputs (in --out-dir, default: this folder)
  mpce_tab_<name>.tex (one table per file, ready to \\input), mpce_supplementary.tex (per-case
  tables), mpce_summary.json (every number quoted in the text), mpce_case_stats.csv,
  mpce_case_tests.csv, mpce_ablation_tests.csv, mpce_best_layouts_maxN.csv; figures (pdf + png)
  in --fig-dir (default ../figures_mpce).

Ranking rule (survivorship-bias correction; used for every case-level ranking below).
  Ranking methods by the mean objective of their *feasible* runs rewards a method that is
  feasible only in a few lucky runs. Within each case a method is therefore ranked normally
  (by the mean feasible objective, higher = better) only if at least half of its runs (15 of
  30; 5 of 10) are feasible. Methods with fewer feasible runs receive the worst ranks, below all
  qualifying methods, ordered by their number of feasible runs (more = better; ties share the
  average rank).
"""
import os, sys, glob, re, json, argparse, time
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker
from scipy.stats import wilcoxon, friedmanchisquare, rankdata, norm, f as f_dist

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

MAIN = ["SSABV", "SSA", "LXSSA", "PSOC", "DE", "BVNS", "SLSQP"]
FOCUS = "SSABV"
ABL = ["SSABV", "LXBV", "RSVNS", "BVNS", "SSA", "LXSSA"]
M9 = ["SSABV", "LXBV", "RSVNS", "BVNS", "SSA", "LXSSA", "PSOC", "DE", "SLSQP"]
LAB = {"SSABV": "SSA-VNS", "SSABV25": "SSA-VNS (25\\%)", "SSABV75": "SSA-VNS (75\\%)", "LXBV": "LX-SSA-VNS",
       "RSVNS": "RS-VNS", "BVNS": "VNS", "SSA": "SSA", "LXSSA": "LX-SSA", "PSOC": "PSO", "PSO": "PSO (old)",
       "DE": "DE", "SLSQP": "MS-SLSQP"}
# Colour follows the method (categorical palette of final_results.py, validated in slot order:
# blue, orange, aqua, yellow, magenta, green, violet, red). RS-VNS is the ninth entity: neutral
# grey with its own marker and dash pattern (composite encoding instead of a generated hue).
COL = {"SSABV": "#2a78d6", "SSA": "#eb6834", "PSOC": "#1baf7a", "PSO": "#1baf7a", "DE": "#eda100",
       "BVNS": "#e87ba4", "SLSQP": "#008300", "LXSSA": "#4a3aa7", "LXBV": "#e34948", "RSVNS": "#8d8b85"}
MRK = {"SSABV": "*", "SSA": "s", "PSOC": "^", "PSO": "^", "DE": "v", "BVNS": "D", "SLSQP": "P", "LXSSA": "o",
       "LXBV": "X", "RSVNS": "h"}
LS = {"SSABV": "-", "SSA": "--", "PSOC": "-.", "PSO": "-.", "DE": ":", "BVNS": (0, (5, 1)),
      "SLSQP": (0, (3, 1, 1, 1)), "LXSSA": (0, (1, 1)), "LXBV": (0, (4, 2, 1, 2)), "RSVNS": (0, (6, 2, 2, 2))}
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#e4e3df"
RADII = {500: 10, 750: 12, 1000: 15}
MID = {500: 6, 750: 8, 1000: 10}
DSN = {"1": "Data Set I", "2": "Data Set II", "HR": "Horns Rev 1"}
CASE = ["Dataset", "Radius", "Turbines"]
KEY = ["Algorithm", "Dataset", "Radius", "Turbines", "Seed", "Budget", "Init"]
SPLITCASES = [(ds, r, n) for ds in ("1", "2") for r, n in ((500, 6), (750, 8), (1000, 10), (500, 10), (750, 12), (1000, 15))]
LARGE = [(ds, r, n) for ds in ("1", "2") for r, n in RADII.items()]
HR16 = ("HR", 0, 16)
IEA_NAME = {16: "IEA37 CS1, 16 turbines", 36: "IEA37 CS1, 36 turbines"}   # text names (both are Case Study 1 scenarios)
BASELINE_IEA = {16: 366941.57116, 36: 737883.09851}    # official example-layout AEP (MWh), used if the published file is absent
BUDGETS = [6030, 30030, 120030]
INSTALLED_HR16 = 139.51          # GWh/yr, recomputed with hornsrev_model below when available
RANK_RULE = ("Within each case a method is ranked by the mean objective of its feasible runs only if at "
             "least half of its runs (15 of 30; 5 of 10) are feasible; methods with fewer feasible runs "
             "receive the worst ranks, below all qualifying methods, ordered by their number of feasible "
             "runs (more = better, ties averaged). This removes the survivorship bias of ranking by the "
             "mean feasible objective.")

plt.rcParams.update({"font.family": "serif", "font.size": 8, "axes.edgecolor": MUTED,
                     "axes.labelcolor": INK, "xtick.color": MUTED, "ytick.color": MUTED,
                     "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.5,
                     "axes.spines.top": False, "axes.spines.right": False,
                     "legend.frameon": False, "lines.linewidth": 1.4, "pdf.fonttype": 42})

LOG = []


def log(msg=""):
    print(msg, flush=True)
    LOG.append(msg)


# ------------------------------------------------------------------ statistics helpers
def holm(p):
    p = np.asarray(p, float); o = np.argsort(p); m = len(p); out = np.empty(m); run = 0
    for i, k in enumerate(o):
        run = max(run, min(1.0, (m - i) * p[k])); out[k] = run
    return out


def fmt_p(p):
    if not np.isfinite(p):
        return "--"
    if p < 1e-3:
        m, e = f"{p:.2e}".split("e")
        return f"${m}\\times10^{{{int(e)}}}$"
    return f"{p:.3f}"


def goodness(df):
    """Seed-level score, higher = better; infeasible runs below every feasible run."""
    return np.where(df.Feasible, df.Objective, -1e12 - np.maximum(0, 1e4 - df.MinSpacing))


def wil(d):
    """Two-sided Wilcoxon signed-rank on paired differences d (zeros dropped) + rank-biserial r."""
    d = np.asarray(d, float); d = d[np.isfinite(d)]
    nz = d[np.abs(d) > 1e-9]
    if len(nz) == 0:
        return 1.0, 0.0, 0
    p = wilcoxon(nz).pvalue
    rk = rankdata(np.abs(nz))
    return float(p), float((rk[nz > 0].sum() - rk[nz < 0].sum()) / rk.sum()), len(nz)


def rank_rule(fmean, nfeas, nruns):
    """Ranks (1 = best) of one case under the survivorship-bias rule (see module docstring)."""
    # means are rounded to 1e-6 (objective units) so that exact ties (e.g. every method at the wake-free
    # optimum for N = 2) are not broken by floating-point noise of a re-evaluation
    fmean = np.round(np.asarray(fmean, float), 6); nfeas = np.asarray(nfeas, float); nruns = np.asarray(nruns, float)
    q = (nfeas >= np.ceil(nruns / 2)) & np.isfinite(fmean)
    r = np.empty(len(fmean))
    if q.any():
        r[q] = rankdata(-fmean[q])
    if (~q).any():
        r[~q] = q.sum() + rankdata(-nfeas[~q])
    return r


def rank_old(fmean, *_):
    """Ranking of final_results.py (mean feasible objective; no feasible run = tied last)."""
    return rankdata(-np.nan_to_num(np.asarray(fmean, float), nan=-np.inf))


def case_stats(D, methods, rank_fn=rank_rule):
    """Per case x method: runs, feasible runs, mean/SD/best objective and mean/SD loss of the feasible runs, rank."""
    D = D[D.Algorithm.isin(methods)]
    g = D.groupby(CASE + ["Algorithm"])
    S = pd.DataFrame({"N": g.size(), "NFeas": g.Feasible.sum()})
    gf = D[D.Feasible].groupby(CASE + ["Algorithm"])
    S["Mean"] = gf.Objective.mean(); S["SD"] = gf.Objective.std(); S["Best"] = gf.Objective.max()
    S["Loss"] = gf.LossPct.mean(); S["LossSD"] = gf.LossPct.std()
    S = S.reset_index()
    S["NFeas"] = S.NFeas.astype(int)
    S["Qualified"] = S.NFeas >= np.ceil(S.N / 2)
    S["Rank"] = np.nan
    for _, idx in S.groupby(CASE).groups.items():
        s = S.loc[idx]
        S.loc[idx, "Rank"] = rank_fn(s.Mean.values, s.NFeas.values, s.N.values)
    return S


def rank_matrix(S, methods):
    R = S.pivot_table(index=CASE, columns="Algorithm", values="Rank").reindex(columns=methods)
    R = R.dropna()
    return R


def friedman_block(R, methods, focus=FOCUS):
    """Case-level Friedman, Iman-Davenport, average ranks, Holm post hoc z tests vs focus and vs best."""
    X = R[methods].values
    n, k = X.shape
    chi, pf = friedmanchisquare(*[X[:, j] for j in range(k)])
    ff = (n - 1) * chi / (n * (k - 1) - chi)
    pff = f_dist.sf(ff, k - 1, (k - 1) * (n - 1))
    avg = X.mean(0); se = np.sqrt(k * (k + 1) / (6 * n))
    fi = methods.index(focus); bi = int(np.argmin(avg))
    zf = {a: (avg[j] - avg[fi]) / se for j, a in enumerate(methods) if j != fi}
    pz = dict(zip(zf, holm([2 * norm.sf(abs(z)) for z in zf.values()])))
    zb = {a: (avg[j] - avg[bi]) / se for j, a in enumerate(methods) if j != bi}
    pb = dict(zip(zb, holm([2 * norm.sf(abs(z)) for z in zb.values()])))
    best = {a: int(((X[:, j] == 1) & ((X == 1).sum(1) == 1)).sum()) for j, a in enumerate(methods)}
    return dict(n_cases=n, chi2=float(chi), p=float(pf), iman_davenport=float(ff), iman_davenport_p=float(pff),
                avg_rank={a: float(avg[j]) for j, a in enumerate(methods)}, best_ranked=methods[bi],
                z_vs_focus={a: float(v) for a, v in zf.items()}, p_holm_vs_focus={a: float(v) for a, v in pz.items()},
                p_holm_vs_best={a: float(v) for a, v in pb.items()}, sole_best_count=best)


def paired_vs(sub, focus, others):
    """Seed-paired Wilcoxon of focus vs each other method (one case), Holm over the comparisons."""
    piv = sub.assign(S=goodness(sub)).pivot_table(index="Seed", columns="Algorithm", values="S")
    rows = []
    for b in others:
        if focus not in piv or b not in piv:
            continue
        d = (piv[focus] - piv[b]).dropna()
        p, rb, nn = wil(d.values)
        rows.append(dict(Baseline=b, P=p, RB=rb, NPairs=len(d)))
    for r, h in zip(rows, holm([r["P"] for r in rows]) if rows else []):
        r["PHolm"] = float(h)
        r["Outcome"] = ("W" if r["RB"] > 0 else "L") if h < 0.05 else "T"
    return rows


def case_mean_wilcoxon(S, focus, others, val="Loss"):
    """Wilcoxon signed-rank on per-case mean wake losses (%), focus vs each method, Holm over methods.

    d = L(other) - L(focus) > 0 means focus better. Cases in which exactly one of the two methods
    qualifies (>= half of its runs feasible) are counted for the qualifying method with a difference
    larger than every observed one (consistent with the ranking rule); cases in which neither
    qualifies are dropped. The plain version (only cases where both qualify) is reported too.
    """
    P = S.pivot_table(index=CASE, columns="Algorithm", values=val)
    Q = S.pivot_table(index=CASE, columns="Algorithm", values="Qualified").astype(float)
    out = {}
    for b in others:
        if b not in P or focus not in P:
            continue
        qa, qb = Q[focus] > 0.5, Q[b] > 0.5
        both = qa & qb
        d = (P[b] - P[focus])[both].values
        big = 10 * (np.nanmax(np.abs(d)) if len(d) and np.isfinite(d).any() else 1.0) + 1.0
        dfull = np.r_[d, np.full(int((qa & ~qb).sum()), big), np.full(int((~qa & qb).sum()), -big)]
        p, rb, nn = wil(dfull)
        p2, rb2, nn2 = wil(d)
        out[b] = dict(p=p, rb=rb, n_cases=len(dfull), only_focus_qualified=int((qa & ~qb).sum()),
                      only_other_qualified=int((~qa & qb).sum()), neither_qualified=int((~qa & ~qb).sum()),
                      mean_dloss_pp=float(np.mean(-d)) if len(d) else np.nan,
                      median_dloss_pp=float(np.median(-d)) if len(d) else np.nan,
                      focus_lower_loss_cases=int((d > 1e-9).sum()), other_lower_loss_cases=int((d < -1e-9).sum()),
                      p_both_qualified=p2, rb_both_qualified=rb2, n_both_qualified=int(both.sum()))
    for key in ("p", "p_both_qualified"):
        bs = list(out)
        for b, h in zip(bs, holm([out[b][key] for b in bs]) if bs else []):
            out[b][key + "_holm"] = float(h)
    return out


def curves(sub):
    return np.array([[float(v) for v in c.split(";")] for c in sub.Curve])


def curve_x(sub, ncp):
    b = int(sub.Budget.iloc[0])
    return np.arange(1, ncp + 1) * max(1, (b - 30) // 200)


def coords(s):
    return np.array([[float(v) for v in p.split()] for p in s.split(";")])


def wtl_str(c):
    return f"{c.get('W', 0)}/{c.get('T', 0)}/{c.get('L', 0)}"


# ------------------------------------------------------------------ figure helpers
class Figs:
    def __init__(self, d):
        self.d = d; self.made = []

    def save(self, fig, name):
        os.makedirs(self.d, exist_ok=True)
        fig.savefig(f"{self.d}/{name}.pdf", bbox_inches="tight")
        fig.savefig(f"{self.d}/{name}.png", dpi=160, bbox_inches="tight")
        plt.close(fig)
        self.made.append(name)


def legend_row(fig, algs, y=1.0, ncol=None):
    h = [plt.Line2D([], [], color=COL[a], ls=LS[a], marker=MRK[a], ms=7 if MRK[a] == "*" else 4, label=LAB[a].replace("\\%", "%"))
         for a in algs]
    fig.legend(handles=h, loc="lower center", ncol=ncol or len(algs), bbox_to_anchor=(0.5, y), fontsize=7.5)


def log_axis(ax):
    ax.set_yscale("log")
    ax.yaxis.set_major_locator(matplotlib.ticker.LogLocator(subs=(1.0, 2.0, 3.0, 5.0)))
    ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:g}"))
    ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())


def median_curve(M, lower_better, q=None):
    """Median (or quantile q) across runs of the best-feasible curve; NaN where < 50% (75% for q3) feasible."""
    fin = np.isfinite(M); frac = fin.mean(0)
    fill = np.inf if lower_better else -np.inf
    X = np.where(fin, M, fill)
    if q is None:
        return np.where(frac >= 0.5, np.median(X, 0), np.nan)
    need = 0.5 if q <= 50 else 0.75
    with np.errstate(invalid="ignore"):
        return np.where(frac >= need, np.percentile(X, q, 0), np.nan)


def conv_panel(ax, sub, algs, loss=True, band=True):
    """Median best-feasible convergence curves (wake loss % if loss else objective) of algs."""
    ideal = sub.Ideal.iloc[0]
    for a in algs:
        s = sub[sub.Algorithm == a]
        if not len(s):
            continue
        M = curves(s)
        x = curve_x(s, M.shape[1])
        V = 100 * (ideal - M) / ideal if loss else M
        med = median_curve(V, loss)
        ax.plot(x, med, color=COL[a], ls=LS[a], lw=1.2)
        if band:
            q1, q3 = median_curve(V, loss, 25), median_curve(V, loss, 75)
            ax.fill_between(x, q1, q3, color=COL[a], alpha=0.12, lw=0)
        kl = np.where(np.isfinite(med))[0]
        if len(kl):
            ax.plot(x[kl[-1]], med[kl[-1]], marker=MRK[a], color=COL[a], ms=7 if MRK[a] == "*" else 4)


# ------------------------------------------------------------------ data loading
def std_cols(df, budget=6030, init="random"):
    df = df.copy()
    if "AEP" in df and "Objective" not in df:
        df = df.rename(columns={"AEP": "Objective", "IdealAEP": "Ideal"})
        df["Dataset"], df["Radius"] = "HR", 0
        df["WakeLoss"] = df.Ideal - df.Objective
    if "Budget" not in df:
        df["Budget"] = budget
    if "Init" not in df:
        df["Init"] = init
    df["Dataset"] = df.Dataset.astype(str).str.replace(r"\.0$", "", regex=True)
    df["Radius"] = pd.to_numeric(df.Radius, errors="coerce").fillna(0).astype(int)
    for c in ("Turbines", "Seed", "Budget"):
        df[c] = df[c].astype(int)
    if "Ideal" not in df:
        df["Ideal"] = np.nan
    if "WakeLoss" not in df:
        df["WakeLoss"] = df.Ideal - df.Objective
    if "MinSpacing" not in df:
        df["MinSpacing"] = np.nan
    df["Feasible"] = df.Feasible.astype(str).str.lower().isin(["true", "1", "1.0"])
    df["LossPct"] = 100 * df.WakeLoss / df.Ideal
    df["Source"] = df.get("Source", "")
    return df


def read_shards(exp, data_dir, partial):
    """Merge mpce_<exp>_s<i>of<k>.csv. Returns (DataFrame or None, status dict)."""
    files = glob.glob(os.path.join(data_dir, f"mpce_{exp}_s*of*.csv"))
    if not files:
        return None, dict(status="missing", shards="0", rows=0)
    by_k = {}
    for fn in files:
        m = re.search(rf"mpce_{exp}_s(\d+)of(\d+)\.csv$", os.path.basename(fn))
        if m:
            by_k.setdefault(int(m.group(2)), {})[int(m.group(1))] = fn
    complete = {k: v for k, v in by_k.items() if set(v) == set(range(k))}
    if complete:
        k = max(complete, key=lambda kk: sum(os.path.getsize(f) for f in complete[kk].values()))
        use, status = list(complete[k].values()), "complete"
        shards = f"{k}/{k}"
    elif partial:
        use = [f for v in by_k.values() for f in v.values()]
        status = "PARTIAL (used, --partial)"
        shards = ", ".join(f"{len(v)}/{k}" for k, v in by_k.items())
    else:
        shards = ", ".join(f"{len(v)}/{k}" for k, v in by_k.items())
        return None, dict(status="incomplete (ignored; use --partial to preview)", shards=shards, rows=0)
    df = pd.concat([pd.read_csv(f) for f in sorted(use)], ignore_index=True)
    df = std_cols(df)
    df = df.drop_duplicates(KEY, keep="last").reset_index(drop=True)
    df["Source"] = f"mpce_{exp}"
    return df, dict(status=status, shards=shards, rows=len(df),
                    methods=sorted(df.Algorithm.unique()), datasets=sorted(df.Dataset.unique()))


def load(data_dir, partial):
    avail, fallbacks = {}, []
    base = []
    for fn, drop in (("fresh_grid.csv", ["VNS"]), ("fresh_vgrid.csv", []), ("fresh_bgrid.csv", []),
                     ("fresh_hr16.csv", ["VNS"]), ("fresh_vhr16.csv", []), ("fresh_bhr16.csv", [])):
        p = os.path.join(HERE, fn)
        if os.path.exists(p):
            d = std_cols(pd.read_csv(p)); d = d[~d.Algorithm.isin(drop)]; d["Source"] = fn
            base.append(d)
            avail[fn] = dict(status="present", rows=len(d), methods=sorted(d.Algorithm.unique()))
        else:
            avail[fn] = dict(status="missing", rows=0)
    B = pd.concat(base, ignore_index=True)
    new = {}
    for exp in ("rsvns", "psoc", "slsqp", "ssasplit", "hr16new", "feas", "b30k", "b120k", "iea16", "iea36"):
        df, st = read_shards(exp, data_dir, partial)
        avail[f"mpce_{exp}"] = st
        if df is not None:
            new[exp] = df
    parts = [B]
    # SLSQP rerun replaces the old-platform SLSQP runs on the 68 cases
    if "slsqp" in new:
        parts[0] = B[~((B.Algorithm == "SLSQP") & B.Dataset.isin(["1", "2"]))]
    else:
        fallbacks.append("SLSQP (68 cases): mpce_slsqp missing -> using the old-platform SLSQP runs of fresh_grid.csv")
        m = (B.Algorithm == "SLSQP") & B.Dataset.isin(["1", "2"])
        B.loc[m, "Source"] = "FALLBACK fresh_grid SLSQP"
    fallbacks.append("SLSQP (Horns Rev 16, 6,030 calls): taken from fresh_hr16.csv (no rerun planned)")
    for exp in ("rsvns", "psoc", "slsqp", "ssasplit", "hr16new"):
        if exp in new:
            parts.append(new[exp])
    A6 = pd.concat(parts, ignore_index=True)
    # PSOC fallback: old-settings PSO relabelled (development only)
    for dom, has in ((["1", "2"], "psoc" in new), (["HR"], "hr16new" in new and (new["hr16new"].Algorithm == "PSOC").any())):
        if not has:
            f = A6[(A6.Algorithm == "PSO") & A6.Dataset.isin(dom)].copy()
            f["Algorithm"] = "PSOC"; f["Source"] = "FALLBACK old PSO"
            A6 = pd.concat([A6, f], ignore_index=True)
            fallbacks.append(f"PSOC ({'68 cases' if dom[0] == '1' else 'Horns Rev 16'}): constriction-PSO runs missing "
                             f"-> using the old-settings PSO runs relabelled as PSOC (DEVELOPMENT ONLY)")
    A6 = A6[A6.Algorithm != "PSO"]
    extra = [new[e] for e in ("feas", "b30k", "b120k", "iea16", "iea36") if e in new]
    ALL = pd.concat([A6] + extra, ignore_index=True)
    ALL = ALL.drop_duplicates(KEY, keep="first").reset_index(drop=True)
    return ALL, avail, fallbacks, new


# ------------------------------------------------------------------ LaTeX helpers
def table(env, caption, label, spec, header, lines, size="\\scriptsize", sep="3pt", resize=False, pos="!t"):
    body = "\n".join(lines)
    tab = f"\\begin{{tabular}}{{{spec}}}\n\\toprule\n{header} \\\\\n\\midrule\n{body}\n\\bottomrule\n\\end{{tabular}}"
    if resize:
        tab = "\\resizebox{\\textwidth}{!}{%\n" + tab + "}"
    return (f"\\begin{{{env}}}[{pos}]\n\\centering\n\\caption{{{caption}}}\n\\label{{{label}}}\n"
            f"{size}\\setlength{{\\tabcolsep}}{{{sep}}}\n{tab}\n\\end{{{env}}}\n")


# ------------------------------------------------------------------ main
def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default=HERE)
    ap.add_argument("--out-dir", default=HERE)
    ap.add_argument("--fig-dir", default=os.path.join(HERE, "..", "figures_mpce"))
    ap.add_argument("--procs", type=int, default=4)
    ap.add_argument("--partial", action="store_true")
    ap.add_argument("--skip-robust", action="store_true")
    ap.add_argument("--only-robust", action="store_true")
    args = ap.parse_args(argv)
    t0 = time.time()
    os.makedirs(args.out_dir, exist_ok=True)
    OUT = lambda fn: os.path.join(args.out_dir, fn)
    FG = Figs(os.path.abspath(args.fig_dir))
    tabs, supp = {}, []

    ALL, avail, fallbacks, new = load(args.data_dir, args.partial)
    log("=" * 78)
    log("DATA AVAILABILITY")
    for k, v in avail.items():
        extra = f"  methods={','.join(v['methods'])}" if v.get("methods") else ""
        log(f"  {k:20s} {v['status']:45s} shards={v.get('shards', '-'):>6s} rows={v['rows']:6d}{extra}")
    for f in fallbacks:
        log(f"  FALLBACK: {f}")
    log("=" * 78)
    summary = dict(generated=time.strftime("%Y-%m-%d %H:%M:%S"), data_availability=avail, fallbacks=fallbacks,
                   ranking_rule=RANK_RULE, labels={k: v.replace("\\%", "%") for k, v in LAB.items()})
    R6 = ALL[(ALL.Budget == 6030) & (ALL.Init == "random")]
    G = R6[R6.Dataset.isin(["1", "2"]) & R6.Algorithm.isin(MAIN)].reset_index(drop=True)
    missing_main = [a for a in MAIN if a not in set(G.Algorithm)]
    if missing_main:
        log(f"  WARNING: main methods without data: {missing_main}")
    MAINP = [a for a in MAIN if a not in missing_main]
    summary["main_methods"] = MAINP

    if args.only_robust:
        robustness(G, MAINP, args, summary, tabs, OUT, only=True)
        return

    # =========================================================== 1. main comparison
    log("\n[1] Main comparison (68 cases)")
    S = case_stats(G, MAINP)
    S.to_csv(OUT("mpce_case_stats.csv"), index=False)
    others = [a for a in MAINP if a != FOCUS]
    crow = []
    for (ds, r, n), sub in G.groupby(CASE):
        for x in paired_vs(sub, FOCUS, others):
            crow.append(dict(Dataset=ds, Radius=r, Turbines=n, **x))
    C = pd.DataFrame(crow)
    C.to_csv(OUT("mpce_case_tests.csv"), index=False)
    R = rank_matrix(S, MAINP)
    FR = friedman_block(R, MAINP)
    FR_old = friedman_block(rank_matrix(case_stats(G, MAINP, rank_old), MAINP), MAINP)
    CW = case_mean_wilcoxon(S, FOCUS, others)
    feas_pct = {a: float(G[G.Algorithm == a].Feasible.mean() * 100) for a in MAINP}
    mean_loss = {a: float(S[(S.Algorithm == a) & S.Qualified].Loss.mean()) for a in MAINP}
    summary["main"] = dict(friedman=FR, friedman_old_rule=dict(avg_rank=FR_old["avg_rank"], chi2=FR_old["chi2"]),
                           case_mean_wilcoxon=CW, feasible_pct=feas_pct,
                           mean_loss_pct_qualified_cases=mean_loss,
                           cases_below_half_feasible={a: int((~S[S.Algorithm == a].Qualified).sum()) for a in MAINP},
                           wtl={b: C[C.Baseline == b].Outcome.value_counts().to_dict() for b in others},
                           rank_biserial_median={b: float(C[C.Baseline == b].RB.median()) for b in others})
    log(f"  Friedman chi2={FR['chi2']:.1f} p={FR['p']:.2e}  ID F={FR['iman_davenport']:.1f}")
    log("  avg ranks: " + ", ".join(f"{LAB[a]} {v:.2f}" for a, v in sorted(FR["avg_rank"].items(), key=lambda t: t[1])))
    log("  W/T/L vs SSA-VNS: " + ", ".join(f"{LAB[b]} {wtl_str(summary['main']['wtl'][b])}" for b in others))

    # --- per-case tables (supplementary)
    for ds in ("1", "2"):
        for r, nmax in RADII.items():
            lines = []
            for n in range(2, nmax + 1):
                s = S[(S.Dataset == ds) & (S.Radius == r) & (S.Turbines == n)].set_index("Algorithm").reindex(MAINP)
                ideal = G[(G.Dataset == ds) & (G.Radius == r) & (G.Turbines == n)].Ideal.iloc[0]
                qm = s.Mean.where(s.Qualified.astype(bool)).round(6)
                best = set(qm.index[qm == qm.max()]) if qm.notna().any() else set()      # ties all bold
                cells = []
                for a in MAINP:
                    if not s.NFeas[a]:
                        cells.append(f"--$^{{0}}$"); continue
                    sd = s.SD[a]
                    c = f"{s.Mean[a]:.1f} ({0 if not np.isfinite(sd) else sd:.1f})"
                    if a in best:
                        c = f"\\textbf{{{c}}}"
                    if s.NFeas[a] < s.N[a]:
                        c += f"$^{{{int(s.NFeas[a])}}}$"
                    cells.append(c)
                lines.append(f"{n} & {ideal:.1f} & " + " & ".join(cells) + " \\\\")
            supp.append(table("table*", f"{DSN[ds]}, {r}-m farm: benchmark objective, mean (SD) over the feasible runs of 30 seed-paired runs at 6,030 objective calls. Bold: highest mean among the methods with at least 15 feasible runs. A superscript gives the number of feasible runs when fewer than 30; ``--'' means no feasible run.",
                              f"tab:res-{ds}-{r}", "c" * (2 + len(MAINP)),
                              "$N$ & Ideal & " + " & ".join(LAB[a] for a in MAINP), lines, resize=True, sep="2pt"))

    # --- W/T/L table
    lines = []
    for ds in ("1", "2"):
        for r in RADII:
            s = C[(C.Dataset == ds) & (C.Radius == r)]
            lines.append(f"{'I' if ds == '1' else 'II'} & {r} & {s.Turbines.nunique()} & " +
                         " & ".join(wtl_str(s[s.Baseline == b].Outcome.value_counts()) for b in others) + " \\\\")
    lines.append("\\midrule\nAll & & %d & " % C[CASE].drop_duplicates().shape[0] +
                 " & ".join(wtl_str(C[C.Baseline == b].Outcome.value_counts()) for b in others) + " \\\\")
    lines.append("$\\tilde r_{\\rm rb}$ & & & " + " & ".join(f"{C[C.Baseline == b].RB.median():+.2f}" for b in others) + " \\\\")
    tabs["wtl"] = table("table", "Pairwise outcome of SSA-VNS against each method: number of cases in which SSA-VNS is significantly better / not significantly different / significantly worse (two-sided Wilcoxon signed-rank test on 30 seed-paired runs, Holm-adjusted over the %d comparisons of each case, $\\alpha=0.05$; infeasible runs rank below all feasible runs). Last row: median over the cases of the rank-biserial correlation $r_{\\rm rb}$ (positive = SSA-VNS better)." % len(others),
                        "tab:wtl", "llc" + "c" * len(others),
                        "DS & $r$ (m) & Cases & " + " & ".join(LAB[b] for b in others), lines, sep="2.5pt")

    # --- Friedman + case-mean Wilcoxon table
    order = sorted(MAINP, key=lambda a: FR["avg_rank"][a])
    lines = []
    for a in order:
        if a == FOCUS:
            pz = pw = rb = dl = "--"
        else:
            pz = fmt_p(FR["p_holm_vs_focus"][a]); pw = fmt_p(CW[a]["p_holm"]); rb = f"{CW[a]['rb']:+.2f}"
            dl = f"${CW[a]['mean_dloss_pp']:+.3f}$"
        lines.append(f"{LAB[a]} & {FR['avg_rank'][a]:.2f} & {FR['sole_best_count'][a]} & {feas_pct[a]:.1f} & "
                     f"{summary['main']['cases_below_half_feasible'][a]} & {pz} & {pw} & {rb} & {dl} \\\\")
    tabs["friedman68"] = table(
        "table*", "Case-level analysis over the %d benchmark cases. Ranks of the mean feasible objective within each case (1 = best); a method with fewer than 15 feasible runs out of 30 in a case is ranked below all other methods, by its number of feasible runs (``$<$15'': number of such cases). Friedman $\\chi^2_F=%.1f$ (%d d.f.), $p=%s$; Iman--Davenport $F_F=%.1f$, $p=%s$. $p_z$: Holm-adjusted $p$ of the average-rank $z$ test against SSA-VNS (best-ranked method: %s). Because mean-rank post hoc tests depend on the pool of compared methods \\cite{Benavoli2016}, pairwise two-sided Wilcoxon signed-rank tests on the per-case mean wake losses (\\%%) are also given: $p_W$ (Holm-adjusted over the %d comparisons), matched-pairs rank-biserial correlation $r_{\\rm rb}$ (positive = SSA-VNS better) and mean difference $\\overline{\\Delta L}$ of the wake loss (percentage points, SSA-VNS minus method, over the cases where both methods have at least 15 feasible runs). ``Best'': cases in which the method alone ranks first; ``Feas.'': percentage of feasible runs."
        % (FR["n_cases"], FR["chi2"], len(MAINP) - 1, fmt_p(FR["p"]).strip("$"), FR["iman_davenport"],
           fmt_p(FR["iman_davenport_p"]).strip("$"), LAB[FR["best_ranked"]], len(others)),
        "tab:friedman68", "lcccccccc",
        "Method & Avg.\\ rank & Best & Feas.\\ (\\%) & $<$15 & $p_z$ & $p_W$ & $r_{\\rm rb}$ & $\\overline{\\Delta L}$ (pp)",
        lines)

    # --- figures: average ranks
    fig, ax = plt.subplots(figsize=(3.4, 2.1))
    for pos, a in enumerate(order[::-1]):
        v = FR["avg_rank"][a]
        ax.plot([1, v], [pos, pos], color=GRID, lw=2, zorder=1)
        ax.scatter(v, pos, s=90 if MRK[a] == "*" else 36, color=COL[a], marker=MRK[a], zorder=3, edgecolor="white", lw=0.6)
        ax.text(v + 0.08, pos, f"{v:.2f}", va="center", fontsize=7, color=INK)
    ax.set_yticks(range(len(order))); ax.set_yticklabels([LAB[a] for a in order[::-1]])
    ax.set_xlim(1, len(MAINP) + 0.3); ax.set_xlabel(f"Average rank over {FR['n_cases']} cases (1 = best)")
    ax.grid(axis="y", visible=False)
    FG.save(fig, "avg_ranks")

    # wake loss vs N (only case means of methods with >= 15 feasible runs are drawn)
    for what, name, ylab in (("Loss", "wakeloss_vs_n", "Mean wake loss (% of ideal)"),
                             ("Feas", "feasibility_vs_n", "Feasible runs (%)")):
        fig, axes = plt.subplots(2, 3, figsize=(7.1, 4.0 if what == "Loss" else 3.6), sharey="row" if what == "Loss" else True)
        for i, ds in enumerate(("1", "2")):
            for j, (r, nmax) in enumerate(RADII.items()):
                ax = axes[i, j]
                s = S[(S.Dataset == ds) & (S.Radius == r)]
                for a in MAINP:
                    m = s[s.Algorithm == a].sort_values("Turbines")
                    y = m.Loss.where(m.Qualified) if what == "Loss" else 100 * m.NFeas / m.N
                    ax.plot(m.Turbines, y, color=COL[a], ls=LS[a], marker=MRK[a], ms=6 if MRK[a] == "*" else 3.5)
                if what == "Feas":
                    ax.set_ylim(-5, 105)
                ax.set_title(f"{DSN[ds]}, $r$ = {r} m", fontsize=8, color=INK)
                ax.set_xticks(range(2, nmax + 1, 2 if nmax > 10 else 1))
                if i == 1: ax.set_xlabel("Number of turbines $N$")
                if j == 0: ax.set_ylabel(ylab)
        legend_row(fig, MAINP)
        fig.tight_layout()
        FG.save(fig, name)

    # convergence at mid and max N
    for tag, lv in (("mid", MID), ("max", RADII)):
        fig, axes = plt.subplots(2, 3, figsize=(7.1, 4.0))
        for i, ds in enumerate(("1", "2")):
            for j, r in enumerate(RADII):
                ax = axes[i, j]; n = lv[r]
                conv_panel(ax, G[(G.Dataset == ds) & (G.Radius == r) & (G.Turbines == n)], MAINP)
                log_axis(ax)
                ax.set_title(f"{DSN[ds]}, $r$ = {r} m, $N$ = {n}", fontsize=8, color=INK)
                if i == 1: ax.set_xlabel("Objective-function calls")
                if j == 0: ax.set_ylabel("Best wake loss (% of ideal)")
        legend_row(fig, MAINP)
        fig.tight_layout()
        FG.save(fig, f"convergence_{tag}")

    # box plots (mid and max N)
    for tag, lv in (("mid", MID), ("max", RADII)):
        fig, axes = plt.subplots(2, 3, figsize=(7.1, 3.8))
        for i, ds in enumerate(("1", "2")):
            for j, r in enumerate(RADII):
                ax = axes[i, j]; n = lv[r]
                sub = G[(G.Dataset == ds) & (G.Radius == r) & (G.Turbines == n) & G.Feasible]
                data = [sub[sub.Algorithm == a].LossPct.values for a in MAINP]
                bp = ax.boxplot(data, widths=0.6, patch_artist=True, showfliers=True,
                                flierprops=dict(marker=".", ms=3, mec=MUTED), medianprops=dict(color=INK, lw=1),
                                whiskerprops=dict(color=MUTED), capprops=dict(color=MUTED))
                for patch, a in zip(bp["boxes"], MAINP):
                    patch.set_facecolor(COL[a]); patch.set_alpha(0.55); patch.set_edgecolor(COL[a])
                cnt = [len(d) for d in data]
                ax.set_xticks(range(1, len(MAINP) + 1))
                ax.set_xticklabels([LAB[a] + ("" if c == 30 else f"\n({c})") for a, c in zip(MAINP, cnt)], rotation=40, fontsize=6)
                ax.set_title(f"{DSN[ds]}, $r$ = {r} m, $N$ = {n}", fontsize=8, color=INK)
                if j == 0: ax.set_ylabel("Wake loss (% of ideal)")
                ax.grid(axis="x", visible=False)
        fig.tight_layout()
        FG.save(fig, f"boxplots_{tag}")

    # best layouts, largest N
    fig, axes = plt.subplots(2, 3, figsize=(7.1, 4.8))
    best_rows = []
    for i, ds in enumerate(("1", "2")):
        for j, (r, n) in enumerate(RADII.items()):
            ax = axes[i, j]
            sub = G[(G.Dataset == ds) & (G.Radius == r) & (G.Turbines == n) & G.Feasible]
            top = sub.sort_values("Objective").iloc[-1]
            fs = sub[sub.Algorithm == FOCUS]
            lx = fs.sort_values("Objective").iloc[-1] if len(fs) else None
            for row, mk, colr in ((top, "s", INK), (lx, "o", COL[FOCUS])):
                if row is None:
                    continue
                xy = coords(row.Coordinates)
                ax.scatter(xy[:, 0], xy[:, 1], marker=mk, s=22 if mk == "s" else 12,
                           facecolor="none" if mk == "s" else colr, edgecolor=colr, lw=1, zorder=3)
                best_rows.append(dict(Dataset=ds, Radius=r, Turbines=n, Algorithm=row.Algorithm, Seed=row.Seed,
                                      Objective=row.Objective, WakeLoss=row.WakeLoss, LossPct=row.LossPct,
                                      Coordinates=row.Coordinates))
            t = np.linspace(0, 2 * np.pi, 200)
            ax.plot(r * np.cos(t), r * np.sin(t), color=MUTED, lw=0.8)
            ax.set_aspect("equal"); ax.set_xlim(-1.08 * r, 1.08 * r); ax.set_ylim(-1.08 * r, 1.08 * r)
            gap = "" if lx is None else f" (best SSA-VNS: +{100 * (top.Objective - lx.Objective) / top.Ideal:.2f} pp loss)"
            ax.set_title(f"{DSN[ds]}, $r$ = {r} m, $N$ = {n}\nbest of all methods: {LAB[top.Algorithm]}{gap if top.Algorithm != FOCUS else ''}",
                         fontsize=8, color=INK)
            ax.tick_params(labelsize=6)
    hl = [plt.Line2D([], [], ls="none", marker="s", ms=6, mfc="none", mec=INK, mew=1, label="best layout of all methods"),
          plt.Line2D([], [], ls="none", marker="o", ms=4.5, color=COL[FOCUS], label="best SSA-VNS layout")]
    fig.legend(handles=hl, loc="lower center", ncol=2, bbox_to_anchor=(0.5, 1.0), fontsize=7.5)
    fig.tight_layout()
    FG.save(fig, "layouts_max")
    BL = pd.DataFrame(best_rows).drop_duplicates()
    BL.to_csv(OUT("mpce_best_layouts_maxN.csv"), index=False)
    summary["main"]["best_layout_maxN"] = [
        dict(case=f"{d}-{r}-{n}", best_method=g.sort_values("Objective").Algorithm.iloc[-1],
             focus_gap_pp=float((g.Objective.max() - g[g.Algorithm == FOCUS].Objective.max()) / g.Ideal.iloc[0] * 100)
             if (g.Algorithm == FOCUS).any() else None)
        for (d, r, n), g in BL.assign(Ideal=BL.Objective + BL.WakeLoss).groupby(CASE)]

    # computational cost
    rt = (G.Seconds / G.Calls * 1000).groupby(G.Algorithm).mean().reindex(MAINP)
    rn = G.groupby(["Algorithm", "Turbines"]).Seconds.mean().unstack().reindex(MAINP)
    src = G.groupby("Algorithm").Source.agg(lambda s: ",".join(sorted(set(s))))
    summary["cost"] = dict(ms_per_call=rt.round(4).to_dict(),
                           sec_per_run_range={a: [float(rn.loc[a].min()), float(rn.loc[a].max())] for a in MAINP},
                           source=src.to_dict(),
                           note="wall-clock times come from the machine that produced each file; runs from different "
                                "files may not be directly comparable")
    lines = [f"{LAB[a]} & {rt[a]:.3f} & {rn.loc[a].min():.1f}--{rn.loc[a].max():.1f} \\\\" for a in MAINP]
    tabs["cost"] = table("table", "Computational cost at 6,030 objective calls: mean wall-clock time per objective call (including the optimizer overhead) and range over $N$ of the mean time per run (single core).",
                         "tab:cost", "lcc", "Method & ms per call & s per run", lines)

    # =========================================================== 2. ablation
    log("\n[2] Ablation")
    GA = R6[R6.Dataset.isin(["1", "2"])]
    ablp = [a for a in ABL if a in set(GA.Algorithm)]
    if "RSVNS" not in ablp:
        log("  RS-VNS runs missing: ablation without the RS-VNS variant/contrast")
    GB = GA[GA.Algorithm.isin(ablp)]
    CONTR = [("SSABV", "SSA", "VNS phase (after SSA)"),
             ("SSABV", "BVNS", "SSA start vs.\\ best initial point"),
             ("SSABV", "RSVNS", "swarm phase vs.\\ random sampling (equal calls)"),
             ("SSABV", "LXBV", "Laplace step inside the hybrid"),
             ("LXSSA", "SSA", "Laplace step alone (no VNS)")]
    CONTR = [c for c in CONTR if c[0] in ablp and c[1] in ablp]
    arows = []
    for (ds, r, n), sub in GB.groupby(CASE):
        piv = sub.assign(S=goodness(sub)).pivot_table(index="Seed", columns="Algorithm", values="S")
        lm = sub[sub.Feasible].groupby("Algorithm").LossPct.mean()
        rr = []
        for a, b, _ in CONTR:
            p, rb, _ = wil((piv[a] - piv[b]).dropna().values)
            rr.append(dict(Dataset=ds, Radius=r, Turbines=n, A=a, B=b, P=p, RB=rb, DLoss=lm.get(a, np.nan) - lm.get(b, np.nan)))
        for x, h in zip(rr, holm([x["P"] for x in rr])):
            x["PHolm"] = h
        arows += rr
    A = pd.DataFrame(arows)
    A["Outcome"] = np.where(A.PHolm < 0.05, np.where(A.RB > 0, "W", "L"), "T")
    A.to_csv(OUT("mpce_ablation_tests.csv"), index=False)
    SA = case_stats(GB, ablp)
    FA = friedman_block(rank_matrix(SA, ablp), ablp)
    comp = {"SSABV": ("SSA", "VNS"), "LXBV": ("LX-SSA", "VNS"), "RSVNS": ("random sampling", "VNS"),
            "BVNS": ("--", "VNS"), "SSA": ("SSA", "--"), "LXSSA": ("LX-SSA", "--")}
    afeas = {a: float(GB[GB.Algorithm == a].Feasible.mean() * 100) for a in ablp}
    rows1 = [f"{LAB[a]} & {comp[a][0]} & {comp[a][1]} & {FA['avg_rank'][a]:.2f} & {afeas[a]:.1f} \\\\"
             for a in sorted(ablp, key=lambda a: FA["avg_rank"][a])]
    rows2 = []
    for a, b, what in CONTR:
        x = A[(A.A == a) & (A.B == b)]
        rows2.append(f"{LAB[a]} vs.\\ {LAB[b]} & {what} & {wtl_str(x.Outcome.value_counts())} & ${x.DLoss.mean():+.3f}$ \\\\")
    tabs["ablation"] = (r"""\begin{table}[!t]
\centering
\caption{Component ablation over the 68 benchmark cases (30 seed-paired runs, 6,030 calls; RS-VNS spends the first half of the budget on random layouts and refines the best one by the same VNS). Top: phase-1 and phase-2 components, average rank among the %d variants (ranking rule of Table~\ref{tab:friedman68}; Friedman $\chi^2_F=%.1f$, %d d.f., $p=%s$) and percentage of feasible runs. Bottom: planned contrasts, number of cases in which the first variant is significantly better / not different / worse (Wilcoxon signed-rank, Holm-adjusted over the %d contrasts of each case, $\alpha=0.05$), and mean difference of the wake loss $\overline{\Delta L}$ (percentage points; negative = first variant better).}
\label{tab:ablation}
\scriptsize\setlength{\tabcolsep}{2.5pt}
\begin{tabular}{lcccc}
\toprule
Variant & Phase 1 & Phase 2 & Avg.\ rank & Feas.\ (\%%) \\
\midrule
""" % (len(ablp), FA["chi2"], len(ablp) - 1, fmt_p(FA["p"]).strip("$"), len(CONTR)) + "\n".join(rows1) + r"""
\bottomrule
\end{tabular}

\smallskip
\begin{tabular}{l>{\raggedright\arraybackslash}p{2.8cm}cc}
\toprule
Contrast & Isolates & W/T/L & $\overline{\Delta L}$ \\
\midrule
""" + "\n".join(rows2) + r"""
\bottomrule
\end{tabular}
\end{table}
""")
    summary["ablation"] = dict(variants=ablp, friedman=FA, feasible_pct=afeas,
                               contrasts={f"{a}-{b}": dict(A[(A.A == a) & (A.B == b)].Outcome.value_counts().to_dict(),
                                                            dloss_pp=float(A[(A.A == a) & (A.B == b)].DLoss.mean()),
                                                            isolates=w.replace("\\", ""))
                                          for a, b, w in CONTR})
    log("  avg ranks: " + ", ".join(f"{LAB[a]} {FA['avg_rank'][a]:.3f}" for a in ablp))
    for k, v in summary["ablation"]["contrasts"].items():
        log(f"  {k}: {wtl_str(v)}  dL={v['dloss_pp']:+.3f}")

    # phase-2 statistic (VNS phase after the swarm): wake loss left at the switch removed by VNS
    ph2 = {}
    for alg, sw in (("SSABV", 3030), ("LXBV", 3030)):
        Hy = GA[(GA.Algorithm == alg) & GA.Feasible]
        if not len(Hy):
            continue
        gains, made_feasible = [], 0
        idx = sw // 30 - 1
        for ideal, obj, cv in zip(Hy.Ideal.values, Hy.Objective.values, Hy.Curve.values):
            c = cv.split(";")[idx]
            at = float(c)
            if not np.isfinite(at):
                made_feasible += 1; continue
            l1, l2 = ideal - at, ideal - obj
            if l1 > 1e-9:
                gains.append(100 * (l1 - l2) / l1)
        nall = int((GA.Algorithm == alg).sum())
        ph2[alg] = dict(mean=float(np.mean(gains)), median=float(np.median(gains)), n_runs_with_loss_at_switch=len(gains),
                        infeasible_at_switch_made_feasible=made_feasible,
                        final_infeasible=int(nall - len(Hy)), switch_call=sw,
                        pct_runs_improved=float(100 * np.mean(np.array(gains) > 1e-9)))
    summary["ablation"]["phase2_loss_reduction_pct"] = ph2
    if "SSABV" in ph2:
        log(f"  phase 2 of SSA-VNS removes {ph2['SSABV']['mean']:.2f}% (median {ph2['SSABV']['median']:.2f}%) of the loss "
            f"left at the switch; {ph2['SSABV']['infeasible_at_switch_made_feasible']} runs infeasible at the switch made feasible")

    # reproduction check against the parent of commit a9779a3 (old rule, five variants, six contrasts)
    OLDV = ["SSABV", "LXBV", "BVNS", "LXSSA", "SSA"]
    if all(a in set(GA.Algorithm) for a in OLDV):
        Go = GA[GA.Algorithm.isin(OLDV)]
        So = case_stats(Go, OLDV, rank_old)
        Fo = friedman_block(rank_matrix(So, OLDV), OLDV)
        OC = [("LXBV", "LXSSA"), ("LXBV", "BVNS"), ("LXBV", "SSABV"), ("SSABV", "SSA"), ("SSABV", "BVNS"), ("LXSSA", "SSA")]
        orows = []
        for _, sub in Go.groupby(CASE):
            piv = sub.assign(S=goodness(sub)).pivot_table(index="Seed", columns="Algorithm", values="S")
            lm = sub[sub.Feasible].groupby("Algorithm").LossPct.mean()
            rr = [dict(A=a, B=b, P=wil((piv[a] - piv[b]).values)[0], RB=wil((piv[a] - piv[b]).values)[1],
                       DL=lm.get(a, np.nan) - lm.get(b, np.nan)) for a, b in OC]
            for x, h in zip(rr, holm([x["P"] for x in rr])):
                x["O"] = ("W" if x["RB"] > 0 else "L") if h < 0.05 else "T"
            orows += rr
        O = pd.DataFrame(orows)
        x = O[(O.A == "LXBV") & (O.B == "SSABV")]
        oc = x.O.value_counts()
        rep = dict(ssabv_avg_rank=Fo["avg_rank"]["SSABV"], expected_ssabv_avg_rank=1.625,
                   lxbv_vs_ssabv_wtl=wtl_str(oc), expected_lxbv_vs_ssabv_wtl="0/66/2",
                   lxbv_minus_ssabv_dloss=float(x.DL.mean()), expected_dloss=0.06306204787105701,
                   phase2_ssabv=ph2.get("SSABV"), expected_phase2_ssabv=dict(mean=36.66372550902314, median=29.194081556313783,
                                                                             infeasible_at_switch=79),
                   old_avg_rank=Fo["avg_rank"])
        rep["pass"] = bool(abs(rep["ssabv_avg_rank"] - 1.625) < 1e-9 and rep["lxbv_vs_ssabv_wtl"] == "0/66/2"
                           and abs(rep["lxbv_minus_ssabv_dloss"] - 0.06306204787105701) < 1e-9
                           and abs(ph2["SSABV"]["mean"] - 36.66372550902314) < 1e-6
                           and ph2["SSABV"]["infeasible_at_switch_made_feasible"] == 79)
        summary["reproduction_check_a9779a3_parent"] = rep
        log(f"  REPRODUCTION CHECK (old rule, 5 variants, Holm over 6 contrasts): SSABV avg rank {rep['ssabv_avg_rank']:.4f} "
            f"(expected 1.625), LXBV vs SSABV {rep['lxbv_vs_ssabv_wtl']} (expected 0/66/2), dL {rep['lxbv_minus_ssabv_dloss']:+.4f} "
            f"(expected +0.0631), phase-2 mean {ph2['SSABV']['mean']:.3f} (expected 36.664) -> {'PASS' if rep['pass'] else 'FAIL'}")

    # ablation convergence (largest N)
    fig, axes = plt.subplots(2, 3, figsize=(7.1, 4.0))
    for i, ds in enumerate(("1", "2")):
        for j, (r, n) in enumerate(RADII.items()):
            ax = axes[i, j]
            conv_panel(ax, GB[(GB.Dataset == ds) & (GB.Radius == r) & (GB.Turbines == n)], ablp, band=False)
            ax.axvline(3030, color=MUTED, lw=0.7, ls=(0, (1, 2)))
            log_axis(ax)
            ax.set_title(f"{DSN[ds]}, $r$ = {r} m, $N$ = {n}", fontsize=8, color=INK)
            if i == 1: ax.set_xlabel("Objective-function calls")
            if j == 0: ax.set_ylabel("Median best wake loss (%)")
    legend_row(fig, ablp)
    fig.tight_layout()
    FG.save(fig, "ablation_convergence")

    # =========================================================== 3. budget split
    log("\n[3] Budget split")
    Sp = R6[R6.Algorithm.isin(["SSABV25", "SSABV", "SSABV75"])]
    Sp = Sp[[(d, r, n) in SPLITCASES for d, r, n in zip(Sp.Dataset, Sp.Radius, Sp.Turbines)]]
    if {"SSABV25", "SSABV75"} <= set(Sp.Algorithm):
        lines, srows = [], []
        for (ds, r, n), s in Sp.groupby(CASE):
            if s.Algorithm.nunique() < 3:
                continue
            piv = s.assign(S=goodness(s)).pivot_table(index="Seed", columns="Algorithm", values="S")
            fm = s[s.Feasible].groupby("Algorithm").Objective.mean(); fc = s.groupby("Algorithm").Feasible.sum()
            lmn = s[s.Feasible].groupby("Algorithm").LossPct.mean()
            ps, rbs = [], []
            for b in ("SSABV25", "SSABV75"):
                p, rb, _ = wil((piv["SSABV"] - piv[b]).dropna().values); ps.append(p); rbs.append(rb)
            ph = holm(ps)
            trio = ("SSABV25", "SSABV", "SSABV75")
            cells = []
            for a in trio:
                v = fm.get(a, np.nan)
                c = "--" if not np.isfinite(v) else f"{v:.1f}"
                if fc.get(a, 0) < 30:
                    c += f"$^{{{int(fc.get(a, 0))}}}$"
                cells.append(c)
            best = max(trio, key=lambda a: fm.get(a, -np.inf) if fc.get(a, 0) >= 15 else -np.inf)
            srows.append(dict(case=f"{ds}-{r}-{n}", best=best, p25_holm=float(ph[0]), p75_holm=float(ph[1]),
                              rb25=rbs[0], rb75=rbs[1], loss={a: float(lmn.get(a, np.nan)) for a in trio}))
            lines.append(f"{'I' if ds == '1' else 'II'} & {r} & {n} & " + " & ".join(cells) + f" & {fmt_p(ph[0])} & {fmt_p(ph[1])} \\\\")
        tabs["split"] = table("table", "Sensitivity of SSA-VNS to the budget split between the SSA and VNS phases (25\\%, 50\\% and 75\\% of the 6,030 calls for SSA): mean benchmark objective of the feasible runs (superscript: feasible runs when fewer than 30) and Holm-adjusted Wilcoxon signed-rank $p$ of the 50\\% split against the 25\\% and 75\\% splits (30 seed-paired runs).",
                              "tab:split", "cccccccc", "DS & $r$ & $N$ & 25\\% & 50\\% & 75\\% & $p$ (25) & $p$ (75)", lines, sep="2.5pt")
        SR = pd.DataFrame(srows)
        summary["split"] = dict(cases=srows, best_count=SR.best.value_counts().to_dict(),
                                sig_vs25=int((SR.p25_holm < 0.05).sum()), sig_vs75=int((SR.p75_holm < 0.05).sum()),
                                sig_vs25_50better=int(((SR.p25_holm < 0.05) & (SR.rb25 > 0)).sum()),
                                sig_vs75_50better=int(((SR.p75_holm < 0.05) & (SR.rb75 > 0)).sum()))
        log(f"  best split counts {summary['split']['best_count']}; significant vs 25%: {summary['split']['sig_vs25']}, vs 75%: {summary['split']['sig_vs75']}")
    else:
        log("  SKIPPED: mpce_ssasplit (SSABV25/SSABV75) not available")
        summary["split"] = None

    # =========================================================== 4. Horns Rev 16
    log("\n[4] Horns Rev 1, 16 turbines")
    try:
        import hornsrev_model as hr
        inst = float(hr.aep_gwh(hr.site(16)[0]))
    except Exception as e:                                  # pragma: no cover
        hr, inst = None, INSTALLED_HR16
        log(f"  hornsrev_model unavailable ({e}); installed AEP = {INSTALLED_HR16}")
    H = R6[(R6.Dataset == "HR") & (R6.Turbines == 16)]
    hm = MAINP + (["RSVNS"] if "RSVNS" in set(H.Algorithm) else [])
    H = H[H.Algorithm.isin(hm)]
    hm = [a for a in hm if a in set(H.Algorithm)]
    if len(H):
        ideal = float(H.Ideal.iloc[0])
        tt = {t["Baseline"]: t for t in paired_vs(H, FOCUS, [a for a in hm if a != FOCUS])}
        piv = H.assign(S=goodness(H)).pivot_table(index="Seed", columns="Algorithm", values="S").dropna()
        pf = friedmanchisquare(*[piv[a] for a in hm])[1] if len(piv) > 1 else np.nan
        SH = case_stats(H, hm)
        SH = SH.set_index("Algorithm")
        lines = [f"Installed layout & {inst:.2f} & {100 * (1 - inst / ideal):.2f} & -- & -- & -- \\\\", "\\midrule"]
        hrsum = {}
        for a in hm:
            s = SH.loc[a]
            if s.NFeas:
                c = [f"{s.Mean:.2f} ({0 if not np.isfinite(s.SD) else s.SD:.2f})", f"{s.Loss:.2f}", f"{int(s.NFeas)}/{int(s.N)}"]
            else:
                c = ["--", "--", f"0/{int(s.N)}"]
            c += ["--", "--"] if a == FOCUS else [fmt_p(tt[a]["PHolm"]), f"{tt[a]['RB']:+.2f}"]
            lines.append(f"{LAB[a]} & " + " & ".join(c) + " \\\\")
            hrsum[a] = dict(mean=float(s.Mean), sd=float(s.SD), best=float(s.Best), loss_pct=float(s.Loss),
                            feasible=int(s.NFeas), runs=int(s.N), rank=float(s.Rank),
                            gain_vs_installed_pct=float(100 * (s.Mean / inst - 1)),
                            p_holm=None if a == FOCUS else tt[a]["PHolm"], rb=None if a == FOCUS else tt[a]["RB"])
        tabs["hr16"] = table("table", "Horns Rev~1 site case, 16-turbine block (Vestas V80 power and thrust curves, measured 12-sector wind climate, installed outline; Jensen wake $k=0.04$; minimum spacing $4D=320$~m; wake-free AEP %.2f~GWh/yr): AEP in GWh/yr, mean (SD) over the feasible runs of 30 seed-paired runs at 6,030 objective calls, mean wake loss (\\%%), feasible runs, Holm-adjusted Wilcoxon signed-rank $p$ of SSA-VNS vs.\\ each method and rank-biserial $r_{\\rm rb}$ (positive = SSA-VNS better); run-level Friedman $p=%s$." % (ideal, fmt_p(pf).strip("$")),
                             "tab:hr-site", "lccccc", "Layout / method & AEP & Loss (\\%) & Feas. & $p_{\\rm Holm}$ & $r_{\\rm rb}$", lines, size="\\footnotesize", sep="3pt")
        summary["hr16"] = dict(installed_aep=inst, ideal_aep=ideal, installed_loss_pct=100 * (1 - inst / ideal),
                               friedman_p=float(pf), methods=hrsum,
                               sources=H.groupby("Algorithm").Source.first().to_dict())
        log("  mean AEP: " + ", ".join(f"{LAB[a]} {hrsum[a]['mean']:.2f}({hrsum[a]['feasible']})" for a in hm) + f"; installed {inst:.2f}")
        # convergence + layouts
        fig, axes = plt.subplots(1, 2, figsize=(7.1, 2.8), gridspec_kw=dict(width_ratios=[1.25, 1]))
        ax = axes[0]
        conv_panel(ax, H, hm, loss=False, band=False)
        ax.axhline(inst, color=INK, lw=0.9, ls=(0, (1, 1)))
        ax.text(1000, inst, "installed layout", va="bottom", fontsize=6.5, color=INK)
        ax.set_title("Horns Rev 1, 16 turbines: median best feasible AEP", fontsize=8, color=INK)
        ax.set_xlabel("Objective-function calls"); ax.set_ylabel("AEP (GWh/yr)")
        lo = np.nanmin([SH.Mean[a] for a in hm if SH.NFeas[a]] + [inst])
        ax.set_ylim(bottom=lo - 2.0)
        ax = axes[1]
        if hr is not None:
            xy0, poly = hr.site(16)
            pp = np.vstack([poly, poly[:1]])
            ax.plot(pp[:, 0], pp[:, 1], color=MUTED, lw=0.8)
            ax.scatter(xy0[:, 0], xy0[:, 1], marker="x", s=14, color=INK, lw=0.8, label=f"installed ({inst:.2f})", zorder=3)
        fz = H[H.Feasible]
        top = fz.sort_values("Objective").iloc[-1]
        ff = fz[fz.Algorithm == FOCUS]
        if len(ff) and top.Algorithm != FOCUS:
            b = ff.sort_values("Objective").iloc[-1]; xy = coords(b.Coordinates)
            ax.scatter(xy[:, 0], xy[:, 1], marker="o", s=12, color=COL[FOCUS], label=f"best SSA-VNS ({b.Objective:.2f})", zorder=4)
        xy = coords(top.Coordinates)
        ax.scatter(xy[:, 0], xy[:, 1], marker="s" if top.Algorithm != FOCUS else "o", s=22 if top.Algorithm != FOCUS else 12,
                   facecolor="none" if top.Algorithm != FOCUS else COL[FOCUS], edgecolor=INK if top.Algorithm != FOCUS else COL[FOCUS],
                   lw=1, label=f"best of all methods: {LAB[top.Algorithm]} ({top.Objective:.2f})", zorder=3)
        summary["hr16"]["best_overall"] = dict(method=top.Algorithm, aep=float(top.Objective), seed=int(top.Seed))
        ax.set_aspect("equal"); ax.tick_params(labelsize=6)
        ax.set_title("Layouts (m); AEP in GWh/yr", fontsize=8, color=INK)
        ax.legend(fontsize=6, loc="upper left", bbox_to_anchor=(0, -0.1), ncol=1)
        legend_row(fig, hm, ncol=len(hm))
        fig.tight_layout()
        FG.save(fig, "hr16")
    else:
        log("  SKIPPED: no Horns Rev data")

    # =========================================================== 5. feasible vs random initialization
    log("\n[5] Feasible initialization")
    FE = ALL[(ALL.Init == "feasible") & (ALL.Budget == 6030)]
    cases7 = LARGE + [HR16]
    if len(FE):
        fm_ = [a for a in M9 if a in set(FE.Algorithm)]
        rows, det = [], []
        for c in cases7:
            f = FE[(FE.Dataset == c[0]) & (FE.Radius == c[1]) & (FE.Turbines == c[2])]
            rr = R6[(R6.Dataset == c[0]) & (R6.Radius == c[1]) & (R6.Turbines == c[2])]
            if not len(f):
                continue
            crow = []
            for a in fm_:
                fa, ra = f[f.Algorithm == a], rr[rr.Algorithm == a]
                if not len(fa) or not len(ra):
                    continue
                pf_ = pd.Series(goodness(fa), index=fa.Seed.values); pr_ = pd.Series(goodness(ra), index=ra.Seed.values)
                d = (pf_ - pr_).dropna()
                p, rb, _ = wil(d.values)
                crow.append(dict(case=f"{c[0]}-{c[1]}-{c[2]}", Algorithm=a, P=p, RB=rb,
                                 LossR=ra[ra.Feasible].LossPct.mean(), LossF=fa[fa.Feasible].LossPct.mean(),
                                 AEPR=ra[ra.Feasible].Objective.mean() if c[0] == "HR" else np.nan,
                                 AEPF=fa[fa.Feasible].Objective.mean() if c[0] == "HR" else np.nan,
                                 FeasR=int(ra.Feasible.sum()), FeasF=int(fa.Feasible.sum()), NR=len(ra), NF=len(fa),
                                 RandomSource=ra.Source.iloc[0]))
            for x, h in zip(crow, holm([x["P"] for x in crow])):
                x["PHolm"] = h
                x["Outcome"] = ("W" if x["RB"] > 0 else "L") if h < 0.05 else "T"
            det += crow
        Dt = pd.DataFrame(det)
        Dt.to_csv(OUT("mpce_feasinit_tests.csv"), index=False)
        lines = []
        for a in fm_:
            x = Dt[Dt.Algorithm == a]
            if not len(x):
                continue
            fb = "$^\\dagger$" if x.RandomSource.str.startswith("FALLBACK").any() else ""
            lines.append(f"{LAB[a]}{fb} & {x.LossR.mean():.3f} & {x.LossF.mean():.3f} & {100 * x.FeasR.sum() / x.NR.sum():.1f} & "
                         f"{100 * x.FeasF.sum() / x.NF.sum():.1f} & {wtl_str(x.Outcome.value_counts())} \\\\")
        tabs["feasinit"] = table("table", "Feasibility-preserving versus uniform random initialization at 6,030 calls on the six largest benchmark cases and the Horns Rev 16-turbine block (%d cases, 30 seed-paired runs each): mean wake loss (\\%%) of the feasible runs averaged over the cases, percentage of feasible runs, and number of cases in which feasible initialization is significantly better / not different / worse (Wilcoxon signed-rank, Holm-adjusted over the methods of each case). Per-case values: Supplementary Table~\\ref{tab:feasinit-cases}.%s" % (Dt.case.nunique(), " $^\\dagger$: random-initialization runs are a development fallback." if "dagger" in "".join(lines) else ""),
                                 "tab:feasinit", "lccccc",
                                 "Method & \\multicolumn{2}{c}{Loss (\\%)} & \\multicolumn{2}{c}{Feas. (\\%)} & W/T/L \\\\\n & random & feasible & random & feasible &", lines)
        dl = []
        for case, x in Dt.groupby("case", sort=False):
            dl.append(f"\\multicolumn{{7}}{{l}}{{\\emph{{{case_name(case)}}}}} \\\\")
            for _, y in x.iterrows():
                isHR = case.startswith("HR")
                vr = f"{y.AEPR:.2f}" if isHR else f"{y.LossR:.3f}"; vf = f"{y.AEPF:.2f}" if isHR else f"{y.LossF:.3f}"
                dl.append(f"{LAB[y.Algorithm]} & {vr} & {vf} & {y.FeasR}/{y.NR} & {y.FeasF}/{y.NF} & {fmt_p(y.PHolm)} & {y.RB:+.2f} \\\\")
        supp.append(table("table", "Feasible versus random initialization per case (6,030 calls): mean wake loss (\\%; Horns Rev: mean AEP in GWh/yr) of the feasible runs, feasible runs, Holm-adjusted Wilcoxon signed-rank $p$ (over the methods of the case) and rank-biserial $r_{\\rm rb}$ (positive = feasible initialization better).",
                          "tab:feasinit-cases", "lcccccc", "Method & Random & Feasible & Feas. (R) & Feas. (F) & $p_{\\rm Holm}$ & $r_{\\rm rb}$", dl, pos="p"))
        summary["feasinit"] = dict(cases=sorted(Dt.case.unique()), per_method={
            a: dict(loss_random=float(x.LossR.mean()), loss_feasible=float(x.LossF.mean()),
                    feas_random_pct=float(100 * x.FeasR.sum() / x.NR.sum()), feas_feasible_pct=float(100 * x.FeasF.sum() / x.NF.sum()),
                    wtl=x.Outcome.value_counts().to_dict(), random_source=sorted(set(x.RandomSource)))
            for a, x in Dt.groupby("Algorithm")})
        # ranks under feasible init
        SF = case_stats(FE, fm_)
        summary["feasinit"]["avg_rank_feasible_init"] = SF.groupby("Algorithm").Rank.mean().to_dict()
        log("  per method (loss random -> feasible): " + ", ".join(
            f"{LAB[a]} {v['loss_random']:.3f}->{v['loss_feasible']:.3f}" for a, v in summary["feasinit"]["per_method"].items()))
    else:
        log("  SKIPPED: mpce_feas not available")
        summary["feasinit"] = None

    # =========================================================== 6. budget scaling
    log("\n[6] Budget scaling")
    BD = ALL[(ALL.Init == "random") & ALL.Algorithm.isin(M9)]
    BD = BD[[(d, r, n) in cases7 for d, r, n in zip(BD.Dataset, BD.Radius, BD.Turbines)]]
    have_b = [b for b in BUDGETS if b != 6030 and (BD.Budget == b).any()]
    if have_b:
        bl = [6030] + have_b
        BD = BD[BD.Budget.isin(bl)]
        bsum, rows_c = {}, []
        for b in bl:
            x = BD[BD.Budget == b]
            ms = [a for a in M9 if a in set(x.Algorithm)]
            Sb = case_stats(x, ms)
            Sb["Budget"] = b
            rows_c.append(Sb)
            Rb = rank_matrix(Sb, ms)
            avg = Rb.mean(axis=0).to_dict()
            grid = Sb[Sb.Dataset != "HR"]
            bsum[b] = dict(n_cases=len(Rb), avg_rank=avg, best=min(avg, key=avg.get),
                           ssabv_rank_position=int(1 + sorted(avg.values()).index(avg[FOCUS])) if FOCUS in avg else None,
                           mean_loss_grid={a: float(grid[(grid.Algorithm == a) & grid.Qualified].Loss.mean()) for a in ms},
                           feas_pct={a: float(100 * x[x.Algorithm == a].Feasible.mean()) for a in ms},
                           hr16_aep={a: float(Sb[(Sb.Dataset == "HR") & (Sb.Algorithm == a)].Mean.iloc[0])
                                     for a in ms if ((Sb.Dataset == "HR") & (Sb.Algorithm == a)).any()},
                           sources={a: sorted(set(x[x.Algorithm == a].Source)) for a in ms})
            # SSA-VNS vs MS-SLSQP per case
            if FOCUS in ms and "SLSQP" in ms:
                vs = {}
                for (d, r, n), sub in x.groupby(CASE):
                    t = paired_vs(sub, FOCUS, ["SLSQP"])
                    ss = Sb[(Sb.Dataset == d) & (Sb.Radius == r) & (Sb.Turbines == n)].set_index("Algorithm")
                    vs[f"{d}-{r}-{n}"] = dict(p=t[0]["P"] if t else None, rb=t[0]["RB"] if t else None,
                                              ssabv=float(ss.Mean.get(FOCUS, np.nan)), slsqp=float(ss.Mean.get("SLSQP", np.nan)),
                                              ssabv_loss=float(ss.Loss.get(FOCUS, np.nan)), slsqp_loss=float(ss.Loss.get("SLSQP", np.nan)))
                bsum[b]["ssabv_vs_slsqp"] = vs
            log(f"  budget {b}: avg ranks " + ", ".join(f"{LAB[a]} {v:.2f}" for a, v in sorted(avg.items(), key=lambda t: t[1])))
        SB = pd.concat(rows_c, ignore_index=True)
        SB.to_csv(OUT("mpce_budget_case_stats.csv"), index=False)
        summary["budget"] = dict(budgets=bl, per_budget=bsum,
                                 ssabv_best_at_all_budgets=all(bsum[b]["best"] == FOCUS for b in bl),
                                 note="HR16 at 120,030 calls uses 10 seeds (qualification: at least 5 feasible runs)")
        ms_all = [a for a in M9 if a in set(BD.Algorithm)]
        hdr = "Method & " + " & ".join(f"\\multicolumn{{3}}{{c}}{{{b:,} calls}}".replace(",", "{,}") for b in bl) + " \\\\\n & " + \
              " & ".join("Rank & Loss & Feas." for _ in bl)
        lines = []
        for a in ms_all:
            c = []
            for b in bl:
                v = bsum[b]
                c += [f"{v['avg_rank'][a]:.2f}" if a in v["avg_rank"] else "--",
                      f"{v['mean_loss_grid'][a]:.3f}" if np.isfinite(v["mean_loss_grid"].get(a, np.nan)) else "--",
                      f"{v['feas_pct'][a]:.0f}" if a in v["feas_pct"] else "--"]
            lines.append(f"{LAB[a]} & " + " & ".join(c) + " \\\\")
        tabs["budget"] = table("table*", "Budget scaling (random initialization) on the six largest benchmark cases and the Horns Rev 16-turbine block: average rank over the %d cases (ranking rule of Table~\\ref{tab:friedman68}), mean wake loss (\\%%) of the feasible runs averaged over the six benchmark cases in which the method has at least half of its runs feasible, and percentage of feasible runs. Per-case values: Table~\\ref{tab:budget-cases}." % len(cases7),
                               "tab:budget", "l" + "ccc" * len(bl), hdr, lines, sep="3pt")
        # per-case table
        dl = []
        for c in cases7:
            sc = SB[(SB.Dataset == c[0]) & (SB.Radius == c[1]) & (SB.Turbines == c[2])]
            if not len(sc):
                continue
            isHR = c[0] == "HR"
            unit = " (AEP, GWh/yr)" if isHR else " (wake loss, \\%)"
            cn = case_name("-".join(map(str, c)))
            dl.append(f"\\multicolumn{{{1 + 2 * len(bl)}}}{{l}}{{\\emph{{{cn}{unit}}}}} \\\\")
            for a in ms_all:
                cc = []
                for b in bl:
                    y = sc[(sc.Budget == b) & (sc.Algorithm == a)]
                    if not len(y):
                        cc += ["--", "--"]; continue
                    y = y.iloc[0]
                    v = y.Mean if isHR else y.Loss
                    cc += ["--" if not np.isfinite(v) else (f"{v:.2f}" if isHR else f"{v:.3f}"), f"{y.NFeas}/{y.N}"]
                dl.append(f"{LAB[a]} & " + " & ".join(cc) + " \\\\")
        tb = table("table", "Budget scaling per case: mean wake loss (\\%; Horns Rev: mean AEP in GWh/yr) of the feasible runs and feasible runs at each budget (random initialization).",
                   "tab:budget-cases", "l" + "cc" * len(bl),
                   "Method & " + " & ".join(f"{b:,}".replace(",", "{,}") + " & Feas." for b in bl), dl, pos="p")
        tabs["budget_cases"] = tb
        supp.append(tb)
        # figure
        fig, axes = plt.subplots(2, 4, figsize=(7.1, 3.9))
        for ax, c in zip(axes.flat, cases7):
            sc = SB[(SB.Dataset == c[0]) & (SB.Radius == c[1]) & (SB.Turbines == c[2])]
            isHR = c[0] == "HR"
            for a in ms_all:
                y = sc[sc.Algorithm == a].sort_values("Budget")
                v = (y.Mean if isHR else y.Loss).where(y.Qualified)
                ax.plot(y.Budget, v, color=COL[a], ls=LS[a], marker=MRK[a], ms=6 if MRK[a] == "*" else 3.5, lw=1.1)
            ax.set_xscale("log")
            ax.set_xticks(bl); ax.set_xticklabels([f"{b // 1000}k" if b % 1000 == 30 else str(b) for b in bl], fontsize=6)
            ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
            ax.set_title(case_name("-".join(map(str, c))).replace("Data Set ", "DS "), fontsize=7.5, color=INK)
            ax.set_ylabel("Mean AEP (GWh/yr)" if isHR else "Mean wake loss (%)", fontsize=7)
        axes.flat[-1].axis("off")
        h = [plt.Line2D([], [], color=COL[a], ls=LS[a], marker=MRK[a], ms=7 if MRK[a] == "*" else 4, label=LAB[a]) for a in ms_all]
        axes.flat[-1].legend(handles=h, loc="center", fontsize=7)
        fig.supxlabel("Objective-function calls (log scale)", fontsize=8)
        fig.tight_layout()
        FG.save(fig, "budget_scaling")
    else:
        log("  SKIPPED: mpce_b30k / mpce_b120k not available")
        summary["budget"] = None

    # =========================================================== 7. IEA37
    log("\n[7] IEA37 Case Study 1 (16- and 36-turbine scenarios)")
    IE = ALL[ALL.Dataset.str.startswith("IEA37")]
    if len(IE):
        pub = load_published(args.data_dir)
        isum, lines, plines, lay = {}, [], [], []
        for n, x in IE.groupby("Turbines"):
            ds = x.Dataset.iloc[0]
            cname = IEA_NAME.get(n, f"IEA37 CS1, {n} turbines")
            ideal = float(x.Ideal.iloc[0]) if x.Ideal.notna().any() else np.nan
            isum[cname] = dict(dataset_id=ds, turbines=int(n), ideal_aep=ideal, budgets={})
            for b, y in x.groupby("Budget"):
                ms = [a for a in M9 if a in set(y.Algorithm)]
                Sy = case_stats(y, ms).set_index("Algorithm")
                tt = {t["Baseline"]: t for t in paired_vs(y, FOCUS, [a for a in ms if a != FOCUS])} if FOCUS in ms else {}
                isum[cname]["budgets"][int(b)] = {
                    a: dict(mean=float(Sy.Mean[a]), best=float(Sy.Best[a]), sd=float(Sy.SD[a]), loss_pct=float(Sy.Loss[a]),
                            feasible=int(Sy.NFeas[a]), runs=int(Sy.N[a]), rank=float(Sy.Rank[a]),
                            p_holm=tt.get(a, {}).get("PHolm"), rb=tt.get(a, {}).get("RB")) for a in ms}
                bs = f"{b:,}".replace(",", "{,}")
                lines.append(f"\\multicolumn{{8}}{{l}}{{\\emph{{{cname}, {bs} calls}}}} \\\\")
                for a in sorted(ms, key=lambda a: Sy.Rank[a]):
                    s = Sy.loc[a]
                    if s.NFeas:
                        v = f"{s.Mean:,.0f} & {s.Best:,.0f} & {0 if not np.isfinite(s.SD) else s.SD:,.0f} & {s.Loss:.2f}".replace(",", "{,}")
                    else:
                        v = "-- & -- & -- & --"
                    lines.append(f"{LAB[a]} & {v} & {s.NFeas}/{s.N} & {s.Rank:g} & "
                                 f"{'--' if a == FOCUS or a not in tt else fmt_p(tt[a]['PHolm'])} \\\\")
            # comparison with the published Case Study 1 results
            ours = x[x.Feasible]
            if not len(ours):
                continue
            top = ours.sort_values("Objective").iloc[-1]
            fo = ours[ours.Algorithm == FOCUS]
            fbest = fo.sort_values("Objective").iloc[-1] if len(fo) else None
            fmean = x[(x.Algorithm == FOCUS) & x.Feasible].groupby("Budget").Objective.mean().to_dict()
            P = pub[pub.Turbines == n] if pub is not None else None
            base = P[P.Baseline] if P is not None else None
            parts = P[~P.Baseline] if P is not None else None
            bl_aep = float(base.AEP.iloc[0]) if base is not None and len(base) else BASELINE_IEA.get(n, np.nan)
            pf = parts[parts.Feasible].sort_values("AEP") if parts is not None else None
            pa = parts.sort_values("AEP") if parts is not None else None
            pf_aep, pf_who = (float(pf.AEP.iloc[-1]), pf.Participant.iloc[-1]) if pf is not None and len(pf) else (np.nan, "--")
            pa_aep, pa_who, pa_feas = (float(pa.AEP.iloc[-1]), pa.Participant.iloc[-1], bool(pa.Feasible.iloc[-1])) \
                if pa is not None and len(pa) else (np.nan, "--", None)
            fb = float(fbest.Objective) if fbest is not None else np.nan
            rank_among = int(1 + (pf.AEP > fb).sum()) if pf is not None and len(pf) else None
            isum[cname]["published"] = dict(
                baseline_aep=bl_aep, best_feasible_published=pf_aep, best_feasible_by=pf_who,
                best_overall_published=pa_aep, best_overall_by=pa_who, best_overall_feasible=pa_feas,
                n_published=int(len(parts)) if parts is not None else 0,
                n_published_feasible=int(len(pf)) if pf is not None else 0,
                our_best=float(top.Objective), our_best_method=top.Algorithm, our_best_budget=int(top.Budget), our_best_seed=int(top.Seed),
                ssabv_best=fb, ssabv_best_budget=int(fbest.Budget) if fbest is not None else None,
                ssabv_mean_by_budget={int(k): float(v) for k, v in fmean.items()},
                ssabv_best_vs_best_feasible_pct=float(100 * (fb / pf_aep - 1)) if np.isfinite(pf_aep) else None,
                ssabv_best_vs_baseline_pct=float(100 * (fb / bl_aep - 1)) if np.isfinite(bl_aep) else None,
                ssabv_mean_vs_baseline_pct={int(k): float(100 * (v / bl_aep - 1)) for k, v in fmean.items()},
                our_best_vs_best_feasible_pct=float(100 * (top.Objective / pf_aep - 1)) if np.isfinite(pf_aep) else None,
                ssabv_best_rank_among_feasible_published=rank_among)
            num = lambda v: "--" if not np.isfinite(v) else f"{v:,.1f}".replace(",", "{,}")
            plines.append(f"{cname} & {num(bl_aep)} & {num(pf_aep)} ({pf_who}) & {num(pa_aep)} ({pa_who}{'' if pa_feas in (None, True) else ', infeas.'}) & "
                          f"{num(top.Objective)} ({LAB[top.Algorithm]}) & {num(fb)} & "
                          + " / ".join(num(fmean[k]) for k in sorted(fmean)) +
                          f" & {'--' if not np.isfinite(pf_aep) else f'{100 * (fb / pf_aep - 1):+.2f}'} & "
                          f"{'--' if not np.isfinite(bl_aep) else f'{100 * (fb / bl_aep - 1):+.2f}'} \\\\")
            lay.append((n, cname, top, fbest, pf_who, x.Radius.iloc[0]))
        bl_txt = " / ".join(f"{b:,}".replace(",", "{,}") for b in sorted(IE.Budget.unique()))
        tabs["iea37"] = table("table", "IEA Wind Task~37 Case Study~1 (16- and 36-turbine scenarios; official AEP model): AEP (MWh) mean, best and SD over the feasible runs, mean wake loss relative to the wake-free AEP (\\%), feasible runs, rank (ranking rule of Table~\\ref{tab:friedman68}) and Holm-adjusted Wilcoxon signed-rank $p$ of SSA-VNS vs.\\ each method, per scenario and budget.",
                              "tab:iea37", "lccccccc", "Method & Mean & Best & SD & Loss (\\%) & Feas. & Rank & $p_{\\rm Holm}$", lines, sep="2.5pt")
        tabs["iea37_published"] = table("table*", "IEA Wind Task~37 Case Study~1: AEP (MWh, official calculator) of the baseline (example) layout, of the best feasible and the best overall participant layout%s, our best feasible layout over all methods and budgets, the best SSA-VNS layout, mean SSA-VNS AEP (budgets %s calls), and difference (\\%%) of the best SSA-VNS layout from the best feasible participant layout and from the baseline. ``infeas.'': violates the boundary or spacing constraint by more than 1~mm." % ("" if pub is not None else " (published results file not available)", bl_txt),
                                        "tab:iea37-pub", "lcccccccc", "Scenario & Baseline & Best feasible publ. & Best publ. & Our best & SSA-VNS best & SSA-VNS mean & $\\Delta_{\\rm publ}$ & $\\Delta_{\\rm base}$", plines, sep="2.5pt", resize=True)
        summary["iea37"] = isum
        summary["iea37_published_file"] = pub is not None
        for k, v in isum.items():
            p = v.get("published", {})
            log(f"  {k}: SSA-VNS best {p.get('ssabv_best', float('nan')):.1f}, best feasible published {p.get('best_feasible_published', float('nan')):.1f} "
                f"({p.get('best_feasible_by')}), baseline {p.get('baseline_aep', float('nan')):.1f}")
        # layouts: best SSA-VNS vs best feasible participant layout (and our best if another method)
        if lay:
            fig, axes = plt.subplots(1, len(lay), figsize=(3.55 * len(lay), 3.3), squeeze=False)
            for ax, (n, cname, top, fbest, pwho, rad) in zip(axes[0], lay):
                t = np.linspace(0, 2 * np.pi, 200)
                ax.plot(rad * np.cos(t), rad * np.sin(t), color=MUTED, lw=0.8)
                try:
                    import iea37_model as iem
                    if pwho not in ("--", None):
                        xy, _ = iem.load_submission(str(pwho).replace("par", ""), n)
                        ax.scatter(xy[:, 0], xy[:, 1], marker="s", s=22, facecolor="none", edgecolor=INK, lw=1,
                                   label=f"best feasible participant ({pwho})", zorder=3)
                except Exception as e:                                  # pragma: no cover
                    log(f"  participant layout not drawn: {e}")
                if fbest is not None:
                    xy = coords(fbest.Coordinates)
                    ax.scatter(xy[:, 0], xy[:, 1], marker="o", s=12, color=COL[FOCUS], label="best SSA-VNS", zorder=4)
                if top.Algorithm != FOCUS:
                    xy = coords(top.Coordinates)
                    ax.scatter(xy[:, 0], xy[:, 1], marker=MRK[top.Algorithm], s=16, facecolor="none", edgecolor=COL[top.Algorithm],
                               lw=1, label=f"our best ({LAB[top.Algorithm]})", zorder=3)
                ax.set_aspect("equal"); ax.tick_params(labelsize=6)
                ax.set_title(f"{cname} (m)", fontsize=8, color=INK)
                ax.legend(fontsize=6, loc="upper left", bbox_to_anchor=(0, -0.08), ncol=1)
            fig.tight_layout()
            FG.save(fig, "iea37_layouts")
        # convergence per scenario (largest budget)
        fig, axes = plt.subplots(1, IE.Turbines.nunique(), figsize=(7.1, 2.6), squeeze=False)
        for ax, (n, x) in zip(axes[0], IE.groupby("Turbines")):
            y = x[x.Budget == x.Budget.max()]
            ms = [a for a in M9 if a in set(y.Algorithm)]
            conv_panel(ax, y, ms, loss=True, band=False)
            ax.set_title(f"{IEA_NAME.get(n, n)}, {int(y.Budget.iloc[0]):,} calls", fontsize=8, color=INK)
            ax.set_xlabel("Objective-function calls"); ax.set_ylabel("Median best wake loss (%)")
        legend_row(fig, [a for a in M9 if a in set(IE.Algorithm)], ncol=5)
        fig.tight_layout()
        FG.save(fig, "iea37_convergence")
    else:
        log("  SKIPPED: mpce_iea16 / mpce_iea36 not available")
        summary["iea37"] = None
    # =========================================================== 8. robustness
    if args.skip_robust:
        log("\n[8] Robustness SKIPPED (--skip-robust)")
    else:
        robustness(G, MAINP, args, summary, tabs, OUT)

    # =========================================================== write
    hdr = "%% generated by mpce_results.py -- do not edit by hand\n"
    for fn in glob.glob(OUT("mpce_tab_*.tex")):             # remove tables of sections skipped in this run
        if os.path.basename(fn)[9:-4] not in tabs and not (args.skip_robust and fn.endswith("_robust.tex")):
            os.remove(fn); log(f"  removed stale {os.path.basename(fn)}")
    for k, v in tabs.items():
        open(OUT(f"mpce_tab_{k}.tex"), "w").write(hdr + v)
    open(OUT("mpce_supplementary.tex"), "w").write(
        hdr + "%% Supplementary tables (per-case results); requires booktabs, graphicx\n" + "\n".join(supp))
    summary["tables"] = [f"mpce_tab_{k}.tex" for k in tabs] + ["mpce_supplementary.tex"]
    summary["figures"] = [f"figures_mpce/{f}.pdf" for f in FG.made]
    summary["log"] = LOG
    json.dump(clean(summary), open(OUT("mpce_summary.json"), "w"), indent=1)
    log(f"\nwrote {len(tabs)} tables, {len(FG.made)} figures, mpce_summary.json in {time.time() - t0:.0f} s")


def robustness(G, methods, args, summary, tabs, OUT, only=False):
    log("\n[8] Robustness re-evaluation")
    import mpce_robustness as mr
    F = G[G.Feasible & G.Algorithm.isin(methods)]
    RV = mr.reevaluate(F, HERE, args.procs, log)
    tex, out = mr.robust_table(RV, methods, LAB, rank_rule)
    tabs["robust"] = tex
    summary["robustness"] = out
    RV.drop(columns=["Coordinates", "Curve", "Key"]).to_csv(OUT("mpce_reevaluation.csv"), index=False)
    for k, v in out.items():
        if k == "rounding_check":
            log(f"  re-evaluation of the rounded coordinates vs recorded objective: {v}"); continue
        log(f"  {k:12s} tau={v['tau']:.2f} same best={v['same_best_pct']:.0f}%  SSA-VNS rank {v['avg_rank'][FOCUS]:.2f}")
    if only:
        open(OUT("mpce_tab_robust.tex"), "w").write("%% generated by mpce_results.py -- do not edit by hand\n" + tex)


def case_name(c):
    p = c.split("-")
    if p[0] == "HR":
        return f"Horns Rev 1, {p[-1]} turbines"
    if p[0] in DSN:
        return f"{DSN[p[0]]}, $r$ = {p[1]} m, $N$ = {p[2]}"
    return c


def load_published(data_dir):
    """iea37_published_results.csv -> DataFrame(Turbines, Participant, AEP, Feasible, Baseline).

    Expected columns (written by the IEA37 agent): Case, Turbines, Participant, AEP_MWh_official_calc,
    Feasible_tol1e-3m (baseline row: Participant starting with "baseline"). Falls back to any
    case/turbines, participant/algorithm and AEP columns if the names differ."""
    for d in (data_dir, HERE, os.path.join(HERE, "iea37_data")):
        fn = os.path.join(d, "iea37_published_results.csv")
        if not os.path.exists(fn):
            continue
        P = pd.read_csv(fn)
        low = {c.lower(): c for c in P.columns}
        pick = lambda *names: next((low[c] for nm in names for c in low if nm in c), None)
        tc = pick("turbines"); cc = pick("case")
        pc = pick("participant", "algorithm", "method", "name")
        ac = pick("aep_mwh_official", "aep_mwh_reported", "aep")
        fc = pick("feasible")
        if ac is None or (tc is None and cc is None):
            log(f"  WARNING: {fn} lacks turbines/AEP columns: {list(P.columns)}")
            return None
        turb = pd.to_numeric(P[tc], errors="coerce") if tc else \
            P[cc].astype(str).str.extract(r"(\d+)", expand=False).astype(float)
        part = P[pc].astype(str) if pc else pd.Series("?", index=P.index)
        feas = P[fc].astype(str).str.lower().isin(["true", "1", "1.0", "yes"]) if fc else pd.Series(True, index=P.index)
        out = pd.DataFrame({"Turbines": turb, "Participant": part, "AEP": pd.to_numeric(P[ac], errors="coerce"),
                            "Feasible": feas, "Baseline": part.str.lower().str.startswith("baseline")})
        out = out.dropna(subset=["AEP", "Turbines"])
        out["Turbines"] = out.Turbines.astype(int)
        log(f"  published IEA37 results: {fn} ({len(out)} rows; AEP column {ac}; feasibility column {fc})")
        return out
    log("  published IEA37 results file not found (comparison table written with placeholders)")
    return None

def clean(o):
    if isinstance(o, dict):
        return {str(k): clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [clean(v) for v in o]
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating, float)):
        return None if not np.isfinite(o) else float(o)
    if isinstance(o, np.bool_):
        return bool(o)
    return o


if __name__ == "__main__":
    main()

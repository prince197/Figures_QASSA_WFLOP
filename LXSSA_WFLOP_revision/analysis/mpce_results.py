"""Results and statistics pipeline for the MPCE resubmission; the proposed ("focus") method is switchable.

Usage:  python mpce_results.py [--focus {PSOBV,PSOC,SSABV}] [--data-dir D] [--out-dir O] [--fig-dir F]
                               [--procs 4] [--partial] [--common-seeds] [--skip-robust] [--only-robust]

  --focus         proposed method (default PSOBV = PSO-VNS; PSOC = PSO with constriction coefficients;
                  SSABV = SSA-VNS). Every section follows it: main comparison (focus first, then the
                  other seven of PSOBV, PSOC, SSABV, SSA, LXSSA, DE, BVNS, SLSQP), Friedman / Holm post
                  hoc, case-mean Wilcoxon and run-level W/T/L against the focus, ablation (built around
                  the focus if it is a two-phase hybrid, otherwise around the hybrid whose phase 1 is the
                  focus), split table (only if split data exist for the focus), figures (focus drawn with
                  a star marker, solid and thicker line, on top), captions and summary keys (generic:
                  summary["focus"], "focus_*").
  --common-seeds  preview aid for incomplete shards: within every case x budget x initialization, keep
                  only the seeds that every method of that group has (so means and ranks are seed-paired).
  The "decision" block (summary["decision"], printed last) compares PSOBV with PSOC and SSABV
  independently of --focus.

Inputs (all runs at 6,030 calls with random initialization unless stated)
  fresh_grid.csv  (LXSSA, SSA, PSO [old settings; development fallback for PSOC], DE,
                   SLSQP [old platform; development fallback for the SLSQP rerun]; VNS ignored)
  fresh_vgrid.csv (BVNS = basic VNS, "VNS"),  fresh_bgrid.csv (LXBV, SSABV)
  fresh_hr16.csv, fresh_vhr16.csv, fresh_bhr16.csv (Horns Rev 1, 16 turbines; AEP/IdealAEP)
  mpce_<exp>_s<i>of<k>.csv from mpce_experiments.py / iea37_experiments.py (read from --data-dir; all
  shards of an experiment are merged; an experiment is used only when all k shards are present, unless
  --partial is given): rsvns, psoc, psobv (68 cases + HR16), slsqp, psosplit (PSOBV25 / PSOBV75 on the
  12 split cases), omega90 (Phase 6: PSOBV90 = PSO-VNS with omega = 0.9 on the 12 split cases), rsdisc (Phase 6:
  RSDVNS = RS-VNS whose Phase-1 samples are uniform in the farm disc, 68 cases; component-analysis control),
  ssasplit (dropped; used if present), hr16new (PSOC, RSVNS on HR16), feas, b30k,
  b120k, iea16, iea36 (nine methods M9) and the PSO-VNS-only arms feasp, b30kp, b120kp, iea16p, iea36p,
  which are merged with the nine-method experiments so that PSO-VNS is a tenth method (M10).
  hrfix (24 shards): ALL Horns Rev 1 runs again after the direction-binning fix (commit 7676da9); once any hrfix
  shard exists, every Horns Rev row of every other source is dropped (see load()).
  iea37_published_results.csv (optional; columns case, participant/algorithm, AEP).

Outputs (in --out-dir, default: this folder)
  mpce_tab_<name>.tex (one table per file, ready to \\input), mpce_supplementary.tex (per-case
  tables), mpce_summary.json (every number quoted in the text), mpce_case_stats.csv,
  mpce_case_tests.csv, mpce_ablation_tests.csv, mpce_best_layouts_maxN.csv; figures (pdf + png)
  in --fig-dir (default ../figures_mpce). At the end, mpce_numbers.py writes mpce_numbers.tex (the "N..." macros
  of the manuscript) and mpce_check_final.py prints PASS / FAIL of its CHECK-FINAL statements. The main-text
  tables (mpce_tab_*.tex) follow the layout of MPCE_PSO_VNS.tex; everything else goes to mpce_supplementary.tex
  (input by MPCE_PSO_VNS_supplement.tex). mpce_tab_baseline.tex (tab:baseline, summary["baseline"]): the previous
  study's method pool (LXBV, SSABV, LXSSA, SSA, PSO, DE, BVNS, SLSQP) ranked with the old and with the constriction
  PSO setting (baseline_setting). Full rebuild incl. LaTeX: sh build_mpce_paper.sh [--partial].

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
from scipy.stats import wilcoxon, friedmanchisquare, rankdata, norm, f as f_dist, t as t_dist

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

FOCUS_CHOICES = ("PSOBV", "PSOC", "SSABV")
FOCUS = "PSOBV"                               # set by --focus (set_focus)
MAIN8 = ["PSOBV", "PSOC", "SSABV", "SSA", "LXSSA", "DE", "BVNS", "SLSQP"]      # main comparison pool
MAIN = list(MAIN8)                            # focus first (set_focus)
M9 = ["SSABV", "LXBV", "RSVNS", "BVNS", "SSA", "LXSSA", "PSOC", "DE", "SLSQP"]  # nine-method experiments
M10 = MAIN8 + ["LXBV", "RSVNS"]              # M9 + PSO-VNS-only arms; focus first (set_focus)
HYBRIDS = ["PSOBV", "SSABV", "LXBV"]          # two-phase hybrids (phase 1 swarm, phase 2 basic VNS)
PHASE1 = {"PSOBV": "PSOC", "SSABV": "SSA", "LXBV": "LXSSA"}
PHASE1_ITER_CALLS = {"PSOC": 30, "SSA": 30, "LXSSA": 60, "RS": 30}  # calls per phase-1 iteration (Np = 30)
DECISION_PAIRS = [("PSOBV", "PSOC"), ("PSOBV", "SSABV")]
LAB = {"PSOBV": "PSO-VNS", "SSABV": "SSA-VNS", "SSABV25": "SSA-VNS (25\\%)", "SSABV75": "SSA-VNS (75\\%)",
       "PSOBV25": "PSO-VNS (25\\%)", "PSOBV75": "PSO-VNS (75\\%)", "PSOBV90": "PSO-VNS (90\\%)", "LXBV": "LX-SSA-VNS",
       "RSVNS": "RS-VNS", "RSDVNS": "RSD-VNS", "BVNS": "VNS", "SSA": "SSA", "LXSSA": "LX-SSA", "PSOC": "PSO", "PSO": "PSO (old)",
       "DE": "DE", "SLSQP": "MS-SLSQP"}
# Colour follows the method, independent of the focus (categorical palette of final_results.py; the
# eight slots blue, orange, aqua, yellow, magenta, green, violet, red go to the eight methods of the
# main comparison). The three focus candidates take violet / blue / aqua, the only trio of the palette
# that passes the all-pairs CVD and normal-vision checks of the dataviz validator (worst pair
# blue-violet: CVD dE 13.0, normal dE 16.3; red-aqua would be CVD dE 6.9 and red-orange a hard fail).
# Hence PSO-VNS = violet, LX-SSA moves to magenta and VNS to red. LX-SSA-VNS and RS-VNS appear only in
# the ablation / nine-method sections: two neutral greys (dE 15.5) with their own markers and dashes
# (composite encoding instead of generated hues). The disc-sampling control RSD-VNS (Phase 6) is a third, lighter
# grey with its own marker (octagon) and dash pattern (same composite encoding). The focus is additionally drawn
# with a star marker, a solid, thicker line and on top (see mk / ls / lw).
COL = {"PSOBV": "#4a3aa7", "SSABV": "#2a78d6", "PSOC": "#1baf7a", "PSO": "#1baf7a", "SSA": "#eb6834",
       "DE": "#eda100", "LXSSA": "#e87ba4", "SLSQP": "#008300", "BVNS": "#e34948", "LXBV": "#5f5e5a",
       "RSVNS": "#8d8b85", "RSDVNS": "#aeaba4"}
MRK = {"PSOBV": "p", "SSABV": "d", "SSA": "s", "PSOC": "^", "PSO": "^", "DE": "v", "BVNS": "D", "SLSQP": "P",
       "LXSSA": "o", "LXBV": "X", "RSVNS": "h", "RSDVNS": "8"}
LS = {"PSOBV": (0, (7, 2)), "SSABV": (0, (3, 1.5)), "SSA": "--", "PSOC": "-.", "PSO": "-.", "DE": ":",
      "BVNS": (0, (5, 1)), "SLSQP": (0, (3, 1, 1, 1)), "LXSSA": (0, (1, 1)), "LXBV": (0, (4, 2, 1, 2)),
      "RSVNS": (0, (6, 2, 2, 2)), "RSDVNS": (0, (2, 1.2))}


def set_focus(f):
    """Make every section follow the focus method f."""
    global FOCUS, MAIN, M10
    assert f in FOCUS_CHOICES, f
    FOCUS = f
    MAIN = [f] + [a for a in MAIN8 if a != f]
    M10 = MAIN + ["LXBV", "RSVNS"]


def mk(a):
    return "*" if a == FOCUS else MRK[a]


def ms(a, base=3.5):
    return base * 1.8 if a == FOCUS else base


def ls(a):
    return "-" if a == FOCUS else LS[a]


def lw(a, base=1.2):
    return base * 1.5 if a == FOCUS else base


def zo(a):
    return 5 if a == FOCUS else 3


def switch_call(alg, budget, np_=30, split=0.5):
    """Number of Phase-1 evaluations of a two-phase variant, i.e. the call after which the VNS phase starts.
    Swarm hybrids (HybridBVNS iteration arithmetic): np_ + round((split B - np_) / per) x per, e.g. 3,030 for
    PSO-VNS at 6,030. Split variants carry the share in their label: PSOBV25 / PSOBV75 / PSOBV90 (omega = 0.25 /
    0.75 / 0.9; PSOBV90 -> 30 + 180 x 30 = 5,430 PSO evaluations at 6,030). RSVNS and the disc-sampling control
    RSDVNS (Phase 6) sample one layout per call: n1 = round(split B) = 3,015 Phase-1 evaluations at 6,030 (rs_vns.RSVNS,
    self.n1; incl. the common initial population of 30). (Until R3-9 of review round 2 RS/RSD used the iteration
    arithmetic too, which gave 3,030 and read the switch 15 VNS evaluations late.)"""
    if alg in ("RSVNS", "RSDVNS"):
        return int(round(split * budget))
    p1 = PHASE1[alg.rstrip("0123456789")]
    if alg[-2:] in ("25", "75", "90"):
        split = int(alg[-2:]) / 100
    per = PHASE1_ITER_CALLS[p1]
    return np_ + max(1, int(round((split * budget - np_) / per))) * per


def switch_index(sw, budget):
    """Index of the last convergence checkpoint at or before call sw (checkpoint k is at call (k + 1) x step,
    step = (B - 30) // 200 = 30 at 6,030): 3,030 -> index 100 (call 3,030), 3,015 -> index 99 (call 3,000; the
    last 15 Phase-1 samples of RS-VNS / RSD-VNS are not observable in the stored curves)."""
    return int(sw // max(1, (budget - 30) // 200)) - 1


HR_ORDER = ["PSOBV", "SSABV", "LXBV", "RSVNS", "BVNS", "PSOC", "SSA", "LXSSA", "DE", "SLSQP"]   # main-text order (tables)
HEAD2 = {a: "\\begin{tabular}{@{}c@{}}%s\\\\%s\\end{tabular}" % tuple(LAB[a].split("-", 1)[0:1] + [LAB[a].split("-", 1)[1]])
         for a in ("SSABV", "LXSSA", "SLSQP", "LXBV", "RSVNS", "PSOBV")}
HEAD2 = {a: v.replace("\\\\", "-\\\\", 1) for a, v in HEAD2.items()}      # two-line column heads, e.g. SSA-/VNS
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
# Horns Rev 1. The direction binning of hornsrev_model.py was fixed in commit 7676da9 (np.round gave the 12
# sectors 7 and 5 five-degree bins alternately, total frequency 1.003); every Horns Rev run is repeated as
# experiment `hrfix` (mpce_hrfix_s<i>of24.csv). The installed and wake-free AEPs are recomputed with
# hornsrev_model; INSTALLED_HR16 is only the fallback if the module cannot be imported (fixed binning).
INSTALLED_HR16 = 139.821         # GWh/yr, installed 16-turbine block, fixed binning (hornsrev_model.aep_gwh)
# The old runs (fresh_hr16/vhr16/bhr16, hr16new, psobv, feas/feasp, b30k(p), b120k(p)) were evaluated with the
# old binning; while they are still loaded (no hrfix file yet), the installed AEP of the same (old) model is
# used so that the comparison stays within one model. Recognized by the wake-free AEP stored with the runs.
LEGACY_HR16 = dict(ideal=149.2632720953485, installed=139.51300511740232)
# PyWake reference for the complete 80-turbine farm (R4-2, D19): read from pywake_check.csv (pywake_check.py,
# PyWake 2.6.20 NOJ(Hornsrev1Site(), V80(), k=0.04) at direction/speed bins IDENTICAL to hornsrev_model.py, row
# Farm=HR80, Bins=ours_5deg_2.5, Model=NOJ_k0.04). The former constant 662.5 GWh/yr (10.96 %, quoted from an older
# text) is not reproduced by PyWake at any bin setting and is no longer used. \NHRPyWakeDiff = relative
# difference (%) of the installed 80-turbine AEP of hornsrev_model.py from that PyWake value; pending if the file
# is missing.
PYWAKE_CHECK_CSV = "pywake_check.csv"
PYWAKE_ROW = dict(Farm="HR80", Bins="ours_5deg_2.5", Model="NOJ_k0.04", Source="PyWake")


def pywake_reference(data_dir=None):
    """PyWake NOJ reference of the installed 80-turbine farm at our bins (dict) or None if pywake_check.csv is absent."""
    for d in ([data_dir] if data_dir else []) + [HERE]:
        fn = os.path.join(d, PYWAKE_CHECK_CSV)
        if os.path.exists(fn):
            P = pd.read_csv(fn)
            r = P[np.logical_and.reduce([P[k].astype(str) == v for k, v in PYWAKE_ROW.items()])]
            o = P[(P.Farm == PYWAKE_ROW["Farm"]) & (P.Bins == PYWAKE_ROW["Bins"]) & (P.Model == "paper_model") & (P.Source == "ours")]
            if len(r):
                r = r.iloc[0]
                return dict(aep=float(r.AEP_GWh), loss_pct=float(r.WakeLossPct), ideal=float(r.IdealAEP_GWh),
                            version=str(r.Version), settings=str(r.Settings), bins=PYWAKE_ROW["Bins"],
                            p_max_abs_diff=float(r.PMaxAbsDiff),
                            ours_csv_aep=float(o.AEP_GWh.iloc[0]) if len(o) else None, file=fn)
    return None
BOOT_N, BOOT_SEED = 10000, 20260928   # bootstrap over cases (percentile 95 % CI), fixed seed
# Practical equivalence of case-mean wake losses (equivalence_block): one margin for EVERY case-mean comparison,
# in percentage points (pp) of wake loss. It was chosen AFTER the primary analysis (post hoc, when "PSO-VNS is
# not significantly different from PSO" had to be turned into a practical-equivalence statement); the paper
# therefore also reports, for every pair, the smallest margin at which equivalence holds (min_margin_pp), so a
# reader can apply any other margin. Do not tune this constant to make a particular pair (non-)equivalent.
EQ_MARGIN = 0.05
# Bayesian signed-rank test (Benavoli et al. 2017, JMLR 18:77): Dirichlet-process prior strength s, pseudo-
# observation z0, Monte Carlo samples of the posterior, fixed seed
BAYES_S, BAYES_Z0, BAYES_N, BAYES_SEED = 0.5, 0.0, 50000, 20260928
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


def fmt_p(p, d=3):
    """p for tables: d significant digits, scientific notation below 0.001."""
    if not np.isfinite(p):
        return "--"
    if p < 1e-3:
        m, e = f"{p:.{d - 1}e}".split("e")
        return f"${m}\\times10^{{{int(e)}}}$"
    return f"{p:.3f}"


def goodness(df):
    """Seed-level score, higher = better; infeasible runs below every feasible run."""
    return np.where(df.Feasible, df.Objective, -1e12 - np.maximum(0, 1e4 - df.MinSpacing))


def wil(d):
    """Two-sided Wilcoxon signed-rank on paired differences d (zeros dropped) + rank-biserial r.

    Zero handling (R3 item 17), exposed as is: differences with |d| <= 1e-9 are dropped before both the test
    and the matched-pairs rank-biserial correlation (Wilcoxon's "wilcox" zero method), so r_rb describes only
    the non-tied pairs; if every difference is zero the test is not run and p = 1, r_rb = 0 (this happens
    for the small-N cases in which all methods reach the same layout). The median r_rb over cases therefore
    mixes such degenerate cases (r_rb = 0) with informative ones; captions say so."""
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


def friedman_block(R, methods, focus=None):
    """Case-level Friedman, Iman-Davenport, average ranks, Holm post hoc z tests vs focus and vs best."""
    focus = focus or (FOCUS if FOCUS in methods else methods[0])
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


def boot_means(x, n=BOOT_N, seed=BOOT_SEED):
    """Bootstrap distribution of the mean of x (n resamples of the cases with replacement, fixed seed)."""
    x = np.asarray(x, float); x = x[np.isfinite(x)]
    if len(x) < 2:
        return None
    rng = np.random.default_rng(seed)
    return x[rng.integers(0, len(x), (n, len(x)))].mean(1)


def boot_ci(x, n=BOOT_N, seed=BOOT_SEED, level=0.95):
    """Percentile bootstrap CI of the mean of x (resampling the cases with replacement, fixed seed)."""
    m = boot_means(x, n, seed)
    if m is None:
        return [None, None]
    a = (1 - level) / 2
    return [float(np.quantile(m, a)), float(np.quantile(m, 1 - a))]


def tost(d, margin=EQ_MARGIN, n=BOOT_N, seed=BOOT_SEED):
    """Equivalence (TOST) of the mean of the case-mean differences d (pp) within (-margin, margin).

    The bootstrap resamples are those of boot_ci (same n, same seed), so ci95 equals the 95 % CI of the \\NCm...
    macros. Decision rule: equivalence holds iff the 90 % percentile bootstrap CI lies strictly inside
    (-margin, margin) (TOST at alpha = 0.05). TOST p-value (primary, "p_tost"): NONPARAMETRIC BOOTSTRAP, the
    inversion of the same percentile interval: p_low = share of bootstrap means <= -margin, p_high = share
    >= +margin, each with the add-one correction (k + 1) / (n + 1), p_tost = max(p_low, p_high); its floor is
    1 / (n + 1) (p_tost_at_floor = True: no bootstrap mean outside the margin, i.e. p <= 1/(n+1)). A paired-t
    TOST on the case means (Schuirmann) is stored for reference as p_tost_t. min_margin_pp = max(|lo90|, |hi90|):
    equivalence holds for every margin larger than this."""
    d = np.asarray(d, float); d = d[np.isfinite(d)]
    m = boot_means(d, n, seed)
    if m is None:
        return None
    lo90, hi90 = float(np.quantile(m, 0.05)), float(np.quantile(m, 0.95))
    k_lo, k_hi = int((m <= -margin).sum()), int((m >= margin).sum())
    p_lo, p_hi = (k_lo + 1) / (n + 1), (k_hi + 1) / (n + 1)
    se = d.std(ddof=1) / np.sqrt(len(d))
    if se > 0:
        pt = max(float(t_dist.sf((d.mean() + margin) / se, len(d) - 1)), float(t_dist.cdf((d.mean() - margin) / se, len(d) - 1)))
    else:
        pt = 0.0 if abs(d.mean()) < margin else 1.0
    return dict(n_cases=int(len(d)), mean_dloss_pp=float(d.mean()),
                ci95_mean_dloss_pp=[float(np.quantile(m, 0.025)), float(np.quantile(m, 0.975))],
                ci90_mean_dloss_pp=[lo90, hi90], margin_pp=margin,
                equivalent=bool(-margin < lo90 and hi90 < margin),
                p_tost=float(max(p_lo, p_hi)), p_tost_low=float(p_lo), p_tost_high=float(p_hi),
                p_tost_at_floor=bool(k_lo == 0 and k_hi == 0), p_tost_t=pt,
                min_margin_pp=float(max(abs(lo90), abs(hi90))))


def bayes_signrank(d, rope=EQ_MARGIN, s=BAYES_S, z0=BAYES_Z0, n=BAYES_N, seed=BAYES_SEED):
    """Bayesian signed-rank test with a region of practical equivalence (Benavoli, Corani, Demsar, Zaffalon 2017,
    "Time for a change", JMLR 18:77, Sec. 5), implemented with numpy.

    d = case-mean differences L(A) - L(B) (pp; negative = A better), rope = (-rope, rope). Prior: Dirichlet
    process with strength s and base measure concentrated on the pseudo-observation z0, so the posterior weights
    of (z0, d_1, ..., d_n) are Dirichlet(s, 1, ..., 1). For every Monte Carlo sample w of the weights
      theta_A   = sum_ij w_i w_j [z_i + z_j < -2 rope]      (A better)
      theta_B   = sum_ij w_i w_j [z_i + z_j > +2 rope]      (B better)
      theta_rope = 1 - theta_A - theta_B                   (sums exactly at +-2 rope count 1/2, as in baycomp)
    (i, j over all n + 1 points incl. i = j and z0; the Walsh averages (z_i + z_j) / 2 compared with the rope).
    Reported: P(A better), P(practically equivalent), P(B better) = share of the posterior samples in which
    theta_A, theta_rope, theta_B is the largest of the three (ties split equally), as in the paper's simplex plots.
    Unit checks (_bayes_selftest, run by equivalence_block): all d = 0 -> P(rope) = 1 (every Walsh average is 0);
    all d = -1 pp with rope 0.05 -> P(A better) ~ 1; all d = +1 pp -> P(B better) ~ 1."""
    d = np.asarray(d, float); d = d[np.isfinite(d)]
    z = np.r_[z0, d]
    ws = np.r_[s, np.ones(len(d))]
    W = np.random.default_rng(seed).dirichlet(ws, n)                    # (n, len(z))
    sums = z[:, None] + z[None, :]
    below = (sums < -2 * rope) + 0.5 * (sums == -2 * rope)
    above = (sums > 2 * rope) + 0.5 * (sums == 2 * rope)
    tA = ((W @ below) * W).sum(1)
    tB = ((W @ above) * W).sum(1)
    T = np.c_[tA, 1 - tA - tB, tB]
    win = (T >= T.max(1, keepdims=True) - 1e-15).astype(float)
    pr = (win / win.sum(1, keepdims=True)).mean(0)
    return dict(p_a_better=float(pr[0]), p_rope=float(pr[1]), p_b_better=float(pr[2]), rope_pp=rope,
                prior_strength=s, prior_z0=z0, mc_samples=n, seed=seed, n_cases=int(len(d)),
                mean_theta=[float(tA.mean()), float((1 - tA - tB).mean()), float(tB.mean())])


def _bayes_selftest():
    """Unit checks of bayes_signrank (see its docstring); cheap (2,000 samples)."""
    r0 = bayes_signrank(np.zeros(20), n=2000)
    rA = bayes_signrank(-np.ones(20), n=2000)
    rB = bayes_signrank(np.ones(20), n=2000)
    ok = r0["p_rope"] > 0.999 and rA["p_a_better"] > 0.99 and rB["p_b_better"] > 0.99
    assert ok, ("bayes_signrank self-test failed", r0, rA, rB)
    return ok


# the component-analysis pairs of the \NCm... macros (ablation.case_mean keys; first minus second): the nine of
# Phase 4 and (Phase 6) the four contrasts of the disc-sampling control RSD-VNS (used only if its runs exist)
EQ_PAIRS_ABL = ["SSABV-RSVNS", "LXBV-RSVNS", "PSOBV-RSVNS", "SSABV-LXBV", "LXSSA-SSA", "SSABV-SSA", "LXBV-LXSSA",
                "PSOBV-PSOC", "PSOBV-BVNS", "SSABV-RSDVNS", "LXBV-RSDVNS", "PSOBV-RSDVNS", "RSDVNS-RSVNS"]


def case_mean_diffs(S, a, b, val="Loss"):
    """Per-case differences val(a) - val(b) over the cases in which BOTH methods qualify (>= half of the runs
    feasible) -- the case set of mean_dloss_pp / ci95_mean_dloss_pp of case_mean_wilcoxon (\\NCm...DL, CI).
    Cases in which only one method qualifies (imputed in the Wilcoxon test only) or neither does are excluded;
    their counts are returned."""
    P = S.pivot_table(index=CASE, columns="Algorithm", values=val)
    Q = S.pivot_table(index=CASE, columns="Algorithm", values="Qualified").astype(float)
    qa, qb = Q[a] > 0.5, Q[b] > 0.5
    both = qa & qb
    return (P[a] - P[b])[both], dict(only_a_qualified=int((qa & ~qb).sum()), only_b_qualified=int((~qa & qb).sum()),
                                     neither_qualified=int((~qa & ~qb).sum()), n_both_qualified=int(both.sum()))


def equivalence_block(S, SA, G, summary, supp):
    """Practical equivalence of the case-mean wake losses at EQ_MARGIN (summary["equivalence"], tab:equivalence)
    and spread / worst run of PSO-VNS vs PSO (summary["spread"]).

    Pairs: the component-analysis pairs of the \\NCm... macros (EQ_PAIRS_ABL) (case stats SA of the ablation variants) and
    the focus vs each other method of the main comparison (case stats S); a pair in both groups is computed from
    both and must agree (same runs). d = L(A) - L(B) (pp, negative = A better) over the cases in which both
    methods qualify (case_mean_diffs). tost(): 90 % bootstrap CI (BOOT_N case resamples, seed BOOT_SEED),
    bootstrap TOST p, minimal margin; bayes_signrank(): Bayesian signed-rank test with ROPE (-EQ_MARGIN,
    EQ_MARGIN)."""
    _bayes_selftest()
    pairs, groups = {}, {}
    main_pairs = [f"{FOCUS}-{b}" for b in summary.get("main_methods", []) if b != FOCUS]
    for grp, keys, Sx in (("ablation", EQ_PAIRS_ABL, SA), ("main", main_pairs, S)):
        for key in keys:
            a, b = key.split("-")
            if Sx is None or a not in set(Sx.Algorithm) or b not in set(Sx.Algorithm):
                continue
            d, cnt = case_mean_diffs(Sx, a, b)
            r = tost(d.values)
            if r is None:
                continue
            r.update(cnt)
            r["bayes"] = bayes_signrank(d.values)
            if key in pairs:                                   # same runs -> identical case means
                assert np.isclose(pairs[key]["mean_dloss_pp"], r["mean_dloss_pp"]) and pairs[key]["n_cases"] == r["n_cases"], key
            else:
                pairs[key] = r
            groups.setdefault(key, []).append(grp)
    for key in pairs:
        pairs[key]["groups"] = groups[key]
    summary["equivalence"] = dict(
        margin_pp=EQ_MARGIN,
        margin_note="one margin for every case-mean comparison, chosen after the primary analysis (post hoc); "
                    "min_margin_pp = smallest margin at which the 90 % CI lies inside (-m, m)",
        rule="equivalent iff the 90 % percentile bootstrap CI of the mean case-mean difference lies strictly inside (-m, m) (TOST, alpha 0.05)",
        p_tost="nonparametric bootstrap TOST p (add-one, floor 1/(BOOT_N+1)); p_tost_t = paired-t TOST on the case means, for reference",
        bootstrap=dict(resamples=BOOT_N, seed=BOOT_SEED, unit="cases"),
        bayes=dict(test="Bayesian signed-rank test (Benavoli et al. 2017, JMLR 18:77)", prior_strength=BAYES_S,
                   prior_z0=BAYES_Z0, mc_samples=BAYES_N, seed=BAYES_SEED, rope=[-EQ_MARGIN, EQ_MARGIN],
                   probabilities="share of posterior samples in which each of theta_A, theta_rope, theta_B is the largest"),
        case_set="cases in which both methods have >= 15 of 30 feasible runs (as mean_dloss_pp / ci95 of the \\NCm macros)",
        sign="d = L(A) - L(B), pp; negative = A (first method) better",
        pairs=pairs)
    for key, r in pairs.items():
        a, b = key.split("-")
        log(f"  EQ {LAB[a]:>10s} vs {LAB[b]:10s} n={r['n_cases']:2d} dL={r['mean_dloss_pp']:+.4f} "
            f"90%CI=[{r['ci90_mean_dloss_pp'][0]:+.4f}, {r['ci90_mean_dloss_pp'][1]:+.4f}] p_TOST={r['p_tost']:.4g} "
            f"(t {r['p_tost_t']:.3g}) m_min={r['min_margin_pp']:.4f} eq@{EQ_MARGIN}={'yes' if r['equivalent'] else 'no'}  "
            f"Bayes P(A)={r['bayes']['p_a_better']:.3f} P(rope)={r['bayes']['p_rope']:.3f} P(B)={r['bayes']['p_b_better']:.3f}")

    # --- supplementary table (tab:equivalence)
    f3 = lambda v: f"{v:.3f}".replace("-", "$-$")
    ci = lambda c: f"[{f3(c[0])}, {f3(c[1])}]"
    up3 = lambda v: f"{np.ceil(v * 1000 - 1e-9) / 1000:.3f}"
    pr2 = lambda v: "$>$0.99" if v > 0.99 else "$<$0.01" if v < 0.01 else f"{v:.2f}"     # never 1.00 / 0.00
    def row(key):
        a, b = key.split("-"); r = pairs[key]; y = r["bayes"]
        pt = ("$\\le" + fmt_p(r["p_tost"], 2).strip("$") + "$") if r["p_tost_at_floor"] else fmt_p(r["p_tost"])
        return (f"{LAB[a]} vs.\\ {LAB[b]} & {r['n_cases']} & ${r['mean_dloss_pp']:+.3f}$ & {ci(r['ci95_mean_dloss_pp'])} & "
                f"{ci(r['ci90_mean_dloss_pp'])} & {pt} & {up3(r['min_margin_pp'])} & {'yes' if r['equivalent'] else 'no'} & "
                f"{pr2(y['p_a_better'])} & {pr2(y['p_rope'])} & {pr2(y['p_b_better'])} \\\\")
    ab_ = [k for k in EQ_PAIRS_ABL if k in pairs]
    mn_ = [k for k in main_pairs if k in pairs and k not in ab_]
    dup = [k for k in main_pairs if k in ab_]
    lines = (["\\multicolumn{11}{l}{\\emph{Component-analysis contrasts}} \\\\"] + [row(k) for k in ab_]
             + ["\\midrule", "\\multicolumn{11}{l}{\\emph{Main comparison, %s vs.\\ each method%s}} \\\\"
                % (LAB[FOCUS], (" (" + " and ".join(LAB[k.split('-')[1]] for k in dup) + ": see above)") if dup else "")]
             + [row(k) for k in mn_])
    nexc = {k: pairs[k]["only_a_qualified"] + pairs[k]["only_b_qualified"] + pairs[k]["neither_qualified"] for k in pairs}
    exc = [f"{LAB[k.split('-')[0]]} vs.\\ {LAB[k.split('-')[1]]} {v}" for k, v in nexc.items() if v]
    supp.append(table(
        "table*", "Practical equivalence of the per-case mean wake losses at 6,030 evaluations. $\\overline{\\Delta L}$: mean "
        "difference of the case-mean wake losses, first minus second method (percentage points, pp; negative = first method "
        "better), over the $n$ cases in which both methods have at least 15 feasible runs out of 30 (cases excluded: %s); "
        "95\\%% and 90\\%% percentile bootstrap CIs (%s resamples of the cases, fixed seed). One equivalence margin "
        "$m=%.2f$~pp is applied to every comparison; it was set after the primary analysis. Two one-sided tests (TOST): "
        "the pair is equivalent at $m$ (``Eq.'') if the 90\\%% CI lies inside $(-m, m)$; $p_{\\rm TOST}$: bootstrap TOST "
        "$p$ (share of bootstrap means beyond $\\mp m$, add-one corrected; $\\le$: no bootstrap mean beyond the margin); "
        "$m_{\\min}$: smallest margin at which equivalence holds, $\\max(|{\\rm lo}_{90}|, |{\\rm hi}_{90}|)$, rounded up. "
        "Bayesian signed-rank test \\cite{Benavoli2017} on the case-mean differences with region of practical equivalence "
        "$(-m, m)$ (Dirichlet-process prior strength $s=%.1f$ at $z_0=0$, %s Monte Carlo samples, fixed seed): posterior "
        "probability that the first method is better ($P_{\\rm A}$), that the two are practically equivalent "
        "($P_{\\rm rope}$) and that the second is better ($P_{\\rm B}$)."
        % ("; ".join(exc) if exc else "none", f"{BOOT_N:,}".replace(",", "{,}"), EQ_MARGIN, BAYES_S,
           f"{BAYES_N:,}".replace(",", "{,}")),
        "tab:equivalence", "lcccccccccc",
        "Pair (first vs.\\ second) & $n$ & $\\overline{\\Delta L}$ (pp) & 95\\% CI & 90\\% CI & $p_{\\rm TOST}$ & $m_{\\min}$ & "
        "Eq. & $P_{\\rm A}$ & $P_{\\rm rope}$ & $P_{\\rm B}$", lines, size="\\scriptsize", sep="2.5pt", pos="!htb"))

    # --- spread and worst run, PSO-VNS vs PSO (68 cases, 6,030 evaluations, random initialization): per case the SD
    # of the wake loss (%) of the feasible runs and the worst (highest-loss) feasible run; cases in which both
    # methods qualify; Wilcoxon signed-rank (wil: |d| <= 1e-9 dropped) on the per-case values, unadjusted
    A_, B_ = "PSOBV", "PSOC"
    if {A_, B_} <= set(S.Algorithm):
        worst = G[G.Feasible].groupby(CASE + ["Algorithm"]).LossPct.max().rename("Worst").reset_index()
        Sw = S.merge(worst, on=CASE + ["Algorithm"], how="left")
        sd_d, cnt = case_mean_diffs(Sw, B_, A_, "LossSD")           # > 0: PSO-VNS has the smaller SD
        wr_d, _ = case_mean_diffs(Sw, B_, A_, "Worst")              # > 0: PSO-VNS has the better (lower-loss) worst run
        Q = Sw[Sw.Qualified.astype(bool)].pivot_table(index=CASE, columns="Algorithm", values=["LossSD", "Worst"])
        both = Q["LossSD"][[A_, B_]].dropna().index
        psd, rsd, _ = wil(sd_d.values)
        pwr, rwr, _ = wil(wr_d.values)
        sp = dict(pair=[A_, B_], n_cases=int(len(sd_d)), **cnt,
                  mean_sd_pp={A_: float(Q["LossSD"].loc[both, A_].mean()), B_: float(Q["LossSD"].loc[both, B_].mean())},
                  median_sd_pp={A_: float(Q["LossSD"].loc[both, A_].median()), B_: float(Q["LossSD"].loc[both, B_].median())},
                  sd_wins=int((sd_d > 1e-9).sum()), sd_losses=int((sd_d < -1e-9).sum()), sd_ties=int((sd_d.abs() <= 1e-9).sum()),
                  sd_p=psd, sd_rb=rsd,
                  mean_worst_loss_pct={A_: float(Q["Worst"].loc[both, A_].mean()), B_: float(Q["Worst"].loc[both, B_].mean())},
                  worst_wins=int((wr_d > 1e-9).sum()), worst_losses=int((wr_d < -1e-9).sum()),
                  worst_ties=int((wr_d.abs() <= 1e-9).sum()), worst_p=pwr, worst_rb=rwr,
                  note="per case: SD of the wake loss (%) of the feasible runs and loss of the worst feasible run, 30 runs, "
                       "6,030 evaluations; wins = cases in which PSO-VNS has the smaller SD / lower worst-run loss "
                       "(|d| <= 1e-9 = tie); Wilcoxon signed-rank on the per-case values (ties dropped), unadjusted; "
                       "rb > 0 = PSO-VNS smaller")
        sp["sd_direction"] = ("PSOBV_smaller" if sp["sd_p"] < 0.05 and sp["sd_rb"] > 0 else
                              "PSOC_smaller" if sp["sd_p"] < 0.05 else "no_significant_difference")
        sp["worst_direction"] = ("PSOBV_better" if sp["worst_p"] < 0.05 and sp["worst_rb"] > 0 else
                                 "PSOC_better" if sp["worst_p"] < 0.05 else "no_significant_difference")
        summary["spread"] = sp
        log(f"  SPREAD PSO-VNS vs PSO ({sp['n_cases']} cases): mean SD {sp['mean_sd_pp'][A_]:.4f} vs {sp['mean_sd_pp'][B_]:.4f} pp, "
            f"smaller SD in {sp['sd_wins']}/{sp['sd_losses']}/{sp['sd_ties']} (PSO-VNS/PSO/tie) cases, p={sp['sd_p']:.3g} "
            f"r_rb={sp['sd_rb']:+.2f}; worst run better in {sp['worst_wins']}/{sp['worst_losses']}/{sp['worst_ties']}, "
            f"p={sp['worst_p']:.3g} r_rb={sp['worst_rb']:+.2f}; mean worst loss {sp['mean_worst_loss_pct'][A_]:.4f} vs "
            f"{sp['mean_worst_loss_pct'][B_]:.4f} %")


def case_mean_wilcoxon(S, focus, others, val="Loss"):
    """Wilcoxon signed-rank on per-case mean wake losses (%), focus vs each method, Holm over methods.

    d = L(other) - L(focus) > 0 means focus better. Cases in which exactly one of the two methods
    qualifies (>= half of its runs feasible) are counted for the qualifying method with a difference
    larger than every observed one (consistent with the ranking rule); cases in which neither
    qualifies are dropped. The plain version (only cases where both qualify) is reported too.
    Also returned: wins / losses = cases entering the test in which the focus has the lower / higher
    mean loss (imputed cases included, exact ties |d| <= 1e-9 excluded), and the 95 % percentile bootstrap
    CI (BOOT_N resamples of the both-qualified cases, fixed seed BOOT_SEED) of mean_dloss_pp =
    mean of L(focus) - L(other) over the both-qualified cases (pp; negative = focus better).
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
                      wins=int((dfull > 1e-9).sum()), losses=int((dfull < -1e-9).sum()),
                      ci95_mean_dloss_pp=boot_ci(-d),
                      p_both_qualified=p2, rb_both_qualified=rb2, n_both_qualified=int(both.sum()))
    for key in ("p", "p_both_qualified"):
        bs = list(out)
        for b, h in zip(bs, holm([out[b][key] for b in bs]) if bs else []):
            out[b][key + "_holm"] = float(h)
    return out


CURVE_TOL = 0.01      # objective units (curves are stored with 3-4 decimals)


def curve_valid(M, ideal):
    """A checkpoint value (best feasible objective so far) must lie in [0, ideal]. Rarely (mostly PSO, whose
    particles sit exactly on a constraint boundary) the tracker of the experiment scripts accepted a layout
    that passes the 1e-6 feasibility tolerance but still carries a (large) penalty, so the recorded value is
    the penalized objective. Such points are treated as unknown (NaN, i.e. like not-yet-feasible)."""
    return (M >= -CURVE_TOL) & (M <= np.asarray(ideal, float).reshape(-1, *([1] * (np.ndim(M) - 1))) + CURVE_TOL)


def curves(sub):
    M = np.array([[float(v) for v in c.split(";")] for c in sub.Curve])
    return np.where(curve_valid(M, sub.Ideal.values), M, np.nan)


def curve_x(sub, ncp):
    b = int(sub.Budget.iloc[0])
    return np.arange(1, ncp + 1) * max(1, (b - 30) // 200)


def coords(s):
    return np.array([[float(v) for v in p.split()] for p in s.split(";")])


def wtl_str(c):
    return f"{c.get('W', 0)}/{c.get('T', 0)}/{c.get('L', 0)}"


def phase2_stat(D, alg, switch=None, validate=True):
    """Share (%) of the wake loss left at the switch call that the VNS phase of hybrid alg removes (feasible
    final runs whose best-feasible curve is finite at the switch), and the number of feasible final runs
    that were still infeasible at the switch (made feasible by phase 2)."""
    X = D[D.Algorithm == alg]
    Hy = X[X.Feasible]
    if not len(Hy):
        return None
    b = int(X.Budget.iloc[0])
    sw = switch or switch_call(alg, b)
    idx = switch_index(sw, b)
    gains, made_feasible, invalid, lsw = [], 0, 0, []
    for ideal, obj, cv in zip(Hy.Ideal.values, Hy.Objective.values, Hy.Curve.values):
        at = float(cv.split(";")[idx])
        if validate and np.isfinite(at) and not curve_valid(np.array([at]), [ideal])[0]:
            invalid += 1; continue                         # penalized value recorded at the switch (see curve_valid)
        if not np.isfinite(at):
            made_feasible += 1; continue
        l1, l2 = ideal - at, ideal - obj
        lsw.append(100 * l1 / ideal)
        if l1 > 1e-9:
            gains.append(100 * (l1 - l2) / l1)
    g = np.array(gains)
    Mall = curves(X)
    feas_sw = float(100 * np.isfinite(Mall[:, idx]).mean())
    return dict(mean=float(g.mean()) if len(g) else np.nan, median=float(np.median(g)) if len(g) else np.nan,
                n_runs_with_loss_at_switch=len(g), infeasible_at_switch_made_feasible=made_feasible,
                excluded_penalized_value_at_switch=invalid,
                final_infeasible=int(len(X) - len(Hy)), runs=int(len(X)), switch_call=int(sw),
                pct_runs_improved=float(100 * np.mean(g > 1e-9)) if len(g) else np.nan,
                mean_loss_at_switch_pct=float(np.mean(lsw)) if lsw else np.nan, feasible_at_switch_pct=feas_sw)


SMIN_M = 8 * 38.5          # minimum spacing 4D = 8R = 308 m (wflop_model.R = 38.5 m), as run_grid of mpce_experiments.py


def phase1_replay(D, n1=3015, npop=30):
    """Geometry-only replay of the Phase 1 of RS-VNS (square sampling) and RSD-VNS (disc sampling) with random
    initialization, 6,030 evaluations, split 0.5: n1 = round(0.5 x 6,030) = 3,015 samples per run, the first 30 of
    them the common initial population (init_pop: uniform in the square [-r, r]^2N). The seeded legacy stream of each
    run is replayed with np.random.RandomState(seed) (identical to the np.random.seed(seed) of RSVNS.__init__):
      RS-VNS : n1 x 2N draws uniform(-r, r) (the population, then one init_pop(1, 2N) layout per sample);
      RSD-VNS: the population (30 x 2N uniform(-r, r)), then per sample N angles uniform(0, 2 pi) followed by N radii
               r sqrt(uniform(0, 1)) (mpce_experiments.RSDVNS.optimize; uniform(0, h) = h x random_sample exactly).
    No objective call. A sample is feasible if every turbine is inside the circle and every pair is >= 4D apart
    (tolerance 1e-6 m, as run_grid). The geometry does not depend on the wind data set, so Data Set II repeats Data
    Set I. Returns, per method, the runs whose Phase-1 samples contain no feasible layout (whole Phase 1 incl. the
    initial population, and the sampled part only), the feasible share of the samples, and a verification against
    the stored convergence curves of D: the checkpoint at call 3,000 (index 99) is finite iff one of the first 3,000
    samples is feasible (the tracker records the best feasible objective)."""
    out = {}
    for alg in ("RSVNS", "RSDVNS"):
        X = D[(D.Algorithm == alg) & D.Dataset.isin(["1", "2"]) & (D.Budget == 6030) & (D.Init == "random")]
        if not len(X):
            continue
        feas_any, feas_samp, nfeas, nsamp, per_case = {}, {}, 0, 0, {}
        first3000, inside_any, spaced_any = {}, {}, {}
        for rad, n in sorted({(int(r), int(n)) for r, n in zip(X.Radius, X.Turbines)}):
            iu = np.triu_indices(n, 1)
            k_case = 0
            for sd in sorted(set(X[(X.Radius == rad) & (X.Turbines == n)].Seed)):
                rs = np.random.RandomState(int(sd))
                P0 = rs.uniform(-rad, rad, (npop, 2 * n))
                if alg == "RSVNS":
                    P1 = rs.uniform(-rad, rad, (n1 - npop, 2 * n))
                else:
                    U = rs.random_sample((n1 - npop, 2, n))
                    ang, rr = 2 * np.pi * U[:, 0], rad * np.sqrt(U[:, 1])
                    P1 = np.stack([rr * np.cos(ang), rr * np.sin(ang)], -1).reshape(n1 - npop, 2 * n)
                XY = np.vstack([P0, P1]).reshape(n1, n, 2)
                inside = np.sqrt((XY ** 2).sum(-1)).max(-1) <= rad + 1e-6
                dd = np.sqrt(((XY[:, :, None] - XY[:, None]) ** 2).sum(-1))[:, iu[0], iu[1]].min(-1)
                fz = inside & (dd >= SMIN_M - 1e-6)
                feas_any[(rad, n, sd)] = bool(fz.any()); feas_samp[(rad, n, sd)] = bool(fz[npop:].any())
                inside_any[(rad, n, sd)] = bool(inside.any()); spaced_any[(rad, n, sd)] = bool((dd >= SMIN_M - 1e-6).any())
                first3000[(rad, n, sd)] = bool(fz[:3000].any())
                nfeas += int(fz.sum()); nsamp += n1
                k_case += int(not fz.any())
            per_case[f"{rad}-{n}"] = k_case
        chk = agree = 0
        for ds, rad, n, sd, cv in zip(X.Dataset, X.Radius, X.Turbines, X.Seed, X.Curve):
            v = cv.split(";")[99]
            fin = v != "nan" and np.isfinite(float(v))
            chk += 1; agree += int(fin == first3000[(int(rad), int(n), sd)])
        runs = {(r_, n_, s_) for r_, n_, s_ in zip(X.Radius.astype(int), X.Turbines.astype(int), X.Seed)}
        nds = X.Dataset.nunique()
        out[alg] = dict(runs=int(len(X)), samples_per_run=n1, initial_population=npop,
                        runs_no_feasible_sample=int(nds * sum(not feas_any[k] for k in runs)),
                        runs_no_feasible_sample_sampled_part=int(nds * sum(not feas_samp[k] for k in runs)),
                        pct_feasible_samples=float(100 * nfeas / nsamp),
                        runs_no_feasible_sample_per_case_dsI=per_case,
                        runs_no_feasible_but_inside_sample=int(nds * sum(not feas_any[k] and inside_any[k] for k in runs)),
                        runs_no_feasible_but_spaced_sample=int(nds * sum(not feas_any[k] and spaced_any[k] for k in runs)),
                        runs_no_inside_sample=int(nds * sum(not inside_any[k] for k in runs)),
                        runs_no_feasible_first3000=int(nds * sum(not first3000[k] for k in runs)),
                        _feas_any={f"{r_}-{n_}-{s_}": feas_any[(r_, n_, s_)] for r_, n_, s_ in runs}, _nds=int(nds),
                        verified_runs=chk, verified_agree=agree, verification_ok=bool(chk > 0 and agree == chk),
                        note="Phase-1 samples = 30 common initial layouts (square) + 2,985 samples (square for RS-VNS, "
                             "disc for RSD-VNS); counts over both data sets (identical geometry and seeds)")
        log(f"  Phase-1 replay {LAB[alg]}: {out[alg]['runs_no_feasible_sample']} of {out[alg]['runs']} runs without a feasible "
            f"sample ({out[alg]['runs_no_feasible_sample_sampled_part']} ignoring the initial population); feasible samples "
            f"{out[alg]['pct_feasible_samples']:.2f}%; curve check {agree}/{chk}")
    # R3-6 (review round 2): how many RS-VNS runs without a feasible sample can the square be blamed for? Paired over
    # the same runs (case, seed; same initial population): disc sampling (RSD-VNS) rescues a run iff it has a feasible
    # sample where square sampling has none; runs without a feasible sample under both are limited by the spacing
    # constraint (packing density), not by the square.
    if "RSVNS" in out and "RSDVNS" in out:
        fa, fb = out["RSVNS"].pop("_feas_any"), out["RSDVNS"].pop("_feas_any")
        nds = out["RSVNS"].pop("_nds"); out["RSDVNS"].pop("_nds")
        common = sorted(set(fa) & set(fb))
        dec = dict(runs=int(nds * len(common)),
                   rs_none_rsd_some=int(nds * sum(not fa[k] and fb[k] for k in common)),
                   rs_none_rsd_none=int(nds * sum(not fa[k] and not fb[k] for k in common)),
                   rs_some_rsd_none=int(nds * sum(fa[k] and not fb[k] for k in common)),
                   note="paired over runs (case, seed, both data sets): rs_none_rsd_some = RS-VNS runs without a feasible "
                        "Phase-1 sample that disc sampling would rescue (the most the square can explain); rs_none_rsd_none "
                        "= runs without a feasible sample even when sampling in the disc (spacing / packing density)")
        out["square_vs_disc"] = dec
        log(f"  Phase-1 replay, square vs disc (paired): of {out['RSVNS']['runs_no_feasible_sample']} RS-VNS runs without a "
            f"feasible sample, {dec['rs_none_rsd_some']} have one with disc sampling and {dec['rs_none_rsd_none']} have none "
            f"either way ({dec['rs_some_rsd_none']} runs feasible with the square only); "
            f"{out['RSVNS']['runs_no_feasible_but_inside_sample']} of them have a sample with all turbines inside the circle")
    else:
        for v in out.values():
            v.pop("_feas_any", None); v.pop("_nds", None)
    return out


def split_section(R6, base, tabs, key, label, primary=True, supp=None):
    """Budget-split table of hybrid `base` (omega = 25 / 50 / 75 % of the calls for phase 1, if its split runs
    exist; omega = 90 % if its runs exist (Phase 6, experiment omega90: PSOBV90, 5,430 PSO evaluations); and
    omega = 100 %, i.e. the phase-1 swarm alone at the same budget and seeds: PSO for PSO-VNS).
    Run-level tests: per case, seed-paired Wilcoxon of 50 % vs 25 %, 50 % vs 75 %, [50 % vs 90 %,] 50 % vs 100 %,
    75 % vs 100 % [and 90 % vs 75 %], Holm-adjusted over these comparisons of the case (one family per case: four
    without omega = 0.9, six with it). Case level: Wilcoxon on the 12 per-case mean losses (case_mean_wilcoxon;
    unadjusted) for the same pairs; the pairs without omega = 0.9 go to out["case_mean"] (keys as before, re-tested
    by mpce_inference_extra.py X01), those with omega = 0.9 to out["case_mean_omega90"]."""
    ids = (f"{base}25", base, f"{base}75")
    i90 = f"{base}90"
    if base not in PHASE1:
        log(f"  SKIPPED: {LAB[base]} is not a two-phase hybrid (no budget split)")
        return None
    p1c = PHASE1[base]                               # omega = 1: the phase-1 swarm alone (same budget, same seeds)
    Sp = R6[R6.Algorithm.isin(ids + (i90, p1c))]
    Sp = Sp[[(d, r, n) in SPLITCASES for d, r, n in zip(Sp.Dataset, Sp.Radius, Sp.Turbines)]]
    if not {ids[0], ids[2]} <= set(Sp.Algorithm):
        log(f"  SKIPPED: no split runs for {LAB[base]} ({ids[0]} / {ids[2]} not in the data)"
            + {"SSABV": " -- mpce_ssasplit not present (experiment dropped)", "PSOBV": " -- mpce_psosplit missing"}.get(base, ""))
        return None
    has90 = i90 in set(Sp.Algorithm)
    has100 = p1c in set(Sp.Algorithm)
    idsN = ids + ((i90,) if has90 else ()) + ((p1c,) if has100 else ())   # in order of omega
    Sp = Sp[Sp.Algorithm.isin(idsN)]
    # (tag, a, b): run-level comparisons of one case, one Holm family; p<tag>_holm / rb<tag> in the per-case rows
    comps = [("25", base, ids[0]), ("75", base, ids[2])]
    comps += [("90", base, i90)] if has90 else []
    comps += [("100", base, p1c), ("75v100", ids[2], p1c)] if has100 else []
    comps += [("90v75", i90, ids[2])] if has90 else []
    shares = {ids[0]: 0.25, base: 0.5, ids[2]: 0.75, i90: 0.9, p1c: 1.0}
    lines, srows = [], []
    for (ds, r, n), s in Sp.groupby(CASE):
        if s.Algorithm.nunique() < len(idsN):
            continue
        piv = s.assign(S=goodness(s)).pivot_table(index="Seed", columns="Algorithm", values="S")
        fm = s[s.Feasible].groupby("Algorithm").Objective.mean(); fc = s.groupby("Algorithm").Feasible.sum()
        nr = s.groupby("Algorithm").size()
        lmn = s[s.Feasible].groupby("Algorithm").LossPct.mean()
        ps, rbs = [], []
        for _, a, b in comps:
            p, rb, _ = wil((piv[a] - piv[b]).dropna().values); ps.append(p); rbs.append(rb)
        ph = holm(ps)
        cells = []
        for a in idsN:
            v = fm.get(a, np.nan)
            c = "--" if not np.isfinite(v) else f"{v:.1f}"
            if fc.get(a, 0) < nr.get(a, 0):
                c += f"$^{{{int(fc.get(a, 0))}}}$"
            cells.append(c)
        best = max(idsN, key=lambda a: fm.get(a, -np.inf) if fc.get(a, 0) >= np.ceil(nr.get(a, 0) / 2) else -np.inf)
        row = dict(case=f"{ds}-{r}-{n}", best=best, loss={a: float(lmn.get(a, np.nan)) for a in idsN})
        for (tag, _, _), h, rb in zip(comps, ph, rbs):
            row[f"p{tag}_holm"] = float(h); row[f"rb{tag}"] = rb
        srows.append(row)
        lines.append(f"{'I' if ds == '1' else 'II'} & {r} & {n} & " + " & ".join(cells) + " & "
                     + " & ".join(fmt_p(x) for x in ph) + " \\\\")
    p1 = LAB[p1c]
    pre = "" if primary else f"Secondary analysis ({LAB[base]}, not the proposed method). "
    SR = pd.DataFrame(srows)
    # compact main-text table: mean loss over the cases, average rank among the splits, W/T/L of 50% vs each
    Ssp = case_stats(Sp, list(idsN))
    Rsp = rank_matrix(Ssp, list(idsN))
    avg_sp = Rsp.mean(axis=0).to_dict()
    mloss = {a: float(np.mean([c["loss"][a] for c in srows])) for a in idsN}
    def wtl_split(tag):
        pcol, rbcol = f"p{tag}_holm", f"rb{tag}"
        if pcol not in SR:
            return None
        w = int(((SR[pcol] < 0.05) & (SR[rbcol] > 0)).sum()); l = int(((SR[pcol] < 0.05) & (SR[rbcol] < 0)).sum())
        return dict(W=w, T=int(len(SR) - w - l), L=l)
    W_ = {tag: wtl_split(tag) for tag, _, _ in comps}
    w25, w75, w90, w100, w75v100, w90v75 = (W_.get(t) for t in ("25", "75", "90", "100", "75v100", "90v75"))
    lower = {a: int(sum(c["loss"][a] < c["loss"][base] - 1e-12 for c in srows)) for a in idsN if a != base}
    cm, cm90 = {}, {}
    for _, a, b in comps:                    # case-mean Wilcoxon (12 cases), a vs b, unadjusted
        x = case_mean_wilcoxon(Ssp, a, [b])[b]
        (cm90 if i90 in (a, b) else cm)[f"{a}-{b}"] = dict(p=x["p"], wins=x["wins"], losses=x["losses"],
                                                          mean_dloss_pp=x["mean_dloss_pp"], n_cases=x["n_cases"])
    wcol = {ids[0]: w25, ids[2]: w75, i90: w90, p1c: w100}
    cl = [f"${shares[a]:g}${' (' + p1 + ' alone)' if a == p1c else ''} & {mloss[a]:.3f} & {avg_sp[a]:.2f} & "
          + ("--" if a == base else wtl_str(wcol[a])) + " \\\\" for a in idsN]
    nsh = {3: "Three", 4: "Four", 5: "Five"}[len(idsN)]
    ncw = {2: "Two", 4: "Four", 3: "Three", 6: "Six"}.get(len(comps), str(len(comps)))
    om = ", ".join(f"${shares[a]:g}$" for a in idsN if a != p1c)
    if primary:
        foot = []
        if has100:
            foot.append("\\multicolumn{4}{l}{$\\omega=0.75$ vs.\\ $\\omega=1$ (%s alone): %s (W/T/L from the $\\omega=0.75$ side)}" % (p1, wtl_str(w75v100)))
        if has90:
            foot.append("\\multicolumn{4}{l}{$\\omega=0.9$ vs.\\ $\\omega=0.75$: %s (W/T/L from the $\\omega=0.9$ side)}" % wtl_str(w90v75))
        foot = "\\\\\n".join(foot) if foot else None
        foot = [foot] if foot else None
        tabs[key] = table("table", f"Budget Split of {LAB[base]} ({len(SR)} Cases, 30 Seeds)",
                          label, "cccc", f"$\\omega$ & Mean loss (\\%) & Avg.\\ rank & $0.5$ vs.\\ $\\omega$ (W/T/L)", cl, sep="4pt", foot=foot,
                          note=f"$\\omega$: share of the 6,030 evaluations given to {p1} (${chr(92)}omega={om[1:]}"
                               f"{'; $' + chr(92) + 'omega=1$: ' + p1 + ' alone, same seeds' if has100 else ''}); cases: middle and "
                               f"largest $N$ of each farm and data set. Avg.\\ rank among the {nsh.lower()} settings. W/T/L: cases in "
                               f"which the default $\\omega=0.5$ is significantly better / not different / worse (run-level "
                               f"Wilcoxon, Holm-adjusted over the {ncw.lower()} comparisons of each case).")
        label = label + "-cases"
    heads = ["$\\omega=0.25$", "$0.5$", "$0.75$"] + (["$0.9$"] if has90 else []) + (["$1$"] if has100 else [])
    lab_p = {"25": "0.5/0.25", "75": "0.5/0.75", "90": "0.5/0.9", "100": "0.5/1", "75v100": "0.75/1", "90v75": "0.9/0.75"}
    pheads = [f"$p$ ({lab_p[t]})" for t, _, _ in comps]
    pct = ", ".join(f"{int(round(100 * shares[a]))}\\%" for a in idsN if a != p1c)
    (supp.append if supp is not None else (lambda t: tabs.__setitem__(key, t)))(table("table*" if has100 else "table", pre + f"Sensitivity of {LAB[base]} to the budget split between the {p1} and VNS phases ({pct} of the 6,030 calls for {p1}{'; 100' + chr(92) + '%: ' + p1 + ' alone, same seeds' if has100 else ''}): mean benchmark objective of the feasible runs (superscript: feasible runs when fewer than 30) and Wilcoxon signed-rank $p$ (30 seed-paired runs), Holm-adjusted over the {len(comps)} comparisons of each case.",
                      label, "ccc" + "c" * (len(heads) + len(pheads)), "DS & $r$ & $N$ & " + " & ".join(heads + pheads), lines, sep="2.5pt" if len(comps) <= 4 else "1.5pt", pos="!htb"))
    out = dict(method=base, primary=primary, cases=srows, best_count=SR.best.value_counts().to_dict(),
               n_cases=int(len(SR)), settings={a: shares[a] for a in idsN}, omega1_method=p1c if has100 else None,
               omega90_method=i90 if has90 else None, n_settings=len(idsN), n_comparisons=len(comps),
               holm_family="per case, over the %d run-level comparisons %s" % (len(comps), ", ".join(f"{a} vs {b}" for _, a, b in comps)),
               avg_rank=avg_sp, mean_loss=mloss, wtl_50_vs_25=w25, wtl_50_vs_75=w75, wtl_50_vs_90=w90,
               wtl_50_vs_100=w100, wtl_75_vs_100=w75v100, wtl_90_vs_75=w90v75, case_mean=cm, case_mean_omega90=cm90,
               n_lower_loss_than_50=lower,
               sig_vs25=int((SR.p25_holm < 0.05).sum()), sig_vs75=int((SR.p75_holm < 0.05).sum()),
               sig_vs25_50better=int(((SR.p25_holm < 0.05) & (SR.rb25 > 0)).sum()),
               sig_vs75_50better=int(((SR.p75_holm < 0.05) & (SR.rb75 > 0)).sum()))
    if has90:
        # omega = 0.9 against omega = 0.75 per case: lower mean loss (feasible runs), and the best setting overall
        out["n_lower_loss_90_than_75"] = int(sum(c["loss"][i90] < c["loss"][ids[2]] - 1e-12 for c in srows))
        out["n_higher_loss_90_than_75"] = int(sum(c["loss"][i90] > c["loss"][ids[2]] + 1e-12 for c in srows))
        out["best_mean_loss"] = min(mloss, key=mloss.get)
        out["best_avg_rank"] = min(avg_sp, key=avg_sp.get)
    log(f"  {LAB[base]}: best split counts {out['best_count']}; significant vs 25%: {out['sig_vs25']}, vs 75%: {out['sig_vs75']}"
        + (f"; 75% vs 100% ({p1}): {wtl_str(w75v100)}" if has100 else "")
        + (f"; 90% vs 75%: {wtl_str(w90v75)}, 50% vs 90%: {wtl_str(w90)}" if has90 else "")
        + f"; avg ranks {', '.join(f'{a} {v:.2f}' for a, v in avg_sp.items())}; mean loss {', '.join(f'{a} {v:.4f}' for a, v in mloss.items())}")
    return out


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
    h = [plt.Line2D([], [], color=COL[a], ls=ls(a), lw=lw(a), marker=mk(a), ms=ms(a, 4), label=LAB[a].replace("\\%", "%"))
         for a in algs]
    fig.legend(handles=h, loc="lower center", ncol=ncol or (len(algs) if len(algs) <= 7 else int(np.ceil(len(algs) / 2))),
               bbox_to_anchor=(0.5, y), fontsize=7.5)


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
        ax.plot(x, med, color=COL[a], ls=ls(a), lw=lw(a), zorder=zo(a))
        if band:
            q1, q3 = median_curve(V, loss, 25), median_curve(V, loss, 75)
            ax.fill_between(x, q1, q3, color=COL[a], alpha=0.12, lw=0)
        kl = np.where(np.isfinite(med))[0]
        if len(kl):
            ax.plot(x[kl[-1]], med[kl[-1]], marker=mk(a), color=COL[a], ms=ms(a, 4), zorder=zo(a))


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
    for exp in ("rsvns", "psoc", "psobv", "slsqp", "psosplit", "omega90", "rsdisc", "ssasplit", "hr16new", "feas", "b30k", "b120k",
                "iea16", "iea36", "feasp", "b30kp", "b120kp", "iea16p", "iea36p"):
        df, st = read_shards(exp, data_dir, partial)
        avail[f"mpce_{exp}"] = st
        if df is not None:
            new[exp] = df
    # Horns Rev 1 rerun with the corrected direction binning (experiment hrfix, 24 shards). As soon as ANY
    # hrfix shard exists, every Horns Rev row (Dataset == "HR") of every other source -- fresh_hr16 / vhr16 /
    # bhr16, hr16new, psobv, feas / feasx / feasp, b30k(p), b120k(p), ... -- is dropped (old binning, not
    # comparable) and replaced by the hrfix rows (run_hr output: Budget, Init, Objective = AEP, Ideal, ...).
    # Incomplete hrfix shards are used as they are (status PARTIAL): mixing the two models is never done, and
    # mpce_numbers.py marks the Horns Rev macros pending until all shards are present.
    hrfix = None
    if glob.glob(os.path.join(data_dir, "mpce_hrfix_s*of*.csv")):
        hrfix, st = read_shards("hrfix", data_dir, True)
        if st["status"] != "complete":
            st["status"] = "PARTIAL (used; old HR runs dropped)"
        avail["mpce_hrfix"] = st
        if hrfix is not None:
            hrfix = hrfix[hrfix.Dataset == "HR"]
    else:
        avail["mpce_hrfix"] = dict(status="missing", shards="0", rows=0)
    if hrfix is not None:
        n_old = int((B.Dataset == "HR").sum()) + sum(int((d.Dataset == "HR").sum()) for d in new.values())
        B = B[B.Dataset != "HR"]
        new = {k: v[v.Dataset != "HR"].copy() for k, v in new.items()}
        fallbacks.append(f"Horns Rev 16: all runs from mpce_hrfix (corrected direction binning, {len(hrfix)} runs); "
                         f"{n_old} runs of the old model dropped")
    parts = [B]
    # SLSQP rerun replaces the old-platform SLSQP runs on the 68 cases
    if "slsqp" in new:
        parts[0] = B[~((B.Algorithm == "SLSQP") & B.Dataset.isin(["1", "2"]))]
    else:
        fallbacks.append("SLSQP (68 cases): mpce_slsqp missing -> using the old-platform SLSQP runs of fresh_grid.csv")
        m = (B.Algorithm == "SLSQP") & B.Dataset.isin(["1", "2"])
        B.loc[m, "Source"] = "FALLBACK fresh_grid SLSQP"
    if hrfix is None:
        fallbacks.append("SLSQP (Horns Rev 16, 6,030 calls): taken from fresh_hr16.csv (old model, replaced by mpce_hrfix)")
        fallbacks.append("Horns Rev 16: mpce_hrfix missing -> old runs (direction binning before commit 7676da9)")
    for exp in ("rsvns", "psoc", "psobv", "slsqp", "psosplit", "omega90", "rsdisc", "ssasplit", "hr16new"):
        if exp in new:
            parts.append(new[exp])
    if hrfix is not None:
        parts.append(hrfix)
    A6 = pd.concat(parts, ignore_index=True)
    # PSOC fallback: old-settings PSO relabelled (development only)
    hr_psoc = (hrfix is not None and (hrfix.Algorithm == "PSOC").any()) or \
              ("hr16new" in new and (new["hr16new"].Algorithm == "PSOC").any())
    for dom, has in ((["1", "2"], "psoc" in new), (["HR"], hr_psoc)):
        if not has:
            f = A6[(A6.Algorithm == "PSO") & A6.Dataset.isin(dom)].copy()
            f["Algorithm"] = "PSOC"; f["Source"] = "FALLBACK old PSO"
            A6 = pd.concat([A6, f], ignore_index=True)
            fallbacks.append(f"PSOC ({'68 cases' if dom[0] == '1' else 'Horns Rev 16'}): constriction-PSO runs missing "
                             f"-> using the old-settings PSO runs relabelled as PSOC (DEVELOPMENT ONLY)")
    new["pso_old"] = A6[(A6.Algorithm == "PSO") & (A6.Source != "FALLBACK old PSO")].copy()   # w = 0.7, c1 = c2 = 2
    A6 = A6[A6.Algorithm != "PSO"]
    # nine-method experiments + their PSO-VNS-only arms (PSO-VNS becomes the tenth method, M10)
    extra = [new[e] for e in ("feas", "b30k", "b120k", "iea16", "iea36", "feasp", "b30kp", "b120kp", "iea16p", "iea36p")
             if e in new]
    ALL = pd.concat([A6] + extra, ignore_index=True)
    ALL = ALL.drop_duplicates(KEY, keep="first").reset_index(drop=True)
    return ALL, avail, fallbacks, new


def common_seeds(ALL):
    """Keep, within every case x budget x initialization, only the seeds that every method of the group has
    (split variants *25 / *75 / *90 are not used to define the intersection)."""
    grp = CASE + ["Budget", "Init"]
    core = ALL[~ALL.Algorithm.str.contains(r"(?:25|75|90)$")]
    sets = core.groupby(grp + ["Algorithm"]).Seed.agg(frozenset)
    keep = sets.groupby(level=list(range(len(grp)))).agg(lambda v: frozenset.intersection(*v))
    k = pd.Series([keep.get(tuple(t), frozenset()) for t in ALL[grp].itertuples(index=False)], index=ALL.index)
    m = np.array([sd in ks for sd, ks in zip(ALL.Seed, k)])
    return ALL[m].reset_index(drop=True), int((~m).sum())


# ------------------------------------------------------------------ LaTeX helpers
def tnote(text, width="\\columnwidth"):
    """Table note below the tabular (R1-20: captions are short noun phrases; definitions, test families and case
    sets go here, in scriptsize, instead of the all-caps IEEE caption)."""
    return "\\par\\vspace{2pt}\\parbox{%s}{\\scriptsize %s}\n" % (width, text)


def table(env, caption, label, spec, header, lines, size="\\scriptsize", sep="3pt", resize=False, pos="!t", foot=None, note=None):
    body = "\n".join(lines)
    ft = "" if not foot else "\n" + "\n".join(foot)
    tab = f"\\begin{{tabular}}{{{spec}}}\n\\toprule\n{header} \\\\\n\\midrule\n{body}\n\\bottomrule{ft}\n\\end{{tabular}}"
    if resize:
        tab = "\\resizebox{\\textwidth}{!}{%\n" + tab + "}"
    nt = "" if not note else tnote(note, "\\textwidth" if env.endswith("*") else "\\columnwidth")
    return (f"\\begin{{{env}}}[{pos}]\n\\centering\n\\caption{{{caption}}}\n\\label{{{label}}}\n"
            f"{size}\\setlength{{\\tabcolsep}}{{{sep}}}\n{tab}\n{nt}\\end{{{env}}}\n")


# ------------------------------------------------------------------ main
def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--focus", choices=FOCUS_CHOICES, default="PSOBV")
    ap.add_argument("--data-dir", default=HERE)
    ap.add_argument("--out-dir", default=HERE)
    ap.add_argument("--fig-dir", default=os.path.join(HERE, "..", "figures_mpce"))
    ap.add_argument("--procs", type=int, default=4)
    ap.add_argument("--partial", action="store_true")
    ap.add_argument("--common-seeds", action="store_true")
    ap.add_argument("--skip-robust", action="store_true")
    ap.add_argument("--only-robust", action="store_true")
    args = ap.parse_args(argv)
    set_focus(args.focus)
    LOG.clear()
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
    dropped = 0
    if args.common_seeds:
        ALL, dropped = common_seeds(ALL)
        log(f"  --common-seeds: {dropped} runs dropped (seeds not shared by every method of a case/budget/init group)")
    fl = LAB[FOCUS]
    log(f"  FOCUS (proposed method): {FOCUS} = {fl}")
    log("=" * 78)
    summary = dict(generated=time.strftime("%Y-%m-%d %H:%M:%S"), focus=FOCUS, focus_label=fl,
                   options=dict(partial=args.partial, common_seeds=args.common_seeds),
                   colours=COL, data_availability=avail, fallbacks=fallbacks,
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
    log(f"  W/T/L vs {fl}: " + ", ".join(f"{LAB[b]} {wtl_str(summary['main']['wtl'][b])}" for b in others))

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
                              "$N$ & Wake-free & " + " & ".join(LAB[a] for a in MAINP), lines, resize=True, sep="2pt"))

    # --- minimum-spacing re-optimization (earlier runs with the original SSA / LX-SSA code; analyze_authors_runs.py)
    #     and constructible turbine counts (packing_capacity.py), copied into the supplement after the per-case tables
    try:
        at = open(os.path.join(HERE, "authors_tables.tex")).read()
        m_ = re.search(r"\\begin\{table\*\}\[!t\]\n(?:(?!\\end\{table\*\}).)*?\\label\{tab:spacing-authors\}.*?\\end\{table\*\}", at, re.S)
        if m_:
            t_ = m_.group(0).replace("Minimum-spacing sensitivity with the authors' LX-SSA and SSA code", "Minimum-spacing re-optimization with the original LX-SSA and SSA code") \
                .replace("DS & Radius &", "DS & $r$ (m) &").replace("[!t]", "[!htb]", 1)
            supp.append(t_)
        cap = pd.read_csv(os.path.join(HERE, "packing_capacity.csv"))
        cl = [f"{r} & " + " & ".join(str(int(cap[(cap.Radius == r) & (cap.Spacing == sp_)].MaxNConstructed.iloc[0])) for sp_ in ("4D", "5D", "6D"))
              + f" & {RADII[r]} \\\\" for r in RADII]
        supp.append(table("table", "Largest turbine count for which a layout satisfying the circular boundary and the minimum spacing was constructed (multi-start packing search; constructive lower bound), compared with the largest tested $N$.",
                          "tab:capacity", "ccccc", "$r$ (m) & $4D$ (308 m) & $5D$ (385 m) & $6D$ (462 m) & Largest tested $N$", cl, size="\\footnotesize", pos="!htb"))
    except Exception as e:                                   # pragma: no cover
        log(f"  spacing tables not copied: {e}")

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
    nword = {6: "Six", 7: "Seven", 8: "Eight", 9: "Nine"}.get(len(others), str(len(others)))
    tabs["wtl"] = table("table", "Run-Level Outcome of %s Against Each Method (W/T/L over the %d Cases)" % (fl, C[CASE].drop_duplicates().shape[0]),
                        "tab:wtl", "llc" + "c" * len(others),
                        "DS & $r$ (m) & Cases & " + " & ".join(HEAD2.get(b, LAB[b]) for b in others), lines, sep="1.6pt",
                        note="W/T/L: cases in which %s is significantly better / not different / worse (two-sided Wilcoxon "
                             "signed-rank test on 30 seed-paired runs, $\\alpha=0.05$, Holm-adjusted over the %s comparisons of "
                             "%s in each case, a different family from Table~\\ref{tab:ablation}; an infeasible run ranks below "
                             "every feasible run). $\\tilde r_{\\rm rb}$: median rank-biserial correlation (positive = %s "
                             "better; zero differences dropped, $r_{\\rm rb}=0$ if all are zero)." % (fl, nword.lower(), fl, fl))
    # W/T/L and mean loss difference by farm size (N >= 10 vs N < 10), for the text
    big = C.Turbines >= 10
    Sq = S[S.Qualified].pivot_table(index=CASE, columns="Algorithm", values="Loss")
    bign = Sq.index.get_level_values("Turbines") >= 10
    by_n = {}
    for b in others:
        x = C[C.Baseline == b]
        d = (Sq[FOCUS] - Sq[b]) if b in Sq else pd.Series(dtype=float)
        by_n[b] = dict(wtl_large=x[x.Turbines >= 10].Outcome.value_counts().to_dict(),
                       wtl_small=x[x.Turbines < 10].Outcome.value_counts().to_dict(),
                       n_large=int(x[x.Turbines >= 10][CASE].drop_duplicates().shape[0]),
                       n_small=int(x[x.Turbines < 10][CASE].drop_duplicates().shape[0]),
                       mean_dloss_pp_large=float(d[bign].mean()) if len(d) else np.nan,
                       mean_dloss_pp_small=float(d[~bign].mean()) if len(d) else np.nan,
                       mean_abs_dloss_pp_small=float(d[~bign].abs().mean()) if len(d) else np.nan,
                       max_abs_dloss_pp_small=float(d[~bign].abs().max()) if len(d) else np.nan)
    summary["main"]["by_n"] = by_n
    # exploratory (post hoc) N >= 10 subgroups of focus vs its phase-1 swarm (PSO for PSO-VNS): case-mean
    # Wilcoxon (same procedure as the main case-mean test, unadjusted) on the subgroup's cases, and significant
    # run-level wins / losses (main-table family: Holm over the 7 comparisons of each case)
    ref1 = PHASE1.get(FOCUS)
    if ref1 in others:
        sub_n = {}
        for tag, dsl in (("Large", ("1", "2")), ("dsILarge", ("1",)), ("dsIILarge", ("2",))):
            Sg = S[S.Dataset.isin(dsl) & (S.Turbines >= 10)]
            x = case_mean_wilcoxon(Sg, FOCUS, [ref1])[ref1]
            cr = C[(C.Baseline == ref1) & C.Dataset.isin(dsl) & (C.Turbines >= 10)].Outcome.value_counts()
            sub_n[tag] = dict(datasets=list(dsl), n_cases=int(Sg[CASE].drop_duplicates().shape[0]), p=x["p"], rb=x["rb"],
                              wins=x["wins"], losses=x["losses"], mean_dloss_pp=x["mean_dloss_pp"],
                              ci95_mean_dloss_pp=x["ci95_mean_dloss_pp"],
                              run_level=dict(W=int(cr.get("W", 0)), T=int(cr.get("T", 0)), L=int(cr.get("L", 0))))
        summary["main"]["subgroup_vs_phase1"] = dict(
            phase1=ref1, groups=sub_n,
            note="exploratory, post hoc split at N >= 10; case-mean Wilcoxon unadjusted; run-level W/T/L from the "
                 "main-table family (Holm over the 7 comparisons of the focus in each case)")
    # imputation in the main case-mean tests (focus vs the seven methods)
    summary["main"]["case_mean_imputation"] = dict(
        imputed={b: int(CW[b]["only_focus_qualified"] + CW[b]["only_other_qualified"]) for b in others},
        dropped={b: int(CW[b]["neither_qualified"]) for b in others},
        imputed_total=int(sum(CW[b]["only_focus_qualified"] + CW[b]["only_other_qualified"] for b in others)),
        dropped_total=int(sum(CW[b]["neither_qualified"] for b in others)),
        note="a case in which exactly one of the two methods has >= 15 feasible runs enters the case-mean test with "
             "a difference larger than every observed one in favour of that method (imputed); cases in which "
             "neither has are dropped; counts summed over the comparisons of the main table")
    # PSO-VNS vs PSO has two different run-level tallies in the paper (R3 item 5): Table wtl counts
    # significant cases with Holm over the 7 comparisons of the focus with the other methods of the main
    # comparison in each case; Table ablation uses Holm over the planned ablation contrasts of each case.
    # The raw per-case p-values are identical; only the Holm family (and thus the adjusted p) differs.
    summary["holm_families"] = dict(
        main_wtl=f"per case: {fl} vs each of the other {len(others)} methods of the main comparison (Table wtl, tab:friedman68 p_W over the {len(others)} methods)",
        ablation="per case: the planned contrasts of Table ablation (run level); case-mean p_W of Table ablation: Holm over the same contrasts",
        baseline="per case: the hybrid vs each of the other 7 methods of the previous study's pool (Table baseline)",
        split="per case: the four run-level comparisons of the split table (50 vs 25, 50 vs 75, 50 vs 100, 75 vs 100)",
        hr16="per setting: the focus vs each other method",
        feasinit="per case: feasible vs random initialization, over the methods of the case",
        why_psovns_vs_pso_differs="The main-table tally (Holm over 7 comparisons) and the ablation tally (Holm over the "
                                  "ablation contrasts) use the same per-case raw p-values of PSO-VNS vs PSO but different "
                                  "Holm families, so the adjusted p and hence W/T/L differ; both are descriptive tallies.")

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
    supp.append(table(
        "table*", "Case-level analysis over the %d benchmark cases. Ranks of the mean feasible benchmark objective within each case (1 = best); a method with fewer than 15 feasible runs out of 30 in a case is ranked below all other methods, by its number of feasible runs (``$<$15'': number of such cases). Friedman $\\chi^2_F=%.1f$ (%d d.f.), $p=%s$; Iman--Davenport $F_F=%.1f$, $p=%s$. $p_z$: Holm-adjusted $p$ of the average-rank $z$ test against %s (best-ranked method: %s). Because mean-rank post hoc tests depend on the pool of compared methods \\cite{Benavoli2016}, pairwise two-sided Wilcoxon signed-rank tests on the per-case mean wake losses (\\%%) are also given: $p_W$ (Holm-adjusted over the %d comparisons), matched-pairs rank-biserial correlation $r_{\\rm rb}$ (positive = %s better) and mean difference $\\overline{\\Delta L}$ of the wake loss (percentage points, %s minus method, over the cases where both methods have at least 15 feasible runs). ``Best'': cases in which the method alone ranks first; ``Feas.'': percentage of feasible runs."
        % (FR["n_cases"], FR["chi2"], len(MAINP) - 1, fmt_p(FR["p"]).strip("$"), FR["iman_davenport"],
           fmt_p(FR["iman_davenport_p"]).strip("$"), fl, LAB[FR["best_ranked"]], len(others), fl, fl),
        "tab:friedman68-full", "lcccccccc",
        "Method & Avg.\\ rank & Best & Feas.\\ (\\%) & $<$15 & $p_z$ & $p_W$ & $r_{\\rm rb}$ & $\\overline{\\Delta L}$ (pp)",
        lines))
    # main-text version (single column; manuscript layout)
    lines = []
    for a in order:
        if a == FOCUS:
            pz = pw = dl = "--"
        else:
            pz = fmt_p(FR["p_holm_vs_focus"][a], 2); pw = fmt_p(CW[a]["p_holm"], 2)
            dl = f"${CW[a]['mean_dloss_pp']:+.3f}$"
        lines.append(f"{LAB[a]} & {FR['avg_rank'][a]:.2f} & {FR['sole_best_count'][a]} & {feas_pct[a]:.1f} & {pz} & {pw} & {dl} \\\\")
    qual_set = (lambda l_: (", ".join(l_) + ", otherwise %d" % FR["n_cases"]) if l_ else "all %d cases" % FR["n_cases"])(
        [f"{CW[b]['n_both_qualified']} for {LAB[b]}" for b in others if CW[b]["n_both_qualified"] < FR["n_cases"]])
    tabs["friedman68"] = table(
        "table", "Case-Level Analysis of the %s Methods over the %d Benchmark Cases" % (
            {7: "Seven", 8: "Eight", 9: "Nine"}.get(len(MAINP), str(len(MAINP))), FR["n_cases"]),
        "tab:friedman68", "lcccccc",
        "Method & Avg.\\ rank & Best & Feas. & $p_z$ & $p_W$ & $\\overline{\\Delta L}$", lines, sep="2pt",
        note="Avg.\\ rank: average rank of the mean feasible objective (1 = best; fewer than 15 feasible runs in a case = "
             "ranked last); Best: cases in which the method alone ranks first; Feas.: feasible runs (\\%%); $p_z$: Holm-adjusted "
             "$p$ of the average-rank test against %s. Because mean-rank tests depend on the pool~\\cite{Benavoli2016}, "
             "$p_W$ gives the two-sided Wilcoxon signed-rank test on the %d per-case mean wake losses (Holm-adjusted over the "
             "%d methods) and $\\overline{\\Delta L}$ the mean wake-loss difference (pp, %s minus method; negative = %s better) "
             "over the cases in which both methods have at least 15 feasible runs (%s). Friedman $\\chi^2_F=%.1f$ (%d d.f.), "
             "$p=%s$; Iman--Davenport $F_F=%.1f$ (%d, %d d.f.), $p=%s$."
             % (fl, FR["n_cases"], len(others), fl, fl, qual_set, FR["chi2"], len(MAINP) - 1, fmt_p(FR["p"], 2).strip("$"),
                FR["iman_davenport"], len(MAINP) - 1, (len(MAINP) - 1) * (FR["n_cases"] - 1),
                fmt_p(FR["iman_davenport_p"], 2).strip("$")))

    # --- further numbers quoted in the text (all written to summary["main"])
    ref = PHASE1.get(FOCUS)                      # phase-1 method of the focus hybrid (PSO for PSO-VNS)
    largest = {}
    for ds in ("1", "2"):
        for r, n in RADII.items():
            s_ = S[(S.Dataset == ds) & (S.Radius == r) & (S.Turbines == n)].set_index("Algorithm")
            largest[f"{ds}-{r}-{n}"] = {a: (float(s_.Loss[a]) if bool(s_.Qualified[a]) else None) for a in MAINP}
    summary["main"]["loss_largest_n"] = largest
    if ref in MAINP:
        dd = [v[ref] - v[FOCUS] for v in largest.values() if v.get(ref) is not None and v.get(FOCUS) is not None]
        summary["main"]["largest_n_loss_reduction_vs_phase1_pp"] = dict(phase1=ref, min=float(min(dd)), max=float(max(dd)),
                                                                        n=len(dd), all_positive=bool(min(dd) > 0))
    sm = S[(S.Turbines <= 3) & S.Qualified]
    summary["main"]["max_loss_n_le_3_pct"] = float(sm.Loss.max())
    summary["main"]["n_runs"] = int(len(G))
    summary["main"]["feasible_runs_by_dataset"] = {a: {ds: int(G[(G.Algorithm == a) & (G.Dataset == ds)].Feasible.sum()) for ds in ("1", "2")}
                                                   for a in MAINP}
    # feasibility of the phase-1 swarm at the switch in the densest cases (500 m, N = 10)
    if FOCUS in HYBRIDS:
        dense = {}
        for ds in ("1", "2"):
            x = G[(G.Algorithm == FOCUS) & (G.Dataset == ds) & (G.Radius == 500) & (G.Turbines == 10)]
            if not len(x):
                continue
            M = curves(x)
            idx = switch_index(switch_call(FOCUS, 6030), 6030)
            ph1 = G[(G.Algorithm == ref) & (G.Dataset == ds) & (G.Radius == 500) & (G.Turbines == 10)] if ref else x.iloc[:0]
            dense[ds] = dict(feasible_at_switch_pct=float(100 * np.isfinite(M[:, idx]).mean()),
                             feasible_final_pct=float(100 * x.Feasible.mean()),
                             phase1_alone_final_feasible_pct=float(100 * ph1.Feasible.mean()) if len(ph1) else None)
        summary["main"]["dense_500_10"] = dense
    # minimum spacing of the feasible final layouts (5D = 385 m, 6D = 462 m)
    Fz = G[G.Feasible]
    mins = np.array([np.min(np.linalg.norm(xy[:, None] - xy[None], axis=2)[np.triu_indices(len(xy), 1)]) if len(xy) > 1 else np.inf
                     for xy in map(coords, Fz.Coordinates)])
    D_ = 77.0
    sp = {}
    for tag, m in (("all", np.ones(len(Fz), bool)), ("n11_15", (Fz.Turbines >= 11).values), ("n6_10", ((Fz.Turbines >= 6) & (Fz.Turbines <= 10)).values)):
        sp[tag] = dict(n_layouts=int(m.sum()), pct_5d=float(100 * np.mean(mins[m] >= 5 * D_ - 1e-6)),
                       pct_6d=float(100 * np.mean(mins[m] >= 6 * D_ - 1e-6)))
    summary["main"]["spacing_share"] = sp
    log(f"  spacing: {sp}")

    # --- earlier PSO setting (w = 0.7, c1 = c2 = 2; fresh_grid.csv) in place of the constriction PSO
    PO = new.get("pso_old")
    if PO is not None and len(PO) and "PSOC" in MAINP:
        PO = PO[PO.Dataset.isin(["1", "2"]) & (PO.Budget == 6030)]
        pool = [("PSO" if a == "PSOC" else a) for a in MAINP]
        GP = pd.concat([G[G.Algorithm != "PSOC"], PO], ignore_index=True)
        SP = case_stats(GP, pool)
        FP = friedman_block(rank_matrix(SP, pool), pool, focus=FOCUS if FOCUS != "PSOC" else "PSO")
        ar = FP["avg_rank"]
        GC = pd.concat([G[G.Algorithm == "PSOC"], PO], ignore_index=True)
        oc = []
        for _, sub in GC.groupby(CASE):
            t = paired_vs(sub, "PSOC", ["PSO"])
            if t:
                oc.append(t[0]["Outcome"])
        SC = case_stats(GC, ["PSOC", "PSO"])
        Lq = SC[SC.Qualified].pivot_table(index=CASE, columns="Algorithm", values="Loss").dropna()
        summary["pso_setting"] = dict(
            old_avg_rank=float(ar["PSO"]), old_rank_position=int(1 + sorted(ar.values()).index(ar["PSO"])),
            n_methods=len(pool), constriction_avg_rank=float(FR["avg_rank"]["PSOC"]),
            constriction_rank_position=int(1 + sorted(FR["avg_rank"].values()).index(FR["avg_rank"]["PSOC"])),
            old_feasible_pct=float(100 * PO.Feasible.mean()), constriction_feasible_pct=float(feas_pct["PSOC"]),
            constriction_vs_old_wtl=pd.Series(oc).value_counts().to_dict(),
            old_minus_constriction_loss_pp=float((Lq["PSO"] - Lq["PSOC"]).mean()),
            n_cases_loss=int(len(Lq)), old_cases_below_half_feasible=int((~SP[SP.Algorithm == "PSO"].Qualified).sum()),
            pool_avg_rank=ar)
        log(f"  old PSO setting: avg rank {ar['PSO']:.2f} (position {summary['pso_setting']['old_rank_position']} of {len(pool)}), "
            f"feasible {summary['pso_setting']['old_feasible_pct']:.1f}%, constriction vs old {summary['pso_setting']['constriction_vs_old_wtl']}")
    else:
        summary["pso_setting"] = None
    summary["baseline"] = baseline_setting(R6, PO, tabs)
    if summary["baseline"]:
        bo, bc = summary["baseline"]["old"]["case_ranks"], summary["baseline"]["constriction"]["case_ranks"]
        bl_ = [f"{'I' if k.split('-')[0] == '1' else 'II'} & {k.split('-')[1]} & {k.split('-')[2]} & "
               + " & ".join(f"{t[k][a]:g}" for t in (bo, bc) for a in ("PSO", "SSABV", "LXBV")) + " \\\\"
               for k in sorted(bo, key=lambda k: (k.split("-")[0], int(k.split("-")[1]), int(k.split("-")[2])))]
        supp.append(table("table", "Per-case ranks (1 = best of 8; feasibility-aware rule) of PSO, SSA-VNS and LX-SSA-VNS in the method pool of Table~\\ref{tab:baseline} with the old PSO setting ($w=0.7$, $c_1=c_2=2$) and with the constriction setting (68 cases, 6,030 evaluations).",
                          "tab:baseline-cases", "ccccccccc",
                          "& & & \\multicolumn{3}{c}{Old setting} & \\multicolumn{3}{c}{Constriction} \\\\\n\\cmidrule(lr){4-6}\\cmidrule(lr){7-9}\nDS & $r$ & $N$ & PSO & SSA-VNS & LX-SSA-VNS & PSO & SSA-VNS & LX-SSA-VNS",
                          bl_, size="\\tiny", sep="2.5pt", pos="p"))

    # --- figures: average ranks
    fig, ax = plt.subplots(figsize=(3.4, 2.1))
    for pos, a in enumerate(order[::-1]):
        v = FR["avg_rank"][a]
        ax.plot([1, v], [pos, pos], color=GRID, lw=2, zorder=1)
        ax.scatter(v, pos, s=90 if a == FOCUS else 36, color=COL[a], marker=mk(a), zorder=zo(a), edgecolor="white", lw=0.6)
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
                    ax.plot(m.Turbines, y, color=COL[a], ls=ls(a), lw=lw(a), marker=mk(a), ms=ms(a), zorder=zo(a))
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
                if i == 1: ax.set_xlabel("Evaluations")
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
            for row, mkr, colr in ((top, "s", INK), (lx, "o", COL[FOCUS])):
                if row is None:
                    continue
                xy = coords(row.Coordinates)
                ax.scatter(xy[:, 0], xy[:, 1], marker=mkr, s=22 if mkr == "s" else 12,
                           facecolor="none" if mkr == "s" else colr, edgecolor=colr, lw=1, zorder=3)
                best_rows.append(dict(Dataset=ds, Radius=r, Turbines=n, Algorithm=row.Algorithm, Seed=row.Seed,
                                      Objective=row.Objective, WakeLoss=row.WakeLoss, LossPct=row.LossPct,
                                      Coordinates=row.Coordinates))
            t = np.linspace(0, 2 * np.pi, 200)
            ax.plot(r * np.cos(t), r * np.sin(t), color=MUTED, lw=0.8)
            ax.set_aspect("equal"); ax.set_xlim(-1.08 * r, 1.08 * r); ax.set_ylim(-1.08 * r, 1.08 * r)
            gap = "" if lx is None else f" (best {fl}: +{100 * (top.Objective - lx.Objective) / top.Ideal:.2f} pp loss)"
            ax.set_title(f"{DSN[ds]}, $r$ = {r} m, $N$ = {n}\nbest of all methods: {LAB[top.Algorithm]}{gap if top.Algorithm != FOCUS else ''}",
                         fontsize=8, color=INK)
            ax.tick_params(labelsize=6)
    hl = [plt.Line2D([], [], ls="none", marker="s", ms=6, mfc="none", mec=INK, mew=1, label="best layout of all methods"),
          plt.Line2D([], [], ls="none", marker="o", ms=4.5, color=COL[FOCUS], label=f"best {fl} layout")]
    fig.legend(handles=hl, loc="lower center", ncol=2, bbox_to_anchor=(0.5, 1.0), fontsize=7.5)
    fig.tight_layout()
    FG.save(fig, "layouts_max")
    BL = pd.DataFrame(best_rows).drop_duplicates()
    BL.to_csv(OUT("mpce_best_layouts_maxN.csv"), index=False)
    cl = []
    for _, b_ in BL.iterrows():
        xy = coords(b_.Coordinates)
        cl.append(f"{'I' if b_.Dataset == '1' else 'II'} & {b_.Radius} & {b_.Turbines} & {LAB[b_.Algorithm]} & {b_.Objective:.1f} & "
                  "\\parbox[t]{9.2cm}{\\raggedright " + "; ".join(f"({x:.1f}, {y:.1f})" for x, y in xy) + "} \\\\")
    supp.append(table("table*", "Coordinates (m, farm center at the origin) of the best layout of all methods and of the best %s layout for the largest $N$ of each farm (6,030 calls); objective in benchmark units." % fl,
                      "tab:coords", "cccccl", "DS & $r$ (m) & $N$ & Method & Objective & Coordinates $(x_i, y_i)$", cl, size="\\tiny", pos="p"))
    summary["main"]["best_layout_maxN"] = [
        dict(case=f"{d}-{r}-{n}", best_method=g.sort_values("Objective").Algorithm.iloc[-1],
             focus_gap_pp=float((g.Objective.max() - g[g.Algorithm == FOCUS].Objective.max()) / g.Ideal.iloc[0] * 100)
             if (g.Algorithm == FOCUS).any() else None)
        for (d, r, n), g in BL.assign(Ideal=BL.Objective + BL.WakeLoss).groupby(CASE)]

    # computational cost
    rt = (G.Seconds / G.Calls * 1000).groupby(G.Algorithm).mean().reindex(MAINP)
    rn = G.groupby(["Algorithm", "Turbines"]).Seconds.mean().unstack().reindex(MAINP)
    src = G.groupby("Algorithm").Source.agg(lambda s: ",".join(sorted(set(s))))
    # time per evaluation (D6): ms per evaluation = Seconds / Calls * 1000 of every run, median over the runs of
    # each method. Data: the 6,030-evaluation runs (random initialization) of the 68 benchmark cases of the main
    # comparison (G: 30 runs x 68 cases per method; Seconds = wall-clock time of the run incl. optimizer
    # overhead, measured by the experiment scripts on the machine that produced each file).
    me = (G.Seconds / G.Calls * 1000).groupby(G.Algorithm).median().reindex(MAINP)
    meta = [a for a in MAINP if a != "SLSQP"]
    summary["cost_per_eval"] = dict(
        median_ms={a: float(me[a]) for a in MAINP}, metaheuristics=meta,
        meta_min_ms=float(me[meta].min()), meta_max_ms=float(me[meta].max()),
        meta_min_method=str(me[meta].idxmin()), meta_max_method=str(me[meta].idxmax()),
        slsqp_ms=float(me["SLSQP"]) if "SLSQP" in me else None,
        slsqp_over_pso=float(me["SLSQP"] / me["PSOC"]) if {"SLSQP", "PSOC"} <= set(MAINP) else None,
        data="6,030-evaluation runs of the 68 cases (random initialization), ms = Seconds / Calls * 1000, median over runs")
    summary["cost"] = dict(ms_per_call=rt.round(4).to_dict(),
                           sec_per_run_range={a: [float(rn.loc[a].min()), float(rn.loc[a].max())] for a in MAINP},
                           source=src.to_dict(),
                           note="wall-clock times come from the machine that produced each file; runs from different "
                                "files may not be directly comparable")
    lines = [f"{LAB[a]} & {rt[a]:.3f} & {rn.loc[a].min():.1f}--{rn.loc[a].max():.1f} \\\\" for a in MAINP]
    supp.append(table("table", "Computational cost at 6,030 objective calls: mean wall-clock time per objective call (including the optimizer overhead) and range over $N$ of the mean time per run (single core). Times come from the machine that produced each result file.",
                      "tab:cost", "lcc", "Method & ms per call & s per run", lines))

    # =========================================================== 2. ablation
    log("\n[2] Ablation")
    GA = R6[R6.Dataset.isin(["1", "2"])]
    have = set(GA.Algorithm)
    # the ablation is built around the focus if it is a two-phase hybrid, otherwise around the hybrid
    # whose phase 1 is the focus (e.g. PSO -> PSO-VNS)
    H0 = FOCUS if FOCUS in HYBRIDS else next(h for h, p1 in PHASE1.items() if p1 == FOCUS)
    if H0 != FOCUS:
        log(f"  focus {fl} is not a two-phase hybrid: ablation built around {LAB[H0]} (phase 1 = {fl})")
    P1 = PHASE1[H0]
    ablv = list(dict.fromkeys([H0] + [h for h in HYBRIDS if h != H0] + ["RSVNS", "RSDVNS", "BVNS", P1] +
                              [PHASE1[h] for h in HYBRIDS if h != H0]))           # table order: hybrids, controls, swarms
    ablp = [a for a in ablv if a in have]
    if [a for a in ablv if a not in have]:
        log(f"  ablation variants without data (dropped with their contrasts): {[a for a in ablv if a not in have]}")
    GB = GA[GA.Algorithm.isin(ablp)]
    CONTR = [(H0, P1, f"VNS phase ({LAB[P1]})")]
    CONTR += [(h, PHASE1[h], f"VNS phase ({LAB[PHASE1[h]]})") for h in HYBRIDS if h != H0]
    CONTR += [(H0, "BVNS", "swarm start vs.\\ best initial point"),
              (H0, "RSVNS", f"{LAB[P1]} vs.\\ random sampling")]
    CONTR += [(h, "RSVNS", f"{LAB[PHASE1[h]]} vs.\\ random sampling") for h in HYBRIDS if h != H0]
    # Phase 6: the disc-sampling control RSD-VNS (experiment rsdisc), a stronger random-sampling Phase 1
    CONTR += [(H0, "RSDVNS", f"{LAB[P1]} vs.\\ disc sampling")]
    CONTR += [(h, "RSDVNS", f"{LAB[PHASE1[h]]} vs.\\ disc sampling") for h in HYBRIDS if h != H0]
    CONTR += [("RSDVNS", "RSVNS", "disc vs.\\ square sampling")]
    CONTR += [(H0, h, f"{LAB[P1]} vs.\\ {LAB[PHASE1[h]]} as Phase~1") for h in HYBRIDS if h != H0]
    CONTR += [("SSABV", "LXBV", "LX-SSA as published (hybrid)")] if H0 != "LXBV" else []   # R1-2 / D14: not "Laplace step"
    CONTR += [("LXSSA", "SSA", "LX-SSA as published (alone)")]
    CONTR = list(dict.fromkeys(CONTR))
    CONTR = [c for c in CONTR if c[0] in ablp and c[1] in ablp]
    arows = []
    for (ds, r, n), sub in GB.groupby(CASE):
        piv = sub.assign(S=goodness(sub)).pivot_table(index="Seed", columns="Algorithm", values="S")
        lm = sub[sub.Feasible].groupby("Algorithm").LossPct.mean()
        rr = []
        for a, b, _ in CONTR:
            if a not in piv or b not in piv:
                continue
            d = (piv[a] - piv[b]).dropna()
            p, rb, _ = wil(d.values)
            rr.append(dict(Dataset=ds, Radius=r, Turbines=n, A=a, B=b, P=p, RB=rb, NPairs=len(d),
                           DLoss=lm.get(a, np.nan) - lm.get(b, np.nan)))
        for x, h in zip(rr, holm([x["P"] for x in rr])):
            x["PHolm"] = h
        arows += rr
    A = pd.DataFrame(arows)
    A["Outcome"] = np.where(A.PHolm < 0.05, np.where(A.RB > 0, "W", "L"), "T")
    A.to_csv(OUT("mpce_ablation_tests.csv"), index=False)
    SA = case_stats(GB, ablp)
    FA = friedman_block(rank_matrix(SA, ablp), ablp, focus=H0)
    # case-mean Wilcoxon of every contrast (same procedure as the main case-mean test, case_mean_wilcoxon):
    # p unadjusted, p_holm = Holm over the contrasts of the table (one family), wins / losses = cases in which
    # the first variant has the lower / higher case-mean loss (imputed cases included), mean_dloss_pp = first
    # minus second (negative = first better) with a 95 % bootstrap CI over the both-qualified cases
    CM = {}
    for a, b, _ in CONTR:
        x = case_mean_wilcoxon(SA, a, [b])[b]
        CM[f"{a}-{b}"] = {k: x[k] for k in ("p", "rb", "n_cases", "wins", "losses", "mean_dloss_pp", "median_dloss_pp",
                                              "ci95_mean_dloss_pp", "only_focus_qualified", "only_other_qualified",
                                              "neither_qualified", "n_both_qualified")}
    for k, h in zip(CM, holm([v["p"] for v in CM.values()])):
        CM[k]["p_holm"] = float(h)
    comp = {"PSOBV": ("PSO", "VNS"), "SSABV": ("SSA", "VNS"), "LXBV": ("LX-SSA", "VNS"), "RSVNS": ("random (square)", "VNS"),
            "RSDVNS": ("random (disc)", "VNS"),
            "BVNS": ("--", "VNS"), "PSOC": ("PSO", "--"), "SSA": ("SSA", "--"), "LXSSA": ("LX-SSA", "--")}
    afeas = {a: float(GB[GB.Algorithm == a].Feasible.mean() * 100) for a in ablp}
    rows1 = [f"{LAB[a]} & {comp[a][0]} & {comp[a][1]} & {FA['avg_rank'][a]:.2f} & {afeas[a]:.1f} \\\\" for a in ablp]
    rows2 = []
    for a, b, what in CONTR:
        x = A[(A.A == a) & (A.B == b)]
        rows2.append(f"{LAB[a]} vs.\\ {LAB[b]} & {what} & {wtl_str(x.Outcome.value_counts())} & "
                     f"{fmt_p(CM[f'{a}-{b}']['p_holm'], 2)} & ${CM[f'{a}-{b}']['mean_dloss_pp']:+.3f}$ \\\\")
    nw = {6: "Six", 7: "Seven", 8: "Eight", 9: "Nine", 10: "Ten", 11: "Eleven", 12: "Twelve", 13: "Thirteen", 14: "Fourteen",
          15: "Fifteen", 16: "Sixteen"}
    ncw = nw.get(len(CONTR), len(CONTR))
    nq = [CM[f"{a}-{b}"]["n_both_qualified"] for a, b, _ in CONTR]
    nqs = (f"{min(nq)}" if min(nq) == max(nq) else f"{min(nq)}--{max(nq)}") + " of %d per contrast" % FA["n_cases"]
    anote = ("Top: average rank among the %s variants (Friedman $\\chi^2_F=%.1f$, %d d.f., $p=%s$) and feasible runs. "
             "Bottom: planned contrasts. W/T/L: cases in which the first variant is significantly better / not different / "
             "worse (run-level Wilcoxon, $\\alpha=0.05$, Holm-adjusted over the %s contrasts of each case, a different family "
             "from Table~\\ref{tab:wtl}); $p_W$: Wilcoxon test on the %d per-case mean wake losses, Holm-adjusted over the %s "
             "contrasts; $\\overline{\\Delta L}$: mean difference of the case-mean wake losses (pp; negative = first variant "
             "better) over the cases in which both variants have at least 15 feasible runs (%s)."
             % (nw.get(len(ablp), len(ablp)).lower(), FA["chi2"], len(ablp) - 1, fmt_p(FA["p"]).strip("$"), ncw.lower(),
                FA["n_cases"], ncw.lower(), nqs))
    tabs["ablation"] = (r"""\begin{table}[!t]
\centering
\caption{Component Analysis over the %d Benchmark Cases (6,030 Calls, 30 Seed-Paired Runs)}
\label{tab:ablation}
\scriptsize\setlength{\tabcolsep}{2.5pt}
\begin{tabular}{lcccc}
\toprule
Variant & Phase 1 & Phase 2 & Avg.\ rank & Feas.\ (\%%) \\
\midrule
""" % FA["n_cases"] + "\n".join(rows1) + r"""
\bottomrule
\end{tabular}

\smallskip
\begin{tabular}{l>{\raggedright\arraybackslash}p{1.95cm}ccc}
\toprule
Contrast & Isolates & W/T/L & $p_W$ & $\overline{\Delta L}$ \\
\midrule
""" + "\n".join(rows2) + r"""
\bottomrule
\end{tabular}
""" + tnote(anote) + r"""\end{table}
""")
    # per-case component analysis (supplement): mean loss of the two variants that have no per-case table
    # elsewhere (LX-SSA-VNS, RS-VNS) and the run-level outcome of every contrast (+ / 0 / -)
    sym = {"W": "$+$", "T": "$\\cdot$", "L": "$-$"}
    Lq = SA.pivot_table(index=CASE, columns="Algorithm", values="Loss")
    Qq = SA.pivot_table(index=CASE, columns="Algorithm", values="Qualified").astype(float)
    extra_v = [a for a in ("LXBV", "RSVNS", "RSDVNS") if a in ablp]
    pl = []
    for (ds, r, n), x in A.groupby(CASE):
        oc = {(y.A, y.B): y.Outcome for y in x.itertuples()}
        lv = [(f"{Lq.loc[(ds, r, n), a]:.3f}" if Qq.loc[(ds, r, n), a] > 0.5 else "--") for a in extra_v]
        pl.append(f"{'I' if ds == '1' else 'II'} & {r} & {n} & " + " & ".join(lv) + " & "
                  + " & ".join(sym.get(oc.get((a, b)), "") for a, b, _ in CONTR) + " \\\\")
    chead = " & ".join(f"C{i + 1}" for i in range(len(CONTR)))
    key_ = "; ".join(f"C{i + 1}: {LAB[a]} vs.\\ {LAB[b]}" for i, (a, b, _) in enumerate(CONTR))
    supp.append(table("table*", "Component analysis per case (68 cases, 6,030 calls, 30 seed-paired runs): mean wake loss (\\%%) of the feasible runs of %s (``--'': fewer than 15 feasible runs), and run-level outcome of each planned contrast ($+$: first variant significantly better, $-$: significantly worse, $\\cdot$: no significant difference; Wilcoxon signed-rank, Holm-adjusted over the %d contrasts of the case). %s." % ((lambda l_: ", ".join(l_[:-1]) + " and " + l_[-1] if len(l_) > 1 else "".join(l_))([LAB[a] for a in extra_v]), len(CONTR), key_),
                      "tab:ablation-cases", "ccc" + "c" * (len(extra_v) + len(CONTR)),
                      "DS & $r$ & $N$ & " + " & ".join(LAB[a] for a in extra_v) + " & " + chead, pl, size="\\tiny", sep="2.2pt", pos="p"))
    summary["ablation"] = dict(hybrid=H0, phase1=P1, variants=ablp, friedman=FA, feasible_pct=afeas, n_contrasts=len(CONTR),
                               holm_family=f"per case, over the {len(CONTR)} contrasts (run level); case-mean p_holm over the same {len(CONTR)} contrasts",
                               contrasts={f"{a}-{b}": dict(A[(A.A == a) & (A.B == b)].Outcome.value_counts().to_dict(),
                                                            dloss_pp=float(A[(A.A == a) & (A.B == b)].DLoss.mean()),
                                                            isolates=w.replace("\\", ""))
                                          for a, b, w in CONTR},
                               case_mean=CM)
    log(f"  ablation hybrid {LAB[H0]}; avg ranks: " + ", ".join(f"{LAB[a]} {FA['avg_rank'][a]:.3f}" for a in ablp))
    for k, v in summary["ablation"]["contrasts"].items():
        log(f"  {k}: {wtl_str(v)}  dL={v['dloss_pp']:+.3f}")

    # phase-2 statistic (VNS phase after the swarm) for every hybrid (and RS-VNS): share of the wake loss
    # left at the switch that the VNS phase removes; runs infeasible at the switch that end feasible
    ph2 = {alg: phase2_stat(GA, alg) for alg in HYBRIDS + ["RSVNS", "RSDVNS"] if alg in have}
    if P1 in have:      # the phase-1 swarm alone over the same second half of the budget (no VNS)
        ph2[P1 + "_continued"] = phase2_stat(GA, P1, switch=switch_call(H0, 6030))
    summary["ablation"]["n_runs"] = int(len(GB))
    summary["ablation"]["n_runs_benchmark_total"] = int(len(GA[GA.Algorithm.isin(set(MAINP) | set(ablp))]))
    summary["ablation"]["phase2_loss_reduction_pct"] = ph2
    # runs whose Phase-1 random samples contain no feasible layout (RS-VNS: square, RSD-VNS: disc), geometry-only replay
    summary["ablation"]["phase1_replay"] = phase1_replay(GA)
    swl = []
    for alg, v in ph2.items():
        if v is None:
            continue
        nm_ = f"{LAB[P1]} continued (no VNS)" if alg.endswith("_continued") else LAB[alg]
        f1 = lambda x, d=1: "--" if x is None or not np.isfinite(x) else f"{x:.{d}f}"
        swl.append(f"{nm_} & {v['switch_call']:,} & {f1(v['feasible_at_switch_pct'])} & {f1(v['mean_loss_at_switch_pct'], 2)} & "
                   f"{f1(v['mean'])} & {f1(v['median'])} & {v['infeasible_at_switch_made_feasible']} \\\\".replace(",", "{,}", 1))
    supp.append(table("table", "Switch point of the two-phase variants (68 cases, 30 runs, 6,030 evaluations): number of Phase-1 evaluations (the VNS phase starts with the next evaluation), runs whose best layout is feasible at the switch (\\%%; read at the last convergence checkpoint at or before the switch, i.e.\\ at 3{,}000 evaluations for RS-VNS and RSD-VNS, whose last 15 Phase-1 samples are not observable in the stored curves), mean wake loss at the switch (\\%%, feasible runs), mean and median share (\\%%) of the wake loss left at the switch that the second phase removes, and runs infeasible at the switch that end feasible. ``%s continued'': %s alone over the same evaluations." % (LAB[P1], LAB[P1]),
                      "tab:switch", "lcccccc", "Variant & Switch & Feas. (\\%) & Loss (\\%) & Mean share & Median share & Made feasible", swl, pos="!htb"))
    for alg, v in ph2.items():
        if v is None:
            continue
        log(f"  phase 2 of {LAB.get(alg, alg)} (switch at call {v['switch_call']:,}) removes {v['mean']:.2f}% (median {v['median']:.2f}%) "
            f"of the loss left at the switch; {v['infeasible_at_switch_made_feasible']} runs infeasible at the switch made feasible")

    # practical equivalence of the case means (EQ_MARGIN) and spread / worst run of PSO-VNS vs PSO
    log("\n[2b] Practical equivalence (margin %.2f pp, set after the primary analysis) and spread" % EQ_MARGIN)
    equivalence_block(S, SA, G, summary, supp)

    # reproduction check against the parent of commit a9779a3 (old rule, five variants, six contrasts);
    # independent of --focus (SSA-VNS numbers of the earlier pipeline)
    OLDV = ["SSABV", "LXBV", "BVNS", "LXSSA", "SSA"]
    if all(a in have for a in OLDV):
        Go = GA[GA.Algorithm.isin(OLDV)]
        So = case_stats(Go, OLDV, rank_old)
        Fo = friedman_block(rank_matrix(So, OLDV), OLDV, focus="SSABV")
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
        p2s = phase2_stat(GA, "SSABV", switch=3030, validate=False)
        rep = dict(ssabv_avg_rank=Fo["avg_rank"]["SSABV"], expected_ssabv_avg_rank=1.625,
                   lxbv_vs_ssabv_wtl=wtl_str(oc), expected_lxbv_vs_ssabv_wtl="0/66/2",
                   lxbv_minus_ssabv_dloss=float(x.DL.mean()), expected_dloss=0.06306204787105701,
                   phase2_ssabv=p2s, expected_phase2_ssabv=dict(mean=36.66372550902314, median=29.194081556313783,
                                                                infeasible_at_switch=79),
                   old_avg_rank=Fo["avg_rank"])
        rep["pass"] = bool(abs(rep["ssabv_avg_rank"] - 1.625) < 1e-9 and rep["lxbv_vs_ssabv_wtl"] == "0/66/2"
                           and abs(rep["lxbv_minus_ssabv_dloss"] - 0.06306204787105701) < 1e-9
                           and abs(p2s["mean"] - 36.66372550902314) < 1e-6
                           and p2s["infeasible_at_switch_made_feasible"] == 79)
        summary["reproduction_check_a9779a3_parent"] = rep
        log(f"  REPRODUCTION CHECK (SSA-VNS, old rule, 5 variants, Holm over 6 contrasts): SSABV avg rank {rep['ssabv_avg_rank']:.4f} "
            f"(expected 1.625), LXBV vs SSABV {rep['lxbv_vs_ssabv_wtl']} (expected 0/66/2), dL {rep['lxbv_minus_ssabv_dloss']:+.4f} "
            f"(expected +0.0631), phase-2 mean {p2s['mean']:.3f} (expected 36.664) -> {'PASS' if rep['pass'] else 'FAIL'}"
            + ("" if not dropped else " [--common-seeds dropped runs: check not applicable]"))

    # ablation convergence (largest N)
    fig, axes = plt.subplots(2, 3, figsize=(7.1, 4.2))
    for i, ds in enumerate(("1", "2")):
        for j, (r, n) in enumerate(RADII.items()):
            ax = axes[i, j]
            conv_panel(ax, GB[(GB.Dataset == ds) & (GB.Radius == r) & (GB.Turbines == n)], ablp, band=False)
            ax.axvline(switch_call(H0, 6030), color=MUTED, lw=0.7, ls=(0, (1, 2)))
            log_axis(ax)
            ax.set_title(f"{DSN[ds]}, $r$ = {r} m, $N$ = {n}", fontsize=8, color=INK)
            if i == 1: ax.set_xlabel("Evaluations")
            if j == 0: ax.set_ylabel("Median best wake loss (%)")
    legend_row(fig, ablp)
    fig.tight_layout()
    FG.save(fig, "ablation_convergence")

    # =========================================================== 3. budget split
    log("\n[3] Budget split")
    summary["split"] = split_section(R6, FOCUS, tabs, "split", "tab:split", primary=True, supp=supp)
    if summary["split"]:
        summary["holm_families"]["split"] = "per case: the %d run-level comparisons of the split table (%s)" % (
            summary["split"]["n_comparisons"], summary["split"]["holm_family"].split("comparisons ", 1)[1])
    if FOCUS != "SSABV":
        summary["split_ssabv"] = split_section(R6, "SSABV", tabs, "split_ssabv", "tab:split-ssabv", primary=False, supp=supp)

    # =========================================================== 4. Horns Rev 16
    log("\n[4] Horns Rev 1, 16 turbines")
    try:
        import hornsrev_model as hr
        inst = float(hr.aep_gwh(hr.site(16)[0]))
        ideal_model = float(hr.aep_gwh(np.zeros((16, 2)), with_wake=False))
        xy80 = hr.site(80)[0]
        a80, i80 = float(hr.aep_gwh(xy80)), float(hr.aep_gwh(xy80, with_wake=False))
        pwr = pywake_reference(args.data_dir)
        summary["hr_validation"] = dict(
            installed80_aep=a80, ideal80_aep=i80, installed80_loss_pct=100 * (1 - a80 / i80),
            pywake_aep=pwr["aep"] if pwr else None, pywake_loss_pct=pwr["loss_pct"] if pwr else None,
            pywake_ideal_aep=pwr["ideal"] if pwr else None, pywake_version=pwr["version"] if pwr else None,
            rel_diff_pct=100 * (a80 / pwr["aep"] - 1) if pwr else None,
            loss_diff_pp=(100 * (1 - a80 / i80) - pwr["loss_pct"]) if pwr else None,
            ours_matches_pywake_check=(abs(pwr["ours_csv_aep"] - a80) < 1e-6) if pwr and pwr["ours_csv_aep"] is not None else None,
            note=("installed 80-turbine Horns Rev 1 farm, hornsrev_model.py (current binning) vs PyWake %s NOJ(k=0.04) at "
                  "identical bins (%s, pywake_check.csv)" % (pwr["version"], pwr["bins"])) if pwr else
                 "pywake_check.csv missing: no PyWake reference (run pywake_check.py with py_wake==2.6.20)")
        if pwr:
            log(f"  80-turbine installed AEP {a80:.2f} GWh/yr vs PyWake {pwr['version']} NOJ {pwr['aep']:.2f} (identical bins) -> "
                f"{summary['hr_validation']['rel_diff_pct']:+.2f}% AEP, {summary['hr_validation']['loss_diff_pp']:+.2f} pp wake loss"
                + ("" if summary["hr_validation"]["ours_matches_pywake_check"] else "  WARNING: our AEP differs from pywake_check.csv"))
        else:
            log("  pywake_check.csv missing: \\NHRPyWakeDiff pending")
    except Exception as e:                                  # pragma: no cover
        hr, inst, ideal_model = None, INSTALLED_HR16, None
        log(f"  hornsrev_model unavailable ({e}); installed AEP = {INSTALLED_HR16}")
    H = R6[(R6.Dataset == "HR") & (R6.Turbines == 16)]
    hr_model = "current"
    if len(H):
        ideal_data = float(H.Ideal.iloc[0])
        if abs(ideal_data - LEGACY_HR16["ideal"]) < 1e-6:      # old runs (binning before 7676da9): same-model installed AEP
            inst, hr_model = LEGACY_HR16["installed"], "legacy (direction binning before commit 7676da9; mpce_hrfix missing)"
            log(f"  Horns Rev runs of the OLD model (wake-free AEP {ideal_data:.4f}); installed AEP of the old model {inst:.3f}")
        elif ideal_model is not None and abs(ideal_data - ideal_model) > 1e-6:
            log(f"  WARNING: Horns Rev wake-free AEP of the runs ({ideal_data:.4f}) differs from hornsrev_model ({ideal_model:.4f})")
    hm = [a for a in HR_ORDER if a in M10]
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
                            runs_above_installed=int(((H.Algorithm == a) & H.Feasible & (H.Objective > inst)).sum()),
                            p_holm=None if a == FOCUS else tt[a]["PHolm"], rb=None if a == FOCUS else tt[a]["RB"])
        supp.append(table("table", "Horns Rev~1 site case, 16-turbine block (Vestas V80 power and thrust curves, measured 12-sector wind climate, installed outline; Jensen wake $k=0.04$; minimum spacing $4D=320$~m; wake-free AEP %.2f~GWh/yr): AEP in GWh/yr, mean (SD) over the feasible runs of 30 seed-paired runs at 6,030 objective calls, mean wake loss (\\%%), feasible runs, Wilcoxon signed-rank $p$ of %s vs.\\ each method (Holm-adjusted over these %d comparisons) and rank-biserial $r_{\\rm rb}$ (positive = %s better; zero differences dropped); run-level Friedman $p=%s$ (over the seeds common to all methods)." % (ideal, fl, len(hm) - 1, fl, fmt_p(pf).strip("$")),
                             "tab:hr-site-full", "lccccc", "Layout / method & AEP & Loss (\\%) & Feas. & $p_{\\rm Holm}$ & $r_{\\rm rb}$", lines, size="\\footnotesize", sep="3pt", pos="!htb"))
        summary["hr16"] = dict(installed_aep=inst, ideal_aep=ideal, installed_loss_pct=100 * (1 - inst / ideal),
                               model=hr_model, friedman_p=float(pf), methods=hrsum,
                               holm_family=f"per setting, {LAB[FOCUS]} vs each of the other {len(hm) - 1} methods",
                               sources=H.groupby("Algorithm").Source.first().to_dict())
        log("  mean AEP: " + ", ".join(f"{LAB[a]} {hrsum[a]['mean']:.2f}({hrsum[a]['feasible']})" for a in hm) + f"; installed {inst:.2f}")
        # convergence + layouts
        fig, axes = plt.subplots(1, 2, figsize=(7.1, 2.8), gridspec_kw=dict(width_ratios=[1.25, 1]))
        ax = axes[0]
        conv_panel(ax, H, hm, loss=False, band=False)
        ax.axhline(inst, color=INK, lw=0.9, ls=(0, (1, 1)))
        ax.text(1000, inst, "installed layout", va="bottom", fontsize=6.5, color=INK)
        ax.set_title("Horns Rev 1, 16 turbines: median best feasible AEP", fontsize=8, color=INK)
        ax.set_xlabel("Evaluations"); ax.set_ylabel("AEP (GWh/yr)")
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
            ax.scatter(xy[:, 0], xy[:, 1], marker="o", s=12, color=COL[FOCUS], label=f"best {fl} ({b.Objective:.2f})", zorder=4)
        xy = coords(top.Coordinates)
        ax.scatter(xy[:, 0], xy[:, 1], marker="s" if top.Algorithm != FOCUS else "o", s=22 if top.Algorithm != FOCUS else 12,
                   facecolor="none" if top.Algorithm != FOCUS else COL[FOCUS], edgecolor=INK if top.Algorithm != FOCUS else COL[FOCUS],
                   lw=1, label=f"best of all methods: {LAB[top.Algorithm]} ({top.Objective:.2f})", zorder=3)
        summary["hr16"]["best_overall"] = dict(method=top.Algorithm, aep=float(top.Objective), seed=int(top.Seed))
        ax.set_aspect("equal"); ax.tick_params(labelsize=6)
        ax.set_title("Layouts (m); AEP in GWh/yr", fontsize=8, color=INK)
        ax.legend(fontsize=6, loc="upper left", bbox_to_anchor=(0, -0.1), ncol=1)
        legend_row(fig, hm, ncol=5)
        fig.tight_layout()
        FG.save(fig, "hr16")
    else:
        log("  SKIPPED: no Horns Rev data")

    # =========================================================== 5. feasible vs random initialization
    log("\n[5] Feasible initialization")
    FE = ALL[(ALL.Init == "feasible") & (ALL.Budget == 6030)]
    cases7 = LARGE + [HR16]
    if len(FE):
        fm_ = [a for a in M10 if a in set(FE.Algorithm)]
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
            lines.append(f"{LAB[a]}{fb} & {x.LossR.mean():.3f} & " + ("--" if not np.isfinite(x.LossF.mean()) else f"{x.LossF.mean():.3f}") + f" & {100 * x.FeasR.sum() / x.NR.sum():.1f} & "
                         f"{100 * x.FeasF.sum() / x.NF.sum():.1f} & {wtl_str(x.Outcome.value_counts())} \\\\")
        supp.append(table("table", "Feasibility-preserving versus uniform random initialization at 6,030 calls on the six largest benchmark cases and the Horns Rev 16-turbine block (%d cases, 30 seed-paired runs each): mean wake loss (\\%%) of the feasible runs averaged over the cases, percentage of feasible runs, and number of cases in which feasible initialization is significantly better / not different / worse (Wilcoxon signed-rank, Holm-adjusted over the methods of each case). Per-case values: Supplementary Table~\\ref{tab:feasinit-cases}.%s" % (Dt.case.nunique(), " $^\\dagger$: random-initialization runs are a development fallback." if "dagger" in "".join(lines) else ""),
                                 "tab:feasinit", "lccccc",
                                 "Method & \\multicolumn{2}{c}{Loss (\\%)} & \\multicolumn{2}{c}{Feas. (\\%)} & W/T/L \\\\\n & random & feasible & random & feasible &", lines, pos="!htb"))
        dl = []
        for case, x in Dt.groupby("case", sort=False):
            dl.append(f"\\multicolumn{{7}}{{l}}{{\\emph{{{case_name(case)}}}}} \\\\")
            for _, y in x.iterrows():
                isHR = case.startswith("HR")
                fmt_ = lambda v, d: "--" if v is None or not np.isfinite(v) else f"{v:.{d}f}"
                vr = fmt_(y.AEPR, 2) if isHR else fmt_(y.LossR, 3); vf = fmt_(y.AEPF, 2) if isHR else fmt_(y.LossF, 3)
                dl.append(f"{LAB[y.Algorithm]} & {vr} & {vf} & {y.FeasR}/{y.NR} & {y.FeasF}/{y.NF} & {fmt_p(y.PHolm)} & {y.RB:+.2f} \\\\")
        supp.append(table("table", "Feasible versus random initialization per case (6,030 calls): mean wake loss (\\%; Horns Rev: mean AEP in GWh/yr) of the feasible runs, feasible runs, Holm-adjusted Wilcoxon signed-rank $p$ (over the methods of the case) and rank-biserial $r_{\\rm rb}$ (positive = feasible initialization better).",
                          "tab:feasinit-cases", "lcccccc", "Method & Random & Feasible & Feas. (R) & Feas. (F) & $p_{\\rm Holm}$ & $r_{\\rm rb}$", dl, pos="p"))
        summary["feasinit"] = dict(cases=sorted(Dt.case.unique()), per_method={
            a: dict(loss_random=float(x.LossR.mean()), loss_feasible=float(x.LossF.mean()),
                    feas_random_pct=float(100 * x.FeasR.sum() / x.NR.sum()), feas_feasible_pct=float(100 * x.FeasF.sum() / x.NF.sum()),
                    wtl=x.Outcome.value_counts().to_dict(), random_source=sorted(set(x.RandomSource)))
            for a, x in Dt.groupby("Algorithm")})
        # runs whose final objective differs from the best initial layout (first convergence checkpoint = best of
        # the 30 initial layouts at 6,030 evaluations; curves are stored with 3-4 decimals -> tolerance 2e-3)
        c0 = pd.to_numeric(FE.Curve.str.split(";").str[0], errors="coerce")
        moved = (FE.Objective - c0).abs() > 2e-3
        summary["feasinit"]["runs_moved_from_best_initial"] = {
            a: dict(moved=int((moved & (FE.Algorithm == a)).sum()), runs=int((FE.Algorithm == a).sum()),
                    moved_grid6=int((moved & (FE.Algorithm == a) & (FE.Dataset != "HR")).sum()),
                    runs_grid6=int(((FE.Algorithm == a) & (FE.Dataset != "HR")).sum()))
            for a in fm_}
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
    BD = ALL[(ALL.Init == "random") & ALL.Algorithm.isin(M10)]
    BD = BD[[(d, r, n) in cases7 for d, r, n in zip(BD.Dataset, BD.Radius, BD.Turbines)]]
    have_b = [b for b in BUDGETS if b != 6030 and (BD.Budget == b).any()]
    if have_b:
        bl = [6030] + have_b
        BD = BD[BD.Budget.isin(bl)]
        bsum, rows_c = {}, []
        for b in bl:
            x = BD[BD.Budget == b]
            mets = [a for a in M10 if a in set(x.Algorithm)]
            Sb = case_stats(x, mets)
            Sb["Budget"] = b
            rows_c.append(Sb)
            Rb = rank_matrix(Sb, mets)
            avg = Rb.mean(axis=0).to_dict()
            grid = Sb[Sb.Dataset != "HR"]
            Rg = rank_matrix(grid, mets)
            avg6 = Rg.mean(axis=0).to_dict()
            bsum[b] = dict(n_cases=len(Rb), avg_rank=avg, best=min(avg, key=avg.get),
                           avg_rank_grid6=avg6, n_cases_grid6=int(len(Rg)), best_grid6=min(avg6, key=avg6.get) if avg6 else None,
                           runs_grid6={a: int(((x.Algorithm == a) & (x.Dataset != "HR")).sum()) for a in mets},
                           runs_hr16={a: int(((x.Algorithm == a) & (x.Dataset == "HR")).sum()) for a in mets},
                           focus_rank_position=int(1 + sorted(avg.values()).index(avg[FOCUS])) if FOCUS in avg else None,
                           mean_loss_grid={a: float(grid[(grid.Algorithm == a) & grid.Qualified].Loss.mean()) for a in mets},
                           feas_pct={a: float(100 * x[x.Algorithm == a].Feasible.mean()) for a in mets},
                           hr16_aep={a: float(Sb[(Sb.Dataset == "HR") & (Sb.Algorithm == a)].Mean.iloc[0])
                                     for a in mets if ((Sb.Dataset == "HR") & (Sb.Algorithm == a)).any()},
                           sources={a: sorted(set(x[x.Algorithm == a].Source)) for a in mets})
            # run-level W/T/L of the focus vs each method over the 7 cases (Holm over the comparisons of each case)
            if FOCUS in mets:
                oc_ = []
                for _, sub in x.groupby(CASE):
                    oc_ += paired_vs(sub, FOCUS, [a for a in mets if a != FOCUS])
                O_ = pd.DataFrame(oc_)
                bsum[b]["wtl_vs_focus"] = {a: {k: int(v) for k, v in O_[O_.Baseline == a].Outcome.value_counts().items()}
                                           for a in mets if a != FOCUS}
            # focus vs MS-SLSQP per case
            if FOCUS in mets and "SLSQP" in mets:
                vs = {}
                for (d, r, n), sub in x.groupby(CASE):
                    t = paired_vs(sub, FOCUS, ["SLSQP"])
                    ss = Sb[(Sb.Dataset == d) & (Sb.Radius == r) & (Sb.Turbines == n)].set_index("Algorithm")
                    vs[f"{d}-{r}-{n}"] = dict(p=t[0]["P"] if t else None, rb=t[0]["RB"] if t else None,
                                              focus=float(ss.Mean.get(FOCUS, np.nan)), slsqp=float(ss.Mean.get("SLSQP", np.nan)),
                                              focus_loss=float(ss.Loss.get(FOCUS, np.nan)), slsqp_loss=float(ss.Loss.get("SLSQP", np.nan)))
                bsum[b]["focus_vs_slsqp"] = vs
            log(f"  budget {b}: avg ranks " + ", ".join(f"{LAB[a]} {v:.2f}" for a, v in sorted(avg.items(), key=lambda t: t[1])))
        SB = pd.concat(rows_c, ignore_index=True)
        SB.to_csv(OUT("mpce_budget_case_stats.csv"), index=False)
        summary["budget"] = dict(budgets=bl, per_budget=bsum,
                                 focus_best_at_all_budgets=all(bsum[b]["best"] == FOCUS for b in bl),
                                 focus_best_grid6_at_all_budgets=all(bsum[b]["best_grid6"] == FOCUS for b in bl),
                                 note="HR16 at 120,030 calls uses 10 seeds (qualification: at least 5 feasible runs)")
        ms_all = [a for a in M10 if a in set(BD.Algorithm)]
        hdr = "Method & " + " & ".join(f"\\multicolumn{{4}}{{c}}{{{b:,} evaluations}}".replace(",", "{,}") for b in bl) + " \\\\\n & " + \
              " & ".join("Rank & Loss & Feas. & W/T/L" for _ in bl)
        lines = []
        for a in ms_all:
            c = []
            for b in bl:
                v = bsum[b]
                c += [f"{v['avg_rank'][a]:.2f}" if a in v["avg_rank"] else "--",
                      f"{v['mean_loss_grid'][a]:.3f}" if np.isfinite(v["mean_loss_grid"].get(a, np.nan)) else "--",
                      f"{v['feas_pct'][a]:.0f}" if a in v["feas_pct"] else "--",
                      "--" if a == FOCUS or a not in v.get("wtl_vs_focus", {}) else wtl_str(v["wtl_vs_focus"][a])]
            lines.append(f"{LAB[a]} & " + " & ".join(c) + " \\\\")
        supp.append(table("table*", "Budget scaling (random initialization) on the six largest benchmark cases and the Horns Rev 16-turbine block: average rank over the %d cases (ranking rule of Table~\\ref{tab:friedman68}), mean wake loss (\\%%) of the feasible runs averaged over the six benchmark cases in which the method has at least half of its runs feasible, percentage of feasible runs, and cases in which %s is significantly better / not different / worse than the method (run-level Wilcoxon signed-rank, Holm-adjusted over the comparisons of %s in each case; exploratory). Per-case values: Table~\\ref{tab:budget-cases}." % (len(cases7), fl, fl),
                               "tab:budget", "l" + "cccc" * len(bl), hdr, lines, sep="2.2pt", pos="!htb"))
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
        supp.append(tb)
        # figure
        fig, axes = plt.subplots(2, 4, figsize=(7.1, 3.9))
        for ax, c in zip(axes.flat, cases7):
            sc = SB[(SB.Dataset == c[0]) & (SB.Radius == c[1]) & (SB.Turbines == c[2])]
            isHR = c[0] == "HR"
            for a in ms_all:
                y = sc[sc.Algorithm == a].sort_values("Budget")
                v = y.Loss.where(y.Qualified)
                ax.plot(y.Budget, v, color=COL[a], ls=ls(a), marker=mk(a), ms=ms(a), lw=lw(a, 1.1), zorder=zo(a))
            ax.set_xscale("log")
            ax.set_xticks(bl); ax.set_xticklabels([f"{b // 1000}k" if b % 1000 == 30 else str(b) for b in bl], fontsize=6)
            ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
            ax.set_title(case_name("-".join(map(str, c))), fontsize=6.5, color=INK)
            ax.set_ylabel("Mean AEP loss (%)" if isHR else "Mean wake loss (%)", fontsize=7)
        axes.flat[-1].axis("off")
        h = [plt.Line2D([], [], color=COL[a], ls=ls(a), lw=lw(a), marker=mk(a), ms=ms(a, 4), label=LAB[a]) for a in ms_all]
        axes.flat[-1].legend(handles=h, loc="center", fontsize=7)
        fig.supxlabel("Evaluations (log scale)", fontsize=8)
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
                mets = [a for a in M10 if a in set(y.Algorithm)]
                Sy = case_stats(y, mets).set_index("Algorithm")
                tt = {t["Baseline"]: t for t in paired_vs(y, FOCUS, [a for a in mets if a != FOCUS])} if FOCUS in mets else {}
                isum[cname]["budgets"][int(b)] = {
                    a: dict(mean=float(Sy.Mean[a]), best=float(Sy.Best[a]), sd=float(Sy.SD[a]), loss_pct=float(Sy.Loss[a]),
                            feasible=int(Sy.NFeas[a]), runs=int(Sy.N[a]), rank=float(Sy.Rank[a]),
                            p_holm=tt.get(a, {}).get("PHolm"), rb=tt.get(a, {}).get("RB")) for a in mets}
                bs = f"{b:,}".replace(",", "{,}")
                lines.append(f"\\multicolumn{{8}}{{l}}{{\\emph{{{cname}, {bs} calls}}}} \\\\")
                for a in sorted(mets, key=lambda a: Sy.Rank[a]):
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
            rank_among = int(1 + (pf.AEP > fb).sum()) if pf is not None and len(pf) and np.isfinite(fb) else None
            isum[cname]["published"] = dict(
                baseline_aep=bl_aep, best_feasible_published=pf_aep, best_feasible_by=pf_who,
                best_overall_published=pa_aep, best_overall_by=pa_who, best_overall_feasible=pa_feas,
                n_published=int(len(parts)) if parts is not None else 0,
                n_published_feasible=int(len(pf)) if pf is not None else 0,
                our_best=float(top.Objective), our_best_method=top.Algorithm, our_best_budget=int(top.Budget), our_best_seed=int(top.Seed),
                focus_best=fb, focus_best_budget=int(fbest.Budget) if fbest is not None else None,
                focus_mean_by_budget={int(k): float(v) for k, v in fmean.items()},
                focus_best_vs_best_feasible_pct=float(100 * (fb / pf_aep - 1)) if np.isfinite(pf_aep) else None,
                focus_best_vs_baseline_pct=float(100 * (fb / bl_aep - 1)) if np.isfinite(bl_aep) else None,
                focus_mean_vs_baseline_pct={int(k): float(100 * (v / bl_aep - 1)) for k, v in fmean.items()},
                our_best_vs_best_feasible_pct=float(100 * (top.Objective / pf_aep - 1)) if np.isfinite(pf_aep) else None,
                focus_best_rank_among_feasible_published=rank_among)
            num = lambda v: "--" if not np.isfinite(v) else f"{v:,.1f}".replace(",", "{,}")
            plines.append(f"{cname} & {num(bl_aep)} & {num(pf_aep)} ({pf_who}) & {num(pa_aep)} ({pa_who}{'' if pa_feas in (None, True) else ', infeas.'}) & "
                          f"{num(top.Objective)} ({LAB[top.Algorithm]}) & {num(fb)} & "
                          + " / ".join(num(fmean[k]) for k in sorted(fmean)) +
                          f" & {'--' if not np.isfinite(pf_aep) else f'{100 * (fb / pf_aep - 1):+.2f}'} & "
                          f"{'--' if not np.isfinite(bl_aep) else f'{100 * (fb / bl_aep - 1):+.2f}'} \\\\")
            lay.append((n, cname, top, fbest, pf_who, x.Radius.iloc[0]))
        bl_txt = " / ".join(f"{b:,}".replace(",", "{,}") for b in sorted(IE.Budget.unique()))
        supp.append(table("table", "IEA Wind Task~37 Case Study~1 (16- and 36-turbine scenarios; official AEP model): AEP (MWh) mean, best and SD over the feasible runs, mean wake loss relative to the wake-free AEP (\\%%), feasible runs, rank (ranking rule of Table~\\ref{tab:friedman68}) and Holm-adjusted Wilcoxon signed-rank $p$ of %s vs.\\ each method, per scenario and budget." % fl,
                              "tab:iea37-detail", "lccccccc", "Method & Mean & Best & SD & Loss (\\%) & Feas. & Rank & $p_{\\rm Holm}$", lines, sep="2.5pt", pos="p"))
        supp.append(table("table*", "IEA Wind Task~37 Case Study~1: AEP (MWh, official calculator) of the baseline (example) layout, of the best feasible and the best overall participant layout%s, our best feasible layout over all methods and budgets, the best %s layout, mean %s AEP (budgets %s calls), and difference (\\%%) of the best %s layout from the best feasible participant layout and from the baseline. ``infeas.'': violates the boundary or spacing constraint by more than 1~mm." % ("" if pub is not None else " (published results file not available)", fl, fl, bl_txt, fl),
                                        "tab:iea37-pub", "lcccccccc", f"Scenario & Baseline & Best feasible publ. & Best publ. & Our best & {fl} best & {fl} mean & $\\Delta_{{\\rm publ}}$ & $\\Delta_{{\\rm base}}$", plines, sep="2.5pt", resize=True, pos="!htb"))
        summary["iea37"] = isum
        summary["iea37_published_file"] = pub is not None
        for k, v in isum.items():
            p = v.get("published", {})
            log(f"  {k}: {fl} best {p.get('focus_best', float('nan')):.1f}, best feasible published {p.get('best_feasible_published', float('nan')):.1f} "
                f"({p.get('best_feasible_by')}), baseline {p.get('baseline_aep', float('nan')):.1f}")
        # layouts: best focus vs best feasible participant layout (and our best if another method)
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
                    ax.scatter(xy[:, 0], xy[:, 1], marker="o", s=12, color=COL[FOCUS], label=f"best {fl}", zorder=4)
                if top.Algorithm != FOCUS:
                    xy = coords(top.Coordinates)
                    ax.scatter(xy[:, 0], xy[:, 1], marker=mk(top.Algorithm), s=16, facecolor="none", edgecolor=COL[top.Algorithm],
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
            mets = [a for a in M10 if a in set(y.Algorithm)]
            conv_panel(ax, y, mets, loss=True, band=False)
            ax.set_title(f"{IEA_NAME.get(n, n)}, {int(y.Budget.iloc[0]):,} evaluations", fontsize=8, color=INK)
            ax.set_xlabel("Evaluations"); ax.set_ylabel("Median best wake loss (%)")
        legend_row(fig, [a for a in M10 if a in set(IE.Algorithm)], ncol=5)
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
        supp.append(tabs.pop("robust").replace("[!t]", "[!htb]", 1).replace("{\\tabcolsep}{3pt}", "{\\tabcolsep}{1.8pt}"))

    # =========================================================== 7b. data checks for CHECK-FINAL (no tables)
    log("\n[7b] Data checks (seed pairing, loss monotonicity, boundary rule, provenance)")
    summary["pairing"] = pairing_check(ALL)
    summary["main"]["loss_monotone"] = loss_monotone(S, FOCUS)
    summary["boundary_rule"] = boundary_rule()
    if summary.get("spread") and (summary.get("hr16") or {}).get("methods"):      # reused: \NHRFeasPSOVNS, \NHRFeasPSO
        hm_ = summary["hr16"]["methods"]
        summary["spread"]["hr16_feasible"] = {a: dict(feasible=hm_[a]["feasible"], runs=hm_[a]["runs"])
                                              for a in ("PSOBV", "PSOC") if a in hm_}
    summary["provenance"] = {f"{a}|{ds}|{b}|{i}": sorted(set(x.Source))
                             for (a, ds, b, i), x in ALL.groupby(["Algorithm", "Dataset", "Budget", "Init"])}
    po_ = new.get("pso_old")
    pv = pd.concat([ALL] + ([po_.assign(Algorithm="PSO")] if po_ is not None and len(po_) else []), ignore_index=True)
    pv = pv.assign(Study=np.where(pv.Dataset.isin(["1", "2"]), "Benchmark",
                                  np.where(pv.Dataset == "HR", "Horns Rev 1", "IEA37")))
    PROV_LAB = {'PSOBV25': 'PSO-VNS ($' + chr(92) + 'omega=0.25$)', 'PSOBV75': 'PSO-VNS ($' + chr(92) + 'omega=0.75$)',
                'PSOBV90': 'PSO-VNS ($' + chr(92) + 'omega=0.9$)'}
    pvl = []
    for (st_, b, i), x in pv.groupby(["Study", "Budget", "Init"]):
        for src, y in x.groupby("Source"):
            pvl.append(f"{st_} & {b:,} & {i} & {', '.join(PROV_LAB.get(a, LAB.get(a, a)) for a in sorted(y.Algorithm.unique()))} & "
                       f"\\texttt{{{src.replace('_', chr(92) + '_')}}} & {len(y)} \\\\".replace(",", "{,}", 1))
    supp.append(table("table*", "Provenance of the per-run results used in this paper: study, budget (evaluations), initialization, methods, source file (\\texttt{mpce\\_<exp>}: all shards \\texttt{mpce\\_<exp>\\_s<i>of<k>.csv} of experiment \\texttt{<exp>} of \\texttt{mpce\\_experiments.py} / \\texttt{iea37\\_experiments.py}; \\texttt{fresh\\_*.csv}: earlier grid runs, see the repository README) and number of runs.",
                      "tab:provenance", "lccp{7.2cm}ll", "Study & Budget & Init. & Methods & Source & Runs",
                      [l_.replace("PSOBV25", "PSO-VNS (0.25)").replace("PSOBV75", "PSO-VNS (0.75)").replace("PSOBV90", "PSO-VNS (0.9)") for l_ in pvl], size="\\scriptsize", pos="p"))

    # =========================================================== 8b. main-text tables and figures in the manuscript layout
    log("\n[8b] Main-text tables (manuscript layout)")
    main_text_tables(ALL, R6, summary, tabs, FG, inst, hr, args)

    # =========================================================== 9. decision block (independent of --focus)
    log("\n[9] Decision block (PSO-VNS vs PSO and vs SSA-VNS; independent of --focus)")
    summary["decision"], dec_lines = decision_block(ALL)

    # =========================================================== write
    hdr = "%% generated by mpce_results.py -- do not edit by hand\n"
    for fn in glob.glob(OUT("mpce_tab_*.tex")):             # remove tables of sections skipped in this run
        if os.path.basename(fn)[9:-4] not in tabs and not (args.skip_robust and fn.endswith(("_robust.tex", "_robust_final.tex"))):
            os.remove(fn); log(f"  removed stale {os.path.basename(fn)}")
    # cross-document references (xr): the main text reads the supplement's labels with prefix "S-" and the
    # supplement (MPCE_PSO_VNS_supplement.tex) reads the main text's labels with prefix "M-"
    supp_lbl = set(re.findall(r"\\label\{([^}]*)\}", "\n".join(supp)))
    main_lbl = (set(re.findall(r"\\label\{([^}]*)\}", "\n".join(tabs.values()))) | {"tab:friedman68"}) - supp_lbl
    xref = lambda txt, lbls, pre: re.sub(r"\\ref\{([^}]*)\}", lambda m_: "\\ref{%s%s}" % (pre if m_.group(1) in lbls else "", m_.group(1)), txt)
    tabs = {k: xref(v, supp_lbl, "S-") for k, v in tabs.items()}
    supp = [xref(v, main_lbl, "M-") for v in supp]
    tabs = {k: evals(v) for k, v in tabs.items()}
    supp = [evals(v) for v in supp]
    for k, v in tabs.items():
        open(OUT(f"mpce_tab_{k}.tex"), "w").write(hdr + v)
    open(OUT("mpce_supplementary.tex"), "w").write(
        hdr + "%% Supplementary tables (per-case results); requires booktabs, graphicx\n" + "\n".join(supp))
    summary["tables"] = [f"mpce_tab_{k}.tex" for k in tabs] + ["mpce_supplementary.tex"]
    summary["figures"] = [f"figures_mpce/{f}.pdf" for f in FG.made]
    summary["log"] = LOG
    json.dump(clean(summary), open(OUT("mpce_summary.json"), "w"), indent=1)
    log(f"\nwrote {len(tabs)} tables, {len(FG.made)} figures, mpce_summary.json in {time.time() - t0:.0f} s")
    print("\n" + "=" * 78)
    for ln in dec_lines:
        print(ln)
    print("=" * 78, flush=True)
    summary["decision"]["printed"] = dec_lines
    json.dump(clean(summary), open(OUT("mpce_summary.json"), "w"), indent=1)
    # numbers quoted in the manuscript (mpce_numbers.tex) and the check of the data-dependent statements
    import mpce_numbers, mpce_check_final
    mpce_numbers.main(["--summary", OUT("mpce_summary.json"), "--out", OUT("mpce_numbers.tex")])
    mpce_check_final.main(["--summary", OUT("mpce_summary.json")])


def pairing_check(ALL):
    """Seed pairing (CHECK-FINAL C43): within every case x budget x initialization, every method has the same
    set of seeds; at 6,030 evaluations the first convergence checkpoint (best feasible objective among the 30
    initial layouts, i.e. the common initial population) is identical for all methods of a case-seed pair
    (compared as stored; checked only where it is finite -- infeasible initial populations are stored as nan).
    Split variants (*25 / *75) are included. Other budgets log their first checkpoint later than call 30."""
    grp = CASE + ["Budget", "Init"]
    sets = ALL.groupby(grp + ["Algorithm"]).Seed.agg(lambda v: tuple(sorted(set(v))))
    nset = sets.groupby(level=list(range(len(grp)))).nunique()
    X = ALL[ALL.Budget == 6030]
    c0 = pd.to_numeric(X.Curve.str.split(";").str[0], errors="coerce")
    X = X.assign(C0=c0)[np.isfinite(c0)]
    spread = X.groupby(CASE + ["Init", "Seed"]).C0.agg(lambda v: float(v.max() - v.min()))
    nm = X.groupby(CASE + ["Init", "Seed"]).Algorithm.nunique()
    out = dict(groups=int(len(nset)), groups_with_different_seed_sets=int((nset > 1).sum()),
               case_seed_pairs_6030=int(len(spread)), pairs_compared=int((nm > 1).sum()),
               pairs_first_checkpoint_differs=int(((spread > 1e-9) & (nm > 1)).sum()),
               max_first_checkpoint_spread=float(spread.max()) if len(spread) else None)
    log(f"  seed pairing: {out}")
    return out


def loss_monotone(S, a):
    """Monotonicity of the case-mean wake loss of method a (CHECK-FINAL C42): violations of (i) loss
    non-decreasing in N within each data set and radius, (ii) loss non-increasing in r at fixed N and data set,
    (iii) Data Set II >= Data Set I at the same (r, N); ties within 1e-6 pp allowed."""
    L = S[(S.Algorithm == a) & S.Qualified].set_index(CASE).Loss.sort_index()
    tol = 1e-6
    vn, vr, vd, npair, dsgt = [], [], [], 0, 0
    for (d, r), x in L.groupby(level=[0, 1]):
        x = x.sort_index(level=2)
        vn += [f"{d}-{r}-{n}" for n, dv in zip(x.index.get_level_values(2)[1:], np.diff(x.values)) if dv < -tol]
    for (d, n), x in L.groupby(level=[0, 2]):
        x = x.sort_index(level=1)
        vr += [f"{d}-{rr}-{n}" for rr, dv in zip(x.index.get_level_values(1)[1:], np.diff(x.values)) if dv > tol]
    for (r, n), x in L.groupby(level=[1, 2]):
        if x.index.get_level_values(0).nunique() == 2:
            npair += 1
            l1, l2 = float(x.xs("1", level=0).iloc[0]), float(x.xs("2", level=0).iloc[0])
            dsgt += int(l2 > l1 + tol)
            if l2 < l1 - tol:
                vd.append(f"{r}-{n}")
    return dict(method=a, violations_n=vn, violations_r=vr, violations_ds=vd, ds_pairs=npair, ds2_higher=dsgt,
                ds_ties=npair - dsgt - len(vd))


def boundary_rule():
    """Boundary-handling check quoted in the robustness section (CHECK-FINAL C41): SSA with radial projection
    onto the circle (ssa_reference.py, ssa_reference_runs.csv: 4D, 3,030 evaluations, non-greedy) versus the
    archived runs of the original SSA with box clipping (../selected_30_run_data.csv, Algorithm SSA,
    EnergyProduction), six cases x 30 runs; two-sided Mann-Whitney U on the feasible runs' objectives."""
    from scipy.stats import mannwhitneyu
    fr, fa = os.path.join(HERE, "ssa_reference_runs.csv"), os.path.join(HERE, "..", "selected_30_run_data.csv")
    if not (os.path.exists(fr) and os.path.exists(fa)):
        log("  boundary rule: input files missing"); return None
    R_, A_ = pd.read_csv(fr), pd.read_csv(fa)
    R_ = R_[(R_.Spacing == "4D") & (R_.Evaluations == 3030) & (~R_.Greedy.astype(bool))]
    rows = []
    for (d, r, n), x in R_.groupby(["Dataset", "Radius", "Turbines"]):
        y = A_[(A_.Algorithm == "SSA") & (A_.Dataset == d) & (A_.Radius == r) & (A_.Turbines == n)].EnergyProduction
        xf = x[x.Feasible.astype(bool)].Objective
        rows.append(dict(case=f"{d}-{r}-{n}", n_proj=int(len(xf)), n_orig=int(len(y)), diff=float(xf.mean() - y.mean()),
                         p=float(mannwhitneyu(xf, y).pvalue)))
    out = dict(cases=rows, n_cases=len(rows), n_higher=int(sum(r_["diff"] > 0 for r_ in rows)),
               max_p=max(r_["p"] for r_ in rows), min_diff=min(r_["diff"] for r_ in rows), max_diff=max(r_["diff"] for r_ in rows))
    log(f"  boundary rule: projection higher in {out['n_higher']} of {out['n_cases']} cases, max p {out['max_p']:.2g}")
    return out


PEND = "\\TBD{}"            # table cell whose data are still missing (the manuscript defines \TBD)
NOFEAS = {"RSVNS"}          # not run with feasible initialization (one packing solve per random sample)
NR = "n/r"                  # table cell: not run


def _complete(x, methods, cases, runs):
    """dict method -> True if x holds `runs` runs of the method in every case of `cases`."""
    out = {}
    for a in methods:
        xa = x[x.Algorithm == a]
        out[a] = all(((xa.Dataset == d) & (xa.Radius == r) & (xa.Turbines == n)).sum() >= runs for d, r, n in cases)
    return out


def main_text_tables(ALL, R6, summary, tabs, FG, inst, hr, args):
    """Tables of the manuscript whose layout differs from the per-section tables: merged feasible-initialization
    + budget table (tab:feasbudget), Horns Rev 16 table with ten methods and loss columns for all settings
    (tab:hr-site), compact IEA37 table (tab:iea37), single-column robustness table (tab:robust-final), and the
    combined IEA37 + Horns Rev layout figure (layouts_iea37_hr16). Missing data -> \TBD{} cells."""
    meth = [a for a in HR_ORDER if a in M10]
    cases6 = LARGE
    # ---------------- feasible initialization + budget (six largest benchmark cases)
    rows, fb = [], dict(methods=meth, random={}, feasible={}, rank={}, complete={})
    sel = lambda X: X[[(d, r, n) in cases6 for d, r, n in zip(X.Dataset, X.Radius, X.Turbines)]]
    X6 = sel(R6[R6.Algorithm.isin(meth)])
    F6 = sel(ALL[(ALL.Init == "feasible") & (ALL.Budget == 6030) & ALL.Algorithm.isin(meth)])
    def loss_feas(X, a):
        xa = X[X.Algorithm == a]
        if not len(xa):
            return None
        Sa = case_stats(xa, [a])
        return dict(loss=float(Sa[Sa.Qualified].Loss.mean()) if Sa.Qualified.any() else None,
                    feas=float(100 * xa.Feasible.mean()), n=int(len(xa)))
    cr = _complete(X6, meth, cases6, 30)
    methF = [a for a in meth if a not in NOFEAS]
    cf = _complete(F6, methF, cases6, 30); cf.update({a: False for a in meth if a in NOFEAS})
    fb["complete"]["random_6030"] = all(cr.values()); fb["complete"]["feasible_6030"] = all(cf[a] for a in methF)
    ranks = {}
    for b in BUDGETS:
        Xb = sel(ALL[(ALL.Init == "random") & (ALL.Budget == b) & ALL.Algorithm.isin(meth)])
        cb = _complete(Xb, meth, cases6, 30)
        fb["complete"][f"rank_{b}"] = all(cb.values())
        if all(cb.values()):
            Sb = case_stats(Xb, meth)
            ranks[b] = rank_matrix(Sb, meth).mean(axis=0).to_dict()
            fb["rank"][b] = ranks[b]
            fb.setdefault("mean_loss", {})[b] = {a: float(Sb[(Sb.Algorithm == a) & Sb.Qualified].Loss.mean()) for a in meth}
            fb.setdefault("feas_pct", {})[b] = {a: float(100 * Xb[Xb.Algorithm == a].Feasible.mean()) for a in meth}
            # cases (of the six) in which each method has the best (lowest) rank; ties for first counted for all
            # tied methods (first_count) or not at all (sole_first_count)
            Rb_ = rank_matrix(Sb, meth)
            top = Rb_.eq(Rb_.min(axis=1), axis=0)
            fb.setdefault("first_count", {})[b] = {a: int(top[a].sum()) for a in meth}
            fb.setdefault("sole_first_count", {})[b] = {a: int((top[a] & (top.sum(axis=1) == 1)).sum()) for a in meth}
            ml = fb["mean_loss"][b]
            close = [a for a in ("BVNS", "SSABV", "RSVNS", "LXBV", "DE") if a in ml and np.isfinite(ml[a])]
            fb.setdefault("close_gap_pp", {})[b] = dict(methods=close, gap={a: float(ml[a] - ml[FOCUS]) for a in close},
                                                       max_gap=float(max(ml[a] - ml[FOCUS] for a in close)) if close else None)
    f3 = lambda v: "--" if v is None or not np.isfinite(v) else f"{v:.3f}"      # data present, no qualified case
    for a in meth:
        r_ = loss_feas(X6, a) if cr[a] else None
        f_ = loss_feas(F6, a) if cf[a] else None
        fb["random"][a] = r_; fb["feasible"][a] = f_
        c = [f3(r_["loss"]) if r_ else PEND, f"{r_['feas']:.1f}" if r_ else PEND,
             f3(f_["loss"]) if f_ else (NR if a in NOFEAS else PEND), f"{f_['feas']:.1f}" if f_ else (NR if a in NOFEAS else PEND)]
        c += [f"{ranks[b][a]:.2f}" if b in ranks else PEND for b in BUDGETS]
        rows.append(f"{LAB[a]} & " + " & ".join(c) + " \\\\")
    if fb["complete"]["feasible_6030"]:
        fb["rank_feasible_init"] = rank_matrix(case_stats(F6[F6.Algorithm.isin(methF)], methF), methF).mean(axis=0).to_dict()
    pend = [k for k, v in fb["complete"].items() if not v]
    short = {"random_6030": "6k", "feasible_6030": "feas", "rank_6030": "6k", "rank_30030": "b30k", "rank_120030": "b120k"}
    foot = ["\\multicolumn{8}{l}{\\TBD{pending: %s}}" % ", ".join(dict.fromkeys(short[k] for k in pend))] if pend else None
    fb["pending"] = pend
    for b in ranks:
        o = sorted(ranks[b], key=ranks[b].get)
        fb.setdefault("best", {})[b] = o[0]
        fb.setdefault("focus_position", {})[b] = 1 + o.index(FOCUS) if FOCUS in o else None
    summary["feasbudget"] = fb
    tabs["feasbudget"] = (
        "\\begin{table}[!t]\n\\centering\n\\caption{Initialization and Budget on the Six Largest Benchmark Cases (30 Seeds)}\n"
        "\\label{tab:feasbudget}\n\\scriptsize\\setlength{\\tabcolsep}{2.4pt}\n\\begin{tabular}{lccccccc}\n\\toprule\n"
        "& \\multicolumn{2}{c}{Random init.} & \\multicolumn{2}{c}{Feasible init.} & \\multicolumn{3}{c}{Avg.\\ rank} \\\\\n"
        "\\cmidrule(lr){2-3}\\cmidrule(lr){4-5}\\cmidrule(lr){6-8}\n"
        "Method & Loss & Feas. & Loss & Feas. & 6,030 & 30,030 & 120,030 \\\\\n\\midrule\n" + "\n".join(rows) +
        "\n\\bottomrule\n" + ("\n".join(foot) + "\n" if foot else "") + "\\end{tabular}\n" +
        tnote("Loss: mean wake loss (\\%) of the feasible runs; Feas.: feasible runs (\\%); both at 6,030 calls with random "
              "and feasibility-preserving initialization. Avg.\\ rank at 6,030, 30,030 and 120,030 calls, random "
              "initialization. n/r: not run.") + "\\end{table}\n")
    log(f"  feasbudget: pending {pend}; ranks {ranks}")

    # ---------------- Horns Rev 16 (ten methods; loss at 6k random / 6k feasible / 30k / 120k)
    hsum = summary.get("hr16") or {}
    H = ALL[(ALL.Dataset == "HR") & (ALL.Turbines == 16)]
    settings = [("6030R", 6030, "random", 30), ("6030F", 6030, "feasible", 30), ("30030R", 30030, "random", 30),
                ("120030R", 120030, "random", 10)]
    hl, hloss = [], {}
    ideal = hsum.get("ideal_aep", np.nan)
    hl.append(f"Installed layout & {inst:.2f} & -- & -- & {100 * (1 - inst / ideal):.2f} & -- & -- & -- \\\\" if np.isfinite(ideal) else
              f"Installed layout & {inst:.2f} & -- & -- & {PEND} & -- & -- & -- \\\\")
    hl.append("\\midrule")
    for a in meth:
        m6 = (hsum.get("methods") or {}).get(a)
        c = [f"{m6['mean']:.2f}" if m6 and m6.get("feasible") else ("--" if m6 else PEND),
             f"{m6['feasible']}/{m6['runs']}" if m6 else PEND,
             "--" if a == FOCUS else (fmt_p(m6["p_holm"], 2) if m6 and m6.get("p_holm") is not None else PEND)]
        hloss[a] = {}
        for tag, b, init, nrun in settings:
            y = H[(H.Algorithm == a) & (H.Budget == b) & (H.Init == init)]
            if init == "feasible" and a in NOFEAS:
                c.append(NR); hloss[a][tag] = None; continue
            if len(y) < nrun:
                c.append(PEND); hloss[a][tag] = None; continue
            nf = int(y.Feasible.sum())
            v = float(y[y.Feasible].LossPct.mean()) if nf else np.nan
            hloss[a][tag] = dict(loss=v if np.isfinite(v) else None, feasible=nf, runs=int(len(y)),
                                 mean_aep=float(y[y.Feasible].Objective.mean()) if nf else None,
                                 runs_above_installed=int((y.Feasible & (y.Objective > inst)).sum()))
            cell = "--" if not nf else f"{v:.2f}"
            if 0 < nf < len(y):
                cell += f"$^{{{nf}}}$"
            c.append(cell)
        hl.append(f"{LAB[a]} & " + " & ".join(c) + " \\\\")
    if hsum:
        hsum["loss_by_setting"] = hloss
    tabs["hr16"] = (
        "\\begin{table}[!t]\n\\centering\n\\caption{Horns Rev~1 16-Turbine Block (5\\textdegree{} Direction Bins): AEP, Feasibility and AEP Loss}\n"
        "\\label{tab:hr-site}\n\\scriptsize\\setlength{\\tabcolsep}{2.2pt}\n\\begin{tabular}{lccccccc}\n\\toprule\n"
        "& \\multicolumn{3}{c}{6,030 calls, R} & \\multicolumn{4}{c}{Loss (\\%)} \\\\\n\\cmidrule(lr){2-4}\\cmidrule(lr){5-8}\n"
        "Layout / method & AEP & Feas. & $p_{\\rm Holm}$ & 6k R & 6k F & 30k & 120k \\\\\n\\midrule\n" + "\n".join(hl) +
        "\n\\bottomrule\n\\end{tabular}\n" +
        tnote("Wake-free AEP %.2f GWh/yr. AEP: mean (GWh/yr) over the feasible runs of 30 seeds at 6,030 calls, random "
              "initialization (R); Feas.: feasible runs; $p_{\\rm Holm}$: run-level Wilcoxon signed-rank $p$ of %s vs.\\ each "
              "method (infeasible runs ranked last), Holm-adjusted over the %d comparisons. Loss: mean AEP loss (\\%%) at "
              "6,030 calls with random (R) and feasibility-preserving (F) initialization and at 30,030 and 120,030 calls "
              "(R; 10 seeds at 120,030); superscript: feasible runs when not all are feasible; n/r: not run."
              % (ideal, LAB[FOCUS], len(meth) - 1)) + "\\end{table}\n")

    # ---------------- IEA37 compact table (best feasible run of each method)
    isum = summary.get("iea37") or {}
    pub = load_published(args.data_dir)
    IE = ALL[ALL.Dataset.str.startswith("IEA37")]
    il, ib = [], {}
    pubv = {}
    for n in (16, 36):
        P = pub[pub.Turbines == n] if pub is not None else None
        base = float(P[P.Baseline].AEP.iloc[0]) if P is not None and P.Baseline.any() else BASELINE_IEA[n]
        pf = P[(~P.Baseline) & P.Feasible] if P is not None else None
        pubv[n] = dict(base=base, best_feasible=float(pf.AEP.max()) if pf is not None and len(pf) else np.nan,
                       best_feasible_by=str(pf.sort_values("AEP").Participant.iloc[-1]) if pf is not None and len(pf) else "--",
                       feasible_aeps=sorted(map(float, pf.AEP)) if pf is not None else [])
    num = lambda v: PEND if v is None or not np.isfinite(v) else f"{v:,.1f}".replace(",", "{,}")
    il.append("Example layout & \\multicolumn{2}{c}{%s} & \\multicolumn{2}{c}{%s} \\\\" % (num(pubv[16]["base"]), num(pubv[36]["base"])))
    il.append("Best published, strict$^{a}$ & \\multicolumn{2}{c}{%s} & \\multicolumn{2}{c}{%s} \\\\" % (num(pubv[16]["best_feasible"]), num(pubv[36]["best_feasible"])))
    # R4-5 / D19: second convention -- published layouts projected radially onto the boundary (iea37_projected.py ->
    # iea37_projected.json, read only; the \NF... macros of mpce_numbers_dir.tex come from the same file). Without
    # the file the table keeps the strict row only and says so in the note.
    proj = None
    try:
        J = json.load(open(os.path.join(HERE, "iea37_projected.json")))
        proj = {n: J["scenarios"][str(n)] for n in (16, 36)}
        for n in (16, 36):
            if abs(proj[n]["strict"]["best_aep"] - pubv[n]["best_feasible"]) > 0.5:
                log(f"  WARNING: iea37_projected.json strict best ({proj[n]['strict']['best_aep']:.1f}) differs from the published file "
                    f"({pubv[n]['best_feasible']:.1f}) for {n} turbines")
            pubv[n]["best_projected"] = float(proj[n]["projected"]["best_aep"])
            pubv[n]["best_projected_by"] = str(proj[n]["projected"]["best_by"])
            pubv[n]["best_projected_excess_m"] = float(proj[n]["projected"]["best_max_excess_m"])
        il.append("Best published, projected$^{b}$ & \\multicolumn{2}{c}{%s} & \\multicolumn{2}{c}{%s} \\\\"
                  % (num(pubv[16]["best_projected"]), num(pubv[36]["best_projected"])))
    except (OSError, KeyError, ValueError) as e:
        log(f"  IEA37: projected published layouts not available ({e}); strict convention only")
        proj = None
    il.append("\\midrule")
    for a in MAIN8:
        c = []
        for n in (16, 36):
            for b in (6030, 30030):
                y = IE[(IE.Turbines == n) & (IE.Budget == b) & (IE.Algorithm == a)]
                if len(y) < 30:
                    c.append(PEND); ib.setdefault(a, {})[f"{n}T_{b}"] = None; continue
                v = float(y[y.Feasible].Objective.max()) if y.Feasible.any() else np.nan
                mv = float(y[y.Feasible].Objective.mean()) if y.Feasible.any() else np.nan
                ib.setdefault(a, {})[f"{n}T_{b}"] = dict(best=v, mean=mv, feasible=int(y.Feasible.sum()), runs=int(len(y)),
                                                        best_rank_among_feasible_published=int(1 + sum(p > v for p in pubv[n]["feasible_aeps"])) if np.isfinite(v) else None,
                                                        mean_rank_among_feasible_published=int(1 + sum(p > mv for p in pubv[n]["feasible_aeps"])) if np.isfinite(mv) else None,
                                                        best_gap_to_best_feasible_published_pct=float(100 * (v / pubv[n]["best_feasible"] - 1)) if np.isfinite(v) else None,
                                                        mean_gap_to_best_feasible_published_pct=float(100 * (mv / pubv[n]["best_feasible"] - 1)) if np.isfinite(mv) else None)
                c.append("--" if not np.isfinite(v) else num(v))
        il.append(f"{LAB[a]} & " + " & ".join(c) + " \\\\")
    summary["iea37_compact"] = dict(published={str(k): {kk: vv for kk, vv in v.items() if kk != "feasible_aeps"} | dict(n_feasible=len(v["feasible_aeps"]))
                                               for k, v in pubv.items()}, methods=ib)
    def whof(key):
        w_ = {pubv[16][key], pubv[36][key]}
        return (("participant~%s in both scenarios" % next(iter(w_)).replace("par", "")) if len(w_) == 1 else
                "participants~%s (16) and %s (36)" % (pubv[16][key].replace("par", ""), pubv[36][key].replace("par", "")))
    def dist(m):
        return f"{m:.1f}~m" if m >= 0.1 else f"{1000 * m:.1f}~mm"
    inote = ("AEP of the best feasible run of each method (30 seeds; --: no feasible run). $^{a}$Best published layout with "
             "every turbine within 1~mm of the boundary (%s). " % whof("best_feasible_by"))
    if proj is not None:
        inote += ("$^{b}$Best published layout after projecting the turbines outside the boundary radially onto it, "
                  "spacing re-checked (%s; up to %s and %s outside before projection)."
                  % (whof("best_projected_by"), dist(pubv[16]["best_projected_excess_m"]), dist(pubv[36]["best_projected_excess_m"])))
    else:
        inote += "Layouts that place turbines outside the boundary by more than 1~mm are not counted."
    tabs["iea37"] = (
        "\\begin{table}[!t]\n\\centering\n\\caption{IEA37 Case Study~1: AEP (MWh) of the Best Run of Each Method and of Published Layouts~\\cite{Baker2019,IEA37repo}}\n"
        "\\label{tab:iea37}\n\\scriptsize\\setlength{\\tabcolsep}{2.5pt}\n\\begin{tabular}{lcccc}\n\\toprule\n"
        "& \\multicolumn{2}{c}{16 turbines ($r=1300$~m)} & \\multicolumn{2}{c}{36 turbines ($r=2000$~m)} \\\\\n\\cmidrule(lr){2-3}\\cmidrule(lr){4-5}\n"
        "Method / source & 6,030 & 30,030 & 6,030 & 30,030 \\\\\n\\midrule\n" + "\n".join(il) + "\n\\bottomrule\n"
        "\\end{tabular}\n" + tnote(inote) + "\\end{table}\n")

    # ---------------- robustness, single column
    rb = summary.get("robustness")
    if rb:
        rl = []
        for col, name in (("Objective", "Benchmark"), ("Cubic", "Cubic curve"), ("CubicCutout", "Cubic, 25 m/s cut-out"), ("Gauss", "Gaussian wake")):
            v = rb.get(col)
            if not v:
                continue
            ar = v["avg_rank"]; best = min(ar, key=ar.get)
            v["best_method"] = best
            v["focus_rank_position"] = int(1 + sorted(ar.values()).index(ar[FOCUS]))
            d1, d2 = v["rel_change_pct"]["1"], v["rel_change_pct"]["2"]
            rl.append(f"{name} & {d1:+.1f} & {d2:+.1f} & {ar[FOCUS]:.2f} & {LAB[best]} & " +
                      ("-- & --" if col == "Objective" else f"{v['tau']:.2f} & {v['same_best_pct']:.0f}") + " \\\\")
        tabs["robust_final"] = (
            "\\begin{table}[!t]\n\\centering\n\\caption{Robustness to the Benchmark Model: Final Layouts Re-Evaluated, Not Re-Optimized}\n"
            "\\label{tab:robust-final}\n\\scriptsize\\setlength{\\tabcolsep}{2.4pt}\n\\begin{tabular}{lcccccc}\n\\toprule\n"
            "& \\multicolumn{2}{c}{$\\Delta$ (\\%)} & & & & \\\\\n\\cmidrule(lr){2-3}\nModel & DS I & DS II & Rank & Best & $\\bar\\tau$ & Same \\\\\n\\midrule\n" +
            "\n".join(rl) + "\n\\bottomrule\n\\end{tabular}\n" +
            tnote("All feasible final layouts of the %d cases. $\\Delta$: mean relative change of the objective (\\%%); Rank: "
                  "average rank of %s; Best: best-ranked method; $\\bar\\tau$: mean Kendall correlation between the benchmark "
                  "and the alternative ordering within a case; Same: cases (\\%%) with unchanged best method. Ranks of all "
                  "methods: supplementary material." % (68, LAB[FOCUS])) + "\\end{table}\n")

    # ---------------- combined layout figure: IEA37 16 / 36 (best focus layout at the largest budget) + Horns Rev 16
    fig, axes = plt.subplots(1, 3, figsize=(7.1, 2.75))
    for ax, n in zip(axes[:2], (16, 36)):
        rad = {16: 1300, 36: 2000}[n]
        t = np.linspace(0, 2 * np.pi, 200)
        ax.plot(rad * np.cos(t), rad * np.sin(t), color=MUTED, lw=0.8)
        who_n = pubv[n]["best_feasible_by"]
        try:
            import iea37_model as iem
            xy, _ = iem.load_submission(str(who_n).replace("par", ""), n)
            ax.scatter(xy[:, 0], xy[:, 1], marker="s", s=20, facecolor="none", edgecolor=INK, lw=0.9,
                       label=f"best feasible published ({num(pubv[n]['best_feasible']).replace('{,}', ',')})", zorder=3)
        except Exception as e:                                  # pragma: no cover
            log(f"  participant layout not drawn: {e}")
        y = IE[(IE.Turbines == n) & (IE.Algorithm == FOCUS) & IE.Feasible]
        if len(y):
            y = y[y.Budget == y.Budget.max()]
            bst = y.sort_values("Objective").iloc[-1]; xy = coords(bst.Coordinates)
            ax.scatter(xy[:, 0], xy[:, 1], marker="o", s=12, color=COL[FOCUS],
                       label=f"best {LAB[FOCUS]}, {int(bst.Budget):,} evaluations ({bst.Objective:,.1f})", zorder=4)
        else:
            ax.text(0, 0, f"{LAB[FOCUS]} runs pending", ha="center", va="center", fontsize=7, color="#c00000")
        ax.set_aspect("equal"); ax.tick_params(labelsize=6)
        ax.set_title(f"IEA37 CS1, {n} turbines (AEP in MWh)", fontsize=7.5, color=INK)
        ax.legend(fontsize=5.8, loc="upper left", bbox_to_anchor=(0, -0.08), ncol=1)
    ax = axes[2]
    H6 = H[(H.Budget == 6030) & (H.Init == "random") & H.Feasible]
    if hr is not None:
        xy0, poly = hr.site(16)
        pp = np.vstack([poly, poly[:1]])
        ax.plot(pp[:, 0], pp[:, 1], color=MUTED, lw=0.8)
        ax.scatter(xy0[:, 0], xy0[:, 1], marker="x", s=14, color=INK, lw=0.8, label=f"installed ({inst:.2f})", zorder=3)
    y = H6[H6.Algorithm == FOCUS]
    if len(y):
        bst = y.sort_values("Objective").iloc[-1]; xy = coords(bst.Coordinates)
        ax.scatter(xy[:, 0], xy[:, 1], marker="o", s=12, color=COL[FOCUS], label=f"best {LAB[FOCUS]}, 6,030 evaluations ({bst.Objective:.2f})", zorder=4)
    ax.set_aspect("equal"); ax.tick_params(labelsize=6)
    ax.set_title("Horns Rev 1, 16 turbines (AEP in GWh/yr)", fontsize=7.5, color=INK)
    ax.legend(fontsize=5.8, loc="upper left", bbox_to_anchor=(0, -0.08), ncol=1)
    fig.tight_layout()
    FG.save(fig, "layouts_iea37_hr16")


BASE_POOL = ["LXBV", "SSABV", "LXSSA", "SSA", "PSO", "DE", "BVNS", "SLSQP"]    # previous study's methods ("PSO" = pool's PSO)
BASE_HYB = ["LXBV", "SSABV"]                                                   # salp-swarm hybrids compared with PSO (W/T/L)


def baseline_setting(R6, PO, tabs):
    """Previous study's method pool (68 cases, 6,030 calls) with the old PSO setting (w = 0.7, c1 = c2 = 2;
    fresh_grid.csv "PSO") and with the constriction setting ("PSOC"): average ranks (feasibility-aware rule),
    best method, PSO rank position, and run-level W/T/L of LX-SSA-VNS and SSA-VNS against PSO (seed-paired
    Wilcoxon, Holm over the seven comparisons of the hybrid in each case, as in the main W/T/L table).
    Writes tabs["baseline"] (tab:baseline); returns summary["baseline"] (None if a method has no data)."""
    if PO is None or not len(PO):
        log("  baseline setting: old-setting PSO runs missing -> section skipped")
        return None
    PO = PO[PO.Dataset.isin(["1", "2"]) & (PO.Budget == 6030)]
    B6 = R6[R6.Dataset.isin(["1", "2"])]
    out = dict(pool=BASE_POOL, hybrids=BASE_HYB, labels={a: LAB["PSOC" if a == "PSO" else a] for a in BASE_POOL},
               note="'PSO' denotes the PSO of the respective pool: old = w 0.7, c1 = c2 = 2 (fresh_grid PSO); "
                    "constriction = Clerc-Kennedy (PSOC). W/T/L from the hybrid's side, Holm over the 7 comparisons "
                    "of the hybrid within each case.")
    for tag, pso in (("old", "PSO"), ("constriction", "PSOC")):
        pool = [pso if a == "PSO" else a for a in BASE_POOL]
        D = pd.concat([B6[B6.Algorithm.isin([a for a in pool if a != "PSO"])], PO] if pso == "PSO"
                      else [B6[B6.Algorithm.isin(pool)]], ignore_index=True)
        miss = [a for a in pool if a not in set(D.Algorithm)]
        if miss:
            log(f"  baseline setting ({tag}): no data for {miss} -> section skipped")
            return None
        S = case_stats(D, pool)
        R = rank_matrix(S, pool)
        FB = friedman_block(R, pool, focus="LXBV")
        case_ranks = {f"{d}-{r}-{n}": {("PSO" if a == pso else a): float(R.loc[(d, r, n), a]) for a in (pso, "SSABV", "LXBV")}
                      for (d, r, n) in R.index}
        ar = {("PSO" if a == pso else a): v for a, v in FB["avg_rank"].items()}
        wt = {}
        for h in BASE_HYB:
            oc = []
            for _, sub in D.groupby(CASE):
                for x in paired_vs(sub, h, [a for a in pool if a != h]):
                    if x["Baseline"] == pso:
                        oc.append(x["Outcome"])
            wt[h] = pd.Series(oc, dtype=object).value_counts().to_dict()
        best = min(ar, key=ar.get)
        out[tag] = dict(pso_code=pso, n_cases=int(FB["n_cases"]), n_runs=int(len(D)), avg_rank=ar, best=best,
                        best_label=out["labels"][best], pso_avg_rank=float(ar["PSO"]),
                        pso_position=int(1 + sorted(ar.values()).index(ar["PSO"])),
                        pso_feasible_pct=float(100 * D[D.Algorithm == pso].Feasible.mean()),
                        wtl_vs_pso={h: dict(W=int(wt[h].get("W", 0)), T=int(wt[h].get("T", 0)), L=int(wt[h].get("L", 0)))
                                    for h in BASE_HYB},
                        friedman_chi2=float(FB["chi2"]), friedman_p=float(FB["p"]), case_ranks=case_ranks)
        log(f"  baseline pool ({tag} PSO): " + ", ".join(f"{out['labels'][a]} {v:.2f}" for a, v in sorted(ar.items(), key=lambda t: t[1]))
            + f"; PSO position {out[tag]['pso_position']}; " + ", ".join(f"{LAB[h]} vs PSO {wtl_str(out[tag]['wtl_vs_pso'][h])}" for h in BASE_HYB))
    o, c = out["old"], out["constriction"]
    lines = []
    for a in BASE_POOL:
        cells = []
        for t in (o, c):
            v = f"{t['avg_rank'][a]:.2f}"
            cells.append(f"\\textbf{{{v}}}" if a == t["best"] else v)
        lines.append(f"{out['labels'][a]} & " + " & ".join(cells) + " \\\\")
    lines.append("\\midrule")
    lines.append("PSO position & %d & %d \\\\" % (o["pso_position"], c["pso_position"]))
    for h in BASE_HYB:
        lines.append(f"{LAB[h]} vs.\\ PSO (W/T/L) & {wtl_str(o['wtl_vs_pso'][h])} & {wtl_str(c['wtl_vs_pso'][h])} \\\\")
    tabs["baseline"] = table(
        "table", "Effect of the PSO Setting on an Eight-Method Pool (the Methods of Table~\\ref{tab:friedman68} with "
        "LX-SSA-VNS Instead of PSO-VNS; %d Cases, 6,030 Evaluations)" % o["n_cases"],
        "tab:baseline", "lcc",
        "& \\multicolumn{2}{c}{PSO setting} \\\\\n\\cmidrule(lr){2-3}\nMethod & Old & Constriction", lines,
        foot=["\\multicolumn{3}{p{0.95\\columnwidth}}{Old: $w=0.7$, $c_1=c_2=2$; constriction: Clerc--Kennedy setting~\\cite{Clerc2002}. "
              "Average rank: 1 = best, feasibility-aware rule (fewer than 15 of 30 feasible runs "
              "in a case = ranked last); bold: best average rank. W/T/L: cases in which the hybrid is significantly "
              "better / not different / worse than PSO (two-sided Wilcoxon signed-rank test, 30 seed-paired runs, "
              "Holm-adjusted over the seven comparisons of the hybrid in each case, $\\alpha=0.05$).}"], sep="4pt",
        size="\\footnotesize")
    return out


def decision_block(ALL):
    """Compact PSO-VNS vs PSO and PSO-VNS vs SSA-VNS comparison (independent of --focus) for choosing the
    proposed method. Pool: the eight methods of the main comparison (MAIN8) wherever they have data.
      * 68 cases (6,030 calls): average ranks of the 8-method comparison, case-mean Wilcoxon on the per-case
        mean wake losses (two-sided p, rank-biserial r, mean dL = L(PSO-VNS) - L(other), pp; negative =
        PSO-VNS better), run-level W/T/L (30 seed-paired runs per case, Holm within the case over the 7
        comparisons of PSO-VNS with the other methods).
      * six largest cases + HR16 at 6,030 / 30,030 / 120,030 calls: the same statistics over these 7 cases.
      * IEA37 CS1 16 / 36 turbines at each budget: mean AEP (MWh) of the feasible runs, rank among the 8,
        run-level Wilcoxon (Holm over the 7 comparisons) of PSO-VNS vs PSO / SSA-VNS.
    """
    A, others = "PSOBV", [a for a in MAIN8 if a != "PSOBV"]
    pairs = [b for _, b in DECISION_PAIRS]
    out, lines = dict(pool=MAIN8, pairs=[f"{a}-{b}" for a, b in DECISION_PAIRS],
                      dloss_sign="dL = L(PSOBV) - L(other) in percentage points of the wake loss; negative = PSO-VNS better",
                      wtl_sign="W/T/L from the PSO-VNS side (W = PSO-VNS significantly better)"), []

    def block(X, label, cases_note):
        mets = [a for a in MAIN8 if a in set(X.Algorithm)]
        if A not in mets:
            return dict(available=False, note=f"no PSO-VNS runs ({label})")
        S = case_stats(X, mets)
        R = rank_matrix(S, mets)
        avg = R.mean(axis=0).to_dict() if len(R) else {}
        oth = [a for a in others if a in mets]
        CW = case_mean_wilcoxon(S, A, oth)
        crow = []
        for _, sub in X.groupby(CASE):
            crow += paired_vs(sub, A, oth)
        C = pd.DataFrame(crow)
        res = dict(available=True, n_cases=int(len(R)), methods=mets, n_comparisons=len(oth),
                   psobv_runs=int((X.Algorithm == A).sum()), psobv_seeds=sorted(map(int, X[X.Algorithm == A].Seed.unique())),
                   avg_rank=avg, rank_position={a: int(1 + sorted(avg.values()).index(avg[a])) for a in avg},
                   mean_loss_pct={a: float(S[(S.Algorithm == a) & S.Qualified & (S.Dataset != "HR")].Loss.mean()) for a in mets},
                   note=cases_note,
                   fallback_sources={a: sorted(set(X[(X.Algorithm == a) & X.Source.str.startswith("FALLBACK")].Source) |
                                               {d for d in X[(X.Algorithm == a) & X.Source.str.startswith("FALLBACK")].Dataset.map(lambda v: f"dataset {v}")})
                                     for a in mets if X[(X.Algorithm == a)].Source.str.startswith("FALLBACK").any()})
        if (X.Dataset == "HR").any():
            h = S[S.Dataset == "HR"].set_index("Algorithm")
            res["hr16_mean_aep"] = {a: float(h.Mean.get(a, np.nan)) for a in mets}
            res["hr16_feasible"] = {a: f"{int(h.NFeas.get(a, 0))}/{int(h.N.get(a, 0))}" for a in mets}
        for b in pairs:
            if b not in mets:
                res[b] = None; continue
            cw = CW.get(b, {})
            wtl = C[C.Baseline == b].Outcome.value_counts() if len(C) else pd.Series(dtype=int)
            res[b] = dict(case_mean_p=cw.get("p"), case_mean_p_holm7=cw.get("p_holm"), case_mean_rb=cw.get("rb"),
                          mean_dloss_pp=cw.get("mean_dloss_pp"), lower_loss_cases=f"{cw.get('focus_lower_loss_cases')}"
                          f" vs {cw.get('other_lower_loss_cases')}", n_cases=cw.get("n_cases"),
                          run_level_wtl=wtl_str(wtl), median_rb_run_level=float(C[C.Baseline == b].RB.median()) if len(C) else None)
        return res

    R6 = ALL[(ALL.Budget == 6030) & (ALL.Init == "random")]
    out["grid68"] = block(R6[R6.Dataset.isin(["1", "2"])], "68 cases", "68 benchmark cases, 6,030 calls")
    cases7 = LARGE + [HR16]
    out["large7"] = {}
    for b in BUDGETS:
        X = ALL[(ALL.Init == "random") & (ALL.Budget == b)]
        X = X[[(d, r, n) in cases7 for d, r, n in zip(X.Dataset, X.Radius, X.Turbines)]]
        out["large7"][b] = block(X, f"{b} calls", "six largest benchmark cases + Horns Rev 16" +
                                 (" (HR16: 10 seeds)" if b == 120030 else "")) if len(X) else dict(available=False, note="no data")
    out["iea37"] = {}
    IE = ALL[ALL.Dataset.str.startswith("IEA37")]
    for (n, b), y in IE.groupby(["Turbines", "Budget"]):
        mets = [a for a in MAIN8 if a in set(y.Algorithm)]
        key = f"{n}T_{b}"
        if A not in mets:
            out["iea37"][key] = dict(available=False, note="no PSO-VNS runs (iea%dp missing)" % n); continue
        Sy = case_stats(y, mets).set_index("Algorithm")
        tt = {t["Baseline"]: t for t in paired_vs(y, A, [a for a in others if a in mets])}
        out["iea37"][key] = dict(available=True, turbines=int(n), budget=int(b),
                                 mean_aep={a: float(Sy.Mean[a]) for a in mets}, best_aep={a: float(Sy.Best[a]) for a in mets},
                                 feasible={a: f"{int(Sy.NFeas[a])}/{int(Sy.N[a])}" for a in mets},
                                 rank={a: float(Sy.Rank[a]) for a in mets},
                                 **{b2: (dict(p_holm7=tt[b2]["PHolm"], rb=tt[b2]["RB"], outcome=tt[b2]["Outcome"],
                                              mean_aep_diff=float(Sy.Mean[A] - Sy.Mean[b2]))
                                         if b2 in tt else None) for b2 in pairs})

    # ---- printable block
    f3 = lambda v: "--" if v is None or not np.isfinite(v) else f"{v:.3f}"
    fp = lambda v: "--" if v is None or not np.isfinite(v) else (f"{v:.1e}" if v < 1e-3 else f"{v:.3f}")
    lines.append("DECISION BLOCK: PSO-VNS (PSOBV) vs PSO (PSOC) and vs SSA-VNS (SSABV); pool = " + ", ".join(LAB[a] for a in MAIN8))
    lines.append("  dL = L(PSO-VNS) - L(other) [pp], negative = PSO-VNS better; W/T/L from the PSO-VNS side (Holm over 7 per case)")

    def show(tag, r):
        if not r or not r.get("available"):
            lines.append(f"  {tag:26s} n/a ({(r or {}).get('note', 'no data')})"); return
        ar = r["avg_rank"]
        lines.append(f"  {tag:26s} cases={r['n_cases']:2d}  PSO-VNS runs={r['psobv_runs']} (seeds {len(r['psobv_seeds'])})  avg rank: "
                     + ", ".join(f"{LAB[a]} {ar[a]:.2f}" for a in ("PSOBV", "PSOC", "SSABV") if a in ar)
                     + f"  (best: {LAB[min(ar, key=ar.get)]} {min(ar.values()):.2f})")
        for b in pairs:
            x = r.get(b)
            if not x:
                lines.append(f"      vs {LAB[b]:8s} n/a"); continue
            lines.append(f"      vs {LAB[b]:8s} case-mean Wilcoxon p={fp(x['case_mean_p'])} (Holm7 {fp(x['case_mean_p_holm7'])}) "
                         f"r_rb={x['case_mean_rb']:+.2f}  mean dL={f3(x['mean_dloss_pp'])} pp  lower-loss cases {x['lower_loss_cases']}  "
                         f"run-level W/T/L {x['run_level_wtl']}")
        for a, v in r.get("fallback_sources", {}).items():
            lines.append(f"      NOTE: {LAB[a]} uses development-fallback runs here ({'; '.join(v)})")
        if "hr16_mean_aep" in r:
            h = r["hr16_mean_aep"]
            lines.append("      HR16 mean AEP (GWh/yr): " + ", ".join(f"{LAB[a]} {h[a]:.2f} ({r['hr16_feasible'][a]})"
                                                             for a in ("PSOBV", "PSOC", "SSABV") if a in h))
    show("68 cases, 6,030 calls", out["grid68"])
    for b in BUDGETS:
        show(f"large-6 + HR16, {b:,}", out["large7"][b])
    for k, r in out["iea37"].items():
        if not r.get("available"):
            lines.append(f"  IEA37 {k:20s} n/a ({r['note']})"); continue
        m = r["mean_aep"]
        s = f"  IEA37 {r['turbines']}T {r['budget']:>7,}   mean AEP (MWh): " + ", ".join(
            f"{LAB[a]} {m[a]:,.0f} [rank {r['rank'][a]:g}, {r['feasible'][a]}]" for a in ("PSOBV", "PSOC", "SSABV") if a in m)
        for b in pairs:
            x = r.get(b)
            if x:
                s += f"; vs {LAB[b]}: dAEP {x['mean_aep_diff']:+,.0f}, p_Holm7 {fp(x['p_holm7'])}, {x['outcome']}"
        lines.append(s)
    return out, lines


def robustness(G, methods, args, summary, tabs, OUT, only=False):
    log("\n[8] Robustness re-evaluation")
    import mpce_robustness as mr
    F = G[G.Feasible & G.Algorithm.isin(methods)]
    nruns = G[G.Algorithm.isin(methods)].groupby(CASE + ["Algorithm"]).size()
    RV = mr.reevaluate(F, HERE, args.procs, log)
    tex, out = mr.robust_table(RV, methods, LAB, rank_rule, n_runs=nruns)
    tabs["robust"] = tex
    summary["robustness"] = out
    RV.drop(columns=["Coordinates", "Curve", "Key"]).to_csv(OUT("mpce_reevaluation.csv"), index=False)
    for k, v in out.items():
        if k == "rounding_check":
            log(f"  re-evaluation of the rounded coordinates vs recorded objective: {v}"); continue
        log(f"  {k:12s} tau={v['tau']:.2f} same best={v['same_best_pct']:.0f}%  {LAB[FOCUS]} rank {v['avg_rank'].get(FOCUS, float('nan')):.2f}")
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

def evals(t):
    """Terminology of the paper in generated tables: "evaluations" for objective calls (R6-12, R6-32)."""
    t = re.sub(r"\b[Oo]bjective(?:-function)? calls\b", "evaluations", t)
    t = re.sub(r"\bper (?:objective )?call\b", "per evaluation", t)
    t = re.sub(r"\bCalls\b", "Evaluations", t)
    return re.sub(r"\bcalls\b", "evaluations", t)


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

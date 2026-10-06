"""Inference-robustness analyses for the MPCE resubmission (Phase 6, workstream W2; answers R3 items 3-7).

Usage (from the repository root or analysis/):  python3 analysis/mpce_inference_extra.py [--data-dir D] [--perm 20000]

Question: do the paper's statistical conclusions depend on analysis choices? Five re-analyses of the SAME
per-run data as mpce_results.py (68 cases = data sets I/II x radii 500/750/1000 m x N = 2..10/12/15; 6,030
evaluations, random initialization, 30 seed-paired runs):

  1. Qualification threshold (survivorship rule). The paper ranks a method in a case by its mean feasible
     objective only if >= 15 of its 30 runs are feasible (otherwise below all qualifying methods, by number of
     feasible runs), and the case-mean Wilcoxon imputes a maximal difference when only one of two methods
     qualifies. Everything is recomputed for thresholds 1, 10, 15 (baseline), 20, 25, 30: the 8-method Friedman
     ranking (tab:friedman68), the best-ranked method, the post hoc set, PSO-VNS vs PSO on the case means and the
     component-analysis contrasts.
  2. Case dependence. The 68 cases form six (data set, radius) clusters with nested N (identified from the data).
     (a) cluster means of the case-mean differences, exact sign test, exact Wilcoxon and an exact cluster
     sign-flip test over the 6 clusters, and a cluster-robust t test (CR1, G-1 d.f.); (b) leave-one-cluster-out
     (LOCO) re-runs of the case-mean Wilcoxon and of the Friedman ranking; (c) cluster (block) bootstrap CIs of the
     mean difference (10,000 resamples of the clusters, fixed seed).
  3. Multiplicity. One Holm correction (and Benjamini-Hochberg as a secondary view) over ALL case-mean Wilcoxon
     tests quoted in the main text (main-table pairs, component-analysis contrasts, N >= 10 subgroups, budget-
     split tests).
  4. Post hoc N >= 10 subgroup. Permutation tests of the interaction (mean PSO-VNS - PSO difference for N >= 10
     minus that for N < 10; unrestricted and stratified by cluster), Spearman correlation of the case difference
     with N (overall, per data set, and a within-cluster stratified trend test), and a paired Data Set II vs I
     test over identical (r, N) cases.
  Phase 6 controls: the four component-analysis contrasts with the disc-sampling control RSD-VNS (experiment rsdisc)
     and the omega = 0.9 budget-split tests (experiment omega90, PSOBV90) enter analyses 1-3 (the split tests, 12 cases
     = 2 per cluster, get their own threshold / cluster / LOCO results: summary keys split_case_level, cluster_split,
     loco_split).
  5. Effect sizes in energy terms: headline differences as % of the benchmark AEP and in MWh/yr of benchmark
     (model) energy. The benchmark objective is 15 x the expected farm power in kW (03_model.tex), so the gross
     benchmark AEP is objective / 15 x 8.76 MWh/yr.
  6. Wake-model uncertainty (to put the +-0.05 pp equivalence margin in context): |L_Gauss - L_Jensen| of the same
     feasible final layouts from mpce_reevaluation.csv (written by mpce_results.py, section 8), and the change of the
     PSO-VNS - PSO case-mean difference between the two wake models (\\NXModelShift... macros, tab:X-modelshift).

Review round 2 (R2 statistics review, lead decision D13; sections 7-11 of this file):
  7. Practical equivalence (margin EQ_MARGIN, read from mpce_summary.json) of the 18 pairs of tab:equivalence at three
     levels of inference: (a) seed level = fixed benchmark (runs resampled within cases, seed-paired), (b) case level
     (bootstrap over cases; reproduces mpce_results.equivalence_block exactly, check X38), (c) cluster level (CR2 variance
     with t(5), Bell-McCaffrey d.f. for reference; restricted wild-cluster bootstrap-t TOST with all 6^6 Webb draws).
     Bayesian signed-rank posterior means of theta (copy of mpce_results.bayes_signrank). Figure equiv_curve.pdf.
  8. Heterogeneity of PSO-VNS - PSO: cases beyond +-margin, equivalence for N < 10 / N >= 10 and for trivial /
     non-trivial wake loss (PSO-VNS case mean < / >= 0.2 pp), cluster means, data set x density interaction
     (density phi = N (l_min / 2r)^2), the margin in benchmark energy per turbine.
  9. Hodges-Lehmann estimates with exact Wilcoxon-inversion 95 % CIs.  10. Exact McNemar test of Horns Rev feasibility.
Tables tab:X-equiv-levels, tab:X-bayes, tab:X-heterogeneity, tab:X-beyond, tab:X-hl, tab:X-mcnemar; checks X38...

The data are loaded and the helpers (Wilcoxon with zero differences dropped, Holm, the ranking rule, the case-
mean Wilcoxon with imputation) are re-implemented here exactly as in mpce_results.py (that file is owned by
another workstream and is not imported); check X01 verifies that the baseline threshold reproduces
mpce_summary.json.

Outputs (analysis/): mpce_summary_extra.json, mpce_numbers_extra.tex (\\NX... macros), mpce_supp_inference.tex
(supplementary tables, labels tab:X-...). Then run mpce_check_extra.py (checks X01...).
"""
import os, re, sys, glob, json, math, time, argparse, itertools
import numpy as np, pandas as pd
from scipy.stats import wilcoxon, friedmanchisquare, rankdata, norm, binomtest, spearmanr, t as t_dist

HERE = os.path.dirname(os.path.abspath(__file__))
CASE = ["Dataset", "Radius", "Turbines"]
KEY = ["Algorithm", "Dataset", "Radius", "Turbines", "Seed", "Budget", "Init"]
MAIN8 = ["PSOBV", "PSOC", "SSABV", "SSA", "LXSSA", "DE", "BVNS", "SLSQP"]
ABL = ["PSOBV", "SSABV", "LXBV", "RSVNS", "RSDVNS", "BVNS", "PSOC", "SSA", "LXSSA"]   # RSDVNS: disc-sampling control (rsdisc)
FOCUS = "PSOBV"
SPLITCASES = [(ds, r, n) for ds in ("1", "2") for r, n in ((500, 6), (750, 8), (1000, 10), (500, 10), (750, 12), (1000, 15))]
THRESHOLDS = [1, 10, 15, 20, 25, 30]
BASE_THR = 15
BOOT_N, BOOT_SEED = 10000, 20260928
PERM_SEED = 20260929
ENERGY_FACTOR = 8.76 / 15.0          # MWh/yr per unit of benchmark objective (objective = 15 x expected power in kW)
LAB = {"PSOBV": "PSO-VNS", "PSOC": "PSO", "SSABV": "SSA-VNS", "SSA": "SSA", "LXSSA": "LX-SSA", "DE": "DE",
       "BVNS": "VNS", "SLSQP": "MS-SLSQP", "LXBV": "LX-SSA-VNS", "RSVNS": "RS-VNS",
       "PSOBV25": "PSO-VNS ($\\omega=0.25$)", "PSOBV75": "PSO-VNS ($\\omega=0.75$)", "PSOBV90": "PSO-VNS ($\\omega=0.9$)",
       "RSDVNS": "RSD-VNS"}
MAC = {"PSOBV": "PSOVNS", "PSOC": "PSO", "SSABV": "SSAVNS", "SSA": "SSA", "LXSSA": "LXSSA", "DE": "DE",
       "BVNS": "VNS", "SLSQP": "MSSLSQP", "LXBV": "LXSSAVNS", "RSVNS": "RSVNS", "RSDVNS": "RSDVNS"}
# main-table pairs (PSO-VNS vs each method) and the fifteen component-analysis contrasts (first vs second; the last
# four, with the disc-sampling control RSD-VNS, were added in Phase 6 -- same family as mpce_results.py)
MAIN_PAIRS = [(FOCUS, b) for b in MAIN8 if b != FOCUS]
ABL_CONTR = [("PSOBV", "PSOC"), ("SSABV", "SSA"), ("LXBV", "LXSSA"), ("PSOBV", "BVNS"), ("PSOBV", "RSVNS"),
             ("SSABV", "RSVNS"), ("LXBV", "RSVNS"), ("PSOBV", "SSABV"), ("PSOBV", "LXBV"), ("SSABV", "LXBV"),
             ("LXSSA", "SSA"), ("SSABV", "RSDVNS"), ("LXBV", "RSDVNS"), ("PSOBV", "RSDVNS"), ("RSDVNS", "RSVNS")]
RSD_CONTR = [("SSABV", "RSDVNS"), ("LXBV", "RSDVNS"), ("PSOBV", "RSDVNS"), ("RSDVNS", "RSVNS")]
ALL_PAIRS = list(dict.fromkeys(MAIN_PAIRS + ABL_CONTR))          # 19 unique pairs (15 before Phase 6 controls)
# the conclusions of the paper that the threshold / LOCO analyses re-test (component analysis: D1 of PHASE4.md)
KEY_CONTR = [("PSOBV", "PSOC"), ("SSABV", "RSVNS"), ("LXBV", "RSVNS"), ("SSABV", "LXBV"), ("PSOBV", "RSVNS")] + RSD_CONTR
EXPECT = {("PSOBV", "PSOC"): "ns", ("SSABV", "RSVNS"): "A", ("LXBV", "RSVNS"): "ns", ("SSABV", "LXBV"): "A",
          ("PSOBV", "RSVNS"): "A",
          # Phase 6 controls (lead decision D11): SSA-VNS and LX-SSA-VNS significantly WORSE than RSD-VNS (B),
          # PSO-VNS better than RSD-VNS, RSD-VNS better than RS-VNS
          ("SSABV", "RSDVNS"): "B", ("LXBV", "RSDVNS"): "B", ("PSOBV", "RSDVNS"): "A", ("RSDVNS", "RSVNS"): "A"}
# A = first method significantly better (p < 0.05, lower mean loss); B = second significantly better; ns = p >= 0.05
# budget split (12 cases, 2 per cluster): the case-mean tests of the split table; omega = 0.9 (PSOBV90, experiment
# omega90) added in Phase 6. Split tests are unadjusted in the paper; EXPECT_SPLIT = the paper's verdicts (raw p).
SPM = ["PSOBV25", "PSOBV", "PSOBV75", "PSOBV90", "PSOC"]
SPLIT_PAIRS = [("PSOBV", "PSOBV25"), ("PSOBV", "PSOBV75"), ("PSOBV", "PSOC"), ("PSOBV75", "PSOC"),
               ("PSOBV", "PSOBV90"), ("PSOBV90", "PSOBV75")]
SPLIT_NAME = {("PSOBV", "PSOBV25"): "SplitFiftyVsTwentyFive", ("PSOBV", "PSOBV75"): "SplitFiftyVsSeventyFive",
              ("PSOBV", "PSOC"): "SplitFiftyVsHundred", ("PSOBV75", "PSOC"): "SplitSeventyFiveVsHundred",
              ("PSOBV", "PSOBV90"): "SplitFiftyVsNinety", ("PSOBV90", "PSOBV75"): "SplitNinetyVsSeventyFive"}
SPLIT_LAB = {"PSOBV25": "$\\omega=0.25$", "PSOBV": "$\\omega=0.5$", "PSOBV75": "$\\omega=0.75$", "PSOBV90": "$\\omega=0.9$",
             "PSOC": "$\\omega=1$ (PSO)"}
EXPECT_SPLIT = {("PSOBV", "PSOBV25"): "A", ("PSOBV", "PSOBV75"): "B", ("PSOBV", "PSOC"): "ns", ("PSOBV75", "PSOC"): "A",
                ("PSOBV", "PSOBV90"): "B", ("PSOBV90", "PSOBV75"): "ns"}
def _eq_margin():
    """The practical-equivalence margin (pp of wake loss) is owned by mpce_results.py (EQ_MARGIN); it is read from
    mpce_summary.json ("equivalence" -> "margin_pp") so that both scripts always use the same value."""
    try:
        return float(json.load(open(os.path.join(HERE, "mpce_summary.json")))["equivalence"]["margin_pp"]), "mpce_summary.json"
    except Exception:                                                   # pragma: no cover
        return 0.05, "default (mpce_summary.json missing)"


EQ_MARGIN, EQ_MARGIN_SRC = _eq_margin()
# ---- Phase 6 / review round 2 (R2 statistics, lead decision D13): equivalence at three inference levels,
# heterogeneity, Bayesian posterior means, Hodges-Lehmann estimates, McNemar (functions in section 7)
# the 18 pairs of tab:equivalence (mpce_results.equivalence_block: component-analysis pairs, then PSO-VNS vs each
# main method); order as in that table
EQ_PAIRS = ["SSABV-RSVNS", "LXBV-RSVNS", "PSOBV-RSVNS", "SSABV-LXBV", "LXSSA-SSA", "SSABV-SSA", "LXBV-LXSSA",
            "PSOBV-PSOC", "PSOBV-BVNS", "SSABV-RSDVNS", "LXBV-RSDVNS", "PSOBV-RSDVNS", "RSDVNS-RSVNS",
            "PSOBV-SSABV", "PSOBV-SSA", "PSOBV-LXSSA", "PSOBV-DE", "PSOBV-SLSQP"]
HL_KEY = [("PSOBV", "PSOC"), ("SSABV", "RSDVNS"), ("LXBV", "RSDVNS"), ("PSOBV", "RSDVNS"), ("SSABV", "RSVNS"), ("SSABV", "LXBV")]
SEED_BOOT_SEED = 20260930            # seed-level (fixed-benchmark) bootstrap: runs resampled within cases
BAYES_S, BAYES_Z0, BAYES_N, BAYES_SEED = 0.5, 0.0, 50000, 20260928   # as mpce_results.py (reproduced, check X38)
WEBB = np.array([-math.sqrt(1.5), -1.0, -math.sqrt(0.5), math.sqrt(0.5), 1.0, math.sqrt(1.5)])  # Webb (2014) 6-point
SMIN_M = 308.0                       # minimum spacing 4D = 308 m (mpce_results.SMIN_M); density phi = N (SMIN_M / 2r)^2
CURVE_GRID = np.round(np.arange(0.0, 0.1201, 0.001), 3)               # margins (pp) of the equivalence curve
TRIVIAL_LOSS_PP = 0.2                # lead request: "trivial" case = case-mean wake loss of PSO-VNS below 0.2 pp


def pk(a, b):
    return f"{a}-{b}"


def plab(a, b):
    return f"{LAB[a]} vs.\\ {LAB[b]}"


def pmac(a, b):
    return f"{MAC[a]}vs{MAC[b]}"


# ------------------------------------------------------------------ data (as mpce_results.load, 68 cases only)
def std_cols(df, budget=6030, init="random"):
    df = df.copy()
    if "Budget" not in df:
        df["Budget"] = budget
    if "Init" not in df:
        df["Init"] = init
    df["Dataset"] = df.Dataset.astype(str).str.replace(r"\.0$", "", regex=True)
    df["Radius"] = pd.to_numeric(df.Radius, errors="coerce").fillna(0).astype(int)
    for c in ("Turbines", "Seed", "Budget"):
        df[c] = df[c].astype(int)
    if "WakeLoss" not in df:
        df["WakeLoss"] = df.Ideal - df.Objective
    df["Feasible"] = df.Feasible.astype(str).str.lower().isin(["true", "1", "1.0"])
    df["LossPct"] = 100 * df.WakeLoss / df.Ideal
    return df


def read_shards(exp, data_dir):
    files = glob.glob(os.path.join(data_dir, f"mpce_{exp}_s*of*.csv"))
    by_k = {}
    for fn in files:
        m = re.search(rf"mpce_{exp}_s(\d+)of(\d+)\.csv$", os.path.basename(fn))
        if m:
            by_k.setdefault(int(m.group(2)), {})[int(m.group(1))] = fn
    complete = {k: v for k, v in by_k.items() if set(v) == set(range(k))}
    if not complete:
        raise SystemExit(f"experiment {exp}: no complete shard set in {data_dir}")
    k = max(complete, key=lambda kk: sum(os.path.getsize(f) for f in complete[kk].values()))
    df = pd.concat([pd.read_csv(f, usecols=lambda c: c not in ("Coordinates", "Curve")) for f in sorted(complete[k].values())],
                   ignore_index=True)
    return std_cols(df).drop_duplicates(KEY, keep="last")


def load(data_dir):
    parts = []
    for fn, drop in (("fresh_grid.csv", ["VNS", "PSO", "SLSQP"]), ("fresh_vgrid.csv", []), ("fresh_bgrid.csv", [])):
        d = std_cols(pd.read_csv(os.path.join(HERE, fn), usecols=lambda c: c not in ("Coordinates", "Curve")))
        parts.append(d[~d.Algorithm.isin(drop)])           # old PSO dropped; SLSQP replaced by mpce_slsqp (as load())
    for exp in ("rsvns", "psoc", "psobv", "slsqp", "psosplit", "omega90", "rsdisc"):   # omega90 / rsdisc: Phase 6 controls
        parts.append(read_shards(exp, data_dir))
    A = pd.concat(parts, ignore_index=True).drop_duplicates(KEY, keep="first")
    A = A[(A.Budget == 6030) & (A.Init == "random") & A.Dataset.isin(["1", "2"])].reset_index(drop=True)
    return A


# ------------------------------------------------------------------ statistics (as in mpce_results.py)
def holm(p):
    p = np.asarray(p, float); o = np.argsort(p); m = len(p); out = np.empty(m); run = 0
    for i, k in enumerate(o):
        run = max(run, min(1.0, (m - i) * p[k])); out[k] = run
    return out


def bh(p):
    p = np.asarray(p, float); m = len(p); o = np.argsort(p); q = np.empty(m); run = 1.0
    for i in range(m - 1, -1, -1):
        k = o[i]; run = min(run, p[k] * m / (i + 1)); q[k] = run
    return q


def wil(d):
    d = np.asarray(d, float); d = d[np.isfinite(d)]
    nz = d[np.abs(d) > 1e-9]
    if len(nz) == 0:
        return 1.0, 0.0, 0
    p = wilcoxon(nz).pvalue
    rk = rankdata(np.abs(nz))
    return float(p), float((rk[nz > 0].sum() - rk[nz < 0].sum()) / rk.sum()), len(nz)


def case_stats(D, methods, thr):
    """Per case x method: runs, feasible runs, mean objective / loss of the feasible runs, Qualified (>= thr)."""
    D = D[D.Algorithm.isin(methods)]
    g = D.groupby(CASE + ["Algorithm"])
    S = pd.DataFrame({"N": g.size(), "NFeas": g.Feasible.sum()})
    gf = D[D.Feasible].groupby(CASE + ["Algorithm"])
    S["Mean"] = gf.Objective.mean(); S["Loss"] = gf.LossPct.mean(); S["Ideal"] = g.Ideal.first()
    S = S.reset_index()
    S["NFeas"] = S.NFeas.astype(int)
    assert (S.N == 30).all(), "every case/method must have 30 runs"
    S["Qualified"] = (S.NFeas >= thr) & S.Mean.notna()
    return S


def ranks_rule(fmean, nfeas, qual):
    fmean = np.round(np.asarray(fmean, float), 6); nfeas = np.asarray(nfeas, float); q = np.asarray(qual, bool)
    r = np.empty(len(fmean))
    if q.any():
        r[q] = rankdata(-fmean[q])
    if (~q).any():
        r[~q] = q.sum() + rankdata(-nfeas[~q])
    return r


def friedman(S, methods, focus=FOCUS):
    S = S[S.Algorithm.isin(methods)]
    rows = []
    for c, s in S.groupby(CASE):
        s = s.set_index("Algorithm").reindex(methods)
        rows.append(ranks_rule(s.Mean.values, s.NFeas.values, s.Qualified.values))
    X = np.array(rows); n, k = X.shape
    chi, pf = friedmanchisquare(*[X[:, j] for j in range(k)])
    avg = X.mean(0); se = np.sqrt(k * (k + 1) / (6 * n)); fi = methods.index(focus)
    zf = {a: (avg[j] - avg[fi]) / se for j, a in enumerate(methods) if j != fi}
    pz = dict(zip(zf, holm([2 * norm.sf(abs(z)) for z in zf.values()])))
    order = sorted(methods, key=lambda a: avg[methods.index(a)])
    return dict(n_cases=n, chi2=float(chi), p=float(pf), avg_rank={a: float(avg[j]) for j, a in enumerate(methods)},
                best_ranked=methods[int(np.argmin(avg))], order=order,
                p_holm_vs_focus={a: float(v) for a, v in pz.items()},
                posthoc_nonsig=sorted(a for a, v in pz.items() if v >= 0.05))


def cm_test(S, a, b, val="Loss"):
    """Case-mean Wilcoxon of a vs b exactly as mpce_results.case_mean_wilcoxon (imputation of maximal differences
    for cases in which only one qualifies; cases in which neither qualifies dropped). d = L(a) - L(b)."""
    P = S.pivot_table(index=CASE, columns="Algorithm", values=val)
    Q = S.pivot_table(index=CASE, columns="Algorithm", values="Qualified").astype(float)
    qa, qb = Q[a] > 0.5, Q[b] > 0.5
    both = qa & qb
    d = (P[b] - P[a])[both].values                       # > 0: a better (as in case_mean_wilcoxon)
    big = 10 * (np.nanmax(np.abs(d)) if len(d) and np.isfinite(d).any() else 1.0) + 1.0
    dfull = np.r_[d, np.full(int((qa & ~qb).sum()), big), np.full(int((~qa & qb).sum()), -big)]
    p, rb, nn = wil(dfull)
    return dict(p=p, rb=rb, n_cases=int(len(dfull)), wins=int((dfull > 1e-9).sum()), losses=int((dfull < -1e-9).sum()),
                mean_dloss_pp=float(np.mean(-d)) if len(d) else float("nan"), n_both=int(both.sum()),
                imputed=int((qa & ~qb).sum() + (~qa & qb).sum()), dropped=int((~qa & ~qb).sum()))


def case_diffs(S, a, b, val="Loss"):
    """Series (index = case) of val(a) - val(b) over the cases in which both qualify."""
    P = S.pivot_table(index=CASE, columns="Algorithm", values=val)
    Q = S.pivot_table(index=CASE, columns="Algorithm", values="Qualified").astype(bool)
    both = Q[a] & Q[b]
    return (P[a] - P[b])[both]


def verdict(x):
    """'A' = first significantly better (p < 0.05 and lower mean loss), 'B' = second significantly better, 'ns'."""
    if x["p"] >= 0.05:
        return "ns"
    return "A" if x["wins"] > x["losses"] else "B"


# ------------------------------------------------------------------ formatting
def num(v, d=2, sign=False):
    if v is None or (isinstance(v, float) and not math.isfinite(v)):
        raise ValueError("missing value")
    s = f"{v:+,.{d}f}" if sign else f"{v:,.{d}f}"
    s = s.replace(",", "{,}")
    return s.replace("-", "\\ensuremath{-}", 1) if s.startswith("-") else s


def pval(p):
    if p < 1e-3:
        m, e = f"{p:.1e}".split("e")
        return f"\\ensuremath{{{m}\\times10^{{{int(e)}}}}}"
    if p >= 0.995:
        return "1.0"
    return f"{p:.2g}" if p < 0.1 else f"{p:.2f}"


def pval_up(p):
    if p >= 0.995:
        return "1.0"
    e = math.floor(math.log10(p))
    if p < 1e-3:
        m = math.ceil(p / 10 ** e * 10 - 1e-9) / 10
        if m >= 10:
            m, e = 1.0, e + 1
        return f"\\ensuremath{{{m:.1f}\\times10^{{{e}}}}}"
    step = 10 ** (e - 1) if p < 0.1 else 0.01
    v = math.ceil(p / step - 1e-9) * step
    d = max(0, -(e - 1)) if p < 0.1 else 2
    return f"{v:.{d}f}"


def pval_down(p):
    """p rounded DOWN (for 'p >= x' statements)."""
    if p >= 0.995:
        return "0.99"
    e = math.floor(math.log10(p))
    if p < 1e-3:
        m = math.floor(p / 10 ** e * 10 + 1e-9) / 10
        return f"\\ensuremath{{{m:.1f}\\times10^{{{e}}}}}"
    step = 10 ** (e - 1) if p < 0.1 else 0.01
    v = math.floor(p / step + 1e-9) * step
    d = max(0, -(e - 1)) if p < 0.1 else 2
    return f"{v:.{d}f}"


def tp(p):
    """p for tables (plain math)."""
    if not np.isfinite(p):
        return "--"
    if p < 1e-3:
        m, e = f"{p:.1e}".split("e")
        return f"${m}\\times10^{{{int(e)}}}$"
    return f"{p:.3f}"


def ci_txt(ci, d=3):
    return "[" + ", ".join(num(v, d) for v in ci) + "]"


def ci_tab(ci, d=3):
    return "[" + ", ".join(f"{v:.{d}f}".replace("-", "$-$") for v in ci) + "]"


def stack(a, b):
    return f"\\begin{{tabular}}{{@{{}}c@{{}}}}{a}\\\\{b}\\end{{tabular}}"


def listing(items):
    items = list(items)
    if len(items) <= 1:
        return "".join(items)
    return ", ".join(items[:-1]) + " and " + items[-1]


WORD = {0: "none", 1: "one", 2: "two", 3: "three", 4: "four", 5: "five", 6: "six", 7: "seven", 8: "eight", 9: "nine",
        10: "ten", 11: "eleven", 12: "twelve", 13: "thirteen", 14: "fourteen", 15: "fifteen", 16: "sixteen",
        17: "seventeen", 18: "eighteen", 19: "nineteen", 20: "twenty"}
WORD.update({20 + i: "twenty-" + WORD[i] for i in range(1, 10)})
WORD[30] = "thirty"


def table(env, caption, label, spec, header, lines, sep="2.5pt", pos="!htb", foot=None, resize=False):
    body = "\n".join(lines)
    ft = "" if not foot else "\n" + "\n".join(foot)
    tab = (f"\\begin{{tabular}}{{{spec}}}\n\\toprule\n{header} \\\\\n"
           f"\\midrule\n{body}\n\\bottomrule{ft}\n\\end{{tabular}}")
    if resize:
        tab = "\\resizebox{\\textwidth}{!}{%\n" + tab + "}"
    return (f"\\begin{{{env}}}[{pos}]\n\\centering\n\\caption{{{caption}}}\n\\label{{{label}}}\n"
            f"\\scriptsize\\setlength{{\\tabcolsep}}{{{sep}}}\n{tab}\n\\end{{{env}}}\n")


# ------------------------------------------------------------------ 1. qualification threshold
def threshold_block(G, Sp):
    """G: runs of the 68 cases; Sp: runs of the 12 budget-split cases (SPM settings)."""
    out = {}
    for thr in THRESHOLDS:
        S = case_stats(G, sorted(set(MAIN8) | set(ABL)), thr)
        SPt = case_stats(Sp, SPM, thr)
        FR = friedman(S, MAIN8)
        FA = friedman(S, ABL)
        main = {b: cm_test(S, FOCUS, b) for b in MAIN8 if b != FOCUS}
        for b, h in zip(main, holm([main[b]["p"] for b in main])):
            main[b]["p_holm"] = float(h)
        abl = {pk(a, b): cm_test(S, a, b) for a, b in ABL_CONTR}
        for k, h in zip(abl, holm([abl[k]["p"] for k in abl])):
            abl[k]["p_holm"] = float(h)
        concl = {
            "best_ranked_PSOBV": FR["best_ranked"] == FOCUS,
            "posthoc_only_PSOC_nonsig": FR["posthoc_nonsig"] == ["PSOC"],
            "PSOBV-PSOC_ns": main["PSOC"]["p_holm"] >= 0.05 and main["PSOC"]["p"] >= 0.05,
        }
        for a, b in KEY_CONTR[1:]:
            x = abl[pk(a, b)]
            ok_raw = verdict(x) == EXPECT[(a, b)]
            ok_holm = verdict(dict(x, p=x["p_holm"])) == EXPECT[(a, b)]
            concl[f"{pk(a, b)}_{EXPECT[(a, b)]}"] = bool(ok_raw and ok_holm)
        split = {pk(a, b): cm_test(SPt, a, b) for a, b in SPLIT_PAIRS}
        for a, b in SPLIT_PAIRS:                                        # unadjusted, as in the split table
            concl[f"split:{pk(a, b)}_{EXPECT_SPLIT[(a, b)]}"] = bool(verdict(split[pk(a, b)]) == EXPECT_SPLIT[(a, b)])
        out[thr] = dict(threshold=thr, friedman=FR, ablation_friedman=dict(avg_rank=FA["avg_rank"], best_ranked=FA["best_ranked"]),
                        main_case_mean=main, ablation_case_mean=abl, split_case_mean=split,
                        imputed_split=int(sum(split[k]["imputed"] for k in split)),
                        imputed_main=int(sum(main[b]["imputed"] for b in main)),
                        dropped_main=int(sum(main[b]["dropped"] for b in main)),
                        imputed_ablation=int(sum(abl[k]["imputed"] for k in abl)),
                        unqualified_cells_main=int((~S[S.Algorithm.isin(MAIN8)].Qualified).sum()),
                        conclusions=concl, all_unchanged=bool(all(concl.values())))
    return out


# ------------------------------------------------------------------ 2. clusters
def cluster_block(S, clusters, pairs=ALL_PAIRS):
    cl_of = {c: (c[0], c[1]) for c in S.set_index(CASE).index.unique()}
    res = {}
    G = len(clusters)
    signs = np.array(list(itertools.product([1, -1], repeat=G)))       # all 2^G sign patterns
    tcrit = t_dist.ppf(0.975, G - 1)
    for a, b in pairs:
        d = case_diffs(S, a, b)                                          # L(a) - L(b), both qualified
        cid = np.array([clusters.index(cl_of[c]) for c in d.index])
        dv = d.values
        n_g = np.array([(cid == g).sum() for g in range(G)]); s_g = np.array([dv[cid == g].sum() for g in range(G)])
        cmeans = s_g / np.where(n_g > 0, n_g, np.nan)
        fav_a = int((cmeans < -1e-9).sum()); fav_b = int((cmeans > 1e-9).sum())
        sign_p = float(binomtest(fav_a, fav_a + fav_b, 0.5).pvalue) if fav_a + fav_b else 1.0
        cm = cmeans[np.isfinite(cmeans)]
        wil_p = float(wilcoxon(cm, method="exact").pvalue) if (np.abs(cm) > 1e-12).any() else 1.0
        # exact cluster sign-flip test of the pooled (case-weighted) mean: flip whole clusters
        tot = np.abs(signs @ s_g); obs = abs(s_g.sum())
        flip_p = float(np.mean(tot >= obs - 1e-12))
        # cluster-robust (CR1) t test of the pooled mean, G-1 d.f.
        est = dv.mean(); e = dv - est
        sg_e = np.array([e[cid == g].sum() for g in range(G)])
        V = G / (G - 1) * (sg_e ** 2).sum() / len(dv) ** 2
        tstat = est / math.sqrt(V) if V > 0 else float("inf")
        cr_p = float(2 * t_dist.sf(abs(tstat), G - 1))
        cr_ci = [float(est - tcrit * math.sqrt(V)), float(est + tcrit * math.sqrt(V))]
        # cluster (block) bootstrap of the pooled mean
        rng = np.random.default_rng(BOOT_SEED)
        idx = rng.integers(0, G, (BOOT_N, G))
        bm = s_g[idx].sum(1) / n_g[idx].sum(1)
        cb_ci = [float(np.quantile(bm, 0.025)), float(np.quantile(bm, 0.975))]
        # case bootstrap (as mpce_results.boot_ci) for comparison
        rng = np.random.default_rng(BOOT_SEED)
        cbm = dv[rng.integers(0, len(dv), (BOOT_N, len(dv)))].mean(1)
        case_ci = [float(np.quantile(cbm, 0.025)), float(np.quantile(cbm, 0.975))]
        res[pk(a, b)] = dict(a=a, b=b, n_cases=int(len(dv)), mean_dloss_pp=float(est),
                             cluster_means_pp={f"{c[0]}-{c[1]}": (None if not np.isfinite(v) else float(v))
                                               for c, v in zip(clusters, cmeans)},
                             cluster_n={f"{c[0]}-{c[1]}": int(v) for c, v in zip(clusters, n_g)},
                             clusters_favour_first=fav_a, clusters_favour_second=fav_b, n_clusters=G,
                             sign_test_p=sign_p, wilcoxon_exact_p=wil_p, cluster_signflip_p=flip_p,
                             cr1_t=float(tstat), cr1_p=cr_p, cr1_ci95=cr_ci,
                             cluster_boot_ci95=cb_ci, case_boot_ci95=case_ci,
                             cluster_boot_excludes_zero=bool(cb_ci[0] > 0 or cb_ci[1] < 0),
                             cr1_excludes_zero=bool(cr_ci[0] > 0 or cr_ci[1] < 0))
    return res


def loco_block(S, clusters, full, pairs=ALL_PAIRS, rank=True):
    out = {}
    per = []
    for c in clusters:
        m = ~((S.Dataset == c[0]) & (S.Radius == c[1]))
        Sx = S[m]
        if rank:
            FR = friedman(Sx, MAIN8)
            row = dict(dropped=f"{c[0]}-{c[1]}", n_cases=FR["n_cases"], best_ranked=FR["best_ranked"],
                       posthoc_nonsig=FR["posthoc_nonsig"], avg_rank_PSOBV=FR["avg_rank"]["PSOBV"],
                       avg_rank_PSOC=FR["avg_rank"]["PSOC"], tests={})
        else:
            row = dict(dropped=f"{c[0]}-{c[1]}", n_cases=int(len(Sx[CASE].drop_duplicates())), tests={})
        for a, b in pairs:
            x = cm_test(Sx, a, b)
            row["tests"][pk(a, b)] = dict(p=x["p"], wins=x["wins"], losses=x["losses"], mean_dloss_pp=x["mean_dloss_pp"],
                                          verdict=verdict(x))
        per.append(row)
    for a, b in pairs:
        k = pk(a, b)
        ps = [r["tests"][k]["p"] for r in per]; mds = [r["tests"][k]["mean_dloss_pp"] for r in per]
        vs = [r["tests"][k]["verdict"] for r in per]
        fsign = np.sign(full[k]["mean_dloss_pp"])
        out[k] = dict(p_min=float(min(ps)), p_max=float(max(ps)), mean_dloss_min=float(min(mds)), mean_dloss_max=float(max(mds)),
                      full_verdict=full[k]["verdict"], verdicts=vs,
                      direction_flips=int(sum(np.sign(v) != fsign for v in mds)),
                      wins_losses_flips=int(sum((r["tests"][k]["wins"] > r["tests"][k]["losses"]) !=
                                                (full[k]["wins"] > full[k]["losses"]) for r in per)),
                      verdict_changes=int(sum(v != full[k]["verdict"] for v in vs)),
                      changed_when_dropping=[r["dropped"] for r, v in zip(per, vs) if v != full[k]["verdict"]])
    if not rank:
        return dict(per_cluster=per, pairs=out)
    return dict(per_cluster=per, pairs=out,
                best_ranked_always_PSOBV=bool(all(r["best_ranked"] == FOCUS for r in per)),
                posthoc_only_PSOC_always=bool(all(r["posthoc_nonsig"] == ["PSOC"] for r in per)))


# ------------------------------------------------------------------ 3. multiplicity
def multiplicity_block(S, G, SP):
    tests = []
    main = {b: cm_test(S, FOCUS, b) for b in MAIN8 if b != FOCUS}
    for b, x in main.items():
        tests.append(dict(key=f"main:{pk(FOCUS, b)}", label=plab(FOCUS, b), family="Main comparison (Table~\\ref{M-tab:friedman68})",
                          n=x["n_cases"], p=x["p"], wins=x["wins"], losses=x["losses"], dl=x["mean_dloss_pp"]))
    for a, b in ABL_CONTR:
        if (a, b) in MAIN_PAIRS:
            continue                                                    # same test as in the main table
        x = cm_test(S, a, b)
        tests.append(dict(key=f"abl:{pk(a, b)}", label=plab(a, b), family="Component analysis (Table~\\ref{M-tab:ablation}; pairs not in the main comparison)",
                          n=x["n_cases"], p=x["p"], wins=x["wins"], losses=x["losses"], dl=x["mean_dloss_pp"]))
    for tag, dsl, txt in (("Large", ("1", "2"), "$N\\ge10$, both data sets"), ("dsILarge", ("1",), "$N\\ge10$, Data Set I"),
                          ("dsIILarge", ("2",), "$N\\ge10$, Data Set II")):
        Sg = S[S.Dataset.isin(dsl) & (S.Turbines >= 10)]
        x = cm_test(Sg, FOCUS, "PSOC")
        tests.append(dict(key=f"sub:{tag}", label=f"PSO-VNS vs.\\ PSO, {txt}", family="Post hoc subgroups of PSO-VNS vs.\\ PSO",
                          n=x["n_cases"], p=x["p"], wins=x["wins"], losses=x["losses"], dl=x["mean_dloss_pp"]))
    for a, b in SPLIT_PAIRS:
        x = cm_test(SP, a, b)
        tests.append(dict(key=f"split:{pk(a, b)}", label=plab(a, b).replace("PSO-VNS vs", "PSO-VNS ($\\omega=0.5$) vs", 1)
                          if a == "PSOBV" else plab(a, b), family="Budget split (Table~\\ref{M-tab:split}, 12 cases)",
                          n=x["n_cases"], p=x["p"], wins=x["wins"], losses=x["losses"], dl=x["mean_dloss_pp"]))
    p = [t["p"] for t in tests]
    for t, h, q in zip(tests, holm(p), bh(p)):
        t["p_holm_all"] = float(h); t["q_bh_all"] = float(q)
        t["sig_raw"] = t["p"] < 0.05; t["sig_holm"] = t["p_holm_all"] < 0.05; t["sig_bh"] = t["q_bh_all"] < 0.05
    lost_holm = [t["key"] for t in tests if t["sig_raw"] and not t["sig_holm"]]
    lost_bh = [t["key"] for t in tests if t["sig_raw"] and not t["sig_bh"]]
    sig_h = [t["p_holm_all"] for t in tests if t["sig_holm"]]
    return dict(n_tests=len(tests), tests=tests, n_sig_raw=int(sum(t["sig_raw"] for t in tests)),
                n_sig_holm=int(sum(t["sig_holm"] for t in tests)), n_sig_bh=int(sum(t["sig_bh"] for t in tests)),
                lost_under_holm=lost_holm, lost_under_bh=lost_bh, max_sig_p_holm=float(max(sig_h)) if sig_h else None,
                note="one family = every case-mean Wilcoxon test quoted in the main text (main-table pairs, component-"
                     "analysis contrasts not already in the main table (incl. the four RSD-VNS contrasts), the three "
                     "N >= 10 subgroups, the six budget-split tests incl. omega = 0.9); Holm (FWER) and Benjamini-"
                     "Hochberg (FDR) over this family")


# ------------------------------------------------------------------ 4. post hoc subgroup
def perm_diff(d, big, strata, B, rng):
    """two-sided permutation p of mean(d[big]) - mean(d[~big]); labels permuted within strata."""
    obs = d[big].mean() - d[~big].mean()
    cnt = 0
    L = np.tile(big, (B, 1))
    for s in np.unique(strata):
        ix = np.where(strata == s)[0]
        if len(ix) > 1:
            L[:, ix] = rng.permuted(L[:, ix], axis=1)
    nb = L.sum(1)
    stat = (L * d).sum(1) / nb - ((~L) * d).sum(1) / (len(d) - nb)
    cnt = int((np.abs(stat) >= abs(obs) - 1e-12).sum())
    return float(obs), float((cnt + 1) / (B + 1))


def strat_spearman(d, n, strata, B, rng):
    """sum over strata of the within-stratum Spearman rho(N, d); permutation of d within strata."""
    us = [s for s in np.unique(strata) if (strata == s).sum() > 2]
    obs, rhos = 0.0, {}
    rank_n = {s: rankdata(n[strata == s]) for s in us}
    for s in us:
        r = spearmanr(n[strata == s], d[strata == s]).correlation
        rhos[s] = float(r); obs += r
    stat = np.zeros(B)
    for s in us:
        ds = d[strata == s]
        P = rng.permuted(np.tile(ds, (B, 1)), axis=1)
        rp = rankdata(P, axis=1)
        rn = rank_n[s]
        rn_c = rn - rn.mean(); rp_c = rp - rp.mean(1, keepdims=True)
        stat += (rp_c @ rn_c) / np.sqrt((rp_c ** 2).sum(1) * (rn_c ** 2).sum())
    p = float(((np.abs(stat) >= abs(obs) - 1e-12).sum() + 1) / (B + 1))
    return float(obs / len(us)), p, rhos


def subgroup_block(S, clusters, B):
    rng = np.random.default_rng(PERM_SEED)
    d = case_diffs(S, FOCUS, "PSOC")                                     # L(PSO-VNS) - L(PSO); negative = PSO-VNS better
    idx = d.index.to_frame(index=False)
    dv = d.values; n = idx.Turbines.values; ds = idx.Dataset.values
    strata = np.array([clusters.index((a, b)) for a, b in zip(idx.Dataset, idx.Radius)])
    big = n >= 10
    out = dict(n_cases=int(len(dv)), mean_large=float(dv[big].mean()), mean_small=float(dv[~big].mean()))
    out["interaction_all"] = dict(zip(("diff_pp", "p_unrestricted"), perm_diff(dv, big, np.zeros_like(strata), B, rng)))
    out["interaction_all"]["p_stratified"] = perm_diff(dv, big, strata, B, rng)[1]
    for tag, dsv in (("dsI", "1"), ("dsII", "2")):
        m = ds == dsv
        o, pu = perm_diff(dv[m], big[m], np.zeros(m.sum()), B, rng)
        ps = perm_diff(dv[m], big[m], strata[m], B, rng)[1]
        out[f"interaction_{tag}"] = dict(diff_pp=o, p_unrestricted=pu, p_stratified=ps,
                                         mean_large=float(dv[m & big].mean()), mean_small=float(dv[m & ~big].mean()),
                                         n_large=int((m & big).sum()), n_small=int((m & ~big).sum()))
    for tag, m in (("all", np.ones(len(dv), bool)), ("dsI", ds == "1"), ("dsII", ds == "2")):
        r = spearmanr(n[m], dv[m])
        mr, sp, rhos = strat_spearman(dv[m], n[m], strata[m], B, rng)
        out[f"spearman_{tag}"] = dict(rho=float(r.correlation), p=float(r.pvalue), n=int(m.sum()),
                                      stratified_mean_rho=mr, stratified_p=sp,
                                      rho_by_cluster={f"{clusters[s][0]}-{clusters[s][1]}": v for s, v in rhos.items()})
    # Data Set II vs Data Set I over identical (r, N): e = d_II - d_I (negative = the PSO-VNS advantage is larger in DS II)
    D2 = d.xs("2", level="Dataset"); D1 = d.xs("1", level="Dataset")
    com = D2.index.intersection(D1.index)
    e = (D2.loc[com] - D1.loc[com])
    nn = np.array([i[1] for i in com])
    for tag, m in (("all", np.ones(len(e), bool)), ("large", nn >= 10), ("small", nn < 10)):
        ev = e.values[m]
        p, rb, k = wil(ev)
        out[f"dsII_vs_dsI_{tag}"] = dict(n_pairs=int(m.sum()), mean_pp=float(ev.mean()), p=p, rb=rb,
                                         dsII_larger_gain=int((ev < -1e-9).sum()), dsI_larger_gain=int((ev > 1e-9).sum()))
    sv = out["spearman_dsII"]; iv = out["interaction_dsII"]; pv = out["dsII_vs_dsI_large"]
    out["pattern_survives"] = bool(iv["diff_pp"] < 0 and iv["p_stratified"] < 0.05 and sv["stratified_mean_rho"] < 0
                                   and sv["stratified_p"] < 0.05 and pv["mean_pp"] < 0 and pv["p"] < 0.05)
    out["overall_interaction_significant"] = bool(out["interaction_all"]["p_stratified"] < 0.05)
    out["overall_trend_significant"] = bool(out["spearman_all"]["stratified_p"] < 0.05)
    out["dsI_trend_negative_significant"] = bool(out["spearman_dsI"]["stratified_mean_rho"] < 0 and out["spearman_dsI"]["stratified_p"] < 0.05)
    out["permutations"] = B
    return out


# ------------------------------------------------------------------ 5. energy
def energy_block(S):
    out = {}
    specs = [("PSOBV", "PSOC", None), ("SSABV", "RSVNS", None), ("PSOBV", "RSVNS", None), ("PSOBV", "PSOC", "dsIILarge")]
    for a, b, grp in specs:
        Sx = S if grp is None else S[(S.Dataset == "2") & (S.Turbines >= 10)]
        dl = case_diffs(Sx, a, b)                                        # pp of wake-free AEP, L(a) - L(b)
        dobj = case_diffs(Sx, a, b, "Mean")                              # benchmark objective, a - b
        ob = Sx.pivot_table(index=CASE, columns="Algorithm", values="Mean")[b].loc[dobj.index]
        ideal = Sx.pivot_table(index=CASE, columns="Algorithm", values="Ideal")[b].loc[dobj.index]
        pct_aep = 100 * dobj / ob                                        # % of the second method's benchmark AEP
        mwh = dobj * ENERGY_FACTOR                                       # MWh/yr of benchmark energy (gain of a)
        nturb = np.array([i[2] for i in dobj.index])
        rng = np.random.default_rng(BOOT_SEED)
        bm = mwh.values[rng.integers(0, len(mwh), (BOOT_N, len(mwh)))].mean(1)
        rng = np.random.default_rng(BOOT_SEED)
        bp = pct_aep.values[rng.integers(0, len(mwh), (BOOT_N, len(mwh)))].mean(1)
        lg = nturb == nturb.max()
        key = pk(a, b) + ("" if grp is None else "_" + grp)
        out[key] = dict(a=a, b=b, group=grp or "all", n_cases=int(len(dl)), mean_dloss_pp=float(dl.mean()),
                        mean_gain_pct_aep=float(pct_aep.mean()), ci95_gain_pct_aep=[float(np.quantile(bp, .025)), float(np.quantile(bp, .975))],
                        max_gain_pct_aep=float(pct_aep.max()), min_gain_pct_aep=float(pct_aep.min()),
                        mean_gain_mwh_yr=float(mwh.mean()), ci95_gain_mwh_yr=[float(np.quantile(bm, .025)), float(np.quantile(bm, .975))],
                        mean_gain_mwh_yr_per_turbine=float((mwh / nturb).mean()),
                        largest_N=int(nturb.max()), mean_gain_mwh_yr_largest_N=float(mwh[lg].mean()),
                        mean_gain_pct_aep_largest_N=float(pct_aep[lg].mean()),
                        mean_aep_second_mwh_yr=float((ob * ENERGY_FACTOR).mean()),
                        mean_wakefree_aep_mwh_yr=float((ideal * ENERGY_FACTOR).mean()),
                        note="benchmark (model) energy: gross AEP = objective / 15 x 8.76 MWh/yr (Kusiak-Song benchmark, "
                             "capacity factor far above real sites); positive gain = first method yields more")
    return out


# ------------------------------------------------------------------ 6. wake-model uncertainty (Jensen vs Gaussian)
def model_shift_block(S, data_dir):
    """Wake-model uncertainty from the paper's own data: mpce_results.py (section 8, mpce_robustness.py) re-evaluates
    every feasible final layout of the eight main methods (68 cases) with the Bastankhah--Porte-Agel Gaussian wake
    (k* = 0.04, linear curve) and the benchmark Jensen model (column Linear: same stored coordinates, rounded to 1 mm)
    and writes them to mpce_reevaluation.csv. Wake loss of a layout under model M: 100 (Ideal - P_M) / Ideal (pp of
    the wake-free power, which does not depend on the wake model). Reported: |L_Gauss - L_Jensen| per layout, and the
    change of the PSO-VNS - PSO case-mean difference between the two models over the cases in which both methods have
    >= 15 feasible runs (qualification rule of the paper; case means over the feasible runs)."""
    fn = os.path.join(data_dir, "mpce_reevaluation.csv")
    if not os.path.exists(fn):
        fn = os.path.join(HERE, "mpce_reevaluation.csv")
    R = pd.read_csv(fn, usecols=["Algorithm", "Dataset", "Radius", "Turbines", "Seed", "Budget", "Init", "Feasible",
                                 "Objective", "Ideal", "Linear", "Gauss"])
    R["Dataset"] = R.Dataset.astype(str).str.replace(r"\.0$", "", regex=True)
    R["Feasible"] = R.Feasible.astype(str).str.lower().isin(["true", "1", "1.0"])
    R = R[(R.Budget == 6030) & (R.Init == "random") & R.Dataset.isin(["1", "2"]) & R.Feasible & R.Algorithm.isin(MAIN8)]
    assert R.Gauss.notna().all() and R.Linear.notna().all()
    R = R.assign(LJ=100 * (R.Ideal - R.Linear) / R.Ideal, LG=100 * (R.Ideal - R.Gauss) / R.Ideal,
                 LRec=100 * (R.Ideal - R.Objective) / R.Ideal)
    ad = (R.LG - R.LJ).abs()
    out = dict(source=os.path.basename(fn), n_layouts=int(len(R)), methods=sorted(R.Algorithm.unique()),
               n_cases=int(len(R[CASE].drop_duplicates())),
               mean_abs_shift_pp=float(ad.mean()), median_abs_shift_pp=float(ad.median()),
               q10_abs_shift_pp=float(ad.quantile(0.1)), q90_abs_shift_pp=float(ad.quantile(0.9)),
               min_abs_shift_pp=float(ad.min()), max_abs_shift_pp=float(ad.max()),
               share_abs_shift_above_margin=float((ad > EQ_MARGIN).mean()),
               mean_signed_shift_pp=float((R.LG - R.LJ).mean()),
               mean_loss_jensen_pp=float(R.LJ.mean()), mean_loss_gauss_pp=float(R.LG.mean()),
               max_abs_rounding_pp=float((R.LJ - R.LRec).abs().max()),
               by_dataset={ds: dict(mean_abs_shift_pp=float(ad[R.Dataset == ds].mean()),
                                    median_abs_shift_pp=float(ad[R.Dataset == ds].median()), n=int((R.Dataset == ds).sum()))
                           for ds in ("1", "2")},
               by_method={a: float(ad[R.Algorithm == a].mean()) for a in MAIN8},
               eq_margin_pp=EQ_MARGIN)
    out["ratio_mean_to_margin"] = out["mean_abs_shift_pp"] / EQ_MARGIN
    out["ratio_median_to_margin"] = out["median_abs_shift_pp"] / EQ_MARGIN
    # case means over the feasible runs, qualification >= BASE_THR feasible runs (as case_stats)
    g = R.groupby(CASE + ["Algorithm"])
    C = pd.DataFrame({"LJ": g.LJ.mean(), "LG": g.LG.mean(), "LRec": g.LRec.mean(), "NF": g.size()})
    pairs = {}
    for b in [x for x in MAIN8 if x != FOCUS]:
        A_, B_ = C.xs(FOCUS, level="Algorithm"), C.xs(b, level="Algorithm")
        com = A_.index.intersection(B_.index)
        ok = (A_.loc[com, "NF"] >= BASE_THR) & (B_.loc[com, "NF"] >= BASE_THR)
        com = com[ok.values]
        dJ = A_.loc[com, "LJ"] - B_.loc[com, "LJ"]; dG = A_.loc[com, "LG"] - B_.loc[com, "LG"]
        dR = A_.loc[com, "LRec"] - B_.loc[com, "LRec"]
        ch = (dG - dJ).abs()
        pairs[pk(FOCUS, b)] = dict(n_cases=int(len(com)), mean_abs_change_pp=float(ch.mean()),
                                   median_abs_change_pp=float(ch.median()), max_abs_change_pp=float(ch.max()),
                                   mean_diff_jensen_pp=float(dJ.mean()), mean_diff_gauss_pp=float(dG.mean()),
                                   mean_diff_recorded_pp=float(dR.mean()),
                                   abs_change_of_mean_pp=float(abs(dG.mean() - dJ.mean())),
                                   sign_changes=int((np.sign(dG) * np.sign(dJ) < 0).sum()),
                                   wilcoxon_p_jensen=wil(dJ.values)[0], wilcoxon_p_gauss=wil(dG.values)[0])
    out["pairs"] = pairs
    # consistency: the recorded-objective PSO-VNS - PSO difference equals the paper's (mpce_results case means)
    ref = case_diffs(S, FOCUS, "PSOC").mean()
    out["pair_recorded_vs_paper_abs_diff"] = float(abs(pairs[pk(FOCUS, "PSOC")]["mean_diff_recorded_pp"] - ref))
    out["note"] = ("|L_Gauss - L_Jensen| of the same feasible final layouts (8 main methods, 68 cases), pp of the wake-free "
                   "power; pair: change of the per-case mean difference PSO-VNS - PSO between the two wake models")
    return out


# ------------------------------------------------------------------ 7. equivalence at three levels (review round 2, D13)
def boot_summary(bm, est, margin=EQ_MARGIN):
    """TOST summary of bootstrap means bm, exactly as mpce_results.tost: 90 % / 95 % percentile CIs, equivalence iff
    the 90 % CI lies strictly inside (-m, m), bootstrap TOST p = max(share <= -m, share >= +m) with add-one correction
    (an inversion of the percentile interval, not a p value calibrated at the null boundary), minimal margin
    max(|lo90|, |hi90|)."""
    B = len(bm)
    lo90, hi90 = (float(v) for v in np.quantile(bm, [0.05, 0.95]))
    lo95, hi95 = (float(v) for v in np.quantile(bm, [0.025, 0.975]))
    k_lo, k_hi = int((bm <= -margin).sum()), int((bm >= margin).sum())
    return dict(mean_dloss_pp=float(est), ci90=[lo90, hi90], ci95=[lo95, hi95],
                p_tost=float(max(k_lo + 1, k_hi + 1) / (B + 1)), p_tost_at_floor=bool(k_lo == 0 and k_hi == 0),
                min_margin_pp=float(max(abs(lo90), abs(hi90))), equivalent=bool(-margin < lo90 and hi90 < margin))


def boot_p_curve(bm, grid):
    """bootstrap TOST p (definition of boot_summary) for every margin of grid."""
    s = np.sort(bm); B = len(s)
    k_lo = np.searchsorted(s, -grid, side="right"); k_hi = B - np.searchsorted(s, grid, side="left")
    return np.maximum(k_lo + 1, k_hi + 1) / (B + 1)


def case_boot_means(d, n=BOOT_N, seed=BOOT_SEED):
    """bootstrap means over the cases -- the same resamples as mpce_results.boot_means / tost (reproduced, X38)."""
    d = np.asarray(d, float)
    return d[np.random.default_rng(seed).integers(0, len(d), (n, len(d)))].mean(1)


def run_arrays(G, a, b, cases):
    """per case: (loss %, feasible) of the 30 seed-paired runs of a and b, ordered by seed (identical seeds asserted)."""
    R = G[G.Algorithm.isin([a, b])].sort_values(CASE + ["Algorithm", "Seed"])
    out = {}
    for c, g in R.groupby(CASE):
        if c not in cases:
            continue
        ga, gb = g[g.Algorithm == a], g[g.Algorithm == b]
        assert (ga.Seed.values == gb.Seed.values).all() and len(ga) == 30, (a, b, c)
        out[c] = (ga.LossPct.values, ga.Feasible.values.astype(bool), gb.LossPct.values, gb.Feasible.values.astype(bool))
    assert len(out) == len(cases)
    return out


def seed_level(G, a, b, d, B=BOOT_N, seed=SEED_BOOT_SEED):
    """(a) Fixed-benchmark (seed-level) inference: the 68 cases are fixed, only the run-to-run (seed) variability is
    random. Bootstrap: within every case the 30 seed-paired runs are resampled with replacement (the SAME seeds for
    both methods), the case mean of each method is the mean over the feasible resampled runs (as the pipeline), and
    the estimate is the mean of the case differences over the cases of d (both methods qualified). Stratified:
    independent seed resamples per case (runs of different cases independent); joint: one seed resample for all cases
    (allows for dependence between cases that share a seed number). If a resample has no feasible run of a method in
    a case (never for >= 15 feasible of 30 in practice) the case mean of the full data is used (counted)."""
    cases = list(d.index)
    RA = run_arrays(G, a, b, set(cases))
    rng_s = np.random.default_rng(seed); Ij = np.random.default_rng(seed + 1).integers(0, 30, (B, 30))
    tot_s = np.zeros(B); tot_j = np.zeros(B); fb = 0; est = 0.0
    all_feas = True; per_seed = np.zeros(30)

    def cmean(L, F, I, full):
        nf = F[I].sum(1); s = np.where(F, L, 0.0)[I].sum(1)
        m = np.where(nf > 0, s / np.maximum(nf, 1), full)
        return m, int((nf == 0).sum())

    for c in cases:
        La, Fa, Lb, Fb = RA[c]
        ma, mb = La[Fa].mean(), Lb[Fb].mean()
        assert abs((ma - mb) - d.loc[c]) < 1e-9, ("case mean mismatch", a, b, c)
        est += ma - mb
        all_feas &= bool(Fa.all() and Fb.all())
        per_seed += La - Lb
        Is = rng_s.integers(0, 30, (B, 30))
        for I, tot in ((Is, tot_s), (Ij, tot_j)):
            xa, fa = cmean(La, Fa, I, ma); xb, fb_ = cmean(Lb, Fb, I, mb)
            tot += xa - xb; fb += fa + fb_
    n = len(cases); bm_s, bm_j = tot_s / n, tot_j / n
    r = boot_summary(bm_s, est / n)
    r.update(n_cases=n, resamples=B, seed=seed, fallback_resamples=fb, all_runs_feasible=all_feas,
             joint=boot_summary(bm_j, est / n))
    if all_feas:                                                       # reviewer's version: benchmark average per seed, t(29)
        ps = per_seed / n; se = ps.std(ddof=1) / math.sqrt(30); q = t_dist.ppf(0.95, 29)
        r["per_seed_t"] = dict(mean=float(ps.mean()), ci90=[float(ps.mean() - q * se), float(ps.mean() + q * se)],
                               se=float(se), n_seeds=30,
                               p_tost=float(max(t_dist.sf((ps.mean() + EQ_MARGIN) / se, 29), t_dist.cdf((ps.mean() - EQ_MARGIN) / se, 29))))
    return r, bm_s


def cr_parts(d, cid, G):
    n = len(d); ng = np.bincount(cid, minlength=G).astype(float)
    return n, ng, 1.0 / np.sqrt(1.0 - ng / n)


def cr2(d, cid, G):
    """(c) cluster-robust inference on the case-weighted mean: CR2 variance (Bell-McCaffrey bias-reduced
    linearization; for the intercept-only model the cluster sums of the residuals are scaled by (1 - n_g/n)^-1/2),
    t with G - 1 = 5 d.f. (primary, as requested) and the Bell-McCaffrey / Imbens-Kolesar d.f. (homoskedastic working
    model) for reference; CR1 (G/(G-1) scaling) for comparison with tab:X-loco."""
    d = np.asarray(d, float); n, ng, adj = cr_parts(d, cid, G)
    est = d.mean(); U = np.bincount(cid, weights=d - est, minlength=G)
    V2 = ((U * adj) ** 2).sum() / n ** 2; V1 = G / (G - 1) * (U ** 2).sum() / n ** 2
    A = np.zeros((n, G))
    for g in range(G):
        A[cid == g, g] = adj[g] / n
    M = np.eye(n) - 1.0 / n
    K = M @ (A @ A.T) @ M
    df_bm = float(np.trace(K) ** 2 / np.trace(K @ K))
    out = dict(mean=float(est), se_cr2=float(math.sqrt(V2)), se_cr1=float(math.sqrt(V1)), df=G - 1, df_bm=df_bm)
    for nm, se, df in (("cr2", math.sqrt(V2), G - 1), ("cr2_bm", math.sqrt(V2), df_bm), ("cr1", math.sqrt(V1), G - 1)):
        q = t_dist.ppf(0.95, df)
        lo, hi = est - q * se, est + q * se
        out[nm] = dict(ci90=[float(lo), float(hi)], min_margin_pp=float(max(abs(lo), abs(hi))),
                       p_tost=float(max(t_dist.sf((est + EQ_MARGIN) / se, df), t_dist.cdf((est - EQ_MARGIN) / se, df))),
                       equivalent=bool(-EQ_MARGIN < lo and hi < EQ_MARGIN),
                       ci95=[float(est - t_dist.ppf(0.975, df) * se), float(est + t_dist.ppf(0.975, df) * se)])
    return out


def webb_all(G):
    return np.array(list(itertools.product(WEBB, repeat=G)))              # 6^G draws (46,656 for G = 6), enumerated


def wild_p(d, cid, G, mu0, W):
    """Restricted wild-cluster bootstrap-t for H0: mean = mu0 (x* = mu0 + w_g (x - mu0)), CR2-studentized, all Webb
    draws W enumerated. Returns (p for H1: mean > mu0, p for H1: mean < mu0)."""
    d = np.asarray(d, float); n, ng, adj = cr_parts(d, cid, G)
    U = np.bincount(cid, weights=d - mu0, minlength=G)
    t_obs = (d.mean() - mu0) / math.sqrt((((U - ng * (d.mean() - mu0)) * adj) ** 2).sum() / n ** 2)
    WU = W * U; S = WU.sum(1)
    E = WU - ng * (S / n)[:, None]
    t = (S / n) / np.sqrt(((E * adj) ** 2).sum(1) / n ** 2)
    return float(np.mean(t >= t_obs - 1e-10)), float(np.mean(t <= t_obs + 1e-10))


def wild_level(d, cid, G, W, grid=None):
    """wild-cluster TOST at EQ_MARGIN, 90 % CI by test inversion (bisection on mu0: lower = boundary where the
    one-sided p for H1: mean > mu0 crosses 0.05, upper likewise) and minimal margin max(-lo, hi); optional p curve."""
    d = np.asarray(d, float); est = d.mean()

    def bound(side):
        f = (lambda mu: wild_p(d, cid, G, mu, W)[0]) if side == "lo" else (lambda mu: wild_p(d, cid, G, mu, W)[1])
        inner, outer = est, est - 2.0 if side == "lo" else est + 2.0
        assert f(inner) > 0.05 and f(outer) <= 0.05, ("bracket", side)
        for _ in range(40):
            mid = 0.5 * (inner + outer)
            if f(mid) > 0.05:
                inner = mid
            else:
                outer = mid
        return 0.5 * (inner + outer)
    lo, hi = bound("lo"), bound("hi")
    pl, ph = wild_p(d, cid, G, -EQ_MARGIN, W)[0], wild_p(d, cid, G, EQ_MARGIN, W)[1]
    out = dict(ci90=[float(lo), float(hi)], min_margin_pp=float(max(-lo, hi)), p_tost=float(max(pl, ph)),
               p_low=pl, p_high=ph, equivalent=bool(max(pl, ph) <= 0.05), draws=int(len(W)), weights="Webb 6-point",
               studentization="CR2")
    if grid is not None:
        out["curve"] = [float(max(wild_p(d, cid, G, -m, W)[0], wild_p(d, cid, G, m, W)[1])) for m in grid]
    return out


def bayes_signrank(d, rope=EQ_MARGIN, s=BAYES_S, z0=BAYES_Z0, n=BAYES_N, seed=BAYES_SEED):
    """Bayesian signed-rank test (Benavoli et al. 2017), a copy of mpce_results.bayes_signrank (same prior, samples and
    seed; reproduced in check X38) extended by 95 % credible intervals of theta. P_* = share of the posterior samples in
    which theta_A / theta_rope / theta_B is the LARGEST of the three, i.e. the posterior probability that each region
    is the most probable one -- NOT the probability that the mean difference lies in the rope. mean_theta = posterior
    means of (theta_A, theta_rope, theta_B): the expected shares of Walsh averages (d_i + d_j)/2 below -rope, inside,
    above +rope."""
    d = np.asarray(d, float); d = d[np.isfinite(d)]
    z = np.r_[z0, d]; W = np.random.default_rng(seed).dirichlet(np.r_[s, np.ones(len(d))], n)
    sums = z[:, None] + z[None, :]
    below = (sums < -2 * rope) + 0.5 * (sums == -2 * rope); above = (sums > 2 * rope) + 0.5 * (sums == 2 * rope)
    tA = ((W @ below) * W).sum(1); tB = ((W @ above) * W).sum(1)
    T = np.c_[tA, 1 - tA - tB, tB]
    win = (T >= T.max(1, keepdims=True) - 1e-15).astype(float)
    pr = (win / win.sum(1, keepdims=True)).mean(0)
    return dict(p_a_better=float(pr[0]), p_rope=float(pr[1]), p_b_better=float(pr[2]), n_cases=int(len(d)),
                mean_theta=[float(v) for v in T.mean(0)],
                ci95_theta=[[float(v) for v in np.quantile(T[:, k], [0.025, 0.975])] for k in range(3)])


def equivalence_levels_block(G, S, clusters, grid=CURVE_GRID):
    """Equivalence of every tab:equivalence pair at three inference levels: (a) seed level (fixed benchmark),
    (b) case level (bootstrap over cases, reproduces mpce_results.equivalence_block), (c) cluster level (CR2, 5 d.f.;
    restricted wild-cluster bootstrap-t TOST with all Webb draws); plus the Bayesian signed-rank posterior means."""
    Gn = len(clusters); W = webb_all(Gn)
    res, curve = {}, {}
    for key in EQ_PAIRS:
        a, b = key.split("-")
        d = case_diffs(S, a, b)
        cid = np.array([clusters.index((c[0], c[1])) for c in d.index])
        dv = d.values
        bm_c = case_boot_means(dv)
        case = boot_summary(bm_c, dv.mean())
        se = dv.std(ddof=1) / math.sqrt(len(dv))
        case["p_tost_t"] = float(max(t_dist.sf((dv.mean() + EQ_MARGIN) / se, len(dv) - 1),
                                     t_dist.cdf((dv.mean() - EQ_MARGIN) / se, len(dv) - 1)))
        seed, bm_s = seed_level(G, a, b, d)
        focus = key == "PSOBV-PSOC"
        cl = cr2(dv, cid, Gn)
        cl["wild"] = wild_level(dv, cid, Gn, W, grid if focus else None)
        cl["n_clusters"] = Gn
        cl["equivalent"] = bool(cl["cr2"]["equivalent"] and cl["wild"]["equivalent"])
        res[key] = dict(a=a, b=b, n_cases=int(len(dv)), mean_dloss_pp=float(dv.mean()), seed=seed, case=case, cluster=cl,
                        bayes=bayes_signrank(dv))
        if focus:
            curve = dict(pair=key, margins=[float(m) for m in grid], seed=[float(v) for v in boot_p_curve(bm_s, grid)],
                         case=[float(v) for v in boot_p_curve(bm_c, grid)],
                         cluster_cr2=[float(max(t_dist.sf((dv.mean() + m) / cl["se_cr2"], Gn - 1),
                                                t_dist.cdf((dv.mean() - m) / cl["se_cr2"], Gn - 1))) for m in grid],
                         cluster_wild=cl["wild"].pop("curve"))
        print(f"EQ3 {key:14s} n={len(dv)} dL={dv.mean():+.4f} | seed 90% [{seed['ci90'][0]:+.4f},{seed['ci90'][1]:+.4f}] "
              f"m={seed['min_margin_pp']:.4f} | case [{case['ci90'][0]:+.4f},{case['ci90'][1]:+.4f}] m={case['min_margin_pp']:.4f} "
              f"p={case['p_tost']:.3g} | CR2 [{cl['cr2']['ci90'][0]:+.4f},{cl['cr2']['ci90'][1]:+.4f}] m={cl['cr2']['min_margin_pp']:.4f} "
              f"(df_BM {cl['df_bm']:.2f}) wild p={cl['wild']['p_tost']:.3g} m={cl['wild']['min_margin_pp']:.4f} | "
              f"theta {np.round(res[key]['bayes']['mean_theta'], 3)}")
    return dict(margin_pp=EQ_MARGIN, margin_source=EQ_MARGIN_SRC, pairs=res, curve=curve,
                levels=dict(seed="fixed benchmark: runs resampled within cases (seed-paired, stratified by case; "
                                 f"{BOOT_N} resamples, seed {SEED_BOOT_SEED}); joint = one seed resample for all cases",
                            case="bootstrap over the cases (as mpce_results.tost; resamples/seed of the pipeline)",
                            cluster="CR2 variance, t with G-1 d.f.; restricted wild-cluster bootstrap-t (Webb 6-point, all "
                                    "6^G draws enumerated, CR2-studentized) TOST and CI by test inversion"))


def reproduce_equivalence(EL):
    """check X38: case level and Bayesian test reproduce mpce_summary.json["equivalence"]."""
    try:
        E = json.load(open(os.path.join(HERE, "mpce_summary.json")))["equivalence"]
        assert list(E["pairs"]) == EQ_PAIRS, f"tab:equivalence pairs differ: {list(E['pairs'])}"
        dev = 0.0
        for k, r in EL["pairs"].items():
            e = E["pairs"][k]; c = r["case"]; y = r["bayes"]
            dev = max(dev, abs(c["mean_dloss_pp"] - e["mean_dloss_pp"]), abs(c["p_tost"] - e["p_tost"]),
                      abs(c["min_margin_pp"] - e["min_margin_pp"]), abs(c["p_tost_t"] - e["p_tost_t"]),
                      *[abs(u - v) for u, v in zip(c["ci90"], e["ci90_mean_dloss_pp"])],
                      *[abs(u - v) for u, v in zip(c["ci95"], e["ci95_mean_dloss_pp"])],
                      abs(y["p_a_better"] - e["bayes"]["p_a_better"]), abs(y["p_rope"] - e["bayes"]["p_rope"]),
                      abs(y["p_b_better"] - e["bayes"]["p_b_better"]),
                      *[abs(u - v) for u, v in zip(y["mean_theta"], e["bayes"]["mean_theta"])])
            assert c["equivalent"] == e["equivalent"] and r["n_cases"] == e["n_cases"], k
        return dict(ok=bool(dev < 1e-9), max_abs_diff=float(dev), margin_pipeline=E["margin_pp"])
    except Exception as ex:                                              # pragma: no cover
        return dict(ok=False, error=repr(ex))


# ------------------------------------------------------------------ 8. heterogeneity of PSO-VNS - PSO
def ols_hc3_parts(X):
    """for a fixed design X: P = (X'X)^-1 X' and leverages h (HC3 variance of coefficient k: sum_i P_ki^2 r_i^2/(1-h_i)^2)."""
    P = np.linalg.solve(X.T @ X, X.T)
    return P, np.einsum("ij,ji->i", X, P)


def studentized_slope(Y, X, P, h, k):
    """coefficient k and its HC3 t statistic for every row of Y (vectorized over permutations)."""
    beta = Y @ P.T; r = Y - beta @ X.T
    se = np.sqrt(((r / (1 - h)) ** 2) @ (P[k] ** 2))
    return beta[:, k], beta[:, k] / se


def heterogeneity_block(S, clusters, B, EL):
    """Per-case heterogeneity of d = L(PSO-VNS) - L(PSO) (pp): cases beyond +-margin, stratified equivalence (N < 10,
    N >= 10), cluster means, and a data set x density interaction.

    Density (stated before this analysis was run, after review round 2): phi = N (l_min / 2r)^2, the nominal share of
    the farm disc covered by the N exclusion discs of diameter l_min = 4D = 308 m (the same ordering as N/r^2).
    Interaction test: the two data sets share the identical 34 (r, N) cases, so e_j = d_II(r,N) - d_I(r,N) is regressed
    on phi_j; the slope is the data set x density interaction (d ~ (r,N) pair + data set + data set:phi). Test:
    studentized (HC3) slope, two-sided, signs of e_j flipped at random (B Monte Carlo draws, fixed seed); the
    studentization keeps the test valid for H0 slope = 0 when the data sets also differ in level (Janssen 1997; Chung
    and Romano 2013). Direction within each data set: d ~ radius fixed effects + phi (34 cases), HC3-studentized slope,
    permutation of d within radius (B draws)."""
    rng = np.random.default_rng(PERM_SEED + 7)
    d = case_diffs(S, FOCUS, "PSOC")
    idx = d.index.to_frame(index=False); dv = d.values
    n = idx.Turbines.values; r = idx.Radius.values.astype(float); ds = idx.Dataset.values
    phi = n * (SMIN_M / (2 * r)) ** 2
    m = EQ_MARGIN
    out = dict(n_cases=int(len(dv)), margin_pp=m,
               beyond_first=int((dv < -m).sum()), beyond_second=int((dv > m).sum()),
               max_first=float(-dv.min()), max_second=float(dv.max()),
               exact_ties=int((np.abs(dv) <= 1e-9).sum()), near_ties=int((np.abs(dv) < 0.005).sum()),
               density_def="phi = N (l_min / 2r)^2, l_min = 308 m")
    out["cases_beyond"] = [dict(dataset=str(a), radius=int(b_), N=int(c), phi=float(p), d=float(v))
                           for a, b_, c, p, v in zip(ds, r, n, phi, dv) if abs(v) > m]
    # stratified by N, and by the size of the wake loss itself (lead request: cases with a trivial wake loss cannot
    # differ by the margin): case-mean wake loss of PSO-VNS < TRIVIAL_LOSS_PP ("trivial") vs >= ("nontrivial")
    lossA = S[S.Algorithm == FOCUS].set_index(CASE).Loss.loc[d.index].values
    out["trivial_loss_pp"] = TRIVIAL_LOSS_PP
    for tag, msk in (("small", n < 10), ("large", n >= 10), ("trivial", lossA < TRIVIAL_LOSS_PP),
                     ("nontrivial", lossA >= TRIVIAL_LOSS_PP)):
        x = dv[msk]
        bs = boot_summary(case_boot_means(x), x.mean())
        bs.update(n_cases=int(msk.sum()), bayes=bayes_signrank(x),
                  wins_first=int((x < -1e-9).sum()), wins_second=int((x > 1e-9).sum()),
                  beyond_first=int((x < -m).sum()), beyond_second=int((x > m).sum()),
                  mean_abs=float(np.abs(x).mean()), loss_range_first=[float(lossA[msk].min()), float(lossA[msk].max())],
                  N_range=[int(n[msk].min()), int(n[msk].max())])
        out[tag] = bs
    # the margin in benchmark energy: margin_pp / 100 x wake-free AEP per turbine (benchmark AEP = objective / 15 x
    # 8.76 MWh/yr, ENERGY_FACTOR, as the \NXEn macros); Ideal = wake-free objective of the case (same for all methods)
    ide = S[S.Algorithm == FOCUS].set_index(CASE).Ideal.loc[d.index].values
    per_t = ide * ENERGY_FACTOR / n
    out["margin_energy"] = {f"ds{ 'I' * int(v) }": dict(wakefree_mwh_yr_per_turbine=float(per_t[ds == v].mean()),
                                                         wakefree_per_turbine_spread=float(per_t[ds == v].max() - per_t[ds == v].min()),
                                                         margin_mwh_yr_per_turbine=float(m / 100 * per_t[ds == v].mean()))
                            for v in ("1", "2")}
    # cluster means
    cm = {}
    for (dsv, rv) in clusters:
        mk = (ds == dsv) & (r == rv)
        cm[f"{dsv}-{rv}"] = float(dv[mk].mean())
    out["cluster_means"] = cm
    out["clusters_beyond_first"] = [k for k, v in cm.items() if v < -m]
    out["clusters_beyond_second"] = [k for k, v in cm.items() if v > m]
    # data set x density interaction on the matched (r, N) pairs
    D1 = d.xs("1", level="Dataset"); D2 = d.xs("2", level="Dataset")
    com = D1.index.intersection(D2.index)
    assert len(com) == len(D1) == len(D2), "data sets must share the (r, N) cases"
    e = (D2.loc[com] - D1.loc[com]).values
    ph = np.array([c[1] * (SMIN_M / (2 * c[0])) ** 2 for c in com])
    X = np.c_[np.ones(len(e)), ph]; P, h = ols_hc3_parts(X)
    b_obs, t_obs = studentized_slope(e[None, :], X, P, h, 1)
    sg = rng.choice([-1.0, 1.0], (B, len(e)))
    _, tb = studentized_slope(sg * e, X, P, h, 1)
    p_int = float(((np.abs(tb) >= abs(t_obs[0]) - 1e-12).sum() + 1) / (B + 1))
    out["interaction"] = dict(n_pairs=int(len(e)), slope_pp_per_unit_phi=float(b_obs[0]), t_hc3=float(t_obs[0]),
                              p_signflip=p_int, perm=B, df=int(len(e) - 2),
                              p_t=float(2 * t_dist.sf(abs(t_obs[0]), len(e) - 2)), phi_range=[float(ph.min()), float(ph.max())],
                              note="slope of d_II - d_I on phi over the 34 identical (r, N) cases; negative = the PSO-VNS "
                                   "advantage grows with density more in Data Set II than in Data Set I")
    for tag, dsv in (("dsI", "1"), ("dsII", "2")):
        mk = ds == dsv; y = dv[mk]; rr = r[mk]; pp = phi[mk]
        rad = sorted(set(rr))
        X = np.c_[np.array([[float(v == q) for q in rad] for v in rr]), pp]; k = X.shape[1] - 1
        P, h = ols_hc3_parts(X)
        b_obs, t_obs = studentized_slope(y[None, :], X, P, h, k)
        Y = np.tile(y, (B, 1))
        for q in rad:
            ix = np.where(rr == q)[0]
            Y[:, ix] = rng.permuted(Y[:, ix], axis=1)
        _, tb = studentized_slope(Y, X, P, h, k)
        pv = float(((np.abs(tb) >= abs(t_obs[0]) - 1e-12).sum() + 1) / (B + 1))
        dense = [c for c in out["cases_beyond"] if c["dataset"] == dsv]
        out[f"slope_{tag}"] = dict(n=int(mk.sum()), slope_pp_per_unit_phi=float(b_obs[0]), t_hc3=float(t_obs[0]),
                                   p_perm=pv, perm=B, beyond_cases=dense,
                                   mean_top_density=float(y[pp >= np.quantile(pp, 0.75)].mean()))
    out["crossover"] = bool(out["slope_dsI"]["slope_pp_per_unit_phi"] > 0 > out["slope_dsII"]["slope_pp_per_unit_phi"])
    print(f"HET beyond +-{m}: {out['beyond_first']} favour PSO-VNS (max {out['max_first']:.3f}), {out['beyond_second']} favour PSO "
          f"(max {out['max_second']:.3f}); ties {out['exact_ties']}")
    print("HET margin energy", out["margin_energy"])
    for tag in ("small", "large", "trivial", "nontrivial"):
        v = out[tag]
        print(f"HET {tag}: n={v['n_cases']} mean {v['mean_dloss_pp']:+.4f} 90% {np.round(v['ci90'], 4)} p {v['p_tost']:.3g} "
              f"eq {v['equivalent']} bayes {v['bayes']['p_a_better']:.3f}/{v['bayes']['p_rope']:.3f}/{v['bayes']['p_b_better']:.3f}")
    print("HET clusters", {k: round(v, 3) for k, v in cm.items()})
    print("HET interaction", out["interaction"], "\n    dsI", {k: v for k, v in out["slope_dsI"].items() if k != "beyond_cases"},
          "\n    dsII", {k: v for k, v in out["slope_dsII"].items() if k != "beyond_cases"})
    return out


# ------------------------------------------------------------------ 9. Hodges-Lehmann estimates
def signrank_cdf(n):
    """exact null CDF of the Wilcoxon signed-rank statistic T+ for n untied non-zero differences."""
    pmf = np.zeros(n * (n + 1) // 2 + 1); pmf[0] = 1.0
    for k in range(1, n + 1):
        new = pmf * 0.5
        new[k:] += 0.5 * pmf[:-k]
        pmf = new
    return np.cumsum(pmf)


def hodges_lehmann(d, alpha=0.05):
    """HL estimate (median of the n(n+1)/2 Walsh averages, i <= j) and the (1 - alpha) CI from the exact inversion of
    the Wilcoxon signed-rank test (as R wilcox.test(conf.int = TRUE)): [W_(k), W_(M+1-k)] with k = the smallest t such
    that P(T+ <= t) >= alpha/2 under H0; exact for continuous data (zero differences / ties make it approximate)."""
    d = np.sort(np.asarray(d, float)); n = len(d); iu = np.triu_indices(n)
    w = np.sort(((d[:, None] + d[None, :]) / 2)[iu]); M = len(w)
    cdf = signrank_cdf(n)
    k = max(int(np.argmax(cdf >= alpha / 2 - 1e-12)), 1)
    return dict(hl=float(np.median(w)), ci95=[float(w[k - 1]), float(w[M - k])], n=int(n), n_walsh=int(M),
                achieved_conf=float(1 - 2 * cdf[k - 1]), median=float(np.median(d)))


def hl_block(S, case_level):
    out = {}
    for a, b in ALL_PAIRS:
        d = case_diffs(S, a, b)
        r = hodges_lehmann(d.values)
        cl = case_level[pk(a, b)]
        rb = case_boot_means(d.values)
        r.update(a=a, b=b, mean=float(d.mean()), ci95_mean=[float(np.quantile(rb, 0.025)), float(np.quantile(rb, 0.975))],
                 wilcoxon_p=cl["p"], wins_first=cl["wins"], wins_second=cl["losses"], n_wilcoxon=cl["n_cases"],
                 key_pair=(a, b) in HL_KEY, excludes_zero=bool(r["ci95"][0] > 0 or r["ci95"][1] < 0))
        out[pk(a, b)] = r
        if (a, b) in HL_KEY:
            print(f"HL {pk(a, b):14s} HL {r['hl']:+.4f} [{r['ci95'][0]:+.4f}, {r['ci95'][1]:+.4f}] mean {r['mean']:+.4f} "
                  f"{np.round(r['ci95_mean'], 4)} Wilcoxon p {cl['p']:.3g}")
    return out


# ------------------------------------------------------------------ 10. McNemar, Horns Rev feasibility
def mcnemar_block(data_dir):
    """Paired exact McNemar test of feasibility, PSO-VNS vs PSO, Horns Rev 16 turbines, 6,030 evaluations, random
    initialization (experiment hrfix, current direction binning), runs paired by seed."""
    H = read_shards("hrfix", data_dir)
    H = H[(H.Dataset == "HR") & (H.Budget == 6030) & (H.Init == "random") & H.Algorithm.isin([FOCUS, "PSOC"])]
    P = H.pivot_table(index="Seed", columns="Algorithm", values="Feasible", aggfunc="first").astype(bool)
    assert len(P) == 30 and P.notna().all().all()
    a, b = P[FOCUS].values, P["PSOC"].values
    n10, n01 = int((a & ~b).sum()), int((~a & b).sum())
    p = float(binomtest(min(n10, n01), n10 + n01, 0.5).pvalue) if n10 + n01 else 1.0
    out = dict(n_seeds=int(len(P)), both=int((a & b).sum()), neither=int((~a & ~b).sum()), only_first=n10, only_second=n01,
               feasible_first=int(a.sum()), feasible_second=int(b.sum()), p_exact=p, source="mpce_hrfix (6,030, random)")
    try:
        sp = json.load(open(os.path.join(HERE, "mpce_summary.json")))["spread"]["hr16_feasible"]
        out["matches_summary"] = bool(sp[FOCUS]["feasible"] == out["feasible_first"] and sp["PSOC"]["feasible"] == out["feasible_second"])
    except Exception as ex:                                              # pragma: no cover
        out["matches_summary"] = False; out["error"] = repr(ex)
    print("McNemar HR", out)
    return out


# ------------------------------------------------------------------ 11. figure: equivalence curve
def equiv_figure(EL, fn):
    """TOST p vs margin for PSO-VNS - PSO at the seed, case and cluster level (figures_mpce/equiv_curve.pdf)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    INK, MUTED, GRID = "#0b0b0b", "#52514e", "#e4e3df"
    # figure style shared by all figures of the paper (= mpce_results.FIG_STYLE; STIX text, drawn at the printed
    # width 0.45 x 16 cm of the main text, so that the fonts print at their nominal size)
    plt.rcParams.update({"font.family": "STIXGeneral", "mathtext.fontset": "stix", "font.size": 8, "axes.titlesize": 8,
                         "axes.labelsize": 8, "xtick.labelsize": 7.5, "ytick.labelsize": 7.5, "legend.fontsize": 7.5,
                         "axes.edgecolor": MUTED, "axes.labelcolor": INK, "text.color": INK, "xtick.color": MUTED,
                         "ytick.color": MUTED, "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
                         "xtick.minor.width": 0.4, "ytick.minor.width": 0.4, "xtick.major.size": 2.5,
                         "ytick.major.size": 2.5, "xtick.minor.size": 1.5, "ytick.minor.size": 1.5, "xtick.major.pad": 2,
                         "ytick.major.pad": 2, "axes.labelpad": 2.5, "axes.titlepad": 3.5, "axes.grid": True,
                         "grid.color": GRID, "grid.linewidth": 0.5, "axes.spines.top": False, "axes.spines.right": False,
                         "legend.frameon": False, "legend.handlelength": 2.2, "legend.columnspacing": 1.2,
                         "legend.handletextpad": 0.5, "legend.borderaxespad": 0.2, "lines.linewidth": 1.2,
                         "savefig.pad_inches": 0.02, "pdf.fonttype": 42})
    C = EL["curve"]; x = np.array(C["margins"]); r = EL["pairs"][C["pair"]]
    # Okabe-Ito colours (colour-blind safe); the two cluster-level curves share green and differ by dash pattern
    lines = [("seed", "#0072b2", "-", "seed level", r["seed"]["min_margin_pp"]),
             ("case", "#d55e00", "-", "case level", r["case"]["min_margin_pp"]),
             ("cluster_cr2", "#009e73", "-", "cluster level, CR2 $t_5$", r["cluster"]["cr2"]["min_margin_pp"]),
             ("cluster_wild", "#009e73", (0, (3, 1.5)), "cluster level, wild bootstrap", r["cluster"]["wild"]["min_margin_pp"])]
    fig, ax = plt.subplots(figsize=(0.45 * 16.0 / 2.54, 2.5))
    floor = 1.0 / (BOOT_N + 1)
    ax.axhline(0.05, color=MUTED, lw=0.8, ls=":", zorder=1)
    ax.axvline(EQ_MARGIN, color=MUTED, lw=0.8, ls=":", zorder=1)
    ax.text(0.001, 0.05 * 0.8, r"$\alpha=0.05$", color=MUTED, ha="left", va="top", fontsize=7.5)
    ax.text(EQ_MARGIN - 0.001, 1.2e-3, f"$m$ = {EQ_MARGIN:g} pp", color=MUTED, ha="right", va="center", fontsize=7.5,
            rotation=90)
    # legend below the axes (keeps the curves free); the dot on the alpha line marks each minimal margin
    hs = []
    for k, col, ls, lab, mm in lines:
        y = np.array(C[k], float)
        if k in ("seed", "case"):                               # bootstrap p has the resolution floor 1/(B+1): stop there
            hit = np.where(y <= floor + 1e-12)[0]
            if len(hit):
                x_, y_ = x[:hit[0] + 1], np.maximum(y[:hit[0] + 1], floor)
            else:
                x_, y_ = x, y
        else:
            x_, y_ = x, y
        h, = ax.plot(x_, y_, color=col, ls=ls, lw=1.5 if k != "cluster_wild" else 1.3, zorder=3,
                     label=f"{lab}, $m_{{\min}}={mm:.3f}$")
        hs.append(h)
        ax.plot([mm], [0.05], marker="o", ms=4.5, color=col, mec="white", mew=0.6, zorder=4)
    ax.set_yscale("log"); ax.set_ylim(8e-5, 1.3); ax.set_xlim(0, 0.12)
    ax.set_xticks([0, 0.02, 0.04, 0.06, 0.08, 0.10, 0.12])
    ax.set_xticklabels(["0", "0.02", "0.04", "0.06", "0.08", "0.10", "0.12"])
    ax.set_xlabel("equivalence margin $m$ (pp of wake loss)")
    ax.set_ylabel(r"TOST $p$ (PSO-VNS $-$ PSO)")
    ax.legend(handles=hs, loc="upper left", bbox_to_anchor=(-0.2, -0.2), ncol=1, handlelength=2.4, labelspacing=0.25)
    fig.tight_layout(pad=0.3)
    fig.savefig(fn, bbox_inches="tight")
    plt.close(fig)
    print("wrote", fn)


# ------------------------------------------------------------------ main
def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default=HERE)
    ap.add_argument("--out-dir", default=HERE)
    ap.add_argument("--perm", type=int, default=20000)
    ap.add_argument("--fig-dir", default=os.path.join(os.path.dirname(HERE), "figures_mpce"))
    args = ap.parse_args(argv)
    t0 = time.time()
    A = load(args.data_dir)
    G = A[A.Algorithm.isin(set(MAIN8) | set(ABL))]
    Sp = A[A.Algorithm.isin(SPM)]
    Sp = Sp[[(d, r, n) in SPLITCASES for d, r, n in zip(Sp.Dataset, Sp.Radius, Sp.Turbines)]]
    assert set(Sp.Algorithm) == set(SPM), f"split settings missing: {set(SPM) - set(Sp.Algorithm)}"
    assert "RSDVNS" in set(G.Algorithm), "RSD-VNS runs (mpce_rsdisc) missing"
    SP = case_stats(Sp, SPM, BASE_THR)
    S = case_stats(G, sorted(set(MAIN8) | set(ABL)), BASE_THR)
    cases = S[CASE].drop_duplicates()
    clusters = sorted({(d, int(r)) for d, r in zip(cases.Dataset, cases.Radius)})
    print(f"loaded {len(A)} runs; {len(cases)} cases in {len(clusters)} clusters {clusters}")

    summ = dict(generated=time.strftime("%Y-%m-%d %H:%M:%S"), script="mpce_inference_extra.py",
                n_cases=int(len(cases)), clusters=[f"{d}-{r}" for d, r in clusters],
                cluster_sizes={f"{d}-{r}": int(((cases.Dataset == d) & (cases.Radius == r)).sum()) for d, r in clusters},
                thresholds=THRESHOLDS, base_threshold=BASE_THR, boot=dict(n=BOOT_N, seed=BOOT_SEED),
                perm=dict(n=args.perm, seed=PERM_SEED), energy_factor_mwh_per_objective=ENERGY_FACTOR)
    # --- 0. reproduction of mpce_summary.json at the baseline threshold
    TH = threshold_block(G, Sp)
    base = TH[BASE_THR]
    rep = dict(ok=None)
    try:
        MS = json.load(open(os.path.join(HERE, "mpce_summary.json")))
        dr = max(abs(base["friedman"]["avg_rank"][a] - MS["main"]["friedman"]["avg_rank"][a]) for a in MAIN8)
        dp = max(abs(math.log10(base["main_case_mean"][b]["p"]) - math.log10(MS["main"]["case_mean_wilcoxon"][b]["p"]))
                 for b in base["main_case_mean"])
        assert set(base["ablation_case_mean"]) == set(MS["ablation"]["case_mean"]), "ablation contrast sets differ"
        dpa = max(max(abs(math.log10(base["ablation_case_mean"][k][f]) - math.log10(MS["ablation"]["case_mean"][k][f]))
                      for f in ("p", "p_holm")) for k in base["ablation_case_mean"])
        dra = max(abs(base["ablation_friedman"]["avg_rank"][a] - MS["ablation"]["friedman"]["avg_rank"][a])
                  for a in MS["ablation"]["variants"])
        assert set(MS["ablation"]["variants"]) == set(ABL), "ablation variants differ"
        sub = MS["main"]["subgroup_vs_phase1"]["groups"]
        dsub = max(abs(math.log10(cm_test(S[S.Dataset.isin(v["datasets"]) & (S.Turbines >= 10)], FOCUS, "PSOC")["p"])
                       - math.log10(v["p"])) for v in sub.values())
        spl = dict(MS["split"]["case_mean"], **(MS["split"].get("case_mean_omega90") or {}))
        assert set(spl) == {pk(a, b) for a, b in SPLIT_PAIRS}, f"split pairs differ: {sorted(spl)}"
        dsp = max(max(abs(math.log10(cm_test(SP, *k.split("-"))["p"]) - math.log10(v["p"])),
                      abs(cm_test(SP, *k.split("-"))["mean_dloss_pp"] - v["mean_dloss_pp"]))
                  for k, v in spl.items())
        dmd = max(abs(base["ablation_case_mean"][k]["mean_dloss_pp"] - MS["ablation"]["case_mean"][k]["mean_dloss_pp"])
                  for k in base["ablation_case_mean"])
        rep = dict(max_abs_diff_avg_rank=dr, max_abs_diff_avg_rank_ablation=dra, max_abs_diff_log10p_main=dp,
                   max_abs_diff_log10p_ablation=dpa, max_abs_diff_mean_dloss_ablation=dmd,
                   max_abs_diff_log10p_subgroups=dsub, max_abs_diff_log10p_split=dsp,
                   n_ablation_contrasts=len(base["ablation_case_mean"]), n_split_tests=len(spl),
                   summary_generated=MS.get("generated"), ok=bool(max(dr, dra, dp, dpa, dmd, dsub, dsp) < 1e-6))
    except Exception as e:                                               # pragma: no cover
        rep = dict(ok=False, error=repr(e))
    summ["reproduction"] = rep
    print("reproduction of mpce_summary.json:", rep)
    summ["threshold"] = {str(k): v for k, v in TH.items()}
    concl_names = list(base["conclusions"])
    summ["threshold_summary"] = dict(
        conclusions=concl_names,
        unchanged={c: bool(all(TH[t]["conclusions"][c] for t in THRESHOLDS)) for c in concl_names},
        all_unchanged=bool(all(TH[t]["all_unchanged"] for t in THRESHOLDS)),
        best_ranked={str(t): TH[t]["friedman"]["best_ranked"] for t in THRESHOLDS},
        rank_PSOBV={str(t): TH[t]["friedman"]["avg_rank"]["PSOBV"] for t in THRESHOLDS},
        rank_PSOC={str(t): TH[t]["friedman"]["avg_rank"]["PSOC"] for t in THRESHOLDS},
        p_PSOBV_PSOC={str(t): TH[t]["main_case_mean"]["PSOC"]["p"] for t in THRESHOLDS},
        order_identical=bool(all(TH[t]["friedman"]["order"] == base["friedman"]["order"] for t in THRESHOLDS)),
        orders={str(t): TH[t]["friedman"]["order"] for t in THRESHOLDS})
    for t in THRESHOLDS:
        print(f"thr {t:2d}: best {TH[t]['friedman']['best_ranked']}, order {TH[t]['friedman']['order']}, "
              f"conclusions {TH[t]['conclusions']}")

    # --- 2. clusters
    CL = cluster_block(S, clusters)
    full = {}
    for a, b in ALL_PAIRS:
        x = cm_test(S, a, b); x["verdict"] = verdict(x); full[pk(a, b)] = x
    LO = loco_block(S, clusters, full)
    summ["case_level"] = full
    summ["cluster"] = CL
    summ["loco"] = LO
    summ["cluster_min_attainable_p"] = float(2 / 2 ** len(clusters))
    # budget split (12 cases = 2 per cluster): the same cluster / LOCO analyses for the six split tests
    fullsp = {}
    for a, b in SPLIT_PAIRS:
        x = cm_test(SP, a, b); x["verdict"] = verdict(x); fullsp[pk(a, b)] = x
    CLS = cluster_block(SP, clusters, SPLIT_PAIRS)
    LOS = loco_block(SP, clusters, fullsp, SPLIT_PAIRS, rank=False)
    summ["split_case_level"] = fullsp
    summ["cluster_split"] = CLS
    summ["loco_split"] = LOS
    for k, v in CLS.items():
        print(f"split cluster {k:16s} fav {v['clusters_favour_first']}/{v['clusters_favour_second']} wil p {v['wilcoxon_exact_p']:.3f} "
              f"CR1 p {v['cr1_p']:.2g} cboot {np.round(v['cluster_boot_ci95'], 3)} | case p {fullsp[k]['p']:.3g} | LOCO p "
              f"{LOS['pairs'][k]['p_min']:.2g}-{LOS['pairs'][k]['p_max']:.2g} changes {LOS['pairs'][k]['verdict_changes']}")
    for k, v in CL.items():
        print(f"cluster {k:14s} fav {v['clusters_favour_first']}/{v['clusters_favour_second']} sign p {v['sign_test_p']:.3f} "
              f"wil p {v['wilcoxon_exact_p']:.3f} flip p {v['cluster_signflip_p']:.3f} CR1 p {v['cr1_p']:.2g} "
              f"cboot {np.round(v['cluster_boot_ci95'], 3)} | LOCO p {LO['pairs'][k]['p_min']:.2g}-{LO['pairs'][k]['p_max']:.2g} "
              f"verdict changes {LO['pairs'][k]['verdict_changes']} flips {LO['pairs'][k]['direction_flips']}")
    # --- 3. multiplicity
    MU = multiplicity_block(S, G, SP)
    summ["multiplicity"] = MU
    for t in MU["tests"]:
        print(f"  {t['key']:28s} n={t['n']:2d} p={t['p']:.3g} holm={t['p_holm_all']:.3g} bh={t['q_bh_all']:.3g}")
    print("lost under Holm:", MU["lost_under_holm"], "lost under BH:", MU["lost_under_bh"])
    # --- 4. subgroup
    SG = subgroup_block(S, clusters, args.perm)
    summ["subgroup"] = SG
    print(json.dumps({k: v for k, v in SG.items()}, indent=1, default=float)[:3000])
    # --- 5. energy
    EN = energy_block(S)
    summ["energy"] = EN
    # --- 6. wake-model uncertainty
    MSH = model_shift_block(S, args.data_dir)
    summ["model_shift"] = MSH
    print(f"model shift: mean |L_G - L_J| {MSH['mean_abs_shift_pp']:.3f} pp, median {MSH['median_abs_shift_pp']:.3f} pp "
          f"({MSH['n_layouts']} layouts; ratio to margin {MSH['ratio_mean_to_margin']:.1f}); PSO-VNS - PSO pair change "
          f"{MSH['pairs']['PSOBV-PSOC']['mean_abs_change_pp']:.4f} pp (J {MSH['pairs']['PSOBV-PSOC']['mean_diff_jensen_pp']:+.4f}, "
          f"G {MSH['pairs']['PSOBV-PSOC']['mean_diff_gauss_pp']:+.4f}); recorded vs paper {MSH['pair_recorded_vs_paper_abs_diff']:.2g}")
    for k, v in EN.items():
        print(f"energy {k}: dL {v['mean_dloss_pp']:+.4f} pp, gain {v['mean_gain_pct_aep']:+.4f} % AEP, "
              f"{v['mean_gain_mwh_yr']:+.1f} MWh/yr (CI {np.round(v['ci95_gain_mwh_yr'], 1)}), largest N {v['mean_gain_mwh_yr_largest_N']:+.1f}")
    # --- 7-11. review round 2 (D13): equivalence at three levels, heterogeneity, HL, McNemar, equivalence curve
    EL = equivalence_levels_block(G, S, clusters)
    EL["reproduction"] = reproduce_equivalence(EL)
    print("reproduction of mpce_summary.json['equivalence']:", EL["reproduction"])
    summ["equivalence_levels"] = EL
    summ["heterogeneity"] = heterogeneity_block(S, clusters, args.perm, EL)
    summ["hodges_lehmann"] = hl_block(S, full)
    summ["hr_mcnemar"] = mcnemar_block(args.data_dir)
    os.makedirs(args.fig_dir, exist_ok=True)
    equiv_figure(EL, os.path.join(args.fig_dir, "equiv_curve.pdf"))
    summ["figure_equiv_curve"] = os.path.join(os.path.basename(args.fig_dir), "equiv_curve.pdf")
    summ["seconds"] = round(time.time() - t0, 1)

    json.dump(clean(summ), open(os.path.join(args.out_dir, "mpce_summary_extra.json"), "w"), indent=1)
    write_macros(summ, os.path.join(args.out_dir, "mpce_numbers_extra.tex"))
    write_tables(summ, os.path.join(args.out_dir, "mpce_supp_inference.tex"))
    print(f"done in {time.time() - t0:.1f} s")


def clean(o):
    if isinstance(o, dict):
        return {str(k): clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [clean(v) for v in o]
    if isinstance(o, (np.bool_,)):
        return bool(o)
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating, float)):
        return None if not np.isfinite(o) else float(o)
    return o


# ------------------------------------------------------------------ macros
def write_macros(s, fn):
    M = []

    def m(name, val, note=""):
        assert re.fullmatch(r"NX[A-Za-z]+", name), name
        M.append((name, val, note))

    TH = s["threshold"]; ts = s["threshold_summary"]; T = [str(t) for t in s["thresholds"]]
    m("NXThrList", listing([str(t) for t in s["thresholds"]]), "qualification thresholds (feasible runs of 30)")
    m("NXThrMin", str(min(s["thresholds"]))); m("NXThrMax", str(max(s["thresholds"])))
    bset = sorted(set(ts["best_ranked"].values()))
    m("NXThrBest", " / ".join(LAB[b] for b in bset), "best-ranked method(s) over all thresholds")
    m("NXThrRankPSOVNSMin", num(min(ts["rank_PSOBV"].values()))); m("NXThrRankPSOVNSMax", num(max(ts["rank_PSOBV"].values())))
    m("NXThrRankPSOMin", num(min(ts["rank_PSOC"].values()))); m("NXThrRankPSOMax", num(max(ts["rank_PSOC"].values())))
    m("NXThrOrderSame", "yes" if ts["order_identical"] else "no", "same 8-method rank order at every threshold")
    m("NXThrPSOPMin", pval_down(min(ts["p_PSOBV_PSOC"].values())), "PSO-VNS vs PSO case-mean p, smallest over thresholds (rounded down)")
    m("NXThrPSOPMax", pval(max(ts["p_PSOBV_PSOC"].values())))
    for (a, b) in KEY_CONTR[1:]:
        ps = [TH[t]["ablation_case_mean"][pk(a, b)]["p"] for t in T]
        if EXPECT[(a, b)] == "ns" or (a, b) in RSD_CONTR:
            m(f"NXThr{pmac(a, b)}PMin", pval_down(min(ps)), f"{a}-{b}: smallest p over thresholds (rounded down)")
        if EXPECT[(a, b)] != "ns" or (a, b) in RSD_CONTR:
            m(f"NXThr{pmac(a, b)}PMax", pval_up(max(ps)), f"{a}-{b}: largest p over thresholds (rounded up)")
    for (a, b) in SPLIT_PAIRS:                                        # budget split (12 cases), unadjusted
        ps = [TH[t]["split_case_mean"][pk(a, b)]["p"] for t in T]
        m(f"NXThr{SPLIT_NAME[(a, b)]}PMin", pval_down(min(ps)), f"split {a}-{b}: smallest p over thresholds (rounded down)")
        m(f"NXThr{SPLIT_NAME[(a, b)]}PMax", pval_up(max(ps)), f"split {a}-{b}: largest p over thresholds (rounded up)")
    nconc = len(ts["conclusions"]); nun = sum(ts["unchanged"].values())
    m("NXThrNConcl", WORD.get(nconc, str(nconc))); m("NXThrNUnchanged", WORD.get(nun, str(nun)))
    changed = [c for c, v in ts["unchanged"].items() if not v]
    m("NXThrChanged", listing(changed).replace("_", " ") if changed else "none")
    m("NXThrChangedAt", listing(sorted({t for t in T if not TH[t]["all_unchanged"]}, key=int)) or "none",
      "thresholds at which at least one conclusion changes")
    imp = [TH[t]["imputed_main"] for t in T]
    m("NXThrImputedMin", str(min(imp))); m("NXThrImputedMax", str(max(imp)))
    # clusters
    CL = s["cluster"]; LO = s["loco"]
    m("NXClusters", WORD[len(s["clusters"])], "number of (data set, radius) clusters")
    m("NXClustersNum", str(len(s["clusters"])))
    m("NXClMinP", pval(s["cluster_min_attainable_p"]), "smallest attainable two-sided exact p over the clusters")
    m("NXBootN", "10{,}000")
    pairs_cl = [(pk(a, b), pmac(a, b), CL, LO) for a, b in ALL_PAIRS] + \
               [(pk(a, b), SPLIT_NAME[(a, b)], s["cluster_split"], s["loco_split"]) for a, b in SPLIT_PAIRS]
    for k, nm, CLx, LOx in pairs_cl:
        c = CLx[k]; lo = LOx["pairs"][k]; a, b = k.split("-")
        m(f"NXCl{nm}Fav", str(c["clusters_favour_first"]), f"clusters in which {a} has the lower mean loss")
        m(f"NXCl{nm}Opp", str(c["clusters_favour_second"]))
        m(f"NXCl{nm}SignP", pval(c["sign_test_p"])); m(f"NXCl{nm}WilP", pval(c["wilcoxon_exact_p"]))
        m(f"NXCl{nm}FlipP", pval(c["cluster_signflip_p"])); m(f"NXCl{nm}CRP", pval(c["cr1_p"]))
        m(f"NXCl{nm}CI", ci_txt(c["cluster_boot_ci95"]), "cluster-bootstrap 95% CI of the mean case difference (pp)")
        m(f"NXCl{nm}CRCI", ci_txt(c["cr1_ci95"]), "cluster-robust t 95% CI (pp)")
        m(f"NXLoco{nm}PMin", pval_down(lo["p_min"])); m(f"NXLoco{nm}PMax", pval_up(lo["p_max"]))
        m(f"NXLoco{nm}Changes", str(lo["verdict_changes"]), "LOCO runs (of 6) in which the significance verdict changes")
        m(f"NXLoco{nm}Flips", str(lo["direction_flips"]), "LOCO runs in which the sign of the mean difference flips")
    m("NXLocoBestAlways", "yes" if LO["best_ranked_always_PSOBV"] else "no")
    nchg = [k for k, v in LO["pairs"].items() if v["verdict_changes"] > 0]
    m("NXLocoNChanged", WORD.get(len(nchg), str(len(nchg))), "pairs whose verdict changes in at least one LOCO run")
    m("NXLocoChanged", listing([plab(*k.split("-")) for k in nchg]) if nchg else "none")
    m("NXLocoNPairs", WORD.get(len(ALL_PAIRS), str(len(ALL_PAIRS))))
    nfl = [k for k, v in LO["pairs"].items() if v["direction_flips"] > 0]
    m("NXLocoFlipped", listing([plab(*k.split("-")) for k in nfl]) if nfl else "none")
    csig = [k for k, v in s["case_level"].items() if v["p"] < 0.05]
    m("NXClNCaseSig", str(len(csig)), "68-case pairs significant on the case means (unadjusted)")
    cun = [k for k in csig if min(CL[k]["clusters_favour_first"], CL[k]["clusters_favour_second"]) == 0]
    m("NXClNCaseSigUnanimous", str(len(cun)), "of these: pairs in which all clusters favour the same method")
    cnu = [k for k in csig if k not in cun]
    m("NXClCaseSigNotUnanimous", listing([plab(*k.split("-")) for k in cnu]) if cnu else "none")
    nsig_cl = [k for k, v in CL.items() if v["wilcoxon_exact_p"] < 0.05]
    m("NXClNSigWil", str(len(nsig_cl)), "pairs significant at the cluster level (exact Wilcoxon over clusters)")
    nall = [k for k, v in CL.items() if min(v["clusters_favour_first"], v["clusters_favour_second"]) == 0]
    m("NXClNUnanimous", str(len(nall)), "pairs in which every cluster favours the same method")
    ncb = [k for k, v in CL.items() if v["cluster_boot_excludes_zero"]]
    m("NXClNBootExcl", str(len(ncb)), "pairs whose cluster-bootstrap CI excludes 0")
    ncr = [k for k, v in CL.items() if v["cr1_p"] < 0.05]
    m("NXClNCRSig", str(len(ncr)), "pairs significant with the cluster-robust t test")
    # multiplicity
    MU = s["multiplicity"]
    m("NXMultN", str(MU["n_tests"])); m("NXMultSigRaw", str(MU["n_sig_raw"]))
    m("NXMultSigHolm", str(MU["n_sig_holm"])); m("NXMultSigBH", str(MU["n_sig_bh"]))
    fams = [t["key"].split(":")[0] for t in MU["tests"]]
    for f_, nm in (("main", "Main"), ("abl", "Abl"), ("sub", "Sub"), ("split", "Split")):
        m(f"NXMultN{nm}", str(fams.count(f_)), f"tests of kind '{f_}' in the all-family Holm")
    m("NXMultNWord", WORD.get(MU["n_tests"], str(MU["n_tests"])))
    lab = {t["key"]: t["label"] for t in MU["tests"]}
    m("NXMultLostHolm", listing([lab[k] for k in MU["lost_under_holm"]]) if MU["lost_under_holm"] else "none")
    m("NXMultLostBH", listing([lab[k] for k in MU["lost_under_bh"]]) if MU["lost_under_bh"] else "none")
    m("NXMultMaxSigHolmP", pval_up(MU["max_sig_p_holm"]), "largest all-family Holm p among the tests still significant (rounded up)")
    tk = {t["key"]: t for t in MU["tests"]}
    for key, nm in (("abl:SSABV-RSVNS", "SSAVNSvsRSVNS"), ("abl:SSABV-LXBV", "SSAVNSvsLXSSAVNS"), ("abl:LXBV-RSVNS", "LXSSAVNSvsRSVNS"),
                    ("abl:PSOBV-RSVNS", "PSOVNSvsRSVNS"), ("main:PSOBV-PSOC", "PSOVNSvsPSO"), ("sub:dsIILarge", "PSOVNSvsPSOdsIILarge"),
                    ("sub:Large", "PSOVNSvsPSOLarge"), ("sub:dsILarge", "PSOVNSvsPSOdsILarge"), ("split:PSOBV-PSOBV75", "SplitFiftyVsSeventyFive"),
                    ("split:PSOBV-PSOBV25", "SplitFiftyVsTwentyFive"), ("split:PSOBV75-PSOC", "SplitSeventyFiveVsHundred"),
                    ("split:PSOBV-PSOC", "SplitFiftyVsHundred"),
                    # Phase 6 controls: RSD-VNS contrasts and omega = 0.9 split tests
                    ("abl:SSABV-RSDVNS", "SSAVNSvsRSDVNS"), ("abl:LXBV-RSDVNS", "LXSSAVNSvsRSDVNS"),
                    ("abl:PSOBV-RSDVNS", "PSOVNSvsRSDVNS"), ("abl:RSDVNS-RSVNS", "RSDVNSvsRSVNS"),
                    ("split:PSOBV-PSOBV90", "SplitFiftyVsNinety"), ("split:PSOBV90-PSOBV75", "SplitNinetyVsSeventyFive")):
        m(f"NXHolm{nm}P", pval_up(tk[key]["p_holm_all"]) if tk[key]["p_holm_all"] < 0.05 else pval(tk[key]["p_holm_all"]),
          "all-family Holm p (rounded up when significant)")
        m(f"NXBH{nm}Q", pval_up(tk[key]["q_bh_all"]) if tk[key]["q_bh_all"] < 0.05 else pval(tk[key]["q_bh_all"]))
    # subgroup
    SG = s["subgroup"]
    ia, i1, i2 = SG["interaction_all"], SG["interaction_dsI"], SG["interaction_dsII"]
    m("NXPermN", f"{s['perm']['n']:,}".replace(",", "{,}"))
    m("NXIntDiff", num(ia["diff_pp"], 3), "mean PSO-VNS - PSO (pp) for N>=10 minus that for N<10 (68 cases)")
    m("NXIntP", pval(ia["p_stratified"]), "permutation p, labels permuted within clusters")
    m("NXIntPUnstr", pval(ia["p_unrestricted"]))
    m("NXIntDsIIDiff", num(i2["diff_pp"], 3)); m("NXIntDsIIP", pval(i2["p_stratified"])); m("NXIntDsIIPUnstr", pval(i2["p_unrestricted"]))
    m("NXIntDsIDiff", num(i1["diff_pp"], 3)); m("NXIntDsIP", pval(i1["p_stratified"])); m("NXIntDsIPUnstr", pval(i1["p_unrestricted"]))
    for tag, nm in (("all", "All"), ("dsI", "DsI"), ("dsII", "DsII")):
        v = SG[f"spearman_{tag}"]
        m(f"NXSpear{nm}", num(v["rho"], 2)); m(f"NXSpear{nm}P", pval(v["p"]))
        m(f"NXSpear{nm}Strat", num(v["stratified_mean_rho"], 2), "mean within-cluster Spearman rho(N, PSO-VNS - PSO)")
        m(f"NXSpear{nm}StratP", pval(v["stratified_p"]), "within-cluster permutation p")
    for tag, nm in (("large", "Large"), ("all", "All"), ("small", "Small")):
        v = SG[f"dsII_vs_dsI_{tag}"]
        m(f"NXDsPair{nm}N", str(v["n_pairs"])); m(f"NXDsPair{nm}Mean", num(v["mean_pp"], 3)); m(f"NXDsPair{nm}P", pval(v["p"]))
        m(f"NXDsPair{nm}IIWins", str(v["dsII_larger_gain"]), "(r,N) pairs in which the PSO-VNS advantage is larger in Data Set II")
    m("NXSubgroupSurvives", "yes" if SG["pattern_survives"] else "no")
    # energy
    EN = s["energy"]
    for k, nm in (("PSOBV-PSOC", "PSOVNSvsPSO"), ("SSABV-RSVNS", "SSAVNSvsRSVNS"), ("PSOBV-RSVNS", "PSOVNSvsRSVNS"),
                  ("PSOBV-PSOC_dsIILarge", "PSOVNSvsPSOdsIILarge")):
        v = EN[k]
        m(f"NXEn{nm}Pct", num(v["mean_gain_pct_aep"], 3), "mean gain of the first method, % of the second's benchmark AEP")
        m(f"NXEn{nm}PctCI", ci_txt(v["ci95_gain_pct_aep"], 3))
        m(f"NXEn{nm}PctMax", num(v["max_gain_pct_aep"], 2)); m(f"NXEn{nm}PctMin", num(v["min_gain_pct_aep"], 2))
        m(f"NXEn{nm}MWh", num(v["mean_gain_mwh_yr"], 0), "mean gain, MWh/yr of benchmark (model) energy per case")
        m(f"NXEn{nm}MWhCI", ci_txt(v["ci95_gain_mwh_yr"], 0))
        m(f"NXEn{nm}MWhTurb", num(v["mean_gain_mwh_yr_per_turbine"], 1), "mean gain per turbine, MWh/yr")
        m(f"NXEn{nm}MWhLargest", num(v["mean_gain_mwh_yr_largest_N"], 0), f"mean gain at N = {v['largest_N']}")
        m(f"NXEn{nm}PctLargest", num(v["mean_gain_pct_aep_largest_N"], 3))
    m("NXEnFactor", "8.76/15", "MWh/yr per unit of benchmark objective")
    # wake-model uncertainty (Jensen vs Gaussian re-evaluation of the same final layouts)
    MS = s["model_shift"]; mp = MS["pairs"]["PSOBV-PSOC"]
    m("NXModelShiftMean", num(MS["mean_abs_shift_pp"], 2), "mean |L_Gauss - L_Jensen| of the same feasible final layouts (pp)")
    m("NXModelShiftMedian", num(MS["median_abs_shift_pp"], 2), "median |L_Gauss - L_Jensen| (pp)")
    m("NXModelShiftPTen", num(MS["q10_abs_shift_pp"], 2)); m("NXModelShiftPNinety", num(MS["q90_abs_shift_pp"], 2))
    m("NXModelShiftSigned", num(MS["mean_signed_shift_pp"], 2, sign=True), "mean L_Gauss - L_Jensen (pp)")
    m("NXModelShiftDsI", num(MS["by_dataset"]["1"]["mean_abs_shift_pp"], 2)); m("NXModelShiftDsII", num(MS["by_dataset"]["2"]["mean_abs_shift_pp"], 2))
    m("NXModelShiftNLayouts", f"{MS['n_layouts']:,}".replace(",", "{,}"), "feasible final layouts re-evaluated (8 main methods)")
    m("NXModelShiftRatio", num(MS["ratio_mean_to_margin"], 1), "mean model shift / equivalence margin")
    m("NXModelShiftMedianRatio", num(MS["ratio_median_to_margin"], 1), "median model shift / equivalence margin")
    m("NXModelShiftAboveMarginPct", num(100 * MS["share_abs_shift_above_margin"], 0), "% of layouts whose shift exceeds the margin")
    m("NXModelShiftPairMean", num(mp["mean_abs_change_pp"], 3), "mean |change| of the per-case PSO-VNS - PSO difference, Gauss vs Jensen (pp)")
    m("NXModelShiftPairMax", num(mp["max_abs_change_pp"], 3))
    m("NXModelShiftPairChangeOfMean", num(mp["abs_change_of_mean_pp"], 3), "|mean_Gauss - mean_Jensen| of PSO-VNS - PSO (pp)")
    m("NXModelShiftPairJensen", num(mp["mean_diff_jensen_pp"], 3), "mean PSO-VNS - PSO (pp), Jensen, re-evaluated coordinates")
    m("NXModelShiftPairGauss", num(mp["mean_diff_gauss_pp"], 3), "mean PSO-VNS - PSO (pp), Gaussian wake")
    m("NXModelShiftPairSignChanges", str(mp["sign_changes"]), "cases in which the sign of PSO-VNS - PSO differs between the models")
    m("NXModelShiftPairN", str(mp["n_cases"]))
    write_macros_r2(s, m)
    names = [x[0] for x in M]
    assert len(names) == len(set(names)), f"duplicate macros: {sorted({x for x in names if names.count(x) > 1})}"
    hdr = ["% generated by mpce_inference_extra.py from the per-run data (" + s["generated"] + ") -- do not edit by hand",
           "% Phase 6, W2 (inference robustness): qualification threshold, cluster-level / leave-one-cluster-out / cluster",
           "% bootstrap, one Holm (and BH) family over all main-text case-mean tests, the post hoc N >= 10 subgroup, and",
           "% effect sizes in benchmark energy. Supplementary tables: mpce_supp_inference.tex (tab:X-...);",
           "% checks: mpce_check_extra.py (X01...). Signs: pair macros 'AvsB' use first minus second (pp of wake loss,",
           "% negative = first better) except the \\NXEn... macros, which give the GAIN of the first method (positive = more",
           "% energy). p values: 'PMax' / 'MaxSig' bounds are rounded up, 'PMin' bounds rounded down.",
           "% reproduction of mpce_summary.json at the 15-run threshold: " + ("ok" if s["reproduction"].get("ok") else "FAILED")]
    lines = hdr + [f"\\newcommand{{\\{n}}}{{{v}}}" + (f"   % {c}" if c else "") for n, v, c in M]
    open(fn, "w").write("\n".join(lines) + "\n")
    print(f"wrote {fn} ({len(M)} macros)")


def up3(v):
    """minimal margins: rounded UP to 3 decimals (as the pipeline's m_min column)."""
    return f"{math.ceil(v * 1000 - 1e-9) / 1000:.3f}"


def ptost(r):
    return ("\\ensuremath{\\le}" + pval(r["p_tost"])) if r.get("p_tost_at_floor") else pval(r["p_tost"])


def bprob(v):
    """Bayesian 'most probable region' probabilities: never 1.00 / 0.00 (as the pipeline)."""
    return "\\ensuremath{>}0.99" if v > 0.99 else "\\ensuremath{<}0.01" if v < 0.01 else f"{v:.2f}"


CL_NAME = {"1-500": "DsIFive", "1-750": "DsISevenFifty", "1-1000": "DsIThousand",
           "2-500": "DsIIFive", "2-750": "DsIISevenFifty", "2-1000": "DsIIThousand"}


def case_list(cases):
    """'500~m: $N=8$ ($+0.160$), ...' grouped by radius (for the \\NXHet...List macros)."""
    out = []
    for rr in sorted({c["radius"] for c in cases}):
        cc = [c for c in cases if c["radius"] == rr]
        out.append(f"{rr}~m: " + ", ".join(f"$N={c['N']}$ (${c['d']:+.2f}$)" for c in cc))
    return "; ".join(out) if out else "none"


def write_macros_r2(s, m):
    """macros of the review-round-2 analyses (sections 7-10): \\NXEq..., \\NXBayMean..., \\NXHet..., \\NXDens...,
    \\NXHL..., \\NXHRMcNemar..., \\NXMargin.... Signs: first minus second (pp of wake loss, negative = first better).
    Bayesian 'Left/Rope/Right': theta_A (first better) / theta_rope / theta_B (second better); the \\NX...Bay{Left,Rope,
    Right} probabilities are the posterior probabilities that each region is the MOST PROBABLE one (not the
    probability that the mean difference lies in the rope); \\NXBayMean... are the posterior means of theta."""
    EL = s["equivalence_levels"]; P_ = EL["pairs"]
    m("NXEqMarginSource", EL["margin_source"].replace("_", "\\_"), "where EQ_MARGIN was read")
    for key, r in P_.items():
        a, b = key.split("-"); nm = pmac(a, b); sd, cs, cl = r["seed"], r["case"], r["cluster"]
        m(f"NXEqMean{nm}", num(r["mean_dloss_pp"], 3), f"{key}: mean case difference (pp), n = {r['n_cases']}")
        m(f"NXEqSeed{nm}CI", ci_txt(sd["ci90"]), "seed level (fixed benchmark): 90% CI, runs resampled within cases")
        m(f"NXEqSeed{nm}CINinetyFive", ci_txt(sd["ci95"]))
        m(f"NXEqSeed{nm}Min", up3(sd["min_margin_pp"]), "minimal equivalence margin (pp, rounded up)")
        m(f"NXEqSeed{nm}P", ptost(sd), "bootstrap TOST p at the margin")
        m(f"NXEqSeed{nm}Holds", "yes" if sd["equivalent"] else "no")
        m(f"NXEqCase{nm}CI", ci_txt(cs["ci90"]), "case level (bootstrap over cases, = pipeline): 90% CI")
        m(f"NXEqCase{nm}CINinetyFive", ci_txt(cs["ci95"]))
        m(f"NXEqCase{nm}Min", up3(cs["min_margin_pp"])); m(f"NXEqCase{nm}P", ptost(cs))
        m(f"NXEqCase{nm}Holds", "yes" if cs["equivalent"] else "no")
        m(f"NXEqClust{nm}CI", ci_txt(cl["cr2"]["ci90"]), "cluster level: CR2 90% CI, t with 5 d.f.")
        m(f"NXEqClust{nm}CINinetyFive", ci_txt(cl["cr2"]["ci95"]))
        m(f"NXEqClust{nm}Min", up3(cl["cr2"]["min_margin_pp"]), "minimal margin, CR2 t(5)")
        m(f"NXEqClust{nm}PCR", pval(cl["cr2"]["p_tost"]), "CR2 t(5) TOST p")
        m(f"NXEqClust{nm}P", pval(cl["wild"]["p_tost"]), "wild-cluster bootstrap-t TOST p (Webb, all draws)")
        m(f"NXEqClust{nm}CIWild", ci_txt(cl["wild"]["ci90"]), "wild-cluster 90% CI (test inversion)")
        m(f"NXEqClust{nm}MinWild", up3(cl["wild"]["min_margin_pp"]))
        m(f"NXEqClust{nm}Holds", "yes" if cl["equivalent"] else "no", "CR2 CI inside +-m AND wild TOST p <= 0.05")
        th = r["bayes"]["mean_theta"]
        for k, v in zip(("Left", "Rope", "Right"), th):
            m(f"NXBayMean{nm}{k}", f"{v:.2f}", "posterior mean of theta (Left = first better)")
    f = P_["PSOBV-PSOC"]
    if "per_seed_t" in f["seed"]:
        m("NXEqSeedJointPSOVNSvsPSOCI", ci_txt(f["seed"]["per_seed_t"]["ci90"]),
          "benchmark-average difference per seed (30 seeds), t(29) 90% CI (the reviewer's version)")
    m("NXEqSeedStratJointPSOVNSvsPSOCI", ci_txt(f["seed"]["joint"]["ci90"]), "seed bootstrap with one seed resample for all cases")
    m("NXEqClustCROnePSOVNSvsPSOCI", ci_txt(f["cluster"]["cr1"]["ci90"]), "CR1 t(5) 90% CI")
    m("NXEqClustDfBM", num(f["cluster"]["df_bm"], 1), "Bell-McCaffrey d.f. of the CR2 variance (PSO-VNS - PSO)")
    m("NXEqClustBMPSOVNSvsPSOCI", ci_txt(f["cluster"]["cr2_bm"]["ci90"]), "CR2 90% CI with Bell-McCaffrey d.f.")
    m("NXEqClustDf", str(f["cluster"]["df"]))
    m("NXEqWildDraws", f"{f['cluster']['wild']['draws']:,}".replace(",", "{,}"), "Webb draws enumerated (6^6)")
    m("NXBayPSOVNSvsPSOCIRope", ci_txt(f["bayes"]["ci95_theta"][1], 2), "95% credible interval of theta_rope")
    m("NXBayPSOVNSvsPSOCILeft", ci_txt(f["bayes"]["ci95_theta"][0], 2))
    nhold = {lv: sum(1 for r in P_.values() if (r[lv]["equivalent"])) for lv in ("seed", "case", "cluster")}
    for lv, nm in (("seed", "Seed"), ("case", "Case"), ("cluster", "Clust")):
        m(f"NXEq{nm}NHolds", str(nhold[lv]), f"tab:equivalence pairs equivalent at the {lv} level")
    # heterogeneity of PSO-VNS - PSO
    H = s["heterogeneity"]
    m("NXHetN", str(H["n_cases"]))
    m("NXHetBeyondPSOVNS", str(H["beyond_first"]), "cases with PSO-VNS - PSO < -margin (PSO-VNS better by more than m)")
    m("NXHetBeyondPSO", str(H["beyond_second"]), "cases with PSO-VNS - PSO > +margin")
    m("NXHetBeyond", str(H["beyond_first"] + H["beyond_second"]))
    m("NXHetMaxPSOVNS", num(H["max_first"], 2), "largest advantage of PSO-VNS in one case (pp)")
    m("NXHetMaxPSO", num(H["max_second"], 2), "largest advantage of PSO in one case (pp)")
    m("NXHetTies", str(H["exact_ties"]), "cases with identical case means"); m("NXHetNearTies", str(H["near_ties"]), "|d| < 0.005 pp")
    for tag, nm in (("small", "Small"), ("large", "Large"), ("trivial", "Trivial"), ("nontrivial", "Nontrivial")):
        v = H[tag]; y = v["bayes"]
        m(f"NXEq{nm}N", str(v["n_cases"]), f"PSO-VNS vs PSO, stratum {tag}")
        m(f"NXEq{nm}Mean", num(v["mean_dloss_pp"], 3)); m(f"NXEq{nm}CI", ci_txt(v["ci90"]), "90% case-bootstrap CI")
        m(f"NXEq{nm}Min", up3(v["min_margin_pp"])); m(f"NXEq{nm}P", ptost(v)); m(f"NXEq{nm}Holds", "yes" if v["equivalent"] else "no")
        m(f"NXEq{nm}MeanAbs", num(v["mean_abs"], 3))
        m(f"NXEq{nm}BeyondPSOVNS", str(v["beyond_first"])); m(f"NXEq{nm}BeyondPSO", str(v["beyond_second"]))
        for k, pv, th in (("Left", y["p_a_better"], y["mean_theta"][0]), ("Rope", y["p_rope"], y["mean_theta"][1]),
                          ("Right", y["p_b_better"], y["mean_theta"][2])):
            m(f"NXEq{nm}Bay{k}", bprob(pv), "P(region most probable); Left = PSO-VNS better")
            m(f"NXEq{nm}BayMean{k}", f"{th:.2f}", "posterior mean of theta")
    m("NXTrivialLoss", f"{H['trivial_loss_pp']:g}", "trivial case: case-mean wake loss of PSO-VNS below this (pp)")
    cm = H["cluster_means"]
    for k, v in cm.items():
        m(f"NXHetCl{CL_NAME[k]}", num(v, 3), f"cluster {k}: mean PSO-VNS - PSO (pp)")
    m("NXHetClBeyond", str(len(H["clusters_beyond_first"]) + len(H["clusters_beyond_second"])), "cluster means beyond +-margin")
    m("NXHetClBeyondPSOVNS", str(len(H["clusters_beyond_first"]))); m("NXHetClBeyondPSO", str(len(H["clusters_beyond_second"])))
    for dsv, dn in (("1", "DsI"), ("2", "DsII")):
        cb = [c for c in H["cases_beyond"] if c["dataset"] == dsv]
        m(f"NXHet{dn}PSOList", case_list([c for c in cb if c["d"] > 0]), "cases beyond +margin (PSO better)")
        m(f"NXHet{dn}PSOVNSList", case_list([c for c in cb if c["d"] < 0]), "cases beyond -margin (PSO-VNS better)")
        m(f"NXHet{dn}NPSO", str(sum(c["d"] > 0 for c in cb))); m(f"NXHet{dn}NPSOVNS", str(sum(c["d"] < 0 for c in cb)))
    it = H["interaction"]
    m("NXDensDef", "\\ensuremath{N\\,(\\ell_{\\min}/2r)^2}", "density: nominal share of the farm disc covered by exclusion discs")
    m("NXDensLmin", f"{SMIN_M:.0f}"); m("NXDensPhiMin", num(it["phi_range"][0], 2)); m("NXDensPhiMax", num(it["phi_range"][1], 2))
    m("NXDensNPairs", str(it["n_pairs"]))
    m("NXDensInterSlope", num(it["slope_pp_per_unit_phi"], 2), "slope of d_II - d_I on phi (pp per unit phi)")
    m("NXDensInterP", pval(it["p_signflip"]), "studentized sign-flip permutation p (two-sided)")
    m("NXDensPerm", f"{it['perm']:,}".replace(",", "{,}"))
    for tag, dn in (("dsI", "DsI"), ("dsII", "DsII")):
        v = H[f"slope_{tag}"]
        m(f"NXDensSlope{dn}", num(v["slope_pp_per_unit_phi"], 2, sign=True), "slope of PSO-VNS - PSO on phi, radius FE (pp per unit phi)")
        m(f"NXDensP{dn}", pval(v["p_perm"]), "studentized within-radius permutation p")
        m(f"NXDensTopMean{dn}", num(v["mean_top_density"], 3, sign=True), "mean d over the densest quarter of the cases")
    for dsn, dn in (("dsI", "DsI"), ("dsII", "DsII")):
        v = H["margin_energy"][dsn]
        m(f"NXMarginMWhTurb{dn}", num(v["margin_mwh_yr_per_turbine"], 1), "margin as MWh/yr per turbine (benchmark energy)")
        m(f"NXWakeFreeMWhTurb{dn}", num(v["wakefree_mwh_yr_per_turbine"], 0), "wake-free benchmark AEP per turbine, MWh/yr")
    # Hodges-Lehmann
    HLr = s["hodges_lehmann"]
    for a, b in HL_KEY:
        v = HLr[pk(a, b)]
        m(f"NXHL{pmac(a, b)}", num(v["hl"], 3, sign=True), "Hodges-Lehmann estimate (median of Walsh averages, pp)")
        m(f"NXHL{pmac(a, b)}CI", ci_txt(v["ci95"]), "95% CI, exact inversion of the Wilcoxon signed-rank test")
    # McNemar
    Mc = s["hr_mcnemar"]
    m("NXHRMcNemarP", pval(Mc["p_exact"]), "exact McNemar p, Horns Rev feasibility PSO-VNS vs PSO (6,030, random starts)")
    m("NXHRMcNemarOnlyPSOVNS", str(Mc["only_first"])); m("NXHRMcNemarOnlyPSO", str(Mc["only_second"]))
    m("NXHRMcNemarBoth", str(Mc["both"])); m("NXHRMcNemarNeither", str(Mc["neither"]))
    m("NXHRMcNemarN", str(Mc["n_seeds"]))


# ------------------------------------------------------------------ tables
def write_tables(s, fn):
    T = []
    TH = s["threshold"]
    cn = s["threshold_summary"]["conclusions"]
    # --- threshold table
    lines = []
    for t in s["thresholds"]:
        x = TH[str(t)]; fr = x["friedman"]; mc = x["main_case_mean"]; ab = x["ablation_case_mean"]
        cells = [f"{t}" + (" (paper)" if t == s["base_threshold"] else ""), str(x["imputed_main"]),
                 LAB[fr["best_ranked"]], f"{fr['avg_rank']['PSOBV']:.2f}", f"{fr['avg_rank']['PSOC']:.2f}",
                 ", ".join(LAB[a] for a in fr["posthoc_nonsig"]) or "--",
                 tp(mc["PSOC"]["p"])] + [tp(ab[pk(a, b)]["p"]) for a, b in KEY_CONTR[1:5]] + \
                ["yes" if x["all_unchanged"] else "\\textbf{no}"]
        lines.append(" & ".join(cells) + " \\\\")
    T.append(table("table*", "Sensitivity of the case-level conclusions to the qualification threshold (a method is ranked in a case by the mean objective of its feasible runs only if at least this many of its 30 runs are feasible; the paper uses 15). Unq.: case--method cells of the eight-method comparison that do not qualify (PSO-VNS always qualifies, so each of them enters a case-mean test of PSO-VNS with an imputed maximal difference). Best-ranked method and average ranks of the Friedman ranking of Table~\\ref{M-tab:friedman68}; $p_z$\\,ns: methods not significantly different from PSO-VNS in the Holm-adjusted average-rank tests. Unadjusted two-sided Wilcoxon $p$ on the per-case mean wake losses (with the same imputation rule at each threshold) for PSO-VNS vs.\\ PSO and the component-analysis contrasts SSA-VNS vs.\\ RS-VNS, LX-SSA-VNS vs.\\ RS-VNS, SSA-VNS vs.\\ LX-SSA-VNS and PSO-VNS vs.\\ RS-VNS. Same: every conclusion of the paper holds at this threshold (best rank of PSO-VNS; only PSO not significantly different in rank; PSO-VNS vs.\\ PSO not significant; SSA-VNS better than RS-VNS; LX-SSA-VNS not different from RS-VNS; SSA-VNS better than LX-SSA-VNS; PSO-VNS better than RS-VNS; each both unadjusted and Holm-adjusted within its table family), and so do those of the controls in Table~\\ref{tab:X-threshold-controls}.",
                   "tab:X-threshold", "c" * 12,
                   "Thresh. & Unq. & Best & " + " & ".join(stack(*h) for h in (("Rank", "PSO-VNS"), ("Rank", "PSO"), ("$p_z$", "ns"),
                   ("PSO-VNS/", "PSO"), ("SSA-VNS/", "RS-VNS"), ("LX-SSA-VNS/", "RS-VNS"), ("SSA-VNS/", "LX-SSA-VNS"), ("PSO-VNS/", "RS-VNS"))) + " & Same", lines))
    # --- threshold table of the Phase 6 controls (RSD-VNS contrasts, budget split incl. omega = 0.9)
    vm = {"A": "\\,$+$", "B": "\\,$-$", "ns": "\\,$\\cdot$"}
    lines = []
    for t in s["thresholds"]:
        x = TH[str(t)]; ab = x["ablation_case_mean"]; sp = x["split_case_mean"]
        cells = [f"{t}" + (" (paper)" if t == s["base_threshold"] else ""), str(x["imputed_split"])]
        cells += [tp(ab[pk(a, b)]["p"]) + vm[verdict(dict(ab[pk(a, b)], p=max(ab[pk(a, b)]["p"], ab[pk(a, b)]["p_holm"])))]
                  for a, b in RSD_CONTR]
        cells += [tp(sp[pk(a, b)]["p"]) + vm[verdict(sp[pk(a, b)])] for a, b in SPLIT_PAIRS]
        lines.append(" & ".join(cells) + " \\\\")
    T.append(table("table*", "Sensitivity of the Phase-6 control comparisons to the qualification threshold (as Table~\\ref{tab:X-threshold}). Left: case-mean Wilcoxon $p$ (68 cases, unadjusted) of the four component-analysis contrasts with the disc-sampling control RSD-VNS; mark: verdict at $\\alpha=0.05$ both unadjusted and Holm-adjusted over the %d component-analysis contrasts ($+$: first method better, $-$: second better, $\\cdot$: not significant). Right: case-mean Wilcoxon $p$ (12 split cases, unadjusted as in Table~\\ref{M-tab:split}) of the budget-split tests of PSO-VNS (share $\\omega$ of the evaluations for PSO; $\\omega=1$: PSO alone); mark: unadjusted verdict. Unq.: split case--setting cells that do not qualify (imputed maximal differences)." % len(ABL_CONTR),
                   "tab:X-threshold-controls", "c" * (2 + len(RSD_CONTR) + len(SPLIT_PAIRS)),
                   "Thresh. & Unq. & " + " & ".join(stack(LAB[a] + "/", LAB[b]) for a, b in RSD_CONTR) + " & "
                   + " & ".join(stack(SPLIT_LAB[a] + "/", SPLIT_LAB[b]) for a, b in SPLIT_PAIRS), lines, resize=True))
    # --- cluster table
    CL = s["cluster"]; LO = s["loco"]
    cls = s["clusters"]
    cname = {c: ("I" if c.split("-")[0] == "1" else "II") + "/" + c.split("-")[1] for c in cls}
    lines = []
    for a, b in ALL_PAIRS:
        c = CL[pk(a, b)]
        cm = [("--" if c["cluster_means_pp"][k] is None else f"{c['cluster_means_pp'][k]:+.3f}".replace("-", "$-$")) for k in cls]
        lines.append(f"{plab(a, b)} & {c['n_cases']} & ${c['mean_dloss_pp']:+.3f}$ & " + " & ".join(cm) +
                     f" & {c['clusters_favour_first']}/{c['clusters_favour_second']} & {tp(c['sign_test_p'])} & {tp(c['wilcoxon_exact_p'])}"
                     f" & {tp(c['cluster_signflip_p'])} & {tp(c['cr1_p'])} & {tp(s['case_level'][pk(a, b)]['p'])} \\\\")
    T.append(table("table*", "Cluster-level analysis. The 68 cases form %d clusters (data set/radius in m) with nested $N$. $\\overline{\\Delta L}$: mean difference of the case-mean wake losses, first minus second method (pp; negative = first better), over the $n$ cases in which both methods have at least 15 feasible runs; cluster columns: the same mean within each cluster. Fav.: clusters in which the first/second method has the lower mean. Tests over the %d cluster means: exact two-sided sign test and exact Wilcoxon signed-rank test; flip: exact cluster sign-flip test of the case-weighted mean (all $2^{%d}$ sign patterns); CR: cluster-robust $t$ test of the case-weighted mean (CR1 variance, %d d.f.). With %d clusters the smallest attainable exact two-sided $p$ is %s. Case: the case-level Wilcoxon $p$ of the paper (68 cases, imputation rule; unadjusted)."
                   % (len(cls), len(cls), len(cls), len(cls) - 1, len(cls), f"{s['cluster_min_attainable_p']:.3f}"),
                   "tab:X-cluster", "lcc" + "c" * len(cls) + "cccccc",
                   "Pair (first vs.\\ second) & $n$ & $\\overline{\\Delta L}$ & " + " & ".join(cname[k] for k in cls) +
                   " & Fav. & $p_{\\rm sign}$ & $p_{\\rm W}$ & $p_{\\rm flip}$ & $p_{\\rm CR}$ & $p_{\\rm case}$", lines, sep="2pt", resize=True))
    lines = []
    vmap = {"A": "first", "B": "second", "ns": "ns"}
    for a, b in ALL_PAIRS:
        c = CL[pk(a, b)]; lo = LO["pairs"][pk(a, b)]
        chg = ", ".join(cname[k] for k in lo["changed_when_dropping"]) or "--"
        lines.append(f"{plab(a, b)} & ${c['mean_dloss_pp']:+.3f}$ & {ci_tab(c['case_boot_ci95'])} & {ci_tab(c['cluster_boot_ci95'])} & "
                     f"{ci_tab(c['cr1_ci95'])} & {vmap[lo['full_verdict']]} & {tp(lo['p_min'])}--{tp(lo['p_max'])} & "
                     f"{lo['direction_flips']} & {chg} \\\\")
    T.append(table("table*", "Dependence of the case-mean comparisons on the six clusters. 95\\%% intervals of the mean case difference $\\overline{\\Delta L}$ (pp, first minus second): percentile bootstrap over the cases (as in the paper), over the clusters (block bootstrap: %s resamples of the six clusters with replacement, case-weighted mean, fixed seed; with six clusters it tends to be too narrow) and cluster-robust $t$ interval (CR1, 5 d.f.). Leave-one-cluster-out (LOCO): the case-mean Wilcoxon test (imputation rule, unadjusted) is repeated six times, each time without one cluster. Verdict: result on all 68 cases (first / second = that method significantly better at $\\alpha=0.05$; ns = not significant); $p$ range over the six LOCO runs; Flips: LOCO runs in which the sign of $\\overline{\\Delta L}$ changes; Changed: clusters whose removal changes the verdict. The best-ranked method of the eight-method Friedman ranking is %s in every LOCO run%s."
                   % ("10{,}000", "PSO-VNS" if LO["best_ranked_always_PSOBV"] else "\\textbf{not} PSO-VNS",
                      ", and PSO is the only method not significantly different from it in rank" if LO["posthoc_only_PSOC_always"] else ""),
                   "tab:X-loco", "lcccccccc",
                   "Pair (first vs.\\ second) & $\\overline{\\Delta L}$ & Case boot. CI & Cluster boot. CI & CR1 CI & Verdict & LOCO $p$ & Flips & Changed",
                   lines, sep="2pt", resize=True))
    # --- budget split: cluster level and LOCO (12 cases = 2 per cluster)
    CLS = s["cluster_split"]; LOS = s["loco_split"]
    lines = []
    for a, b in SPLIT_PAIRS:
        c = CLS[pk(a, b)]; lo = LOS["pairs"][pk(a, b)]
        cm = [("--" if c["cluster_means_pp"][k] is None else f"{c['cluster_means_pp'][k]:+.3f}".replace("-", "$-$")) for k in cls]
        chg = ", ".join(cname[k] for k in lo["changed_when_dropping"]) or "--"
        lines.append(f"{SPLIT_LAB[a]} vs.\\ {SPLIT_LAB[b]} & {c['n_cases']} & ${c['mean_dloss_pp']:+.3f}$ & " + " & ".join(cm) +
                     f" & {c['clusters_favour_first']}/{c['clusters_favour_second']} & {tp(c['wilcoxon_exact_p'])} & {tp(c['cr1_p'])} & "
                     f"{ci_tab(c['cluster_boot_ci95'])} & {tp(s['split_case_level'][pk(a, b)]['p'])} & {vmap[lo['full_verdict']]} & "
                     f"{tp(lo['p_min'])}--{tp(lo['p_max'])} & {chg} \\\\")
    T.append(table("table*", "Budget-split tests of PSO-VNS (Table~\\ref{M-tab:split}; share $\\omega$ of the 6,030 evaluations for PSO, $\\omega=1$: PSO alone) at the cluster level. The 12 split cases are two per (data set, radius) cluster. Columns as in Tables~\\ref{tab:X-cluster} and \\ref{tab:X-loco}: mean case difference (pp, first minus second; negative = first better) overall and per cluster; Fav.: clusters favouring the first/second setting; $p_{\\rm W}$: exact Wilcoxon over the six cluster means; $p_{\\rm CR}$: cluster-robust $t$ (CR1, 5 d.f.); cluster-bootstrap 95\\% CI; $p_{\\rm case}$: case-mean Wilcoxon over the 12 cases (unadjusted, as in the paper); Verdict and LOCO $p$ range (10 cases each) as in Table~\\ref{tab:X-loco}; Changed: clusters whose removal changes the verdict.",
                   "tab:X-cluster-split", "lcc" + "c" * len(cls) + "cccccccc",
                   "Split test & $n$ & $\\overline{\\Delta L}$ & " + " & ".join(cname[k] for k in cls) +
                   " & Fav. & $p_{\\rm W}$ & $p_{\\rm CR}$ & Cluster boot. CI & $p_{\\rm case}$ & Verdict & LOCO $p$ & Changed", lines, sep="2pt", resize=True))
    # --- multiplicity
    MU = s["multiplicity"]
    lines, fam = [], None
    for t in MU["tests"]:
        if t["family"] != fam:
            fam = t["family"]; lines.append(f"\\multicolumn{{8}}{{l}}{{\\emph{{{fam}}}}} \\\\")
        lines.append(f"{t['label']} & {t['n']} & {t['wins']}/{t['losses']} & ${t['dl']:+.3f}$ & {tp(t['p'])} & {tp(t['p_holm_all'])} & {tp(t['q_bh_all'])} & "
                     + ("yes" if t["sig_holm"] else ("\\textbf{no}" if t["sig_raw"] else "--")) + " \\\\")
    T.append(table("table*", "One multiplicity family over all %d case-mean Wilcoxon tests quoted in the main text (tests that appear in two tables are counted once). $n$: cases entering the test (imputation rule); W/L: cases in which the first method has the lower/higher case-mean wake loss; $\\overline{\\Delta L}$: mean difference (pp, first minus second, both-qualified cases); $p$: unadjusted two-sided Wilcoxon $p$; $p_{\\rm Holm}$: Holm-adjusted over all %d tests (family-wise error); $q_{\\rm BH}$: Benjamini--Hochberg adjusted (false discovery rate). Holm sig.: significant after the all-family Holm correction (``--'': not significant even unadjusted; \\textbf{no}: significant only unadjusted). %d of %d tests are significant unadjusted, %d after Holm and %d after Benjamini--Hochberg."
                   % (MU["n_tests"], MU["n_tests"], MU["n_sig_raw"], MU["n_tests"], MU["n_sig_holm"], MU["n_sig_bh"]),
                   "tab:X-multiplicity", "lccccccc",
                   "Test (first vs.\\ second) & $n$ & W/L & $\\overline{\\Delta L}$ & $p$ & $p_{\\rm Holm}$ & $q_{\\rm BH}$ & Holm sig.", lines))
    # --- subgroup
    SG = s["subgroup"]
    ia, i1, i2 = SG["interaction_all"], SG["interaction_dsI"], SG["interaction_dsII"]
    lines = [
        f"Interaction $N\\ge10$ vs.\\ $N<10$, all cases & {SG['n_cases']} & ${ia['diff_pp']:+.3f}$ & {tp(ia['p_unrestricted'])} & {tp(ia['p_stratified'])} \\\\",
        f"Interaction $N\\ge10$ vs.\\ $N<10$, Data Set I & {i1['n_large'] + i1['n_small']} & ${i1['diff_pp']:+.3f}$ & {tp(i1['p_unrestricted'])} & {tp(i1['p_stratified'])} \\\\",
        f"Interaction $N\\ge10$ vs.\\ $N<10$, Data Set II & {i2['n_large'] + i2['n_small']} & ${i2['diff_pp']:+.3f}$ & {tp(i2['p_unrestricted'])} & {tp(i2['p_stratified'])} \\\\"]
    for tag, txt in (("all", "all cases"), ("dsI", "Data Set I"), ("dsII", "Data Set II")):
        v = SG[f"spearman_{tag}"]
        lines.append(f"Spearman $\\rho(N, \\Delta L)$, {txt} & {v['n']} & ${v['rho']:+.2f}$ / ${v['stratified_mean_rho']:+.2f}$ & {tp(v['p'])} & {tp(v['stratified_p'])} \\\\")
    for tag, txt in (("large", "$N\\ge10$"), ("small", "$N<10$"), ("all", "all $N$")):
        v = SG[f"dsII_vs_dsI_{tag}"]
        lines.append(f"Data Set II minus I, same $(r,N)$, {txt} ({v['dsII_larger_gain']}/{v['dsI_larger_gain']}) & {v['n_pairs']} & ${v['mean_pp']:+.3f}$ & {tp(v['p'])} & -- \\\\")
    T.append(table("table", "Tests of the post hoc $N\\ge10$ pattern of PSO-VNS vs.\\ PSO. $\\Delta L$: case-mean wake loss of PSO-VNS minus that of PSO (pp; negative = PSO-VNS better). Interaction: mean $\\Delta L$ of the $N\\ge10$ cases minus that of the $N<10$ cases, two-sided permutation test (%s permutations, fixed seed) of the group labels over all cases (unrestr.) or within each (data set, radius) cluster (strat.). Spearman: correlation of $\\Delta L$ with $N$ (value: pooled / mean within-cluster $\\rho$; $p$: asymptotic / within-cluster permutation of $\\Delta L$, which uses $N$ only within a farm and so is not confounded with the radius). Data Set II minus I: difference of $\\Delta L$ between the two data sets over the identical $(r, N)$ cases (in parentheses: pairs in which the PSO-VNS advantage is larger in Data Set II / in Data Set I), Wilcoxon signed-rank test."
                   % ("{:,}".format(SG["permutations"]).replace(",", "{,}")),
                   "tab:X-subgroup", "lccc c", "Test & $n$ & Estimate & $p$ (unrestr.) & $p$ (strat.)", lines, sep="3pt"))
    # --- energy
    EN = s["energy"]
    lines = []
    for k, txt in (("PSOBV-PSOC", "PSO-VNS vs.\\ PSO"), ("PSOBV-PSOC_dsIILarge", "PSO-VNS vs.\\ PSO (Data Set II, $N\\ge10$)"),
                   ("SSABV-RSVNS", "SSA-VNS vs.\\ RS-VNS"), ("PSOBV-RSVNS", "PSO-VNS vs.\\ RS-VNS")):
        v = EN[k]
        lines.append(f"{txt} & {v['n_cases']} & ${v['mean_dloss_pp']:+.3f}$ & ${v['mean_gain_pct_aep']:+.3f}$ {ci_tab(v['ci95_gain_pct_aep'])} & "
                     f"${v['mean_gain_mwh_yr']:+.0f}$ {ci_tab(v['ci95_gain_mwh_yr'], 0)} & ${v['mean_gain_mwh_yr_per_turbine']:+.1f}$ & "
                     f"${v['mean_gain_mwh_yr_largest_N']:+.0f}$ ({v['largest_N']}) \\\\")
    T.append(table("table*", "Headline differences in benchmark energy. $\\overline{\\Delta L}$: mean difference of the case-mean wake losses (pp of the wake-free AEP, first minus second). Gain: mean over the cases of the first method's benchmark AEP minus the second's, in \\% of the second's and in MWh/yr of benchmark (model) energy per case, with 95\\% percentile bootstrap CIs over the cases (10{,}000 resamples, fixed seed); per turbine: gain divided by $N$; largest $N$: mean gain over the cases with the largest $N$ of the set. Benchmark AEP $=$ objective$/15\\times8.76$~MWh/yr (Section~\\ref{M-sec:powermodel}); the Kusiak--Song benchmark has a capacity factor far above real sites (62\\% wake-free for Data Set~I), so these values indicate the order of magnitude of the differences, not the energy of a real farm.",
                   "tab:X-energy", "lcccccc",
                   "Pair (first vs.\\ second) & Cases & $\\overline{\\Delta L}$ (pp) & Gain (\\% AEP) & Gain (MWh/yr) & per turbine & largest $N$", lines))
    # --- wake-model uncertainty
    MS = s["model_shift"]
    lines = [f"All feasible final layouts & {MS['n_layouts']:,} & {MS['mean_abs_shift_pp']:.3f} & {MS['median_abs_shift_pp']:.3f} & "
             f"[{MS['q10_abs_shift_pp']:.3f}, {MS['q90_abs_shift_pp']:.3f}] & ${MS['mean_signed_shift_pp']:+.3f}$ \\\\".replace(",", "{,}", 1)]
    for ds, txt in (("1", "Data Set I"), ("2", "Data Set II")):
        v = MS["by_dataset"][ds]
        lines.append(f"\\quad {txt} & {v['n']:,} & {v['mean_abs_shift_pp']:.3f} & {v['median_abs_shift_pp']:.3f} & -- & -- \\\\".replace(",", "{,}", 1))
    lines.append("\\midrule")
    lines.append("\\multicolumn{6}{l}{\\emph{Case-mean difference PSO-VNS minus other method (pp)}} \\\\")
    lines.append("Pair & Cases & Jensen & Gaussian & mean $|$change$|$ & sign changes \\\\")
    for b in [x for x in MAIN8 if x != FOCUS]:
        v = MS["pairs"][pk(FOCUS, b)]
        lines.append(f"{plab(FOCUS, b)} & {v['n_cases']} & ${v['mean_diff_jensen_pp']:+.3f}$ & ${v['mean_diff_gauss_pp']:+.3f}$ & "
                     f"{v['mean_abs_change_pp']:.3f} & {v['sign_changes']} \\\\")
    T.append(table("table", "Wake-model uncertainty from the study's own layouts. Every feasible final layout of the eight main methods (68 cases, 6{,}030 evaluations) is re-evaluated (not re-optimized) with the benchmark Jensen model and with the Gaussian wake of Bastankhah and Port\\'e-Agel ($k^*=0.04$, Table~\\ref{tab:robust}); wake loss in pp of the wake-free power. Top: absolute difference $|L_{\\rm Gauss}-L_{\\rm Jensen}|$ of the same layout (mean, median, 10--90\\%% quantiles) and mean signed difference, to be compared with the equivalence margin of $\\pm%g$~pp. Bottom: mean over the cases (both methods $\\ge15$ feasible runs) of the case-mean difference under each model, the mean absolute change of the per-case difference between the models, and the cases in which its sign differs." % EQ_MARGIN,
                   "tab:X-modelshift", "lccccc", "Layouts & $n$ & mean & median & 10--90\\% & signed", lines, sep="3pt"))
    write_tables_r2(s, T)
    hdr = ("% generated by mpce_inference_extra.py (" + s["generated"] + ") -- do not edit by hand\n"
           "% Supplementary tables of the inference-robustness analyses (Phase 6, W2); requires booktabs\n\n")
    open(fn, "w").write(hdr + "\n".join(T))
    print(f"wrote {fn} ({len(T)} tables)")


def write_tables_r2(s, T):
    """supplementary tables of the review-round-2 analyses: tab:X-equiv-levels, tab:X-bayes, tab:X-heterogeneity,
    tab:X-beyond, tab:X-hl, tab:X-mcnemar."""
    EL = s["equivalence_levels"]; P_ = EL["pairs"]; mg = EL["margin_pp"]
    f3 = lambda v: f"{v:+.3f}".replace("-", "$-$")
    yn = lambda v: "yes" if v else "no"
    tpt = lambda r: ("$\\le$" + tp(r["p_tost"]).strip("$") if r.get("p_tost_at_floor") else tp(r["p_tost"]))
    abl = [k for k in EQ_PAIRS if k in P_][:13]; mn = [k for k in EQ_PAIRS if k in P_][13:]
    # --- equivalence at three levels
    def row(k):
        r = P_[k]; a, b = k.split("-"); sd, cs, cl = r["seed"], r["case"], r["cluster"]
        return (f"{plab(a, b)} & {r['n_cases']} & ${r['mean_dloss_pp']:+.3f}$ & {ci_tab(sd['ci90'])} & {up3(sd['min_margin_pp'])} & "
                f"{yn(sd['equivalent'])} & {ci_tab(cs['ci90'])} & {up3(cs['min_margin_pp'])} & {tpt(cs)} & {yn(cs['equivalent'])} & "
                f"{ci_tab(cl['cr2']['ci90'])} & {up3(cl['cr2']['min_margin_pp'])} & {tp(cl['wild']['p_tost'])} & "
                f"{up3(cl['wild']['min_margin_pp'])} & {yn(cl['equivalent'])} \\\\")
    lines = (["\\multicolumn{15}{l}{\\emph{Component-analysis contrasts}} \\\\"] + [row(k) for k in abl] +
             ["\\midrule", "\\multicolumn{15}{l}{\\emph{Main comparison, PSO-VNS vs.\\ each method (PSO and VNS: see above)}} \\\\"] +
             [row(k) for k in mn])
    f = P_["PSOBV-PSOC"]["cluster"]
    T.append(table("table*", "Practical equivalence (margin $m=%g$~pp, set after the primary analysis) of the pairs of Table~\\ref{tab:equivalence} at three levels of inference. $\\overline{\\Delta L}$: mean over the $n$ cases of the case-mean wake-loss difference, first minus second method (pp; negative = first better; the same estimate at every level). \\emph{Seed level} (fixed benchmark: the %d cases are fixed and only the run-to-run variability is random): %s bootstrap resamples of the 30 seed-paired runs within every case (same seeds for both methods; case means over the feasible resampled runs), percentile 90\\%% CI. \\emph{Case level} (generalization to exchangeable cases, as Table~\\ref{tab:equivalence}): percentile bootstrap over the cases; $p_{\\rm TOST}$ as there. \\emph{Cluster level} (generalization to new farms; %d (data set, radius) clusters): cluster-robust 90\\%% CI with the CR2 (bias-reduced) variance and a $t$ distribution with %d d.f. (Bell--McCaffrey d.f.\\ %.1f for PSO-VNS vs.\\ PSO), and a restricted wild-cluster bootstrap-$t$ TOST (CR2-studentized, all %s draws of the Webb six-point weights enumerated; $m_{\\min}$ by test inversion). $m_{\\min}$: smallest margin at which equivalence holds (rounded up). Eq.: equivalent at $\\pm m$ (cluster level: CR2 CI inside $\\pm m$ and wild $p_{\\rm TOST}\\le0.05$)."
                   % (mg, s["n_cases"], f"{BOOT_N:,}".replace(",", "{,}"), f["n_clusters"], f["df"], f["df_bm"],
                      f"{f['wild']['draws']:,}".replace(",", "{,}")),
                   "tab:X-equiv-levels", "lcccccccccccccc",
                   "& & & \\multicolumn{3}{c}{Seed level} & \\multicolumn{4}{c}{Case level} & \\multicolumn{5}{c}{Cluster level} \\\\\n"
                   "\\cmidrule(lr){4-6}\\cmidrule(lr){7-10}\\cmidrule(lr){11-15}\n"
                   "Pair (first vs.\\ second) & $n$ & $\\overline{\\Delta L}$ & 90\\% CI & $m_{\\min}$ & Eq. & 90\\% CI & $m_{\\min}$ & "
                   "$p_{\\rm TOST}$ & Eq. & CR2 90\\% CI & $m_{\\min}$ & wild $p_{\\rm TOST}$ & wild $m_{\\min}$ & Eq.",
                   lines, sep="2pt", resize=True))
    # --- Bayesian posterior means
    pr2 = lambda v: "$>$0.99" if v > 0.99 else "$<$0.01" if v < 0.01 else f"{v:.2f}"
    def brow(lab, n, y):
        return (f"{lab} & {n} & {pr2(y['p_a_better'])} & {pr2(y['p_rope'])} & {pr2(y['p_b_better'])} & "
                + " & ".join(f"{t:.2f} [{c[0]:.2f}, {c[1]:.2f}]" for t, c in zip(y["mean_theta"], y["ci95_theta"])) + " \\\\")
    H = s["heterogeneity"]
    lines = [brow(plab(*k.split("-")), P_[k]["n_cases"], P_[k]["bayes"]) for k in EQ_PAIRS if k in P_]
    lines += ["\\midrule", "\\multicolumn{8}{l}{\\emph{PSO-VNS vs.\\ PSO, subsets of the cases}} \\\\"]
    for tag, txt in (("small", "$N<10$"), ("large", "$N\\ge10$"), ("trivial", "wake loss of PSO-VNS $<%g$~pp" % H["trivial_loss_pp"]),
                     ("nontrivial", "wake loss of PSO-VNS $\\ge%g$~pp" % H["trivial_loss_pp"])):
        lines.append(brow(txt, H[tag]["n_cases"], H[tag]["bayes"]))
    T.append(table("table*", "Bayesian signed-rank test \\cite{Benavoli2017} on the case-mean differences (first minus second; region of practical equivalence $(-%g, %g)$~pp; prior strength $s=%g$ at $z_0=0$; %s posterior samples, fixed seed; as Table~\\ref{tab:equivalence}). $\\theta_{\\rm A}$, $\\theta_{\\rm rope}$, $\\theta_{\\rm B}$: posterior shares of the Walsh averages $(d_i+d_j)/2$ below $-m$ (first better), inside and above $+m$ (second better) -- the chance that a new pair of cases favours the first method by more than $m$, is within $\\pm m$, or favours the second. $P_{\\rm A}$, $P_{\\rm rope}$, $P_{\\rm B}$: posterior probability that each region is the \\emph{most probable} of the three (the quantity of Table~\\ref{tab:equivalence}); this is not the probability that the mean difference lies within $\\pm m$. Posterior means of $\\theta$ with 95\\%% credible intervals."
                   % (mg, mg, BAYES_S, f"{BAYES_N:,}".replace(",", "{,}")),
                   "tab:X-bayes", "lccccccc",
                   "Pair (first vs.\\ second) & $n$ & $P_{\\rm A}$ & $P_{\\rm rope}$ & $P_{\\rm B}$ & $\\theta_{\\rm A}$ & $\\theta_{\\rm rope}$ & $\\theta_{\\rm B}$",
                   lines, sep="2.5pt"))
    # --- heterogeneity of PSO-VNS - PSO
    cls = s["clusters"]; cname = {c: ("I" if c.split("-")[0] == "1" else "II") + "/" + c.split("-")[1] for c in cls}
    lines = [f"All cases & {H['n_cases']} & {H['beyond_first']}/{H['beyond_second']} & ${P_['PSOBV-PSOC']['mean_dloss_pp']:+.3f}$ & "
             f"{ci_tab(P_['PSOBV-PSOC']['case']['ci90'])} & {up3(P_['PSOBV-PSOC']['case']['min_margin_pp'])} & "
             f"{yn(P_['PSOBV-PSOC']['case']['equivalent'])} & "
             + " / ".join(pr2(v) for v in (P_['PSOBV-PSOC']['bayes']['p_a_better'], P_['PSOBV-PSOC']['bayes']['p_rope'],
                                           P_['PSOBV-PSOC']['bayes']['p_b_better'])) + " \\\\"]
    for tag, txt in (("small", "$N<10$"), ("large", "$N\\ge10$"), ("trivial", "loss $<%g$~pp" % H["trivial_loss_pp"]),
                     ("nontrivial", "loss $\\ge%g$~pp" % H["trivial_loss_pp"])):
        v = H[tag]; y = v["bayes"]
        lines.append(f"{txt} & {v['n_cases']} & {v['beyond_first']}/{v['beyond_second']} & ${v['mean_dloss_pp']:+.3f}$ & {ci_tab(v['ci90'])} & "
                     f"{up3(v['min_margin_pp'])} & {yn(v['equivalent'])} & "
                     + " / ".join(pr2(q) for q in (y["p_a_better"], y["p_rope"], y["p_b_better"])) + " \\\\")
    lines += ["\\midrule", "\\multicolumn{8}{l}{\\emph{Cluster means (data set/radius in m), pp:} " +
              ", ".join(f"{cname[k]} {f3(H['cluster_means'][k])}" for k in cls) +
              f"; beyond $\\pm m$: {len(H['clusters_beyond_first'])} favour PSO-VNS, {len(H['clusters_beyond_second'])} favour PSO}} \\\\"]
    it = H["interaction"]; s1, s2 = H["slope_dsI"], H["slope_dsII"]
    lines += ["\\midrule",
              "\\multicolumn{8}{l}{\\emph{Data set $\\times$ density $\\phi=N(\\ell_{\\min}/2r)^2$ (slope, pp per unit $\\phi$; studentized permutation $p$)}} \\\\",
              f"\\multicolumn{{8}}{{l}}{{Interaction, $\\Delta L_{{\\rm II}}-\\Delta L_{{\\rm I}}$ on $\\phi$ over the {it['n_pairs']} identical $(r,N)$ cases: "
              f"slope {f3(it['slope_pp_per_unit_phi'])}, $p={tp(it['p_signflip']).strip('$')}$ (sign flips)}} \\\\",
              f"\\multicolumn{{8}}{{l}}{{Data Set I: slope {f3(s1['slope_pp_per_unit_phi'])} (PSO better at higher density), $p={tp(s1['p_perm']).strip('$')}$; "
              f"Data Set II: slope {f3(s2['slope_pp_per_unit_phi'])} (PSO-VNS better), $p={tp(s2['p_perm']).strip('$')}$ (within-radius permutation)}} \\\\"]
    T.append(table("table*", "Heterogeneity of PSO-VNS vs.\\ PSO ($\\Delta L$: case-mean wake loss of PSO-VNS minus that of PSO, pp; negative = PSO-VNS better). Beyond: cases with $\\Delta L<-m$ / $\\Delta L>+m$ ($m=%g$~pp); 90\\%% CI, $m_{\\min}$, Eq.\\ as the case level of Table~\\ref{tab:X-equiv-levels}; Bayes: $P_{\\rm A}/P_{\\rm rope}/P_{\\rm B}$ (region most probable; A = PSO-VNS better) of the Bayesian signed-rank test (Table~\\ref{tab:X-bayes}). Subsets: turbine number $N<10$ / $N\\ge10$, and cases in which the case-mean wake loss of PSO-VNS is below / at least %g~pp (a case with a trivial wake loss cannot differ by $m$). Density: $\\phi=N(\\ell_{\\min}/2r)^2$ with $\\ell_{\\min}=%d$~m, the nominal share of the farm disc covered by the $N$ exclusion discs (range %.2f--%.2f); stated after review, before this test was run. Interaction test: slope of the difference between the data sets over the identical $(r,N)$ cases, HC3-studentized, %s random sign flips of the pairs (valid when the data sets also differ in level); per data set: regression on $\\phi$ with radius fixed effects, HC3-studentized, permutation of $\\Delta L$ within radius. The cases beyond $\\pm m$ are listed in Table~\\ref{tab:X-beyond}."
                   % (mg, H["trivial_loss_pp"], SMIN_M, it["phi_range"][0], it["phi_range"][1], f"{it['perm']:,}".replace(",", "{,}")),
                   "tab:X-heterogeneity", "lccccccc",
                   "Cases & $n$ & Beyond & $\\overline{\\Delta L}$ & 90\\% CI & $m_{\\min}$ & Eq. & Bayes", lines, sep="3pt"))
    # --- the cases beyond +-m
    cb = sorted(H["cases_beyond"], key=lambda c: (c["dataset"], c["radius"], c["N"]))
    lines = [f"{'I' if c['dataset'] == '1' else 'II'} & {c['radius']} & {c['N']} & {c['phi']:.2f} & ${c['d']:+.3f}$ & "
             f"{'PSO-VNS' if c['d'] < 0 else 'PSO'} \\\\" for c in cb]
    T.append(table("table", "The %d of %d cases in which PSO-VNS and PSO differ by more than the margin ($|\\Delta L|>%g$~pp; $\\Delta L$ = PSO-VNS minus PSO, pp): %d favour PSO-VNS (up to %.2f~pp) and %d favour PSO (up to %.2f~pp). $\\phi=N(\\ell_{\\min}/2r)^2$."
                   % (len(cb), H["n_cases"], mg, H["beyond_first"], H["max_first"], H["beyond_second"], H["max_second"]),
                   "tab:X-beyond", "cccccc", "Data set & $r$ (m) & $N$ & $\\phi$ & $\\Delta L$ & Better", lines, sep="4pt"))
    # --- Hodges-Lehmann
    HLr = s["hodges_lehmann"]; lines = []
    for a, b in ALL_PAIRS:
        v = HLr[pk(a, b)]
        lab = ("\\textbf{" + plab(a, b) + "}") if v["key_pair"] else plab(a, b)
        lines.append(f"{lab} & {v['n']} & ${v['mean']:+.3f}$ & {ci_tab(v['ci95_mean'])} & ${v['hl']:+.3f}$ & {ci_tab(v['ci95'])} & "
                     f"{v['wins_first']}/{v['wins_second']} & {tp(v['wilcoxon_p'])} \\\\")
    T.append(table("table*", "Estimates that match the tests. Mean: mean case-mean difference (first minus second, pp) with the 95\\% percentile bootstrap CI over the cases (as in the paper). HL: Hodges--Lehmann estimate (median of the $n(n+1)/2$ Walsh averages of the case-mean differences), the location estimate that belongs to the Wilcoxon signed-rank test, with the 95\\% CI from the exact inversion of that test (exact for untied differences; zero differences are kept). W/L and $p_{\\rm W}$: cases in which the first/second method has the lower case mean and the case-mean Wilcoxon $p$ of the paper (zero differences dropped; imputation when only one method qualifies, so $n$ of the test can exceed the $n$ of the estimates). Bold: the pairs discussed in the text."
                   , "tab:X-hl", "lccccccc",
                   "Pair (first vs.\\ second) & $n$ & Mean & 95\\% CI & HL & 95\\% CI & W/L & $p_{\\rm W}$", lines, sep="3pt", resize=True))
    # --- McNemar
    Mc = s["hr_mcnemar"]
    lines = [f"PSO-VNS feasible & {Mc['both']} & {Mc['only_first']} \\\\", f"PSO-VNS infeasible & {Mc['only_second']} & {Mc['neither']} \\\\"]
    T.append(table("table", "Feasibility of the final layouts on the Horns Rev~1 16-turbine block (6,030 evaluations, random starts), PSO-VNS vs.\\ PSO, %d runs paired by seed. Exact McNemar test on the discordant pairs: $p=%s$ (two-sided). One site and one initialization; on the %d benchmark cases both methods are always feasible."
                   % (Mc["n_seeds"], tp(Mc["p_exact"]).strip("$"), s["n_cases"]),
                   "tab:X-mcnemar", "lcc", "& PSO feasible & PSO infeasible", lines, sep="4pt"))


if __name__ == "__main__":
    main()

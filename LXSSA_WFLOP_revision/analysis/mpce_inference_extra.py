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
  5. Effect sizes in energy terms: headline differences as % of the benchmark AEP and in MWh/yr of benchmark
     (model) energy. The benchmark objective is 15 x the expected farm power in kW (03_model.tex), so the gross
     benchmark AEP is objective / 15 x 8.76 MWh/yr.

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
ABL = ["PSOBV", "SSABV", "LXBV", "RSVNS", "BVNS", "PSOC", "SSA", "LXSSA"]
FOCUS = "PSOBV"
SPLITCASES = [(ds, r, n) for ds in ("1", "2") for r, n in ((500, 6), (750, 8), (1000, 10), (500, 10), (750, 12), (1000, 15))]
THRESHOLDS = [1, 10, 15, 20, 25, 30]
BASE_THR = 15
BOOT_N, BOOT_SEED = 10000, 20260928
PERM_SEED = 20260929
ENERGY_FACTOR = 8.76 / 15.0          # MWh/yr per unit of benchmark objective (objective = 15 x expected power in kW)
LAB = {"PSOBV": "PSO-VNS", "PSOC": "PSO", "SSABV": "SSA-VNS", "SSA": "SSA", "LXSSA": "LX-SSA", "DE": "DE",
       "BVNS": "VNS", "SLSQP": "MS-SLSQP", "LXBV": "LX-SSA-VNS", "RSVNS": "RS-VNS",
       "PSOBV25": "PSO-VNS ($\\omega=0.25$)", "PSOBV75": "PSO-VNS ($\\omega=0.75$)"}
MAC = {"PSOBV": "PSOVNS", "PSOC": "PSO", "SSABV": "SSAVNS", "SSA": "SSA", "LXSSA": "LXSSA", "DE": "DE",
       "BVNS": "VNS", "SLSQP": "MSSLSQP", "LXBV": "LXSSAVNS", "RSVNS": "RSVNS"}
# main-table pairs (PSO-VNS vs each method) and the eleven component-analysis contrasts (first vs second)
MAIN_PAIRS = [(FOCUS, b) for b in MAIN8 if b != FOCUS]
ABL_CONTR = [("PSOBV", "PSOC"), ("SSABV", "SSA"), ("LXBV", "LXSSA"), ("PSOBV", "BVNS"), ("PSOBV", "RSVNS"),
             ("SSABV", "RSVNS"), ("LXBV", "RSVNS"), ("PSOBV", "SSABV"), ("PSOBV", "LXBV"), ("SSABV", "LXBV"),
             ("LXSSA", "SSA")]
ALL_PAIRS = list(dict.fromkeys(MAIN_PAIRS + ABL_CONTR))          # 15 unique pairs
# the conclusions of the paper that the threshold / LOCO analyses re-test (component analysis: D1 of PHASE4.md)
KEY_CONTR = [("PSOBV", "PSOC"), ("SSABV", "RSVNS"), ("LXBV", "RSVNS"), ("SSABV", "LXBV"), ("PSOBV", "RSVNS")]
EXPECT = {("PSOBV", "PSOC"): "ns", ("SSABV", "RSVNS"): "A", ("LXBV", "RSVNS"): "ns", ("SSABV", "LXBV"): "A",
          ("PSOBV", "RSVNS"): "A"}     # A = first method significantly better (p < 0.05, lower mean loss); ns = p >= 0.05


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
    for exp in ("rsvns", "psoc", "psobv", "slsqp", "psosplit"):
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


def listing(items):
    items = list(items)
    if len(items) <= 1:
        return "".join(items)
    return ", ".join(items[:-1]) + " and " + items[-1]


WORD = {0: "none", 1: "one", 2: "two", 3: "three", 4: "four", 5: "five", 6: "six", 7: "seven", 8: "eight", 9: "nine",
        10: "ten", 11: "eleven", 12: "twelve"}


def table(env, caption, label, spec, header, lines, sep="2.5pt", pos="!htb", foot=None):
    body = "\n".join(lines)
    ft = "" if not foot else "\n" + "\n".join(foot)
    return (f"\\begin{{{env}}}[{pos}]\n\\centering\n\\caption{{{caption}}}\n\\label{{{label}}}\n"
            f"\\scriptsize\\setlength{{\\tabcolsep}}{{{sep}}}\n\\begin{{tabular}}{{{spec}}}\n\\toprule\n{header} \\\\\n"
            f"\\midrule\n{body}\n\\bottomrule{ft}\n\\end{{tabular}}\n\\end{{{env}}}\n")


# ------------------------------------------------------------------ 1. qualification threshold
def threshold_block(G, SP):
    out = {}
    for thr in THRESHOLDS:
        S = case_stats(G, sorted(set(MAIN8) | set(ABL)), thr)
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
        out[thr] = dict(threshold=thr, friedman=FR, ablation_friedman=dict(avg_rank=FA["avg_rank"], best_ranked=FA["best_ranked"]),
                        main_case_mean=main, ablation_case_mean=abl,
                        imputed_main=int(sum(main[b]["imputed"] for b in main)),
                        dropped_main=int(sum(main[b]["dropped"] for b in main)),
                        imputed_ablation=int(sum(abl[k]["imputed"] for k in abl)),
                        unqualified_cells_main=int((~S[S.Algorithm.isin(MAIN8)].Qualified).sum()),
                        conclusions=concl, all_unchanged=bool(all(concl.values())))
    return out


# ------------------------------------------------------------------ 2. clusters
def cluster_block(S, clusters, perm_seed):
    cl_of = {c: (c[0], c[1]) for c in S.set_index(CASE).index.unique()}
    res = {}
    G = len(clusters)
    signs = np.array(list(itertools.product([1, -1], repeat=G)))       # all 2^G sign patterns
    tcrit = t_dist.ppf(0.975, G - 1)
    for a, b in ALL_PAIRS:
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


def loco_block(S, clusters, full):
    out = {}
    per = []
    for c in clusters:
        m = ~((S.Dataset == c[0]) & (S.Radius == c[1]))
        Sx = S[m]
        FR = friedman(Sx, MAIN8)
        row = dict(dropped=f"{c[0]}-{c[1]}", n_cases=FR["n_cases"], best_ranked=FR["best_ranked"],
                   posthoc_nonsig=FR["posthoc_nonsig"], avg_rank_PSOBV=FR["avg_rank"]["PSOBV"],
                   avg_rank_PSOC=FR["avg_rank"]["PSOC"], tests={})
        for a, b in ALL_PAIRS:
            x = cm_test(Sx, a, b)
            row["tests"][pk(a, b)] = dict(p=x["p"], wins=x["wins"], losses=x["losses"], mean_dloss_pp=x["mean_dloss_pp"],
                                          verdict=verdict(x))
        per.append(row)
    for a, b in ALL_PAIRS:
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
    return dict(per_cluster=per, pairs=out,
                best_ranked_always_PSOBV=bool(all(r["best_ranked"] == FOCUS for r in per)),
                posthoc_only_PSOC_always=bool(all(r["posthoc_nonsig"] == ["PSOC"] for r in per)))


# ------------------------------------------------------------------ 3. multiplicity
def multiplicity_block(S, G, SP):
    tests = []
    main = {b: cm_test(S, FOCUS, b) for b in MAIN8 if b != FOCUS}
    for b, x in main.items():
        tests.append(dict(key=f"main:{pk(FOCUS, b)}", label=plab(FOCUS, b), family="main table (tab:friedman68)",
                          n=x["n_cases"], p=x["p"], wins=x["wins"], losses=x["losses"], dl=x["mean_dloss_pp"]))
    for a, b in ABL_CONTR:
        if (a, b) in MAIN_PAIRS:
            continue                                                    # same test as in the main table
        x = cm_test(S, a, b)
        tests.append(dict(key=f"abl:{pk(a, b)}", label=plab(a, b), family="component analysis (tab:ablation)",
                          n=x["n_cases"], p=x["p"], wins=x["wins"], losses=x["losses"], dl=x["mean_dloss_pp"]))
    for tag, dsl, txt in (("Large", ("1", "2"), "$N\\ge10$, both data sets"), ("dsILarge", ("1",), "$N\\ge10$, Data Set I"),
                          ("dsIILarge", ("2",), "$N\\ge10$, Data Set II")):
        Sg = S[S.Dataset.isin(dsl) & (S.Turbines >= 10)]
        x = cm_test(Sg, FOCUS, "PSOC")
        tests.append(dict(key=f"sub:{tag}", label=f"PSO-VNS vs.\\ PSO, {txt}", family="post hoc subgroup (Sec. results)",
                          n=x["n_cases"], p=x["p"], wins=x["wins"], losses=x["losses"], dl=x["mean_dloss_pp"]))
    for a, b in (("PSOBV", "PSOBV25"), ("PSOBV", "PSOBV75"), ("PSOBV", "PSOC"), ("PSOBV75", "PSOC")):
        x = cm_test(SP, a, b)
        tests.append(dict(key=f"split:{pk(a, b)}", label=plab(a, b).replace("PSO-VNS vs", "PSO-VNS ($\\omega=0.5$) vs", 1)
                          if a == "PSOBV" else plab(a, b), family="budget split (tab:split, 12 cases)",
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
                     "analysis contrasts not already in the main table, the three N >= 10 subgroups, the four budget-"
                     "split tests); Holm (FWER) and Benjamini-Hochberg (FDR) over this family")


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


# ------------------------------------------------------------------ main
def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default=HERE)
    ap.add_argument("--out-dir", default=HERE)
    ap.add_argument("--perm", type=int, default=20000)
    args = ap.parse_args(argv)
    t0 = time.time()
    A = load(args.data_dir)
    G = A[A.Algorithm.isin(set(MAIN8) | set(ABL))]
    SPm = ["PSOBV25", "PSOBV", "PSOBV75", "PSOC"]
    Sp = A[A.Algorithm.isin(SPm)]
    Sp = Sp[[(d, r, n) in SPLITCASES for d, r, n in zip(Sp.Dataset, Sp.Radius, Sp.Turbines)]]
    SP = case_stats(Sp, SPm, BASE_THR)
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
    TH = threshold_block(G, SP)
    base = TH[BASE_THR]
    rep = dict(ok=None)
    try:
        MS = json.load(open(os.path.join(HERE, "mpce_summary.json")))
        dr = max(abs(base["friedman"]["avg_rank"][a] - MS["main"]["friedman"]["avg_rank"][a]) for a in MAIN8)
        dp = max(abs(math.log10(base["main_case_mean"][b]["p"]) - math.log10(MS["main"]["case_mean_wilcoxon"][b]["p"]))
                 for b in base["main_case_mean"])
        dpa = max(abs(math.log10(base["ablation_case_mean"][k]["p"]) - math.log10(MS["ablation"]["case_mean"][k]["p"]))
                  for k in base["ablation_case_mean"])
        sub = MS["main"]["subgroup_vs_phase1"]["groups"]
        dsub = max(abs(math.log10(cm_test(S[S.Dataset.isin(v["datasets"]) & (S.Turbines >= 10)], FOCUS, "PSOC")["p"])
                       - math.log10(v["p"])) for v in sub.values())
        dsp = max(abs(math.log10(cm_test(SP, *k.split("-"))["p"]) - math.log10(v["p"])) for k, v in MS["split"]["case_mean"].items())
        rep = dict(max_abs_diff_avg_rank=dr, max_abs_diff_log10p_main=dp, max_abs_diff_log10p_ablation=dpa,
                   max_abs_diff_log10p_subgroups=dsub, max_abs_diff_log10p_split=dsp,
                   summary_generated=MS.get("generated"), ok=bool(max(dr, dp, dpa, dsub, dsp) < 1e-6))
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
    CL = cluster_block(S, clusters, PERM_SEED)
    full = {}
    for a, b in ALL_PAIRS:
        x = cm_test(S, a, b); x["verdict"] = verdict(x); full[pk(a, b)] = x
    LO = loco_block(S, clusters, full)
    summ["case_level"] = full
    summ["cluster"] = CL
    summ["loco"] = LO
    summ["cluster_min_attainable_p"] = float(2 / 2 ** len(clusters))
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
    for k, v in EN.items():
        print(f"energy {k}: dL {v['mean_dloss_pp']:+.4f} pp, gain {v['mean_gain_pct_aep']:+.4f} % AEP, "
              f"{v['mean_gain_mwh_yr']:+.1f} MWh/yr (CI {np.round(v['ci95_gain_mwh_yr'], 1)}), largest N {v['mean_gain_mwh_yr_largest_N']:+.1f}")
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
        if EXPECT[(a, b)] == "ns":
            m(f"NXThr{pmac(a, b)}PMin", pval_down(min(ps)), f"{a}-{b}: smallest p over thresholds (rounded down)")
        else:
            m(f"NXThr{pmac(a, b)}PMax", pval_up(max(ps)), f"{a}-{b}: largest p over thresholds (rounded up)")
    nconc = len(ts["conclusions"]); nun = sum(ts["unchanged"].values())
    m("NXThrNConcl", WORD.get(nconc, str(nconc))); m("NXThrNUnchanged", WORD.get(nun, str(nun)))
    changed = [c for c, v in ts["unchanged"].items() if not v]
    m("NXThrChanged", listing(changed).replace("_", " ") if changed else "none")
    imp = [TH[t]["imputed_main"] for t in T]
    m("NXThrImputedMin", str(min(imp))); m("NXThrImputedMax", str(max(imp)))
    # clusters
    CL = s["cluster"]; LO = s["loco"]
    m("NXClusters", WORD[len(s["clusters"])], "number of (data set, radius) clusters")
    m("NXClustersNum", str(len(s["clusters"])))
    m("NXClMinP", pval(s["cluster_min_attainable_p"]), "smallest attainable two-sided exact p over the clusters")
    m("NXBootN", "10{,}000")
    for a, b in ALL_PAIRS:
        k = pk(a, b); c = CL[k]; lo = LO["pairs"][k]; nm = pmac(a, b)
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
    lab = {t["key"]: t["label"] for t in MU["tests"]}
    m("NXMultLostHolm", listing([lab[k] for k in MU["lost_under_holm"]]) if MU["lost_under_holm"] else "none")
    m("NXMultLostBH", listing([lab[k] for k in MU["lost_under_bh"]]) if MU["lost_under_bh"] else "none")
    m("NXMultMaxSigHolmP", pval_up(MU["max_sig_p_holm"]), "largest all-family Holm p among the tests still significant (rounded up)")
    tk = {t["key"]: t for t in MU["tests"]}
    for key, nm in (("abl:SSABV-RSVNS", "SSAVNSvsRSVNS"), ("abl:SSABV-LXBV", "SSAVNSvsLXSSAVNS"), ("abl:LXBV-RSVNS", "LXSSAVNSvsRSVNS"),
                    ("abl:PSOBV-RSVNS", "PSOVNSvsRSVNS"), ("main:PSOBV-PSOC", "PSOVNSvsPSO"), ("sub:dsIILarge", "PSOVNSvsPSOdsIILarge"),
                    ("sub:Large", "PSOVNSvsPSOLarge"), ("sub:dsILarge", "PSOVNSvsPSOdsILarge"), ("split:PSOBV-PSOBV75", "SplitFiftyVsSeventyFive"),
                    ("split:PSOBV-PSOBV25", "SplitFiftyVsTwentyFive"), ("split:PSOBV75-PSOC", "SplitSeventyFiveVsHundred"),
                    ("split:PSOBV-PSOC", "SplitFiftyVsHundred")):
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
    hdr = ["% generated by mpce_inference_extra.py from the per-run data (%s) -- do not edit by hand" % s["generated"],
           "% Phase 6, W2 (inference robustness): qualification threshold, cluster-level / leave-one-cluster-out / cluster",
           "% bootstrap, one Holm (and BH) family over all main-text case-mean tests, the post hoc N >= 10 subgroup, and",
           "% effect sizes in benchmark energy. Supplementary tables: mpce_supp_inference.tex (tab:X-...);",
           "% checks: mpce_check_extra.py (X01...). Signs: pair macros 'AvsB' use first minus second (pp of wake loss,",
           "% negative = first better) except the \\NXEn... macros, which give the GAIN of the first method (positive = more",
           "% energy). p values: 'PMax' / 'MaxSig' bounds are rounded up, 'PMin' bounds rounded down.",
           "% reproduction of mpce_summary.json at the 15-run threshold: %s" % ("ok" if s["reproduction"].get("ok") else "FAILED")]
    lines = hdr + [f"\\newcommand{{\\{n}}}{{{v}}}" + (f"   % {c}" if c else "") for n, v, c in M]
    open(fn, "w").write("\n".join(lines) + "\n")
    print(f"wrote {fn} ({len(M)} macros)")


# ------------------------------------------------------------------ tables
def write_tables(s, fn):
    T = []
    TH = s["threshold"]
    cn = s["threshold_summary"]["conclusions"]
    # --- threshold table
    lines = []
    for t in s["thresholds"]:
        x = TH[str(t)]; fr = x["friedman"]; mc = x["main_case_mean"]; ab = x["ablation_case_mean"]
        cells = [f"{t}" + (" (paper)" if t == s["base_threshold"] else ""), str(x["unqualified_cells_main"]), str(x["imputed_main"]),
                 LAB[fr["best_ranked"]], f"{fr['avg_rank']['PSOBV']:.2f}", f"{fr['avg_rank']['PSOC']:.2f}",
                 ", ".join(LAB[a] for a in fr["posthoc_nonsig"]) or "--",
                 tp(mc["PSOC"]["p"])] + [tp(ab[pk(a, b)]["p"]) for a, b in KEY_CONTR[1:]] + \
                ["yes" if x["all_unchanged"] else "\\textbf{no}"]
        lines.append(" & ".join(cells) + " \\\\")
    T.append(table("table*", "Sensitivity of the case-level conclusions to the qualification threshold (a method is ranked in a case by the mean objective of its feasible runs only if at least this many of its 30 runs are feasible; the paper uses 15). Unq.: case--method cells of the eight-method comparison that do not qualify; Imp.: cases entering the seven case-mean tests of PSO-VNS with an imputed maximal difference (only one method qualifies). Best-ranked method and average ranks of the Friedman ranking of Table~\\ref{M-tab:friedman68}; $p_z$\\,ns: methods not significantly different from PSO-VNS in the Holm-adjusted average-rank tests. Unadjusted two-sided Wilcoxon $p$ on the per-case mean wake losses (with the same imputation rule at each threshold) for PSO-VNS vs.\\ PSO and the component-analysis contrasts SSA-VNS vs.\\ RS-VNS, LX-SSA-VNS vs.\\ RS-VNS, SSA-VNS vs.\\ LX-SSA-VNS and PSO-VNS vs.\\ RS-VNS. Same: every conclusion of the paper holds at this threshold (best rank of PSO-VNS; only PSO not significantly different in rank; PSO-VNS vs.\\ PSO not significant; SSA-VNS better than RS-VNS; LX-SSA-VNS not different from RS-VNS; SSA-VNS better than LX-SSA-VNS; PSO-VNS better than RS-VNS; each both unadjusted and Holm-adjusted within its table family).",
                   "tab:X-threshold", "cccccccccccc",
                   "Threshold & Unq. & Imp. & Best & PSO-VNS & PSO & $p_z$\\,ns & PSO-VNS/PSO & SSA-VNS/RS-VNS & LX-SSA-VNS/RS-VNS & SSA-VNS/LX-SSA-VNS & PSO-VNS/RS-VNS & Same", lines))
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
                   " & Fav. & $p_{\\rm sign}$ & $p_{\\rm W}$ & $p_{\\rm flip}$ & $p_{\\rm CR}$ & $p_{\\rm case}$", lines, sep="2pt"))
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
                   lines, sep="2pt"))
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
    hdr = ("%% generated by mpce_inference_extra.py (%s) -- do not edit by hand\n"
           "%% Supplementary tables of the inference-robustness analyses (Phase 6, W2); requires booktabs\n\n" % s["generated"])
    open(fn, "w").write(hdr + "\n".join(T))
    print(f"wrote {fn} ({len(T)} tables)")


if __name__ == "__main__":
    main()

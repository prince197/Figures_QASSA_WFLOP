"""Per-wind-farm results and statistics of the nine methods, from the STORED run records only (no optimization rerun).

Usage (from analysis/):  python3 rev4_farm_tables.py [--out-dir D]
Outputs (default: this folder): rev4_farm_tables.tex (five booktabs tables, ready to \\input), rev4_farm_tables.json
(every number of the tables, the ranking order of each table, the data status and the cross-checks against the
paper's generated summaries and printed tables), rev4_farm_tables.log (copy of stdout).

Pool (display order of every table):  PSO-VNS [PSOBV], PSO [PSOC], SSA-VNS [SSABV], SSA [SSA], LX-SSA [LXSSA],
DE [DE], GA [GA], VNS [BVNS], MS-SLSQP [SLSQP]. Not used: the old PSO setting (PSO), LXBV, RSVNS, RSDVNS, SLSQPX,
PSOSLSQPX. All runs: random initialization, 30 seed-paired runs (seeds 1-30), stored feasibility labels.

Farms and records
  I    Kusiak-Song benchmark, 68 cases (Data Sets I/II x r = 500/750/1000 m x N), 6,030 evaluations: the eight methods
       exactly as in the paper (rev2_analysis.main_runs = mpce_inference_extra.load: fresh_grid [SSA, LXSSA, DE],
       fresh_vgrid [BVNS], fresh_bgrid [SSABV], mpce_psoc, mpce_psobv, mpce_slsqp) + GA (rev2_ga_s*of10).
       Objective = benchmark objective; reported as wake loss (%) = 100 WakeLoss / Ideal.
  II   Horns Rev 1 16-turbine block: mpce_hrfix_s*of24 (Init random; 6,030 and 30,030) + GA rev2_gahr_s*of2.
       AEP in GWh/yr (5-deg bins, the evaluator used in optimization), AEP loss (%) relative to the wake-free AEP.
  III  Lillgrund 4x4 block: rev2_lg16_s*of4 (6,030) and rev2_lg16b_s*of30 (30,030); GA from rev4_galg_s*of*.csv
       (6,030) and rev4_galgb_s*of*.csv (30,030) if present and complete (30 seeds of that budget); otherwise the GA
       row is marked "pending" and the statistics are those of the eight available methods (rerun this script once
       the GA files exist).
  IV   IEA37 Case Study 1, 16 and 36 turbines, 6,030 and 30,030: mpce_iea16/36 + PSO-VNS arms mpce_iea16p/36p +
       GA rev2_gaiea_s*of2 (split by Turbines). AEP in MWh (official IEA37 model).

Statistics (functions of mpce_results.py = MR, rev2_analysis.py = R2 and rev3_sites.py = R3, imported unchanged)
  * qualification and ranking: MR.case_stats / MR.rank_rule (a method with fewer than 15 of 30 feasible runs is
    ranked below every qualified method, by its number of feasible runs); ranks are within the nine-method pool.
  * benchmark: MR.rank_matrix + MR.friedman_block (Friedman, Iman-Davenport, Holm z tests vs PSO-VNS);
    MR.case_mean_wilcoxon (Wilcoxon signed-rank on the 68 case-mean wake losses, PSO-VNS vs each method, imputation
    where only one method qualifies, Holm over the 8 comparisons); R2.wtl_block (run-level W/T/L: seed-paired
    Wilcoxon per case on MR.goodness, Holm over the 8 comparisons of PSO-VNS in each case).
  * sites: R3.pool_analysis (pairing of seeds, case stats, run-level Friedman over seeds on MR.goodness, seed-paired
    Wilcoxon MR.paired_vs of PSO-VNS vs each other method, Holm within that family) -- the procedure of the paper's
    site tables (tab:hr-site / tab:iea-means and Tables S-sa-hr-pool / S-sa-iea-pool): feasibility first (every
    infeasible run below every feasible run, infeasible runs ordered by minimum spacing), then the objective.
    Additionally (JSON only): the all-run paired endpoint of Section archive-reliability (feasible beats infeasible,
    better objective wins with ties within 1e-9, two infeasible runs tie; rule of rev3_inference.outcome), score
    (W + T/2)/n, exact two-sided sign test over the discordant pairs, Holm over the same family.
No existing file is modified.
"""
import os, sys, glob, re, json, time, argparse, warnings
import numpy as np, pandas as pd
from scipy.stats import binomtest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import mpce_results as MR          # noqa: E402  (functions and constants only)
import rev2_analysis as R2         # noqa: E402
import rev3_sites as R3            # noqa: E402

warnings.filterwarnings("ignore", category=RuntimeWarning)
POOL = ["PSOBV", "PSOC", "SSABV", "SSA", "LXSSA", "DE", "GA", "BVNS", "SLSQP"]
MAIN8 = list(MR.MAIN8)
FOCUS = "PSOBV"
LAB = {"PSOBV": "PSO-VNS", "PSOC": "PSO", "SSABV": "SSA-VNS", "SSA": "SSA", "LXSSA": "LX-SSA", "DE": "DE", "GA": "GA",
       "BVNS": "VNS", "SLSQP": "MS-SLSQP"}
CASE = MR.CASE
CLUSTERS = [(ds, r) for ds in ("1", "2") for r in (500, 750, 1000)]
BUDS = (6030, 30030)
OBJ_TIE = 1e-9
REV2_SUMMARY = os.path.join(HERE, "..", "SWEVO_rev2", "experiments_docs", "rev2_summary.json")
LOG = []


def log(msg=""):
    print(msg, flush=True)
    LOG.append(str(msg))


# ====================================================================== formatting
def fp2(p):
    return "--" if p is None or not np.isfinite(p) else MR.fmt_p(float(p), 2)


def num(v, d, big=False):
    if v is None or not np.isfinite(v):
        return "--"
    s = f"{v:,.{d}f}" if big else f"{v:.{d}f}"
    return s.replace(",", "{,}").replace("-", "$-$")


def bold(s, on):
    return f"\\textbf{{{s}}}" if on else s


def sgn(v, d=3):
    return "--" if v is None or not np.isfinite(v) else f"${v:+.{d}f}$"


# ====================================================================== data
def read_rev4(exp):
    """Merge rev4_<exp>_s<i>of<k>.csv (Coordinates / Curve dropped), as rev2_analysis.read_rev2."""
    by_k = {}
    for fn in glob.glob(os.path.join(HERE, f"rev4_{exp}_s*of*.csv")):
        m = re.search(rf"rev4_{exp}_s(\d+)of(\d+)\.csv$", os.path.basename(fn))
        if m:
            by_k.setdefault(int(m.group(2)), {})[int(m.group(1))] = fn
    if not by_k:
        return None, dict(status="missing", rows=0)
    complete = {k: v for k, v in by_k.items() if set(v) == set(range(k))}
    if complete:
        k = max(complete)
        use, status = list(complete[k].values()), "complete shard set"
    else:
        use, status = [f for v in by_k.values() for f in v.values()], "PARTIAL shard set"
    parts = []
    for f in sorted(use):
        try:
            parts.append(pd.read_csv(f, usecols=lambda c: c not in ("Coordinates", "Curve")))
        except pd.errors.EmptyDataError:
            pass
    if not parts:
        return None, dict(status="empty", rows=0, files=[os.path.basename(f) for f in use])
    df = MR.std_cols(pd.concat(parts, ignore_index=True)).drop_duplicates(MR.KEY, keep="last").reset_index(drop=True)
    df["Source"] = f"rev4_{exp}"
    return df, dict(status=status, rows=int(len(df)), files=[os.path.basename(f) for f in sorted(use)],
                    methods=sorted(df.Algorithm.unique()))


def load_benchmark(status):
    A = R2.main_runs()
    Ng, st = R2.read_rev2("ga", HERE)
    status["benchmark_8"] = dict(source="rev2_analysis.main_runs (mpce_inference_extra.load)", rows=int(len(A[A.Algorithm.isin(MAIN8)])))
    status["rev2_ga"] = st
    D = pd.concat([A[A.Algorithm.isin(MAIN8)], Ng[Ng.Algorithm == "GA"]], ignore_index=True)
    return D


def load_lg(status):
    parts = []
    for exp in ("lg16", "lg16b"):
        d, st = R2.read_rev2(exp, HERE)
        status[f"rev2_{exp}"] = st
        if d is not None:
            parts.append(d)
    L = pd.concat(parts, ignore_index=True).drop_duplicates(MR.KEY, keep="last")
    L = L[L.Dataset == "LG"]
    for exp in ("galg", "galgb"):
        d, st = read_rev4(exp)
        status[f"rev4_{exp}"] = st
        if d is not None:
            d = d[(d.Algorithm == "GA") & (d.Dataset == "LG")]
            if len(d) and abs(float(d.Ideal.iloc[0]) - float(L.Ideal.iloc[0])) > 1e-6:
                status[f"rev4_{exp}"]["warning"] = "wake-free AEP differs from rev2_lg16: model mismatch, GA not used"
                continue
            L = pd.concat([L, d], ignore_index=True)
    return L.drop_duplicates(MR.KEY, keep="last")


def ga_ready(y):
    """GA usable in a site x budget group only with the full set of seeds of PSO-VNS (else: pending / partial)."""
    g = set(y[y.Algorithm == "GA"].Seed)
    ref = set(y[y.Algorithm == FOCUS].Seed)
    if not g:
        return False, "pending"
    if not ref <= g:
        return False, f"partial ({len(g)}/{len(ref)} seeds; not used)"
    return True, "used"


# ====================================================================== all-run paired endpoint (JSON only)
def allrun(y, a, b):
    P = y.pivot_table(index="Seed", columns="Algorithm", values=["Objective", "Feasible"], aggfunc="first")
    fa, fb = P["Feasible"][a].astype(bool).values, P["Feasible"][b].astype(bool).values
    d = (P["Objective"][a] - P["Objective"][b]).values
    O = np.full(len(fa), 0.5)
    O[fa & ~fb] = 1.0; O[~fa & fb] = 0.0
    both = fa & fb
    O[both & (d > OBJ_TIE)] = 1.0; O[both & (d < -OBJ_TIE)] = 0.0
    w, l = int((O == 1).sum()), int((O == 0).sum()); t = int(len(O) - w - l)
    p = float(binomtest(w, w + l).pvalue) if w + l else 1.0
    return dict(wins=w, ties=t, losses=l, n=int(len(O)), score=float(O.mean()), p_sign=p)


# ====================================================================== T-I (benchmark)
def benchmark(out, tex):
    status = out["data_status"]
    D = load_benchmark(status)
    G, drop = R2.pair_seeds(D, CASE, POOL)
    comp = R2.completeness(G, POOL, CASE)
    S = MR.case_stats(G, POOL)
    R = MR.rank_matrix(S, POOL)
    FR = MR.friedman_block(R, POOL, FOCUS)
    others = [a for a in POOL if a != FOCUS]
    CW = MR.case_mean_wilcoxon(S, FOCUS, others)
    WT = R2.wtl_block(G, FOCUS, others)
    n_cases = int(FR["n_cases"])
    res = dict(n_cases=n_cases, n_runs=int(len(G)), dropped_unpaired_runs=int(drop), completeness=comp, methods={})
    cl_key = lambda ds, r: f"DS{'I' if ds == '1' else 'II'}_r{r}"
    for a in POOL:
        s = S[S.Algorithm == a]
        m = dict(feasible_pct=float(100 * G[G.Algorithm == a].Feasible.mean()), avg_rank=float(FR["avg_rank"][a]),
                 sole_best_cases=int(FR["sole_best_count"][a]), n_qualified=int(s.Qualified.sum()),
                 mean_loss_pct_qualified=float(s[s.Qualified].Loss.mean()), clusters={})
        for ds, r in CLUSTERS:
            c = s[(s.Dataset == ds) & (s.Radius == r)]
            Rc = R[(R.index.get_level_values("Dataset") == ds) & (R.index.get_level_values("Radius") == r)]
            m["clusters"][cl_key(ds, r)] = dict(n_cases=int(len(c)), n_qualified=int(c.Qualified.sum()),
                                                mean_loss_pct_qualified=float(c[c.Qualified].Loss.mean()) if c.Qualified.any() else None,
                                                feasible_pct=float(100 * G[(G.Algorithm == a) & (G.Dataset == ds) & (G.Radius == r)].Feasible.mean()),
                                                avg_rank_9pool=float(Rc[a].mean()))
        res["methods"][a] = m
    order = sorted(POOL, key=lambda a: FR["avg_rank"][a])
    res["ranking_order"] = order
    out["tables"]["tab:farm-I-res"] = res
    log("[I] 9-method avg ranks: " + ", ".join(f"{LAB[a]} {FR['avg_rank'][a]:.2f}" for a in order))

    # ---- T-I-a
    keys = [cl_key(ds, r) for ds, r in CLUSTERS]
    best = {}
    for k in keys + ["all"]:
        vals = {}
        for a in POOL:
            c = res["methods"][a]["clusters"][k] if k != "all" else dict(n_cases=n_cases, n_qualified=res["methods"][a]["n_qualified"],
                                                                        mean_loss_pct_qualified=res["methods"][a]["mean_loss_pct_qualified"])
            if c["n_qualified"] == c["n_cases"]:
                vals[a] = round(c["mean_loss_pct_qualified"], 3)
        best[k] = min(vals.values()) if vals else None
    bf = max(round(res["methods"][a]["feasible_pct"], 1) for a in POOL)
    br = min(round(FR["avg_rank"][a], 2) for a in POOL)
    lines = []
    for a in POOL:
        m = res["methods"][a]
        cells = []
        for k in keys + ["all"]:
            c = m["clusters"][k] if k != "all" else dict(n_cases=n_cases, n_qualified=m["n_qualified"], mean_loss_pct_qualified=m["mean_loss_pct_qualified"])
            v = c["mean_loss_pct_qualified"]
            s_ = bold(num(v, 3), c["n_qualified"] == c["n_cases"] and round(v, 3) == best[k])
            if c["n_qualified"] < c["n_cases"]:
                s_ += f"$^{{{c['n_qualified']}}}$"
            cells.append(s_)
        cells.append(bold(f"{m['feasible_pct']:.1f}", round(m["feasible_pct"], 1) == bf))
        cells.append(bold(f"{m['avg_rank']:.2f}", round(m["avg_rank"], 2) == br))
        lines.append(f"{LAB[a]} & " + " & ".join(cells) + " \\\\")
    ncl = {k: res["methods"][FOCUS]["clusters"][k]["n_cases"] for k in keys}
    head = ("& \\multicolumn{3}{c}{Data Set I} & \\multicolumn{3}{c}{Data Set II} & & & \\\\\n"
            "\\cmidrule(lr){2-4}\\cmidrule(lr){5-7}\n"
            "Method & " + " & ".join(f"{r}~m ({ncl[cl_key(ds, r)]})" for ds, r in CLUSTERS)
            + f" & All ({n_cases}) & Feas.\\ (\\%) & Avg.\\ rank")
    tex.append(MR.table(
        "table", "Type~I farms (Kusiak--Song benchmark): mean wake loss (\\%) by data set and farm radius, feasible runs and "
        "average rank of the nine methods.", "tab:farm-I-res", "l" + "c" * 9, head, lines, size="\\scriptsize", sep="2.6pt",
        pos="!htb",
        note="%d benchmark cases (number of cases per cluster in parentheses), 30 seed-paired runs per case and method, 6{,}030 "
             "evaluations, random initialization. Wake loss: per case the mean wake loss of the feasible runs, averaged over the "
             "cases of the cluster in which the method has at least 15 of 30 feasible runs (superscript: number of such cases when "
             "fewer than all; such means are over easier case subsets and are not comparable). Bold: lowest mean among the methods "
             "qualified in every case of the column; highest share of feasible runs; best average rank. Avg.\\ rank: "
             "feasibility-aware rank of the mean feasible objective within the nine-method pool (rule of "
             "Table~\\ref{tab:friedman68}), averaged over the %d cases." % (n_cases, n_cases)))

    # ---- T-I-b
    stats = dict(friedman=dict(chi2=FR["chi2"], df=len(POOL) - 1, p=FR["p"], iman_davenport=FR["iman_davenport"],
                               iman_davenport_df=[len(POOL) - 1, (len(POOL) - 1) * (n_cases - 1)],
                               iman_davenport_p=FR["iman_davenport_p"], n_cases=n_cases, best_ranked=FR["best_ranked"]),
                 avg_rank=FR["avg_rank"], vs_psovns={}, ranking_order=order,
                 notes="p_z: Holm z test of average ranks vs PSO-VNS (8 comparisons); case-mean Wilcoxon: MR.case_mean_wilcoxon "
                       "(Holm over 8); W/T/L: run-level, R2.wtl_block (Holm over the 8 comparisons of PSO-VNS per case); "
                       "mean_dloss_pp = PSO-VNS minus method over the cases where both qualify (negative = PSO-VNS better)")
    for b in others:
        c = CW[b]
        ties = n_cases - c["wins"] - c["losses"]
        stats["vs_psovns"][b] = dict(p_z_holm=FR["p_holm_vs_focus"][b], z=FR["z_vs_focus"][b],
                                     case_mean_better=c["wins"], case_mean_equal_or_neither=ties, case_mean_worse=c["losses"],
                                     p_wilcoxon=c["p"], p_wilcoxon_holm=c["p_holm"], rb=c["rb"],
                                     mean_dloss_pp=c["mean_dloss_pp"], ci95_mean_dloss_pp=c["ci95_mean_dloss_pp"],
                                     n_both_qualified=c["n_both_qualified"], imputed=c["only_focus_qualified"] + c["only_other_qualified"],
                                     run_level_wtl=dict(W=WT[b]["W"], T=WT[b]["T"], L=WT[b]["L"]), median_rb_run_level=WT[b]["median_rb"])
    out["tables"]["tab:farm-I-stats"] = stats
    lines = []
    for a in POOL:
        if a == FOCUS:
            lines.append(f"{LAB[a]} & {FR['avg_rank'][a]:.2f} & -- & -- & -- & -- & -- \\\\")
            continue
        v = stats["vs_psovns"][a]
        lines.append(f"{LAB[a]} & {FR['avg_rank'][a]:.2f} & {fp2(v['p_z_holm'])} & "
                     f"{v['case_mean_better']}/{v['case_mean_equal_or_neither']}/{v['case_mean_worse']} & {fp2(v['p_wilcoxon_holm'])} & "
                     f"{sgn(v['mean_dloss_pp'])} & {v['run_level_wtl']['W']}/{v['run_level_wtl']['T']}/{v['run_level_wtl']['L']} \\\\")
    qn = [f"{v['n_both_qualified']} for {LAB[b]}" for b, v in stats["vs_psovns"].items() if v["n_both_qualified"] < n_cases]
    head = ("Method & Avg.\\ rank & $p_z$ & Case means & $p_W$ & $\\overline{\\Delta L}$ (pp) & W/T/L")
    f = stats["friedman"]
    tex.append(MR.table(
        "table", "Type~I farms (Kusiak--Song benchmark): Friedman test of the nine methods and comparisons of PSO-VNS with each "
        "method.", "tab:farm-I-stats", "lcccccc", head, lines, size="\\footnotesize", sep="4pt", pos="!htb",
        note="%d cases, 30 seed-paired runs, 6{,}030 evaluations. Friedman test on the feasibility-aware case ranks of the nine "
             "methods: $\\chi^2_F=%.1f$ (%d d.f.), $p=%s$; Iman--Davenport $F_F=%.1f$ (%d, %d d.f.), $p=%s$. $p_z$: Holm-adjusted "
             "$p$ of the average-rank $z$ test against PSO-VNS (8 comparisons). Case means: cases in which PSO-VNS has the lower / "
             "equal (or neither method qualifies) / higher mean wake loss (a case in which only one method has at least 15 "
             "feasible runs counts for that method). $p_W$: two-sided Wilcoxon signed-rank test on the %d per-case mean wake "
             "losses, PSO-VNS vs.\\ method (cases in which only one method qualifies imputed in its favour), Holm-adjusted over "
             "the 8 comparisons. $\\overline{\\Delta L}$: mean difference of the case-mean wake losses, PSO-VNS minus method "
             "(percentage points; negative = PSO-VNS better), over the cases in which both qualify (%s). W/T/L: cases in which "
             "PSO-VNS is significantly better / not different / worse in the seed-paired Wilcoxon test of the 30 runs (infeasible "
             "runs below every feasible run; Holm over the 8 comparisons of PSO-VNS in each case, $\\alpha=0.05$)."
             % (n_cases, f["chi2"], f["df"], fp2(f["p"]).strip("$"), f["iman_davenport"], f["iman_davenport_df"][0],
                f["iman_davenport_df"][1], fp2(f["iman_davenport_p"]).strip("$"), n_cases,
                (", ".join(qn) + ", otherwise %d" % n_cases) if qn else "all %d cases" % n_cases)))
    log(f"[I] Friedman chi2={f['chi2']:.1f} p={f['p']:.3g}; ID F={f['iman_davenport']:.1f}")
    log("[I] PSO-VNS vs: " + ", ".join(f"{LAB[b]} pW={v['p_wilcoxon_holm']:.3g} dL={v['mean_dloss_pp']:+.3f} "
                                       f"WTL={v['run_level_wtl']['W']}/{v['run_level_wtl']['T']}/{v['run_level_wtl']['L']}"
                                       for b, v in stats["vs_psovns"].items()))
    return G


# ====================================================================== sites (T-II, T-III, T-IV)
def site_setting(y, unit_big):
    """One site x budget (x scenario): R3.pool_analysis on the available methods of the pool + all-run endpoint."""
    ok, ga_status = ga_ready(y)
    pool = [a for a in POOL if a != "GA" or ok]
    y = y[y.Algorithm.isin(pool)]
    r = R3.pool_analysis(y, pool, {"psovns_vs_all": (FOCUS, None)})
    fam = r["families"]["psovns_vs_all"]
    yy, _ = R2.pair_seeds(y, CASE, [a for a in pool if a in set(y.Algorithm)])
    ar = {b: allrun(yy, FOCUS, b) for b in fam["tests"]}
    for b, h in zip(list(ar), MR.holm([ar[b]["p_sign"] for b in ar]) if ar else []):
        ar[b]["p_sign_holm"] = float(h)
    out = dict(pool=pool, ga_status=ga_status, n_seeds=r["n_seeds"], dropped_unpaired_runs=r["dropped_unpaired_runs"],
               friedman_run_level=r["friedman"], best=r["best"], methods={})
    for a in POOL:
        if a not in r["methods"]:
            out["methods"][a] = dict(status=ga_status if a == "GA" else "missing")
            continue
        v = dict(r["methods"][a])
        v.pop("median_seconds", None)
        if a != FOCUS:
            t = fam["tests"][a]
            v["psovns_vs"] = dict(p=t["p"], p_holm=t["p_holm"], rb=t["rb"], outcome=t["outcome"], n_pairs=t["n_pairs"],
                                  seed_wins=t["wins"], seed_losses=t["losses"], seed_ties=t["ties"],
                                  family_size=fam["comparisons"], all_run_endpoint=ar[a])
        out["methods"][a] = v
    out["ranking_order"] = sorted([a for a in POOL if a in r["methods"]],
                                  key=lambda a: (r["methods"][a]["rank"], -(r["methods"][a]["mean"] if r["methods"][a]["mean"] is not None else -1e18)))
    return out


def site_rows(sets, big, dmean, dsd):
    """Rows of one site panel; sets = list of setting dicts (one per budget)."""
    lines = []
    for a in POOL:
        c = []
        for st in sets:
            v = st["methods"][a]
            if "rank" not in v:
                c.append(f"\\multicolumn{{5}}{{c}}{{{v['status']}}}")
                continue
            q = [w for w in st["methods"].values() if "rank" in w]
            mx_best = max([w["best"] for w in q if w["best"] is not None], default=None)
            isb = a == st["best"]
            mean = "--" if v["mean"] is None else (num(v["mean"], dmean, big) + "\\,$\\pm$\\," + num(v["sd"] if v["sd"] is not None else 0.0, dsd, big))
            c += [f"{v['feasible']}/{v['runs']}",
                  bold(mean, isb and v["mean"] is not None),
                  bold(num(v["best"], dmean, big), v["best"] is not None and v["best"] == mx_best),
                  bold(f"{v['rank']:g}", isb)]
            if a == FOCUS:
                c.append("--")
            else:
                t = v["psovns_vs"]
                mark = "" if t["outcome"] == "T" else ("$^{+}$" if t["outcome"] == "W" else "$^{-}$")
                c.append(fp2(t["p_holm"]) + mark)
        lines.append(f"{LAB[a]} & " + " & ".join(c) + " \\\\")
    return lines


SITE_HEAD = ("& \\multicolumn{5}{c}{6{,}030 evaluations} & \\multicolumn{5}{c}{30{,}030 evaluations} \\\\\n"
             "\\cmidrule(lr){2-6}\\cmidrule(lr){7-11}\n"
             "Method & Feas. & %s & Best & Rank & $p_{\\rm Holm}$ & Feas. & %s & Best & Rank & $p_{\\rm Holm}$")
SITE_NOTE = ("Random initialization, 30 seed-paired runs (seeds 1--30) per method and budget; stored feasibility labels. Feas.: "
             "feasible runs; Mean\\,$\\pm$\\,SD: mean AEP (%s) and its standard deviation over the feasible runs (mean AEP losses "
             "relative to the wake-free AEP: rev4\\_farm\\_tables.json); Best: best feasible run; --: no feasible run. Rank: feasibility-aware "
             "rank of the mean AEP within the methods of the pool (rule of Table~\\ref{tab:friedman68}: a method with fewer than "
             "15 feasible runs ranks below the qualified methods, by its number of feasible runs). Bold: rank-1 method (mean, "
             "rank), highest best run. $p_{\\rm Holm}$: two-sided seed-paired Wilcoxon signed-rank test "
             "of PSO-VNS against the method (feasibility first: every infeasible run below every feasible run, infeasible runs "
             "ordered by minimum spacing), Holm-adjusted over the comparisons of PSO-VNS at that budget%s; $^{+}$/$^{-}$: "
             "PSO-VNS significantly better/worse ($\\alpha=0.05$).")


def friedman_txt(st):
    f = st["friedman_run_level"]
    return f"$\\chi^2_F={f['chi2']:.1f}$ ({f['df']} d.f.), $p={fp2(f['p']).strip('$')}$"


def block_hr(out, tex):
    HH = R3.load_hr()
    out["data_status"]["mpce_hrfix+rev2_gahr"] = dict(rows=int(len(HH)), methods=sorted(HH.Algorithm.unique()))
    import hornsrev_model as hrm
    inst = float(hrm.aep_gwh(hrm.site(16)[0]))
    ideal = float(HH.Ideal.iloc[0])
    res = dict(site="Horns Rev 1, 16-turbine block (5-deg bins)", unit="GWh/yr", installed_aep=inst, ideal_aep=ideal,
               installed_loss_pct=100 * (1 - inst / ideal), budgets={})
    for b in BUDS:
        y = HH[(HH.Budget == b) & (HH.Init == "random")]
        res["budgets"][str(b)] = site_setting(y, False)
        log_site("II HR", b, res["budgets"][str(b)])
    res["ranking_order"] = {k: v["ranking_order"] for k, v in res["budgets"].items()}
    out["tables"]["tab:farm-II"] = res
    S6, S30 = res["budgets"]["6030"], res["budgets"]["30030"]
    lines = [f"Installed block & -- & \\multicolumn{{2}}{{c}}{{{inst:.2f}}} & -- & -- & "
             f"-- & \\multicolumn{{2}}{{c}}{{{inst:.2f}}} & -- & -- \\\\", "\\midrule"]
    lines += site_rows([S6, S30], False, 2, 2)
    fam = len(S6["pool"]) - 1
    tex.append(MR.table(
        "table*", "Type~II farm (Horns Rev~1, 16-turbine block): AEP, feasibility, rank and seed-paired tests of PSO-VNS against "
        "the other eight methods.", "tab:farm-II", "l" + "c" * 10, SITE_HEAD % ("Mean\\,$\\pm$\\,SD", "Mean\\,$\\pm$\\,SD"),
        lines, size="\\scriptsize", sep="2.2pt", pos="!htb",
        note=(SITE_NOTE % ("GWh/yr", " (%d comparisons)" % fam))
        + " 5$^\\circ$ direction bins (the evaluator used in optimization), minimum spacing $4D$; wake-free AEP %.2f~GWh/yr. "
          "Installed block: the $7D$ grid of the farm (AEP loss %.2f\\%%). Run-level Friedman test over the seeds: 6{,}030: %s; 30{,}030: %s."
        % (ideal, res["installed_loss_pct"], friedman_txt(S6), friedman_txt(S30))))


def block_lg(out, tex):
    L = load_lg(out["data_status"])
    import rev2_site_model as lg
    xy0, _ = lg.site(16)
    inst = float(lg.aep_gwh(xy0))
    ideal = float(L.Ideal.iloc[0])
    res = dict(site="Lillgrund, 4x4 block (5-deg bins, minimum spacing 3D)", unit="GWh/yr", installed_aep=inst, ideal_aep=ideal,
               installed_loss_pct=100 * (1 - inst / ideal), budgets={})
    for b in BUDS:
        y = L[(L.Budget == b) & (L.Init == "random") & L.Algorithm.isin(POOL)]
        res["budgets"][str(b)] = site_setting(y, False)
        log_site("III LG", b, res["budgets"][str(b)])
    res["ranking_order"] = {k: v["ranking_order"] for k, v in res["budgets"].items()}
    res["ga_status"] = {k: v["ga_status"] for k, v in res["budgets"].items()}
    out["tables"]["tab:farm-III"] = res
    S6, S30 = res["budgets"]["6030"], res["budgets"]["30030"]
    lines = [f"Installed block & -- & \\multicolumn{{2}}{{c}}{{{inst:.2f}}} & -- & -- & "
             f"-- & \\multicolumn{{2}}{{c}}{{{inst:.2f}}} & -- & -- \\\\", "\\midrule"]
    lines += site_rows([S6, S30], False, 2, 2)
    pend = [b for b in BUDS if res["budgets"][str(b)]["ga_status"] != "used"]
    fams = "; ".join(f"{b:,}".replace(",", "{,}") + f": {len(res['budgets'][str(b)]['pool']) - 1}" for b in BUDS)
    ga_note = ("" if not pend else
               " GA: %s at %s evaluations (record %s not yet available or incomplete); ranks and tests at that budget are those "
               "of the %s other methods. [GA PENDING -- rerun rev4\\_farm\\_tables.py]"
               % ("pending" if all(res["budgets"][str(b)]["ga_status"] == "pending" for b in pend) else "incomplete",
                  " and ".join(f"{b:,}".replace(",", "{,}") for b in pend),
                  " / ".join({6030: "rev4\\_galg", 30030: "rev4\\_galgb"}[b] for b in pend),
                  " / ".join(str(len(res["budgets"][str(b)]["pool"])) for b in pend)))
    tex.append(MR.table(
        "table*", "Type~III farm (Lillgrund, $4\\times4$ block): AEP, feasibility, rank and seed-paired tests of PSO-VNS against the "
        "other methods.", "tab:farm-III", "l" + "c" * 10, SITE_HEAD % ("Mean\\,$\\pm$\\,SD", "Mean\\,$\\pm$\\,SD"),
        lines, size="\\scriptsize", sep="2.2pt", pos="!htb",
        note=(SITE_NOTE % ("GWh/yr", " (comparisons: %s)" % fams))
        + " Siemens SWT-2.3-93, 12-sector Weibull climate of the Lillgrund met mast, Jensen wake $k=0.04$, 5$^\\circ$ direction "
          "bins, minimum spacing $3D$; wake-free AEP %.2f~GWh/yr. Installed block: the $4\\times4$ block of the farm (AEP loss %.2f\\%%). Run-level "
          "Friedman test over the seeds: 6{,}030: %s; 30{,}030: %s.%s"
        % (ideal, res["installed_loss_pct"], friedman_txt(S6), friedman_txt(S30), ga_note)))


def block_iea(out, tex):
    parts = [R2.read_mpce(e) for e in ("iea16", "iea36", "iea16p", "iea36p")]
    g1, st = R2.read_rev2("gaiea", HERE)
    out["data_status"]["rev2_gaiea"] = st
    IE = pd.concat([p for p in parts if p is not None] + [g1[g1.Algorithm == "GA"]], ignore_index=True)
    IE = IE[IE.Dataset.str.startswith("IEA37") & IE.Algorithm.isin(POOL)]
    out["data_status"]["mpce_iea16/36(+p)"] = dict(rows=int(len(IE)), methods=sorted(IE.Algorithm.unique()))
    pv = R3.published_iea()
    res = dict(site="IEA37 Case Study 1", unit="MWh", published={str(n): v for n, v in pv.items()}, settings={})
    for n in (16, 36):
        for b in BUDS:
            y = IE[(IE.Turbines == n) & (IE.Budget == b) & (IE.Init == "random")]
            st_ = site_setting(y, True)
            for a, v in st_["methods"].items():
                if v.get("best") is not None:
                    v["best_vs_published_strict_pct"] = 100 * (v["best"] / pv[n]["strict"] - 1)
            res["settings"][f"{n}T_{b}"] = st_
            log_site(f"IV IEA{n}", b, st_)
    res["ranking_order"] = {k: v["ranking_order"] for k, v in res["settings"].items()}
    out["tables"]["tab:farm-IV"] = res
    lines = []
    for i, n in enumerate((16, 36)):
        if i:
            lines.append("\\midrule")
        lines.append(f"\\multicolumn{{11}}{{l}}{{\\emph{{{n} turbines ($r={1300 if n == 16 else 2000}$~m)}}; published: example "
                     f"{num(pv[n]['example'], 0, True)}, best strict {num(pv[n]['strict'], 0, True)}, "
                     f"best projected {num(pv[n]['projected'], 0, True)}}} \\\\")
        lines += site_rows([res["settings"][f"{n}T_{b}"] for b in BUDS], True, 0, 0)
    frs = "; ".join(f"{n} turbines, " + f"{b:,}".replace(",", "{,}") + ": " + friedman_txt(res["settings"][f"{n}T_{b}"])
                    for n in (16, 36) for b in BUDS)
    tex.append(MR.table(
        "table*", "Type~IV farms (IEA37 Case Study~1, 16 and 36 turbines): AEP, feasibility, rank and seed-paired tests of PSO-VNS "
        "against the other eight methods.", "tab:farm-IV", "l" + "c" * 10, SITE_HEAD % ("Mean\\,$\\pm$\\,SD", "Mean\\,$\\pm$\\,SD"),
        lines, size="\\scriptsize", sep="1.6pt", pos="!htb",
        note=(SITE_NOTE % ("MWh, rounded to integers", " (8 comparisons)"))
        + " Official IEA37 Gaussian-wake model and wind rose (16 directions, one wind speed), circular boundary, minimum spacing "
          "$2D$. Published layouts~\\cite{Baker2019,IEA37repo}: example layout, best submitted layout with every turbine within "
          "1~mm of the boundary (strict) and after radial projection onto it. Run-level Friedman tests over the seeds: %s." % frs))


def log_site(tag, b, st):
    log(f"[{tag}] {b}: GA {st['ga_status']}; best {st['best']}; order " + " > ".join(
        f"{LAB[a]}({st['methods'][a]['feasible']}, {st['methods'][a]['mean'] if st['methods'][a]['mean'] is None else round(st['methods'][a]['mean'], 2)})"
        for a in st["ranking_order"]))
    log("      PSO-VNS vs: " + ", ".join(f"{LAB[a]} {v['psovns_vs']['p_holm']:.2g}{v['psovns_vs']['outcome']}"
                                         for a, v in st["methods"].items() if "psovns_vs" in v))


# ====================================================================== cross-checks
PRINTED = {   # values printed in the manuscript (SWEVO_rev2/latex_source/optA/sw/06_results.tex, 08_beyond.tex)
    "tab:friedman68 avg rank (8 methods)": {"PSOBV": 1.82, "PSOC": 1.97, "SSABV": 3.45, "BVNS": 4.40, "SSA": 5.38,
                                             "SLSQP": 5.88, "LXSSA": 6.06, "DE": 7.04},
    "tab:ga-main avg rank (9 methods)": {"PSOBV": 2.18, "PSOC": 2.40, "GA": 2.74, "SSABV": 4.32, "BVNS": 5.21, "SSA": 6.34,
                                         "SLSQP": 6.85, "LXSSA": 7.01, "DE": 7.96},
    "tab:ga-main feasible % (9 methods)": {"PSOBV": 100.0, "PSOC": 100.0, "GA": 98.6, "SSABV": 99.6, "BVNS": 100.0, "SSA": 97.6,
                                           "SLSQP": 100.0, "LXSSA": 96.9, "DE": 71.4},
    "tab:hr-site AEP 6k R": {"GA": 138.12, "PSOBV": 138.29, "SSABV": 137.59, "BVNS": 137.35, "PSOC": 137.79, "SSA": 136.25,
                             "LXSSA": 136.00, "SLSQP": 136.31},
    "tab:hr-site loss 6k R": {"GA": 7.18, "PSOBV": 7.07, "SSABV": 7.54, "BVNS": 7.70, "PSOC": 7.40, "SSA": 8.44, "LXSSA": 8.61,
                              "SLSQP": 8.40},
    "tab:hr-site loss 30k": {"GA": 6.16, "PSOBV": 6.52, "SSABV": 7.19, "BVNS": 6.85, "PSOC": 6.78, "SSA": 8.17, "LXSSA": 8.50,
                             "DE": 7.40, "SLSQP": 7.82},
    "tab:iea-means 16T_6030": {"PSOBV": 395031.8, "PSOC": 390540.2, "GA": 391037.6, "SSABV": 388562.9, "SSA": 378073.1,
                               "LXSSA": 373977.4, "DE": 338610.6, "BVNS": 387637.3, "SLSQP": 382605.2},
    "tab:iea-means 16T_30030": {"PSOBV": 402583.1, "PSOC": 401347.0, "GA": 401766.1, "SSABV": 396920.5, "SSA": 392875.6,
                                "LXSSA": 390267.4, "DE": 359952.3, "BVNS": 396192.0, "SLSQP": 394567.0},
    "tab:iea-means 36T_6030": {"PSOBV": 736832.4, "PSOC": 729188.5, "GA": 774971.9, "SSABV": 728870.3, "SSA": 717258.5,
                               "LXSSA": 709660.5, "BVNS": 741565.6, "SLSQP": 733342.0},
    "tab:iea-means 36T_30030": {"PSOBV": 807446.0, "PSOC": 802252.3, "GA": 823061.5, "SSABV": 791404.4, "SSA": 771298.0,
                                "LXSSA": 752030.3, "BVNS": 801797.4, "SLSQP": 798913.0},
    "tab:iea37 best 16T_6030": {"PSOBV": 402782.4, "PSOC": 403753.1, "SSABV": 403251.3, "SSA": 395150.2, "LXSSA": 393284.6,
                                "DE": 357391.3, "BVNS": 400169.3, "SLSQP": 405433.8, "GA": 399691.3},
    "Sec. lillgrund text": {"PSOBV feasible 6030": 17, "PSOBV feasible 30030": 30, "PSOC feasible 6030": 0,
                            "PSOC feasible 30030": 19, "PSOBV mean 6030": 113.07, "PSOBV mean 30030": 114.59},
}


def chk(name, ours, ref, tol):
    bad = {k: (ours.get(k), v) for k, v in ref.items() if ours.get(k) is None or abs(ours[k] - v) > tol}
    r = dict(check=name, n=len(ref), tolerance=tol, ok=not bad, mismatches={k: dict(ours=o, ref=v) for k, (o, v) in bad.items()})
    log(f"  [{'OK' if r['ok'] else 'MISMATCH'}] {name}" + ("" if r["ok"] else f": {r['mismatches']}"))
    return r


def cross_checks(out, G9):
    log("\n==== cross-checks")
    C = []
    T = out["tables"]
    # ---- benchmark, 8-method pool (the paper's main comparison), recomputed from the same runs
    G8 = G9[G9.Algorithm.isin(MAIN8)]
    S8 = MR.case_stats(G8, MAIN8)
    FR8 = MR.friedman_block(MR.rank_matrix(S8, MAIN8), MAIN8, FOCUS)
    o8 = [a for a in MAIN8 if a != FOCUS]
    CW8 = MR.case_mean_wilcoxon(S8, FOCUS, o8)
    WT8 = R2.wtl_block(G8, FOCUS, o8)
    eight = dict(avg_rank=FR8["avg_rank"], chi2=FR8["chi2"], p=FR8["p"], iman_davenport=FR8["iman_davenport"],
                 p_wilcoxon_holm={b: CW8[b]["p_holm"] for b in o8}, mean_dloss_pp={b: CW8[b]["mean_dloss_pp"] for b in o8},
                 wtl={b: {k: WT8[b][k] for k in "WTL"} for b in o8},
                 mean_loss_pct_qualified={a: float(S8[(S8.Algorithm == a) & S8.Qualified].Loss.mean()) for a in MAIN8})
    out["eight_method_pool_benchmark"] = eight
    M = json.load(open(os.path.join(HERE, "mpce_summary.json")))["main"]
    C.append(chk("8-method avg ranks = mpce_summary.json main.friedman.avg_rank", FR8["avg_rank"], M["friedman"]["avg_rank"], 1e-9))
    C.append(chk("8-method Friedman chi2 / Iman-Davenport = mpce_summary.json", dict(chi2=FR8["chi2"], ff=FR8["iman_davenport"]),
                 dict(chi2=M["friedman"]["chi2"], ff=M["friedman"]["iman_davenport"]), 1e-6))
    C.append(chk("8-method case-mean Wilcoxon p_holm = mpce_summary.json", eight["p_wilcoxon_holm"],
                 {b: v["p_holm"] for b, v in M["case_mean_wilcoxon"].items()}, 1e-12))
    C.append(chk("8-method mean dL (pp) = mpce_summary.json", eight["mean_dloss_pp"],
                 {b: v["mean_dloss_pp"] for b, v in M["case_mean_wilcoxon"].items()}, 1e-9))
    C.append(chk("8-method run-level W/T/L = mpce_summary.json main.wtl (tab:wtl)",
                 {f"{b}_{k}": float(v[k]) for b, v in eight["wtl"].items() for k in "WTL"},
                 {f"{b}_{k}": float(v.get(k, 0)) for b, v in M["wtl"].items() for k in "WTL"}, 0))
    C.append(chk("8-method mean loss (qualified cases) = mpce_summary.json", eight["mean_loss_pct_qualified"],
                 M["mean_loss_pct_qualified_cases"], 1e-9))
    C.append(chk("printed tab:friedman68 avg ranks (2 d.p.)", {a: round(v, 2) for a, v in FR8["avg_rank"].items()},
                 PRINTED["tab:friedman68 avg rank (8 methods)"], 1e-9))
    # ---- benchmark, 9-method pool
    f9 = T["tab:farm-I-stats"]
    if os.path.exists(REV2_SUMMARY):
        g = json.load(open(REV2_SUMMARY))["ga68"]
        C.append(chk("9-method avg ranks = rev2_summary.json ga68.friedman.avg_rank", f9["avg_rank"], g["friedman"]["avg_rank"], 1e-9))
        C.append(chk("9-method Friedman chi2 = rev2_summary.json", dict(chi2=f9["friedman"]["chi2"]), dict(chi2=g["friedman"]["chi2"]), 1e-6))
        C.append(chk("9-method feasible % = rev2_summary.json", {a: v["feasible_pct"] for a, v in T["tab:farm-I-res"]["methods"].items()},
                     g["feasible_pct"], 1e-9))
        C.append(chk("9-method mean loss (qualified) = rev2_summary.json",
                     {a: v["mean_loss_pct_qualified"] for a, v in T["tab:farm-I-res"]["methods"].items()}, g["mean_loss_pct_qualified"], 1e-9))
    C.append(chk("printed tab:ga-main avg ranks (2 d.p.)", {a: round(v, 2) for a, v in f9["avg_rank"].items()},
                 PRINTED["tab:ga-main avg rank (9 methods)"], 1e-9))
    C.append(chk("printed tab:ga-main feasible % (1 d.p.)", {a: round(v["feasible_pct"], 1) for a, v in T["tab:farm-I-res"]["methods"].items()},
                 PRINTED["tab:ga-main feasible % (9 methods)"], 1e-9))
    # ---- Horns Rev 1: rev3_sites.json hr_pool (11-method pool; per-method statistics are pool-independent, raw p identical)
    J = json.load(open(os.path.join(HERE, "rev3_sites.json")))
    hr = T["tab:farm-II"]["budgets"]
    for b in BUDS:
        ref = J["hr_pool"]["budgets"][str(b)]
        for key in ("mean", "sd", "best", "loss_pct", "feasible"):
            C.append(chk(f"HR {b} {key} = rev3_sites.json hr_pool", {a: hr[str(b)]["methods"][a][key] for a in POOL
                                                                      if hr[str(b)]["methods"][a].get(key) is not None},
                         {a: ref["methods"][a][key] for a in POOL if ref["methods"][a][key] is not None}, 1e-9))
        C.append(chk(f"HR {b} raw paired p (PSO-VNS vs) = rev3_sites.json hr_pool psovns_vs_all",
                     {a: hr[str(b)]["methods"][a]["psovns_vs"]["p"] for a in POOL if a != FOCUS},
                     {a: ref["families"]["psovns_vs_all"]["tests"][a]["p"] for a in POOL if a != FOCUS}, 1e-12))
    C.append(chk("printed tab:hr-site AEP 6k (2 d.p.)", {a: round(v["mean"], 2) for a, v in hr["6030"]["methods"].items() if v["mean"]},
                 PRINTED["tab:hr-site AEP 6k R"], 1e-9))
    C.append(chk("printed tab:hr-site loss 6k (2 d.p.)", {a: round(v["loss_pct"], 2) for a, v in hr["6030"]["methods"].items() if v["loss_pct"]},
                 PRINTED["tab:hr-site loss 6k R"], 1e-9))
    C.append(chk("printed tab:hr-site loss 30k (2 d.p.)", {a: round(v["loss_pct"], 2) for a, v in hr["30030"]["methods"].items() if v["loss_pct"]},
                 PRINTED["tab:hr-site loss 30k"], 1e-9))
    # ---- IEA37: rev3_sites.json iea_pool (13-method pool)
    ie = T["tab:farm-IV"]["settings"]
    for k, st in ie.items():
        ref = J["iea_pool"]["settings"][k]
        for key in ("mean", "best", "feasible"):
            C.append(chk(f"IEA {k} {key} = rev3_sites.json iea_pool", {a: v[key] for a, v in st["methods"].items() if v.get(key) is not None},
                         {a: ref["methods"][a][key] for a in POOL if ref["methods"][a][key] is not None}, 1e-6))
        C.append(chk(f"IEA {k} raw paired p = rev3_sites.json iea_pool psovns_vs_all",
                     {a: st["methods"][a]["psovns_vs"]["p"] for a in POOL if a != FOCUS},
                     {a: ref["families"]["psovns_vs_all"]["tests"][a]["p"] for a in POOL if a != FOCUS}, 1e-12))
        C.append(chk(f"printed tab:iea-means {k} (1 d.p.)", {a: round(v["mean"], 1) for a, v in st["methods"].items() if v.get("mean")},
                     PRINTED[f"tab:iea-means {k}"], 1e-6))
    C.append(chk("printed tab:iea37 best 16T_6030 (1 d.p.)", {a: round(v["best"], 1) for a, v in ie["16T_6030"]["methods"].items() if v.get("best")},
                 PRINTED["tab:iea37 best 16T_6030"], 1e-6))
    # ---- Lillgrund: rev2_summary.json lillgrund (pool SP9 incl. RSD-VNS); p_holm reproduced with that family
    lgt = T["tab:farm-III"]["budgets"]
    if os.path.exists(REV2_SUMMARY):
        ref = json.load(open(REV2_SUMMARY))["lillgrund"]["budgets"]
        L = load_lg({})
        for b in BUDS:
            rm = ref[str(b)]["methods"]
            for key in ("mean", "sd", "best", "loss_pct", "feasible"):
                C.append(chk(f"LG {b} {key} = rev2_summary.json lillgrund", {a: v[key] for a, v in lgt[str(b)]["methods"].items()
                                                                              if v.get(key) is not None and a in rm},
                             {a: rm[a][key] for a in POOL if a in rm and rm[a][key] is not None}, 1e-9))
            y = L[(L.Budget == b) & (L.Init == "random") & L.Algorithm.isin(R2.SP9)]
            y, _ = R2.pair_seeds(y, CASE, R2.SP9)
            tt = {t["Baseline"]: t["PHolm"] for t in MR.paired_vs(y, FOCUS, [a for a in R2.SP9 if a != FOCUS])}
            C.append(chk(f"LG {b} p_holm with the rev2 family (8 comparisons incl. RSD-VNS) = rev2_summary.json", tt,
                         {a: v["p_holm_psovns_vs"] for a, v in rm.items() if v["p_holm_psovns_vs"] is not None}, 1e-12))
    o = {"PSOBV feasible 6030": lgt["6030"]["methods"]["PSOBV"]["feasible"], "PSOBV feasible 30030": lgt["30030"]["methods"]["PSOBV"]["feasible"],
         "PSOC feasible 6030": lgt["6030"]["methods"]["PSOC"]["feasible"], "PSOC feasible 30030": lgt["30030"]["methods"]["PSOC"]["feasible"],
         "PSOBV mean 6030": round(lgt["6030"]["methods"]["PSOBV"]["mean"], 2), "PSOBV mean 30030": round(lgt["30030"]["methods"]["PSOBV"]["mean"], 2)}
    C.append(chk("printed Sec. lillgrund numbers", o, PRINTED["Sec. lillgrund text"], 1e-9))
    sc = {f"{a} {b}": round(lgt[str(b)]["methods"][a]["psovns_vs"]["all_run_endpoint"]["score"], 3) for a in ("PSOC", "SLSQP") for b in BUDS}
    C.append(chk("printed Sec. archive-reliability Lillgrund all-run scores (3 d.p.)", sc,
                 {"PSOC 6030": 0.783, "PSOC 30030": 0.933, "SLSQP 6030": 0.333, "SLSQP 30030": 0.933}, 1e-9))
    out["cross_checks"] = dict(all_ok=all(c["ok"] for c in C), n_checks=len(C), n_failed=sum(not c["ok"] for c in C), checks=C)
    log(f"  cross-checks: {len(C) - out['cross_checks']['n_failed']}/{len(C)} OK")


# ====================================================================== main
def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", default=HERE)
    args = ap.parse_args(argv)
    os.makedirs(args.out_dir, exist_ok=True)
    t0 = time.time()
    out = dict(generated=time.strftime("%Y-%m-%d %H:%M:%S"), script="rev4_farm_tables.py",
               pool=POOL, labels=LAB, focus=FOCUS,
               protocol=dict(ranking=MR.RANK_RULE, run_level_score="mpce_results.goodness (objective if feasible; every "
                             "infeasible run below every feasible run, infeasible runs ordered by minimum spacing)",
                             paired_test="mpce_results.paired_vs (two-sided Wilcoxon signed-rank on seed-paired score "
                                         "differences, |d| <= 1e-9 dropped), Holm within the family of PSO-VNS comparisons",
                             all_run_endpoint="feasible beats infeasible; both feasible: larger objective wins (tie if "
                                              "|d| <= 1e-9); both infeasible: tie; score (W + T/2)/n; exact sign test on W vs L",
                             seeds="1-30, stored runs only; no optimization rerun", initialization="random"),
               data_status={}, tables={})
    tex = ["%% generated by rev4_farm_tables.py from the stored run records -- do not edit by hand\n"
           "%% requires booktabs; labels tab:farm-I-res, tab:farm-I-stats, tab:farm-II, tab:farm-III, tab:farm-IV\n"]
    log("==== Type I: benchmark")
    G9 = benchmark(out, tex)
    log("\n==== Type II: Horns Rev 1")
    block_hr(out, tex)
    log("\n==== Type III: Lillgrund")
    block_lg(out, tex)
    log("\n==== Type IV: IEA37")
    block_iea(out, tex)
    cross_checks(out, G9)
    out["runtime_s"] = time.time() - t0
    json.dump(R3.clean(out), open(os.path.join(args.out_dir, "rev4_farm_tables.json"), "w"), indent=1)
    open(os.path.join(args.out_dir, "rev4_farm_tables.tex"), "w").write("\n".join(tex))
    open(os.path.join(args.out_dir, "rev4_farm_tables.log"), "w").write("\n".join(LOG) + "\n")
    log(f"done in {time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()

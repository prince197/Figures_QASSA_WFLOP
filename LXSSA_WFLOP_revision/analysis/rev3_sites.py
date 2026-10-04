"""Revision 3, topic `sites` (items B2, B6, C3 of the remaining-work report; reviewer items R3, R8, R9).

Usage (from analysis/):  python3 rev3_sites.py [--only hr,iea,eval,lg,time,init] [--recompute-eval] [--recompute-init]
                                               [--out-dir D]
Outputs (default: this folder): rev3_sites.json, rev3_sites_tables.tex, rev3_sites_reeval.csv (per-layout
re-evaluations, cached), rev3_sites_inittime.csv (initialization timings, cached), rev3_sites.log (stdout copy).

Blocks
  hr    B2(a)  Horns Rev 1 16-turbine block, ONE pool = the ten methods of tab:hr-site (mpce_hrfix) + GA (rev2_gahr),
               random starts, 6,030 and 30,030 evaluations, seeds 1-30. Per budget: feasible runs, mean / best AEP and
               AEP loss of the feasible runs, feasibility-aware rank of the mean AEP (mpce_results.rank_rule: a method
               with fewer than 15 feasible runs ranks below the qualified ones, by its feasible runs), run-level
               Friedman test over the 30 seeds (blocks = seeds; score = mpce_results.goodness: AEP of a feasible run,
               every infeasible run below every feasible run, infeasible runs ordered by minimum spacing) and the mean
               seed-level rank of each method, seed-paired Wilcoxon signed-rank tests (mpce_results.paired_vs, same
               score). Families (Holm within each, two-sided alpha = 0.05):
                 F_best(b) = the best method (feasibility-aware rank 1) vs each of the other 10 methods at budget b;
                 F_PV(b)   = PSO-VNS vs each of the other 10 methods at budget b (identical to F_best(b) if PSO-VNS is best).
  iea   B2(b)  IEA37 Case Study 1, 16 and 36 turbines, 6,030 and 30,030 evaluations: ONE pool of 13 methods = the ten
               methods (mpce_iea16/36 + PSO-VNS arms mpce_iea16p/36p) + GA (rev2_gaiea) + MS-SLSQP and PSO-SLSQP with the
               exact analytic gradient (rev2_grad). Same statistics per scenario x budget; families (Holm within each):
               F_best(n, b) = best method vs the other 12; F_PV(n, b) = PSO-VNS vs the other 12 (reported in the JSON);
               F_free(n, b) = best gradient-free method vs the other gradient-free methods (10 comparisons; JSON).
               Feasibility = the stored labels (rounded coordinates cannot re-verify them; supplement S11).
  eval  C3     Horns Rev 1: every stored layout of mpce_hrfix (970) and rev2_gahr (60) re-evaluated (not re-optimized)
               with the paper's 5-deg bins (replay of the stored AEP), 1-deg bins centred at 0.5, 1.5, ... (primary, as
               mpce_direction.py), the phase checks 1 deg at integer centres and 5 deg at 0, 5, ..., and PyWake 2.6.20
               NOJ (k = 0.04, defaults) at the 1-deg bins (as pywake_check.py --layouts). Tests as Table S44 / tab:F-hr /
               tab:F-hrpair: mean differences on jointly feasible seeds with the unadjusted Wilcoxon p, and the paper's
               run-level test (infeasible below feasible) with Holm over the 10 comparisons of GA.
  lg    C3     Lillgrund 16-turbine block (rev2_lg16 / rev2_lg16b, 540 layouts): replay at 5 deg, 1-deg bins
               (rev2_site_model.bins(1, 0.5)) and PyWake NOJ at 1 deg; PSO-VNS vs each method (paper's run-level test,
               Holm over the 8 comparisons per budget) and the highest mean AEP under each evaluator.
  time  B6     Elapsed time from the Seconds column of the original records (benchmark main comparison and controls,
               budget, feasibility-preserving initialization, Horns Rev 1, IEA37) and of the revision records: median
               seconds per run by method and N / budget, ratios to PSO-VNS within a study, ms per evaluation (median and
               mean of Seconds / Calls; the latter is the statistic of Table S-cost), total CPU time, shard wall times
               from the run logs.
  init  B6     Wall time of the feasibility-preserving initialization (feasible_init.make_generator: L-BFGS-B packing,
               30 layouts per run, exactly as drawn by the optimizers after np.random.seed(seed)) for the six largest
               benchmark cases and the Horns Rev 1 block, seeds 1-30, measured now on this machine (one process).
               Check: the regenerated population contains the stored final layout of the feasible-start PSO run of
               the same seed (PSO never leaves the best initial layout from feasible starts).

All statistics reuse mpce_results.py (MR) and rev2_analysis.py (R2) functions; no existing file is modified.
"""
import os, sys, glob, re, json, time, math, argparse, warnings
import numpy as np, pandas as pd
from scipy.stats import friedmanchisquare, rankdata, wilcoxon, f as f_dist

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import mpce_results as MR          # noqa: E402
import rev2_analysis as R2         # noqa: E402
import mpce_inference_extra as MX  # noqa: E402

warnings.filterwarnings("ignore", category=RuntimeWarning)
CASE = MR.CASE
HR10 = list(MR.HR_ORDER)                     # the ten methods of tab:hr-site
HR11 = HR10 + ["GA"]
GRADM = ["SLSQPX", "PSOSLSQPX"]
IEA13 = HR10 + ["GA"] + GRADM
FREE = HR10 + ["GA"]                         # gradient-free (MS-SLSQP with forward differences counted here)
LAB = dict(R2.LAB)
LAB.update(SLSQPX="MS-SLSQP, exact grad.", PSOSLSQPX="PSO-SLSQP, exact grad.", RSDVNS="RSD-VNS")
LAB_TXT = {k: v.replace("\\%", "%").replace("$", "") for k, v in LAB.items()}
BUDS = (6030, 30030)
LOG = []


def log(msg=""):
    print(msg, flush=True)
    LOG.append(str(msg))


def fp(p):
    return "--" if p is None or not np.isfinite(p) else MR.fmt_p(float(p), 2)


def big(v, d=1):
    return "--" if v is None or not np.isfinite(v) else f"{v:,.{d}f}".replace(",", "{,}")


# ====================================================================== run-level pool analysis
def pool_analysis(y, methods, families):
    """Run-level analysis of ONE pool on one site x budget (single case).
    families: dict name -> (focus or callable(best) -> focus, list of comparators or None = all others)."""
    y = y[y.Algorithm.isin(methods)]
    mets = [a for a in methods if a in set(y.Algorithm)]
    y, drop = R2.pair_seeds(y, CASE, mets)
    S = MR.case_stats(y, mets).set_index("Algorithm")
    piv = y.assign(Sc=MR.goodness(y)).pivot_table(index="Seed", columns="Algorithm", values="Sc")[mets].dropna()
    n, k = piv.shape
    chi, pfr = friedmanchisquare(*[piv[a].values for a in mets])
    ff = (n - 1) * chi / (n * (k - 1) - chi)
    pff = float(f_dist.sf(ff, k - 1, (k - 1) * (n - 1)))
    srank = pd.DataFrame(np.vstack([rankdata(-r) for r in piv.values]), index=piv.index, columns=mets).mean()
    # best method: feasibility-aware rank 1 (ties -> higher mean AEP)
    best = sorted(mets, key=lambda a: (S.Rank[a], -(S.Mean[a] if np.isfinite(S.Mean[a]) else -1e18)))[0]
    res = dict(methods={}, n_seeds=int(n), n_methods=int(k), dropped_unpaired_runs=int(drop), best=best,
               friedman=dict(chi2=float(chi), df=int(k - 1), p=float(pfr), iman_davenport=float(ff), iman_davenport_p=pff,
                             blocks="seeds (run level)", score="mpce_results.goodness (infeasible below feasible)"),
               families={})
    for a in mets:
        nf = int(S.NFeas[a])
        res["methods"][a] = dict(runs=int(S.N[a]), feasible=nf,
                                 mean=float(S.Mean[a]) if nf else None, sd=float(S.SD[a]) if nf > 1 else None,
                                 best=float(S.Best[a]) if nf else None, loss_pct=float(S.Loss[a]) if nf else None,
                                 rank=float(S.Rank[a]), qualified=bool(S.Qualified[a]), seed_rank=float(srank[a]),
                                 median_seconds=float(y[y.Algorithm == a].Seconds.median()) if "Seconds" in y else None)
    for name, (foc, comp) in families.items():
        f = foc(best) if callable(foc) else foc
        if f not in mets:
            continue
        others = [b for b in (comp(best) if callable(comp) else (comp or mets)) if b != f and b in mets]
        rows = MR.paired_vs(y, f, others)
        fam = dict(focus=f, comparisons=len(rows), correction="Holm within this family", tests={})
        for r in rows:
            d = (piv[f] - piv[r["Baseline"]]).values
            fam["tests"][r["Baseline"]] = dict(p=r["P"], p_holm=r["PHolm"], rb=r["RB"], outcome=r["Outcome"],
                                               n_pairs=r["NPairs"], wins=int((d > 1e-9).sum()), losses=int((d < -1e-9).sum()),
                                               ties=int((np.abs(d) <= 1e-9).sum()))
        res["families"][name] = fam
    return res


# ====================================================================== loaders
def raw_shards(pattern, cols=None):
    fs = sorted(glob.glob(os.path.join(HERE, pattern)))
    return pd.concat([pd.read_csv(f, usecols=cols) for f in fs], ignore_index=True) if fs else None


def load_hr():
    H = R2.read_mpce("hrfix")
    H = H[(H.Dataset == "HR") & (H.Turbines == 16)]
    G, _ = R2.read_rev2("gahr", HERE)
    G = G[(G.Algorithm == "GA") & (G.Dataset == "HR")]
    return pd.concat([H, G], ignore_index=True)


def load_iea():
    parts = [R2.read_mpce(e) for e in ("iea16", "iea36", "iea16p", "iea36p")]
    g1, _ = R2.read_rev2("gaiea", HERE)
    g2, _ = R2.read_rev2("grad", HERE)
    IE = pd.concat([p for p in parts if p is not None] + [g1[g1.Algorithm == "GA"], g2[g2.Algorithm.isin(GRADM)]],
                   ignore_index=True)
    return IE[IE.Dataset.str.startswith("IEA37")]


# ====================================================================== B2(a): Horns Rev 1 pool
def block_hr(summ, tex):
    HH = load_hr()
    import hornsrev_model as hrm
    inst = float(hrm.aep_gwh(hrm.site(16)[0]))
    ideal = float(HH.Ideal.iloc[0])
    out = dict(pool=HR11, pool_note="ten methods of tab:hr-site (mpce_hrfix) + GA (rev2_gahr); random starts; seeds 1-30",
               installed_aep=inst, ideal_aep=ideal, budgets={})
    for b in BUDS:
        y = HH[(HH.Budget == b) & (HH.Init == "random")]
        r = pool_analysis(y, HR11, {"best_vs_all": (lambda best: best, None), "psovns_vs_all": ("PSOBV", None)})
        for a, v in r["methods"].items():
            yy = y[(y.Algorithm == a) & y.Feasible]
            v["runs_above_installed"] = int((yy.Objective > inst).sum())
        out["budgets"][b] = r
        log(f"[hr] {b}: best {r['best']}; Friedman p={r['friedman']['p']:.3g}; " + ", ".join(
            f"{a} {v['mean'] if v['mean'] is None else round(v['mean'], 2)} ({v['feasible']}) r{v['rank']:g}/{v['seed_rank']:.2f}"
            for a, v in sorted(r["methods"].items(), key=lambda t: t[1]["rank"])))
        for fn, fam in r["families"].items():
            log(f"      {fn} ({fam['focus']}): " + ", ".join(f"{k} {t['p_holm']:.2g}{t['outcome']}" for k, t in fam["tests"].items()))
    summ["hr_pool"] = out
    # ---- LaTeX: rows ordered by the feasibility-aware rank at 30,030 evaluations
    B6, B30 = out["budgets"][6030], out["budgets"][30030]
    order = sorted(HR11, key=lambda a: (B30["methods"][a]["rank"], B6["methods"][a]["rank"]))
    two6 = B6["best"] != "PSOBV"
    two30 = B30["best"] != "PSOBV"
    il = f"{100 * (1 - inst / ideal):.2f}"
    lines = ["Installed block & " + " & ".join(["--", f"{inst:.2f}", il, "--", "--", "--"] + (["--"] if two6 else [])) + " & & "
             + " & ".join(["--", f"{inst:.2f}", il, "--", "--", "--"] + (["--"] if two30 else [])) + " \\\\", "\\midrule"]

    def cells(r, a, two):
        v = r["methods"][a]
        best, fb, fpv = r["best"], r["families"]["best_vs_all"], r["families"]["psovns_vs_all"]
        mean = "--" if v["mean"] is None else f"{v['mean']:.2f}"
        loss = "--" if v["loss_pct"] is None else f"{v['loss_pct']:.2f}"
        if a == best:
            mean, loss = f"\\textbf{{{mean}}}", f"\\textbf{{{loss}}}"
        pb = "best" if a == best else (fp(fb["tests"][a]["p_holm"]) + ("" if fb["tests"][a]["outcome"] == "T" else
                                                                        "$^{\\ast}$"))
        c = [f"{v['feasible']}/{v['runs']}", mean, loss, f"{v['rank']:g}", f"{v['seed_rank']:.2f}", pb]
        if two:
            c.append("--" if a == "PSOBV" else fp(fpv["tests"][a]["p_holm"]))
        return c
    for a in order:
        lines.append(f"{LAB[a]} & " + " & ".join(cells(B6, a, two6)) + " & & " + " & ".join(cells(B30, a, two30)) + " \\\\")
    n6, n30 = (7 if two6 else 6), (7 if two30 else 6)
    hdr6 = "Feas. & AEP & Loss & Rank & $\\bar R_s$ & $p_{\\rm best}$" + (" & $p_{\\rm PV}$" if two6 else "")
    hdr30 = "Feas. & AEP & Loss & Rank & $\\bar R_s$ & $p_{\\rm best}$" + (" & $p_{\\rm PV}$" if two30 else "")
    head = (f"& \\multicolumn{{{n6}}}{{c}}{{6{{,}}030 evaluations}} & & \\multicolumn{{{n30}}}{{c}}{{30{{,}}030 evaluations}} \\\\\n"
            f"\\cmidrule(lr){{2-{1 + n6}}}\\cmidrule(lr){{{3 + n6}-{2 + n6 + n30}}}\n"
            f"Method & {hdr6} & & {hdr30}")
    fr = lambda r: f"$\\chi^2_F={r['friedman']['chi2']:.1f}$ ({r['friedman']['df']} d.f.), $p={fp(r['friedman']['p']).strip('$')}$"
    note = ("Pool: the ten methods of Table~\\ref{M-tab:hr-site} and GA (11 methods), random initialization, 30 seed-paired runs "
            "(seeds 1--30) per method and budget, 5$^\\circ$ direction bins (the evaluator used in optimization). "
            f"Wake-free AEP {ideal:.2f}~GWh/yr. Feas.: feasible runs; AEP: mean (GWh/yr) and Loss: mean AEP loss (\\%) of the "
            "feasible runs (bold: best method); Rank: feasibility-aware rank of the mean AEP within this pool (rule of "
            "Table~\\ref{M-tab:friedman68}: fewer than 15 feasible runs ranks below the qualified methods, by feasible runs); "
            "$\\bar R_s$: mean seed-level rank (run-level Friedman ranks over the 30 seeds, infeasible runs below every feasible "
            "run, two infeasible runs ordered by minimum spacing; 1 = best). $p_{\\rm best}$: seed-paired Wilcoxon signed-rank "
            "test (same ordering) of the best method (Rank 1 at that budget) against the method, Holm-adjusted within the family of "
            "its 10 comparisons at that budget; $^{\\ast}$: significant at 0.05 (in every such case the best method is better)"
            + ("; $p_{\\rm PV}$: the same for PSO-VNS against the method (second family of 10 comparisons at that budget)" if (two6 or two30) else "")
            + f". Run-level Friedman test: 6{{,}}030: {fr(B6)}; 30{{,}}030: {fr(B30)}.")
    tex.append(MR.table("table*", "Horns Rev~1 16-turbine block, recomputed pool with GA (11 methods; 5$^\\circ$ bins, random "
                        "starts, 30 seeds, 6{,}030 and 30{,}030 evaluations): feasibility, mean AEP of the feasible runs, "
                        "feasibility-aware and run-level ranks, and seed-paired tests of the best method.",
                        "tab:S-r3-hr-pool", "l" + "c" * n6 + "c" + "c" * n30, head, lines, sep="1.7pt", pos="!htb", note=note))


# ====================================================================== B2(b): IEA37 pool
def published_iea():
    pub = MR.load_published(HERE)
    pv = {}
    for n in (16, 36):
        P = pub[pub.Turbines == n]
        pf = P[(~P.Baseline) & P.Feasible]
        pv[n] = dict(example=float(P[P.Baseline].AEP.iloc[0]), strict=float(pf.AEP.max()))
    J = json.load(open(os.path.join(HERE, "iea37_projected.json")))
    for n in (16, 36):
        pv[n]["projected"] = float(J["scenarios"][str(n)]["projected"]["best_aep"])
    return pv


def block_iea(summ, tex):
    IE = load_iea()
    pv = published_iea()
    out = dict(pool=IEA13, pool_note="ten methods of tab:iea37-detail (mpce_iea16/36 + mpce_iea16p/36p) + GA (rev2_gaiea) "
               "+ MS-SLSQP and PSO-SLSQP with exact analytic gradient (rev2_grad, c_g = 3); random starts; seeds 1-30",
               published={str(n): v for n, v in pv.items()}, settings={})
    for n in (16, 36):
        for b in BUDS:
            y = IE[(IE.Turbines == n) & (IE.Budget == b) & (IE.Init == "random")]
            r = pool_analysis(y, IEA13, {
                "best_vs_all": (lambda best: best, None), "psovns_vs_all": ("PSOBV", None),
                "best_free_vs_free": (None, None)})
            # best gradient-free method and its family (computed separately: needs the ranks)
            free = [a for a in FREE if a in r["methods"]]
            bf = sorted(free, key=lambda a: (r["methods"][a]["rank"], -(r["methods"][a]["mean"] or -1e18)))[0]
            r2 = pool_analysis(y, IEA13, {"best_free_vs_free": (bf, free)})
            r["families"]["best_free_vs_free"] = r2["families"]["best_free_vs_free"]
            r["best_gradient_free"] = bf
            for a, v in r["methods"].items():
                if v["best"] is not None:
                    v["best_vs_strict_pct"] = 100 * (v["best"] / pv[n]["strict"] - 1)
                    v["best_vs_projected_pct"] = 100 * (v["best"] / pv[n]["projected"] - 1)
            out["settings"][f"{n}T_{b}"] = r
            log(f"[iea] {n}T {b}: best {r['best']}, best free {bf}; Friedman p={r['friedman']['p']:.3g}; " + ", ".join(
                f"{a} {v['mean'] if v['mean'] is None else round(v['mean'], 1)} ({v['feasible']}) r{v['rank']:g}"
                for a, v in sorted(r["methods"].items(), key=lambda t: t[1]["rank"])))
            for fn, fam in r["families"].items():
                log(f"      {fn} ({fam['focus']}): " + ", ".join(f"{k} {t['p_holm']:.2g}{t['outcome']}" for k, t in fam["tests"].items()))
    summ["iea_pool"] = out
    # ---- LaTeX: two panels (16 / 36 turbines), columns per budget: Feas., Mean, Best, Rank, p_best
    lines = []
    for n in (16, 36):
        S6, S30 = out["settings"][f"{n}T_6030"], out["settings"][f"{n}T_30030"]
        lines.append(f"\\multicolumn{{11}}{{l}}{{\\emph{{{n} turbines ($r={1300 if n == 16 else 2000}$~m)}}}} \\\\")
        lines.append(f"Published & \\multicolumn{{10}}{{l}}{{example layout {big(pv[n]['example'], 0)}; best feasible, strict "
                     f"{big(pv[n]['strict'], 0)}; best after projection {big(pv[n]['projected'], 0)}}} \\\\")
        order = sorted(IEA13, key=lambda a: (S30["methods"][a]["rank"], S6["methods"][a]["rank"]))
        for a in order:
            c = []
            for r in (S6, S30):
                v = r["methods"][a]
                mean, bst = big(v["mean"], 0), big(v["best"], 0)
                if a == r["best"]:
                    mean = f"\\textbf{{{mean}}}"
                fam = r["families"]["best_vs_all"]
                pb = "best" if a == r["best"] else fp(fam["tests"][a]["p_holm"]) + ("" if fam["tests"][a]["outcome"] == "T" else "$^{\\ast}$")
                c += [f"{v['feasible']}/{v['runs']}", mean, bst, f"{v['rank']:g}", pb]
            lines.append(f"{LAB[a]} & " + " & ".join(c) + " \\\\")
        if n == 16:
            lines.append("\\midrule")
    head = ("& \\multicolumn{5}{c}{6{,}030 evaluations} & \\multicolumn{5}{c}{30{,}030 evaluations} \\\\\n"
            "\\cmidrule(lr){2-6}\\cmidrule(lr){7-11}\n"
            "Method & Feas. & Mean AEP & Best AEP & Rank & $p_{\\rm best}$ & Feas. & Mean AEP & Best AEP & Rank & $p_{\\rm best}$")
    frs = "; ".join(f"{n} turbines, {b:,}".replace(",", "{,}") + f": $p={fp(out['settings'][f'{n}T_{b}']['friedman']['p']).strip('$')}$"
                    for n in (16, 36) for b in BUDS)
    note = ("Pool (13 methods, one pool per scenario and budget): the ten methods of Table~\\ref{tab:iea37-detail}, GA, and "
            "MS-SLSQP and PSO-SLSQP with the exact analytic gradient (each gradient charged $c_g=3$ evaluations); random "
            "initialization, 30 seed-paired runs (seeds 1--30). AEP in MWh, rounded to integers (official IEA37 Gaussian-wake model); Mean: over "
            "the feasible runs (bold: best method); Best: best feasible run; Feas.: feasible runs as originally labelled (the "
            "rounded stored coordinates cannot re-verify the strict $10^{-6}$~m labels; Section~\\ref{sec:S-archive}). Rank: "
            "feasibility-aware rank of the mean AEP within this pool (rule of Table~\\ref{M-tab:friedman68}). $p_{\\rm best}$: "
            "seed-paired Wilcoxon signed-rank test of the best method (Rank 1) against the method (infeasible runs below every "
            "feasible run), Holm-adjusted within the family of its 12 comparisons in that scenario and budget; $^{\\ast}$: "
            "significant at 0.05 (the best method better in every such case). Published rows: example layout, best feasible "
            "submitted layout with every turbine within 1~mm of the boundary (strict) and after radial projection onto it "
            f"(projected)~\\cite{{Baker2019,IEA37repo}}. Run-level Friedman tests: {frs}.")
    tex.append(MR.table("table*", "IEA37 Case Study~1, recomputed pool with GA and the exact-gradient methods (13 methods; "
                        "16 and 36 turbines, random starts, 30 seeds, 6{,}030 and 30{,}030 evaluations): feasibility, mean and "
                        "best AEP, feasibility-aware rank and seed-paired tests of the best method.",
                        "tab:S-r3-iea-pool", "l" + "ccccc" * 2, head, lines, sep="1.5pt", pos="!htb", note=note))


# ====================================================================== C3: Horns Rev 1 alternate evaluators
HR_EV = ["5deg_2.5", "1deg_0.5", "1deg_0", "5deg_0", "PyWakeNOJ_1deg"]
EV_LAB = {"5deg_2.5": "5$^\\circ$", "1deg_0.5": "1$^\\circ$", "PyWakeNOJ_1deg": "PyWake NOJ 1$^\\circ$",
          "1deg_0": "1$^\\circ$ (int.)", "5deg_0": "5$^\\circ$ at 0$^\\circ$"}


def parse_xy(c):
    return np.array([[float(v) for v in p.split()] for p in c.split(";")])


def reeval_layouts(recompute, cache):
    """Per-layout re-evaluation of Horns Rev 1 (hrfix + GA) and Lillgrund (lg16, lg16b); cached CSV."""
    if os.path.exists(cache) and not recompute:
        return pd.read_csv(cache)
    import mpce_direction as D
    import hornsrev_model as hrm
    import rev2_site_model as lg
    import py_wake
    from py_wake import NOJ
    from py_wake.examples.data.hornsrev1 import Hornsrev1Site, V80, wt16_x, wt16_y
    from py_wake.examples.data import lillgrund as L
    rows = []
    t0 = time.time()
    # ---- Horns Rev 1
    cols = ["Algorithm", "Dataset", "Turbines", "Seed", "Budget", "Init", "Feasible", "Objective", "MinSpacing", "Seconds", "Coordinates"]
    X = pd.concat([raw_shards("mpce_hrfix_s*of24.csv", cols).assign(Source="mpce_hrfix"),
                   raw_shards("rev2_gahr_s*of*.csv", cols).assign(Source="rev2_gahr")], ignore_index=True)
    X = X[X.Dataset.astype(str) == "HR"].drop_duplicates(["Algorithm", "Seed", "Budget", "Init"], keep="last")
    B = {k: D.hr_bins(*v) for k, v in D.HR_BINS.items()}
    mdl = NOJ(Hornsrev1Site(), V80(), k=0.04)
    cx, cy = np.mean(wt16_x), np.mean(wt16_y)
    b1 = D.hr_bins(1.0, 0.5)
    xyi, _ = hrm.site(16)
    inst = {k: D.hr_aep(xyi, b) for k, b in B.items()}
    inst["PyWakeNOJ_1deg"] = float(mdl(np.asarray(wt16_x, float), np.asarray(wt16_y, float), wd=b1["wd"], ws=hrm.WS).aep().sum())
    rows.append(dict(Site="HR", Algorithm="INSTALLED", Seed=0, Budget=0, Init="--", Feasible=True, Objective=hrm.aep_gwh(xyi),
                     Source="hornsrev_model.site(16)", **inst))
    for i, r in enumerate(X.itertuples(index=False)):
        xy = parse_xy(r.Coordinates)
        e = {k: D.hr_aep(xy, b) for k, b in B.items()}
        e["PyWakeNOJ_1deg"] = float(mdl(xy[:, 0] + cx, xy[:, 1] + cy, wd=b1["wd"], ws=hrm.WS).aep().sum())
        rows.append(dict(Site="HR", Algorithm=r.Algorithm, Seed=r.Seed, Budget=r.Budget, Init=r.Init, Feasible=r.Feasible,
                         Objective=r.Objective, Source=r.Source, **e))
        if i % 200 == 0:
            log(f"  [eval] HR {i}/{len(X)} {time.time() - t0:.0f} s")
    # ---- Lillgrund
    Y = pd.concat([raw_shards("rev2_lg16_s*of*.csv", cols).assign(Source="rev2_lg16"),
                   raw_shards("rev2_lg16b_s*of*.csv", cols).assign(Source="rev2_lg16b")], ignore_index=True)
    Y = Y[Y.Dataset.astype(str) == "LG"].drop_duplicates(["Algorithm", "Seed", "Budget", "Init"], keep="last")
    wd1, f1, p1 = lg.bins(1.0, 0.5)
    lm = NOJ(L.LillgrundSite(), L.SWT23(), k=0.04)
    gx, gy = lg.WT_X[lg.I16].mean(), lg.WT_Y[lg.I16].mean()
    xyl, _ = lg.site(16)
    li = dict(**{"5deg_2.5": lg.aep_gwh(xyl), "1deg_0.5": lg.aep_gwh(xyl, True, wd1, f1, p1),
                 "PyWakeNOJ_1deg": float(lm(lg.WT_X[lg.I16], lg.WT_Y[lg.I16], wd=wd1, ws=lg.WS).aep().sum())})
    rows.append(dict(Site="LG", Algorithm="INSTALLED", Seed=0, Budget=0, Init="--", Feasible=True, Objective=li["5deg_2.5"],
                     Source="rev2_site_model.site(16)", **li))
    for i, r in enumerate(Y.itertuples(index=False)):
        xy = parse_xy(r.Coordinates)
        e = {"5deg_2.5": lg.aep_gwh(xy), "1deg_0.5": lg.aep_gwh(xy, True, wd1, f1, p1),
             "PyWakeNOJ_1deg": float(lm(xy[:, 0] + gx, xy[:, 1] + gy, wd=wd1, ws=lg.WS).aep().sum())}
        rows.append(dict(Site="LG", Algorithm=r.Algorithm, Seed=r.Seed, Budget=r.Budget, Init=r.Init, Feasible=r.Feasible,
                         Objective=r.Objective, Source=r.Source, **e))
        if i % 200 == 0:
            log(f"  [eval] LG {i}/{len(Y)} {time.time() - t0:.0f} s")
    out = pd.DataFrame(rows)
    out["Feasible"] = out.Feasible.astype(str).str.lower().isin(["true", "1", "1.0"])
    out.insert(0, "PyWakeVersion", py_wake.__version__)
    out.to_csv(cache, index=False, float_format="%.10g")
    log(f"  [eval] wrote {cache} ({len(out)} rows) in {time.time() - t0:.0f} s")
    return out


def jf_pair(s, a, b, c):
    """a - b on jointly feasible seeds (as tab:F-hrpair): n, mean difference, unadjusted Wilcoxon p, wins."""
    A_ = s[(s.Algorithm == a) & s.Feasible].set_index("Seed")[c]
    B_ = s[(s.Algorithm == b) & s.Feasible].set_index("Seed")[c]
    j = A_.index.intersection(B_.index)
    if len(j) < 2:
        return dict(n=int(len(j)))
    dd = (A_[j] - B_[j]).values
    p = float(wilcoxon(dd).pvalue) if (np.abs(dd) > 1e-12).any() else 1.0
    return dict(n=int(len(j)), mean_diff_gwh=float(dd.mean()), p=p, first_better=int((dd > 0).sum()), second_better=int((dd < 0).sum()))


def paper_tests(s, focus, others, c):
    Y = s.copy()
    Y["Objective"] = Y[c]
    Y["MinSpacing"] = Y.get("MinSpacing", pd.Series(np.nan, index=Y.index)).fillna(0.0)
    rows = MR.paired_vs(Y, focus, [b for b in others if b != focus and b in set(Y.Algorithm)])
    return {r["Baseline"]: dict(p_holm=r["PHolm"], rb=r["RB"], outcome=r["Outcome"], n_pairs=r["NPairs"]) for r in rows}


def add_minspacing(E, site):
    """Attach MinSpacing (needed by the paper's ordering of infeasible runs) from the stored records."""
    if site == "HR":
        Hm = load_hr()[["Algorithm", "Seed", "Budget", "Init", "MinSpacing"]]
    else:
        parts = [R2.read_rev2(e, HERE)[0] for e in ("lg16", "lg16b")]
        Hm = pd.concat(parts, ignore_index=True)[["Algorithm", "Seed", "Budget", "Init", "MinSpacing"]]
    return E.merge(Hm, on=["Algorithm", "Seed", "Budget", "Init"], how="left")


def block_eval(summ, tex, recompute):
    E = reeval_layouts(recompute, os.path.join(HERE, "rev3_sites_reeval.csv"))
    H = E[(E.Site == "HR")].copy()
    inst = H[H.Algorithm == "INSTALLED"].iloc[0]
    X = add_minspacing(H[H.Algorithm != "INSTALLED"].copy(), "HR")
    X["Seed"] = X.Seed.astype(int); X["Budget"] = X.Budget.astype(int)
    out = dict(installed={c: float(inst[c]) for c in HR_EV}, n_layouts=int(len(X)),
               note="5deg_2.5 = the paper's evaluator replayed on the stored (0.01-m rounded) coordinates")
    # reproduction of the stored AEP (5-deg replay) and of the published tab:F-hr values (ten methods)
    rep = {}
    for src, g in X[X.Feasible].groupby("Source"):
        d = g["5deg_2.5"] - g.Objective
        rep[src] = dict(n_feasible=int(len(g)), max_abs_dev_gwh=float(d.abs().max()), median_abs_dev_gwh=float(d.abs().median()),
                        n_abs_dev_gt_0_01=int((d.abs() > 0.01).sum()), n_higher_gt_0_01=int((d > 0.01).sum()),
                        mean_signed_dev_gwh=float(d.mean()))
        for b in BUDS:
            gb = g[(g.Budget == b) & (g.Init == "random")]
            if len(gb):
                db = gb["5deg_2.5"] - gb.Objective
                rep[src][f"{b}"] = dict(n=int(len(gb)), max_abs_dev_gwh=float(db.abs().max()), mean_signed_dev_gwh=float(db.mean()),
                                        mean_stored=float(gb.Objective.mean()), mean_replay=float(gb["5deg_2.5"].mean()))
    out["replay_5deg"] = rep
    pw_old = pd.read_csv(os.path.join(HERE, "pywake_check_hr16runs.csv"))
    m = X.merge(pw_old, on=["Algorithm", "Seed", "Budget", "Init"], suffixes=("", "_old"))
    out["pywake_vs_stored_pywake_check"] = dict(n=int(len(m)), max_abs_diff_gwh=float((m.PyWakeNOJ_1deg - m.PyWakeNOJ_1deg_old).abs().max()))
    log(f"[eval] HR replay: " + json.dumps({k: (v['n_feasible'], round(v['max_abs_dev_gwh'], 4), v['n_abs_dev_gt_0_01']) for k, v in rep.items()})
        + f"; PyWake vs pywake_check_hr16runs.csv max |diff| {out['pywake_vs_stored_pywake_check']['max_abs_diff_gwh']:.2e}")
    meth = {}
    for b in BUDS:
        s = X[(X.Budget == b) & (X.Init == "random")]
        for a in HR11:
            f = s[(s.Algorithm == a) & s.Feasible]
            e = dict(runs=int((s.Algorithm == a).sum()), feasible=int(len(f)))
            for c in HR_EV:
                e[f"mean_{c}"] = float(f[c].mean()) if len(f) else None
                e[f"best_{c}"] = float(f[c].max()) if len(f) else None
                e[f"above_{c}"] = int((f[c] > out["installed"][c]).sum())
            e["mean_drop_5to1_gwh"] = float((f["5deg_2.5"] - f["1deg_0.5"]).mean()) if len(f) else None
            meth[f"{b}_{a}"] = e
    out["methods"] = meth
    hi, pairs, ptest = {}, {}, {}
    for b in BUDS:
        s = X[(X.Budget == b) & (X.Init == "random")]
        for c in HR_EV:
            q = {a: meth[f"{b}_{a}"][f"mean_{c}"] for a in HR11
                 if meth[f"{b}_{a}"]["feasible"] >= 15 and meth[f"{b}_{a}"][f"mean_{c}"] is not None}
            oq = sorted(q, key=lambda a: -q[a])
            hi[f"{b}_{c}"] = dict(order_qualified=oq, best=oq[0], best_aep=q[oq[0]], second=oq[1], second_aep=q[oq[1]],
                                  next_best_excluding_ga=[a for a in oq if a != "GA"][0])
            # GA vs PSO-VNS and vs the next best method (excluding GA) under this evaluator, jointly feasible seeds
            for ref in ("PSOBV", hi[f"{b}_{c}"]["next_best_excluding_ga"]):
                pairs[f"{b}_GA_{ref}_{c}"] = jf_pair(s, "GA", ref, c)
            for ref in HR10:
                pairs.setdefault(f"{b}_GA_{ref}_{c}", jf_pair(s, "GA", ref, c))
            ptest[f"{b}_{c}"] = paper_tests(s, "GA", HR10, c)
    out["highest_mean"] = hi
    out["pairs_joint_feasible_GA_minus"] = pairs
    out["paper_test_GA_vs"] = ptest
    for b in BUDS:
        for c in ("5deg_2.5", "1deg_0.5", "PyWakeNOJ_1deg"):
            h = hi[f"{b}_{c}"]
            nb = h["next_best_excluding_ga"]
            q1, q2 = pairs[f"{b}_GA_PSOBV_{c}"], pairs[f"{b}_GA_{nb}_{c}"]
            log(f"[eval] HR {b} {c}: order {', '.join(h['order_qualified'][:4])} ({h['best_aep']:.2f}, {h['second_aep']:.2f}); "
                f"GA-PSO-VNS {q1['mean_diff_gwh']:+.3f} p={q1['p']:.2g} ({q1['first_better']}/{q1['second_better']}); "
                f"GA-{nb} {q2['mean_diff_gwh']:+.3f} p={q2['p']:.2g}; paper test GA vs PSO-VNS pHolm={ptest[f'{b}_{c}']['PSOBV']['p_holm']:.2g}"
                f" max pHolm={max(t['p_holm'] for t in ptest[f'{b}_{c}'].values()):.2g}")
    out["above_installed_GA"] = {f"{b}_{c}": meth[f"{b}_GA"][f"above_{c}"] for b in BUDS for c in HR_EV}
    out["mean_drop_5to1_gwh"] = {f"{b}_{a}": meth[f"{b}_{a}"]["mean_drop_5to1_gwh"] for b in BUDS for a in HR11}
    summ["hr_eval"] = out
    # ---- LaTeX: mean AEP under three evaluators (eleven methods) + GA tests
    lines = [f"Installed block & -- & " + " & ".join(f"{out['installed'][c]:.2f}" for c in ("5deg_2.5", "1deg_0.5", "PyWakeNOJ_1deg"))
             + " & & -- & " + " & ".join(f"{out['installed'][c]:.2f}" for c in ("5deg_2.5", "1deg_0.5", "PyWakeNOJ_1deg")) + " \\\\",
             "\\midrule"]
    order = sorted(HR11, key=lambda a: -(meth[f"30030_{a}"]["mean_1deg_0.5"] or 0) if meth[f"30030_{a}"]["feasible"] >= 15 else 1e9)
    for a in order:
        c = []
        for b in BUDS:
            e = meth[f"{b}_{a}"]
            c.append(f"{e['feasible']}/{e['runs']}")
            for cc in ("5deg_2.5", "1deg_0.5", "PyWakeNOJ_1deg"):
                v = e[f"mean_{cc}"]
                t = "--" if v is None else f"{v:.2f}"
                if v is not None and hi[f"{b}_{cc}"]["best"] == a:
                    t = f"\\textbf{{{t}}}"
                c.append(t)
            if b == 6030:
                c.append("")
        lines.append(f"{LAB[a]} & " + " & ".join(c) + " \\\\")
    lines.append("\\midrule")
    lines.append("\\multicolumn{10}{l}{\\emph{GA minus method, jointly feasible seeds (number in Feas.): mean $\\Delta$AEP (GWh/yr), Wilcoxon $p$}} \\\\")
    refs = []
    for b in BUDS:
        for c in ("5deg_2.5", "1deg_0.5", "PyWakeNOJ_1deg"):
            nb = hi[f"{b}_{c}"]["next_best_excluding_ga"]
            if nb not in refs:
                refs.append(nb)
    if "PSOBV" in refs:
        refs.remove("PSOBV")
    for ref in ["PSOBV"] + refs:
        c, cp_ = [], []
        for b in BUDS:
            q0 = pairs[f"{b}_GA_{ref}_5deg_2.5"]
            c.append(f"{q0.get('n', 0)}"); cp_.append("")
            for cc in ("5deg_2.5", "1deg_0.5", "PyWakeNOJ_1deg"):
                q = pairs[f"{b}_GA_{ref}_{cc}"]
                c.append("--" if "p" not in q else (f"{q['mean_diff_gwh']:+.2f}").replace("-", "$-$"))
                cp_.append("" if "p" not in q else fp(q["p"]))
            if b == 6030:
                c.append(""); cp_.append("")
        lines.append(f"vs.\\ {LAB[ref]}: $\\Delta$AEP & " + " & ".join(c) + " \\\\")
        lines.append(f"\\quad Wilcoxon $p$ & " + " & ".join(cp_) + " \\\\")
    # paper test row: GA vs PSO-VNS, Holm over 10
    c = []
    for b in BUDS:
        c.append("")
        for cc in ("5deg_2.5", "1deg_0.5", "PyWakeNOJ_1deg"):
            t = ptest[f"{b}_{cc}"]["PSOBV"]
            c.append(f"{fp(t['p_holm'])} ({t['outcome']})")
        if b == 6030:
            c.append("")
    lines.append("vs.\\ PSO-VNS: $p_{\\rm Holm}$ & " + " & ".join(c) + " \\\\")
    head = ("& \\multicolumn{4}{c}{6{,}030 evaluations} & & \\multicolumn{4}{c}{30{,}030 evaluations} \\\\\n"
            "\\cmidrule(lr){2-5}\\cmidrule(lr){7-10}\n"
            "Layout / method & Feas. & 5$^\\circ$ & 1$^\\circ$ & PyWake & & Feas. & 5$^\\circ$ & 1$^\\circ$ & PyWake")
    repga = rep.get("rev2_gahr", {})
    note = ("Every stored final layout re-evaluated (not re-optimized), as in Table~\\ref{tab:F-hr}: 5$^\\circ$: the paper's model "
            "replayed on the stored coordinates (rounded to 1~cm); 1$^\\circ$: bins centred at 0.5$^\\circ$, 1.5$^\\circ$, \\ldots "
            "(each 30$^\\circ$ sector keeps its Weibull $A$, $k$, frequency spread over its 30 bins); PyWake: NOJ of PyWake~2.6.20 "
            "($k=0.04$, defaults) at the same 1$^\\circ$ bins. Pool: the ten methods of Table~\\ref{M-tab:hr-site} and GA, random "
            "initialization, 30 seeds. Mean AEP (GWh/yr) of the feasible runs (bold: highest among methods with at least 15 "
            "feasible runs). Rows ordered by the 1$^\\circ$ mean at 30{,}030 evaluations. Lower block: GA minus PSO-VNS and minus "
            "the next best method under some evaluator, on the jointly feasible seeds (count in the Feas.\\ column; positive = GA "
            "better), with the unadjusted Wilcoxon signed-rank $p$; last row: the paper's run-level test of GA against PSO-VNS "
            "(infeasible runs below every feasible run), Holm-adjusted over the 10 comparisons of GA under the same evaluator "
            f"(W/T/L from GA's side). GA replay of the stored 5$^\\circ$ AEP: maximum absolute difference "
            f"{repga.get('max_abs_dev_gwh', float('nan')):.3f}~GWh/yr, {repga.get('n_abs_dev_gt_0_01', 0)} of "
            f"{repga.get('n_feasible', 0)} feasible layouts differ by more than 0.01~GWh/yr.")
    tex.append(MR.table("table*", "Horns Rev~1 16-turbine block: the GA layouts and those of the ten methods of the main study "
                        "re-evaluated with 1$^\\circ$ direction bins and with PyWake NOJ (random starts, 30 seeds, 6{,}030 and "
                        "30{,}030 evaluations).", "tab:S-r3-hr-eval", "lcccccccccc".replace("lcccccccccc", "lccccc" + "cccc"),
                        head, lines, sep="2pt", pos="!htb", note=note))


# ====================================================================== C3: Lillgrund 1-deg
SP9 = list(R2.SP9)


def block_lg(summ, tex, recompute):
    E = reeval_layouts(recompute, os.path.join(HERE, "rev3_sites_reeval.csv"))
    L = E[E.Site == "LG"].copy()
    inst = L[L.Algorithm == "INSTALLED"].iloc[0]
    X = add_minspacing(L[L.Algorithm != "INSTALLED"].copy(), "LG")
    X["Seed"] = X.Seed.astype(int); X["Budget"] = X.Budget.astype(int)
    cols = ["5deg_2.5", "1deg_0.5", "PyWakeNOJ_1deg"]
    out = dict(installed={c: float(inst[c]) for c in cols}, n_layouts=int(len(X)),
               installed_drop_5to1_pct=float(100 * (1 - inst["1deg_0.5"] / inst["5deg_2.5"])), budgets={})
    d = X[X.Feasible]["5deg_2.5"] - X[X.Feasible].Objective
    out["replay_5deg"] = dict(n_feasible=int(X.Feasible.sum()), max_abs_dev_gwh=float(d.abs().max()),
                              n_abs_dev_gt_0_01=int((d.abs() > 0.01).sum()), mean_signed_dev_gwh=float(d.mean()))
    for b in BUDS:
        s = X[X.Budget == b]
        bo = dict(methods={}, highest_mean={}, psovns_tests={})
        for a in SP9:
            f = s[(s.Algorithm == a) & s.Feasible]
            e = dict(runs=int((s.Algorithm == a).sum()), feasible=int(len(f)))
            for c in cols:
                e[f"mean_{c}"] = float(f[c].mean()) if len(f) else None
                e[f"above_{c}"] = int((f[c] > out["installed"][c]).sum())
            e["mean_drop_5to1_pct"] = float(100 * (1 - f["1deg_0.5"].mean() / f["5deg_2.5"].mean())) if len(f) else None
            bo["methods"][a] = e
        for c in cols:
            q = {a: e[f"mean_{c}"] for a, e in bo["methods"].items() if e["feasible"] >= 15 and e[f"mean_{c}"] is not None}
            oq = sorted(q, key=lambda a: -q[a])
            bo["highest_mean"][c] = dict(order_qualified=oq, values={a: q[a] for a in oq})
            bo["psovns_tests"][c] = paper_tests(s, "PSOBV", SP9, c)
        out["budgets"][b] = bo
        log(f"[lg] {b}: " + "; ".join(f"{c}: " + ", ".join(f"{a} {bo['highest_mean'][c]['values'][a]:.2f}" for a in bo['highest_mean'][c]['order_qualified'][:3])
                                     for c in cols))
        log(f"      PSO-VNS paper test (Holm 8) 1deg: " + ", ".join(f"{k} {t['p_holm']:.2g}{t['outcome']}" for k, t in bo["psovns_tests"]["1deg_0.5"].items()))
    out["above_installed_total"] = {c: int(sum(out["budgets"][b]["methods"][a][f"above_{c}"] for b in BUDS for a in SP9)) for c in cols}
    F = X[X.Feasible]
    out["mean_drop_5to1_pct_all_feasible"] = float(100 * ((F["5deg_2.5"] - F["1deg_0.5"]) / F["5deg_2.5"]).mean())
    summ["lg_eval"] = out
    lines = [f"Installed block & -- & " + " & ".join(f"{out['installed'][c]:.2f}" for c in cols) + " & & -- & "
             + " & ".join(f"{out['installed'][c]:.2f}" for c in cols) + " & -- \\\\", "\\midrule"]
    B30 = out["budgets"][30030]
    order = sorted(SP9, key=lambda a: -(B30["methods"][a]["mean_1deg_0.5"] or 0) if B30["methods"][a]["feasible"] >= 15 else 1e9 - B30["methods"][a]["feasible"])
    for a in order:
        c = []
        for b in BUDS:
            bo = out["budgets"][b]
            e = bo["methods"][a]
            c.append(f"{e['feasible']}/{e['runs']}")
            for cc in cols:
                v = e[f"mean_{cc}"]
                t = "--" if v is None else f"{v:.2f}"
                if v is not None and bo["highest_mean"][cc]["order_qualified"] and bo["highest_mean"][cc]["order_qualified"][0] == a:
                    t = f"\\textbf{{{t}}}"
                c.append(t)
            if b == 6030:
                c.append("")
            else:
                t = bo["psovns_tests"]["1deg_0.5"].get(a)
                c.append("--" if t is None else f"{fp(t['p_holm'])} ({t['outcome']})")
        lines.append(f"{LAB[a]} & " + " & ".join(c) + " \\\\")
    head = ("& \\multicolumn{4}{c}{6{,}030 evaluations} & & \\multicolumn{5}{c}{30{,}030 evaluations} \\\\\n"
            "\\cmidrule(lr){2-5}\\cmidrule(lr){7-11}\n"
            "Layout / method & Feas. & 5$^\\circ$ & 1$^\\circ$ & PyWake & & Feas. & 5$^\\circ$ & 1$^\\circ$ & PyWake & PSO-VNS vs., 1$^\\circ$")
    note = ("Stored final layouts of Table~\\ref{tab:S-rev-site} re-evaluated (not re-optimized): 5$^\\circ$: the model used in "
            "optimization replayed on the stored coordinates (rounded to 1~cm; maximum absolute difference to the stored AEP "
            f"{out['replay_5deg']['max_abs_dev_gwh']:.3f}~GWh/yr, {out['replay_5deg']['n_abs_dev_gt_0_01']} of "
            f"{out['replay_5deg']['n_feasible']} feasible layouts differ by more than 0.01); 1$^\\circ$: bins centred at 0.5$^\\circ$, "
            "1.5$^\\circ$, \\ldots{} with the sector's Weibull parameters; PyWake: NOJ of PyWake~2.6.20 ($k=0.04$, defaults) at the "
            "same 1$^\\circ$ bins. Mean AEP (GWh/yr) of the feasible runs (stored labels; bold: highest among methods with at least "
            "15 feasible runs). Last column: the paper's run-level test of PSO-VNS against the method with the 1$^\\circ$ AEP "
            "(infeasible runs below every feasible run), Holm over its 8 comparisons, W/T/L from PSO-VNS's side.")
    tex.append(MR.table("table*", "Lillgrund 16-turbine block (minimum spacing $3D$): final layouts re-evaluated with 1$^\\circ$ "
                        "direction bins and PyWake NOJ (nine methods, random starts, 30 seeds, 6{,}030 and 30{,}030 evaluations).",
                        "tab:S-r3-lg-eval", "lccccccccccc", head, lines, sep="2pt", pos="!htb", note=note))


# ====================================================================== B6: elapsed time
def read_orig(exp, cols=("Algorithm", "Dataset", "Radius", "Turbines", "Seed", "Budget", "Init", "Calls", "Seconds", "Feasible", "Objective", "Ideal", "WakeLoss")):
    fs = sorted(glob.glob(os.path.join(HERE, f"mpce_{exp}_s*of*.csv")))
    if not fs:
        return None
    d = pd.concat([pd.read_csv(f, usecols=lambda c: c in cols) for f in fs], ignore_index=True)
    d = MR.std_cols(d)
    return d.drop_duplicates(MR.KEY, keep="last").reset_index(drop=True).assign(Source=f"mpce_{exp}")


def shard_walltimes():
    """'<exp> <file> <n> runs in <s> s' lines of run_*.log."""
    rows = []
    for fn in sorted(glob.glob(os.path.join(HERE, "run_*.log"))):
        seen = set()
        for line in open(fn, errors="ignore"):
            m = re.match(r"(\w+) (\S+\.csv) (\d+) runs in (\d+) s", line.strip())
            if m and m.group(2) not in seen:
                seen.add(m.group(2))
                rows.append(dict(log=os.path.basename(fn), exp=m.group(1), file=m.group(2), runs=int(m.group(3)), wall_s=int(m.group(4))))
    return pd.DataFrame(rows)


def block_time(summ, tex):
    out = {}
    # ---------------- benchmark main comparison + controls (68 cases, 6,030, random): MX.load
    A = MX.load(HERE)
    A = A.assign(ms=1000 * A.Seconds / A.Calls)
    main = list(MR.MAIN8) + ["LXBV", "RSVNS", "RSDVNS"]
    src = {"PSOBV": "mpce_psobv", "PSOC": "mpce_psoc", "SLSQP": "mpce_slsqp", "RSVNS": "mpce_rsvns", "RSDVNS": "mpce_rsdisc",
           "SSA": "fresh_grid", "LXSSA": "fresh_grid", "DE": "fresh_grid", "BVNS": "fresh_vgrid", "SSABV": "fresh_bgrid",
           "LXBV": "fresh_bgrid", "PSOBV25": "mpce_psosplit", "PSOBV75": "mpce_psosplit", "PSOBV90": "mpce_omega90"}
    bm = {}
    for a in main + ["PSOBV25", "PSOBV75", "PSOBV90"]:
        x = A[A.Algorithm == a]
        if not len(x):
            continue
        byN = x.groupby("Turbines").Seconds.median()
        msN = x.groupby("Turbines").ms.median()
        bm[a] = dict(source=src.get(a), runs=int(len(x)), median_s_per_run=float(x.Seconds.median()),
                     median_s_by_N={int(k): float(v) for k, v in byN.items()},
                     median_ms_per_eval=float(x.ms.median()), mean_ms_per_eval=float(x.ms.mean()),
                     pooled_ms_per_eval=float(1000 * x.Seconds.sum() / x.Calls.sum()),
                     median_ms_by_N={int(k): float(v) for k, v in msN.items()},
                     total_cpu_h=float(x.Seconds.sum() / 3600))
    pv = bm["PSOBV"]
    for a, v in bm.items():
        v["ratio_median_s_to_psovns"] = v["median_s_per_run"] / pv["median_s_per_run"]
        v["ratio_by_N_to_psovns"] = {n: v["median_s_by_N"][n] / pv["median_s_by_N"][n] for n in v["median_s_by_N"]}
    # the per-evaluation statement of the manuscript (Section 9.7 / S-cost) -- which statistic gives which number
    M8 = list(MR.MAIN8)
    meta = [a for a in M8 if a != "SLSQP"]
    ver = dict(
        median_of_runs=dict(meta_min=min(bm[a]["median_ms_per_eval"] for a in meta), meta_max=max(bm[a]["median_ms_per_eval"] for a in meta),
                            slsqp=bm["SLSQP"]["median_ms_per_eval"], pso=bm["PSOC"]["median_ms_per_eval"],
                            slsqp_over_pso=bm["SLSQP"]["median_ms_per_eval"] / bm["PSOC"]["median_ms_per_eval"],
                            source="mpce_results.py summary['cost_per_eval'] -> mpce_numbers.tex \\NMsEvalMin/Max/SLSQP, \\NSLSQPSlowdown"),
        mean_of_runs=dict(meta_min=min(bm[a]["mean_ms_per_eval"] for a in meta), meta_max=max(bm[a]["mean_ms_per_eval"] for a in meta),
                          slsqp=bm["SLSQP"]["mean_ms_per_eval"], pso=bm["PSOC"]["mean_ms_per_eval"],
                          slsqp_over_pso=bm["SLSQP"]["mean_ms_per_eval"] / bm["PSOC"]["mean_ms_per_eval"],
                          source="mpce_results.py summary['cost'] (Table S-cost, 'ms per call') -> \\NMsPerCallMin/Max"),
        pooled=dict(slsqp=bm["SLSQP"]["pooled_ms_per_eval"], pso=bm["PSOC"]["pooled_ms_per_eval"],
                    slsqp_over_pso=bm["SLSQP"]["pooled_ms_per_eval"] / bm["PSOC"]["pooled_ms_per_eval"]),
        slsqp_over_pso_median_by_N={n: bm["SLSQP"]["median_ms_by_N"][n] / bm["PSOC"]["median_ms_by_N"][n] for n in bm["PSOC"]["median_ms_by_N"]},
        slsqp_over_psovns_median_s_by_N=bm["SLSQP"]["ratio_by_N_to_psovns"])
    fo = MR.std_cols(pd.read_csv(os.path.join(HERE, "fresh_grid.csv"), usecols=lambda c: c not in ("Coordinates", "Curve")))
    fo = fo[fo.Algorithm.isin(["SLSQP", "PSO"])].assign(ms=lambda d: 1000 * d.Seconds / d.Calls)
    earlier = {}
    for a, g in fo.groupby("Algorithm"):
        earlier[a] = dict(source="fresh_grid.csv (earlier study / platform)", runs=int(len(g)),
                          median_ms_per_eval=float(g.ms.median()), mean_ms_per_eval=float(g.ms.mean()),
                          median_s_by_N={int(k): float(v) for k, v in g.groupby("Turbines").Seconds.median().items()})
    ver["slsqp_rerun_small_vs_large_N"] = dict(
        median_ms_N_le_7=float(A[(A.Algorithm == "SLSQP") & (A.Turbines <= 7)].ms.median()),
        median_ms_N_ge_8=float(A[(A.Algorithm == "SLSQP") & (A.Turbines >= 8)].ms.median()),
        pso_median_ms_N_le_7=float(A[(A.Algorithm == "PSOC") & (A.Turbines <= 7)].ms.median()),
        pso_median_ms_N_ge_8=float(A[(A.Algorithm == "PSOC") & (A.Turbines >= 8)].ms.median()))
    ver["earlier_platform"] = dict(slsqp_median_ms=earlier["SLSQP"]["median_ms_per_eval"], slsqp_mean_ms=earlier["SLSQP"]["mean_ms_per_eval"],
                                   pso_old_median_ms=earlier["PSO"]["median_ms_per_eval"],
                                   slsqp_over_psoc_median=earlier["SLSQP"]["median_ms_per_eval"] / bm["PSOC"]["median_ms_per_eval"])
    out["benchmark"] = dict(methods=bm, per_eval_verification=ver, slsqp_earlier_platform=earlier["SLSQP"], pso_old_setting=earlier["PSO"],
                            data="MX.load: 68 cases, 6,030 evaluations, random init.; Seconds = wall-clock time of the run "
                                 "(perf_counter around the optimizer, incl. its overhead, excl. writing), single process per run")
    log(f"[time] benchmark ms/eval median: " + ", ".join(f"{a} {v['median_ms_per_eval']:.3f}" for a, v in bm.items()))
    log(f"[time] benchmark ms/eval mean:   " + ", ".join(f"{a} {v['mean_ms_per_eval']:.3f}" for a, v in bm.items()))
    log(f"[time] SLSQP/PSO per eval: median {ver['median_of_runs']['slsqp_over_pso']:.2f}, mean {ver['mean_of_runs']['slsqp_over_pso']:.2f}, "
        f"pooled {ver['pooled']['slsqp_over_pso']:.2f}; by N (median): " + ", ".join(f"{n}:{v:.1f}" for n, v in ver["slsqp_over_pso_median_by_N"].items()))
    # ---------------- budget study and feasible init (six largest cases)
    LARGE = [(ds, r, n) for ds in ("1", "2") for r, n in ((500, 10), (750, 12), (1000, 15))]
    isL = lambda d: pd.Series([(a, b, c) in LARGE for a, b, c in zip(d.Dataset, d.Radius, d.Turbines)], index=d.index)
    bud = pd.concat([d for d in (read_orig(e) for e in ("b30k", "b30kp", "b120k", "b120kp", "feas", "feasp")) if d is not None],
                    ignore_index=True)
    bud = bud[bud.Dataset.isin(["1", "2"])]
    r6 = A[isL(A)].assign(Source="main 6,030")
    BL = pd.concat([r6, bud], ignore_index=True)
    BL = BL[isL(BL)]
    m9 = ["PSOBV", "SSABV", "LXBV", "RSVNS", "BVNS", "PSOC", "SSA", "LXSSA", "DE", "SLSQP"]
    bs = {}
    for (b, ini), g in BL.groupby(["Budget", "Init"]):
        k = f"{b}_{ini}"
        med = g.groupby("Algorithm").Seconds.median()
        bs[k] = dict(median_s={a: float(med[a]) for a in m9 if a in med},
                     ratio_to_psovns={a: float(med[a] / med["PSOBV"]) for a in m9 if a in med and "PSOBV" in med},
                     median_ms_per_eval={a: float((1000 * g[g.Algorithm == a].Seconds / g[g.Algorithm == a].Calls).median()) for a in m9 if a in med},
                     runs=int(len(g)), total_cpu_h=float(g.Seconds.sum() / 3600))
    # feasible minus random median time on the same six cases at 6,030 (rough initialization overhead in the runs)
    fr = {}
    for a in m9:
        x1 = BL[(BL.Budget == 6030) & (BL.Init == "feasible") & (BL.Algorithm == a)]
        x0 = BL[(BL.Budget == 6030) & (BL.Init == "random") & (BL.Algorithm == a)]
        if len(x1) and len(x0):
            m = x1.merge(x0, on=["Dataset", "Radius", "Turbines", "Seed"], suffixes=("_f", "_r"))
            fr[a] = dict(median_diff_s=float((m.Seconds_f - m.Seconds_r).median()), n=int(len(m)))
    out["budget_feasinit"] = dict(settings=bs, feasible_minus_random_s=fr, cases="six largest benchmark cases, 30 seeds")
    log("[time] budget/feas median s: " + "; ".join(f"{k}: " + ", ".join(f"{a} {v:.1f}" for a, v in s_["median_s"].items()) for k, s_ in bs.items()))
    # ---------------- Horns Rev 1 (hrfix + GA)
    HH = load_hr()
    hs = {}
    for (b, ini), g in HH.groupby(["Budget", "Init"]):
        med = g.groupby("Algorithm").Seconds.median()
        hs[f"{b}_{ini}"] = dict(median_s={a: float(med[a]) for a in HR11 if a in med},
                                ratio_to_psovns={a: float(med[a] / med["PSOBV"]) for a in HR11 if a in med and "PSOBV" in med},
                                median_ms_per_eval={a: float((1000 * g[g.Algorithm == a].Seconds / g[g.Algorithm == a].Calls).median()) for a in HR11 if a in med},
                                runs=int(len(g)))
    out["hr"] = dict(settings=hs, note="mpce_hrfix (one batch, 24 shards) and rev2_gahr (GA, separate revision batch)")
    log("[time] HR median s: " + "; ".join(f"{k}: " + ", ".join(f"{a} {v:.1f}" for a, v in s_["median_s"].items()) for k, s_ in hs.items()))
    # ---------------- IEA37
    IE = load_iea()
    ie = {}
    for (n, b), g in IE.groupby(["Turbines", "Budget"]):
        med = g.groupby("Algorithm").Seconds.median()
        ie[f"{n}T_{b}"] = dict(median_s={a: float(med[a]) for a in IEA13 if a in med},
                               ratio_to_psovns={a: float(med[a] / med["PSOBV"]) for a in IEA13 if a in med},
                               ratio_to_ga={a: float(med[a] / med["GA"]) for a in IEA13 if a in med and "GA" in med},
                               mean_s={a: float(g[g.Algorithm == a].Seconds.mean()) for a in IEA13 if a in med},
                               median_ms_per_eval={a: float((1000 * g[g.Algorithm == a].Seconds / g[g.Algorithm == a].Calls).median()) for a in IEA13 if a in med})
    out["iea"] = dict(settings=ie, note="mpce_iea16/36 (nine methods), mpce_iea16p/36p (PSO-VNS arm, separate batch), rev2_gaiea, "
                                        "rev2_grad (revision batches); ratios across batches indicative")
    log("[time] IEA median s: " + "; ".join(f"{k}: " + ", ".join(f"{a} {v:.1f}" for a, v in s_["median_s"].items()) for k, s_ in ie.items()))
    # ---------------- revision studies of Table S-time (verification)
    rv = {}
    sp, _ = R2.read_rev2("spacing", HERE)
    if sp is not None:
        med = sp.groupby("Algorithm").Seconds.median()
        rv["spacing"] = dict(median_s=med.to_dict(), ratio_to_psovns=(med / med["PSOBV"]).to_dict(), runs=int(len(sp)),
                             cases=sorted(set(map(str, zip(sp.Dataset, sp.Radius, sp.Turbines)))).__len__())
    for e in ("lg16", "lg16b"):
        d, _ = R2.read_rev2(e, HERE)
        med = d.groupby("Algorithm").Seconds.median()
        rv[e] = dict(median_s=med.to_dict(), ratio_to_psovns=(med / med["PSOBV"]).to_dict())
    out["revision_studies"] = rv
    # ---------------- total compute: the 36,190 original records (composition of 05_setup.tex) + revision records
    def secs(df):
        return float(df.Seconds.sum())
    comp = {}
    comp["main_and_controls_68cases (11 methods x 2,040)"] = (int(A[A.Algorithm.isin(main)].shape[0]), secs(A[A.Algorithm.isin(main)]))
    old = MR.std_cols(pd.read_csv(os.path.join(HERE, "fresh_grid.csv"), usecols=lambda c: c not in ("Coordinates", "Curve")))
    old = old[old.Algorithm == "PSO"]
    comp["old PSO setting (fresh_grid PSO, 2,040)"] = (int(len(old)), secs(old))
    spl = A[A.Algorithm.isin(["PSOBV25", "PSOBV75", "PSOBV90"])]
    comp["budget split (psosplit + omega90)"] = (int(len(spl)), secs(spl))
    cs = read_orig("csweep")
    comp["coefficient sweep (csweep)"] = (int(len(cs)), secs(cs))
    bb = pd.concat([read_orig(e) for e in ("b30k", "b30kp", "b120k", "b120kp")], ignore_index=True)
    bb = bb[bb.Dataset.isin(["1", "2"])]
    comp["budget 30,030 / 120,030 (six largest cases)"] = (int(len(bb)), secs(bb))
    ff = pd.concat([read_orig(e) for e in ("feas", "feasp")], ignore_index=True)
    ff = ff[ff.Dataset.isin(["1", "2"])]
    comp["feasible initialization (six largest cases)"] = (int(len(ff)), secs(ff))
    hf = read_orig("hrfix")
    comp["Horns Rev 1 (hrfix)"] = (int(len(hf)), secs(hf))
    ia = pd.concat([read_orig(e) for e in ("iea16", "iea36", "iea16p", "iea36p")], ignore_index=True)
    comp["IEA37 (10 methods)"] = (int(len(ia)), secs(ia))
    tot_n = sum(v[0] for v in comp.values()); tot_s = sum(v[1] for v in comp.values())
    rev = {}
    for e in ("ga", "gahr", "gaiea", "grad", "laplace", "spacing", "lg16", "lg16b"):
        d, _ = R2.read_rev2(e, HERE)
        rev[e] = (int(len(d)), secs(d))
    rev_n = sum(v[0] for v in rev.values()); rev_s = sum(v[1] for v in rev.values())
    # superseded but run (not used in the paper): HR rows of b30k/b120k/feas/psobv, fresh HR, hr16new
    sw = shard_walltimes()
    wt = {}
    if len(sw):
        for exp, g in sw.groupby("exp"):
            files = list(g.file)
            ss = 0.0
            for f in files:
                p = os.path.join(HERE, f)
                if os.path.exists(p):
                    ss += float(pd.read_csv(p, usecols=["Seconds"]).Seconds.sum())
            wt[exp] = dict(shards_with_log=int(len(g)), runs=int(g.runs.sum()), sum_shard_wall_h=float(g.wall_s.sum() / 3600),
                           max_shard_wall_h=float(g.wall_s.max() / 3600), sum_run_seconds_h=float(ss / 3600),
                           ratio_cpu_to_shard_wall=float(ss / g.wall_s.sum()) if g.wall_s.sum() else None)
    out["compute"] = dict(original_components={k: dict(records=v[0], cpu_h=v[1] / 3600) for k, v in comp.items()},
                          original_records=tot_n, original_cpu_h=tot_s / 3600,
                          revision_components={k: dict(records=v[0], cpu_h=v[1] / 3600) for k, v in rev.items()},
                          revision_records=rev_n, revision_cpu_h=rev_s / 3600,
                          shard_logs=wt,
                          note="Seconds of a run = serial wall-clock time of that run in one process; the sum over the "
                               "records is the total serial-equivalent compute (CPU-seconds). Runs were executed in "
                               "parallel worker processes (multiprocessing Pool, several shards on a 4-core machine), so the "
                               "elapsed time of a study is much shorter than its total; ratio_cpu_to_shard_wall = sum of the "
                               "run Seconds of a logged shard / its logged wall time (= average number of concurrently busy "
                               "workers).")
    log(f"[time] original records {tot_n} = {tot_s / 3600:.1f} CPU-h; revision records {rev_n} = {rev_s / 3600:.1f} CPU-h")
    log("[time] components: " + "; ".join(f"{k} {v[0]} {v[1] / 3600:.1f} h" for k, v in comp.items()))
    if wt:
        log("[time] shard logs: " + "; ".join(f"{k}: {v['shards_with_log']} shards, sum wall {v['sum_shard_wall_h']:.2f} h, "
                                            f"run-sum {v['sum_run_seconds_h']:.2f} h, x{v['ratio_cpu_to_shard_wall']:.2f}" for k, v in wt.items()))
    summ["time"] = out


def block_init(summ, recompute):
    """Wall time of the feasibility-preserving initialization (30 layouts per run), measured now."""
    cache = os.path.join(HERE, "rev3_sites_inittime.csv")
    import feasible_init as FI
    import hornsrev_model as hrm
    from wflop_model import R
    cases = [(ds, r, n) for ds in ("1", "2") for r, n in ((500, 10), (750, 12), (1000, 15))] + [("HR", 0, 16)]
    if os.path.exists(cache) and not recompute:
        T = pd.read_csv(cache, dtype={"Dataset": str})
    else:
        rows = []
        orig_uniform = FI._uniform_in_site
        cnt = [0]

        def counting(n, circle_r, poly):
            cnt[0] += 1
            return orig_uniform(n, circle_r, poly)
        FI._uniform_in_site = counting
        # reference: stored final layouts of the feasible-start PSO runs (PSO never leaves the best initial layout)
        st = pd.concat([pd.read_csv(f, usecols=["Algorithm", "Dataset", "Radius", "Turbines", "Seed", "Init", "Coordinates"])
                        for f in sorted(glob.glob(os.path.join(HERE, "mpce_feasx_s*of8.csv")))
                        + sorted(glob.glob(os.path.join(HERE, "mpce_hrfix_s*of24.csv")))], ignore_index=True)
        st = st[(st.Algorithm == "PSOC") & (st.Init == "feasible")]
        st["Dataset"] = st.Dataset.astype(str).str.replace(r"\.0$", "", regex=True)
        st["Radius"] = pd.to_numeric(st.Radius, errors="coerce").fillna(0).astype(int)
        t_all = time.time()
        for ds, r, n in cases:
            if ds == "HR":
                _, poly = hrm.site(16)
                gen = FI.make_generator(4 * hrm.D, poly=poly)
            else:
                gen = FI.make_generator(8 * R, circle_r=r)
            for seed in range(1, 31):
                np.random.seed(seed)
                cnt[0] = 0
                t, tc = time.perf_counter(), time.process_time()
                P = gen(30, 2 * n)
                sec, cpu = time.perf_counter() - t, time.process_time() - tc
                ref = st[(st.Dataset == ds) & (st.Radius == r) & (st.Turbines == n) & (st.Seed == seed)]
                dmin = np.nan
                if len(ref):
                    q = parse_xy(ref.Coordinates.iloc[0]).ravel()
                    dmin = float(np.abs(P - q[None]).max(1).min())
                rows.append(dict(Dataset=ds, Radius=r, Turbines=n, Seed=seed, Layouts=30, Seconds=sec, CPUSeconds=cpu, Tries=cnt[0], LoadAvg1=os.getloadavg()[0],
                                 PSOFeasibleStartMatchMaxAbsM=dmin))
            log(f"  [init] {ds}-{r}-{n}: median {np.median([x['Seconds'] for x in rows[-30:]]):.2f} s "
                f"(CPU {np.median([x['CPUSeconds'] for x in rows[-30:]]):.2f} s, load {np.mean([x['LoadAvg1'] for x in rows[-30:]]):.1f}), "
                f"max {max(x['Seconds'] for x in rows[-30:]):.2f} s, tries {sum(x['Tries'] for x in rows[-30:])}/900, "
                f"PSO match max {np.nanmax([x['PSOFeasibleStartMatchMaxAbsM'] for x in rows[-30:]]):.4f} m ({time.time() - t_all:.0f} s)")
        FI._uniform_in_site = orig_uniform
        T = pd.DataFrame(rows)
        T.insert(0, "Machine", f"{os.cpu_count()} cores; measured {time.strftime('%Y-%m-%d %H:%M')}, one process")
        T.to_csv(cache, index=False, float_format="%.6g")
    out = dict(design="feasible_init.make_generator, 30 layouts per run (N_p = 30) after np.random.seed(seed), seeds 1-30; "
                      "spacing 8R = 308 m (benchmark) / 4D = 320 m (Horns Rev 1)", cases={})
    for (ds, r, n), g in T.groupby(["Dataset", "Radius", "Turbines"], sort=False):
        out["cases"][f"{ds}-{r}-{n}"] = dict(median_s=float(g.Seconds.median()), max_s=float(g.Seconds.max()),
                                             median_cpu_s=float(g.CPUSeconds.median()), max_cpu_s=float(g.CPUSeconds.max()),
                                             mean_loadavg1=float(g.LoadAvg1.mean()),
                                             min_s=float(g.Seconds.min()), mean_tries_per_layout=float(g.Tries.sum() / g.Layouts.sum()),
                                             max_tries_run=int(g.Tries.max()),
                                             pso_match_max_abs_m=float(g.PSOFeasibleStartMatchMaxAbsM.max()),
                                             pso_matched=int((g.PSOFeasibleStartMatchMaxAbsM <= 0.0051).sum()))
    bench = T[T.Dataset != "HR"]
    out["benchmark_six_largest"] = dict(median_s=float(bench.Seconds.median()), max_s=float(bench.Seconds.max()),
                                        median_cpu_s=float(bench.CPUSeconds.median()), max_cpu_s=float(bench.CPUSeconds.max()))
    out["note"] = ("Seconds = wall time, CPUSeconds = process CPU time of the initialization of one run (30 layouts); the machine "
                   "(4 cores) was shared with other jobs (1-min load average in LoadAvg1), so wall times are inflated relative "
                   "to CPU times; CPU time is the better estimate of the cost on an idle core")
    out["machine"] = str(T.Machine.iloc[0]) if "Machine" in T else None
    summ["init_time"] = out
    log("[init] " + "; ".join(f"{k}: med {v['median_s']:.2f} max {v['max_s']:.2f} s, tries/layout {v['mean_tries_per_layout']:.2f}, "
                              f"PSO match {v['pso_matched']}/30" for k, v in out["cases"].items()))


PROBE_CODE = r'''
import sys, os, time, json
sys.path.insert(0, HEREDIR)
os.chdir(HEREDIR)
import mpce_experiments as E, pandas as pd
st = pd.read_csv("mpce_slsqp_s0of1.csv", usecols=["Dataset", "Radius", "Turbines", "Seed", "Objective", "Seconds"])
for n in (7, 8, 10, 15):
    t = time.process_time(); r = E.run_grid(("SLSQP", 1, 1000 if n == 15 else 750, n, 1, 6030, "random")); c = time.process_time() - t
    s = st[(st.Dataset == 1) & (st.Radius == (1000 if n == 15 else 750)) & (st.Turbines == n) & (st.Seed == 1)].iloc[0]
    print(json.dumps(dict(Dataset=1, Radius=1000 if n == 15 else 750, Turbines=n, Seed=1, Seconds=r["Seconds"], CPUSeconds=c,
                          Objective=r["Objective"], StoredObjective=float(s.Objective), StoredSeconds=float(s.Seconds),
                          Calls=r["Calls"], BLASThreads=os.environ.get("OPENBLAS_NUM_THREADS"),
                          LoadAvg1=os.getloadavg()[0])), flush=True)
'''


def block_probe(summ, recompute):
    """Re-time single MS-SLSQP runs (data set I, seed 1, N = 7, 8, 10, 15) with single-threaded BLAS, to test whether
    the step in the stored mpce_slsqp times between N = 7 and N = 8 is a property of the method."""
    import subprocess
    cache = os.path.join(HERE, "rev3_sites_slsqp_probe.csv")
    if recompute or not os.path.exists(cache):
        env = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1")
        code = PROBE_CODE.replace("HEREDIR", repr(HERE))
        p = subprocess.run([sys.executable, "-c", code], env=env, capture_output=True, text=True, timeout=1800)
        rows = [json.loads(l) for l in p.stdout.splitlines() if l.startswith("{")]
        pd.DataFrame(rows).to_csv(cache, index=False, float_format="%.6g")
    P = pd.read_csv(cache)
    summ["slsqp_probe"] = dict(rows=P.to_dict("records"),
                               note="MS-SLSQP (mpce_experiments.run_grid, unchanged) rerun now with OPENBLAS/OMP/MKL threads = 1, "
                                    "one process, shared machine; StoredSeconds/StoredObjective from mpce_slsqp. SLSQP results "
                                    "are platform dependent, so the objective need not match the stored one.")
    log("[probe] " + "; ".join(f"N={int(r.Turbines)}: {r.Seconds:.2f} s (CPU {r.CPUSeconds:.2f}) vs stored {r.StoredSeconds:.2f} s; "
                              f"obj {r.Objective:.1f} vs stored {r.StoredObjective:.1f}" for r in P.itertuples()))


def time_table(summ, tex):
    """tab:S-r3-time (benchmark, budget, feasible init., initialization timing) and tab:S-r3-time-sites."""
    T = summ["time"]
    bm = T["benchmark"]["methods"]
    old = T["benchmark"].get("slsqp_earlier_platform", {})
    bf = T["budget_feasinit"]["settings"]
    order = ["PSOBV", "PSOC", "SSABV", "SSA", "LXSSA", "DE", "BVNS", "SLSQP", "LXBV", "RSVNS", "RSDVNS"]

    def rc(med, ratio):
        if med is None or not np.isfinite(med):
            return "--"
        m = f"{med:.1f}" if med >= 10 else f"{med:.2f}"
        return m if ratio is None or not np.isfinite(ratio) else f"{m} ({ratio:.2f})"
    lines = []
    for a in order:
        v = bm[a]
        c = [f"{v['median_ms_per_eval']:.2f}", f"{v['mean_ms_per_eval']:.2f}",
             rc(v["median_s_by_N"].get(5), None), rc(v["median_s_by_N"].get(10), None),
             rc(v["median_s_by_N"].get(15), v["ratio_by_N_to_psovns"].get(15))]
        for k in ("6030_feasible", "30030_random", "120030_random"):
            s_ = bf.get(k, {})
            c.append(rc(s_.get("median_s", {}).get(a), s_.get("ratio_to_psovns", {}).get(a)) if a in s_.get("median_s", {}) else "n/r")
        lines.append(f"{LAB[a]} & " + " & ".join(c) + " \\\\")
    if old:
        lines.append(f"MS-SLSQP, earlier$^{{a}}$ & {old['median_ms_per_eval']:.2f} & {old['mean_ms_per_eval']:.2f} & "
                     f"{rc(old['median_s_by_N'].get(5), None)} & {rc(old['median_s_by_N'].get(10), None)} & "
                     f"{rc(old['median_s_by_N'].get(15), None)} & & & \\\\")
    it = summ.get("init_time")
    if it:
        cs = it["cases"]; b6 = it["benchmark_six_largest"]; hr = cs.get("HR-0-16", {})
        lines.append("\\midrule")
        lines.append("\\multicolumn{9}{l}{\\emph{Initialization alone (30 layouts per run, L-BFGS-B packing, not charged), "
                     "re-run now, seeds 1--30: wall [CPU] time}} \\\\")
        lines.append("\\multicolumn{3}{l}{Six largest cases} & \\multicolumn{6}{l}{median %.2f~s [%.2f~s], maximum %.2f~s [%.2f~s]} \\\\"
                     % (b6["median_s"], b6["median_cpu_s"], b6["max_s"], b6["max_cpu_s"]))
        lines.append("\\multicolumn{3}{l}{Horns Rev~1 block} & \\multicolumn{6}{l}{median %.2f~s [%.2f~s], maximum %.2f~s [%.2f~s]} \\\\"
                     % (hr.get("median_s", np.nan), hr.get("median_cpu_s", np.nan), hr.get("max_s", np.nan), hr.get("max_cpu_s", np.nan)))
    cp = T["compute"]
    head = ("& \\multicolumn{5}{c}{68 cases, 6{,}030 evaluations, random init.} & \\multicolumn{3}{c}{Six largest cases, s per run} \\\\\n"
            "\\cmidrule(lr){2-6}\\cmidrule(lr){7-9}\n"
            "& \\multicolumn{2}{c}{ms per evaluation} & \\multicolumn{3}{c}{s per run} & 6{,}030 & 30{,}030 & 120{,}030 \\\\\n"
            "\\cmidrule(lr){2-3}\\cmidrule(lr){4-6}\n"
            "Method & median & mean & $N=5$ & $N=10$ & $N=15$ & feas.\\ init. & random & random")
    pr = summ.get("slsqp_probe", {}).get("rows", [])
    prs = ", ".join(f"$N={int(r['Turbines'])}$: {r['Seconds']:.1f}~s (stored {r['StoredSeconds']:.1f}~s)" for r in pr)
    note = ("Elapsed (wall-clock) seconds of a run as stored in the original records (Seconds: time of the optimizer call in one "
            "worker process, incl. optimizer overhead and, for MS-SLSQP, the SLSQP subproblems and constraint evaluations, which are "
            "not charged). ms per evaluation: Seconds/Evaluations$\\times10^3$ of each run, median and mean over the 2{,}040 runs "
            "of the method (the mean is the statistic of Table~\\ref{tab:cost}); s per run: median over the runs with the given $N$ "
            "(60 runs); six largest cases: median over 180 runs (6{,}030 feasible init.\\ includes the packing). In parentheses: "
            "ratio to PSO-VNS. Methods come from different run batches (SSA, LX-SSA, DE, VNS, SSA-VNS, LX-SSA-VNS: files of the "
            "earlier study), so ratios are indicative. MS-SLSQP (rerun batch mpce\\_slsqp): the time per run jumps from 3.1~s at "
            "$N=7$ to 16--22~s at $N\\ge8$, independent of $r$; re-timed now with single-threaded BLAS: " + prs + ". "
            "$^{a}$Same method, runs of the earlier study (fresh\\_grid). n/r: not run. Total serial compute (sum of Seconds): "
            "%.1f CPU-hours for the %s original records and %.1f for the %s revision records; the runs were executed in four "
            "parallel worker processes, so the elapsed time of a study was about a quarter of its total."
            % (cp["original_cpu_h"], f"{cp['original_records']:,}".replace(",", "{,}"), cp["revision_cpu_h"],
               f"{cp['revision_records']:,}".replace(",", "{,}")))
    tex.append(MR.table("table*", "Elapsed time per evaluation and per run in the original benchmark records (median seconds; ratio "
                        "to PSO-VNS), and wall time of the feasibility-preserving initialization.", "tab:S-r3-time",
                        "lcccccccc", head, lines, sep="1.7pt", pos="!htb", note=note))
    # ---- sites
    hs = T["hr"]["settings"]; ie = T["iea"]["settings"]
    lines = []
    for a in ["PSOBV", "PSOC", "SSABV", "SSA", "LXSSA", "DE", "BVNS", "SLSQP", "LXBV", "RSVNS", "GA", "SLSQPX", "PSOSLSQPX"]:
        c = []
        for k in ("6030_random", "6030_feasible", "30030_random", "120030_random"):
            s_ = hs.get(k, {})
            c.append(rc(s_["median_s"].get(a), s_["ratio_to_psovns"].get(a)) if a in s_.get("median_s", {}) else "n/r")
        for k in ("16T_6030", "16T_30030", "36T_6030", "36T_30030"):
            s_ = ie.get(k, {})
            c.append(rc(s_["median_s"].get(a), s_["ratio_to_psovns"].get(a)) if a in s_.get("median_s", {}) else "n/r")
        lines.append(f"{LAB[a]} & " + " & ".join(c) + " \\\\")
    head = ("& \\multicolumn{4}{c}{Horns Rev~1 (16 turbines)} & \\multicolumn{4}{c}{IEA37 Case Study~1} \\\\\n"
            "\\cmidrule(lr){2-5}\\cmidrule(lr){6-9}\n"
            "& 6{,}030 R & 6{,}030 F & 30{,}030 & 120{,}030 & 16, 6{,}030 & 16, 30{,}030 & 36, 6{,}030 & 36, 30{,}030")
    note = ("Median elapsed seconds per run in the stored records (30 seeds; 10 at 120{,}030), in parentheses the ratio to PSO-VNS "
            "in the same setting. R/F: random / feasibility-preserving initialization (F includes the packing, about 3.5~s of CPU "
            "time per run at Horns Rev~1). Batches: Horns Rev~1 mpce\\_hrfix (ten methods, one batch) and rev2\\_gahr (GA); IEA37 "
            "mpce\\_iea16/36 (nine methods), mpce\\_iea16p/36p (PSO-VNS), rev2\\_gaiea (GA) and rev2\\_grad (exact-gradient "
            "methods). Ratios across batches are indicative. n/r: not run.")
    tex.append(MR.table("table*", "Elapsed time per run on Horns Rev~1 and IEA37 Case Study~1 in the original and revision "
                        "records (median seconds; ratio to PSO-VNS).", "tab:S-r3-time-sites", "lcccccccc", head, lines,
                        sep="2pt", pos="!htb", note=note))


# ====================================================================== main
def clean(o):
    if isinstance(o, dict):
        return {str(k): clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [clean(v) for v in o]
    if isinstance(o, (np.floating, float)):
        return None if not np.isfinite(o) else float(o)
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, np.bool_):
        return bool(o)
    return o


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", default="hr,iea,eval,lg,time,init,probe")
    ap.add_argument("--recompute-probe", action="store_true")
    ap.add_argument("--recompute-eval", action="store_true")
    ap.add_argument("--recompute-init", action="store_true")
    ap.add_argument("--out-dir", default=HERE)
    args = ap.parse_args(argv)
    t0 = time.time()
    blocks = [b.strip() for b in args.only.split(",") if b.strip()]
    jf = os.path.join(args.out_dir, "rev3_sites.json")
    summ = json.load(open(jf)) if os.path.exists(jf) and set(blocks) != {"hr", "iea", "eval", "lg", "time", "init", "probe"} else {}
    summ.update(generated=time.strftime("%Y-%m-%d %H:%M:%S"), script="rev3_sites.py",
                protocol=dict(ranking=MR.RANK_RULE, run_level_score="mpce_results.goodness: AEP if feasible, else "
                              "-1e12 - max(0, 1e4 - MinSpacing) (every infeasible run below every feasible run)",
                              paired_test="mpce_results.paired_vs: two-sided Wilcoxon signed-rank on seed-paired score "
                                          "differences (|d| <= 1e-9 dropped), Holm within the stated family, alpha 0.05",
                              friedman="scipy friedmanchisquare over seeds (blocks) x methods of one pool, run-level score",
                              feasibility="stored labels (10^-6 m tolerance at run time)", seeds="1-30 (stored runs)",
                              rng="no random resampling in this script; the initialization timing uses np.random.seed(seed), seeds 1-30"))
    tex = ["%% generated by rev3_sites.py -- do not edit by hand\n"]
    if "hr" in blocks:
        log("\n==== hr"); block_hr(summ, tex)
    if "iea" in blocks:
        log("\n==== iea"); block_iea(summ, tex)
    if "eval" in blocks:
        log("\n==== eval"); block_eval(summ, tex, args.recompute_eval)
    if "lg" in blocks:
        log("\n==== lg"); block_lg(summ, tex, args.recompute_eval)
    if "time" in blocks:
        log("\n==== time"); block_time(summ, tex)
    if "init" in blocks:
        log("\n==== init"); block_init(summ, args.recompute_init)
    if "probe" in blocks:
        log("\n==== probe"); block_probe(summ, args.recompute_probe)
    if "time" in summ:
        time_table(summ, tex)
    summ["runtime_s"] = time.time() - t0
    json.dump(clean(summ), open(jf, "w"), indent=1)
    if set(blocks) >= {"hr", "iea", "eval", "lg", "time"}:
        open(os.path.join(args.out_dir, "rev3_sites_tables.tex"), "w").write("\n".join(tex))
    open(os.path.join(args.out_dir, "rev3_sites.log"), "w").write("\n".join(LOG) + "\n")
    log(f"done in {time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()

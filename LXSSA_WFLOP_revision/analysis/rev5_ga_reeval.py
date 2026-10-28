"""rev5: GA added to the re-evaluation studies of the stored final layouts (power curve / wake model and direction
resolution). No optimization is run; only the stored final layouts of GA are re-evaluated, with exactly the procedure
already applied to the other methods, and the generated tables / numbers are recomputed with GA in the pool.

Usage (from analysis/):  python3 rev5_ga_reeval.py [--procs 2] [--only robust|direction]
Then:                    python3 mpce_check_dir.py --summary rev5_ga_summary_dir.json --macros rev5_ga_numbers_dir.tex

GA records: rev2_ga_s0of10.csv ... rev2_ga_s9of10.csv (real-coded GA, 68 cases, seeds 1-30, 6,030 evaluations, random
initialization; 2,040 runs), read with mpce_results.std_cols as rev2_analysis.read_rev2 (duplicates on KEY: last kept).
Coordinates: the stored run-record coordinates (3 decimals = 1 mm), exactly as for the other methods, whose re-evaluations
(mpce_robustness.reevaluate, mpce_direction.bench_reevaluate) use the Coordinates column of their run records
(fresh_*.csv / mpce_*.csv), not the full-precision rerun captures rev3_fullprec_*.csv. For GA, full-precision captures
exist for 3 records only (rev3_fullprec_rev2_ga.csv); they are re-evaluated here as a side check (JSON only).

(a) Power curve and wake model (mpce_results.py section 8 = mpce_robustness.py; Table tab:robust-final, Table S-tab:robust,
    \\N... robustness macros of mpce_numbers.py, wake-model shift of mpce_inference_extra.model_shift_block = \\NXModelShift...
    macros and Table S-tab:X-modelshift).
    Pool: PSO-VNS, PSO, SSA-VNS, SSA, LX-SSA, DE, GA, VNS, MS-SLSQP (display order of the nine-method tables).
    Every feasible final layout is re-evaluated with mpce_robustness._reval (benchmark linear curve + Jensen, cubic curve,
    cubic curve with 25 m/s cut-out, Gaussian wake k* = 0.04 with the linear curve). The eight methods' values come from
    mpce_reevaluation_cache.csv (read only, never written); GA's new layouts are evaluated with the same function (Pool of
    --procs processes, chunksize 25 as in reevaluate) and cached in rev5_ga_reevaluation_cache.csv. Statistics:
    mpce_robustness.robust_table with mpce_results.rank_rule (unchanged), runs per case and method counted as in
    mpce_results.robustness. The eight-method pool is recomputed first from the same cache and must reproduce
    mpce_summary.json["robustness"] exactly (check recorded in the JSON).
    Outputs: rev5_ga_reevaluation.csv (as mpce_reevaluation.csv, nine methods), rev5_ga_tab_robust_final.tex (generator
    format of mpce_tab_robust_final.tex), rev5_ga_tab_robust_supp.tex (S-tab:robust as in mpce_supplementary.tex, and
    S-tab:X-modelshift as in mpce_supp_inference.tex), rev5_ga_numbers_robust.tex (the robustness \\N... and \\NXModelShift...
    macros with GA in the pool), rev5_ga_summary_robust.json.
(b) Direction resolution (mpce_direction.py, benchmark block): every feasible final layout of the eleven methods + GA
    re-evaluated with 1, 3 and 15 sub-bins per 15-deg bin (Jensen) and 1, 15 sub-bins (Gaussian), mpce_direction.bench_
    reevaluate / bench_block unchanged; main pool MAIN8 + GA (statistics of tab:F-bench / tab:F-equiv, \\NF... macros),
    component-analysis pool ABL9 unchanged (GA is not a component variant). Horns Rev, PyWake and IEA37 blocks are recomputed
    unchanged (GA's Horns Rev runs are handled in rev3_sites.py). The eight-method pool is evaluated first on the same
    re-evaluations and must reproduce mpce_summary_dir.json["benchmark"] (check recorded in the JSON).
    Outputs: rev5_ga_summary_dir.json, rev5_ga_numbers_dir.tex, rev5_ga_supp_direction.tex, rev5_ga_direction_layouts.csv.gz.
No existing file is modified.
"""
import os, sys, glob, re, json, time, argparse, tempfile, shutil, math
import numpy as np, pandas as pd
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import mpce_results as M            # noqa: E402  (functions / constants only)
import mpce_robustness as MRB       # noqa: E402

POOL9 = ["PSOBV", "PSOC", "SSABV", "SSA", "LXSSA", "DE", "GA", "BVNS", "SLSQP"]
FOCUS = "PSOBV"
LAB9 = dict(M.LAB, GA="GA")
LOG = []


def log(s=""):
    print(s, flush=True)
    LOG.append(str(s))


def jclean(o):
    if isinstance(o, dict):
        return {str(k): jclean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [jclean(v) for v in o]
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating, float)):
        return None if not np.isfinite(o) else float(o)
    if isinstance(o, np.bool_):
        return bool(o)
    return o


def load_ga(with_coords=True):
    files = {}
    for fn in glob.glob(os.path.join(HERE, "rev2_ga_s*of*.csv")):
        m = re.search(r"rev2_ga_s(\d+)of(\d+)\.csv$", os.path.basename(fn))
        if m:
            files.setdefault(int(m.group(2)), {})[int(m.group(1))] = fn
    k = 10
    assert k in files and set(files[k]) == set(range(k)), files.keys()
    df = pd.concat([pd.read_csv(files[k][i]) for i in range(k)], ignore_index=True)
    df = M.std_cols(df).drop_duplicates(M.KEY, keep="last").reset_index(drop=True)
    df["Source"] = "rev2_ga"
    assert len(df) == 2040 and set(df.Algorithm) == {"GA"} and (df.Budget == 6030).all() and (df.Init == "random").all()
    return df


def benchmark_runs():
    M.set_focus(FOCUS)
    ALL, avail, fb, new = M.load(HERE, False)
    return ALL, avail, fb


def maxdiff(a, b, path="", out=None):
    """Largest absolute difference between two nested JSON-like structures (numbers only; keys of a)."""
    if out is None:
        out = dict(max=0.0, where="", mismatch=[])
    if isinstance(a, dict):
        for k, v in a.items():
            if k in ("per_case",):
                pass
            if not isinstance(b, dict) or str(k) not in {str(x) for x in b}:
                out["mismatch"].append(f"{path}/{k}: missing")
                continue
            bb = b[k] if k in b else b[str(k)]
            maxdiff(v, bb, f"{path}/{k}", out)
    elif isinstance(a, (list, tuple)):
        if not isinstance(b, (list, tuple)) or len(a) != len(b):
            out["mismatch"].append(f"{path}: list")
        else:
            for i, (x, y) in enumerate(zip(a, b)):
                maxdiff(x, y, f"{path}[{i}]", out)
    elif isinstance(a, (bool, np.bool_)) or isinstance(a, str) or a is None:
        if a != b:
            out["mismatch"].append(f"{path}: {a!r} != {b!r}")
    else:
        if b is None or not np.isfinite(float(a)):
            if not (b is None and (a is None or not np.isfinite(float(a)))):
                out["mismatch"].append(f"{path}: {a!r} != {b!r}")
        else:
            d = abs(float(a) - float(b))
            if d > out["max"]:
                out["max"], out["where"] = d, path
    return out


# ====================================================================== (a) power curve and wake model
def robust_final_tex(rb, focus=FOCUS, lab=LAB9):
    """tab:robust-final exactly as mpce_results.main_text_tables (section 'robustness, single column')."""
    rl = []
    for col, name in (("Objective", "Benchmark"), ("Cubic", "Cubic curve"), ("CubicCutout", "Cubic, 25 m/s cut-out"), ("Gauss", "Gaussian wake")):
        v = rb.get(col)
        if not v:
            continue
        ar = v["avg_rank"]; best = min(ar, key=ar.get)
        v["best_method"] = best
        v["focus_rank_position"] = int(1 + sorted(ar.values()).index(ar[focus]))
        d1, d2 = v["rel_change_pct"]["1"], v["rel_change_pct"]["2"]
        rl.append(f"{name} & {d1:+.1f} & {d2:+.1f} & {ar[focus]:.2f} & {lab[best]} & " +
                  ("-- & --" if col == "Objective" else f"{v['tau']:.2f} & {v['same_best_pct']:.0f}") + " \\\\")
    return ("\\begin{table}[!t]\n\\centering\n\\caption{Robustness to the Benchmark Model: Final Layouts Re-Evaluated, Not Re-Optimized}\n"
            "\\label{tab:robust-final}\n\\scriptsize\\setlength{\\tabcolsep}{2.4pt}\n\\begin{tabular}{lcccccc}\n\\toprule\n"
            "& \\multicolumn{2}{c}{$\\Delta$ (\\%)} & & & & \\\\\n\\cmidrule(lr){2-3}\nModel & DS I & DS II & Rank & Best & $\\bar\\tau$ & Same \\\\\n\\midrule\n" +
            "\n".join(rl) + "\n\\bottomrule\n\\end{tabular}\n" +
            M.tnote("All feasible final layouts of the %d cases. $\\Delta$: mean relative change of the objective (\\%%); Rank: "
                    "average rank of %s; Best: best-ranked method; $\\bar\\tau$: mean Kendall correlation between the benchmark "
                    "and the alternative ordering within a case; Same: cases (\\%%) with unchanged best method. Ranks of all "
                    "methods: supplementary material." % (68, lab[focus])) + "\\end{table}\n")


def supp_robust_tex(tex):
    """S-tab:robust as written into mpce_supplementary.tex by mpce_results.main (placement, tabcolsep, xref, evals)."""
    t = tex.replace("[!t]", "[!htb]", 1).replace("{\\tabcolsep}{3pt}", "{\\tabcolsep}{1.8pt}")
    t = re.sub(r"\\ref\{(tab:friedman68)\}", r"\\ref{M-\1}", t)
    return M.evals(t)


def robust_macros(rb, focus=FOCUS, lab=LAB9):
    """The robustness macros of mpce_numbers.py (section 'robustness'), same formulas."""
    import mpce_numbers as MN
    out = []
    P = lambda n, v: out.append(f"\\newcommand{{\\{n}}}{{{v}}}")
    P("NCubicOverI", MN.num(-rb["Cubic"]["rel_change_pct"]["1"], 1))
    P("NCubicOverII", MN.num(-rb["Cubic"]["rel_change_pct"]["2"], 1))
    for col, nm in (("Cubic", "Cubic"), ("CubicCutout", "Cutout"), ("Gauss", "Gauss")):
        P(f"NTau{nm}", MN.num(rb[col]["tau"], 2))
        P(f"NSame{nm}", MN.num(rb[col]["same_best_pct"], 0))
        P(f"NRank{nm}PSOVNS", MN.num(rb[col]["avg_rank"][focus], 2))
        P(f"NRank{nm}PSO", MN.num(rb[col]["avg_rank"]["PSOC"], 2))
        P(f"NPos{nm}PSOVNS", MN.ORD[1 + sorted(rb[col]["avg_rank"].values()).index(rb[col]["avg_rank"][focus])])
        P(f"NBest{nm}", lab[min(rb[col]["avg_rank"], key=rb[col]["avg_rank"].get)])
    P("NGaussDeltaI", MN.num(rb["Gauss"]["rel_change_pct"]["1"], 1, sign=True))
    return out


def model_shift(RVcsv_df, S, pool):
    """mpce_inference_extra.model_shift_block on a re-evaluation file of the given pool (function unchanged; its module
    constant MAIN8 is set to the pool for the call and restored)."""
    import mpce_inference_extra as MX
    tmp = tempfile.mkdtemp(prefix="rev5_ms_")
    try:
        RVcsv_df.to_csv(os.path.join(tmp, "mpce_reevaluation.csv"), index=False)
        old, oldlab = MX.MAIN8, dict(MX.LAB)
        MX.MAIN8 = list(pool); MX.LAB["GA"] = "GA"
        try:
            ms = MX.model_shift_block(S, tmp)
        finally:
            MX.MAIN8 = old; MX.LAB.clear(); MX.LAB.update(oldlab)
    finally:
        shutil.rmtree(tmp)
    ms["source"] = "rev5_ga_reevaluation.csv" if "GA" in pool else "mpce_reevaluation.csv"
    return ms


def model_shift_macros(MS):
    """\\NXModelShift... macros, formulas of mpce_inference_extra.write_macros."""
    import mpce_inference_extra as MX
    num = MX.num
    out = []
    m = lambda n, v, c="": out.append(f"\\newcommand{{\\{n}}}{{{v}}}" + (f"   % {c}" if c else ""))
    mp = MS["pairs"]["PSOBV-PSOC"]
    m("NXModelShiftMean", num(MS["mean_abs_shift_pp"], 2), "mean |L_Gauss - L_Jensen| of the same feasible final layouts (pp)")
    m("NXModelShiftMedian", num(MS["median_abs_shift_pp"], 2), "median |L_Gauss - L_Jensen| (pp)")
    m("NXModelShiftPTen", num(MS["q10_abs_shift_pp"], 2)); m("NXModelShiftPNinety", num(MS["q90_abs_shift_pp"], 2))
    m("NXModelShiftSigned", num(MS["mean_signed_shift_pp"], 2, sign=True), "mean L_Gauss - L_Jensen (pp)")
    m("NXModelShiftDsI", num(MS["by_dataset"]["1"]["mean_abs_shift_pp"], 2)); m("NXModelShiftDsII", num(MS["by_dataset"]["2"]["mean_abs_shift_pp"], 2))
    m("NXModelShiftNLayouts", f"{MS['n_layouts']:,}".replace(",", "{,}"), "feasible final layouts re-evaluated (9 methods incl. GA)")
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
    return out


def model_shift_tex(MS, pool):
    """S-tab:X-modelshift as mpce_inference_extra.write_tables (wake-model uncertainty), pool size in the caption."""
    import mpce_inference_extra as MX
    LABX = dict(MX.LAB, GA="GA")
    plab = lambda a, b: f"{LABX[a]} vs.\\ {LABX[b]}"
    lines = [f"All feasible final layouts & {MS['n_layouts']:,} & {MS['mean_abs_shift_pp']:.3f} & {MS['median_abs_shift_pp']:.3f} & "
             f"[{MS['q10_abs_shift_pp']:.3f}, {MS['q90_abs_shift_pp']:.3f}] & ${MS['mean_signed_shift_pp']:+.3f}$ \\\\".replace(",", "{,}", 1)]
    for ds, txt in (("1", "Data Set I"), ("2", "Data Set II")):
        v = MS["by_dataset"][ds]
        lines.append(f"\\quad {txt} & {v['n']:,} & {v['mean_abs_shift_pp']:.3f} & {v['median_abs_shift_pp']:.3f} & -- & -- \\\\".replace(",", "{,}", 1))
    lines.append("\\midrule")
    lines.append("\\multicolumn{6}{l}{\\emph{Case-mean difference PSO-VNS minus other method (pp)}} \\\\")
    lines.append("Pair & Cases & Jensen & Gaussian & mean $|$change$|$ & sign changes \\\\")
    for b in [x for x in pool if x != FOCUS]:
        v = MS["pairs"][f"{FOCUS}-{b}"]
        lines.append(f"{plab(FOCUS, b)} & {v['n_cases']} & ${v['mean_diff_jensen_pp']:+.3f}$ & ${v['mean_diff_gauss_pp']:+.3f}$ & "
                     f"{v['mean_abs_change_pp']:.3f} & {v['sign_changes']} \\\\")
    return MX.table("table", "Wake-model uncertainty from the study's own layouts. Every feasible final layout of the nine methods (68 cases, 6{,}030 evaluations) is re-evaluated (not re-optimized) with the benchmark Jensen model and with the Gaussian wake of Bastankhah and Port\\'e-Agel ($k^*=0.04$, Table~\\ref{tab:robust}); wake loss in pp of the wake-free power. Top: absolute difference $|L_{\\rm Gauss}-L_{\\rm Jensen}|$ of the same layout (mean, median, 10--90\\%% quantiles) and mean signed difference, to be compared with the equivalence margin of $\\pm%g$~pp. Bottom: mean over the cases (both methods $\\ge15$ feasible runs) of the case-mean difference under each model, the mean absolute change of the per-case difference between the models, and the cases in which its sign differs." % MX.EQ_MARGIN,
                    "tab:X-modelshift", "lccccc", "Layouts & $n$ & mean & median & 10--90\\% & signed", lines, sep="3pt")


def part_robust(args):
    t0 = time.time()
    log("=" * 78 + "\n(a) power curve and wake model: GA added to the re-evaluation of mpce_robustness.py")
    ALL, avail, fb = benchmark_runs()
    R6 = ALL[(ALL.Budget == 6030) & (ALL.Init == "random")]
    G8 = R6[R6.Dataset.isin(["1", "2"]) & R6.Algorithm.isin(M.MAIN)].reset_index(drop=True)
    MAINP = [a for a in M.MAIN if a in set(G8.Algorithm)]
    assert MAINP == M.MAIN8, MAINP
    GA = load_ga()
    G9 = pd.concat([G8, GA], ignore_index=True)
    out = dict(generated=time.strftime("%Y-%m-%d %H:%M:%S"), script="rev5_ga_reeval.py", pool=POOL9,
               ga_records=dict(files=[f"rev2_ga_s{i}of10.csv" for i in range(10)], runs=int(len(GA)),
                               feasible=int(GA.Feasible.sum()), coordinates="stored run-record coordinates (3 decimals, 1 mm), "
                               "as for the other methods (mpce_robustness.reevaluate uses the run records' Coordinates column)"))
    # ---- re-evaluation (cache of the eight methods read only; GA cache separate)
    cache = pd.read_csv(os.path.join(HERE, "mpce_reevaluation_cache.csv"))
    gfn = os.path.join(HERE, "rev5_ga_reevaluation_cache.csv")
    gcache = pd.read_csv(gfn) if os.path.exists(gfn) else pd.DataFrame(columns=["Key"] + MRB.VALS)
    def feas(G, methods):
        F = G[G.Feasible & G.Algorithm.isin(methods)].reset_index(drop=True).copy()
        F["Key"] = [MRB._key(str(d), c) for d, c in zip(F.Dataset, F.Coordinates)]
        return F
    F8, F9 = feas(G8, MAINP), feas(G9, POOL9)
    miss8 = int((~F8.Key.isin(set(cache.Key))).sum())
    out["cache_check"] = dict(n_layouts_8=int(len(F8)), n_missing_in_cache_8=miss8, cache_rows=int(len(cache)),
                              cache_duplicate_keys=int(cache.Key.duplicated().sum()))
    assert miss8 == 0, f"{miss8} eight-method layouts not in mpce_reevaluation_cache.csv"
    have = set(cache.Key) | set(gcache.Key)
    todo = F9[~F9.Key.isin(have)].drop_duplicates("Key")
    assert set(todo.Algorithm) <= {"GA"}, set(todo.Algorithm)
    out["cache_check"]["ga_layouts"] = int((F9.Algorithm == "GA").sum())
    out["cache_check"]["ga_keys_already_in_8method_cache"] = int(F9[F9.Algorithm == "GA"].Key.isin(set(cache.Key)).sum())
    if len(todo):
        log(f"  re-evaluating {len(todo)} new GA layouts on {args.procs} processes (mpce_robustness._reval) ...")
        t = time.time()
        with Pool(args.procs) as pool:
            res = pool.map(MRB._reval, list(zip(todo.Dataset.astype(str), todo.Coordinates)), chunksize=25)
        new = pd.DataFrame(res); new.insert(0, "Key", todo.Key.values)
        gcache = pd.concat([gcache, new], ignore_index=True)
        gcache.to_csv(gfn, index=False)
        log(f"  done in {time.time() - t:.0f} s; GA cache {len(gcache)} layouts -> rev5_ga_reevaluation_cache.csv")
    # spot check: the eight-method cache equals a fresh evaluation with the current model (random sample)
    rs = np.random.RandomState(20261028)
    smp = F8.drop_duplicates("Key").iloc[rs.choice(F8.Key.nunique(), size=min(args.spot, F8.Key.nunique()), replace=False)]
    with Pool(args.procs) as pool:
        chk = pd.DataFrame(pool.map(MRB._reval, list(zip(smp.Dataset.astype(str), smp.Coordinates)), chunksize=25))
    cv = cache.drop_duplicates("Key").set_index("Key").loc[smp.Key.values, MRB.VALS].values
    out["cache_check"]["spot_check"] = dict(n=int(len(smp)), max_rel_diff=float(np.max(np.abs(chk[MRB.VALS].values - cv) / np.abs(cv))))
    log(f"  cache spot check ({len(smp)} eight-method layouts re-evaluated now): max rel diff {out['cache_check']['spot_check']['max_rel_diff']:.2e}")
    C = pd.concat([cache, gcache], ignore_index=True)
    C[MRB.VALS] = C[MRB.VALS].astype(float)
    # ---- eight-method reproduction (exactly mpce_results.robustness)
    nr8 = G8[G8.Algorithm.isin(MAINP)].groupby(M.CASE + ["Algorithm"]).size()
    RV8 = F8.merge(cache, on="Key", how="left")
    tex8, out8 = MRB.robust_table(RV8, MAINP, M.LAB, M.rank_rule, n_runs=nr8)
    ref = json.load(open(os.path.join(HERE, "mpce_summary.json")))["robustness"]
    d8 = maxdiff(jclean(out8), ref)
    old_csv = pd.read_csv(os.path.join(HERE, "mpce_reevaluation.csv"))
    k = ["Algorithm", "Dataset", "Radius", "Turbines", "Seed"]
    old_csv["Dataset"] = old_csv.Dataset.astype(str)
    mm = RV8.merge(old_csv[k + MRB.VALS], on=k, suffixes=("", "_old"))
    out["reproduction_8"] = dict(max_abs_diff_vs_mpce_summary=d8["max"], where=d8["where"], mismatches=d8["mismatch"],
                                 rows_vs_mpce_reevaluation_csv=int(len(mm)), rows_8=int(len(RV8)),
                                 max_abs_value_diff_vs_mpce_reevaluation_csv=float(np.nanmax(np.abs(
                                     mm[MRB.VALS].values - mm[[v + "_old" for v in MRB.VALS]].values))))
    log(f"  eight-method pool vs mpce_summary.json robustness: max |diff| {d8['max']:.2e} ({d8['where']}), mismatches {d8['mismatch']}")
    assert d8["max"] < 1e-9 and not d8["mismatch"]
    # ---- nine-method pool
    nr9 = G9[G9.Algorithm.isin(POOL9)].groupby(M.CASE + ["Algorithm"]).size()
    RV9 = F9.merge(C.drop_duplicates("Key"), on="Key", how="left")
    assert RV9[MRB.VALS].notna().all().all() and len(RV9) == len(F9)
    tex9, out9 = MRB.robust_table(RV9, POOL9, LAB9, M.rank_rule, n_runs=nr9)
    out["robustness_8"] = jclean(out8)
    out["robustness_9"] = jclean(out9)
    rf9 = robust_final_tex(out["robustness_9"]); rf8 = robust_final_tex(out["robustness_8"], lab=M.LAB)
    # the eight-method single-column table must equal the stored generator output
    stored = open(os.path.join(HERE, "mpce_tab_robust_final.tex")).read().split("\n", 1)[1]
    out["reproduction_8"]["robust_final_tex_identical"] = bool(rf8 == stored)
    log(f"  eight-method tab:robust-final identical to mpce_tab_robust_final.tex: {rf8 == stored}")
    hdr = "%% generated by rev5_ga_reeval.py (pipeline of mpce_results.py section 8, GA added to the pool) -- do not edit by hand\n"
    open(os.path.join(HERE, "rev5_ga_tab_robust_final.tex"), "w").write(hdr + rf9)
    # GA's own results
    ga = RV9[RV9.Algorithm == "GA"]
    oth = RV9[RV9.Algorithm != "GA"]
    gres = dict(n_layouts=int(len(ga)), n_layouts_by_ds={ds: int((ga.Dataset == ds).sum()) for ds in ("1", "2")})
    for col in ("Cubic", "CubicCutout", "Gauss"):
        gres[f"rel_change_pct_{col}"] = {ds: float(100 * (ga[ga.Dataset == ds][col] / ga[ga.Dataset == ds]["Linear"] - 1).mean()) for ds in ("1", "2")}
        gres[f"rel_change_pct_{col}_other8"] = {ds: float(100 * (oth[oth.Dataset == ds][col] / oth[oth.Dataset == ds]["Linear"] - 1).mean()) for ds in ("1", "2")}
    for col, nm in (("Objective", "Benchmark"), ("Cubic", "Cubic"), ("CubicCutout", "CubicCutout"), ("Gauss", "Gauss")):
        ar = out9[col]["avg_rank"]
        gres[f"avg_rank_{nm}"] = ar["GA"]
        gres[f"rank_position_{nm}"] = int(1 + sorted(ar.values()).index(ar["GA"]))
        gres[f"order_{nm}"] = sorted(ar, key=ar.get)
    dev = (ga.Linear - ga.Objective).abs() / ga.Objective
    gres["rounding_check"] = dict(median_rel_dev=float(dev.median()), max_rel_dev=float(dev.max()), n_rel_dev_above_1e_4=int((dev > 1e-4).sum()))
    gres["mean_loss_jensen_pct"] = float((100 * (ga.Ideal - ga.Linear) / ga.Ideal).mean())
    gres["mean_loss_gauss_pct"] = float((100 * (ga.Ideal - ga.Gauss) / ga.Ideal).mean())
    out["ga"] = gres
    # ---- wake-model shift (mpce_inference_extra.model_shift_block), 8 (reproduction) and 9 methods
    S9 = M.case_stats(G9, POOL9)
    cols = [c for c in RV9.columns if c not in ("Coordinates", "Curve", "Key")]
    RV9[cols].to_csv(os.path.join(HERE, "rev5_ga_reevaluation.csv"), index=False)
    ms8 = model_shift(RV8[[c for c in RV8.columns if c not in ("Coordinates", "Curve", "Key")]], M.case_stats(G8, MAINP), MAINP)
    ref_ms = json.load(open(os.path.join(HERE, "mpce_summary_extra.json")))["model_shift"]
    dms = maxdiff({k_: v for k_, v in jclean(ms8).items() if k_ not in ("source", "note")}, ref_ms)
    out["reproduction_8"]["model_shift_max_abs_diff_vs_mpce_summary_extra"] = dms["max"]
    out["reproduction_8"]["model_shift_mismatches"] = dms["mismatch"]
    log(f"  eight-method model shift vs mpce_summary_extra.json: max |diff| {dms['max']:.2e}, mismatches {dms['mismatch']}")
    ms9 = model_shift(RV9[cols], S9, POOL9)
    out["model_shift_8"] = jclean(ms8); out["model_shift_9"] = jclean(ms9)
    # ---- CHECK-FINAL C21 / C22 and CHECK-EXTRA X29 / X30 conditions with the nine-method pool
    import mpce_check_final as CF
    import mpce_check_extra as CX
    cf = {c[0]: c for c in CF.CONDITIONS}
    cx = {c[0]: c for c in CX.CHECKS}
    CX.EQ_MARGIN = json.load(open(os.path.join(HERE, "mpce_summary_extra.json")))["equivalence_levels"]["margin_pp"]
    chk = {}
    for tag, rob, ms in (("8", out["robustness_8"], out["model_shift_8"]), ("9", out["robustness_9"], out["model_shift_9"])):
        for cid in ("C21", "C22"):
            try:
                chk[f"{cid}_{tag}"] = bool(cf[cid][3]({"robustness": rob}))
            except Exception as e:
                chk[f"{cid}_{tag}"] = f"error {e!r}"
        for cid in ("X29", "X30"):
            chk[f"{cid}_{tag}"] = bool(cx[cid][2]({"model_shift": ms}))
    out["checks"] = chk
    log("  checks (8 -> 9 methods): " + ", ".join(f"{c} {chk[c + '_8']} -> {chk[c + '_9']}" for c in ("C21", "C22", "X29", "X30")))
    # ---- supplementary tables and macros
    supp = ["%% generated by rev5_ga_reeval.py -- do not edit by hand",
            "%% Table S-tab:robust (mpce_supplementary.tex) and Table S-tab:X-modelshift (mpce_supp_inference.tex) with GA in the pool",
            "", supp_robust_tex(tex9), model_shift_tex(ms9, POOL9)]
    open(os.path.join(HERE, "rev5_ga_tab_robust_supp.tex"), "w").write("\n".join(supp))
    mac = ["% generated by rev5_ga_reeval.py -- do not edit by hand",
           "% robustness macros of mpce_numbers.py (section 'robustness') and \\NXModelShift... macros of mpce_inference_extra.py,",
           "% recomputed with GA in the pool (nine methods). Same names as the macros they replace."]
    mac += robust_macros(out["robustness_9"]) + model_shift_macros(ms9)
    open(os.path.join(HERE, "rev5_ga_numbers_robust.tex"), "w").write("\n".join(mac) + "\n")
    # ---- full-precision side check (3 GA records with rev3 full-precision captures)
    fp = pd.read_csv(os.path.join(HERE, "rev3_fullprec_rev2_ga.csv"))
    rows = []
    for _, r in fp.iterrows():
        xy_full = ";".join(" ".join(p.split()) for p in r.CoordinatesFull.split(";"))
        vf = MRB._reval((str(r.Dataset), xy_full))
        g = RV9[(RV9.Algorithm == "GA") & (RV9.Dataset == str(r.Dataset)) & (RV9.Radius == int(r.Radius)) &
                (RV9.Turbines == int(r.Turbines)) & (RV9.Seed == int(r.Seed))]
        e = dict(Dataset=str(r.Dataset), Radius=int(r.Radius), Turbines=int(r.Turbines), Seed=int(r.Seed),
                 in_pool=bool(len(g)), full={k_: float(v) for k_, v in vf.items()})
        if len(g):
            e["rounded"] = {k_: float(g.iloc[0][k_]) for k_ in MRB.VALS}
            e["max_rel_diff"] = float(max(abs(vf[k_] - g.iloc[0][k_]) / abs(g.iloc[0][k_]) for k_ in MRB.VALS))
        rows.append(e)
    out["fullprec_side_check"] = rows
    out["runtime_s"] = time.time() - t0
    out["log"] = LOG
    json.dump(jclean(out), open(os.path.join(HERE, "rev5_ga_summary_robust.json"), "w"), indent=1)
    log(f"  nine-method pool: " + "; ".join(f"{c}: PSO-VNS {out9[c]['avg_rank']['PSOBV']:.2f}, best "
                                            f"{LAB9[min(out9[c]['avg_rank'], key=out9[c]['avg_rank'].get)]}, tau {out9[c]['tau']:.2f}, "
                                            f"same {out9[c]['same_best_pct']:.0f}" for c in ("Objective", "Cubic", "CubicCutout", "Gauss")))
    log(f"  wrote rev5_ga_tab_robust_final.tex, rev5_ga_tab_robust_supp.tex, rev5_ga_numbers_robust.tex, rev5_ga_reevaluation.csv, "
        f"rev5_ga_summary_robust.json in {time.time() - t0:.0f} s")


# ====================================================================== (b) direction resolution
def part_direction(args):
    t0 = time.time()
    log("=" * 78 + "\n(b) direction resolution: GA added to mpce_direction.py (benchmark block)")
    import mpce_direction as D
    LOG.clear()
    ALL, avail, fb = benchmark_runs()
    M.set_focus(D.FOCUS)
    D.BM_MARGIN = M.EQ_MARGIN
    stamp = time.strftime("%Y-%m-%d %H:%M:%S")
    ALLM8 = list(D.ALLM)
    G8 = ALL[(ALL.Budget == 6030) & (ALL.Init == "random") & ALL.Dataset.isin(["1", "2"]) & ALL.Algorithm.isin(ALLM8)].reset_index(drop=True)
    GA = load_ga()
    G = pd.concat([G8, GA], ignore_index=True)
    log(f"[1] Benchmark: {len(G)} runs ({G.Feasible.sum()} feasible), {G.groupby(['Dataset', 'Radius', 'Turbines']).ngroups} cases "
        f"(GA: {len(GA)} runs, {GA.Feasible.sum()} feasible)")
    F = D.bench_reevaluate(G, args.procs, log)
    # eight-method reproduction on the same re-evaluations (module pools unchanged)
    BM8, _ = D.bench_block(M, G8, F[F.Algorithm != "GA"].reset_index(drop=True), log)
    ref = json.load(open(os.path.join(HERE, "mpce_summary_dir.json")))["benchmark"]
    d8 = maxdiff(jclean(BM8), ref)
    log(f"  eight-method pool vs mpce_summary_dir.json benchmark: max |diff| {d8['max']:.2e} ({d8['where']}), mismatches {d8['mismatch'][:5]}")
    # nine-method main pool (GA added); component-analysis pool unchanged
    D.MAIN8 = list(POOL9)
    D.ALLM = list(POOL9) + ["LXBV", "RSVNS", "RSDVNS"]
    D.LAB["GA"] = "GA"; D.MAC["GA"] = "GA"
    log("  --- nine-method pool ---")
    BM, Fm = D.bench_block(M, G, F, log)
    log(f"  reproduction (S = 1 vs recorded objective): {BM['reproduction']}")
    log("  mean wake loss 15 -> 5 -> 1 deg: " + ", ".join(
        f"{D.LAB[a]} {BM['res']['rec']['mean_loss_qualified'][a]:.3f}/{BM['res']['J3']['mean_loss_qualified'][a]:.3f}/"
        f"{BM['res']['J15']['mean_loss_qualified'][a]:.3f}" for a in D.MAIN8))
    # GA's own direction-resolution numbers
    fg = Fm[Fm.Algorithm == "GA"]
    BM["ga"] = dict(n_layouts=int(len(fg)), layout_loss_rise_pp=BM["layout_loss_rise_pp"]["GA"],
                    layout_loss_rise_rel_pct=BM["layout_loss_rise_rel_pct"]["GA"],
                    reproduction=dict(median_rel_dev=float(((fg.J1 - fg.Objective).abs() / fg.Objective).median()),
                                      max_rel_dev=float(((fg.J1 - fg.Objective).abs() / fg.Objective).max())))
    log("\n[2] Horns Rev 16 (unchanged; GA's Horns Rev runs: rev3_sites.py)")
    HR, HX = D.hr_block(M, ALL, log)
    pwf = os.path.join(HERE, "pywake_check.csv")
    PW = pd.read_csv(pwf) if os.path.exists(pwf) else None
    import iea37_projected as IPm
    tmp = tempfile.mkdtemp(prefix="rev5_iea_")
    try:
        IP = IPm.run(M, ALL, log, out_dir=tmp)
    finally:
        shutil.rmtree(tmp)
    summary = dict(generated=stamp, script="rev5_ga_reeval.py (mpce_direction.py functions; GA added to the main pool)",
                   eq_margin_pp=M.EQ_MARGIN, benchmark=BM, horns_rev=HR,
                   pywake=None if PW is None else PW.to_dict(orient="records"), iea37_projected=IP,
                   data_availability={k: v.get("status") for k, v in avail.items()}, fallbacks=fb,
                   ga_records=dict(files=[f"rev2_ga_s{i}of10.csv" for i in range(10)], runs=int(len(GA)), feasible=int(GA.Feasible.sum()),
                                   coordinates="stored run-record coordinates (3 decimals, 1 mm), as for the other methods"),
                   reproduction_8=dict(max_abs_diff_vs_mpce_summary_dir=d8["max"], where=d8["where"], mismatches=d8["mismatch"]),
                   benchmark_8=BM8, log=LOG)
    ref_all = json.load(open(os.path.join(HERE, "mpce_summary_dir.json")))
    for blk in ("horns_rev", "iea37_projected"):
        dd = maxdiff(jclean(summary[blk]), ref_all[blk])
        summary["reproduction_8"][f"{blk}_max_abs_diff"] = dd["max"]; summary["reproduction_8"][f"{blk}_mismatches"] = dd["mismatch"]
        log(f"  {blk} vs mpce_summary_dir.json: max |diff| {dd['max']:.2e}, mismatches {dd['mismatch'][:5]}")
    json.dump(D.jclean(summary), open(os.path.join(HERE, "rev5_ga_summary_dir.json"), "w"), indent=1)
    keep = ["Algorithm", "Dataset", "Radius", "Turbines", "Seed", "Objective", "Ideal", "J1", "J3", "J15", "I1", "G1", "G15"]
    Fm[keep].to_csv(os.path.join(HERE, "rev5_ga_direction_layouts.csv.gz"), index=False, float_format="%.8g")
    fn_mac = os.path.join(HERE, "rev5_ga_numbers_dir.tex")
    n = D.write_macros(BM, HR, PW, IP, fn_mac, stamp)
    t = open(fn_mac).read().replace("% generated by mpce_direction.py", "% generated by rev5_ga_reeval.py with the functions of mpce_direction.py, GA added to the main pool,", 1)
    t = t.replace("(68 cases, 11 methods)", "(68 cases, 12 methods incl. GA)").replace("over the 8 methods", "over the 9 methods").replace(
        "Kendall tau of the 8 average ranks", "Kendall tau of the 9 average ranks")
    open(fn_mac, "w").write(t)
    fn_sup = os.path.join(HERE, "rev5_ga_supp_direction.tex")
    D.write_supp(BM, HR, PW, IP, fn_sup, stamp)
    t = open(fn_sup).read().replace("% generated by mpce_direction.py", "% generated by rev5_ga_reeval.py with the functions of mpce_direction.py, GA added to the main pool,", 1)
    open(fn_sup, "w").write(t)
    log(f"\nwrote rev5_ga_summary_dir.json, rev5_ga_numbers_dir.tex ({n} macros), rev5_ga_supp_direction.tex, "
        f"rev5_ga_direction_layouts.csv.gz in {time.time() - t0:.0f} s")
    summary["log"] = LOG
    json.dump(D.jclean(summary), open(os.path.join(HERE, "rev5_ga_summary_dir.json"), "w"), indent=1)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--procs", type=int, default=2)
    ap.add_argument("--spot", type=int, default=300, help="eight-method cached layouts re-evaluated now as a cache check")
    ap.add_argument("--only", choices=("robust", "direction"))
    args = ap.parse_args(argv)
    if args.only in (None, "robust"):
        part_robust(args)
    if args.only in (None, "direction"):
        part_direction(args)


if __name__ == "__main__":
    main()

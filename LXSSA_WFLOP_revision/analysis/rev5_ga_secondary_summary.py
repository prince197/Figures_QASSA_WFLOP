"""GA results of the two secondary benchmark studies (rev5_ga_secondary.py) for the report.

Usage:  python3 rev5_ga_secondary_summary.py      (after rev5_ga_secondary.py, mpce_results.py and rev2_analysis.py)
Output: rev5_ga_secondary_summary.json (and a printed digest)

Reads the run records (rev5_ga_*_s0of1.csv, rev2_ga_s*of10.csv, the stored runs of the other methods) and the
regenerated outputs mpce_summary.json, mpce_feasinit_tests.csv, mpce_budget_case_stats.csv and
../SWEVO_rev2/experiments_docs/rev2_summary.json. New computations (functions of mpce_results.py / rev2_analysis.py,
unchanged):
  * record checks: runs, seeds, Calls == Budget, Seconds;
  * seed pairing at 6,030 feasible: first convergence checkpoint (best feasible of the 30 initial layouts) of GA equal
    to that of every other method, per case and seed (as mpce_results.pairing_check);
  * budget study, run-level PSO-VNS vs GA on the SIX benchmark cases (mpce_results.paired_vs per case with the
    comparison family of the budget section: PSO-VNS vs each of the other ten methods, Holm within the case);
  * runs whose final objective differs from the best initial layout at 6,030 feasible (as summary feasinit);
  * spacing study: all-run paired outcome of PSO-VNS vs GA at 5D / 6D (10 compared split cases x 30 seeds; feasible
    beats infeasible, two feasible runs by objective with ties within 1e-9, two infeasible runs tie; score
    (W + T/2) / n; 95 % percentile case bootstrap, 20,000 resamples of the 10 cases, seed BOOT_SEED -- the seed of
    the hand-typed table S-tab:S-allrun-spacing is not stored, so its intervals are recomputed for every row with
    this seed and both are reported); the PSO-VNS - GA case-mean contrast (rev2_analysis.contrast), reported only
    (not part of the Holm family of the spacing contrast table).
"""
import os, sys, glob, json
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import mpce_results as MR          # noqa: E402
import rev2_analysis as RA         # noqa: E402

BOOT_N, BOOT_SEED = 20000, 20260928
CASE = MR.CASE
EXP_DOCS = os.path.join(HERE, "..", "SWEVO_rev2", "experiments_docs")


def inc(D, cases, neg=False):
    """Boolean array: row's case (Dataset, Radius, Turbines) in cases (not in cases if neg)."""
    cs = set(cases)
    return np.array([(c in cs) != neg for c in zip(D.Dataset, D.Radius, D.Turbines)], bool)


def rd(pat):
    fs = sorted(glob.glob(os.path.join(HERE, pat)))
    return MR.std_cols(pd.concat([pd.read_csv(f) for f in fs], ignore_index=True)), [os.path.basename(f) for f in fs]


def record_checks(D, budget, init, seeds, cases):
    x = D[(D.Budget == budget) & (D.Init == init)]
    per = x.groupby(CASE).Seed.apply(lambda s: sorted(s) == list(seeds))
    return dict(runs=int(len(x)), cases=int(x[CASE].drop_duplicates().shape[0]), expected_cases=len(cases),
                seeds_complete=bool(per.all()), calls_equal_budget=bool((x.Calls == budget).all()),
                feasible=int(x.Feasible.sum()), seconds_median=float(x.Seconds.median()),
                seconds_mean=float(x.Seconds.mean()), seconds_sum=float(x.Seconds.sum()))


def allrun(sp, first, second, spc, rng_seed=BOOT_SEED):
    key = ["Dataset", "Radius", "Turbines", "Seed"]
    A = sp[(sp.Algorithm == first) & (sp.Spacing == spc)].set_index(key)
    B = sp[(sp.Algorithm == second) & (sp.Spacing == spc)].set_index(key)
    j = A.join(B, lsuffix="_a", rsuffix="_b", how="inner")
    fa, fb = j.Feasible_a.astype(bool), j.Feasible_b.astype(bool)
    w = (fa & ~fb) | (fa & fb & (j.Objective_a > j.Objective_b + 1e-9))
    l_ = (~fa & fb) | (fa & fb & (j.Objective_b > j.Objective_a + 1e-9))
    t = ~(w | l_)
    j = j.assign(W=w.astype(int), L=l_.astype(int), T=t.astype(int)).reset_index()
    per = j.groupby(["Dataset", "Radius", "Turbines"])[["W", "T", "L"]].sum()
    n = per.sum(1).values
    sc = (per.W.values + per["T"].values / 2)
    rng = np.random.default_rng(rng_seed)
    idx = rng.integers(0, len(per), size=(BOOT_N, len(per)))
    bs = sc[idx].sum(1) / n[idx].sum(1)
    W, T, L = int(per.W.sum()), int(per["T"].sum()), int(per.L.sum())
    return dict(n_pairs=int(n.sum()), n_cases=int(len(per)), W=W, T=T, L=L, score=(W + T / 2) / n.sum(),
                ci95=[float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))],
                feas_second_pct=float(100 * fb.mean()), feas_first_pct=float(100 * fa.mean()))


def main():
    out = dict(note=__doc__.split("\n")[0])
    ms = json.load(open(os.path.join(HERE, "mpce_summary.json")))
    # ------------------------------------------------------------ study 1: records
    G1, f1 = rd("rev5_ga_b30k_s*of*.csv"); G2, f2 = rd("rev5_ga_b120k_s*of*.csv"); G3, f3 = rd("rev5_ga_feas_s*of*.csv")
    G4, f4 = rd("rev2_ga_s*of10.csv")
    G4 = G4[inc(G4, MR.LARGE)]
    rec = {"6030 random (rev2_ga)": record_checks(G4, 6030, "random", range(1, 31), MR.LARGE),
           "30030 random": record_checks(G1, 30030, "random", range(1, 31), MR.LARGE),
           "120030 random": record_checks(G2, 120030, "random", range(1, 31), MR.LARGE),
           "6030 feasible": record_checks(G3, 6030, "feasible", range(1, 31), MR.LARGE)}
    hrf = glob.glob(os.path.join(HERE, "rev5_ga_hr_s*of*.csv"))
    if hrf:
        GH, fh = rd("rev5_ga_hr_s*of*.csv")
        rec["HR16 120030 random"] = record_checks(GH, 120030, "random", range(1, 11), [MR.HR16])
        rec["HR16 6030 feasible"] = record_checks(GH, 6030, "feasible", range(1, 31), [MR.HR16])
    out["budget_feasinit_records"] = rec
    # seed pairing at 6,030 feasible (first checkpoint = best feasible of the 30 initial layouts)
    FE = pd.concat([pd.read_csv(os.path.join(HERE, "mpce_feas_s0of1.csv")), pd.read_csv(os.path.join(HERE, "mpce_feasp_s0of1.csv"))]
                   + ([pd.read_csv(f) for f in hrf] if hrf else []) + [pd.read_csv(f) for f in sorted(glob.glob(os.path.join(HERE, "rev5_ga_feas_s*of*.csv")))],
                   ignore_index=True)
    FE = MR.std_cols(FE)
    FE = FE[(FE.Init == "feasible").values & inc(FE, MR.LARGE)]
    FE["C0"] = pd.to_numeric(FE.Curve.str.split(";").str[0], errors="coerce")
    g = FE[FE.Algorithm == "GA"].set_index(CASE + ["Seed"]).C0
    mism, comp = 0, 0
    for a in sorted(set(FE.Algorithm) - {"GA"}):
        o = FE[FE.Algorithm == a].set_index(CASE + ["Seed"]).C0
        j = pd.concat([g, o], axis=1, join="inner").dropna()
        comp += len(j); mism += int((j.iloc[:, 0] != j.iloc[:, 1]).sum())
    out["feasible_init_pairing"] = dict(compared_case_seed_method_pairs=comp, mismatches=mism,
                                        note="first stored convergence checkpoint (3 decimals) of GA vs each other method, 6,030 feasible, six largest cases")
    # runs whose final objective differs from the best initial layout (Curve = best feasible objective; computed in
    # mpce_results section 5, summary feasinit)
    out["feasinit_moved_from_best_initial"] = (ms.get("feasinit") or {}).get("runs_moved_from_best_initial")
    # ------------------------------------------------------------ study 1: generated numbers
    fb = ms["feasbudget"]
    out["feasbudget"] = dict(methods=fb["methods"], rank=fb["rank"], mean_loss=fb["mean_loss"], feas_pct=fb["feas_pct"],
                             first_count=fb["first_count"], sole_first_count=fb["sole_first_count"],
                             random_6030=fb["random"], feasible_6030=fb["feasible"],
                             rank_feasible_init=fb.get("rank_feasible_init"), best=fb.get("best"),
                             focus_position=fb.get("focus_position"), close_gap_pp=fb.get("close_gap_pp"), pending=fb["pending"])
    bud = ms["budget"]["per_budget"]
    out["budget_supp"] = {b: dict(avg_rank_7cases=v["avg_rank"], best=v["best"], avg_rank_grid6=v["avg_rank_grid6"],
                                  wtl_vs_PSOBV_7cases=v.get("wtl_vs_focus", {}).get("GA"),
                                  mean_loss_grid_GA=v["mean_loss_grid"].get("GA"), feas_pct_GA=v["feas_pct"].get("GA"),
                                  hr16_aep=v.get("hr16_aep"))
                          for b, v in bud.items()}
    # run-level PSO-VNS vs GA on the six benchmark cases (family of the budget section: all other methods of the case)
    SEC = MR.with_ga(MR.M10)
    ALLB = []
    for pat in ("mpce_b30k_s*of3.csv", "mpce_b30kp_s0of1.csv", "mpce_b120k_s*of8.csv", "mpce_b120kp_s*of4.csv",
                "rev5_ga_b30k_s*of*.csv", "rev5_ga_b120k_s*of*.csv"):
        ALLB.append(rd(pat)[0])
    B = pd.concat(ALLB, ignore_index=True)
    B = B[inc(B, MR.LARGE)]
    w6 = {}
    for b in (30030, 120030):
        x = B[B.Budget == b]
        oc = []
        for c, sub in x.groupby(CASE):
            t = MR.paired_vs(sub, "PSOBV", [a for a in SEC if a != "PSOBV" and a in set(sub.Algorithm)])
            oc += [dict(case="-".join(map(str, c)), **r) for r in t if r["Baseline"] == "GA"]
        w6[b] = dict(W=sum(r["Outcome"] == "W" for r in oc), T=sum(r["Outcome"] == "T" for r in oc),
                     L=sum(r["Outcome"] == "L" for r in oc),
                     per_case={r["case"]: dict(p_holm=r["PHolm"], rb=r["RB"], outcome=r["Outcome"]) for r in oc})
    out["budget_PSOBV_vs_GA_run_level_six_cases"] = w6
    # feasible vs random per case for GA (mpce_feasinit_tests.csv; Holm over the methods of the case)
    ft = pd.read_csv(os.path.join(HERE, "mpce_feasinit_tests.csv"))
    out["feasinit_tests_GA"] = ft[ft.Algorithm == "GA"].drop(columns=["RandomSource"]).to_dict("records")
    bc = pd.read_csv(os.path.join(HERE, "mpce_budget_case_stats.csv"))
    out["budget_case_stats_GA"] = bc[bc.Algorithm == "GA"][["Dataset", "Radius", "Turbines", "Budget", "N", "NFeas", "Loss", "Mean", "Rank"]].to_dict("records")
    # ------------------------------------------------------------ study 2: spacing
    S5, fs = rd("rev5_ga_spacing_s*of*.csv")
    excl = [("1", 500, 10), ("2", 500, 10)]
    out["spacing_records"] = {sp: record_checks(S5[S5.Spacing == sp], 6030, "random", range(1, 31), MR.SPLITCASES)
                              for sp in ("5D", "6D")}
    sp_old, _ = rd("rev2_spacing_s*of40.csv")
    sp_all = pd.concat([sp_old, S5], ignore_index=True)
    sp_all = sp_all[inc(sp_all, excl, neg=True)]
    ar = {}
    for comp in ["PSOC", "SSABV", "SSA", "LXSSA", "DE", "GA", "BVNS", "SLSQP", "RSDVNS", ("SSABV", "RSDVNS")]:
        first, second = comp if isinstance(comp, tuple) else ("PSOBV", comp)
        ar[f"{first}-{second}"] = {spc: allrun(sp_all, first, second, spc) for spc in ("5D", "6D")}
    out["spacing_allrun"] = ar
    out["spacing_seconds_10cases_5D6D_pooled_median"] = {a: float(sp_all[sp_all.Algorithm == a].Seconds.median())
                                                         for a in sorted(set(sp_all.Algorithm))}
    rs = json.load(open(os.path.join(EXP_DOCS, "rev2_summary.json")))["spacing"]
    out["spacing_generated"] = {sp: dict(friedman_excl={k: v for k, v in rs["spacing"][sp]["excl"]["friedman"].items()
                                                        if k in ("chi2", "p", "iman_davenport", "iman_davenport_p", "avg_rank",
                                                                 "p_holm_vs_focus", "sole_best_count", "n_cases")},
                                         avg_rank_all=rs["spacing"][sp]["all"]["friedman"]["avg_rank"],
                                         feasible_pct=rs["spacing"][sp]["excl"]["feasible_pct"],
                                         n_qualified=rs["spacing"][sp]["excl"]["n_qualified"],
                                         mean_loss_pct_qualified=rs["spacing"][sp]["excl"]["mean_loss_pct_qualified"],
                                         contrasts={k: (dict(n=r["n_cases"], mean=r["mean_dloss_pp"], p=r["wilcoxon"]["p"])
                                                        if r else None) for k, r in rs["spacing"][sp]["excl"]["contrasts"].items()})
                                for sp in ("4D", "5D", "6D")}
    out["spacing_holm_over_contrasts_excl"] = rs.get("holm_over_contrasts_excl")
    # PSO-VNS - GA case-mean contrast (informational; not in the Holm family of the contrast table)
    M4 = RA.main_runs()
    M4 = M4[M4.Algorithm.isin(["PSOBV"]).values & inc(M4, MR.SPLITCASES)].assign(Spacing="4D")
    g4, _ = RA.read_rev2("ga", HERE)
    g4 = g4[(g4.Algorithm == "GA").values & inc(g4, MR.SPLITCASES)].assign(Spacing="4D")
    D = pd.concat([M4, g4, sp_old[sp_old.Algorithm == "PSOBV"], S5], ignore_index=True)
    D = D[inc(D, excl, neg=True)]
    con = {}
    for spc in ("4D", "5D", "6D"):
        Gs = D[D.Spacing == spc].reset_index(drop=True)
        S = MR.case_stats(Gs, ["PSOBV", "GA"])
        r = RA.contrast(Gs, S, "PSOBV", "GA")
        con[spc] = dict(n_cases=r["n_cases"], mean_dloss_pp=r["mean_dloss_pp"],
                        case_ci90=r["case_level"]["ci90_mean_dloss_pp"] if r["case_level"] else None,
                        seed_ci90=r["seed_level"]["ci90"] if r.get("seed_level") else None,
                        wilcoxon_p=r["wilcoxon"]["p"], wins_PSOBV=r["wilcoxon"]["wins_first"], wins_GA=r["wilcoxon"]["wins_second"],
                        imputed=r["wilcoxon"]["imputed"], dropped=r["wilcoxon"]["dropped"],
                        hl=r["hodges_lehmann"]["hl"] if r.get("hodges_lehmann") else None)
    out["spacing_contrast_PSOBV_minus_GA"] = con
    json.dump(RA.clean(out), open(os.path.join(HERE, "rev5_ga_secondary_summary.json"), "w"), indent=1)
    # ------------------------------------------------------------ digest
    print(json.dumps(RA.clean({k: out[k] for k in ("budget_feasinit_records", "feasible_init_pairing", "spacing_records")}), indent=1))
    for b in ("6030", "30030", "120030"):
        print("feasbudget", b, "ranks", {a: round(v, 2) for a, v in fb["rank"][b].items()})
    print("written rev5_ga_secondary_summary.json")


if __name__ == "__main__":
    main()

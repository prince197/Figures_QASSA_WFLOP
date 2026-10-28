"""GA in the all-run paired outcome (tab:allrun, S-tab:S-sa-allrun, S-tab:S-sa-dep-bench, S-tab:archive-paired) and
coordinate-precision audit of the GA records.

Usage (from any directory):  python3 analysis/rev5_ga_allrun.py [--out-dir DIR]
Outputs (analysis/ by default): rev5_ga_allrun.json, rev5_ga_allrun.log (stdout, if redirected). The report
rev5_ga_allrun_report.md is written by hand from the JSON. Nothing existing is modified; no optimizer is run.

Procedure: the functions of the generating scripts are imported, not re-implemented or modified:
  benchmark, pooled   rev3_inference.run_matrix / outcome / boot_scores, scipy binomtest (as rev3_inference.b1_block),
                      mpce_inference_extra.holm, case_stats / case_diffs           (tab:S-sa-allrun; tab:allrun pooled)
  benchmark, seeds    rev3_dependence.seed_inference / group_inference / wtl       (tab:S-sa-dep-bench; tab:allrun seeds)
  Lillgrund           audit_archive.reliability (score, bootstrap 20,000, seed 20261004 + budget, exact sign test,
                      Holm over all PSO-VNS comparisons in the table), rev3_dependence.cp (Clopper-Pearson)
                                                                                    (S-tab:archive-paired; tab:allrun)
  precision audit     rev3_precision_audit.audit_file / summarize                   (S-tab:S-sa-precision)
Step 0 (reproduction): every existing row of the three benchmark families and of the Lillgrund table is recomputed
with these functions and asserted to equal the stored JSON / CSV outputs (rev3_inference.json, rev3_dependence.json,
SWEVO_rev2/validation/paired_outcomes.csv) exactly.
Step 1 (GA): GA records rev2_ga_s*of10.csv (benchmark: 68 cases, seeds 1-30, 6,030 evaluations, random starts; loaded
with mpce_inference_extra.std_cols as all other records), rev4_galg_s0of1.csv / rev4_galgb_s0of1.csv (Lillgrund, 6,030 /
30,030). GA is appended as the 13th benchmark pair (index j = 12): pooled bootstrap seed B1_SEED + 10 j (+1 for the
violation variant), seed-level bootstrap seed BOOT_SEED + 12 (both offsets unused by the existing pairs). Holm families
enlarged: benchmark PSO-VNS comparisons 7 -> 8, Lillgrund 16 -> 18 (nine methods at two budgets); the 5 sampling-control
contrasts are unaffected.
Step 2 (precision): the GA benchmark records are part of the existing audit (study rev2_ga in rev3_precision_audit.py;
class-(iii) reruns in rev3_fullprec_rev2_ga*.csv); this is re-checked by re-running audit_file on them. The Lillgrund GA
records (rev4_galg*, written after the audit) are audited here with the same audit_file.
"""
import os, sys, glob, json, math, time, argparse
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np, pandas as pd
from scipy.stats import binomtest

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import mpce_inference_extra as X            # noqa: E402
import rev3_inference as RI                 # noqa: E402
import rev3_dependence as RD                # noqa: E402
import audit_archive as AA                  # noqa: E402
import rev3_precision_audit as PA           # noqa: E402

VALID = os.path.join(ROOT, "SWEVO_rev2", "validation")
GA_J = len(RI.B1_PAIRS)                      # 12: next pair index
MAIN8_GA = RI.MAIN7 + ["GA"]


def clean(o):
    return RD.clean(o)


# ------------------------------------------------------------------ benchmark data
def load_bench():
    A = X.load(HERE)
    G = A[A.Algorithm.isin(set(X.MAIN8) | set(X.ABL))].copy()
    files = sorted(glob.glob(os.path.join(HERE, "rev2_ga_s*of10.csv")))
    assert len(files) == 10
    ga = pd.concat([pd.read_csv(f, usecols=lambda c: c != "Curve") for f in files], ignore_index=True)
    ga = X.std_cols(ga)
    assert not ga.duplicated(X.KEY).any()
    ga = ga[(ga.Budget == 6030) & (ga.Init == "random") & ga.Dataset.isin(["1", "2"])]
    assert len(ga) == 2040 and set(ga.Algorithm) == {"GA"}
    # coordinates -> V (as rev3_inference.main)
    C = RI.load_coords(HERE)
    G = G.merge(C[X.KEY + ["Coordinates", "Objective"]].rename(columns={"Objective": "ObjC"}), on=X.KEY, how="left",
                validate="1:1")
    assert G.Coordinates.notna().all() and np.allclose(G.ObjC, G.Objective, rtol=0, atol=0)
    G = G.drop(columns="ObjC")
    G = pd.concat([G, ga[[c for c in G.columns if c in ga.columns]]], ignore_index=True)
    G["V"] = 0.0
    inf = ~G.Feasible
    G.loc[inf, "V"] = [RI.violation(c, r)[0] for c, r in zip(G.Coordinates[inf], G.Radius[inf])]
    S0 = X.case_stats(G[G.Algorithm != "GA"], sorted(set(X.MAIN8) | set(X.ABL)), X.BASE_THR)
    S = X.case_stats(G, sorted(set(X.MAIN8) | set(X.ABL) | {"GA"}), X.BASE_THR)
    cases = [tuple(c) for c in S0[X.CASE].drop_duplicates().sort_values(X.CASE).itertuples(index=False)]
    assert len(cases) == 68
    # identical case set and ideal value for GA
    I = S.pivot_table(index=X.CASE, columns="Algorithm", values="Ideal")
    ideal_maxrel = float(((I["GA"] - I["PSOBV"]).abs() / I["PSOBV"]).max())
    return G, S0, S, cases, ideal_maxrel


def pooled_row(RM, S, a, b, j, B):
    """rev3_inference.b1_block, one pair (same order of operations and RNG seeds)."""
    Ra, Rb = RM[a], RM[b]
    row = dict(a=a, b=b, n_pairs=int(Ra["feas"].size),
               feas_a_pct=float(100 * Ra["feas"].mean()), feas_b_pct=float(100 * Rb["feas"].mean()))
    d = X.case_diffs(S, a, b)
    row["cond_mean_dloss_pp"] = float(d.mean()); row["cond_n_cases"] = int(len(d))
    row["both_infeasible_pairs"] = int((~Ra["feas"] & ~Rb["feas"]).sum())
    for var, viol in (("main", False), ("violation", True)):
        O = RI.outcome(Ra, Rb, viol)
        W, T, L = int((O == 1).sum()), int((O == 0.5).sum()), int((O == 0).sum())
        r = dict(W=W, T=T, L=L, score=float(O.mean()),
                 p_sign=float(binomtest(W, W + L, 0.5).pvalue) if W + L > 0 else 1.0,
                 case_score_gt_half=int((O.mean(1) > 0.5 + 1e-12).sum()),
                 case_score_lt_half=int((O.mean(1) < 0.5 - 1e-12).sum()))
        r["ci95"] = RI.boot_scores(O, B, RI.B1_SEED + 10 * j + (1 if viol else 0))
        row[var] = r
    return row


def seed_row(RM, cases, a, b, j):
    """rev3_dependence.bench_block, one pair."""
    O = RI.outcome(RM[a], RM[b])
    W, T, L = RD.wtl(O)
    return dict(a=a, b=b, W=W, T=T, L=L, n=int(O.size), seed=RD.seed_inference(O, len(cases), RD.BOOT_SEED + j),
                group=RD.group_inference(O, cases))


def same(x, y, path=""):
    """exact equality of nested JSON-like structures (floats compared exactly)."""
    if isinstance(x, dict):
        for k in x:
            if k in y:
                same(x[k], y[k], path + "/" + k)
        return
    if isinstance(x, (list, tuple)):
        assert len(x) == len(y), path
        for i, (u, v) in enumerate(zip(x, y)):
            same(u, v, f"{path}[{i}]")
        return
    assert x == y, (path, x, y)


def bench_block(B):
    G, S0, S, cases, ideal_maxrel = load_bench()
    algs = sorted({x for p in RI.B1_PAIRS for x in p} | {"GA"})
    RM = {a: RI.run_matrix(G, a, cases) for a in algs}
    REF1 = json.load(open(os.path.join(HERE, "rev3_inference.json")))["B1"]["pairs"]
    REF2 = json.load(open(os.path.join(HERE, "rev3_dependence.json")))["benchmark"]["pairs"]
    pooled, seeds = {}, {}
    t0 = time.time()
    for j, (a, b) in enumerate(RI.B1_PAIRS):                       # step 0: reproduction of all 12 existing rows
        k = RI.pk(a, b)
        pooled[k] = pooled_row(RM, S0, a, b, j, B)
        seeds[k] = seed_row(RM, cases, a, b, j)
        ref1 = {kk: vv for kk, vv in REF1[k].items()}
        for var in ("main", "violation"):
            ref1[var] = {kk: vv for kk, vv in ref1[var].items() if not kk.startswith("p_sign_holm")}
        same(pooled[k], ref1, "B1/" + k)
        same(seeds[k]["seed"], {kk: v for kk, v in REF2[k]["seed"].items() if not kk.endswith("_holm")}, "seed/" + k)
        same(seeds[k]["group"], {kk: v for kk, v in REF2[k]["group"].items() if not kk.endswith("_holm")}, "group/" + k)
        assert (seeds[k]["W"], seeds[k]["T"], seeds[k]["L"]) == (REF2[k]["W"], REF2[k]["T"], REF2[k]["L"])
        print(f"reproduced {k:14s} W/T/L {pooled[k]['main']['W']}/{pooled[k]['main']['T']}/{pooled[k]['main']['L']} "
              f"caseCI {np.round(pooled[k]['main']['ci95']['case'], 3)} seedCI {np.round(seeds[k]['seed']['ci95_joint_seed'], 3)}"
              f" ({time.time() - t0:.0f} s)", flush=True)
    # Holm, old families (reproduction of the adjusted values)
    old = {}
    fam7 = [RI.pk("PSOBV", b) for b in RI.MAIN7]
    ctrl = [RI.pk(a, b) for a, b in RI.CTRL5]
    for fam, keys in (("main7", fam7), ("controls5", ctrl)):
        for var in ("main", "violation"):
            for k, h in zip(keys, X.holm([pooled[k][var]["p_sign"] for k in keys])):
                assert float(h) == REF1[k][var]["p_sign_holm_" + fam], (k, var)
        for key, sub, out in (("p_signflip", "seed", "p_signflip_holm"), ("p_wilcoxon", "seed", "p_wilcoxon_holm"),
                              ("p_signflip_groups", "group", "p_signflip_groups_holm")):
            for k, h in zip(keys, X.holm([seeds[k][sub][key] for k in keys])):
                assert float(h) == REF2[k][sub][out], (k, out)
    for k in fam7:
        old[k] = dict(pooled_sign_holm=REF1[k]["main"]["p_sign_holm_main7"],
                      pooled_sign_holm_violation=REF1[k]["violation"]["p_sign_holm_main7"],
                      signflip_holm=REF2[k]["seed"]["p_signflip_holm"], wilcoxon_holm=REF2[k]["seed"]["p_wilcoxon_holm"],
                      groups_holm=REF2[k]["group"]["p_signflip_groups_holm"])
    # step 1: GA
    k = RI.pk("PSOBV", "GA")
    pooled[k] = pooled_row(RM, S, "PSOBV", "GA", GA_J, B)
    seeds[k] = seed_row(RM, cases, "PSOBV", "GA", GA_J)
    fam8 = [RI.pk("PSOBV", b) for b in MAIN8_GA]
    new = {kk: {} for kk in fam8}
    for var, nm in (("main", "pooled_sign_holm"), ("violation", "pooled_sign_holm_violation")):
        for kk, h in zip(fam8, X.holm([pooled[kk][var]["p_sign"] for kk in fam8])):
            new[kk][nm] = float(h); pooled[kk][var]["p_sign_holm_main8"] = float(h)
    for key, sub, out, nm in (("p_signflip", "seed", "p_signflip_holm8", "signflip_holm"),
                              ("p_wilcoxon", "seed", "p_wilcoxon_holm8", "wilcoxon_holm"),
                              ("p_signflip_groups", "group", "p_signflip_groups_holm8", "groups_holm")):
        for kk, h in zip(fam8, X.holm([seeds[kk][sub][key] for kk in fam8])):
            new[kk][nm] = float(h); seeds[kk][sub][out] = float(h)
    changes = {kk: {nm: dict(old=old[kk][nm], new=new[kk][nm]) for nm in old[kk]} for kk in fam7}
    sig = {kk: dict(pooled=new[kk]["pooled_sign_holm"] < RD.ALPHA, seed=new[kk]["signflip_holm"] < RD.ALPHA) for kk in fam8}
    g, sd = pooled[k], seeds[k]["seed"]
    print(f"GA  {k}: feas {g['feas_a_pct']:.1f}/{g['feas_b_pct']:.1f} dL {g['cond_mean_dloss_pp']:+.3f} ({g['cond_n_cases']} cases) "
          f"W/T/L {g['main']['W']}/{g['main']['T']}/{g['main']['L']} score {g['main']['score']:.4f} "
          f"case {np.round(g['main']['ci95']['case'], 3)} 2st {np.round(g['main']['ci95']['two_stage'], 3)} p {g['main']['p_sign']:.3g} "
          f"Holm8 {new[k]['pooled_sign_holm']:.3g} | viol score {g['violation']['score']:.4f} | seedCI "
          f"{np.round(sd['ci95_joint_seed'], 3)} seeds +{sd['seeds_above']}/-{sd['seeds_below']} p_flip {sd['p_signflip']:.3g} "
          f"Holm8 {new[k]['signflip_holm']:.3g} p_W {sd['p_wilcoxon']:.3g} groups +{seeds[k]['group']['groups_above']}/"
          f"-{seeds[k]['group']['groups_below']} p_grp {seeds[k]['group']['p_signflip_groups']:.3g}", flush=True)
    for kk in fam7:
        print("  Holm 7->8", kk, {nm: f"{v['old']:.3g} -> {v['new']:.3g}" for nm, v in changes[kk].items()}, flush=True)
    return dict(pairs_reproduced=[RI.pk(a, b) for a, b in RI.B1_PAIRS], reproduction="exact (all fields, all 12 pairs, "
                "both variants, Holm main7/controls5 of rev3_inference.json and rev3_dependence.json)",
                ga_pair=dict(pooled=pooled[k], seed=seeds[k]), holm_family_main8=fam8, holm_changes_main7_to_main8=changes,
                adjusted_main8=new, significant_main8=sig, ga_ideal_max_rel_diff_vs_psovns=ideal_maxrel,
                seeds_used=dict(pooled_case_boot=RI.B1_SEED + 10 * GA_J, pooled_violation=RI.B1_SEED + 10 * GA_J + 1,
                                seed_level=RD.BOOT_SEED + GA_J))


# ------------------------------------------------------------------ Lillgrund
def lillgrund_block():
    fr = {}
    for fam, pat in (("lg16", "rev2_lg16_s*of4.csv"), ("lg16b", "rev2_lg16b_s*of30.csv")):
        files = sorted(glob.glob(os.path.join(HERE, pat)))
        fr[fam] = pd.concat([pd.read_csv(f) for f in files], ignore_index=True)
        assert len(fr[fam]) == AA.FAMILIES[fam]
    rows0, pairs0 = AA.reliability(fr)
    ref = pd.read_csv(os.path.join(VALID, "paired_outcomes.csv"), float_precision="round_trip")
    assert len(pairs0) == len(ref) == 16
    for p in pairs0:
        q = ref[(ref.budget == p["budget"]) & (ref.second == p["second"])].iloc[0]
        for c in ("wins", "ties", "losses", "superiority", "bootstrap95_lower", "bootstrap95_upper", "exact_sign_p", "holm16_p"):
            assert p[c] == q[c], (p, c, q[c])
    print("reproduced Lillgrund paired_outcomes.csv (16 rows, exact)", flush=True)
    ga = {"lg16": pd.read_csv(os.path.join(HERE, "rev4_galg_s0of1.csv")),
          "lg16b": pd.read_csv(os.path.join(HERE, "rev4_galgb_s0of1.csv"))}
    for fam in fr:
        g = ga[fam]
        assert len(g) == 30 and set(g.Algorithm) == {"GA"} and sorted(g.Seed) == list(range(1, 31))
        assert int(g.Budget.iloc[0]) == int(fr[fam].Budget.iloc[0]) and (g.Calls == g.Budget).all()
        assert np.allclose(g.Ideal, fr[fam].Ideal.iloc[0], rtol=0, atol=1e-9)
    fr2 = {fam: pd.concat([fr[fam], ga[fam][fr[fam].columns]], ignore_index=True) for fam in fr}
    AA.LABELS.setdefault("GA", "GA")                                # label only (in memory)
    rows, pairs = AA.reliability(fr2)
    assert len(pairs) == 18
    adj = AA.holm([p["exact_sign_p"] for p in pairs])
    for p, v in zip(pairs, adj):
        p["holm18_p"] = float(v); del p["holm16_p"]                 # reliability() names the Holm column holm16_p
        p["cp95_conditional"] = RD.cp(p["wins"], p["wins"] + p["losses"])
    old = {(p["budget"], p["second"]): p for p in pairs0}
    changes = []
    for p in pairs:
        o = old.get((p["budget"], p["second"]))
        if o is not None:
            assert (o["wins"], o["ties"], o["losses"], o["superiority"], o["bootstrap95_lower"], o["bootstrap95_upper"]) == \
                   (p["wins"], p["ties"], p["losses"], p["superiority"], p["bootstrap95_lower"], p["bootstrap95_upper"])
            changes.append(dict(budget=p["budget"], second=p["second"], label=p["label"], holm16_p=o["holm16_p"],
                                holm18_p=p["holm18_p"], sig16=o["holm16_p"] < 0.05, sig18=p["holm18_p"] < 0.05))
    gap = [p for p in pairs if p["second"] == "GA"]
    garows = [r for r in rows if r["algorithm"] in ("GA", "PSOBV")]
    for p in gap:
        print(f"LG GA {p['budget']}: {p['wins']}/{p['ties']}/{p['losses']} score {p['superiority']:.4f} "
              f"[{p['bootstrap95_lower']:.3f}, {p['bootstrap95_upper']:.3f}] p {p['exact_sign_p']:.4g} Holm18 {p['holm18_p']:.4g} "
              f"CP {np.round(p['cp95_conditional'], 3)}", flush=True)
    for c in changes:
        print(f"  Holm 16->18 {c['budget']} {c['label']:8s} {c['holm16_p']:.4g} -> {c['holm18_p']:.4g}", flush=True)
    return dict(reproduction="audit_archive.reliability on rev2_lg16*/rev2_lg16b* equals SWEVO_rev2/validation/"
                "paired_outcomes.csv exactly (16 rows)", ga_pairs=gap, reliability_rows=garows, holm16_to_holm18=changes,
                all_pairs_holm18=pairs)


# ------------------------------------------------------------------ precision audit
def precision_block():
    out = {}
    # GA benchmark records: part of the existing audit; re-check by recomputation
    rec = pd.read_csv(os.path.join(HERE, "rev3_precision_audit_records.csv"))
    st = rec[rec.Study == "rev2_ga"]
    re_rows = []
    for f in sorted(glob.glob(os.path.join(HERE, "rev2_ga_s*of10.csv"))):
        re_rows += PA.audit_file(f, "rev2_ga", "used")
    R = pd.DataFrame(re_rows)
    m = st.merge(R[["File", "Row", "Class", "ClassHybrid", "Feasible"]], on=["File", "Row"], suffixes=("", "_re"), validate="1:1")
    assert len(m) == 2040 and (m.Class == m.Class_re).all() and (m.Feasible == m.Feasible_re).all()
    fp = pd.read_csv(os.path.join(HERE, "rev3_fullprec_rev2_ga.csv"))
    fpl = pd.read_csv(os.path.join(HERE, "rev3_fullprec_rev2_ga_local.csv"))
    iii = st[st.Class == "iii"][["File", "Row"]].sort_values(["File", "Row"]).values.tolist()
    assert sorted(fp[["File", "Row"]].values.tolist()) == iii == sorted(fpl[["File", "Row"]].values.tolist())
    aj = json.load(open(os.path.join(HERE, "rev3_precision_audit.json")))
    out["benchmark_rev2_ga"] = dict(
        in_existing_audit=True, study="rev2_ga (rev3_precision_audit.STUDIES)", records=int(len(st)),
        labelled_feasible=int(st.Feasible.sum()), class_i=int((st.Class == "i").sum()), class_ii=int((st.Class == "ii").sum()),
        class_iii=int((st.Class == "iii").sum()), recomputed_classes_equal=True,
        json_summary={c: aj["studies"]["rev2_ga"][c] for c in ("records", "class_i", "class_ii", "class_iii")},
        class_iii_reruns=[dict(File=r.File, Row=int(r.Row), Dataset=str(r.Dataset), Radius=int(r.Radius), Turbines=int(r.Turbines),
                               Seed=int(r.Seed), StoredLabel=bool(r.StoredLabel), RerunLabel=bool(r.RerunLabel),
                               LabelChanged=bool(r.LabelChanged), LabelDecided=bool(r.LabelDecided), Reproduction=r.Reproduction,
                               Margin_m=float(r.Margin)) for r in fp.itertuples()],
        class_iii_all_decided_unchanged=bool(fp.LabelDecided.all() and not fp.LabelChanged.any()),
        local_rerun_all_match=bool(fpl.AllMatch.all()))
    # Lillgrund GA records (not in the existing audit): same audit_file
    lg = {}
    rows = []
    for st_id, fn in (("rev4_galg", "rev4_galg_s0of1.csv"), ("rev4_galgb", "rev4_galgb_s0of1.csv")):
        r = PA.audit_file(os.path.join(HERE, fn), st_id, "used")
        rows += r
    D = pd.DataFrame(rows)
    summ = PA.summarize(D)
    for st_id, s in summ.items():
        g = D[D.Study == st_id]
        s["min_decimals"] = int(g.Decimals.min())
        s["min_abs_slack_spacing_m"] = float(g.SlackSpacing.abs().min())
        s["min_abs_slack_boundary_m"] = float(g.SlackBoundary.abs().min())
        s["max_bound_spacing_m"] = float(g.BoundSpacing.max())
        s["max_bound_boundary_m"] = float(g.BoundBoundary.max())
        lg[st_id] = s
    out["lillgrund_rev4_galg"] = dict(in_existing_audit=False, studies=lg,
                                      records=int(len(D)), class_i=int((D.Class == "i").sum()),
                                      class_ii=int((D.Class == "ii").sum()), class_iii=int((D.Class == "iii").sum()),
                                      reruns_needed=int((D.Class == "iii").sum()))
    tot = aj["total"]
    out["audit_totals_existing"] = dict(records=tot["records"], class_i=tot["class_i"], class_iii=tot["class_iii"])
    out["audit_totals_with_lillgrund_ga"] = dict(records=tot["records"] + int(len(D)),
                                                 class_i=tot["class_i"] + int((D.Class == "i").sum()),
                                                 class_ii=tot["class_ii"] + int((D.Class == "ii").sum()),
                                                 class_iii=tot["class_iii"] + int((D.Class == "iii").sum()))
    print("precision: rev2_ga", {k: out["benchmark_rev2_ga"][k] for k in ("records", "class_i", "class_ii", "class_iii")},
          "| rev4_galg*", {k: out["lillgrund_rev4_galg"][k] for k in ("records", "class_i", "class_ii", "class_iii")}, flush=True)
    for st_id, s in lg.items():
        print(f"  {st_id}: decimals {s['decimals']} feasible {s['labelled_feasible']} classes {s['class_i']}/{s['class_ii']}/"
              f"{s['class_iii']} min|slack| sp {s['min_abs_slack_spacing_m']:.3g} bd {s['min_abs_slack_boundary_m']:.3g} m; "
              f"max replay diff {s['max_abs_replay_diff']:.3g} GWh ({s['max_abs_replay_diff_pp']:.3g} pp)", flush=True)
    return out


# ------------------------------------------------------------------ LaTeX rows (tab:allrun columns and supplement)
def lat(R):
    b = R["benchmark"]["ga_pair"]; m = b["pooled"]["main"]; s = b["seed"]["seed"]; g = b["seed"]["group"]
    adj = R["benchmark"]["adjusted_main8"]["PSOBV-GA"]
    L = {}
    L["tab_allrun_benchmark"] = (f"PSO-VNS vs.\\ GA & {m['W']}/{m['T']}/{m['L']} & {m['score']:.3f} & {RI.ci(m['ci95']['case'])} & "
                                 f"{RI.ptex(adj['pooled_sign_holm'])} & {RD.fci(s['ci95_joint_seed'])} & {RD.ptex(adj['signflip_holm'])} \\\\")
    def p3(p):
        if p >= 0.001:
            return f"{p:.3f}"
        m_, e = f"{p:.2e}".split("e")
        return f"${m_}\\times10^{{{int(e)}}}$"
    for p in R["lillgrund"]["ga_pairs"]:
        L[f"tab_allrun_lillgrund_{p['budget']}"] = (
            f"PSO-VNS vs.\\ GA, {p['budget']:,} & {p['wins']}/{p['ties']}/{p['losses']} & {p['superiority']:.3f} & -- & -- & "
            f"[{p['bootstrap95_lower']:.3f}, {p['bootstrap95_upper']:.3f}] & {p3(p['holm18_p'])} \\\\")
        L[f"S_archive_paired_{p['budget']}"] = (
            f"{p['budget']:,} & GA & {p['wins']}/{p['ties']}/{p['losses']} & {p['superiority']:.3f} & "
            f"[{p['bootstrap95_lower']:.3f}, {p['bootstrap95_upper']:.3f}] & "
            f"[{p['cp95_conditional'][0]:.3f}, {p['cp95_conditional'][1]:.3f}] & {p['holm18_p']:.3g} \\\\")
    v = b["pooled"]["violation"]
    L["S_sa_allrun"] = (f"PSO-VNS vs.\\ GA & {b['pooled']['feas_a_pct']:.1f}/{b['pooled']['feas_b_pct']:.1f} & "
                        f"{RI.f3(b['pooled']['cond_mean_dloss_pp'], True)} & {m['W']}/{m['T']}/{m['L']} & {m['score']:.3f} & "
                        f"{RI.ci(m['ci95']['case'])} & {RI.ci(m['ci95']['two_stage'])} & {RI.ptex(adj['pooled_sign_holm'])} & "
                        f"{v['score']:.3f} \\\\")
    L["S_sa_dep_bench"] = (f"PSO-VNS vs.\\ GA & {b['seed']['W']}/{b['seed']['T']}/{b['seed']['L']} & {s['score']:.3f} & "
                           f"{RD.fci(s['ci95_joint_seed'])} & {s['seeds_above']}/{s['seeds_below']} & {RD.ptex(adj['signflip_holm'])} & "
                           f"{RD.ptex(adj['wilcoxon_holm'])} & {g['groups_above']}/{g['groups_below']} & {RD.ptex(adj['groups_holm'])} \\\\")
    return L


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", default=HERE)
    ap.add_argument("--b1-boot", type=int, default=RI.B1_BOOT)
    a = ap.parse_args(argv)
    t0 = time.time()
    R = dict(generated=time.strftime("%Y-%m-%d %H:%M:%S"), script="rev5_ga_allrun.py",
             rule="feasible beats infeasible; two feasible runs by the objective (tie |d| <= 1e-9); two infeasible runs tie; "
                  "score (W + T/2)/n",
             sources=dict(benchmark_ga="rev2_ga_s*of10.csv", lillgrund_ga=["rev4_galg_s0of1.csv", "rev4_galgb_s0of1.csv"]))
    R["benchmark"] = bench_block(a.b1_boot)
    R["lillgrund"] = lillgrund_block()
    R["precision"] = precision_block()
    R["latex_rows"] = lat(R)
    R["seconds"] = round(time.time() - t0, 1)
    json.dump(clean(R), open(os.path.join(a.out_dir, "rev5_ga_allrun.json"), "w"), indent=1)
    for k, v in R["latex_rows"].items():
        print(f"{k}: {v}")
    print(f"done in {R['seconds']} s", flush=True)


if __name__ == "__main__":
    main()

"""Check the outcome-dependent statements of MPCE_PSO_VNS.tex against mpce_summary.json.

Usage (from analysis/):  python3 mpce_check_final.py [--summary mpce_summary.json] [--tex ../MPCE_PSO_VNS.tex]

Every qualitative statement in the manuscript that was written from the data (e.g. "significantly better",
"never worse", "ranks second under the Gaussian model") carries a comment  % CHECK-FINAL [Cnn]: <condition>.
This script evaluates the same conditions (CONDITIONS below) and prints PASS / FAIL / PENDING (data missing)
for each; it also reports IDs used in the manuscript but not defined here, and vice versa. A FAIL means the
sentence next to the comment must be rewritten. Numbers themselves come from mpce_numbers.tex and need no check.
mpce_results.py runs this script at the end of every run. Exit code 1 if any condition fails.

Reference scan (not evaluated here): the comments % CHECK-EXTRA [Xnn] (mpce_check_extra.py), % CHECK-DIAG [Dnn]
(mpce_check_diag.py) and % CHECK-THEORY [Tnn] (make_theory_figures.py --check) of the manuscript and optA/*.tex are
listed with the ids they reference, and ids that their script does not define are flagged (the ids a script
defines are read from its source text; the scripts themselves evaluate the conditions). A tag covers every
bracketed id of its letter up to the next CHECK- keyword of the comment, incl. lists ([C19]/[C30]) and ranges
([X02]-[X09]).
"""
import os, re, sys, json, argparse

HERE = os.path.dirname(os.path.abspath(__file__))


def g(d, *path):
    """Nested lookup; raises KeyError (-> PENDING) when a key or section is missing."""
    for p in path:
        if d is None:
            raise KeyError("/".join(map(str, path)))
        d = d[p]
    if d is None:
        raise KeyError("/".join(map(str, path)))
    return d


def wtl(d):
    return int(d.get("W", 0)), int(d.get("T", 0)), int(d.get("L", 0))


def rankpos(ar, a):
    return 1 + sorted(ar.values()).index(ar[a])


# (id, statement in the text, keys, condition)
CONDITIONS = [
    ("C01", "PSO-VNS has the best average rank of the eight methods (abstract, Sec. VI-A, conclusion)",
     "main.friedman.best_ranked", lambda s: g(s, "main", "friedman", "best_ranked") == "PSOBV"),
    ("C02", "the Friedman test rejects equal performance", "main.friedman.p", lambda s: g(s, "main", "friedman", "p") < 0.05),
    ("C03", "post hoc: PSO-VNS significantly better than every method except PSO (grammar of the sentence assumes exactly PSO is non-significant)",
     "main.friedman.p_holm_vs_focus", lambda s: all((p >= 0.05) == (a == "PSOC") for a, p in g(s, "main", "friedman", "p_holm_vs_focus").items())),
    ("C04", "vs PSO: at least as good (W >= L, lower mean loss) but not significant over all cases (case-mean Wilcoxon p >= 0.05)",
     "main.wtl.PSOC, main.mean_loss_pct_qualified_cases, main.case_mean_wilcoxon.PSOC.p",
     lambda s: wtl(g(s, "main", "wtl", "PSOC"))[0] >= wtl(g(s, "main", "wtl", "PSOC"))[2]
     and g(s, "main", "mean_loss_pct_qualified_cases", "PSOBV") < g(s, "main", "mean_loss_pct_qualified_cases", "PSOC")
     and g(s, "main", "case_mean_wilcoxon", "PSOC", "p") >= 0.05),
    ("C05", "vs PSO: gain concentrates at N >= 10 (W > L there, mean gain > 0) and nearly equal loss for N < 10 (mean |dL| < 0.05 pp)",
     "main.by_n.PSOC", lambda s: (lambda b: wtl(b["wtl_large"])[0] > wtl(b["wtl_large"])[2] and b["mean_dloss_pp_large"] < 0
                                  and b["mean_abs_dloss_pp_small"] < 0.05)(g(s, "main", "by_n", "PSOC"))),
    ("C06", "vs SSA-VNS: never significantly worse, significant on the case means",
     "main.wtl.SSABV.L, main.case_mean_wilcoxon.SSABV.p_holm",
     lambda s: wtl(g(s, "main", "wtl", "SSABV"))[2] == 0 and g(s, "main", "case_mean_wilcoxon", "SSABV", "p_holm") < 0.05),
    ("C07", "06 sentence 'Against SSA, LX-SSA and DE, it is significantly better in ... cases, and against VNS and MS-SLSQP in ...' "
            "(re-pointed; attach [C07] there): each W >= 40, L <= 2, case-mean p_holm < 0.05",
     "main.wtl.*, main.case_mean_wilcoxon.*.p_holm",
     lambda s: all(wtl(g(s, "main", "wtl", a))[0] >= 40 and wtl(g(s, "main", "wtl", a))[2] <= 2
                   and g(s, "main", "case_mean_wilcoxon", a, "p_holm") < 0.05 for a in ("SSA", "LXSSA", "DE", "BVNS", "SLSQP"))),
    ("C08", "PSO-VNS returns a feasible layout in every run (68 cases and Horns Rev 16)",
     "main.feasible_pct.PSOBV, hr16.methods.PSOBV.feasible",
     lambda s: g(s, "main", "feasible_pct", "PSOBV") == 100 and g(s, "hr16", "methods", "PSOBV", "feasible") == g(s, "hr16", "methods", "PSOBV", "runs")),
    ("C09", "feasible-run counts of PSO, SSA, LX-SSA and DE are identical for both wind data sets",
     "main.feasible_runs_by_dataset", lambda s: all(v["1"] == v["2"] for a, v in g(s, "main", "feasible_runs_by_dataset").items()
                                                   if a in ("PSOC", "SSA", "LXSSA", "DE"))),
    ("C10", "VNS phase removes (slightly) more of the loss left at the switch than continuing PSO, and both keep improving",
     "ablation.phase2_loss_reduction_pct.PSOBV/PSOC_continued.mean",
     lambda s: g(s, "ablation", "phase2_loss_reduction_pct", "PSOBV", "mean") >= g(s, "ablation", "phase2_loss_reduction_pct", "PSOC_continued", "mean") > 0),
    ("C11", "VNS phases of SSA-VNS / RS-VNS remove larger shares but start from worse or more often infeasible layouts; VNS alone ranks below every hybrid and RS-VNS",
     "ablation.phase2_loss_reduction_pct.*, ablation.friedman.avg_rank",
     lambda s: (lambda p, ar: p["SSABV"]["mean"] > p["PSOBV"]["mean"] and p["RSVNS"]["mean"] > p["PSOBV"]["mean"]
                and p["PSOBV"]["feasible_at_switch_pct"] > p["SSABV"]["feasible_at_switch_pct"] > p["RSVNS"]["feasible_at_switch_pct"]
                and p["PSOBV"]["mean_loss_at_switch_pct"] < p["SSABV"]["mean_loss_at_switch_pct"]
                and all(ar["BVNS"] > ar[h] for h in ("PSOBV", "SSABV", "LXBV", "RSVNS")))(
         g(s, "ablation", "phase2_loss_reduction_pct"), g(s, "ablation", "friedman", "avg_rank"))),
    ("C12", "old PSO setting (w=0.7, c1=c2=2) understates PSO: constriction PSO never significantly worse, better in most cases, "
            "moves PSO up to second place, ahead of SSA, LX-SSA and SSA-VNS",
     "pso_setting.*, main.friedman.avg_rank",
     lambda s: (lambda p, ar: wtl(p["constriction_vs_old_wtl"])[2] == 0 and wtl(p["constriction_vs_old_wtl"])[0] >= 34
                and p["constriction_rank_position"] == 2 and p["old_rank_position"] > p["constriction_rank_position"]
                and all(ar["PSOC"] < ar[a] for a in ("SSA", "LXSSA", "SSABV")))(g(s, "pso_setting"), g(s, "main", "friedman", "avg_rank"))),
    ("C13", "ablation: VNS phase gives a small gain (PSO-VNS vs PSO W > L)",
     "ablation.contrasts.PSOBV-PSOC", lambda s: wtl(g(s, "ablation", "contrasts", "PSOBV-PSOC"))[0] > wtl(g(s, "ablation", "contrasts", "PSOBV-PSOC"))[2]),
    ("C14", "the PSO start beats the best initial point (VNS) and random sampling (RS-VNS): W >= 40, L <= 2 / L == 0",
     "ablation.contrasts.PSOBV-BVNS, PSOBV-RSVNS",
     lambda s: wtl(g(s, "ablation", "contrasts", "PSOBV-BVNS"))[0] >= 40 and wtl(g(s, "ablation", "contrasts", "PSOBV-BVNS"))[2] <= 2
     and wtl(g(s, "ablation", "contrasts", "PSOBV-RSVNS"))[0] >= 40 and wtl(g(s, "ablation", "contrasts", "PSOBV-RSVNS"))[2] == 0),
    ("C15", "SSA phase gives only a small gain over random sampling (D1): SSA-VNS vs RS-VNS case-mean p < 0.05 with SSA-VNS lower "
            "(dL < 0, |dL| < 0.1 pp and less than half of |dL| of PSO-VNS vs RS-VNS), tied at run level in most cases (T > 34 of 68)",
     "ablation.case_mean.SSABV-RSVNS / PSOBV-RSVNS, ablation.contrasts.SSABV-RSVNS",
     lambda s: (lambda c, cp, w: c["p"] < 0.05 and c["mean_dloss_pp"] < 0 and abs(c["mean_dloss_pp"]) < 0.1
                and abs(c["mean_dloss_pp"]) < 0.5 * abs(cp["mean_dloss_pp"]) and w[1] > 34)(
         g(s, "ablation", "case_mean", "SSABV-RSVNS"), g(s, "ablation", "case_mean", "PSOBV-RSVNS"),
         wtl(g(s, "ablation", "contrasts", "SSABV-RSVNS")))),
    ("C16", "PSO is the better first phase than SSA / LX-SSA (never worse, W >= 34) and its best layout at the switch is more often feasible and has lower loss",
     "ablation.contrasts.PSOBV-SSABV/PSOBV-LXBV, ablation.phase2_loss_reduction_pct",
     lambda s: all(wtl(g(s, "ablation", "contrasts", k))[2] == 0 and wtl(g(s, "ablation", "contrasts", k))[0] >= 34 for k in ("PSOBV-SSABV", "PSOBV-LXBV"))
     and (lambda p: all(p["PSOBV"]["feasible_at_switch_pct"] > p[h]["feasible_at_switch_pct"] and p["PSOBV"]["mean_loss_at_switch_pct"] < p[h]["mean_loss_at_switch_pct"]
                        for h in ("SSABV", "LXBV")))(g(s, "ablation", "phase2_loss_reduction_pct"))),
    ("C17", "LX-SSA as published is worse (D1, D14; Table IV rows 'LX-SSA as published (hybrid)' / '(alone)'): LX-SSA never "
            "significantly better than SSA and LX-SSA-VNS never significantly better than SSA-VNS (run level); on the case means "
            "SSA-VNS beats LX-SSA-VNS and SSA beats LX-SSA (p < 0.05)",
     "ablation.contrasts.LXSSA-SSA.W, SSABV-LXBV.L, ablation.case_mean.SSABV-LXBV / LXSSA-SSA",
     lambda s: wtl(g(s, "ablation", "contrasts", "LXSSA-SSA"))[0] == 0 and wtl(g(s, "ablation", "contrasts", "SSABV-LXBV"))[2] == 0
     and g(s, "ablation", "case_mean", "SSABV-LXBV", "p") < 0.05 and g(s, "ablation", "case_mean", "SSABV-LXBV", "mean_dloss_pp") < 0
     and g(s, "ablation", "case_mean", "LXSSA-SSA", "p") < 0.05 and g(s, "ablation", "case_mean", "LXSSA-SSA", "mean_dloss_pp") > 0),
    ("C18", "split: 75% has lower loss than 50% in >= 10 of 12 cases, 50% never significantly better than 75%, 25% never significantly better than 50%, ranks 75 < 50 < 25",
     "split.n_lower_loss_than_50, split.wtl_50_vs_75, split.wtl_50_vs_25, split.avg_rank",
     lambda s: (lambda p: p["n_lower_loss_than_50"]["PSOBV75"] >= 10 and p["wtl_50_vs_75"]["W"] == 0 and p["wtl_50_vs_25"]["L"] == 0
                and p["avg_rank"]["PSOBV75"] < p["avg_rank"]["PSOBV"] < p["avg_rank"]["PSOBV25"])(g(s, "split"))),
    ("C19", "Horns Rev 16 (6,030, random starts, the paper's 5-deg bins): PSO-VNS has the highest mean AEP and is significantly "
            "better than every other method in the run-level test (DATA condition; D19: the AEP edge vanishes with 1-deg bins "
            "and PyWake, see mpce_check_dir.py F24/F25, so the text states the 6,030 result as a feasibility result)",
     "hr16.methods.*.mean/p_holm",
     lambda s: (lambda h: all(h["PSOBV"]["mean"] > v["mean"] for a, v in h.items() if a != "PSOBV" and v.get("mean") is not None)
                and all(v["p_holm"] < 0.05 for a, v in h.items() if a != "PSOBV"))(g(s, "hr16", "methods"))),
    ("C20", "Horns Rev 16: PSO finds feasible layouts in only part of its random-start runs (fewer than 30 of 30)",
     "hr16.methods.PSOC.feasible", lambda s: g(s, "hr16", "methods", "PSOC", "feasible") < g(s, "hr16", "methods", "PSOC", "runs")),
    ("C21", "ranking stable under the cubic curves: tau >= 0.9 and PSO-VNS keeps the best rank",
     "robustness.Cubic/CubicCutout.tau, avg_rank",
     lambda s: all(g(s, "robustness", m, "tau") >= 0.9 and rankpos(g(s, "robustness", m, "avg_rank"), "PSOBV") == 1 for m in ("Cubic", "CubicCutout"))),
    ("C22", "Gaussian wake: ranking less stable (tau < 0.8, same best < 75%), PSO ranks first and PSO-VNS second (both ahead of the rest)",
     "robustness.Gauss.tau, same_best_pct, avg_rank",
     lambda s: (lambda r: r["tau"] < 0.8 and r["same_best_pct"] < 75 and rankpos(r["avg_rank"], "PSOC") == 1 and rankpos(r["avg_rank"], "PSOBV") == 2)(g(s, "robustness", "Gauss"))),
    ("C23", "optimized layouts often lie on the 4D constraint (fewer than 10% of the N = 11-15 layouts satisfy 5D)",
     "main.spacing_share.n11_15.pct_5d", lambda s: g(s, "main", "spacing_share", "n11_15", "pct_5d") < 10),
    ("C24", "densest cases: the VNS phase repairs the PSO global best that is often infeasible at the switch (feasible at switch < final)",
     "main.dense_500_10", lambda s: all(v["feasible_at_switch_pct"] < v["feasible_final_pct"] for v in g(s, "main", "dense_500_10").values())),
    ("C25", "Horns Rev 16: larger budgets lower the mean AEP loss of PSO-VNS (30,030 and 120,030 below 6,030)",
     "hr16.loss_by_setting.PSOBV", lambda s: (lambda h: h["30030R"]["loss"] < h["6030R"]["loss"] and h["120030R"]["loss"] < h["6030R"]["loss"])(
         g(s, "hr16", "loss_by_setting", "PSOBV"))),
    ("C26", "budget 6,030 on the six largest cases: PSO-VNS ranks first (column 6,030 of Table feasbudget)",
     "feasbudget.rank.6030", lambda s: min(g(s, "feasbudget", "rank", "6030"), key=g(s, "feasbudget", "rank", "6030").get) == "PSOBV"),
    ("C27", "previous study's method pool (Table baseline): under the old PSO setting BOTH salp-swarm hybrids (SSA-VNS and "
            "LX-SSA-VNS) rank ahead of PSO; under the corrected setting PSO ranks ahead of every salp-swarm method (SSA, LX-SSA, SSA-VNS, LX-SSA-VNS)",
     "baseline.old.avg_rank, baseline.constriction.avg_rank",
     lambda s: (lambda o, c: all(o[h] < o["PSO"] for h in ("SSABV", "LXBV"))
                and all(c["PSO"] < c[a] for a in ("SSA", "LXSSA", "SSABV", "LXBV")))(
         g(s, "baseline", "old", "avg_rank"), g(s, "baseline", "constriction", "avg_rank"))),
    ("C28", "feasibility-preserving initialization, six largest cases: VNS alone (BVNS) has the best average rank and PSO-VNS ranks second",
     "feasbudget.rank_feasible_init",
     lambda s: (lambda r: rankpos(r, "BVNS") == 1 and min(r, key=r.get) == "BVNS" and sum(v <= r["PSOBV"] for v in r.values()) == 2)(
         g(s, "feasbudget", "rank_feasible_init"))),
    ("C29", "six largest cases: PSO-VNS has the best average rank at all three budgets (6,030, 30,030, 120,030)",
     "feasbudget.rank.{6030,30030,120030}",
     lambda s: all((lambda r: all(r["PSOBV"] < v for a, v in r.items() if a != "PSOBV"))(g(s, "feasbudget", "rank", b))
                   for b in ("6030", "30030", "120030"))),
    ("C30", "Horns Rev 1 at 30,030 evaluations (5-deg bins): PSO-VNS has the highest mean AEP (reworded, D16/D19: the former "
            "second part 'at least one of its runs exceeds the installed layout' is dropped -- no run does with 1-deg bins, F21 -- "
            "and is no longer evaluated; see C62; survival at 1 deg / PyWake: F26)",
     "hr16.loss_by_setting.*.30030R.mean_aep",
     lambda s: (lambda L: all(L["PSOBV"]["30030R"]["mean_aep"] > v["30030R"]["mean_aep"] for a, v in L.items()
                              if a != "PSOBV" and v and v.get("30030R") and v["30030R"].get("mean_aep") is not None))(
         g(s, "hr16", "loss_by_setting"))),
]



# Phase-4 conditions (optA/PHASE4.md D1-D7; R5 proposals C31-C43, with the meanings of C37-C40 set by the 08
# section and of C43 by the lead). Conditions whose sentences are being rewritten evaluate DATA only.
def _cm(s, k):
    return g(s, "ablation", "case_mean", k)


def _wtl_abl(s, k):
    return wtl(g(s, "ablation", "contrasts", k))


def _hr_mean(s, a, tag):
    return g(s, "hr16", "loss_by_setting", a, tag, "mean_aep")


CONDITIONS += [
    ("C31", "Horns Rev 1, 6,030 evaluations (5-deg bins): the PSO-VNS mean AEP is below the installed block (DATA condition; "
            "neither the mean nor any run is claimed to exceed it, C62)",
     "hr16.methods.PSOBV.mean, hr16.installed_aep",
     lambda s: g(s, "hr16", "methods", "PSOBV", "mean") < g(s, "hr16", "installed_aep")),
    ("C32", "PSO-VNS vs PSO (D2): not significant over all 68 cases (case-mean p >= 0.05); Data Set II, N >= 10: PSO-VNS lower "
            "case-mean loss in more cases than PSO (wins > losses, mean dL < 0); Data Set I, N >= 10: PSO lower in at least half",
     "main.case_mean_wilcoxon.PSOC.p, main.subgroup_vs_phase1.groups.dsIILarge/dsILarge",
     lambda s: (lambda G2, G1: g(s, "main", "case_mean_wilcoxon", "PSOC", "p") >= 0.05
                and G2["wins"] > G2["losses"] and G2["mean_dloss_pp"] < 0 and G1["losses"] >= G1["n_cases"] / 2)(
         g(s, "main", "subgroup_vs_phase1", "groups", "dsIILarge"), g(s, "main", "subgroup_vs_phase1", "groups", "dsILarge"))),
    ("C33", "VNS phase (D1): significant for SSA and LX-SSA (run level W > L with L = 0; case-mean p < 0.05, dL < 0), "
            "not significant for PSO (case-mean p >= 0.05)",
     "ablation.contrasts/case_mean SSABV-SSA, LXBV-LXSSA, PSOBV-PSOC",
     lambda s: all(_wtl_abl(s, k)[0] > _wtl_abl(s, k)[2] == 0 and _cm(s, k)["p"] < 0.05 and _cm(s, k)["mean_dloss_pp"] < 0
                   for k in ("SSABV-SSA", "LXBV-LXSSA")) and _cm(s, "PSOBV-PSOC")["p"] >= 0.05),
    ("C34", "budget (D2): PSO-VNS has the best average rank on the six largest cases at all three budgets; first-place counts are "
            "generated for every budget and PSO-VNS is first in at least one case at 120,030",
     "feasbudget.rank, feasbudget.first_count",
     lambda s: all((lambda r: all(r["PSOBV"] < v for a, v in r.items() if a != "PSOBV"))(g(s, "feasbudget", "rank", b))
                   for b in ("6030", "30030", "120030"))
     and all(g(s, "feasbudget", "first_count", b, "PSOBV") >= 0 for b in ("6030", "30030", "120030"))
     and g(s, "feasbudget", "first_count", "120030", "PSOBV") >= 1),
    ("C35", "LX-SSA phase gives no gain over random sampling (D1): LX-SSA-VNS vs RS-VNS case-mean p >= 0.05, run level T >= 60, W and L <= 5",
     "ablation.case_mean.LXBV-RSVNS.p, ablation.contrasts.LXBV-RSVNS",
     lambda s: _cm(s, "LXBV-RSVNS")["p"] >= 0.05 and (lambda w: w[1] >= 60 and w[0] <= 5 and w[2] <= 5)(_wtl_abl(s, "LXBV-RSVNS"))),
    ("C36", "component analysis: the Friedman test rejects equal average ranks (\\NAblFriedVerb)",
     "ablation.friedman.p", lambda s: g(s, "ablation", "friedman", "p") < 0.05),
    ("C37", "budget, 120,030 evaluations (08): SSA-VNS ranks ahead of PSO on the six largest cases, and VNS, SSA-VNS, RS-VNS, "
            "LX-SSA-VNS and DE are all within 0.2 pp (\\NBudCloseGapOneTwentyK) of the PSO-VNS mean wake loss",
     "feasbudget.rank.120030, feasbudget.close_gap_pp.120030",
     lambda s: (lambda r, c: r["SSABV"] < r["PSOC"] and set(c["methods"]) == {"BVNS", "SSABV", "RSVNS", "LXBV", "DE"}
                and 0 <= c["max_gap"] < 0.2)(g(s, "feasbudget", "rank", "120030"), g(s, "feasbudget", "close_gap_pp", "120030"))),
    ("C38", "feasible starts (08): VNS alone has the lowest mean wake loss and the best average rank on the six largest cases; "
            "PSO and DE never move from their best initial layout (final objective = first checkpoint in every run)",
     "feasbudget.feasible.*.loss, feasbudget.rank_feasible_init, feasinit.runs_moved_from_best_initial",
     lambda s: (lambda f, r, m: min((a for a in f if f[a] and f[a].get("loss") is not None), key=lambda a: f[a]["loss"]) == "BVNS"
                and min(r, key=r.get) == "BVNS" and m["PSOC"]["moved_grid6"] == 0 and m["DE"]["moved_grid6"] == 0
                and m["PSOC"]["runs_grid6"] > 0 and m["DE"]["runs_grid6"] > 0)(
         g(s, "feasbudget", "feasible"), g(s, "feasbudget", "rank_feasible_init"), g(s, "feasinit", "runs_moved_from_best_initial"))),
    ("C39", "Horns Rev 1 (08): PSO-VNS has the highest mean AEP of all methods at 6,030 and 30,030 evaluations, but not at 120,030",
     "hr16.loss_by_setting.*.{6030R,30030R,120030R}.mean_aep",
     lambda s: (lambda lb: all(max((a for a in lb if lb[a] and lb[a].get(t) and lb[a][t].get("mean_aep") is not None),
                                   key=lambda a: lb[a][t]["mean_aep"]) == "PSOBV" for t in ("6030R", "30030R"))
                and max((a for a in lb if lb[a] and lb[a].get("120030R") and lb[a]["120030R"].get("mean_aep") is not None),
                        key=lambda a: lb[a]["120030R"]["mean_aep"]) != "PSOBV")(g(s, "hr16", "loss_by_setting"))),
    ("C40", "IEA37 CS1 (08): the best layout of all methods (both budgets) is below the best feasible published layout in both "
            "scenarios, and the gap is larger for 36 turbines",
     "iea37.*.published.our_best_vs_best_feasible_pct",
     lambda s: (lambda i: all(v["published"]["our_best_vs_best_feasible_pct"] < 0 for v in i.values())
                and abs(i["IEA37 CS1, 36 turbines"]["published"]["our_best_vs_best_feasible_pct"])
                > abs(i["IEA37 CS1, 16 turbines"]["published"]["our_best_vs_best_feasible_pct"]))(g(s, "iea37"))),
    ("C41", "boundary rule (09): SSA with radial projection has the higher mean objective than the original SSA in all six "
            "cases (Mann-Whitney p < 0.01 in every case)",
     "boundary_rule", lambda s: (lambda b: b["n_higher"] == b["n_cases"] == 6 and b["max_p"] < 0.01)(g(s, "boundary_rule"))),
    ("C42", "quality (06): the PSO-VNS case-mean wake loss grows with N, falls as the farm becomes larger and is higher for "
            "Data Set II than for Data Set I (ties allowed)",
     "main.loss_monotone",
     lambda s: (lambda m: not m["violations_n"] and not m["violations_r"] and not m["violations_ds"] and m["ds2_higher"] > 0)(
         g(s, "main", "loss_monotone"))),
    ("C43", "seed pairing: every method has the same seed set in each case x budget x initialization, and the first logged "
            "(best initial) objective is identical across all methods for every case-seed pair at 6,030 evaluations where finite",
     "pairing", lambda s: (lambda p: p["groups_with_different_seed_sets"] == 0 and p["pairs_first_checkpoint_differs"] == 0
                           and p["pairs_compared"] > 0)(g(s, "pairing"))),
    ("C44", "PSO phase (D1): a PSO phase gives a clearly larger gain over random sampling than an SSA phase (PSO-VNS vs RS-VNS "
            "case-mean p < 0.05 and dL at least twice that of SSA-VNS vs RS-VNS)",
     "ablation.case_mean.PSOBV-RSVNS, SSABV-RSVNS",
     lambda s: _cm(s, "PSOBV-RSVNS")["p"] < 0.05 and _cm(s, "PSOBV-RSVNS")["mean_dloss_pp"] < 2 * _cm(s, "SSABV-RSVNS")["mean_dloss_pp"] < 0),
    ("C45", "budget (D3): at 120,030 evaluations SSA-VNS ranks ahead of PSO on the six largest cases",
     "feasbudget.rank.120030", lambda s: g(s, "feasbudget", "rank", "120030", "SSABV") < g(s, "feasbudget", "rank", "120030", "PSOC")),
    ("C46", "cost (D6): the metaheuristics need a similar time per evaluation (max/min < 2.5) and MS-SLSQP needs more than every "
            "metaheuristic (slowdown vs PSO > 1)",
     "cost_per_eval", lambda s: (lambda c: c["meta_max_ms"] / c["meta_min_ms"] < 2.5 and c["slsqp_ms"] > c["meta_max_ms"]
                                 and c["slsqp_over_pso"] > 1)(g(s, "cost_per_eval"))),
    ("C47", "budget split (D7): omega = 0.75 is the best of the tested settings (lowest mean loss and average rank of 0.25, 0.5, "
            "0.75, 0.9 [Phase 6, if present] and omega = 1 = PSO alone), and omega = 1 has a higher mean loss than omega = 0.75",
     "split.mean_loss, split.avg_rank",
     lambda s: (lambda p: p["mean_loss"]["PSOC"] > p["mean_loss"]["PSOBV75"] and min(p["mean_loss"], key=p["mean_loss"].get) == "PSOBV75"
                and min(p["avg_rank"], key=p["avg_rank"].get) == "PSOBV75")(g(s, "split"))),
    ("C48", "IEA37 CS1, 36 turbines, 30,030 evaluations: the best layout of all methods comes from the gradient-based MS-SLSQP",
     "iea37_compact.methods.*.36T_30030.best",
     lambda s: (lambda m: max((a for a in m if m[a].get("36T_30030") and m[a]["36T_30030"].get("best") is not None),
                              key=lambda a: m[a]["36T_30030"]["best"]) == "SLSQP")(g(s, "iea37_compact", "methods"))),
    # practical equivalence at the single margin equivalence.margin_pp (= EQ_MARGIN of mpce_results.py, set after the
    # primary analysis; never tune it to make these pass)
    ("C49", "equivalence (PSO-VNS vs PSO): the case means are practically equivalent at the margin (90 % bootstrap CI inside "
            "(-m, m)) and the Bayesian signed-rank test gives P(rope) > 0.99 (printed as \\NBayPSOVNSvsPSORope = '>0.99')",
     "equivalence.pairs.PSOBV-PSOC.{equivalent, ci90_mean_dloss_pp, bayes.p_rope}, equivalence.margin_pp",
     lambda s: (lambda r, m: r["equivalent"] and -m < r["ci90_mean_dloss_pp"][0] and r["ci90_mean_dloss_pp"][1] < m
                and r["bayes"]["p_rope"] > 0.99)(g(s, "equivalence", "pairs", "PSOBV-PSOC"), g(s, "equivalence", "margin_pp"))),
    ("C50", "equivalence (LX-SSA-VNS vs RS-VNS): practically equivalent at the margin (90 % bootstrap CI inside (-m, m))",
     "equivalence.pairs.LXBV-RSVNS.{equivalent, ci90_mean_dloss_pp}, equivalence.margin_pp",
     lambda s: (lambda r, m: r["equivalent"] and -m < r["ci90_mean_dloss_pp"][0] and r["ci90_mean_dloss_pp"][1] < m)(
         g(s, "equivalence", "pairs", "LXBV-RSVNS"), g(s, "equivalence", "margin_pp"))),
    ("C51", "equivalence (SSA-VNS vs RS-VNS): NOT equivalent at the margin (90 % bootstrap CI not inside (-m, m)); "
            "SSA-VNS has the lower mean loss and a mean gain larger than m (the CI itself does not exclude gains below m)",
     "equivalence.pairs.SSABV-RSVNS.{equivalent, ci90_mean_dloss_pp, mean_dloss_pp}, equivalence.margin_pp",
     lambda s: (lambda r, m: not r["equivalent"] and not (-m < r["ci90_mean_dloss_pp"][0] and r["ci90_mean_dloss_pp"][1] < m)
                and r["mean_dloss_pp"] < -m)(g(s, "equivalence", "pairs", "SSABV-RSVNS"), g(s, "equivalence", "margin_pp"))),
    ("C52", "spread (PSO-VNS vs PSO, 68 cases, 6,030): no significant difference in the per-case SD of the feasible-run wake "
            "loss nor in the worst feasible run (Wilcoxon p >= 0.05 for both)",
     "spread.sd_p, spread.worst_p",
     lambda s: g(s, "spread", "sd_p") >= 0.05 and g(s, "spread", "worst_p") >= 0.05),
]


# Phase-6 conditions: the disc-sampling control RSD-VNS (experiment rsdisc; component analysis, Holm over the grown
# family of contrasts) and the split omega = 0.9 (experiment omega90). DATA conditions; the sentences follow them.
def _eq(s, k):
    return g(s, "equivalence", "pairs", k)


CONDITIONS += [
    ("C53", "SSA phase vs disc sampling: SSA-VNS is NOT better than RSD-VNS -- never significantly better at run level (W = 0), "
            "higher case mean in more cases than lower (losses > wins), case-mean Wilcoxon p_Holm < 0.05 AGAINST SSA-VNS "
            "(r_rb < 0) and mean dL > 0 (\\NCmSSAVNSvsRSDVNSDL)",
     "ablation.contrasts.SSABV-RSDVNS, ablation.case_mean.SSABV-RSDVNS",
     lambda s: (lambda w, c: w[0] == 0 and c["losses"] > c["wins"] and c["p_holm"] < 0.05 and c["rb"] < 0 and c["mean_dloss_pp"] > 0)(
         _wtl_abl(s, "SSABV-RSDVNS"), _cm(s, "SSABV-RSDVNS"))),
    ("C54", "SSA-VNS vs RSD-VNS, size of the difference: the 95 % bootstrap CI of the mean case-mean difference contains 0, "
            "practical equivalence at +-m is NOT shown (90 % CI not inside (-m, m)), and the Bayesian test gives P(SSA-VNS "
            "better) < 0.01 and P(rope) > 0.5 (\\NBaySSAVNSvsRSDVNSRope)",
     "ablation.case_mean.SSABV-RSDVNS.ci95_mean_dloss_pp, equivalence.pairs.SSABV-RSDVNS",
     lambda s: (lambda c, r, m: c["ci95_mean_dloss_pp"][0] < 0 < c["ci95_mean_dloss_pp"][1] and not r["equivalent"]
                and not (-m < r["ci90_mean_dloss_pp"][0] and r["ci90_mean_dloss_pp"][1] < m)
                and r["bayes"]["p_a_better"] < 0.01 and r["bayes"]["p_rope"] > 0.5)(
         _cm(s, "SSABV-RSDVNS"), _eq(s, "SSABV-RSDVNS"), g(s, "equivalence", "margin_pp"))),
    ("C55", "LX-SSA phase vs disc sampling: LX-SSA-VNS is worse than RSD-VNS -- run level L > W, case-mean p_Holm < 0.05 with "
            "dL > 0 and its 95 % CI above 0, not equivalent at +-m",
     "ablation.contrasts/case_mean LXBV-RSDVNS, equivalence.pairs.LXBV-RSDVNS",
     lambda s: (lambda w, c, r: w[2] > w[0] and c["p_holm"] < 0.05 and c["mean_dloss_pp"] > 0 and c["ci95_mean_dloss_pp"][0] > 0
                and not r["equivalent"])(_wtl_abl(s, "LXBV-RSDVNS"), _cm(s, "LXBV-RSDVNS"), _eq(s, "LXBV-RSDVNS"))),
    ("C56", "PSO phase vs disc sampling: PSO-VNS is better than RSD-VNS -- run level W >= 40 and L = 0, case-mean p_Holm < 0.05 "
            "with dL < 0, and the whole 90 % CI below -m (the gain exceeds the equivalence margin)",
     "ablation.contrasts/case_mean PSOBV-RSDVNS, equivalence.pairs.PSOBV-RSDVNS",
     lambda s: (lambda w, c, r, m: w[0] >= 40 and w[2] == 0 and c["p_holm"] < 0.05 and c["mean_dloss_pp"] < 0
                and r["ci90_mean_dloss_pp"][1] < -m)(
         _wtl_abl(s, "PSOBV-RSDVNS"), _cm(s, "PSOBV-RSDVNS"), _eq(s, "PSOBV-RSDVNS"), g(s, "equivalence", "margin_pp"))),
    ("C57", "disc sampling is the stronger control: RSD-VNS vs RS-VNS never significantly worse at run level (L = 0, W > 0), "
            "case-mean p_Holm < 0.05 with dL < 0; its Phase 1 is feasible at the switch more often and fewer of its runs have "
            "no feasible Phase-1 sample (geometry replay, verified against the stored curves)",
     "ablation.contrasts/case_mean RSDVNS-RSVNS, ablation.phase2_loss_reduction_pct, ablation.phase1_replay",
     lambda s: (lambda w, c, p, rp: w[2] == 0 and w[0] > 0 and c["p_holm"] < 0.05 and c["mean_dloss_pp"] < 0
                and p["RSDVNS"]["feasible_at_switch_pct"] > p["RSVNS"]["feasible_at_switch_pct"]
                and rp["RSDVNS"]["verification_ok"] and rp["RSVNS"]["verification_ok"]
                and rp["RSDVNS"]["runs_no_feasible_sample"] < rp["RSVNS"]["runs_no_feasible_sample"])(
         _wtl_abl(s, "RSDVNS-RSVNS"), _cm(s, "RSDVNS-RSVNS"), g(s, "ablation", "phase2_loss_reduction_pct"),
         g(s, "ablation", "phase1_replay"))),
    ("C58", "component-analysis ranks with RSD-VNS: PSO-VNS first, PSO second, RSD-VNS third, i.e. the disc-sampling control "
            "ranks ahead of SSA-VNS, LX-SSA-VNS and RS-VNS (\\NRankRSDVNS, \\NAblOrder)",
     "ablation.friedman.avg_rank",
     lambda s: (lambda ar: rankpos(ar, "PSOBV") == 1 and rankpos(ar, "PSOC") == 2 and rankpos(ar, "RSDVNS") == 3
                and all(ar["RSDVNS"] < ar[a] for a in ("SSABV", "LXBV", "RSVNS")))(g(s, "ablation", "friedman", "avg_rank"))),
    ("C59", "budget split, omega = 0.9: no significant run-level difference from omega = 0.75 in any case (W = L = 0 from either "
            "side), case-mean p >= 0.05; omega = 0.9 has a higher mean loss and a worse average rank than 0.75 (so 0.75 stays "
            "the best tested setting, see C47) but a lower mean loss than omega = 1 (PSO alone) and than omega = 0.5",
     "split.wtl_90_vs_75, split.case_mean_omega90.PSOBV90-PSOBV75, split.mean_loss, split.avg_rank",
     lambda s: (lambda p: wtl(p["wtl_90_vs_75"])[0] == 0 and wtl(p["wtl_90_vs_75"])[2] == 0
                and p["case_mean_omega90"]["PSOBV90-PSOBV75"]["p"] >= 0.05
                and p["mean_loss"]["PSOBV90"] > p["mean_loss"]["PSOBV75"] and p["avg_rank"]["PSOBV90"] > p["avg_rank"]["PSOBV75"]
                and p["mean_loss"]["PSOBV90"] < p["mean_loss"]["PSOC"] and p["mean_loss"]["PSOBV90"] < p["mean_loss"]["PSOBV"])(g(s, "split"))),
    ("C60", "RS-VNS weakness is mainly spacing, not the square (R3-6): paired over the same runs, fewer of the RS-VNS runs "
            "without a feasible Phase-1 sample are rescued by disc sampling (\\NRsSquareExplains) than stay without one "
            "under disc sampling too (\\NRsSpacingDominates); the counts add up to the RS-VNS total and the RSD-VNS total "
            "(\\NRSRunsNoFeasSample, \\NRSDRunsNoFeasSample); both replays reproduce the stored curves; and every RS-VNS run "
            "without a feasible sample does draw samples wholly inside the circle (\\NRsNoFeasInside = \\NRSRunsNoFeasSample), "
            "so spacing, not the boundary, is binding in all of them",
     "ablation.phase1_replay.square_vs_disc, ablation.phase1_replay.RSVNS/RSDVNS (incl. RSVNS.runs_no_feasible_but_inside_sample)",
     lambda s: (lambda rp: (lambda d: d["rs_none_rsd_some"] < d["rs_none_rsd_none"]
                            and d["rs_none_rsd_some"] + d["rs_none_rsd_none"] == rp["RSVNS"]["runs_no_feasible_sample"]
                            and d["rs_none_rsd_none"] + d["rs_some_rsd_none"] == rp["RSDVNS"]["runs_no_feasible_sample"]
                            and rp["RSVNS"]["verification_ok"] and rp["RSDVNS"]["verification_ok"]
                            and rp["RSVNS"]["runs_no_feasible_but_inside_sample"] == rp["RSVNS"]["runs_no_feasible_sample"])(
         rp["square_vs_disc"]))(g(s, "ablation", "phase1_replay"))),
    ("C61", "switch point (R3-9): RS-VNS and RSD-VNS switch after round(0.5 x 6,030) = 3,015 Phase-1 evaluations, the swarm "
            "hybrids after 3,030; the RS-VNS / RSD-VNS runs infeasible at the switch (read at the last checkpoint before it, "
            "call 3,000) equal the replayed runs without a feasible sample among the first 3,000 samples, and for RS-VNS also "
            "among all 3,015 (\\NRSRunsNoFeasSample; RSD-VNS: at most 15 samples later, a few runs more are feasible)",
     "ablation.phase2_loss_reduction_pct.*.switch_call/feasible_at_switch_pct/runs, ablation.phase1_replay",
     lambda s: (lambda p, rp: p["RSVNS"]["switch_call"] == 3015 and p["RSDVNS"]["switch_call"] == 3015
                and all(p[a]["switch_call"] == 3030 for a in ("PSOBV", "SSABV", "LXBV"))
                and all(round(p[a]["runs"] * (1 - p[a]["feasible_at_switch_pct"] / 100)) == rp[a]["runs_no_feasible_first3000"]
                        >= rp[a]["runs_no_feasible_sample"] for a in ("RSVNS", "RSDVNS"))
                and rp["RSVNS"]["runs_no_feasible_first3000"] == rp["RSVNS"]["runs_no_feasible_sample"])(
         g(s, "ablation", "phase2_loss_reduction_pct"), g(s, "ablation", "phase1_replay"))),
]


# Text conditions: evaluated on the manuscript text (comments stripped), not on the summary. D16/D19: with 1-deg
# bins (F21) and under PyWake no optimized Horns Rev run exceeds the installed block, so the paper must not claim
# runs above the installed layout: none of the \NHR...Above... macros (runs_above_installed) may be used.
ABOVE_MACRO = re.compile(r"\\NHR[A-Za-z]*Above[A-Za-z]*")


def _uncommented(line):
    """Line without its LaTeX comment (first % not preceded by a backslash)."""
    m = re.search(r"(?<!\\)%", line)
    return line[:m.start()] if m else line


def _above_uses(T):
    return sorted({f"{fn}: {m}" for fn, lines in T.items() for ln in lines
                   if not re.search(r"\\(?:provide|new|renew)command", ln) for m in ABOVE_MACRO.findall(ln)})


TEXT_CONDITIONS = [
    ("C62", "Horns Rev (D16/D19): no claim of optimized runs above the installed layout -- no \\NHR...Above... macro "
            "(runs_above_installed) is used in the text of the manuscript, optA/*.tex or the supplement",
     "manuscript text (comments stripped)", lambda T: not _above_uses(T)),
]


# ------------------------------------------------------------------ reference scan of the CHECK- comments
KIND_LETTER = {"FINAL": "C", "EXTRA": "X", "DIAG": "D", "THEORY": "T", "DIR": "F", "CSWEEP": "S"}
KIND_SOURCE = {"EXTRA": ("mpce_check_extra.py", r"\(\s*\"(X\d+)\""), "DIAG": ("mpce_check_diag.py", r"\(\s*\"(D\d+)\""),
               "THEORY": ("make_theory_figures.py", r"(?:has|append)\(\s*\(?\s*\"(T\d+)\""),
               "DIR": ("mpce_check_dir.py", r"\(\s*\"(F\d+)\""), "CSWEEP": ("mpce_check_csweep.py", r"\(\s*\"(S\d+)\"")}
MALFORMED = []


def scan_refs(files):
    """{kind: {id: set(files)}} of the ids referenced by the % CHECK-<kind> comments of the given .tex files.
    Malformed tags (unknown kind, or no id of the kind's letter such as "CHECK-FINAL [NEW, D2]") go to MALFORMED."""
    refs = {k: {} for k in KIND_LETTER}
    for fn in files:
        for no, line in enumerate(open(fn, encoding="utf-8", errors="replace"), 1):
            if "CHECK-" not in line or "%" not in line:
                continue
            com = line[line.index("%"):]
            parts = re.split(r"CHECK-([A-Za-z]*)", com)
            for kind, scope in zip(parts[1::2], parts[2::2]):
                where = f"{os.path.basename(fn)}:{no}: CHECK-{kind}{scope.rstrip()[:60]}"
                if kind not in KIND_LETTER:
                    MALFORMED.append(f"{where}  (unknown kind 'CHECK-{kind}'; known: "
                                     + ", ".join("CHECK-" + k for k in KIND_LETTER) + ")")
                    continue
                L = KIND_LETTER[kind]
                if re.match(r"\s*:?\s*\[[A-Z]nn\]", scope):          # documentation of the convention ("[Cnn]")
                    continue
                ids = set()
                for m in re.finditer(r"\[(%s)(\d{2,3})[^\]]*\]\s*(?:--|-|–|—)\s*\[%s(\d{2,3})" % (L, L), scope):
                    ids |= {f"{L}{i:02d}" for i in range(int(m.group(2)), int(m.group(3)) + 1)}
                for br in re.findall(r"\[([^\]]*)\]", scope):
                    ids |= set(re.findall(r"\b(%s\d{2,3})\b" % L, br))
                for i in ids:
                    refs[kind].setdefault(i, set()).add(os.path.basename(fn))
                if not ids:                                           # R5-12: reported and counted as an error
                    MALFORMED.append(f"{where}  (no [{L}nn] id)")
    return refs


def defined_ids(kind):
    fn, pat = KIND_SOURCE[kind]
    p = os.path.join(HERE, fn)
    return set(re.findall(pat, open(p).read())) if os.path.exists(p) else None


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--summary", default=os.path.join(HERE, "mpce_summary.json"))
    ap.add_argument("--tex", default=os.path.join(HERE, "..", "MPCE_PSO_VNS.tex"))
    a = ap.parse_args(argv)
    s = json.load(open(a.summary))
    del MALFORMED[:]
    texs, T = [], None
    if os.path.exists(a.tex):
        import glob as _g
        texs = [a.tex] + sorted(_g.glob(os.path.join(os.path.dirname(a.tex), "optA", "*.tex")))
        sup = os.path.join(os.path.dirname(a.tex), "MPCE_PSO_VNS_supplement.tex")
        texs += [sup] if os.path.exists(sup) else []
        T = {os.path.basename(f): [_uncommented(l) for l in open(f, encoding="utf-8", errors="replace")] for f in texs}
    res = {}
    print(f"CHECK-FINAL conditions (summary generated {s.get('generated', '?')})")
    for cid, text, keys, fn in CONDITIONS:
        try:
            r = "PASS" if fn(s) else "FAIL"
        except (KeyError, TypeError, ValueError, IndexError, AttributeError, ZeroDivisionError):
            r = "PENDING"
        res[cid] = r
        print(f"  [{cid}] {r:7s} {text}\n            keys: {keys}")
    for cid, text, keys, fn in TEXT_CONDITIONS:
        r = "PENDING" if T is None else ("PASS" if fn(T) else "FAIL")
        res[cid] = r
        print(f"  [{cid}] {r:7s} {text}\n            keys: {keys}")
        if cid == "C62" and r == "FAIL":
            for u in _above_uses(T):
                print(f"            used: {u}")
    errors = 0
    refs = None
    if texs:
        refs = scan_refs(texs)
        used = set(refs["FINAL"])
        for cid in sorted(used - set(res)):
            print(f"  [{cid}] UNDEFINED: used in the manuscript but not defined in mpce_check_final.py")
            errors += 1
        for cid in sorted(set(res) - used - {c[0] for c in TEXT_CONDITIONS}):
            print(f"  [{cid}] (not referenced in the manuscript)")
        for m in MALFORMED:
            print(f"  ERROR: malformed CHECK tag (not evaluated): {m}")
        errors += len(MALFORMED)
    n = {k: sum(v == k for v in res.values()) for k in ("PASS", "FAIL", "PENDING")}
    if refs is not None:
        print("  referenced checks of the other scripts (not evaluated here; run the script named):")
        for kind in ("EXTRA", "DIAG", "DIR", "CSWEEP", "THEORY"):
            ids = sorted(refs[kind])
            dfn = defined_ids(kind)
            und = [i for i in ids if dfn is not None and i not in dfn]
            errors += len(und)
            print(f"    CHECK-{kind:6s} ({KIND_SOURCE[kind][0]}): {len(ids)} ids referenced: {', '.join(ids) if ids else 'none'}"
                  + (f"; UNDEFINED in {KIND_SOURCE[kind][0]}: {', '.join(und)}" if und else "")
                  + (f"; defined but not referenced: {', '.join(sorted(dfn - set(ids)))}" if dfn else "")
                  + ("" if dfn is not None else f"; {KIND_SOURCE[kind][0]} not found"))
    fb = [f for f in s.get("fallbacks", []) if "->" in f]
    for f in fb:
        print(f"  NOTE (data): {f}")
    print(f"  summary: {n['PASS']} PASS, {n['FAIL']} FAIL, {n['PENDING']} PENDING, {errors} TAG ERRORS "
          f"({len(MALFORMED)} malformed, {errors - len(MALFORMED)} undefined ids)")
    return 1 if n["FAIL"] or errors else 0


if __name__ == "__main__":
    sys.exit(main())

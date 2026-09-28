"""Checks X01... of the inference-robustness analyses (Phase 6, W2) against mpce_summary_extra.json.

Usage:  python3 analysis/mpce_check_extra.py [--summary analysis/mpce_summary_extra.json]

Each check is ONE sentence that the paper (main text or supplement) may state; the condition is evaluated on the
output of mpce_inference_extra.py and printed as PASS / FAIL. A FAIL means the sentence must not be written (or
must be rewritten). Exit code 1 if any check fails. The sentences use the \\NX... macros of mpce_numbers_extra.tex
for every number.
"""
import os, re, sys, json, argparse
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
SIG_PAIRS_EXPECTED_NS = {"PSOBV-PSOC", "LXBV-RSVNS"}          # the two pairs the paper calls not significant
RSD = ["SSABV-RSDVNS", "LXBV-RSDVNS", "PSOBV-RSDVNS", "RSDVNS-RSVNS"]  # Phase 6 disc-sampling control contrasts
EQ_MARGIN = 0.05                     # overwritten in main() by mpce_summary_extra.json["equivalence_levels"]["margin_pp"]


def thr(s):
    return [s["threshold"][str(t)] for t in s["thresholds"]]


def case_sig(s):
    return {k for k, v in s["case_level"].items() if v["p"] < 0.05}


def mt(s, key):
    return next(t for t in s["multiplicity"]["tests"] if t["key"] == key)


# ---- review round 2 (D13) helpers
def eq(s, k):
    return s["equivalence_levels"]["pairs"][k]


def m_(s):
    return s["equivalence_levels"]["margin_pp"]


def inside(ci, m):
    return -m < ci[0] and ci[1] < m


def het(s):
    return s["heterogeneity"]


def fig_ok(s):
    f = os.path.join(os.path.dirname(HERE), s["figure_equiv_curve"])
    return os.path.exists(f) and abs(os.path.getmtime(f) - os.path.getmtime(os.path.join(HERE, "mpce_summary_extra.json"))) < 600


CHECKS = [
    ("X01", "The re-implementation reproduces the pipeline exactly at the paper's 15-run threshold (average ranks of the main "
            "and component-analysis Friedman rankings, and every case-mean p of the main table, of the fifteen component-analysis "
            "contrasts incl. RSD-VNS (raw and Holm), of the N>=10 subgroups and of the six budget-split tests incl. omega = 0.9 "
            "equal mpce_summary.json).",
     lambda s: s["reproduction"]["ok"] is True),
    ("X02", "The best-ranked method of the eight-method comparison is PSO-VNS for every qualification threshold from 1 to 30 "
            "feasible runs (\\NXThrList).",
     lambda s: all(x["friedman"]["best_ranked"] == "PSOBV" for x in thr(s))),
    ("X03", "The complete average-rank order of the eight methods is the same at every threshold 1-30.",
     lambda s: s["threshold_summary"]["order_identical"]),
    ("X04", "At every threshold, PSO is the only method whose average rank does not differ significantly from that of PSO-VNS "
            "(Holm post hoc).",
     lambda s: all(x["friedman"]["posthoc_nonsig"] == ["PSOC"] for x in thr(s))),
    ("X05", "At every threshold, PSO-VNS and PSO do not differ significantly on the case means (p >= 0.05, p between "
            "\\NXThrPSOPMin and \\NXThrPSOPMax).",
     lambda s: all(x["main_case_mean"]["PSOC"]["p"] >= 0.05 and x["main_case_mean"]["PSOC"]["p_holm"] >= 0.05 for x in thr(s))),
    ("X06", "At every threshold, SSA-VNS has a significantly lower case-mean wake loss than RS-VNS (unadjusted and Holm over "
            "the component-analysis contrasts; p <= \\NXThrSSAVNSvsRSVNSPMax).",
     lambda s: all(x["conclusions"]["SSABV-RSVNS_A"] for x in thr(s))),
    ("X07", "At every threshold, LX-SSA-VNS and RS-VNS do not differ significantly (p >= \\NXThrLXSSAVNSvsRSVNSPMin).",
     lambda s: all(x["conclusions"]["LXBV-RSVNS_ns"] for x in thr(s))),
    ("X08", "At every threshold, SSA-VNS has a significantly lower case-mean wake loss than LX-SSA-VNS.",
     lambda s: all(x["conclusions"]["SSABV-LXBV_A"] for x in thr(s))),
    ("X09", "At every threshold, PSO-VNS has a significantly lower case-mean wake loss than RS-VNS.",
     lambda s: all(x["conclusions"]["PSOBV-RSVNS_A"] for x in thr(s))),
    ("X10", "The 68 cases form six (data set, radius) clusters with nested N (9, 11 and 14 cases per radius in each data set).",
     lambda s: s["n_cases"] == 68 and len(s["clusters"]) == 6
     and sorted(s["cluster_sizes"].values()) == [9, 9, 11, 11, 14, 14]),
    ("X11", "Every comparison that is significant on the 68 case means favours the same method in all six clusters (the "
            "strongest result attainable with six clusters, exact p = \\NXClMinP), except the two contrasts of the SSA hybrids "
            "with the disc-sampling control (\\NXClCaseSigNotUnanimous: \\NXClSSAVNSvsRSDVNSFav/\\NXClSSAVNSvsRSDVNSOpp and "
            "\\NXClLXSSAVNSvsRSDVNSFav/\\NXClLXSSAVNSvsRSDVNSOpp clusters).",
     lambda s: all(min(s["cluster"][k]["clusters_favour_first"], s["cluster"][k]["clusters_favour_second"]) == 0
                   and s["cluster"][k]["wilcoxon_exact_p"] <= s["cluster_min_attainable_p"] + 1e-12
                   for k in case_sig(s) - {"SSABV-RSDVNS", "LXBV-RSDVNS"})
     and all(min(s["cluster"][k]["clusters_favour_first"], s["cluster"][k]["clusters_favour_second"]) > 0
             for k in ("SSABV-RSDVNS", "LXBV-RSDVNS"))
     and case_sig(s) == set(s["case_level"]) - SIG_PAIRS_EXPECTED_NS),
    ("X12", "At the cluster level PSO-VNS and PSO do not differ either: PSO-VNS has the lower mean in \\NXClPSOVNSvsPSOFav of the "
            "six clusters (all cluster-level p >= 0.05) and the cluster-bootstrap CI \\NXClPSOVNSvsPSOCI contains zero.",
     lambda s: (lambda c: c["clusters_favour_first"] > 0 and c["clusters_favour_second"] > 0
                and min(c["sign_test_p"], c["wilcoxon_exact_p"], c["cluster_signflip_p"], c["cr1_p"]) >= 0.05
                and not c["cluster_boot_excludes_zero"])(s["cluster"]["PSOBV-PSOC"])),
    ("X13", "LX-SSA-VNS and RS-VNS split the clusters (\\NXClLXSSAVNSvsRSVNSFav/\\NXClLXSSAVNSvsRSVNSOpp) and do not differ at the "
            "cluster level.",
     lambda s: (lambda c: c["clusters_favour_first"] > 0 and c["clusters_favour_second"] > 0
                and min(c["sign_test_p"], c["wilcoxon_exact_p"], c["cluster_signflip_p"], c["cr1_p"]) >= 0.05)(s["cluster"]["LXBV-RSVNS"])),
    ("X14", "Leaving out any one of the six clusters changes none of the \\NXLocoNPairs case-mean verdicts (significant or not, "
            "and which method is better when significant) except that of SSA-VNS vs RSD-VNS (\\NXLocoSSAVNSvsRSDVNSChanges of six "
            "runs, p up to \\NXLocoSSAVNSvsRSDVNSPMax); PSO-VNS keeps the best average rank, and PSO remains the only method not "
            "significantly behind it.",
     lambda s: {k for k, v in s["loco"]["pairs"].items() if v["verdict_changes"] > 0} == {"SSABV-RSDVNS"}
     and s["loco"]["best_ranked_always_PSOBV"] and s["loco"]["posthoc_only_PSOC_always"]),
    ("X15", "In the leave-one-cluster-out runs the sign of the mean difference changes only for the two non-significant pairs "
            "(PSO-VNS vs PSO, LX-SSA-VNS vs RS-VNS) and for SSA-VNS vs RSD-VNS; PSO-VNS vs PSO stays non-significant "
            "(p >= \\NXLocoPSOVNSvsPSOPMin).",
     lambda s: {k for k, v in s["loco"]["pairs"].items() if v["direction_flips"] > 0} <= SIG_PAIRS_EXPECTED_NS | {"SSABV-RSDVNS"}
     and s["loco"]["pairs"]["PSOBV-PSOC"]["p_min"] >= 0.05),
    ("X16", "The cluster-bootstrap 95% CI excludes zero for every comparison that is significant on the case means except "
            "SSA-VNS vs RSD-VNS (\\NXClSSAVNSvsRSDVNSCI), and contains zero for the two that are not significant.",
     lambda s: all(s["cluster"][k]["cluster_boot_excludes_zero"] == (k in case_sig(s) - {"SSABV-RSDVNS"}) for k in s["cluster"])),
    ("X17", "CAVEAT: three case-level-significant comparisons are not significant with the cluster-robust t test (5 d.f.): "
            "SSA-VNS vs RS-VNS (p = \\NXClSSAVNSvsRSVNSCRP, although all six clusters favour SSA-VNS; consistent but small), "
            "SSA-VNS vs RSD-VNS (p = \\NXClSSAVNSvsRSDVNSCRP) and LX-SSA-VNS vs RSD-VNS (p = \\NXClLXSSAVNSvsRSDVNSCRP); all other "
            "significant comparisons are.",
     lambda s: {k for k in case_sig(s) if s["cluster"][k]["cr1_p"] >= 0.05} == {"SSABV-RSVNS", "SSABV-RSDVNS", "LXBV-RSDVNS"}
     and s["cluster"]["SSABV-RSVNS"]["clusters_favour_first"] == 6),
    ("X18", "With one Holm correction over all \\NXMultN case-mean tests quoted in the main text (\\NXMultNMain main-table, "
            "\\NXMultNAbl component-analysis, \\NXMultNSub subgroup and \\NXMultNSplit budget-split tests), every test significant "
            "unadjusted stays significant except the post hoc N>=10 subgroup over both data sets and omega = 0.5 vs omega = 0.9 "
            "(\\NXMultLostHolm); under Benjamini-Hochberg only the former is lost (\\NXMultLostBH).",
     lambda s: sorted(s["multiplicity"]["lost_under_holm"]) == ["split:PSOBV-PSOBV90", "sub:Large"]
     and s["multiplicity"]["lost_under_bh"] == ["sub:Large"]),
    ("X19", "The Data Set II, N>=10 subgroup of PSO-VNS vs PSO remains significant after the all-family Holm correction "
            "(p_Holm <= \\NXHolmPSOVNSvsPSOdsIILargeP).",
     lambda s: mt(s, "sub:dsIILarge")["sig_holm"] and mt(s, "sub:dsIILarge")["wins"] > mt(s, "sub:dsIILarge")["losses"]),
    ("X20", "SSA-VNS vs RS-VNS and SSA-VNS vs LX-SSA-VNS remain significant after the all-family Holm correction; PSO-VNS vs PSO "
            "and LX-SSA-VNS vs RS-VNS are not significant with or without it.",
     lambda s: mt(s, "abl:SSABV-RSVNS")["sig_holm"] and mt(s, "abl:SSABV-LXBV")["sig_holm"]
     and not mt(s, "main:PSOBV-PSOC")["sig_raw"] and not mt(s, "abl:LXBV-RSVNS")["sig_raw"]),
    ("X21", "The budget-split case-mean tests omega = 0.5 vs 0.25, 0.5 vs 0.75 and 0.75 vs 1 (plain PSO) remain significant after "
            "the all-family Holm correction; omega = 0.5 vs 1 is not significant.",
     lambda s: all(mt(s, k)["sig_holm"] for k in ("split:PSOBV-PSOBV25", "split:PSOBV-PSOBV75", "split:PSOBV75-PSOC"))
     and not mt(s, "split:PSOBV-PSOC")["sig_raw"]),
    ("X22", "In Data Set II the PSO-VNS advantage over PSO is significantly larger for N>=10 than for N<10 (permutation within "
            "clusters, p = \\NXIntDsIIP); in Data Set I there is no such difference (p = \\NXIntDsIP, opposite sign).",
     lambda s: s["subgroup"]["interaction_dsII"]["diff_pp"] < 0 and s["subgroup"]["interaction_dsII"]["p_stratified"] < 0.05
     and s["subgroup"]["interaction_dsI"]["diff_pp"] > 0 and s["subgroup"]["interaction_dsI"]["p_stratified"] >= 0.05),
    ("X23", "Without any threshold, the PSO-VNS advantage grows with N within every Data Set II farm (mean within-farm Spearman "
            "\\NXSpearDsIIStrat, p = \\NXSpearDsIIStratP) but not over all cases (p = \\NXSpearAllStratP), because in Data Set I "
            "the mean within-farm trend has the opposite sign (\\NXSpearDsIStrat, p = \\NXSpearDsIStratP); 'the gain grows with N' must therefore not be written without "
            "'in Data Set II'.",
     lambda s: s["subgroup"]["spearman_dsII"]["stratified_mean_rho"] < 0 and s["subgroup"]["spearman_dsII"]["stratified_p"] < 0.05
     and all(v < 0 for v in s["subgroup"]["spearman_dsII"]["rho_by_cluster"].values())
     and s["subgroup"]["spearman_all"]["stratified_p"] >= 0.05 and s["subgroup"]["spearman_dsI"]["stratified_mean_rho"] > 0),
    ("X24", "For the same farm and N>=10, the PSO-VNS advantage over PSO is larger in Data Set II than in Data Set I in "
            "\\NXDsPairLargeIIWins of \\NXDsPairLargeN cases (Wilcoxon p = \\NXDsPairLargeP), but not for N<10 (p = \\NXDsPairSmallP).",
     lambda s: s["subgroup"]["dsII_vs_dsI_large"]["p"] < 0.05 and s["subgroup"]["dsII_vs_dsI_large"]["mean_pp"] < 0
     and s["subgroup"]["dsII_vs_dsI_small"]["p"] >= 0.05),
    ("X25", "Over the 68 cases the benchmark-energy gain of PSO-VNS over PSO is below 0.05% of AEP (\\NXEnPSOVNSvsPSOPct%, "
            "\\NXEnPSOVNSvsPSOMWh MWh/yr per case) and its 95% CI includes zero.",
     lambda s: (lambda e: abs(e["mean_gain_pct_aep"]) < 0.05 and e["ci95_gain_pct_aep"][0] < 0 < e["ci95_gain_pct_aep"][1])(s["energy"]["PSOBV-PSOC"])),
    ("X26", "The benchmark-energy gain of SSA-VNS over RS-VNS is positive (CI excludes zero) but below 0.1% of AEP and less than a "
            "fifth of that of PSO-VNS over RS-VNS.",
     lambda s: (lambda a, b: 0 < a["ci95_gain_pct_aep"][0] and a["mean_gain_pct_aep"] < 0.1
                and a["mean_gain_pct_aep"] < 0.2 * b["mean_gain_pct_aep"])(s["energy"]["SSABV-RSVNS"], s["energy"]["PSOBV-RSVNS"])),
    ("X27", "In Data Set II with N>=10 the gain of PSO-VNS over PSO is about 0.2% of benchmark AEP (\\NXEnPSOVNSvsPSOdsIILargePct%, "
            "CI \\NXEnPSOVNSvsPSOdsIILargePctCI) and positive in every case.",
     lambda s: (lambda e: 0.1 <= e["mean_gain_pct_aep"] < 0.3 and e["ci95_gain_pct_aep"][0] > 0
                and e["min_gain_pct_aep"] > 0)(s["energy"]["PSOBV-PSOC_dsIILarge"])),
    ("X29", "The wake-model uncertainty is several times the equivalence margin: the wake loss of the same final layouts "
            "differs between the Jensen and the Gaussian wake model by \\NXModelShiftMean pp on average "
            "(\\NXModelShiftRatio times the margin; > 5 x 0.05 pp). [Replaces the earlier '> 10 x' version, which the data do not support.]",
     lambda s: s["model_shift"]["mean_abs_shift_pp"] > 5 * EQ_MARGIN),
    ("X30", "The wake loss of the same final layouts differs between the Jensen and Gaussian wake models by \\NXModelShiftMean pp "
            "on average (median \\NXModelShiftMedian pp; \\NXModelShiftRatio times the equivalence margin; above the margin for "
            "\\NXModelShiftAboveMarginPct% of the \\NXModelShiftNLayouts layouts), whereas the per-case PSO-VNS - PSO difference "
            "changes by only \\NXModelShiftPairMean pp on average between the models (below the margin; mean \\NXModelShiftPairJensen "
            "vs \\NXModelShiftPairGauss pp).",
     lambda s: (lambda m, p: m["mean_abs_shift_pp"] > 5 * EQ_MARGIN and m["median_abs_shift_pp"] > 5 * EQ_MARGIN
                and m["share_abs_shift_above_margin"] > 0.5 and p["mean_abs_change_pp"] < EQ_MARGIN
                and abs(p["mean_diff_jensen_pp"]) < EQ_MARGIN and abs(p["mean_diff_gauss_pp"]) < EQ_MARGIN
                and s["model_shift"]["pair_recorded_vs_paper_abs_diff"] < 1e-9)(s["model_shift"], s["model_shift"]["pairs"]["PSOBV-PSOC"])),
    ("X31", "At every qualification threshold 1-30 the verdicts of the four RSD-VNS contrasts (SSA-VNS and LX-SSA-VNS worse, "
            "PSO-VNS better than RSD-VNS; RSD-VNS better than RS-VNS; unadjusted and Holm over the component-analysis contrasts) "
            "and of the six budget-split tests (incl. omega = 0.9 better than 0.5 unadjusted, 0.9 vs 0.75 not significant) are unchanged "
            "(SSA-VNS vs RSD-VNS p between \\NXThrSSAVNSvsRSDVNSPMin and \\NXThrSSAVNSvsRSDVNSPMax).",
     lambda s: all(v for x in thr(s) for c, v in x["conclusions"].items() if "RSDVNS" in c or c.startswith("split:"))
     and sum(1 for c in thr(s)[0]["conclusions"] if "RSDVNS" in c or c.startswith("split:")) == 10),
    ("X32", "PSO-VNS vs RSD-VNS and RSD-VNS vs RS-VNS are robust to case dependence: all six clusters favour PSO-VNS and RSD-VNS "
            "respectively (exact p = \\NXClMinP), the cluster-robust t tests are significant (p = \\NXClPSOVNSvsRSDVNSCRP, "
            "\\NXClRSDVNSvsRSVNSCRP), the cluster-bootstrap CIs exclude zero (\\NXClPSOVNSvsRSDVNSCI, \\NXClRSDVNSvsRSVNSCI) and no "
            "leave-one-cluster-out run changes the verdict.",
     lambda s: all(s["cluster"][k]["clusters_favour_first"] == 6 and s["cluster"][k]["cr1_p"] < 0.05
                   and s["cluster"][k]["cluster_boot_excludes_zero"] and s["loco"]["pairs"][k]["verdict_changes"] == 0
                   and s["case_level"][k]["verdict"] == "A" for k in ("PSOBV-RSDVNS", "RSDVNS-RSVNS"))),
    ("X33", "CAVEAT: the case-level result that SSA-VNS is worse than RSD-VNS (Wilcoxon, also after the all-family Holm, p_Holm <= "
            "\\NXHolmSSAVNSvsRSDVNSP) does NOT survive the cluster analyses: RSD-VNS has the lower mean in only \\NXClSSAVNSvsRSDVNSOpp of "
            "six clusters (all cluster-level p >= 0.05, cluster-robust p = \\NXClSSAVNSvsRSDVNSCRP), the cluster-bootstrap CI "
            "\\NXClSSAVNSvsRSDVNSCI contains zero and one leave-one-cluster-out run loses significance (p up to "
            "\\NXLocoSSAVNSvsRSDVNSPMax). Write 'the SSA phase gives no gain over disc sampling', not 'SSA-VNS is significantly worse'.",
     lambda s: (lambda c, t: s["case_level"]["SSABV-RSDVNS"]["verdict"] == "B" and t["sig_holm"]
                and c["clusters_favour_first"] > 0 and min(c["sign_test_p"], c["wilcoxon_exact_p"], c["cluster_signflip_p"], c["cr1_p"]) >= 0.05
                and not c["cluster_boot_excludes_zero"] and s["loco"]["pairs"]["SSABV-RSDVNS"]["verdict_changes"] >= 1)(
         s["cluster"]["SSABV-RSDVNS"], mt(s, "abl:SSABV-RSDVNS"))),
    ("X34", "LX-SSA-VNS is worse than RSD-VNS on the case means (all-family Holm p <= \\NXHolmLXSSAVNSvsRSDVNSP) in \\NXClLXSSAVNSvsRSDVNSOpp "
            "of six clusters, with a cluster-bootstrap CI above zero (\\NXClLXSSAVNSvsRSDVNSCI) and no leave-one-cluster-out change; but "
            "the cluster-robust t test is not significant (p = \\NXClLXSSAVNSvsRSDVNSCRP) -> 'worse', not 'clearly worse at the farm level'.",
     lambda s: (lambda c: s["case_level"]["LXBV-RSDVNS"]["verdict"] == "B" and mt(s, "abl:LXBV-RSDVNS")["sig_holm"]
                and c["clusters_favour_second"] >= 5 and c["cluster_boot_ci95"][0] > 0 and c["cr1_p"] >= 0.05
                and s["loco"]["pairs"]["LXBV-RSDVNS"]["verdict_changes"] == 0)(s["cluster"]["LXBV-RSDVNS"])),
    ("X35", "All four RSD-VNS contrasts remain significant after the all-family Holm correction over \\NXMultN tests "
            "(SSA-VNS vs RSD-VNS p_Holm <= \\NXHolmSSAVNSvsRSDVNSP, LX-SSA-VNS vs RSD-VNS \\NXHolmLXSSAVNSvsRSDVNSP, PSO-VNS vs RSD-VNS "
            "\\NXHolmPSOVNSvsRSDVNSP, RSD-VNS vs RS-VNS \\NXHolmRSDVNSvsRSVNSP).",
     lambda s: all(mt(s, f"abl:{k}")["sig_holm"] for k in RSD)),
    ("X36", "omega = 0.9 does not beat omega = 0.75 at any level (case means p = \\NXThrSplitNinetyVsSeventyFivePMin-"
            "\\NXThrSplitNinetyVsSeventyFivePMax over thresholds, all-family Holm p = \\NXHolmSplitNinetyVsSeventyFiveP, "
            "\\NXClSplitNinetyVsSeventyFiveFav/\\NXClSplitNinetyVsSeventyFiveOpp clusters, cluster-bootstrap CI "
            "\\NXClSplitNinetyVsSeventyFiveCI); omega = 0.9 is better than 0.5 unadjusted and in all six clusters, but not after the "
            "all-family Holm (p = \\NXHolmSplitFiftyVsNinetyP) and not in \\NXLocoSplitFiftyVsNinetyChanges of six LOCO runs.",
     lambda s: (lambda c9, c5: s["split_case_level"]["PSOBV90-PSOBV75"]["p"] >= 0.05 and not mt(s, "split:PSOBV90-PSOBV75")["sig_holm"]
                and not c9["cluster_boot_excludes_zero"] and c9["cr1_p"] >= 0.05
                and s["split_case_level"]["PSOBV-PSOBV90"]["verdict"] == "B" and c5["clusters_favour_second"] == 6
                and not mt(s, "split:PSOBV-PSOBV90")["sig_holm"] and s["loco_split"]["pairs"]["PSOBV-PSOBV90"]["verdict_changes"] >= 1)(
         s["cluster_split"]["PSOBV90-PSOBV75"], s["cluster_split"]["PSOBV-PSOBV90"])),
    ("X37", "The budget-split conclusions omega = 0.5 better than 0.25, 0.75 better than 0.5 and 0.75 better than 1 (PSO) hold in "
            "all six clusters (cluster-bootstrap CIs exclude zero, no LOCO change); 0.5 vs 1 splits the clusters.",
     lambda s: all(min(s["cluster_split"][k]["clusters_favour_first"], s["cluster_split"][k]["clusters_favour_second"]) == 0
                   and s["cluster_split"][k]["cluster_boot_excludes_zero"] and s["loco_split"]["pairs"][k]["verdict_changes"] == 0
                   for k in ("PSOBV-PSOBV25", "PSOBV-PSOBV75", "PSOBV75-PSOC"))
     and min(s["cluster_split"]["PSOBV-PSOC"]["clusters_favour_first"], s["cluster_split"]["PSOBV-PSOC"]["clusters_favour_second"]) > 0),
    # ---------------- review round 2 (R2 statistics, lead decision D13)
    ("X38", "The case-level equivalence results (mean, 90%/95% CIs, bootstrap TOST p, t-TOST p, minimal margin, verdict) and "
            "the Bayesian signed-rank probabilities and posterior means of all 18 pairs of tab:equivalence are reproduced "
            "exactly from the per-run data (margin \\NXEqMarginSource).",
     lambda s: s["equivalence_levels"]["reproduction"]["ok"] is True and s["equivalence_levels"]["margin_pp"] == s["equivalence_levels"]["reproduction"]["margin_pipeline"]),
    ("X39", "Equivalence of PSO-VNS and PSO at +-margin holds at the seed level (fixed benchmark, 90% CI \\NXEqSeedPSOVNSvsPSOCI) "
            "and at the case level (\\NXEqCasePSOVNSvsPSOCI, p_TOST = \\NXEqCasePSOVNSvsPSOP) but NOT at the cluster level "
            "(CR2 90% CI \\NXEqClustPSOVNSvsPSOCI, wild-cluster TOST p = \\NXEqClustPSOVNSvsPSOP).",
     lambda s: eq(s, "PSOBV-PSOC")["seed"]["equivalent"] and eq(s, "PSOBV-PSOC")["case"]["equivalent"]
     and not eq(s, "PSOBV-PSOC")["cluster"]["cr2"]["equivalent"] and not eq(s, "PSOBV-PSOC")["cluster"]["wild"]["equivalent"]
     and eq(s, "PSOBV-PSOC")["cluster"]["wild"]["p_tost"] > 0.05),
    ("X40", "On the fixed benchmark PSO-VNS has a small but statistically significant average advantage over PSO that lies "
            "inside the margin: seed-level 90% CI \\NXEqSeedPSOVNSvsPSOCI (95%: \\NXEqSeedPSOVNSvsPSOCINinetyFive; one seed "
            "resample for all cases: \\NXEqSeedStratJointPSOVNSvsPSOCI; benchmark average per seed, t(29): "
            "\\NXEqSeedJointPSOVNSvsPSOCI).",
     lambda s: (lambda r: all(-m_(s) < c[0] and c[1] < 0 for c in (r["ci90"], r["ci95"], r["joint"]["ci90"], r["per_seed_t"]["ci90"]))
                and r["fallback_resamples"] == 0)(eq(s, "PSOBV-PSOC")["seed"])),
    ("X41", "The minimal equivalence margin of PSO-VNS vs PSO grows with the level of generalization: \\NXEqSeedPSOVNSvsPSOMin pp "
            "(seed) < \\NXEqCasePSOVNSvsPSOMin pp (case) < margin < \\NXEqClustPSOVNSvsPSOMin pp (cluster, CR2; wild bootstrap "
            "\\NXEqClustPSOVNSvsPSOMinWild pp).",
     lambda s: (lambda r: r["seed"]["min_margin_pp"] < r["case"]["min_margin_pp"] < m_(s) < min(r["cluster"]["cr2"]["min_margin_pp"],
                                                                                            r["cluster"]["wild"]["min_margin_pp"]))(eq(s, "PSOBV-PSOC"))),
    ("X42", "The cluster-level verdict for PSO-VNS vs PSO does not depend on the variant: CR1 (\\NXEqClustCROnePSOVNSvsPSOCI), CR2 with "
            "5 d.f. and CR2 with Bell-McCaffrey d.f. (\\NXEqClustDfBM; \\NXEqClustBMPSOVNSvsPSOCI) all extend beyond -margin.",
     lambda s: (lambda c: 4.0 < c["df_bm"] <= 5.5 and all(c[v]["ci90"][0] < -m_(s) and not c[v]["equivalent"] for v in ("cr1", "cr2", "cr2_bm")))(
         eq(s, "PSOBV-PSOC")["cluster"])),
    ("X43", "Per case, \\NXHetBeyond of the \\NXHetN cases differ by more than +-margin, in both directions: \\NXHetBeyondPSOVNS favour "
            "PSO-VNS (by up to \\NXHetMaxPSOVNS pp) and \\NXHetBeyondPSO favour PSO (by up to \\NXHetMaxPSO pp); \\NXHetTies cases have "
            "identical case means.",
     lambda s: (lambda h: h["n_cases"] == 68 and h["beyond_first"] > 0 and h["beyond_second"] > 0
                and h["beyond_first"] + h["beyond_second"] == len(h["cases_beyond"]) and h["max_first"] > m_(s) and h["max_second"] > m_(s))(het(s))),
    ("X44", "The average equivalence rests on the small layouts: for N < 10 (\\NXEqSmallN cases) PSO-VNS and PSO are equivalent "
            "(90% CI \\NXEqSmallCI; Bayesian probability that equivalence is the most probable outcome \\NXEqSmallBayRope); for "
            "N >= 10 (\\NXEqLargeN cases) they are not (mean \\NXEqLargeMean pp, 90% CI \\NXEqLargeCI, which excludes zero but "
            "is NOT entirely beyond -margin; P(PSO-VNS better most probable) = \\NXEqLargeBayLeft).",
     lambda s: (lambda sm, lg: sm["equivalent"] and sm["bayes"]["p_rope"] > 0.99 and not lg["equivalent"] and lg["ci90"][1] < 0
                and lg["ci90"][1] > -m_(s) and lg["bayes"]["p_a_better"] > 0.5 and sm["n_cases"] + lg["n_cases"] == 68)(het(s)["small"], het(s)["large"])),
    ("X45", "Cases with a trivial wake loss dominate the average: in the \\NXEqTrivialN cases in which the case-mean wake loss of "
            "PSO-VNS is below \\NXTrivialLoss pp the two are equivalent (90% CI \\NXEqTrivialCI), in the other \\NXEqNontrivialN cases "
            "equivalence cannot be shown (90% CI \\NXEqNontrivialCI, minimal margin \\NXEqNontrivialMin pp; P(rope most probable) "
            "\\NXEqNontrivialBayRope), and every case beyond the margin is non-trivial.",
     lambda s: (lambda t, n_: t["equivalent"] and not n_["equivalent"] and t["n_cases"] + n_["n_cases"] == 68
                and n_["beyond_first"] + n_["beyond_second"] == het(s)["beyond_first"] + het(s)["beyond_second"])(het(s)["trivial"], het(s)["nontrivial"])),
    ("X46", "\\NXHetClBeyond of the six cluster means of PSO-VNS - PSO exceed the margin, in opposite directions (\\NXHetClBeyondPSOVNS "
            "favour PSO-VNS, \\NXHetClBeyondPSO favour PSO; Data Set I 500 m \\NXHetClDsIFive, Data Set II 1000 m \\NXHetClDsIIThousand pp).",
     lambda s: len(het(s)["clusters_beyond_first"]) > 0 and len(het(s)["clusters_beyond_second"]) > 0),
    ("X47", "Data set x density crossover (density \\NXDensDef): the difference between the data sets changes with density (slope "
            "\\NXDensInterSlope pp per unit density, studentized sign-flip p = \\NXDensInterP); within Data Set II PSO-VNS gains "
            "with density (slope \\NXDensSlopeDsII, p = \\NXDensPDsII), within Data Set I PSO gains (slope \\NXDensSlopeDsI, p = "
            "\\NXDensPDsI). All \\NXHetDsINPSO Data Set I cases beyond the margin favour PSO and are dense (\\NXHetDsIPSOList); the "
            "Data Set II cases beyond the margin favour PSO-VNS (\\NXHetDsIINPSOVNS) except \\NXHetDsIINPSO (\\NXHetDsIIPSOList).",
     lambda s: (lambda h: h["interaction"]["p_signflip"] < 0.05 and h["interaction"]["slope_pp_per_unit_phi"] < 0
                and h["slope_dsII"]["slope_pp_per_unit_phi"] < 0 and h["slope_dsII"]["p_perm"] < 0.05
                and h["slope_dsI"]["slope_pp_per_unit_phi"] > 0 and h["slope_dsI"]["p_perm"] < 0.05 and h["crossover"]
                and all(c["d"] > 0 and c["phi"] > 0.35 for c in h["cases_beyond"] if c["dataset"] == "1")
                and sum(c["d"] > 0 for c in h["cases_beyond"] if c["dataset"] == "2") <= 1)(het(s))),
    ("X48", "Bayesian wording: 'above 0.99' is the posterior probability that practical equivalence is the MOST PROBABLE of the three "
            "outcomes; the posterior mean shares are theta_PSO-VNS = \\NXBayMeanPSOVNSvsPSOLeft, theta_rope = \\NXBayMeanPSOVNSvsPSORope "
            "(95% credible interval \\NXBayPSOVNSvsPSOCIRope), theta_PSO = \\NXBayMeanPSOVNSvsPSORight, i.e. about one in four pairs of "
            "cases favours PSO-VNS by more than the margin.",
     lambda s: (lambda y: y["p_rope"] > 0.99 and y["mean_theta"][1] < 0.75 and 0.2 < y["mean_theta"][0] < 0.3)(eq(s, "PSOBV-PSOC")["bayes"])),
    ("X49", "The Hodges-Lehmann estimate of SSA-VNS - RSD-VNS is positive (\\NXHLSSAVNSvsRSDVNS pp, 95% CI \\NXHLSSAVNSvsRSDVNSCI, "
            "excludes zero, consistent with the Wilcoxon test), whereas the 95% CI of the mean contains zero.",
     lambda s: (lambda v: v["hl"] > 0 and v["ci95"][0] > 0 and v["ci95_mean"][0] < 0 < v["ci95_mean"][1] and v["wilcoxon_p"] < 0.05)(
         s["hodges_lehmann"]["SSABV-RSDVNS"])),
    ("X50", "The Hodges-Lehmann estimate of PSO-VNS - PSO is close to zero (\\NXHLPSOVNSvsPSO pp, 95% CI \\NXHLPSOVNSvsPSOCI inside "
            "+-margin and containing zero); the mean (\\NXEqMeanPSOVNSvsPSO pp) is larger in magnitude because a few large cases dominate it.",
     lambda s: (lambda v: v["ci95"][0] < 0 < v["ci95"][1] and inside(v["ci95"], m_(s)) and v["mean"] < v["hl"] < 0)(s["hodges_lehmann"]["PSOBV-PSOC"])),
    ("X51", "The Hodges-Lehmann estimates of the other key pairs agree in sign with the Wilcoxon verdicts and their 95% CIs exclude "
            "zero: LX-SSA-VNS worse than RSD-VNS (\\NXHLLXSSAVNSvsRSDVNS, \\NXHLLXSSAVNSvsRSDVNSCI), PSO-VNS better than RSD-VNS "
            "(\\NXHLPSOVNSvsRSDVNS), SSA-VNS better than RS-VNS (\\NXHLSSAVNSvsRSVNS) and than LX-SSA-VNS (\\NXHLSSAVNSvsLXSSAVNS).",
     lambda s: all((lambda v: v["excludes_zero"] and np.sign(v["hl"]) == sg and v["wilcoxon_p"] < 0.05)(s["hodges_lehmann"][k])
                   for k, sg in (("LXBV-RSDVNS", 1), ("PSOBV-RSDVNS", -1), ("SSABV-RSVNS", -1), ("SSABV-LXBV", -1)))),
    ("X52", "Feasibility on Horns Rev (one site, 16 turbines, random starts, 6,030 evaluations): PSO-VNS is feasible in every run, PSO "
            "not; all \\NXHRMcNemarOnlyPSOVNS discordant seeds favour PSO-VNS, none favours PSO (\\NXHRMcNemarOnlyPSO), exact McNemar "
            "p = \\NXHRMcNemarP.",
     lambda s: (lambda c: c["only_second"] == 0 and c["only_first"] > 0 and c["p_exact"] < 0.01 and c["matches_summary"]
                and c["feasible_first"] == c["n_seeds"])(s["hr_mcnemar"])),
    ("X53", "D14 (one-sided): a gain of the SSA phase over sampling in the disc larger than the margin is ruled out at the seed and case "
            "levels (lower 90% limits of SSA-VNS - RSD-VNS above -margin: \\NXEqSeedSSAVNSvsRSDVNSCI, \\NXEqCaseSSAVNSvsRSDVNSCI) but "
            "NOT at the cluster level, where the lower limit lies just beyond it (CR2 \\NXEqClustSSAVNSvsRSDVNSCI, wild "
            "\\NXEqClustSSAVNSvsRSDVNSCIWild).",
     lambda s: (lambda r: r["seed"]["ci90"][0] > -m_(s) and r["case"]["ci90"][0] > -m_(s) and r["cluster"]["cr2"]["ci90"][0] < -m_(s)
                and r["cluster"]["wild"]["ci90"][0] < -m_(s))(eq(s, "SSABV-RSDVNS"))),
    ("X54", "SSA-VNS vs RSD-VNS on the fixed benchmark: equivalent at +-margin (seed-level minimal margin \\NXEqSeedSSAVNSvsRSDVNSMin pp) "
            "with the 90% CI above zero (\\NXEqSeedSSAVNSvsRSDVNSCI) but the 95% CI not (\\NXEqSeedSSAVNSvsRSDVNSCINinetyFive); at the "
            "case and cluster levels equivalence is not shown (\\NXEqCaseSSAVNSvsRSDVNSMin, \\NXEqClustSSAVNSvsRSDVNSMin pp).",
     lambda s: (lambda r: r["seed"]["equivalent"] and r["seed"]["ci90"][0] > 0 and r["seed"]["ci95"][0] < 0
                and not r["case"]["equivalent"] and not r["cluster"]["equivalent"])(eq(s, "SSABV-RSDVNS"))),
    ("X55", "LX-SSA-VNS and RS-VNS are equivalent at the seed and case levels (minimal margins \\NXEqSeedLXSSAVNSvsRSVNSMin, "
            "\\NXEqCaseLXSSAVNSvsRSVNSMin pp) but not at the cluster level (\\NXEqClustLXSSAVNSvsRSVNSMin pp); no pair of "
            "tab:equivalence is equivalent at the cluster level (\\NXEqClustNHolds of 18; seed \\NXEqSeedNHolds, case \\NXEqCaseNHolds).",
     lambda s: (lambda r: r["seed"]["equivalent"] and r["case"]["equivalent"] and not r["cluster"]["equivalent"])(eq(s, "LXBV-RSVNS"))
     and sum(r["cluster"]["equivalent"] for r in s["equivalence_levels"]["pairs"].values()) == 0
     and len(s["equivalence_levels"]["pairs"]) == 18),
    ("X56", "In benchmark energy the margin corresponds to \\NXMarginMWhTurbDsI MWh/yr per turbine in Data Set I and \\NXMarginMWhTurbDsII "
            "MWh/yr in Data Set II (wake-free benchmark AEP \\NXWakeFreeMWhTurbDsI and \\NXWakeFreeMWhTurbDsII MWh/yr per turbine); the "
            "mean PSO-VNS - PSO gain is \\NXEnPSOVNSvsPSOMWhTurb MWh/yr per turbine, far below both.",
     lambda s: all(abs(v["margin_mwh_yr_per_turbine"] - m_(s) / 100 * v["wakefree_mwh_yr_per_turbine"]) < 1e-9
                   and v["wakefree_per_turbine_spread"] < 1e-6 for v in het(s)["margin_energy"].values())
     and abs(s["energy"]["PSOBV-PSOC"]["mean_gain_mwh_yr_per_turbine"]) < 0.5 * min(v["margin_mwh_yr_per_turbine"] for v in het(s)["margin_energy"].values())),
    ("X57", "Figure equiv_curve.pdf (TOST p vs margin, PSO-VNS vs PSO, seed / case / cluster level) was generated in the same run.",
     lambda s: fig_ok(s) and len(s["equivalence_levels"]["curve"]["margins"]) == len(s["equivalence_levels"]["curve"]["cluster_wild"])),
    ("X28", "mpce_numbers_extra.tex and mpce_supp_inference.tex were generated in the same run as mpce_summary_extra.json.",
     lambda s: all(s["generated"] in open(os.path.join(HERE, f)).readline()
                   for f in ("mpce_numbers_extra.tex", "mpce_supp_inference.tex"))),
]


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--summary", default=os.path.join(HERE, "mpce_summary_extra.json"))
    args = ap.parse_args(argv)
    s = json.load(open(args.summary))
    global EQ_MARGIN
    EQ_MARGIN = s["equivalence_levels"]["margin_pp"]
    macros = set(re.findall(r"\\newcommand\{\\(NX[A-Za-z]+)\}", open(os.path.join(HERE, "mpce_numbers_extra.tex")).read()))
    nf = 0
    for cid, sentence, cond in CHECKS:
        try:
            ok = bool(cond(s))
        except Exception as e:                                    # pragma: no cover
            ok = False; sentence += f"  [error: {e!r}]"
        missing = [m for m in re.findall(r"\\(NX[A-Za-z]+)", sentence) if m not in macros]
        if missing:
            ok = False; sentence += f"  [undefined macros: {missing}]"
        nf += not ok
        print(f"{'PASS' if ok else 'FAIL'}  {cid}  {sentence}")
    print(f"\n{len(CHECKS) - nf} PASS / {nf} FAIL")
    return 1 if nf else 0


if __name__ == "__main__":
    sys.exit(main())

"""Checks X01... of the inference-robustness analyses (Phase 6, W2) against mpce_summary_extra.json.

Usage:  python3 analysis/mpce_check_extra.py [--summary analysis/mpce_summary_extra.json]

Each check is ONE sentence that the paper (main text or supplement) may state; the condition is evaluated on the
output of mpce_inference_extra.py and printed as PASS / FAIL. A FAIL means the sentence must not be written (or
must be rewritten). Exit code 1 if any check fails. The sentences use the \\NX... macros of mpce_numbers_extra.tex
for every number.
"""
import os, re, sys, json, argparse

HERE = os.path.dirname(os.path.abspath(__file__))
SIG_PAIRS_EXPECTED_NS = {"PSOBV-PSOC", "LXBV-RSVNS"}          # the two pairs the paper calls not significant


def thr(s):
    return [s["threshold"][str(t)] for t in s["thresholds"]]


def case_sig(s):
    return {k for k, v in s["case_level"].items() if v["p"] < 0.05}


def mt(s, key):
    return next(t for t in s["multiplicity"]["tests"] if t["key"] == key)


CHECKS = [
    ("X01", "The re-implementation reproduces the pipeline exactly at the paper's 15-run threshold (average ranks and every "
            "case-mean p of the main table, component analysis, N>=10 subgroups and budget split equal mpce_summary.json).",
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
            "strongest result attainable with six clusters, exact p = \\NXClMinP).",
     lambda s: all(min(s["cluster"][k]["clusters_favour_first"], s["cluster"][k]["clusters_favour_second"]) == 0
                   and s["cluster"][k]["wilcoxon_exact_p"] <= s["cluster_min_attainable_p"] + 1e-12 for k in case_sig(s))
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
    ("X14", "Leaving out any one of the six clusters changes none of the fifteen case-mean verdicts (significant or not, "
            "and which method is better when significant), PSO-VNS keeps the best average rank, and PSO remains the only method not significantly behind it.",
     lambda s: all(v["verdict_changes"] == 0 for v in s["loco"]["pairs"].values())
     and s["loco"]["best_ranked_always_PSOBV"] and s["loco"]["posthoc_only_PSOC_always"]),
    ("X15", "In the leave-one-cluster-out runs the sign of the mean difference changes only for the two non-significant pairs "
            "(PSO-VNS vs PSO, LX-SSA-VNS vs RS-VNS); PSO-VNS vs PSO stays non-significant (p >= \\NXLocoPSOVNSvsPSOPMin).",
     lambda s: {k for k, v in s["loco"]["pairs"].items() if v["direction_flips"] > 0} <= SIG_PAIRS_EXPECTED_NS
     and s["loco"]["pairs"]["PSOBV-PSOC"]["p_min"] >= 0.05),
    ("X16", "The cluster-bootstrap 95% CI excludes zero for every comparison that is significant on the case means and contains "
            "zero for the two that are not.",
     lambda s: all(s["cluster"][k]["cluster_boot_excludes_zero"] == (k in case_sig(s)) for k in s["cluster"])),
    ("X17", "CAVEAT: SSA-VNS vs RS-VNS is the only case-level-significant comparison that is not significant with the "
            "cluster-robust t test (5 d.f., p = \\NXClSSAVNSvsRSVNSCRP; CI \\NXClSSAVNSvsRSVNSCRCI), although all six clusters "
            "favour SSA-VNS; the SSA gain over random sampling is consistent but small.",
     lambda s: {k for k in case_sig(s) if s["cluster"][k]["cr1_p"] >= 0.05} == {"SSABV-RSVNS"}
     and s["cluster"]["SSABV-RSVNS"]["clusters_favour_first"] == 6),
    ("X18", "With one Holm correction over all \\NXMultN case-mean tests quoted in the main text, every test significant "
            "unadjusted stays significant except the post hoc N>=10 subgroup over both data sets (same under Benjamini-Hochberg).",
     lambda s: s["multiplicity"]["lost_under_holm"] == ["sub:Large"] and s["multiplicity"]["lost_under_bh"] == ["sub:Large"]),
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
    ("X28", "mpce_numbers_extra.tex and mpce_supp_inference.tex were generated in the same run as mpce_summary_extra.json.",
     lambda s: all(s["generated"] in open(os.path.join(HERE, f)).readline()
                   for f in ("mpce_numbers_extra.tex", "mpce_supp_inference.tex"))),
]


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--summary", default=os.path.join(HERE, "mpce_summary_extra.json"))
    args = ap.parse_args(argv)
    s = json.load(open(args.summary))
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

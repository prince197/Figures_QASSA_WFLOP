"""Checks F01... of the direction-resolution, PyWake and IEA37-projection analyses (sensitivity analyses).

Usage:  python3 analysis/mpce_check_dir.py [--summary analysis/mpce_summary_dir.json]

Each check is ONE sentence that the paper (main text or supplement) may write, with the \\NF... macros of
mpce_numbers_dir.tex for every number; the condition is evaluated on the outputs of mpce_direction.py,
pywake_check.py and iea37_projected.py and printed as PASS / FAIL. A FAIL means the sentence must not be written
(or must be rewritten). Exit code 1 if any check fails.
"""
import os, re, sys, json, argparse
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
M = 0.05                                                   # equivalence margin (pp), as mpce_results.EQ_MARGIN
VNS_BASED = ["PSOBV", "SSABV", "LXBV", "RSVNS", "BVNS"]
NON_VNS = ["PSOC", "SSA", "LXSSA", "DE", "SLSQP"]


def R(s, k):
    return s["benchmark"]["res"][k]


def pair(s, k):
    return R(s, k)["pso_pair"]


def hr(s):
    return s["horns_rev"]


def hp(s, bud, b, c):
    return hr(s)["pairs_joint_feasible"][f"{bud}_random_{b}_{c}"]


def pw(P, farm, bins, model, src="PyWake"):
    r = P[(P.Farm == farm) & (P.Bins == bins) & (P.Model == model) & (P.Source == src)]
    return r.iloc[0]


def macro(fn, name):
    t = open(fn).read()
    m = re.search(r"\\newcommand\{\\%s\}\{(.*?)\}(?:\s|$)" % name, t)
    return m.group(1) if m else None


def paper_num(name):
    v = macro(os.path.join(HERE, "mpce_numbers.tex"), name)
    return float(v.replace("\\ensuremath{-}", "-").replace("{,}", "")) if v else None


CHECKS = [
    # ------------------------------------------------------------------ benchmark: method and reproduction
    ("F01", "Splitting every 15-deg bin into one sub-bin reproduces the benchmark objective (median relative deviation < 1e-6; "
            "deviations only from the 1-mm rounding of the stored coordinates), the wake-free objective is identical for 1, 3 and 15 "
            "sub-bins, and the case-level statistics of the re-evaluation equal those of the recorded objective (PSO-VNS - PSO "
            "within 0.001 pp, no sign change, same best-ranked method).",
     lambda s, P, I: s["benchmark"]["reproduction"]["median_rel_dev"] < 1e-6 and s["benchmark"]["reproduction"]["ideal_equal_all_subs"]
     and s["benchmark"]["reproduction"]["ideal_matches_recorded"]
     and abs(pair(s, "J1")["all"]["mean_dloss_pp"] - pair(s, "rec")["all"]["mean_dloss_pp"]) < 1e-3
     and pair(s, "J1")["sign_changes_vs_recorded"] == 0 and R(s, "J1")["friedman"]["best_ranked"] == "PSOBV"),
    ("F02", "The 15-deg statistics of this analysis reproduce mpce_summary.json (average ranks of the eight methods, PSO-VNS - PSO "
            "case-mean difference and p, the MS-SLSQP gap).",
     lambda s, P, I: (lambda m: all(abs(R(s, "rec")["friedman"]["avg_rank"][a] - m["main"]["friedman"]["avg_rank"][a]) < 1e-9
                                    for a in m["main"]["friedman"]["avg_rank"])
                      and abs(pair(s, "rec")["all"]["mean_dloss_pp"] - m["main"]["case_mean_wilcoxon"]["PSOC"]["mean_dloss_pp"]) < 1e-9
                      and abs(pair(s, "rec")["all_wilcoxon"]["p"] - m["main"]["case_mean_wilcoxon"]["PSOC"]["p"]) < 1e-9
                      and abs(R(s, "rec")["case_mean_wilcoxon"]["SLSQP"]["mean_dloss_pp"] - m["main"]["case_mean_wilcoxon"]["SLSQP"]["mean_dloss_pp"]) < 1e-9)
     (json.load(open(os.path.join(HERE, "mpce_summary.json"))))),
    # ------------------------------------------------------------------ benchmark: level shift
    ("F03", "With 1-deg sub-bins the mean wake loss of every one of the eight methods rises, by \\NFLossRiseMin to \\NFLossRiseMax pp "
            "(\\NFLossRiseRelMin to \\NFLossRiseRelMax % of its benchmark value); for PSO-VNS from \\NFLossPSOVNSFifteen % to "
            "\\NFLossPSOVNSOne %, about \\NFLossRiseMarginRatio times the +-0.05 pp equivalence margin.",
     lambda s, P, I: all(R(s, "J15")["mean_loss_qualified"][a] > R(s, "rec")["mean_loss_qualified"][a] for a in R(s, "rec")["mean_loss_qualified"])
     and (R(s, "J15")["mean_loss_qualified"]["PSOBV"] - R(s, "rec")["mean_loss_qualified"]["PSOBV"]) / M >= 15),
    ("F04", "Most of the change appears already with 5-deg sub-bins (for PSO-VNS more than 70 % of the 15-to-1-deg increase).",
     lambda s, P, I: (R(s, "J3")["mean_loss_qualified"]["PSOBV"] - R(s, "rec")["mean_loss_qualified"]["PSOBV"])
     / (R(s, "J15")["mean_loss_qualified"]["PSOBV"] - R(s, "rec")["mean_loss_qualified"]["PSOBV"]) > 0.7),
    ("F05", "The finer rose penalizes the VNS-based methods and PSO most and DE and MS-SLSQP least (\\NFLossRiseMinBy smallest, "
            "increase of DE and MS-SLSQP below 0.75 pp, of PSO-VNS, PSO and VNS above 1.0 pp).",
     lambda s, P, I: (lambda r: min(r, key=r.get) == "DE" and r["DE"] < 0.75 and r["SLSQP"] < 0.75
                      and min(r["PSOBV"], r["PSOC"], r["BVNS"]) > 1.0)
     ({a: R(s, "J15")["mean_loss_qualified"][a] - R(s, "rec")["mean_loss_qualified"][a] for a in R(s, "rec")["mean_loss_qualified"]})),
    # ------------------------------------------------------------------ benchmark: ranking
    ("F06", "PSO-VNS remains the best-ranked of the eight methods with 5-deg and with 1-deg sub-bins (average rank \\NFRankPSOVNSOne vs "
            "PSO \\NFRankPSOOne at 1 deg), and PSO remains the only method whose average rank is not significantly different "
            "(Holm).",
     lambda s, P, I: R(s, "J3")["friedman"]["best_ranked"] == "PSOBV" and R(s, "J15")["friedman"]["best_ranked"] == "PSOBV"
     and [a for a, p in R(s, "J15")["friedman"]["p_holm_vs_focus"].items() if p >= 0.05] == ["PSOC"]),
    ("F07", "The rest of the order changes: MS-SLSQP moves from sixth to third (average rank \\NFRankMSSLSQPFifteen to "
            "\\NFRankMSSLSQPOne) and VNS from fourth to sixth; the within-case ordering agrees only moderately with the benchmark "
            "(mean Kendall tau \\NFTauOne; same best method in \\NFSameBestOne % of the cases).",
     lambda s, P, I: R(s, "rec")["rank_order"].index("SLSQP") == 5 and R(s, "J15")["rank_order"].index("SLSQP") == 2
     and R(s, "rec")["rank_order"].index("BVNS") == 3 and R(s, "J15")["rank_order"].index("BVNS") == 5
     and R(s, "J15")["tau_within_case"] < 0.7 and R(s, "J15")["same_best_pct"] < 50),
    ("F08", "The gap of MS-SLSQP to PSO-VNS shrinks from \\NFSLSQPGapFifteen to \\NFSLSQPGapOne pp (to less than a third) but "
            "remains significant (Holm).",
     lambda s, P, I: (lambda a, b: b < a / 3 and b > 0 and R(s, "J15")["case_mean_wilcoxon"]["SLSQP"]["p_holm"] < 0.05)
     (-R(s, "rec")["case_mean_wilcoxon"]["SLSQP"]["mean_dloss_pp"], -R(s, "J15")["case_mean_wilcoxon"]["SLSQP"]["mean_dloss_pp"])),
    ("F09", "Under the Gaussian wake the best-ranked method is PSO with the 15-deg bins but PSO-VNS with 1-deg sub-bins, so the PSO "
            "lead of the Gaussian re-evaluation is itself a property of the coarse rose.",
     lambda s, P, I: R(s, "G1")["friedman"]["best_ranked"] == "PSOC" and R(s, "G15")["friedman"]["best_ranked"] == "PSOBV"),
    # ------------------------------------------------------------------ benchmark: PSO-VNS vs PSO
    ("F10", "With 1-deg sub-bins, PSO-VNS - PSO is \\NFPairMeanOne pp (90 % CI \\NFPairCIOne): the average difference is still "
            "smaller than the margin, but the CI extends beyond -0.05 pp, so equivalence at +-0.05 pp can no longer be shown "
            "(smallest margin \\NFPairMinMarginOne pp), and PSO-VNS is significantly better on the case means (p = \\NFPairPOne).",
     lambda s, P, I: (lambda e, p: -M < e["mean_dloss_pp"] < 0 and not e["equivalent"] and e["ci90_mean_dloss_pp"][0] < -M
                      and e["ci90_mean_dloss_pp"][1] < 0 and p < 0.05)(pair(s, "J15")["all"], pair(s, "J15")["all_wilcoxon"]["p"])),
    ("F11", "The same holds with 5-deg sub-bins (\\NFPairMeanFive pp, 90 % CI \\NFPairCIFive; not equivalent; p = \\NFPairPFive).",
     lambda s, P, I: (lambda e: not e["equivalent"] and e["ci90_mean_dloss_pp"][1] < 0)(pair(s, "J3")["all"])
     and pair(s, "J3")["all_wilcoxon"]["p"] < 0.05),
    ("F12", "The per-case sign of PSO-VNS - PSO changes in \\NFSignChangesOne of the \\NFSignChangesBase cases with a nonzero "
            "difference under both resolutions (more than a third).",
     lambda s, P, I: pair(s, "J15")["sign_changes_vs_recorded"] * 3 > pair(s, "J15")["n_cases_nonzero_both"]),
    ("F13", "On the \\NFNonTrivN non-trivial cases (PSO-VNS wake loss >= \\NFNonTrivThr %), equivalence is shown neither with the "
            "benchmark bins (90 % CI \\NFNonTrivCIFifteen) nor with 1-deg sub-bins (\\NFNonTrivCIOne, where PSO-VNS is "
            "significantly better, p = \\NFNonTrivPOne).",
     lambda s, P, I: not pair(s, "rec")["nontrivial"]["equivalent"] and not pair(s, "J15")["nontrivial"]["equivalent"]
     and pair(s, "J15")["nontrivial"]["ci90_mean_dloss_pp"][1] < 0 and pair(s, "J15")["nontrivial_wilcoxon_p"] < 0.05
     and s["benchmark"]["n_nontrivial"] == 44),
    ("F14", "The Data Set II N >= 10 advantage of PSO-VNS survives the finer rose: lower mean wake loss in \\NFDSIINTenWinsOne of "
            "\\NFDSIINTenN cases (\\NFDSIINTenMeanOne pp, p = \\NFDSIINTenPOne), as with 15-deg bins (\\NFDSIINTenMeanFifteen pp).",
     lambda s, P, I: (lambda d: d["focus_better"] == d["n"] == 10 and d["wilcoxon_p"] < 0.05 and d["mean"] < 0)(pair(s, "J15")["ds2_n10"])
     and pair(s, "rec")["ds2_n10"]["focus_better"] == 10),
    ("F15", "With 1-deg sub-bins the largest per-case advantage of PSO over PSO-VNS falls from \\NFPairMaxPSOBetterFifteen to "
            "\\NFPairMaxPSOBetterOne pp and PSO-VNS has the lower loss in \\NFPairWinsOne of 68 cases (\\NFPairLossesOne for PSO).",
     lambda s, P, I: pair(s, "J15")["max_other_better"] < pair(s, "rec")["max_other_better"]
     and pair(s, "J15")["focus_better"] > pair(s, "rec")["focus_better"]),
    # ------------------------------------------------------------------ benchmark: component analysis
    ("F16", "The component-analysis conclusions hold with 1-deg sub-bins: SSA-VNS and LX-SSA-VNS have a higher case-mean wake loss "
            "than the disc-sampling control RSD-VNS (p = \\NFAblSSAVNSvsRSDVNSOneP, \\NFAblLXSSAVNSvsRSDVNSOneP), PSO-VNS is better "
            "than RSD-VNS, RS-VNS and VNS, and PSO-VNS keeps the best average rank of the nine variants.",
     lambda s, P, I: (lambda c: c["SSABV-RSDVNS"]["mean_dloss_pp"] > 0 and c["SSABV-RSDVNS"]["p"] < 0.05
                      and c["LXBV-RSDVNS"]["mean_dloss_pp"] > 0 and c["LXBV-RSDVNS"]["p"] < 0.05
                      and all(c[k]["mean_dloss_pp"] < 0 and c[k]["p"] < 0.05 for k in ("PSOBV-RSDVNS", "PSOBV-RSVNS", "PSOBV-BVNS")))
     (R(s, "J15")["ablation_contrasts"]) and R(s, "J15")["ablation_friedman"]["best_ranked"] == "PSOBV"),
    ("F17", "The small gain of SSA-VNS over the square-sampling control RS-VNS is not significant with 1-deg sub-bins "
            "(\\NFAblSSAVNSvsRSVNSOne pp, p = \\NFAblSSAVNSvsRSVNSOneP).",
     lambda s, P, I: R(s, "J15")["ablation_contrasts"]["SSABV-RSVNS"]["p"] >= 0.05
     and R(s, "rec")["ablation_contrasts"]["SSABV-RSVNS"]["p"] < 0.05),
    # ------------------------------------------------------------------ Horns Rev
    ("F20", "The 1-deg Horns Rev bins give every 30-deg sector exactly 30 bins with total frequency 1; the paper's 5-deg setting of "
            "this script reproduces hornsrev_model.py exactly for the installed block and, for the stored runs, up to the 1-cm "
            "rounding of the stored coordinates (median deviation < 1e-4 GWh/yr; fewer than 12 % of the 901 feasible runs lower by "
            "more than 0.01 GWh/yr, none higher; the count of runs above the installed block, 16, is the same for recorded and "
            "re-evaluated values). The 5-deg vs 1-deg comparisons therefore use re-evaluated values of the same coordinates.",
     lambda s, P, I: hr(s)["bins"]["1deg_0.5"]["bins_per_sector"] == [30] and abs(hr(s)["bins"]["1deg_0.5"]["freq_sum"] - 1) < 1e-12
     and hr(s)["bins"]["5deg_2.5"]["bins_per_sector"] == [6] and hr(s)["installed_reproduces_hornsrev_model"]
     and hr(s)["reproduction"]["median_abs_dev_gwh"] < 1e-4 and hr(s)["reproduction"]["n_abs_dev_gt_0_01"] < 0.12 * 901
     and hr(s)["reproduction"]["n_higher_than_recorded_gt_0_01"] <= 1
     and hr(s)["reproduction"]["runs_above_installed_recorded"] == hr(s)["above_installed_total"]["5deg_2.5"] == 16
     and hr(s)["n_runs"] == 970 and hr(s)["n_feasible"] == 901),
    ("F21", "With 1-deg bins no feasible run of any method at any budget exceeds the installed block (\\NFHRAboveOneTotal of "
            "\\NFHRNFeas, against \\NFHRAboveFiveTotal with the paper's 5-deg bins); the same holds with 1-deg bins at integer "
            "centres and with 5-deg bins centred at 0 deg.",
     lambda s, P, I: hr(s)["above_installed_total"]["1deg_0.5"] == 0 and hr(s)["above_installed_total"]["1deg_0"] == 0
     and hr(s)["above_installed_total"]["5deg_0"] == 0 and hr(s)["above_installed_total"]["5deg_2.5"] > 0),
    ("F22", "Under PyWake's NOJ model (1-deg bins) only \\NFHRAbovePyWakeTotal run exceeds the installed block, a DE run at "
            "120,030 evaluations; no PSO-VNS run does.",
     lambda s, P, I: hr(s)["above_installed_total"].get("PyWakeNOJ_1deg") == 1
     and hr(s)["above_installed_by_method"]["PyWakeNOJ_1deg"]["DE"] == 1
     and hr(s)["methods"]["120030_random_DE"]["above_PyWakeNOJ_1deg"] == 1),
    ("F23", "The installed block loses \\NFHRInstDrop GWh/yr under the finer rose, the optimized layouts on average \\NFHRMeanDrop "
            "GWh/yr; the VNS-based methods lose most (every VNS-based method more than every other method).",
     lambda s, P, I: hr(s)["mean_drop_gwh_all"] > 2 * hr(s)["installed_drop_gwh"]
     and min(hr(s)["mean_drop_gwh_by_method"][a] for a in VNS_BASED) > max(hr(s)["mean_drop_gwh_by_method"][a] for a in NON_VNS)),
    ("F24", "At 6,030 evaluations (random starts) the AEP advantage of PSO-VNS over PSO on the \\NFHRPairSixKOneN jointly feasible "
            "seeds is significant only with the paper's 5-deg bins (\\NFHRPairSixKFive GWh/yr, p = \\NFHRPairSixKFiveP); with "
            "1-deg bins (\\NFHRPairSixKOne, p = \\NFHRPairSixKOneP) and with PyWake NOJ (\\NFHRPairSixKPyWake, p = "
            "\\NFHRPairSixKPyWakeP) PSO has the (non-significantly) higher mean; the robust 6,030 result is feasibility.",
     lambda s, P, I: hp(s, 6030, "PSOC", "5deg_2.5")["p"] < 0.05 and hp(s, 6030, "PSOC", "5deg_2.5")["mean_diff_gwh"] > 0
     and hp(s, 6030, "PSOC", "1deg_0.5")["p"] >= 0.05 and hp(s, 6030, "PSOC", "1deg_0.5")["mean_diff_gwh"] < 0
     and hp(s, 6030, "PSOC", "PyWakeNOJ_1deg")["p"] >= 0.05 and hp(s, 6030, "PSOC", "PyWakeNOJ_1deg")["mean_diff_gwh"] < 0
     and hp(s, 6030, "PSOC", "1deg_0.5")["n"] == 19),
    ("F25", "With 1-deg bins PSO-VNS no longer has the highest mean AEP at 6,030 evaluations (random starts): PSO "
            "(\\NFHRSixKRBestOneAEP GWh/yr on its \\NFHRSixKRPSOFeas feasible runs) is marginally ahead of PSO-VNS "
            "(\\NFHRSixKRSecondOneAEP).",
     lambda s, P, I: hr(s)["highest_mean"]["6030_random_1deg_0.5"]["best_qualified"] == "PSOC"
     and hr(s)["highest_mean"]["6030_random_1deg_0.5"]["second_qualified"] == "PSOBV"
     and hr(s)["highest_mean"]["6030_random_5deg_2.5"]["best_qualified"] == "PSOBV"),
    ("F26", "At 30,030 evaluations PSO-VNS keeps the highest mean AEP with 1-deg bins (\\NFHRThirtyKBestOneAEP vs "
            "\\NFHRThirtyKSecondOne \\NFHRThirtyKSecondOneAEP GWh/yr) and under PyWake NOJ (\\NFHRThirtyKBestPyWakeAEP), and its "
            "advantage over PSO on the \\NFHRPairThirtyKOneN jointly feasible seeds is significant under both (p = "
            "\\NFHRPairThirtyKOneP, \\NFHRPairThirtyKPyWakeP).",
     lambda s, P, I: all(hr(s)["highest_mean"][f"30030_random_{c}"]["best_qualified"] == "PSOBV"
                         and hr(s)["highest_mean"][f"30030_random_{c}"]["best_any"] == "PSOBV"
                         for c in ("5deg_2.5", "1deg_0.5", "1deg_0", "PyWakeNOJ_1deg"))
     and all(hp(s, 30030, "PSOC", c)["p"] < 0.05 and hp(s, 30030, "PSOC", c)["mean_diff_gwh"] > 0 for c in ("1deg_0.5", "PyWakeNOJ_1deg"))),
    ("F27", "At 30,030 evaluations PSO-VNS is significantly better than every other method in the paper's run-level test with "
            "1-deg bins (Holm p <= \\NFHRThirtyKMaxPOne), but on the jointly feasible seeds the comparison with DE rests on 8 runs "
            "and is not significant.",
     lambda s, P, I: all(v["outcome"] == "W" for v in hr(s)["paper_test"]["30030_random_1deg_0.5"].values())
     and hp(s, 30030, "DE", "1deg_0.5")["n"] == 8 and hp(s, 30030, "DE", "1deg_0.5")["p"] >= 0.05),
    ("F28", "With 5-deg bins centred at 0 deg instead of 2.5 deg, the PSO-VNS advantage over PSO at 30,030 is not significant.",
     lambda s, P, I: hp(s, 30030, "PSOC", "5deg_0")["p"] >= 0.05),
    ("F29", "At 120,030 evaluations DE has the highest mean AEP under 5-deg bins, 1-deg bins and PyWake NOJ.",
     lambda s, P, I: all(hr(s)["highest_mean"][f"120030_random_{c}"]["best_qualified"] == "DE" for c in ("5deg_2.5", "1deg_0.5", "PyWakeNOJ_1deg"))),
    ("F30", "The wake loss of the installed block depends on the bin phase: \\NFHRInstLossFive % with 5-deg bins centred at 2.5 deg, "
            "\\NFHRInstLossFivePhaseZero % centred at 0 deg, \\NFHRInstLossOne % with 1-deg bins.",
     lambda s, P, I: hr(s)["installed_loss_pct"]["5deg_0"] - hr(s)["installed_loss_pct"]["5deg_2.5"] > 0.3),
    # ------------------------------------------------------------------ PyWake
    ("F31", "PyWake \\NFPyWakeVersion was run at bins identical to ours (bin probabilities equal to \\NFPyWakePMaxDiff) and gives the "
            "same wake-free AEP (\\NFPyWakeIdealAEP GWh/yr).",
     lambda s, P, I: P is not None and P[P.Source == "PyWake"].PMaxAbsDiff.max() < 1e-12
     and abs(pw(P, "HR80", "ours_5deg_2.5", "NOJ_k0.04").IdealAEP_GWh - pw(P, "HR80", "ours_5deg_2.5", "paper_model", "ours").IdealAEP_GWh) < 1e-6
     and str(pw(P, "HR80", "ours_5deg_2.5", "NOJ_k0.04").Version) == "2.6.20"),
    ("F32", "For the installed 80-turbine farm, PyWake NOJ (k = 0.04) gives \\NFPyWakeAEP GWh/yr (wake loss \\NFPyWakeLossPct %) and "
            "our model \\NFPyWakeOurAEP GWh/yr (\\NFPyWakeOurLossPct %): our AEP is \\NFPyWakeDiffPct % higher (less than 1 %), "
            "our wake loss \\NFPyWakeLossDiffPP pp lower.",
     lambda s, P, I: (lambda v, o: 0 < 100 * (o.AEP_GWh / v.AEP_GWh - 1) < 1.0 and o.WakeLossPct < v.WakeLossPct)
     (pw(P, "HR80", "ours_5deg_2.5", "NOJ_k0.04"), pw(P, "HR80", "ours_5deg_2.5", "paper_model", "ours"))),
    ("F33", "Rotor averaging does not explain the difference: PyWake's hub-centre variant (\\NFPyWakeRotorCenterAEP GWh/yr) is further "
            "from our value than the default NOJ.",
     lambda s, P, I: abs(pw(P, "HR80", "ours_5deg_2.5", "NOJ_k0.04_RotorCenter").AEP_GWh - pw(P, "HR80", "ours_5deg_2.5", "paper_model", "ours").AEP_GWh)
     > abs(pw(P, "HR80", "ours_5deg_2.5", "NOJ_k0.04").AEP_GWh - pw(P, "HR80", "ours_5deg_2.5", "paper_model", "ours").AEP_GWh)),
    ("F34", "The difference comes from the thrust coefficient taken at the free-stream speed: with C_T at the local (waked) speed our "
            "model reproduces PyWake's NOJ with hub-centre test and momentum-theory induction exactly (\\NFOurLocalCTAEP GWh/yr; "
            "all bin settings and both farms agree to 1e-6 GWh/yr), and the remaining gap to PyWake's default NOJ "
            "(< 0.25 % of AEP) is rotor averaging and the Madsen a(C_T) polynomial.",
     lambda s, P, I: all(abs(pw(P, f, b, "NOJ_k0.04_RotorCenter_momentum").AEP_GWh - pw(P, f, b, "paper_model_localCT", "ours").AEP_GWh) < 1e-6
                         for f in ("HR80", "HR16") for b in ("ours_5deg_2.5", "1deg_0.5", "pywake_default_1deg_0"))
     and abs(pw(P, "HR80", "ours_5deg_2.5", "NOJ_k0.04").AEP_GWh / pw(P, "HR80", "ours_5deg_2.5", "paper_model_localCT", "ours").AEP_GWh - 1) < 0.0025),
    ("F35", "For the installed 16-turbine block PyWake NOJ gives \\NFPyWakeSixteenAEP GWh/yr against our \\NFPyWakeSixteenOurAEP "
            "(\\NFPyWakeSixteenDiffPct %).",
     lambda s, P, I: 0 < pw(P, "HR16", "ours_5deg_2.5", "paper_model", "ours").AEP_GWh / pw(P, "HR16", "ours_5deg_2.5", "NOJ_k0.04").AEP_GWh - 1 < 0.01),
    ("F36", "The PyWake reference constant of mpce_results.py (662.5 GWh/yr, 10.96 %) is not reproduced by PyWake at our bins "
            "(it lies between the 5-deg and PyWake-default values and differs from both by more than 0.5 GWh/yr).",
     lambda s, P, I: abs(pw(P, "HR80", "ours_5deg_2.5", "NOJ_k0.04").AEP_GWh - 662.5) > 0.5
     and abs(pw(P, "HR80", "pywake_default_1deg_0", "NOJ_k0.04").AEP_GWh - 662.5) > 0.5),
    # ------------------------------------------------------------------ IEA37
    ("F40", "The projection script reproduces iea37_published_results.csv (AEP of every published layout, the 1-mm feasibility flags "
            "and the boundary-projected AEPs).",
     lambda s, P, I: all(I["consistency"][k] for k in ("aep_matches_file", "projected_aep_matches_file", "strict_matches_file"))
     and I["consistency"]["n_matched"] == 24),
    ("F41", "\\NFIEANChanged published layouts violate only the boundary and become feasible after radial projection "
            "(\\NFIEAChanged); all but one violate it by less than 1 cm; the two layouts that violate the 2D spacing "
            "(participants 5 and 7, 36 turbines) stay infeasible.",
     lambda s, P, I: sorted((c["participant"], c["turbines"]) for c in I["changed_status"])
     == sorted([("par8", 16), ("par11", 16), ("par12", 16), ("par8", 36), ("par12", 36)])
     and sum(c["max_excess_m"] >= 0.01 for c in I["changed_status"]) == 1
     and sorted((c["participant"], c["turbines"]) for c in I["still_infeasible"]) == [("par5", 36), ("par7", 36)]),
    ("F42", "With projection the best feasible published layout is that of participant 12 in both scenarios (\\NFIEABestProjSixteen "
            "and \\NFIEABestProjThirtySix MWh; minimum spacing \\NFIEAProjSpacingSixteen and \\NFIEAProjSpacingThirtySix m), "
            "instead of participant 4 (\\NFIEABestStrictSixteen, \\NFIEABestStrictThirtySix).",
     lambda s, P, I: all(I["scenarios"][n]["projected"]["best_by"] == "par12" and I["scenarios"][n]["strict"]["best_by"] == "par4"
                         and I["scenarios"][n]["projected"]["best_min_spacing_m"] >= 260 for n in ("16", "36"))),
    ("F43", "The strict-convention gaps reproduce the paper's \\NIEA... macros (best and mean PSO-VNS gap, rank, at both budgets and "
            "both farm sizes).",
     lambda s, P, I: all(abs(I["scenarios"][n]["ours"][str(b)]["strict"][k] - paper_num(f"NIEA{m}{nm}{bn}")) < 0.006
                         for n, nm in (("16", "Sixteen"), ("36", "ThirtySix")) for b, bn in ((6030, "SixK"), (30030, "ThirtyK"))
                         for k, m in (("best_gap_pct", "Gap"), ("mean_gap_pct", "MeanGap")) if paper_num(f"NIEA{m}{nm}{bn}") is not None)),
    ("F44", "With projected published layouts the gap of the best PSO-VNS layout at 30,030 evaluations widens from "
            "\\NFIEAGapSixteenThirtyKStrict % to \\NFIEAGapSixteenThirtyKProj % (16 turbines) and from \\NFIEAGapThirtySixThirtyKStrict % "
            "to \\NFIEAGapThirtySixThirtyKProj % (36 turbines), i.e. \\NFIEAGapLossSixteenThirtyKProj and "
            "\\NFIEAGapLossThirtySixThirtyKProj pp more wake loss; it ranks \\NFIEARankSixteenThirtyKProj and "
            "\\NFIEARankThirtySixThirtyKProj among the published layouts.",
     lambda s, P, I: all(I["scenarios"][n]["ours"]["30030"]["projected"]["best_gap_pct"] < I["scenarios"][n]["ours"]["30030"]["strict"]["best_gap_pct"] < 0
                         and I["scenarios"][n]["ours"]["30030"]["projected"]["best_gap_loss_pp"] > 0 for n in ("16", "36"))
     and I["scenarios"]["16"]["ours"]["30030"]["projected"]["best_rank"] == 4 and I["scenarios"]["36"]["ours"]["30030"]["projected"]["best_rank"] == 9),
    ("F45", "The best of all our layouts is \\NFIEAOurBestGapThirtySixThirtyKProj % below the best feasible published layout at 36 "
            "turbines (\\NFIEAOurBestByThirtySixThirtyK, 30,030) and \\NFIEAOurBestGapSixteenThirtyKProj % at 16 turbines "
            "(\\NFIEAOurBestBySixteenThirtyK).",
     lambda s, P, I: I["scenarios"]["36"]["ours"]["30030"]["our_best_method"] == "SLSQP"
     and I["scenarios"]["16"]["ours"]["30030"]["our_best_method"] == "PSOBV"
     and I["scenarios"]["36"]["ours"]["30030"]["projected"]["our_best_gap_pct"] < 0),
    ("F46", "In wake-loss terms, the best PSO-VNS layout at 30,030 has \\NFIEAPSOVNSLossSixteenThirtyK % (16 turbines) and "
            "\\NFIEAPSOVNSLossThirtySixThirtyK % (36) against \\NFIEABestProjLossSixteen % and \\NFIEABestProjLossThirtySix % for "
            "the best feasible (projected) published layouts.",
     lambda s, P, I: all(I["scenarios"][n]["ours"]["30030"]["psovns_best_loss_pct"] > I["scenarios"][n]["projected"]["best_loss_pct"]
                         for n in ("16", "36"))),
]


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--summary", default=os.path.join(HERE, "mpce_summary_dir.json"))
    args = ap.parse_args(argv)
    s = json.load(open(args.summary))
    pf = os.path.join(HERE, "pywake_check.csv")
    P = pd.read_csv(pf) if os.path.exists(pf) else None
    I = json.load(open(os.path.join(HERE, "iea37_projected.json")))
    macros = open(os.path.join(HERE, "mpce_numbers_dir.tex")).read()
    nfail = 0
    for cid, text, fn in CHECKS:
        try:
            ok = bool(fn(s, P, I))
        except Exception as e:                     # missing data = FAIL
            ok, text = False, text + f"  [error: {type(e).__name__}: {e}]"
        used = [m for m in re.findall(r"\\(NF[A-Za-z]+)", text) if f"\\newcommand{{\\{m}}}" not in macros]
        if used:
            ok, text = False, text + f"  [undefined macros: {', '.join(used)}]"
        nfail += not ok
        print(f"[{cid}] {'PASS' if ok else 'FAIL'}  {text}")
    print(f"\n{len(CHECKS) - nfail}/{len(CHECKS)} checks passed")
    sys.exit(1 if nfail else 0)


if __name__ == "__main__":
    main()

"""Write analysis/mpce_numbers.tex: one \\newcommand{\\N...}{...} per number (or data-dependent phrase) quoted in
MPCE_PSO_VNS.tex, computed from mpce_summary.json (written by mpce_results.py).

Usage (from analysis/):  python3 mpce_numbers.py [--summary mpce_summary.json] [--out mpce_numbers.tex] [--allow-partial]

mpce_results.py calls this script at the end of every run, so normally nothing has to be run by hand.

Pending data. Every macro lists the experiments (mpce_<exp>_s<i>of<k>.csv) it depends on. If one of them is not
complete (missing, or only some shards present), the macro expands to \\TBD{pending: <what>} so that the paper
still compiles; --allow-partial prints the values computed from incomplete shards instead (preview only).
The main comparison uses the old-platform MS-SLSQP runs until mpce_slsqp exists (flagged by mpce_results.py);
these macros are NOT marked pending, because the SLSQP rerun only replaces one baseline.

Macro names contain letters only. Method codes: PSOVNS (PSOBV), PSO (PSOC, constriction coefficients),
SSAVNS (SSABV), SSA, LXSSA, DE, VNS (BVNS), SLSQP (MS-SLSQP), LXSSAVNS (LXBV), RSVNS, RSDVNS (RSD-VNS, the
disc-sampling control, experiment rsdisc). Split settings: TwentyFive / Fifty / SeventyFive / Ninety (PSOBV90,
experiment omega90) / Hundred (PSO alone).
"""
import os, json, argparse, math

HERE = os.path.dirname(os.path.abspath(__file__))
CODE = {"PSOBV": "PSOVNS", "PSOC": "PSO", "SSABV": "SSAVNS", "SSA": "SSA", "LXSSA": "LXSSA", "DE": "DE",
        "BVNS": "VNS", "SLSQP": "SLSQP", "LXBV": "LXSSAVNS", "RSVNS": "RSVNS", "RSDVNS": "RSDVNS"}
LAB = {"PSOBV": "PSO-VNS", "PSOC": "PSO", "SSABV": "SSA-VNS", "SSA": "SSA", "LXSSA": "LX-SSA", "DE": "DE",
       "BVNS": "VNS", "SLSQP": "MS-SLSQP", "LXBV": "LX-SSA-VNS", "RSVNS": "RS-VNS", "RSDVNS": "RSD-VNS"}
MAIN8 = ["PSOBV", "PSOC", "SSABV", "SSA", "LXSSA", "DE", "BVNS", "SLSQP"]
M10 = ["PSOBV", "SSABV", "LXBV", "RSVNS", "BVNS", "PSOC", "SSA", "LXSSA", "DE", "SLSQP"]
ORD = {1: "first", 2: "second", 3: "third", 4: "fourth", 5: "fifth", 6: "sixth", 7: "seventh", 8: "eighth",
       9: "ninth", 10: "tenth"}
WORD = {0: "no", 1: "one", 2: "two", 3: "three", 4: "four", 5: "five", 6: "six", 7: "seven", 8: "eight",
        9: "nine", 10: "ten", 11: "eleven", 12: "twelve", 13: "thirteen", 14: "fourteen", 15: "fifteen", 16: "sixteen"}
BUDGETS = (6030, 30030, 120030)
BNAME = {6030: "SixK", 30030: "ThirtyK", 120030: "OneTwentyK"}
BTXT = {6030: "6{,}030", 30030: "30{,}030", 120030: "120{,}030"}

# experiments each group of macros needs (see the module docstring)
REQ_MAIN = ("psobv", "psoc")
REQ_ABL = ("psobv", "psoc", "rsvns", "rsdisc")    # RSD-VNS is a variant of the component analysis (ranks, Holm family)
REQ_SPLIT = ("psosplit", "omega90")                # omega = 0.9 is a setting of the split table (ranks, Holm family)
REQ_HR = ("psobv", "hr16new")
REQ_FEAS = ("feas", "feasp")
REQ_B30 = ("b30k", "b30kp")
REQ_B120 = ("b120k", "b120kp")
REQ_IEA = ("iea16", "iea36", "iea16p", "iea36p")


# ------------------------------------------------------------------ formatting
def num(v, d=2, sign=False):
    if v is None or (isinstance(v, float) and not math.isfinite(v)):
        raise ValueError("missing value")
    s = f"{v:+,.{d}f}" if sign else f"{v:,.{d}f}"
    s = s.replace(",", "{,}")
    return s.replace("-", "\\ensuremath{-}", 1) if s.startswith("-") else s      # works in and outside math


def numd(v, d=2):
    """like num, but '--' when the value does not exist (e.g. a method without feasible runs)."""
    return "--" if v is None else num(v, d)


def pval(p):
    """p-value for print: two significant digits; scientific notation below 0.001, e.g. 1.9x10^-11 (\\ensuremath, so
    the macro works in text and in math)."""
    if p is None or not math.isfinite(p):
        raise ValueError("missing p")
    if p < 1e-3:
        m, e = f"{p:.1e}".split("e")
        return f"\\ensuremath{{{m}\\times10^{{{int(e)}}}}}"
    if p >= 0.995:
        return "1.0"
    return f"{p:.2g}" if p < 0.1 else f"{p:.2f}"


def pval_up(p):
    """like pval, but rounded UP (for bounds such as "p_Holm <= x"): the printed value is never
    smaller than p."""
    if p is None or not math.isfinite(p):
        raise ValueError("missing p")
    if p >= 0.995:
        return "1.0"
    e = math.floor(math.log10(p))
    if p < 1e-3:
        m = math.ceil(p / 10 ** e * 10 - 1e-9) / 10
        if m >= 10:
            m, e = 1.0, e + 1
        return f"\\ensuremath{{{m:.1f}\\times10^{{{e}}}}}"
    step = 10 ** (e - 1) if p < 0.1 else 0.01            # two significant digits (0.1 <= p: two decimals)
    v = math.ceil(p / step - 1e-9) * step
    d = max(0, -(e - 1)) if p < 0.1 else 2
    return f"{v:.{d}f}"


def prob(v):
    """posterior probability for print (the \\NBay... macros): two decimals, but ">0.99" above 0.99 and "<0.01" below
    0.01 (never 1.00 / 0.00); \\ensuremath, so the macro works in text and in math."""
    if v is None or not math.isfinite(v):
        raise ValueError("missing probability")
    if v > 0.99:
        return "\\ensuremath{>0.99}"
    if v < 0.01:
        return "\\ensuremath{<0.01}"
    return num(v, 2)


def ceil_num(v, d=2):
    """number rounded UP to d decimals (for bounds such as "within x")."""
    return num(math.ceil(v * 10 ** d - 1e-9) / 10 ** d, d)


def listing(items):
    items = list(items)
    if len(items) <= 1:
        return "".join(items)
    return ", ".join(items[:-1]) + " and " + items[-1]


def wtl(d):
    d = d or {}
    return int(d.get("W", 0)), int(d.get("T", 0)), int(d.get("L", 0))


class Macros:
    def __init__(self, summary, allow_partial):
        self.s = summary
        self.allow = allow_partial
        self.out, self.pending, self.errors = [], [], []
        av = summary.get("data_availability", {})
        self.status = {k.replace("mpce_", ""): v.get("status", "missing") for k, v in av.items() if k.startswith("mpce_")}

    def missing(self, req):
        bad = []
        for e in req:
            st = self.status.get(e, "missing")
            if st == "complete" or (self.allow and st.startswith("PARTIAL")):
                continue
            bad.append(e)
        return bad

    def put(self, name, fn, req=(), what=None):
        assert name.isalpha(), name
        bad = self.missing(req)
        val = None
        if not bad:
            try:
                val = fn()
            except (KeyError, TypeError, ValueError, IndexError, ZeroDivisionError, AttributeError) as e:
                self.errors.append(f"{name}: {type(e).__name__} {e}")
                bad = ["value missing in summary"]
        if val is None:
            txt = f"\\TBD{{pending: {', '.join(bad)}}}"
            self.pending.append(name)
        else:
            txt = str(val)
        self.out.append(f"\\newcommand{{\\{name}}}{{{txt}}}")


def build(s, allow_partial=False):
    M = Macros(s, allow_partial)
    P = M.put
    # Horns Rev: once any mpce_hrfix shard exists, every Horns Rev number comes from it (mpce_results.py drops
    # the old-model runs) and the Horns Rev macros are pending until all hrfix shards are present
    hrfix = M.status.get("hrfix", "missing") != "missing"
    REQ_HRX = ("hrfix",) if hrfix else REQ_HR
    hr_req = lambda req: REQ_HRX if hrfix else req + REQ_HR
    main = s.get("main") or {}
    fr = main.get("friedman") or {}
    cw = main.get("case_mean_wilcoxon") or {}
    focus = s.get("focus", "PSOBV")
    others = [a for a in MAIN8 if a != focus]

    # ---------------- main comparison (68 cases)
    P("NRunsMain", lambda: num(main["n_runs"], 0), REQ_MAIN, "runs of the main comparison")
    P("NRunsBenchmark", lambda: num(s["ablation"]["n_runs_benchmark_total"], 0), REQ_ABL, "runs on the 68 cases")
    P("NFriedChi", lambda: num(fr["chi2"], 1), REQ_MAIN)
    P("NFriedDf", lambda: str(len(fr["avg_rank"]) - 1), REQ_MAIN)
    P("NFriedP", lambda: pval(fr["p"]), REQ_MAIN)
    P("NImanF", lambda: num(fr["iman_davenport"], 1), REQ_MAIN)
    P("NFriedVerb", lambda: "rejects" if fr["p"] < 0.05 else "does not reject", REQ_MAIN)
    for a in MAIN8:
        c = CODE[a]
        P(f"NRank{c}", lambda a=a: num(fr["avg_rank"][a], 2), REQ_MAIN)
        P(f"NSoleBest{c}", lambda a=a: str(fr["sole_best_count"][a]), REQ_MAIN)
        P(f"NLoss{c}", lambda a=a: num(main["mean_loss_pct_qualified_cases"][a], 3), REQ_MAIN)
        P(f"NFeas{c}", lambda a=a: num(main["feasible_pct"][a], 1), REQ_MAIN)
    order = sorted(fr.get("avg_rank", {}), key=lambda a: fr["avg_rank"][a])
    P("NRankFirst", lambda: LAB[order[0]], REQ_MAIN)
    P("NRankOrderOthers", lambda: listing(f"{LAB[a]} ({num(fr['avg_rank'][a], 2)})" for a in order[1:]), REQ_MAIN)
    ph = fr.get("p_holm_vs_focus", {})
    sig = [a for a in order if a in ph and ph[a] < 0.05]
    ns = [a for a in order if a in ph and ph[a] >= 0.05]
    P("NPostHocSigList", lambda: listing(LAB[a] for a in sig) if sig else "none of the other methods", REQ_MAIN)
    P("NPostHocSigMaxP", lambda: pval_up(max(ph[a] for a in sig)), REQ_MAIN)       # used as "p_Holm <= x"
    P("NPostHocNonSigList", lambda: listing(LAB[a] for a in ns) if ns else "any method", REQ_MAIN)
    P("NPostHocNonSigP", lambda: pval(min(ph[a] for a in ns)) if ns else "--", REQ_MAIN)
    for b in others:
        c = CODE[b]
        W, T, L = wtl(main["wtl"].get(b)) if main.get("wtl") else (None,) * 3
        P(f"NWtl{c}", lambda W=W, T=T, L=L: f"{W}/{T}/{L}", REQ_MAIN)
        P(f"NWtl{c}W", lambda W=W: str(W), REQ_MAIN)
        P(f"NWtl{c}L", lambda L=L: str(L), REQ_MAIN)
        P(f"NWtl{c}T", lambda T=T: str(T), REQ_MAIN)
        P(f"NRrb{c}", lambda b=b: num(main["rank_biserial_median"][b], 2, sign=True), REQ_MAIN)
        P(f"NPW{c}", lambda b=b: pval(cw[b]["p_holm"]), REQ_MAIN)            # Holm over the 7 methods
        P(f"NPWraw{c}", lambda b=b: pval(cw[b]["p"]), REQ_MAIN)             # unadjusted
        P(f"NDL{c}", lambda b=b: num(cw[b]["mean_dloss_pp"], 3, sign=True), REQ_MAIN)
        P(f"NDLabs{c}", lambda b=b: num(-cw[b]["mean_dloss_pp"], 2), REQ_MAIN)
        P(f"NLowerLoss{c}", lambda b=b: str(cw[b]["focus_lower_loss_cases"]), REQ_MAIN)
        P(f"NHigherLoss{c}", lambda b=b: str(cw[b]["other_lower_loss_cases"]), REQ_MAIN)
    P("NWtlPSOSig", lambda: str(sum(wtl(main["wtl"]["PSOC"])[::2])), REQ_MAIN)
    byn = main.get("by_n", {})
    P("NNLarge", lambda: str(byn["PSOC"]["n_large"]), REQ_MAIN)
    P("NWtlPSOLargeW", lambda: str(wtl(byn["PSOC"]["wtl_large"])[0]), REQ_MAIN)
    P("NWtlPSOLargeL", lambda: str(wtl(byn["PSOC"]["wtl_large"])[2]), REQ_MAIN)
    P("NWtlPSOSmallW", lambda: str(wtl(byn["PSOC"]["wtl_small"])[0]), REQ_MAIN)
    P("NWtlPSOSmallL", lambda: str(wtl(byn["PSOC"]["wtl_small"])[2]), REQ_MAIN)
    P("NGainPSOLarge", lambda: num(-byn["PSOC"]["mean_dloss_pp_large"], 2), REQ_MAIN)
    P("NAbsDiffPSOSmall", lambda: num(byn["PSOC"]["mean_abs_dloss_pp_small"], 2), REQ_MAIN)
    # losses of PSO-VNS against the other methods (phrase)
    def losses_phrase():
        lost = [(b, wtl(main["wtl"][b])[2]) for b in others if b != "PSOC" and wtl(main["wtl"][b])[2] > 0]
        never = [b for b in others if b != "PSOC" and wtl(main["wtl"][b])[2] == 0]
        parts = []
        if never:
            parts.append(f"it is never significantly worse than {listing(LAB[b] for b in never)}")
        if lost:
            parts.append("significantly worse than " + listing(f"{LAB[b]} in {WORD.get(n, n)} case{'s' if n > 1 else ''}" for b, n in lost))
        return " and ".join(parts)
    P("NLossesPhrase", losses_phrase, REQ_MAIN)
    P("NLossMaxSmallN", lambda: num(main["max_loss_n_le_3_pct"], 2), REQ_MAIN)
    big = (main.get("loss_largest_n") or {}).get("1-1000-15", {})
    for a in ("PSOBV", "PSOC", "SSABV", "SSA"):
        P(f"NLossBigI{CODE[a]}", lambda a=a: num(big[a], 2), REQ_MAIN, "loss DS I, 1000 m, N = 15")
    g = main.get("largest_n_loss_reduction_vs_phase1_pp") or {}
    P("NBigGainMin", lambda: num(g["min"], 2, sign=True), REQ_MAIN)
    P("NBigGainMax", lambda: num(g["max"], 2, sign=True), REQ_MAIN)
    P("NBigGainPos", lambda: WORD[sum(1 for v in main["loss_largest_n"].values() if v["PSOC"] is not None and v["PSOBV"] is not None and v["PSOC"] > v["PSOBV"])], REQ_MAIN)
    dn = main.get("dense_500_10") or {}
    P("NDenseSwitchFeas", lambda: num(sum(v["feasible_at_switch_pct"] for v in dn.values()) / len(dn), 1), REQ_MAIN)
    P("NDenseFinalFeas", lambda: num(sum(v["feasible_final_pct"] for v in dn.values()) / len(dn), 1), REQ_MAIN)
    P("NDensePSOFinalFeas", lambda: num(sum(v["phase1_alone_final_feasible_pct"] for v in dn.values()) / len(dn), 1), REQ_MAIN)
    fbd = main.get("feasible_runs_by_dataset") or {}

    # ---------------- phase 2 / convergence
    p2 = (s.get("ablation") or {}).get("phase2_loss_reduction_pct") or {}
    P("NPhTwoMean", lambda: num(p2["PSOBV"]["mean"], 1), REQ_ABL)
    P("NPhTwoMedian", lambda: num(p2["PSOBV"]["median"], 1), REQ_ABL)
    P("NPhTwoMadeFeas", lambda: str(p2["PSOBV"]["infeasible_at_switch_made_feasible"]), REQ_ABL)
    P("NPsoContMean", lambda: num(p2["PSOC_continued"]["mean"], 1), REQ_ABL)
    P("NPsoContMedian", lambda: num(p2["PSOC_continued"]["median"], 1), REQ_ABL)
    for a in ("PSOBV", "SSABV", "LXBV", "RSVNS", "RSDVNS"):
        P(f"NPhTwoMean{CODE[a]}", lambda a=a: num(p2[a]["mean"], 1), REQ_ABL)
        P(f"NSwitchFeas{CODE[a]}", lambda a=a: num(p2[a]["feasible_at_switch_pct"], 1), REQ_ABL)
        P(f"NSwitchLoss{CODE[a]}", lambda a=a: num(p2[a]["mean_loss_at_switch_pct"], 2), REQ_ABL)

    # ---------------- earlier PSO setting
    ps = s.get("pso_setting") or {}
    P("NOldPSORank", lambda: num(ps["old_avg_rank"], 2), REQ_MAIN)
    P("NOldPSOPos", lambda: ORD[ps["old_rank_position"]], REQ_MAIN)
    P("NPSOPos", lambda: ORD[ps["constriction_rank_position"]], REQ_MAIN)
    P("NOldPSOFeas", lambda: num(ps["old_feasible_pct"], 1), REQ_MAIN)
    P("NOldPSOWorse", lambda: str(wtl(ps["constriction_vs_old_wtl"])[0]), REQ_MAIN)
    P("NOldPSOBetter", lambda: WORD.get(wtl(ps["constriction_vs_old_wtl"])[2], str(wtl(ps["constriction_vs_old_wtl"])[2])), REQ_MAIN)
    P("NOldPSODL", lambda: num(ps["old_minus_constriction_loss_pp"], 2), REQ_MAIN)
    P("NOldPSOCases", lambda: str(ps["n_cases_loss"]), REQ_MAIN)
    P("NOldPSOBelowHalf", lambda: str(ps["old_cases_below_half_feasible"]), REQ_MAIN)

    # ---------------- previous study's method pool with the old / constriction PSO setting (tab:baseline)
    bs = s.get("baseline") or {}
    for tag, nm in (("old", "Old"), ("constriction", "New")):
        P(f"NBase{nm}Best", lambda tag=tag: bs[tag]["best_label"], REQ_MAIN)
        P(f"NBase{nm}PSORank", lambda tag=tag: num(bs[tag]["avg_rank"]["PSO"], 2), REQ_MAIN)
        P(f"NBase{nm}LXBVRank", lambda tag=tag: num(bs[tag]["avg_rank"]["LXBV"], 2), REQ_MAIN)
        P(f"NBase{nm}SSAVNSRank", lambda tag=tag: num(bs[tag]["avg_rank"]["SSABV"], 2), REQ_MAIN)
        P(f"NBase{nm}PSOPos", lambda tag=tag: str(int(bs[tag]["pso_position"])), REQ_MAIN)
        P(f"NBase{nm}LXBVvsPSO", lambda tag=tag: "%d/%d/%d" % wtl(bs[tag]["wtl_vs_pso"]["LXBV"]), REQ_MAIN)
        P(f"NBase{nm}SSAVNSvsPSO", lambda tag=tag: "%d/%d/%d" % wtl(bs[tag]["wtl_vs_pso"]["SSABV"]), REQ_MAIN)

    # ---------------- ablation
    ab = s.get("ablation") or {}
    fa = ab.get("friedman") or {}
    P("NAblChi", lambda: num(fa["chi2"], 1), REQ_ABL)
    P("NAblP", lambda: pval(fa["p"]), REQ_ABL)
    P("NAblNVar", lambda: WORD[len(fa["avg_rank"])], REQ_ABL)
    P("NAblOrder", lambda: listing(f"{LAB[a]} ({num(v, 2)})" for a, v in sorted(fa["avg_rank"].items(), key=lambda t: t[1])), REQ_ABL)
    ct = ab.get("contrasts") or {}
    for key in ("PSOBV-PSOC", "PSOBV-BVNS", "PSOBV-RSVNS", "PSOBV-SSABV", "PSOBV-LXBV", "SSABV-RSVNS", "SSABV-LXBV", "LXSSA-SSA",
                "SSABV-SSA", "LXBV-LXSSA", "LXBV-RSVNS", "PSOBV-RSDVNS", "SSABV-RSDVNS", "LXBV-RSDVNS", "RSDVNS-RSVNS"):
        a, b = key.split("-")
        nm = f"NAbl{CODE[a]}vs{CODE[b]}"
        P(nm, lambda key=key: "%d/%d/%d" % wtl(ct[key]), REQ_ABL)
        P(nm + "W", lambda key=key: str(wtl(ct[key])[0]), REQ_ABL)
        P(nm + "T", lambda key=key: str(wtl(ct[key])[1]), REQ_ABL)
        P(nm + "L", lambda key=key: str(wtl(ct[key])[2]), REQ_ABL)
        P(nm + "DL", lambda key=key: num(ct[key]["dloss_pp"], 3, sign=True), REQ_ABL)

    # ---------------- budget split
    sp = s.get("split") or {}
    P("NSplitCases", lambda: WORD.get(sp["n_cases"], str(sp["n_cases"])), REQ_SPLIT)
    for w, k in ((25, "PSOBV25"), (50, "PSOBV"), (75, "PSOBV75")):
        nm = {25: "TwentyFive", 50: "Fifty", 75: "SeventyFive"}[w]
        P(f"NSplitLoss{nm}", lambda k=k: num(sp["mean_loss"][k], 3), REQ_SPLIT)
        P(f"NSplitRank{nm}", lambda k=k: num(sp["avg_rank"][k], 2), REQ_SPLIT)
        P(f"NSplitBest{nm}", lambda k=k: WORD.get(sp["best_count"].get(k, 0), str(sp["best_count"].get(k, 0))), REQ_SPLIT)
    P("NSplitSeventyFiveLower", lambda: WORD.get(sp["n_lower_loss_than_50"]["PSOBV75"], str(sp["n_lower_loss_than_50"]["PSOBV75"])), REQ_SPLIT)
    P("NSplitSeventyFiveSig", lambda: WORD.get(sp["wtl_50_vs_75"]["L"], str(sp["wtl_50_vs_75"]["L"])), REQ_SPLIT)
    P("NSplitFiftyBetterSeventyFive", lambda: WORD.get(sp["wtl_50_vs_75"]["W"], str(sp["wtl_50_vs_75"]["W"])), REQ_SPLIT)
    P("NSplitTwentyFiveSig", lambda: WORD.get(sp["wtl_50_vs_25"]["W"], str(sp["wtl_50_vs_25"]["W"])), REQ_SPLIT)
    P("NSplitTwentyFiveBetter", lambda: WORD.get(sp["wtl_50_vs_25"]["L"], str(sp["wtl_50_vs_25"]["L"])), REQ_SPLIT)
    P("NSplitGainSeventyFive", lambda: num(sp["mean_loss"]["PSOBV"] - sp["mean_loss"]["PSOBV75"], 2), REQ_SPLIT)

    # ---------------- Horns Rev 16
    hr = s.get("hr16") or {}
    hm = hr.get("methods") or {}
    P("NHRInstalled", lambda: num(hr["installed_aep"], 2), REQ_HRX)
    P("NHRInstalledLoss", lambda: num(hr["installed_loss_pct"], 2), REQ_HRX)
    P("NHRIdeal", lambda: num(hr["ideal_aep"], 2), REQ_HRX)
    for a in M10:
        c = CODE[a]
        P(f"NHRMean{c}", lambda a=a: numd(hm[a]["mean"], 2), REQ_HRX)
        P(f"NHRBest{c}", lambda a=a: numd(hm[a]["best"], 2), REQ_HRX)
        P(f"NHRFeas{c}", lambda a=a: f"{hm[a]['feasible']}/{hm[a]['runs']}", REQ_HRX)
        P(f"NHRAbove{c}", lambda a=a: str(hm[a]["runs_above_installed"]), REQ_HRX)
    P("NHRSDPSOVNS", lambda: num(hm["PSOBV"]["sd"], 2), REQ_HRX)
    P("NHRLossSixKPSOVNS", lambda: num(hm["PSOBV"]["loss_pct"], 2), REQ_HRX)

    def above_phrase():
        o = [a for a in M10 if a != focus and a in hm and hm[a]["runs_above_installed"] > 0]
        if not o:
            return "no other method reaches it"
        runs = lambda k: f"{WORD.get(k, k)} run{'s' if k != 1 else ''}"
        return ("of the other methods, only " + listing(f"{LAB[a]} ({runs(hm[a]['runs_above_installed'])})" for a in o)
                + (" does so" if len(o) == 1 else " do so"))
    P("NHRAboveOthersPhrase", above_phrase, REQ_HRX)
    P("NHRMaxP", lambda: pval_up(max(v["p_holm"] for a, v in hm.items() if a != focus and v.get("p_holm") is not None)), REQ_HRX)   # "p_Holm <= x"
    P("NHRAboveOthers", lambda: (lambda o: "no other method" if not o else listing(f"{LAB[a]} ({hm[a]['runs_above_installed']})" for a in o))(
        [a for a in M10 if a != focus and a in hm and hm[a]["runs_above_installed"] > 0]), REQ_HRX)
    lb = hr.get("loss_by_setting") or {}
    for tag, nm, req in (("6030F", "FeasInit", REQ_FEAS), ("30030R", "ThirtyK", REQ_B30), ("120030R", "OneTwentyK", REQ_B120)):
        P(f"NHRLoss{nm}PSOVNS", lambda tag=tag: num(lb["PSOBV"][tag]["loss"], 2), hr_req(req))
        P(f"NHRMean{nm}PSOVNS", lambda tag=tag: num(lb["PSOBV"][tag]["mean_aep"], 2), hr_req(req))
        P(f"NHRAbove{nm}PSOVNS", lambda tag=tag: str(lb["PSOBV"][tag]["runs_above_installed"]), hr_req(req))
        P(f"NHRFeas{nm}PSO", lambda tag=tag: f"{lb['PSOC'][tag]['feasible']}/{lb['PSOC'][tag]['runs']}", hr_req(req))
    # Horns Rev 1 block at 30,030 evaluations (PSO-VNS): mean AEP and runs above the installed layout ('x of y')
    P("NHRThirtyKMean", lambda: num(lb["PSOBV"]["30030R"]["mean_aep"], 2), hr_req(REQ_B30))
    P("NHRThirtyKAbove", lambda: f"{lb['PSOBV']['30030R']['runs_above_installed']} of {lb['PSOBV']['30030R']['runs']}", hr_req(REQ_B30))

    # ---------------- feasible initialization and budget (six largest benchmark cases)
    fb = s.get("feasbudget") or {}
    ran, fea = fb.get("random") or {}, fb.get("feasible") or {}
    for a in ("PSOBV", "PSOC", "SSABV", "SSA", "LXSSA", "DE", "RSVNS"):
        c = CODE[a]
        P(f"NFbRandFeas{c}", lambda a=a: num(ran[a]["feas"], 1), REQ_MAIN)
        P(f"NFbRandLoss{c}", lambda a=a: numd(ran[a]["loss"], 2), REQ_MAIN)
        P(f"NFbFeasFeas{c}", lambda a=a: "n/r" if a == "RSVNS" and not fea.get(a) else num(fea[a]["feas"], 1), REQ_FEAS)
        P(f"NFbFeasLoss{c}", lambda a=a: "n/r" if a == "RSVNS" and not fea.get(a) else numd(fea[a]["loss"], 2), REQ_FEAS)
    P("NFbFeasBestLoss", lambda: LAB[min((a for a in fea if fea[a] and fea[a]["loss"] is not None), key=lambda a: fea[a]["loss"])], REQ_FEAS)
    P("NFbFeasBestRank", lambda: LAB[min(fb["rank_feasible_init"], key=fb["rank_feasible_init"].get)], REQ_FEAS)
    P("NFbFeasRankPSOVNS", lambda: num(fb["rank_feasible_init"]["PSOBV"], 2), REQ_FEAS)
    rk = fb.get("rank") or {}
    for b in BUDGETS:
        req = {6030: REQ_MAIN, 30030: REQ_B30, 120030: REQ_B120}[b]
        for a in ("PSOBV", "PSOC", "BVNS", "RSVNS", "SLSQP", "SSABV"):
            P(f"NBudRank{CODE[a]}{BNAME[b]}", lambda a=a, b=b: num(rk[str(b)][a], 2), req)
            P(f"NBudLoss{CODE[a]}{BNAME[b]}", lambda a=a, b=b: num(fb["mean_loss"][str(b)][a], 2), req)
        P(f"NBudBest{BNAME[b]}", lambda b=b: LAB[min(rk[str(b)], key=rk[str(b)].get)], req)
        P(f"NBudGainPSO{BNAME[b]}", lambda b=b: num(fb["mean_loss"][str(b)]["PSOC"] - fb["mean_loss"][str(b)]["PSOBV"], 2, sign=True), req)

    def budget_phrase():
        best = {b: min(rk[str(b)], key=rk[str(b)].get) for b in BUDGETS}
        if all(best[b] == focus for b in BUDGETS):
            return f"{LAB[focus]} keeps the best average rank at all three budgets"
        lost = [b for b in BUDGETS if best[b] != focus]
        return (f"{LAB[focus]} has the best average rank at " + listing(BTXT[b] for b in BUDGETS if best[b] == focus) +
                " calls, and " + listing(f"{LAB[best[b]]} at {BTXT[b]}" for b in lost) + " calls")
    P("NBudgetPhrase", budget_phrase, REQ_MAIN + REQ_B30 + REQ_B120, "budget ranking")

    def feas_phrase():
        # which methods gain feasibility, and does the best method change
        gain = [a for a in M10 if ran.get(a) and fea.get(a) and fea[a]["feas"] - ran[a]["feas"] >= 5]
        txt = ("the share of feasible runs rises for " + listing(f"{LAB[a]} ({num(ran[a]['feas'], 1)}\\% to {num(fea[a]['feas'], 1)}\\%)" for a in gain)
               if gain else "the share of feasible runs changes by less than five percentage points for every method")
        return txt
    P("NFeasInitPhrase", feas_phrase, REQ_MAIN + REQ_FEAS, "feasible-initialization summary")

    # ---------------- IEA37
    ic = (s.get("iea37_compact") or {}).get("methods", {})
    pub = (s.get("iea37_compact") or {}).get("published", {})
    for n, nm in ((16, "Sixteen"), (36, "ThirtySix")):
        req_n = (f"iea{n}", f"iea{n}p")
        P(f"NIEAPubFeas{nm}", lambda n=n: str(pub[str(n)]["n_feasible"]), ())
        for b in (6030, 30030):
            k = f"{n}T_{b}"
            P(f"NIEABest{nm}{BNAME[b]}", lambda k=k: num(ic["PSOBV"][k]["best"], 1), req_n)
            P(f"NIEAGap{nm}{BNAME[b]}", lambda k=k: num(ic["PSOBV"][k]["best_gap_to_best_feasible_published_pct"], 2, sign=True), req_n)
            P(f"NIEAMeanGap{nm}{BNAME[b]}", lambda k=k: num(ic["PSOBV"][k]["mean_gap_to_best_feasible_published_pct"], 2, sign=True), req_n)
            P(f"NIEARank{nm}{BNAME[b]}", lambda k=k: ORD[ic["PSOBV"][k]["best_rank_among_feasible_published"]], req_n)
            P(f"NIEAMeanRank{nm}{BNAME[b]}", lambda k=k: ORD[ic["PSOBV"][k]["mean_rank_among_feasible_published"]], req_n)
            P(f"NIEABestOf{nm}{BNAME[b]}", lambda k=k: LAB[max((a for a in ic if ic[a].get(k) and ic[a][k]["best"] is not None),
                                                                 key=lambda a: ic[a][k]["best"])], req_n)

    P("NIEAPhrase", lambda: "reach " + listing(
        f"{num(100 + ic['PSOBV'][f'{n}T_30030']['best_gap_to_best_feasible_published_pct'], 1)}\\% ({n} turbines)" for n in (16, 36))
        + " of the AEP of the best feasible published layouts", REQ_IEA, "IEA37 summary")

    # ---------------- cost
    co = s.get("cost") or {}
    P("NMsPerCallMin", lambda: num(min(co["ms_per_call"].values()), 2), REQ_MAIN)
    P("NMsPerCallMax", lambda: num(max(co["ms_per_call"].values()), 2), REQ_MAIN)
    P("NSecPerRunMin", lambda: num(co["sec_per_run_range"]["PSOBV"][0], 1), REQ_MAIN)
    P("NSecPerRunMax", lambda: num(co["sec_per_run_range"]["PSOBV"][1], 1), REQ_MAIN)

    # ---------------- robustness
    rb = s.get("robustness") or {}
    P("NCubicOverI", lambda: num(-rb["Cubic"]["rel_change_pct"]["1"], 1), REQ_MAIN)
    P("NCubicOverII", lambda: num(-rb["Cubic"]["rel_change_pct"]["2"], 1), REQ_MAIN)
    for col, nm in (("Cubic", "Cubic"), ("CubicCutout", "Cutout"), ("Gauss", "Gauss")):
        P(f"NTau{nm}", lambda col=col: num(rb[col]["tau"], 2), REQ_MAIN)
        P(f"NSame{nm}", lambda col=col: num(rb[col]["same_best_pct"], 0), REQ_MAIN)
        P(f"NRank{nm}PSOVNS", lambda col=col: num(rb[col]["avg_rank"][focus], 2), REQ_MAIN)
        P(f"NRank{nm}PSO", lambda col=col: num(rb[col]["avg_rank"]["PSOC"], 2), REQ_MAIN)
        P(f"NPos{nm}PSOVNS", lambda col=col: ORD[1 + sorted(rb[col]["avg_rank"].values()).index(rb[col]["avg_rank"][focus])], REQ_MAIN)
        P(f"NBest{nm}", lambda col=col: LAB[min(rb[col]["avg_rank"], key=rb[col]["avg_rank"].get)], REQ_MAIN)
    P("NGaussDeltaI", lambda: num(rb["Gauss"]["rel_change_pct"]["1"], 1, sign=True), REQ_MAIN)

    # ---------------- minimum spacing of the final layouts
    sh = main.get("spacing_share") or {}
    P("NSpFiveAll", lambda: num(sh["all"]["pct_5d"], 0), REQ_MAIN)
    P("NSpSixAll", lambda: num(sh["all"]["pct_6d"], 0), REQ_MAIN)
    P("NSpFiveLarge", lambda: num(sh["n11_15"]["pct_5d"], 1), REQ_MAIN)
    P("NSpSixLarge", lambda: num(sh["n11_15"]["pct_6d"], 1), REQ_MAIN)
    # ---------------- macros of the component analysis, budgets, Horns Rev 1 and IEA37
    # case-mean Wilcoxon of the component-analysis contrasts (summary ablation.case_mean; procedure of
    # case_mean_wilcoxon, the same as \NPW...): P unadjusted two-sided p (PHolm: Holm over the contrasts of
    # Table ablation), Wins / Losses = cases in which A has the lower / higher case-mean loss (cases in which
    # only one method has >= 15 feasible runs counted for that method), DL = mean of L(A) - L(B) over the
    # cases in which both qualify (pp, signed, negative = A better), CI = 95 % percentile bootstrap CI of DL
    # (10,000 resamples of the cases, fixed seed) as "[a, b]".
    cm = ab.get("case_mean") or {}
    CMP = {"SSAVNSvsRSVNS": "SSABV-RSVNS", "LXSSAVNSvsRSVNS": "LXBV-RSVNS", "PSOVNSvsRSVNS": "PSOBV-RSVNS",
           "SSAVNSvsLXSSAVNS": "SSABV-LXBV", "LXSSAvsSSA": "LXSSA-SSA", "SSAVNSvsSSA": "SSABV-SSA",
           "LXSSAVNSvsLXSSA": "LXBV-LXSSA", "PSOVNSvsPSO": "PSOBV-PSOC", "PSOVNSvsVNS": "PSOBV-BVNS",
           "SSAVNSvsRSDVNS": "SSABV-RSDVNS", "LXSSAVNSvsRSDVNS": "LXBV-RSDVNS", "PSOVNSvsRSDVNS": "PSOBV-RSDVNS",
           "RSDVNSvsRSVNS": "RSDVNS-RSVNS"}
    for nm, key in CMP.items():
        P(f"NCm{nm}P", lambda key=key: pval(cm[key]["p"]), REQ_ABL)
        P(f"NCm{nm}PHolm", lambda key=key: pval(cm[key]["p_holm"]), REQ_ABL)
        P(f"NCm{nm}Wins", lambda key=key: str(cm[key]["wins"]), REQ_ABL)
        P(f"NCm{nm}Losses", lambda key=key: str(cm[key]["losses"]), REQ_ABL)
        P(f"NCm{nm}DL", lambda key=key: num(cm[key]["mean_dloss_pp"], 3, sign=True), REQ_ABL)
        P(f"NCm{nm}CI", lambda key=key: "[%s, %s]" % tuple(num(v, 3) for v in cm[key]["ci95_mean_dloss_pp"]), REQ_ABL)
    P("NAblFriedVerb", lambda: "rejects" if fa["p"] < 0.05 else "does not reject", REQ_ABL)
    # practical equivalence of the case means (summary equivalence; mpce_results.equivalence_block): one margin
    # EQ_MARGIN (pp) for every pair, set after the primary analysis. Same case set as \NCm...DL / CI (both methods
    # >= 15 feasible runs). CI = 90 % percentile bootstrap CI of the mean L(A) - L(B) (10,000 case resamples,
    # fixed seed) as "[a, b]"; P = bootstrap TOST p (add-one; at its floor 1/10001 it is a bound, see
    # p_tost_at_floor); Min = smallest margin at which equivalence holds, rounded UP to 3 decimals; Holds =
    # "yes" iff the 90 % CI lies inside (-m, m). Bay...: Bayesian signed-rank test (Benavoli et al. 2017) with
    # ROPE (-m, m): Better = P(A better), Rope = P(practically equivalent), Worse = P(B better), 2 decimals (prob():
    # "\\ensuremath{>0.99}" above 0.99, "\\ensuremath{<0.01}" below 0.01).
    eq = s.get("equivalence") or {}
    eqp = eq.get("pairs") or {}
    P("NEqMargin", lambda: num(eq["margin_pp"], 2), REQ_ABL)
    for nm, key in (("PSOVNSvsPSO", "PSOBV-PSOC"), ("SSAVNSvsRSVNS", "SSABV-RSVNS"), ("LXSSAVNSvsRSVNS", "LXBV-RSVNS"),
                    ("PSOVNSvsRSVNS", "PSOBV-RSVNS"), ("SSAVNSvsLXSSAVNS", "SSABV-LXBV"),
                    ("SSAVNSvsRSDVNS", "SSABV-RSDVNS"), ("LXSSAVNSvsRSDVNS", "LXBV-RSDVNS"), ("PSOVNSvsRSDVNS", "PSOBV-RSDVNS"),
                    ("RSDVNSvsRSVNS", "RSDVNS-RSVNS")):
        P(f"NEq{nm}CI", lambda key=key: "[%s, %s]" % tuple(num(v, 3) for v in eqp[key]["ci90_mean_dloss_pp"]), REQ_ABL)
        P(f"NEq{nm}P", lambda key=key: pval(eqp[key]["p_tost"]), REQ_ABL)
        P(f"NEq{nm}Min", lambda key=key: ceil_num(eqp[key]["min_margin_pp"], 3), REQ_ABL)
        P(f"NEq{nm}Holds", lambda key=key: "yes" if eqp[key]["equivalent"] else "no", REQ_ABL)
        P(f"NBay{nm}Better", lambda key=key: prob(eqp[key]["bayes"]["p_a_better"]), REQ_ABL)
        P(f"NBay{nm}Rope", lambda key=key: prob(eqp[key]["bayes"]["p_rope"]), REQ_ABL)
        P(f"NBay{nm}Worse", lambda key=key: prob(eqp[key]["bayes"]["p_b_better"]), REQ_ABL)
    # spread and worst run, PSO-VNS vs PSO (summary spread): mean over the cases of the per-case SD of the
    # feasible-run wake loss (pp, 3 decimals); cases in which PSO-VNS has the smaller / larger SD and the better /
    # worse worst feasible run (ties excluded); Wilcoxon p on the per-case SDs / worst runs (unadjusted)
    spr = s.get("spread") or {}
    P("NSdPSOVNS", lambda: num(spr["mean_sd_pp"]["PSOBV"], 3), REQ_MAIN)
    P("NSdPSO", lambda: num(spr["mean_sd_pp"]["PSOC"], 3), REQ_MAIN)
    P("NSdWinsPSOVNS", lambda: str(spr["sd_wins"]), REQ_MAIN)
    P("NSdLossesPSOVNS", lambda: str(spr["sd_losses"]), REQ_MAIN)
    P("NSdP", lambda: pval(spr["sd_p"]), REQ_MAIN)
    P("NWorstWinsPSOVNS", lambda: str(spr["worst_wins"]), REQ_MAIN)
    P("NWorstLossesPSOVNS", lambda: str(spr["worst_losses"]), REQ_MAIN)
    P("NWorstP", lambda: pval(spr["worst_p"]), REQ_MAIN)
    # exploratory N >= 10 subgroups, PSO-VNS vs PSO (case-mean Wilcoxon unadjusted; run level: main-table family)
    sg = (main.get("subgroup_vs_phase1") or {}).get("groups") or {}
    for g_ in ("Large", "dsILarge", "dsIILarge"):
        P(f"NCmPSOVNSvsPSO{g_}P", lambda g_=g_: pval(sg[g_]["p"]), REQ_MAIN)
        P(f"NCmPSOVNSvsPSO{g_}Wins", lambda g_=g_: str(sg[g_]["wins"]), REQ_MAIN)
        P(f"NCmPSOVNSvsPSO{g_}Losses", lambda g_=g_: str(sg[g_]["losses"]), REQ_MAIN)
        P(f"NCmPSOVNSvsPSO{g_}N", lambda g_=g_: str(sg[g_]["n_cases"]), REQ_MAIN)
        P(f"NCmPSOVNSvsPSO{g_}DL", lambda g_=g_: num(sg[g_]["mean_dloss_pp"], 3, sign=True), REQ_MAIN)
    for g_ in ("dsILarge", "dsIILarge"):
        P(f"NWtlPSO{g_}W", lambda g_=g_: str(sg[g_]["run_level"]["W"]), REQ_MAIN)
        P(f"NWtlPSO{g_}L", lambda g_=g_: str(sg[g_]["run_level"]["L"]), REQ_MAIN)
    P("NImanP", lambda: pval(fr["iman_davenport_p"]), REQ_MAIN)
    ci_ = main.get("case_mean_imputation") or {}
    P("NCmImputed", lambda: str(ci_["imputed_total"]), REQ_MAIN)      # plain count, summed over the 7 comparisons
    P("NCmDropped", lambda: str(ci_["dropped_total"]), REQ_MAIN)
    P("NCmImputedPhrase", lambda: listing(f"{v} for {LAB[b]}" for b, v in ci_["imputed"].items() if v) or "none", REQ_MAIN)
    # time per evaluation (median over the 6,030-evaluation runs of the 68 cases; see mpce_results.py)
    ce = s.get("cost_per_eval") or {}
    P("NMsEvalMin", lambda: num(ce["meta_min_ms"], 2), REQ_MAIN)
    P("NMsEvalMax", lambda: num(ce["meta_max_ms"], 2), REQ_MAIN)
    P("NMsEvalSLSQP", lambda: num(ce["slsqp_ms"], 2), REQ_MAIN)
    P("NMsEvalPSO", lambda: num(ce["median_ms"]["PSOC"], 2), REQ_MAIN)
    P("NSLSQPSlowdown", lambda: num(ce["slsqp_over_pso"], 1), REQ_MAIN)       # plain number (text adds "times")
    # budget split: omega = 1 (PSO alone, same seeds and budget) as the end point
    P("NSplitRankHundred", lambda: num(sp["avg_rank"]["PSOC"], 2), REQ_SPLIT + ("psoc",))
    P("NSplitLossHundred", lambda: num(sp["mean_loss"]["PSOC"], 3), REQ_SPLIT + ("psoc",))
    P("NSplitHundredVsSeventyFive", lambda: "%d/%d/%d" % wtl(sp["wtl_75_vs_100"]), REQ_SPLIT + ("psoc",))
    P("NSplitCmP", lambda: pval(sp["case_mean"]["PSOBV-PSOBV75"]["p"]), REQ_SPLIT)     # case-mean 0.75 vs 0.5
    P("NSplitCmHundredP", lambda: pval(sp["case_mean"]["PSOBV75-PSOC"]["p"]), REQ_SPLIT + ("psoc",))
    P("NSplitCmHundredWins", lambda: str(sp["case_mean"]["PSOBV75-PSOC"]["wins"]), REQ_SPLIT + ("psoc",))
    # first places of PSO-VNS among the six largest cases per budget (ties for first included)
    for b in BUDGETS:
        req = {6030: REQ_MAIN, 30030: REQ_B30, 120030: REQ_B120}[b]
        P(f"NBudFirstPSOVNS{BNAME[b]}", lambda b=b: str(fb["first_count"][str(b)]["PSOBV"]), req)
    P("NBudCloseGapOneTwentyK", lambda: ceil_num(fb["close_gap_pp"]["120030"]["max_gap"], 2), REQ_B120)  # rounded up ("within")
    P("NBudCloseListOneTwentyK", lambda: listing(LAB[a] for a in fb["close_gap_pp"]["120030"]["methods"]), REQ_B120)
    # Horns Rev: validation of the 80-turbine farm against PyWake, methods below PSO-VNS at 120,030
    # the reference is PyWake 2.6.20 NOJ(k=0.04) at bins IDENTICAL to ours (pywake_check.csv, row
    # HR80 / ours_5deg_2.5 / NOJ_k0.04), not the constant 662.5; same value as \NFPyWakeDiffPct
    hv = s.get("hr_validation") or {}
    P("NHRPyWakeDiff", lambda: num(abs(hv["rel_diff_pct"]), 2), ())                 # plain number (text adds \%); our AEP higher
    P("NHRPyWakeOurs", lambda: num(hv["installed80_aep"], 1), ())
    P("NHRPyWakeRef", lambda: num(hv["pywake_aep"], 1), ())
    P("NHRPyWakeLossDiff", lambda: num(abs(hv["loss_diff_pp"]), 2), ())             # pp; our wake loss LOWER than PyWake's
    P("NHRLowerOneTwentyK", lambda: (lambda o: listing(LAB[a] for a in o) if o else "no other method")(
        [a for a in M10 if a != "PSOBV" and (lb.get(a) or {}).get("120030R") and lb[a]["120030R"]["loss"] is not None
         and lb[a]["120030R"]["loss"] < lb["PSOBV"]["120030R"]["loss"]]), hr_req(REQ_B120))
    def hr_best(tag):
        c = [a for a in M10 if (lb.get(a) or {}).get(tag) and lb[a][tag].get("mean_aep") is not None]
        return LAB[max(c, key=lambda a: lb[a][tag]["mean_aep"])]
    P("NHRBestHundredTwentyK", lambda: hr_best("120030R"), hr_req(REQ_B120))      # method with the highest mean AEP
    P("NHRBestThirtyK", lambda: hr_best("30030R"), hr_req(REQ_B30))
    P("NHRBestSixK", lambda: hr_best("6030R"), REQ_HRX)
    P("NHRBestFeasInit", lambda: hr_best("6030F"), hr_req(REQ_FEAS))
    # boundary rule (robustness section): radial projection vs box clipping for SSA
    br = s.get("boundary_rule") or {}
    P("NBoundCases", lambda: WORD.get(br["n_cases"], str(br["n_cases"])), ())
    P("NBoundHigher", lambda: WORD.get(br["n_higher"], str(br["n_higher"])), ())
    P("NBoundMaxP", lambda: pval_up(br["max_p"]), ())
    P("NBoundDiffMin", lambda: num(br["min_diff"], 0), ())
    P("NBoundDiffMax", lambda: num(br["max_diff"], 0), ())

    # ---------------- macros of the additional controls: omega = 0.9 and the disc-sampling control
    # budget split with omega = 0.9 (PSOBV90, 5,430 PSO evaluations): mean loss (%), average rank among the settings of
    # Table split, and the run-level W/T/L of omega = 0.9 against omega = 0.75 FROM THE omega = 0.9 SIDE (W = 0.9
    # significantly better; same Holm family as the table, over the comparisons of each case). SeventyFiveVsNinety:
    # the same tally from the 0.75 side; FiftyVsNinety: 0.5 against 0.9 (0.5 side, the table column).
    P("NSplitLossNinety", lambda: num(sp["mean_loss"]["PSOBV90"], 3), REQ_SPLIT)
    P("NSplitRankNinety", lambda: num(sp["avg_rank"]["PSOBV90"], 2), REQ_SPLIT)
    P("NSplitBestNinety", lambda: WORD.get(sp["best_count"].get("PSOBV90", 0), str(sp["best_count"].get("PSOBV90", 0))), REQ_SPLIT)
    P("NSplitNinetyVsSeventyFive", lambda: "%d/%d/%d" % wtl(sp["wtl_90_vs_75"]), REQ_SPLIT)
    P("NSplitSeventyFiveVsNinety", lambda: "%d/%d/%d" % wtl(sp["wtl_90_vs_75"])[::-1], REQ_SPLIT)
    P("NSplitFiftyVsNinety", lambda: "%d/%d/%d" % wtl(sp["wtl_50_vs_90"]), REQ_SPLIT)
    P("NSplitNinetyLowerThanSeventyFive", lambda: WORD.get(sp["n_lower_loss_90_than_75"], str(sp["n_lower_loss_90_than_75"])), REQ_SPLIT)
    P("NSplitCmNinetyP", lambda: pval(sp["case_mean_omega90"]["PSOBV90-PSOBV75"]["p"]), REQ_SPLIT)     # case-mean 0.9 vs 0.75, unadjusted
    P("NSplitBestSetting", lambda: "\\ensuremath{\\omega=%g}" % {"PSOBV25": 0.25, "PSOBV": 0.5, "PSOBV75": 0.75, "PSOBV90": 0.9, "PSOC": 1}[sp["best_mean_loss"]], REQ_SPLIT)
    P("NSplitNSettings", lambda: WORD[sp["n_settings"]], REQ_SPLIT)                    # "five"
    P("NSplitNComparisons", lambda: WORD[sp["n_comparisons"]], REQ_SPLIT)             # Holm family per case ("six")
    # component analysis with RSD-VNS: family size, average rank in the pool of Table ablation, feasibility, switch
    P("NAblNContrasts", lambda: WORD[ab["n_contrasts"]], REQ_ABL)                     # "fifteen" (Holm family per case)
    P("NAblDf", lambda: str(len(fa["avg_rank"]) - 1), REQ_ABL)                          # Friedman d.f. of Table ablation ("8")
    for a in (ab.get("variants") or []):
        P(f"NAblRank{CODE[a]}", lambda a=a: num(fa["avg_rank"][a], 2), REQ_ABL)
    P("NRankRSDVNS", lambda: num(fa["avg_rank"]["RSDVNS"], 2), REQ_ABL)
    P("NFeasRSDVNS", lambda: num(ab["feasible_pct"]["RSDVNS"], 1), REQ_ABL)
    P("NFeasRSVNS", lambda: num(ab["feasible_pct"]["RSVNS"], 1), REQ_ABL)
    # runs (of 2,040) whose 3,015 Phase-1 samples (30 common initial layouts + 2,985 samples) contain no feasible layout;
    # geometry-only replay of the seeded streams, verified against the stored curves (mpce_results.phase1_replay)
    rp = ab.get("phase1_replay") or {}
    def replay_runs(a):
        if not rp[a]["verification_ok"]:
            raise ValueError(f"Phase-1 replay of {a} does not reproduce the stored curves")
        return num(rp[a]["runs_no_feasible_sample"], 0)
    P("NRSDRunsNoFeasSample", lambda: replay_runs("RSDVNS"), REQ_ABL)
    P("NRSRunsNoFeasSample", lambda: replay_runs("RSVNS"), REQ_ABL)
    P("NRSDPctFeasSamples", lambda: num(rp["RSDVNS"]["pct_feasible_samples"], 1), REQ_ABL)
    P("NRSPctFeasSamples", lambda: num(rp["RSVNS"]["pct_feasible_samples"], 1), REQ_ABL)
    # square vs disc, paired over the same runs. \NRsSquareExplains RS-VNS runs without a feasible
    # sample would have one with disc sampling (the most the square can explain); \NRsSpacingDominates have none either
    # way (spacing / packing density); \NRsNoFeasInside of the RS-VNS runs without a feasible sample do sample layouts
    # with all turbines inside the circle, but none of them is spaced
    sq = rp.get("square_vs_disc") or {}
    P("NRsSquareExplains", lambda: num(sq["rs_none_rsd_some"], 0), REQ_ABL)
    P("NRsSpacingDominates", lambda: num(sq["rs_none_rsd_none"], 0), REQ_ABL)
    P("NRsSquareOnly", lambda: num(sq["rs_some_rsd_none"], 0), REQ_ABL)
    P("NRsNoFeasInside", lambda: num(rp["RSVNS"]["runs_no_feasible_but_inside_sample"], 0), REQ_ABL)
    P("NRsNoInside", lambda: num(rp["RSVNS"]["runs_no_inside_sample"], 0), REQ_ABL)
    return M



def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--summary", default=os.path.join(HERE, "mpce_summary.json"))
    ap.add_argument("--out", default=os.path.join(HERE, "mpce_numbers.tex"))
    ap.add_argument("--allow-partial", action="store_true")
    a = ap.parse_args(argv)
    s = json.load(open(a.summary))
    M = build(s, a.allow_partial)
    fb = [f for f in s.get("fallbacks", []) if "FALLBACK" in f or "fallback" in f or "old-platform" in f]
    hdr = ["% generated by mpce_numbers.py from mpce_summary.json (" + s.get("generated", "?") + ") -- do not edit by hand",
           "% pending macros (expand to \\TBD{pending: ...}): " + (", ".join(M.pending) if M.pending else "none")]
    hdr += ["% data note: " + f for f in s.get("fallbacks", [])]
    open(a.out, "w").write("\n".join(hdr + M.out) + "\n")
    print(f"mpce_numbers.py: wrote {len(M.out)} macros to {os.path.relpath(a.out)} ({len(M.pending)} pending)")
    for e in M.errors:
        print("  value missing:", e)
    return M


if __name__ == "__main__":
    main()

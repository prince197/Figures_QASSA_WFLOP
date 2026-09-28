"""Check the outcome-dependent statements of MPCE_PSO_VNS.tex against mpce_summary.json.

Usage (from analysis/):  python3 mpce_check_final.py [--summary mpce_summary.json] [--tex ../MPCE_PSO_VNS.tex]

Every qualitative statement in the manuscript that was written from the data (e.g. "significantly better",
"never worse", "ranks second under the Gaussian model") carries a comment  % CHECK-FINAL [Cnn]: <condition>.
This script evaluates the same conditions (CONDITIONS below) and prints PASS / FAIL / PENDING (data missing)
for each; it also reports IDs used in the manuscript but not defined here, and vice versa. A FAIL means the
sentence next to the comment must be rewritten. Numbers themselves come from mpce_numbers.tex and need no check.
mpce_results.py runs this script at the end of every run. Exit code 1 if any condition fails.
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
    ("C07", "clearly outperforms SSA, LX-SSA, DE, VNS and MS-SLSQP (each: W >= 40, L <= 2, case-mean p_holm < 0.05)",
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
    ("C15", "SSA-VNS and RS-VNS are statistically tied (T >= 60 of 68, W and L <= 5)",
     "ablation.contrasts.SSABV-RSVNS",
     lambda s: (lambda w: w[1] >= 60 and w[0] <= 5 and w[2] <= 5)(wtl(g(s, "ablation", "contrasts", "SSABV-RSVNS")))),
    ("C16", "PSO is the better first phase than SSA / LX-SSA (never worse, W >= 34) and its best layout at the switch is more often feasible and has lower loss",
     "ablation.contrasts.PSOBV-SSABV/PSOBV-LXBV, ablation.phase2_loss_reduction_pct",
     lambda s: all(wtl(g(s, "ablation", "contrasts", k))[2] == 0 and wtl(g(s, "ablation", "contrasts", k))[0] >= 34 for k in ("PSOBV-SSABV", "PSOBV-LXBV"))
     and (lambda p: all(p["PSOBV"]["feasible_at_switch_pct"] > p[h]["feasible_at_switch_pct"] and p["PSOBV"]["mean_loss_at_switch_pct"] < p[h]["mean_loss_at_switch_pct"]
                        for h in ("SSABV", "LXBV")))(g(s, "ablation", "phase2_loss_reduction_pct"))),
    ("C17", "Laplace step gives no gain: LX-SSA never significantly better than SSA; LX-SSA-VNS never significantly better than SSA-VNS",
     "ablation.contrasts.LXSSA-SSA.W, SSABV-LXBV.L",
     lambda s: wtl(g(s, "ablation", "contrasts", "LXSSA-SSA"))[0] == 0 and wtl(g(s, "ablation", "contrasts", "SSABV-LXBV"))[2] == 0),
    ("C18", "split: 75% has lower loss than 50% in >= 10 of 12 cases, 50% never significantly better than 75%, 25% never significantly better than 50%, ranks 75 < 50 < 25",
     "split.n_lower_loss_than_50, split.wtl_50_vs_75, split.wtl_50_vs_25, split.avg_rank",
     lambda s: (lambda p: p["n_lower_loss_than_50"]["PSOBV75"] >= 10 and p["wtl_50_vs_75"]["W"] == 0 and p["wtl_50_vs_25"]["L"] == 0
                and p["avg_rank"]["PSOBV75"] < p["avg_rank"]["PSOBV"] < p["avg_rank"]["PSOBV25"])(g(s, "split"))),
    ("C19", "Horns Rev 16: PSO-VNS has the highest mean AEP, is significantly better than every other method, and its best run exceeds the installed layout",
     "hr16.methods.*.mean/p_holm/best, hr16.installed_aep",
     lambda s: (lambda h, inst: all(h["PSOBV"]["mean"] > v["mean"] for a, v in h.items() if a != "PSOBV" and v.get("mean") is not None)
                and all(v["p_holm"] < 0.05 for a, v in h.items() if a != "PSOBV") and h["PSOBV"]["best"] > inst)(g(s, "hr16", "methods"), g(s, "hr16", "installed_aep"))),
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
    ("C27", "previous study's method pool (Table baseline): under the old PSO setting at least one salp-swarm hybrid (SSA-VNS or "
            "LX-SSA-VNS) ranks ahead of PSO; under the corrected setting PSO ranks ahead of every salp-swarm method (SSA, LX-SSA, SSA-VNS, LX-SSA-VNS)",
     "baseline.old.avg_rank, baseline.constriction.avg_rank",
     lambda s: (lambda o, c: any(o[h] < o["PSO"] for h in ("SSABV", "LXBV"))
                and all(c["PSO"] < c[a] for a in ("SSA", "LXSSA", "SSABV", "LXBV")))(
         g(s, "baseline", "old", "avg_rank"), g(s, "baseline", "constriction", "avg_rank"))),
]


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--summary", default=os.path.join(HERE, "mpce_summary.json"))
    ap.add_argument("--tex", default=os.path.join(HERE, "..", "MPCE_PSO_VNS.tex"))
    a = ap.parse_args(argv)
    s = json.load(open(a.summary))
    res = {}
    print(f"CHECK-FINAL conditions (summary generated {s.get('generated', '?')})")
    for cid, text, keys, fn in CONDITIONS:
        try:
            r = "PASS" if fn(s) else "FAIL"
        except (KeyError, TypeError, ValueError, IndexError):
            r = "PENDING"
        res[cid] = r
        print(f"  [{cid}] {r:7s} {text}\n            keys: {keys}")
    used = set()
    if os.path.exists(a.tex):
        import glob as _g
        txt = open(a.tex).read() + "".join(open(f).read() for f in sorted(_g.glob(os.path.join(os.path.dirname(a.tex), "optA", "*.tex"))))
        used = set(re.findall(r"CHECK-FINAL \[(C\d+)\]", txt))
        for cid in sorted(used - set(res)):
            print(f"  [{cid}] UNDEFINED: used in the manuscript but not defined in mpce_check_final.py")
        for cid in sorted(set(res) - used):
            print(f"  [{cid}] (not referenced in the manuscript)")
    n = {k: sum(v == k for v in res.values()) for k in ("PASS", "FAIL", "PENDING")}
    print(f"  summary: {n['PASS']} PASS, {n['FAIL']} FAIL, {n['PENDING']} PENDING")
    fb = [f for f in s.get("fallbacks", []) if "->" in f]
    for f in fb:
        print(f"  NOTE (data): {f}")
    return 1 if n["FAIL"] else 0


if __name__ == "__main__":
    sys.exit(main())

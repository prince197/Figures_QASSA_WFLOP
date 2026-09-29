"""Check the statements drawn from the PSO coefficient sweep / bound-handling experiment (mpce_csweep.py) against
mpce_summary_csweep.json.

Usage (from analysis/):  python3 mpce_check_csweep.py [--summary mpce_summary_csweep.json]

Every sentence of the paper or supplement that rests on this experiment should carry a comment
% CHECK-CSWEEP [Snn]: <condition>. This script evaluates the conditions below and prints PASS / FAIL / PENDING
(data missing). A FAIL means the sentence next to the comment must be rewritten. Exit code 1 if any FAIL.
"""
import os, sys, json, argparse

HERE = os.path.dirname(os.path.abspath(__file__))
SWEEP = ["PSOW07C12", "PSOW07C14", "PSOW07C16", "PSOW07C17", "PSOW07C18", "PSOW07C19", "PSOW07C20"]


def g(d, *path):
    for p in path:
        if d is None:
            raise KeyError("/".join(map(str, path)))
        d = d[str(p)] if isinstance(d, dict) and str(p) in d else d[p]
    if d is None:
        raise KeyError("/".join(map(str, path)))
    return d


def st(s, a, k):
    return g(s, "settings", a, k)


BELOW = ["PSOW07C12", "PSOW07C14", "PSOW07C16", "PSOW07C17"]         # order-2 stable at w = 0.7 (c < 1.7486)
ABOVE = ["PSOW07C18", "PSOW07C19", "PSOW07C20"]                      # not order-2 stable


def pr(s, key, what, k):
    return g(s, "paired", key, what, k)


def step(s, key, c_from):
    return next(x for x in g(s, "response_steps", key) if abs(x["c_from"] - c_from) < 1e-9)["change_per_0p1"]


def fdrop(s, c_from):
    return next(x for x in g(s, "steps") if abs(x["c_from"] - c_from) < 1e-9)["drop_pp_per_0p1"]


CONDITIONS = [
    ("S01", "the c = 2.0 sweep point reproduces all 360 stored runs of the old setting (fresh_grid PSO) bit for bit",
     lambda s: g(s, "verification", "old_runs_identical") == g(s, "verification", "old_runs_compared") == 360),
    ("S02", "order-2 stable side: for every c <= 1.7 at least 98 % of the final layouts are feasible and all 12 cases "
            "are qualified (>= 15 of 30 feasible)",
     lambda s: all(st(s, a, "feas_pct") >= 98 and st(s, a, "cases_qualified") == 12 for a in BELOW)),
    ("S03", "beyond the boundary feasibility falls monotonically with c (c = 1.8 > 1.9 > 2.0) and is below 95 % for "
            "every c > c*; below 60 % at c = 2.0",
     lambda s: st(s, "PSOW07C18", "feas_pct") > st(s, "PSOW07C19", "feas_pct") > st(s, "PSOW07C20", "feas_pct")
     and all(st(s, a, "feas_pct") < 95 for a in ABOVE) and st(s, "PSOW07C20", "feas_pct") < 60),
    ("S04", "the loss of feasibility starts in the interval (1.7, 1.8) that contains c* = 1.7486: every step below it "
            "loses < 1 pp per 0.1 in c, the step 1.7 -> 1.8 loses >= 5 pp, and c = 1.7 vs 1.8 differs in feasibility "
            "(McNemar p < 0.001)",
     lambda s: all(fdrop(s, c) < 1 for c in (1.2, 1.4, 1.6)) and fdrop(s, 1.7) >= 5
     and next(x for x in g(s, "steps") if x["c_from"] == 1.7)["contains_c_star"]
     and pr(s, "PSOW07C17_vs_PSOW07C18", "feas", "p") < 0.001),
    ("S05", "no step at c*: the decline beyond c* is progressive and accelerating (drop per 0.1 in c: 1.9 -> 2.0 > "
            "1.8 -> 1.9 > 1.7 -> 1.8), i.e. the largest drop is not the step that contains c*",
     lambda s: fdrop(s, 1.9) > fdrop(s, 1.8) > fdrop(s, 1.7) and not g(s, "largest_drop", "contains_c_star")),
    ("S06", "the swarm dynamics change smoothly through c*: the final spread grows monotonically with c, already by a "
            "factor >= 5 from c = 1.4 to 1.7 (stable side), and the largest decrease of the feasible-evaluation share "
            "(iterations 101-200) per 0.1 in c lies below c* (1.6 -> 1.7)",
     lambda s: all(st(s, SWEEP[i], "spread_rel_median") < st(s, SWEEP[i + 1], "spread_rel_median") for i in range(6))
     and g(s, "spread_ratio_c17_over_c14") >= 5
     and min(g(s, "response_steps", "feas_eval_late_pct_mean"), key=lambda x: x["change_per_0p1"])["c_to"] <= 1.7),
    ("S07", "final quality degrades already on the stable side: mean loss (common cases) is non-decreasing from c = 1.4 "
            "to 2.0, and c = 1.6 beats c = 1.7 in >= 10 of 12 cases (Wilcoxon p < 0.05)",
     lambda s: all(st(s, SWEEP[i], "loss_common") <= st(s, SWEEP[i + 1], "loss_common") + 1e-12 for i in range(1, 6))
     and pr(s, "PSOW07C16_vs_PSOW07C17", "loss", "a_better") >= 10 and pr(s, "PSOW07C16_vs_PSOW07C17", "loss", "p") < 0.05),
    ("S08", "crossing from c = 1.7 to 1.8 worsens the loss in all 12 cases (Wilcoxon p < 0.001) and the average rank "
            "change 1.7 -> 1.8 is the largest along the sweep",
     lambda s: pr(s, "PSOW07C17_vs_PSOW07C18", "loss", "a_better") == 12 and pr(s, "PSOW07C17_vs_PSOW07C18", "loss", "p") < 0.001
     and max(g(s, "response_steps", "avg_rank"), key=lambda x: x["change_per_0p1"])["c_from"] == 1.7),
    ("S09", "zeroing the clipped velocity components does NOT rescue the old setting: feasibility within 3 pp of the "
            "implemented rule (McNemar p > 0.05), no more qualified cases, loss not better",
     lambda s: abs(st(s, "PSOOLD_VZERO", "feas_pct") - st(s, "PSOW07C20", "feas_pct")) <= 3
     and pr(s, "PSOOLD_VZERO_vs_PSOW07C20", "feas", "p") > 0.05
     and st(s, "PSOOLD_VZERO", "cases_qualified") <= st(s, "PSOW07C20", "cases_qualified")
     and pr(s, "PSOOLD_VZERO_vs_PSOW07C20", "loss", "mean_diff_pp") >= 0),
    ("S10", "velocity clamping (vmax = 0.2 (ub - lb)) largely rescues feasibility of the old setting: >= 90 % feasible "
            "(+ >= 30 pp, McNemar p < 0.001, no run lost), all 12 cases qualified, lower loss in every common case",
     lambda s: st(s, "PSOOLD_VMAX", "feas_pct") >= 90
     and st(s, "PSOOLD_VMAX", "feas_pct") - st(s, "PSOW07C20", "feas_pct") >= 30
     and pr(s, "PSOOLD_VMAX_vs_PSOW07C20", "feas", "p") < 0.001 and pr(s, "PSOOLD_VMAX_vs_PSOW07C20", "feas", "b_only") == 0
     and st(s, "PSOOLD_VMAX", "cases_qualified") == 12
     and pr(s, "PSOOLD_VMAX_vs_PSOW07C20", "loss", "a_better") == pr(s, "PSOOLD_VMAX_vs_PSOW07C20", "loss", "cases") > 0),
    ("S11", "... but not fully: with vmax the old setting stays worse than the constriction setting and than c = 1.7 in "
            "feasibility (McNemar p < 0.01) and in loss in all 12 cases",
     lambda s: all(pr(s, k, "feas", "p") < 0.01 and pr(s, k, "loss", "b_better") == 12
                   for k in ("PSOOLD_VMAX_vs_PSOC", "PSOOLD_VMAX_vs_PSOW07C17"))),
    ("S12", "the old setting with vmax is statistically indistinguishable from c = 1.8 without vmax (feasibility McNemar "
            "p > 0.05, loss Wilcoxon p > 0.05)",
     lambda s: pr(s, "PSOOLD_VMAX_vs_PSOW07C18", "feas", "p") > 0.05 and pr(s, "PSOOLD_VMAX_vs_PSOW07C18", "loss", "p") > 0.05),
    ("S13", "the constriction setting has the best average rank of the ten settings; c = 1.4 (w = 0.7) is 100 % feasible "
            "and not significantly different from it in loss (Wilcoxon p > 0.05)",
     lambda s: g(s, "best_avg_rank") == "PSOC" and st(s, "PSOW07C14", "feas_pct") == 100
     and pr(s, "PSOW07C14_vs_PSOC", "loss", "p") > 0.05),
    ("S14", "too small coefficients also hurt: c = 1.2 has a higher mean loss (common cases) and a worse average rank than "
            "c = 1.4 and c = 1.6",
     lambda s: st(s, "PSOW07C12", "loss_common") > max(st(s, "PSOW07C14", "loss_common"), st(s, "PSOW07C16", "loss_common"))
     and st(s, "PSOW07C12", "avg_rank") > max(st(s, "PSOW07C14", "avg_rank"), st(s, "PSOW07C16", "avg_rank"))),
    ("S15", "feasibility of every run is identical in Data Sets I and II (same geometry), so feasibility statistics rest "
            "on 6 geometries x 30 seeds (McNemar on 180 pairs)",
     lambda s: g(s, "feas_identical_across_datasets") is True),
]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--summary", default=os.path.join(HERE, "mpce_summary_csweep.json"))
    a = ap.parse_args()
    if not os.path.exists(a.summary):
        print(f"PENDING: {a.summary} missing (run mpce_csweep.py)"); sys.exit(0)
    s = json.load(open(a.summary))
    n = dict(PASS=0, FAIL=0, PENDING=0)
    for cid, text, cond in CONDITIONS:
        try:
            r = "PASS" if cond(s) else "FAIL"
        except (KeyError, TypeError, IndexError, ValueError):
            r = "PENDING"
        n[r] += 1
        print(f"{cid} {r:7s} {text}")
    print(f"\n{n['PASS']} PASS / {n['FAIL']} FAIL / {n['PENDING']} PENDING")
    sys.exit(1 if n["FAIL"] else 0)


if __name__ == "__main__":
    main()

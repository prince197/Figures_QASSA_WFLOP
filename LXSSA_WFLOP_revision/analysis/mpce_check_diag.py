"""Check the statements drawn from the Phase-6 diagnostics (mpce_diagnostics.py) against mpce_summary_diag.json.

Usage (from analysis/):  python3 mpce_check_diag.py [--summary mpce_summary_diag.json]

Every sentence of the paper or supplement that rests on a diagnostic should carry a comment
% CHECK-DIAG [Dnn]: <condition>. This script evaluates the conditions below and prints PASS / FAIL / PENDING
(data missing). A FAIL means the sentence next to the comment must be rewritten. Exit code 1 if any FAIL.
"""
import os, sys, json, argparse

HERE = os.path.dirname(os.path.abspath(__file__))


def g(d, *path):
    for p in path:
        if d is None:
            raise KeyError("/".join(map(str, path)))
        d = d[str(p)] if isinstance(d, dict) and str(p) in d else d[p]
    if d is None:
        raise KeyError("/".join(map(str, path)))
    return d


V3 = ("PSOBV", "SSABV", "RSVNS")

CONDITIONS = [
    ("D01", "instrumented runs reproduce the stored runs of the paper (WakeLoss, coordinates, curve, calls)",
     lambda s: g(s, "verification", "stored_match") == g(s, "verification", "instrumented_runs")
     and not g(s, "verification", "stored_mismatch") and not g(s, "verification", "stored_missing")),
    ("D02", "instrumentation does not change behaviour: seed-1 reruns with the original classes are bit-identical",
     lambda s: g(s, "verification", "bitident_runs") > 0
     and g(s, "verification", "bitident_match") == g(s, "verification", "bitident_runs")),
    ("D03", "old PSO setting keeps a much wider swarm: final spread and |V| at least 3x the constriction values "
            "(median paired ratio) and larger in every paired run",
     lambda s: g(s, "T1", "paired", "spread_ratio_median") >= 3 and g(s, "T1", "paired", "vel_ratio_median") >= 3
     and g(s, "T1", "paired", "spread_ratio_min") > 1 and g(s, "T1", "paired", "vel_ratio_min") > 1),
    ("D04", "old setting clips more coordinates: mean clipped share over iterations 1-200 at least 2x constriction, and "
            "higher in every paired run over iterations 101-200",
     lambda s: g(s, "T1", "PSO", "clip_pct_all") >= 2 * g(s, "T1", "PSOC", "clip_pct_all")
     and g(s, "T1", "paired", "n_old_clip_higher_late") == g(s, "T1", "runs_per_setting")),
    ("D05", "old setting rarely samples feasible layouts (< 10 % of particle evaluations in iterations 101-200 vs "
            "> 25 % for constriction) and returns fewer feasible global bests (constriction: all)",
     lambda s: g(s, "T1", "PSO", "feas_particles_pct_late") < 10 and g(s, "T1", "PSOC", "feas_particles_pct_late") > 25
     and g(s, "T1", "PSO", "gbest_feasible_final") < g(s, "T1", "PSOC", "gbest_feasible_final") == g(s, "T1", "runs_per_setting")),
    ("D06", "stored runs, classified with the paper's 1e-6 m rule (R3-7, R5-F2): the infeasible finals of the old setting "
            "are exactly the runs without the Feasible flag (runs - feasible, 2,040 - 1,754 = 286, i.e. 86.0 % feasible), "
            "boundary-only + spacing-only + both add up to them, none is undecidable at the 1-mm storage precision after the "
            "reproduced runs (all reproduced runs match the stored ones) and none contradicts the flag; most violate ONLY "
            "the spacing constraint (> 50 %; boundary-only < spacing-only); the constriction setting has none",
     lambda s: (lambda v, w, ver: v["infeasible"] == v["runs"] - v["feasible_flag"]
                and v["boundary_only"] + v["spacing_only"] + v["both"] == v["infeasible"]
                and v["undecided"] == 0 and v["flag_contradictions"] == 0
                and ver["t1b_match"] == ver["t1b_runs"] == ver["t1b_expected"]
                and v["spacing_only"] > 0.5 * v["infeasible"] and v["boundary_only"] < v["spacing_only"]
                and w["infeasible"] == 0 and w["runs"] == w["feasible_flag"])(
         g(s, "T1_stored", "PSO"), g(s, "T1_stored", "PSOC"), g(s, "verification"))),
    ("D07", "old setting: most infeasible particle evaluations violate both constraints, i.e. the swarm samples mostly "
            "outside the circle (turbine outside the circle in > 50 % of evaluations, iterations 101-200)",
     lambda s: g(s, "T1", "PSO", "infeasible_both_pct") > 50 and g(s, "T1", "PSO", "bviol_particles_pct_late") > 50),
    ("D08", "the VNS phase is essentially one compass-search descent: for PSO-VNS, SSA-VNS and RS-VNS the median number "
            "of shaking steps is <= 2, shaking costs < 1 % of the Phase-2 evaluations, >= 85 % of the Phase-2 improvement "
            "comes from the first descent, and >= 80 % of the runs have no accepted shaking cycle",
     lambda s: all(g(s, "T2", a, "shakes_median") <= 2 and g(s, "T2", a, "shake_evals_pct_of_phase2") < 1
                   and g(s, "T2", a, "imp_descent_share_pct") >= 85
                   and g(s, "T2", a, "runs_no_accepted_cycle") >= 0.8 * g(s, "T2", a, "runs") for a in V3)),
    ("D09", "for N >= 12 the median run has no shaking step: the first descent does not finish within the Phase-2 budget "
            "in at least 90 % of the runs with N = 12 or 15 (and in all runs with N = 15)",
     lambda s: all(g(s, "T2", a, "shakes_median_by_n", n) == 0 for a in V3 for n in (12, 15))
     and all(g(s, "T2", a, "descent_not_finished_by_n", 12) + g(s, "T2", a, "descent_not_finished_by_n", 15)
             >= 0.9 * 4 * g(s, "T2", "seeds")                   # N = 12 and 15, both data sets
             and g(s, "T2", a, "descent_not_finished_by_n", 15) == 2 * g(s, "T2", "seeds") for a in V3)),
    ("D10", "accepted shaking cycles occur only in the smallest neighbourhood k = 1; k_max = 5 is never reached",
     lambda s: all(sum(g(s, "T2", a, "accepted_by_k")[1:]) == 0 and g(s, "T2", a, "runs_reaching_kmax") == 0 for a in V3)),
    ("D11", "VNS never turns a feasible incumbent infeasible; infeasible switch points are repaired only by the first "
            "descent (never by a later cycle)",
     lambda s: all(g(s, "T2", a, "lost_feasibility") == 0 and g(s, "T2", a, "made_feasible_cycle") == 0 for a in V3)),
    ("D12", "the start handed to VNS is the best feasible Phase-1 layout whenever Phase 1 evaluated one "
            "(so an 'rsbest' control would equal RS-VNS)",
     lambda s: all(g(s, "T2", a, "start_not_best_feasible") == 0 for a in V3)),
    ("D13", "feasible starts: PSO and DE never improve on the best initial layout (instrumented runs and all stored runs)",
     lambda s: all(g(s, "T3", a, "improved_runs") == 0 and g(s, "T3", "stored", a, "runs_final_differs_from_initial_best") == 0
                   for a in ("PSOC", "DE"))),
    ("D14", "feasible starts: every feasible PSO candidate is the G-particle (zero velocity) re-evaluating G; no personal "
            "best is ever updated; DE produces no feasible trial and accepts none",
     lambda s: g(s, "T3", "PSOC", "feasible_candidates_not_G_reevaluation") == 0
     and g(s, "T3", "PSOC", "iterations_exactly_one_feasible") == g(s, "T3", "PSOC", "iterations")
     and g(s, "T3", "PSOC", "runs_gbest_particle_zero_velocity_throughout") == g(s, "T3", "PSOC", "runs")
     and g(s, "T3", "PSOC", "pbest_updates") == 0
     and g(s, "T3", "DE", "feasible_candidates") == 0 and g(s, "T3", "DE", "accepted_trials") == 0),
    ("D15", "label symmetry alone does not explain the stall: after optimal relabelling the first PSO move is still "
            "feasible in < 1 % of draws; only a common weight (convex combination) is feasible more often",
     lambda s: g(s, "T3", "PSOC", "mixture", "raw", "feasible_pct") < 1
     and g(s, "T3", "PSOC", "mixture", "matched", "feasible_pct") < 1
     and g(s, "T3", "PSOC", "mixture", "matched_scalar", "feasible_pct") > 5 * max(0.01, g(s, "T3", "PSOC", "mixture", "matched", "feasible_pct"))),
    ("D16", "the feasible initial layouts are distinct, not permutations of each other (matched RMS distance to G > 50 m)",
     lambda s: g(s, "T3", "PSOC", "init_matched_rms_min_m") > 50),
    ("D17", "RS-VNS Phase-1 replay agrees with the runs, the inside-circle rate matches (pi/4)^N (within 5 % relative "
            "for N <= 10), and for N = 15 in the 1000-m farm none of the 90,450 box samples is feasible",
     lambda s: g(s, "verification", "rs_replay_agree") == g(s, "verification", "rs_replay_checked") > 0
     and all(abs(v["p_inside"] / v["p_inside_theory"] - 1) < 0.05 for v in g(s, "rs_replay").values() if v["n"] <= 10)
     and g(s, "rs_replay", "1000-15", "feasible") == 0),
    ("D18", "share of runs with a feasible layout by the switch (3,030 calls; RS-VNS 3,015, read at the checkpoint at "
            "3,000): PSO > SSA > RS (68 cases)",
     lambda s: g(s, "first_feasible", "PSOBV", "pct_by_switch") > g(s, "first_feasible", "SSABV", "pct_by_switch")
     > g(s, "first_feasible", "RSVNS", "pct_by_switch")),
    ("D19", "the planned cloud experiments are small: omega90 and rsdisc each < 2 CPU-hours (stored run times)",
     lambda s: g(s, "T4", "omega90", "cpu_hours") < 2 and g(s, "T4", "rsdisc", "cpu_hours") < 2),
    ("D20", "VNS cycles (R3-2): the shake-plus-search cycles take a sizeable share of the Phase-2 evaluations but give "
            "little of the Phase-2 gain -- for every method their evaluation share (runs with a feasible switch point, "
            "pooled) is more than twice their gain share (accepted cycles plus the cycle cut off by the budget); for "
            "PSO-VNS the cycles take >= 15 % of the evaluations and give < 5 % of the gain (accepted cycles < 1 %)",
     lambda s: all(g(s, "T2", a, "cycle_evals_pct_of_phase2_imp_runs") > 2 * g(s, "T2", a, "imp_cycles_all_share_pct") for a in V3)
     and g(s, "T2", "PSOBV", "cycle_evals_pct_of_phase2") >= 15 and g(s, "T2", "PSOBV", "imp_cycles_all_share_pct") < 5
     and g(s, "T2", "PSOBV", "imp_cycles_share_pct") < 1),
    ("D21", "old-setting boundary violations are box-tangent pinning (R3-7): at least 90 % of the stored old-setting "
            "finals that violate the boundary (1e-6 m rule) have an outside turbine with a coordinate on the box bound",
     lambda s: g(s, "T1_stored", "PSO", "tangent") >= 0.9 * g(s, "T1_stored", "PSO", "any_boundary") > 0),
    ("D22", "mpce_numbers_diag.tex is in sync with mpce_summary_diag.json (R5 reproducibility 2): the macros regenerated "
            "from the stored summary with mpce_diagnostics.macros equal the file",
     lambda s: _macros_in_sync(s)),
]


def _macros_in_sync(s):
    import re
    sys.path.insert(0, HERE)
    import mpce_diagnostics as MD            # post-processing only (no run is started on import)
    want = {k: str(v) for k, v in MD.macros(s).items()}
    have = dict(re.findall(r"\\newcommand\{\\(\w+)\}\{(.*)\}\s*$", open(os.path.join(HERE, "mpce_numbers_diag.tex")).read(), re.M))
    bad = sorted(k for k in set(want) | set(have) if want.get(k) != have.get(k))
    if bad:
        print("  D22 differences:", ", ".join(f"{k}: file {have.get(k)} / summary {want.get(k)}" for k in bad[:10]))
    return not bad


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--summary", default=os.path.join(HERE, "mpce_summary_diag.json"))
    a = ap.parse_args()
    s = json.load(open(a.summary))
    n = dict(PASS=0, FAIL=0, PENDING=0)
    for cid, text, cond in CONDITIONS:
        try:
            r = "PASS" if cond(s) else "FAIL"
        except (KeyError, TypeError, IndexError, ImportError, FileNotFoundError) as e:
            r = f"PENDING"
        n[r] += 1
        print(f"{cid} {r:7s} {text}")
    print(f"\n{n['PASS']} PASS / {n['FAIL']} FAIL / {n['PENDING']} PENDING")
    sys.exit(1 if n["FAIL"] else 0)


if __name__ == "__main__":
    main()

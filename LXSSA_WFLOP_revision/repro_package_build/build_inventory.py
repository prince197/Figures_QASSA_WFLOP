#!/usr/bin/env python3
"""Inventory of the SWEVO reproducibility package: classify every candidate file of the repository.

Usage (from anywhere):
    python3 repro_package_build/build_inventory.py [--root REPO] [--out repro_package_build]

Writes, in --out:
  inventory.csv       one row per candidate file: path, category, decision (include / exclude / pending), bytes,
                      note. Every file below analysis/ (recursively, caches excluded), figures_mpce/,
                      figures_final/, SWEVO_rev2/experiments_docs/ and SWEVO_rev2/validation/ is classified;
                      a file that no rule matches is reported as "unclassified" (excluded) so that nothing is
                      silently added or lost.
  package_files.txt   the included paths (relative to the repository root), one per line, in the order used by
                      make_package.sh.
Rules are first-match; analysis/rev3_* files (revision 3, produced by other work packages) are included
automatically when present, and the expected-but-absent ones are listed with decision "pending".
The script only reads the repository; it never modifies a file outside --out.
"""
import argparse, csv, fnmatch, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))

# (pattern relative to the repository root, category, decision, note). First match wins.
RULES = [
    # ---- caches and editor files: never packaged -------------------------------------------------------------
    ("*/__pycache__/*", "cache", "exclude", "Python byte-code cache"),
    ("*.pyc", "cache", "exclude", "Python byte-code cache"),
    # ---- legacy / superseded material (other papers or superseded versions) ----------------------------------
    ("analysis/superseded_linear_penalty/*", "legacy", "exclude", "earlier runs with a reconstructed linear penalty (superseded)"),
    ("analysis/build_*_manuscript.py", "legacy", "exclude", "manuscript builders of earlier papers (LX-SSA versions)"),
    ("analysis/build_section7.py", "legacy", "exclude", "manuscript builder of an earlier paper"),
    ("analysis/final_frontmatter.py", "legacy", "exclude", "manuscript builder of an earlier paper"),
    ("analysis/hybrid_frontmatter.py", "legacy", "exclude", "manuscript builder of an earlier paper"),
    ("analysis/final_prose.py", "legacy", "exclude", "manuscript builder of an earlier paper"),
    ("analysis/flatten_swevo.py", "manuscript-tool", "exclude", "LaTeX flattening of the SWEVO manuscript (manuscript sources are not part of the data package)"),
    ("analysis/swevo_bib.py", "manuscript-tool", "exclude", "bibliography tool of the manuscript"),
    ("analysis/build_swevo.sh", "manuscript-tool", "exclude", "pdflatex of the manuscript; needs the manuscript sources (not in the data package)"),
    ("analysis/build_mpce_paper.sh", "manuscript-tool", "exclude", "MPCE-version build (pdflatex of MPCE_PSO_VNS.tex); the Python stages are in regenerate_all.sh"),
    ("analysis/final_results.py", "legacy", "exclude", "earlier LX-SSA-VNS paper (final_* outputs)"),
    ("analysis/final_tables.tex", "legacy", "exclude", "earlier paper"),
    ("analysis/final_summary.json", "legacy", "exclude", "earlier paper"),
    ("analysis/final_case_tests.csv", "legacy", "exclude", "earlier paper"),
    ("analysis/final_ablation_tests.csv", "legacy", "exclude", "earlier paper"),
    ("analysis/final_best_layouts_maxN.csv", "input", "include", "best layouts of an earlier analysis; read by make_model_figures.py (example layout of Fig. 1a)"),
    ("analysis/final_robust_*", "legacy", "exclude", "earlier paper"),
    ("analysis/robustness_final.py", "legacy", "exclude", "earlier paper; its output final_reevaluation.csv is kept as the cache seed of mpce_robustness.py"),
    ("analysis/final_reevaluation.csv", "input", "include", "cache seed read by mpce_robustness.py (re-evaluations of the fresh_* layouts)"),
    ("analysis/fresh_results.py", "legacy", "exclude", "six-method study of an earlier paper"),
    ("analysis/fresh_tables.tex", "legacy", "exclude", "earlier paper"),
    ("analysis/fresh_summary.json", "legacy", "exclude", "earlier paper"),
    ("analysis/fresh_case_tests.csv", "legacy", "exclude", "earlier paper"),
    ("analysis/fresh_best_layouts_maxN.csv", "legacy", "exclude", "earlier paper"),
    ("analysis/fresh_robust_*", "legacy", "exclude", "earlier paper"),
    ("analysis/fresh_reevaluation.csv", "legacy", "exclude", "earlier paper"),
    ("analysis/robustness_fresh.py", "legacy", "exclude", "earlier paper"),
    ("analysis/hybrid_lxssa_bvns.py", "optimizer", "include", "two-phase hybrids PSO-VNS / SSA-VNS / LX-SSA-VNS"),
    ("analysis/hybrid_lxssa_vns.py", "model", "include", "superseded hybrid with the modified VNS; kept because full_grid_experiments.py imports it (rows VNS/LXVNS are not used)"),
    ("analysis/hybrid_*", "legacy", "exclude", "superseded intermediate hybrid version"),
    ("analysis/robustness_hybrid.py", "legacy", "exclude", "superseded intermediate hybrid version"),
    ("analysis/bvns_*", "legacy", "exclude", "superseded intermediate version"),
    ("analysis/fresh_hgrid.csv", "runs-audited", "include", "superseded / earlier-paper runs; not used by the SWEVO results but part of the 56,540 records of the precision audit (rev3_precision_audit.py, C4)"),
    ("analysis/fresh_hsplit.csv", "runs-audited", "include", "superseded / earlier-paper runs; not used by the SWEVO results but part of the 56,540 records of the precision audit (rev3_precision_audit.py, C4)"),
    ("analysis/fresh_hhr*.csv", "runs-audited", "include", "superseded / earlier-paper runs; not used by the SWEVO results but part of the 56,540 records of the precision audit (rev3_precision_audit.py, C4)"),
    ("analysis/fresh_bsplit.csv", "runs-audited", "include", "superseded / earlier-paper runs; not used by the SWEVO results but part of the 56,540 records of the precision audit (rev3_precision_audit.py, C4)"),
    ("analysis/fresh_hr80.csv", "runs-audited", "include", "superseded / earlier-paper runs; not used by the SWEVO results but part of the 56,540 records of the precision audit (rev3_precision_audit.py, C4)"),
    ("analysis/fresh_vhr80.csv", "runs-audited", "include", "superseded / earlier-paper runs; not used by the SWEVO results but part of the 56,540 records of the precision audit (rev3_precision_audit.py, C4)"),
    ("analysis/fresh_bhr80.csv", "runs-audited", "include", "superseded / earlier-paper runs; not used by the SWEVO results but part of the 56,540 records of the precision audit (rev3_precision_audit.py, C4)"),
    ("analysis/followup_statistics.py", "legacy", "exclude", "earlier analyses (superseded)"),
    ("analysis/reevaluate_layouts.py", "legacy", "exclude", "earlier analyses (superseded)"),
    ("analysis/make_robustness_tables.py", "legacy", "exclude", "earlier analyses (superseded)"),
    ("analysis/robustness_tables.tex", "legacy", "exclude", "earlier analyses (superseded)"),
    ("analysis/robustness_text.tex", "legacy", "exclude", "earlier draft text"),
    ("analysis/robustness_summary.csv", "legacy", "exclude", "earlier analyses (superseded)"),
    ("analysis/layout_reevaluation.csv", "legacy", "exclude", "earlier analyses (superseded)"),
    ("analysis/bootstrap_ci.csv", "legacy", "exclude", "earlier analyses (superseded)"),
    ("analysis/friedman_case_level.csv", "legacy", "exclude", "earlier analyses (superseded)"),
    ("analysis/batch_runtime_per_run.csv", "legacy", "exclude", "earlier analyses (superseded)"),
    ("analysis/best_layouts_selected_cases.csv", "legacy", "exclude", "earlier analyses (superseded)"),
    ("analysis/extra_experiments.py", "legacy", "exclude", "earlier VNS/SLSQP/Horns Rev runs (superseded by mpce_* and hrfix)"),
    ("analysis/extra_experiments.log", "legacy", "exclude", "earlier paper"),
    ("analysis/extra_runs_*.csv", "legacy", "exclude", "earlier paper (superseded)"),
    ("analysis/extra_tables.tex", "legacy", "exclude", "earlier paper"),
    ("analysis/analyze_extra_runs.py", "legacy", "exclude", "earlier paper"),
    ("analysis/six_method_average_ranks.csv", "legacy", "exclude", "earlier paper"),
    ("analysis/sec_*.tex", "legacy", "exclude", "manuscript text of an earlier paper"),
    ("analysis/hornsrev_site_text.tex", "legacy", "exclude", "manuscript text of an earlier paper"),
    ("analysis/README.md", "legacy-doc", "exclude", "history of all analysis versions; superseded for SWEVO by README_reproduce.md and the package README"),
    # ---- revision 3 (other work packages; included automatically when present) -------------------------------
    ("analysis/rev3_precision_results/*.jsonl", "rev3-data", "include", "raw full-precision rerun results (C4), one JSON line per record; merged by rev3_precision_rerun.py collect"),
    ("analysis/rev3p_*.log", "rev3-log", "include", "worker logs of the 35 cloud shards of the full-precision reruns (C4)"),
    ("analysis/rev3_review_report.md", "rev3-editorial", "exclude", "internal editorial review notes on the manuscript text"),
    ("analysis/rev3_suggested_text.md", "rev3-editorial", "exclude", "internal drafting notes for the manuscript text"),
    ("analysis/rev3_refs_*", "rev3-editorial", "exclude", "internal reference-list check of the manuscript"),
    ("analysis/rev3_*.py", "rev3-code", "include", "revision-3 script"),
    ("analysis/rev3_*.csv", "rev3-data", "include", "revision-3 run records / outputs"),
    ("analysis/rev3_*.csv.gz", "rev3-data", "include", "revision-3 run records / outputs"),
    ("analysis/rev3_*.json", "rev3-output", "include", "revision-3 output"),
    ("analysis/rev3_*.tex", "rev3-output", "include", "revision-3 table"),
    ("analysis/rev3_*.log", "rev3-log", "include", "revision-3 log"),
    ("analysis/rev3_*.md", "rev3-doc", "include", "revision-3 manifest / prespecification"),
    ("analysis/rev3_*", "unclassified", "exclude", "revision-3 file of an unexpected type (check by hand)"),
    # ---- models, objectives and site data ---------------------------------------------------------------------
    ("analysis/wflop_model.py", "model", "include", "Kusiak-Song benchmark evaluator (Jensen wake, Weibull rose), cubic curve, Gaussian wake"),
    ("analysis/authors_objective.py", "model", "include", "penalized objective (exact quadratic penalty of the original objective.py)"),
    ("analysis/objective_original.py", "model", "include", "the original authors' objective.py, verbatim (reference for validate_evaluator.py)"),
    ("analysis/hornsrev_model.py", "model", "include", "Horns Rev 1 site model (V80, measured 12-sector climate from PyWake 2.6.20)"),
    ("analysis/iea37_model.py", "model", "include", "IEA37 Case Study 1 evaluator"),
    ("analysis/rev2_site_model.py", "model", "include", "Lillgrund block model (PyWake inputs; wind climate after Gocmen and Giebel 2016)"),
    ("analysis/iea37_data/*", "input", "include", "IEA Wind Task 37 case-study files (repo iea37-wflo-casestudies, cs1-2, commit 267f6e5), unmodified"),
    ("analysis/record_io.py", "model", "include", "17-significant-digit coordinate writer / reader"),
    # ---- optimizers ---------------------------------------------------------------------------------------------
    ("analysis/authors_optimizers.py", "optimizer", "include", "GA/PSO/DE/SSA/LX-SSA of the original code (logic and random-call order unchanged)"),
    ("analysis/extra_baselines.py", "optimizer", "include", "MS-SLSQP (and the modified VNS of an earlier version)"),
    ("analysis/original_vns.py", "optimizer", "include", "basic VNS"),
    ("analysis/hybrid_lxssa_bvns.py", "optimizer", "include", "two-phase hybrids PSO-VNS / SSA-VNS / LX-SSA-VNS"),
    ("analysis/rs_vns.py", "optimizer", "include", "RS-VNS (random-sampling phase 1 + VNS)"),
    ("analysis/feasible_init.py", "optimizer", "include", "feasibility-preserving initialization"),
    ("analysis/init_hook.py", "optimizer", "include", "initial-population hook (random / feasible)"),
    # ---- experiment drivers -------------------------------------------------------------------------------------
    ("analysis/full_grid_experiments.py", "driver", "include", "runs fresh_grid / fresh_vgrid / fresh_bgrid (+ old Horns Rev files)"),
    ("analysis/mpce_experiments.py", "driver", "include", "runs mpce_<exp>_s<i>of<k>.csv (PSO, PSO-VNS, RS-VNS, RSD-VNS, MS-SLSQP, split, init., budget, hrfix, csweep)"),
    ("analysis/iea37_experiments.py", "driver", "include", "runs mpce_iea16/iea36(+p)"),
    ("analysis/rev2_ga.py", "driver", "include", "revision 2: real-coded GA (ga, gahr, gaiea)"),
    ("analysis/rev2_gradient.py", "driver", "include", "revision 2: exact-gradient MS-SLSQP / PSO-SLSQP on IEA37"),
    ("analysis/rev2_laplace.py", "driver", "include", "revision 2: Laplace ablation"),
    ("analysis/rev2_site.py", "driver", "include", "revision 2: Lillgrund block runs"),
    ("analysis/rev2_spacing.py", "driver", "include", "revision 2: 5D / 6D spacing re-optimization"),
    ("analysis/authors_experiments.py", "driver", "include", "runs authors_runs_*.csv (spacing table tab:spacing-authors)"),
    ("analysis/ssa_reference.py", "driver", "include", "SSA with radial projection (boundary-rule check)"),
    ("analysis/mpce_diagnostics.py", "driver", "include", "instrumented diagnostic reruns (Section S-diag)"),
    ("analysis/packing_capacity.py", "driver", "include", "packing bounds (tab:capacity)"),
    ("analysis/pywake_check.py", "driver", "include", "Horns Rev 1 PyWake NOJ cross-check (needs py_wake 2.6.20)"),
    ("analysis/iea37_projected.py", "driver", "include", "projected IEA37 layouts"),
    ("analysis/calibrate_authors_code.py", "check", "include", "re-run vs recorded distributions of the original code"),
    ("analysis/validate_evaluator.py", "check", "include", "evaluator vs 720 archived objectives"),
    # ---- analysis and checks ------------------------------------------------------------------------------------
    ("analysis/mpce_results.py", "analysis", "include", "main tables, figures, mpce_summary.json"),
    ("analysis/mpce_robustness.py", "analysis", "include", "cubic-curve / Gaussian-wake re-evaluation (cached)"),
    ("analysis/mpce_numbers.py", "analysis", "include", "number macros mpce_numbers.tex"),
    ("analysis/mpce_inference_extra.py", "analysis", "include", "equivalence at three levels, inference robustness"),
    ("analysis/mpce_direction.py", "analysis", "include", "direction resolution, Horns Rev 1 re-evaluation"),
    ("analysis/mpce_csweep.py", "analysis", "include", "PSO coefficient sweep analysis"),
    ("analysis/make_theory_figures.py", "analysis", "include", "theory figure and checks T01-T13"),
    ("analysis/make_model_figures.py", "analysis", "include", "model schematics (figures_final/)"),
    ("analysis/analyze_authors_runs.py", "analysis", "include", "authors_tables.tex (spacing table) and calibration table"),
    ("analysis/rev2_analysis.py", "analysis", "include", "revision-2 statistics and tables (rev2_summary.json, rev2_tables.tex)"),
    ("analysis/audit_archive.py", "check", "include", "archive audit of the revision-2 records (expects the experiments/ + validation/ layout)"),
    ("analysis/check_revision.py", "check", "include", "focused checks of the revision-2 package"),
    ("analysis/mpce_check_*.py", "check", "include", "check scripts (C/X/D/F/S ids)"),
    ("analysis/README_reproduce.md", "doc", "include", "repository reproduction notes (MPCE/SWEVO pipeline)"),
    ("analysis/requirements.txt", "doc", "include", "environment of the original runs (superseded by the package requirements.txt)"),
    ("analysis/literature_values.*", "input", "include", "hand-collected literature values (Section S-literature)"),
    # ---- run records ------------------------------------------------------------------------------------------
    ("analysis/fresh_grid.csv", "runs-original", "include", "68 cases: SSA, LX-SSA, DE, old-setting PSO (+ unused SLSQP/VNS rows)"),
    ("analysis/fresh_vgrid.csv", "runs-original", "include", "68 cases: VNS"),
    ("analysis/fresh_bgrid.csv", "runs-original", "include", "68 cases: SSA-VNS, LX-SSA-VNS"),
    ("analysis/fresh_hr16.csv", "runs-superseded-read", "include", "old Horns Rev binning; read by mpce_results.py but every HR row is replaced by hrfix"),
    ("analysis/fresh_vhr16.csv", "runs-superseded-read", "include", "old Horns Rev binning; read, then replaced by hrfix"),
    ("analysis/fresh_bhr16.csv", "runs-superseded-read", "include", "old Horns Rev binning; read, then replaced by hrfix"),
    ("analysis/mpce_feasx_s*of8.csv", "runs-original-shards", "include", "the 8 shards concatenated into mpce_feas_s0of1.csv (provenance; not read by the analysis)"),
    ("analysis/mpce_*_s*of*.csv", "runs-original", "include", "per-run records of mpce_experiments.py / iea37_experiments.py"),
    ("analysis/rev2_*_s*of*.csv", "runs-rev2", "include", "revision-2 per-run records"),
    ("analysis/authors_runs_*.csv", "runs-earlier", "include", "re-runs of the original code (inputs of analyze_authors_runs.py: spacing table)"),
    ("analysis/ssa_reference_runs.csv", "runs-earlier", "include", "SSA radial-projection runs (read by mpce_results.py)"),
    # ---- inputs / caches / reference values ----------------------------------------------------------------
    ("analysis/iea37_published_results.csv", "input", "include", "published IEA37 CS1 layouts evaluated with the official calculator"),
    ("analysis/pywake_check*.csv", "input", "include", "PyWake NOJ reference values (read by mpce_results.py / mpce_direction.py)"),
    ("analysis/pywake_check*.json", "input", "include", "PyWake check of the stored Horns Rev runs"),
    ("analysis/mpce_reevaluation_cache.csv", "input", "include", "re-evaluation cache (key md5 of data set + coordinates)"),
    ("analysis/packing_capacity.csv", "output", "include", "packing bounds (input of mpce_results.py)"),
    ("analysis/calibration_vs_recorded.csv", "output", "include", "Mann-Whitney calibration table"),
    ("analysis/iea37_projected.*", "output", "include", "projected IEA37 layouts"),
    # ---- generated outputs -----------------------------------------------------------------------------------
    ("analysis/mpce_numbers*.tex", "output", "include", "generated number macros"),
    ("analysis/mpce_tab_*.tex", "output", "include", "generated main-text tables"),
    ("analysis/mpce_supp*.tex", "output", "include", "generated supplementary tables"),
    ("analysis/mpce_summary*.json", "output", "include", "generated summaries"),
    ("analysis/mpce_*.csv", "output", "include", "generated statistics"),
    ("analysis/mpce_direction_layouts.csv.gz", "output", "include", "1-degree re-evaluation of all final layouts"),
    ("analysis/authors_tables.tex", "output", "include", "spacing / calibration tables"),
    ("analysis/authors_spacing_all.csv", "output", "include", "output of analyze_authors_runs.py"),
    ("analysis/authors_runtime_6030.csv", "output", "include", "output of analyze_authors_runs.py"),
    ("analysis/equal_budget_tests.csv", "output", "include", "output of analyze_authors_runs.py"),
    ("analysis/hornsrev_tests.csv", "output", "include", "output of analyze_authors_runs.py"),
    # ---- logs ------------------------------------------------------------------------------------------------
    ("analysis/check_*.log", "log", "include", "stored outputs of the check scripts"),
    ("analysis/csweep.log", "log", "include", "log of mpce_csweep.py"),
    ("analysis/direction.log", "log", "include", "log of mpce_direction.py"),
    ("analysis/run_*.log", "log", "include", "launcher logs of the original experiment shards"),
    ("analysis/authors_experiments.log", "log", "include", "launcher log of authors_experiments.py"),
    # ---- figures -----------------------------------------------------------------------------------------------
    ("figures_mpce/*.pdf", "figure", "include", "generated figures (mpce_results.py, mpce_inference_extra.py, mpce_csweep.py, mpce_diagnostics.py, make_theory_figures.py)"),
    ("figures_mpce/*.png", "figure", "include", "PNG previews of the generated figures"),
    ("figures_final/fig_wind_farm.pdf", "figure", "include", "model schematic (make_model_figures.py)"),
    ("figures_final/fig_wake_model.pdf", "figure", "include", "model schematic (make_model_figures.py)"),
    ("figures_final/fig_half_cone.pdf", "figure", "include", "model schematic (make_model_figures.py)"),
    ("figures_final/*", "legacy", "exclude", "figures of an earlier paper"),
    # ---- revision-2 package documents ------------------------------------------------------------------------
    ("SWEVO_rev2/experiments_docs/rev2_summary.json", "output", "include", "revision-2 summary (rev2_analysis.py)"),
    ("SWEVO_rev2/experiments_docs/rev2_tables.tex", "output", "include", "revision-2 tables (rev2_analysis.py)"),
    ("SWEVO_rev2/experiments_docs/run_rev2_*.sh", "driver", "include", "revision-2 launch scripts (exact shard commands)"),
    ("SWEVO_rev2/validation/*", "output", "include", "revision-2 audit outputs (audit_archive.py, check_revision.py)"),
    # ---- root-level inputs -------------------------------------------------------------------------------------
    ("selected_30_run_data.csv", "runs-archived", "include", "720 archived runs of the original code (evaluator validation, calibration, boundary-rule check)"),
]

SCAN = ["analysis", "figures_mpce", "figures_final", "SWEVO_rev2/experiments_docs", "SWEVO_rev2/validation"]
EXTRA_FILES = ["selected_30_run_data.csv"]

# Revision-3 outputs other work packages were asked to produce (Section D / brief). Listed as "pending" until a
# matching file exists. Pattern -> what it is.
PENDING_HOOKS = [
    ("analysis/rev3_fullprec_*.csv", "full-precision (17 significant digits) side files of the legacy coordinates (C4)"),
    ("analysis/rev3_*_manifest.md", "prespecification of every revision-3 run set"),
    ("analysis/rev3_*.py", "revision-3 drivers and analyses"),
    ("analysis/rev3_*.csv", "revision-3 run records (seeds 31-60 unless stated)"),
    ("analysis/rev3_*.tex", "revision-3 tables (tab:S-r3-*)"),
    ("analysis/rev3_*.json", "revision-3 summaries"),
]


def classify(rel):
    for pat, cat, dec, note in RULES:
        if fnmatch.fnmatch(rel, pat):
            return cat, dec, note
    return "unclassified", "exclude", "no rule matched (check by hand)"


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--root", default=os.path.dirname(HERE))
    ap.add_argument("--out", default=HERE)
    a = ap.parse_args(argv)
    root = os.path.abspath(a.root)
    rows = []
    paths = []
    for d in SCAN:
        for dp, dns, fns in os.walk(os.path.join(root, d)):
            dns[:] = sorted(x for x in dns if x != "__pycache__")
            for fn in sorted(fns):
                paths.append(os.path.relpath(os.path.join(dp, fn), root))
    paths += [p for p in EXTRA_FILES if os.path.exists(os.path.join(root, p))]
    for rel in sorted(set(paths)):
        cat, dec, note = classify(rel)
        rows.append(dict(path=rel, category=cat, decision=dec, bytes=os.path.getsize(os.path.join(root, rel)), note=note))
    for pat, what in PENDING_HOOKS:
        if not any(fnmatch.fnmatch(r["path"], pat) for r in rows):
            rows.append(dict(path=pat, category="rev3-pending", decision="pending", bytes=0,
                             note=f"not present when the inventory was built: {what}"))
    os.makedirs(a.out, exist_ok=True)
    with open(os.path.join(a.out, "inventory.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["path", "category", "decision", "bytes", "note"])
        w.writeheader(); w.writerows(rows)
    inc = [r["path"] for r in rows if r["decision"] == "include"]
    with open(os.path.join(a.out, "package_files.txt"), "w") as fh:
        fh.write("\n".join(inc) + "\n")
    from collections import Counter
    c = Counter((r["decision"], r["category"]) for r in rows)
    for (dec, cat), n in sorted(c.items()):
        b = sum(r["bytes"] for r in rows if r["decision"] == dec and r["category"] == cat)
        print(f"{dec:8s} {cat:24s} {n:5d} files {b / 1e6:9.2f} MB")
    un = [r["path"] for r in rows if r["category"] == "unclassified"]
    if un:
        print("UNCLASSIFIED (excluded):", ", ".join(un))
    print(f"included: {len(inc)} files, {sum(r['bytes'] for r in rows if r['decision'] == 'include') / 1e6:.1f} MB")
    return 0


if __name__ == "__main__":
    sys.exit(main())

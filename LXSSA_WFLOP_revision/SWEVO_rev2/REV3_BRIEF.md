# Revision 3 of the SWEVO manuscript — brief for all agents

Basis: the authors' package `SWEVO_rev2/` (canonical LaTeX sources in `SWEVO_rev2/latex_source/`, revision notes
`SWEVO_rev2/REVISION_CHANGELOG.md`, `AUTHOR_CONFIRMATION.md`) and the remaining-work report
(text: /tmp/claude-0/-home-user-Figures-QASSA-WFLOP/83ceea2c-3e62-519f-bca6-caf988ef6cfa/scratchpad/rev_work.txt;
items R1–R10, A1–A8, B1–B7, C1–C4, D, E1–E6).

KEY FACT: the revision-2 authors believed the 36,190 original run records and the core modules were missing. They are
NOT missing: they are in `analysis/` of this repository (fresh_*.csv, mpce_*_s*of*.csv, wflop_model.py,
authors_objective.py, authors_optimizers.py, extra_baselines.py, feasible_init.py, hornsrev_model.py,
hybrid_lxssa_bvns.py, iea37_experiments.py, iea37_model.py, init_hook.py, mpce_experiments.py,
mpce_inference_extra.py, mpce_results.py, original_vns.py, rs_vns.py, check scripts, ...). The 11,640 revision-2 records
(analysis/rev2_*_s*of*.csv) and their drivers (analysis/rev2_*.py, record_io.py, audit_archive.py, check_revision.py)
were copied in as well. `python3 analysis/rev2_analysis.py --out-dir <scratch>` reproduces SWEVO_rev2's
rev2_summary.json exactly (4,670 values), so this environment matches the one of the runs.
Loaders: `mpce_inference_extra.load(HERE)` (main benchmark, 6,030 evaluations, random init, all 12 methods incl.
controls, without coordinates), `mpce_results.py` (load / std_cols / case_stats / rank_rule / paired_vs / tost /
case_mean_wilcoxon / HR and IEA37 loaders), `rev2_analysis.py` (read_rev2 and the blocks for GA, sites, gradients,
spacing). Reuse these functions; do not modify existing scripts (write new ones).

## Rules
- Write only NEW files named `analysis/rev3_<topic>*.{py,csv,json,tex,log}` (your topic is in your task), plus
  scratch files in /tmp/claude-0/-home-user-Figures-QASSA-WFLOP/83ceea2c-3e62-519f-bca6-caf988ef6cfa/scratchpad/rev3_<topic>/.
  Do NOT edit manuscript/supplement sources, existing scripts, existing CSVs. No git commit/push (the lead commits).
- CPU: the machine has 4 cores shared by several agents. Use at most the number of worker processes stated in your
  task; long jobs in the background with progress logs; never block on > 20 min without a progress check.
- New optimization runs: prespecify them first in `analysis/rev3_<topic>_manifest.md` (question, cases, methods,
  seeds, budget accounting, comparison family, margin, outcomes) before running; use NEW seeds 31–60 unless the task
  says otherwise; write coordinates with 17 significant digits (`record_io.encode_coordinates`); same CSV columns as
  the existing run files (Algorithm, Dataset, Radius, Turbines, Seed, Budget, Init, Calls, Objective, Ideal, WakeLoss,
  Feasible, MinSpacing, Seconds, Coordinates, Curve, + any extra columns at the end).
- Report ALL outcomes, including unfavorable ones. Never fabricate; if something cannot be done, say exactly why.
- Statistics conventions of the paper: wake loss in % of ideal, ΔL = first minus second (pp, negative = first
  better), margin m = 0.05 pp, Holm within a stated family, two-sided α = 0.05, 90% intervals for equivalence,
  fixed RNG seeds recorded in the output.
- LaTeX tables you produce: booktabs, \scriptsize or \footnotesize, fit \textwidth of the SWEVO preprint (A4, 2.5 cm
  margins, 12pt; check with a quick standalone compile in your scratch dir: 0 overfull), label `tab:S-r3-<name>`,
  caption that states cases, seeds, budget, endpoint, direction of differences and the method pool.
- Final report to the lead: files written, what each table/number shows (with the key numbers), the exact sentences
  you suggest for the manuscript (and where), limitations, and run-times.

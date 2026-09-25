# Analysis supporting the revised manuscript

All scripts run with Python 3.11, NumPy, SciPy and pandas, from inside this folder.
Input: `../selected_30_run_data.csv` (720 stored follow-up runs with coordinates).

| Script | Manuscript item | Output |
|---|---|---|
| `wflop_model.py` | Benchmark evaluator (Jensen/Weibull, Kusiak–Song discretization), cubic power curve, Gaussian wake, Horns Rev 1 wind rose (dataset 3) | library |
| `validate_evaluator.py` | Sec. VIII-A: reproduces all 720 recorded objectives (max abs. error 8.7e-11) | stdout |
| `authors_optimizers.py` | Authors' GA, PSO, DE, SSA, LX-SSA code (logic and random-call order unchanged) | library |
| `authors_objective.py` | Penalized minimization objective, vectorized; default `penalty="authors"` is the exact quadratic penalty of `objective.py` | library |
| `calibrate_authors_code.py` | Sec. VII-A: reruns vs. recorded follow-up runs (count vs. magnitude penalty) | `calibration_*.csv` |
| `authors_experiments.py` | Sec. VII and VIII-D: `budget` (equal budgets 3,030/6,030), `crowded`, `spacing` (5D/6D), `hornsrev`; about 30 min on 4 cores | `authors_runs_*.csv` |
| `analyze_authors_runs.py` | Sec. VII: calibration, seed-paired Friedman + Holm–Wilcoxon, feasibility, runtime, LaTeX tables | `authors_tables.tex`, `equal_budget_tests.csv`, `hornsrev_tests.csv`, `calibration_vs_recorded.csv`, `authors_runtime_6030.csv` |
| `followup_statistics.py` | Sec. VI-D: case-level Friedman + Holm, bootstrap CIs, A12, batch runtimes, best layouts | `friedman_case_level.csv`, `bootstrap_ci.csv`, `batch_runtime_per_run.csv`, `best_layouts_selected_cases.csv` |
| `reevaluate_layouts.py` | Sec. VIII-B/C: cubic power curve, Gaussian wake, 5D/6D check of stored layouts | `layout_reevaluation.csv`, `robustness_summary.csv` |
| `packing_capacity.py` | Sec. VIII-D: constructible turbine count per radius/spacing | `packing_capacity.csv` |
| `make_robustness_tables.py` | LaTeX tables of Sec. VIII-B/C/D (capacity) | `robustness_tables.tex` |
| `objective_original.py` | The authors' `objective.py` (verbatim); `authors_objective.py` matches it to < 1e-6 relative | library |
| `extra_baselines.py` | VNS and multistart SLSQP (R3-3) | library |
| `hornsrev_model.py` | Horns Rev 1 site model: V80 power/Ct curves, measured wind climate, real outline (cross-checked against PyWake NOJ within 0.3%) | library |
| `extra_experiments.py` | `extra` (VNS/SLSQP, ~11 min), `hr16` (~5 min), `hr80` (~25 min) | `extra_runs_*.csv` |
| `analyze_extra_runs.py` | Sec. VII-C/D/F tables | `extra_tables.tex`, `six_method_average_ranks.csv` |
| `build_section7.py` | Rebuilds Section VII of the manuscript from the generated tables | – |
| `ssa_reference.py` | Sec. VII-H only: SSA with radial projection (boundary-handling comparison) | `ssa_reference_runs.csv` |

`robustness_text.tex` holds the original draft text of Section VIII; the current text is in the main manuscript.

Caveats:
- LX-SSA runs use the authors' code with phi = 0, chi = 1. At the recorded budgets the reruns reproduce the recorded follow-up distributions (22 of 24 Mann–Whitney p ≥ 0.05).
- All controlled runs use the authors' exact penalty (`penalty="authors"`); earlier runs with a reconstructed linear penalty are kept in `superseded_linear_penalty/` for reference only.
- Horns Rev 1 wind rose (12 sectors, Weibull A/k/f) is taken from DTU PyWake 2.6.20 (MIT licence), `py_wake/examples/data/hornsrev1.py`.

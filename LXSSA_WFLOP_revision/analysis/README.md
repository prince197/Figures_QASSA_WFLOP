# Analysis supporting the revised manuscript

## Hybrid LX-SSA-VNS study (current manuscript, Sections VI–XI)

| Script | Purpose | Output |
|---|---|---|
| `hybrid_lxssa_vns.py` | The proposed hybrid (LX-SSA phase + VNS phase, split rho) | library |
| `full_grid_experiments.py hgrid` | 68 cases × LX-SSA-VNS × 30 seeds (paired with `fresh_grid.csv`) | `fresh_hgrid.csv` |
| `full_grid_experiments.py hsplit` | Split ablation: rho = 25% / 75% on 12 cases | `fresh_hsplit.csv` |
| `full_grid_experiments.py hhr16` / `hhr80` | Horns Rev 1 with the hybrid | `fresh_hhr16.csv`, `fresh_hhr80.csv` |
| `hybrid_results.py` | Eight-method tables, statistics, ablation and figures | `hybrid_tables.tex`, `hybrid_summary.json`, `hybrid_case_tests.csv`, `hybrid_best_layouts_maxN.csv`, `../figures_hybrid/*.pdf` |
| `robustness_hybrid.py` | Cubic power curve / Gaussian wake re-evaluation, eight methods | `hybrid_robust_table.tex`, `hybrid_robust_summary.json`, `hybrid_reevaluation.csv` |
| `original_vns.py` | Original basic VNS (BVNS; Mladenovic & Hansen 1997, continuous form of Mladenovic et al. 2008): l_inf-shell shaking of the whole layout + complete best-improvement compass local search | library |
| `full_grid_experiments.py vgrid` / `vhr16` / `vhr80` | BVNS on the 68 cases and Horns Rev (paired seeds) | `fresh_vgrid.csv`, `fresh_vhr16.csv`, `fresh_vhr80.csv` |
| `bvns_compare.py` | BVNS vs each method (Holm over seven per case) | `bvns_tables.tex`, `bvns_summary.json`, `bvns_case_tests.csv` |
| `sec_hybrid.tex` | Section VI (algorithm) | – |
| `build_hybrid_manuscript.py` | Writes Sections VI–XII and the appendix from the outputs above | – |
| `hybrid_frontmatter.py` | Title, abstract, introduction and motivation section | – |

## Six-method study (reference runs reused by the hybrid study)

| Script | Purpose | Output |
|---|---|---|
| `full_grid_experiments.py grid` | 68 cases × 6 methods × 30 seeds, 6,030 calls, convergence recording (~2.2 h on 4 cores) | `fresh_grid.csv` |
| `full_grid_experiments.py hr16` / `hr80` | Horns Rev 1 site case, 16 and 80 turbines | `fresh_hr16.csv`, `fresh_hr80.csv` |
| `fresh_results.py` | All result tables, statistics and figures | `fresh_tables.tex`, `fresh_summary.json`, `fresh_case_tests.csv`, `fresh_best_layouts_maxN.csv`, `../figures_fresh/*.pdf` |
| `robustness_fresh.py` | Cubic power curve / Gaussian wake re-evaluation of all fresh layouts | `fresh_robust_table.tex`, `fresh_reevaluation.csv` |
| `build_fresh_manuscript.py` | Writes Sections VI–XI and the appendix of the manuscript from the outputs above | – |
| `sec_setup.tex` | Text of the Experimental Setup section | – |

Libraries used by the fresh study: `authors_optimizers.py`, `authors_objective.py` (exact `objective.py` penalty), `objective_original.py`, `extra_baselines.py` (VNS, MS-SLSQP), `wflop_model.py`, `hornsrev_model.py`. The spacing table uses `authors_runs_spacing.csv` (from `authors_experiments.py spacing`) and the capacity table `packing_capacity.csv`.

## Earlier analyses (superseded; not used in the current manuscript)

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

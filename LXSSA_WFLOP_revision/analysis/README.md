# Analysis supporting the revised manuscript

## MPCE resubmission (current): proposed method PSO-VNS (`PSOBV`)

| Script | Purpose | Output |
|---|---|---|
| `mpce_experiments.py`, `iea37_experiments.py` | Runs (see their docstrings for experiments and method ids) | `mpce_<exp>_s<i>of<k>.csv` |
| `mpce_results.py [--focus PSOBV\|PSOC\|SSABV] [--partial] [--common-seeds] [--skip-robust]` | All tables, statistics and figures; every section follows `--focus` (default `PSOBV`); the decision block (`summary["decision"]`, printed last) compares PSO-VNS with PSO and SSA-VNS regardless of the focus | `mpce_tab_*.tex`, `mpce_supplementary.tex`, `mpce_summary.json`, `mpce_*_tests.csv`, `../figures_mpce/*` |
| `mpce_robustness.py` | Cubic power curve / Gaussian wake re-evaluation of the main-comparison layouts (cached in `mpce_reevaluation_cache.csv`) | `mpce_tab_robust.tex`, `mpce_reevaluation.csv` |

| `mpce_numbers.py` | Every number / data-dependent phrase quoted in `MPCE_PSO_VNS.tex` as `\N...` macros (letters only); a macro whose experiments are missing or incomplete expands to `\TBD{pending: <experiments>}`; run automatically at the end of `mpce_results.py` | `mpce_numbers.tex` |
| `mpce_check_final.py` | Evaluates the conditions behind the outcome-dependent statements of the manuscript (comments `% CHECK-FINAL [Cnn]` in `MPCE_PSO_VNS.tex`) and prints PASS / FAIL / PENDING; run automatically at the end of `mpce_results.py` | stdout |
| `build_mpce_paper.sh [--partial]` | **One command after new shards arrive**: `mpce_results.py` (+ numbers + checks), then `pdflatex` of `MPCE_PSO_VNS_supplement.tex` and `MPCE_PSO_VNS.tex` in alternation (cross-references via `xr`: the main text reads the supplement's labels with prefix `S-`, the supplement the main text's with prefix `M-`) | `../MPCE_PSO_VNS.pdf`, `../MPCE_PSO_VNS_supplement.pdf` |

Main-text tables written by `mpce_results.py` in the manuscript layout: `mpce_tab_friedman68.tex` (single column), `mpce_tab_wtl.tex`, `mpce_tab_ablation.tex` (eight variants, eight contrasts), `mpce_tab_split.tex` (compact), `mpce_tab_feasbudget.tex` (six largest cases: random vs feasible initialization, ranks at 6,030 / 30,030 / 120,030 calls), `mpce_tab_hr16.tex` (ten methods, loss at 6k R / 6k F / 30k / 120k), `mpce_tab_iea37.tex` (best runs vs published), `mpce_tab_robust_final.tex` (single column). A column whose runs are incomplete prints `\TBD{}` cells. All other tables (per-case results S1--S6, spacing re-optimization S7, full Friedman table, cost, per-case split / initialization / budget / Horns Rev / IEA37 / robustness tables, best-layout coordinates) go to `mpce_supplementary.tex`, which `../MPCE_PSO_VNS_supplement.tex` inputs. Main-text figures: `avg_ranks`, `wakeloss_vs_n`, `ablation_convergence`, `budget_scaling`, `layouts_iea37_hr16` (IEA37 16/36 + Horns Rev 16 layouts).

Method ids: `PSOBV` = PSO-VNS (constriction-PSO phase 1 for half the budget, then basic VNS), `PSOC` = PSO with Clerc–Kennedy constriction coefficients, `SSABV` = SSA-VNS, `LXBV` = LX-SSA-VNS, `RSVNS` = random-sampling phase 1 + VNS, `BVNS` = basic VNS, `SLSQP` = MS-SLSQP; split variants `PSOBV25` / `PSOBV75`.

## Final study (previous LXSSA manuscript): hybrid LX-SSA-VNS = LX-SSA + original basic VNS

| Script | Purpose | Output |
|---|---|---|
| `original_vns.py` | Original basic VNS (Mladenovic & Hansen 1997; continuous form of Mladenovic et al. 2008) | library |
| `hybrid_lxssa_bvns.py` | Final hybrid (phase 1 LX-SSA or SSA, phase 2 original VNS) | library |
| `full_grid_experiments.py vgrid` / `vhr16` / `vhr80` | Original VNS on the 68 cases and Horns Rev | `fresh_vgrid.csv`, `fresh_vhr16.csv`, `fresh_vhr80.csv` |
| `full_grid_experiments.py bgrid` | LX-SSA-VNS (`LXBV`) and SSA-VNS (`SSABV`, not reported in the manuscript) on the 68 cases | `fresh_bgrid.csv` |
| `full_grid_experiments.py bsplit` | Budget split 25% / 75% on 12 cases | `fresh_bsplit.csv` |
| `full_grid_experiments.py bhr16` / `bhr80` | Hybrids on Horns Rev | `fresh_bhr16.csv`, `fresh_bhr80.csv` |
| `final_results.py` | Seven-method tables, statistics, component ablation (SSA, LX-SSA, VNS, LX-SSA-VNS), split table, figures | `final_tables.tex`, `final_summary.json`, `final_case_tests.csv`, `final_ablation_tests.csv`, `final_best_layouts_maxN.csv`, `../figures_final/*.pdf` |
| `robustness_final.py` | Cubic power curve / Gaussian wake re-evaluation, seven methods | `final_robust_table.tex`, `final_robust_summary.json`, `final_reevaluation.csv` |
| `sec_hybrid.tex`, `sec_setup.tex` | Section VI (algorithm) and the Experimental Setup source | – |
| `build_final_manuscript.py`, `final_frontmatter.py` | Write Sections VI–XII, appendix, abstract/contributions/motivation | – |
| `final_prose.py` | Final prose pass in the authors' style (run once after the two builders); switches Figs. 1–3 to the redrawn figures | – |
| `make_model_figures.py` | Redrawn model figures (circular farm, Jensen wake, wake half cone) | `../figures_final/fig_wind_farm.pdf`, `fig_wake_model.pdf`, `fig_half_cone.pdf` |

Method ids in the CSV files: `LXBV` = LX-SSA-VNS (final), `SSABV` = SSA-VNS (run, not reported), `BVNS` = original VNS (shown as "VNS"), `LXSSA`, `SSA`, `PSO`, `DE`, `SLSQP` = MS-SLSQP. The ids `VNS` (modified VNS) and `LXVNS` (earlier hybrid with the modified VNS) in `fresh_grid.csv` / `fresh_hgrid.csv` belong to superseded versions and are not used in the manuscript.

## Superseded intermediate versions (kept for reference, not used in the manuscript)

`hybrid_lxssa_vns.py`, `hybrid_results.py`, `robustness_hybrid.py`, `build_hybrid_manuscript.py`, `hybrid_frontmatter.py`, `bvns_compare.py` and their outputs (`hybrid_*`, `bvns_*`, `fresh_hgrid.csv`, `fresh_hsplit.csv`, `fresh_hhr*.csv`).

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

# Reproducing the results of MPCE_PSO_VNS (PSO-VNS paper)

All commands run from `analysis/` with the environment pinned in `requirements.txt`
(Python 3.11.15, NumPy 2.4.6, SciPy 1.17.1, pandas 3.0.6, Matplotlib 3.11.2). MS-SLSQP (SciPy SLSQP) is
platform dependent; all other methods reproduce bitwise for a given seed.

## 1. Per-run data: which command produced which CSV

Every run file has one row per run (method, case, seed, budget, initialization, evaluations `Calls`, final
objective, feasibility, minimum spacing, wall-clock `Seconds`, final coordinates, convergence curve). The
`fresh_*.csv` files have no `Budget` / `Init` columns (all 6,030 evaluations, random initialization; filled in by
`mpce_results.py`); their Horns Rev files use `AEP` / `IdealAEP` instead of `Objective` / `Ideal`.
Convergence checkpoints: every (B-30)//200 evaluations (201 values at 6,030, 200 at 30,030 / 120,030).

| CSV | Command | Content |
|---|---|---|
| `fresh_grid.csv` | `python3 full_grid_experiments.py grid` | 68 cases, SSA, LX-SSA, DE, old-setting PSO (`PSO`, w = 0.7, c1 = c2 = 2); its `SLSQP` and `VNS` rows are not used |
| `fresh_vgrid.csv` | `python3 full_grid_experiments.py vgrid` | 68 cases, basic VNS (`BVNS`) |
| `fresh_bgrid.csv` | `python3 full_grid_experiments.py bgrid` | 68 cases, LX-SSA-VNS (`LXBV`), SSA-VNS (`SSABV`) |
| `fresh_hr16.csv`, `fresh_vhr16.csv`, `fresh_bhr16.csv` | `full_grid_experiments.py hr16 / vhr16 / bhr16` | Horns Rev 1, 16 turbines, old direction binning; **superseded by `hrfix`** (see below) |
| `mpce_rsvns_s0of1.csv` | `python3 mpce_experiments.py rsvns` | 68 cases, RS-VNS |
| `mpce_psoc_s0of1.csv` | `python3 mpce_experiments.py psoc` | 68 cases, PSO (constriction) |
| `mpce_psobv_s<i>of2.csv` | `python3 mpce_experiments.py psobv <i> 2` | 68 cases + Horns Rev 16 (HR rows superseded by `hrfix`), PSO-VNS |
| `mpce_slsqp_s0of1.csv` | `python3 mpce_experiments.py slsqp` | 68 cases, MS-SLSQP (rerun on the current platform) |
| `mpce_psosplit_s0of1.csv` | `python3 mpce_experiments.py psosplit` | 12 split cases, PSO-VNS with omega = 0.25 / 0.75 (omega = 1 is `mpce_psoc`) |
| `mpce_hr16new_s0of1.csv` | `python3 mpce_experiments.py hr16new` | Horns Rev 16, PSO and RS-VNS (superseded by `hrfix`) |
| `mpce_feas_s0of1.csv` | `python3 mpce_experiments.py feasx <i> 8` for i = 0..7, then the 8 shards `mpce_feasx_s<i>of8.csv` concatenated into `mpce_feas_s0of1.csv` (1,680 rows, identical objectives / coordinates) | feasible initialization, six largest cases + Horns Rev 16, the 8 methods of M9 without RS-VNS. The documented `feas` experiment (9 methods incl. RS-VNS) was **not** run as such: RS-VNS would need one packing solve per random sample. `mpce_results.py` reads only `mpce_feas_s*`, never the `feasx` shards |
| `mpce_feasp_s0of1.csv` | `python3 mpce_experiments.py feasp` | PSO-VNS arm of the feasible-initialization study |
| `mpce_b30k_s<i>of3.csv`, `mpce_b30kp_s0of1.csv` | `mpce_experiments.py b30k <i> 3`, `b30kp` | 30,030 evaluations, six largest cases + Horns Rev 16 (M9, PSO-VNS arm) |
| `mpce_b120k_s<i>of8.csv`, `mpce_b120kp_s<i>of4.csv` | `mpce_experiments.py b120k <i> 8`, `b120kp <i> 4` | 120,030 evaluations (Horns Rev 16: 10 seeds) |
| `mpce_iea16_s0of1.csv`, `mpce_iea36_s0of1.csv`, `mpce_iea16p_s0of1.csv`, `mpce_iea36p_s0of1.csv` | `python3 iea37_experiments.py iea16` / `iea36` / `iea16p` / `iea36p` | IEA37 Case Study 1, 16 and 36 turbines, 6,030 and 30,030 evaluations |
| `mpce_hrfix_s<i>of24.csv` | `python3 mpce_experiments.py hrfix <i> 24` for i = 0..23 | **All** Horns Rev 1 16-turbine runs after the direction-binning fix of `hornsrev_model.py` (commit 7676da9): 10 methods at 6,030 random / 6,030 feasible (no RS-VNS) / 30,030 (30 seeds) / 120,030 (10 seeds) |
| `iea37_published_results.csv` | built from the IEA37 repository (`iea37_data/`) | published CS1 layouts evaluated with the official calculator |
| `ssa_reference_runs.csv` | `python3 ssa_reference.py` | SSA with radial boundary projection (boundary-rule check, robustness section) |
| `../selected_30_run_data.csv` | archived runs of the earlier study | original SSA / LX-SSA / PSO / DE runs (evaluator validation, calibration, boundary-rule check) |

Horns Rev loading rule (`mpce_results.py`, `load`): as soon as any `mpce_hrfix_s*.csv` exists, every Horns Rev
row of every other file is dropped and all Horns Rev results come from `hrfix` (incomplete shards are used, and
`mpce_numbers.py` prints the Horns Rev macros as `\TBD{pending: hrfix}` until all 24 shards are present). Without
any `hrfix` file the old runs are used together with the installed AEP of the old model (139.51 GWh/yr), so that
runs and reference come from the same model.

## 2. Tables, figures and numbers: which script produces what

`sh build_mpce_paper.sh [--partial]` runs `python3 mpce_results.py` (which calls `mpce_numbers.py` and
`mpce_check_final.py`) and then pdflatex of the supplement and the main text.

| Output | Produced by |
|---|---|
| `mpce_tab_baseline.tex` (Table baseline), `mpce_tab_friedman68.tex`, `mpce_tab_wtl.tex`, `mpce_tab_ablation.tex`, `mpce_tab_split.tex`, `mpce_tab_feasbudget.tex`, `mpce_tab_hr16.tex`, `mpce_tab_iea37.tex`, `mpce_tab_robust_final.tex` | `mpce_results.py` |
| `mpce_supplementary.tex` (per-case results, full Friedman table, cost, component analysis per case, switch point, split per case, Horns Rev, feasible initialization, budget, IEA37, robustness, provenance, baseline per case; copies of the spacing and packing tables) | `mpce_results.py` |
| `mpce_summary.json` (every quoted number), `mpce_case_stats.csv`, `mpce_case_tests.csv`, `mpce_ablation_tests.csv`, `mpce_feasinit_tests.csv`, `mpce_budget_case_stats.csv`, `mpce_best_layouts_maxN.csv`, `mpce_reevaluation.csv` | `mpce_results.py` |
| `../figures_mpce/*.pdf` (average ranks, wake loss / feasibility vs N, convergence, box plots, layouts, ablation convergence, budget scaling, Horns Rev, IEA37) | `mpce_results.py` |
| robustness re-evaluation (cubic power curve, Gaussian wake) | `mpce_robustness.py` (called by `mpce_results.py`; cached in `mpce_reevaluation_cache.csv` -- delete the cache to recompute) |
| `mpce_numbers.tex` (`\N...` macros) | `mpce_numbers.py` from `mpce_summary.json` |
| CHECK-FINAL PASS / FAIL list | `mpce_check_final.py` (reads `mpce_summary.json` and the `% CHECK-FINAL [Cnn]` comments of `../MPCE_PSO_VNS.tex` and `../optA/*.tex`) |

Statistics conventions (see the docstrings of `mpce_results.py`): run-level Wilcoxon signed-rank tests on 30
seed-paired runs, Holm-adjusted within each case over the comparisons of that case (the family differs between
tables, e.g. PSO-VNS vs PSO: 7 comparisons in Table wtl, the component-analysis contrasts in Table ablation);
case-mean Wilcoxon tests on the per-case mean wake losses (a case where only one method has at least 15 feasible
runs counts as a maximal difference in its favour, cases where neither has are dropped); 95 % bootstrap CIs:
10,000 resamples of the cases, seed 20260928; zero differences are dropped in the Wilcoxon test and the
rank-biserial correlation.

## 3. Not regenerated by `build_mpce_paper.sh`

| Item | Script / source |
|---|---|
| Model schematics (`../figures_final/fig_wind_farm.pdf`, `fig_wake_model.pdf`, `fig_half_cone.pdf`) | `python3 make_model_figures.py` |
| Minimum-spacing re-optimization table (`tab:spacing-authors`, copied from `authors_tables.tex`) | `python3 authors_experiments.py spacing` -> `authors_runs_spacing.csv`, then `python3 analyze_authors_runs.py` -> `authors_tables.tex` |
| Packing bounds (`tab:capacity`, `packing_capacity.csv`) | `python3 packing_capacity.py` |
| Evaluator checks (720 archived objectives within 8.7e-11 absolute; 22 of 24 Mann-Whitney p >= 0.05) | `python3 validate_evaluator.py`; `python3 calibrate_authors_code.py SSA,LXSSA <penalty>` -> `calibration_vs_recorded.csv` |
| Horns Rev 1 validation against PyWake (80 turbines; `\NHRPyWakeDiff`) | the model side is recomputed by `mpce_results.py`; the PyWake reference (662.5 GWh/yr) is a constant (PyWake is not installed here) |
| Literature comparison table | `literature_values.csv` / `literature_values.md` (collected by hand) |
| Per-run data (Section 1) | the experiment scripts above (hours of CPU time) |

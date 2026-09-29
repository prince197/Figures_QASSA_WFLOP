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
| `mpce_omega90_s0of1.csv` | `python3 mpce_experiments.py omega90 0 1` | Phase 6: 12 split cases, PSO-VNS with omega = 0.9 (label `PSOBV90`; 5,430 PSO evaluations, then VNS), 30 seeds |
| `mpce_rsdisc_s0of1.csv` | `python3 mpce_experiments.py rsdisc 0 1` | Phase 6: 68 cases, RSD-VNS (label `RSDVNS`): RS-VNS whose 2,985 Phase-1 samples after the common initial population are uniform in the farm disc (3,015 Phase-1 evaluations), 30 seeds |
| `mpce_hr16new_s0of1.csv` | `python3 mpce_experiments.py hr16new` | Horns Rev 16, PSO and RS-VNS (superseded by `hrfix`) |
| `mpce_feas_s0of1.csv` | `python3 mpce_experiments.py feasx <i> 8` for i = 0..7, then the 8 shards `mpce_feasx_s<i>of8.csv` concatenated in shard order with `python3 -c "import glob, pandas as pd; pd.concat([pd.read_csv(f) for f in sorted(glob.glob('mpce_feasx_s*of8.csv'))], ignore_index=True).to_csv('mpce_feas_s0of1.csv', index=False)"` (1,680 rows; this command reproduces the stored file byte for byte -- the default pandas float parser drops the last digit of some floats, below 1e-15 relative) | feasible initialization, six largest cases + Horns Rev 16, the 8 methods of M9 without RS-VNS. The documented `feas` experiment (9 methods incl. RS-VNS) was **not** run as such: RS-VNS would need one packing solve per random sample. `mpce_results.py` reads only `mpce_feas_s*`, never the `feasx` shards |
| `mpce_feasp_s0of1.csv` | `python3 mpce_experiments.py feasp` | PSO-VNS arm of the feasible-initialization study |
| `mpce_b30k_s<i>of3.csv`, `mpce_b30kp_s0of1.csv` | `mpce_experiments.py b30k <i> 3`, `b30kp` | 30,030 evaluations, six largest cases + Horns Rev 16 (M9, PSO-VNS arm) |
| `mpce_b120k_s<i>of8.csv`, `mpce_b120kp_s<i>of4.csv` | `mpce_experiments.py b120k <i> 8`, `b120kp <i> 4` | 120,030 evaluations (Horns Rev 16: 10 seeds) |
| `mpce_iea16_s0of1.csv`, `mpce_iea36_s0of1.csv`, `mpce_iea16p_s0of1.csv`, `mpce_iea36p_s0of1.csv` | `python3 iea37_experiments.py iea16` / `iea36` / `iea16p` / `iea36p` | IEA37 Case Study 1, 16 and 36 turbines, 6,030 and 30,030 evaluations |
| `mpce_hrfix_s<i>of24.csv` | `python3 mpce_experiments.py hrfix <i> 24` for i = 0..23 | **All** Horns Rev 1 16-turbine runs after the direction-binning fix of `hornsrev_model.py` (commit 7676da9): 10 methods at 6,030 random / 6,030 feasible (no RS-VNS) / 30,030 (30 seeds) / 120,030 (10 seeds) |
| `iea37_published_results.csv` | built from the official files in `iea37_data/` (Section 1a) | published CS1 layouts (baseline + participants 1-12, 16 and 36 turbines) evaluated with the official calculator; optional for `mpce_results.py` (without it C40 and the published-layout comparisons are PENDING) |
| `ssa_reference_runs.csv` | `python3 ssa_reference.py` | SSA with radial boundary projection (boundary-rule check, robustness section) |
| `../selected_30_run_data.csv` | archived runs of the earlier study | original SSA / LX-SSA / PSO / DE runs (evaluator validation, calibration, boundary-rule check) |

### 1a. IEA37 published results (`iea37_published_results.csv`)

Source: the IEA Wind Task 37 case-study repository (https://github.com/IEAWindTask37/iea37-wflo-casestudies, folder
`cs1-2`, commit 267f6e5), copied unmodified into `iea37_data/` (`iea37-ex16/36.yaml`, `iea37-cs1-results/iea37-par<p>-opt<n>.yaml`,
`iea37-aepcalc.py`, `iea37-335mw.yaml`, `iea37-windrose.yaml`). The file was assembled once by the IEA37 agent; the
numeric columns (and the columns `mpce_results.load_published` reads: `Turbines`, `Participant`,
`AEP_MWh_official_calc`, `Feasible_tol1e-3m`) are recomputed and compared with the stored file by

```
python3 - <<'PY'
import importlib.util, os, numpy as np, pandas as pd
import iea37_model as M
sp = importlib.util.spec_from_file_location("calc", os.path.join(M.DATA, "iea37-aepcalc.py")); C = importlib.util.module_from_spec(sp); sp.loader.exec_module(C)
wd, wf, ws = C.getWindRoseYAML(os.path.join(M.DATA, "iea37-windrose.yaml"))
ci, co, rws, rp, dia = C.getTurbAtrbtYAML(os.path.join(M.DATA, "iea37-335mw.yaml"))
rows = []
for n in (16, 36):
    for par, fn in [("baseline (example layout)", os.path.join(M.DATA, f"iea37-ex{n}.yaml"))] + \
                   [(f"par{p}", os.path.join(M.DATA, "iea37-cs1-results", f"iea37-par{p}-opt{n}.yaml")) for p in range(1, 13)]:
        xy, rep, _ = M.load_layout(fn); tc, _, _ = C.getTurbLocYAML(fn)
        ex = float(np.sqrt((xy ** 2).sum(1)).max() - M.RADIUS[n])
        rows.append(dict(Turbines=n, Participant=par, AEP_MWh_reported=rep,
                         AEP_MWh_official_calc=round(float(C.calcAEP(tc, wf, ws, wd, dia, ci, co, rws, rp).sum()), 5),
                         AEP_MWh_iea37_model=M.aep(xy), MinSpacing_m=M.min_spacing(xy), MaxBoundaryExcess_m=ex,
                         **{"Feasible_tol1e-3m": bool(M.min_spacing(xy) >= M.SMIN - 1e-3 and ex <= 1e-3)}))
N, O = pd.DataFrame(rows), pd.read_csv("iea37_published_results.csv")
for c in N.columns[2:]:
    print(c, "max rel. diff", np.max(np.abs(N[c].astype(float).values - O[c].astype(float).values) / np.maximum(1, np.abs(O[c].astype(float).values))))
PY
```

(all differences <= 2e-16 on 2026-09-28). The descriptive columns (`Case`, `BoundaryRadius_m`, `Algorithm`, `Source`
URL, `AEP_MWh_boundary_projected`, `RankAmongParticipants`, `Notes`) are not used by the analysis.

Horns Rev loading rule (`mpce_results.py`, `load`): as soon as any `mpce_hrfix_s*.csv` exists, every Horns Rev
row of every other file is dropped and all Horns Rev results come from `hrfix` (incomplete shards are used, and
`mpce_numbers.py` prints the Horns Rev macros as `\TBD{pending: hrfix}` until all 24 shards are present). Without
any `hrfix` file the old runs are used together with the installed AEP of the old model (139.51 GWh/yr), so that
runs and reference come from the same model.

## 2. Tables, figures and numbers: which script produces what

`sh build_mpce_paper.sh [--partial]` runs `python3 mpce_results.py` (which calls `mpce_numbers.py` and
`mpce_check_final.py`), then `python3 mpce_inference_extra.py`, `python3 mpce_check_extra.py` and
`python3 mpce_check_diag.py` (logs `check_extra.log`, `check_diag.log`, `check_final.log`; a FAIL is reported at the end
but does not stop the build), and then pdflatex of the supplement and the main text. The diagnostics and the theory
checks are separate steps (Section 3).

| Output | Produced by |
|---|---|
| `mpce_tab_baseline.tex` (Table baseline), `mpce_tab_friedman68.tex`, `mpce_tab_wtl.tex`, `mpce_tab_ablation.tex`, `mpce_tab_split.tex`, `mpce_tab_feasbudget.tex`, `mpce_tab_hr16.tex`, `mpce_tab_iea37.tex`, `mpce_tab_robust_final.tex` | `mpce_results.py` |
| `mpce_supplementary.tex` (per-case results, full Friedman table, cost, component analysis per case, switch point, split per case, Horns Rev, feasible initialization, budget, IEA37, robustness, provenance, baseline per case; copies of the spacing and packing tables) | `mpce_results.py` |
| `mpce_summary.json` (every quoted number), `mpce_case_stats.csv`, `mpce_case_tests.csv`, `mpce_ablation_tests.csv`, `mpce_feasinit_tests.csv`, `mpce_budget_case_stats.csv`, `mpce_best_layouts_maxN.csv`, `mpce_reevaluation.csv` | `mpce_results.py` |
| `../figures_mpce/*.pdf` (average ranks, wake loss / feasibility vs N, convergence, box plots, layouts, ablation convergence, budget scaling, Horns Rev, IEA37) | `mpce_results.py` |
| robustness re-evaluation (cubic power curve, Gaussian wake) | `mpce_robustness.py` (called by `mpce_results.py`; cached in `mpce_reevaluation_cache.csv` -- delete the cache to recompute) |
| `mpce_numbers.tex` (`\N...` macros) | `mpce_numbers.py` from `mpce_summary.json` |
| CHECK-FINAL PASS / FAIL list | `mpce_check_final.py` (reads `mpce_summary.json` and the `% CHECK-FINAL [Cnn]` comments of `../MPCE_PSO_VNS.tex`, `../MPCE_PSO_VNS_supplement.tex` and `../optA/*.tex`; it also lists the ids referenced by `% CHECK-EXTRA [Xnn]`, `% CHECK-DIAG [Dnn]`, `% CHECK-THEORY [Tnn]` and flags ids their script does not define, without evaluating them) |
| `mpce_numbers_extra.tex` (`\NX...`), `mpce_supp_inference.tex`, `mpce_summary_extra.json` | `mpce_inference_extra.py` (reads the per-run CSVs and `mpce_summary.json`; check X01 = it reproduces `mpce_summary.json`) |
| CHECK-EXTRA PASS / FAIL list (X01...) | `mpce_check_extra.py` (reads `mpce_summary_extra.json`) |
| CHECK-DIAG PASS / FAIL list (D01...) | `mpce_check_diag.py` (reads the stored `mpce_summary_diag.json`; the diagnostics are not rerun by the build) |

Switch point of the two-phase variants (`mpce_results.switch_call`): swarm hybrids switch after
30 + round((0.5 x 6,030 - 30) / c) x c Phase-1 evaluations (3,030 for PSO-VNS, SSA-VNS, LX-SSA-VNS); RS-VNS and RSD-VNS after
round(0.5 x 6,030) = 3,015 (`rs_vns.RSVNS.n1`). Feasibility and loss at the switch are read at the last convergence
checkpoint at or before it (call 3,030, resp. 3,000 for RS/RSD; check C61). Before review round 2 RS/RSD were read at
3,030 (R3-9).

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
| Evaluator checks (720 archived objectives within 8.7e-11 absolute; 22 of 24 Mann-Whitney p >= 0.05) | `python3 validate_evaluator.py` (reads `../selected_30_run_data.csv`; prints `records 720`, largest absolute objective difference 8.73e-11); `python3 calibrate_authors_code.py SSA,LXSSA <penalty>` -> `calibration_vs_recorded.csv` (Mann-Whitney tests of re-run vs. recorded objectives). The numbers are typed in 05_setup.tex, not generated as macros |
| IEA37 calculator check (relative AEP difference < 1e-11 against the official example layouts) | `python3 iea37_model.py` (prints `rel.err` -8.8e-12 / -2.8e-12 / 2.9e-12 for ex16 / ex36 / ex64 and the largest per-direction difference, 5e-6 MWh); the published layouts: Section 1a |
| Horns Rev 1 validation against PyWake (80 turbines and 16-turbine block; `\NHRPyWakeDiff`, `\NFPyWake...`) | `python3 pywake_check.py` (requires `pip install py_wake==2.6.20`, optional; not in requirements.txt) runs PyWake NOJ (k = 0.04) at bins identical to `hornsrev_model.py` and writes `pywake_check.csv` (+ `pywake_check_hr16runs.csv/.json` with `--layouts`); the CSV is committed, and `mpce_results.py` (`pywake_reference()`) and `mpce_direction.py` read it. Result: our AEP is 0.81 % higher (C_T evaluated at the free-stream speed; with local C_T our model reproduces PyWake hub-centre NOJ exactly). The former constant 662.5 GWh/yr (quoted from the previous study) is no longer used |
| Direction resolution (1° re-evaluation of all final layouts; Horns Rev 1°; projected IEA37 layouts; `\NF...`, `mpce_supp_direction.tex`) | `python3 mpce_direction.py --procs 2` (~75 s; also runs `iea37_projected.py`), then `python3 mpce_check_dir.py` (F01...) — run by the build |
| PSO coefficient sweep and bound handling (`mpce_csweep_s0of1.csv`; `\NS...`, `mpce_supp_csweep.tex`, `../figures_mpce/csweep.pdf`) | data: `python3 mpce_experiments.py csweep 0 1 --procs=4` (~2.6 CPU-h); analysis: `python3 mpce_csweep.py`, `python3 mpce_check_csweep.py` (S01...) — analysis run by the build |
| Robustness re-evaluation (cubic power curve, Gaussian wake; `mpce_reevaluation.csv`, `mpce_tab_robust_final.tex`) | `mpce_robustness.py` via `mpce_results.py` reads the cache `mpce_reevaluation_cache.csv` (key: md5 of data set + coordinate string) and evaluates only layouts not yet in it. Full recomputation: `rm mpce_reevaluation_cache.csv && python3 mpce_results.py` (the cache is then re-seeded from `final_reevaluation.csv` of the previous version, and every other layout is re-evaluated) |
| Literature comparison table | `literature_values.csv` / `literature_values.md` (collected by hand) |
| Per-run data (Section 1) | the experiment scripts above (hours of CPU time) |
| Diagnostics (instrumented runs; `mpce_numbers_diag.tex` (`\ND...`), `mpce_supp_diag.tex`, `mpce_summary_diag.json`, `../figures_mpce/diag_*.pdf`) | `python3 mpce_diagnostics.py [--procs=2] [--cache=DIR]` (~18 min on 2 cores; with `--cache` the raw runs are stored in `DIR/diag_raw.pkl` and reused, so a re-run only post-processes (~1.5 min); `--rerun=T1,T2,T3,T1b` recomputes chosen studies, and studies missing from the cache are run). Study T1b (review round 2, R3-7): the 23 stored old-setting PSO runs whose final boundary status cannot be decided from the 1-mm coordinates (and that are not among the instrumented T1 runs) are re-run uninstrumented (~40 s) and verified against the stored runs; the final layouts are then classified with the paper's 1e-6 m rule (infeasible = stored `Feasible` flag, 286 of 2,040). Then `python3 mpce_check_diag.py` (also run by the build on the stored summary; D22 checks that `mpce_numbers_diag.tex` equals the macros regenerated from `mpce_summary_diag.json`) |
| Theory figures and checks T01... (`../figures_mpce/theory_*.pdf`) | `python3 make_theory_figures.py --check` (~3 min) |

# Reproducibility package — SWEVO paper on baseline configuration and random-sampling controls in WFLO

**Version 1.0.0-rc (release candidate, 2026-10-04).** Package for the article
*Baseline Configuration and Random-Sampling Controls in Metaheuristic Comparisons for Wind Farm Layout Optimization*
(P. Solanki, P. Dwivedi, V. Garg, V. Shukla; submitted to *Swarm and Evolutionary Computation*).

**Status of the public deposit: pending author action.** No DOI exists yet. The authors create the versioned
public deposit (e.g. a GitHub release archived by Zenodo), then enter the DOI in `CITATION.cff` (commented
placeholder) and in the Data availability statement of the article. The licence is not yet chosen; see
`LICENSE_SUGGESTION.md` (suggestion only: MIT for code, CC-BY-4.0 for data — **to be confirmed by the authors**).

## 1. Purpose and contents

The package contains everything needed to regenerate every table, figure and number of the article and its
supplement from the stored per-run records, and to rerun every optimization study:

| Part | Where | Count |
|---|---|---|
| Per-run records of the original study (one row per run: method, case, seed, budget, initialization, evaluations, objective, feasibility, minimum spacing, seconds, final coordinates, convergence curve) | `analysis/fresh_{grid,vgrid,bgrid}.csv`, `analysis/mpce_*_s<i>of<k>.csv` | **36,190** runs (`count_records.py`) + the superseded Horns Rev rows that the analysis drops + the 8 `feasx` shards merged into `mpce_feas_s0of1.csv` |
| Per-run records of revision 2 (GA, exact gradients, Laplace ablation, 5D/6D spacing, Lillgrund) | `analysis/rev2_*_s<i>of<k>.csv` | **11,640** runs |
| Records and outputs of revision 3 | `analysis/rev3_*` | the files present when `make_package.sh` ran (listed in `inventory.csv`; see Section 7) |
| Models and objectives | `analysis/wflop_model.py`, `authors_objective.py`, `objective_original.py`, `hornsrev_model.py`, `iea37_model.py` (+ `iea37_data/`), `rev2_site_model.py` (Lillgrund), `record_io.py` | |
| Optimizers | `analysis/authors_optimizers.py`, `extra_baselines.py`, `original_vns.py`, `hybrid_lxssa_bvns.py`, `rs_vns.py`, `feasible_init.py`, `init_hook.py` (+ `hybrid_lxssa_vns.py`, imported by `full_grid_experiments.py`) | |
| Experiment drivers | `analysis/full_grid_experiments.py`, `mpce_experiments.py`, `iea37_experiments.py`, `rev2_{ga,gradient,laplace,site,spacing}.py`, `authors_experiments.py`, `ssa_reference.py`, `mpce_diagnostics.py`, `packing_capacity.py`, `pywake_check.py`, `iea37_projected.py`; launch scripts `SWEVO_rev2/experiments_docs/run_rev2_*.sh` | |
| Analysis and checks | `analysis/mpce_results.py`, `mpce_robustness.py`, `mpce_numbers.py`, `mpce_inference_extra.py`, `mpce_direction.py`, `mpce_csweep.py`, `make_theory_figures.py`, `make_model_figures.py`, `analyze_authors_runs.py`, `rev2_analysis.py`, `audit_archive.py`, `check_revision.py`, `validate_evaluator.py`, `calibrate_authors_code.py`, `mpce_check_{final,extra,diag,dir,csweep}.py` | |
| Generated outputs (as used by the article) | `analysis/mpce_tab_*.tex`, `mpce_supp*.tex`, `mpce_supplementary.tex`, `mpce_numbers*.tex`, `mpce_summary*.json`, statistics CSVs, `figures_mpce/`, `figures_final/fig_*.pdf`, `SWEVO_rev2/experiments_docs/rev2_{summary.json,tables.tex}`, `SWEVO_rev2/validation/` | |
| Package tools | `regenerate_all.sh`, `compare_outputs.py`, `verify_determinism.py`, `count_records.py`, `build_inventory.py`, `make_package.sh`, `inventory.csv`, `MANIFEST.sha256` | |
| Verification logs of this release | `verification/` | |

Files of the repository that belong to earlier papers or superseded versions (manuscript builders of earlier
papers, `superseded_linear_penalty/`, `final_*`, `hybrid_*`, `bvns_*`, `fresh_results.py`, `extra_*`, the
earlier `figures_final/` content, ...) are excluded; `inventory.csv` lists every candidate file with its decision
(include / exclude / pending) and the reason.

## 2. Layout

The package keeps the repository paths, so that every script runs unchanged:

```
SWEVO_repro_<version>/
  README.md  CITATION.cff  LICENSE_SUGGESTION.md  requirements.txt  requirements-pywake.txt  pip-freeze-full.txt
  Dockerfile  MANIFEST.sha256  inventory.csv  regenerate_all.sh  compare_outputs.py  verify_determinism.py
  count_records.py  build_inventory.py  make_package.sh
  analysis/                 code, run records, generated tables / numbers / summaries, iea37_data/
  figures_mpce/             generated figures (PDF + PNG previews)
  figures_final/            model schematics fig_wind_farm / fig_wake_model / fig_half_cone (PDF)
  selected_30_run_data.csv  720 archived runs of the original code (evaluator validation, calibration)
  SWEVO_rev2/experiments_docs/   revision-2 summary/tables and launch scripts
  SWEVO_rev2/validation/         revision-2 archive audit outputs
  verification/             logs of the verification of this release
```

All scripts are run from `analysis/` (they find their inputs relative to their own location).

## 3. Environment

Verified with CPython **3.11.15** on Linux x86_64 (Ubuntu 24.04), NumPy **2.4.6**, SciPy **1.17.1**, pandas
**3.0.6**, Matplotlib **3.11.2**, PyYAML **6.0.1** (pip wheels; NumPy/SciPy with their bundled OpenBLAS 0.3.31) —
the environment of the runs. Optional: PyWake **2.6.20** (+ xarray, netCDF4, autograd, ...) for the Horns Rev 1 and
Lillgrund cross-checks only.

```sh
python3.11 -m venv .venv && . .venv/bin/activate
python -m pip install -r requirements.txt               # everything except the PyWake cross-checks
python -m pip install -r requirements-pywake.txt        # optional
```
or with Docker: `docker build -t swevo-repro . && docker run --rm -v "$PWD":/pkg -w /pkg swevo-repro`.
`pip-freeze-full.txt` is the complete `pip freeze` of the verification container (reference only).
Bitwise reproduction of MS-SLSQP runs (SciPy SLSQP) needs the same NumPy/SciPy build and BLAS; all other
optimizers reproduce bit for bit for a given seed (Section 6).

## 4. Verify the checksums

```sh
cd SWEVO_repro_<version>
sha256sum -c MANIFEST.sha256          # every file of the package except MANIFEST.sha256 itself
```
The sha256 of the zip file is in `SWEVO_repro_<version>.zip.sha256` (next to the zip) and, after the deposit, on
the deposit page. `SWEVO_rev2/validation/raw_data_checksums.json` holds the checksums of the 140 revision-2 run
files as delivered in revision 2; `python3 analysis/check_revision.py` (run by `regenerate_all.sh`) verifies them.

## 5. Regenerate the tables, figures and numbers

One command regenerates everything that is computed from the stored run records, in a scratch copy, and compares
each regenerated file with the stored copy:

```sh
sh regenerate_all.sh /tmp/regen                     # ~6 min on one core
WITH_SLOW=1 WITH_DETERMINISM=1 sh regenerate_all.sh /tmp/regen     # + theory figure, evaluator, packing, schematics, reruns
WITH_PYWAKE=1 sh regenerate_all.sh /tmp/regen       # + PyWake cross-checks (requirements-pywake.txt)
MANUSCRIPT=/path/to/latex_source WITH_SLOW=1 sh regenerate_all.sh /tmp/regen   # + CHECK-tag scan and theory text checks
```
Results: `/tmp/regen/logs/steps.txt` (exit code and seconds of every step), `/tmp/regen/logs/compare.txt`
(IDENTICAL / IDENTICAL-NORMALIZED / DIFFERENT per file). The diagnostics outputs (`mpce_*diag*`, `diag_pso_dynamics.*`) are not regenerated (instrumented reruns, `mpce_diagnostics.py`); their stored copies are checked by D01–D22. IDENTICAL-NORMALIZED means identical after removing
only time stamps of the generation, PDF creation dates, gzip header times, absolute paths, run-time statements
and environment version strings (see `compare_outputs.py`). The package files themselves are never written.

### 5.1 Map: manuscript item -> script -> command -> input files -> runtime

Commands are run in `analysis/`; runtimes are wall-clock seconds on one core of the verification machine
(shared, loaded); "stored" means the article uses a stored output whose generation needs a long rerun.

| Manuscript item | Script | Command | Input files | Runtime |
|---|---|---|---|---|
| Table `tab:baseline` (PSO setting, eight-method pool) | `mpce_results.py` | `python3 mpce_results.py --procs 1` | `fresh_grid.csv`, `fresh_vgrid.csv`, `fresh_bgrid.csv`, `mpce_psoc_s0of1.csv`, `mpce_slsqp_s0of1.csv` | 56–89 s (all `mpce_results.py` outputs together) |
| Tables `tab:friedman68`, `tab:wtl`; Figs. `fig:wakeloss` (`wakeloss_vs_n.pdf`), `fig:ablation` (`ablation_convergence.pdf`); supplementary per-case tables, cost (`tab:cost`), provenance, `mpce_supplementary.tex`; `avg_ranks`, `feasibility_vs_n`, `convergence_*`, `boxplots_*`, `layouts_max` figures | `mpce_results.py` | same run | the 13 benchmark run files: `fresh_{grid,vgrid,bgrid}.csv`, `mpce_{psoc,psobv,rsvns,rsdisc,slsqp}_s*of*.csv` | (same run) |
| Table `tab:ablation` (component analysis, 15 contrasts) | `mpce_results.py` | same run | + `mpce_rsvns_s0of1.csv`, `mpce_rsdisc_s0of1.csv` | (same run) |
| Split table (`mpce_tab_split.tex`, supplement) | `mpce_results.py` | same run | `mpce_psosplit_s0of1.csv`, `mpce_omega90_s0of1.csv` | (same run) |
| Table `tab:feasbudget`, Fig. `fig:budget` (`budget_scaling.pdf`) | `mpce_results.py` | same run | `mpce_feas_s0of1.csv`, `mpce_feasp_s0of1.csv`, `mpce_b30k*_s*.csv`, `mpce_b120k*_s*.csv`, `mpce_hrfix_s*of24.csv` | (same run) |
| Table `tab:hr-site` (Horns Rev 1), `hr16.pdf`, `layouts_iea37_hr16.pdf` | `mpce_results.py` (`mpce_tab_hr16.tex`); GA row: `rev2_analysis.py` (`rev2_summary.json` -> `ga_hr16`) | as above; `python3 rev2_analysis.py --out-dir OUT` | `mpce_hrfix_s*of24.csv`, `pywake_check.csv`; `rev2_gahr_s*of2.csv` | 56 s; 14 s |
| Table `tab:iea37` (best runs vs published layouts), `iea37_*.pdf` | `mpce_results.py` (`mpce_tab_iea37.tex`); GA and analytic-gradient rows: `rev2_analysis.py` (`ga_iea37`, `grad_iea37`) | as above | `mpce_iea16*/iea36*_s0of1.csv`, `iea37_published_results.csv`; `rev2_gaiea_*.csv`, `rev2_grad_*.csv` | 56 s; 14 s |
| Tables `tab:iea-means`, `tab:ga-main` (typed from summaries) | `rev2_analysis.py`, `mpce_results.py` | as above | `rev2_ga_*.csv`, `rev2_gaiea_*.csv`, `rev2_grad_*.csv` + main-comparison files | 14 s |
| Table `tab:robust-final` (cubic curve, Gaussian wake) | `mpce_results.py` -> `mpce_robustness.py` | same run (cache `mpce_reevaluation_cache.csv`; full recomputation: delete the cache; not timed here) | final layouts of the run files, `final_reevaluation.csv` (cache seed) | (same run) |
| Table `tab:equiv-main`, Fig. `fig:equiv-curve` (`equiv_curve.pdf`), supplementary Section S-inference (`mpce_supp_inference.tex`), macros `\NX...` (`mpce_numbers_extra.tex`) | `mpce_inference_extra.py` | `python3 mpce_inference_extra.py` | benchmark run files, `mpce_summary.json`, `mpce_reevaluation.csv` | 31 s |
| Fig. `fig:S-csweep` (`csweep.pdf`), `mpce_supp_csweep.tex`, `\NS...` (`mpce_numbers_csweep.tex`) | `mpce_csweep.py` | `python3 mpce_csweep.py` | `mpce_csweep_s0of1.csv`, `mpce_psoc_s0of1.csv`, `fresh_grid.csv` | 4 s |
| Direction resolution (1° bins), Horns Rev 1 PyWake comparison, projected IEA37 layouts: `mpce_supp_direction.tex`, `\NF...` (`mpce_numbers_dir.tex`) | `mpce_direction.py` (+ `iea37_projected.py`) | `python3 mpce_direction.py --procs 1` | all final layouts, `pywake_check.csv`, `pywake_check_hr16runs.csv`, `iea37_published_results.csv` | 136–185 s |
| Fig. 1 (`fig:model`: `figures_final/fig_wind_farm.pdf`, `fig_wake_model.pdf`, `fig_half_cone.pdf`) | `make_model_figures.py` | `python3 make_model_figures.py` | `final_best_layouts_maxN.csv` (example layout) | 1 s |
| Fig. `fig:S-th-stability` (`theory_stability.pdf`), Section S-theory numbers (checks T01–T13) | `make_theory_figures.py` | `python3 make_theory_figures.py` (figure); `--check --swevo` with the manuscript sources (Section 5) | none (Monte Carlo, seed 20260928) | 210–290 s (with or without `--check`) |
| Number macros `\N...` (`mpce_numbers.tex`) | `mpce_numbers.py` (called by `mpce_results.py`) | `python3 mpce_numbers.py` | `mpce_summary.json` | < 1 s |
| Outcome-dependent statements (CHECK-FINAL C01–C62) | `mpce_check_final.py` (called by `mpce_results.py`) | `python3 mpce_check_final.py [--tex .../SWEVO_manuscript.tex]` | `mpce_summary.json` (+ manuscript for the tag scan) | < 1 s |
| Checks X01–X57, D01–D22, F01–F46, S01–S15 | `mpce_check_extra.py`, `mpce_check_diag.py`, `mpce_check_dir.py`, `mpce_check_csweep.py` | `python3 mpce_check_<x>.py` | `mpce_summary_{extra,diag,dir,csweep}.json` | 1 s each |
| Algorithm diagnostics (Section S-diag, `mpce_supp_diag.tex`, `\ND...`, `diag_pso_dynamics.pdf`) | `mpce_diagnostics.py` | `python3 mpce_diagnostics.py --procs=2` (instrumented reruns; stored summary `mpce_summary_diag.json` checked by `mpce_check_diag.py`) | `fresh_grid.csv`, `fresh_bgrid.csv`, `mpce_psoc/psobv/rsvns/feas/psosplit/hrfix` | stored; ~18 min on 2 cores (not rerun here) |
| Revision-2 tables (Laplace ablation, GA, 5D/6D spacing, exact gradients, Lillgrund; `rev2_tables.tex`, `rev2_summary.json`) | `rev2_analysis.py` | `python3 rev2_analysis.py --out-dir OUT` | `rev2_*_s*of*.csv` + main-comparison and IEA37/Horns Rev files | 14 s |
| Lillgrund reliability tables (supplement `archive_reliability_tables.tex` = `SWEVO_rev2/validation/reliability_tables.tex`), archive audit | `audit_archive.py`, `check_revision.py` | run in a folder with `experiments/` (the `rev2_*` CSVs + the two scripts + `record_io.py`, `rev2_site_model.py`) and `validation/`; done by `regenerate_all.sh` | `rev2_*_s*of*.csv` | 3 s |
| Spacing re-optimization of SSA / LX-SSA (`tab:spacing-authors`, `authors_tables.tex`) | `analyze_authors_runs.py` | `python3 analyze_authors_runs.py` | `authors_runs_{budget,crowded,spacing,hornsrev}.csv`, `../selected_30_run_data.csv` | 1 s |
| Packing capacity (`tab:capacity`) | `packing_capacity.py` | `python3 packing_capacity.py` | none | 39 s |
| Evaluator verification (720 archived objectives, max error 8.7e-11) | `validate_evaluator.py` | `python3 validate_evaluator.py` | `../selected_30_run_data.csv` | 1 s |
| IEA37 calculator check (relative AEP difference < 1e-11) | `iea37_model.py` | `python3 iea37_model.py` | `iea37_data/` | < 1 s |
| Calibration of the reused code (2 of 24 Mann–Whitney tests reject) | `calibrate_authors_code.py` | `python3 calibrate_authors_code.py SSA,LXSSA authors` | `../selected_30_run_data.csv` | reruns, minutes (not run here) |
| Horns Rev 1 / Lillgrund vs PyWake | `pywake_check.py`, `rev2_site_model.py` | `python3 pywake_check.py --layouts --out-dir OUT`; `python3 rev2_site_model.py --pywake` | `mpce_hrfix_s*of24.csv` | 46 s; 4 s |
| Typed tables `tab:rq-matrix`, `tab:params`, `tab:design`, `tab:recommendations`; run counts | design constants; `count_records.py` | `python3 ../count_records.py --analysis-dir .` | all run files | 5 s |
| Supplementary `tab:S-eqclus`, `tab:S-allrun-spacing`, `tab:S-time` (added in revision 2) | typed by the revision-2 authors from stored tables and records; no generating script in the revision-2 archive. Revision-3 scripts (`rev3_inference.py`, `rev3_numbers_audit.py`, if present) recompute them | — | stored tables, `rev2_spacing_*.csv` | — |

Note on the LaTeX sources: the article's sections paste the generated tables between `% >>>>> begin
analysis/<file>` / `% <<<<< end` markers and the revision-2 authors edited the pasted copies (sentence-case
captions, typography, the added GA / analytic-gradient rows of `tab:hr-site` and `tab:iea37`, whose values are in
`rev2_summary.json`). The generated files of this package are the unedited pipeline outputs.

## 6. Rerun the studies (optimization runs)

Every run file is produced by one command (run in `analysis/`; `--procs P` sets the worker processes; shards
`i k` split a study into k parts). Expected CPU time = sum of the stored per-run `Seconds` (single-core time of the
original machines; divide by the number of workers).

| Run files | Command | Runs | CPU time |
|---|---|---|---|
| `fresh_grid.csv` (SSA, LX-SSA, DE, old PSO; + unused SLSQP/VNS rows) | `python3 full_grid_experiments.py grid` | 12,240 | 8.8 h |
| `fresh_vgrid.csv` (VNS) | `python3 full_grid_experiments.py vgrid` | 2,040 | 1.3 h |
| `fresh_bgrid.csv` (SSA-VNS, LX-SSA-VNS) | `python3 full_grid_experiments.py bgrid` | 4,080 | 2.7 h |
| `fresh_{hr16,vhr16,bhr16}.csv` (old HR binning, superseded by `hrfix`) | `full_grid_experiments.py hr16 / vhr16 / bhr16` | 270 | 0.4 h |
| `mpce_psoc_s0of1.csv` (PSO) | `python3 mpce_experiments.py psoc` | 2,040 | 0.8 h |
| `mpce_psobv_s<i>of2.csv` (PSO-VNS) | `python3 mpce_experiments.py psobv <i> 2` | 2,070 | 1.2 h |
| `mpce_rsvns_s0of1.csv` (RS-VNS) | `python3 mpce_experiments.py rsvns` | 2,040 | 1.0 h |
| `mpce_rsdisc_s0of1.csv` (RSD-VNS) | `python3 mpce_experiments.py rsdisc 0 1` | 2,040 | 1.3 h |
| `mpce_slsqp_s0of1.csv` (MS-SLSQP) | `python3 mpce_experiments.py slsqp` | 2,040 | 6.2 h |
| `mpce_psosplit_s0of1.csv`, `mpce_omega90_s0of1.csv` (split) | `mpce_experiments.py psosplit`; `omega90 0 1` | 1,080 | 0.8 h |
| `mpce_csweep_s0of1.csv` (PSO coefficient sweep) | `python3 mpce_experiments.py csweep 0 1 --procs=4` | 3,240 | 2.6 h |
| `mpce_feasx_s<i>of8.csv` -> `mpce_feas_s0of1.csv`; `mpce_feasp_s0of1.csv` (feasible init.) | `mpce_experiments.py feasx <i> 8` (i = 0..7), then the concatenation command of `analysis/README_reproduce.md`; `feasp` | 1,890 | 5.1 h |
| `mpce_b30k_s<i>of3.csv`, `mpce_b30kp_s0of1.csv` | `mpce_experiments.py b30k <i> 3`; `b30kp` | 2,100 | 17.0 h |
| `mpce_b120k_s<i>of8.csv`, `mpce_b120kp_s<i>of4.csv` | `mpce_experiments.py b120k <i> 8`; `b120kp <i> 4` | 1,900 | 51.6 h |
| `mpce_hrfix_s<i>of24.csv` (all Horns Rev 1 runs) | `python3 mpce_experiments.py hrfix <i> 24` (i = 0..23) | 970 | 9.6 h |
| `mpce_hr16new_s0of1.csv` (superseded) | `python3 mpce_experiments.py hr16new` | 60 | 0.2 h |
| `mpce_iea16/iea36(+p)_s0of1.csv` (IEA37) | `python3 iea37_experiments.py iea16 / iea36 / iea16p / iea36p` | 1,200 | 3.2 h |
| `rev2_laplace_s<i>of12.csv` | `python3 rev2_laplace.py laplace <i> 12 --procs 2` | 2,160 | 1.3 h |
| `rev2_ga_s<i>of10.csv`, `rev2_gahr_s<i>of2.csv`, `rev2_gaiea_s<i>of2.csv` | `python3 rev2_ga.py ga <i> 10` / `gahr <i> 2` / `gaiea <i> 2` | 2,220 | 1.6 h |
| `rev2_spacing_s<i>of40.csv` | `python3 rev2_spacing.py spacing <i> 40 --procs 2` | 6,480 | 4.7 h |
| `rev2_lg16_s<i>of4.csv`, `rev2_lg16b_s<i>of30.csv` (Lillgrund) | `python3 rev2_site.py lg16 <i> 4 --procs=2`; `lg16b <i> 30` | 540 | 2.9 h |
| `rev2_grad_s<i>of40.csv` (exact gradients) | `python3 rev2_gradient.py grad <i> 40 --procs=2` | 240 | 2.0 h |
| `authors_runs_*.csv` | `python3 authors_experiments.py budget / crowded / spacing / hornsrev` | — | ~30 min on 4 cores |
| revision-3 studies | see the prespecification files `analysis/rev3_*_manifest.md` and the docstrings of `analysis/rev3_*.py` | — | — |

Total stored CPU time of the original and revision-2 runs: about 126 CPU-hours. The revision-2 launchers with the
exact shard counts are `SWEVO_rev2/experiments_docs/run_rev2_resumable.sh` (resumable) and `run_rev2_hpc.sh`.
Determinism: a rerun of a stored run reproduces all fields except `Seconds`; `verify_determinism.py` checks seed 1
of two benchmark cases for PSO-VNS and SSA-VNS against the stored records (Section 8).

## 7. Known limits

- **Legacy coordinate rounding.** The original and revision-2 run files store final coordinates with 3 decimals
  (benchmark, IEA37: 1 mm), 2 decimals (Horns Rev 1, Lillgrund: 1 cm) and convergence curves with 3–4 decimals;
  objective, wake loss, feasibility label and minimum spacing are stored at full double precision from the
  unrounded layout. The 10^-6 m feasibility tolerance of the paper therefore cannot be re-verified from the
  stored coordinates alone (e.g. 236 of 240 labelled-feasible exact-gradient runs fail a strict recheck of the
  rounded coordinates; `SWEVO_rev2/validation/archive_audit.json`). The drivers are deterministic for a given
  seed (MS-SLSQP only on the same NumPy/SciPy/BLAS build; reruns in Section 8 are bit-identical to the records),
  so the full-precision layouts can be recovered by rerunning the stored runs; revision 3 does this and writes them with 17 significant digits
  (`record_io.encode_coordinates`) as side files `analysis/rev3_fullprec_*.csv` / `analysis/rev3_precision_*`
  (included when present; see `inventory.csv` for the status at packaging time). New runs use the 17-digit writer.
- **MS-SLSQP** results depend on the NumPy/SciPy build and BLAS; `mpce_slsqp_s0of1.csv` was rerun on the
  platform of the other runs; the old-platform SLSQP rows of `fresh_grid.csv` are not used.
- **Diagnostics** (`mpce_diagnostics.py`, instrumented reruns) and the full robustness re-evaluation (delete
  `mpce_reevaluation_cache.csv`) are not rerun by `regenerate_all.sh`; their stored summaries are checked
  (D01–D22) and their cached values reused.
- **Elapsed times** (`Seconds`) come from the machines and batches that produced each file and are indicative only.
- The article's LaTeX sources are not part of this data package by default (`MANUSCRIPT_DIR` in `make_package.sh`
  adds them); the CHECK-tag scan and the theory text checks need them.
- Third-party inputs (IEA37 files, PyWake-derived site data) keep their own terms (`LICENSE_SUGGESTION.md`).

## 8. Verification of this release (2026-10-04)

See `verification/verification_summary.txt` and the logs next to it. In short: all regenerated tables, number
macros, summaries, statistics files and figures are byte-identical to the stored copies or differ only in
generation time stamps / PDF creation dates / gzip header times / recorded environment versions; all check
scripts pass (C01–C62, X01–X57, D01–D22, F01–F46, S01–S15, T01–T13); the revision-2 summary (`rev2_summary.json`)
and tables are reproduced; reruns of stored PSO-VNS and SSA-VNS runs are bit-identical to the records.

## 9. How to cite

See `CITATION.cff` (DOI to be added after the deposit).

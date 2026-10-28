# GA in the budget / feasible-initialization, spacing and Horns Rev 1 studies (rev5_ga_secondary) — reproduction

The real-coded GA (`rev2_ga.RCGA`, label `GA`) was missing from the budget-scaling / feasibility-preserving
initialization study, from the 5D / 6D spacing re-optimization and from the Horns Rev 1 settings 120,030 (random) and
6,030 (feasible). These runs add it with the design, runner, seeds, case set, budgets and record columns of the other
methods of each study. No existing record is modified.

## 1. Runs

Run from this `analysis/` directory (Python 3.11.15, NumPy 2.4.6, SciPy 1.17.1, pandas 3.0.6; 4 cores).

```sh
python3 rev5_ga_secondary.py verify                      # seed pairing + runner checks  -> rev5_ga_secondary_verify.log
python3 rev5_ga_secondary.py b120k 0 1 --procs=4         # -> rev5_ga_b120k_s0of1.csv  (180 runs)
python3 rev5_ga_secondary.py b30k  0 1 --procs=4         # -> rev5_ga_b30k_s0of1.csv   (180 runs)
python3 rev5_ga_secondary.py hr    0 1 --procs=4         # -> rev5_ga_hr_s0of1.csv     ( 40 runs)
python3 rev5_ga_secondary.py feas  0 1 --procs=4         # -> rev5_ga_feas_s0of1.csv   (180 runs)
python3 rev5_ga_secondary.py spacing 0 1 --procs=4       # -> rev5_ga_spacing_s0of1.csv (720 runs)
```

| Mode | Design (as the other methods) | Runner | Runs |
|---|---|---|---|
| `b30k` | six largest benchmark cases (`mpce_experiments.LARGE`), seeds 1–30, 30,030 calls, random init. (as `mpce_experiments.py b30k` / `b30kp`) | `mpce_experiments.run_grid` via `rev2_ga._with_ga` | 180 |
| `b120k` | six largest cases, seeds 1–30, 120,030 calls, random (as `b120k` / `b120kp`) | same | 180 |
| `feas` | six largest cases, seeds 1–30, 6,030 calls, feasibility-preserving init. (as `feas` / `feasx` / `feasp`; `init_hook.GEN = feasible_init.make_generator(4D, circle)`) | same | 180 |
| `hr` | Horns Rev 1 16-turbine block: 120,030 random, seeds 1–10, and 6,030 feasible, seeds 1–30 (as `hrfix`) | `mpce_experiments.run_hr` via `rev2_ga._with_ga` (as `rev2_ga.py gahr`) | 40 |
| `spacing` | 12 split cases (`mpce_experiments.SPLITCASES`), ell_min = 5D and 6D, seeds 1–30, 6,030 calls, random (as `rev2_spacing.py spacing`) | `rev2_spacing.run_grid_sp` with the dispatcher `mpce_experiments.run_method_x` swapped for `_run_method_x_ga` (GA -> `rev2_ga.run_method_ga`, every other label unchanged) | 720 |

GA at 6,030 random on the benchmark is the stored `rev2_ga_s*of10.csv` (68 cases, seeds 1–30) and on Horns Rev 1 at
6,030 / 30,030 random the stored `rev2_gahr_s*of2.csv`; they are reused (6,030 column of the budget study, 4D column of
the spacing study), as the stored runs of the other methods are.

GA settings (unchanged `rev2_ga.RCGA`): N_p = 30, binary tournament, SBX eta_c = 15, p_c = 0.9 (per-variable
probability 0.5), polynomial mutation eta_m = 20, p_m = 1/n, clipping to the box, generational replacement with
elitism; T = (B - 30) // 30 generations, i.e. exactly B calls. Tracker, penalty, feasibility test and convergence
checkpoints ((B - 30) // 200 calls) are those of the runner. Seeding: `np.random.seed(seed)` in the constructor, then
the common `init_pop` draw (random or feasible generator), so GA's initial population is that of the other methods.

Record columns: those of the other methods of the study (`mpce_*`: Algorithm ... Seconds, Coordinates, Curve; spacing:
plus Spacing, SMin). Coordinates are written with 17 significant digits (`record_io.encode_coordinates`, as the current
`rev2_ga.py` / `rev2_spacing.py`); the archived `mpce_b30k/b120k/feas` records of the other methods store 3 decimals.

## 2. Verification (`rev5_ga_secondary_verify.log`)

* Seed pairing: for seeds 1 and 2 the first 30 evaluated layouts of GA are identical (`np.array_equal`) to those of
  PSO and PSO-VNS in each design: 30,030 random (DS I, 500 m, N = 10), 120,030 random (DS II, 1000 m, N = 15), 6,030
  feasible (DS I, 750 m, N = 12), spacing 5D (DS II, 750 m, N = 8), spacing 6D (DS I, 1000 m, N = 10), Horns Rev 1
  120,030 random and 6,030 feasible: 28 / 28 comparisons identical.
* Runner check: the same code path reproduces stored records of another method bit for bit (Objective, Calls,
  Feasible): PSO, DS I, 500 m, N = 10, seed 1, 6,030 feasible (`mpce_feas_s0of1.csv`); PSO, DS II, 750 m, N = 8,
  seed 1, 5D (`rev2_spacing`); PSO, Horns Rev 1, seed 1, 6,030 feasible (`mpce_hrfix`).
* At 6,030 feasible, the first stored convergence checkpoint (best feasible of the 30 initial layouts) of every GA run
  equals that of every other method of the same case and seed (`rev5_ga_secondary_summary.json`,
  `feasible_init_pairing`).

## 3. Analysis (GA added to the pools)

```sh
python3 mpce_results.py                 # GA in sections 1 (Type I supplement), 5, 6 and tab:feasbudget; --no-ga-secondary: as before
python3 rev2_analysis.py --out-dir ../SWEVO_rev2/experiments_docs   # GA in the spacing pool and the Horns Rev GA cells; --no-ga-spacing: as before
python3 rev5_ga_secondary_summary.py    # -> rev5_ga_secondary_summary.json (GA results, all-run spacing rows, checks)
```

With `--no-ga-secondary` / `--no-ga-spacing` both generators reproduce the stored outputs (all `mpce_*` tables,
supplement, numbers, CSVs and figures byte-identical or identical after time-stamp normalization; `mpce_summary.json`
differs only by the GA entries of its `labels` / `colours` dictionaries; `rev2_summary.json` / `rev2_tables.tex`
identical), checked with `../repro_package_build/compare_outputs.py`. Changes of the outputs: see
`rev5_ga_secondary_report.md`.

## 4. Timings

See `rev5_ga_secondary_run.log` (wall time per mode on 4 worker processes; the machine was shared with another
4-process-equivalent job during most of the runs, so the per-run `Seconds` are inflated and are not comparable with the
timings of other batches).

| Mode | Runs | Wall time (4 processes) | Sum of per-run Seconds | Median s per run |
|---|---|---|---|---|
| b120k | 180 | 2,282 s (38 min) | 9,083 s | 49.5 |
| b30k | 180 | 504 s | 2,005 s | 10.8 |
| hr | 40 | 307 s | 1,205 s | 97.6 (120,030); 6.7 (6,030 feasible) |
| feas | 180 | 283 s | 1,128 s | 5.8 (incl. the packing initialization) |
| spacing | 720 | 382 s | 1,527 s | 2.06 (5D), 2.07 (6D) |
| total | 1,300 | 3,758 s (63 min) | 14,948 s (4.2 CPU-h) | |

Runs made 2026-10-28 03:43-04:46 UTC; `verify` 2026-10-28 (about 1 min). Analysis: `mpce_results.py` about 2 min,
`rev2_analysis.py` about 1 min, `rev5_ga_secondary_summary.py` about 1 min.

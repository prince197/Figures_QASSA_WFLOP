# rev3_fine -- direct optimization with 1-degree direction bins (revision 3, item C2 / reviewer concern R8)

Status: PRESPECIFIED (written 2026-10-04, before any run of this experiment). Driver `rev3_fine.py`, analysis
`rev3_fine_analysis.py`. Experiment id `rev3_fine` (primary arm) and `rev3_fine15` (optional control arm); the
legacy experiments (`mpce_*`, `rev2_*`) are not touched and their file names are not reused.

## Question
Section "Direction Resolution" (09_robust.tex) only RE-EVALUATES stored layouts, optimized with the 24 x 15-deg
benchmark rose, with 5- and 1-deg sub-bins, and states: "Which method would lead if the search itself used 1-deg
bins was not tested." Question: which method leads when the SEARCH itself uses the 1-deg rose, and does the leader
(and the PSO-VNS vs. PSO verdict) differ from the 1-deg re-evaluation of the same methods' 15-deg layouts on the same
cases?

## Objective (the only change to the benchmark protocol)
The benchmark objective with S = 15 sub-bins per 15-deg bin, i.e. 360 bins of 1 deg centred at 0.5, 1.5, ... deg
(mathematical angle as wflop_model.THETA): sub-bin m of bin j has centre theta_j - 7.5 + (m + 0.5) deg, the Weibull
parameters (psi_j, k_j) of bin j and frequency w_j / 15. This is exactly `mpce_direction.bench_objective(xy, ds, 15)`
(the evaluator of the paper's 1-deg re-evaluation); `rev3_fine.FineObjective` is a faster implementation of the same
function (identical floating-point operations on the in-wake pairs; validated bit for bit, see below). The wake-free
value (Ideal) is the same for every S. Penalty (authors' exact penalty, 1e10 scale, `authors_objective.make_objective`
with farm_objective replaced by the 1-deg objective), box repair (clipping to [-r, r]), feasibility tolerance (1e-6 m;
spacing 8R = 308 m = 4D, inside the circle of radius r), Tracker (every call counts, feasible-best curve with 201
checkpoints) and output columns are those of `mpce_experiments.run_grid`. Wake model: Jensen (Kusiak-Song) as in the
paper.

## Cases (fixed; 8 = both data sets x low and high density)
Per data set (I and II): (r = 500 m, N = 10), (r = 750 m, N = 6), (r = 750 m, N = 12), (r = 1000 m, N = 15).
N = 10 at 500 m, 12 at 750 m and 15 at 1000 m are the largest N of the benchmark (high density); 750 m, N = 6 is a
low-density case.

## Methods (as in the paper, Table 2; N_p = 30; labels of the stored data)
PSO-VNS (PSOBV; constriction PSO phase 1, omega = 0.5, then basic VNS), PSO (PSOC; w = 0.7298, c1 = c2 = 1.49618),
GA (rev2_ga.RCGA; SBX eta_c = 15, p_c = 0.9, polynomial mutation eta_m = 20, p_m = 1/n, tournament, elitism),
MS-SLSQP (SLSQP; SciPy SLSQP, maxiter 200 per start, ftol 1e-9, forward differences with step 1 m; minimizes the
1-deg wake loss under the scaled constraints), RSD-VNS (RSDVNS; initial population, then disc samples, VNS).
Code paths: `mpce_experiments.run_method` / `run_method_x` and `rev2_ga.run_method_ga`, unchanged.

## Seeds, budget, initialization
NEW seeds 31-60 (30 seed-paired runs per method and case: equal seeds share the initial population). Budget
B = 6,030 objective evaluations per run; each evaluation of the 1-deg objective counts once (initial population,
MS-SLSQP finite-difference calls and repeated evaluations included; constraint evaluations not charged). Random
initialization (uniform in [-r, r]^{2N}, `init_hook.GEN = None`).
Primary arm: 8 cases x 5 methods x 30 seeds = 1,200 runs (`rev3_fine_s<i>of<k>.csv`, column Rose = "1deg").
Optional control arm (only if cloud time is left; analysed if present, never required): the same 1,200 runs with the
15-deg objective (`rev3_fine15_s<i>of<k>.csv`, Rose = "15deg"), giving seed-paired direct-vs-re-evaluated
differences.

## Outputs per run
Standard columns (Algorithm, Dataset, Radius, Turbines, Seed, Budget, Init, Calls, Objective, Ideal, WakeLoss,
Feasible, MinSpacing, Seconds, Coordinates [17 significant digits, record_io.encode_coordinates], Curve), then
Rose, SubBins, Obj15, Ideal15, ObjFine, IdealFine (the final layout evaluated with BOTH roses: 15-deg benchmark
objective and 1-deg objective), Experiment. Objective / Ideal / WakeLoss are those of the run's own rose.

## Reference arm (no new runs)
The stored main-comparison runs of the same 5 methods, cases and seeds 1-30 (mpce_psobv, mpce_psoc, mpce_slsqp,
mpce_rsdisc, rev2_ga), re-evaluated with `mpce_direction.bench_objective(xy, ds, 15)` (feasible rows; infeasible rows
keep their label and never enter a mean), i.e. the procedure of the paper's 1-deg re-evaluation restricted to these
8 cases and 5 methods. Also reported: the recorded 15-deg values of the same runs.

## Endpoints and analysis (all exploratory; conventions of the paper)
Endpoint: wake loss in % of the wake-free objective, L = 100 (Ideal - Objective) / Ideal, under the 1-deg objective.
Delta L = first minus second (pp; negative = first better). Margin m = 0.05 pp. Two-sided alpha = 0.05.
1. PRIMARY comparison: PSO-VNS - PSO on the direct-1-deg arm: mean of the per-case mean differences over the cases
   in which both qualify (>= 15 of 30 runs feasible); 90% seed-level interval (within-case seed-paired bootstrap,
   mpce_inference_extra.seed_level, stratified, 10,000 resamples, seed 20260930) and 90% case-level interval
   (bootstrap over the cases, mpce_results.tost, 10,000 resamples, seed 20260928); TOST equivalence verdict at
   m = 0.05 pp at each level; case-mean Wilcoxon (mpce_inference_extra.cm_test; with 8 cases the smallest attainable
   two-sided p is 0.0078); Hodges-Lehmann. Exploratory: not part of the paper's primary equivalence questions, no
   multiplicity adjustment.
2. Feasibility-aware ranks of the five methods (mpce_results.case_stats / rank_rule: qualified = >= 15 of 30
   feasible runs, ranked by the mean objective of the feasible runs rounded to 1e-6, others below by number of
   feasible runs), average ranks over the 8 cases, Friedman / Iman-Davenport, Holm z tests vs. the best method
   (family of 4) and vs. PSO-VNS (family of 4); feasibility rates; mean wake loss of qualified methods.
   LEADER = method with the best (lowest) average rank (ties reported as ties).
3. All-run paired scores: PSO-VNS against each other method, all 240 seed pairs (8 cases x 30 seeds): a feasible run
   beats an infeasible run; two feasible runs: the higher objective wins (tie if |difference| <= 1e-6 objective
   units); two infeasible runs tie. Score = (W + T/2) / 240 with 95% case-bootstrap percentile interval (20,000
   resamples of the 8 cases, all seed pairs of a case kept together, seed 20261004); exact two-sided sign test over
   discordant pairs, Holm over the 4 comparisons.
4. Comparison with the 1-deg RE-EVALUATION (reference arm; same cases and methods, seeds 1-30): items 1-3 repeated on
   the reference arm; does the leader change (leader of direct arm vs. leader of reference arm)? Kendall tau between
   the two average-rank orders; share of cases with the same best method; per method, the per-case mean 1-deg wake
   loss of direct optimization minus that of the re-evaluated 15-deg layouts (pp; unpaired seed sets 31-60 vs 1-30,
   case-level 90% bootstrap interval) -- the gain from optimizing on the fine rose. If the optional control arm
   exists: the same difference seed-paired (seeds 31-60 both arms).
Interpretation rules fixed in advance: "PSO-VNS leads under direct 1-deg optimization" iff it has the best average
rank; "PSO-VNS better than PSO" iff both 90% intervals (seed and case level) lie below 0; "equivalent" iff the 90%
interval of a level lies inside (-0.05, 0.05) pp (stated per level); otherwise inconclusive. "Leader changes" iff
the best-average-rank method differs between the direct and the reference arm. With 8 cases, case-level statements
describe these 8 cases, not the 68-case benchmark.

## Validation (before the study; results in rev3_fine_validate.log)
(1) With the 15-deg objective swapped back, the driver reproduces stored main-comparison runs bit for bit (Objective,
    Calls, Feasible, Curve, 3-decimal coordinates) for seed 1 of two cases and each method (PSO-VNS, PSO, MS-SLSQP,
    RSD-VNS from mpce_*; GA from rev2_ga).
(2) The fast 1-deg objective of stored layouts equals mpce_direction.bench_objective(xy, ds, 15) bit for bit (and the
    stored J15 column of mpce_direction_layouts.csv.gz to its 3-decimal rounding).

## Reproducibility notes
Each shard writes its rows as they finish (resumable: keys already in the shard file are skipped). SLSQP results
depend on the SciPy / BLAS platform; NumPy's SIMD exp / pow kernels depend on the CPU, so runs are exactly
reproducible only on the same software and CPU class (the log records Python, NumPy, SciPy and the CPU flags).
Reference environment: Python 3.11, NumPy 2.4.6, SciPy 1.17.1, pandas 3.0.6.

## Execution plan (added after validation and timing; no design change)
Validation (rev3_fine_validate.log): PASSED -- (1) 10/10 stored seed-1 runs reproduced bit for bit with the 15-deg
objective (cases I-750-6 and II-1000-15; PSO-VNS, PSO, GA, MS-SLSQP, RSD-VNS); (2) 1,200/1,200 stored layouts give
bit-identical 1-deg objective values with FineObjective and mpce_direction.bench_objective(., ., 15).
Timing (rev3_fine_timing.log; measured on the shared 4-core machine at load average 10-13, so wall times are
inflated, probably by ~1.5x): one 1-deg run at N = 15 takes 13-15 s (PSO-VNS, PSO, GA, RSD-VNS) and 35-37 s
(MS-SLSQP); estimated total 3.1-4.6 CPU-hours for the 1,200 runs.
Shards: 5 (tasks[i::5]; exactly 6 of the 30 seeds of every case x method per shard, 240 runs, est. 37-55 min each):
    cd analysis && python3 rev3_fine.py <i> 5 --procs 1 > rev3_fine_s<i>of5.log 2>&1      (i = 0..4)
Optional 15-deg control arm, 3 shards (est. 2 CPU-hours in total):
    cd analysis && python3 rev3_fine.py <i> 3 --procs 1 --rose 15 > rev3_fine15_s<i>of3.log 2>&1   (i = 0..2)
Analysis after the runs: cd analysis && python3 rev3_fine_analysis.py   (-> rev3_fine.json, rev3_fine_tables.tex)

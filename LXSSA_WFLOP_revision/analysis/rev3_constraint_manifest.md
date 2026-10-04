# Revision 3, item C1 (reviewer concern R5): alternative constraint handling -- pre-specification

Written 2026-10-04 BEFORE any study run (only a smoke test of the driver, 6 single runs on seed 1, was executed to
check that the code works and to time it; those runs are not part of the study and are not stored).
Driver: `analysis/rev3_constraint.py`; analysis: `analysis/rev3_constraint_analysis.py`.

## Question
Do the method conclusions survive a constraint handling that (a) does not add boundary violations in m^2 to spacing
violations in m and (b) does not depend on the penalty factor mu = 1e10 (manuscript Section 3.4, Eq. (6))?

## Variants (only the constraint handling changes)
| id | label | objective seen by the optimizer | repair |
|----|-------|----------------------------------|--------|
| (i) | `pen` | F_p of Eq. (6): L + sum_{g_b>0} (1 + mu g_b)^2 + sum_{g_s>0} (1 + mu g_s)^2, mu = 1e10, g_b = x^2+y^2-r^2 (m^2), g_s = l_min - l_ij (m) (`authors_objective.make_objective`, unchanged) | box clipping to [-r, r] after every move, inside the operators (velocity kept in PSO) |
| (ii) | `deb` | Deb's feasibility rules on dimensionless violations, v = sum_i max(0, g_b_i)/r^2 + sum_{i<j} max(0, g_s_ij)/l_min, realized by the scalar surrogate F_deb = L if feasible, BIG + v otherwise (BIG = 2^18 = 262,144) | box clipping, unchanged |
| (iii) | `proj` | F_deb as in (ii) | box clipping inside the operators, then radial projection p <- r p/|p| of every turbine outside the circle, applied in the objective wrapper in place (written back into the optimizer's arrays: Lamarckian repair) |

"Feasible" is the study's rule throughout (identical code to `mpce_experiments.run_grid`): l_ij >= l_min - 1e-6 m
for all pairs and sqrt(x_i^2+y_i^2) <= r + 1e-6 m for all turbines, l_min = 8 R = 308 m. A layout feasible under
this tolerance gets F_deb = L even if a violation below 1e-6 m remains (v is then reported but not used).

**Order equivalence of F_deb.** Deb's rules (Deb 2000, CMAME 186:311-338): (1) feasible before infeasible;
(2) two feasible by the objective L; (3) two infeasible by the total (normalized) violation v. Every comparison of
two solutions in the five methods -- PSO personal/global best update (`<`), GA (RCGA) binary tournament (`<`),
elitism (argmin/argmax), SSA food source and population sorting (argsort, `<`), DE one-to-one replacement (`<`),
basic VNS acceptance, shaking acceptance and best-improvement local search (`<`), the hybrids' hand-over of the
Phase-1 best -- uses objective values only through comparisons; no method uses F arithmetically (no
fitness-proportional selection, averaging or differences of F). Hence any strictly order-preserving scalar
realization of Deb's order leaves every algorithm's behaviour identical to a lexicographic comparison. With
0 <= L <= ideal <= 2.107e5 < BIG and v > 0 for infeasible layouts: rule (1) holds exactly; rule (2) is exact
(F = L); rule (3) is realized by the floating-point sum BIG + v, which is monotone non-decreasing in v, so no pair is
ever reversed; two violations closer than ulp(BIG)/2 = 2^-35 ~ 2.9e-11 (dimensionless) become a tie, and a tie keeps
the incumbent (strict `<`), which is Deb's tie treatment. 2.9e-11 is far below the 1e-6 m tolerance (1e-6/l_min =
3.2e-9). `rev3_constraint.py validate` checks this on 4,000 sampled layouts per case (all pairs: reversed pairs must
be 0; ties introduced are counted). The driver asserts L < BIG and v > 0 for every infeasible evaluation.

## Methods (5)
PSO-VNS (`PSOBV`, omega = 0.5, constriction PSO phase), PSO (`PSOC`, constriction w = 0.7298, c1 = c2 = 1.49618),
GA (`GA`, `rev2_ga.RCGA`: binary tournament, SBX eta_c = 15, p_c = 0.9, polynomial mutation eta_m = 20, p_m = 1/n,
elitism), SSA-VNS (`SSABV`, omega = 0.5), DE (`DE`, DE/rand/1/bin, F = 0.5, CR = 0.9; the poorly feasible
baseline). All settings as in the paper (Table 2); N_p = 30. Code paths: `mpce_experiments.run_grid` /
`run_method` (and `rev2_ga.run_method_ga` for GA) unchanged; only `make_objective` (the objective/penalty) is
swapped, and for (iii) the projection is done inside the swapped objective.

## Cases (6, fixed)
Both data sets; dense and less dense; largest N at r = 500 and 1000 m, and one middle N:
(I, 500 m, N = 10), (I, 750 m, N = 6), (I, 1000 m, N = 15), (II, 500 m, N = 10), (II, 750 m, N = 6),
(II, 1000 m, N = 15). Density phi = N (l_min/2r)^2: 0.95 (500/10), 0.25 (750/6), 0.36 (1000/15).

## Seeds and the reference
Seeds 31-60 (new) for ALL three variants, including the reference (i). Justification: the comparison between
variants is then paired by seed and by initial population (all methods draw the same seeded initial population of
30), so variant effects are not confounded with seed sets, and the reviewer asked for new seeds 31-60 for new
studies. The stored runs of seeds 1-30 are used (a) for the bit-for-bit validation of the driver (variant (i),
seeds 1-2, all 5 methods x 6 cases = 60 runs: objective, wake loss, feasibility, calls, curve and coordinates at
the stored precision must be identical), and (b) as a descriptive replication check of (i) (seeds 1-30 stored vs
31-60 new: feasibility rates and conditional means). They do not enter the inferential comparisons.

## Budget accounting
6,030 objective evaluations per run, every evaluation counted (the existing Tracker and the optimizers' own call
counters; HybridBVNS raises at the budget), random starts (uniform in the box, `init_hook.GEN = None`); 201
convergence checkpoints of the best feasible wake loss. The surrogate and the projection cost no extra evaluations.
The driver asserts Calls = 6,030 evaluations counted by the wrapper. Size: 3 variants x 5 methods x 6 cases x
30 seeds = 2,700 runs (+ 60 validation runs), 2 worker processes.

## Outcomes (per variant; wake loss in % of the ideal, pp differences = PSO-VNS minus comparator)
1. Feasibility rate of the final layouts per method x case and pooled (180 runs; Wilson 95% interval).
2. Conditional mean wake loss (feasible runs only) per method x case; pooled as the mean over cases (cases where
   the method has no feasible run are listed as missing).
3. **Primary: all-run paired score** of PSO-VNS against each of the other four methods, per case (30 seed pairs) and
   pooled (180 pairs): a feasible run beats an infeasible one; two feasible runs are compared by wake loss (tie if
   the objectives differ by <= 1e-9 legacy units); two infeasible runs tie. Score = (W + T/2)/n (0.5 = no
   difference). 95% interval: two-stage percentile bootstrap over cases and seeds (resample the 6 cases, then the 30
   seed pairs within each resampled case; 20,000 resamples, seed 20261004); per case: seed bootstrap. Test: exact
   two-sided sign test on W vs L (ties dropped).
   **Comparison family**: for each variant separately, the 4 comparators x (6 cases + pooled) = 28 sign tests,
   Holm-adjusted within the family, alpha = 0.05 two-sided.
4. Secondary (margin): the same score with two feasible runs counted as a tie when |Delta L| <= m = 0.05 pp; and the
   mean Delta L (pp) over seed pairs in which both runs are feasible, with 95% and 90% two-stage bootstrap intervals;
   "equivalent" if the 90% interval lies within [-0.05, 0.05], "PSO-VNS better beyond the margin" if the 95%
   interval lies below -0.05. Wilcoxon signed-rank on these Delta L (Holm within the same 28-test structure).
5. Feasibility-aware average rank per variant: in each (case, seed) block the 5 methods are ranked, feasible runs
   by wake loss (1 = best), infeasible runs after all feasible ones and tied among themselves (mid-rank); average
   over the 180 blocks, with two-stage bootstrap 95% intervals; and the paper's case-level rule
   (`mpce_results.rank_rule`: mean feasible objective if >= 15 of 30 runs feasible, otherwise below, by number of
   feasible runs) averaged over the 6 cases.
6. Rank-order change: the order of the five methods by the run-level average rank in (ii) and (iii) compared with
   (i) (identical / which pairs swap; Kendall tau), and the same for the case-level rule.
7. Exploratory, not in the family: per method, change between variants on identical seeds (pooled feasibility:
   exact McNemar; Delta L of seed pairs feasible under both variants: mean with 95% bootstrap interval);
   evaluation-level diagnostics (share of feasible evaluations, first feasible evaluation, projections).

## Decision rule (stated before running)
The method conclusions are said to **survive** a variant if, relative to (i) on the same seeds, (a) the Holm
category (significantly > 0.5 / not significant / significantly < 0.5) of the pooled all-run score of PSO-VNS
against each of the four comparators is unchanged, and (b) the run-level feasibility-aware rank order of the five
methods is unchanged; a swap of two adjacent methods whose 95% rank intervals overlap is reported as "order not
resolved" rather than as a reversal. Any change in a per-case category is reported. All outcomes are reported,
including unfavourable ones.

## Outputs
`rev3_constraint_{pen,deb,proj}.csv` (columns of `run_grid` + Variant, Violation, MaxBoundExcess,
MaxSpacingDeficit, FeasEvalPct, FirstFeasCall, Projections; coordinates with 17 significant digits),
`rev3_constraint_validate.{csv,json}`, `rev3_constraint.json`, `rev3_constraint_tables.tex`
(`tab:S-r3-constraint`), `rev3_constraint_*.log`.

# Reviewer status — fresh numerical study

The numerical sections of `LXSSA_WFLOP_reviewer_revised.tex` (21 pages) were rebuilt from scratch. The earlier result tables (historical EA/ACO/PF/BBO/SSA/LX-SSA/QA-SSA values and the recorded follow-up data) are removed. All results now come from one controlled study with the original optimizer code and the original `objective.py`:
- 68 benchmark cases (Data Sets I and II; 500 m with N = 2–10, 750 m with N = 2–12, 1000 m with N = 2–15);
- six methods (LX-SSA, SSA, PSO, DE, VNS, MS-SLSQP), 30 seed-paired runs each, 6,030 objective calls per run, 12,240 runs in total;
- the Horns Rev 1 site case.

| ID | Status | Where |
|---|---|---|
| R1-1 Wake/Weibull | Addressed | Sec. III-A |
| R1-2 4D spacing | Addressed | Sec. VIII-B: capacity (Table XII), 5D/6D compliance of all fresh layouts, LX-SSA/SSA re-optimization at 4D/5D/6D (Table XIII) |
| R1-3 Power curve | Addressed | Sec. VIII-A, Table XI (cubic curve ± cut-out on all fresh layouts) |
| R1-4 Leader/follower | Addressed | Algorithms 1–2 match the code |
| R1-5 Laplace mechanism / ablation | Addressed | φ = 0, χ = 1; LX-SSA vs SSA seed-paired over 68 cases (Table III) |
| R1-6 Constraints | Addressed | Sec. III-C matches `objective.py` (box clip + quadratic penalty) |
| R1-7 Baselines | Addressed | SSA, PSO, DE, VNS, MS-SLSQP at equal budget |
| R1-8 Statistics | Addressed | Case-level Friedman (Table II), per-case Friedman + Holm–Wilcoxon + effect sizes (Table III, CSV), box plots, convergence curves |
| R1-9 Real wind | Addressed | Sec. VII-H: Horns Rev 1 (real V80, wind climate, outline; 16 and 80 turbines vs installed layout; PyWake cross-check within 0.3%) |
| R1-10 Turbine count | Addressed | Setup section |
| R1-11 Coordinates/cost | Addressed | Appendix Table XIV + CSVs for all runs; timing per run (Sec. VII-G) |
| R1-12 Presentation | Addressed | Compiles cleanly |
| R2-1 Motivation | Addressed | Contributions point to Sections VII–VIII |
| R2-2 Model novelty | Addressed | – |
| R2-3 Laplace parameters | Addressed | φ = 0, χ = 1, not tuned |
| R2-4 Significance/runtime/boxplots | Addressed | Tables II–III, Figs. 4–9 |
| R2-5 Direct benchmark | Addressed | Sec. VII |
| R3-1 Contribution | Addressed as far as evidence allows | Honest controlled comparison |
| R3-2 Jensen rationale | Addressed | Gaussian re-evaluation (Table XI) |
| R3-3 Modern baselines | Addressed | VNS and gradient-based MS-SLSQP; surrogate methods discussed only |
| R3-4 Literature | Addressed; **verify DOIs** | – |
| R3-5 Table typo | Moot | Old tables removed |

## Main results

1. **Overall ranking** over 68 cases (Friedman p = 1.8e-27): VNS 1.60, SSA 3.15, PSO 3.33, MS-SLSQP 3.87, **LX-SSA 3.97**, DE 5.08.
2. **LX-SSA vs SSA:** 0 significantly better, 59 not different, 9 significantly worse.
3. **LX-SSA vs VNS:** 0 / 27 / 41. **vs PSO:** 18 / 39 / 11. **vs MS-SLSQP:** 18 / 31 / 19. **vs DE:** 46 / 18 / 4.
4. **Feasibility:** MS-SLSQP 99.9%, VNS 99.6%, SSA 97.6%, LX-SSA 96.9%, PSO 86.0%, DE 71.4%.
5. **Horns Rev 1:** no method beats the installed layout. For 80 turbines only MS-SLSQP finds feasible layouts (1.6% below installed).
6. **Robustness:** the algorithm ordering is robust to the power curve (τ = 0.97) but not to the wake model (Gaussian τ = 0.58). VNS stays best; LX-SSA stays fifth.

## Still for the authors

- Confirm φ = 0, χ = 1.
- Verify all reference DOIs.
- Update the cover letter and response letter to the new results.

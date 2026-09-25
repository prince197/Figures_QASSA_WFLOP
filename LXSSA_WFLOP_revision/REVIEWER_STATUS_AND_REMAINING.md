# Reviewer status after integrating the authors' optimizer code

Checked against *MPCE_Reviewer_Suggestions_Incorporated_vs_Remaining.pdf*. The main source is `LXSSA_WFLOP_reviewer_revised.tex` (28 pages); it compiles with IEEEtran with no errors and no undefined references.

| ID | Status | Where in the manuscript |
|---|---|---|
| R1-1 Wake/Weibull | Addressed | Sec. III-A, Eq. (3) |
| R1-2 4D spacing | Addressed | Sec. VIII-D: geometric capacity (Table XVII), stored-layout 5D/6D compliance, LX-SSA/SSA re-optimization at 4D/5D/6D with the authors' code (Table XVIII) |
| R1-3 Linear power curve | Addressed (re-evaluation) | Sec. VIII-B, Table XV; objective → kW → AEP in Sec. III-B |
| R1-4 Leader/follower | Addressed | Algorithms 1–2 now match the code |
| R1-5 Laplace mechanism / ablation | **Addressed** | φ = 0, χ = 1 stated (Eq. 11 corrected to the code's form); seed-paired equal-budget SSA vs LX-SSA ablation, Sec. VII-B |
| R1-6 Constraints | **Addressed (corrected)** | Sec. III-C rewritten to match `objective.py` exactly: box clipping + quadratic penalty (1 + 10^10 g)^2 (was: radial projection) |
| R1-7 PSO/DE baselines | **Addressed** | Sec. VII: equal budgets 3,030 and 6,030, seed-paired |
| R1-8 Statistics | Addressed | Run-level Friedman, Holm–Wilcoxon signed-rank, rank-biserial (Table XII); recorded-data tests (Tables IX–X) |
| R1-9 Real wind | **Addressed** | Sec. VII-E (Horns Rev 1 wind rose, benchmark turbine) and Sec. VII-F (Horns Rev 1 site: real V80 power/Ct curves, measured wind climate, real outline, 16- and 80-turbine cases vs installed layout; evaluator cross-checked against PyWake within 0.3%) |
| R1-10 Turbine count | Addressed | – |
| R1-11 Coordinates/cost | Addressed | Coordinates (Table XIX + CSVs); individually timed runs with hardware/software (Sec. VII-E) |
| R1-12 Presentation | Addressed | – |
| R2-1 Motivation | Addressed | Contribution bullets point to evidence |
| R2-2 Model novelty | Addressed | – |
| R2-3 Laplace parameters | **Addressed** | φ = 0, χ = 1, interpretation, not tuned |
| R2-4 Significance/runtime | Addressed | Sec. VII-B/E |
| R2-5 Direct benchmark | Addressed | Sec. VII |
| R3-1 Contribution | Addressed as far as evidence allows | Honest equal-budget result; feasibility advantage of SSA family |
| R3-2 Jensen rationale | Addressed (+ Gaussian re-evaluation) | Sec. VIII-C |
| R3-3 Modern baselines | **Addressed** | Sec. VII-C: VNS and multistart SLSQP (gradient-based) at equal budgets, seed-paired; surrogate methods discussed only |
| R3-4 Literature | Addressed; **verify DOIs** | 8 new references incl. PyWake, Mladenović & Hansen (VNS), Kraft (SLSQP) |
| R3-5 Table typo | Addressed | – |

## Key results the authors must read

All controlled runs use the authors' exact `objective.py` (quadratic penalty) and optimizer code; calibration against the recorded runs: 22/24 comparisons p ≥ 0.05, means within 0.21%.

1. **At equal budgets, LX-SSA does not outperform SSA or PSO.** SSA has the higher mean in all 6 cases at both budgets; LX-SSA is never significantly better than SSA. LX-SSA beats DE consistently. Its lead in the recorded data comes from its doubled evaluation budget.
2. **VNS beats LX-SSA** significantly in 5 of 6 cases at both budgets and has the best average rank of all six methods. Multistart SLSQP ranks 5th of 6 (the Jensen/24-bin objective is non-smooth).
3. **Crowded farms:** LX-SSA/SSA are feasible far more often than PSO/DE, but VNS and SLSQP are more reliable still (27–30 of 30).
4. **Horns Rev 1 site:** no method beats the installed layout at 6,030 calls (16-turbine block: VNS 137.97 vs installed 139.51 GWh/yr). For the full 80-turbine farm, none of LX-SSA, SSA, PSO, DE or VNS finds a feasible layout; only SLSQP does (1.6% below installed).
5. **Runtime:** ~0.41 ms per call for LX-SSA and SSA; the recorded 4–15× LX-SSA slowdown is not reproduced.

## To confirm or supply

- **φ = 0, χ = 1** were used for the recorded runs (strongly supported by calibration, but please confirm).
- **Dataset II, 500 m, N = 3** table entry (flagged with †).
- DOIs of all references.
- Cover letter (not yet updated for these results).

# Reviewer status after integrating the authors' optimizer code

Checked against *MPCE_Reviewer_Suggestions_Incorporated_vs_Remaining.pdf*. The main source is `LXSSA_WFLOP_reviewer_revised.tex` (26 pages); it compiles with IEEEtran with no errors and no undefined references.

| ID | Status | Where in the manuscript |
|---|---|---|
| R1-1 Wake/Weibull | Addressed | Sec. III-A, Eq. (3) |
| R1-2 4D spacing | Addressed | Sec. VIII-D: geometric capacity (Table XVII), stored-layout 5D/6D compliance, LX-SSA/SSA re-optimization at 4D/5D/6D with the authors' code (Table XVIII) |
| R1-3 Linear power curve | Addressed (re-evaluation) | Sec. VIII-B, Table XV; objective → kW → AEP in Sec. III-B |
| R1-4 Leader/follower | Addressed | Algorithms 1–2 now match the code |
| R1-5 Laplace mechanism / ablation | **Addressed** | φ = 0, χ = 1 stated (Eq. 11 corrected to the code's form); seed-paired equal-budget SSA vs LX-SSA ablation, Sec. VII-B |
| R1-6 Constraints | **Addressed (corrected)** | Sec. III-C rewritten to match the code: box clipping + 10^10 penalty (was: radial projection) |
| R1-7 PSO/DE baselines | **Addressed** | Sec. VII: equal budgets 3,030 and 6,030, seed-paired |
| R1-8 Statistics | Addressed | Run-level Friedman, Holm–Wilcoxon signed-rank, rank-biserial (Table XII); recorded-data tests (Tables IX–X) |
| R1-9 Real wind | **Addressed (partly)** | Sec. VII-D: Horns Rev 1 measured wind rose; benchmark turbine and circular farm kept |
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
| R3-3 Modern baselines | Partly | PSO/DE benchmarked; VNS/gradient/surrogate only discussed |
| R3-4 Literature | Addressed; **verify DOIs** | 6 new references incl. PyWake |
| R3-5 Table typo | Addressed | – |

## Key results the authors must read

1. **At equal budgets, LX-SSA does not outperform SSA or PSO.** At 3,030 calls it ties SSA in all 6 cases; at 6,030 calls SSA is significantly better in 2 cases and PSO in 3. LX-SSA beats DE consistently. Its lead in the recorded data comes from its doubled evaluation budget. The abstract, discussion and conclusion now say this.
2. **Clear positive result:** in the most crowded cases, LX-SSA/SSA return feasible layouts in 22–30 of 30 runs; PSO and DE in 0–2.
3. **Horns Rev 1:** same pattern (SSA best mean in 3 of 4 cases; LX-SSA tied with SSA except 1000 m/15).
4. **Runtime:** at equal calls, all four algorithms take ~0.4 ms per call; the recorded 4–15× LX-SSA slowdown is not reproduced.

## To confirm or supply

- **Penalty form in `objective.py`.** I assumed wake loss + 10^10 × (total spacing shortfall + total boundary excess). If the code uses a different form, correct Eq. (7) and re-run `authors_experiments.py`.
- **φ = 0, χ = 1** were used for the recorded runs (strongly supported by calibration, but please confirm).
- **Dataset II, 500 m, N = 3** table entry (flagged with †).
- DOIs of all references.
- Cover letter (not yet updated for these results).

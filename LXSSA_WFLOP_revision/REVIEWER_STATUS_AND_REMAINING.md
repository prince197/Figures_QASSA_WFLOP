# Reviewer status — hybrid LX-SSA-VNS manuscript

`LXSSA_WFLOP_reviewer_revised.tex` (23 pages) now proposes **LX-SSA-VNS**, a two-phase hybrid: the published LX-SSA (Solanki & Deep 2023, φ = 0, χ = 1) explores for the first 50% of the budget (3,030 calls), then variable neighbourhood search intensifies the LX-SSA food source for the remaining 3,000 calls (Sec. VI, Algorithm 3).

All results come from one controlled study with the original optimizer code and the original `objective.py`:
- 68 benchmark cases (Data Sets I and II; 500 m with N = 2–10, 750 m with N = 2–12, 1000 m with N = 2–15);
- eight methods (LX-SSA-VNS, LX-SSA, SSA, PSO, DE, modified VNS, MS-SLSQP, and the **original basic VNS = BVNS**), 30 seed-paired runs each, 6,030 objective calls per run, 16,320 runs in total. Earlier methods were **not** re-run: new methods reuse the same seeds, so they are paired with the existing runs;
- a budget-split ablation (25% / 50% / 75% to LX-SSA; 12 cases, 720 extra runs);
- the Horns Rev 1 site case (16 turbines × 30 seeds, 80 turbines × 10 seeds).

| ID | Status | Where |
|---|---|---|
| R1-1 Wake/Weibull | Addressed | Sec. III-A |
| R1-2 4D spacing | Addressed | Sec. IX-B: capacity (Table XIV), 5D/6D compliance, LX-SSA/SSA re-optimization at 4D/5D/6D (Table XV) |
| R1-3 Power curve | Addressed | Sec. IX-A, Table XIII (cubic curve ± cut-out, all eight methods) |
| R1-4 Leader/follower | Addressed | Algorithms 1–3 match the code |
| R1-5 Laplace mechanism / ablation | Addressed | φ = 0, χ = 1; ablation of the hybrid (Sec. VIII-G, Table XI); original vs modified VNS (Sec. VIII-F, Table X) |
| R1-6 Constraints | Addressed | Sec. III-C matches `objective.py`; both hybrid phases use it |
| R1-7 Baselines | Addressed | LX-SSA, SSA, PSO, DE, original BVNS, modified VNS, MS-SLSQP at equal budget |
| R1-8 Statistics | Addressed | Case-level Friedman (Table II), per-case Friedman + Holm–Wilcoxon + effect sizes (Table III, CSV), box plots, convergence curves |
| R1-9 Real wind | Addressed | Sec. VIII-J: Horns Rev 1 (real V80, wind climate, outline; 16 and 80 turbines vs installed layout) |
| R1-10 Turbine count | Addressed | Setup section |
| R1-11 Coordinates/cost | Addressed | Appendix Table XVI + CSVs for all runs; timing (Sec. VIII-I) |
| R1-12 Presentation | Addressed | Compiles cleanly |
| R2-1 Motivation | Addressed | Sec. II (motivation for the hybrid) |
| R2-2 Model novelty | Addressed | Novelty is the hybrid algorithm, not the wake model |
| R2-3 Laplace parameters | Addressed | φ = 0, χ = 1, not tuned; split ρ examined |
| R2-4 Significance/runtime/boxplots | Addressed | Tables II–III, Figs. 4–12 |
| R2-5 Direct benchmark | Addressed | Sec. VIII |
| R3-1 Contribution | Addressed | New hybrid + honest controlled comparison |
| R3-2 Jensen rationale | Addressed | Gaussian re-evaluation (Table XIII) |
| R3-3 Modern baselines | Addressed | Original BVNS, modified VNS and gradient-based MS-SLSQP; surrogate methods discussed only |
| R3-4 Literature | Addressed; **verify DOIs** | – |
| R3-5 Table typo | Moot | Old tables removed |

## Main results

1. **Overall ranking** over 68 cases (Friedman χ² = 212.5, p = 2.6e-42): modified VNS 2.09, **LX-SSA-VNS 2.74**, **BVNS (original VNS) 3.62**, SSA 4.79, PSO 4.84, MS-SLSQP 5.38, LX-SSA 5.68, DE 6.86. Hybrid vs VNS (Holm p = 0.123) and hybrid vs BVNS (p = 0.068) not significant at case level; VNS significantly better than BVNS (p = 5.1e-4); hybrid significantly better than the other five (p ≤ 2.9e-6).
2. **Hybrid run-level W/T/L (68 cases, Holm over 7):** vs LX-SSA 35/33/0, vs SSA 27/41/0, vs PSO 32/36/0, vs DE 52/15/1, vs MS-SLSQP 34/29/5, vs VNS 0/53/15, **vs BVNS 12/49/7**.
3. **Original BVNS (Table X):** better than LX-SSA in 32 cases, SSA 24, PSO 30, DE 53, MS-SLSQP 24; vs modified VNS 8 better (all N = 3–4) / 27 worse (all N ≥ 6); only method feasible in all 2,040 runs.
4. **Feasibility:** BVNS 100%, MS-SLSQP 99.9%, VNS 99.6%, LX-SSA-VNS 98.1%, SSA 97.6%, LX-SSA 96.9%, PSO 86.0%, DE 71.4%.
5. **Ablation:** the VNS phase removes 36.4% (median 31.1%) of the wake loss left by the LX-SSA phase and repairs 82 infeasible runs. More budget to LX-SSA (75%) is worse in 11/12 cases (4 significant); more to VNS (25%) slightly better in 9/12 (1 significant).
6. **Horns Rev 1 (16 turbines):** VNS 137.97, BVNS 137.67, hybrid 137.59 GWh/yr; hybrid better than LX-SSA and MS-SLSQP, not significantly different from SSA, VNS, BVNS; installed layout 139.51. For 80 turbines only MS-SLSQP is feasible.
7. **Robustness:** ordering robust to the power curve (τ = 0.96), not to the wake model (Gaussian τ = 0.57); VNS and LX-SSA-VNS keep the two best ranks under all models; BVNS drops to 4.62 under the Gaussian model.

## Honest limitation to keep in the paper

The modified VNS is still the best method at this budget, and the hybrid is not significantly better than the original BVNS at case level. The hybrid clearly improves LX-SSA, but its gain comes from the VNS phase. Reviewers may ask why LX-SSA is needed at all; the paper says this openly (Secs. VIII-F, VIII-G, XI, XII).

## Still for the authors

- Confirm φ = 0, χ = 1.
- Verify all reference DOIs.
- Update the cover letter and response letter to the new framing (proposed hybrid).

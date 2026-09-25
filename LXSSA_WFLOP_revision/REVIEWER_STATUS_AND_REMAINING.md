# Reviewer status — hybrid LX-SSA-VNS manuscript

`LXSSA_WFLOP_reviewer_revised.tex` (22 pages) now proposes **LX-SSA-VNS**, a two-phase hybrid: the published LX-SSA (Solanki & Deep 2023, φ = 0, χ = 1) explores for the first 50% of the budget (3,030 calls), then variable neighbourhood search intensifies the LX-SSA food source for the remaining 3,000 calls (Sec. VI, Algorithm 3).

All results come from one controlled study with the original optimizer code and the original `objective.py`:
- 68 benchmark cases (Data Sets I and II; 500 m with N = 2–10, 750 m with N = 2–12, 1000 m with N = 2–15);
- seven methods (LX-SSA-VNS, LX-SSA, SSA, PSO, DE, VNS, MS-SLSQP), 30 seed-paired runs each, 6,030 objective calls per run, 14,280 runs in total. The six reference methods were **not** re-run: the hybrid runs reuse the same seeds, so they are paired with the existing runs;
- a budget-split ablation (25% / 50% / 75% to LX-SSA; 12 cases, 720 extra runs);
- the Horns Rev 1 site case (16 turbines × 30 seeds, 80 turbines × 10 seeds).

| ID | Status | Where |
|---|---|---|
| R1-1 Wake/Weibull | Addressed | Sec. III-A |
| R1-2 4D spacing | Addressed | Sec. IX-B: capacity (Table XIII), 5D/6D compliance, LX-SSA/SSA re-optimization at 4D/5D/6D (Table XIV) |
| R1-3 Power curve | Addressed | Sec. IX-A, Table XII (cubic curve ± cut-out, all seven methods) |
| R1-4 Leader/follower | Addressed | Algorithms 1–3 match the code |
| R1-5 Laplace mechanism / ablation | Addressed | φ = 0, χ = 1; ablation of the hybrid (Sec. VIII-F, Table X) |
| R1-6 Constraints | Addressed | Sec. III-C matches `objective.py`; both hybrid phases use it |
| R1-7 Baselines | Addressed | LX-SSA, SSA, PSO, DE, VNS, MS-SLSQP at equal budget |
| R1-8 Statistics | Addressed | Case-level Friedman (Table II), per-case Friedman + Holm–Wilcoxon + effect sizes (Table III, CSV), box plots, convergence curves |
| R1-9 Real wind | Addressed | Sec. VIII-I: Horns Rev 1 (real V80, wind climate, outline; 16 and 80 turbines vs installed layout) |
| R1-10 Turbine count | Addressed | Setup section |
| R1-11 Coordinates/cost | Addressed | Appendix Table XV + CSVs for all runs; timing (Sec. VIII-H) |
| R1-12 Presentation | Addressed | Compiles cleanly |
| R2-1 Motivation | Addressed | Sec. II (motivation for the hybrid) |
| R2-2 Model novelty | Addressed | Novelty is the hybrid algorithm, not the wake model |
| R2-3 Laplace parameters | Addressed | φ = 0, χ = 1, not tuned; split ρ examined |
| R2-4 Significance/runtime/boxplots | Addressed | Tables II–III, Figs. 4–12 |
| R2-5 Direct benchmark | Addressed | Sec. VIII |
| R3-1 Contribution | Addressed | New hybrid + honest controlled comparison |
| R3-2 Jensen rationale | Addressed | Gaussian re-evaluation (Table XII) |
| R3-3 Modern baselines | Addressed | VNS and gradient-based MS-SLSQP; surrogate methods discussed only |
| R3-4 Literature | Addressed; **verify DOIs** | – |
| R3-5 Table typo | Moot | Old tables removed |

## Main results

1. **Overall ranking** over 68 cases (Friedman χ² = 190.9, p = 1.7e-38): VNS 1.85, **LX-SSA-VNS 2.38**, SSA 4.07, PSO 4.12, MS-SLSQP 4.71, LX-SSA 4.91, DE 5.95. Hybrid vs VNS not significant at case level (Holm p = 0.153); hybrid significantly better than the other five (Holm p ≤ 1e-5).
2. **Hybrid run-level W/T/L (68 cases):** vs LX-SSA 36/32/0, vs SSA 27/41/0, vs PSO 33/35/0, vs DE 52/15/1, vs MS-SLSQP 34/29/5, **vs VNS 0/52/16** (all 16 losses at N ≥ 7).
3. **Feasibility:** MS-SLSQP 99.9%, VNS 99.6%, LX-SSA-VNS 98.1%, SSA 97.6%, LX-SSA 96.9%, PSO 86.0%, DE 71.4%.
4. **Ablation:** the VNS phase removes 36.4% (median 31.1%) of the wake loss left by the LX-SSA phase and repairs 82 infeasible runs. More budget to LX-SSA (75%) is worse in 11/12 cases (4 significant); more to VNS (25%) slightly better in 9/12 (1 significant).
5. **Horns Rev 1 (16 turbines):** hybrid 137.59 GWh/yr — better than LX-SSA, SSA, MS-SLSQP; worse than VNS (137.97); installed layout 139.51. For 80 turbines only MS-SLSQP is feasible.
6. **Robustness:** ordering robust to the power curve (τ = 0.96), not to the wake model (Gaussian τ = 0.58); VNS and LX-SSA-VNS keep the two best ranks under all models.

## Honest limitation to keep in the paper

Stand-alone VNS is still the best method at this budget. The hybrid clearly improves LX-SSA, but its gain comes from the VNS phase. Reviewers may ask why LX-SSA is needed at all; the paper says this openly (Secs. VIII-F, XI, XII).

## Still for the authors

- Confirm φ = 0, χ = 1.
- Verify all reference DOIs.
- Update the cover letter and response letter to the new framing (proposed hybrid).

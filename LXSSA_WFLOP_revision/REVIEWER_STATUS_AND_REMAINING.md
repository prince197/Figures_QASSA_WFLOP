# Reviewer status — final manuscript (hybrid LX-SSA-VNS with the original VNS)

`LXSSA_WFLOP_reviewer_revised.tex` (21 pages) proposes **LX-SSA-VNS**: the published LX-SSA (Solanki & Deep 2023, φ = 0, χ = 1) explores for the first 50% of the budget (3,030 calls), then the **original basic VNS** (Mladenović & Hansen 1997; continuous form of Mladenović et al. 2008: l∞-shell shaking of the whole layout + complete best-improvement compass local search) intensifies the LX-SSA food source for the remaining 3,000 calls (Sec. VI, Algorithm 3). The earlier modified VNS is no longer used anywhere in the paper.

All results come from one controlled study with the original optimizer code and the original `objective.py`:
- 68 benchmark cases (Data Sets I and II; 500 m with N = 2–10, 750 m with N = 2–12, 1000 m with N = 2–15);
- seven methods (LX-SSA-VNS, LX-SSA, SSA, PSO, DE, VNS = original basic VNS, MS-SLSQP), 30 seed-paired runs each, 6,030 calls per run (14,280 runs);
- component ablation: SSA, LX-SSA, VNS, LX-SSA-VNS on all 68 cases, and budget split 25/50/75% on 12 cases (+720 runs);
- Horns Rev 1 site case (16 turbines × 30 seeds, 80 turbines × 10 seeds).

| ID | Status | Where |
|---|---|---|
| R1-1 Wake/Weibull | Addressed | Sec. III-A |
| R1-2 4D spacing | Addressed | Sec. IX-B: capacity (Table XIV), 5D/6D compliance, LX-SSA/SSA re-optimization at 4D/5D/6D (Table XV) |
| R1-3 Power curve | Addressed | Sec. IX-A, Table XIII (cubic curve ± cut-out, all seven methods) |
| R1-4 Leader/follower | Addressed | Algorithms 1–3 match the code |
| R1-5 Laplace mechanism / ablation | Addressed | φ = 0, χ = 1; component ablation SSA / LX-SSA / VNS / LX-SSA-VNS (Table X, Fig. 10) and budget split (Table XI) |
| R1-6 Constraints | Addressed | Sec. III-C matches `objective.py`; both hybrid phases use it |
| R1-7 Baselines | Addressed | LX-SSA, SSA, PSO, DE, original VNS, MS-SLSQP at equal budget |
| R1-8 Statistics | Addressed | Case-level Friedman (Table II), per-case Friedman + Holm–Wilcoxon + effect sizes (Table III, CSV), box plots, convergence curves |
| R1-9 Real wind | Addressed | Sec. VIII-I: Horns Rev 1 (real V80, wind climate, outline; 16 and 80 turbines vs installed layout) |
| R1-10 Turbine count | Addressed | Setup section |
| R1-11 Coordinates/cost | Addressed | Appendix Table XVI + CSVs for all runs; timing (Sec. VIII-H) |
| R1-12 Presentation | Addressed | Compiles cleanly |
| R2-1 Motivation | Addressed | Sec. II (motivation for the hybrid) |
| R2-2 Model novelty | Addressed | Novelty is the hybrid algorithm, not the wake model |
| R2-3 Laplace parameters | Addressed | φ = 0, χ = 1 (confirmed by the authors as the original code values), not tuned; split ρ examined |
| R2-4 Significance/runtime/boxplots | Addressed | Tables II–III, Figs. 4–13 |
| R2-5 Direct benchmark | Addressed | Sec. VIII |
| R3-1 Contribution | Addressed | New hybrid + honest controlled comparison |
| R3-2 Jensen rationale | Addressed | Gaussian re-evaluation (Table XIII) |
| R3-3 Modern baselines | Addressed | Original basic VNS and gradient-based MS-SLSQP; surrogate methods discussed only |
| R3-4 Literature | Addressed; **verify DOIs** | – |
| R3-5 Table typo | Moot | Old tables removed |

## Main results

1. **Overall ranking** over 68 cases (Friedman χ² = 168.4, p = 9.9e-34): **LX-SSA-VNS 1.89**, VNS 2.85, SSA 3.85, PSO 4.14, MS-SLSQP 4.44, LX-SSA 4.76, DE 6.07. The hybrid is significantly better than every method at case level (vs VNS Holm p = 0.010; others p ≤ 2.3e-7).
2. **Hybrid run-level W/T/L (Holm over 6):** vs LX-SSA 39/29/0, SSA 32/36/0, PSO 36/32/0, DE 52/16/0, VNS 6/60/2, MS-SLSQP 36/28/4.
3. **Feasibility:** VNS 100%, LX-SSA-VNS 99.9%, MS-SLSQP 99.9%, SSA 97.6%, LX-SSA 96.9%, PSO 86.0%, DE 71.4%.
4. **Ablation (Table X):** ranks LX-SSA-VNS 1.50, VNS 2.15, SSA 2.81, LX-SSA 3.54 (Friedman χ² = 103.7, p = 2.5e-22). VNS phase: LX-SSA-VNS vs LX-SSA 43/25/0. LX-SSA start vs best initial point: LX-SSA-VNS vs VNS 10/56/2. Laplace step alone: LX-SSA vs SSA 0/59/9. VNS phase removes 39.0% (median 31.9%) of the remaining wake loss and repairs 118 infeasible runs. Budget split: insensitive (one significant difference in 12 cases).
5. **Horns Rev 1 (16 turbines):** LX-SSA-VNS 137.77 GWh/yr (highest mean, 30/30 feasible), VNS 137.67 (n.s.), better than LX-SSA, SSA, MS-SLSQP; installed layout 139.51. For 80 turbines only MS-SLSQP is feasible.
6. **Robustness:** ordering robust to the power curve (τ = 0.97), not to the wake model (Gaussian τ = 0.54); LX-SSA-VNS keeps the best average rank under all models.

## Honest limitation to keep in the paper

The Laplace step on its own does not improve on SSA (LX-SSA vs SSA 0/59/9); the paper says so. SSA-VNS runs exist (`SSABV` in `fresh_bgrid.csv`) but, at the authors' request, are not reported in the manuscript; a reviewer may ask for an SSA + VNS comparison.

## Author confirmations and remaining steps

- φ = 0, χ = 1 confirmed by the authors as the values of the original LX-SSA code.
- The paper is submitted as a new manuscript; no response-to-reviewers letter is prepared. New cover letter: `Cover_letter.tex`; competing-interest declaration: `conflict_of_interest.tex`.
- References: checked against publisher records (see git history for the corrections).
- Before submission: one final proofread by the authors and, if required by the journal, a statement on AI-assisted editing.

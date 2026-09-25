# Reviewer status after the September 2026 integration pass

Checked against *MPCE_Reviewer_Suggestions_Incorporated_vs_Remaining.pdf*. The main source is `LXSSA_WFLOP_reviewer_revised.tex`; it compiles with IEEEtran with no undefined references.

| ID | Status now | What this pass added (main .tex) | Still outstanding (needs author data or new runs) |
|---|---|---|---|
| R1-1 Wake/Weibull | Addressed | – | Final check against code |
| R1-2 4D spacing | **Addressed (new experiment)** | Sec. VII-D: geometric capacity at 4D/5D/6D (Table XIII), 5D/6D compliance of 720 stored layouts, reference-SSA re-optimization at 4D/5D/6D (Table XIV) | Re-run with the original LX-SSA code if desired |
| R1-3 Linear power curve | **Addressed (re-evaluation)** | Sec. VII-B, Table XI: cubic curve ± cut-out, % change in expected power, algorithm ordering; objective→kW→AEP conversion in Sec. III-B | Manufacturer curve (optional) |
| R1-4 Leader/follower | Addressed | Lemma 1 rewritten compactly (identical formula) | Final check against code |
| R1-5 Laplace mechanism | Partial | Text now states the published phi, chi are retained without retuning and no sensitivity was done | **Numerical phi, chi**; equal-budget matched SSA vs LX-SSA ablation |
| R1-6 Constraints | Addressed | Hard-coded "Section III-B" replaced by `\ref` | – |
| R1-7 PSO/DE baselines | Partial (data integrated) | Runtime column, budget-effect quantification (Sec. VII-D) | **Equal-budget rerun** |
| R1-8 Statistics | **Addressed (exploratory)** | Case-level Friedman + Holm post hoc (Table X), bootstrap 95% CIs, rank-biserial/A12, list of non-significant cases | Historical-table raw vectors |
| R1-9 Real wind | **Pending** | – | Measured wind rose (e.g. Horns Rev 1) with LX-SSA; needs phi, chi |
| R1-10 Turbine count | Addressed | – | – |
| R1-11 Coordinates/cost | **Mostly addressed** | Appendix B coordinate table (Table XV) + CSVs; batch runtimes; timed reference runs with hardware/software | Hardware of archived runs; historical coordinates |
| R1-12 Presentation | Addressed | Proofreading: ω indexing, Table I refs, section numbering, Huang/Alba typos, wrong intro physics claim, overfull equations | Final author read |
| R2-1 Motivation | Addressed | Contribution bullets now point to tables/sections | – |
| R2-2 Model novelty | Addressed | – | – |
| R2-3 Laplace parameters | Partial | Fallback statement added (published settings retained, no retuning) | **Numerical phi, chi** |
| R2-4 Significance/runtime | **Addressed (exploratory)** | CIs, effect sizes, Friedman, runtime column, explicit tie cases | Independent timings of archived runs |
| R2-5 Direct benchmark | Partial | as R1-7 | Equal-budget rerun |
| R3-1 Contribution | Partial → improved | Robustness section + statistics | Real wind, equal budget |
| R3-2 Jensen rationale | **Addressed (+ optional check done)** | Sec. VII-C, Table XII: Gaussian (Bastankhah–Porté-Agel) re-evaluation, k* = 0.04 / 0.022 | Re-optimization under Gaussian (optional) |
| R3-3 Modern baselines | Partial | as R1-7 | VNS/gradient discussed only |
| R3-4 Literature | Addressed | 5 new references (Carrillo 2013, Bastankhah 2014, Demšar 2006, Holm 1979, Vargha 2000) | **Verify all DOIs** (Crossref/publisher sites were not reachable from the build environment) |
| R3-5 Table typo | Addressed | – | – |

## Findings the authors must review before resubmission

1. **Calibration gap.** An independent SSA written strictly from Eqs. (8)/(9) scores significantly higher than the *recorded* SSA **and** LX-SSA at the same 3,030-evaluation budget in all six cases (e.g. +339 and +543 above the recorded LX-SSA mean at 750 m/8). The manuscript reports this in Sec. VII-D. Please check the original implementation (initialization, bound handling, leader step). If you find the cause and re-run, this paragraph should be revised.
2. **Dataset II, 500 m, N = 3.** The entries 42137.21 exceed the ideal 21947.06. They are flagged with † in Table III, not changed. Replace them with the original record (very likely 21947.06 / 0, as in the 750 m table).
3. **LX-SSA runtime.** Recorded batch time per run is 0.14 s for LX-SSA versus ~0.01–0.04 s for the others. That is 4–15× slower, more than the 2× evaluation count explains. It is reported descriptively.
4. **Population size.** Evaluation counts (3030/6030) imply a population of 30 in the follow-up runs, not N_p = 5·2N. This is stated in Sec. VI-D.

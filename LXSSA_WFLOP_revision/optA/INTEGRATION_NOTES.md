# Integration notes (lead) — items for Phase 2 / Phase 4

- 01_front Nomenclature: E_ideal should read "wake-free benchmark objective of one turbine" (not expected power) — from 03_model agent.
- 03_model says a feasible layout always beats an infeasible one for violations > ~1e-7; 05_setup's final-feasibility tolerance is 1e-6 m — make wording consistent.
- 03_model's last sentence references \label{sec:setup} — 05_setup must keep it.
- 09_robust: \TBD{macro: ...} for boundary-handling stats → pipeline agent adds \NBoundDiffMin/\NBoundDiffMax/\NBoundMaxP/\NBoundCases; replace in 09 during Phase 2.
- 07_ablation is a top-level \section (sec:ablation); 06_results told to drop "and the ablation" and trim C11 overlap.
- 07 optional: run-level contrasts SSA-VNS vs SSA, LX-SSA-VNS vs LX-SSA (not needed).
- Cover letter pulls macros from analysis/mpce_numbers.tex; must stay in LXSSA_WFLOP_revision/.
- New refs cited by 02_intro: Sorensen2015, BartzBeielstein2020 → merge via optA/merge_newrefs.py.

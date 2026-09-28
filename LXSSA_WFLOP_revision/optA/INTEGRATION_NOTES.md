# Integration notes (lead) — items for Phase 2 / Phase 4

- 01_front Nomenclature: E_ideal should read "wake-free benchmark objective of one turbine" (not expected power) — from 03_model agent.
- 03_model says a feasible layout always beats an infeasible one for violations > ~1e-7; 05_setup's final-feasibility tolerance is 1e-6 m — make wording consistent.
- 03_model's last sentence references \label{sec:setup} — 05_setup must keep it.
- 09_robust: \TBD{macro: ...} for boundary-handling stats → pipeline agent adds \NBoundDiffMin/\NBoundDiffMax/\NBoundMaxP/\NBoundCases; replace in 09 during Phase 2.
- 07_ablation is a top-level \section (sec:ablation); 06_results told to drop "and the ablation" and trim C11 overlap.
- 07 optional: run-level contrasts SSA-VNS vs SSA, LX-SSA-VNS vs LX-SSA (not needed).
- Cover letter pulls macros from analysis/mpce_numbers.tex; must stay in LXSSA_WFLOP_revision/.
- New refs cited by 02_intro: Sorensen2015, BartzBeielstein2020 → merge via optA/merge_newrefs.py.
- 02_intro: claims without CHECK ids: "VNS alone ranks first with feasible starts" (\NFbFeasBestRank) and "best rank at all three budgets" (\NBudgetPhrase) → consider adding C28 (feasible-start best = BVNS) and C29 (PSOBV best rank at 30,030 and 120,030) in mpce_check_final.py; outline/contribution 4 must point to the discussion label if 10 creates sec:discussion; Horns Rev claim wording must match (C19 is at 6,030; 30,030 needs its own check).
- Overfull box in 04_methods (equation ~line 37) — check in integration.
- 10_limits_concl: new label sec:discussion → 02_intro outline/contribution 4 should reference sec:discussion for the recommendations.
- 10: "best published IEA37 layouts obtained with gradient-based methods" — only participant 4 (SNOPT+WEC per Baker 2019, unverified) is known; reword to "the best feasible published layout (reported as a gradient-based method, \TBD{verify})" or drop "gradient-based".
- 10: "behind salp-swarm hybrids" plural — fine if both LXBV and SSABV rank ahead of old PSO in tab:baseline (check C27 numbers).
- Pipeline asked for C28 (feasible-start: BVNS 1st, PSOBV 2nd), C29 (PSOBV best at all budgets), C30 (HR16 30k mean > installed; macros \NHRThirtyKMean, \NHRThirtyKAbove).
- 04_methods: new labels sec:vns, sec:template; constriction factor symbol is κ (χ = Laplace scale); references sec:psosetting (06) — 06 must keep it; possible duplication with 05_setup (iteration counts, RS-VNS/LX-SSA-VNS details) → trim in 05 during integration.

# Phase 4 — fixing the eight reviews (all agents must follow this, together with optA/BRIEF.md)

Reviews: optA/reviews/R1_editor.md, R2_metaheuristics.md, R3_statistics.md, R4_windenergy.md,
R5_claims_audit.md, R6_copyedit.md, R7_reproducibility.md, R8_supplement.md. Earlier lead notes:
optA/INTEGRATION_NOTES.md. Read ALL eight reviews and fix every item that concerns YOUR file(s).
BRIEF.md hard rules still apply (own files only, macros only for numbers, no pdflatex in the shared folder,
use `sh optA/compile_copy.sh <private scratch dir>`, no commit/push, ≤ 12 pages overall).

## Lead decisions (these override BRIEF.md where they differ)

D1 — M2 is restated; the old wording was contradicted by the paper's own case-mean test.
  Old: "salp-swarm phases add nothing beyond random sampling; the Laplace step gives no gain; VNS helps every swarm".
  New: "Against a random-sampling control (RS-VNS), an SSA phase gives only a small gain (tied at run level,
  \NAblSSAVNSvsRSVNS; case-mean difference \NCmSSAVNSvsRSVNSDL pp, p = \NCmSSAVNSvsRSVNSP), an LX-SSA phase gives
  no gain (p = \NCmLXSSAVNSvsRSVNSP), and the Laplace step makes the hybrid worse (SSA-VNS better than LX-SSA-VNS,
  p = \NCmSSAVNSvsLXSSAVNSP). A PSO phase gives a clearly larger gain. The VNS phase significantly improves both
  salp-swarm methods (\NAblSSAVNSvsSSA, \NAblLXSSAVNSvsLXSSA); for PSO its gain is small and not significant over
  all 68 cases (p = \NPWPSO)."
  Also say plainly that RS-VNS is a weak, "no-search" control: uniform samples in the bounding square are rarely
  feasible (all N turbines inside the circle with probability (π/4)^N), so "better than RS-VNS" is a minimum bar.
  Thesis sentence (abstract/intro/conclusion): "... while salp-swarm phases add little beyond random sampling and the
  Laplace step makes the hybrid worse."
D2 — M3 limits restated. Do NOT write "at least as good as PSO". Write: PSO-VNS has the best average rank and is
  feasible in every run, but over the 68 cases it is not significantly different from PSO (p = \NPWPSO,
  run level \NWtlPSO, i.e. PSO is significantly better in \NWtlPSOL cases). Its advantage is concentrated in
  Data Set II with N ≥ 10 (\NCmPSOVNSvsPSOdsIILargeWins of \NCmPSOVNSvsPSOdsIILargeN cases lower, p =
  \NCmPSOVNSvsPSOdsIILargeP); in Data Set I with N ≥ 10 PSO has the lower loss in \NCmPSOVNSvsPSOdsILargeLosses of
  \NCmPSOVNSvsPSOdsILargeN cases. Call the N ≥ 10 split exploratory (post hoc). Do not write "gain grows with N".
  Budget claim: "best average rank over the six largest cases at all three budgets (exploratory, no test; first
  in \NBudFirstPSOVNSOneTwentyK of 6 cases at 120,030 evaluations)".
D3 — M1 wording. The old setting w = 0.7, c1 = c2 = 2 satisfies the order-1 (mean) stability condition but violates
  the order-2 (variance) condition of Poli (φ = c1 + c2 = 4 > 24(1−w²)/(7−5w) = 3.50); write "violates the
  order-2 stability condition" (not "non-convergent" as a bare statement; "does not converge in variance" is fine).
  Note that these conditions assume stagnation and unbounded search, and that in the code velocities are kept when
  positions are clipped to the box, which together plausibly explains the low feasibility. State that the ranking
  reversal is shown at 6,030 evaluations; at 120,030 evaluations SSA-VNS ranks ahead of PSO (\NBudRankSSAVNSOneTwentyK
  vs \NBudRankPSOOneTwentyK) — report it, do not hide it. Use "changes the ranking" rather than "reverses".
D4 — Statistics (R3). State in 05: the 68 cases are not independent (6 farm radii with nested N; two data sets), so
  p-values are descriptive of this benchmark, not of a population of farms; Holm correction for run-level tests is
  within each case (family = the comparisons of that case), W/T/L are descriptive tallies; the case-mean test treats
  a case in which only one method has ≥ 15 feasible runs as a maximal difference in its favour (\NCmImputed cases)
  and drops cases where neither has (\NCmDropped); rank-biserial excludes zero differences. Report Iman–Davenport
  p (\NImanP). Report effect sizes with 95 % bootstrap CIs over cases (\NCm...CI) where a comparison is a headline
  claim. The two different PSO-VNS vs PSO tallies (main table vs ablation table) come from different Holm families
  — the caption/text must say so. "p_Holm ≤ x" statements must use macros rounded upward.
D5 — Horns Rev (R4 found a direction-binning bug; the model is fixed and ALL Horns Rev runs are being rerun).
  Use only the existing \NHR... macros; they will be regenerated from the new runs. Installed 16-turbine block AEP
  and wake-free AEP change (\NHRInstalled, \NHRIdeal). Do not state any Horns Rev result as a hard claim except
  through macros with a CHECK-FINAL condition; the lead rewrites Horns Rev sentences after the reruns. Validation vs
  PyWake for the 80-turbine farm: \NHRPyWakeDiff (about 0.6 %, no longer 0.3 %). The abstract must not say
  PSO-VNS "exceeds" the installed layout unless the mean does; phrase: "its best run (and \NHRAbovePSOVNS of 30
  runs) exceed the installed block".
D6 — Cost (R5 false claim): replace "all methods need a similar time" by the per-evaluation times
  (\NMsEvalMin–\NMsEvalMax ms for the metaheuristics; MS-SLSQP \NMsEvalSLSQP ms, about \NSLSQPSlowdown× PSO).
D7 — Budget split: ω = 1 (plain PSO at the same budget) is added to the split table as the end point
  (\NSplitRankHundred, \NSplitLossHundred); 0.75 remains "best of the tested values on 12 cases", 0.5 remains the
  pre-specified default. No new runs.
D8 — Reproducibility (R7): fix the mismatches (MS-SLSQP minimizes the wake loss with explicit constraints, not F_p;
  8.7e-11 is an absolute difference; RS-VNS uses 3,015 samples vs 3,030 swarm calls — say "about the same";
  feasibility tolerance consistent = 1e-6 m everywhere), add MS-SLSQP settings (maxiter 200, ftol 1e-9, scaled
  constraints, restarts from new random layouts after the 30 initial ones), Horns Rev settings (boundary +0.1 %,
  distance-outside penalty, box half-width 943 m) and feasible-initialization parameters — main text only if short,
  otherwise supplement. Data availability must not promise "one command regenerates everything"; say which script
  produces which table (list lives in the supplement / repository README).
D9 — New references (R2): only these, and only with details you can confirm in optA/reviews/R2_metaheuristics.md:
  Cleghorn & Engelbrecht 2018 (PSO stability), Bratton & Kennedy 2007 (standard PSO), Castelli et al. 2022 (SSA
  critique), Zaharie 2002 (DE critical values), Hansen & Ostermeier 2001 (CMA-ES, for the limitation "no CMA-ES /
  L-SHADE baseline"), Benavoli et al. 2017 (Bayesian comparison). Put new \bibitem entries in
  optA/newrefs/<yourfile>.tex (the lead merges them). Max 6 new refs overall; prefer the first four.
D10 — Keep the page budget: the whole main PDF ≤ 12 pages. Every addition must be paid for by a cut in the same file.
  Move detail to the supplement (tell the supplement agent via your report, the lead coordinates).

## Macro contract (defined by the pipeline agent in analysis/mpce_numbers.tex; placeholders in
## optA/phase4_macros_stub.tex print [TBD] until then — use these exact names)
- \NCm<A>vs<B>{P,Wins,Losses,DL,CI}: case-mean Wilcoxon p (same procedure as \NPW...), number of cases where A has
  the lower / higher mean loss, mean difference A − B of the case-mean wake losses in pp (signed, negative = A
  better), 95 % bootstrap CI of that mean as "[a, b]". Pairs <A>vs<B>: SSAVNSvsRSVNS, LXSSAVNSvsRSVNS,
  PSOVNSvsRSVNS, SSAVNSvsLXSSAVNS, LXSSAvsSSA, SSAVNSvsSSA, LXSSAVNSvsLXSSA, PSOVNSvsPSO, PSOVNSvsVNS.
- \NAbl<A>vs<B>{,W,T,L,DL} (run-level, ablation family) for new pairs SSAVNSvsSSA, LXSSAVNSvsLXSSA, LXSSAVNSvsRSVNS.
- \NCmPSOVNSvsPSO<G>{P,Wins,Losses,N} for G = Large (N ≥ 10, both data sets), dsILarge, dsIILarge;
  \NWtlPSOdsILargeW/L, \NWtlPSOdsIILargeW/L (significant run-level wins/losses of PSO-VNS vs PSO in the subgroup).
- \NImanP; \NMsEvalMin, \NMsEvalMax, \NMsEvalSLSQP, \NSLSQPSlowdown; \NSplitRankHundred, \NSplitLossHundred,
  \NSplitHundredVsSeventyFive (run-level W/T/L of ω = 0.75 vs ω = 1); \NBudFirstPSOVNSSixK/ThirtyK/OneTwentyK;
  \NHRPyWakeDiff; \NCmImputed, \NCmDropped.

## Report
List: review items fixed (by reviewer and number), items deliberately not fixed (and why), labels
added/removed, figures/tables moved, new refs, macros you needed that are not in the contract, compile result
(pages) of your private compile.

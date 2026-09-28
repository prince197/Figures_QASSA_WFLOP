# Option A brief — reframing of the MPCE manuscript (all agents must follow this)

The manuscript (`MPCE_PSO_VNS.tex`) is split into section files `optA/01_front.tex` … `optA/11_back.tex`.
Each agent owns exactly ONE file (named in its task) and must not edit any other file.
No new experiments: every result already exists. Numbers come ONLY from macros in
`analysis/mpce_numbers.tex` (generated from `analysis/mpce_summary.json`). Never type a result number by hand;
if a number you need has no macro, write `\TBD{macro: <what>}` and list it in your report.

## Why we reframe
Rigor is strong, but a paper whose headline is "new algorithm PSO-VNS wins" is weak: PSO-VNS is not significantly
better than a correctly configured PSO over all 68 cases (case-mean Wilcoxon p = \NPWPSO), VNS alone ranks first
with feasible starts, PSO ranks first under the Gaussian re-evaluation, and on IEA37 the methods stay below the
gradient-based published optima. The reframed paper claims exactly what the evidence supports.

## New thesis (one sentence)
In continuous wind farm layout optimization, the conclusions of metaheuristic comparisons depend strongly on how the
baselines are configured and how the comparison is controlled; under a controlled protocol, a simple two-phase design —
a convergent PSO followed by basic VNS (PSO-VNS) — is the most consistent performer, while salp-swarm phases and the
Laplace step add nothing beyond random sampling.

## Title (use in 01_front; other agents refer to "this paper", never to the title)
"Revisiting Metaheuristic Comparisons for Continuous Wind Farm Layout Optimization: Baseline Configuration,
Component Analysis and a Two-Phase PSO–VNS Design"
(short running head: "SOLANKI et al.: REVISITING METAHEURISTIC COMPARISONS FOR WIND FARM LAYOUT OPTIMIZATION")

## Three messages (in this order everywhere: abstract, introduction, results order, conclusion)
M1 — Baseline configuration changes the conclusions. The PSO setting w = 0.7, c1 = c2 = 2, common in WFLOP
comparisons, violates the particle convergence condition. With it, PSO ranks \NOldPSOPos-th (feasible in \NOldPSOFeas\%
of runs) and salp-swarm hybrids appear superior; with the Clerc–Kennedy constriction setting PSO ranks \NPSOPos-nd and
outperforms all salp-swarm methods. A new table `\input{analysis/mpce_tab_baseline.tex}` (label `tab:baseline`,
produced by the pipeline agent) shows the ranking of the previous study's method pool (LX-SSA-VNS, SSA-VNS, LX-SSA,
SSA, PSO, DE, VNS, MS-SLSQP) under the old and the corrected PSO setting. Its macros (defined by the pipeline agent):
`\NBaseOldBest`, `\NBaseNewBest` (name of best method in that pool), `\NBaseOldPSORank`, `\NBaseNewPSORank`,
`\NBaseOldLXBVRank`, `\NBaseNewLXBVRank`, `\NBaseOldSSAVNSRank`, `\NBaseNewSSAVNSRank`,
`\NBaseOldLXBVvsPSO`, `\NBaseNewLXBVvsPSO`, `\NBaseOldSSAVNSvsPSO`, `\NBaseNewSSAVNSvsPSO` (run-level W/T/L strings).
State this carefully: as a property of the setting, never as criticism of specific papers.
M2 — What matters in two-phase hybrids (component analysis). The VNS phase helps every swarm; the swarm phase helps
only if the swarm is effective: SSA and LX-SSA phases are statistically tied with random sampling (RS-VNS), the Laplace
step gives no gain, a PSO phase beats random sampling, VNS alone and the salp phases.
M3 — PSO-VNS as a simple, strong reference method, with its limits stated plainly: best average rank (\NRankPSOVNS),
at least as good as PSO overall and better for N >= 10, best rank at all three budgets, feasible in every run, beats the
installed Horns Rev 1 block at 30,030 evaluations; BUT not significantly better than PSO over all cases, VNS alone is
best with feasible starts (PSO and DE stall there — explain why), PSO ranks first under the Gaussian re-evaluation, and
on IEA37 the gap to gradient-based published optima remains.

## Contributions (for 02_intro; exactly these four)
1. A controlled evaluation protocol for continuous WFLOP (equal budgets counting every objective and gradient call,
   seed-paired runs, feasibility-aware ranking, case-mean and run-level tests with Holm correction, ablation with a
   random-sampling control, robustness re-evaluation, public code and data) applied to 8 methods, 68 benchmark cases,
   a Horns Rev 1 block and IEA37 Case Study 1.
2. Evidence that the PSO baseline setting reverses the ranking of methods (M1).
3. A component analysis of two-phase swarm–VNS hybrids (M2).
4. PSO-VNS as a simple reference method, with an honest account of where it wins and where it does not (M3), and
   practical recommendations for benchmarking WFLOP methods.

## Structure (files) and page budget — the whole PDF must stay ≤ 12 pages, main text ending on p. 11 at the latest;
aim for 11 pages total. Cut rather than add. Budgets are for the typeset two-column text incl. floats.
- 01_front: title, authors, abstract (≤ 230 words, no citations, the three messages with a few key macros), index
  terms (5–6: Wind farm layout optimization; Particle swarm optimization; Variable neighborhood search; Benchmarking;
  Metaheuristics; Wake effect), Nomenclature (keep compact; drop symbols no longer used in the main text).
- 02_intro (≤ 1.0 page): motivation for MPCE readers (layout decisions fix energy yield; many metaheuristics proposed,
  comparisons hard to trust), short focused literature (keep MPCE refs), the gap (baselines, budgets, feasibility,
  statistics), the thesis, the four contributions, paper outline.
- 03_model (≤ 1.2 pages): keep the model; tighten wording; move anything derivational that is not needed to the
  supplement only if already there (do not delete content that has no other home).
- 04_methods (≤ 1.3 pages): PSO with constriction (+ convergence condition, M1 hook), SSA/LX-SSA briefly, basic VNS,
  PSO-VNS algorithm box; frame SSA/LX-SSA/RS-VNS as the phase-1 alternatives compared in the component analysis.
- 05_setup (≤ 0.9 page): the protocol as a first-class element ("Evaluation protocol"): budgets, seeds/pairing,
  feasibility-aware ranking rule, statistics, ablation design incl. RS-VNS control, beyond-benchmark tests.
- 06_results (≤ 2.0 pages): order = (A) Effect of the baseline setting (M1, new tab:baseline, existing
  "Effect of the PSO Coefficients" content merged here), (B) overall comparison (existing tab:friedman68, tab:wtl,
  figures), (C) solution quality/feasibility/convergence (condense).
- 07_ablation (≤ 0.9 page): component analysis (M2) + budget split (50% pre-specified, 75% better — report honestly).
- 08_beyond (≤ 1.3 pages): budget scaling, feasible initialization (keep the stall explanation), Horns Rev 1, IEA37
  (state the gap to published optima plainly), cost.
- 09_robust (≤ 0.5 page): power curve / wake model / spacing / boundary.
- 10_limits_concl (≤ 0.9 page): rename to "Discussion and Recommendations" + "Limitations" + "Conclusion". Add a
  concise list "Recommendations for benchmarking WFLOP methods" (5–7 items: convergent baseline settings; equal budgets
  incl. gradient calls; seed pairing; feasibility-aware ranking; case-mean + run-level tests with multiplicity
  correction; random-sampling control in ablations; report feasible-start and model-robustness checks; open code).
- 11_back: acknowledgment, data availability, bibliography (IEEE order), biographies.

## Hard rules
- Keep every `\input{analysis/...}` table and every figure unless your section budget forces a move; if you move a
  table/figure to the supplement, say so in your report (the lead handles the supplement).
- Keep every label (`\label{...}`) that other files may reference; if you remove one, list it in your report.
- Keep every `% CHECK-FINAL [Cxx]` comment attached to its sentence; if you rewrite the sentence, keep the claim within
  the condition. You may reuse existing IDs; a new ID C27 (baseline table) is defined by the pipeline agent:
  "C27: under the old PSO setting at least one salp-swarm hybrid ranks ahead of PSO; under the corrected setting PSO
  ranks ahead of every salp-swarm method (SSA, LX-SSA, SSA-VNS, LX-SSA-VNS)".
- New citations: do NOT edit 11_back.tex; put the complete IEEE-style `\bibitem{Key} ...` into
  `optA/newrefs/<yourfile>.tex` and cite `Key`. Only cite works you are sure exist (title, venue, year, DOI).
- Terminology: PSO-VNS, PSO (always the constriction setting unless "old setting" is said), SSA, LX-SSA, SSA-VNS,
  LX-SSA-VNS, RS-VNS, VNS, MS-SLSQP; "wake loss"; "benchmark objective"; "evaluations" for objective calls;
  "feasibility-aware ranking".
- Tone: measured, no hype words (novel, superior, outstanding, significantly unless statistically), IEEE style,
  authors' plain style, present tense for findings.
- Do not commit or push. To test-compile, run from `LXSSA_WFLOP_revision/`:
  `sh optA/compile_copy.sh <your private scratch dir>` (never pdflatex in the shared folder).
- Report: what you changed, labels removed/added, figures/tables moved, new refs, missing macros, compile result.

# R8: Supplementary material and main/supplement consistency

Scope: MPCE_PSO_VNS_supplement.pdf (18 pp.) and .tex, analysis/mpce_supplementary.tex (generated tables), MPCE_PSO_VNS.pdf, optA/*.tex, main-text generated tables analysis/mpce_tab_*.tex. Both logs: no undefined or multiply defined references. All 21 `\ref{S-...}` in the main text resolve to the right content (S1 wind data, S-I.B SSA pseudocode, S-II.A PSO setting, S2-S7 per-case, Fig. S3, S-II.D, S10 split, Fig. S9/S14 budget, S12 feasinit, S-VII literature, S18 cost, S-V.A Gauss, S20/S21 spacing, S19 via the Table IX caption, Figs. S12-S14, S-VI.A power). No "Further Generated Tables" section is printed, so every generated table is placed. There are no leftovers of the older framing (no "proposed" LX-SSA-VNS/SSA-VNS, no Horns Rev 80 turbines), and "PSO" always means the constriction setting.

## Issues

1. [major] analysis/mpce_supplementary.tex (Table S12) — "DE & nan & 3.252 & 0/30 & 30/30" — raw `nan` is printed in 5 cells (Table S12, p. 10) wherever DE has no feasible random-start run. — Print "--", as Tables S13/S14/S15 already do, in the feasinit-cases writer of mpce_results.py.

2. [major] MPCE_PSO_VNS_supplement.tex S-VII — "[TBD: author: literature table pending verification of the published values against the original papers" — the literature section is only a one-sentence placeholder with a red TBD. The main text (08_beyond.tex l. 47) points to it and carries its own \TBD. — Add the verified table, or drop S-VII and the main-text sentence before submission. Never submit with a TBD.

3. [major] analysis/mpce_supplementary.tex (Tables S11, S13) vs analysis/mpce_tab_feasbudget.tex (Table VI) — "average rank over the 7 cases" — the supplement summaries cover the 6 largest cases plus Horns Rev, but main Table VI and the text cover the 6 largest cases only. As a result, the same quantities differ: PSO-VNS ranks 1.43/1.14/2.43 (S13) vs 1.50/1.17/2.17 (VI), PSO feasibility is 95% vs 100.0%, and PSO-VNS random-start loss is 4.509 (S11) vs 4.084 (VI). — Use the same case set as Table VI, or drop S11/S13 (neither is cited from the main text) and keep only the per-case Tables S12/S14.

4. [major] MPCE_PSO_VNS_supplement.tex S-II.A — "Table~\ref{M-tab:baseline} ranks the method pool of the previous version" — the headline baseline result (old PSO feasible in 86.0% of runs, worse in 56/68 cases, +0.61 pp) has no per-case support in the supplement. S-II.A describes only the protocol. — Add a per-case table for old-setting vs constriction PSO (means, feasible runs, Wilcoxon p), or add an "old PSO" column to Tables S2-S7.

5. [major] MPCE_PSO_VNS_supplement.tex S-III — "uses the same 68 cases and runs as the main comparison, plus RS-VNS and LX-SSA-VNS" — "Component Analysis: Details" holds only the split table. The per-case values of LX-SSA-VNS and RS-VNS, the per-case contrasts behind Table IV, and the switch-point feasibility/loss statistics quoted in Section VII are all missing. — Add a per-case ablation table (the 2 extra variants plus contrast W/T/L per case) and the switch statistics.

6. [major] analysis/mpce_supplementary.tex (Table S21) — "Minimum-spacing re-optimization with the original LX-SSA and SSA code ... $p$ for LX-SSA vs.\ SSA" — this is a leftover of the LX-SSA framing: the only spacing re-optimization tests LX-SSA against SSA, a question the new paper does not ask. PSO-VNS and PSO are absent. — Re-run the six cases with PSO-VNS and PSO (at least), or relabel the table as an archived check and say why it is limited to SSA/LX-SSA.

7. [major] optA/09_robust.tex l. 15 — "an otherwise identical SSA that projects every turbine radially ... in all six representative cases" — this boundary-rule claim (Mann-Whitney, p<0.01, 3,030 evaluations) has no table or figure in the supplement. It is also again SSA-centred and uses another test and budget. — Add a small table (6 cases: clip vs radial projection, means, p), or cut the claim.

8. [minor] Supplement, unreferenced items — Tables S8, S9, S15, S16, S17 and Figs. S1, S2, S4-S8, S10, S11 are never cited from the main text. Only "Section S-II.D" and Fig. S3 are cited. — Cite Table S8 in VI-B (p_W, r_rb), Figs. S5-S7 in VI-C, Table S15/Fig. S10 in VIII-B, Tables S16-S17/Fig. S11 in VIII-C, and Table S9/Fig. S8 where layouts are discussed. Otherwise drop the items.

9. [minor] analysis/mpce_supplementary.tex (Table S18) and 08_beyond.tex l. 47 — "all methods need a similar time (0.24--1.82~ms per evaluation" — a 7.5x range is not "similar". The caption adds "Times come from the machine that produced each result file", so the ms/call values are not comparable across methods. — Re-time all methods on one machine, or state that the timings are indicative and say "within a factor of about 8 (MS-SLSQP slowest)".

10. [minor] MPCE_PSO_VNS_supplement.tex Fig. S3/S4 captions — "median (line) and interquartile range (band) of the best feasible wake loss over 30 runs" — DE is missing from all Fig. S3 panels, and curves start late, but neither is explained. The log y-axis is not stated. PSO (teal) and MS-SLSQP (green) are hard to tell apart. — Add "curves start once at least 15 runs are feasible; DE has no feasible median at the largest N; log scale", and use more distinct colours.

11. [minor] figures_mpce/layouts_max.pdf (Fig. S8) — "best of all methods: MS-SLSQP (best PSO-VNS: +0.21 pp loss)" — the two-line panel titles overlap the neighbouring titles in the top row (p. 8) and are unreadable. — Shorten the titles (e.g., "Best: MS-SLSQP (+0.21 pp)") or add horizontal spacing.

12. [minor] figures_mpce/budget_scaling.pdf (Fig. S9) — `\pipefig{0.7\textwidth}` — 7 panels plus a legend at 0.7 text width. Tick labels ("6k 30k 120k") and axis labels are far too small to read in print. — Use full \textwidth and a 2x4 grid with larger fonts.

13. [minor] analysis/mpce_supplementary.tex (Table S12) — "PSO & 138.06 & 136.74 & 19/30 & 30/30 & 0.560 & +0.23" — for Horns Rev, feasible starts give PSO a *lower* AEP, yet r_rb is positive ("positive = feasible initialization better"). The caption does not say how infeasible random runs enter the paired test. — State the pairing rule for infeasible runs (e.g., "infeasible runs rank last"), or check the sign.

14. [minor] analysis/mpce_supplementary.tex (Tables S11/S12) — "PSO ... 5.955" and "DE ... 5.955" — PSO and DE have identical feasible-start values in every case. The main text explains why (neither moves from the best initial layout), but the tables do not. — Add a footnote: "PSO and DE never improve on the best initial layout with feasible starts (Section VIII-A); their values coincide."

15. [minor] MPCE_PSO_VNS_supplement.tex Alg. S1 — "evaluate both candidates and keep the better one ... Apply box repair, evaluate $F_p$ for the population" — the pseudocode re-evaluates the followers after they were already evaluated. The "2N_p evaluations" per LX-SSA iteration only works out with this double count, and T (200 for SSA vs 100 for LX-SSA at 6,030) is not stated. — Make the evaluation count explicit per line and give T for both methods.

16. [minor] analysis/mpce_supplementary.tex (Table S9) and Fig. S8 caption — "of the best layout of all methods and of the best PSO-VNS layout" — in 5 of 6 farms only one row or marker set appears, because PSO-VNS is itself best. — Add "(one entry when PSO-VNS is the best method)".

17. [minor] MPCE_PSO_VNS_supplement.tex Fig. S7 caption — "Percentage of feasible runs versus the number of turbines." — the caption is not self-contained: it omits the methods, the 30 runs, the 6,030 evaluations and the random initialization. The Fig. S2 caption does not say that DE curves stop once DE has too few feasible runs. — Complete both captions in the style of Fig. S3.

18. [minor] Terminology, supplement tables — "6,030 objective calls" / "calls" (19 times) vs "evaluations" in the supplement prose and main text. Axes read "Objective-function calls" while captions say "versus evaluations". — Use "evaluations" throughout, as the intro defines, and relabel the figure axes.

19. [minor] analysis/mpce_supplementary.tex (Table S19) — "Gaussian wake ($k^*=0.04$)" — the symbol is K_G in eq. (S1) and in the main text, and Table S15 uses "Jensen wake $k=0.04$" for the benchmark K. — Use K_G and K consistently.

20. [minor] analysis/mpce_supplementary.tex (Tables S21, S2-S7) — "DS & $r$ (m) ... 1 & 500" and column "Ideal" / "ideal objective" — data sets are coded 1/2 here but I/II elsewhere (S9, S10, S12). "Ideal" is used where the main text says "wake-free benchmark objective". — Use I/II and "wake-free" everywhere.

21. [minor] analysis/mpce_supplementary.tex (Tables S10, S11-S16) — "(25\%, 50\% and 75\% of the 6,030 calls for PSO)" / "$p$ (25)" — the split is given in % while the main text uses ω=0.25/0.5/0.75. The method order also differs: LX-SSA-VNS and RS-VNS are appended last in S11-S16 but grouped with the hybrids in main Table VI. — Label the split columns ω=0.25/0.5/0.75 and use the main-text method order.

22. [minor] MPCE_PSO_VNS_supplement.tex S-V heading — "Robustness Checks: Details" — the main Section IX is titled "Sensitivity to the Modeling Assumptions". Main sections are numbered "VI-A", but the supplement uses "S-II.A". — Rename to "Sensitivity to the Modeling Assumptions: Details", and consider "S-II-A" to match the main style.

23. [minor] MPCE_PSO_VNS_supplement.tex Fig. S12 — "The points show an optimized layout for 12 turbines in the 750-m farm (Wind Data Set I)" — the figure comes from the old pipeline (figures_final/), and neither the method nor the version that produced the layout is named. — Name the method, or regenerate it from the PSO-VNS layout of Table S9.

24. [minor] analysis/mpce_supplementary.tex (Table S11 caption) — "Per-case values: Supplementary Table~\ref{tab:feasinit-cases}" — the supplement refers to itself as "Supplementary Table". The Table S17 `\resizebox` also shrinks the text below scriptsize and is barely legible. — Say "Table S12". Split S17 into two rows or use landscape instead of resizing.

25. [minor] Layout and length — at 18 pages, the length is reasonable. However, Fig. S1 (8 average ranks at full column width) repeats main Table II, Fig. S10 (right) repeats main Fig. 2 (right), and p. 11 is half empty. — Drop Fig. S1 (or shrink it to 0.5\columnwidth), and check float placement around Table S13/Fig. S9.

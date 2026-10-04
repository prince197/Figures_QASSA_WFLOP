# Revision 3, topic `numbers`: traceability of the numbers in the revision-2 manuscript

Read-only audit of `SWEVO_rev2/latex_source/` (main text: `optA/swevo_front.tex`, `optA/sw/02_intro.tex` ...
`10_limits_concl.tex`; supplement: `SWEVO_supplement.tex`, `optA/sw/supp_theory.tex`) against the repository's
generated files and data. No source, script or CSV was modified.

Files: `analysis/rev3_numbers_audit.py` (the audit), `analysis/rev3_numbers_audit.json` (all items),
`analysis/rev3_numbers_theory_check.py` (copy of `make_theory_figures.py` that checks the revision-2 theory text and
writes its figure to a scratch directory), `analysis/rev3_numbers_checks.log` (check-suite and verification logs).
Reproduce: `python3 analysis/rev3_numbers_audit.py [--rev2-out DIR]` (about 15 s, plus 15 s for
`rev2_analysis.py` if DIR holds no regenerated `rev2_summary.json`).

## 1. Mismatches and wrong numbers (most important first)

| # | Where (rev-2 file:line) | Text | Correct value / source | Assessment |
|---|---|---|---|---|
| 1 | `optA/sw/08_beyond.tex:271` and `SWEVO_supplement.tex:443` (prints the hand-edited `\NMsEval*` macros of `SWEVO_rev2/latex_source/analysis/mpce_numbers.tex`) | "the metaheuristics need 0.24--0.47 ms and MS-SLSQP 1.82 ms, about 7.5 times as long as PSO" | Generated macros (median of Seconds/Calls over the 2,040 runs per method, `mpce_results.py` `cost_per_eval`): **0.23--0.45 ms, 0.53 ms, 2.3 times**. The revision-2 values are the **means** (Table S-cost, `summary["cost"]["ms_per_call"]`: PSO 0.2444, DE 0.4668, MS-SLSQP 1.8202); the mean ratio is 1.8202/0.2444 = **7.45, i.e. "7.4"** (7.5 only arises from the rounded 1.820/0.244 = 7.459) | Hand edit, not generated. The two statistics differ because the MS-SLSQP time per evaluation is bimodal: median 0.42--0.51 ms for N <= 7 and 2.65--3.57 ms for N >= 8 (2.1--2.5 times PSO for N <= 7, 9--14 times for N >= 8). Both are correct statistics of the same records; the text must name one. Either restore the generated median values (0.23--0.45, 0.53, 2.3; then the text should say "median" and not point to Table S-cost, which gives means), or keep the means, say "mean (Table S-cost)" and write "about 7.4 times"; with the means, the generator `mpce_numbers.py` no longer reproduces `mpce_numbers.tex` (rerunning it reverts the edit). Saying the N dependence ("2.1--2.5 times PSO for N <= 7, 9--14 times for N >= 8") would be the most informative. |
| 2 | `optA/sw/08_beyond.tex:271` | "MS-SLSQP ... takes ... 2.3 and 1.6 times as long on the Lillgrund block at 6,030 and 30,030 evaluations" | median 12.788 / 5.684 s = **2.2498 -> 2.2** (and 1.596 -> 1.6) | Double rounding (Table S-time prints 2.25, rounded again to 2.3). Write 2.2 (or 2.25). |
| 3 | `optA/sw/08_beyond.tex:130` | "PyWake ... gives 661.39 GWh/yr and ours 666.75 GWh/yr (0.81% more; wake loss 11.11% vs. 10.39%)" | ours 10.39%, PyWake 11.11% (`\NFPyWakeOurLossPct`, `\NFPyWakeLossPct`; the generated text read "wake loss 10.39% vs 11.11%") | Revision 2 swapped the order; inside the parenthesis that starts with "ours ... 0.81% more" it now reads as ours = 11.11%. Values correct, attribution ambiguous; write "wake loss 10.39% vs. 11.11% (PyWake)". The same construction is in `SWEVO_supplement.tex:845` (Lillgrund: "our model 116.10 GWh/yr (2.38% more; AEP loss 18.65% vs. 16.72%)" = PyWake vs. ours; values verified with PyWake 2.6.20, see 2c). |
| 4 | `SWEVO_supplement.tex` tab:S-eqclus (cluster-mean column) | RSD-VNS vs. RS-VNS: -0.071, -0.042; PSO-VNS vs. RSD-VNS: -0.233 | from the run records: **-0.070 (-0.07035), -0.043 (-0.04250), -0.234 (-0.23353)** | Third-decimal slips (the table was built from rounded per-case table values); all estimates, CIs and verdicts of the table reproduce exactly. |
| 5 | `SWEVO_supplement.tex:841` | "parallelogram (area 1.081 km^2, perimeter 4.24 km)", Oler bound 15.7 | exact corner parallelogram: **1.079 km^2**, 4.24 km, 15.71; boundary actually used (enlarged by 0.2%): 1.084 km^2, 4.25 km, 15.75 (`rev2_site_model.site`) | Area matches neither polygon; conclusion (< 16 turbines at 4D) unchanged. Write 1.079. |
| 6 | `SWEVO_supplement.tex:953` (text before tab:S-time) | "The records of PSO-VNS on IEA37 belong to the missing original archive" | The records exist (`analysis/mpce_iea16p_s0of1.csv`, `mpce_iea36p_s0of1.csv`); `SWEVO_supplement.tex:775` itself quotes their mean times (1.3, 5.7, 3.2, 17.1 s; reproduced exactly) | Wrong statement. Medians from those records (PSO-VNS 16T: 1.31 / 5.64 s; 36T: 3.17 / 17.20 s at 6,030 / 30,030) could be added to tab:S-time. |
| 7 | `optA/sw/08_beyond.tex:154` (CHECK-DIR [F23]) | comment refers to `\NFHRInstDrop` 0.32 and `\NFHRMeanDrop` 1.08 | values are printed correctly in l. 148 | The GA paragraph inserted at l. 150 separates the comment from its paragraph (the only CHECK comment detached by revision 2). Move the comment block of l. 151--154 after l. 148. |

Checked and correct although they look suspicious:
* `optA/sw/10_limits_concl.tex:41` "lie 1.93% (16 turbines) and 6.23% (36 turbines) below": the generated text printed
  "-1.93% ... below" (double negative); the revision-2 wording is the correct fix.
* `optA/sw/08_beyond.tex:200` (gradient paragraph) "about 3--6 times those of PSO-VNS": ratios of the mean times are
  3.3--6.1 (medians 2.3--5.6), so "about 3--6" holds on the basis the supplement uses (l. 775).
* "1.82" in `08_beyond.tex:271` is not the average rank 1.82 of Table 7 (coincidence of values).

No other number of the main text or the supplement differs from its source (sections 2 and 2c).

## 2. Part 1: CHECK comments and the macros they name

208 CHECK comments in the 13 files; 57 of them name 130 macros (with or without backslash). Values: repository
`analysis/mpce_numbers*.tex` (2,782 macros; the revision-2 copies differ only in the five `\NMsEval*`/`\NSLSQPSlowdown`
values and the wording of `\NXMultLostHolm`).

* 72 numeric macros: value printed in the comment's paragraph (same format); 6 text macros found verbatim.
* 3 "undefined": `\NCm...`, `\NIEA...`, `\NHR...` are wildcard prefixes in the comments, not macros.
* 2 "no paragraph": `06_results.tex:168` (`\NRankPSOVNS ... \NRankDE`) documents a commented-out figure.
* 6 text macros not in the paragraph: `\NXClCaseSigNotUnanimous`, `\NXLocoChanged` (05_setup:87, 06_results:159),
  `\NAblNContrasts`/`\NAblNVar` (07_ablation:66, 83): the paragraphs state them qualitatively ("qualify only the
  significant contrasts of the salp-swarm hybrids with the controls", "leaving out one cluster changes none of them",
  "nine variants" in the table note); consistent.
* 41 numeric "not in paragraph" (none is a wrong number):
  * detached comment (item 7 above): `\NFHRInstDrop`, `\NFHRMeanDrop` (08_beyond:154 -> printed at l. 148).
  * printed in an adjacent paragraph: `\NXMarginMWhTurbDsI/DsII` 4.1/2.1 (05_setup:83 -> l. 79), `\NDShakesPSOVNS`
    (04_methods:174), `\NAblDf` 8 (07_ablation:66 -> table note l. 61), `\NEqMargin` 0.05 (07_ablation:72),
    `\NDRsFeasNFifteen` 0 (07_ablation:87 -> l. 81).
  * printed elsewhere in the same file: `\NSBoundC` 1.7486 (04_methods:73 -> l. 59), `\NDVnsRunsPSOVNS` 120
    (05_setup:51 -> l. 104), `\NRsSpacingDominates` 488 (07_ablation:82 -> l. 86), `\NSwitchFeasSSAVNS` 95.7
    (07_ablation:91 -> l. 104), `\NSFeasVMax` 92.2 (10_limits:73 -> l. 8), `\NBudCloseGapOneTwentyK` 0.17 (10_limits:80 -> l. 18).
  * not printed in the file (the comment documents the value behind a qualitative sentence; the repository macro
    version did not print them either): 02_intro:22 (`\NDDescentSharePSOVNS` 97.5, `\NDCycleEvalsPctPSOVNS` 21.5,
    `\NDCycleSharePSOVNS` 0.5), 02_intro:24 (`\NFPairCIOne`), 02_intro:25 (`\NXHetBeyondPSOVNS` 11, `\NXHetBeyondPSO` 8),
    02_intro:26 (`\NXHRMcNemarP`, `\NCmPSOVNSvsPSOdsIILarge{Wins,N}`, `\NXHolmPSOVNSvsPSOdsIILargeP`), 04_methods:73
    (`\NSFeasCOneEight/Nine`, `\NSVmaxFrac`), 04_methods:152 (`\NRsSpacingDominates`), 04_methods:174
    (`\NDDescentSharePSOVNS`, `\NDCycleAllSharePSOVNS`), 05_setup:51 (`\NDDynRuns`, `\NDFsRunsPSO`, `\NDFsRunsDE`),
    05_setup:87 (`\NXClSSAVNSvsRSVNSCRP`), 07_ablation:99 (`\NFAblSSAVNSvsRSVNSOneP`), 07_ablation:110
    (`\NBaySSAVNSvsRSVNSBetter`), 08_beyond:66 (`\NDFsCandDE`), 10_limits:10 (`\NSFeasBelow`, `\NSCBelow`), 10_limits:19
    (`\NFRankOrderOne`), 10_limits:33 (`\NFeasPSO`). The qualitative sentences were read and agree with the values.

## 3. Part 2: all literal numbers

Extraction from the non-comment text (bibliography and author block skipped), 2,623 numbers in the main text and
3,382 in the supplement. The revision-2 files were aligned token by token with the repository's macro versions of the
same files (macros expanded with the repository values), so every number standing where a macro stood is compared with
that macro ("anchored").

| Class | Main | Supplement |
|---|---|---|
| anchored to a macro, value equal (incl. macros still used in the supplement) | 489 | 344 |
| inlined generated tables (mpce_tab_baseline, friedman68, wtl, ablation, feasbudget, robust_final), numbers identical to `analysis/*.tex` | 486 | -- |
| inlined generated tables with rows added in rev. 2 (hr16: GA row; iea37: GA and exact-gradient rows); added values match `rev2_tables.tex` / `rev2_summary.json` | 157 | -- |
| traced to a printed source or json value, specific (>= 3 significant digits, few neighbouring values in the pool) | 293 | 753 |
| traced, but the pool also holds neighbouring values (not specific) | 97 | 469 |
| traced, <= 2 significant digits (weak: 30, 68, 0.05, 1.5, ...) | 1,076 | 1,750 |
| only in older repository outputs (secondary) | 10 | 22 |
| structural / untraced | 15 | 44 |

* Anchored but different: only the four per-evaluation times (item 1), in the main text and in the supplement. Text
  macros replaced by different wording: `\NXMultLostHolm` (05_setup:86, supplement:355; same content), `\NXClusters`
  ("six clusters" -> "six designed groups"), the abstract's Horns Rev feasibility "30/30 vs 19/30" and McNemar p
  (removed with the rewritten sentence), `\NFIEAGap...Proj` (sign fixed, see above), `\NFPyWake...LossPct` (item 3).
* Typed revision-2 tables in the supplement vs. the regenerated `rev2_tables.tex`: all rows identical for
  tab:S-rev-laplace, -ga, -ga-hr, -ga-iea, -spacing-rank, -spacing, -grad, -grad-calls, -site (differences only in
  captions/notes; the added note values 3.149% (LXREP, `/laplace/mean_loss_pct_all_qualified`) and "PSO qualifies in
  8 of 10 cases, SSA-VNS and RSD-VNS in 8" (`/spacing/spacing/6D/excl/n_qualified`) are correct).
  `archive_reliability_tables.tex` = `validation/reliability_tables.tex`, and a rerun of `audit_archive.py` (scratch
  copy) reproduces `reliability.csv`, `paired_outcomes.csv` and `reliability_tables.tex` byte for byte.
* Prose of the five revision-2 supplement sections against a section-scoped pool (matching subtree of the regenerated
  `rev2_summary.json` + the section's regenerated tables): 221 of 276 non-weak numbers found; the other 55 are design
  constants (2,160, 6,480, 540 runs, 3,030/15,030, 308/385/462 m), derived quantities, or values not stored in any
  output; all of these were recomputed (2c) or checked by hand.
* The 59 structural/untraced numbers are derived constants, each checked: R/K = 38.5/0.075 = 513 m; 770 m = 10D,
  21.6% / 8.8% deficit (C_T = 0.8); 1,472 kW; 936.38 kW = 14,045.74/15 (62% capacity factor); 4.6e-8 = (459.1-1)/1e10;
  0.0015625 = 0.05/2^5; design-table run counts 6,120 = 3 x 2,040, 1,080 = 3 x 360, 3,240 = 9 x 360, 2,220 =
  2,040+60+120 (row counts of the run files; total 47,830); "more than 40,000 runs"; 1/10,001 (add-one bootstrap);
  4,530/5,430/2,985/15,030 evaluations; 943 m = largest |coordinate| of the Horns Rev parallelogram
  (`hornsrev_model.site`); 1.0005 (`feasible_init.py`); 1.003 (`hornsrev_model.py` docstring); theory constants
  2.99236, 3.3475, 3.4971, 601 x 601, table tab:S-th-box (T13 of the theory check).
* Main-text tables typed by hand: tab:equiv-main (06_results) reproduces tab:X-equiv-levels and the `\NBay...Rope`
  macros exactly for all six pairs (the 1-degree row = `\NFPairMeanOne`, `\NFPairCIOne`); tab:ga-main (08_beyond)
  = tab:R2-ga; tab:iea-means (08_beyond) values = tab:R2-grad / tab:R2-ga-iea; tab:design run counts check.

## 4. Part 2c: targeted recomputations (from run records; 63 OK, 4 different = items 2, 4, 5)

* tab:S-time: all 27 spacing / Lillgrund medians with ratios and the 6 IEA37 entries reproduce exactly (spacing:
  10 split cases, 5D and 6D pooled). Main text 08_beyond:271: 1.8 (1.81) OK, 6.5--6.7 (6.48--6.68) and 3.9--4.3
  (3.88--4.29) OK, Lillgrund 2.3 -> 2.2 (item 2).
* Exact-gradient accounting (rev2_grad): FunCalls + 3 GradCalls + Unused = Budget in all 240 runs; means 1,666.6 /
  1,646.8 objective and 1,454.2 / 1,460.9 gradient evaluations (text 1,667/1,647/1,454/1,461); max(FunCalls + 19
  GradCalls) = 29,486 / 29,550 (c_g = 20: 30,952 / 31,020, so "any c_g <= 19" is exact); gradient share 72.3--73.0%.
  SDs 1,801--6,101 vs. PSO-VNS 4,587--19,268 MWh exact. Mean times per run (sup:775) exact for all four methods.
  `rev2_gradient.py check`: normwise FD error 7.6e-10 / 1.34e-9 ("at most 1.3e-9" OK); the cost ratio is machine
  dependent (2.13 / 2.29 here vs. 2.29 / 2.26 in the text; c_g = 3 either way); the JAX comparison (6e-15) could not be
  rerun (JAX not installed).
* Spacing all-run table (tab:S-allrun-spacing): all Feas., W/T/L and scores of the 9 rows exact (bootstrap CIs not
  recomputed: seed not recorded).
* tab:S-eqclus: all case-weighted means, equal-cluster means and t5 CIs exact; three cluster means off by 0.001 (item 4).
* Lillgrund: PyWake check rerun (PyWake 2.6.20): 113.40 / 116.10 GWh/yr, 2.38%, 18.65% vs. 16.72%, local-C_T 114.19,
  1-degree 0.50% (114.18 vs. 113.62), 1.65% drop; Horns Rev 0.53% / 0.23% from `pywake_check.csv`; best PSO-VNS run at
  30,030 116.16 GWh/yr (+0.05%), minimum spacing 3.00D--3.08D; McNemar p 1.5e-5 / 9.8e-4; 19 jointly feasible seeds,
  17 higher, 114.60 vs. 113.62: all exact. Area 1.081 km^2: item 5.
* IEA37 gaps: GA 3.14% / 5.20%, PSO-SLSQP 0.04 / 1.84% strict and 0.64 / 3.92% projected: exact. All-run reliability
  paragraph (08_beyond:265) = `validation/paired_outcomes.csv`; record counts 11,640 exact.

## 5. Part 3: hand-edited generated files (`SWEVO_rev2/latex_source/analysis/*.tex` vs. `analysis/*.tex`)

Only one numeric change: `mpce_numbers.tex` l. 651--655, `\NMsEvalMin` 0.23 -> 0.24, `\NMsEvalMax` 0.45 -> 0.47,
`\NMsEvalSLSQP` 0.53 -> 1.82, `\NMsEvalPSO` 0.23 -> 0.24, `\NSLSQPSlowdown` 2.3 -> 7.5. Recomputed from the 2,040
6,030-evaluation runs per method (Calls = 6,030 in every run; files fresh_grid/bgrid/vgrid, mpce_psobv, mpce_psoc,
mpce_slsqp):

| ms per evaluation | PSO | PSO-VNS | SSA-VNS | SSA | LX-SSA | DE | VNS | MS-SLSQP | SLSQP / PSO |
|---|---|---|---|---|---|---|---|---|---|
| median over runs (generator definition) | 0.2348 | 0.3407 | 0.3780 | 0.3856 | 0.3821 | 0.4468 | 0.3583 | 0.5314 | 2.26 |
| mean over runs (= Table S-cost) | 0.2444 | 0.3490 | 0.3879 | 0.3986 | 0.3956 | 0.4668 | 0.3711 | 1.8202 | 7.45 |

So the generated macros are correct for the definition in `mpce_results.py` (median), the revision-2 values are
correct means except that the ratio rounds to 7.4, not 7.5 (item 1). The other differences are wording only:
`\NXMultLostHolm` text, sentence-case captions (mpce_tab_*), "ranked last" -> "ranked below the qualified methods, by
number of feasible runs" in the ranking notes (matches `rank_rule`, which orders non-qualified methods by their
feasible runs), math-mode minus signs in mpce_tab_robust_final, caption wording in mpce_supp_inference and
mpce_supp_diag ("relabelling" -> "relabeling"), "~\cite" spacing. `archive_reliability_tables.tex` is new (identical
to `validation/reliability_tables.tex`).

## 6. Part 4: check suites on the repository data

| Suite | Result |
|---|---|
| `mpce_check_final.py --tex SWEVO_rev2/latex_source/SWEVO_manuscript.tex` | 62 PASS, 0 FAIL, 0 PENDING, 0 TAG ERRORS (0 malformed, 0 undefined ids); not referenced: C05, C09, C11, C19, C30, C39, C43 (same as for the repository macro version) |
| `mpce_check_extra.py` | 57 PASS / 0 FAIL |
| `mpce_check_diag.py` | 22 PASS / 0 FAIL / 0 PENDING |
| `mpce_check_dir.py` | 41/41 passed |
| `mpce_check_csweep.py` | 15 PASS / 0 FAIL / 0 PENDING |
| `rev3_numbers_theory_check.py --check` (rev-2 supp_theory.tex + 04_methods.tex) | 57 PASS, 0 FAIL |

These are the counts quoted in the comment at `05_setup.tex:102`. Limitation: the suites test data conditions and the
CHECK-id coverage of the comments; since revision 2 typed the numbers, they no longer test that the printed numbers
equal the macros. That is what parts 1--3 of this audit do.

## 7. How much is traced

Of the 6,005 numbers, 5,946 (99.0%) are traced to a source (anchored macro, identical generated table, specific or
non-specific pool value, or recomputed in 2c); 2,826 of them are only weakly traceable by value (<= 2 significant
digits, mostly design constants); the remaining 59 are derived constants checked by hand (section 3). Wrong or
ambiguous values: items 1--5 of section 1 (main text: the per-evaluation times and their ratio, the Lillgrund time
ratio 2.3, the PyWake wake-loss order; supplement: the same per-evaluation times via the hand-edited macros, three
tab:S-eqclus cluster means, the Lillgrund area, the Lillgrund wake-loss order).

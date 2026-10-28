# GA in the inference-sensitivity, imputation and timing tables (rev5_ga_infer_time)

Script: `analysis/rev5_ga_infer_time.py` (one process, about 15 s, no optimizer run, no existing file modified,
nothing under `SWEVO_rev2/` written). Every number below comes from `rev5_ga_infer_time.json` or `rev5_ga_infer_time.log`.

## 1. Outputs (new files in analysis/)

| File | Contents |
|---|---|
| `rev5_ga_inference_extra_tables.tex` | complete replacements of `tab:X-threshold` and `tab:X-multiplicity` (same labels), nine-method pool |
| `rev5_ga_inference_extra_macros.tex` | `\renewcommand` for the 35 `\NX...` macros of `mpce_numbers_extra.tex` whose value changes (old value in a comment) |
| `rev5_ga_imputation_tables.tex` | `tab:S-sa-imputation` with GA (main family PSO-VNS vs. 8 methods, Holm over 8), latex_source label and wording |
| `rev5_ga_time_tables.tex` | `tab:S-sa-time` (GA row added) and `tab:S-sa-time-sites` (GA HR 6,030 F and 120,030 filled), latex_source labels and wording |
| `rev5_ga_infer_time.json` / `.log` | old and new blocks, all cell changes, changed macros, LOCO ranking check |

Integration, not done because `latex_source` was not to be edited (tested on a scratch copy of `latex_source`: the
supplement compiles with no errors and no overfull boxes, and all four tables appear with the GA content):
```
\input{analysis/rev5_ga_inference_extra_macros.tex}        % after \input{analysis/rev5_ga_numbers_robust.tex} (line 45)
\CaptureSuppTables{analysis/rev5_ga_inference_extra_tables.tex}   % after line 486 (rev5_ga_tab_robust_supp.tex)
\CaptureSuppTables{analysis/rev5_ga_imputation_tables.tex}        % after line 1057 (rev5_ga_inference_tables.tex)
\CaptureSuppTables{analysis/rev5_ga_time_tables.tex}              % after line 1059 (rev5_ga_sites_tables.tex)
```
The files must also be copied to `latex_source/analysis/`.

## 2. Procedure and reproduction (step 0)

GA records: `rev2_ga_s*of10.csv`, 2,040 runs = 68 cases x seeds 1-30, 6,030 evaluations, random starts, loaded with
`mpce_inference_extra.std_cols` like all the other records. GA is placed after DE (order of `tab:friedman68-full`).
The identical case set and ideal values are checked.

Before GA was added, every existing output was regenerated with the same functions and asserted to be identical:
`mpce_summary_extra.json` threshold / threshold_summary / multiplicity (recomputed from the run records),
`mpce_supp_inference.tex` (all 15 tables), `mpce_numbers_extra.tex` (all 1,026 macros), `rev3_inference.json` B5 (recomputed)
and `rev3_inference_tables.tex`, `rev3_sites.json` time (recomputed with `block_time`) and the two time tables of
`rev3_sites_tables.tex`. After the label and wording edits of the copies, the latex_source versions of `X-threshold`,
`X-multiplicity`, `S-sa-imputation`, `S-sa-time` and `S-sa-time-sites` were also reproduced exactly.

How each table was regenerated with GA:
* X-threshold / X-multiplicity: `mpce_inference_extra.threshold_block` / `multiplicity_block`. Its module constant `MAIN8`
  (with `MAIN_PAIRS`) is set to the nine methods for the call. The table text comes from `write_tables` and the macros from `write_macros`.
* S-sa-imputation: `rev3_inference.imputation_block` on the case statistics with GA (main family of 8 pairs), and
  `rev3_inference.tables`. The ablation family is unchanged, which was asserted.
* S-sa-time: `block_time` and `time_table` keep their method lists in local variables, so they cannot add GA.
  The GA row therefore uses copies of their statements (`bm_entry`, `bud_entry`, `time_row`, `rc`). Those copies were
  first checked: they reproduce every stored benchmark value, budget value and table row exactly.
* S-sa-time-sites: `block_time` with `rev3_sites.load_hr` extended by `rev5_ga_hr_s0of1.csv`. All non-GA cells were
  asserted unchanged.

## 3. Significance verdicts

Only one verdict changes. With GA in the pool it is the one point where the threshold analysis no longer says "no change":
* **X-threshold, threshold 30.** PSO-VNS vs. GA (case-mean Wilcoxon) gives p = 0.0136 (table: 0.014) and Holm over 8
  p = 0.0272, with W/L 40/19. GA has fewer than 30 feasible runs in 4 cases, and these enter as imputed maximal
  differences in favour of PSO-VNS. So the main-text conclusion "PSO-VNS not significantly different from GA" fails at
  that threshold. At the other thresholds the test is not significant:
  * thresholds 1, 10 and 15: p = 0.0506, Holm-8 p = 0.101, 0 imputed cases;
  * thresholds 20 and 25: p = 0.0506, 2 imputed cases.

  The "Same" column of row 30 becomes **no** (yes -> \textbf{no}).
* **X-threshold, everything else.** PSO-VNS is best-ranked at every threshold. The nine-method rank order is the same
  at all thresholds: PSO-VNS, PSO, GA, SSA-VNS, VNS, SSA, MS-SLSQP, LX-SSA, DE. The rank tests find exactly {GA, PSO}
  not significantly different from PSO-VNS at every threshold. PSO-VNS vs. PSO stays not significant (Holm over 8).
  All component and control conclusions are unchanged.
* **Leave-one-cluster-out ranking (checked, table not regenerated).** PSO-VNS is best in all 6 LOCO runs, and
  {GA, PSO} are the methods not significantly different from it in each run.
* **X-multiplicity (one Holm and BH family, 28 -> 29 tests).** PSO-VNS vs. GA is added: n 68, W/L 39/20, ΔL −0.026,
  p 0.051, p_Holm 0.339, q_BH 0.061, not significant even unadjusted. The counts stay 23 significant unadjusted,
  21 after Holm and 22 after BH. The tests lost under Holm and under BH are the same as before.
* **S-sa-imputation.** The new PSO-VNS vs. GA row is "0/0 & 9 & 0.101 & 0.070 & 68 & 0.101 & 0.070 & ns/ns/ns/ns".
  These are Holm-8 values; Pratt unadjusted p is 0.035, but it is not significant after Holm. No verdict changes in
  the main family (before and after: none), and none in the ablation family.

## 4. Changed numbers (old -> new)

### tab:X-threshold
Columns **Rank GA** (after Rank PSO) and **PSO-VNS/GA** p (after PSO-VNS/PSO) were added (14 columns). The tabular is
scaled with `\resizebox{\textwidth}`, as the other wide X tables are; without it the table is 52 pt too wide.

| Thresh. | Unq. | Rank PSO-VNS | Rank PSO | Rank GA (new) | p_z ns | PSO-VNS/GA p (new) | Same |
|---|---|---|---|---|---|---|---|
| 1 | 8 -> 8 | 1.82 -> 2.18 | 1.97 -> 2.40 | 2.74 | PSO -> GA, PSO | 0.051 | yes |
| 10 | 18 -> 18 | 1.82 -> 2.18 | 1.97 -> 2.40 | 2.74 | PSO -> GA, PSO | 0.051 | yes |
| 15 (paper) | 22 -> 22 | 1.82 -> 2.18 | 1.97 -> 2.40 | 2.74 | PSO -> GA, PSO | 0.051 | yes |
| 20 | 26 -> 28 | 1.82 -> 2.18 | 1.97 -> 2.40 | 2.76 | PSO -> GA, PSO | 0.051 | yes |
| 25 | 30 -> 32 | 1.82 -> 2.18 | 1.97 -> 2.40 | 2.76 | PSO -> GA, PSO | 0.051 | yes |
| 30 | 47 -> 51 | 1.82 -> 2.16 | 1.97 -> 2.37 | 2.90 | PSO -> GA, PSO | 0.014 | yes -> **no** |

These columns are unchanged: Best (PSO-VNS), PSO-VNS/PSO p (0.296), and the four component-contrast p values.

Caption changes:
* "eight-method comparison" -> "nine-method comparison (with GA)";
* "nine-method Friedman ranking";
* PSO-VNS vs. GA added to the list of tests;
* the conclusion list now reads "only PSO and GA not significantly different in rank; PSO-VNS vs. PSO and PSO-VNS vs. GA not significant";
* a new final sentence: "With GA in the pool, the only change is at threshold 30: PSO-VNS vs. GA becomes significant
  ($p=0.014$; Holm over 8: $p=0.027$), because GA has fewer than 30 feasible runs in 4 cases, which enter the test
  with an imputed maximal difference in favour of PSO-VNS."

Definition of the conclusions with GA (documented in the script). "posthoc_only_PSOC_nonsig" becomes
"posthoc_only_PSOC_GA_nonsig", and "PSOBV-GA_ns" (unadjusted and Holm over 8) is added. This gives 18 conclusions
instead of 17, matching the main-text statement for the nine-method pool. The JSON also keeps the 8-method definition
evaluated on the nine-method ranking (`threshold_new_raw_conclusions`), which fails at every threshold because the
not-significant set is now {GA, PSO}.

### tab:X-multiplicity
* Caption: "all 28" -> "all 29" (3 times). The counts sentence becomes "23 of 29 tests are significant unadjusted, 21 after Holm and 22 after BH".
* New row: `PSO-VNS vs.\ GA & 68 & 39/20 & $-0.026$ & 0.051 & 0.339 & 0.061 & --`.
* p_Holm changes:
  * PSO-VNS vs. SSA-VNS 4.2e-10 -> 4.4e-10
  * vs. SSA 2.0e-10 -> 2.1e-10
  * vs. LX-SSA 2.0e-10 -> 2.1e-10
  * vs. DE 3.4e-10 -> 3.5e-10
  * vs. VNS 1.1e-9 -> 1.2e-9
  * LX-SSA-VNS vs. LX-SSA 2.0e-10 -> 2.1e-10
  * PSO-VNS vs. RS-VNS 5.6e-10 -> 5.8e-10
  * PSO-VNS vs. LX-SSA-VNS 3.7e-10 -> 3.9e-10
  * SSA-VNS vs. LX-SSA-VNS 7.5e-5 -> 8.1e-5
  * LX-SSA vs. SSA 3.0e-6 -> 3.2e-6
  * SSA-VNS vs. RSD-VNS 0.037 -> 0.042
  * LX-SSA-VNS vs. RSD-VNS 2.1e-5 -> 2.2e-5
  * PSO-VNS vs. RSD-VNS 7.2e-10 -> 7.6e-10
  * RSD-VNS vs. RS-VNS 9.0e-7 -> 9.5e-7
  * N≥10 both data sets 0.291 -> 0.339
  * N≥10 Data Set II 0.018 -> 0.020
  * ω 0.5 vs 0.25 0.016 -> 0.018
  * ω 0.5 vs 0.75 0.016 -> 0.018
  * ω 0.5 vs 0.9 0.188 -> 0.215
* q_BH changes:
  * PSO-VNS vs. PSO 0.331 -> 0.330
  * vs. SSA-VNS 7.6e-11 -> 7.9e-11
  * vs. SSA 5.3e-11 -> 5.5e-11
  * vs. LX-SSA 5.3e-11 -> 5.5e-11
  * vs. DE 7.6e-11 -> 7.9e-11
  * vs. VNS 1.6e-10 -> 1.7e-10
  * vs. MS-SLSQP 5.3e-11 -> 5.5e-11
  * LX-SSA-VNS vs. LX-SSA 5.3e-11 -> 5.5e-11
  * PSO-VNS vs. RS-VNS 9.3e-11 -> 9.6e-11
  * SSA-VNS vs. RS-VNS 7.8e-4 -> 8.1e-4
  * PSO-VNS vs. LX-SSA-VNS 7.6e-11 -> 7.9e-11
  * LX-SSA vs. SSA 4.0e-7 -> 4.2e-7
  * LX-SSA-VNS vs. RSD-VNS 2.7e-6 -> 2.8e-6
  * PSO-VNS vs. RSD-VNS 1.1e-10 -> 1.2e-10
  * RSD-VNS vs. RS-VNS 1.2e-7 -> 1.3e-7
  * N≥10 both 0.059 -> 0.061
  * N≥10 DS I 0.347 -> 0.346
  * ω 0.5 vs PSO 0.440 -> 0.439
  * ω 0.75 vs PSO 8.0e-4 -> 8.3e-4
  * ω 0.5 vs 0.9 0.034 -> 0.035
  * ω 0.9 vs 0.75 0.177 -> 0.176
* Every "Holm sig." entry is unchanged.

### tab:S-sa-imputation
* Family title: "Holm over 7" -> "Holm over 8".
* New row: PSO-VNS vs. GA (above).
* The p values of the other 6 PSO-VNS rows change. Order of the values: imputed Wilcoxon / imputed Pratt / dropped Wilcoxon / dropped Pratt. All verdicts stay A/A/A/A.
  * vs. SSA-VNS: 5.7e-11 -> 7.6e-11, 1.0e-11 -> 1.3e-11, 8.1e-11 -> 9.8e-11, 1.8e-11 -> 2.1e-11
  * vs. SSA and vs. LX-SSA: 4.5e-11 -> 5.3e-11, 8.3e-12 -> 9.5e-12, 8.1e-11 -> 9.8e-11, 1.8e-11 -> 2.1e-11
  * vs. DE: 5.6e-11 -> 7.0e-11, 8.3e-12 -> 9.5e-12, 3.2e-12 -> 3.6e-12, 4.1e-9 -> 6.1e-9
  * vs. VNS: 1.3e-10 -> 1.9e-10, 4.1e-11 -> 6.1e-11, 1.3e-10 -> 1.9e-10, 6.1e-11 -> 8.2e-11
  * vs. MS-SLSQP: 4.2e-11 -> 4.8e-11, 1.8e-11 -> 2.4e-11, 3.6e-11 -> 4.2e-11, 2.4e-11 -> 3.0e-11
* PSO-VNS vs. PSO and the component rows are unchanged.

### tab:S-sa-time
* New row after DE: `GA$^{b}$ & 0.32 & 0.33 & 1.65 & 2.23 & 3.01 (0.87) & 5.80 (0.73) & 10.8 (0.49) & 49.5 (0.87)`.
  * Columns: ms/eval median, ms/eval mean, s per run at N=5/10/15, then the six largest cases (feasible init., 30,030, 120,030).
  * The benchmark columns come from rev2_ga.
  * The last three cells come from the rev5_ga batch (`rev5_ga_feas`, `rev5_ga_b30k`, `rev5_ga_b120k`, 180 runs each).
* New note "$^{b}$": "GA: 68-case runs of a separate run batch (rev2\_ga); its six-largest-case runs with feasible
  initialization, 30,030 and 120,030 evaluations (rev5\_ga) were made later on a machine shared with another job, so
  these three times are inflated and indicative only."
* Worker-hours sentence: "... 12.4 for the 11,640 records of the additional experiments" -> "... 12.4 for the 11,640
  records of the additional experiments and 4.9 for the 1,360 later GA runs (rev4\_galg, rev5\_ga)".
  * The 11,640 / 12.4 are kept separate on purpose. Supplement line 1197 lists these 1,360 GA runs separately as not
    covered by the audit, and line 1221 quotes "12 for the 11,640 runs".
  * If the caller prefers the merged figure: 13,000 records, 17.3 h.
* Every existing cell is unchanged.

### tab:S-sa-time-sites
* GA, HR 6,030 F: n/r -> 6.73 (0.61)$^{a}$.
* GA, HR 120,030: n/r -> 97.6 (0.80)$^{a}$.
* Batch note: "rev2\_gahr (GA);" -> "rev2\_gahr (GA, random initialization, 6,030 and 30,030) and rev5\_ga\_hr (GA,
  6,030 F and 120,030; $^{a}$run on a machine shared with another job, so these two times are inflated and indicative only);".
* All other cells are unchanged.

### Changed \NX macros (rev5_ga_inference_extra_macros.tex)

Threshold macros:

| Macro | Old -> new |
|---|---|
| NXThrRankPSOVNSMin | 1.82 -> 2.16 |
| NXThrRankPSOVNSMax | 1.82 -> 2.18 |
| NXThrRankPSOMin | 1.97 -> 2.37 |
| NXThrRankPSOMax | 1.97 -> 2.40 |
| NXThrNConcl | seventeen -> eighteen |
| NXThrChanged | none -> "PSOBV-GA ns" (raw key, not for prose) |
| NXThrChangedAt | none -> 30 |
| NXThrImputedMax | 47 -> 51 |

Multiplicity macros:

| Macro | Old -> new |
|---|---|
| NXMultN | 28 -> 29 |
| NXMultNMain | 7 -> 8 |
| NXMultNWord | twenty-eight -> twenty-nine |
| NXMultMaxSigHolmP | 0.038 -> 0.042 |

Holm p and BH q macros:

| Macro | Old -> new |
|---|---|
| NXHolmSSAVNSvsRSVNSP | 0.0058 -> 0.0063 |
| NXBHSSAVNSvsRSVNSQ | 7.8e-4 -> 8.1e-4 |
| NXHolmSSAVNSvsLXSSAVNSP | 7.6e-5 -> 8.1e-5 |
| NXHolmPSOVNSvsRSVNSP | 5.6e-10 -> 5.9e-10 |
| NXBHPSOVNSvsRSVNSQ | 9.3e-11 -> 9.7e-11 |
| NXHolmPSOVNSvsPSOdsIILargeP | 0.018 -> 0.020 |
| NXBHPSOVNSvsPSOdsIILargeQ | 0.0028 -> 0.0029 |
| NXHolmPSOVNSvsPSOLargeP | 0.29 -> 0.34 |
| NXBHPSOVNSvsPSOLargeQ | 0.059 -> 0.061 |
| NXHolmSplitFiftyVsSeventyFiveP | 0.017 -> 0.018 |
| NXBHSplitFiftyVsSeventyFiveQ | 0.0022 -> 0.0023 |
| NXHolmSplitFiftyVsTwentyFiveP | 0.017 -> 0.018 |
| NXBHSplitFiftyVsTwentyFiveQ | 0.0022 -> 0.0023 |
| NXHolmSplitSeventyFiveVsHundredP | 0.0059 -> 0.0064 |
| NXBHSplitSeventyFiveVsHundredQ | 8.1e-4 -> 8.4e-4 |
| NXHolmSSAVNSvsRSDVNSP | 0.038 -> 0.042 |
| NXBHSSAVNSvsRSDVNSQ | 0.0062 -> 0.0065 |
| NXHolmLXSSAVNSvsRSDVNSP | 2.1e-5 -> 2.2e-5 |
| NXBHLXSSAVNSvsRSDVNSQ | 2.8e-6 -> 2.9e-6 |
| NXHolmPSOVNSvsRSDVNSP | 7.3e-10 -> 7.6e-10 |
| NXHolmRSDVNSvsRSVNSP | 9.0e-7 -> 9.5e-7 |
| NXHolmSplitFiftyVsNinetyP | 0.19 -> 0.21 |
| NXBHSplitFiftyVsNinetyQ | 0.035 -> 0.036 |

These macros are unchanged: NXThrList, NXThrBest, NXThrOrderSame (yes), NXThrImputedMin, NXMultSigRaw (23),
NXMultSigHolm (21), NXMultSigBH (22), NXMultLostHolm and NXMultLostBH.

Macros the prose uses:
* Supplement: NXThrNConcl, NXMultN, NXHolmPSOVNSvsPSOdsIILargeP, NXHolmSSAVNSvsRSDVNSP and NXHolmSplitFiftyVsNinetyP.
* Main text: NXMultN, in a comment only (05_setup.tex:155).

## 5. Sentences that must change

**Supplement line 467 (qualification threshold).** This sentence becomes false as written: with the new macros it
would say that all eighteen conclusions are unchanged.

Current text:
> with each of the six thresholds \NXThrList{} feasible runs of 30 (the paper uses 15), PSO-VNS is best-ranked, the
> complete rank order of the eight methods (pool without GA) is the same, and all \NXThrNConcl{} conclusions listed in
> the tables are unchanged.

Replacement:
> with each of the six thresholds \NXThrList{} feasible runs of 30 (the paper uses 15), PSO-VNS is best-ranked, the
> complete rank order of the nine methods is the same, and PSO and GA are the only methods not significantly different
> from PSO-VNS in rank. Of the \NXThrNConcl{} conclusions listed in the tables, all but one hold at every threshold: when
> all 30 runs must be feasible, PSO-VNS is significantly better than GA in the case-mean test ($p=0.014$; Holm over 8:
> $p=0.027$), because the four cases in which GA has fewer than 30 feasible runs enter as maximal differences in favor of
> PSO-VNS; at the other thresholds the test is not significant ($p=0.051$).

The `% CHECK-EXTRA [X02]-[X09], [X31]` comment refers to mpce_check_extra.py checks on the eight-method JSON. Those
checks were not re-run for GA.

**Supplement line 472 (multiplicity).** The values are carried by macros: 29, 21 of 23, p_Holm ≤ 0.020 / ≤ 0.042 /
= 0.21. No verdict changes. Only the parenthesis must change:
"(main comparison in the pool without GA; not over the other tests of the paper)" ->
"(main comparison: PSO-VNS vs.\ each of the other eight methods, GA included; not over the other tests of the paper)".
Optional addition: "PSO-VNS vs.\ GA is not significant even unadjusted ($p=0.051$)."

**Supplement line 488, `\SuppTableNote{tab:X-threshold}`.** Delete it, or replace it with:
"Note: ranks and tests of the main comparison in the nine-method pool with GA, as in Table~\ref{tab:friedman68-full}."

**Supplement line 494, `\SuppTableNote{tab:X-multiplicity}`.** Delete it, or replace it with:
"Note: the main-comparison rows are the eight comparisons of PSO-VNS in the nine-method pool (GA included)."

**Supplement line 1106 (imputation).** Change
"changes no verdict of the main comparison (pool without GA) or of the component contrasts with imputed cases" to
"changes no verdict of the main comparison (PSO-VNS vs.\ each of the other eight methods, GA included; Holm over 8) or
of the component contrasts with imputed cases". The statement stays true: no verdict changes.

**Supplement line 1107, `\SuppTableNote{tab:S-sa-imputation}`.** Delete it, or replace it with:
"Note: the main-comparison rows are the eight comparisons of PSO-VNS in the nine-method pool (Holm over 8)."

**Supplement line 1132 (time).** No quoted number changes:
* The metaheuristic range 0.23–0.45 ms per evaluation (median) still holds; GA is 0.32 ms (median) and 0.33 ms (mean).
* MS-SLSQP 0.53 ms, 2.3x PSO, 101 worker-hours and 36,190 runs are unaffected.

Suggested addition: "GA, from a separate run batch, needs 0.32~ms per evaluation; its runs of the later batch (six
largest cases with feasible starts and at 30,030 and 120,030 evaluations; Horns Rev~1 with feasible starts and at
120,030) were made on a machine shared with another job, so their times are indicative only."

Line 563 (`tab:cost`) has the same 0.23–0.45 ms sentence and also stays true. `tab:cost` itself is still eight-method.

**Main text, optA/sw/05_setup.tex:154.** This sentence is not in the supplement, but it quotes these tables and states
the verdict that now changes.

Current text:
> Six qualification thresholds (1--30 feasible runs), dropping the cases with an imputed maximal difference and keeping
> zero differences change no verdict; ... one Holm family over all 28 case-mean tests makes only ... nonsignificant

Replacement:
> Six qualification thresholds (1--30 feasible runs), dropping the cases with an imputed maximal difference and keeping
> zero differences change no verdict, except that PSO-VNS becomes significantly better than GA when all 30 runs must be
> feasible (four cases with fewer than 30 feasible GA runs enter as maximal differences); ... one Holm family over all
> 29 case-mean tests makes only [unchanged] nonsignificant

The CHECK-EXTRA comment on line 155 ("six thresholds (NXThrList) change no verdict") is no longer true with GA.
`SWEVO_manuscript_full.tex` is the flattened copy and must be regenerated.

**Caption of tab:X-loco (generated by mpce_inference_extra.py; not regenerated here).** Its sentence "The best-ranked
method of the eight-method Friedman ranking is PSO-VNS in every LOCO run, and PSO is the only method not significantly
different from it in rank" becomes, with GA (checked): "... of the nine-method Friedman ranking is PSO-VNS in every LOCO
run, and PSO and GA are the only methods not significantly different from it in rank."

## 6. Caveats and remaining eight-method tables

* The X-threshold conclusion set had to be adapted to the nine-method pool, because the old rule
  "only PSO not significantly different" fails for a definitional reason. The adaptation and the two added columns
  are my choices; the raw evaluation is kept in the JSON.
* The X-threshold verdict change at threshold 30 depends on the imputation rule: 4 imputed cases.
* These tables still have eight-method main-comparison rows (no PSO-VNS vs. GA), which conflicts with the co-author rule.
  They were outside this task and were not regenerated:
  * `tab:X-cluster`, `tab:X-loco` and `tab:X-hl`: rows of ALL_PAIRS;
  * `tab:X-equiv-levels` and `tab:X-bayes`: EQ_PAIRS, the main-comparison rows without GA;
  * `tab:cost`.
* The GA timings of the rev5_ga batch, and the rev2_ga timings, come from batches different from those of the other
  methods. Ratios across batches are indicative, as the tables already state.

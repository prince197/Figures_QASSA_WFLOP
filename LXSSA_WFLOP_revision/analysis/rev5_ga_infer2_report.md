# GA in the cluster, LOCO, Hodges-Lehmann, equivalence-level, Bayesian and cost tables (rev5_ga_infer2)

Script: `analysis/rev5_ga_infer2.py`. One process, about 45 s. No optimizer is run, no existing file is modified, and
nothing under `SWEVO_rev2/` is written. Every number below comes from `rev5_ga_infer2.json` or `rev5_ga_infer2.log`.

## 1. Outputs (new files in analysis/)

| File | Contents |
|---|---|
| `rev5_ga_inference_extra2_tables.tex` | complete replacements of `tab:X-cluster`, `tab:X-loco`, `tab:X-equiv-levels`, `tab:X-bayes` and `tab:X-hl`, using the same labels and the latex_source wording, with the pair PSO-VNS vs. GA added |
| `rev5_ga_inference_extra2_macros.tex` | `\renewcommand` for the 4 `\NX...` macros whose value changes (old value in a comment), plus `\newcommand` for the 34 new PSO-VNS-vs-GA macros |
| `rev5_ga_cost_tables.tex` | `tab:cost` with a GA row |
| `rev5_ga_infer2.json` / `.log` | GA blocks, LOCO rankings, cell and text diffs, macro changes, cost values |

Integration was tested on a scratch copy of `latex_source` but not applied, because `latex_source` was not to be edited.
In that test the supplement compiled with no errors and no overfull boxes, the same as the unmodified copy. The GA rows
and the new captions appear in the PDF. The lines to add are:
```
\input{analysis/rev5_ga_inference_extra2_macros.tex}               % after line 47 (rev5_ga_inference_extra_macros.tex)
\CaptureSuppTables{analysis/rev5_ga_cost_tables.tex}               % after line 112 (mpce_supplementary.tex)
\CaptureSuppTables{analysis/rev5_ga_inference_extra2_tables.tex}   % after line 487
```
The three files must also be copied to `latex_source/analysis/`. A later `\CaptureSuppTables` overrides the stored body of
the same label (global `\let`), so the old tables need not be removed.

## 2. Procedure and reproduction (step 0)

**Data.** The GA records are `rev2_ga_s*of10.csv`: 2,040 runs, i.e. 68 cases x seeds 1-30, at 6,030 evaluations with
random starts. They are loaded with `mpce_inference_extra.std_cols`, like all the other records. The script checks that
the case set and the ideal values are identical to those of the other methods. GA is placed after DE, the order of
`tab:friedman68-full`:
* `ALL_PAIRS`: 19 -> 20 pairs, with PSO-VNS vs. GA after PSO-VNS vs. DE.
* `EQ_PAIRS`: 18 -> 19 pairs, with PSO-VNS vs. GA in the main-comparison block after PSO-VNS vs. DE.
* The Friedman ranking of the LOCO runs uses the nine-method pool, with module constant `MAIN8` set to the nine methods
  for that call.

**Step 0: every existing output was regenerated and asserted identical** (exact float equality of the JSON,
character equality of the text) before GA was added:
* `mpce_summary_extra.json`: the `case_level`, `cluster`, `loco`, `equivalence_levels` and `hodges_lehmann` blocks,
  recomputed from the run records with `cluster_block`, `loco_block`, `equivalence_levels_block` and `hl_block`;
* `mpce_supp_inference.tex` (all 15 tables), from `write_tables`;
* `mpce_numbers_extra.tex` (all 1,026 macros), from `write_macros`;
* the latex_source copies of the five tables, after its wording edits:
  * X-equiv-levels: "sensitivity to exchangeable-case assumptions" and "dependence within the six designed benchmark groups";
  * X-bayes: `~\cite`;
* `tab:cost`: the statements of `mpce_results.main` were copied into `cost_rows` and reproduce the table in
  `analysis/` and in `latex_source/analysis/mpce_supplementary.tex`, as well as `mpce_summary.json` `cost.ms_per_call`
  and `cost_per_eval.median_ms`.

**Holm families.** None of these six tables has a multiplicity adjustment. The cluster, LOCO and Hodges-Lehmann p values
are unadjusted. The TOST, bootstrap and Bayesian results of each pair use their own fixed seeds. There is therefore no
Holm family to grow. The script asserts that every value of the existing pairs is unchanged with GA in the pool, so only
the GA rows are new.

## 3. Verdicts and findings with GA

1. **LOCO: verdict change for PSO-VNS vs. GA.** On all 68 cases the test is not significant: p = 0.051, W/L 39/20,
   mean ΔL = −0.026 pp. Leaving out one cluster makes it significant in favor of PSO-VNS (unadjusted) in 3 of the 6 runs:
   * without I/750: p = 0.0115;
   * without I/1000: p = 0.0204;
   * without II/750: p = 0.0368.

   In the other three runs it stays not significant: without I/500 p = 0.103, without II/500 p = 0.227, without II/1000
   p = 0.352. The direction never flips. As a result, `\NXLocoChanged` becomes "PSO-VNS vs.\ GA and SSA-VNS vs.\ RSD-VNS".
2. **LOCO ranking.** In all 6 LOCO runs PSO-VNS is best-ranked, and {GA, PSO} are exactly the methods not significantly
   different from it in rank (Holm average-rank tests). This confirms the caption wording "PSO and GA are the only
   methods not significantly different".
3. **Cluster level (X-cluster).** PSO-VNS vs. GA is not significant on any test. PSO-VNS has the lower mean in 4 of 6
   clusters. The p values are: sign 0.688, exact Wilcoxon 0.438, flip 0.375, CR1 0.370, case 0.051. No count macro
   changes, because the pair is not significant on the case means.
4. **Equivalence (X-equiv-levels): PSO-VNS vs. GA is equivalent at the seed level.**

   | Level | 90% CI (pp) | m_min (pp) | p_TOST | Equivalent |
   |---|---|---|---|---|
   | Seed | [−0.043, −0.009] | 0.043 | 0.011 | yes; also under joint seed resampling, [−0.043, −0.009] |
   | Case | [−0.062, 0.006] | 0.063 (rounded up; raw 0.0623) | 0.130 | no |
   | Cluster | CR2 [−0.078, 0.027] | 0.079 | wild 0.183, wild m_min 0.081 | no |

   The number of pairs equivalent at the seed level (`\NXEqSeedNHolds`) goes from 3 to 4. Still no pair is equivalent at
   the cluster level.
5. **Bayesian (X-bayes).** P_A < 0.01, P_rope > 0.99 (0.994) and P_B < 0.01. The posterior mean θ values are
   0.26 [0.15, 0.40] / 0.57 [0.45, 0.69] / 0.16 [0.07, 0.28].
6. **Hodges-Lehmann (X-hl).** HL = −0.012 pp, with 95% CI [−0.026, −0.000]. The upper limit is −0.00018, so the
   interval just excludes 0, while the case-mean Wilcoxon test gives p = 0.051. The two disagree because the HL CI keeps
   the 9 zero differences and the Wilcoxon test drops them, as the caption already explains. The mean is −0.026, with
   95% CI [−0.071, 0.012].
7. **Cost (tab:cost).** GA needs 0.326 ms per evaluation (mean over 2,040 runs) and 1.4–3.0 s per run. Its median is
   0.316 ms per evaluation. The metaheuristic median range stays 0.23–0.45 ms, and the mean range is 0.24–0.47 ms. The
   GA times come from a separate batch (rev2_ga), so the row is marked `GA$^{a}$` and the caption has the note.

## 4. Changed numbers (old -> new)

No existing cell of any of the six tables changes; each table gains one row.

| Table | New row |
|---|---|
| X-cluster | `PSO-VNS vs.\ GA & 68 & $-0.026$ & $-$0.027 & +0.046 & +0.022 & $-$0.148 & $-$0.007 & $-$0.065 & 4/2 & 0.688 & 0.438 & 0.375 & 0.370 & 0.051` (cluster columns I/500, I/750, I/1000, II/500, II/750, II/1000) |
| X-loco | `PSO-VNS vs.\ GA & $-0.026$ & [$-$0.071, 0.012] & [$-$0.077, 0.017] & [$-$0.093, 0.041] & ns & 0.012--0.352 & 0 & I/750, I/1000, II/750` |
| X-equiv-levels | `PSO-VNS vs.\ GA & 68 & $-0.026$ & [$-$0.043, $-$0.009] & 0.043 & yes & [$-$0.062, 0.006] & 0.063 & 0.130 & no & [$-$0.078, 0.027] & 0.079 & 0.183 & 0.081 & no` |
| X-bayes | `PSO-VNS vs.\ GA & 68 & $<$0.01 & $>$0.99 & $<$0.01 & 0.26 [0.15, 0.40] & 0.57 [0.45, 0.69] & 0.16 [0.07, 0.28]` |
| X-hl | `PSO-VNS vs.\ GA & 68 & $-0.026$ & [$-$0.071, 0.012] & $-0.012$ & [$-$0.026, $-$0.000] & 39/20 & 0.051` (not bold) |
| cost | `GA$^{a}$ & 0.326 & 1.4--3.0` (after DE) |

Caption changes:
* **X-loco:** "The best-ranked method of the eight-method Friedman ranking is PSO-VNS in every LOCO run, and PSO is the
  only method not significantly different from it in rank." becomes "The best-ranked method of the nine-method Friedman
  ranking (with GA) is PSO-VNS in every LOCO run, and PSO and GA are the only methods not significantly different from
  it in rank." The generator's own clause tests "only PSO", which is false with GA, so the script sets the GA clause
  after asserting {GA, PSO} in all six runs.
* **X-equiv-levels:** "of the pairs of Table~\ref{tab:equivalence} at three levels" becomes "of the pairs of
  Table~\ref{tab:equivalence} and of PSO-VNS vs.\ GA at three levels".
* **X-bayes:** "as Table~\ref{tab:equivalence})." becomes "as Table~\ref{tab:equivalence}; with PSO-VNS vs.\ GA added)."
* **cost:** the caption gains "$^{a}$GA: runs of a separate batch (rev2\_ga), so its times are indicative only."

Changed macros:

| Macro | Old -> new |
|---|---|
| NXLocoNPairs | nineteen -> twenty |
| NXLocoNChanged | one -> two |
| NXLocoChanged | SSA-VNS vs.\ RSD-VNS -> PSO-VNS vs.\ GA and SSA-VNS vs.\ RSD-VNS |
| NXEqSeedNHolds | 3 -> 4 (not used in the text) |

These macros are unchanged:
* NXClNCaseSig, NXClNCaseSigUnanimous, NXClCaseSigNotUnanimous, NXClNSigWil, NXClNUnanimous, NXClNBootExcl, NXClNCRSig;
* NXLocoFlipped, NXLocoBestAlways (yes);
* NXEqCaseNHolds, NXEqClustNHolds.

New macros (34) cover PSO-VNS vs. GA:
* `\NXClPSOVNSvsGA{Fav,Opp,SignP,WilP,FlipP,CRP,CI,CRCI}`;
* `\NXLocoPSOVNSvsGA{PMin,PMax,Changes,Flips}` = 0.011 / 0.36 / 3 / 0;
* `\NXEq{Mean,Seed*,Case*,Clust*}PSOVNSvsGA`;
* `\NXBayMeanPSOVNSvsGA{Left,Rope,Right}`.

Values are in the macros file.

## 5. Sentences that must change

**Supplement line 470 (dependence between cases).** With the new macros the sentence becomes ungrammatical ("only that
of PSO-VNS vs. GA and SSA-VNS vs. RSD-VNS"), and it hides the direction of the change.

Current:
> ... and leaving out any one cluster changes, of the \NXLocoNPairs{} verdicts, only that of \NXLocoChanged{}; PSO-VNS
> is best-ranked in every leave-one-cluster-out run.

Replacement:
> ... and leaving out any one cluster changes, of the \NXLocoNPairs{} verdicts, only those of \NXLocoChanged{}: PSO-VNS
> vs.\ GA, not significant on all 68 cases ($p=0.051$), becomes significant in favor of PSO-VNS (unadjusted) when one of
> the clusters I/750, I/1000 or II/750~m is left out ($p=0.012$--$0.037$); PSO-VNS is best-ranked in every
> leave-one-cluster-out run, and PSO and GA are the only methods not significantly different from it in rank.

The rest of the paragraph stays true: "equivalence at the cluster level is not shown for any pair".

**Supplement lines 451–454 (practical equivalence).** No statement becomes false. "No pair of Table~\ref{tab:equivalence}
is equivalent at the cluster level" still holds, and it also holds for GA. Suggested addition after line 454:
> PSO-VNS vs.\ GA (Table~\ref{tab:X-equiv-levels}) is equivalent at the seed level (90\% CI \NXEqSeedPSOVNSvsGACI~pp,
> $m_{\min}=\NXEqSeedPSOVNSvsGAMin$~pp) but not at the case level ($m_{\min}=\NXEqCasePSOVNSvsGAMin$~pp) or the cluster
> level ($m_{\min}=\NXEqClustPSOVNSvsGAMin$~pp); practical equivalence is the most probable outcome of the Bayesian test
> ($P_{\rm rope}>0.99$).

**Supplement line 775 (S-xp-ga) and main text optA/sw/06_type1.tex:74.** "Equivalence of GA and PSO-VNS is not shown"
is a case-level statement, and it stays true at that level. 06_type1 itself defines the seed level as "the level at which
verdicts are decided" (tab:equiv-main caption), so the sentence is incomplete. Suggested wording for 06_type1:74:
> ... whereas equivalence of GA and PSO-VNS is not shown at the case level (90\% intervals of GA minus PSO and GA minus
> PSO-VNS: [$-$0.022, 0.041] and [$-$0.006, 0.062]~pp; Table~\ref{S-tab:S-xp-ga}); on the fixed benchmark (seed level)
> PSO-VNS and GA are equivalent (90\% CI of PSO-VNS minus GA [$-$0.043, $-$0.009]~pp), not at the cluster level
> (Table~\ref{S-tab:X-equiv-levels}).

Line 775 needs the same addition. Note on m_min for PSO-VNS vs. GA at the case level: line 775 gives 0.062, while
X-equiv-levels gives 0.063. The raw value is 0.0623; X-equiv-levels and tab:equivalence round m_min up.

**Main-text table caption, optA/sw/06_type1.tex:128 (tab:equiv-main).**
"All 18 pairs, TOST $p$ values and minimal margins: Tables~\ref{S-tab:X-equiv-levels} and~\ref{S-tab:equivalence}." becomes
"All pairs (18 in Table~\ref{S-tab:equivalence}; with PSO-VNS vs.\ GA, 19 in Table~\ref{S-tab:X-equiv-levels}), TOST $p$ values and minimal margins: ..."

**Comment optA/sw/06_type1.tex:149.** "CHECK-EXTRA [X55] ... main-comparison pairs other than PSO: whole 90% CI beyond -m
at all three levels" is false for GA. Suggested comment: "... other than PSO and GA; PSO-VNS vs GA equivalent at the
seed level only, rev5_ga_infer2_report.md".

**Main text, optA/sw/05_setup.tex:154.** "Cluster analyses over the six (data set, radius) clusters qualify only the
significant contrasts of the salp-swarm hybrids with the controls" stays literally true: no significant result is
weakened beyond those. But the LOCO analysis now changes a verdict in the other direction. Suggested:
> ... qualify only the significant contrasts of the salp-swarm hybrids with the controls, and leaving out one cluster
> makes PSO-VNS vs.\ GA significant (unadjusted) in three of six runs; ...

The `[X14]` comment on line 155 uses `\NXLocoChanged`, so it updates automatically.

**tab:cost, supplement lines 562–563.**
* Line 562 needs no change: the metaheuristic median range is still 0.23–0.45 ms, and GA's median is 0.32 ms.
* Note 563, "means over the 2,040 runs of each method", still holds for GA. Suggested addition: "GA: runs of a separate
  batch (rev2\_ga); its times are indicative only."
* Main text optA/sw/05_setup.tex `tab:runtime` (0.23–0.45 ms) is unchanged.

**Supplement line 482 (Hodges–Lehmann).** No change is needed: the GA row is not a discussed (bold) pair.

**Lines that refer to tab:equivalence's 18 pairs.** These stay correct, because `tab:equivalence` itself was not given GA:
* supplement line 461 (comment);
* optA/sw/07_ablation.tex:139 and :141 ("for all 18 pairs");
* the tab:S-sa-tostaudit sentence at supplement line 1091 ("all 18 pairs at all five levels").

If `tab:equivalence` (mpce_results.equivalence_block) were also given GA, its row would carry the case-level values
above. `equivalence_levels_block` reproduces that function exactly (check X38):
* ΔL −0.026, 95% CI [−0.071, 0.012], 90% CI [−0.062, 0.006];
* p_TOST 0.130, m_min 0.063, not equivalent;
* P_A / P_rope / P_B = <0.01 / >0.99 / <0.01.

## 6. Caveats

* The "PSO-VNS vs. GA" LOCO verdict change is unadjusted, as is every test of tab:X-loco. On all 68 cases the pair is
  not significant even unadjusted (p = 0.051).
* The X-equiv-levels and X-bayes caption edits, and the `GA$^{a}$` marker with its note in tab:cost, are my wording.
* Leftover from the previous task, not in scope here: optA/sw/06_type1.tex still says "$p_{\rm Holm}\le0.018$"
  (line 118, table row) and "Holm over all 28 case-mean tests" (line 122, the note below the table). With GA in the multiplicity family
  (rev5_ga_infer_time_report.md), these values are ≤0.020 and 29.

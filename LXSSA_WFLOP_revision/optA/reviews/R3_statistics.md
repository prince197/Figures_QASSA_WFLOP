# R3 — Statistics review (MPCE_PSO_VNS, optA sources)

Reviewed: optA/01_front, 05_setup, 06_results, 07_ablation, 08_beyond, 09_robust, 10_limits_concl; analysis/mpce_results.py (holm, wil, rank_rule, case_stats, friedman_block, paired_vs, case_mean_wilcoxon, baseline_setting, split_section, HR block); analysis/mpce_summary.json; analysis/mpce_numbers.tex; generated tables (friedman68, wtl, split, baseline). The PDFs could not be rendered in this container (no poppler), so the review is based on the sources and generated tables.

## Overall assessment

The design is careful and better than most in this field: equal budgets that count gradient calls, 30 seed-paired runs from common initial populations, Holm correction, a case-mean Wilcoxon added because of the Benavoli critique of mean-rank post hoc tests, run-level W/T/L with rank-biserial effect sizes, and a feasibility-aware ranking rule. The code does what the text describes (Holm step-down is correct, Iman–Davenport is correct, post hoc z = Δrank/sqrt(k(k+1)/6n) is correct). Every number I checked matches mpce_summary.json and mpce_numbers.tex. Checked: χ²_F = 319.3, p = 4.6e-65, F_ID = 136.5; the ablation χ² = 326.8; p_z(PSO) = 0.713; case-mean p(PSO) = 0.296; W/T/L 14/50/4, 49/19/0 and so on; median r_rb; the baseline ranks; the old-PSO numbers 5.12 / 86.0% / 56 W; HR p_Holm max = 6.7e-4 (SSA-VNS).

The main problems are in inference, not arithmetic:
1. **Tests are applied unevenly.** The case-mean test is the headline test for PSO-VNS vs PSO, where it is non-significant. It is ignored for SSA-VNS vs RS-VNS, where it *is* significant: p = 4.5e-4, which I recomputed from mpce_ablation_tests.csv.
2. **Non-significance is read as equivalence or non-inferiority.** Examples: "tied", "no better than random sampling", "at least as good".
3. **The N ≥ 10 subgroup claim is post hoc** and only borderline.
4. **"VNS helps every swarm" is contradicted for PSO** by the paper's own tests.
5. **Effect sizes have no uncertainty intervals**, and the 68 cases are treated as independent even though they form 6 nested clusters.

None of these needs new experiments. All the fixes are re-analyses of existing CSVs or rewording.

## Issues

1. **[major] 07_ablation.tex** — "SSA-VNS and RS-VNS are statistically tied in \NAblSSAVNSvsRSVNST{} of the 68 cases … As start generators for VNS, SSA and LX-SSA are no better than random sampling". The same claim appears in the abstract ("salp swarm phases … add nothing") and the Conclusion ("statistically tied with random sampling").
   - **Problem:** 66 non-significant per-case run-level tests (30 runs, Holm within case) are not evidence of a tie. The paper's own case-level test, the Wilcoxon on per-case mean losses used as the primary case-level test in Sec. 5, gives SSA-VNS lower loss than RS-VNS in 42 of 60 non-tied cases (mean −0.068 pp, p = 4.5×10⁻⁴ raw, still < 0.01 after Holm over the 8 contrasts), computed from DLoss in mpce_ablation_tests.csv. For LX-SSA-VNS vs RS-VNS the claim holds: 27 vs 33 cases, p = 0.72.
   - **Fix:** Add a case-mean Wilcoxon column (p_W, Holm over the 8 contrasts) to Table ablation. Replace the text with: "The salp-swarm phases add little beyond random sampling: SSA-VNS is significantly better than RS-VNS in only 2 of 68 cases at the run level, and its advantage over the 68 case means is small (0.07 pp) although significant (Wilcoxon, p = 4.5×10⁻⁴). LX-SSA-VNS and RS-VNS do not differ (p = 0.72). Both are clearly worse than PSO-VNS (0.33 and 0.39 pp)." Abstract: "salp swarm phases add little beyond random sampling".

2. **[major] 01_front.tex / 07_ablation.tex** — "a variable neighborhood search (VNS) phase helps every swarm" / "Adding VNS helps every swarm: each hybrid ranks ahead of the same swarm run alone".
   - **Problem:** For PSO the evidence is 8/58/2 at the run level, post hoc p_Holm = 0.58 (ablation.friedman.p_holm_vs_focus.PSOC) and case-mean p = 0.30. Ranking ahead by 0.24 average rank is not "helps". The claim is significant only for SSA and LX-SSA.
   - **Fix:** "A VNS phase significantly improves SSA and LX-SSA. For a convergent PSO, switching to VNS for the second half of the budget gives a small gain that is not significant over all cases (p = 0.30)."

3. **[major] 06_results.tex** — "at least as good as a correctly configured PSO and better for $N\ge10$". The same claim appears in the Conclusion ("is better than PSO for N≥10") and the Discussion ("appears mainly for N≥10").
   - **Problem, "at least as good":** This is a non-inferiority claim, and no non-inferiority test was run. PSO-VNS has 4 significant run-level losses, PSO has the lower mean loss in 26 of 68 cases, and the median ΔL is −0.0002 pp.
   - **Problem, "better for N≥10":** This is a post hoc subgroup with a threshold that was not pre-specified. The subgroup is confounded with radius (6 of the 10 cases per data set are r = 1000 m). Its evidence is 9 W vs 2 L (sign test p = 0.065) or a case-mean Wilcoxon on 20 cases with p = 0.048 (my computation, unadjusted).
   - **Problem, cluster means:** The (DS, r) cluster means of ΔL(PSO − PSO-VNS) are −0.064, −0.043, +0.008, +0.063, +0.031 and +0.090 pp, so PSO is better on average in the two denser Data Set I farms.
   - **Fix:** "not significantly different from PSO over all cases (p = 0.30; W/T/L 14/50/4). In an exploratory, post hoc split, PSO-VNS has the lower mean loss in 14 of the 20 cases with N ≥ 10 (mean 0.08 pp; Wilcoxon p = 0.048, unadjusted)." Add a pre-specifiable trend test, e.g. the Spearman correlation of ΔL with N stratified by (DS, r), or a stratified sign-flip test.

4. **[major] 05_setup.tex** — "Case ranks: the Friedman test … Case means: … Wilcoxon signed-rank tests … on the 68 per-case mean wake losses".
   - **Problem:** Every case-level test treats the 68 cases as independent exchangeable blocks. They are 6 clusters (2 wind data sets × 3 radii) with nested N = 2…15 that share site and wind rose. The 6 cases with N = 2 are ties for all methods. The cases are a fixed design, not a sample, so p-values such as 4.6×10⁻⁶⁵ or 1.9×10⁻¹¹ convey false precision and cannot generalize beyond this benchmark.
   - **Fix:** Add to the Statistics paragraph: "The 68 cases are not independent draws: they form six site–wind clusters with nested N. The case-level p-values therefore describe this benchmark, not a population of WFLOP instances, and p-values below 10⁻¹⁰ are reported as such." Add a cluster-level sensitivity analysis, either the 6 cluster means of ΔL per comparison or a within-cluster sign-flip permutation test.

5. **[major] 05_setup.tex / 06_results.tex** — "adjusted over the comparisons of the case and counted as significant wins, ties and significant losses (W/T/L)".
   - **Problem, multiplicity:** Holm controls error only within a case (7 comparisons), not across the 68 cases or across the more than 15 test families in the paper (main, ablation, baseline, old-vs-new, split, budgets, feasible-init, HR, IEA37, robustness, 24 + 6 Mann–Whitney tests). Under the null, about 3 of 68 cases per column can be "significant" by chance. The 4 losses vs PSO, the 1 loss vs VNS and the 2 SSA-VNS wins vs RS-VNS are within that range.
   - **Problem, family definition:** The same pair gets different counts in different tables purely because the Holm family changes: PSO-VNS vs PSO is 14/50/4 in Table wtl and 8/58/2 in Table ablation.
   - **Fix:** In 05_setup, declare the primary endpoint ("PSO-VNS vs each method, case-mean Wilcoxon, Holm over 7") and label all other analyses as secondary or exploratory. Add: "W/T/L counts are descriptive: Holm controls the error only within a case, so a few significant cases per column can arise by chance." In 07_ablation, state explicitly that 8/58/2 and 14/50/4 differ only in the Holm family.

6. **[major] 06_results.tex / Table friedman68** — "its mean wake loss is \NLossPSOVNS\% against \NLossPSO\%, and the difference of the per-case means is not significant".
   - **Problem:** No effect size has an uncertainty interval: not ΔL, not the average ranks, not r_rb. For PSO the mean ΔL (−0.018 pp) and the median (−0.0002 pp) disagree, so the mean is driven by a few cases. The case-level r_rb (0.155 for PSO, from case_mean_wilcoxon) is computed but not reported.
   - **Fix:** Add 95% CIs to Table friedman68 (bootstrap over cases, preferably a cluster bootstrap over (DS, r)) for ΔL, plus a Hodges–Lehmann estimate and the case-level r_rb. Text: "mean difference −0.018 pp (95% CI [a, b]; Hodges–Lehmann −0.000x pp; r_rb = 0.16; p = 0.30)".

7. **[major] 05_setup.tex** — "so that a few favorable feasible runs cannot earn a good rank". The summary's ranking_rule also says "This removes the survivorship bias".
   - **Problem:** The rule only caps survivorship bias. A method with 15–29 feasible runs is still ranked by the mean of its survivors, and a method with 14 feasible runs is ranked below every qualifying method regardless of quality. The 50% threshold is arbitrary.
   - **Problem, ΔL coverage:** ΔL in Table friedman68 is averaged only over cases in which both methods qualify (50 for DE, 66 for SSA and LX-SSA), and the caption does not say so.
   - **Problem, infeasible runs:** In the run-level score (goodness()), infeasible runs are ordered by minimum spacing only, so boundary violations are ignored.
   - **Fix:** Reword to "limits (but does not remove) survivorship bias". Add a sensitivity analysis: thresholds of 80% and 100%, plus a ranking by the per-case median of the run-level score (infeasible = worst). Report that the old rule gives nearly identical ranks (main.friedman_old_rule: SSA 5.35, LX-SSA 6.10, DE 7.03 vs 5.38, 6.06, 7.04), which is reassuring. Add to the caption of Table friedman68: "ΔL over the cases in which both methods have ≥ 15 feasible runs (50 for DE, 66 for SSA and LX-SSA)". Make the caption's "Ranked Last" consistent with the code: "ranked below all qualifying methods, by number of feasible runs".

8. **[minor] 05_setup.tex** — "Wilcoxon signed-rank tests of PSO-VNS against each method on the 68 per-case mean wake losses".
   - **Problem:** case_mean_wilcoxon() imputes a difference larger than any observed one when only one method qualifies (18 cases for DE, 2 each for SSA and LX-SSA). This is reasonable but undisclosed.
   - **Fix:** Append: "a case in which only one of the two methods has at least 15 feasible runs counts as the largest difference in its favor (18 cases for DE, 2 for SSA and LX-SSA)".

9. **[minor] 06_results.tex / Table friedman68** — "The Friedman test \NFriedVerb{} the hypothesis … ($\chi^2_F=\NFriedChi$, $p=\NFriedP$)" and "($p_{\rm Holm}\le\NPostHocSigMaxP$)".
   - **Problem, Iman–Davenport:** The Setup announces the Iman–Davenport correction, but the text reports the χ² p, and the table gives F_F without d.f. or p. The code has iman_davenport_p = 6.7e-109.
   - **Problem, post hoc SE:** The post hoc SE sqrt(k(k+1)/6n) ignores ties, which are frequent for small N. This makes the tests conservative.
   - **Problem, rounding:** The largest significant p_Holm is 2.04×10⁻⁴, so "≤ 2.0×10⁻⁴" is false after rounding.
   - **Fix:** "(Iman–Davenport F_F = 136.5 with 7 and 469 d.f., p < 10⁻¹⁰)". Change to "p_Holm ≤ 2.1×10⁻⁴", or format the bound with round-up. Note that tied ranks make the z tests conservative.

10. **[minor] 06_results.tex** — "PSO-VNS is significantly better in \NWtlSSAVNSW{} cases and never worse (… $p=\NPWrawSSAVNS$)".
    - **Problem:** The text quotes raw p (1.9×10⁻¹¹), Table friedman68 gives Holm p for the same test (5.7×10⁻¹¹), and the abstract uses \NPWPSO (Holm). It is also unclear to the reader which p is meant.
    - **Fix:** Use \NPWSSAVNS / \NPWPSO throughout and write "p_Holm = …".

11. **[major] 08_beyond.tex / 06_results.tex / 10_limits_concl.tex** — "Only the first place of PSO-VNS is stable across budgets" / "has … the best rank at all three budgets".
    - **Problem:** These average ranks come from 10 methods on only 6 cases, with no test. The Nemenyi CD for k = 10, n = 6 is about 5.5 ranks. At 120,030 evaluations the ranks are PSO-VNS 2.17, VNS 3.17 and DE 3.33, which are indistinguishable. On Horns Rev 1 at 120,030, VNS, RS-VNS and DE have lower loss than PSO-VNS.
    - **Fix:** "On the six largest cases PSO-VNS has the lowest average rank at all three budgets (1.50, 1.17, 2.17), but with six cases these differences are not statistically resolvable, and on Horns Rev 1 at 120,030 evaluations VNS, RS-VNS and DE reach lower losses." Add run-level W/T/L of PSO-VNS vs each method per budget to Table feasbudget.

12. **[minor] 07_ablation.tex / 10_limits_concl.tex** — "A longer PSO phase is better" / "PSO-VNS (preferably with $\omega=0.75$".
    - **Problem:** This is a post hoc finding on 12 of the test cases with 3 levels. There is no across-case test. The recommendation is selected in-sample. The natural fourth level, ω = 1 (plain PSO, with runs already available on the same seeds), is omitted. The monotone "more PSO is better" trend would predict PSO ≥ PSO-VNS, which bears on the main claim. Across-case evidence exists but is not reported: case-mean Wilcoxon for 75% vs 50% gives 11 of 12 lower and exact p = 0.0015, and 50% vs 25% gives the same.
    - **Fix:** Add ω = 1 (PSO) as a row in Table split and report the case-mean Wilcoxon. Replace with: "On these twelve cases, a longer PSO phase gave lower losses (11 of 12 cases; Wilcoxon on case means p = 0.0015; run-level 5 significant). As this was assessed on the test cases, ω = 0.75 is a candidate setting to be validated, not a recommendation." Drop "preferably" from the Discussion.

13. **[minor] 09_robust.tex / 01_front.tex** — "\NBestGauss{} ranks first (\NRankGaussPSO) and PSO-VNS \NPosGaussPSOVNS{} … The small advantage of the VNS phase over PSO is thus specific to the wake model".
    - **Problem:** The rank difference is 1.91 vs 1.97 and untested. It is treated as a reversal, and a causal attribution is drawn from it. The abstract lists "PSO ranks first under a Gaussian wake model" as a finding.
    - **Fix:** "Under the Gaussian re-evaluation PSO and PSO-VNS are indistinguishable (average ranks 1.91 and 1.97), both ahead of all other methods, so the small Jensen-model advantage of the VNS phase does not carry over." Apply the same wording in the abstract and the Limitations.

14. **[minor] 01_front.tex / 10_limits_concl.tex / 08_beyond.tex** — "exceeds the energy yield of the installed Horns Rev~1 block" / "its best Horns Rev~1 layout exceeds the installed one" / "significantly higher than that of every other method".
    - **Problem, best-of-30:** At 6,030 evaluations this is a best-of-30 selection. The mean (138.73) is below the installed layout (139.51), and 6 of 30 runs exceed it (95% CI ≈ 8–39%).
    - **Problem, test description:** The HR test is the seed-paired run-level Wilcoxon with infeasible runs ranked last (PSO has 11 infeasible runs, DE has 30), not a test of mean AEP. DE has no mean.
    - **Fix:** Abstract: "at 30,030 evaluations its mean AEP exceeds that of the installed Horns Rev 1 block (at 6,030, 6 of 30 runs do)". 08: "significantly better than every other method in the seed-paired run-level test (infeasible runs ranked last; p_Holm ≤ 6.7×10⁻⁴)".

15. **[minor] 01_front.tex / 06_results.tex (baseline)** — "with it, PSO ranks \NOldPSOPos{} … and salp swarm hybrids appear stronger, whereas with the constriction setting PSO ranks \NPSOPos{}".
    - **Problem, mixed pools:** The sentence mixes two pools. "Fifth / second" refers to the main pool with the old PSO swapped in (5.12). "Hybrids appear stronger" refers to the earlier study's pool (Table baseline, PSO position 5 → 1). Between the two columns of Table baseline, the ranks of all other methods change only mechanically, because PSO moves.
    - **What is sound:** The reversal claim itself holds, because it rests on run-level W/T/L (35/33/0 → 0/21/47), not on the pool-dependent ranks.
    - **Problem, pairing:** The old-vs-constriction Wilcoxon ("significantly worse in 56 … cases") assumes that runs from the earlier code version (fresh_grid.csv) share initial populations with the new runs by seed. This is asserted, not shown.
    - **Fix:** Abstract: "with it, PSO would rank fifth of the eight methods of the main comparison (second with the constriction setting), and in the earlier study's method pool both salp-swarm hybrids rank ahead of it". Add to Table baseline: "All other methods use identical runs in both columns; their rank changes arise only from the change of PSO." State how seed pairing across code versions was verified (e.g., identical initial F_p per seed); if it cannot be verified, use Mann–Whitney.

16. **[minor] 07_ablation.tex** — "The Laplace step of LX-SSA gives no gain" and "LX-SSA-VNS ranks even behind RS-VNS on average".
    - **Problem, Laplace step:** "No gain" is understated. On case means, LX-SSA is worse than SSA in 53 of 62 non-tied cases (+0.166 pp, p = 8×10⁻⁸), and LX-SSA-VNS is worse than SSA-VNS in 47 of 60 (+0.063 pp, p = 5×10⁻⁶).
    - **Problem, "ranks even behind":** This rests on a 0.16 average-rank difference. LX-SSA-VNS actually has a 0.005 pp *lower* mean loss than RS-VNS (p = 0.72).
    - **Fix:** "The Laplace step is detrimental: LX-SSA is significantly worse than SSA in 7 cases and never better, and worse on the case means (p = 8×10⁻⁸); inside the hybrid, SSA-VNS is 3/65/0 against LX-SSA-VNS (p = 5×10⁻⁶)". Replace "ranks even behind RS-VNS" with "does not differ from RS-VNS (p = 0.72)".

17. **[minor] Table wtl caption / 05_setup.tex** — "Last Row: Median Rank-Biserial Correlation".
    - **Problem:** wil() drops zero differences before computing r_rb, so in small-N cases r_rb describes only the few non-tied pairs. Cases in which all differences are zero get r_rb = 0. The median over 68 cases therefore mixes degenerate and informative cases.
    - **Fix:** Add to the caption: "zero differences are dropped (Wilcoxon's method); cases in which all 30 differences are zero have r_rb = 0". Alternatively use Pratt's treatment, and report the case-level r_rb (item 6) as the headline effect size.

18. **[minor] 05_setup.tex** — "reproduces the archived result distributions of the earlier version (22 of 24 Mann--Whitney tests with $p\ge0.05$)".
    - **Problem:** Non-rejection is not evidence of reproduction, because power is unknown and about 1.2 rejections are expected by chance.
    - **Fix:** "…reproduces the archived distributions: across 24 comparisons, the median difference is at most X% of the objective (Vargha–Delaney A12 between a and b); 2 of 24 Mann–Whitney tests reject at α = 0.05, close to the 1.2 expected by chance." Alternatively, use a TOST equivalence test with a stated margin.

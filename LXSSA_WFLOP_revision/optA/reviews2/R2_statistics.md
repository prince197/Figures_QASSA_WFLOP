# R2 -- Statistical review (nonparametric tests, equivalence, Bayesian comparison, multiplicity, clustered data)

Scope: main text `MPCE_PSO_VNS.tex` + `optA/0*.tex`, `optA/10_limits_concl.tex`; macros `analysis/mpce_numbers*.tex`;
supplement `MPCE_PSO_VNS_supplement.tex`, `optA/supp_theory.tex`, `analysis/mpce_supp_*.tex`; code
`analysis/mpce_results.py`, `analysis/mpce_inference_extra.py`. Severity: MAJOR / MODERATE / MINOR.

## What was recomputed independently (from the per-run CSVs)

The recomputation used only `mpce_psoc_s*`, `mpce_psobv_s*`, `mpce_rsdisc_s*` and `fresh_bgrid.csv` (6,030 evaluations, random init, 68 cases). The scripts are in the session scratchpad; the repository is unchanged.

| Quantity | Recomputed | Macro / JSON | Match |
|---|---|---|---|
| PSO-VNS - PSO, n, mean | 68, -0.0180 pp | `\NCmPSOVNSvsPSODL` -0.018 | yes |
| lower / higher / exact ties | 34 / 26 / 8 | 34 / 26 | yes |
| case-mean Wilcoxon p | 0.296 | `\NCmPSOVNSvsPSOP` 0.30 | yes |
| 90 % percentile bootstrap CI | [-0.0407, 0.0040] | `\NEqPSOVNSvsPSOCI` [-0.041, 0.004] | yes |
| 95 % percentile bootstrap CI | [-0.0446, 0.0078] | `\NCmPSOVNSvsPSOCI` [-0.045, 0.008] | yes |
| bootstrap p_TOST | 0.0107 | `\NEqPSOVNSvsPSOP` 0.011 | yes |
| Bayes signed-rank P(rope) | 0.997 (s = 0.5, z0 = 0) | ">0.99" | yes |
| DS II, N >= 10 | 10/10 lower, exact p = 0.00195 | `\NCmPSOVNSvsPSOdsIILargeP` 0.002 | yes |
| mean abs(d), N < 10 | 0.027 | `\NAbsDiffPSOSmall` 0.03 | yes |
| SSA-VNS vs RSD-VNS | 18/41, p = 0.0046, +0.024 [-0.033, 0.068] | same | yes |
| LX-SSA-VNS vs RSD-VNS | 10/50, p = 1.4e-6, +0.087 [0.039, 0.133] | same | yes |
| PSO-VNS vs RSD-VNS | 57/2, p = 3.6e-11, -0.306 [-0.396, -0.227] | same | yes |
| cluster CR1 95 % CI, PSO-VNS - PSO | [-0.0798, 0.0437] | `mpce_summary_extra.json` | yes |

**New numbers not in the paper** (they drive issues 1–5):

- **Cluster-robust 90 % CIs for PSO-VNS - PSO** (the equivalence interval):
  - CR1, t(5): [-0.067, 0.030]
  - CR2, t(5): [-0.067, 0.031]
  - cluster (block) bootstrap: [-0.0525, 0.020]
  - wild-cluster restricted bootstrap-t (Webb 6-point weights, all 6^6 draws enumerated): TOST p = 0.104
  - **Equivalence at ±0.05 pp fails under every cluster-aware method.**
- **Cluster means of d = L(PSO-VNS) - L(PSO), in pp:**
  - DS I: 500 m +0.064; 750 m +0.043; 1000 m -0.008
  - DS II: 500 m -0.063; 750 m -0.031; 1000 m -0.090
  - Three of the six exceed the margin, in opposite directions.
- **Per-case d, dense layouts:**
  - DS I, 500 m, N = 8, 9, 10: +0.160, +0.274, +0.168 (PSO better)
  - DS I, 750 m, N = 11, 12: +0.183, +0.178 (PSO better)
  - DS II, 1000 m, N = 11–15: -0.19 to -0.38 (PSO-VNS better)
  - Overall, 19 of 68 cases have abs(d) > 0.05: 11 favour PSO-VNS and 8 favour PSO.
- **Fixed-benchmark (seed-level) analysis.** Averaging d over the 68 cases within each seed gives 30 paired replicates. Mean -0.018, 90 % t-CI [-0.0265, -0.0096]. Conditional on this fixed benchmark, the difference is statistically nonzero and inside the margin.
- **Bayesian signed-rank test, detail:**
  - Posterior means (theta_A, theta_rope, theta_B) = (0.25, 0.62, 0.12).
  - N < 10 only: P(rope) = 1.00.
  - N >= 10 only (n = 20): P(rope) = 0.05, P(PSO-VNS better) = 0.93. Mean -0.079, 90 % CI [-0.143, -0.017].
  - Comparing sums with r instead of 2r gives P(rope) = 0.97. Prior strength s = 1 gives 0.997.
- **Hodges–Lehmann (HL) estimates, with 95 % bootstrap CIs:**
  - SSA-VNS - RSD-VNS: +0.025 [0.008, 0.060]
  - PSO-VNS - PSO: -0.003 [-0.011, 0.003]; median -0.0002

No mismatch was found between the text and the macros for the numbers checked. The hard-coded "above 0.99" / "> 0.99" is consistent with P(rope) = 0.997.

**Implementation checks that passed:**

- `holm`, `bh` and the Friedman/Iman–Davenport formulas are correct.
- `bayes_signrank` (`mpce_results.py:372`) matches Benavoli et al. (2017) and baycomp `SignedRankTest`:
  - Dirichlet(s, 1, ..., 1) weights on (z0, d_1..d_n), with prior 0.5 at z0 = 0
  - sums over all i, j including i = j and z0
  - threshold +-2·rope, with ties counted as 1/2
  - the reported value is the share of samples in which each theta is the largest
- `tost` (`mpce_results.py:340`): equivalence holds iff the 90 % percentile CI lies inside ±m. This is internally consistent.

---

## Issues

### 1. MAJOR — The practical-equivalence conclusion does not survive the paper's own case-dependence analysis, yet the paper says it does

- **Where the equivalence claim is made:**
  - `optA/01_front.tex:13` ("practically equivalent ... These conclusions hold across qualification thresholds, farm clusters and one Holm family")
  - `optA/02_intro.tex:14`, `:32`; `optA/05_setup.tex:15`, `:19`; `optA/06_results.tex:30`
  - `optA/10_limits_concl.tex:12`, `:38`; supplement `MPCE_PSO_VNS_supplement.tex:303`
- **Problem.** The TOST and the Bayesian test resample or weight the 68 cases as if they were independent. The paper states that the cases form six nested clusters. The cluster analyses (`mpce_inference_extra.py:389`) only re-check significance verdicts; they never re-check the equivalence verdict.
  - The CR1 95 % interval for PSO-VNS - PSO, printed in the supplement's X-loco table, is [-0.080, 0.044]. It lies far outside ±0.05.
  - The recomputed cluster-aware 90 % intervals are CR1/CR2 [-0.067, 0.030] and cluster bootstrap [-0.0525, 0.020]. The wild-cluster TOST p is 0.10. At the cluster level, equivalence cannot be claimed.
  - "Hold across farm clusters" is therefore true only for the non-significance of the difference, not for equivalence.
- **Fix:**
  1. Add a cluster-aware TOST, primary: CR2 with small-sample df, or wild-cluster restricted bootstrap with Webb weights. Add 90 % cluster-aware CIs to `tab:equivalence`.
  2. Remove "farm clusters" from the list of robustness checks that the equivalence statement survives. Alternatively, restate it as: "equivalent at ±0.05 pp when cases are treated as exchangeable; with cluster-robust inference the 90 % interval [-0.067, 0.030] pp is too wide to show equivalence."
  3. Add a CHECK-EXTRA item that tests equivalence with the cluster-robust interval.

### 2. MAJOR — The target of inference is undefined (fixed benchmark vs. population of cases), and the chosen procedure fits neither consistently

- **Where:** `optA/10_limits_concl.tex:30` ("the p-values describe this benchmark"); `optA/05_setup.tex:15` (bootstrap over cases).
- **Problem.** The two readings call for different procedures:
  - If the inference is conditional on the 68 fixed cases (as the Limitations say), the only sampling variability is over seeds. Resampling cases then answers a different question: generalization to other cases.
  - A seed-level analysis of the benchmark-average difference gives 90 % CI [-0.027, -0.010]. That is significantly nonzero and within the margin ("equivalent and different").
  - If the inference is about new cases, the case bootstrap ignores the cluster structure (issue 1).
  - The paper mixes the two. The main-text case-mean Wilcoxon ("no significant difference", p = 0.30) and the run-level 14/50/4 answer yet other questions.
- **Fix.** State the estimand and population explicitly, and report two layers:
  1. **Fixed benchmark:** seed-level (or run-level) paired inference on the benchmark-average difference, e.g. -0.018 pp, 90 % CI [-0.027, -0.010].
  2. **Generalization:** cluster-aware inference.

  A single coherent option is the Bayesian hierarchical correlated t-test (Corani, Benavoli, Demšar, Mangili, Zaffalon 2017, *Machine Learning* 106:1817). It uses the 30 runs per case and models between-case heterogeneity.

### 3. MAJOR — "Practically equivalent on average" hides large, opposite-signed case- and cluster-level differences (case-mix dilution)

- **Where:** `optA/01_front.tex:13`, `optA/02_intro.tex:9`, `:32`, `optA/06_results.tex:30`, `:39`, `optA/10_limits_concl.tex:12`, `:38`.
- **Problem.** The mean over 68 cases is pulled toward zero by the many small-N cases in which both methods reach the same layout: 8 exact ties and 19 cases with abs(d) < 0.005.
  - Case level: 19 of 68 cases exceed the margin (11 favour PSO-VNS by up to 0.41 pp, 8 favour PSO by up to 0.27 pp).
  - Cluster level: three of six cluster means exceed the margin.
  - The Bayesian result is entirely driven by the N < 10 cases. N < 10 alone gives P(rope) = 1.00. N >= 10 alone gives P(rope) = 0.05 and P(PSO-VNS better) = 0.93; its 90 % CI [-0.143, -0.017] lies wholly outside the rope on the negative side.
  - Consequently, "equivalent" in abstract-level sentences is a property of this case mix, not of the two algorithms.
- **Fix:**
  1. Report the distribution of per-case differences: share beyond ±m in each direction, and a quantile plot.
  2. Report equivalence stratified by N (or by the density N·D²/r²).
  3. Qualify every "practically equivalent" as "on average over the 68-case mix, which is dominated by small layouts."
  4. Add a sensitivity analysis with a relative margin, e.g. a percentage of the case's wake loss or per turbine. Cases with near-zero wake loss cannot differ by 0.05 pp.

### 4. MAJOR — The post hoc margin, and its model-shift justification, are not sound

- **Where:**
  - `optA/05_setup.tex:15` ("Jensen and Gaussian wake losses ... differ by \NXModelShiftMean pp (\NXModelShiftRatio times the margin), but the PSO-VNS-PSO difference changes by only \NXModelShiftPairMean pp")
  - `MPCE_PSO_VNS_supplement.tex:302` ("The margin is small against wake-model uncertainty")
  - `analysis/mpce_inference_extra.py:627–661`
  - `analysis/mpce_results.py:172–178`: the comment says the margin was chosen when "not significantly different ... had to be turned into a practical-equivalence statement"
- **Problems:**
  1. The margin was set after seeing the data. The minimal margin, 0.041, is just below the chosen 0.05, so the TOST error guarantee does not hold.
  2. The level shift abs(L_Gauss - L_Jensen) = 0.38 pp is a common-mode bias of the same layouts. It largely cancels in paired comparisons and says nothing about which *between-method* difference is practically irrelevant.
  3. The relevant quantity, the change of the paired difference, is below the margin: 0.031 pp per case and 0.009 pp for the 68-case mean. That shows differences of the margin's size are reproducible across wake models, which argues against dismissing them as within model error. (18 sign changes are in cases with near-zero d.)
  4. The Gaussian model is uncalibrated (k* = 0.04), so the size of the "uncertainty" is arbitrary.
- **Fix:**
  1. Remove the model-shift sentence as a justification of the margin. Keep Table X-modelshift only as evidence that paired differences are model-robust.
  2. Justify ±0.05 pp by an external, data-independent criterion: e.g. 0.05 % of AEP against typical AEP uncertainty or the economic value per farm. State it as a convention that was fixed post hoc.
  3. Make the equivalence curve (p_TOST vs m, both case-level and cluster-robust) the primary display.
  4. In the abstract, write "within ±0.05 pp (post hoc margin; minimal margin 0.041 pp, 0.067 pp cluster-robust)" or drop "practically equivalent" from the abstract.

### 5. MODERATE — Bayesian signed-rank test: correct implementation, overstated interpretation

- **Where:** `optA/01_front.tex:13` ("posterior probability above 0.99"); `optA/02_intro.tex:32` ("Bayesian probability of practical equivalence above 0.99"); `optA/06_results.tex:30`; the caption of `tab:equivalence` generated at `mpce_results.py:505–517`.
- **Problems:**
  - The reported number is P(theta_rope is the largest of the three). It is not the probability that the methods are equivalent, and not that the mean difference lies in the rope.
  - The posterior means of theta are (0.25, 0.62, 0.12). A new pair of cases has an estimated 25 % chance that its Walsh average favours PSO-VNS by more than 0.05 pp.
  - The estimand (Walsh averages, i.e. a pseudo-median functional) differs from the TOST estimand (the mean).
  - The test also treats cases as iid and ignores run-level uncertainty.
- **Fix:**
  1. Rephrase as "posterior probability that practical equivalence is the most probable of the three outcomes."
  2. Add posterior means and 95 % credible intervals of theta_A, theta_rope and theta_B.
  3. Report the N-stratified result.
  4. Prefer, or add, the hierarchical model (issue 2).
  5. State the prior choice (s = 0.5, z0 = 0, as in baycomp) and that the result is insensitive to s (s = 1 gives 0.997).

### 6. MODERATE — Wilcoxon tests are paired with mean effect sizes (estimand mismatch)

- **Where:**
  - `optA/02_intro.tex:22` ("the case means, not the clusters, favor RSD-VNS: \NCmSSAVNSvsRSDVNSDL pp")
  - `optA/07_ablation.tex:14`
  - `optA/05_setup.tex:15` (Wilcoxon primary, means with bootstrap CIs)
  - supplement `:291`
- **Problem.** For SSA-VNS vs RSD-VNS, the Wilcoxon test is significant (p = 0.0046, p_Holm = 0.014; 41 vs 18 cases), but the reported mean +0.024 pp has a 95 % CI [-0.033, 0.068] that contains 0. The Wilcoxon test is about symmetry / pseudo-median. Its matching estimate is Hodges–Lehmann: +0.025, 95 % CI [0.008, 0.060].
  - For PSO-VNS vs PSO the mean (-0.018) and the HL estimate (-0.003) differ six-fold.
  - Readers are shown a significant test next to a non-significant effect size.
- **Fix:**
  1. Report the HL estimate with its CI next to every Wilcoxon p, or use a mean-based test (bootstrap/t) for the headline effects.
  2. For "SSA adds nothing beyond random sampling", use a one-sided non-inferiority statement: "the SSA phase's gain over RSD-VNS is at most 0.022 pp (upper end of the 90 % CI)". The data support this directly.

### 7. MODERATE — The cluster-level tools are weak and are applied asymmetrically

- **Where:** `analysis/mpce_inference_extra.py:389–435`; `analysis/mpce_supp_inference.tex:45`, `:78`; supplement `:320`; `optA/05_setup.tex:19`.
- **Problems:**
  1. CR1 with G = 6 and t(5) over-rejects when clusters are few and heterogeneous, and the variance here sits in the large-N cases (Cameron–Gelbach–Miller 2008; MacKinnon–Webb 2017; Imbens–Kolesár 2016).
  2. The block bootstrap with 6 clusters is too narrow (the paper acknowledges this).
  3. The exact sign, Wilcoxon and sign-flip tests have minimum p = 2/64 = 0.031. A single dissenting cluster pushes p to at least 0.22, so they have almost no power.
  4. The (data set, radius) clusters ignore crossed dependence. DS I and DS II share the identical (r, N) geometry and seeds, and nested N makes neighbouring cases similar.
  5. The cluster analysis is used to downgrade unfavourable claims (SSA-VNS worse than RSD-VNS; LX-SSA-VNS "worse", not "clearly worse"). It is not applied to the favourable ones: equivalence (issue 1) and the DS II N >= 10 subgroup, which rests on effectively 3 clusters, one containing a single N >= 10 case, so a cluster sign test cannot go below p = 0.25.
- **Fix:**
  1. Replace CR1 with CR2 and Bell–McCaffrey / Imbens–Kolesár degrees of freedom, or with a wild-cluster restricted bootstrap-t using Webb 6-point weights (all 46,656 draws enumerable).
  2. Model the design explicitly: d ~ dataset × radius + f(N), or a mixed model with cluster and (r, N)-geometry effects.
  3. Present cluster results as sensitivity analyses and apply the same criteria to every headline claim.

### 8. MODERATE — The post hoc subgroup (DS II, N >= 10) is promoted to a finding; the opposite DS I pattern is under-reported

- **Where:** `optA/01_front.tex:13`; `optA/02_intro.tex:32`; `optA/06_results.tex:30` ("Data Set I shows no such effect"); `optA/10_limits_concl.tex:12`, `:38`; `mpce_inference_extra.py:515–605`.
- **Problems:**
  - The cut-point (N >= 10) and the data set were chosen after inspecting the data. Holm over 28 tests does not adjust for this selection among implicit alternatives (other cut-points, radii, data sets).
  - In DS I the trend goes the other way: stratified rho = +0.34, p = 0.065. The DS I N >= 10 mean is +0.052 pp (above the margin, favouring PSO). The 500 m and 750 m N >= 10 cluster means are +0.17 and +0.14. With a cut at N >= 8, DS I 500 m adds +0.16 and +0.27.
  - The honest summary is a data set × density crossover interaction, not "gains in the larger DS II layouts".
  - The interaction permutation test compares raw mean differences under strong heteroscedasticity (d is about 0 for small N). It therefore tests exchangeability rather than equal means. A studentized (Welch) statistic is needed.
- **Fix:**
  1. Report both directions symmetrically, with one pre-stated model (d ~ dataset × density, cluster fixed effects).
  2. Label the subgroup "exploratory" in the abstract, or remove it there.
  3. Use a studentized permutation statistic.
  4. Ideally, confirm on held-out cases defined before running: e.g. additional N or radii for both data sets, since run costs are small.

### 9. MODERATE — The multiplicity claim is overstated

- **Where:** `optA/02_intro.tex:14` ("one Holm family over all main-text tests"); `optA/01_front.tex:13`; `optA/05_setup.tex:19` is correctly worded.
- **Problem.** The family (`mpce_inference_extra.py:475`) holds only the 28 case-mean Wilcoxon tests. It excludes:
  - the Friedman post hoc tests and the run-level tests
  - the 20 TOSTs
  - the interaction and Spearman permutation tests (quoted in the main text with p = 1e-4 and 3e-4)
  - the SD and worst-run tests, and the Horns Rev tests
- **Fix.** Say "all 28 case-mean Wilcoxon tests". Either add the subgroup interaction and trend tests to the family, or label them exploratory.

### 10. MODERATE — The causal attribution "the Laplace step is harmful" is not identified by the design

- **Where:** `optA/01_front.tex:13`; `optA/02_intro.tex:9`, `:22`; `optA/10_limits_concl.tex:38`. It contradicts `optA/07_ablation.tex:24`: "these contrasts do not separate the Laplace step from this cost". LX-SSA uses 2N_p evaluations per iteration, so it runs half the iterations.
- **Fix.** Write "LX-SSA (Laplace step, at twice the cost per iteration) is worse". Alternatively, add an equal-iteration or equal-cost-per-iteration control.

### 11. MODERATE — "Feasibility reliability" rests on one site under one initialization and is untested

- **Where:** `optA/01_front.tex:13`, `optA/02_intro.tex:32`, `optA/06_results.tex:39`, `optA/10_limits_concl.tex:12`, `:38`; macros `\NHRFeasPSOVNS` 30/30, `\NHRFeasPSO` 19/30, `\NHRFeasFeasInitPSO` 30/30.
- **Problem.** On all 68 benchmark cases both methods are 100 % feasible. On Horns Rev with feasible starts PSO is also 30/30.
- **Fix:**
  1. Report the paired exact McNemar test for the Horns Rev block (11 discordant seeds, p ≈ 0.001).
  2. Qualify the claim as "one 16-turbine site with random starts".
  3. Do not generalize it in the abstract.

### 12. MINOR — TOST details

- **Where:** `mpce_results.py:340–369`; `optA/supp_theory.tex:214–223`; the `tab:equivalence` caption.
- **Problems:**
  1. The "bootstrap TOST p" is the share of bootstrap means beyond ±m: an inversion of the percentile interval (an achieved significance level). It is not a p-value calibrated under the null boundary, and the add-one correction is ad hoc.
  2. Lemma S-tost-ci proves the t-based equivalence, but the decision uses the percentile bootstrap. Percentile intervals under-cover for n about 68 with heavy tails. Here the t interval ([-0.0406, 0.0045]) and the BCa interval ([-0.043, 0.002]) give the same verdict.
  3. Twenty equivalence tests are reported, but only the PSO-VNS vs PSO test is confirmatory.
- **Fix.** Declare one primary TOST (t-based, as in the lemma, or BCa). Relabel the bootstrap quantity. Mark the other equivalence rows as descriptive.

### 13. MINOR — The same pair gets two W/T/L counts and two Holm p-values

- **Where:** `analysis/mpce_tab_wtl.tex` (PSO-VNS vs PSO 14/50/4, Holm over 7) vs `analysis/mpce_tab_ablation.tex:28` (8/58/2, Holm over 15; p_W 0.592 vs 0.30 raw); `optA/06_results.tex:30`, `optA/07_ablation.tex:5`, `:8`.
- **Problem.** The caption notes the different family, but W/T/L then depends on how many comparisons share the case.
- **Fix.** Use one family for this pair, or report unadjusted per-case outcomes plus one benchmark-level adjustment.

### 14. MINOR — Imputation and the ranking of infeasible runs

- **Where:** `mpce_results.py:560–600` (`big = 10·max abs(d) + 1`, line 581); `optA/05_setup.tex:15`; `goodness`, `mpce_results.py:220`.
- **Problems:**
  - The rule is valid for a rank test, but for DE (18 of 68 cases imputed) the case-mean Wilcoxon becomes largely a feasibility test.
  - The mean ΔL excludes those cases (survivorship; this is stated).
  - Infeasible runs are ordered only by minimum spacing; boundary violations are ignored.
- **Fix.** Name the rule as a prioritized composite endpoint (feasibility first, then loss). Put `p_both_qualified` next to it in the table. Run a sensitivity analysis that orders infeasible runs by total violation.

### 15. MINOR — Handling of zero differences

- **Where:** `wil`, `mpce_results.py:225`; and the same helper in `mpce_inference_extra.py`.
- **Problem.** Zeros are dropped (the "wilcox" method); there are 8 exact zeros for PSO-VNS vs PSO.
- **Fix.** Add a Pratt-method sensitivity check, at least for the headline pairs.

### 16. MINOR — Wording inconsistent with the statistics

- `optA/07_ablation.tex:8`: "switching to VNS gives a small gain (... p = 0.30) that is practically equivalent to zero". Better: "no practically relevant gain".
- `optA/06_results.tex:30`: "For N < 10 both reach nearly the same wake loss (mean absolute difference 0.03 pp)". The mean is correct, but N < 10 includes per-case differences up to 0.27 pp (DS I 500 m, N = 9) and -0.17 pp (DS II 500 m, N = 9). Add "with exceptions at N = 8–9 in the 500-m farms".
- `optA/06_results.tex:30`: "the gain grows with N within every Data Set II farm". The sign holds in all three clusters (rho = -0.60, -0.70, -0.67), but the trend is not monotone: PSO is better at N = 4–7 in DS II. "Is larger for the largest N" is more accurate.

### 17. MINOR — Selection in the split study

- **Where:** `optA/07_ablation.tex:32`.
- **Problem.** ω = 0.75 is "best" among 5 settings on the same 12 cases; the p-values are unadjusted except in the all-family table.
- **Fix.** This is already hedged ("a candidate to be confirmed"). Also state that the Holm-adjusted ω = 0.9 vs 0.5 comparison is non-significant.

---

## Conclusions that ARE well supported

These hold at the case, run, threshold, Holm (all-family) and cluster levels:
- PSO-VNS beats RSD-VNS, RS-VNS, VNS, SSA-VNS, LX-SSA-VNS, SSA, LX-SSA, DE and MS-SLSQP: 6/6 clusters, cluster-robust p <= 0.006, whole CIs far from 0.
- The VNS phase improves both salp-swarm methods.
- RSD-VNS beats RS-VNS.
- The old PSO setting is worse than constriction: 56 of 68 cases, never better.

The ranking conclusions are reported with appropriate hedges ("average ranks not tested", "exploratory").

## Overall verdict

**Major revision (statistics).** The pipeline is careful, transparent and reproducible:
- I reproduced every checked number from the run files.
- The Bayesian signed-rank implementation matches Benavoli et al. and baycomp.
- Holm, BH and Friedman are correct.
- The TOST decision is internally consistent.

The weak point is the headline claim "PSO-VNS and PSO are practically equivalent (posterior probability above 0.99) ... and this holds across farm clusters":
- The margin was chosen post hoc, just above the observed minimal margin.
- The model-shift justification does not support it.
- The claim fails under cluster-aware inference (90 % CI about [-0.067, 0.030]).
- It averages over opposite-signed differences that often exceed the margin.
- The Bayesian probability is misread.

The favourable post hoc subgroup is promoted to the abstract while its mirror image in DS I is under-reported.

## Top 5 fixes

1. **Make equivalence cluster-aware and define the estimand** (issues 1, 2). Report a fixed-benchmark seed-level CI ([-0.027, -0.010]) and a cluster-robust TOST (CR2 / wild-cluster bootstrap). Drop "hold across farm clusters" for the equivalence claim.
2. **Rework the margin** (issue 4). Remove the model-shift justification. Anchor ±0.05 pp externally and present the equivalence curve, i.e. m_min both case-level (0.041) and cluster-robust (about 0.067). Flag the margin as post hoc in the abstract.
3. **Report heterogeneity honestly** (issues 3, 8). Show the per-case distribution (19/68 beyond ±m, both directions) and stratified equivalence by N or density. Present one pre-stated data set × density interaction model with both DS I and DS II directions. Demote the DS II N >= 10 gain to exploratory, or confirm it on held-out cases.
4. **Correct the Bayesian reporting** (issue 5). Rephrase "probability of equivalence". Add the theta posterior means (0.25/0.62/0.12) and the N-stratified result. Prefer the hierarchical Bayesian model on run-level data.
5. **Align estimands and remove overclaims** (issues 6, 9, 10, 11):
   - add HL estimates with CIs next to the Wilcoxon p-values
   - "all 28 case-mean Wilcoxon tests", not "all main-text tests"
   - "LX-SSA is worse", not "the Laplace step is harmful"
   - feasibility reliability qualified as a single-site result, with a McNemar test

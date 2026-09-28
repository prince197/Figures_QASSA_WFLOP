# Review R3: metaheuristics and PSO theory

Manuscript: "Baseline Configuration and Controls in Metaheuristic Comparisons for Wind Farm Layout Optimization, with a PSO-VNS Reference Method". Reviewed files: `MPCE_PSO_VNS.tex` with `optA/01-11`, macros in `analysis/mpce_numbers*.tex`, supplement theory in `optA/supp_theory.tex`, diagnostics in `analysis/mpce_supp_diag.tex` and `analysis/mpce_supplementary.tex`, and code in `analysis/authors_optimizers.py`, `hybrid_lxssa_bvns.py`, `original_vns.py`, `rs_vns.py`, `mpce_experiments.py` (RSDVNS, PSOBV90), `mpce_diagnostics.py`, `mpce_results.py` and `make_theory_figures.py`.

I re-derived the propositions by hand. I recomputed every typed constant independently of `make_theory_figures.py` (NumPy/SciPy, my own Monte Carlo) and checked the code paths behind the propositions.

---

## 0. What I verified (no action needed)

**Stability proposition (Prop. 1 / S1): correct.**
- The per-coordinate recursion matches the code. `r1 = np.random.rand(dim)` means ρ is drawn per coordinate and per move.
- Mean: μ = (c1+c2)/2. Variance: σ² = (c1²+c2²)/12.
- The mean recursion and its Schur-Cohn conditions give the order-1 condition c1+c2 < 4(1+w).
- The second-moment matrix M, its characteristic polynomial and χ(1) all check out: χ(1) = (1+w)(1−w²) − (1+w)σ² − (1−w)m² = 2μ(1−w²) − (1−w)μ² − (1+w)σ². So does χ(−1).
- Jury's cubic conditions check out, and so does the proof that (iv) follows from (i), including the w < 0 branch.
- With c1 = c2 = c, χ(1) = c[12(1−w²) − c(7−5w)]/6. This is Poli's bound.
- The claim that c1 ≠ c2 is stricter holds, because σ² is minimised at c1 = c2 for fixed c1+c2.
- Cov(φ, η) = c1c2(c1−c2)(p−g)/(12(c1+c2)), which vanishes only if c1 = c2.
- V∞ = (1+w)q/χ(1), with q = c1²c2²(p−g)²/(6(c1+c2)²). For c1 = c2 this reduces to Poli's variance factor.
- A Monte Carlo check with c1 ≠ c2 (w = 0.6, c1 = 1.0, c2 = 1.8, 2·10⁵ particles × 400 iterations) gave mean 0.357 (predicted 0.357) and variance 0.2495 (predicted 0.2490).

**Typed constants: all correct to the printed digits.**

| Constant | Typed | Recomputed |
|---|---|---|
| 24(1−0.7²)/(7−3.5) | 3.50 | 3.4971 |
| 24(1−0.7298²)/(7−5·0.7298) | 3.35 | 3.3475 |
| 4(1.7) | 6.8 | 6.8 |
| 4(1.7298) | 6.92 | 6.9192 |
| 2·1.49618 | 2.99 | 2.99236 |
| Margin, constriction | +0.36 | +0.355 |
| Margin, old setting | −0.50 | −0.503 |
| Relative margins | 10.6% / 14.4% | 10.61% / 14.38% |
| ρ1 = √w | 0.84 / 0.85 | 0.8367 / 0.8543 |
| ρ2 (spectral radius of M) | 1.12 / 0.94 | 1.1237 / 0.9442 |
| Left eigenvector u (old setting) | (0.90, 0.21, 0.39) | (0.897, 0.207, 0.391) |
| V∞/(p−g)² (constriction) | 1.09 | 1.0874 |
| SD factor | 1.04 | 1.0428 |
| κ·φ = 0.72984·2.05 | 1.49618 | 1.49618 |
| Zaharie √((1−0.45)/30) | 0.14 | 0.1354 |
| ξ1(T/2) = 2e⁻⁴ | 0.0366 | 0.03663 |
| E K = 1+0.9·29 | 27.1 | 27.1 |
| h̄ = 0.05/2⁵ | 0.0015625 r | six step sizes, 0.05/2⁶ < 10⁻³ |
| (π/4)^N and the bound (b) for N = 2, 10, 12, 15, 16, 20, 36 | 0.62/0.52, 0.089/6.5e−3, 0.055/0.012, 0.027/7.0e−3, 0.021/0.011, 8.0e−3/7.0e−4, 1.7e−4/4.3e−5 | identical |
| Exact P2 (r = 500, s = 308) | 0.443 | 0.7188·(π/4)² = 0.4434 |
| Monte Carlo P16 (IEA37-16) | 1.5e−4 | 1.48e−4 (my own 10⁶ draws) |
| Clopper-Pearson "0 of 10⁷" ends and "expected feasible of 3,015" | as typed | consistent |

**Box proposition (S3):**
- The lens formula at centre distance d = r is correct.
- Monotonicity in the centre distance is correct: the derivative equals minus the chord length.
- The disjoint-half-disc argument for |U_k| ≥ kA(s/2) is correct, and so is the chain rule.

**Prefix proposition (S2): holds in the code.**
- `HybridBVNS.__init__` seeds once, and `PSO(seed=None)` does not reseed.
- `init_pop` and the per-particle `rand(dim)` pairs are drawn in the same order as in stand-alone `PSOC`.
- `BVNS.search` initialises `best = (G, F(G))`.
- Part (c) holds because SSA's ξ1 depends on t/T.

**Remaining lemmas and propositions:**
- Clip lemma: correct. `np.clip` leaves V unchanged, and p, g ≤ r because they are evaluated positions.
- Feasible-start lemma and Proposition S-frozen: correct.
- DE move structure: matches the code (`rand() < CR or j == j_rand`).
- VNS local-search proposition: correct. The visited coordinates lie on a union of three lattices intersected with [−r, r], and the last sweep at h̄ has no strict improvement.

I found **no mathematical error**. The problems are minor proof gaps (issues 12-13) and, more importantly, how far the theory and the diagnostics support the causal claims in the abstract and conclusions.

---

## 1. Numbered issues

Severity scale: **Major**, **Moderate** or **Minor**.

### 1. [Major] Two leaps in the causal claim "violates order-2, so the swarm does not contract and most evaluations are infeasible"

**Location:** Abstract; Intro contribution 2; Sec. 6.1; Discussion; Conclusion ("which violates the order-2 stability condition and rarely samples feasible layouts").

**Problem.** Proposition 1 assumes stagnation and no bounds. The code clips positions and keeps the velocity, which is absorbing without velocity reset (Lemma S-clip). The supplement itself says the moment results "no longer describe the trajectory" once a bound is hit. Yet the old setting clips 8.7% of the coordinates per iteration. With n = 20-30 coordinates, that is about 1 − 0.913^30 ≈ 93% of layouts per iteration with a clipped coordinate. By Lemma S-clip(a), each such layout is (almost surely) outside the disc. The diagnostics agree: 86.3% of the old setting's infeasible evaluations violate both constraints.

The observed non-contraction and infeasibility therefore come from at least two confounded factors:
- larger random steps (c = 2 against 1.496, and w = 0.7 against 0.7298);
- bound handling that pins coordinates at the box, the effect studied by Helwig, Branke and Mostaghim (IEEE TEC 2013).

A two-point comparison of settings cannot attribute the effect to crossing the order-2 boundary. Moreover, the historical c = 2 PSOs (Kennedy and Eberhart 1995; Shi and Eberhart 1998) used velocity clamping (Vmax), which the code lacks. The "old setting" as run is "w = 0.7, c = 2, no Vmax, absorbing clip without velocity reset".

**Fix. (a)** Add a cheap dose-response run on the six diagnostic cases × 10 seeds: w = 0.7 with c1 = c2 ∈ {1.5, 1.6, 1.7, 1.75, 1.8, 1.9, 2.0}. The order-2 boundary is at c = 1.7486. Report spread, clipped share, feasible share and final loss against c. A jump at about 1.75 would support the order-2 story; a smooth trend would not.

**(b)** Run the old setting with (i) velocity set to zero on clip, (ii) reflection and (iii) Vmax = 0.2·(ub−lb). This separates order-2 instability from bound handling.

**(c)** Until then, replace "so" and "which violates ... and rarely samples" with "is consistent with". State in Sec. 4.1 that the old setting is run without Vmax.

### 2. [Major] Phase 2 is a truncated best-improvement compass search; "VNS" conclusions are conclusions about this local search

**Location:** Sec. 4.3 and 5.3 (PSO-VNS paragraph); Sec. 7 ("VNS phase"); Abstract ("VNS acts as a local-search descent"); Prop. S-ls; Table D-vns.

**Problem.**
- For N ≥ 12, the first descent is unfinished in 39 of 40 runs, and the median run never shakes.
- Best-improvement polling costs 2n evaluations per move: 60 at N = 15. The 3,000 Phase-2 evaluations therefore allow at most about 50 moves.
- The ℓ∞ shells (0.1kr to 0.5r on all 2N coordinates) are enormous next to the 308-m spacing, which explains why only k = 1 is ever accepted.
- The metric "shaking uses 0.02% of Phase-2 evaluations" is misleading, because a shake costs one evaluation. The relevant figure is the cycles (shake plus local search): about 21.5% of the Phase-2 evaluations (100 − 78.5) for 0.5% of the gain.
- Prop. S-ls gives a certificate only after a completed descent, so it says nothing about the large cases where the claimed Data Set II gains live.

Consequently the findings "VNS gives PSO only a small gain", "PSO-VNS ≈ PSO" and "neighborhoods matter only for small N" describe this particular inefficient local search. They are not findings about VNS as such.

**Fix.**
- Report "cycles use X% of the Phase-2 evaluations and give Y% of the gain" instead of the 0.02% shaking share.
- Say in the main text that Prop. S-ls is vacuous for N ≥ 12 at B = 6,030.
- Either rename the hybrid (for example "PSO + compass search", with VNS as the nominal framework) or qualify every VNS claim.
- Add one control with an efficient local search on the 12 split cases, such as opportunistic (first-improvement) randomized-order polling or turbine-wise k-neighborhoods (Cazzaro and Pisinger 2022). This shows whether "PSO-VNS ≈ PSO" survives a better Phase 2.

### 3. [Major] Comparison fairness: only PSO's configuration was examined

**Location:** Sec. 5.1; Limitations; Abstract ("We compare eight methods ..."); Intro contribution 4 ("best average rank").

**Problem.** The paper's thesis is that baseline configuration changes conclusions, yet:
- DE is run untuned. With F = 0.5 and CR = 0.9, 27 of 30 coordinates come from inter-layout combinations. DE is feasible in 1.1% of runs on the largest cases and 71.4% overall, and 18 cases are imputed.
- MS-SLSQP uses forward differences on a piecewise-constant (top-hat) objective, so its gradients are zero or undefined on plateaus.
- CMA-ES, the de facto standard for continuous black-box optimization at about 200n evaluations, is absent. So is L-SHADE.

Limitations mentions the missing CMA-ES, but the ranking claims ("ranks first of eight", "PSO-VNS ... best average rank") are not qualified.

**Fix.**
- Add CMA-ES (pycma defaults: λ = 4 + ⌊3 ln n⌋, σ0 = 0.3r, the same F_p, box via `BoundTransform`) on the 68 cases × 30 seeds. This costs about as much as one existing method.
- Add a small DE sensitivity run on the 12 split cases (CR ∈ {0.1, 0.5, 0.9} × F ∈ {0.5, 0.8}), and/or DE with Deb's feasibility rules.
- Otherwise, state in the abstract that the ranking is among untuned default baselines without CMA-ES.

### 4. [Moderate] The PSO phase's advantage is conditional on random (infeasible) starts

**Location:** Abstract ("PSO gives a clear gain"); Conclusion, second finding; Discussion ("a constriction PSO is a sound first choice and PSO-VNS a more reliable one").

**Problem.** From feasible starts, PSO never moves (0 of 210 stored runs), and PSO-VNS (3.94%) is worse than VNS alone (3.81%). That start costs no objective evaluations. Under random starts, PSO-VNS is feasible at the switch in 98.7% of runs against 79.2% for RSD-VNS, and the mean loss of feasible switch points is 1.63% against 1.59%.

The data thus suggest that the PSO phase mainly buys feasibility. The paper says this is "open" in Sec. 7, but the Conclusion and Recommendations present the PSO phase as generally better.

**Fix.**
- Qualify the second finding and the recommendation with "from random starts".
- Add the obvious practical corollary: with a cheap feasibility-preserving initialization, VNS alone is best.
- As a small analysis, compare PSO-VNS with RSD-VNS only on case-seed pairs where both are feasible at the switch. This separates "feasibility" from "better start".

### 5. [Moderate] "Salp-swarm phases add nothing beyond random sampling" and "the Laplace step is harmful" are over-generalised

**Location:** Abstract; Intro (thesis paragraph and contribution 3); Conclusion.

**Problems.**
- **(a) Budget and initialization.** The statement is measured at 6,030 evaluations, ω = 0.5, random starts and the local search of issue 2. At 120,030 evaluations SSA-VNS ranks ahead of PSO, and RSD-VNS was not run at 30k or 120k. From feasible starts, SSA-VNS (3.96%) essentially ties PSO-VNS (3.94%).
- **(b) Laplace step.** The abstract's "(the Laplace step is harmful)" and the Intro's "makes the hybrid worse" contradict Sec. 7: "these contrasts do not separate the Laplace step from this cost".
  - Per iteration, LX-SSA spends 30 evaluations on 15 followers × 2 candidates.
  - It then spends 30 more on re-evaluating the population, 15 of which are duplicates of already-known follower values.
  - So a quarter of its budget is pure duplication, and the published implementation, not the operator, may be what is harmful.

**Fix.**
- Add "at 6,030 evaluations from random starts" to the abstract and conclusion statements.
- Replace "the Laplace step is harmful" with "LX-SSA as published (2N_p evaluations per iteration) is worse".
- Alternatively, run a cost-matched ablation: LX-SSA with cached fitness (1.5N_p per iteration), or SSA whose followers use only the Laplace candidate (N_p per iteration).
- Also state plainly that on the case means RSD-VNS is Wilcoxon-better than SSA-VNS (18 against 41 cases, p_Holm = 0.014) while the mean-difference CI contains 0. "Adds nothing" is then the conservative reading.

### 6. [Moderate] The RS-VNS weakness is attributed to (π/4)^N, but the data show the spacing constraint dominates

**Location:** Sec. 7 ("a square sample has all N turbines inside the circle only with probability (π/4)^N, **so** 730 of its 2,040 runs draw no feasible sample"); Intro contribution 3; Sec. 4.1 after Prop. 2; Discussion.

**Problem.** The "so" is a non-sequitur, as the supplement's Table D-rsreplay shows:
- At r = 500 m and N = 8, 14.2% of square samples (about 430 per run) have all turbines inside, yet 0 are feasible. Disc sampling also gives 0%.
- At N = 15, r = 1000 m, (π/4)^N · 3,015 ≈ 80 samples per run lie fully inside.
- RSD-VNS, which never violates the boundary, still has 494 runs without a feasible sample.

At most 730 − 494 = 236 of the 730 runs can be attributed to square rather than disc sampling. The rest is packing density.

**Fix.** Reword along these lines: "a square sample is feasible with probability at most (π/4)^N, and far less because of the spacing constraint (Table S-th-box). 730 RS-VNS runs and 494 RSD-VNS runs contain no feasible sample." Make the same change in the Intro and the Discussion.

### 7. [Moderate] "Spacing-only" classification of the old setting's infeasible finals uses a tolerance inconsistent with the penalty

**Location:** Sec. 6.1 ("mostly violate only the spacing constraint (246 of 282)"); Sec. 4.1; `mpce_diagnostics.stored_violation_summary`.

**Problem.**
- The stored layouts are classified with a 1e−3 m radial tolerance, because coordinates are stored to 1 mm.
- For a turbine pinned at x = ±r by clipping, g^b = y². A radial excess of 1e−3 m then allows |y| up to √(2r·10⁻³) ≈ 1.0-1.4 m, that is g^b up to about 1-2 m² and a penalty up to about 10²⁰. In the optimizer's metric this is a large boundary violation.
- The table note reports 56 infeasible old finals with a coordinate on the box, but only 36 (6 + 30) are counted as boundary-violating. So at least 20 of the 246 "spacing-only" finals have a turbine pinned at a box tangent point.

The m² scaling makes these tangent points cheap: the penalty is second order in y. This is a penalty-scaling artefact, not a property of the swarm, and it should not appear alongside the spread and feasibility evidence as part of the "mechanism".

**Fix.**
- Reclassify using the instrumented reruns, which hold unrounded coordinates, with the penalty-consistent threshold of about 5·10⁻⁸ in g units.
- Report "pinned at a box tangent point" as its own category.
- Keep the m²-versus-m explanation, but present it as a consequence of penalty scaling. Suggest that future work normalises violations, for example with (‖x_i‖ − r)/r and (ℓ_min − ℓ_ij)/ℓ_min.

### 8. [Moderate] "Practically equivalent on average" rests on a post-hoc margin and on the easy cases

**Location:** Abstract; Sec. 6.2; Sec. 7; Table S-equivalence.

**Problem.**
- The margin of 0.05 pp was fixed after the primary analysis, and the minimal margin is 0.041 pp.
- 48 of the 68 cases have N < 10, where the mean |ΔL| is 0.03 pp and both methods reach nearly the same layout. These cases pull the mean toward 0.
- For the 20 cases with N ≥ 10, the mean difference is −0.079 pp. That exceeds the margin (p = 0.048 unadjusted, 0.29 in the all-family Holm).
- In Data Set II with N ≥ 10, the difference is −0.21 pp, four times the margin.

"Equivalent on average" is therefore true for this case mix, not for the problems on which an optimizer matters.

**Fix.** Report TOST and m_min for N ≥ 10 alongside the 68-case result, state the heterogeneity in the abstract sentence itself, and avoid "PSO-VNS and a well-configured PSO are practically equivalent" without "on this case mix".

### 9. [Minor] Switch point of RS-VNS and RSD-VNS is misreported (3,030 instead of 3,015)

**Location:** `mpce_results.switch_call` (maps "RS" to 30 + round(99.5)·30 = 3,030); `mpce_supplementary.tex`, Table tab:switch; macros `NSwitchFeasRSVNS` (66.7) and `NSwitchFeasRSDVNS`.

**Problem.** RS-VNS switches at call 3,015, but the feasibility is read at the checkpoint for call 3,030, after 15 compass-search evaluations. I checked this on `mpce_rsvns_s0of1.csv`:
- At call 3,000 (index 99), 730 runs are infeasible, so 64.2% are feasible. This matches the replay.
- At call 3,030 (index 100), 680 runs are infeasible, so 66.7% are feasible.

The "made feasible" count (668) and the Phase-2 shares are affected in the same way. The main text (Sec. 4.2) correctly says 3,015.

**Fix.** Use the switch at 3,015 (checkpoint index 99) for RS and RSD, and regenerate the table and macros.

### 10. [Minor] The Zaharie remark is framed as a defect

**Location:** Sec. 4.1, last sentence ("Likewise, the DE scale factor 0.5 exceeds Zaharie's critical value ...").

**Problem.** F > F_crit means the expected population variance does not collapse under mutation and crossover. That is the recommended side of the threshold, not an instability analogous to the old PSO setting. "Likewise" invites the wrong reading.

**Fix.** Rephrase: "F = 0.5 exceeds Zaharie's critical value 0.14, so DE's weak results are not due to premature variance loss." Optionally note that for binomial crossover the effective probability is CR(1 − 1/n) + 1/n. The numerical change is negligible.

### 11. [Minor] Penalty threshold and wake-loss bounds are stated inconsistently

**Location:** Sec. 3.3 and S-th-budget ("about 10⁻⁷"; "0 ≤ L ≤ N E_ideal"); Lemma S-feas-improve ("plus the small negative part of the power curve"); main-text Proposition prop:feasible-start ("only if it is feasible").

**Problem.**
- The threshold is (√2.2·10⁵ − 1)/10¹⁰ = 4.6·10⁻⁸, not about 10⁻⁷.
- The power curve is negative just above cut-in (f(3.5) = −6.99 kW), so L can slightly exceed N E_ideal.
- A trial can replace G if its violations are below this threshold, so "only if it is feasible" is too strong.
- The reporting tolerance (1e−6 m) is far looser than the penalty's effective tolerance. `curve_valid` in `mpce_results.py` already has to discard penalized values accepted as "feasible".
- For F_p above about 10²¹, the wake loss is below the float64 ulp, so infeasible comparisons are penalty-only in the strict numerical sense.

**Fix.**
- Write "about 5·10⁻⁸".
- Write "L ≤ N E_ideal + ε".
- Use "feasible up to violations below 5·10⁻⁸" in Prop. prop:feasible-start.
- Add one sentence on the tolerance mismatch.

### 12. [Minor] Proof gap in Prop. S1(b), necessity

**Problem.** Order-2 stability is defined through E x_t² only. The necessity argument assumes that the whole vector z_t, including E y_t y_{t−1}, converges. This follows from row 1 of M when w·m ≠ 0; the edge cases w = 0 or m = 0 need one line.

**Fix.** Either define order-2 stability through the full second-moment vector, or add the row-1 argument and treat w·m = 0 separately.

### 13. [Minor] Proof gap in Prop. S-settings: unboundedness is shown, not E x_t² → ∞

**Problem.** u^T z_t → ∞ only gives E y_t² + E y_{t−1}² → ∞.

**Fix.** Use the right eigenvector of the simple, real dominant eigenvalue, v = (ρ, mρ/(ρ+w), 1) with v1 = ρ > 0. The other eigenvalues have modulus 0.55 < ρ, so z_t/ρ^t → c·v with c > 0, hence E y_t² ~ cρ^{t+1} → ∞. Also, in the main text, cite Table S-th-settings rather than the figure for "grow by the factor 1.12 per iteration", and add "asymptotically".

### 14. [Minor] Literature attribution for the order-2 region

**Location:** Sec. 4.1 ("the region also holds under a weaker, non-stagnant assumption [Cleghorn2018]"); S-th-pso.

**Problem.** Please check whether Cleghorn and Engelbrecht (2018) state the region for general c1, c2 or only in terms of c1 + c2. If they state it generally in terms of c1 + c2, it conflicts with the (correct) claim here that c1 ≠ c2 is stricter under stagnation.

**Fix.** Say "for c1 = c2". Consider citing Jiang, Luo and Yang (Inf. Process. Lett. 2007) for the same region, and Liu (Evol. Comput. 2015) or Bonyadi and Michalewicz (IEEE TEC 2016) for c1 ≠ c2. The general condition in the paper is likely known, so present it as a restatement rather than a new result.

### 15. [Minor] Box-repaired shaking is not uniform on the ℓ∞ shell

**Location:** Sec. 4.3 ("draws a point uniformly from N_k(X) (box-repaired)"); S-th-vns.

**Problem.** Clipping moves probability mass onto the box faces and can put the point inside δ_{k−1}.

**Fix.** Say "uniform on the shell, then clipped (so neither uniform nor necessarily in the shell)".

### 16. [Minor] Rounding of T1 is Python's round-half-to-even

**Location:** Sec. 4.3 ("T1 = 100 (99.5 rounded up)").

**Problem.** Python's `round` rounds half to even. It happens to give 100 here; for ω = 0.25, 0.75 and 0.9 the values are 49, 150 and 180.

**Fix.** Write "rounded half to even (Python `round`)".

### 17. [Minor] Prop. S-frozen(c) and S-settings describe a regime the runs never reach

**Problem.** The limiting SD is about 1.04|P−G|, and the matched RMS distance of initial layouts to G is 291 m (Table D-feasstart). The resulting samples leave the box almost surely, so assumption (A2) fails immediately.

**Fix.** Add a sentence that (c) is illustrative only.

### 18. [Minor] Mixed statistics in Table D-vns

**Problem.** "First descent: evaluations (median) 2,860" sits next to "share 78.5%", but the share is a pooled ratio over runs.

**Fix.** Label it "pooled share", or give the median share.

### 19. [Minor] "Most evaluations are infeasible" does not distinguish the two settings

**Location:** Abstract.

**Problem.** The constriction setting also has 68.4% infeasible evaluations over iterations 1-200. What distinguishes the settings is feasible particle evaluations in the second half: 2.6% for the old setting against 48.4% for the constriction setting.

**Fix.** Say that directly.

---

## 2. Do the diagnostics support the stated mechanisms?

- **Old-setting spread and feasibility.** The descriptive differences are real and paired: spread 8.9×, |V| 13×, clipped 8.7% against 1.6%, feasible late evaluations 2.6% against 48.4%. The attribution to order-2 instability is not shown (issue 1), and the spacing-only argument is partly a tolerance and penalty-scaling artefact (issue 7).
- **Compass-search descent dominating Phase 2.** Well supported (97.5% of the gain from the first descent; unfinished in 39 of 40 runs with N ≥ 12). This finding undercuts the "VNS" framing (issue 2).
- **Feasible-start stall.** Well supported descriptively: 0 of 21 instrumented and 0 of 210 stored runs move. The relabelling test sensibly rules out the label-mismatch explanation. Proposition S-frozen and the lemma are correct and correctly marked as not implying the stall.
- **RS-VNS weakness.** Supported, but misattributed to the square rather than to packing density (issue 6). There is also an off-by-15 switch reading (issue 9).

## 3. Do the conclusions follow?

- **"PSO setting changes the ranking".** Yes, at 6,030 evaluations. The "because order-2" explanation is plausible but not established (issue 1).
- **"Salp-swarm phases add nothing beyond random sampling; PSO gives a clear gain".** Yes for this budget, split, local search and random starts. The PSO gain is plausibly a feasibility gain (issues 4 and 5). The Laplace-step wording contradicts the paper's own caveat (issue 5b).
- **"PSO-VNS ≈ PSO".** Holds for the 68-case average at a post-hoc margin. It does not hold for N ≥ 10, and it is specific to a truncated best-improvement compass search (issues 2 and 8). The prefix proposition makes this contrast clean, which is a genuine strength.

## 4. Overall verdict

The theory supplement is mathematically sound, carefully scoped ("what the conditions do not cover") and numerically exact; I found no errors in any proposition or typed constant. The weaknesses lie in how the theory and the diagnostics are used. The abstract and conclusions turn a stagnation, no-bounds stability result into a causal explanation of bounded, clipped behaviour, and generalise controls run at one budget, one split and one inefficient local search into statements about "VNS", "salp-swarm phases" and "the Laplace step". The comparison pool also omits the standard continuous black-box baseline (CMA-ES) and leaves the other baselines untuned, although the paper's own thesis is that configuration matters.

**Recommendation: major revision.** The theory section itself needs only minor revision.

## 5. Top 5 fixes

1. **Establish or soften the order-2 causality (issue 1).** Run a c-sweep across the order-2 boundary (w = 0.7, c from 1.5 to 2.0, boundary 1.749) and old-setting runs with velocity reset, reflection or Vmax on the six diagnostic cases. Otherwise replace "so"/"which violates" with "is consistent with", and state that the old setting ran without Vmax.
2. **Qualify the Phase-2 "VNS" and the PSO-VNS ≈ PSO claims (issues 2 and 8).** Report the cycles' evaluation and gain shares, say that Prop. S-ls is vacuous for N ≥ 12, and add one efficient local-search control (opportunistic polling or turbine-wise neighborhoods). Report equivalence for N ≥ 10, where the mean difference of −0.079 pp exceeds the margin.
3. **Scope the component-analysis conclusions (issues 4 and 5).** Add "at 6,030 evaluations, from random starts" to the salp-swarm and PSO-phase findings. Make the PSO gain explicitly a gain over infeasible random starts (VNS alone wins from feasible starts). Replace "the Laplace step is harmful" with "LX-SSA as published is worse", or run a cost-matched ablation.
4. **Correct the control and violation bookkeeping (issues 6, 7 and 9).** RS-VNS weakness is mainly spacing, not (π/4)^N. Reclassify "spacing-only" finals with the penalty-consistent threshold and report box-tangent pinning. Read the RS and RSD switch at 3,015, not 3,030.
5. **Strengthen or qualify the baseline pool (issue 3).** Add CMA-ES (pycma defaults, same penalty, bound transform) and a small DE sensitivity check, or state in the abstract that the ranking is among untuned default baselines without CMA-ES.

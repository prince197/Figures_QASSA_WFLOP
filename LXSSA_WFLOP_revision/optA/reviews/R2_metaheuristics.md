# Reviewer 2: metaheuristics, swarm intelligence and benchmarking

Manuscript: *Revisiting Metaheuristic Comparisons for Continuous Wind Farm Layout Optimization: Baseline Configuration, Component Analysis and a Two-Phase PSO-VNS Design* (MPCE_PSO_VNS.pdf, 12 pp.; sources optA/01-11; code in analysis/).

## Overall assessment

This revision is much better than a typical "new hybrid beats old baselines" paper. The authors report equal budgets that count every call (including gradient calls), seed-paired runs, feasibility-aware ranking, Holm-corrected tests, a random-sampling control and a candid limitations section. I checked the PSO description against the code and the theory, and it is correct:

- The inertia-weight form is equivalent to Clerc–Kennedy constriction: kappa = 0.72984 and c = kappa * 2.05 = 1.49618.
- Poli's bound 24(1 - w^2)/(7 - 5w) gives 3.497 ≈ 3.50 at w = 0.7 and 3.3475 ≈ 3.35 at w = 0.7298.
- The asynchronous gbest update, zero initial velocity, no velocity clamping and clipping match `authors_optimizers.PSO`.
- PSO-VNS really does reproduce PSO up to the switch. It uses the same class, the same seeded global RNG stream and T1 = 100, which gives 3,030 calls.

Three things still need work:

1. **The instability claim is imprecise.** The old setting is order-1 stable and order-2 unstable, and this only holds under stagnation with no bounds. In the code, clipping to the square interacts with that instability, and that interaction is probably what produces the infeasibility.
2. **The component-analysis wording overreaches.** "Add nothing beyond random sampling" rests on non-significant run-level tests. The paper's own case-mean test contradicts it: my recomputation from `analysis/mpce_ablation_tests.csv` gives SSA-VNS better than RS-VNS with Wilcoxon p = 4.5e-4. The RS-VNS control is also a weak one: uniform sampling in the bounding square.
3. **Several headline statements hold only at the single small budget B = 6,030.** At 120,030 calls PSO ranks behind SSA-VNS.

The literature on standard PSO, memetic and PSO+VNS hybrids, and benchmarking guidelines is thin. Recommendation: **major revision**. The core findings are probably robust, but several claims need to be stated precisely.

## Issues

1. **[major] 02_intro.tex / 04_methods.tex**
   - **Quote:** "the particles of PSO converge only if the inertia weight $w$ and the acceleration coefficients $c_1$, $c_2$ lie in a known region"
   - **Problem:** The statement does not say which kind of convergence is meant. Order-1 (mean) stability needs |w| < 1 and 0 < c1 + c2 < 4(1 + w). That bound is 6.8 at w = 0.7, so the old setting (c1 + c2 = 4) *is* order-1 stable: the expected position converges. Only order-2 (variance) stability, Poli's condition (9), is violated. The condition is also derived under stagnation (fixed P and G) and without bounds.
   - **Fix:** Replace with: "the expected particle positions converge (order-1 stability) if |w| < 1 and c1 + c2 < 4(1 + w), and their variance converges as well (order-2 stability) only if, in addition, (9) holds [Poli 2009; Cleghorn & Engelbrecht 2018]. The setting w = 0.7, c1 = c2 = 2 is order-1 stable but order-2 unstable." In Section III-A, add that (9) was derived under stagnation and has been shown to hold under a weaker non-stagnant assumption.
   - **Reference to add:** C. W. Cleghorn and A. P. Engelbrecht, "Particle swarm stability: a theoretical extension using the non-stagnate distribution assumption," *Swarm Intelligence*, vol. 12, no. 1, pp. 1–22, 2018, doi:10.1007/s11721-017-0141-x.

2. **[major] 04_methods.tex**
   - **Quote:** "The variance of the positions then grows even when the bests no longer change, and the swarm keeps sampling far from its best layouts"
   - **Problem:** In the implementation the velocity is unclamped and kept after clipping to [-r, r]. The position variance therefore cannot grow without bound. Instead, coordinates pile up on ±r, and the corners of the square lie outside the circular farm. This clipping/instability interaction is a plausible cause of the 86.0% feasibility and 8 "ranked-last" cases of the old setting. The mechanism is asserted but not measured.
   - **Fix:**
     - Report, for both settings, the share of coordinates clipped per iteration and the share of infeasible runs that violate the boundary versus the spacing constraint.
     - Replace the sentence with: "Under stagnation the variance of the sampling distribution diverges; with box clipping and unclamped velocities this means that particles oscillate with large amplitude and are frequently clipped to the square, whose corners lie outside the farm (Table S.x)."
     - State explicitly that velocities are kept after clipping. This boundary-handling choice is not neutral for PSO.

3. **[minor] 04_methods.tex**
   - **Quote:** "We use the standard global-best PSO~\cite{Kennedy1995} with unmodified update rules."
   - **Problem:** Two corrections are needed:
     - "Standard PSO" usually means Bratton & Kennedy's 2007 definition. That definition uses constriction but a *local ring* topology, so the gbest variant used here is canonical, not standard.
     - Kennedy & Eberhart 1995 has no inertia weight. The equivalence of the inertia and constriction forms, which the paper states, is the subject of Eberhart & Shi (2000) and should be cited there.
   - **Fix:** Write "the canonical global-best PSO [Kennedy1995] in the inertia-weight form [Shi1998] ... equivalent to the constriction form [Clerc2002; Eberhart & Shi 2000]; note that the standard PSO of [Bratton & Kennedy 2007] uses a ring topology."
   - **References to add:**
     - R. C. Eberhart and Y. Shi, "Comparing inertia weights and constriction factors in particle swarm optimization," in *Proc. 2000 Congr. Evol. Comput. (CEC)*, vol. 1, pp. 84–88, 2000.
     - D. Bratton and J. Kennedy, "Defining a standard for particle swarm optimization," in *Proc. IEEE Swarm Intell. Symp. (SIS)*, pp. 120–127, 2007, doi:10.1109/SIS.2007.368035.

4. **[minor] 04_methods.tex**
   - **Quote:** "the positions are clipped to the bounding square after each move (box repair)"
   - **Problem:** The update order is otherwise described correctly: asynchronous, with P and G updated after each evaluation, which matches the code. Two points are missing:
     - The velocity is not reset or reflected at the bound.
     - Clipping is to the square [-r, r]^2, not to the circular site, so every boundary violation must be removed by the penalty.
     Both choices affect PSO more than the population-replacement methods.
   - **Fix:** Add "the velocity is left unchanged when a coordinate is clipped (no reflection or damping)". Cite a boundary-handling study, or state this as a limitation.

5. **[major] 01_front.tex (abstract), 02_intro.tex, 10_limits_concl.tex**
   - **Quote:** "with the constriction setting PSO ranks \NPSOPos{} and outperforms all salp swarm methods"
   - **Problem:** This is shown only at B = 6,030 (200 PSO iterations for up to 30 variables). Table VI shows the reverse at larger budgets. At 120,030 calls PSO has average rank 5.83, SSA-VNS 4.17 and LX-SSA-VNS 5.83. The old setting was run only at 6,030 calls, so its behavior at larger budgets is unknown. A more explorative, order-2-unstable swarm may do relatively better when given more time.
   - **Fix:**
     - Add "at 6,030 evaluations" to every occurrence (abstract, contribution 2, Section VI-A and the conclusion).
     - Add a sentence in Section VIII-A: "at 120,030 evaluations PSO alone falls behind SSA-VNS on the six largest cases (average rank 5.83 vs. 4.17)."
     - Ideally, run the old setting on the six largest cases at 30,030 and 120,030 calls.

6. **[major] 06_results.tex**
   - **Quote:** "Changing the coefficients of one baseline thus reverses the order of PSO and the salp-swarm hybrids."
   - **Problem:** Only PSO's parameters were revisited, which creates a tuning asymmetry. The change is theory-driven, not tuned on the test cases, so this is not a fairness violation in itself. But the paper recommends that every baseline be configured inside its working range, and it does not apply that check to the other methods:
     - DE: F = 0.5, CR = 0.9.
     - VNS: delta_k, h0, h_min.
     - MS-SLSQP: the 1-m finite-difference step.
     - SSA has no parameters.

     For DE the check is easy and favorable: Zaharie's critical value is F_crit = sqrt((1 - CR/2)/N_p) ≈ 0.135 < 0.5, so DE does not collapse prematurely.
   - **Fix:**
     - State that the constriction coefficients were fixed a priori from the literature and not selected on these cases.
     - Report the Zaharie check for DE.
     - Add a small one-factor sensitivity study for DE (CR ∈ {0.1, 0.9}), VNS (delta_k ∈ {0.05k·r, 0.1k·r}) and the SLSQP step on the six largest cases. Alternatively, state explicitly that all non-PSO settings are untuned defaults and that the rankings are conditional on them, as recommended by LaTorre et al. (2021).
   - **References to add:**
     - D. Zaharie, "Critical values for the control parameters of differential evolution algorithms," in *Proc. MENDEL 2002, 8th Int. Conf. Soft Computing*, Brno, pp. 62–67, 2002.
     - A. LaTorre, D. Molina, E. Osaba, J. Poyatos, J. Del Ser, F. Herrera, "A prescription of methodological guidelines for comparing bio-inspired optimization algorithms," *Swarm Evol. Comput.*, vol. 67, 100973, 2021, doi:10.1016/j.swevo.2021.100973.

7. **[major] 07_ablation.tex**
   - **Quote:** "The salp-swarm phases do not: SSA-VNS and RS-VNS are statistically tied in \NAblSSAVNSvsRSVNST{} of the 68 cases"
   - **Problem:** A non-significant result is not evidence of equivalence. The run-level tests (30 runs, Holm over 8 contrasts) have low power for differences of about 0.07 pp. The paper's own case-mean statistic (Section V-B, item 2), recomputed from `mpce_ablation_tests.csv`, gives:
     - SSA-VNS has lower mean loss than RS-VNS in 42 cases and higher in 18. The two-sided Wilcoxon p is 4.5e-4 and ΔL = −0.068 pp.
     - SSA-VNS's average rank is also better (3.79 vs. 4.46).

     The statistics are also applied asymmetrically. For PSO-VNS vs. PSO the paper stresses the non-significant case-mean test (p = 0.30). For SSA-VNS vs. RS-VNS it reports only the run-level ties, and here the case-mean test *is* significant.

     The same applies to "Among the Phase-1 options tested, only PSO is effective."
   - **Fix:**
     - Report the case-mean Wilcoxon test for every contrast in Table IV.
     - Replace "statistically tied ... no better than random sampling" with: "SSA gives VNS a slightly better start than random sampling (ΔL = −0.07 pp; case-mean p = 4.5 × 10⁻⁴, significant in 2 of 68 cases at run level), an order of magnitude less than PSO (−0.40 pp)."
     - If "adds nothing" is to be kept, use an equivalence test (TOST) or a Bayesian signed-rank test with a region of practical equivalence, e.g. ±0.1 pp.
   - **Reference to add:** A. Benavoli, G. Corani, J. Demšar, M. Zaffalon, "Time for a change: a tutorial for comparing multiple classifiers through Bayesian analysis," *J. Mach. Learn. Res.*, vol. 18, no. 77, pp. 1–36, 2017.

8. **[major] 04_methods.tex / 05_setup.tex**
   - **Quote:** "pure random sampling of $\operatorname{round}(\omega B)$ layouts, uniform in the bounding square and including the initial population"
   - **Problem:** This control is weak by construction:
     - A uniform point in the square lies inside the circle with probability π/4, so all N turbines lie inside with probability (π/4)^N, about 0.03 for N = 15, before spacing is even considered.
     - Under the μ = 10¹⁰ penalty, "best of 3,015" is therefore the *least-violating* layout, not a low-wake one. This fits the 66.7% feasibility at the switch.
     - "PSO beats random sampling" thus largely means "a penalty-guided swarm finds feasibility faster than box sampling".
     - Minor: RS uses 3,015 Phase-1 calls, while the swarms use N_p(T1 + 1) = 3,030. The template says "the same number of calls".
   - **Fix:** Add at least one stronger control on the six largest cases:
     - (a) random sampling uniform in the disc (polar sampling);
     - (b) random multistart of the same compass search with an equal budget;
     - (c) random sampling with a cheap feasibility repair.

     Qualify the conclusion as "beyond uniform box sampling" and correct the call count. Relate the control to the random-search literature for WFLOP.
   - **Reference to add:** J. Feng and W. Z. Shen, "Solving the wind farm layout optimization problem using random search algorithm," *Renew. Energy*, vol. 78, pp. 182–192, 2015, doi:10.1016/j.renene.2015.01.005.

9. **[major] 07_ablation.tex / 04_methods.tex**
   - **Quote:** "whereas on larger layouts the single-turbine moves of the compass search pay off"
   - **Problem:** The Phase-2 "VNS" gets 3,000 calls.
     - One best-improvement sweep of the compass search costs 2n calls: 60 for N = 15.
     - h must be halved at least 6 times, from 0.05r to below 10⁻³r.
     - So the initial descent from G alone costs at least about 360 calls, and typically much more.
     - The largest shaking neighborhood (δ5 = 0.5r in ℓ∞ on all 2N coordinates) is close to a restart.

     For large N, Phase 2 may therefore perform only a few shaking steps. In that case the component analysis measures "PSO + compass-search refinement", not "PSO + VNS". The paper's own explanation points to the local search, not the neighborhood structure.
   - **Fix:**
     - Report per run the number of shaking steps, the number of improvements and the maximum k reached, by N.
     - Add a "PSO + local search only (restart the compass search from G with smaller h)" variant to Table IV.
     - If shaking rarely happens, rename the claim ("a local-search phase helps every swarm").

10. **[major] 06_results.tex**
    - **Quote:** "feasible in every run, at least as good as a correctly configured PSO and better for $N\ge10$"
    - **Problem:** Table III shows PSO-VNS significantly *worse* than PSO in 4 cases (14/50/4). For N ≥ 10 it wins 9 of 20 cases and loses 2, and the case-mean difference is not significant (p = 0.30). Neither "at least as good" nor "better for N ≥ 10" is supported.
    - **Fix:** Replace with: "feasible in every run and not significantly different from a correctly configured PSO over all cases (p = 0.30); significantly better in 14 and worse in 4 cases, with the wins concentrated at N ≥ 10 (9 wins, 2 losses in 20 cases)." Make the same change in the conclusion ("is better than PSO for N ≥ 10").

11. **[minor] 07_ablation.tex**
    - **Quote:** "Adding VNS helps every swarm: each hybrid ranks ahead of the same swarm run alone for the full budget"
    - **Problem:**
      - For PSO the contrast is clean. The code confirms that Phase 1 reproduces the first 3,030 calls of the PSO run.
      - For SSA and LX-SSA it is not clean. ξ1 = 2·exp(−(4t/T)²) uses T = T1 in the hybrid, so Phase 1 is a *compressed* SSA with a different schedule, not the first half of the stand-alone run. The SSA-VNS vs. SSA contrast therefore confounds the VNS phase with the schedule.
      - The PSO contrast itself is only 8/58/2 (Table IV).
    - **Fix:** Add "for PSO this contrast isolates the second phase; for SSA and LX-SSA it also changes the ξ1 schedule". Optionally add an SSA-with-T1-schedule-alone control. Tone the sentence down to "the VNS phase improves on SSA and LX-SSA and is not worse than continuing PSO".

12. **[minor] 04_methods.tex**
    - **Quote:** "Comparing SSA with LX-SSA isolates the Laplace step."
    - **Problem:** In `authors_optimizers.LXSSA`, each follower evaluates both the candidate and the midpoint follower. Then the whole population is re-evaluated, so the chosen follower is evaluated twice. That makes N_p/2 redundant calls per iteration, a cost of 2N_p per iteration and half as many iterations as SSA (100 vs. 200). The comparison therefore confounds the Laplace step with the greedy selection, the halved iteration count and the wasted evaluations.
    - **Fix:** Write "compares LX-SSA as published (including its 2N_p calls per iteration) with SSA". Either remove the redundant re-evaluation in a variant, or state this as a limitation of the Laplace-step conclusion.

13. **[minor] 07_ablation.tex**
    - **Quote:** "The Laplace step of LX-SSA gives no gain"
    - **Problem:** This understates the result:
      - Alone, LX-SSA is significantly *worse* than SSA in 7 cases (0/61/7, +0.166 pp).
      - In the hybrid, SSA-VNS beats LX-SSA-VNS on the case means (my recomputation: 47 of 60 non-zero cases lower, Wilcoxon p = 5.4e-6).
    - **Fix:** "The Laplace step, as implemented, is harmful: LX-SSA is significantly worse than SSA in 7 cases and never better, and LX-SSA-VNS is worse than SSA-VNS on the case means (p = ...)". Link the explanation to issue 12.

14. **[major] 02_intro.tex / 01_front.tex / 10_limits_concl.tex**
    - **Quote:** "the Laplace step of the Laplacian SSA (LX-SSA)~\cite{Solanki2023}, proposed earlier by two of the authors, add nothing beyond random sampling"
    - **Problem:**
      - Beyond issue 7, the conclusion is stated generically, but it holds for these implementations, this budget, penalty-based constraint handling and random starts.
      - It is also not placed in the literature. Castelli et al. (2022) showed conceptual and mathematical flaws in SSA. In particular, the leader update is a uniform perturbation of the food source with a shrinking radius, and followers only average. This explains *why* SSA behaves like random search around the incumbent, and it strengthens the paper's message.
    - **Fix:** Qualify: "on this benchmark and at 6,030 evaluations, SSA phases give at most a marginal gain over random sampling". Cite and discuss the review.
    - **Reference to add:** M. Castelli, L. Manzoni, L. Mariot, M. S. Nobile, A. Tangherloni, "Salp Swarm Optimization: a critical review," *Expert Syst. Appl.*, vol. 189, 116029, 2022, doi:10.1016/j.eswa.2021.116029.

15. **[major] 08_beyond.tex**
    - **Quote:** "PSO and DE never improve on the best initial layout in these runs"
    - **Problem:** This is a striking result: in all 210 feasible-start runs, the swarm's global best never moves. PSO's search dynamics therefore contribute nothing once feasibility is given. This implies that PSO's advantage as a Phase-1 generator from random starts comes mainly from penalty-guided feasibility attainment, not from better exploration of the wake-loss landscape. The component analysis cannot separate the two, yet the paper attributes it to a "convergent PSO" that "contracts around a good region".
    - **Fix:**
      - Decompose the Phase-1 advantage: compare the switch points *conditional on feasibility*, the wake loss of feasible switch points, and the calls to first feasibility for PSO, SSA and RS.
      - State in Sections VII and X: "our evidence shows that a convergent PSO is an effective feasibility-seeking start generator from random starts; it does not show better exploration of the wake-loss landscape".
      - Mention label symmetry (permutation invariance of turbines) as a structural reason why recombination-type moves fail.

16. **[minor] 04_methods.tex**
    - **Quote:** "the swarm contracts around a good region within its share of the budget"
    - **Problem:** Order-2 stability guarantees that particles contract around (φ1·P + φ2·G)/(φ1 + φ2). It does not guarantee that this region is good or locally optimal: the canonical gbest PSO is not a guaranteed local optimizer. The sentence presents an empirical observation as a consequence of (9).
    - **Fix:** "with coefficients that satisfy (9) the particles contract around their personal and global bests instead of oscillating across the site; empirically (Section VII) this gives VNS a better start".

17. **[minor] 02_intro.tex**
    - **Quote:** "under a controlled protocol a simple two-phase design, a convergent PSO followed by the basic VNS~\cite{Mladenovic1997,Hansen2001} (PSO-VNS), is the most consistent performer"
    - **Problem:** The pool has no state-of-the-art continuous black-box optimizer: no CMA-ES/IPOP-CMA-ES and no L-SHADE. Most members are weak (SSA variants, untuned DE with 71% feasibility). "Most consistent performer" is relative to this pool.
    - **Fix:**
      - Add "among the eight methods compared".
      - Preferably add CMA-ES with restarts on at least the six largest cases.
      - Cite the WFLOP comparison of eight optimizers by Thomas et al. (2023), which found that initialization and multistarts matter as much as the algorithm.
    - **References to add:**
      - N. Hansen and A. Ostermeier, "Completely derandomized self-adaptation in evolution strategies," *Evol. Comput.*, vol. 9, no. 2, pp. 159–195, 2001, doi:10.1162/106365601750190398.
      - J. J. Thomas, N. F. Baker, P. Malisani, E. Quaeghebeur, S. Sanchez Perez-Moreno, J. Jasa, C. Bay, F. Tilli, D. Bieniek, N. Robinson, A. P. J. Stanley, W. Holt, A. Ning, "A comparison of eight optimization methods applied to a wind farm layout optimization problem," *Wind Energ. Sci.*, vol. 8, pp. 865–891, 2023, doi:10.5194/wes-8-865-2023.

18. **[minor] 05_setup.tex**
    - **Quote:** "SciPy's SLSQP~\cite{Kraft1988} with explicit constraints and forward-difference gradients (1-m step)"
    - **Problem:** The top-hat Jensen objective is discontinuous (Section III-C says so). Forward differences with a 1-m step are therefore either zero or dominated by jumps, so MS-SLSQP is handicapped by design. The IEA37 winner (Table VIII footnote) used wake expansion continuation precisely to get around this.
    - **Fix:** State this handicap where MS-SLSQP results are discussed, and do not read its low rank as evidence against gradient methods. Cite the continuation approach.
    - **Reference to add:** J. J. Thomas and A. Ning, "A method for reducing multi-modality in the wind farm layout optimization problem," *J. Phys.: Conf. Ser.*, vol. 1037, 042012, 2018, doi:10.1088/1742-6596/1037/4/042012.

19. **[minor] 04_methods.tex**
    - **Quote:** "PSO-VNS is not a new PSO variant: Phase~1 is the plain PSO of Section~\ref{sec:pso}"
    - **Problem:** PSO followed by, or embedded with, VNS or local search is an established memetic design, and even the name "PSO-VNS" has been used before. The paper cites no memetic or hybrid-metaheuristic literature. That makes the design look more novel than it is and misses results on when to switch and how to split the budget.
    - **Fix:** Add a short paragraph in the introduction or Section IV positioning PSO-VNS as a sequential (relay) memetic hybrid. Cite continuous VNS variants for the Phase-2 design choices.
    - **References to add:**
      - Y. G. Petalas, K. E. Parsopoulos, M. N. Vrahatis, "Memetic particle swarm optimization," *Ann. Oper. Res.*, vol. 156, no. 1, pp. 99–127, 2007.
      - M. F. Tasgetiren, Y.-C. Liang, M. Sevkli, G. Gencyilmaz, "A particle swarm optimization algorithm for makespan and total flowtime minimization in the permutation flowshop sequencing problem," *Eur. J. Oper. Res.*, vol. 177, pp. 1930–1947, 2007.
      - E. Carrizosa, M. Dražić, Z. Dražić, N. Mladenović, "Gaussian variable neighborhood search for continuous optimization," *Comput. Oper. Res.*, vol. 39, no. 9, pp. 2206–2213, 2012.

20. **[minor] 02_intro.tex**
    - **Quote:** "The benchmarking literature has long warned that such comparisons are hard to interpret~\cite{Sorensen2015,BartzBeielstein2020}"
    - **Problem:** The seven recommendations in Section X largely restate existing guidelines (equal budgets, parameter reporting, paired statistics, ablation against controls). They should be credited so that the contribution is framed as *applying and quantifying* them for WFLOP rather than proposing them.
    - **Fix:** Add LaTorre et al. (2021) [issue 6] and Benavoli et al. (2017) [issue 7] here and in Section X. Write "we follow and quantify for WFLOP the guidelines of [...]".

21. **[minor] 06_results.tex**
    - **Quote:** "a PSO baseline outside the convergence region understates what PSO achieves on the WFLOP"
    - **Problem:** This generalizes from one budget, one benchmark, one boundary-handling rule and one alternative setting. The order-2-unstable setting could be competitive at larger budgets or with reflection or velocity damping at the bounds (issue 2).
    - **Fix:** "on this benchmark, at 6,030 evaluations and with clipping, a PSO baseline outside the order-2 stability region performed markedly worse than the constriction setting, so a method compared with it can appear stronger than it is".

22. **[minor] 07_ablation.tex**
    - **Quote:** "In line with the component analysis, most of the benefit comes from the PSO phase, and VNS serves as a final refinement."
    - **Problem:** The split study is monotone in ω: 0.25 < 0.5 < 0.75. ω = 1 is exactly PSO, and PSO-VNS vs. PSO at ω = 0.5 is not significant. So the optimum may lie between 0.75 and 1, and the VNS contribution may be negligible. The recommendation ω = 0.75 is an extrapolation at the edge of the tested range.
    - **Fix:** Add ω ∈ {0.9, 1.0 (= PSO)} to Table V, which is cheap (12 cases). Report whether the ω = 0.75 hybrid is significantly better than PSO alone. Adjust the recommendation in Section X accordingly.

23. **[minor] 07_ablation.tex / 06_results.tex**
    - **Quote:** "under the Holm adjustment over the contrasts of Table~\ref{tab:ablation}"
    - **Problem:** The same pair (PSO-VNS vs. PSO, same runs) is reported as 14/50/4 in Table III and 8/58/2 in Table IV. The difference comes from different Holm families (7 vs. 8 comparisons per case), but readers will see an inconsistency.
    - **Fix:** Add a footnote to Table IV: "W/T/L differ from Table III because the Holm family differs (8 planned contrasts instead of 7 comparisons against PSO-VNS); the unadjusted counts are ..."

24. **[minor] 05_setup.tex**
    - **Quote:** "Every run of the main comparison and the component analysis has exactly $B=6{,}030$ objective calls"
    - **Problem:** All 68-case conclusions use one small budget: 200 PSO iterations for up to 30 variables. The budget study (6 cases) shows large rank changes: PSO 1.83 → 5.83, VNS 6.50 → 3.17, DE 10.00 → 3.33. Fixed-budget ranks at one budget are a single cut through the anytime behavior.
    - **Fix:**
      - Report anytime results for all 68 cases from the existing convergence logs, e.g. ECDFs of the target-hitting calls or average ranks at 1,000, 3,000 and 6,030 calls.
      - State in the abstract that the 68-case results are at 6,030 evaluations.

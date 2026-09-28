# R5 - Hostile claims audit (prose vs. data)

Scope: 01_front (abstract), 02_intro, 06_results, 07_ablation, 08_beyond, 09_robust, 10_limits_concl.
Data: analysis/mpce_summary.json (generated 2026-09-28 10:31:27), mpce_numbers.tex, mpce_tab_*.tex,
mpce_case_stats.csv, mpce_case_tests.csv, mpce_ablation_tests.csv, mpce_budget_case_stats.csv, raw run CSVs
(mpce_rsvns_s0of1.csv, fresh_bgrid.csv, fresh_vgrid.csv). mpce_check_final.py: 30 PASS / 0 FAIL / 0 PENDING
(C07 is defined but no longer referenced anywhere).

The case-mean statistics I computed below use the paper's own rule: the mean wake loss of the feasible runs,
qualified only with at least 15 feasible runs, and a two-sided Wilcoxon test on the 68 case means (Sec. V, test 2).

Format: [category] file - quote - evidence - fix.

## A. False, or contradicted by the data

1. [false] 07_ablation / 01_front / 02_intro / 10_limits_concl - "The salp-swarm phases do not ... SSA and LX-SSA are no better than random sampling"; abstract "salp swarm phases ... add nothing"; conclusion "the SSA and LX-SSA phases are statistically tied with random sampling".
   Evidence: these claims hold for LX-SSA but not for SSA. SSA-VNS vs RS-VNS on the 68 case means: SSA-VNS has the lower loss in 42 cases, RS-VNS in 18, mean -0.068 pp, and the Wilcoxon test gives p = 4.5e-4. This is the same case-mean test the paper uses to call PSO-VNS vs PSO "not significant" (p = 0.30). Ablation average ranks: SSA-VNS 3.79 vs RS-VNS 4.46. Only the per-case run-level test is tied (2/66/0). For LX-SSA-VNS vs RS-VNS: p = 0.72, LX-SSA-VNS lower in 27 cases vs 33, so "tied" holds there.
   Condition C15 checks only the run-level W/T/L, so it passes although the case-mean test contradicts the claim.
   Fix (07): "The salp-swarm phases add little: SSA-VNS and RS-VNS are tied in \NAblSSAVNSvsRSVNST{} of the 68 cases at run level (\NAblSSAVNSvsRSVNS), and SSA-VNS's small mean advantage (\NAblSSAVNSvsRSVNSDL{} pp) is far below that of PSO-VNS (\NAblPSOVNSvsRSVNSDL{} pp); LX-SSA-VNS is not better than RS-VNS (case-mean p = \NAblLXBVvsRSVNSP)." In the abstract, replace "salp swarm phases and the Laplace step add nothing" with "salp swarm phases add little beyond random sampling and the Laplace step adds nothing". Use the same wording in the intro (contribution 3) and the conclusion.
   New condition: add ablation.case_mean_wilcoxon for SSABV-RSVNS and LXBV-RSVNS to the summary. Wording "tied/no better" requires p >= 0.05 for both; otherwise the generated phrase must say "slightly better".

2. [false] 08_beyond (Horns Rev) - "This mean is below the installed layout, which only the best run (\NHRBestPSOVNS~GWh/yr) and \NHRAbovePSOVNS{} of the 30 runs exceed".
   Evidence: runs_above_installed = 6. "Only the best run" contradicts "6 of the 30 runs".
   Fix: "This mean is below the installed layout, which \NHRAbovePSOVNS{} of the 30 runs exceed (best run \NHRBestPSOVNS~GWh/yr)".
   New condition C31: hr16.methods.PSOBV.mean < hr16.installed_aep and runs_above_installed >= 1. The sentence depends on it and is currently covered only by a "verified" comment.

3. [false / overstated] 08_beyond (cost) - "At equal numbers of evaluations, all methods need a similar time (\NMsPerCallMin--\NMsPerCallMax~ms per evaluation".
   Evidence: cost.ms_per_call ranges from 0.24 ms (PSO) to 1.82 ms (MS-SLSQP), a factor of 7.5. The seven metaheuristics range from 0.24 to 0.47 ms. MS-SLSQP runs take 2.6-21.7 s, against 1.0-4.3 s for the others.
   Fix: "the metaheuristics need a similar time (0.24--0.47 ms per evaluation), MS-SLSQP about \NMsPerCallMax~ms because of its QP subproblems (...)". Generate the metaheuristic range as new macros (e.g. \NMsPerCallMetaMax).

## B. Overstated

4. [overstated] 06_results / 10_limits_concl / 02_intro - "at least as good as a correctly configured PSO and better for $N\ge10$"; conclusion "is better than PSO for $N\ge10$"; discussion "appears mainly for $N\ge10$".
   Evidence (mpce_case_stats.csv, mpce_case_tests.csv): for N >= 10, the Wilcoxon test on the 20 case means gives p = 0.048, unadjusted and borderline. The whole gain comes from Data Set II: all 10 DS II cases with N >= 10 favour PSO-VNS (-0.006 to -0.41 pp), with 8 significant wins at run level. In Data Set I with N >= 10, PSO has the lower mean loss in 6 of 10 cases (e.g. 1-750-11 +0.18, 1-750-12 +0.18, 1-500-10 +0.17 pp). There PSO-VNS is significantly better in 1 case and significantly worse in 2 (1-750-11, 1-750-12). "At least as good": PSO-VNS is significantly worse in 4 cases (1-500-8, 1-750-11, 1-750-12, 2-500-4).
   Fix (06): "... feasible in every run, not worse than a correctly configured PSO on aggregate (\NWtlPSO), and better for $N\ge10$ in Data Set~II, ...". Conclusion: "is better than PSO for the larger layouts of Data Set~II". Discussion: "appears mainly for $N\ge10$ with the broader wind rose (Data Set~II)".
   New condition C32: by_n.PSOC per data set. The claim "for N >= 10" requires W > L in both data sets. Otherwise the text must name the data set, which today means W_DSI(N>=10) = 1 and L = 2.

5. [overstated] 07_ablation - "switching to VNS instead of continuing the swarm gives a small gain that grows with $N$".
   Evidence: mean dL (PSO-VNS - PSO) by N is +0.005 to +0.025 pp for N = 4-9 (PSO-VNS slightly worse) and -0.04 to -0.20 pp only for N >= 10. The Spearman correlation between N and dL is -0.25 (p = 0.038). The gain comes from DS II (see 4).
   Fix: "gives a small gain that appears only for $N\ge10$, mainly in Data Set~II".
   Covered by C32 as proposed.

6. [overstated] 01_front / 02_intro / 07_ablation / 10_limits_concl - "a variable neighborhood search (VNS) phase helps every swarm"; conclusion "the VNS phase improves every swarm".
   Evidence: this holds for SSA (avg rank 6.56 to 3.79, loss 1.805 to 1.653 %) and LX-SSA (7.29 to 4.62). For PSO, however, the effect is 8/58/2, not significant on the case means (p = 0.30) or on the average ranks (p_Holm = 0.58 in the ablation, 0.71 in the main comparison). The paper itself calls it "not significant". No condition checks the SSA and LX-SSA part (C13 only checks PSO W > L).
   Fix (abstract): "a VNS phase clearly improves the salp swarm methods and gives PSO a small, non-significant gain". Same in the conclusion.
   New condition C33: ablation.friedman.avg_rank PSOBV < PSOC, SSABV < SSA, LXBV < LXSSA.

7. [overstated] 01_front - "PSO-VNS ... is feasible in every run and exceeds the energy yield of the installed Horns Rev~1 block".
   Evidence: at 6,030 calls the mean is 138.73 against 139.51 installed, and only 6/30 runs exceed it. The mean exceeds it only at 30,030 calls (139.62). C19 checks only the best run.
   Fix: "... is feasible in every run, and its best layout exceeds the energy yield of the installed Horns Rev~1 block". This matches the intro and the conclusion.

8. [overstated] 06_results, 02_intro, 10_limits_concl, 01_front context - "with up to 20 times the budget, \NBudgetPhrase"; "keeps the best rank at all three budgets"; 08 "Only the first place of PSO-VNS is stable across budgets."
   Evidence: the budget study covers only the six largest cases, and there is no post hoc or pairwise test. At 120,030 calls PSO-VNS ranks first in only 3 of 6 cases (DE first in 2, MS-SLSQP in 1). In 1-500-10 it ranks 5th. Its mean loss (3.14 %) is 0.03 pp below VNS (3.17 %) and 0.08 pp below SSA-VNS (3.23 %). On Horns Rev at 120,030 calls, VNS, RS-VNS and DE are better.
   Fix (06): "with up to 20 times the budget, it keeps the best average rank on the six largest cases, although at 120,030 calls its lead is small and untested". Intro/conclusion: "keeps the best average rank on the six largest cases at all three budgets". 08: "Only PSO-VNS keeps its (first) place in average rank at all three budgets; at 120,030 calls it ranks first in only \NBudFirstCountOneTwentyK{} of the six cases".
   New condition C34: add feasbudget.first_count per budget, generate the number, and require >= 3 for the wording "first place".

9. [overstated] 10_limits_concl (discussion) - "which method appears best depends as much on the comparison protocol as on the methods".
   Evidence: nothing quantifies "as much as". The data show a reversal for one baseline setting and for feasible starts, not an equal effect size.
   Fix: "which method appears best can depend on the comparison protocol as strongly as on the methods themselves".

10. [overstated] 02_intro (contribution 2) - "Evidence that the PSO baseline setting reverses the ranking of the methods".
    Evidence: only the order of PSO relative to the salp-swarm methods reverses. In both pools SSA-VNS stays ahead of LX-SSA-VNS, SSA stays ahead of LX-SSA, and DE stays last (Table baseline).
    Fix: "... reverses the order of PSO and the salp-swarm methods".

11. [overstated] 01_front - "a gap to the published IEA Wind Task~37 optima remains".
    Evidence: the reference is the best feasible published layout (par4), not a proven optimum, and an infeasible published layout is higher.
    Fix: "a gap to the best published IEA Wind Task~37 layouts remains".

12. [overstated/ambiguous] 07_ablation (split) - "$\omega=0.25$ is significantly worse than the default in \NSplitTwentyFiveSig{} cases and better in \NSplitTwentyFiveBetter{} case".
    Evidence: \NSplitTwentyFiveBetter is the count of significant wins (wtl_50_vs_25.L = 0). But ω = 0.25 has the lower mean loss in 1 of 12 cases (split.n_lower_loss_than_50.PSOBV25 = 1). The parallel ω = 0.75 clause uses unqualified "lower" to mean the raw mean, so readers will read "better in no case" as raw means, which is false.
    Fix: "... and significantly better in \NSplitTwentyFiveBetter{} case".

## C. Unsupported, or hand-typed results outside the pipeline

13. [unsupported] 02_intro (contribution 3), 10_limits_concl - "The SSA and LX-SSA phases are statistically tied with random sampling (RS-VNS)".
    Evidence: there is no LX-SSA-VNS vs RS-VNS contrast in Table ablation or in the summary. By my computation it is tied (case-mean p = 0.72; run level about 1/67/0 with an approximate Bonferroni correction), but the paper does not report it. For SSA see finding 1.
    Fix: add the contrast LXBV-RSVNS to ablation.contrasts and Table ablation, with a condition (T >= 60, W and L <= 5).

14. [uncovered] 07_ablation - "and LX-SSA-VNS ranks even behind RS-VNS on average".
    Evidence: 4.62 vs 4.46 (true, but close). No condition checks it.
    New condition C35: ablation.friedman.avg_rank.LXBV > avg_rank.RSVNS.

15. [uncovered] 07_ablation - "The Friedman test rejects equal average ranks". "Rejects" is hand-typed.
    Evidence: p = 1.1e-66 (true).
    Fix: use a generated verb, \NAblFriedVerb, as for the main test.
    New condition C36: ablation.friedman.p < 0.05.

16. [uncovered] 08_beyond (budget) - "PSO falls back from an average rank of ... to ..., VNS (...) and DE move up, and at 120,030 evaluations VNS and SSA-VNS come close to PSO-VNS".
    Evidence: PSO 1.83 to 5.83, VNS 6.50 to 3.17, DE 10.00 to 3.33 (true). Covered only by a "verified" comment. "Come close" is also true of RS-VNS (3.25), DE (3.31) and LX-SSA-VNS (3.31), so singling out VNS and SSA-VNS is selective.
    New condition C37: feasbudget.rank PSOC[120030] > PSOC[6030] + 2, and BVNS and DE [120030] < [6030] - 2. Also require loss[120030] - loss_PSOBV[120030] < 0.15 pp for every method named as "close".
    Fix: "VNS, SSA-VNS, RS-VNS and DE come within 0.2 percentage points of PSO-VNS" (generated).

17. [uncovered] 08_beyond (feasible starts) - "VNS alone is then best (\NFbFeasBestLoss{} has the lowest mean wake loss ...)" and "PSO and DE never improve on the best initial layout in these runs"; limitations "PSO and DE do not move from the best initial layout".
    Evidence: feasible-init losses are BVNS 3.810 and PSOBV 3.938 (true). PSOC and DE are both 5.5488 %, identical, which is consistent with no movement. Covered only by a "verified" comment. C28 checks the rank only. The phrase "VNS alone is then best" is hand-typed next to the macros and would silently conflict if the macros changed.
    New condition C38: argmin feasbudget.feasible.*.loss == BVNS, and |loss PSOC - loss DE| < 1e-9 (both equal to the initial best). Also add a summary field counting feasible-init PSO/DE runs whose final objective differs from the initial best, and require 0.
    "almost always make turbines overlap" is a mechanism that nobody measured. Either log the rejection rate or soften to "typically".

18. [uncovered] 08_beyond (Horns Rev) - "at 120,030 evaluations the loss is \NHRLossOneTwentyKPSOVNS\%, but VNS, RS-VNS and DE are lower".
    Evidence: 6.12, 6.24 and 6.34 < 6.45 (true; SSA-VNS 6.48 is not lower). Covered only by a comment. Based on 10 seeds with no test.
    New condition C39: hr16.loss_by_setting BVNS/RSVNS/DE [120030R].loss < PSOBV [120030R].loss. Generate the list of methods below PSO-VNS as a phrase macro, and add "(10 seeds)".

19. [uncovered] 08_beyond (IEA37) - "None of our methods reaches the best feasible published layouts"; "The gap is larger for 36 turbines, where our best layout ... comes from the gradient-based \NIEABestOfThirtySixThirtyK"; 06 "the best layouts remain below the best published ones"; conclusion "on IEA37 it stays below the best published layouts".
    Evidence: our best 16T is 413,332 < 418,924 and our best 36T is 833,365 (MS-SLSQP) < 863,676 (true). The gap is -1.33 % vs -4.20 % (true). No condition covers these statements. "Gradient-based" is hand-typed and becomes false if the macro changes to a metaheuristic.
    New condition C40: iea37.*.published.our_best < best_feasible_published for 16 and 36 turbines; |gap36| > |gap16|; and NIEABestOfThirtySixThirtyK == SLSQP. Alternatively, generate the adjective.

20. [unsupported, hand-typed results] 09_robust - "attains higher mean benchmark objectives than the original SSA in all six representative cases at 3,030 evaluations (Mann--Whitney test, $p<0.01$ in every case)".
    Evidence: the numbers come from ../selected_30_run_data.csv (archived), which is outside mpce_summary.json and the macro pipeline. The only support is a comment ("largest p = 0.0035"). "All six" and "p < 0.01" are hand-typed result numbers.
    Fix: compute them in mpce_results.py (summary key boundary_rule), generate \NBoundaryCases and \NBoundaryMaxP, and add condition C41: all cases have a higher mean and max p < 0.01.

21. [uncovered] 06_results (quality) - "The wake loss grows with $N$, falls as the farm becomes larger, and is higher for the broader wind rose of Data Set~II".
    Evidence (PSO-VNS case means): the loss is monotone in N apart from 0/-0 ties at N = 2 and monotone in r. DS II > DS I in 31 of 34 pairs, with 3 ties at zero loss. The claim is true, but no condition covers it.
    New condition C42: these three monotonicity checks on the PSO-VNS case means, with ties allowed.

22. [uncovered] 10_limits_concl (discussion) - "\NBaseNewBest{} ranks first in that pool, and PSO is ahead of every salp-swarm method in both".
    Evidence: baseline.constriction.best = PSO (true). C27 covers "ahead", but the grammar assumes that BaseNewBest is PSO and that NBaseOldBest is a salp-swarm hybrid.
    New condition: extend C27 so that baseline.old.best is in {SSABV, LXBV} and baseline.constriction.best == PSO.

23. [macro hygiene] 02_intro (contribution 4) - "VNS alone ranks first with feasible starts, PSO ranks first under a Gaussian wake model".
    Evidence: the values are correct (C22, C28), but the method names are hand-typed. The abstract, 06 and 10 use \NFbFeasBestRank and \NBestGauss.
    Fix: "\NFbFeasBestRank{} alone ranks first with feasible starts, \NBestGauss{} ranks first under a Gaussian wake model". The same applies to 08 "VNS alone is then best" (see 17).

24. [scope] 01_front / 02_intro - "30 seed-paired runs" and "start from the same initial population".
    Evidence: the summary has options.common_seeds = False, so the pipeline never checks that all methods share seeds 1-30 per case. Sec. V asserts the pairing but no check backs it. I could not verify identical initial populations across old (fresh_*) and new (mpce_*) files within the time budget: the first Curve value is NaN for infeasible starts.
    New condition C43: per case, budget and init, the set of seeds is identical for all methods. Add as well a check that the initial best objective or initial-population hash is equal across methods for a sample of seeds.

## D. Checked and supported (no action)

Old PSO: fifth, 86.0 %, 56 worse / 0 better, 8 cases below half feasible (C12). Baseline reversal (C27). Friedman and post hoc (C01-C03). PSO W/T/L 14/50/4, p = 0.30 (C04). SSA-VNS never worse (C06). Feasibility claims (C08, C09, C24). Phase-2 shares (C10, C11). PSO start vs VNS/RS-VNS, and PSO vs SSA/LX-SSA as Phase 1 (C14, C16). Laplace step (C17; LX-SSA is in fact worse in 7 cases). Split (C18). HR16 mean, p-values and feasibility (C19, C20, C25, C30). Cubic and Gaussian robustness (C21, C22). Spacing (C23). Budget-best (C26, C29). Feasible-start rank (C28). All other numbers are macros. Only design constants are typed by hand (68, 30 seeds, 6,030/3,030/3,015 calls, 7 d.f., 15-run threshold); the exception is item 20.

Housekeeping: C07 ("clearly outperforms SSA, LX-SSA, DE, VNS, MS-SLSQP") is defined but no longer referenced. Re-attach it to the 06 sentence "Against SSA, LX-SSA and DE, it is significantly better in ...", or delete it.

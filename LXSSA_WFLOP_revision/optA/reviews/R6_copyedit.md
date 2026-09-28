# R6 — Copy edit (IEEE/MPCE style) of MPCE_PSO_VNS.pdf (12 pp.) and optA/01–11

Scope: grammar, clarity, duplication across sections, terminology, cross-references, IEEE conventions, Nomenclature, abstract, headings, references, flow. Quotes are from the .tex sources (macros shown as in source) unless marked [PDF].
Abstract length: about 238 words in the compiled PDF (target ≤ 250: OK, but only a small margin; see #3).

Format: `N. [severity] file — "quote" — problem — replacement`

---

## A. Duplicated sentences and explanations (grouped)

1. [major] 02_intro / 04_methods / 06_results / 10_limits_concl — "This is a property of the setting, not of any particular study" (06); "We report this as a property of the setting, not as a finding about any particular study" (02); "This is a property of the setting, not of PSO" (04); "is a property of the PSO parameters, not of any particular study" (10) — the same disclaimer appears four times. In 10 it is also illogical, because a setting cannot be "a property of the PSO parameters". — Keep it once, in 06 §VI-A. In 02, delete the sentence. In 04, end the paragraph with "Section~\ref{sec:psosetting} shows that this changes the ranking of the compared methods." In 10, write "The old setting violates \eqref{eq:poli}; with it, PSO ranks …".

2. [major] 02_intro / 04_methods / 06_results — "common in PSO baselines for the WFLOP and used in the earlier version of this study" — this appositive appears almost verbatim three times (Intro ¶3, §III-A, §VI-A), and each time with "violates the convergence condition". — In 06, open §VI-A with: "The old setting (Section~\ref{sec:pso}) violates~\eqref{eq:poli}. We ran PSO with it on the same cases with the same budget and seeds." In 04, keep the full statement. In 02, keep the short mention.

3. [major] 01_front / 02_intro / 06_results / 10_limits_concl — "it is not significantly better than PSO over all cases …, VNS alone ranks first … PSO ranks first under a Gaussian wake model" — the four-part caveat (n.s. vs. PSO, VNS best from feasible starts, Gaussian reversal, IEA37 gap) appears five times: Abstract, Contribution 4, §VI-B "These results support…", §XI, §XII. — Keep it in the Abstract, §VI-B and §XII. Cut Contribution 4 after "…exceeds the energy production of the installed one" and add "; its limits are analyzed in Sections~\ref{sec:beyond} and~\ref{sec:robustness}." §XI already lists the Gaussian and feasible-start limits, so there delete "and when the final layouts are re-evaluated with a Gaussian wake, \NBestGauss{} ranks ahead of PSO-VNS".

4. [major] 08_beyond / 10_limits_concl — "velocity and difference moves between feasible layouts with arbitrary turbine labels make turbines overlap and are rejected by the feasibility-first comparison" — this repeats the explanation in §VIII-A almost word for word. — In §XI, write: "With feasibility-preserving initialization, \NFbFeasBestRank{} ranks first and PSO and DE stall (Section~\ref{sec:feasinit}); our statements on PSO-VNS therefore refer to random starts."

5. [major] 02_intro / 10_limits_concl — "the VNS phase helps every swarm, but the swarm phase helps only if the swarm is effective. The SSA and LX-SSA phases are statistically tied with random sampling" — Contribution 3 and the second finding of §XII are near-verbatim copies. — Rephrase §XII: "Second, VNS improves every swarm, yet only a convergent PSO phase gives VNS a better start than random sampling; the SSA and LX-SSA phases and the Laplace step add nothing."

6. [minor] 04_methods / 05_setup — "every optimizer draws its initial population first, so runs with the same seed start from the \emph{same} initial population for all methods and are paired" (05) — this duplicates the last sentence of §IV-B: "Phase~1 draws the initial population first, so runs with the same seed start from the same population as all other methods and are paired." — Delete the sentence in §IV-B, since §V-B is the right place for it.

7. [minor] 04_methods / 07_ablation / 10_limits_concl — "A swarm phase is useful only if it gives VNS a better start than random sampling with the same number of calls" (04); "The first phase of a hybrid should therefore be judged against a random-sampling control of equal cost" (07); "a swarm phase in a hybrid earns its place only if it beats a random-sampling control" (10) — the same principle is stated three times, plus Recommendation 6. — Keep it in §IV-A and as Recommendation 6. Delete the last sentence of §VII "Choice of Phase 1", and in §X shorten it to "…and a swarm phase must beat a random-sampling control (Section~\ref{sec:ablation})".

8. [minor] 04_methods / 07_ablation — "PSO-VNS and PSO differ only in how they spend the rest of the budget, and comparing them isolates the use of its second part" (04) versus "PSO-VNS and PSO share their first 3,030 calls, so their contrast measures only the use of the second half" (07) — the same explanation appears twice, and "its" is ambiguous. — In 07, write "Because PSO-VNS and PSO share their first 3,030 calls (Section~\ref{sec:hybrid}), switching to VNS…". In 04, write "…isolates the use of the second half of the budget."

9. [minor] 04_methods / 05_setup / 07_ablation — "which evaluates $\operatorname{round}(\omega B)=3{,}015$ uniform random layouts (including the initial population) and refines the best one by the same VNS" — RS-VNS is fully defined in §IV-A and then defined again in §V-A and §VII. — In §V-A write "and the random-sampling control RS-VNS (Section~\ref{sec:template})". In §VII keep only "RS-VNS spends Phase~1 on random layouts".

10. [minor] 06_results / 08_beyond — "with up to 20 times the budget, \NBudgetPhrase{}" (06) and "Only the first place of PSO-VNS is stable across budgets" (08, after "\NBudgetPhrase{}" in the same paragraph) — the same claim is made three times. The rendered text of 06 also doubles the idea ("with up to 20 times the budget, PSO-VNS keeps the best average rank at all three budgets"). — In 06, write "it also keeps the best average rank at 5 and 20 times the budget (Section~\ref{sec:budget})". In 08, delete the final sentence "Only the first place…".

11. [minor] 09_robust / 10_limits_concl — "our energy values are benchmark values, not AEP predictions" (09); "Since the energy values are benchmark values" (10 §X); "its energies are benchmark values" (10 §XI) — the same point is made three times. — Keep §IX. In §X, write "Because the benchmark models are simplified, a layout should be re-evaluated…". In §XI, delete "its energies are benchmark values,".

## B. Terminology consistency

12. [major] 05_setup / 07 / 08 / 09 / tables — "6,030 calls" (§V, §VII, Table IV) vs "6,030 evaluations" (§VIII, §IX, Table I caption) vs "objective calls" vs "Objective-function calls" (Fig. 1 axis) — there are four names for the budget unit. — Use "evaluations" throughout: "objective evaluations" at first use in §IV-A, "evaluations" after that. This also matches $B$ in the Nomenclature ("number of evaluations"). Relabel the Fig. 1 axis "Objective evaluations".

13. [minor] 01_front — "salp swarm hybrids appear stronger … outperforms all salp swarm methods … salp swarm phases" — as a compound adjective this is hyphenated everywhere else ("salp-swarm hybrids"). The body also uses "salp phases" (§X, §XII). — Abstract: "salp-swarm hybrids", "salp-swarm methods", "salp-swarm phases". §X and §XII: "salp-swarm phases".

14. [minor] 02_intro / 06 / 10 / 11 — "the setting $w=0.7$, $c_1=c_2=2$" is restated in full after §III-A defines the *old setting* (06 ¶1, 10 §X, §XI, §XII, 05 §V-A "the PSO runs with $w=0.7$, $c_1=c_2=2$"); 11 also says "the earlier PSO setting". — After the definition in §III-A, always write "the old setting". Repeat the values only in §XII if you want them there.

15. [minor] 05_setup / 06_results / 10 / Table I — "the earlier version of this study" / "the earlier study" / "the previous study" / "PREVIOUS STUDY'S METHODS" — one referent has three names. — Use "the earlier version of this study" once, then "the earlier method pool" (Table I caption: "…OF THE EARLIER METHOD POOL…").

16. [minor] 02_intro / 06 fig caption / 07 — "an ablation with a random-sampling control" (02); "eight ablation variants" (Fig. 1); TABLE IV "ABLATION" — these alternate with "component analysis", which is the term the section heading and §V-A define. — Use "component analysis" in the text and "component-analysis variants" in the Fig. 1 caption, or define it once: "component analysis (ablation)".

17. [minor] 03_model / 06 / Table II — "feasibility-aware ranking" (final ranking with the ≥15-feasible-runs rule) vs "feasibility-first rule/comparison" (the $F_p$ order). §II-C says "The comparison of final runs follows the same order (feasibility-aware ranking…)", but §V-B adds the 15/30 threshold, so it is not "the same order". — "The ranking of final runs (feasibility-aware ranking, Section~\ref{sec:setup}) extends this order with a minimum number of feasible runs."

18. [minor] 02_intro / 08 / Fig. 2 — "IEA Wind Task~37 Case Study~1 (IEA37)" then "IEA37 Case Study~1" (08, 10) and "IEA37 CS1" (Fig. 2 panels) — if IEA37 stands for the case study, then "IEA37 Case Study 1" is redundant. — Define it as "IEA Wind Task~37 (IEA37) Case Study~1" and write "IEA37 Case Study~1" everywhere, including the Fig. 2 panel titles.

19. [minor] 03_model / 05_setup — "under Wind Data Set~I" (03) vs "Data Set~I" everywhere else — this is inconsistent. — "under Data Set~I".

## C. Grammar, clarity, sentence length

20. [major] 10_limits_concl [PDF] — "with the constriction coefficients, PSO ranks second in the main comparison, PSO ranks first in that pool" — the macro \NBaseNewBest expands to "PSO", which creates a stuttering repetition. — "with the constriction coefficients, PSO ranks second in the main comparison and first in that pool, ahead of every salp-swarm method in both (Table~\ref{tab:baseline})."

21. [major] 08_beyond [PDF] — "VNS alone is then best (VNS has the lowest mean wake loss, VNS the best average rank)" — the macros expand to a triple repetition. "PSO-VNS ranks 2.83 on average" is also ungrammatical. — "VNS alone is then best, with the lowest mean wake loss and the best average rank; PSO-VNS has an average rank of \NFbFeasRankPSOVNS{} (…)".

22. [major] 08_beyond — "This mean is below the installed layout, which only the best run (\NHRBestPSOVNS~GWh/yr) and \NHRAbovePSOVNS{} of the 30 runs exceed" — the best run is one of the 6, so the sentence reads as if 7 runs exceed. The mean is also compared with a layout, not with its AEP. — "This mean is below the AEP of the installed layout, which \NHRAbovePSOVNS{} of the 30 runs exceed (best run \NHRBestPSOVNS~GWh/yr); \NHRAboveOthersPhrase."

23. [major] 09_robust — "attains higher mean benchmark objectives than the original SSA in all six representative cases … We therefore apply the rule of Section~\ref{sec:constraints} to all methods." — non sequitur: the projection rule performs better, but the paper keeps clipping "therefore". — "Because the boundary rule alone changes the results, we apply the same clipping rule (Section~\ref{sec:constraints}) to all methods."

24. [minor] 02_intro — "PSO ranks \NPSOPos{} and ahead of every salp-swarm method, including an SSA-VNS hybrid that we had developed first" — "ranks second and ahead" is ungrammatical, and "developed first" is unclear. — "PSO ranks \NPSOPos{}, ahead of every salp-swarm method, including the SSA-VNS hybrid that we originally proposed."

25. [minor] 02_intro — "Our thesis is that the conclusions of such comparisons depend strongly on how the baselines are configured" — one sentence of about 75 words carries three claims. — Split it after "controlled": "Our thesis is that such conclusions depend strongly on how the baselines are configured and how the comparison is controlled. Under a controlled protocol, a simple two-phase design … is the most consistent performer, whereas …".

26. [minor] 06_results — "it would place PSO \NOldPSOPos{} (average rank \NOldPSORank) instead of \NPSOPos, where PSO ranks ahead of SSA, LX-SSA and SSA-VNS" — "where" is ambiguous. — "…instead of \NPSOPos{}; with the constriction setting, PSO ranks ahead of SSA, LX-SSA and SSA-VNS."

27. [minor] 06_results / 07_ablation — "PSO-VNS is significantly better in \NWtlPSOW{} cases, worse in \NWtlPSOL{} and tied in \NWtlPSOT" (14/50/4, Table III) vs "\NAblPSOVNSvsPSO{} under the Holm adjustment over the contrasts" (8/58/2, Table IV) — (a) the order breaks the W/T/L convention; (b) readers see two different W/T/L results for the same pair with no warning. — (a) "better in 14 cases, tied in 50 and worse in 4". (b) In §VII add: "(8/58/2 here versus 14/50/4 in Table~\ref{tab:wtl}, because the Holm family differs: eight contrasts instead of seven comparisons)".

28. [minor] 05_setup — "and it varies $\omega$ on \NSplitCases{} cases" — this is tacked onto a sentence of about 55 words, and "it" is far from its antecedent. — End the sentence after "by the same VNS." Then add: "A split study varies $\omega$ on \NSplitCases{} cases (Section~\ref{sec:split})."

29. [minor] 05_setup — "and with the official IEA37 calculator within a relative $10^{-11}$" — a noun is missing, and "(720 archived final objectives within $8.7\times10^{-11}$)" contradicts "$10^{-6}$". — "…within a relative difference of $10^{-6}$ (maximum $8.7\times10^{-11}$ over 720 archived final objectives), …, and with the official IEA37 calculator within a relative difference of $10^{-11}$."

30. [minor] 05_setup — "only if at least 15 of its 30 runs (5 of 10) are feasible" — "(5 of 10)" is unexplained at this point. — "(5 of 10 in the 10-seed Horns Rev~1 runs)".

31. [minor] 06_results — "and makes \NPhTwoMadeFeas{} runs that were infeasible at the switch feasible" — the object and its complement are too far apart. — "and turns \NPhTwoMadeFeas{} runs that were infeasible at the switch into feasible ones".

32. [minor] 08_beyond — "the best layout would rank \NIEARankSixteenThirtyK{} of \NIEAPubFeasSixteen{}" → [PDF] "third of 9 and seventh of 8" — ordinals in words are mixed with counts in numerals. Also, "at 120,030 evaluations the loss is …, but VNS, RS-VNS and DE are lower" compares methods with a loss. — Write "third of nine and seventh of eight" (use a word-form macro), and "…but those of VNS, RS-VNS and DE are lower".

33. [minor] 10_limits_concl — "The three findings share one implication: in continuous WFLOP" — this transition is abrupt: no "three findings" are listed in §IX/§X, and the Introduction lists four contributions. The article is also missing ("the continuous WFLOP", as in 02). — "The three main findings (Sections~\ref{sec:results}--\ref{sec:robustness}) share one implication: in the continuous WFLOP, …". Also write "the continuous WFLOP" in §XII.

34. [minor] 10_limits_concl — "\item rank feasibility-aware: report the share of feasible runs" — an adjective is used as an adverb. — "\item use feasibility-aware ranking: report the share of feasible runs…". In §XI, change "still differ by \NIEAGap…\% and …\%" (which renders negative values) to "still fall short of … by 1.33\% and 4.20\%".

## D. Structure, headings, flow

35. [minor] 07_ablation / 09_robust — Run-in heads use "\emph{Laplace step.}" (period) in §V/§VII/§VIII but "\emph{Power curve and wake model:}" (colon) in §IX. §VII has a lone subsection "A. Budget Split" (no B), and §IX starts with no lead-in sentence, unlike §VI and §VIII. — Use the period form throughout. Make "Budget split" a run-in paragraph, or promote the §VII run-ins to subsections. Open §IX with: "This section checks whether the conclusions of Section~\ref{sec:results} depend on the benchmark models."

36. [minor] 08_beyond — "\subsection{Published Kusiak--Song Results and Computational Cost}" — this two-sentence subsection mixes two unrelated topics and ends the section abruptly. The runtime sentence belongs with the protocol. — Move the runtime sentence to §V-C "Verification and code". Keep the literature pointer as one sentence at the end of §VIII-C, or put it in §VI.

## E. IEEE conventions, math, numbers, references

37. [minor] 03_model / 09 table — "$E_{\rm ideal}=14045.74$" — the paper uses thousands separators elsewhere (6,030; 366,941.6), so this is inconsistent. Table IX [PDF] prints "-17.7 -36.9 … -21.5" with hyphens instead of minus signs. The units "GWh/yr" and "MWh per year" are also mixed. — "$E_{\rm ideal}=14{,}045.74$". In the table generator, use "$-17.7$". Use "GWh/yr" (or "MWh/yr") consistently.

38. [minor] 01_front (Nomenclature) — Symbols used but not listed: $n$ (=2N, undefined in §III), $\mathbf H$, $\xi_1,\xi_2$, $T$, $lb_m,ub_m$, $\gamma,\phi,\chi$, $\mathcal N_k,\delta_k$, $h$, $\mathbf e_m$ (undefined in text), $u_{ij},v_{ij}$, $f(s)$, $s_{\text{cut-in}},s_{\text{rated}},P_{\text{rated}},\lambda_0,\lambda_1$, $p_\theta,p_s$, $\theta_{\rm met}$, $\boldsymbol\rho_{1,2}$, $K_{\rm G}$, $r_{\rm rb}$, $\bar\tau$, $\Delta L$. Also, $\alpha$ means both the half-cone angle and the significance level ("two-sided with $\alpha=0.05$"). — Add these entries, define "$n=2N$" and "$\mathbf e_m$, the $m$th unit vector" in §III. Replace "with $\alpha=0.05$" with "at the 0.05 significance level" (and likewise in the table captions). Reword "$\ell_{ij}$ Distance of $i$ and $j$" as "Distance between turbines $i$ and $j$", and "$d_{ij}$" as "Downstream distance of turbine $i$ from turbine $j$".

39. [minor] 01_front / 05_setup — Acronyms are not expanded at first use: MS-SLSQP/SLSQP, AEP (first in §II-B), SD (§VIII-B), IEA. The Index Terms are not in alphabetical order (IEEE). The Abstract also overstates Horns Rev: "exceeds the energy yield of the installed Horns Rev~1 block" is true only of the best run at 6,030 evaluations. — "MS-SLSQP, a multistart sequential least-squares programming method (SciPy's SLSQP~\cite{Kraft1988}) …"; "annual energy production (AEP)"; "standard deviation (SD)". Index Terms: "Benchmarking, metaheuristics, particle swarm optimization, variable neighborhood search, wake effect, wind farm layout optimization". Abstract: "…and its best layout exceeds the energy yield of the installed Horns Rev~1 block".

40. [minor] 11_back / 05 / 08 — References and supplement pointers. (a) [30] Solanki2023 has no volume or pages: add them, or add "early access". (b) [40] PyWake and [41] IEA37repo have no year and no "Accessed: <date>" (IEEE online format). (c) [22] should follow IEEE arXiv form: "…open issues," 2020, *arXiv:2007.03488*. [22] is acceptable, but consider "arXiv:2007.03488, 2020." (d) Supplement pointers are inconsistent: "of the supplement" (§VIII-A, §VIII-D) vs "of the supplementary material", and some have neither ("Table S14", "Table S12", "Table S18"). Main-text first mentions are out of order: Figs. S12–S14 come before Fig. S3, and Table S14 before S12. Supplement sections render as "Section S-VI.A" while main sections render as "VI-A". — Use "of the supplementary material" at the first pointer of each section. Renumber the supplement so that first mentions are in order, or cite in order. Set `\thesubsection` in the supplement to `\thesection-\Alph{subsection}` so it renders as "Section S-VI-A" or "S-VI-A".

---

### Items checked and found OK
- Reference list numbering follows order of first citation ([1]–[43]); en-dashes in page ranges; "et al." italicized.
- Equation/table/figure numbering: Tables I–IX and Figs. 1–2 are first cited in numerical order. (Table I is first cited in the Introduction; Fig. 1 in §VI-C, with its caption pointing forward to Table IV, which is acceptable.)
- en-dashes in numeric ranges ("1–30", "N = 11–15", "0.24–1.82 ms") and nonbreaking spaces before units are correct.
- Unresolved placeholders remain (author-side): \TBD in the thanks footnotes, Acknowledgment, Data Availability DOI, §VIII-D, the Table VIII footnote, §XI, and the biographies.

# W5 — Literature and positioning (Phase 6)

Owner: W5. Files: `optA/drafts/positioning.md` (this file) and `optA/newrefs/phase6.tex` (12 verified `\bibitem`s).
All references were verified by web search on 2026-09-28. Publisher sites, Crossref and dblp could not be fetched
from the container, so each entry was checked against at least one indexed record (publisher record, JMLR/arXiv,
RePEc, OSTI/ADS, institutional repository). The source is in the `% VERIFY-REF` line of each entry. Months and
issue numbers that could not be confirmed are left out.

Note on comment placement: each verification comment is the first line *after* its `\bibitem{Key}` line, not
above it. `optA/merge_newrefs.py` attaches any text before a `\bibitem` to the previous entry, so a comment placed
above would travel with the wrong reference after reordering. This follows the `% VERIFY-REF` convention of
`11_back.tex`.

## 0. Verified references in `newrefs/phase6.tex` (12; none duplicates an existing key)

| Group | Key | Short reference | Where to cite |
|---|---|---|---|
| (a) | `Benavoli2017` | Benavoli, Corani, Demšar, Zaffalon, JMLR 18(77):1–36, 2017 | 05_setup statistics (Bayesian signed-rank + ROPE), tab:equivalence caption |
| (a) | `LaTorre2021` | LaTorre et al., Swarm Evol. Comput. 67:100973, 2021 | intro (benchmarking lit.), contribution 1, Sec. X recommendations, 05_setup ("non-PSO settings are untuned defaults") |
| (a) | `Derrac2011` | Derrac, García, Molina, Herrera, Swarm Evol. Comput. 1(1):3–18, 2011 | intro next to Demsar2006; 05_setup (Friedman/Wilcoxon/Holm) |
| (a) | `Hooker1995` | Hooker, J. Heuristics 1(1):33–42, 1995 | intro: controlled vs. competitive testing; component analysis ("why, not only which") |
| (a) | `Campelo2019` | Campelo & Takahashi, J. Heuristics 25(2):305–338, 2019 | 05_setup (30 runs and 68 instances; power); Limitations (power of run-level tests) |
| (a) | `Lakens2017` | Lakens, Soc. Psychol. Personal. Sci. 8(4):355–362, 2017 | 05_setup / tab:equivalence (TOST, smallest effect size of interest = margin) |
| (b) | `CamachoVillalon2023` | Camacho-Villalón, Dorigo, Stützle, ITOR 30(6):2945–2971, 2023 | intro next to Sorensen2015/Castelli2022 |
| (c) | `Bonyadi2017` | Bonyadi & Michalewicz, Evol. Comput. 25(1):1–54, 2017 | 04_methods PSO paragraph (stability analyses, inertia/constriction forms) next to Poli2009/Cleghorn2018 |
| (d) | `Thomas2023` | Thomas et al., Wind Energy Sci. 8:865–891, 2023 | intro WFLOP paragraph next to Baker2019; Limitations (no CMA-ES etc.) |
| (d) | `Wilson2018` | Wilson et al., Renew. Energy 126:681–691, 2018 (GECCO 2015 WFLO competition) | intro WFLOP comparison work |
| (d) | `Azlan2021` | Azlan, Kurnia, Tan, Ismadi, Renew. Sustain. Energy Rev. 135:110047, 2021 | intro (review of methods on the classical benchmark wind scenarios) |
| (d) | `Feng2015` | Feng & Shen, Renew. Energy 78:182–192, 2015 | 04_methods / 07_ablation: RS-VNS control, relation to random search in the WFLOP |

What each reference says (only what the verification sources confirm):
- **Benavoli2017**: argues for Bayesian analysis over null-hypothesis tests when comparing algorithms on
  several data sets, including the Bayesian signed-rank test with a region of practical equivalence (ROPE). The
  W1 pipeline already cites it.
- **LaTorre2021**: gives four guidelines for proposing and comparing bio-inspired algorithms: (1) choice of
  benchmarks and reference algorithms, (2) statistical validation with the correct tests, (3) component analysis
  and parameter tuning of the proposal, and (4) why the algorithm is useful, meaning a real advance over the state
  of the art and not only novelty.
- **Derrac2011**: a tutorial on nonparametric pairwise and multiple-comparison procedures for evolutionary and
  swarm algorithms over multiple problems.
- **Hooker1995**: competitive testing tells which algorithm is faster, not why. It argues for controlled
  experiments.
- **Campelo2019**: a method for choosing the number of problem instances and runs per instance for a desired power
  and accuracy.
- **Lakens2017**: equivalence tests (TOST) with bounds set by the smallest effect size of interest. A
  non-significant difference is not evidence of absence.
- **CamachoVillalon2023**: shows that six metaphor-based algorithms are misleading, in that they reuse known search
  ideas under new metaphors.
- **Bonyadi2017**: a review of PSO for continuous single-objective problems, covering its analyses (stability and
  convergence) and its modifications.
- **Thomas2023**: eight methods (gradient-based, gradient-free, hybrid), each run or directed by researchers
  experienced with it, on one common objective (81 turbines, four discrete regions, concave boundaries). The
  layouts were all different but had similar AEP.
- **Wilson2018**: results of the second GECCO WFLO competition (2015). A common evaluation API was used, and the top
  four of eight teams are compared.
- **Azlan2021**: reviews layout optimization methods under the three classical wind scenarios, which are used as a
  performance benchmark.
- **Feng2015**: a random-search algorithm that iteratively improves a feasible layout inside the feasible space. It
  shows that random search is a recognised WFLOP method in its own right, so a random-sampling control is a natural
  baseline.

### Optional extras (verified the same way, NOT in phase6.tex because of the ≤ 12 limit; add if space allows)

```latex
\bibitem{Eberhart2000}
% VERIFY-REF (2026-09-28): Semantic Scholar record, Google Scholar lookup by DOI -- authors, title, Proc. CEC 2000 vol. 1 pp. 84-88, DOI confirmed
R. C. Eberhart and Y. Shi, ``Comparing inertia weights and constriction factors in particle swarm optimization,'' in \textit{Proc. 2000 Congr. Evol. Comput. (CEC)}, 2000, vol. 1, pp. 84--88, doi: 10.1109/CEC.2000.870279.

\bibitem{ThomasNing2018}
% VERIFY-REF (2026-09-28): IOPscience record 10.1088/1742-6596/1037/4/042012, BYU ScholarsArchive -- authors, title, J. Phys.: Conf. Ser. vol. 1037, art. 042012, 2018, DOI confirmed (wake expansion continuation)
J. J. Thomas and A. Ning, ``A method for reducing multi-modality in the wind farm layout optimization problem,'' \textit{J. Phys.: Conf. Ser.}, vol. 1037, no. 4, Art. no. 042012, 2018, doi: 10.1088/1742-6596/1037/4/042012.

\bibitem{Aranha2022}
% VERIFY-REF (2026-09-28): Springer record 10.1007/s11721-021-00202-9 (PDF header "Swarm Intelligence (2022) 16:1-6"), Aston research portal -- 8 authors, title, vol. 16, pp. 1-6, 2022, DOI confirmed
C. Aranha \textit{et al.}, ``Metaphor-based metaheuristics, a call for action: The elephant in the room,'' \textit{Swarm Intell.}, vol. 16, no. 1, pp. 1--6, 2022, doi: 10.1007/s11721-021-00202-9.

\bibitem{Shakoor2016}
% VERIFY-REF (2026-09-28): IDEAS/RePEc v58y2016icp1048-1059, NASA ADS 2016RSERv..58.1048S -- authors, title, vol. 58, pp. 1048-1059, 2016, DOI confirmed
R. Shakoor, M. Y. Hassan, A. Raheem, and Y.-K. Wu, ``Wake effect modeling: A review of wind farm layout optimization using Jensen's model,'' \textit{Renew. Sustain. Energy Rev.}, vol. 58, pp. 1048--1059, 2016, doi: 10.1016/j.rser.2015.12.229.

\bibitem{Wei2023}
% VERIFY-REF (2026-09-28): IEEE Xplore document 10113802 (title), reference list of J. Electr. Eng. Technol. 10.1007/s42835-026-02821-w (vol. 11, no. 4, pp. 1784-1794, DOI) -- MPCE anchor for R1 issue 24
S. Wei, H. Wang, Y. Fu, F. Li, and L. Huang, ``Electrical system planning of large-scale offshore wind farm based on $N+$ design considering optimization of upper power limits of wind turbines,'' \textit{J. Mod. Power Syst. Clean Energy}, vol. 11, no. 4, pp. 1784--1794, 2023, doi: 10.35833/MPCE.2022.000656.
```

Priorities if the pool can grow:
- `Eberhart2000`: R2 issue 3 asks for it for the equivalence of the inertia and constriction forms, which
  04_methods states.
- `ThomasNing2018`: R2 issue 18. It gives the reason MS-SLSQP is handicapped and names the IEA37 winner's
  continuation method.
- `Wei2023`: only for the MPCE version (R1 issue 24: a recent MPCE wind-farm planning paper).
- `Aranha2022` and `Shakoor2016` are optional.

Not verified, so do not cite (R1 issue 3): I could not find a verifiable WFLOP paper that states the setting
w = 0.7, c1 = c2 = 2. Checking this needs the full texts, and the publisher sites are blocked here. The authors
must check the full texts of candidate papers themselves, e.g., their own earlier work or the PSO-based WFLOP papers
already cited. Otherwise use R1's softened wording: "a setting that is easy to adopt from the original PSO
(c1 = c2 = 2~\cite{Kennedy1995}) and lies outside the order-2 stability region". The current 02_intro already
takes this route ("keeps the acceleration coefficients of the original PSO"), and the abstract no longer says
"common".

## (i) Related-work paragraph for the introduction (drop-in LaTeX; replaces/extends paragraph 3 of 02_intro)

Six sentences. Numbers are macros from W1 (`\NEqMargin`, `\NEqPSOVNSvsPSOCI`, `\NEqPSOVNSvsPSOHolds` = yes at the
time of writing). The last sentence may only stay if `\NEqPSOVNSvsPSOHolds` is still "yes" in the final build
(needs a CHECK-FINAL id from W1).

```latex
Guidelines for comparing metaheuristics are well established: controlled rather than competitive experiments~\cite{Hooker1995}, equal budgets and reported parameter settings, component analysis of the proposed method, and paired nonparametric tests with multiplicity correction~\cite{Demsar2006,Derrac2011,BartzBeielstein2020,LaTorre2021}, with sample sizes chosen for a stated power~\cite{Campelo2019} and equivalence or Bayesian analyses when the question is whether two methods differ by a practically relevant amount~\cite{Benavoli2017,Lakens2017}. Critical studies have shown that several metaphor-based algorithms, SSA among them, reuse known search operators~\cite{Sorensen2015,CamachoVillalon2023,Castelli2022}, and that the behavior of PSO depends on whether its coefficients lie in the stability region~\cite{Poli2009,Bonyadi2017,Cleghorn2018}. In the WFLOP, controlled comparisons with a common objective exist: the GECCO layout competition~\cite{Wilson2018} and the IEA Wind Task~37 studies~\cite{Baker2019,Thomas2023}, in which optimizers applied by experienced users reached different layouts of similar energy yield, while reviews of the classical benchmark wind scenarios list many metaheuristics compared under differing budgets and settings~\cite{Azlan2021}, and random search is itself a WFLOP method~\cite{Feng2015}. These studies compare methods as submitted or as published; they do not examine how the configuration of a baseline or the choice of control changes which method appears best. This paper applies the established practice to the continuous Kusiak--Song WFLOP and adds three elements that are, to our knowledge, not reported for this problem: evidence that one baseline setting, a PSO outside the order-2 stability region, reverses the order of PSO and the salp-swarm hybrids at 6,030 evaluations; a component analysis of two-phase swarm--VNS hybrids against a random-sampling control of equal cost; and an equivalence analysis with a stated margin. With a margin of $\pm\NEqMargin$ percentage points of wake loss, a well-configured PSO and PSO-VNS are practically equivalent over the 68 cases (90\% interval of the mean difference \NEqPSOVNSvsPSOCI), so PSO-VNS is proposed as a consistent reference method, not as a superior one.
% CHECK-FINAL [W1 id for tab:equivalence]: \NEqPSOVNSvsPSOHolds == yes (TOST, margin \NEqMargin, case means)
```

Integration notes for the lead:
- The paragraph takes over the role of the current paragraph 3 opener ("Many of these studies introduce a new or
  modified metaheuristic …"). Keep the four WFLOP-specific issues (baseline, budgets, feasibility, statistics)
  that follow it. If the intro exceeds 1.0 page, drop the `Campelo2019` clause and the `Feng2015` clause first. Both
  can move to 05_setup and 04_methods.
- Contribution 1 should then say "the application of established practice~\cite{Derrac2011,LaTorre2021}…", as R1
  issue 5 and R2 issue 20 request.
- Check the "90\% interval" wording against W1. TOST at α = 0.05 corresponds to a 90% CI. If W1 reports a 95%
  interval, change the label.
- The margin is fixed after the primary analysis (PHASE6 decision). State this and report the minimal equivalence
  margin (`\NEqPSOVNSvsPSOMin`), as Lakens2017 recommends for a transparent smallest effect size of interest.
- Do not add "first" or "novel". "To our knowledge, not reported for this problem" is the strongest claim the
  search supports. I found no WFLOP study with a random-sampling ablation control or an equivalence test, but this
  was a web search, not a systematic review.

## (ii) Literature-related objections per review and how the new references answer them

**R1 (editor, MPCE)**
- Issue 5: the protocol is established practice, not a contribution. Cite `Derrac2011`, `LaTorre2021`,
  `BartzBeielstein2020` and `Demsar2006` in contribution 1. Frame the contribution as application to the WFLOP plus
  two WFLOP-specific elements (feasibility-aware ranking, random-sampling control) plus the equivalence analysis.
  The related-work paragraph above does this.
- Issue 3: "common in WFLOP comparisons" is uncited. No verifiable source was found. Keep the softened wording that
  ties the setting to `Kennedy1995` (already done in the current 02_intro).
- Issue 13: the title. See (iii); no "Revisiting" and no "Design".
- Issue 15: the tone of the Sörensen citation. Pairing `Sorensen2015` with `CamachoVillalon2023` and `LaTorre2021`
  makes it a methodological point, not a verdict on the WFLOP literature. The paragraph says the WFLOP studies
  "compare methods as submitted or as published". That describes their design and does not criticise them.
- Issue 24: MPCE anchoring. The optional `Wei2023` (verified MPCE, 2023, offshore wind-farm electrical planning) can
  join `Srikakulapu2018`/`Hou2019` in paragraph 1 of the MPCE version. `Solanki2023` volume/pages are still open
  (owner: authors).
- Issue 2 (scope): see (iv).

**R2 (metaheuristics)**
- Issue 1: order-1/order-2 stability. `Poli2009` and `Cleghorn2018` exist; `Bonyadi2017` adds the review-level
  reference for PSO stability analyses.
- Issue 3: inertia vs. constriction equivalence and "standard PSO". `Bratton2007` exists; `Eberhart2000` is
  verified but optional (see above). If it is not added, cite `Clerc2002` and `Bonyadi2017` for the equivalence.
- Issue 6: tuning asymmetry. `LaTorre2021` (guideline 3: parameter tuning) supports the sentence "all non-PSO
  settings are untuned defaults and the rankings are conditional on them". `Zaharie2002` exists for the DE check.
- Issue 7: non-significance is not equivalence. `Lakens2017` (TOST) and `Benavoli2017` (Bayesian signed-rank with
  ROPE) justify tab:equivalence. The W1 results support "LX-SSA-VNS equivalent to RS-VNS" and "SSA-VNS a small but
  non-equivalent gain", which is the precise wording R2 and R3 ask for.
- Issue 8: RS-VNS is a weak control. `Feng2015` places the control in the WFLOP random-search literature. It also
  shows that a stronger random search, starting feasible and staying in the feasible space, exists. Cite it in the
  Limitations sentence on stronger controls.
- Issue 14: the SSA critique. `Castelli2022` exists; `CamachoVillalon2023` generalises the point to other
  metaphor-based methods.
- Issue 17: pool without CMA-ES/L-SHADE, and Thomas et al. `Thomas2023` is now verified. Cite it in the intro and in
  Limitations: a comparison of eight methods, each applied by experienced users on a common objective, gave
  different layouts of similar AEP. Do not attribute to it the claim that "initialization and multistarts matter as
  much as the algorithm" (R2's paraphrase); I could not verify that from the abstract.
- Issue 18: MS-SLSQP handicap. The optional `ThomasNing2018` (wake expansion continuation) is verified.
- Issue 19: memetic/PSO+VNS literature (Petalas 2007, Tasgetiren 2007, Carrizosa 2012). Not verified here because
  of the ≤ 12 limit and the priority list given to W5. `Mladenovic2008` (continuous VNS) exists and covers part of
  it. If the lead wants a memetic citation, verify Petalas et al., Ann. Oper. Res. 156:99–127, 2007 separately.
- Issue 20: credit the guidelines. Use `LaTorre2021` and `Benavoli2017` in the intro and in Sec. X ("we follow and
  quantify for the WFLOP the guidelines of~\cite{LaTorre2021,BartzBeielstein2020}").

**R3 (statistics)**
- Issues 1–3: non-significance read as a tie or as non-inferiority. Answered by the TOST (`Lakens2017`) and the
  Bayesian ROPE (`Benavoli2017`) in tab:equivalence.
- Issue 4: 68 non-independent cases. `Campelo2019` separates instance-level from run-level sample size. Use it to
  say that the conclusions are conditional on the 68 instances (six clusters) and that the 30 runs control accuracy
  within an instance. The Limitations sentence already says the p-values describe this benchmark.
- Issue 5: multiplicity across families. `Derrac2011` covers Holm and the family-wise procedures. No new
  literature is needed beyond defining the families.
- Issue 18: reproduction of archived distributions. TOST (`Lakens2017`) with a stated margin is the recommended
  replacement for "22 of 24 Mann–Whitney tests with p ≥ 0.05".

**R4 (wind energy)**
- Issue 6: practical relevance and scale. `Thomas2023` (81 turbines, gradient-based and gradient-free methods on a
  common objective) and `Baker2019` give the context for large farms. Suggested sentence for Sec. X: "Our results
  do not show how PSO-VNS scales to farms of 50–150 turbines, where comparisons with a common objective include
  gradient-based methods~\cite{Baker2019,Thomas2023}."
- Issues 10/12/13 (K = 0.04 offshore value, IEA37 attribution and discretisation). These need wake-model or IEA37
  references that are outside W5's priority list. For the IEA37 participant, check `Baker2019` Table 4 (author
  task). If an offshore K reference is wanted, the verified review `Shakoor2016` covers Jensen-model practice in
  WFLOP. It is not a calibration source, so do not cite it for the value 0.04.
- The wake-model review `Shakoor2016` and the benchmark review `Azlan2021` show that the Jensen/Kusiak–Song-type
  benchmark is the common test bed. This supports R1 issue 10's justification: "we keep the benchmark because most
  published metaheuristic comparisons use it".

**R5 (claims audit)** — no literature objection. The audit concerns claim-to-evidence links. The new references
support the reworded claims ("equivalent within ±margin" via `Lakens2017`, and "no detectable benefit" replaced by a
TOST statement).

**R6 (copy-edit)** — reference-format items only (Solanki2023 volume/pages, PyWake/IEA37 access dates, arXiv
form). All 12 new entries follow the IEEE journal form used in 11_back.tex: abbreviated journal names, "Art. no."
for article numbers, en-dash page ranges, DOI last, and "et al." for more than six authors (Thomas2023, Wilson2018).
The index terms in (iii) are in alphabetical order, as IEEE requires.

**R7 (reproducibility)** — no literature objection. `LaTorre2021` and `BartzBeielstein2020` can be cited for making
code and data public (recommendation 7), but the Zenodo DOI remains an author task.

**R8 (supplement)** — Issue 2: the S-VII literature table of published Kusiak–Song results is a placeholder. No
new reference solves this. Either the authors verify the published values or the table and its main-text sentence
are dropped. `Azlan2021` could be cited in the one remaining sentence as a review of results on the classical
benchmark scenarios.

## (iii) Title (without "Design") and keywords

Options (≤ 18 words, no "Revisiting", en dash in PSO–VNS):
1. **Baseline Configuration and Controls in Metaheuristic Comparisons for Wind Farm Layout Optimization, with a
   PSO–VNS Reference Method** (preferred: names M1, M2 and M3 and avoids "novel" and "design")
2. Controlled Comparison of Metaheuristics for Continuous Wind Farm Layout Optimization: Baseline Configuration,
   Component Analysis and a PSO–VNS Reference Method
3. How Baseline Configuration Changes Metaheuristic Comparisons in Wind Farm Layout Optimization: A Controlled
   Study with a PSO–VNS Reference Method

Running head for option 1: "SOLANKI et al.: BASELINE CONFIGURATION IN METAHEURISTIC COMPARISONS FOR WIND FARM
LAYOUT OPTIMIZATION".

Keywords (6, IEEE alphabetical): Benchmarking, equivalence testing, particle swarm optimization, variable
neighborhood search, wake effect, wind farm layout optimization.
For MPCE, "equivalence testing" may be replaced by "metaheuristics", which is the current set.

## (iv) Cover-letter scope arguments

**Swarm and Evolutionary Computation (SEC)** — the better fit for the core of the paper:
- The protocol implements the two methodological papers that SEC itself published (`Derrac2011`, `LaTorre2021`).
  In particular it follows LaTorre's guideline on component analysis and parameter tuning, and it quantifies what
  ignoring that guideline does to a real-world comparison.
- M1 connects PSO stability theory (order-2 region; `Poli2009`, `Cleghorn2018`) to an empirical reversal of a
  method ranking. This is a theory-to-practice result on a swarm algorithm.
- M2 is a controlled component analysis with a random-sampling control, together with equivalence testing. It
  answers the metaphor critique (`Sorensen2015`, `CamachoVillalon2023`, `Castelli2022`) with evidence instead of
  argument, including a candid negative result on the authors' own LX-SSA.
- The WFLOP is a constrained real-world application with public benchmarks (Kusiak–Song, IEA37) and public code
  and data.

**MPCE (Journal of Modern Power Systems and Clean Energy)**:
- The layout fixes a farm's lifetime energy yield. Give one sentence in AEP terms, from macros (R1 issue 4).
- The paper uses applied cases that planners recognise: the Horns Rev 1 block and IEA Wind Task 37 Case Study 1.
- It continues MPCE's own line of work on wind-farm layout and collector design (`Srikakulapu2018`, `Hou2019`,
  `Asaah2021`; optionally `Wei2023`). Its message for practitioners is that published optimizer rankings used to
  choose layout tools depend on how baselines are configured.
- The recommendations section gives planners and tool developers a checklist for evaluating layout optimizers.
- Be frank in the letter: the core is methodological, on a stylized benchmark. R1 rates the scope risk as
  moderate to high, and the letter should not overstate the power-system content.

Recommendation: if the authors are free to choose, SEC matches the contribution. MPCE requires the AEP translation
and the practical framing above to be in the paper, not only in the letter.

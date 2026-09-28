# R1 (round 2): Handling-editor assessment

Manuscript: MPCE_PSO_VNS.tex + optA/01–11 (built PDF: 13 pages; the target is 12). Supplement: 38 pages.
I read the built PDF, so the numbers quoted below are the typeset macro values.
Line numbers refer to the optA/*.tex source files. Each paragraph there is a single long line.

Space accounting used below: one column holds about 57 text lines at a 12-pt pitch, a page about 114 lines, and half a page about 57 lines.
Page 13 currently holds refs [55]–[59] (about 13 lines) and four biography stubs. Real biographies (about 5–6 lines each), the funding line and the acknowledgment will add about 25–30 lines.
**To end on p. 12 with some safety, the main text must lose at least 45 lines. I recommend cutting about 75–100 lines (0.65–0.9 page).**

---

## A. Assessment in brief

- **Contribution and significance.** The paper is now honest and technically careful. Its real contribution is methodological, a benchmarking case study on continuous WFLOP. Three findings carry it:
  - M1: a non-convergent PSO baseline reverses the ranking.
  - M2: salp-swarm first phases do no better than random sampling in the farm.
  - M3: PSO-VNS ≈ a well-configured PSO.

  M1 and M2 are useful to the WFLOP community. They are not new to the EC community: PSO order-2 theory and "metaphor" critiques already exist. Significance for SWEVO rests on two things:
  - showing that the configuration error actually occurs in WFLOP practice (not yet shown; issue 5);
  - the WFLOP-specific mechanism: clipping plus a penalty with mismatched units produces feasibility collapse.
- **Framing.** M1–M3 are presented coherently and in the same order everywhere. However, the title and the "reference method" label still sell PSO-VNS, while the evidence says:
  - PSO-VNS ≈ PSO;
  - the "VNS" is a single compass-search descent;
  - its advantages are conditional (random starts, one Horns Rev block, a post hoc subgroup).

  The paper would be stronger if the title and abstract led with configuration and controls, and PSO-VNS served as the vehicle of the component analysis.
- **Overclaims that remain:**
  - "the Laplace step is harmful" (confounded with the doubled evaluation cost);
  - the abstract's robustness sentence (the cluster analysis qualifies precisely the M2 contrasts);
  - "practically equivalent on average" (driven by 48 of 68 easy cases with N < 10; the post hoc margin sits just above the observed bound);
  - "feasibility reliability" (one site, random starts only).
- **Underexplained:**
  - RSD-VNS, the control that carries M2, is not defined in the Methods or Protocol sections;
  - the text and Table IV report different ΔL values for the same contrasts;
  - the IEA37 table has 8 methods but the text says ten.
- **Length.** Content is repeated four times (abstract → contributions → results → discussion/conclusion). The contributions list alone is about 60 lines, and nine tables plus one full-width figure compete for space. The cuts in Section D remove about 1 page of candidates with no loss of evidence, because everything moved already has a home in the supplement.

---

## B. Numbered issues

### Major

**1. [major] Abstract overclaims and ambiguities** — `optA/01_front.tex:13`

Five problems:
- (a) "with the constriction setting it ranks first" silently switches to the Table I pool (no PSO-VNS). In the main comparison PSO ranks second.
- (b) "(the Laplace step is harmful)" contradicts `07_ablation.tex:24`: the contrasts "do not separate the Laplace step from this cost".
- (c) "These conclusions hold across qualification thresholds, farm clusters and one Holm family" is false for M2. At cluster level, exactly SSA-VNS vs RSD-VNS and LX-SSA-VNS vs RSD-VNS are the exceptions (`05_setup.tex:19`).
- (d) Two things are unqualified:
  - The equivalence margin was fixed after the primary analysis.
  - "gains in larger Data Set II layouts" is a post hoc subgroup of 10 cases; Data Set I shows the opposite tendency (PSO lower in 6/10).
- (e) The conditions under which PSO-VNS leads are missing: random starts and 6,030 evaluations.

**Fix.** Replace lines 13 with the following (≈225 words, macros kept):

> Wind farm layouts fix a plant's lifetime energy yield, yet comparisons of layout optimizers are hard to trust. We compare eight methods for continuous wind farm layout optimization under a controlled protocol: 68 benchmark cases, equal budgets of 6,030 evaluations from random initial layouts, seed-paired runs, feasibility-aware ranking, and Holm-corrected and equivalence tests. First, baseline configuration changes the conclusions. The particle swarm optimization (PSO) setting $w=0.7$, $c_1=c_2=2$ violates the order-2 (variance) convergence condition; its swarm does not contract, few of its evaluations are feasible, and in an eight-method pool PSO ranks behind two salp-swarm hybrids. With the standard constriction setting, the same pool ranks PSO first. Second, a component analysis with equal-cost random-sampling controls shows that, as a first phase before local search, a salp swarm (SSA) adds nothing beyond random sampling in the farm, the Laplacian SSA as published is worse, and PSO gives a clear gain; at this budget the second phase, variable neighborhood search (VNS), reduces to one compass-search descent. Third, PSO followed by VNS (PSO-VNS) has the best average rank but is practically equivalent to the well-configured PSO (post hoc margin $\pm\NEqMargin$ percentage points of wake loss); it is more reliably feasible on a Horns Rev~1 block (\NHRFeasPSOVNS{} vs \NHRFeasPSO{} runs) and, post hoc, better on the larger Data Set~II layouts. From feasible starts VNS alone ranks first, and the best published IEA Wind Task~37 layouts remain better. We derive seven recommendations for benchmarking layout optimizers.

Update the CHECK comments (C27: use "the same pool"; drop the [X11]/[X14] clause).

**2. [major] "The Laplace step is harmful": claim versus stated confound**

Locations:
- `01_front.tex:13`
- `02_intro.tex:9` ("the Laplace step … makes the hybrid worse")
- `02_intro.tex:22` ("The Laplace step is harmful ($p=\NCmSSAVNSvsLXSSAVNSP$)")
- `07_ablation.tex:24` ("The Laplace step makes the salp swarm worse.")
- `10_limits_concl.tex` (the Conclusion does not repeat it, which is good)

Problem. LX-SSA as published spends 2N_p evaluations per iteration and runs half as many iterations (`04_methods.tex:48`). The contrasts therefore test "Laplace step + half the iterations", which the paper itself states. A SWEVO reviewer will call this an overclaim. It also reads as the authors disowning their own earlier method on confounded evidence.

Fix:
- Intro:9: "…while salp-swarm first phases add nothing beyond random sampling and the Laplacian SSA (LX-SSA)~\cite{Solanki2023}, proposed earlier by two of the authors, is worse as published."
- Intro:22: "LX-SSA as published (Laplace step and twice the evaluations per iteration) is worse than SSA ($p=\NCmSSAVNSvsLXSSAVNSP$)."
- 07:24: replace the first sentence with "LX-SSA as published makes the salp swarm worse." Then move the last sentence ("As LX-SSA costs …") to directly after it, so that the caveat precedes the numbers.
- Rename the Table IV contrast labels "Laplace step in the hybrid/alone" to "LX-SSA as published (hybrid/alone)". This changes the pipeline in `analysis/mpce_results.py`.

**3. [major] "VNS" in name only; the method name, keywords and recommendations still sell VNS**

Locations: `04_methods.tex:89`, `07_ablation.tex:8`, `01_front.tex:1,31`, `10_limits_concl.tex:16`.

Problem. At 6,030 evaluations, shaking uses 0.02% of the Phase-2 evaluations and the first descent yields 97.5% of the Phase-2 gain. The descent is unfinished in 39 of 40 runs with N ≥ 12. The tested design is therefore "PSO + one unfinished compass-search descent", and the VNS neighborhoods are never exercised where it matters.

The text now admits this (good). Yet the title, keywords ("variable neighborhood search") and the Discussion ("PSO-VNS a more reliable one") present VNS as the active ingredient. The same caveat applies to "VNS alone ranks first from feasible starts": that is a compass search from the best feasible initial layout. This is arguably the most practically relevant finding of the paper, and it is buried in Section VIII-A.

Fix:
- (i) Keywords: replace "variable neighborhood search" with "local search" or "hybrid metaheuristics".
- (ii) Add to 04_methods:89 after "…matter only for small $N$": "PSO-VNS at this budget is therefore best read as PSO followed by a compass-search descent; we keep the VNS name because the implemented algorithm is the basic VNS."
- (iii) In the Discussion (10:16) write: "…from feasible starts, a single local search (VNS alone) is competitive."

**4. [major] "Practically equivalent on average" depends on benchmark composition, and the post hoc margin sits just above the observed bound**

Locations: `05_setup.tex:15`, `06_results.tex:30`, `02_intro.tex:32`, `10_limits_concl.tex:12,30`.

Problem. 48 of the 68 cases have N < 10. There, PSO and PSO-VNS coincide (mean |ΔL| = 0.03 pp). The 68-case equivalence is therefore largely a statement about easy cases.

For the 20 cases with N ≥ 10 the picture differs (supplement Table S-X-multiplicity; macro `\NCmPSOVNSvsPSOLargeDL`):
- pooled ΔL = −0.079 pp, which exceeds the ±0.05 margin;
- unadjusted p = 0.048, Holm p = 0.29.

So on the practically relevant layouts the comparison is **inconclusive**: neither significant nor shown equivalent. In addition, the margin (0.05) was set after the primary analysis and lies just above the minimal margin (0.041). A reviewer will read this as tailored, despite the transparency.

Fix, using existing data only:
- (a) Report the 90% CI and the minimal margin of PSO-VNS − PSO for N ≥ 10, pooled and per data set. Add them to Table S-equivalence and add a macro. Then add one sentence at 06_results:30, after "(Table S13)": "Over the 20 cases with $N\ge10$, the mean difference is \NCmPSOVNSvsPSOLargeDL~pp (90\% CI \TBD{macro: N≥10 CI}), neither significant after Holm correction nor equivalent at the margin; the 68-case equivalence is thus carried mainly by the \TBD{macro: 48} cases with $N<10$, in which both methods reach nearly the same wake loss."
- (b) Justify the margin in energy terms where it is introduced (05:15): "0.05 pp is 0.05% of AEP, well below the wake-model uncertainty (\NXModelShiftMean~pp between Jensen and Gaussian)." Then call equivalence a secondary, post hoc analysis in the abstract (done in the issue 1 text) and in the contributions.

**5. [major] Generality of M1: the rank reversal is demonstrated only on the authors' own earlier pool, and its mechanism is confounded with the boundary handling**

Locations: `02_intro.tex:7`, `04_methods.tex:19`, `06_results.tex:9–16`, Table I caption `analysis/mpce_tab_baseline.tex:4`, `05_setup.tex:10,31`.

Problem (a): no evidence of prevalence. The intro gives no WFLOP paper that uses $w=0.7$, $c_1=c_2=2$. The Table I pool is "the Method Pool of the Old-Setting Comparison", and `05_setup.tex:10,31` twice mention "the earlier version". Readers cannot see that earlier version, and M1 reads as a correction of the authors' own previous study. That is still a legitimate case study, but it must be framed that way or backed by evidence of use.

Problem (b): the mechanism is confounded. Clipped coordinates keep their velocity: no reflection or damping (`04_methods.tex:19`). The penalty mixes m² and m (`03_model.tex:38`). Section IX itself shows that the boundary rule alone changes SSA's results significantly in 6/6 cases. The observed collapse (8.7% of coordinates clipped per iteration, 2.6% feasible evaluations) is therefore caused jointly by order-2 instability *and* this boundary handling, but the text attributes it to the stability condition alone ("violates … so the swarm does not contract").

Fix:
- (i) Cite two to four verified WFLOP studies that use $c_1=c_2=2$ with $w\approx0.7$, without comment. If none can be verified, write at 02_intro:7: "…a setting that is easy to adopt from the original PSO papers and that was used in our own earlier comparison, lies in the order-1 … region".
- (ii) Delete "(where logged, the best initial objectives agree, also with the earlier version)" at 05:10. At 05:31 write "reproduces the result distributions of our earlier implementation".
- (iii) Change the Table I caption to "…on an Eight-Method Pool (the Methods of Table~II with LX-SSA-VNS Instead of PSO-VNS)…".
- (iv) At 06_results:9, after "Instrumented reruns show the mechanism", add "(with our boundary handling, in which a clipped coordinate keeps its velocity)". Add one sentence to the Limitations: "The effect of the old setting was measured with one boundary rule (clipping without velocity reset) and one penalty scaling; with velocity damping or reflection it may be smaller."

**6. [major] The key control RSD-VNS is not defined in Methods or Protocol; RS-VNS, a control the paper calls weak, gets a proposition instead**

Locations: `04_methods.tex:59–63`, `05_setup.tex:7`, `07_ablation.tex:14`.

Problem. The abstract, the intro and M2 rest on "random sampling in the farm" (RSD-VNS). The template (04:59) and the protocol (05:7: "adds LX-SSA-VNS and the no-search control RS-VNS") mention only RS-VNS. RSD-VNS first appears in contribution 3 and is defined only in Section VII. The supplement (line 163) holds an important detail: the first 30 RSD-VNS layouts are the common square-sampled initial population. That detail is not in the main text.

Meanwhile Proposition 2 and about 8 lines of IV-A are spent proving that the weak control is weak. RSD-VNS was also not run in the budget, feasible-start, Horns Rev and IEA37 studies, and this is not stated.

Fix. In 04:59, replace "…in the control \emph{RS-VNS}. VNS alone…" with:

> "…in the controls \emph{RSD-VNS} (uniform in the farm disc) and \emph{RS-VNS} (uniform in the bounding square), each after the common initial population. VNS alone…"

Then shorten 04:60–63 (Proposition 2 and the sentence after it) to:

> "A square sample lies inside the farm only with probability $(\pi/4)^N$ (Proposition~\ref{S-prop:S-box}), so RS-VNS is a weak, no-search minimum bar and RSD-VNS the relevant control."

Move the proposition to the supplement; the contribution 3 and Section VII references then point to S-prop:S-box. In 05:7, write "…adds LX-SSA-VNS and the controls RSD-VNS and RS-VNS (Section~\ref{sec:template}), whose 3,015 random layouts…". In 05:25, add "(RSD-VNS: benchmark only)".

**7. [major] Text and Table IV report different numbers for the same contrasts**

Location: `07_ablation.tex:8,14,24` vs `analysis/mpce_tab_ablation.tex`.

| Contrast | ΔL, text | ΔL, Table IV |
|---|---|---|
| SSA-VNS vs SSA | −0.328 | −0.294 |
| LX-SSA-VNS vs LX-SSA | −0.431 | −0.394 |
| LX-SSA vs SSA | +0.166 | +0.153 |

The text uses the `\NAbl…DL` macros; the table matches Table S-X-multiplicity (n = 66). In addition, the text quotes unadjusted p (e.g., 1.4×10⁻⁶ for LX-SSA-VNS vs RSD-VNS) next to a table showing Holm p (8.2×10⁻⁶), and the intro quotes the unadjusted p = 5.4×10⁻⁶ as the headline for the Laplace claim. Reviewers will read this as an error.

Fix. Use one ΔL definition throughout: the case means over cases where both variants qualify, as in the table. Either switch the three text macros to the `\NCm…DL` versions or state in the Table IV caption how the two differ.

Also quote Holm-adjusted p in the text and intro to match the table: change 07:5 to "Here $p$ is the Holm-adjusted $p_W$ of Table~\ref{tab:ablation}". Then use `\NCmLXSSAVNSvsRSDVNSPHolm`, `\NCmSSAVNSvsLXSSAVNSPHolm` (if not present, add macros) and `\NCmLXSSAvsSSAPHolm`.

**8. [major] Practical scope: the ranking is a small-budget phenomenon, and one Horns Rev fact is hidden**

Locations: `08_beyond.tex:8,35`, `10_limits_concl.tex:30,38`.

Problem. Three results limit the practical reach of the 6,030-evaluation ranking:
- At 120,030 evaluations five methods come within 0.17 pp of PSO-VNS, and SSA-VNS overtakes PSO (so M1's "ahead of every salp-swarm method" is also budget-specific).
- On Horns Rev at 6,030 evaluations, **0 of 30 PSO-VNS runs beat the installed 7D grid**, although the optimizers may use 4D spacing (and have no cable or loads constraints). At 120k, DE has the best mean.
- On IEA37 36-turbine, the best PSO-VNS layout would rank 7th of 8 feasible submissions.

These results are reported, but only as local sentences, and the Conclusion does not mention the budget at all.

Fix:
- Conclusion (10:38), first sentence: "At a budget of 6,030 evaluations and from random starts, three findings emerge."
- Add to the Limitations (10:30): "At 6,030 evaluations no optimized Horns Rev~1 layout exceeds the installed $7D$ grid despite the looser $4D$ spacing, and at 120,030 evaluations the leading methods lie within \NBudCloseGapOneTwentyK~pp of each other: the rankings reported here are those of a small-budget regime."

**9. [major] Submission blockers still present**

Locations:
- `01_front.tex:4,5` (funding, corresponding author)
- `11_back.tex:2` (acknowledgment)
- `11_back.tex:5` (Zenodo DOI)
- `11_back.tex:144` (Solanki2023 volume/pages)
- `11_back.tex:203–215` (all four biographies)
- `analysis/mpce_tab_iea37.tex:25` ("[TBD: verify against Baker et al.]")
- `MPCE_PSO_VNS.tex:35` (journal header)

Red placeholders in the PDF mean desk return. The header says *Journal of Modern Power Systems and Clean Energy*: if the target is SWEVO or an IEEE energy journal, the class, header and nomenclature policy change. Biographies are not used in SWEVO, which would itself solve most of the page-13 overflow.

Fix. Resolve every `\TBD` and fix the target journal before the final page count. Budget about 25–30 lines for real MPCE/IEEE biographies.

### Minor

**10. [minor] Title** — `01_front.tex:1`; running head `MPCE_PSO_VNS.tex:36`.

The title is 20 words with an appended clause ("…, with a PSO–VNS Reference Method"). That clause invites "where is the novelty?" and conflicts with M3 and issue 3. The running head ("BASELINES AND CONTROLS …") does not match the title.

Fix, preferred (13 words): **"How Baseline Configuration and Random-Sampling Controls Change Metaheuristic Rankings in Wind Farm Layout Optimization"**. Alternative: "Baseline Configuration and Controls in Metaheuristic Comparisons for Wind Farm Layout Optimization". Running head: "SOLANKI et al.: BASELINE CONFIGURATION AND CONTROLS IN WIND FARM LAYOUT OPTIMIZATION".

**11. [minor] Keywords** — `01_front.tex:31`.

"Wake effect" adds little. "Variable neighborhood search" oversells (issue 3). In an energy journal, "stability" means power-system stability, so the paper should say "convergence condition" rather than "stability condition" in the abstract and Discussion.

Fix: "Benchmarking, equivalence testing, hybrid metaheuristics, particle swarm optimization, statistical comparison, wind farm layout optimization." Keep "order-2 stability" in Section III-A only, glossed once as "(variance convergence)".

**12. [minor] "Reference method" and "feasibility reliability" rest on one site**

Locations: `02_intro.tex:32`, `06_results.tex:39`, `10_limits_concl.tex:12,16,38`.

On the 68 cases PSO is also 100% feasible. The 30/30 vs 19/30 advantage is one 16-turbine parallelogram block with random starts (PSO: 28/30 at 30k; 30/30 with feasible starts).

Fix. Qualify each instance as "(Horns Rev~1 block, random starts)". In 10:16, change "PSO-VNS a more reliable one" to "PSO-VNS, which repairs residual infeasibility at no extra cost, a slightly more reliable one".

**13. [minor] IEA37: text says ten methods, table shows eight** — `05_setup.tex:29` ("all ten methods run on them") vs `analysis/mpce_tab_iea37.tex` (no LX-SSA-VNS, RS-VNS).

Fix: add the two rows or write "the eight methods of the main comparison".

**14. [minor] Duplicate mechanism text and a confusing counterfactual rank**

- `04_methods.tex:28`, the sentence from "In our runs (Fig.~\ref{S-fig:D-psodyn}), the old swarm…" to "…(Section~\ref{sec:psosetting})", repeats `06_results.tex:9` almost number for number.
- `06_results.tex:9` ends with "In Table~\ref{tab:friedman68}, the old setting would place PSO fifth (average rank 5.12)". This introduces a second "fifth" (5.12) next to Table I's 5.04. Readers will ask which is right.

Fix: delete both (these deletions are cut C4 and C5).

**15. [minor] Misleading DE sentence** — `04_methods.tex:28` ("Likewise, the DE scale factor 0.5 … exceeds Zaharie's critical value …").

"Likewise" suggests that DE is also misconfigured. In fact $F$ above Zaharie's critical value means the population variance does not collapse, i.e., DE is on the explorative side. That is exactly what a reviewer will connect with DE's 71.4% feasibility and last rank. It is the paper's own M1 logic applied to another baseline, and the paper leaves it unexamined.

Fix: "The DE setting ($F=0.5$, $\mathrm{CR}=0.9$) lies above Zaharie's critical value $0.14$~\cite{Zaharie2002}, i.e., on the explorative side; it is a standard, untuned default." Add to the Limitations: "DE's low feasibility (Table~\ref{tab:friedman68}) may reflect the same penalty-dominated search as the old PSO setting; we did not tune it."

**16. [minor] Flow: Fig. 1 and Section VI-C belong to the component analysis** — `06_results.tex:45–63`.

Fig. 1 shows the nine ablation variants (RSD-VNS, RS-VNS, LX-SSA-VNS) before they are introduced in Section VII. The VI-C "after the switch" paragraph is a VNS-phase result.

Fix: move the `figure*` block and 06:63 into 07 after the "VNS phase" paragraph. Keep 06:56–58 (quality/feasibility) as the end of VI-B, without a subsection heading.

**17. [minor] Tension between two statements on what PSO contributes**

- `07_ablation.tex:28`: "Whether PSO helps mainly by faster feasibility or by better exploration remains open."
- `08_beyond.tex:19`: "From random starts, PSO may thus help PSO-VNS mainly by reaching feasibility…"

Fix: at 07:28, write "…remains open here; the feasible-start study (Section~\ref{sec:feasinit}) suggests that feasibility is the main contribution."

**18. [minor] Proposition 3 formalizes a triviality** — `08_beyond.tex:16–18`.

"A trial replaces G only if it is feasible and better" follows directly from the penalty. Presenting it as a proposition invites "padding" comments, and the paper itself says the stall "is observed, not implied".

Fix: demote it to the supplement (see cut C7).

**19. [minor] Section structure** — `10_limits_concl.tex:1,28,36`; `04_methods.tex:1,54`.

The paper has twelve top-level sections. Section III ("Preliminaries") and Section IV ("Two-Phase Design") are each short, and Limitations is a separate section.

Fix:
- merge III and IV into "Methods" (saves one heading);
- make Limitations "X-B";
- fold VIII-D (two sentences, `08_beyond.tex:54–57`) into V-C "Verification and code".

**20. [minor] Captions too long for IEEE all-caps typesetting**

Location: `analysis/mpce_tab_friedman68.tex` (12 caps lines), `mpce_tab_ablation.tex` (11), `mpce_tab_wtl.tex` (8), `mpce_tab_hr16.tex` (8), `mpce_tab_baseline.tex` (7).

Fix. Keep captions to a noun phrase of one or two lines (e.g., "Case-Level Analysis of the Eight Methods over the 68 Benchmark Cases"). Move the definitions to a `\multicolumn` note in `\scriptsize` below the table, as Table I already does. This is a pipeline change in `mpce_results.py` (C12 below).

**21. [minor] Nomenclature** — `01_front.tex:34–50`.

It has twelve entries (about 20 lines across the two columns of p. 1). Most symbols are defined at first use, and the section pushes the introduction down.

Fix: delete it if the target journal does not require one (SWEVO and IEEE Trans. do not). Otherwise reduce it to $B$, $F_p$, $N$, $N_p$, $w,c_1,c_2$, $\omega$.

**22. [minor] Intro literature paragraph** — `02_intro.tex:5`.

This is a list of methods with no stated relevance. It is acceptable for MPCE but weak for SWEVO.

Fix (optional): cut to "…has been solved with many metaheuristics~\cite{…}, with gradient-based and pseudo-gradient methods~\cite{…} and with surrogates~\cite{…}". Keep one surrogate citation instead of three; this drops two references (saves about 6 reference lines).

---

## C. Figures and tables

- Keep in the main text:
  - Table I: carries M1;
  - Table II: the overall comparison;
  - Table IV: carries M2, but shorten the caption;
  - Table VI: budget/feasible starts;
  - Table VII: Horns Rev;
  - Table VIII: IEA37;
  - Algorithm 1;
  - Fig. 1, reduced as described below.
- Move to the supplement:
  - Table V (split): it tests a secondary parameter, and every number needed is in the text;
  - Table IX (robust-final): the text already states τ̄, "same best" and the Gaussian flip, and the supplement already has the full Table S-robust.
- Merge Table III into Table II. The per-cluster W/T/L rows are robustness detail that already lives in Tables S-X-cluster/S-X-loco. Add two columns to Table II: "W/T/L" (All row) and "$\tilde r_{\rm rb}$".
- Fig. 1: one row (the three Data Set II panels, where the differences show), or height reduced to about 60%. Report that the full six panels are in the supplement.
- M1 has no figure. Do not add one given the page limit; Table I is clear.

---

## D. Cuts to reach 12 pages (main text → supplement), with locations and savings

Tier 1 alone gives about 0.7 page; with Tier 2 the total is about 1.6 pages of candidates. **Minimum recommended: all of Tier 1 plus C6.** That is about 105 lines, about 0.9 page. It leaves room for real biographies, the funding line and float slack, and may bring the paper to 11 pages plus references.

### Tier 1: pure de-duplication, no evidence lost (about 80 lines)

**C1. Contributions list** — `02_intro.tex:12–39`. It currently spans p. 2 col. 1 l. 419 to col. 2 l. 465, about 60 lines, with about 25 numbers that are repeated in Sections VI–X. Replace lines 12–39 (keep the CHECK comments or move them) with the text below. **Saves about 35 lines.**

```latex
Our contributions are:
\begin{enumerate}
\item \emph{Protocol.} Established benchmarking practice~\cite{Derrac2011,BartzBeielstein2020,LaTorre2021} applied to the continuous WFLOP (equal budgets counting gradient evaluations, seed-paired runs, case-mean and run-level tests with Holm correction, public code and data), with three elements not reported for this problem to our knowledge: feasibility-aware ranking, equal-cost random-sampling controls for hybrids, and a (post hoc) equivalence analysis~\cite{Lakens2017,Benavoli2017}; applied to 68 benchmark cases, a Horns Rev~1 block and IEA37 Case Study~1 (Sections~\ref{sec:setup}, \ref{sec:beyond}).
\item \emph{Baseline configuration.} The PSO setting $w=0.7$, $c_1=c_2=2$ violates the order-2 condition (Proposition~\ref{prop:pso-stability}); with our boundary handling its swarm does not contract and rarely samples feasible layouts, and it reverses the ranking of PSO and the salp-swarm hybrids (Table~\ref{tab:baseline}; Section~\ref{sec:psosetting}).
\item \emph{Component analysis.} Against random sampling in the farm (RSD-VNS), salp-swarm first phases add nothing and LX-SSA as published is worse, whereas a PSO phase gives a clear gain; at 6,030 evaluations the VNS phase is one compass-search descent (Section~\ref{sec:ablation}).
\item \emph{PSO-VNS versus PSO.} PSO-VNS has the best average rank and is feasible in every run, but is practically equivalent to a well-configured PSO on average; we report where it gains (Horns Rev~1 feasibility, larger Data Set~II layouts) and where it does not (feasible starts, Gaussian re-evaluation, IEA37), and derive seven benchmarking recommendations (Sections~\ref{sec:overall}, \ref{sec:beyond}--\ref{sec:discussion}).
\end{enumerate}
```

**C2. Table IX → supplement** — `09_robust.tex:5` (`\input{analysis/mpce_tab_robust_final.tex}`). About 14 lines, p. 10 col. 1 y 270–430. Move the `\input` into the supplement next to Table S-robust. In 09:7, change "(Table~\ref{tab:robust-final})" to "(Table~\ref{S-tab:robust})". Add "(\NCubicOverI\%/\NCubicOverII\% with the linear curve; similar with a 25~m/s cut-out)" if wanted. **Saves about 14 lines.**

**C3. Discussion paragraphs 1–3** — `10_limits_concl.tex:3–14`. About 32 lines (p. 10 col. 2 l. 445 to p. 11 col. 1 l. 140). They restate VI-A, VII and VI-B with the same numbers. Replace 10:3, 10:7 and 10:12 with this one paragraph. **Saves about 22 lines.**

```latex
In this benchmark, which method appears best depends on the protocol as much as on the methods. An advantage over PSO is informative only if the PSO baseline is convergent (Section~\ref{sec:psosetting}); an advantage of a hybrid's first phase only if it beats random sampling \emph{in the site} at equal cost, as a square-sampling control is too weak (Section~\ref{sec:ablation}); and the absence of a significant difference is not equivalence. At 6,030 evaluations and from random starts, a constriction PSO and PSO-VNS are interchangeable on average, PSO-VNS repairing residual infeasibility; from feasible starts a single local search (VNS alone) is competitive. For WFLOP comparisons we recommend~\cite{BartzBeielstein2020,LaTorre2021}:
```

The current 10:16 lead-in sentence is folded into the paragraph above; delete it.

**C4. Duplicate diagnostics in III-A** — `04_methods.tex:28`. Delete the passage from "In our runs (Fig.~\ref{S-fig:D-psodyn}), the old swarm does not contract…" to "…the setting changes the ranking (Section~\ref{sec:psosetting})." Keep the preceding sentences. Replace the "Likewise, the DE…" sentence with the issue-15 wording. **Saves about 6 lines.**

**C5. Counterfactual rank sentence** — `06_results.tex:9`, the last sentence: "In Table~\ref{tab:friedman68}, the old setting would place PSO \NOldPSOPos{} (average rank \NOldPSORank), not \NPSOPos." Delete it (issue 14). **Saves about 2 lines.**

### Tier 2: move secondary material (about 100 lines available; take C6 at least)

**C6. Table V (budget split) → supplement.** Location: `07_ablation.tex:31–32`. The table is about 15 lines (p. 8 col. 2 y 357–540); the paragraph is about 17 lines. Move `\input{analysis/mpce_tab_split.tex}` to the supplement, before Table S-split-cases, and replace 07:32 with the text below. **Saves about 26 lines.**

```latex
\emph{Budget split.}\label{sec:split} On \NSplitCases{} cases (middle and largest $N$ of each farm; Tables~\ref{S-tab:split} and~\ref{S-tab:split-cases}), $\omega=0.75$ has a lower mean wake loss than the preset $\omega=0.5$ in \NSplitSeventyFiveLower{} cases (significantly in \NSplitSeventyFiveSig; \NSplitGainSeventyFive~pp on average) and is never significantly worse; $\omega=0.25$ is worse in \NSplitTwentyFiveSig{} cases, and plain PSO ($\omega=1$) is worse than $\omega=0.75$ (W/T/L of $\omega=0.75$: \NSplitHundredVsSeventyFive), so the best tested split is not the PSO end point. We keep the preset $\omega=0.5$, which may slightly understate PSO-VNS; a larger split should be confirmed on the target problem.
```

**C7. Proposition 3 and the stall mechanism** — `08_beyond.tex:14–19`. About 22 lines, p. 9 col. 1 y 406–700. Move Proposition 3 and the relabelling detail to S-sec:S-th-feas, and replace 08:14–19 with the text below. **Saves about 12 lines.**

```latex
\emph{Feasible starts.} With feasibility-preserving initialization, \NFeasInitPhrase. \NFbFeasBestRank{} alone is then best in mean wake loss and average rank; PSO-VNS follows (rank \NFbFeasRankPSOVNS, loss \NFbFeasLossPSOVNS\%; \NFbRandLossPSOVNS\% with random starts; Table~\ref{S-tab:feasinit-cases}). PSO and DE never improve on the best initial layout (\NDFsStoredRunsPSO{} runs each): their moves change all coordinates at once and, between distinct dense layouts, make turbines collide (median minimum spacing after the first PSO move \NDFsFirstMinSpPSO~m, limit \NDSmin~m), so no trial is both feasible and better and the global best never moves (Section~\ref{S-sec:S-th-feas}, Table~\ref{S-tab:D-feasstart}). VNS, which moves one coordinate at a time, is not affected. The advantage of PSO-VNS thus depends on the initialization: from random starts, PSO may help mainly by reaching feasibility.
```

**C8. Table III → two columns of Table II.** Location: `06_results.tex:21` (`\input{analysis/mpce_tab_wtl.tex}`). Table III is about 20 lines (p. 6 col. 2 y 320–600). Add "W/T/L" (All row) and "$\tilde r_{\rm rb}$" to `mpce_tab_friedman68.tex`, and move the per-cluster table to the supplement next to S-X-cluster. This is a pipeline change in `mpce_results.py`. Update the references at 06:30 and in the Table IV caption ("a different family from Table~\ref{tab:friedman68}"). **Saves about 14 lines.**

**C9. Nomenclature** — `01_front.tex:34–50`. Delete it, or reduce it to six entries (issue 21). **Saves 12–20 lines.**

**C10. Post hoc Data Set II exploration** — `06_results.tex:30`. Replace the passage from "The advantage of PSO-VNS is concentrated…" to "…($p=\NXHolmPSOVNSvsPSOLargeP$)." with the text below. The permutation, Spearman and Data Set I detail go to Table S-X-subgroup. **Saves about 6 lines, and folds in the issue-4 sentence.**

```latex
In a post hoc split, the advantage of PSO-VNS is concentrated in the larger Data Set~II layouts ($N\ge10$: lower mean wake loss in \NCmPSOVNSvsPSOdsIILargeWins{} of \NCmPSOVNSvsPSOdsIILargeN{} cases, $p_{\rm Holm}=\NXHolmPSOVNSvsPSOdsIILargeP$ over all \NXMultN{} case-mean tests; Table~\ref{S-tab:X-subgroup}); Data Set~I shows no such effect, and over all 20 cases with $N\ge10$ the difference (\NCmPSOVNSvsPSOLargeDL~pp) is neither significant after Holm correction nor shown equivalent.
```

**C11. Verification and code** — `05_setup.tex:31`. Keep the first sentence (objective and IEA37 calculator match) and the sentence on public code. Move the Mann–Whitney re-validation, the bit-match count and the hardware/software versions to Section S-I or the Reproducibility Map. **Saves about 5 lines.**

**C12. Caption slimming** (issue 20): Tables II, IV, VII, I. **Saves about 12 lines.**

**C13. Fig. 1 to one row** — `06_results.tex:49`. Show the Data Set II panels only, and point to the supplement for all six. The figure is about 270 pt at full width; one row is about 150 pt. **Saves about 20 column lines.**

**C14. Structure merges** (issue 19): III + IV into one section; Limitations as a subsection; VIII-D folded into V-C. **Saves about 6 lines.**

**Suggested package:**

| Cut | Lines saved |
|---|---|
| C1 | 35 |
| C2 | 14 |
| C3 | 22 |
| C4 | 6 |
| C5 | 2 |
| C6 | 26 |
| **Total** | **≈105 (≈0.9 page)** |

If floats still push the references past p. 12, add C7 (12) and C8 (14). None of these cuts removes a claim; every moved item already has a numbered home in the supplement.

---

## E. Overall verdict

**Major revision before external review.** The science is careful and the tone is honest. The WFLOP-specific lessons are worth publishing:
- a non-convergent PSO baseline plus clipping collapses feasibility and flips rankings;
- salp-swarm first phases add nothing beyond in-site random sampling;
- equivalence, not superiority, is the right claim for PSO-VNS vs PSO.

But the headline text still claims more than the evidence supports:
- the Laplace step (issue 2);
- cluster robustness (issue 1);
- equivalence on a benchmark dominated by easy cases (issue 4);
- a "VNS" that is a single descent (issue 3).

The key control is also missing from the Methods (issue 6), the text and table numbers conflict (issue 7), and the paper is one page too long through repetition, not through evidence. The fixes are textual or use existing data; no new runs are needed.

## The 5 most important fixes

1. **Rewrite the abstract** (issue 1). Use the same-pool wording; replace "Laplace step harmful" with "LX-SSA as published"; drop "farm clusters" from the robustness sentence; mark the margin and the Data Set II subgroup as post hoc; state "6,030 evaluations, random starts".
2. **Qualify the PSO-VNS ≈ PSO equivalence** (issue 4). Report the N ≥ 10 CI and minimal margin (pooled ΔL −0.079 pp is outside ±0.05 pp, and Holm p = 0.29, so the result is inconclusive there). Justify the margin in AEP terms, and label the equivalence analysis post hoc/secondary.
3. **Remove the Laplace-step causal claim everywhere** (issue 2) and **define RSD-VNS in Sections IV-A and V-A** (issue 6). Demote Proposition 2 about the weak square control.
4. **Harmonize the numbers in Section VII and Table IV** (issue 7): the same ΔL definition, and Holm-adjusted p in the text and intro. Fix the IEA37 "ten methods" statement (issue 13).
5. **Cut about 105 lines**: C1 (contributions), C2 (Table IX), C3 (Discussion restatement), C4/C5 (duplicates), C6 (Table V). Resolve all `\TBD` and fix the target journal (issue 9) before the final page count. Retitle without the PSO–VNS clause (issue 10).

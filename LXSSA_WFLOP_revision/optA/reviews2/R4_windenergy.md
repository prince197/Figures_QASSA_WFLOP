# R4 (round 2): Wind-energy review of the manuscript (wake modeling, WFLOP, Horns Rev, IEA Task 37)

Manuscript: MPCE_PSO_VNS.tex + optA/*.tex (macros analysis/mpce_numbers*.tex), supplement MPCE_PSO_VNS_supplement.tex.
Code checked: analysis/wflop_model.py and analysis/authors_objective.py / objective_original.py (Kusiak–Song objective),
analysis/hornsrev_model.py (after commit 7676da9, the direction-binning fix), analysis/iea37_model.py,
analysis/mpce_results.py (Gaussian re-evaluation, HR validation), analysis/mpce_inference_extra.py (energy, model shift).

## 0. What I checked

The reviewer did not edit any repository file. All scripts ran in a scratch directory.

1. **Arithmetic and energy translation.** E_ideal = 14,045.74 gives 936.4 kW, a capacity factor (CF) of 62.4 % and 8,202.7 MWh/yr for Data Set I. For Data Set II the values are 7,315.38, 487.7 kW, CF 32.5 % and 4,272 MWh/yr. The factor 8.76/15 is correct. The code wake loss matches the macro (PSO-VNS 1.322 %).
2. **IEA37.** `iea37_model.aep` reproduces every published AEP in `iea37_published_results.csv` to the MWh, so the calculator is correct.
3. **Horns Rev (HR) model re-evaluation.** I re-evaluated all 970 hrfix runs, i.e. all feasible final layouts plus the installed block, under four settings:
   - the paper's 5° bins centred at 2.5°;
   - 5° bins centred at 0°;
   - 1° bins, both at half-degree and at integer-degree centres;
   - **PyWake 2.6.20**, which I installed in scratch, with three models: `NOJ(k=0.04)` (default rotor averaging), `BastankhahGaussian(k=0.0324555)` and TurbOPark.
4. **Benchmark model re-evaluation.** I re-evaluated the 68-case final layouts of PSO-VNS, PSO, VNS, SSA-VNS and MS-SLSQP with each 15° Kusiak–Song bin split into 15 sub-bins of 1°. Each sub-bin keeps its bin's Weibull parameters and gets 1/15 of its frequency, i.e. a piecewise-constant rose.

The verification results drive most of the issues below.

---

## 1. Numbered issues

### Issue 1 (CRITICAL): the Horns Rev "runs above installed" and "highest mean AEP at 6,030" results come from the 5° direction discretization

**Location**
- optA/08_beyond.tex:35. "PSO-VNS has the highest mean AEP … significantly better than every other method … 3 of its 30 runs exceed the installed layout" at 30,030 evaluations; 2 runs at 120,030 via \NHRAbove….
- analysis/mpce_tab_hr16.tex.
- analysis/hornsrev_model.py:33–37.

**Finding.** The optimized layouts gain about 1 GWh/yr, i.e. about 0.7 pp of AEP, from the position of the 5° direction bins. The installed grid gains only 0.32 GWh/yr from it.

Mean AEP in GWh/yr of feasible runs. Settings: random start. Counts in brackets are runs that exceed the installed block under the same model.

| Setting (budget) | This paper: 5° @2.5° | 1° bins | PyWake NOJ k=0.04 | PyWake Bastankhah | PyWake TurbOPark |
|---|---|---|---|---|---|
| Installed block | 139.82 | 139.50 | 139.17 | 140.90 | 136.23 |
| PSO-VNS (6,030) | 138.29 [0] | 136.99 [0] | 137.22 [0] | 138.60 [0] | 133.79 [0] |
| PSO (6,030, 19 feasible) | 137.79 [0] | **137.04** [0] | **137.31** [0] | **138.71** [0] | **133.89** [0] |
| PSO-VNS (30,030) | 139.11 [**3**] | 137.71 [0] | 137.85 [0] | 139.30 [0] | 134.54 [0] |
| DE (120,030) | 139.59 [**4**] | 138.23 [0] | 138.33 [1] | 139.86 [0] | 135.06 [0] |

**Consequences**
- **Runs above installed.** No PSO-VNS run at any budget beats the installed layout under a 1° rose or any PyWake model. Of all 970 feasible optimized layouts, only one DE run under PyWake NOJ beats it.
- **Highest mean AEP at 6,030.** On the 19 seeds where both methods are feasible, the paper's model gives PSO-VNS an advantage over PSO of +0.54 GWh/yr (Wilcoxon p = 0.016). The advantage reverses under the other models:

  | Model | PSO-VNS − PSO (GWh/yr) | p |
  |---|---|---|
  | 1° bins | −0.06 | 0.65 |
  | PyWake NOJ | −0.12 | 0.40 |
  | PyWake Bastankhah | −0.14 | 0.37 |
  | PyWake TurbOPark | −0.12 | 0.44 |

- **What survives.** Under every model the only robust result at 6,030 is feasibility: 30/30 vs 19/30. The claim of the highest mean AEP at 30,030 survives all models (PyWake NOJ: PSO-VNS 137.85, PSO 137.65, DE 137.55).
- **Why.** The VNS-based methods (compass search) lose the most under the finer rose: about 1.4 GWh/yr against 0.5–0.8 for MS-SLSQP, PSO and DE. They exploit the coarse rose best.
- **The bin phase alone moves the reference.** With the same 5° bins centred at 0°, the installed block's wake loss is 6.41 %, against 6.04 % at 2.5°. The rows of the grid are aligned with 270°. So the "fix" did two things: it normalized the frequencies and it shifted the bin phase, and the shift is not neutral for the reference layout.

**Fix**
1. Re-evaluate every HR layout and the installed block at 1° resolution. The best option is PyWake itself, since the data come from it.
2. Remove "N of its 30 runs exceed the installed layout" at all budgets, or report it only together with the finer-rose result (0/30).
3. Restate the 6,030 result as a feasibility result.
4. Ideally, re-run HR with 1° bins or randomly phased bins so that optimizers cannot tune to the bin centres. HR costs are small.
5. Add one sentence to Limitations: under a 1° rose, no optimized layout at any budget reaches the installed grid.

### Issue 2 (MAJOR): the PyWake validation of the HR model uses a stale constant, and the stated cause of the difference is wrong

**Location**
- optA/08_beyond.tex:35 ("our AEP agrees with PyWake within \NHRPyWakeDiff\% = 0.64 %").
- MPCE_PSO_VNS_supplement.tex:349 ("… with the same turbine curves … The remaining difference is expected from … hub-centre in-wake test instead of a rotor-averaged deficit").
- analysis/mpce_results.py:166–171 (`PYWAKE_HR80_AEP = 662.5`, "PyWake is not installed in this container, so the value is a constant").

**Finding**
- **The reference value is not reproduced.** The 662.5 GWh/yr / 10.96 % constant comes from an earlier text (hornsrev_site_text.tex), whose own model gave 10.95 % before the fix. PyWake 2.6.20 `NOJ(site, V80, k=0.04)` for the 80-turbine farm gives:

  | Direction bins | PyWake AEP (GWh/yr) | PyWake wake loss | Ours (GWh/yr) | Our wake loss |
  |---|---|---|---|---|
  | PyWake default, 1° | 661.93 | 11.04 % | — | — |
  | Same 5° bins as ours | 661.39 | 11.11 % | 666.75 | 10.39 % |

- **Rotor averaging does not explain it.** The PyWake hub-centre variant (`RotorCenter`) at the same bins gives 659.16 GWh/yr (11.41 %). It is further from our model, not closer.
- **Likely causes.** The remaining 0.7–1.0 pp of wake loss (7–9 % of the wake loss) comes from:
  - the thrust coefficient taken at the free-stream speed (hornsrev_model.py:42, 96); V80 C_T drops sharply between 10 and 13 m/s;
  - the bin phase: our model with 1° integer-degree bins gives 663.33 GWh/yr.
- **Wake-free AEP.** It matches PyWake exactly (744.04 GWh/yr). All the discrepancy is in the wake loss, so quoting an AEP difference (0.64 %) understates it: it is about 6 % of the wake loss.

**Fix**
- Run PyWake, pinned in requirements.txt, and report the wake loss and AEP at identical direction bins.
- Correct the attribution in the supplement.
- Either compute C_T at the effective inflow speed (propagate downwind), as PyWake does, or state it as a known bias that lowers the wake loss.
- Remove "with the same settings", or make it true.

### Issue 3 (MAJOR): the 15° Kusiak–Song rose has blind wake gaps, and the optimized layouts exploit them by more than 20 times the equivalence margin

**Location**
- optA/10_limits_concl.tex:30 ("coarse direction bins that optimized layouts can exploit"; stated but not quantified).
- optA/09_robust.tex:7.
- optA/05_setup.tex:15.

**Finding: the blind gaps**
- With K = 0.075 and the cone vertex R/K = 513 m upstream, a turbine 7.5° off every bin centre is never waked in any of the 24 bins once it is at least 686 m (8.9 D) away. I checked this with `wflop_model.jensen_deficits`.
- In the 750 m and 1000 m farms, many turbine pairs can sit in these gaps.

**Finding: re-evaluation with a 1° piecewise-constant rose (no re-optimization)**
- **Wake losses rise sharply.** Mean wake loss rises from 1.32 % to 2.34 % for PSO-VNS and from 1.34 % to 2.38 % for PSO. That is +1.0 pp, or +77 %, about 20 times the ±0.05 pp margin. By data set it rises from 0.73 % to 1.55 % (DS I) and from 2.22 % to 3.24 % (DS II).
- **Sign changes.** PSO-VNS − PSO changes sign in 26 of 68 cases.
- **Equivalence no longer holds.** The case-mean difference becomes −0.043 pp. The 90 % bootstrap CI is [−0.063, −0.025] and the Wilcoxon p is below 0.001. PSO-VNS becomes significantly better, and the ±0.05 pp equivalence is lost.
- **The MS-SLSQP gap was mostly exploitation.** MS-SLSQP exploits the bins least. Its gap to PSO-VNS shrinks from 0.44 pp to 0.13 pp, and its average rank among the five methods I re-evaluated improves from 4.49 to 3.21.
- **Robust findings.** The Data Set II N ≥ 10 advantage of PSO-VNS survives (10 of 10 cases, −0.16 pp). PSO-VNS keeps the best average rank (1.75 against 2.10 for PSO), which is to the authors' credit.

**Fix**
- Add a "direction-resolution" row to Table tab:robust-final: 15 sub-bins per sector, re-evaluated.
- Report the TOST result and the PSO-VNS − PSO sign changes under it.
- Move the statement from Limitations into Section 9 with numbers.
- State that the ranking is a property of the discretized benchmark, not of a continuous rose.

### Issue 4 (MAJOR): the ±0.05 pp equivalence margin is an absolute margin diluted by trivial cases, and is smaller than every model-form uncertainty the paper itself measures

**Location**
- optA/01_front.tex:13 (abstract).
- optA/02_intro.tex:14 and :32.
- optA/05_setup.tex:15 (the margin, and "the PSO-VNS−PSO difference changes by **only** \NXModelShiftPairMean = 0.031 pp").
- optA/06_results.tex:30.
- optA/10_limits_concl.tex:12.

**Finding: dilution by trivial cases**
- 15 of the 68 cases have a PSO-VNS wake loss below 0.05 pp, and 24 below 0.2 pp. These are small N in large discs, where any two methods are trivially "equivalent".
- On the 44 cases with a PSO-VNS wake loss of at least 0.2 pp, the 90 % bootstrap CI of PSO-VNS − PSO is [−0.060, +0.008] pp. That is outside the margin, so equivalence is not shown where it matters.
- Per case, the difference reaches 0.41 pp (DS II, 500 m, N = 10).

**Finding: model sensitivity**
- The Jensen-to-Gaussian shift of the per-case PSO-VNS − PSO difference is 0.031 pp. That is 62 % of the margin and 1.7 times the headline mean difference (−0.018). It also flips the sign in \NXModelShiftPairSignChanges = 18 of 68 cases. "Only" is not warranted.
- Level shifts are 0.38 pp from Jensen to Gaussian and about 1 pp from 15° to 1° bins (Issue 3). Either one is larger than the margin by a factor of 7 to 20.

**Fix**
1. Pre-specify, or at least justify, the margin in energy and uncertainty terms. For example, use a relative margin of 2–5 % of the case's wake loss, or 0.05 % of AEP, which is ≈4.1 MWh/turbine/yr for DS I and ≈2.1 for DS II.
2. Report TOST on the non-trivial cases, and under the 1° and Gaussian re-evaluations.
3. Drop "only" at 05_setup.tex:15 and report the 18 sign changes in the main text.
4. Rephrase "practically equivalent on average" as "equivalent within ±0.05 pp on the 68-case benchmark average, a margin below the model-form uncertainty".

### Issue 5 (MAJOR): the IEA37 gap statements depend on a 1 mm feasibility tolerance and understate the gap

**Location**
- optA/08_beyond.tex:50.
- analysis/mpce_tab_iea37.tex:25 (footnote, with a `\TBD{verify against Baker et al. (2019)}` still present).
- optA/10_limits_concl.tex:34 (\NIEAPhrase "98.7 % / 95.8 %").
- analysis/iea37_published_results.csv (column `Feasible_tol1e-3m`).

**Finding: the tolerance decides which published layout is "best feasible"**
- Published layouts are judged at 1 mm, our own at 10⁻⁶ m.
- The 36-turbine submission par12 (882,383 MWh) lies only **4.9 mm** outside a 2,000 m circle. After radial projection it is feasible: minimum spacing 596 m, AEP **882,382.8 MWh**.
- The 16-turbine par12 lies 3.5 m outside. Projected, it is feasible with AEP **421,451.1 MWh**, which is above par4 (418,924.4).
- I verified both with iea37_model.

**Finding: the gaps are larger with the projected layouts**

| | Paper | With the projected layouts |
|---|---|---|
| Gap of best PSO-VNS, 16 turbines | −1.33 % | −1.93 % |
| Gap of best PSO-VNS, 36 turbines | −4.20 % | −6.23 % |
| Gap of best of our methods (MS-SLSQP), 36 turbines | — | −5.56 % |
| Rank of best PSO-VNS, 16 turbines | third of 9 | about fourth |
| Rank of best PSO-VNS, 36 turbines | seventh of 8 | about ninth |

**Finding: AEP ratios hide the wake-loss difference**
- 16 turbines: 11.97 % wake loss for PSO-VNS against 10.78 % (par4) or 10.24 % (par12 projected).
- 36 turbines: 21.68 % against 18.25 % or 16.48 %. That is 3.4–5.2 pp, or 19–32 % more wake loss.
- Task 37 readers compare wake loss, not AEP ratios.

**Fix**
1. Report both conventions: strict, and tolerance-projected (radial projection plus a spacing check).
2. Give the gaps in AEP % and in wake-loss pp.
3. Resolve the TBD from Baker et al. (2019) and the attribution "Participant 4 = SNOPT+WEC", which is currently "assumed".
4. State the budget context from Thomas et al. (2023), who report optimizer effort on a closely related IEA37 problem.

### Issue 6 (MAJOR): the energy translation is arithmetically correct but lacks practical context, and one intro sentence overstates it

**Location**
- optA/02_intro.tex:3 ("a difference of a few tenths of a percentage point is a lifetime energy difference of that size").
- optA/03_model.tex:25.
- optA/06_results.tex:30 (\NXEnPSOVNSvsPSOPct = 0.020 % of AEP).
- MPCE_PSO_VNS_supplement.tex:326.

**Finding**
- The 68-case gain of PSO-VNS over PSO is 6 MWh/yr per case, or 0.2 MWh per turbine per year. At €50/MWh that is about €10 per turbine per year: economically nil.
- The DS II N ≥ 10 gain (0.22 % of AEP, 111 MWh/yr per case) is modest but real within the model, and robust to the finer rose (Issue 3).
- The paper states the DS I capacity factor (62 %) but not DS II (32.5 %, which is realistic).
- It does not relate these magnitudes to the wake-loss uncertainty of engineering models, which is typically tens of percent of the wake loss, nor to its own measured model shifts.
- The intro sentence is true only inside the model.

**Fix**
- Add one short paragraph or table row with the following:
  - the margin and headline gains in MWh per turbine per year for both data sets;
  - a € figure with a stated price;
  - the model-form shifts measured in the paper: 0.38 pp from Jensen to Gaussian, about 1 pp from the direction resolution, 0.7–1.0 pp for the HR model against PyWake.
- Qualify 02_intro.tex:3 with "in a given wake model".
- Add the DS II capacity factor.

### Issue 7 (MAJOR): the MS-SLSQP baseline says little about gradient-based layout optimization, which is the practitioner standard

**Location**
- optA/05_setup.tex:7.
- optA/10_limits_concl.tex:32.
- optA/06_results.tex:30 (PSO-VNS "significantly better than … MS-SLSQP in most cases").

**Finding**
- The benchmark objective, a top-hat hub-centre Jensen model on a 15° rose, is piecewise constant in the lateral positions. Finite-difference gradients are zero or undefined almost everywhere and cost 2N+1 evaluations each.
- On the smooth IEA37 objective, MS-SLSQP is the best of all methods at 16 turbines / 6,030 and 36 turbines / 30,030. That is the pattern a wind-energy reader expects: TOPFARM, FLORIS/FLOWFarm and SNOPT+WEC all use smooth models with analytic or AD gradients.
- Issue 3 shows that much of MS-SLSQP's benchmark deficit is exploitation by the other methods.

**Fix**
- Add one AD-gradient SLSQP or SNOPT run on IEA37 (JAX or autograd), with wake-expansion continuation if feasible and gradient cost accounted for.
- Otherwise, restrict the MS-SLSQP conclusions explicitly to "finite-difference SLSQP on a discontinuous objective".
- Drop any implication about gradient methods in general.

### Issue 8 (MODERATE): the 4D vs 7D caveat on Horns Rev is framed backwards, and the HR problem has no practical counterpart

**Location**
- optA/08_beyond.tex:35 ("may use the benchmark spacing 4D, whereas the installed grid keeps 7D (no run with 7D was made)").
- optA/05_setup.tex:27.

**Finding**
- All feasible HR layouts have a minimum spacing between 4.00 D and 4.75 D; 100 % are below 5 D.
- The relaxed 4D constraint enlarges the feasible set, so it should make beating the 7D grid easier. That the optimizers still do not beat it (Issue 1) is informative, not a caveat in their favour.
- At 7D inside the corner-turbine parallelogram, the problem is nearly degenerate: the installed grid is close to the only packing.
- At 4D, added turbulence and fatigue loads would be prohibitive offshore. The Jensen model does not see them.

**Fix**
- Reframe the caveat: the optimized layouts use about 4–4.75 D, a spacing that would not be accepted offshore for load reasons, and still do not exceed the 7D grid under a finer rose.
- Consider a 5D/6D run or an enlarged boundary. The paper has the infrastructure (Table S-capacity).

### Issue 9 (MODERATE): the Gaussian re-evaluation is a layout-robustness test, not an optimizer-robustness test

**Location**
- optA/09_robust.tex:7.
- MPCE_PSO_VNS_supplement.tex:398.
- optA/06_results.tex:39.

**Finding**
- Jensen-optimized layouts are re-scored with Bastankhah (K_G = 0.04, hub centre, same 24 × 15° bins, not calibrated to K = 0.075).
- The low τ̄ = 0.66 and the rank flip mix three effects:
  - the different wake shape;
  - the uncalibrated expansion;
  - bin-gap exploitation that a smooth profile partly penalizes, since the 15° bins are unchanged.
- Nothing in it says which optimizer is better under a Gaussian model.

**Fix**
- Re-optimize at least the six largest cases with PSO-VNS, PSO, VNS and MS-SLSQP under the Gaussian model with 1–5° bins.
- Calibrate K_G to the Jensen K, or state explicitly that it is not calibrated.
- If that is not done, qualify "PSO ranks narrowly ahead" as a statement about re-scored Jensen layouts.

### Issue 10 (MODERATE): the HR run-level "significantly better than every other method" merges feasibility with AEP

**Location.** optA/08_beyond.tex:35 ($p_{\rm Holm}\le\NHRMaxP$), and analysis/mpce_tab_hr16.tex (column $p_{\rm Holm}$).

**Finding**
- Against PSO the test ranks PSO's 11 infeasible runs last, so it is largely a feasibility test.
- On the 19 jointly feasible seeds, the AEP advantage is not robust to the model (Issue 1).
- Against DE (0/30 feasible) the test is purely about feasibility.

**Fix.** Report feasibility (count or Fisher test) and AEP on the jointly feasible seeds separately, as the benchmark does with W/T/L plus case means.

### Issue 11 (MODERATE): simplifications of the HR model should be stated where the model is described

**Location.** analysis/hornsrev_model.py:42, :93–96; optA/05_setup.tex:27; MPCE_PSO_VNS_supplement.tex:195.

**Finding.** The model has four simplifications that affect the wake loss and the HR reference comparison:
- C_T at the free-stream speed, not the local speed;
- hub-centre top-hat;
- no turbulence dependence;
- a single-block evaluation. The block is the north-west corner, indices 0–3 of columns 0–3. It is upwind for the prevailing W–SW sectors, but waked by the rest of the farm for E/SE winds, about 20 % frequency.

The supplement mentions C_T and the hub-centre test. The main text only says "the Jensen wake with K = 0.04".

**Fix**
- In the main text, add: "C_T at free-stream speed; underestimates the wake loss by ~0.7–1 pp relative to PyWake's NOJ at the same bins".
- State that the block is the NW corner and what its exclusion of external wakes implies.
- Better, evaluate with PyWake directly.

### Issue 12 (MINOR): the non-physical upstream deficit should be quantified

**Location.** optA/03_model.tex:12.

**Finding**
- In 315 of 4,080 PSO-VNS and PSO final layouts, a turbine receives the upstream (between the vertex and the rotor) deficit.
- The mean contribution is 0.007 pp (14 % of the margin), with a maximum of 1.48 pp in one layout.
- For a pair aligned at 4D, both turbines get the same deficit. The benchmark therefore penalizes close aligned pairs twice.

**Fix.** Report these numbers in one sentence. It shows that "kept for comparability" is harmless on average.

### Issue 13 (MINOR): benchmark power and thrust assumptions and the high-wind regime

**Location.** optA/03_model.tex:16 and :20.

**Finding**
- In DS I (ψ = 13 m/s, k = 2), 31 % of the time is above rated and 2.5 % above 25 m/s.
- Constant C_T = 0.8 overstates the deficits there, and the missing cut-out adds about 2.5 % of time at rated power.
- The cubic and cut-out re-evaluation covers the power curve but not the C_T.

**Fix.** Give these two percentages. Note that DS I results weight a regime in which real turbines have C_T of about 0.1–0.3.

### Issue 14 (MINOR): the highlights overclaim

**Location.** Highlights_PSO_VNS.txt:4 ("Constriction PSO followed by basic VNS is the most consistent layout optimizer").

**Finding.** The claim contradicts the paper itself:
- VNS is first from feasible starts;
- PSO is first under the Gaussian re-evaluation;
- DE is best on HR at 120k;
- the best published IEA37 layouts are better.

**Fix.** "…is a consistent reference method on the Kusiak–Song benchmark with random starts".

### Issue 15 (MINOR): leftover markers and wording

**Location and fix**
- optA/05_setup.tex:28 (`% CHECK (lead, D5): … confirm this description`): resolve it.
- mpce_tab_iea37.tex:25 `\TBD{verify…}`: resolve it.
- Supplement :349, "with the same settings": see Issue 2.
- optA/08_beyond.tex:50: add "(at a feasibility tolerance of 1 mm for published layouts; our layouts 10⁻⁶ m)".

---

## 2. Positive points a wind-energy reviewer would acknowledge

- The IEA37 AEP implementation is exact against the official calculator and the published files.
- The 16-direction, single-speed, Gaussian nature of the case is stated correctly.
- The benchmark energy conversion (objective/15 × 8.76) and the 62 % capacity-factor warning are correct and candid.
- The Horns Rev data (V80 curves, 12-sector Weibull, layout) match PyWake 2.6.20 exactly. The wake-free AEP is identical: 744.04 GWh/yr for 80 turbines and 148.81 for 16.
- The frequency-normalization fix was needed and correctly implemented.
- The authors already avoid claiming that the installed HR layout can be improved, and they list turbulence, loads, cables and LCOE as out of scope.
- The Data Set II N ≥ 10 advantage of PSO-VNS is robust to the finer rose.

## 3. Overall verdict

**Major revision.**

The benchmarking protocol is careful. As a statement about metaheuristics on the discretized Kusiak–Song benchmark, the paper is largely sound.

The wind-energy-facing claims, however, rest on model discretization artifacts that are larger than the effects claimed:
- the Horns Rev "highest mean AEP at 6,030" and "runs above installed";
- the practical meaning of the ±0.05 pp equivalence;
- the IEA37 gap sizes.

Specifically:
- The optimizers, the VNS-based ones especially, gain about 1 pp of wake loss by exploiting 5° and 15° bins.
- The HR model validation rests on a non-reproduced constant.
- The IEA37 "best feasible" reference hinges on a 1 mm tolerance.

None of this invalidates the paper's methodological message, but the HR and practitioner statements must be corrected or re-scoped. The fixes are cheap: re-evaluation only, plus PyWake, which is pip-installable.

## 4. Top 5 fixes

1. **Horns Rev.** Re-evaluate all HR layouts and the installed block with 1° bins and with PyWake (NOJ, Bastankhah, TurbOPark).
   - Delete or qualify "N runs exceed the installed layout", which is 0 under every finer model.
   - Restate the 6,030 result as feasibility only.
   - Keep "highest mean AEP at 30,030", which is robust.
2. **PyWake validation.** Replace the stale 662.5 GWh/yr constant with an actual PyWake 2.6.20 run at identical bins: NOJ gives 661.39 GWh/yr (11.11 %) against our 666.75 (10.39 %).
   - Report the wake-loss difference.
   - Attribute it correctly: C_T at free-stream speed and bin phase, not rotor averaging.
3. **Direction resolution.** Add a direction-resolution sensitivity (15 × 1° sub-bins) to Section 9 and Table tab:robust-final.
   - Report the equivalence test under it: CI [−0.063, −0.025], equivalence lost.
   - Report it on non-trivial cases: 44 cases, CI [−0.060, 0.008].
   - Justify the margin in MWh per turbine per year against the paper's own model-form shifts, and drop "only 0.031 pp".
4. **IEA37 gaps.** Recompute with boundary-projected published layouts (par12: 421,451 and 882,383 MWh), in both AEP % and wake-loss pp.
   - Report both tolerance conventions.
   - Resolve the Baker et al. (2019) TBD and the participant attribution.
5. **Practitioner scope.** Either add a smooth-model test (Gaussian re-optimization with fine bins, plus one AD-gradient optimizer on IEA37), or re-scope the practitioner-facing claims and the highlight to "the Kusiak–Song benchmark with random starts".

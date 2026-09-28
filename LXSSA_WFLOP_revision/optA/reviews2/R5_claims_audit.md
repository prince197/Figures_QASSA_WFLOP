# R5 (round 2): claims-vs-data and reproducibility audit

Scope: main text `optA/01_front.tex` to `optA/10_limits_concl.tex`; the prose of `MPCE_PSO_VNS_supplement.tex` and `optA/supp_theory.tex`; the generated table captions `analysis/mpce_tab_*.tex`; `analysis/README_reproduce.md` and `analysis/build_mpce_paper.sh`.

Method: I expanded every `\N...` macro to its value from the three `mpce_numbers*.tex` files, read every sentence against the CHECK comment next to it, read the condition code of each cited check, and looked up the underlying values in `mpce_summary.json`, `mpce_summary_extra.json` and `mpce_summary_diag.json`. I did not change any repository file except this report.

## 0. Check runs (2026-09-28)

| Script | Result | Notes |
|---|---|---|
| `mpce_check_final.py` | 59 PASS / 0 FAIL / 0 PENDING | 10 are defined but never referenced: C11, C30, C31, C33, C34, C35, C36, C43, C46, C48. Several of these cover sentences that are in the text (see U-items below). |
| `mpce_check_extra.py` | 37 PASS / 0 FAIL | X01, X12, X13, X15 and X28 are defined but not referenced. |
| `mpce_check_diag.py` | 19 PASS / 0 FAIL / 0 PENDING | D18 and D19 are not referenced. |
| `make_theory_figures.py --check` | 57 PASS / 0 FAIL | Run on a scratch copy so that `figures_mpce/theory_stability.pdf` was not overwritten. T13 is not referenced. The T-checks are substring tests: they confirm that a number appears somewhere in `supp_theory.tex` or `04_methods.tex`, not that it appears in the right place. |
| LaTeX logs (existing) | 0 undefined references | Main text 13 pp, supplement 38 pp. The PDFs are newer than all sources. |

Tags that the scanner cannot parse: `CHECK-FINAL [NEW, D2]` appears in 01:24 and 02:35, and `CHECK-FINAL [NEW, D3]` in 02:21. They contain no `Cnn`, so the scanner silently ignores them. They should be `[C32]` and `[C45]`.

Legend: **FALSE** = the text contradicts the data or another part of the paper. **OVERSTATED** = the wording is stronger than the data or the check support. **UNCHECKED** = the claim has no working check, or a hand-typed number has no macro. **STALE** = old wording, comment or number.

## 1. FALSE (number or statement contradicts the data or the paper itself)

| # | Location | Finding | Fix |
|---|---|---|---|
| F1 | 07:8 and 07:24 (`\NAbl...DL`) vs Table `tab:ablation` (`\NCm...DL`) | The same symbol ΔL has different values in the text and in the table beside it. SSA-VNS vs SSA: text −0.328, table −0.294. LX-SSA-VNS vs LX-SSA: text −0.431, table −0.394. LX-SSA vs SSA: text +0.166, table +0.153. The text macros use `ablation.contrasts[*].dloss_pp`, which averages over a different set of cases than `case_mean.mean_dloss_pp`. 07:5 defines ΔL as the mean difference of the case means, which is the table's value. | Use `\NCmSSAVNSvsSSADL`, `\NCmLXSSAVNSvsLXSSADL` and `\NCmLXSSAvsSSADL` in 07. Also rename or drop the `\NAbl*DL` macros so they cannot be reused. |
| F2 | 04:28 ("the 282 infeasible final layouts") vs 06:9 and 02:16 ("feasible in 86.0% of the runs") | 86.0% feasible means 2040 − 1754 = **286** infeasible runs under the paper's 10⁻⁶ m rule (03:38). `stored_violation_summary()` in `mpce_diagnostics.py` uses a 1 mm tolerance, so it reports 282. The paper therefore gives two different infeasible counts. | Either classify with the stored `Feasible` flag (286), or write "282 of the 286 infeasible final layouts violate a constraint by more than 1 mm". Make D06 assert `infeasible == runs − feasible_flag`. |
| F3 | Supplement §S-inference:317 ("with any threshold from 1 to 30 feasible runs") and the X02/X05 descriptions ("every threshold 1–30") | Only six thresholds were run: `thresholds = [1, 10, 15, 20, 25, 30]`. Intermediate values were never tested. | Write "with thresholds of 1, 10, 15, 20, 25 and 30" (`\NXThrList`), as 05:19 and 06:23 already do, and fix the check texts. |
| F4 | Internal contradiction: `supp_theory.tex` (equivalence table: "Only the left column supports statements such as 'adds nothing'") vs 01 abstract, 02:9/22, 07:14, 10:7/38 and supplement 291/295 ("SSA adds nothing" / "gains nothing" beyond disc sampling) | For SSA-VNS vs RSD-VNS, the 95% CI [−0.033, 0.068] contains 0 and the 90% CI [−0.022, 0.062] is **not** inside ±0.05. By the paper's own four-cell table, this is "inconclusive". The claim is defensible only one-sided: the 90% lower bound −0.022 > −0.05 rules out a gain larger than the margin (non-superiority). No check tests this: C54 does not check `ci90[0] > −m`. | Say "no gain of practical size (a gain larger than ±m is excluded, 90% CI lower end −0.022 pp)". Add a one-sided non-superiority row or sentence to the supp_theory table. Add `r["ci90_mean_dloss_pp"][0] > -m` to C54. |
| F5 | 08:41 comment "(C30 no longer asserted here ...)" vs 08:35 text | The text asserts exactly C30: "At 30,030 evaluations, PSO-VNS still has the highest mean AEP ... and 3 of its 30 runs exceed the installed layout". C30 passes but is unreferenced. | Delete the comment and attach `% CHECK-FINAL [C30]` (and C39) to the sentence. |

## 2. OVERSTATED

| # | Location | Finding | Fix |
|---|---|---|---|
| O1 | Abstract ("the Laplace step is harmful"); 02:9 ("the Laplace step ... makes the hybrid worse"); 02:22 ("The Laplace step is harmful (p = 5.4e-6)"); 07:24 first sentence; supplement 295 | 07:24 itself says these contrasts "do not separate the Laplace step from this cost": LX-SSA spends 2N_p evaluations per iteration and so runs half as many iterations. The data show that **LX-SSA as published** is worse, not that the Laplace step is. | Write "LX-SSA as published (Laplace step plus doubled evaluation cost) is worse" everywhere. Keep "harmful" only with that qualifier. |
| O2 | Abstract ("the swarm does not contract"); 04:28 ("the old swarm does not contract") | The diagnostics show that it does contract: the median spread goes from 773 m to 319 m (0.43 r), mostly in the first iterations, and then plateaus. D03 compares only with the constriction swarm (8.9×). No check compares the old swarm with its initial spread. | "the swarm stops contracting at about 0.4 r (8.9 times the constriction spread)". Optionally add a D-check on `curves_median.PSO.spread`. |
| O3 | Abstract ("so ... most evaluations are infeasible") | This does not distinguish the two settings: most evaluations are infeasible under **both**, 98.3% (old) vs 68.4% (constriction) over the whole run. The contrast is only clear late in the run (97.4% vs 51.6%). | "almost all evaluations late in the run are infeasible (97% vs 52% with constriction)", using the `\NDFeasPart...` macros. |
| O4 | Conclusion 10:38 ("puts it ahead of every salp-swarm method") | There is no budget qualifier, yet C45 shows that SSA-VNS ranks ahead of PSO at 120,030 evaluations (06:13, 08:8 and 02:16 say this). | Add "at 6,030 evaluations". |
| O5 | 02:9 thesis ("salp-swarm phases add nothing beyond random sampling"); Conclusion 10:38 (same, unqualified) | Against **square** sampling (RS-VNS), the SSA phase does give a significant gain (C15: p = 4.5e-4, 2/66/0). Only the in-farm (disc) control shows no gain. | Write "beyond random sampling in the farm", as the abstract does. |
| O6 | Abstract ("with the constriction setting it ranks first") | This refers to the previous study's pool (Table `tab:baseline`). A reader of the abstract will take it to mean first of the eight methods, which would contradict "PSO-VNS ... best average rank". | "ranks first in the previous study's method pool". |
| O7 | Abstract and 02:14 ("These conclusions hold across qualification thresholds, farm clusters and one Holm family") | "LX-SSA is worse" is not significant with the cluster-robust t test (X34, p = 0.087). SSA-VNS vs RS-VNS is not cluster-robust either (X17, p = 0.061). The "larger Data Set II layouts" subgroup has no cluster-level analysis (only a within-cluster permutation, X22). | "The main conclusions hold across thresholds and one Holm family; at the farm-cluster level all but the two salp-vs-disc contrasts are unanimous". Alternatively, name the exceptions. |
| O8 | 06:39 ("it adds feasibility in every run"); 10:16 ("PSO-VNS a more reliable one") | PSO is also 100% feasible on all 68 cases and on the six largest cases at every budget. The reliability advantage rests on a single site (Horns Rev 1: 30/30 vs 19/30), and the spread does not differ (C52). | "is feasible in every run, whereas PSO is not on the Horns Rev 1 block (19/30)" and "more reliable in reaching feasibility on the Horns Rev 1 block". |
| O9 | 07:8 ("switching to VNS gives a small gain (−0.018, p = 0.30) that is practically equivalent to zero") | C13 checks only run-level W > L. The case-mean difference is not significant and is equivalent to zero, so "gives a gain" reads stronger than the data allow. | "a small, nonsignificant difference ... practically equivalent to zero". |
| O10 | 02:22 ("the VNS phase is one compass-search descent"); 10:38 ("here one compass-search descent") | D08 gives "essentially": 97.5% of the gain comes from the first descent, the median run has 1 shaking step, and 6/120 runs have an accepted cycle. | Add "essentially", as 04:89, 07:8 and the supplement do. |
| O11 | Supplement 297 and 427 ("the neighbourhood structure of VNS matters only for small N") | The main text (07:8) says "can matter only". Nothing shows that it matters for small N. | "can matter only for small N". |
| O12 | 08:19 ("PSO may thus help PSO-VNS mainly by reaching feasibility") vs 07:28 ("Whether PSO helps mainly by faster feasibility or by better exploration remains open") | The two sections lean in different directions. | Harmonize by keeping "remains open" in both, or cite D18 in both. |
| O13 | Notation of the Data Set II Holm p: 06:30 "p = 0.018"; 10:12 "p = 0.018"; 02:32 and the supplement "p_Holm ≤ 0.018" | `\NXHolmPSOVNSvsPSOdsIILargeP` is an upper bound rounded up. The same test is also said to be Holm-adjusted over "all 28 case-mean tests of the paper" (05, 06), "over all tests" (10) and "quoted in the main text" (supplement, X18). | Use "≤" everywhere and write "over all 28 case-mean tests of the main text". |

## 3. UNCHECKED (claims without a working check, or hand-typed data)

| # | Location | Claim | Fix |
|---|---|---|---|
| U1 | 02:22 ("significantly improves both salp-swarm methods but PSO only slightly"); 07:8 ("The VNS phase significantly improves both salp-swarm methods") | The attached C13 tests only PSO-VNS vs PSO W > L. The right check, **C33**, is unreferenced, and 07:12 has only a "verified" comment. | Attach `[C33]`. |
| U2 | 01, 02:28, 06:23/30/32, 07:10 ("posterior probability above 0.99", "P(ROPE) > 0.99") | This is hand-typed. `\NBayPSOVNSvsPSORope` (= "> 0.99") exists, and the supplement uses it. | Use the macro. |
| U3 | 05:31 | "largest difference 8.7×10⁻¹¹ over 720 archived final objectives", "relative difference < 10⁻¹¹" (IEA37 calculator), "2 of 24 Mann–Whitney tests reject; 1.2 expected" are hand-typed and have no macro or check. They come from `validate_evaluator.py`, `iea37_model.py` (`__main__`) and `calibrate_authors_code.py`, none of which the build runs. | Write these into `mpce_summary.json` (or a small `mpce_numbers_valid.tex`) and add a check. At minimum, list `python3 iea37_model.py` in README §3. |
| U4 | 08:35 (IEA37: "the best of our layouts comes from MS-SLSQP") | **C48** exists but is unreferenced. | Attach `[C48]`. |
| U5 | 08:57 (cost: 0.23–0.45 ms, MS-SLSQP 2.3 × PSO) | **C46** exists but is unreferenced. | Attach `[C46]`. |
| U6 | 08:8 ("ranks first in 4, 5 and 3 of the six cases") | **C34** exists but is unreferenced. | Attach `[C34]`. |
| U7 | 05:10 (seed pairing; "the best initial objectives agree") | **C43** exists but is unreferenced. There is only a "verified" comment (792 pairs; C43 compares 1,002). | Attach `[C43]`. |
| U8 | 07:5 ("their average ranks differ (Friedman p = 6.0e-76)") | **C36** exists but is unreferenced. | Attach `[C36]`. |
| U9 | Supplement 291 ("an LX-SSA phase none" vs RS-VNS) | **C35** exists but is unreferenced. | Attach `[C35]`. |
| U10 | 01:27 and 02:38 (`CHECK (manual)` for the IEA37 gaps) | C40 automates this. | Replace the manual tag with `[C40]`. |
| U11 | 08:8 ("while VNS (6.50 to 3.17) and DE move up") | The implemented C37 does not test the rank moves that the 08:11 comment proposes. | Extend C37, or cite the generated ranks only. |
| U12 | 08:19 ("most trials also place a turbine outside the site") | D14 does not test the 96.4% / 99.4% outside-site shares. They appear only in the comment. | Add them to D14. |
| U13 | 08:14 ("PSO-VNS follows" in mean wake loss and rank, feasible starts) | C28 checks the rank only. The loss is 3.938 vs 3.952 (LX-SSA-VNS) and 3.955 (SSA-VNS), i.e. second by a margin of 0.014 pp, and unchecked. | Extend C28/C38 to the loss, or say "follows in rank". |
| U14 | Supplement 305 ("the Bayesian analysis finds it more likely practically negligible, P_rope 0.70 vs P_A 0.30") | C51 does not test P_rope > P_A. | Add this to C51. |
| U15 | 05:28 (`CHECK (lead, D5)` Horns Rev description) and 05:33 (`CHECK (manual)` CPU) | These are unresolved manual confirmations. | Resolve them and delete the tags. |
| U16 | 09:15 ("the largest 500-m cases may be infeasible" at 5D/6D) | This comes from `packing_capacity.py`, which is not in the build, and has no check. | Add a macro and check, or list it explicitly as manual. |
| U17 | 06:39 and 09:7 ("PSO ranks narrowly ahead" under the Gaussian model) | C22 does not bound the gap (1.91 vs 1.97). | Acceptable. Optionally add `rank(PSOBV) − rank(PSOC) < 0.2` to C22. |
| U18 | Check-code weaknesses | C27 uses `any()`, but the text says "behind salp-swarm hybrids" / "both". The data support it (1.99 and 2.65 < 5.04), but the check does not. The D04 description says "iterations 101–200" but the check tests `clip_pct_all`. D09 hard-codes 36 and 20. | Use `all()` in C27, fix the D04 text, and derive the D09 constants. |

## 4. STALE

| # | Location | Finding | Fix |
|---|---|---|---|
| S1 | 01:23 and 02:33 comments ("C49 tests only > 0.5 ...") | C49 now tests P_rope > 0.99. | Update the comments. |
| S2 | 01:24, 02:21, 02:35 (`[NEW, D2]`, `[NEW, D3]`) | The scanner cannot parse these tags. | Replace with `[C32]` and `[C45]`. |
| S3 | 08:11, 08:29, 08:40, 08:52 ("C37/C38/C39/C40 proposed in R5-xx, to be defined in mpce_check_final.py") | All four are now defined. The 08:11 condition text differs from the implemented C37. | Replace with the implemented condition texts. |
| S4 | 08:42 ("verified ... (pre-rerun data): at 120,030 BVNS 6.12 ... PSOBV 6.45; PSOBV lowest at 6k R (7.05) and 30k (6.46)") | These are pre-hrfix values. Current values: PSO-VNS 7.07 / 6.52 / 6.40, and DE has the highest mean AEP at 120,030. | Delete the comment. |
| S5 | 08:51 ("verified ... 36T/30k ... gap ... −4.20 % ('larger for 36 turbines')") | −4.20% is the best PSO-VNS gap. The best-of-all gap (the C40 quantity) is −3.51%. The text is correct; the comment conflates the two. | Correct the comment. |
| S6 | 07:26 ("verified ... LXSSA-SSA case means 9 lower / 53 higher, p = 8.3e-8") | The data give 9/51 and p = 1.9e-7 (the text macro is correct). | Update or delete the comment. |
| S7 | Supplement Table `tab:S-repro` ("checks X01–X28") | The checks now run X01–X37. | Change to X01–X37. |
| S8 | `supp_theory.tex` header ("checks T01-T12") | T13 exists. | Change to T01–T13. |
| S9 | 11_back ("\TBD{Zenodo DOI}"); 01 (funding, corresponding-author TBDs); Acknowledgment TBD | These are author placeholders. | For the authors. |
| S10 | Searches for "eleven contrasts", LOCO "no verdict changes" and "ten times" | None left in the prose. The LOCO wording in 06:23 ("changes none of them" = comparisons **with PSO-VNS**) is correct per X14. | None needed. |

## 5. Consistency across abstract, introduction, results, discussion, conclusion and supplement

- The core message is consistent: constriction vs old PSO; salp-swarm phases give nothing beyond disc sampling; PSO phase > both controls; PSO-VNS ≈ PSO; reference method, not superior. The inconsistent points are:
  - "Laplace step" (O1).
  - "random sampling" with and without "in the farm" (O5).
  - "is one descent" vs "essentially" (O10).
  - "neighbourhoods matter" vs "can matter" (O11).
  - The feasibility-vs-exploration interpretation (O12).
  - The Holm notation and family wording (O13).
  - The ΔL values in 07 vs Table ablation (F1).
  - The 282 vs 286 infeasible runs (F2).
- Numbers repeated across sections (e.g. +0.024, +0.087, −0.306, 730/494, 30/30 vs 19/30, 8.9×, 2.6% vs 48.4%) all come from the same macros and agree.
- The W/T/L of PSO-VNS vs SSA-VNS differs between 06 (49 wins) and 07 (44 wins) because of the different Holm families. The table captions explain this; the text does not. Add a half-sentence in 07.

## 6. Reproducibility: README_reproduce.md + build_mpce_paper.sh

The build regenerates `mpce_numbers.tex`, `mpce_numbers_extra.tex`, all main tables and the figures from the per-run CSVs, runs the three check scripts and compiles both PDFs. The following gaps remain:

1. **Theory checks are not in the build, and `make_theory_figures.py` has no check-only mode.** `--check` always rewrites `figures_mpce/theory_stability.pdf`, and the T-checks are substring tests. Add a `--check-only` flag (no figure output) and call it from the build; it needs about 3 minutes, or less with a cached Monte Carlo.
2. **`mpce_numbers_diag.tex` is never checked against `mpce_summary_diag.json`.** No analogue of X28 exists, and the build does not rerun `mpce_diagnostics.py` (about 18 minutes). Add a D-check comparing the generation stamps.
3. **The IEA37 published-results file has no generating command.** README §1 says `iea37_published_results.csv` is "built from the IEA37 repository", but no script is named. `mpce_results.py` treats the file as optional, so without it C40 and all published-layout comparisons go PENDING.
4. **The IEA37 calculator check is not documented.** The "< 10⁻¹¹ relative" figure comes from `python3 iea37_model.py`, which README §3 does not list.
5. **PyYAML is missing from `requirements.txt`.** `iea37_model._yaml` needs it for all IEA37 loading.
6. **The PyWake reference is a constant.** `PYWAKE_HR80_AEP = 662.5`: PyWake is not installed and not pinned, so `\NHRPyWakeDiff` (0.64%) cannot be reproduced. Pin `py_wake==2.6.20` and give a script or command that computes 662.5.
7. **Evaluator validation is not in the build** (`validate_evaluator.py`, `calibrate_authors_code.py`), and its numbers are hand-typed (U3).
8. **`mpce_feas_s0of1.csv` is concatenated by hand from the 8 `feasx` shards.** No command is given. Add a one-line pandas command or script.
9. **The robustness re-evaluation reads a cache by default** (`mpce_reevaluation_cache.csv`). The build never recomputes it. Document `rm mpce_reevaluation_cache.csv` as the full-rebuild step, or add a `--no-cache` flag.
10. **Supplement Table S-repro says "grid scripts of the previous version (`full_grid_experiments.py` and successors)".** "Successors" is vague. README gives the exact ids (`grid`, `vgrid`, `bgrid`); use them.
11. **The hardware line in 05:31 is still a manual CHECK** (U15), and the timings come from different machines (supplement §S-cost).
12. **The scanner in `mpce_check_final.py` accepts malformed tags silently** (`CHECK-FINAL [NEW, D2]`). Make it report any CHECK-FINAL comment that contains no `Cnn` id.

## 7. Top 10 fixes (by priority)

1. **F1**: in 07:8/24, replace `\NAblSSAVNSvsSSADL`, `\NAblLXSSAVNSvsLXSSADL` and `\NAblLXSSAvsSSADL` by the `\NCm...DL` macros, so the text matches Table ablation (−0.294 / −0.394 / +0.153).
2. **O1**: qualify "the Laplace step is harmful" in the abstract, 02:9, 02:22, 07:24 and the supplement as "LX-SSA as published (Laplace step and doubled evaluation cost)". 07 itself says the two are confounded.
3. **F4**: make "SSA adds nothing beyond disc sampling" consistent with the supp_theory equivalence rule. State it as non-superiority (a gain > m is excluded, 90% lower end −0.022), add that case to the table, and add `ci90[0] > −m` to C54.
4. **F3**: in the supplement, replace "any threshold from 1 to 30" by the six thresholds tested, and fix the wording of X02–X09 and X31.
5. **F2**: reconcile 282 (1 mm rule, diagnostics) with 286 (10⁻⁶ m rule, 86.0% feasible). Classify by the stored `Feasible` flag, or state the tolerance.
6. **O2/O3**: in the abstract and 04:28, replace "does not contract" and "most evaluations are infeasible" with the measured contrast: it plateaus at about 0.4 r, 8.9× wider; 97% vs 52% of evaluations infeasible late in the run.
7. **O4/O5/O6**: in the conclusion add "at 6,030 evaluations"; in the thesis and conclusion add "in the farm" to "random sampling"; in the abstract say "ranks first in the previous study's pool".
8. **U1–U10 and S2**: attach the existing but unreferenced checks (C30, C33, C34, C35, C36, C43, C46, C48, C40), replace `[NEW, D2]`/`[NEW, D3]` by `[C32]`/`[C45]`, and use `\NBayPSOVNSvsPSORope` instead of the hand-typed "0.99".
9. **S3–S8**: clean up stale comments: "proposed/to be defined" C37–C40, pre-hrfix Horns Rev values, the "C49 tests only > 0.5" note, the LX-SSA vs SSA "verified" p, X01–X28 → X01–X37, T01–T12 → T01–T13.
10. **Reproducibility items 1–6**:
    - Add a `--check-only` theory check and a diagnostics sync check to the build.
    - Add `pyyaml` (and a pinned `py_wake`, or a note) to `requirements.txt`.
    - Name the command that builds `iea37_published_results.csv`, and list `iea37_model.py` among the validation scripts.
    - Move the evaluator-validation numbers of 05:31 into generated macros with a check.

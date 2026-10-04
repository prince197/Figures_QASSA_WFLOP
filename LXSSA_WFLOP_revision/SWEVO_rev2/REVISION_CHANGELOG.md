# Disposition of the twelve ZIP-audit areas

**Overall status: substantially improved; further author work required before submission.** Five areas are fully addressed at the scope stated below; seven are partially addressed. “Fully addressed” concerns the specified correction, not validation of every original result or journal acceptance. No original optimizer was rerun because required modules are missing.

| Area | Status | Incorporated correction | Remaining work |
|---|---|---|---|
| A1. Conflicting manuscript versions | Fully addressed | Recovered the fuller later scientific source into every modular section; generated the single-file manuscript from that canonical source; rebuilt both PDFs and refreshed supplementary references. | Future edits must follow the canonical-source workflow. |
| A2. Incomplete reproduction package | Partially addressed | Added a working standalone archive audit, focused checks, build script, environment pins, coverage inventory and honest reproduction instructions. Replaced unsupported “all checks pass here” wording with historical attribution. | Supply 36,190 original optimization records, diagnostic records, core implementation, original analyses and environment. |
| A3. Coordinate precision / strict feasibility | Partially addressed | Quantified rounded-coordinate failures and Lillgrund AEP replay differences; retained original labels/objectives. Five future writers now preserve 17 significant digits, with binary round-trip checks. | Recover original unrounded final positions or rerun affected experiments. Strict archived feasibility is not independently repaired. |
| A4. PSO variance proposition | Fully addressed | Main and supplementary statements now say zero limit variance iff c1*c2*(p−g)=0; p=g is equivalent only when both coefficients are positive. Added proof explanation and checked the zero-coefficient counterexample and both PSO settings numerically. | No experimental rerun needed for this algebraic correction. |
| A5. Conditional ranking versus reliability | Partially addressed | Explicitly named qualification-then-conditional quality; added Lillgrund success counts/Wilson intervals and paired all-run superiority with bootstrap intervals, exact sign tests and a 16-comparison Holm family. Explained the 17/30 versus 30/30 reversal. | Apply an agreed all-run endpoint to all original benchmark/site cases after supplying the missing records; recover precise feasibility evidence. |
| A6. Inference estimands and calibration | Partially addressed | Removed new-farm generalization from group sensitivity, described case exchangeability as hypothetical, clarified existing case weights and cross-case seed dependence, and limited the intersection–union claim to calibrated constituent tests. Bootstrap equivalence is described as approximate interval inclusion. | Original inference code/data are needed for joint seed resampling, null calibration and any preregistered or population-level claims. |
| A7. Visibility of expanded baselines | Fully addressed | Added GA to the main Horns Rev table, GA and both analytic-gradient variants to main IEA37 best-run results, an IEA37 means/feasibility table, and the nine-method benchmark summary. Historical eight-method ranks remain explicitly separate. | Stronger DE/CMA-ES experiments remain an extension, not a result claimed here. |
| A8. Analytic-gradient interpretation | Fully addressed | Explained analytic objective AND constraint derivatives, different hybrid restart choices, piecewise smoothness, uncharged constraint/solver costs and limitations of batch timings. Added the c_g=19 accounting bounds as saved-prefix arithmetic, not a new run. | An objective-gradient-only ablation would be new research; the revised text does not claim it was done. |
| A9. Evaluator transfer versus re-optimization | Fully addressed | Identified finer-rose/Gaussian/power-curve studies as re-evaluation of saved layouts; identified spacing and gradient studies as direct re-optimization; retained the reversal at 6D. | Optimization under alternative evaluators is an optional new experiment and is not implied by the current results. |
| A10. Constraint handling / baseline fairness | Partially addressed | Added implementation-specific scope, penalty-unit and repair differences, and limits of the single DE configuration. | Run normalized-penalty/repair controls and stronger or tuned baselines before claiming algorithm-family superiority. |
| A11. Declarations and author facts | Partially addressed | Synchronized manuscript and separate declaration text with the supplied roles/funding/interest statements; disclosed ChatGPT assistance; replaced unsupported open-repository promises with the archive's actual status; updated cover letter/highlights. | Authors must confirm all facts, affiliations, roles, AI wording, manuscript approval and submission eligibility. No public DOI has been invented. |
| A12. Presentation and submission package | Partially addressed | Rebuilt PDFs, refreshed references, added legible compact main tables, concise highlights, a shorter cover letter, a starting guide and this implementation record. Archived original loose notes separately. | The main manuscript remains long and dense; authors should choose any further shortening. Journal acceptance and editorial formatting are not guaranteed. |

## New Lillgrund analysis

The endpoint uses ORIGINAL feasibility labels: feasible beats infeasible; feasible pairs compare AEP; both infeasible tie. It does not reproduce the older spacing-ordered signed-rank endpoint, and its p-values must not be substituted into that analysis. It also does not repair label uncertainty from rounded saved layouts.

- At 6,030 evaluations, PSO-VNS versus PSO: 17 wins / 13 ties / 0 losses, superiority 0.7833, 95% paired-bootstrap CI [0.7000, 0.8667], Holm sign-test p=0.00016785.
- At 6,030, versus MS-SLSQP: 10 / 0 / 20, superiority 0.3333, CI [0.1667, 0.5000], Holm p=0.39495. The sign test does not establish a corrected significant difference; the previously reported signed-rank p-value is a different statistic.
- At 30,030, versus PSO and versus MS-SLSQP: 28 / 0 / 2, superiority 0.9333, CI [0.8333, 1.0000], Holm p=0.000012152.

Bootstrap: 20,000 resamples of 30 matched seed pairs, NumPy default_rng seed 20261004+budget; percentile interval. Holm covers eight opponents at two budgets. Wilson intervals are marginal 95% success-rate intervals. These uncertainties describe seeds at the Lillgrund block; they are not site-population intervals.

## Precision findings preserved as open issues

At 10^-6 m tolerance, rounded positions contradict feasibility labels in 236 exact-gradient, 78 Laplace, 229 spacing, 17 Lillgrund 6,030 and 73 Lillgrund 30,030 records. The Horns Rev boundary cannot be rechecked without its missing model. The maximum rounded-coordinate AEP replay discrepancy among labelled-feasible large-budget Lillgrund runs is 0.225072 GWh (PSO seed 18). No such record was relabelled or numerically overwritten.

## Verification completed for this revision

See `validation/archive_audit.json`, `focused_checks.json` and `revision_verification.json` for machine-readable evidence. The supported scope is stored-record and document verification. Historical 200+ claim checks, diagnostic reruns, external-model checks and full-study statistics were not rerun from unavailable code/data.

---

# Revision 3 (4 October 2026) — the "second revision"

Naming in the submitted documents: "first revision" = internal revision 2 (GA, Laplace, spacing, exact-gradient and
Lillgrund experiments); "this revision" = internal revision 3.

This section is appended to the record above; the earlier sections are kept unchanged as history.

**Basis.** This revision answers the remaining-work report on revision 2 (reviewer items R1–R10; work items A1–A8,
B1–B7, C1–C4, D, E1–E6). Triage: `REV3_TRIAGE.md`. Decisions D1–D14: `REV3_INTEGRATION.md`.

**Premise corrected.** Revision 2 assumed that the 36,190 original run records and fifteen core modules were
missing. They were not: they are in `analysis/`, and `rev2_analysis.py` reproduces the revision-2 summary exactly
(4,670 values). The statements built on that premise were deleted from the manuscript, the supplement, the Data
availability statement, the cover letter, `START_HERE.md`, `AUTHOR_CONFIRMATION.md` and
`latex_source/analysis/README_reproduce.md`. These statements had said the material was missing, "historical",
"author-reported" or "cannot be rerun". Areas A2, A3, A5, A6 and A10 of the table above are therefore closed by this
revision, as listed below.

**New runs.** 5,100 optimization runs, all with new seeds 31–60, prespecified in `analysis/rev3_*_manifest.md`,
with coordinates written to 17 significant digits:
- 2,700 runs of the constraint-handling study;
- 1,200 runs of direct 1-degree optimization;
- 1,200 seed-paired 15-degree control runs.

The deterministic full-precision reruns of the stored records come on top of these (C4, count
`\TBD{C4 rerun count}`). All other revision-3 results are re-analyses of stored records. The supplementary tables
`tab:S-r3-*` come from `analysis/rev3_*_tables.tex`. Section numbers below refer to the revision-3 build.

## Reviewer items

| Item | What was done | Where |
|---|---|---|
| R1 / A1 | Declarations complete. The generative-AI statement now follows the Elsevier template and names all tools and their uses (D14). Facts still need author confirmation. | Back matter; `declarations_content.tex`; `AUTHOR_CONFIRMATION.md` items 5–6 |
| R2 / A2 / D | Deposit-ready reproducibility package built and verified: README map, `MANIFEST.sha256`, requirements, `Dockerfile`, `CITATION.cff`. The complete pipeline was re-executed: 91 outputs identical apart from time stamps; all checks pass (62 C, 57 X, 22 D, 41 F, 15 S, 57 T); reruns of stored gradient-free runs are bit-identical except elapsed time (SLSQP-based runs are not bit-reproducible across CPU/BLAS builds). Data availability rewritten (D1): release 1.0.0, to be created from the release candidate 1.0.0-rc after a rebuild with the revision-3 records. Public DOI: pending author action (`\TBD{DOI}`). | `repro_package_build/` (+ `verification/`); Data availability; §6.4; supplement Reproducibility Map and archive section; `latex_source/analysis/README_reproduce.md` |
| R3 / B2 | Site pools recomputed with GA and exact gradients (D9). Horns Rev 1, pool of 11: PSO-VNS best at 6,030 evaluations, GA best at 30,030. IEA37, pool of 13: an exact-gradient method ranks first in every scenario and budget. | §9.3, §9.5; `tab:S-r3-hr-pool`, `tab:S-r3-iea-pool` |
| R4 / B1, B5 | All-run endpoint on the 68 cases (D4). PSO-VNS vs PSO: 950/348/742 of 2,040 seed pairs, score 0.551 [0.516, 0.587]. Against the other six main methods: 0.73–0.90. Imputation sensitivity: no verdict changes (D7). | §6.2, §7.3, §8; `tab:S-r3-allrun`, `tab:S-r3-imputation` |
| R5 / C1 | Constraint-handling study, 2,700 runs (D10). The conclusions depend partly on constraint handling. With projection repair the rank order is unchanged, but PSO-VNS no longer beats PSO significantly; with Deb's rules and box clipping, PSO-VNS's all-run advantage over GA and SSA-VNS is lost. | §3.4, §10.4, §11, §12, §13; `tab:S-r3-constraint`, `tab:S-r3-constraint-cases` |
| R6 | Proposition 1(c), corrected in revision 2. Verified by the theory checks (T01–T13). | — |
| R7 / B3, B4, B7 | TOST calibration by null simulation: both primary claims hold (D5). Synchronized-seed resampling: SSA-VNS vs RSD-VNS seed-level equivalence is not robust; its non-superiority is (D6; footnote b of `tab:equiv-main`, now Table 8). Equal-cluster target with a wild-cluster bootstrap: no verdict changes (D8). | §6.2; `tab:S-r3-tostaudit`, `tab:S-r3-syncseed`, `tab:S-r3-eqclus` |
| R8 / C2, C3 | Direct 1-degree optimization, 1,200 runs plus a 1,200-run 15-degree control arm (D11). The leader PSO-VNS is unchanged, and its advantage over PSO grows beyond the margin. GA (Horns Rev 1) and Lillgrund layouts re-evaluated with 1-degree bins and PyWake NOJ: GA's Horns Rev lead holds; at Lillgrund, PSO-VNS's lead over the runners-up (SSA-VNS, VNS, MS-SLSQP) is not significant at 1 degree, MS-SLSQP wins at 6,030, and under PyWake MS-SLSQP has the highest mean at both budgets (D9). | §9.3, §9.4, §10.5; `tab:S-r3-fine`, `tab:S-r3-hr-eval`, `tab:S-r3-lg-eval` |
| R9 / B6 | Per-evaluation cost given as medians. The MS-SLSQP time jump at N >= 8 is attributed to thread contention (re-timed single-threaded). Initialization time and total CPU time added (D2). | §9.7; `tab:S-r3-time`, `tab:S-r3-time-sites` |
| R10 / E1–E4 | Length and repetition pass: new results in one or two sentences with a pointer to the supplement, repeated caveats removed. Number fixes (D3). Main bibliography corrected and renumbered in order of first citation (D13). | All sections; `analysis/rev3_numbers_report.md`; `analysis/rev3_bib_order.py` |

## Work items A–E

- **A1, A2, A4.** Author confirmations remain; see `AUTHOR_CONFIRMATION.md` items 5–7. The AI declaration is
  identical in the manuscript back matter and `declarations_content.tex`, with `\TBD{authors: confirm the list of
  tools}`.
- **A3.** The primary equivalence questions are unchanged and their team agreement is still pending
  (`AUTHOR_CONFIRMATION.md` item 4).
- **A5.** Solanki and Deep (2023) now reads "published online May 18, 2023 (Online First)". Volume, issue and pages
  have not been assigned yet.
- **A6.** Reference audit: 85 entries OK, Thomas2023 corrected (no issue number), Solanki2023 corrected, IEA37repo
  comment corrected (`analysis/rev3_refs_report.md`, `rev3_refs_corrections.tex`). The main list is in
  first-citation order: 88 entries, all cited, none missing. The supplement list is ordered by I5, with uncited
  Benavoli2016 and Cleghorn2018 deleted.
- **A7.** SWEVO guide check: abstract ≤ 250 words, highlights, keywords, declarations and reference order. The
  length limit and the live guide could not be verified (`rev3_refs_report.md` Section 6).
- **A8.** `AUTHOR_CONFIRMATION.md` items 1 and 3 are resolved. Item 2 is resolved by audit and rerun, with the
  final count `\TBD{C4 rerun count}`.
- **B1–B7, C1–C3.** See the reviewer items above. All outcomes are reported, including unfavourable ones: the
  constraint-handling sensitivity, the Lillgrund 1-degree result, the 15.0% case-level size of the SSA-VNS
  non-superiority test before calibration, and the non-robust seed-level equivalence of SSA-VNS vs RSD-VNS.
- **C4.** All 56,540 stored records were audited (D12):
  - 50,826 labels confirmed from the rounded coordinates;
  - 0 contradicted;
  - 5,714 undecidable within the rounding bound.

  These undecidable records are rerun deterministically at full precision. Result (template of
  `analysis/rev3_precision_notes.md`): 1,160 + B_c bit-identical, 412 + F_c float noise, 1,572 + C_c confirmed,
  0 + X_c changed, 77 + U_c SLSQP-based undecided; the cloud counts are marked
  `\TBD{C4: enter the cloud counts $B_c$, $F_c$, $C_c$, $X_c$, $U_c$}`. The revision-2 statement that the Horns Rev 1
  boundary could not be audited is withdrawn; `hornsrev_model.py` is present.
- **D.** See R2. Remaining: the public deposit, the DOI and the licence choice (author action).
- **E1–E4.** See R10. Documents rebuilt with `latex_source/analysis/build_swevo.sh`.
- **E5–E6.** Highlights and cover letter are revised against the final text (`latex_source/SWEVO_highlights.txt`,
  `SWEVO_cover_letter.tex`). The point-by-point response to the reviewers should map R1–R10 to the rows of the table
  above and should state that only the author-action items are disclosed rather than resolved.

## Status of the revision-2 audit areas after revision 3

- **A2, A3, A5, A6, A10:** closed by revision 3, at the scope stated above.
- **A11:** closed apart from author confirmation and the DOI.
- **A12:** open only for the authors' final read of the rebuilt PDFs and any further shortening.
- **Still open, author action only:** confirmations, the public DOI, the final C4 rerun count.

## Fixes after the adversarial review (`analysis/rev3_review_report.md`, 4 October 2026)

- Reproducibility claims qualified: bit-identical reruns only for the gradient-free methods; SLSQP-based runs are not
  bit-reproducible across CPU/BLAS builds (§6.4, §9.6, S9, S11.11, S12.1, cover letter, response letter). The C4
  placeholders follow the template of `analysis/rev3_precision_notes.md` with one `\TBD` per place (§9.6, S11.11,
  response C4 row, `README_reproduce.md`); S12.2 now points to S11.11.
- Package version: release 1.0.0, to be created from 1.0.0-rc; the response (R2) keeps
  `\TBD{lead: confirm after package rebuild}` for the rebuilt package.
- §6.4: the calibration of the original code is recomputed from stored outputs, not re-executed.
- Lillgrund at 1 degree (§9.4, S10.5, S11.7, §11.2, §13): not significant against SSA-VNS (0.16), MS-SLSQP (0.58)
  and, as at 5 degrees, VNS; MS-SLSQP wins at 6,030 (0.0014); PyWake: MS-SLSQP highest at both budgets, 16 layouts
  above the installed block.
- Overclaiming: "(eight cases)" for direct 1-degree claims, "(six cases, five methods)" for the constraint study, the
  projection result rephrased; Table 1 "on the benchmark".
- Precision: 56,540 vs 52,930 reconciled; replay error and its effect on re-evaluations quantified (§9.6, S11.11).
- Timing: exact-gradient ratio 2.3–5.6 (medians, Table S76); initialization CPU median 3.85 s (wall 4.8 s);
  `latex_source/analysis/mpce_numbers.tex` restored to the pipeline values.
- Stale text: 52,930 runs (§3.2), 13-method IEA37 pool (§6.3), GA pointers (§6.1), "favor"; 68-case all-run
  numbers only in §9.6; Threats count "the last two on six and eight cases".
- Formatting: Table "Expanded nine-method benchmark pool" moved to §7.2 and its best method bolded; the IEA37 means
  table without the feasibility columns (counts in its note).

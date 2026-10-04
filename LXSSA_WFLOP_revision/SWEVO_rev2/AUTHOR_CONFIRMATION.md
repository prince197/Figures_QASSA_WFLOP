# Items the authors must resolve before submission

Updated for revision 3 (4 October 2026). Items 1–3 are resolved by this revision; item 2 still needs the final rerun
count. Items 4 and 8 are resolved as far as analysis can resolve them, and the decisions that remain are the
authors'. Items 5–7 need facts or actions that only the authors can supply. Search the sources for `\TBD` to find
every open value; each one prints in red.

1. **Recover the complete research archive. — RESOLVED (revision 3).** Nothing was missing. The 36,190 original
   records, the diagnostic outputs, all model, optimizer and driver modules and all analysis and check scripts are in
   `analysis/`. The complete analysis pipeline was re-executed in the environment of the runs (CPython 3.11.15,
   NumPy 2.4.6, SciPy 1.17.1, pandas 3.0.6):
   - all 91 regenerated or compared outputs are identical apart from time stamps;
   - all checks pass (62 C, 57 X, 22 D, 41 F, 15 S, 57 T);
   - `rev2_analysis.py` reproduces its 4,670 values;
   - reruns of stored runs are bit-identical except for the elapsed time.

   The deposit-ready package is `repro_package_build/` (README, `MANIFEST.sha256`, requirements, `Dockerfile`,
   `CITATION.cff`); see `repro_package_build/verification/verification_summary.txt`. The public deposit is still
   open, see item 7.
2. **Recover full-precision coordinates. — RESOLVED BY AUDIT AND RERUN; final count pending.** All 56,540 stored
   records were audited with the correct site geometry, including Horns Rev 1 (`analysis/rev3_precision_audit.json`,
   supplementary table `tab:S-r3-precision`):
   - 50,826 labels are confirmed from the rounded coordinates;
   - 0 labels are contradicted;
   - 5,714 are undecidable within the rounding bound.

   The undecidable records are rerun deterministically with the 17-digit writer (`analysis/rev3_precision_rerun.py`;
   prespecified in `rev3_precision_manifest.md`). No coordinates were inferred by projecting the rounded layouts.
   **Open:** the lead enters the final result in place of `\TBD{C4 rerun: N bit-identical, labels confirmed/changed}`
   and `\TBD{C4 rerun count}` (Data availability).
3. **Complete benchmark reliability analysis. — RESOLVED (revision 3).** The all-run endpoint (feasible beats
   infeasible, then the objective) is applied to all 68 benchmark cases (`tab:S-r3-allrun`):
   - PSO-VNS vs PSO: 950/348/742 of 2,040 seed pairs, score 0.551, 95% CI [0.516, 0.587];
   - PSO-VNS vs the other six main methods: 0.73–0.90;
   - SSA-VNS vs RSD-VNS: 0.476.

   Feasibility rates and conditional means are reported separately. The Lillgrund and 5D/6D spacing endpoints are
   unchanged.
4. **Review inferential claims. — analyses done; decision remains with the authors.** Revision 3 adds:
   - a null-simulation calibration of the bootstrap TOST (`tab:S-r3-tostaudit`): both primary claims hold after
     calibration;
   - joint (synchronized) seed resampling (`tab:S-r3-syncseed`): the seed-level equivalence of SSA-VNS vs RSD-VNS is
     not robust, but its non-superiority is;
   - imputation sensitivity (`tab:S-r3-imputation`): no verdict changes;
   - a wild-cluster bootstrap for the equal-cluster target (`tab:S-r3-eqclus`): no verdict changes.

   The authors still decide whether to keep the approximate interval-based equivalence sensitivities. As a team they
   should also confirm the two primary equivalence questions (PSO-VNS vs PSO; SSA-VNS vs RSD-VNS non-superiority) and
   that every other equivalence result is labelled exploratory (item A3 of the remaining-work report).
5. **Confirm author details.** Confirm the names, affiliations and corresponding-author e-mail addresses that were
   kept in the manuscript. Do not substitute affiliations from unrelated earlier projects without author approval.
6. **Confirm declarations.** Every author confirms the CRediT roles, the no-grant funding statement and the
   no-conflict declaration. The generative-AI declaration now follows the Elsevier template and names Claude
   (Anthropic) and ChatGPT (OpenAI) and their uses. Confirm that this list of tools is complete and remove
   `\TBD{authors: confirm the list of tools}`. Elsevier asks that AI used in the research process itself (for
   example code for the analyses) also be described in the Methods; decide whether one sentence in Section 5 ("Code,
   data and environment") is needed. The same wording appears in `optA/swevo_back.tex` and
   `declarations_content.tex`.
7. **Confirm submission eligibility and make the public deposit.**
   - Confirm author approval, originality, concurrent-submission status, third-party and data permissions, and any
     statements the submission system requires.
   - Create the versioned public deposit of `repro_package_build/` (for example a GitHub release archived by Zenodo)
     and choose the licence (`LICENSE_SUGGESTION.md`: MIT for code, CC-BY-4.0 for data, suggestion only).
   - Enter the DOI in `CITATION.cff` and replace `\TBD{DOI}` in the Data availability statement (both copies).
   - Supply the volume, issue and pages of Solanki and Deep (2023) if they have been assigned. The entry now reads
     "published online May 18, 2023 (Online First)".
   - Check the live SWEVO Guide for Authors for limits that could not be verified (length, current reference
     style). See `analysis/rev3_refs_report.md` Section 6.
8. **Choose any new experiments or shortening. — mostly done.** Revision 3 ran:
   - the constraint-handling study (2,700 runs; Deb's rules with box clipping or with radial projection;
     `tab:S-r3-constraint`). Conclusions depend partly on constraint handling, as now stated in the text;
   - direct 1-degree optimization (1,200 runs + 1,200 seed-paired 15-degree control runs; `tab:S-r3-fine`);
   - the TOST null calibration;
   - the re-evaluation of the GA and Lillgrund layouts under alternative evaluators.

   Not run, and not claimed: tuned DE, CMA-ES or L-SHADE baselines and larger or full farms. Revision 3 adds the new
   results in short sentences that point to supplementary tables and removes repeated caveats (item E1). The
   authors decide on any further shortening against the journal's current instructions.

**Readiness decision:** the scientific archive is complete and verified. What remains before submission is author
work:
- confirmations (items 4–7);
- the public DOI;
- the final full-precision rerun count;
- a last read of the rebuilt PDFs.

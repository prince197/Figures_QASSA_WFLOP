# Revision 3 — integration brief (manuscript, supplement, back matter)

Read first: SWEVO_rev2/REV3_BRIEF.md (context), SWEVO_rev2/REV3_TRIAGE.md, the remaining-work report text
(/tmp/claude-0/-home-user-Figures-QASSA-WFLOP/83ceea2c-3e62-519f-bca6-caf988ef6cfa/scratchpad/rev_work.txt) and the
analysis reports listed below (they are in LXSSA_WFLOP_revision/analysis/ unless stated).

Canonical sources: SWEVO_rev2/latex_source/ (SWEVO_manuscript.tex -> manuscript_preamble.tex, optA/swevo_front.tex,
optA/sw/*.tex, optA/swevo_back.tex; SWEVO_supplement.tex + optA/sw/supp_theory.tex). Numbers are typed literally in
these files (no \N macros); keep the "% CHECK..." comments attached to their paragraphs (move them with the text,
update the numbers they mention if the text changes; never break the tag format).
Cross-references: the main text cites supplement items as \ref{S-<label>} (resolved through
xref/supplement_labels.tex, refreshed by analysis/build_swevo.sh); the supplement cites main items as \ref{M-<label>}.
New supplementary tables (label tab:S-r3-...) are added to the supplement by agent I5; in the main text cite them as
Table~\ref{S-tab:S-r3-...}. They will resolve after the build.

## File ownership (edit ONLY your files)
- I1: optA/sw/05_setup.tex, optA/sw/06_results.tex, optA/sw/07_ablation.tex
- I2: optA/sw/08_beyond.tex
- I3: optA/sw/03_model.tex, optA/sw/09_robust.tex, optA/sw/10_limits_concl.tex
- I4: optA/swevo_back.tex, declarations_content.tex, SWEVO_declarations.tex, analysis/README_reproduce.md (of latex_source),
      SWEVO_rev2/START_HERE.md, SWEVO_rev2/AUTHOR_CONFIRMATION.md, SWEVO_rev2/REVISION_CHANGELOG.md (append a
      "Revision 3" section; do not rewrite history)
- I5: SWEVO_supplement.tex, optA/sw/supp_theory.tex, latex_source/analysis/ (copy the rev3 *_tables.tex files there)
- I6: optA/swevo_front.tex (abstract, keywords, nomenclature), optA/sw/02_intro.tex, optA/sw/02b_related.tex,
      SWEVO_highlights.txt, SWEVO_cover_letter.tex
Do not edit manuscript_preamble.tex, capture.tex, 04_methods.tex or other files without asking the lead (send a note in
your report instead). Compile only in a private scratch copy of latex_source
(/tmp/claude-0/-home-user-Figures-QASSA-WFLOP/83ceea2c-3e62-519f-bca6-caf988ef6cfa/scratchpad/int_<you>/, e.g.
`cp -r SWEVO_rev2/latex_source <dir> && cd <dir> && sh analysis/build_swevo.sh`); undefined S-tab:S-r3-* references are
expected until I5 has finished — everything else must compile with 0 errors and no new overfull boxes. No git.

## Length (reviewer E1; the authors' own earlier limit was 40 pages)
The manuscript is 50 pages now. Every agent must make its files NET SHORTER or at most equal in length: put new
results in one or two sentences with a pointer to the supplementary table; remove repeated full statements (the
PSO-VNS vs PSO equivalence verdict, the spacing caveat, the 1-degree caveat, "added after the primary analysis",
"exploratory") where they repeat what an earlier section already says; keep numbers in Results, implications in
Discussion, limitations in Threats, short bounded findings in Conclusions. Report your before/after word counts.

## Decisions (apply them consistently)
D1 Archive / reproducibility (R2, D, A2, A8): the complete archive EXISTS (36,190 original + 11,640 revision records +
   the rev3 records listed below, all modules and analysis scripts) and was re-executed for this revision: every
   regenerated table/macro/summary/figure is identical apart from time stamps (91 files), all checks pass (62 C, 57 X,
   22 D, 41 F, 15 S, 57 T), rev2_analysis reproduces its 4,670 values, reruns of stored runs are bit-identical except
   elapsed time. Delete every statement that records/modules/historical checks are missing, "author-reported",
   "historical", "cannot be rerun". Deposit-ready package: repro_package_build/ (README, MANIFEST.sha256,
   requirements, Dockerfile, CITATION.cff); public DOI pending author action — say exactly that. Source text:
   repro agent report = section "Suggested replacement text" in repro_package_build/README.md and
   repro_package_build/verification/verification_summary.txt.
D2 Per-evaluation cost (numbers audit item 1; sites report B6): quote the MEDIAN per evaluation (0.23–0.45 ms for the
   metaheuristics, MS-SLSQP 0.53 ms, 2.3x PSO) and add that the stored MS-SLSQP times jump at N >= 8 (16–22 s per run,
   3.4 ms/eval) because of thread contention in that batch (re-timed single-threaded: 4–7 s), earlier runs of the same
   method 3.3–4.0 s; Table S-cost reports means (say "mean"). Initialization (uncharged): median 4.8 s (max 15.5 s)
   per run on the six largest cases, more than a 6,030-evaluation PSO-VNS run (2.8 s). 101 CPU-hours for the 36,190
   original runs, run as four parallel processes. (analysis/rev3_sites.json; tables tab:S-r3-time, tab:S-r3-time-sites)
D3 Number fixes (analysis/rev3_numbers_report.md): Lillgrund time ratio 2.2 (not 2.3); Horns Rev wake-loss order
   "10.39% (ours) vs. 11.11% (PyWake)"; tab:S-eqclus cluster means -0.070/-0.043/-0.234 (superseded by tab:S-r3-eqclus);
   Lillgrund area 1.084 km^2 (enlarged boundary used); "PSO-VNS IEA37 records missing" is false; move CHECK-DIR [F23]
   comment back after its paragraph in 08_beyond.tex.
D4 All-run endpoint on the 68 cases (B1; analysis/rev3_inference.json, tab:S-r3-allrun): feasible beats infeasible,
   then objective; PSO-VNS vs PSO 950/348/742 of 2,040 seed pairs, score 0.551 (95% case-bootstrap CI [0.516, 0.587],
   two-stage [0.511, 0.591]; Holm sign test p = 4.7e-7): better more often than not, by amounts within the margin —
   consistent with the seed-level interval; the conditional case-mean test is not significant (p = 0.30). Against the
   other six main methods: scores 0.73–0.90, all intervals above 0.5; rank order unchanged except SSA and MS-SLSQP swap.
   SSA-VNS vs RSD-VNS 0.476 [0.452, 0.500] (two-stage [0.445, 0.507]). The qualification paragraph of 6.2 must now say
   the all-run endpoint is reported for the 68 cases (not only for Lillgrund).
D5 TOST calibration (B3; tab:S-r3-tostaudit): bootstrap p_TOST and the 90%-interval rule agree for all 18 pairs at all
   levels (0 discordant). Null simulation at the margin: empirical sizes 4.4–7.4% for PSO-VNS vs PSO (case and seed),
   4.7–6.0% at the seed level for SSA-VNS vs RSD-VNS, but 15.0% at the case level for its non-superiority test (one
   outlying case, Data Set I, r = 500 m, N = 10, d = -1.33 pp); null-calibrated p: 0.015 (PSO-VNS vs PSO, case) and
   0.007 (SSA-VNS non-superiority, case) — both primary claims hold after calibration. Replace "its nominal 5% error
   rate has not been verified by null simulation" accordingly.
D6 Synchronized seeds (B4; tab:S-r3-syncseed): seed k gives identical initial populations across cases with equal N
   (rescaled by r; identical for the two data sets). Joint seed resampling widens seed-level intervals by up to 1.9x;
   PSO-VNS vs PSO unchanged [-0.026, -0.010]; SSA-VNS vs RSD-VNS [-0.002, 0.050]: its seed-level equivalence and
   seed-level significance are not robust, its non-superiority is. Update Table 7 (tab:equiv-main) footnote b.
D7 Imputation (B5; tab:S-r3-imputation): dropping imputed cases or Pratt zero handling changes no verdict.
D8 Equal-cluster (B7; tab:S-r3-eqclus): recomputed from records (agrees with the rev2 reconstruction within 0.0005);
   wild-cluster bootstrap for the equal-cluster target: PSO-VNS - PSO [-0.062, 0.036] (p_TOST 0.108), SSA-VNS - RSD-VNS
   [-0.067, 0.082]; no verdict changes.
D9 Site pools (B2; tab:S-r3-hr-pool, tab:S-r3-iea-pool) and alternate evaluators (C3; tab:S-r3-hr-eval,
   tab:S-r3-lg-eval): use the sentences of the sites report (analysis/rev3_sites.json; the "Suggested sentences"
   reproduced in analysis/rev3_suggested_text.md, which also has the inference, constraint and reproducibility texts):
   HR pool of 11: PSO-VNS best at 6,030 (GA second, p_Holm 0.53), GA best at 30,030 (beats all ten, p_Holm <= 1.1e-4);
   GA's lead holds at 1 deg (+0.48 GWh/yr, p = 2.8e-4) and under PyWake NOJ (+0.43, p = 4.2e-4); with 1-deg bins none of
   the 961 feasible optimized layouts exceeds the installed block (22 with 5-deg bins). IEA37 pool of 13: an
   exact-gradient method first in every scenario/budget, better than every gradient-free method in all 30 seed pairs
   (p_Holm 2.2e-8); among gradient-free methods PSO-VNS leads with 16 turbines, GA with 36. Lillgrund at 1 deg: installed
   block -1.65%, optimized layouts -1.14% on average (PSO-VNS at 30,030: -1.69%); PSO-VNS keeps the highest mean at 30,030
   (112.64 GWh/yr) but is no longer significantly better than VNS or MS-SLSQP (p_Holm 0.58); at 6,030 MS-SLSQP has the
   higher mean (and wins the paired test, p_Holm 0.0014). The Lillgrund feasibility advantage of PSO-VNS over PSO (17/30
   vs 0/30; 30/30 vs 19/30) is unaffected (feasibility does not depend on the evaluator).
D10 Constraint handling (C1, R5; tab:S-r3-constraint, tab:S-r3-constraint-cases; analysis/rev3_constraint.json and
   rev3_constraint_manifest.md): six cases, PSO-VNS, PSO, GA, SSA-VNS, DE, seeds 31–60, 2,700 runs; variants (i) current
   penalty + box clipping, (ii) Deb's feasibility rules on dimensionless violations (g_b/r^2, g_s/l_min) + box clipping,
   (iii) Deb's rules + radial projection onto the circle. Variant (i) reproduces 60/60 stored runs bit for bit. (iii):
   rank order unchanged (PSO-VNS < PSO < GA < SSA-VNS < DE), all methods 0.46–0.80 pp better, DE feasibility 33.3% -> 80.0%,
   PSO-VNS keeps significant all-run advantages over GA, SSA-VNS, DE (0.683, 0.850, 0.844) but not over PSO (0.561, ns;
   conditional ΔL -0.095 pp [-0.162, -0.037]). (ii): PSO feasibility 98.9% -> 73.3% (mostly boundary violations: a turbine
   clipped to the square near a tangent point is cheap under normalized violations), PSO-VNS conditional loss +0.23 pp,
   all-run advantage over GA (0.517) and SSA-VNS (0.583) no longer significant; GA first by case-level rank (1.83 vs 2.17),
   PSO-VNS first by run-level rank (2.19 vs 2.28). Pre-registered survival rule: fails for (ii), fails only on PSO for
   (iii). Message: the comparison depends partly on constraint handling — on how boundary and spacing violations are
   traded off under box clipping; conclusions are for the implemented penalty/repair, and are robust to the projection
   repair except for PSO-VNS vs PSO. Replace every "no normalized-penalty/alternative-repair/lexicographic experiment was
   run/supplied".
D11 Direct 1-degree optimization (C2, R8; tab:S-r3-fine; analysis/rev3_fine.json, rev3_fine_manifest.md): 8 cases (both
   data sets; r=500 N=10; r=750 N=6, 12; r=1000 N=15), PSO-VNS, PSO, GA, MS-SLSQP, RSD-VNS, seeds 31–60, 6,030 evaluations
   of the 1-degree objective (1,200 runs) + a seed-paired 15-degree control arm (1,200 runs, same seeds). Results (rerun
   rev3_fine_analysis.py when analysis/rev3_fine15_s0of1.csv has 1,200 rows, and take the final numbers from
   rev3_fine.json; preview): direct 1-deg optimization: average ranks PSO-VNS 1.12, PSO 2.62, GA 3.38, RSD-VNS 3.62,
   MS-SLSQP 4.25; PSO-VNS - PSO -0.194 pp (seed 90% [-0.234, -0.152], case 90% [-0.272, -0.114], p_W = 0.0078 over 8
   cases) — beyond the margin; all-run scores of PSO-VNS 0.70–0.87 vs every method (all p_Holm < 1e-8). Direct 1-deg
   optimization lowers the 1-deg wake loss relative to re-evaluated 15-deg layouts by 0.3 (MS-SLSQP) to 1.2 pp (RSD-VNS).
   So the leader (PSO-VNS) does not change when the search itself uses 1-deg bins, and its advantage over PSO grows.
   Replace "Which method would lead if the search itself used 1-deg bins was not tested" and related caveats.
D12 Coordinate precision (C4; tab:S-r3-precision; analysis/rev3_precision_audit.json, rev3_fullprec_*.csv when present):
   all 56,540 stored records audited with the correct site geometry (incl. Horns Rev 1 via hornsrev_model.py): 50,826
   labels confirmed from the rounded coordinates with slack larger than the rounding bound, 0 contradicted, 5,714
   undecidable within the rounding bound (constraints active at the optimum). Full-precision deterministic reruns of the
   undecidable records are in progress: write the sentence with the placeholder \TBD{C4 rerun: N bit-identical, labels
   confirmed/changed} — the lead fills it in.
D13 References (analysis/rev3_refs_report.md, rev3_refs_corrections.tex): apply the corrected bibitems (Solanki2023:
   "published online May 18, 2023 (Online First)"; Thomas2023 without "no. 5"; IEA37repo comment; supplement PyWake
   order; delete uncited Benavoli2016 and Cleghorn2018 from the supplement list); renumber both bibliographies in order
   of first citation (I4: main, I5: supplement) — for the main list use analysis/swevo_bib.py logic adapted to the
   SWEVO_rev2 paths (compile once, read \citation order from the .aux).
D14 AI declaration (A4): "During the preparation of this work the authors used Claude (Anthropic) to assist with code
   for the analyses, the analyses and experiments of the revisions, consistency checks and the drafting and editing of
   the text, and ChatGPT (OpenAI) to assist with the manuscript revision, archive audit, document compilation and
   package preparation. After using these tools, the authors reviewed and edited the content as needed and take full
   responsibility for the content of the published article. \TBD{authors: confirm the list of tools}" (same text in
   declarations_content.tex / SWEVO_declarations.tex).

## Report to the lead
Files changed; for each reviewer item (R1–R10, B1–B7, C1–C4, E1–E6) what you changed and where; before/after word
counts; compile result of your scratch build (errors, undefined refs other than the expected S-tab:S-r3-*, overfull).

## Update (final C2 numbers, control arm complete; analysis/rev3_fine.json, rev3_fine_tables.tex)
Direct 1-deg optimization (seeds 31–60): ranks PSO-VNS 1.12, PSO 2.62, GA 3.38, RSD-VNS 3.62, MS-SLSQP 4.25; PSO-VNS - PSO
-0.194 pp (seed 90% [-0.234, -0.152], case 90% [-0.272, -0.114], p_W = 0.0078, 8 cases). Seed-paired 15-deg control arm
(same seeds, optimized at 15 deg): ranks PSO-VNS 1.88, PSO 2.25, GA 2.62, MS-SLSQP 4.00, RSD-VNS 4.25; re-evaluated at 1 deg:
PSO-VNS 2.25, MS-SLSQP 2.62, PSO 2.88, GA 3.12, RSD-VNS 4.12 (leader unchanged, as for the stored seeds 1–30). Optimizing at
1 deg instead of 15 deg lowers the 1-deg wake loss on the same seeds by 1.12 pp (PSO-VNS), 1.01 (PSO), 1.09 (RSD-VNS),
0.84 (GA) and 0.32 pp (MS-SLSQP). The leader does not change; PSO-VNS's advantage over PSO is about twice that of the
re-evaluated layouts (-0.086 pp, control arm at 1 deg) and lies beyond the 0.05 pp margin.

## Lead's final-pass list (collected from agent reports)
- 0.33 -> 0.32 pp (MS-SLSQP gain from 1-deg optimization) wherever it appears.
- 7.3/8 vs 9.6: the 68-case all-run numbers appear in both; keep them in one place (9.6 or 7.3) and point from the other.
- tab:ga-main (Table 13) is not cited: cite it in 06_results (GA paragraph).
- 05_setup: check qualification paragraph mentions the 68-case all-run endpoint; 6.3 IEA37 sentence vs 13-method pool.
- 10_limits_concl: "Lillgrund ... 5-deg bins only" and "No dimensionless-violation ... experiment was run" (I3).
- 02b_related Table 1: all-run outcome now on the 68 cases (I6).
- Supplement: tab:S-r3-iea-pool caption (generated rev3_sites_tables.tex) "cannot re-verify" -> align with D12 after C4.
- Precision text: explain 56,540 audited records vs 52,930 runs in Table 3 (superseded / archived rows).
- Fill the C4 \TBD placeholders (08_beyond 9.6, supplement S11.11 and S-archive, swevo_back, declarations_content,
  README_reproduce, AUTHOR_CONFIRMATION, REVISION_CHANGELOG).
- build_swevo.sh line 2 stale comment; refresh xref/*_labels.tex in canonical; regenerate SWEVO_manuscript_full.tex.
- Bibliography order rerun at the end (analysis/rev3_bib_order.py) for main and supplement.
- Data availability: "release 1.0.0" vs package 1.0.0-rc — align.

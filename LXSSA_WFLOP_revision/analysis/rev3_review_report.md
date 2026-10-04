# Adversarial review of the integrated revision-3 SWEVO manuscript (4 Oct 2026)

Build: scratch copy of `SWEVO_rev2/latex_source`, `sh analysis/build_swevo.sh`. Result: 0 errors, 0 undefined references, no "??", 0 overfull boxes. Page counts: manuscript 50, supplement 69, response 5, cover letter 1.
Sources checked: `analysis/rev3_{inference,sites,constraint,fine,precision_audit}.json`, `rev3_precision_rerun_summary.json` (C4 reruns, still in progress, 22:11), `rev3_*_analysis.log`, `rev2_summary.json` (regenerated in scratch; it reproduces exactly), `repro_package_build/` and `verification/verification_summary.txt`.
Page and line numbers below refer to the scratch-built PDFs. M = manuscript, S = supplement, R = response letter, CL = cover letter.

Nearly all numbers that revision 3 added match their JSON sources: all-run scores, TOST sizes, calibrated p values, synchronized-seed intervals, equal-cluster intervals, Horns Rev and IEA37 pools, the GA re-evaluation, the constraint study, the direct 1° study, the precision counts, and 52,930 / 16,740 / 5,100. The problems lie in the claims and wording around those numbers. They are listed most severe first.

---

## A. Severe: claims that the evidence contradicts

1. **"Reruns of stored runs are bit-identical except for elapsed time" does not hold for SLSQP-based runs.**
   - Where the claim appears: M 6.4 p22 l.722–723; S9 p45 ("reproduce the records bit for bit in every field except the elapsed time"); S12.1 p66; R p1 (premise paragraph); CL ¶4.
   - What the C4 reruns show so far (`rev3_precision_rerun_summary.json`): of 1,679 reruns, 1,183 are bit-identical, 419 match only to float noise, 11 are "near" and 66 diverged.
   - All 66 diverged reruns are SLSQP-based: MS-SLSQP rows of `fresh_hr16` (27), `mpce_iea` (29) and `rev2_grad` (10). Relative objective differences reach 2.3%.
   - One IEA37 36-turbine MS-SLSQP rerun (seed 22) ends infeasible where the stored run is feasible (`strict_changed_undecided: 1`).
   - The `rev3_sites` MS-SLSQP probe also gives different objectives (e.g. 97,714 vs. stored 97,810.9).
   - The package verification demonstrated bit-identity only for 4 PSO-VNS and SSA-VNS reruns (plus 60/60 constraint-driver runs and 62 diagnostic runs).
   - **Fix:** write "bit-identical for the gradient-free methods; SLSQP-based runs are platform dependent (BLAS/SciPy build) and reproduce only to float noise or diverge."
   - Extend the C4 placeholder (M p37 l.1177; S11.11 p64; S12.2 p66; R p4) to cover this outcome, e.g.: "N bit-identical, K float noise, J not reproducible (SLSQP, platform dependent; labels kept as stored); labels confirmed X / changed 0."

2. **The "deposit-ready" package does not contain the revision-3 data.**
   - `repro_package_build/dist/SWEVO_repro_1.0.0-rc.zip` was built at 19:40.
   - It contains no revision-3 run records or outputs: `rev3_constraint_{pen,deb,proj}.csv`, `rev3_fine_s*`/`rev3_fine15_*` CSVs, `rev3_{inference,constraint,fine}.json`, the `rev3_*_tables.tex` files, the `rev3_fullprec_*` CSVs and `rev3_precision_rerun_summary.json` are all missing. Only scripts, manifests and the early JSONs are included.
   - R2 (p1) nevertheless says it "contains all run records (… 5,100 new runs, plus the full-precision reruns) … verified from a clean unpack".
   - Version mismatch: Data availability (M p45 l.1485; declarations) and R2 say "release 1.0.0", but the package is `1.0.0-rc` (`CITATION.cff`, README).
   - **Fix:** rebuild and re-verify the package after C4 finishes, then use one version string everywhere.

3. **6.4 says the original-code calibration was re-executed; S12 says it was not.**
   - M 6.4 p22 l.731–737 is headed "Evaluator and code validation (re-executed for this revision)" and includes "reproduces the recorded result distributions of the original code (2 of 24 Mann–Whitney tests reject)".
   - S12.1 p66 says "Not re-executed were … the calibration reruns of the original code". The verification summary confirms that only `calibration_vs_recorded.csv` was regenerated from stored outputs.
   - **Fix:** move that sentence out from under the "re-executed" heading, or say "recomputed from the stored calibration runs".

4. **The Lillgrund 1° sentence is wrong as worded.**
   - Where: M 9.4 p35 l.1100–1101 and S11.7 p61: "no longer significantly better than VNS or MS-SLSQP (pHolm = 0.58)".
   - At 5° and 30,030 evaluations, PSO-VNS was already not significantly better than VNS (pHolm 0.058; `rev2_summary.json` lillgrund/30030/BVNS). So "no longer" is wrong for VNS.
   - At 1° it also loses significance against SSA-VNS: pHolm 0.16, compared with 0.016 at 5° (`rev3_sites.log` [lg]).
   - **Correct text:** "…but is no longer significantly better than SSA-VNS or MS-SLSQP (pHolm = 0.16, 0.58) and, as at 5°, not than VNS".
   - The main text also omits that at 6,030 MS-SLSQP wins the paired test at 1° (pHolm = 0.0014; S11.7 states it). Add it.
   - Omitted PyWake results: under PyWake NOJ at 1°, MS-SLSQP has the highest mean at both budgets (112.61 vs. 112.42 and 113.01 vs. 112.93 GWh/yr). Also, 16 layouts exceed the installed block, not 1 (`lg_eval.above_installed_total`). Yet M p35 l.1097 reads as if only one layout exceeds it under every evaluator.
   - M 11.2 p42 l.1361–1362 and M 13 p45 l.1459–1460 say "the Lillgrund AEP lead of PSO-VNS is not significant with 1° bins". PSO-VNS still beats PSO, SSA, LX-SSA, DE and RSD-VNS at 1°. Write "its lead over the runners-up (SSA-VNS, VNS, MS-SLSQP)".

5. **Threats p44 l.1435–1437 miscounts the case coverage.**
   - Text: "The Lillgrund, spacing, GA, exact-gradient, constraint-handling and direct 1° studies were added …, the last three on six to eight cases."
   - The exact-gradient study ran on IEA37 (two scenarios), not on six to eight cases.
   - **Fix:** "…the last two on six and eight cases".

6. **The response letter makes statements that are now false.**
   - R1 p1: "The only remaining [TBD] markers are … DOI, AI tool list, cover-letter confirmation". The C4 placeholders are also still open: M p37, M p45 (Data availability), S p64, S p66, R p4 C4 and the declarations.
   - R10 and A6 claim the bibliographies are in order of first citation. In the main list, Solanki2023 is cited first (p3 l.48) but numbered [26], while Mladenović 1997 and Hansen 2001, cited on l.49, are [24] and [25] (checked in `SWEVO_manuscript.aux`). Rerun `rev3_bib_order.py`. The supplement list (23 references) is in correct order.
   - A7 says the abstract has 239 words; a count from the source gives about 242. Still within the 250 limit.
   - Revision naming is inconsistent:
     - R calls this "Revision 3".
     - The supplement calls S10 "the first revision" and S11 "the second revision" (S intro, S11 p56).
     - The package README calls the 11,640 runs "revision 2".
     - Data availability says "first revision … this revision".
     - Pick one scheme.

## B. Overclaiming

7. **The 8-case direct 1° study is generalized.**
   - Abstract p1 l.23–26: "ranks first, also when optimizing with 1° bins … but reaches 0.19 pp with 1° bins".
   - The same generalization appears in intro item 4 (p4 l.84–86), Table 1 rows 4–5 (p5), 11.1 (p41 l.1329–1330), 13 (p44 l.1450–1455), and in CL ¶3 ("PSO-VNS stays first when the search itself uses 1° bins").
   - The study covers 8 cases, 5 methods and seeds 31–60, and is exploratory per `rev3_fine_manifest.md`.
   - "Reaches 0.19 pp with 1° bins" also blurs the 1° re-evaluation result (−0.043 pp, Table 7) with direct optimization (−0.194 pp).
   - **Fix:** write "when optimized directly with 1° bins (eight cases)" in each place.

8. **The 6-case constraint study is generalized, and "except PSO-VNS versus PSO" is misleading.**
   - Where: abstract p1 l.31–33, intro item 5 (p4 l.91–94), Table 1 row 3, Threats p43 l.1395–1398, and the highlight "Conclusions depend partly on…" (acceptable as is).
   - Add "(six cases, five methods)".
   - Under projection, the run-level order is unchanged, including PSO-VNS ahead of PSO (conditional ΔL −0.095 pp [−0.162, −0.037]). What fails is the all-run sign test (score 0.561, pHolm 1.0).
   - Also omitted: under projection PSO-VNS is feasible less often than PSO (98.9% vs. 100%, Table S77), and the case-level rule order changes (DE 3.83 ahead of SSA-VNS 4.67).
   - **Suggested wording:** "radial projection keeps the ranking, but PSO-VNS's all-run advantage over PSO is no longer significant". M 10.4 and 13 are already close to this.

9. **Table 1, last row ("the leader holds with 1° bins (Sections 9 and 10)") overgeneralizes.**
   - It is false for Horns Rev 1 at 6,030 evaluations: at 1°, GA has 137.13 GWh/yr, PSO 137.04 and PSO-VNS 136.99 (no significant differences).
   - It is false for Lillgrund at 6,030: MS-SLSQP 112.10 vs. PSO-VNS 111.97, MS-SLSQP better (pHolm 0.0014).
   - **Fix:** limit the statement to "on the benchmark".

10. **Several other sentences are imprecise.**
    - Abstract p1 l.28–29: "GA performs better … and with 36 IEA37 turbines". Add "among gradient-free methods"; the exact-gradient methods beat GA.
    - Intro item 4 p4 l.89 says PSO-VNS "does not lead … on IEA37", but 9.5 shows it leading the gradient-free methods with 16 turbines. Qualify.
    - CL ¶2 says the conclusions are tested "on Horns Rev 1, Lillgrund and IEA37 cases with a real-coded GA and exact-gradient baselines". GA was not run on Lillgrund (S10.5), and the exact-gradient methods ran only on IEA37.
    - CL ¶3 lists "a full-precision verification of the feasibility labels of the stored runs" as done. It is pending, and cannot be complete for the SLSQP runs (item 1).
    - Threats p44 l.1410: "Except in the direct 1° study … the search used discretized roses". A 1° rose is also discretized; write "coarse (15°/5°) roses".

## C. Smaller inconsistencies and stale text

11. **Precision audit (M 9.6 p37–38; S11.11).**
    - The text does not reconcile the 56,540 audited records with the 52,930 runs in Table 3.
    - The audit covers 44,900 original-era records (including 3,210 superseded and Horns Rev pre-fix rows) plus 11,640 revision records. The 5,100 revision-3 runs are not audited because they are stored with 17 digits. Add one sentence saying so.
    - The replay error is understated. M p38 l.1179 cites only "Lillgrund: up to 0.225 GWh/yr". Per `rev3_precision_audit.json` and the per-record CSV:
      - maximum 0.94 pp on the main benchmark (fresh_grid, LX-SSA) and 2.70 pp at 120,030 evaluations;
      - 0.53 GWh/yr on Horns Rev;
      - 905 of 51,810 feasible records (250 of 27,569 main-comparison layouts, mostly SSA, LX-SSA and MS-SLSQP) replay off by more than 0.05 pp.
    - Every re-evaluation (1°/5°, Gaussian, cubic, PyWake) uses these rounded coordinates. State the caveat in the Section 10 intro.
    - Align leftover wording with C4 and D1:
      - Table S72 footnote: "cannot re-verify".
      - Table S81 title: "from the supplied run records".
      - Table S81 note: "Rounded archival coordinates do not independently reproduce all labels".

12. **Timing.**
    - M 9.5 p36 l.1143–1144 says the exact-gradient methods take "about 3–6 times as long as PSO-VNS (Table S76)". The medians in Table S76 give 3.3–5.6 (MS-SLSQP) and 2.3–3.5 (PSO-SLSQP), as S12.4 p67 states; "3–6" is a ratio of means. Use 2.3–5.6, or say "means".
    - M 9.7 l.1192 gives the same comparison against GA. Keep one of the two.
    - The initialization median of 4.8 s (M 9.7; S11.8) is wall time on a shared machine (load average 12.7 on 4 cores, `rev3_sites.json` init_time). The CPU median is 3.85 s (maximum 9.7 s). The conclusion (longer than a 2.8 s PSO-VNS run) still holds, but quote CPU time or say "wall time under load".
    - "Earlier runs of the same method 3.3–4.0 s" matches N ≥ 10. For N ≥ 8 the range is 2.95–4.02 s (fresh_grid medians).
    - `latex_source/analysis/mpce_numbers.tex` contains hand-edited timing macros (0.24/0.47/1.82/7.5) where the pipeline gives 0.23/0.45/0.53/2.3. Several `mpce_tab_*.tex` copies are also hand-edited (captions, typography; verification summary item a). The macros are unused, but this conflicts with "identical to those used in the manuscript". Regenerate them or note it.

13. **Stale counts and pointers.**
    - M 3.2 p10 l.321: "more than 40,000 optimization runs". The total is 52,930.
    - M 6.3 p21 l.711: "the ten methods of the budget study run on them". The IEA37 pool has 13 methods.
    - M 6.1 p17 l.559: GA "analyzed as a ninth method in Section S10.2". GA is now also in Table 13 and Tables 10–12.
    - M 6.2 p19 l.587 points to "Section 7.2" for the 68-case all-run endpoint, while 7.3 points to 9.6.
    - Spelling: M 9.6 "favour" vs. "favor" elsewhere.
    - Check counts: S9 and Table S55 list "checks T01–T13", while M 6.4 says "57 theory checks". Both are true (13 IDs, 57 conditions), but say it once consistently.

## D. Repetition that could be cut (E1). The 10 largest passages; keep the first-named location

1. **68-case all-run endpoint:** 7.2 p25 l.817–819; 7.3 p28 l.854–856; 8.3 p31 l.972–974 (0.476 with both CIs); 9.6 p36–37 l.1159–1165 (repeats 0.551, 0.73–0.90, the swap and 0.476); also S11.1. Keep only 9.6, and point to it from 7.2.
2. **Constraint-handling verdict:** abstract; intro item 5 p4; Table 1; 10.4 p39–40; 11 intro p40 l.1297–1299; Threats p43 l.1394–1400; 13 p45 l.1461–1464; Table 15 item 9. Threats and Conclusions repeat 10.4 almost word for word.
3. **Direct 1° result (−0.19 pp, leader unchanged):** abstract; intro item 4 p4 l.84–86; Table 1; 7.3 p28 l.866–869; 10.5 p40; 11.1 p41 l.1328–1332; 11.3 p42 l.1364–1367; 13 p44 l.1450–1455.
4. **PSO-VNS vs. PSO equivalence verdict:** intro item 4; Table 1; 6.2 p20 l.652–663; 7.3 p26–28; 8.1 p28 l.889–891; 8.3 p31 l.962; 11.1 p41 l.1321–1327; 13 p44 l.1451–1454.
5. **GA's Horns Rev lead at 30,030 (also at 1° and PyWake):** abstract; intro item 4; Table 1; 7.2 p25 l.806–808; 9.3 p34; 11.2 p42 l.1360–1361; Threats p44 l.1421–1422; 13 p45 l.1457–1459.
6. **"VNS phase is essentially one compass-search descent"** (39 of 40 runs, 21.5% / 4.7% / 0.5% repeated): 5.2 p16 l.525–529 and 8.1 p28 l.891–904. Also intro item 3 p3 l.75–76; 11.2 p41–42 l.1347–1350; Threats p43 l.1404–1406.
7. **SSA phase: no gain at 4D/5D, gain at 6D:** intro item 3; Table 1; 8.3 p30 l.927–929 and p31 l.975–976; 10.3 p39 l.1235–1236; 11.1 p41 l.1317–1318; 11.2 p41 l.1343–1346; 13 p44 l.1446–1448.
8. **Feasibility advantage on the two blocks:** 7.2 p25 l.819–821; 9.3; 9.4; 11.2 p42 l.1355–1356; Threats p44 l.1416–1417; 13 p45 l.1456.
9. **Reproducibility statement:** intro item 1 p3 l.63–65; 6.4 p22; Threats p42 l.1389–1390; Data availability p45; S9; S12.1. Keep it in 6.4 and Data availability.
10. **Old PSO setting, clamping and large steps:** 4.1 p13 l.427–430; 7.1 p23 l.769–781; 11.1 p41 l.1307–1313; Threats p42–43 l.1393–1394.

Also duplicated: the 1.322% → 2.335% numbers (3.1 p9 l.292–294 and 10.5 p40 l.1267), and the gradient-cost argument (2.1 p6 l.141–145, 164–167 and 9.5 p36 l.1118–1125).

Length: 50 pages against the authors' own 40-page limit. R E1 claims "the manuscript stays at 50 pages".

## E. Formatting

- **[TBD] markers.**
  - Expected C4 placeholders: M p37 l.1177; M p45 l.1487 (Data availability); S11.11 p64; S12.2 p66; R p4 (C4 row); declarations ("C4 rerun count").
  - Author-only markers: M p45 l.1491 (DOI); M p46 l.1501 and declarations (AI tool list); CL ("confirm not under consideration…").
- **"??" and undefined references:** none. Overfull boxes: none. The supplement log warns "Text page 48 contains only floats".
- **Bibliography:** main references [24]–[26] are out of first-citation order (item 6). The supplement is in order.
- **Highlights:** 64–76 characters each (≤ 85). **Abstract:** about 242 words (≤ 250). **Keywords:** 7.
- **Table legibility.**
  - Table 12 (p37): columns run together ("395,031.8 30/30 3"). Increase `\tabcolsep`, or move "Feas." to a footnote, since all values are 30/30 except 6.
  - Table 13 is first cited in 7.2 (p25) but printed on p37 in 9.6; move it next to 7.2.
  - Table 13 does not bold the best method, although R E2 says the tables mark it.
  - Table 7 (p27) has no row for the direct 1° result, which the abstract headlines; optionally add it as an exploratory row.
- **Odd line breaks:** none beyond normal justification in the PDF text.

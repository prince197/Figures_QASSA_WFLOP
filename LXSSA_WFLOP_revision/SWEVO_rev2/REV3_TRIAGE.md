# Triage of "SWEVO revision 2 — remaining work required before submission" (4 Oct 2026)

Verdict per item: G = genuine and addressed here completely; GP = genuine, addressed here as far as possible
(remainder needs the authors); A = genuine but author-only (cannot be done by analysis or editing);
N = not genuine / based on a wrong premise (explained). Status is filled in when the work is done.

Premise check. The report (and the revision-2 changelog, START_HERE, AUTHOR_CONFIRMATION, Section 6.4, Data
availability, Section S11, analysis/README_reproduce.md) assumes that the 36,190 original run records and fifteen core
modules are missing. This premise is wrong for this project: all of them are in the repository (analysis/), and the
revision-2 analysis script reproduces its own summary exactly (4,670 values) when run on them here. Every item of
Section B and C4 is therefore doable, and the "historical checks cannot be rerun" statements must be corrected.

| Item | Concern | Verdict | Work done (filled in) |
|---|---|---|---|
| R1 / A1 | declarations, author facts | A | declarations complete (Elsevier AI template); facts need author confirmation |
| R2 / A2 / D | public availability, reproducibility package | GP | release 1.0.0 built with all revision-3 records and verified from a clean unpack (checksums, full regeneration, determinism); public DOI and licence = author action |
| R3 / B2 | GA and exact gradients in principal site tables, recomputed pools | G | Horns Rev 1 pool of 11, IEA37 pool of 13 recomputed (`tab:S-r3-hr-pool`, `-iea-pool`) |
| R4 / B1, B5 | all-run outcome on the 68 cases; imputation sensitivity | G | all-run endpoint on 68 cases (`tab:S-r3-allrun`); imputation changes no verdict |
| R5 / C1 | penalty/repair confounding | G | 2,700-run constraint-handling study; partial dependence reported |
| R6 | Proposition 1(c) | done in rev2 | verified (theory checks T01–T13 pass) |
| R7 / B3, B4, B7 | TOST calibration, synchronized seeds, equal-cluster wild bootstrap | G | null-calibrated TOST, joint seed resampling, wild-cluster bootstrap; one equivalence found non-robust and reported |
| R8 / C2, C3 | re-evaluation vs re-optimization | G | 1,200 direct 1° runs + 1,200-run 15° control; GA and Lillgrund layouts re-evaluated at 1° and with PyWake |
| R9 / B6 | charged evaluations vs elapsed time | G | per-evaluation medians, initialization CPU time, total CPU time |
| R10 / E1–E4 | presentation | G | repetition pass, number fixes, bibliography in first-citation order |
| C4 | full-precision coordinates | G | 56,540 records audited; 5,714 rerun at full precision: 0 labels changed, 920 SLSQP-based undecidable |
| A3 | primary equivalence questions agreed | A | open (team agreement) |
| A4 | AI declaration accuracy | A (wording updated to name all tools used) | authors confirm the tool list |
| A5 | Solanki & Deep 2023 volume/pages | GP | entry corrected to Online First; volume/pages not yet assigned |
| A6 | reference audit | GP | 88 entries checked, 3 corrected; list in first-citation order |
| A7 | SWEVO Guide for Authors | GP | checkable items pass; live length limit to be checked by authors |
| A8 | AUTHOR_CONFIRMATION items | A / partly resolved (archive and precision items) | items 1–3 resolved; 4–7 author-only |
| E5, E6 | highlights/cover letter; response letter | G | highlights and cover letter revised; 5-page response letter written |

Summary: 18 rows covering 36 item codes (R1–R10, A1–A8, B1–B7, C1–C4, D, E1–E6). Genuine and addressed completely: 10 rows (R3–R5, R7–R10, C4, E5/E6, and R6 already
done in revision 2 and verified). Genuine, addressed as far as possible: 4 rows (R2/A2/D, A5, A6, A7; the remainder is
the DOI, licence, unassigned volume/pages and the live guide). Author-only: 4 rows (R1/A1, A3, A4, A8 items 4–7).
Not genuine as stated: the shared premise that the original records and core modules were missing (they are in
`analysis/`); every item that depended on it was done.

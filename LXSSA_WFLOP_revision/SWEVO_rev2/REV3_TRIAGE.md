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
| R1 / A1 | declarations, author facts | A | text already complete; author confirmation needed |
| R2 / A2 / D | public availability, reproducibility package | GP | complete package built and verified; public DOI = author action |
| R3 / B2 | GA and exact gradients in principal site tables, recomputed pools | G | |
| R4 / B1, B5 | all-run outcome on the 68 cases; imputation sensitivity | G | |
| R5 / C1 | penalty/repair confounding | G | targeted constraint-handling study run |
| R6 | Proposition 1(c) | done in rev2 | verified (theory check) |
| R7 / B3, B4, B7 | TOST calibration, synchronized seeds, equal-cluster wild bootstrap | G | |
| R8 / C2, C3 | re-evaluation vs re-optimization | G | direct 1° optimization; GA (and Lillgrund) layouts re-evaluated |
| R9 / B6 | charged evaluations vs elapsed time | G | original Seconds + initialization timing |
| R10 / E1–E4 | presentation | G | |
| C4 | full-precision coordinates | G | deterministic reruns with 17-digit writer |
| A3 | primary equivalence questions agreed | A | |
| A4 | AI declaration accuracy | A (wording updated to name all tools used) | |
| A5 | Solanki & Deep 2023 volume/pages | GP | |
| A6 | reference audit | GP | |
| A7 | SWEVO Guide for Authors | GP | |
| A8 | AUTHOR_CONFIRMATION items | A / partly resolved (archive and precision items) | |
| E5, E6 | highlights/cover letter; response letter | G | |

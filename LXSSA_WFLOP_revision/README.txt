LXSSA-WFLOP FOLLOW-UP REVIEWER REVISION
====================================

MAIN MANUSCRIPT SOURCE:
  LXSSA_WFLOP_reviewer_revised.tex

ADDED FOLLOW-UP ANALYSIS:
  new_analysis_fragment.tex
  selected_30_run_data.csv
  selected_descriptive_stats.csv
  selected_omnibus_stats.csv
  selected_pairwise_stats.csv
  selected_power_boxplots.pdf and .png
  REVIEWER_STATUS_AND_REMAINING.md

COMPILATION:
  Updated TeX requires IEEEtran.cls and algorithm.sty. LXSSA_WFLOP_reviewer_revised.pdf
  is compiled from the current source (September 2026 integration pass).
  new_analysis_fragment.tex is superseded: its content is inlined and extended in the main .tex.

ROBUSTNESS / SENSITIVITY ANALYSIS (new):
  analysis/  -- scripts, CSV outputs and README (evaluator validation, power-curve and
  Gaussian-wake re-evaluation, 4D/5D/6D spacing study, Friedman/bootstrap statistics).

The original tables and the DS1/DS2 workbooks are distinct experiments.
The follow-up uses 6030 calls for LX-SSA and 3030 for each baseline;
its significance tests are exploratory.

REVIEW TRACKING:
  Reviewer_Comment_Status.md
  AUTHOR_ACTION_REQUIRED.txt

IMPORTANT FRAMING:
  LX-SSA was introduced by Solanki and Deep in 2023.
  This manuscript EMPLOYS/APPLIES the existing LX-SSA to WFLOP; it does not propose or redesign LX-SSA.

The package intentionally excludes older duplicate manuscript sources (main.tex, Improved LXSSA_WFLOP.tex, MPCE_jrnl.tex, etc.) to prevent stale claims and formatting errors from re-entering the submission.

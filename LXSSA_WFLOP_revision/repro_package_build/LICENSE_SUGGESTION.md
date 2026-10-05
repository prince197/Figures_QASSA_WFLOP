# Licence — SUGGESTION, to be confirmed by the authors

**No licence has been chosen.** The authors decide on the licence before the public deposit. Until then, this
package is provided to the editor and reviewers for the review of the article only, and all rights are reserved by
the authors.

Suggested (common for research code and data, compatible with Zenodo and GitHub):

| Part | Suggested licence | Files |
|---|---|---|
| Code | MIT License | `analysis/*.py`, `*.sh`, the package scripts (`regenerate_all.sh`, `compare_outputs.py`, `verify_determinism.py`, `count_records.py`, `build_inventory.py`, `make_package.sh`) |
| Data and generated outputs | Creative Commons Attribution 4.0 International (CC-BY-4.0) | run records (`analysis/fresh_*.csv`, `analysis/mpce_*.csv`, `analysis/rev2_*.csv`, `analysis/rev3_*.csv`, `selected_30_run_data.csv`, ...), generated tables, number macros, summaries, figures |

Third-party material keeps its own terms and must be checked before the deposit:

- `analysis/iea37_data/`: files of the IEA Wind Task 37 case-study repository
  (https://github.com/IEAWindTask37/iea37-wflo-casestudies, folder `cs1-2`, commit 267f6e5), copied unmodified;
  the terms of that repository apply (check it before redistribution, or replace the folder by a download script).
- `analysis/hornsrev_model.py` and `analysis/rev2_site_model.py` contain wind-climate and turbine data taken from
  DTU PyWake 2.6.20 (MIT licence; Lillgrund wind climate after Göçmen and Giebel, 2016); keep the attribution in
  the file headers.
- `analysis/objective_original.py` and `analysis/authors_optimizers.py` are the original code of the earlier study;
  the authors confirm that they may be redistributed under the chosen licence.

If the suggestion is accepted: add `LICENSE` (MIT text, copyright holders = the four authors, year 2026) and
`LICENSE-DATA` (CC-BY-4.0 text or link), set `license:` in `CITATION.cff` (e.g. `license: MIT` and state the data
licence in the README), and delete this file.

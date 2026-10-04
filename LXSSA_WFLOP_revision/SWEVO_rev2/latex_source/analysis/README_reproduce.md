# Reproducing the results of the SWEVO manuscript

This folder (`latex_source/analysis/`) holds the copies of the generated tables and number macros that the LaTeX
sources input, the archive-reliability table and the document build script. It does not regenerate results. All
results are regenerated from the **reproducibility package**, which contains the complete archive.

## Where the code and data are

- **Package:** `repro_package_build/` of the project repository (README, `MANIFEST.sha256`, `requirements.txt`,
  `requirements-pywake.txt`, `Dockerfile`, `CITATION.cff`, `regenerate_all.sh`, `verify_determinism.py`,
  `count_records.py`, `inventory.csv`, verification logs in `verification/`). The release candidate archive is
  `repro_package_build/dist/SWEVO_repro_1.0.0-rc.zip` (+ `.sha256`).
- **Contents:** the per-run records of every optimization run with final coordinates and convergence curves
  (36,190 runs of the original study, 11,640 runs of revision 2 and the revision-3 studies `analysis/rev3_*`), all
  models and evaluators (benchmark, Horns Rev 1, IEA37, Lillgrund), all optimizers, the experiment drivers and the
  analysis and check scripts. The package README maps every table, figure and number to its script, command, input
  files and run time, and every run file to the command that produced it.
- **Public deposit:** pending author action. No DOI exists yet; the authors create the versioned public deposit and
  enter the DOI in `CITATION.cff` and in the Data availability statement (`optA/swevo_back.tex`,
  `declarations_content.tex`).

## Environment

CPython 3.11.15, NumPy 2.4.6, SciPy 1.17.1, pandas 3.0.6, Matplotlib 3.11.2 (the environment of the runs);
optional PyWake 2.6.20 for the Horns Rev 1 and Lillgrund cross-checks. Install with
`python -m pip install -r requirements.txt` in the unpacked package, or build the `Dockerfile`.

## Regenerate and verify

In the unpacked package:

```sh
sha256sum -c MANIFEST.sha256                       # integrity of every file
sh regenerate_all.sh /tmp/regen                    # all tables, macros, summaries, figures + all checks (~6 min)
WITH_SLOW=1 WITH_DETERMINISM=1 WITH_PYWAKE=1 MANUSCRIPT=/path/to/latex_source sh regenerate_all.sh /tmp/regen
```

`/tmp/regen/logs/compare.txt` lists each regenerated file as IDENTICAL, IDENTICAL-NORMALIZED (differs only in
generation time stamps, PDF creation dates or recorded versions) or DIFFERENT.

Verification of this revision (`repro_package_build/verification/verification_summary.txt`): 91 regenerated or
compared outputs, 0 different; all checks pass (CHECK-FINAL 62, CHECK-EXTRA 57, CHECK-DIAG 22, CHECK-DIR 41,
CHECK-CSWEEP 15, theory checks 57); `rev2_analysis.py` reproduces its summary (4,670 values) and tables;
`audit_archive.py` and `check_revision.py` pass; reruns of stored runs with the archived drivers are bit-identical to
the records in every field except the elapsed time. Not re-executed: the optimization studies themselves (about 126
CPU-hours; commands in Section 6 of the package README), the instrumented diagnostic reruns (stored outputs checked
by D01–D22) and the calibration reruns of the reused code.

## Revision-3 studies

Prespecified in `analysis/rev3_*_manifest.md`; drivers, records, analyses and tables in `analysis/rev3_*`:
constraint handling (`rev3_constraint*`, 2,700 runs), direct 1-degree optimization with a seed-paired 15-degree arm
(`rev3_fine*`, 1,200 + 1,200 runs), statistical audits (`rev3_inference*`), site pools, alternate evaluators and
timing (`rev3_sites*`), numbers audit (`rev3_numbers*`) and coordinate precision (`rev3_precision*`).

## Coordinate precision

Run files written before revision 3 store final coordinates rounded to 1 mm (benchmark, IEA37) or 1 cm (Horns Rev 1,
Lillgrund); objective, feasibility label and minimum spacing are stored at full precision from the unrounded layout.
The audit of all 56,540 stored records (`rev3_precision_audit.json`) confirms 50,826 labels from the rounded
coordinates, contradicts none and leaves 5,714 undecidable within the rounding bound. These are rerun
deterministically with the 17-digit writer (`record_io.encode_coordinates`; `rev3_precision_rerun.py`, side files
`rev3_fullprec_*.csv`); result: [TBD: C4 rerun — N bit-identical, labels confirmed/changed]. New runs store 17
significant digits.

## Build the documents

```sh
sh latex_source/analysis/build_swevo.sh
```

compiles the manuscript, the single-file manuscript, the supplement, the cover letter and the declarations and
refreshes the cross-document labels (TeX Live with `elsarticle`, `algorithms`, `algorithmicx`, `lmodern`, `lineno`,
`environ`, `placeins`). The main bibliography is a manual `thebibliography` in order of first citation; after
citations change, run, from `LXSSA_WFLOP_revision/`, `python3 analysis/rev3_bib_order.py
SWEVO_rev2/latex_source/SWEVO_manuscript.aux SWEVO_rev2/latex_source/optA/swevo_back.tex` and rebuild.

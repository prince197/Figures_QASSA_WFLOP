# Revised SWEVO package — 4 October 2026

This is an evidence-supported revision of the uploaded `SWEVO_revision2_full.zip`. It is **not yet a submission-ready or fully reproducible research archive**. Missing original data/code, legacy coordinate precision and author confirmations remain open.

Read `REVISION_CHANGELOG.md` for the disposition of the 12 audit areas and `AUTHOR_CONFIRMATION.md` before submitting. Five areas are closed at the stated scope; seven are partially addressed. The detailed implementation report is in `review/`.

## Documents

- `latex_source/SWEVO_manuscript.pdf`: revised main manuscript.
- `latex_source/SWEVO_manuscript_full.pdf`: identical scientific content compiled from generated single-file source.
- `latex_source/SWEVO_supplement.pdf`: revised supplement, including the executable archive audit and new reliability tables.
- `latex_source/SWEVO_cover_letter.pdf`: revised draft cover letter; complete archive limitations are disclosed.
- `latex_source/SWEVO_declarations.pdf`: consistent declaration wording carried from supplied author statements, with the present AI-assisted revision disclosed. Author verification is still required.
- `latex_source/SWEVO_highlights.txt`: five concise highlights.

## What changed

The fuller manuscript is now canonical in modular form; the shorter source no longer drops later experiments. The PSO zero-variance statement covers zero acceleration coefficients. Inference language distinguishes fixed-benchmark estimates from assumptions about exchangeable cases and group dependence. GA and analytic-gradient results appear in the main tables, with IEA37 conditional means and feasibility. New Lillgrund paired all-run estimates complement the historical conditional-quality rankings. Gradient, constraint handling and alternative-evaluator claims have narrower, accurate scope. Data availability and historical verification statements now disclose the supplied archive's limits.

All 140 original CSV files are unchanged. Existing results have not been fabricated, silently relabelled or presented as freshly rerun. Five experiment drivers preserve 17 significant digits in future coordinate outputs; this cannot restore precision lost in legacy records.

## Run the supported checks

From this folder, with Python 3.12 and an available package installer:

```sh
python3 -m pip install -r requirements-audit.txt
python3 experiments/audit_archive.py
python3 experiments/check_revision.py
```

The first command regenerates derived records in `validation/`; the second checks binary float round-trips, statistical reference values, the PSO correction and raw CSV hashes. Regenerating the audit does not regenerate the original research results or optimize new layouts.

To include regenerated reliability tables in the supplement, copy `validation/reliability_tables.tex` to `latex_source/analysis/archive_reliability_tables.tex` before building.

## Build the documents

Use a TeX Live installation with `elsarticle`, `algorithms`, `algorithmicx`, `lmodern`, `lineno`, `environ`, `placeins` and the standard AMS/graphics packages:

```sh
sh latex_source/analysis/build_swevo.sh
```

No BibTeX step is needed: references are in the supplied `thebibliography`. The build refreshes cross-references and generates the full manuscript from the modular source. Edit `latex_source/optA/swevo_front.tex`, `optA/sw/*.tex`, `optA/swevo_back.tex` or `manuscript_preamble.tex`, then rebuild. **Do not edit the generated full file as a separate version.** Included frozen label files support single-document compilation but should be refreshed after structural changes.

Original optimization launchers and `rev2_analysis.py` remain blocked by missing project modules. Do not treat the historical reproducibility map as a list of runnable included scripts; see `latex_source/analysis/README_reproduce.md`.

#!/bin/sh
# Regenerates the results part of MPCE_PSO_VNS.tex (all mpce_tab_*.tex, mpce_supplementary.tex, figures_mpce/,
# mpce_summary.json, mpce_numbers.tex, CHECK-FINAL list) from the per-run CSVs, then compiles both PDFs.
# NOT regenerated here (separate scripts, see README_reproduce.md section 3): model schematics
# (make_model_figures.py), spacing re-optimization table (analyze_authors_runs.py -> authors_tables.tex),
# packing bounds (packing_capacity.py), evaluator checks (validate_evaluator.py, calibrate_authors_code.py),
# literature table (literature_values.csv), and the per-run data themselves (mpce_experiments.py etc.).
# Python stages only (no LaTeX):  cd analysis && python3 mpce_results.py
# Usage once new shards have arrived:
#   sh analysis/build_mpce_paper.sh            (all experiments complete)
#   sh analysis/build_mpce_paper.sh --partial  (preview with incomplete shards)
# Steps: mpce_results.py (tables, figures, mpce_summary.json; it also runs mpce_numbers.py -> mpce_numbers.tex
# and mpce_check_final.py -> PASS/FAIL of the CHECK-FINAL statements), then pdflatex of the supplement and the
# main text in alternation (they cross-reference each other through xr).
set -e
HERE="$(cd "$(dirname "$0")" && pwd)"
cd "$HERE"
python3 mpce_results.py "$@"
cd "$HERE/.."
for i in 1 2 3; do
  pdflatex -interaction=nonstopmode -halt-on-error MPCE_PSO_VNS_supplement.tex > /dev/null || { echo "supplement: LaTeX error, see MPCE_PSO_VNS_supplement.log"; exit 1; }
  pdflatex -interaction=nonstopmode -halt-on-error MPCE_PSO_VNS.tex > /dev/null || { echo "main text: LaTeX error, see MPCE_PSO_VNS.log"; exit 1; }
done
for f in MPCE_PSO_VNS MPCE_PSO_VNS_supplement; do
  echo "$f.pdf: $(grep -o 'Output written.*' $f.log | sed 's/.*(\([0-9]* pages\).*/\1/')," \
       "undefined refs: $(grep -c 'undefined' $f.log || true), overfull boxes: $(grep -c 'Overfull' $f.log || true)"
done
echo "pending: $(grep -c 'TBD{pending' analysis/mpce_numbers.tex) macros in analysis/mpce_numbers.tex (listed in its 2nd line)," \
     "$(cat analysis/mpce_tab_*.tex | grep -o 'TBD{}' | wc -l) empty table cells; author placeholders in MPCE_PSO_VNS.tex: $(cat MPCE_PSO_VNS.tex optA/*.tex | grep -o 'TBD{[^p}]' | wc -l)"

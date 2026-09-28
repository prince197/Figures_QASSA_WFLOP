#!/bin/sh
# One command to regenerate everything for MPCE_PSO_VNS.tex once new shards have arrived:
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
echo "TBD markers printed in the main text: $(grep -o 'TBD' MPCE_PSO_VNS.tex analysis/mpce_numbers.tex analysis/mpce_tab_*.tex | wc -l) in the sources (pending macros listed at the top of analysis/mpce_numbers.tex)"

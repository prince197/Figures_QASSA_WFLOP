#!/bin/sh
# Compiles the Swarm and Evolutionary Computation (SWEVO) package from the current number macros and tables:
# SWEVO_supplement.tex and SWEVO_manuscript.tex alternately (three rounds; they cross-reference each other
# through xr), then SWEVO_cover_letter.tex and SWEVO_declarations.tex if present (twice each).
# The Python stages of the pipeline are NOT rerun here; to regenerate tables, figures and \N... macros first,
# run  sh analysis/build_mpce_paper.sh  (it also rebuilds the MPCE version) or the scripts listed in
# analysis/README_reproduce.md. Brief: optA/SWEVO.md.
# Usage (from anywhere):  sh analysis/build_swevo.sh
set -e
HERE="$(cd "$(dirname "$0")" && pwd)"
cd "$HERE/.."
for f in SWEVO_supplement SWEVO_manuscript; do
  [ -f "$f.tex" ] || { echo "$f.tex not found in $(pwd)"; exit 1; }
done
for i in 1 2 3; do
  pdflatex -interaction=nonstopmode -halt-on-error SWEVO_supplement.tex > /dev/null || { echo "supplement: LaTeX error, see SWEVO_supplement.log"; grep -A5 '^!' SWEVO_supplement.log | head -20; exit 1; }
  pdflatex -interaction=nonstopmode -halt-on-error SWEVO_manuscript.tex > /dev/null || { echo "manuscript: LaTeX error, see SWEVO_manuscript.log"; grep -A5 '^!' SWEVO_manuscript.log | head -20; exit 1; }
done
DOCS="SWEVO_manuscript SWEVO_supplement"
for f in SWEVO_cover_letter SWEVO_declarations; do
  if [ -f "$f.tex" ]; then
    for i in 1 2; do
      pdflatex -interaction=nonstopmode -halt-on-error "$f.tex" > /dev/null || { echo "$f: LaTeX error, see $f.log"; grep -A5 '^!' "$f.log" | head -20; exit 1; }
    done
    DOCS="$DOCS $f"
  else
    echo "$f.tex not present (skipped)"
  fi
done
for f in $DOCS; do
  echo "$f.pdf: $(grep -o 'Output written.*' $f.log | sed 's/.*(\([0-9]* pages*\).*/\1/')," \
       "undefined refs: $(grep -c 'undefined' $f.log || true), overfull boxes: $(grep -c 'Overfull' $f.log || true)"
done
echo "pending: $(grep -c 'TBD{pending' analysis/mpce_numbers.tex) macros in analysis/mpce_numbers.tex (listed in its 2nd line);" \
     "author placeholders in the SWEVO manuscript files: $(cat SWEVO_manuscript.tex optA/swevo_front.tex optA/sw/*.tex $(ls optA/swevo_back.tex 2>/dev/null) | grep -o 'TBD{[^p}]' | wc -l)"

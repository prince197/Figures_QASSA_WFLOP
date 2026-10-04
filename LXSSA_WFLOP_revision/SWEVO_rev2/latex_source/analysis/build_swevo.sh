#!/bin/sh
# Build documents only (the results themselves are regenerated from the reproducibility package; see analysis/README_reproduce.md).
set -eu
cd "$(dirname "$0")/.."
command -v pdflatex >/dev/null || { echo 'Install TeX Live with elsarticle, algorithms and algorithmicx.' >&2; exit 1; }
for pass in 1 2 3; do
  pdflatex -interaction=nonstopmode -halt-on-error SWEVO_supplement.tex > build_supplement.log 2>&1
  python3 analysis/refresh_references.py
  pdflatex -interaction=nonstopmode -halt-on-error SWEVO_manuscript.tex > build_manuscript.log 2>&1
  python3 analysis/refresh_references.py
done
python3 analysis/flatten_manuscript.py
pdflatex -interaction=nonstopmode -halt-on-error SWEVO_manuscript_full.tex > build_full.log 2>&1
pdflatex -interaction=nonstopmode -halt-on-error SWEVO_manuscript_full.tex >> build_full.log 2>&1
pdflatex -interaction=nonstopmode -halt-on-error SWEVO_manuscript_full.tex >> build_full.log 2>&1
for doc in SWEVO_cover_letter SWEVO_declarations SWEVO_response_to_reviewers; do
  pdflatex -interaction=nonstopmode -halt-on-error "$doc.tex" > "build_${doc}.log" 2>&1
  pdflatex -interaction=nonstopmode -halt-on-error "$doc.tex" >> "build_${doc}.log" 2>&1
done
python3 - <<'PY'
from pathlib import Path
for name in ('SWEVO_manuscript', 'SWEVO_manuscript_full', 'SWEVO_supplement', 'SWEVO_cover_letter', 'SWEVO_declarations', 'SWEVO_response_to_reviewers'):
    path = Path(name + '.pdf')
    if not path.read_bytes().rstrip().endswith(b'%%EOF'):
        raise RuntimeError(f'Incomplete PDF: {path}; inspect its build log and recompile.')
PY
echo 'Built manuscript, full manuscript, supplement, cover letter and declarations PDFs.'

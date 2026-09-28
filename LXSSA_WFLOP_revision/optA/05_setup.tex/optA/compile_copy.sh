#!/bin/sh
# Compile the manuscript in a private copy so parallel agents do not clash on .aux/.log files.
# Usage: sh optA/compile_copy.sh <scratch_dir>   (run from LXSSA_WFLOP_revision/)
set -e
D="$1"; [ -n "$D" ] || { echo "usage: sh optA/compile_copy.sh <scratch_dir>"; exit 1; }
mkdir -p "$D/analysis"
cp MPCE_PSO_VNS.tex MPCE_PSO_VNS_supplement.tex MPCE_PSO_VNS_supplement.aux "$D/" 2>/dev/null || true
cp -r optA figures_mpce "$D/"
[ -d figures_final ] && cp -r figures_final "$D/"
cp analysis/*.tex "$D/analysis/"
cd "$D"
for i in 1 2; do pdflatex -interaction=nonstopmode -halt-on-error MPCE_PSO_VNS.tex > build.log 2>&1 || { echo "LaTeX ERROR:"; grep -A5 '^!' MPCE_PSO_VNS.log | head -20; exit 1; }; done
python3 - <<'PY'
import re
log = open("MPCE_PSO_VNS.log", errors="ignore").read()
print("pages:", re.findall(r"Output written on .*?\((\d+) pages", log))
print("overfull boxes:", log.count("Overfull"), " undefined refs:", len(re.findall(r"Reference .* undefined", log)))
PY

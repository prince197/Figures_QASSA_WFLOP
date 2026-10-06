#!/bin/sh
# Regenerate every table, figure, number macro, summary and check of the SWEVO paper from the stored run records,
# in a scratch copy, and compare the results with the stored copies of this package.
#
# Usage (from the unpacked package folder):   sh regenerate_all.sh [WORKDIR]
#   WORKDIR   scratch folder (default: ./regen_work); it is DELETED and recreated. The package files are never
#             written: everything runs in WORKDIR, a copy of analysis/, figures_mpce/, figures_final/ and the inputs.
# Environment variables:
#   PROCS=1            worker processes for mpce_results.py / mpce_direction.py (default 1)
#   WITH_SLOW=1        also run the slower optional steps: theory Monte Carlo check (~3-6 min), packing bounds,
#                      evaluator validation, IEA37 calculator check, authors-run tables, model schematics,
#                      graphical abstract,
#                      inference sensitivity analyses with null calibrations (rev3_inference.py, ~50 min)
#   WITH_PYWAKE=1      also run the PyWake cross-checks (needs requirements-pywake.txt)
#   WITH_DETERMINISM=1 also rerun seed 1 of two benchmark cases for PSO-VNS and SSA-VNS (verify_determinism.py, ~30 s)
#   MANUSCRIPT=DIR     a folder with the SWEVO LaTeX sources (SWEVO_manuscript.tex, SWEVO_supplement.tex, optA/sw/*.tex);
#                      enables the CHECK-tag scan of mpce_check_final.py and make_theory_figures.py --check --swevo
# Output: WORKDIR/logs/<step>.log, WORKDIR/logs/steps.txt (exit code and seconds per step) and
#         WORKDIR/logs/compare.txt (IDENTICAL / IDENTICAL-NORMALIZED / DIFFERENT per regenerated file).
# Exit code: 0 if every step succeeded and every compared file is identical (after removing time stamps only).
set -u
PKG=$(cd "$(dirname "$0")" && pwd)
WORK=${1:-$PKG/regen_work}
PROCS=${PROCS:-1}
case "$WORK" in /*) ;; *) WORK="$(pwd)/$WORK" ;; esac
[ "$WORK" = "$PKG" ] && { echo "WORKDIR must not be the package folder"; exit 2; }
rm -rf "$WORK"; mkdir -p "$WORK/logs"
cp -R "$PKG/analysis" "$PKG/figures_mpce" "$PKG/figures_final" "$PKG/selected_30_run_data.csv" "$WORK"/
mkdir -p "$WORK/SWEVO_rev2" && cp -R "$PKG/SWEVO_rev2/validation" "$PKG/SWEVO_rev2/experiments_docs" "$WORK/SWEVO_rev2/"   # reference files read by some analyses
LOG="$WORK/logs"; FAIL=0
step() { name=$1; shift; t0=$(date +%s); ( "$@" ) > "$LOG/$name.log" 2>&1; rc=$?
         [ $rc -eq 0 ] || FAIL=1; echo "$name: exit $rc, $(( $(date +%s) - t0 )) s" | tee -a "$LOG/steps.txt"; }
export MPLBACKEND=Agg OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
A="$WORK/analysis"
cd "$A"
python3 -c "import sys, numpy, scipy, pandas, matplotlib; print('python', sys.version.split()[0], 'numpy', numpy.__version__, 'scipy', scipy.__version__, 'pandas', pandas.__version__, 'matplotlib', matplotlib.__version__)" | tee "$LOG/environment.txt"

# 0. record counts (36,190 main study + 11,640 additional experiments + sensitivity-study files)
step count_records python3 "$PKG/count_records.py" --analysis-dir "$A"
# 1. main pipeline: tables, figures, mpce_summary.json, mpce_numbers.tex, CHECK-FINAL C01-C62   (~1 min)
TEXARG=""
if [ -n "${MANUSCRIPT:-}" ]; then mkdir -p "$WORK/ms"; cp -R "$MANUSCRIPT"/. "$WORK/ms/"; TEXARG="$WORK/ms/SWEVO_manuscript.tex"; fi
step mpce_results python3 mpce_results.py --procs "$PROCS"
[ -n "$TEXARG" ] && step mpce_check_final_tags python3 mpce_check_final.py --tex "$TEXARG"
# 2. inference robustness, equivalence at three levels (mpce_numbers_extra.tex, mpce_supp_inference.tex)   (~30 s)
step mpce_inference_extra python3 mpce_inference_extra.py
# 3. direction resolution, Horns Rev 1 re-evaluation, projected IEA37 layouts   (~2.5 min on 1 process)
step mpce_direction python3 mpce_direction.py --procs "$PROCS"
# 4. PSO coefficient sweep analysis   (~5 s)
step mpce_csweep python3 mpce_csweep.py
# 5. checks on the regenerated summaries (X01-X57, D01-D22 on the stored diagnostics summary, F01-F46, S01-S15)
step mpce_check_extra python3 mpce_check_extra.py
step mpce_check_diag python3 mpce_check_diag.py
step mpce_check_dir python3 mpce_check_dir.py
step mpce_check_csweep python3 mpce_check_csweep.py
# 6. statistics and tables of the additional experiments   (~15 s)
step rev2_analysis python3 rev2_analysis.py --out-dir "$WORK/rev2_out"
# 6b. sensitivity studies and analyses (constraint handling, direct 1-degree optimization, sites/budgets with their
#     stored caches, precision audit and collection of the full-precision reruns)   (~1.5 min)
mkdir -p "$WORK/rev3_out"
step rev3_constraint python3 rev3_constraint_analysis.py --out-dir "$WORK/rev3_out"
step rev3_fine python3 rev3_fine_analysis.py --out-dir "$WORK/rev3_out"
step rev3_sites python3 rev3_sites.py --out-dir "$WORK/rev3_out"
step rev3_precision_audit python3 rev3_precision_audit.py --out-dir "$WORK/rev3_out"
step rev3_precision_collect python3 rev3_precision_rerun.py collect      # rewrites rev3_fullprec_<study>.csv in WORKDIR
step rev3_precision_table python3 rev3_precision_rerun.py table
step rev3_feasbounds python3 rev3_feasbounds.py --skip-transfer --out-dir "$WORK/rev3_out"   # bounds of undecided labels (transfer block reused)
# 7. archive audit and focused checks of the additional experiments, in the layout they expect (experiments/ + validation/)   (~5 s)
mkdir -p "$WORK/rev2_audit/experiments" "$WORK/rev2_audit/validation"
cp "$A"/rev2_*_s*of*.csv "$A"/audit_archive.py "$A"/check_revision.py "$A"/record_io.py "$A"/rev2_site_model.py "$WORK/rev2_audit/experiments/"
cp "$PKG"/SWEVO_rev2/validation/* "$WORK/rev2_audit/validation/"
cd "$WORK/rev2_audit"
step audit_archive python3 experiments/audit_archive.py
step check_revision python3 experiments/check_revision.py
cd "$A"
# 8. optional slower steps
if [ "${WITH_SLOW:-0}" = "1" ]; then
  if [ -n "$TEXARG" ]; then   # figure + numbers + checks T01-T13 against the manuscript text (~4 min)
    mkdir -p "$WORK/optA/sw"; cp "$WORK/ms/optA/sw/supp_theory.tex" "$WORK/optA/supp_theory.tex"; cp "$WORK/ms/optA/sw/04_methods.tex" "$WORK/optA/sw/"
    step theory_check python3 make_theory_figures.py --check --swevo
  else                         # figure + numbers only (~4 min)
    step theory_figure python3 make_theory_figures.py
  fi
  step validate_evaluator python3 validate_evaluator.py
  step iea37_calculator python3 iea37_model.py
  step packing_capacity python3 packing_capacity.py
  step analyze_authors_runs python3 analyze_authors_runs.py
  step model_figures python3 make_model_figures.py
  step graphical_abstract python3 make_graphical_abstract.py
  mkdir -p "$WORK/rev3_slow"     # inference sensitivity analyses with the null calibrations (~50 min)
  step rev3_inference python3 rev3_inference.py --out-dir "$WORK/rev3_slow"
  step rev3_dependence python3 rev3_dependence.py --out-dir "$WORK/rev3_slow"   # seed-level all-run inference + joint calibration (~12 min)
fi
if [ "${WITH_PYWAKE:-0}" = "1" ]; then
  mkdir -p "$WORK/pywake_out"
  step pywake_check python3 pywake_check.py --layouts --out-dir "$WORK/pywake_out"
  step lillgrund_pywake python3 rev2_site_model.py --pywake
fi
if [ "${WITH_PYWAKE:-0}" = "1" ]; then mkdir -p "$WORK/rev3_pw"; step rev3_feasbounds_transfer python3 rev3_feasbounds.py --out-dir "$WORK/rev3_pw"; fi
[ "${WITH_DETERMINISM:-0}" = "1" ] && step determinism python3 "$PKG/verify_determinism.py" --analysis-dir "$A"

# 9. compare with the stored copies
C="$LOG/compare.txt"; : > "$C"
cmpf() { python3 "$PKG/compare_outputs.py" "$@" >> "$C" 2>&1 || FAIL=1; }
cmpf "$PKG/analysis" "$A" 'mpce_numbers*.tex' 'mpce_tab_*.tex' 'mpce_supp*.tex' 'mpce_summary*.json' mpce_case_stats.csv \
     mpce_case_tests.csv mpce_ablation_tests.csv mpce_feasinit_tests.csv mpce_budget_case_stats.csv mpce_best_layouts_maxN.csv \
     mpce_reevaluation.csv mpce_reevaluation_cache.csv mpce_direction_layouts.csv.gz mpce_direction_hr16.csv \
     iea37_projected.csv iea37_projected.json
cmpf "$PKG/figures_mpce" "$WORK/figures_mpce" '*.pdf' '*.png'
cmpf "$PKG/SWEVO_rev2/experiments_docs" "$WORK/rev2_out" rev2_summary.json rev2_tables.tex
cmpf "$PKG/SWEVO_rev2/validation" "$WORK/rev2_audit/validation" archive_audit.json focused_checks.json paired_outcomes.csv \
     raw_data_checksums.json reliability.csv reliability_tables.tex
if [ "${WITH_SLOW:-0}" = "1" ]; then
  cmpf "$PKG/analysis" "$A" packing_capacity.csv authors_tables.tex calibration_vs_recorded.csv equal_budget_tests.csv \
       hornsrev_tests.csv authors_spacing_all.csv authors_runtime_6030.csv
  cmpf "$PKG/figures_final" "$WORK/figures_final" fig_wind_farm.pdf fig_wake_model.pdf fig_half_cone.pdf graphical_abstract.pdf graphical_abstract.png graphical_abstract.tif
fi
cmpf "$PKG/analysis" "$WORK/rev3_out" rev3_constraint.json rev3_constraint_tables.tex rev3_fine.json rev3_fine_tables.tex \
     rev3_sites.json rev3_sites_tables.tex rev3_precision_audit.json rev3_precision_audit_records.csv \
     rev3_feasbounds.json rev3_feasbounds_tables.tex
cmpf "$PKG/analysis" "$A" rev3_precision_rerun_summary.json rev3_precision_tables.tex 'rev3_fullprec_*.csv'
[ "${WITH_SLOW:-0}" = "1" ] && cmpf "$PKG/analysis" "$WORK/rev3_slow" rev3_inference.json rev3_inference_tables.tex rev3_dependence.json rev3_dependence_tables.tex
[ "${WITH_PYWAKE:-0}" = "1" ] && cmpf "$PKG/analysis" "$WORK/rev3_pw" rev3_feasbounds.json rev3_feasbounds_tables.tex
[ "${WITH_PYWAKE:-0}" = "1" ] && cmpf "$PKG/analysis" "$WORK/pywake_out" pywake_check.csv pywake_check_hr16runs.csv pywake_check_hr16runs.json
grep -h "files identical" "$C"
echo "steps: $LOG/steps.txt; comparison: $C"
grep -h "^DIFFERENT\|^MISSING" "$C" || echo "no DIFFERENT / MISSING files"
exit $FAIL

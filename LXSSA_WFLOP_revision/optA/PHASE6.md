# Phase 6 — strengthening round (all agents read this, plus optA/BRIEF.md and optA/PHASE4.md)

Goal: make every claim scientifically, theoretically, numerically and experimentally robust so that the
reviews (optA/reviews/R1–R8) leave no open objection. NO new optimization algorithms. Controls, diagnostics,
instrumentation, statistics and theory are allowed.

Decided by the user: practical-equivalence margin ±0.05 pp of wake loss (EQ_MARGIN), stated as fixed after
the primary analysis, with the minimal equivalence margin reported so readers can apply their own.

## Workstreams (Phase A, parallel) — each owns ONLY the files named; nothing else
- W1 pipeline equivalence (analysis/mpce_results.py, mpce_numbers.py, mpce_check_final.py): TOST, Bayesian
  signed-rank with ROPE, spread/worst-run, tab:equivalence.  ONLY W1 edits these three files.
- W2 inference robustness: new analysis/mpce_inference_extra.py → analysis/mpce_numbers_extra.tex (macros
  prefixed \NX...), analysis/mpce_supp_inference.tex (supp tables), analysis/mpce_summary_extra.json,
  analysis/mpce_check_extra.py (checks X01…).
- W3 diagnostics & controls: new analysis/mpce_diagnostics.py; new experiments appended to
  analysis/mpce_experiments.py (only new `if exp == ...` blocks + helpers; never change existing behaviour);
  outputs analysis/mpce_numbers_diag.tex (\ND...), analysis/mpce_supp_diag.tex, figures_mpce/diag_*.pdf,
  analysis/mpce_summary_diag.json, analysis/mpce_check_diag.py (checks D01…).
- W4 theory: optA/drafts/theory_main.tex (≤ 0.5 page for the main text) and optA/supp_theory.tex (supplement
  section, \input by the lead), figures_mpce/theory_*.pdf via analysis/make_theory_figures.py.
- W5 literature & positioning: optA/drafts/positioning.md + optA/newrefs/phase6.tex (verified \bibitem only).

Phase B (lead): integrate drafts into optA/*.tex and the supplement, adversarial re-review, final build.

## Rules
- Numbers only via macros generated from data (never typed by hand in .tex). Every new claim gets a check.
- Do not run pdflatex in the shared folder; to test LaTeX, copy what you need into your private scratch dir.
- Do not commit or push. Use private scratch dirs under
  /tmp/claude-0/-home-user-Figures-QASSA-WFLOP/83ceea2c-3e62-519f-bca6-caf988ef6cfa/scratchpad/p6_<you>/ .
- Local CPU is 4 cores shared by several agents: use at most 2 processes each for local runs; anything
  larger than ~20 CPU-minutes must be prepared as an experiment in mpce_experiments.py and reported to the
  lead, who runs it on cloud workers.
- Report: files created, macros (with values), checks (PASS/FAIL), and what the lead must integrate where.

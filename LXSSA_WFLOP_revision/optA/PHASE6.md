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

## W1 results (equivalence, ±0.05 pp) — notes for Phase B text
- PSO-VNS vs PSO: ΔL −0.018, 90% CI [−0.041, 0.004], TOST p 0.011, minimal margin 0.041, EQUIVALENT; Bayes P(rope) 0.997 (write "> 0.99"). C49.
- LX-SSA-VNS vs RS-VNS: EQUIVALENT (min margin 0.035; P(rope) 0.99). C50.
- SSA-VNS vs RS-VNS: NOT shown equivalent (90% CI [−0.101, −0.039]; min margin 0.102); Bayes 0.30 better / 0.70 rope → write "a small gain that could not be shown to be equivalent to zero at ±0.05 pp; the Bayesian analysis finds it more likely practically negligible (0.70)". Never "exceeds the margin". C51.
- SSA-VNS vs LX-SSA-VNS: not equivalent; Bayes 0.53 better / 0.47 rope.
- Spread: no significant difference (SD p 0.55, worst p 0.80) → do NOT claim PSO-VNS is more robust in spread; the robustness advantage is feasibility (Horns Rev 30/30 vs 19/30). C52.
- Supplement: add \SuppTableOpt{tab:equivalence} after tab:switch; add Benavoli2017 bibitem.

## W2 results (inference robustness) — notes for Phase B text
- Threshold 1–30: all 7 conclusions and the full rank order unchanged (X02–X09).
- Clusters (6): all 13 significant pairs favour the same method in 6/6 clusters; LOCO changes no verdict, PSO-VNS always best.
- CAVEAT: SSA-VNS vs RS-VNS: 6/6 clusters and cluster-bootstrap CI excludes 0, but cluster-robust t p = 0.061 → "consistent but small".
- "Gain grows with N" NOT supported overall; only within Data Set II (ρ = −0.66, p = 3e-4). Always write "in Data Set II".
- Multiplicity (22 tests, one Holm): only the pooled N ≥ 10 subgroup is lost; DS II N ≥ 10 survives (p_Holm 0.012).
- Energy: PSO-VNS vs PSO ≈ 0.020 % AEP (CI contains 0); DS II N ≥ 10 ≈ 0.22 %; SSA vs RS 0.071 %; PSO-VNS vs RS 0.41 % (benchmark energy).
- Integration: \input{analysis/mpce_numbers_extra.tex} in both preambles; \input{analysis/mpce_supp_inference.tex} in supplement; build: add mpce_inference_extra.py + mpce_check_extra.py after mpce_results.py.

## W4 results (theory) — notes for Phase B
- optA/drafts/theory_main.tex: Block A → 04 sec:pso (prop:pso-stability; Poli bound "for c1 = c2"; order-1 uses E[φ]); Block B → 04 sec:template (prop:box, (π/4)^N); Block C → 08 feasible starts (prop:feasible-start; stall is OBSERVED, not implied). ~+0.3 page net → cut elsewhere.
- optA/supp_theory.tex → \input in supplement (sec:S-theory); replace its two \TBD with W3 section refs; cite Benavoli2017 in sec:S-th-equiv.
- prop:S-prefix: PSO-VNS final value never worse than PSO's value at B1 (exact prefix property) — useful for the equivalence discussion.
- checks T01–T13: python analysis/make_theory_figures.py --check (3 min) → add to reproducibility map (not to every build).

## W3 results (diagnostics) — notes for Phase B (CONTRADICTIONS to fix in text)
- At 6,030 evaluations the "VNS phase" is essentially ONE compass-search descent (shaking 0.02 % of evals, 0.5 % of gain; no shaking for N ≥ 12; k ≥ 3 never tried). Text must say so honestly: the component analysis measures the value of a local-search (compass-search) phase; the VNS neighbourhood structure matters only for small N. (D08–D10)
- Old PSO setting: swarm spread 8.9×, |V| 13× the constriction setting; clipped 8.7 % vs 1.6 %; feasible particle evals 2.6 % vs 48.4 %. BUT stored infeasible final layouts are mostly SPACING-only (246/282): boundary penalty is in m² vs spacing in m → the best infeasible layout lies inside the circle with turbines too close. Replace "corners explain" in 04/06. (D03–D07)
- Feasible starts: label symmetry does NOT explain the stall (relabelling: 0.1 % feasible). Cause: independent per-coordinate weights with overshoot + spacing violations between distinct dense layouts; G-holder frozen (V = 0); no pbest ever updates; DE 0 of 126,000 trials feasible. Replace 08 wording. (D13–D16)
- RS-VNS: 730 of 2,040 runs have no feasible sample at all; replay agrees with (π/4)^N. (D17)
- Cloud experiments: omega90 (4 shards), rsdisc (8 shards). Pipeline must learn PSOBV90 and RSDVNS.

## Phase B (integration) — rules for section agents
- Read ALL notes above (W1–W4 results, including the CONTRADICTIONS) and the drafts: optA/drafts/theory_main.tex,
  optA/drafts/positioning.md, optA/newrefs/phase6.tex; data macros: analysis/mpce_numbers.tex (\N…),
  analysis/mpce_numbers_extra.tex (\NX…), analysis/mpce_numbers_diag.tex (\ND…); pending omega90/rsdisc macros
  are in optA/phase6_macros_stub.tex (names final; values arrive later).
- Each agent owns only its files. Page budget: the main PDF is exactly 12 pages now and must stay ≤ 12 —
  every file must end NO LONGER than now (04 may grow ≤ 0.15 page if 06/08 shrink); move details to the
  supplement (the supplement agent adds the new supplementary sections; point to them with \ref{S-...}).
- Comment conventions: keep `% CHECK-FINAL [Cnn]`; new claims from W2 get `% CHECK-EXTRA [Xnn]`, from W3
  `% CHECK-DIAG [Dnn]`, from W4 `% CHECK-THEORY [Tnn]`; equivalence claims `% CHECK-FINAL [C49]–[C52]`.
- New refs: only keys in optA/newrefs/*.tex or 11_back.tex (Benavoli2017, LaTorre2021, Derrac2011, Hooker1995,
  Campelo2019, Lakens2017, CamachoVillalon2023, Bonyadi2017, Thomas2023, Wilson2018, Azlan2021, Feng2015).
- Title (lead decision): "Baseline Configuration and Controls in Metaheuristic Comparisons for Wind Farm Layout
  Optimization, with a PSO–VNS Reference Method"; running head: "SOLANKI et al.: BASELINES AND CONTROLS IN
  METAHEURISTIC COMPARISONS FOR WIND FARM LAYOUT OPTIMIZATION".
- Headline M3 now: "a well-configured PSO and PSO-VNS are practically equivalent on average (±\NEqMargin pp);
  PSO-VNS adds feasibility reliability (Horns Rev \NHRFeasPSOVNS{} vs \NHRFeasPSO{} of 30) and gains in the
  larger Data Set II layouts; at 6,030 evaluations its second phase acts as one compass-search descent".

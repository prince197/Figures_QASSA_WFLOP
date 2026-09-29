# SWEVO expansion (main manuscript up to 40 pages incl. references) — brief for all agents

The SWEVO version now uses its own section copies optA/sw/*.tex (the MPCE files optA/0*.tex, 10_*.tex stay
unchanged — never edit them). Current SWEVO_manuscript.pdf: 23 pages (elsarticle preprint 12pt, a4, 2.5 cm margins,
line numbers). Target: 37–40 pages INCLUDING references (~5 pp) and back matter (~1.5 pp). Each agent has a page
budget for ADDED material; stay within it (measure with a private compile).

Goal: a complete, self-contained journal article for Swarm and Evolutionary Computation readers — bring the most
important evidence back from the supplement into the main text and deepen the explanations; no new claims beyond the
data; no new experiments; no new algorithms.

Tools available in the SWEVO main file:
- \SuppTable{<label>} places a generated supplementary table (captured from analysis/mpce_supplementary.tex) in the
  main text; standalone generated files can be \input directly (analysis/mpce_tab_split.tex,
  analysis/mpce_tab_robust_final.tex, analysis/mpce_supp_csweep.tex, analysis/mpce_supp_direction.tex,
  analysis/mpce_supp_inference.tex, analysis/mpce_supp_diag.tex — the latter four contain several floats; if you need
  only one of them, copy nothing: instead use \CaptureSuppTables{<file>} in your section file and then \SuppTable{label}
  for the one you want (the capture is global; do it once near the top of your file)).
- Figures: figures_mpce/*.pdf (avg_ranks, wakeloss_vs_n, boxplots_max/mid, convergence_max/mid, feasibility_vs_n,
  budget_scaling, hr16, iea37_convergence, iea37_layouts, layouts_iea37_hr16, layouts_max, ablation_convergence,
  csweep, equiv_curve, diag_pso_dynamics, theory_stability) and figures_final/*.pdf (fig_wind_farm, fig_wake_model,
  fig_half_cone). Use \pipefig{<width>}{<file>}{<height>} or \includegraphics; widths in \textwidth fractions.
- Theorem environments: proposition, lemma, definition, remark, proof (amsthm).
- Text sources to adapt: MPCE_PSO_VNS_supplement.tex / SWEVO_supplement.tex and optA/supp_theory.tex (full proofs,
  implementation details, algorithm boxes), optA/reviews2/*.md (what reviewers wanted explained).

Rules
- Numbers only via existing macros (analysis/mpce_numbers*.tex); keep and add CHECK tags with valid ids.
- Every label you bring into the main text that also exists in the supplement: append a line to optA/sw/moved.md
  ("<label> -> main, <your file>"); the supplement agent will remove the duplicate and re-point references.
  Do NOT \ref{S-<label>} for something you placed in the main text; use \ref{<label>}.
- New citations: only keys already in optA/swevo_back.tex, except the related-work agent (see its task).
- Compile privately: copy SWEVO_manuscript.tex, SWEVO_supplement.tex, SWEVO_supplement.aux, optA/, analysis/,
  figures_mpce/, figures_final/ into /tmp/claude-0/-home-user-Figures-QASSA-WFLOP/83ceea2c-3e62-519f-bca6-caf988ef6cfa/scratchpad/sx_<you>/
  and run pdflatex twice; 0 errors, 0 overfull boxes, no "??" from your section (S- refs to moved items will be fixed
  by the supplement agent). No commit/push. Report: what you added, pages added, labels moved, CHECK ids.

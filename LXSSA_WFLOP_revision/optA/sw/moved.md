fig:S-csweep -> main, optA/sw/06_results.tex (label unchanged; the whole generated file analysis/mpce_supp_csweep.tex is \input in the main text, so the supplement must no longer \input it; refer to it as M-fig:S-csweep)
tab:S-csweep -> main, optA/sw/06_results.tex (label unchanged; same \input of analysis/mpce_supp_csweep.tex)
fig:D-psodyn -> main, optA/sw/06_results.tex (label unchanged; captured from analysis/mpce_supp_diag.tex and placed as figure*; tab:D-psodyn, tab:D-vns, tab:D-feasstart, tab:D-rsreplay stay in the supplement: the supplement's \input of mpce_supp_diag.tex must skip only this figure)
fig:S-ranks -> main, optA/sw/06_results.tex (new main label fig:ranks)
fig:S-wakeloss -> main, optA/sw/06_results.tex (new main label fig:wakeloss)
fig:S-equiv-curve -> main, optA/sw/06_results.tex (new main label fig:equiv-curve)
tab:X-equiv-levels -> main, optA/sw/06_results.tex (label unchanged; captured from analysis/mpce_supp_inference.tex via \CaptureSuppTables + \SuppTable)
tab:X-heterogeneity -> main, optA/sw/06_results.tex (label unchanged; captured from analysis/mpce_supp_inference.tex; notes wrapped as in the supplement)
def:S-equiv -> main, optA/sw/05_setup.tex (label unchanged; Definition "Practical equivalence at three levels of inference" adapted from optA/supp_theory.tex sec:S-th-equiv; the supplement keeps Lemma lem:S-tost-ci, the four-outcome table and the Bayesian discussion, so it should cite M-def:S-equiv instead of defining it again)
(note, no label moved) sec:S-stats content (ranking rule, run scores, zero handling, Holm families) is now stated in main sec:setup-stats (optA/sw/05_setup.tex); the supplement subsection may be shortened to a pointer to Section M-sec:setup-stats
tab:split -> main, optA/sw/07_ablation.tex (from analysis/mpce_tab_split.tex; supplement alias M-tab:split now resolves to the main table)
tab:D-vns -> main, optA/sw/07_ablation.tex (captured from analysis/mpce_supp_diag.tex)
tab:D-rsreplay -> main, optA/sw/07_ablation.tex (captured from analysis/mpce_supp_diag.tex)
tab:switch -> main, optA/sw/07_ablation.tex
tab:equivalence -> main, optA/sw/07_ablation.tex
(note) tab:D-vns and tab:D-rsreplay are now placed in the main text by optA/sw/07_ablation.tex (this supersedes the remark in the fig:D-psodyn line that they stay in the supplement); tab:D-psodyn and tab:D-feasstart are not placed by 07
fig:S-budget -> main, optA/sw/08_beyond.tex (new main label fig:budget; figures_mpce/budget_scaling.pdf; the supplement should drop its figure and cite M-fig:budget)
tab:D-feasstart -> main, optA/sw/08_beyond.tex (label unchanged; captured from analysis/mpce_supp_diag.tex; the supplement's Section S-diag-feas and supp_theory.tex (S-th-feas) should cite M-tab:D-feasstart)
lem:S-feas-improve -> main, optA/sw/08_beyond.tex (new main label lem:feas-improve, with a short proof; the supplement may keep its full statement or point to M-lem:feas-improve)
prop:S-frozen -> main (parts (a) and (b) only), optA/sw/08_beyond.tex (new main label prop:feasible-start, with a short proof; part (c) and prop:S-de stay in the supplement; the DE crossover count 27.1 of 30 is quoted in the main text)
tab:robust-final -> main, optA/sw/09_robust.tex (label unchanged; \input{analysis/mpce_tab_robust_final.tex}; the supplement must no longer \input it in sec:S-robust-tables and should drop the sentence "Table tab:robust-final (moved here from the main text) ...", citing M-tab:robust-final instead; its full Table tab:robust stays)
tab:F-bench -> main, optA/sw/09_robust.tex (label unchanged; captured from analysis/mpce_supp_direction.tex in optA/sw/08_beyond.tex; the supplement's \input of mpce_supp_direction.tex must skip it)
(note, not moved: over the +3-page budget) tab:F-equiv, tab:F-hr, tab:F-pywake, tab:F-iea-gap, fig:S-hr16, fig:S-layouts-sites and fig:S-iea37-conv stay in the supplement and are cited from 08/09 as S-...; optA/sw/08_beyond.tex also captures analysis/mpce_supp_direction.tex (global) and defines r@M-tab:friedman68 / r@M-tab:hr-site as aliases so that the M- references inside captured captions resolve in the main text

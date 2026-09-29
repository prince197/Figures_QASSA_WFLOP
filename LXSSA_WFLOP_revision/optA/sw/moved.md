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

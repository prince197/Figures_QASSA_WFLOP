# Suggested manuscript text from the revision-3 analysis agents (starting points; adapt, shorten, keep numbers)

## Sites (B2, B6, C3) — rev3_sites.json
§9.3 (replace the GA paragraph): "In one pool of the ten methods and GA (Table~\ref{S-tab:S-r3-hr-pool}), PSO-VNS is best
at 6,030 evaluations (138.29 GWh/yr; GA 138.12, $p_{\rm Holm}=0.53$; better than the other nine,
$p_{\rm Holm}\le5.1\times10^{-4}$) and GA is best at 30,030 (139.64 GWh/yr; better than all ten,
$p_{\rm Holm}\le1.1\times10^{-4}$, largest against PSO-VNS, which is second and better than the remaining nine,
$p_{\rm Holm}\le0.009$). GA's lead at 30,030 holds when its layouts are re-evaluated with 1$^\circ$ bins
(+0.48 GWh/yr over PSO-VNS, $p=2.8\times10^{-4}$) and with PyWake NOJ (+0.43, $p=4.2\times10^{-4}$); at 6,030 no
difference between GA, PSO and PSO-VNS is significant under any evaluator (Table~\ref{S-tab:S-r3-hr-eval}). With
1$^\circ$ bins, none of the 961 feasible optimized layouts, GA included, exceeds the installed block (22 with 5$^\circ$
bins)." Delete "its layouts were not re-evaluated ... only under the 5° evaluator".
§9.4 (replace "the layouts were not re-evaluated with 1° bins..."): "Re-evaluated with 1$^\circ$ bins, the installed
block loses 1.65\% and the optimized layouts 1.14\% on average (PSO-VNS at 30,030: 1.69\%); PSO-VNS keeps the highest
mean AEP at 30,030 (112.64 GWh/yr) but is then no longer significantly better than VNS or MS-SLSQP
($p_{\rm Holm}=0.58$), and at 6,030 MS-SLSQP has the higher mean (Table~\ref{S-tab:S-r3-lg-eval})."
§9.5 (add): "In one pool of all 13 methods (Table~\ref{S-tab:S-r3-iea-pool}), an exact-gradient method ranks first in
every scenario and budget and is better than every gradient-free method in all 30 seed pairs
($p_{\rm Holm}=2.2\times10^{-8}$); the two exact-gradient methods differ significantly only at 30,030 evaluations.
Among the gradient-free methods, PSO-VNS ranks first with 16 turbines (not different from GA at 30,030,
$p_{\rm Holm}=0.26$) and GA with 36 turbines, where PSO-VNS ranks fourth or fifth and does not differ from VNS."
§9.7 (replace the per-evaluation sentence): "Per evaluation, the metaheuristics need 0.23–0.45 ms (median over runs) and
MS-SLSQP 0.53 ms, about twice as long as PSO when its runs are not slowed by thread contention: in the MS-SLSQP batch,
runs with $N\ge8$ took 16–22 s (3.4 ms per evaluation), whereas the same runs re-timed with single-threaded linear
algebra take 4–7 s and the earlier runs of the same method 3.3–4.0 s (Table~\ref{S-tab:S-r3-time}). The
feasibility-preserving initialization, which is not charged, takes a median of 4.8 s (up to 15.5 s) for the 30 layouts
of a run on the six largest cases, more than a 6,030-evaluation PSO-VNS run itself (2.8 s). The 36,190 original runs
took 101 CPU-hours in total and ran in four parallel processes."
Supplement: Table S-cost caption "mean wall-clock time per evaluation" + note on the N>=8 step and the re-timing;
Table S-time text: replace "PSO-VNS IEA37 records are missing" with ratios to PSO-VNS (tab:S-r3-time-sites); Section
S-rev-ga: replace "GA layouts not re-evaluated" with a pointer to tab:S-r3-hr-eval.
Other numbers: exact-gradient methods vs PSO-VNS time ratios: exact MS-SLSQP 3.3/5.6/4.9/4.9x, PSO-SLSQP 2.3/3.3/3.5/3.3x
(16T 6k, 16T 30k, 36T 6k, 36T 30k) — "about 3–6 times" holds. IEA37 PSO-VNS median run time 1.31/5.64/3.17/17.2 s.
Feasible-init run time PSO-VNS 7.9 s vs 2.8 s random (six largest). Revision records 12.4 CPU-h.

## Inference (B1, B3–B5, B7) — rev3_inference.json
§6.2 (TOST): "A null simulation at the margin (Table~\ref{S-tab:S-r3-tostaudit}) gives empirical sizes of 4.4–7.4\%
for PSO-VNS vs.\ PSO (case and seed level) and 4.7–6.0\% at the seed level for SSA-VNS vs.\ RSD-VNS, but 15.0\% at
the case level for the latter's non-superiority test, because of one outlying case (Data Set~I, $r=500$~m, $N=10$,
$d=-1.33$~pp); calibrated against this null distribution, the case-level lower limits still reject ($p=0.015$ and
$p=0.007$). The bootstrap $p_{\rm TOST}$ and the 90\% interval rule agree for all 18 pairs at every level."
§6.2 (seed dependence): "Seed $k$ seeds the same random stream in every case (identical initial populations for the two
data sets at equal $r$ and $N$). Resampling one seed vector jointly for all cases widens the seed-level intervals by up
to 1.9$\times$; it leaves PSO-VNS vs.\ PSO unchanged ([$-0.026$, $-0.010$]) but moves SSA-VNS vs.\ RSD-VNS to
[$-0.002$, 0.050], so its seed-level equivalence is not robust, whereas non-superiority is
(Table~\ref{S-tab:S-r3-syncseed})."
§6.2 (imputation): "Dropping these cases, or keeping zero differences (Pratt), changes no verdict of Tables 5 and 8
(Table~\ref{S-tab:S-r3-imputation})."
§7.3: "On an all-run paired endpoint (feasible beats infeasible, then the objective), PSO-VNS wins 950, ties 348 and
loses 742 of 2,040 seed pairs against PSO (score 0.551, 95\% case-bootstrap CI [0.516, 0.587], two-stage [0.511,
0.591]; Holm sign test $p=4.7\times10^{-7}$; Table~\ref{S-tab:S-r3-allrun}): it is better more often than not,
consistent with the seed-level interval, but by amounts within the margin."
§8: "On the all-run endpoint SSA-VNS scores 0.476 against RSD-VNS (case CI [0.452, 0.500], two-stage [0.445, 0.507]),
in the direction of the case means." / "On the 68 benchmark cases the same endpoint gives PSO-VNS scores of 0.73–0.90
against SSA-VNS, VNS, MS-SLSQP, SSA, LX-SSA and DE, all intervals above 0.5, and leaves the rank order unchanged except
that MS-SLSQP and SSA swap."
S-addl: "Recomputed from the run records (Table~\ref{tab:S-r3-eqclus}), the values agree with the reconstruction within
0.0005; the restricted wild-cluster bootstrap-t for the equal-cluster mean (all $6^6$ Webb draws) gives PSO-VNS $-$ PSO
[$-0.062$, 0.036] ($p_{\rm TOST}=0.108$) and SSA-VNS $-$ RSD-VNS [$-0.067$, 0.082]; no verdict changes."
Footnote b of tab:equiv-main: add "(not under joint seed resampling, Table~\ref{S-tab:S-r3-syncseed})".

## Constraint handling (C1) — rev3_constraint.json
§3.4 (after "...ranked below feasible ones"): "Because this penalty mixes units and depends on $\mu$, a targeted study
(Table~\ref{S-tab:S-r3-constraint}) reran PSO-VNS, PSO, GA, SSA-VNS and DE on six cases with new seeds 31–60 under two
alternatives: Deb's feasibility rules on dimensionless violations $g^{\rm b}/r^2$ and $g^{\rm s}/\ell_{\min}$ (feasible
before infeasible, then wake loss, then total violation), with box clipping, and the same rules with radial projection
of turbines onto the circle (Section~\ref{sec:boundary})."
§10.4 (replace last sentence): see analysis report: "The rule matters for all methods. In the constraint-handling study
(Table~\ref{S-tab:S-r3-constraint}; six cases, seeds 31–60, 6,030 evaluations), Deb's rules with radial projection
lowered the conditional wake loss of every method by 0.46–0.80 pp and raised the feasibility of DE from 33.3\% to
80.0\%; the run-level rank order PSO-VNS, PSO, GA, SSA-VNS, DE was unchanged, and PSO-VNS remained significantly better
than GA, SSA-VNS and DE in the all-run paired score (0.683, 0.850, 0.844), but not than PSO (0.561; conditional
$\Delta L=-0.095$~pp [$-0.162$, $-0.037$]). With Deb's rules and box clipping, by contrast, a turbine clipped to the
square near a tangent point is cheap: PSO's feasibility fell from 98.9\% to 73.3\% (mostly boundary violations),
PSO-VNS's conditional loss rose by 0.23~pp, and its all-run advantage over GA (0.517) and SSA-VNS (0.583) was no longer
significant; GA had the best case-level rank (1.83 vs.\ 2.17), PSO-VNS the best run-level rank (2.19 vs.\ 2.28)."
§11 / §12: see REV3_INTEGRATION.md D10.

## Reproducibility (D) — repro_package_build/
Data availability (replace): "All code and data of this study are contained in a versioned reproducibility package
(release 1.0.0): the per-run records of every optimization run, with final layouts and convergence curves (36,190 runs
of the original study, 11,640 runs of the first revision and [N3] runs added in this revision), the benchmark, Horns
Rev~1, IEA37 and Lillgrund evaluators, the implementations of all compared methods, the experiment drivers, the
analysis and check scripts that generate every table, figure and number of this article and its supplement, an
environment specification and SHA-256 checksums of all files. The package will be deposited publicly at \TBD{DOI}
before publication; until then it is available to the editor and reviewers from the corresponding author. Horns Rev~1
and Lillgrund inputs derive from PyWake~\cite{PyWake} (Lillgrund wind climate:~\cite{Gocmen2016}), and IEA37 inputs
from~\cite{IEA37repo}."
§6.4 (replace "Archived and independently executable checks"): "\emph{Generated results and automated checks.} Every
table, figure and number of the article and the supplement is generated from the stored per-run records by the analysis
scripts of the reproducibility package (Data availability), and every outcome-dependent statement is linked to an
automated condition. For this revision the complete analysis pipeline was re-executed from the package (36,190
original and 11,640 revision records; CPython 3.11.15, NumPy 2.4.6, SciPy 1.17.1, pandas 3.0.6): all regenerated
tables, number macros, summaries and figures are identical to those used in the manuscript apart from generation time
stamps, and all 62 main-result, 57 statistical, 22 diagnostic, 41 direction-resolution, 15 coefficient-sweep and 57
theory checks pass; the analysis of the revision experiments reproduces its summary (4,670 values) and tables exactly,
and the audit of the revision records (complete shards, unique keys, evaluation counts, objective identities) passes.
Reruns of stored runs with the archived drivers reproduce the records bit for bit in every field except the elapsed
time." Then "Historically reported instrumented diagnostics" -> "Instrumented diagnostics"; "Historically reported
evaluator and code validation" -> "Evaluator and code validation" (+ "(both re-executed for this revision; the Horns
Rev~1 PyWake comparison reproduces the stored values)").
Supplement S-repro (retitle "Reproducibility Map"; first paragraph): "This map lists the script that produces each part
of the results. All scripts, run records and generated outputs are contained in the reproducibility package (Data
availability); its README gives, for every table, figure and number of the article, the command, the input files and
the run time, and for every optimization study the exact command and shard layout. For this revision the complete
analysis pipeline was re-executed from the package in the environment of the runs (CPython 3.11.15, NumPy 2.4.6, SciPy
1.17.1, pandas 3.0.6, Matplotlib 3.11.2): all 91 regenerated or compared outputs are byte-identical to those used in
the manuscript or differ only in generation time stamps, PDF creation dates and recorded software versions, and all
automated checks pass (Section~\ref{sec:S-archive}). The public versioned deposit (DOI) is pending."
S-archive "Coverage and limits" (replace): package content (36,190 + 11,640 + rev3 records, all modules, generated
outputs); audit_archive.py checks pass; complete analysis pipeline re-executed, all checks pass; reruns bit-identical
except elapsed time; NOT re-executed: the optimization runs themselves (about 126 CPU-hours), the instrumented
diagnostic reruns (stored outputs pass D01–D22) and the calibration reruns of the original code; public DOI pending.
"Serialization precision": delete "The Horns Rev~1 boundary could not be audited without its missing model"; add the
C4 audit/rerun results (REV3_INTEGRATION.md D12).

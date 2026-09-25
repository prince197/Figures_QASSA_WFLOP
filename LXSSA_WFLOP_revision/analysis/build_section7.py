"""Rebuild Section VII (controlled re-runs) of the manuscript from the generated tables.
Numbers quoted in the text were checked against authors_runs_*.csv / extra_runs_*.csv."""
import sys

MS = "../LXSSA_WFLOP_reviewer_revised.tex"


def blocks(path):
    tabs = open(path).read()
    return {b.split("\\label{")[1].split("}")[0]: ("\\begin{table" + b).rstrip()
            for b in tabs.split("\\begin{table") if "\\label{" in b}


def main(hr_text):
    A = blocks("authors_tables.tex")
    X = blocks("extra_tables.tex")
    sec = r"""\section{Controlled Re-runs with the Original Implementation}
\label{sec:equalbudget}
The follow-up data of Section~\ref{sec:additional-runlevel} compare algorithms at unequal budgets. To remove this confound, the original implementations of SSA, LX-SSA, PSO and DE, i.e., the code that produced the follow-up data, were re-run with the original objective function~\eqref{eq:penalty}. The objective was re-implemented in vectorized form and checked against the original code (relative difference below $10^{-6}$), and its wake and energy terms reproduce the recorded follow-up objective values exactly (Section~\ref{sec:validation}).

\subsection{Protocol and Calibration}
\label{sec:rerun-protocol}
\emph{Calibration.} At the recorded budgets (6,030 calls for LX-SSA, 3,030 for the others) the re-runs reproduce the recorded follow-up distributions. In 22 of the 24 algorithm--case comparisons a Mann--Whitney test does not distinguish re-run from recorded objective values at the 5\% level (smallest $p=0.0061$), and the largest relative difference between re-run and recorded means is 0.21\%. The counted objective calls are exactly 3,030 and 6,030, as recorded. This confirms that the reported LX-SSA results correspond to $\phi=0$ and $\chi=1$. The only notable difference is that 26 of the 30 DE re-runs are feasible in each 750-m eight-turbine case, whereas all recorded DE runs there are feasible.

\emph{Design.} Each algorithm uses a population of 30 and seeds 1--30. Every optimizer seeds the random-number generator and draws its initial population first, so runs with the same seed start from the \emph{same} initial population for all algorithms (verified numerically). The runs are therefore paired by seed, which justifies the Wilcoxon signed-rank and run-level Friedman tests that the independent-sample analysis of Section~\ref{sec:additional-runlevel} could not use. Two budgets are matched in counted objective calls: 3,030 (LX-SSA 50 iterations; SSA, PSO and DE 100) and 6,030 (LX-SSA 100 iterations; the others 200). For ranking, a feasible run beats an infeasible one, and infeasible runs are compared by their spacing shortfall.

\subsection{Equal-Budget Comparison}
""" + A["tab:equalbudget"] + "\n" + A["tab:equalbudget-tests"] + r"""
Tables~\ref{tab:equalbudget} and~\ref{tab:equalbudget-tests} show that the advantage of LX-SSA seen at unequal budgets does not persist at equal cost:
\begin{itemize}
\item \emph{LX-SSA vs.\ SSA (mechanism-level ablation).} SSA has the higher mean objective in all six cases at both budgets. At 3,030 calls it is significantly better in two cases (Data Set~I at 500~m/4, $p_{\rm Holm}=0.012$, and at 750~m/8, $p_{\rm Holm}=0.020$); at 6,030 calls no difference is significant. Under an equal number of objective calls, the Laplace follower candidate therefore does not improve on the SSA follower update for this problem. A plausible reason is that LX-SSA spends two calls per follower per iteration, so at a fixed budget it performs only half as many iterations.
\item \emph{LX-SSA vs.\ PSO.} PSO is significantly better in Data Set~I at 500~m/4 at 3,030 calls and in both 500~m/4 cases at 6,030 calls; LX-SSA is better only in Data Set~II at 1000~m/8 at 3,030 calls ($p_{\rm Holm}=0.050$). All other comparisons are not significant.
\item \emph{LX-SSA vs.\ DE.} LX-SSA is significantly better in five of six cases at 3,030 calls and in all four eight-turbine cases at 6,030 calls, with large effects ($r_{rb}\ge0.69$).
\end{itemize}
Averaged over the six cases, the Friedman rank of LX-SSA is 2.28 at 3,030 calls (SSA 1.88, PSO 2.20, DE 3.64) and 2.36 at 6,030 calls (SSA 2.13, PSO 2.15, DE 3.36). LX-SSA is thus clearly better than DE and close to, but not better than, SSA and PSO.

\subsection{Additional Baselines: VNS and Gradient-Based Optimization}
\label{sec:vns-slsqp}
""" + X["tab:vns-slsqp"] + r"""
Reviewer~3 asked for comparison with variable neighbourhood search and gradient-based methods. Two such baselines were implemented with the same objective, bounds, seeds, initial populations and budgets:
\begin{itemize}
\item \emph{VNS}~\cite{Mladenovic1997,Cazzaro2022}: basic variable neighbourhood search adapted to continuous coordinates. It starts from the best member of the shared initial population. Shaking moves $k\in\{1,2,3\}$ randomly chosen turbines by Gaussian steps with standard deviation $0.05r$, $0.15r$ and $0.40r$. Each shake is followed by a ten-call first-improvement local search of single-turbine Gaussian moves with an adaptive step (enlarged by 1.5 after a success and reduced by 0.8 after a failure). An improvement is accepted and resets $k=1$; otherwise $k$ increases cyclically.
\item \emph{Multistart SLSQP (MS-SLSQP)}: SciPy's sequential least-squares quadratic programming solver~\cite{Kraft1988} minimizes the wake loss with the spacing and circular-boundary constraints imposed explicitly. Gradients are forward differences with a 1-m step, and every objective evaluation, including those for gradients, counts against the budget. Starts are taken from the shared initial population in order of $F_p$, and the solver is restarted until the budget is exhausted.
\end{itemize}
Table~\ref{tab:vns-slsqp} shows that VNS is the strongest method in the representative cases. It is significantly better than LX-SSA in five of the six cases at both budgets (all except Data Set~II at 500~m/4), with large effects ($r_{rb}$ between $-0.49$ and $-0.84$). Over the six cases it has the best average rank among all six methods (2.09 at 3,030 calls and 2.28 at 6,030 calls; LX-SSA 3.52 and 3.46). MS-SLSQP performs worse: LX-SSA is significantly better in two cases at each budget and not different in the others, and MS-SLSQP ranks fifth of six. This is consistent with the non-smooth objective: the top-hat Jensen wake and the 24 discrete direction bins make finite-difference gradients zero or discontinuous over much of the design space.

\subsection{Crowded Cases and Feasibility}
""" + X["tab:crowded-all"] + r"""
For the largest turbine count of each historical table, Table~\ref{tab:crowded-all} reports how often each method returns a feasible layout at 6,030 calls. The two salp-swarm algorithms are far more reliable than PSO and DE: LX-SSA returns feasible layouts in 13, 28 and 29 of 30 runs for the 500-m, 750-m and 1000-m cases, and SSA in 12, 30 and 30, whereas PSO succeeds in at most 9 and DE in at most one run. VNS and MS-SLSQP are more reliable still (27--30 feasible runs in every case), and one of them attains the highest mean objective in every crowded case. Among the salp-swarm algorithms, SSA attains the higher mean objective in all six cases.

\subsection{Measured Wind Climate with the Benchmark Turbine}
\label{sec:hornsrev}
""" + A["tab:hornsrev"] + r"""
As a first step towards measured wind data, the four algorithms were applied to the measured wind climate of the Horns Rev~1 offshore wind farm (12 directional sectors with Weibull scale, shape and frequency per sector, from the open-source PyWake distribution~\cite{PyWake}), keeping the benchmark turbine, wake model and circular farm. Table~\ref{tab:hornsrev} shows the same pattern as the synthetic data. SSA has the highest mean in the three larger cases. LX-SSA does not differ significantly from SSA except for 1000~m/15, where SSA is better ($p_{\rm Holm}=0.005$; mean wake losses 4.88\% vs.\ 4.25\% of the ideal, AEP 98.7 vs.\ 99.3~GWh). LX-SSA is significantly better than DE in the three larger cases and better than PSO at 750~m/8 and 1000~m/15; at 1000~m/15, PSO finds a feasible layout in only 9 of 30 runs and DE in none.

""" + hr_text + r"""
\subsection{Measured Run Time}
Each run was timed individually. At equal numbers of objective calls the population-based methods take essentially the same time: 0.41~ms per call for LX-SSA and SSA, 0.42~ms for PSO, 0.43~ms for VNS and 0.48~ms for DE, i.e., about 1.8--3.7~s per 6,030-call run depending on $N$. MS-SLSQP needs 0.51~ms per call because of its quadratic subproblems. The run time is dominated by the objective evaluation, and the Laplace step adds no measurable overhead. The large differences in the recorded batch times of Table~\ref{tab:followup-results} (up to 15 times longer for LX-SSA) are therefore not reproduced and are attributable to the undocumented recording conditions rather than to the algorithms. Environment: Linux container, Intel Xeon processor at 2.80~GHz, 4 virtual cores with four runs executed in parallel, Python~3.11.15, NumPy~2.4.6, SciPy~1.17.1.

\subsection{Effect of Boundary Handling}
The implementation clips coordinates to the bounding square and penalizes points outside the circle (Section~\ref{sec:constraints}). An otherwise identical SSA that instead projects every turbine radially onto the circle and ranks candidates by the feasibility-first rules attains higher mean objectives than the SSA re-runs in all six representative cases at 3,030 calls, by 26 to 649 objective units (Mann--Whitney $p<10^{-4}$ in every case; for example 111774.0 vs.\ 111367.1 for Data Set~I at 750~m/8). The boundary-handling rule therefore affects results by more than the differences between the salp-swarm algorithms, and it must be held fixed, as here, in any comparison.

"""
    t = open(MS).read()
    i = t.index(r"\section{Controlled Re-runs with the Original Implementation}")
    j = t.index(r"\section{Robustness of the Benchmark Assumptions}")
    t = t[:i] + sec + t[j:]
    open(MS, "w").write(t)
    print("Section VII rebuilt")


if __name__ == "__main__":
    main(open(sys.argv[1]).read() if len(sys.argv) > 1 else "")

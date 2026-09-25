"""Front matter (title, abstract, introduction, motivation) for the hybrid LX-SSA-VNS manuscript."""
MS = "../LXSSA_WFLOP_reviewer_revised.tex"
t = open(MS).read()


def between(t, start, end, new):
    i = t.index(start); j = t.index(end, i)
    return t[:i] + new + t[j:]


t = t.replace(r"\title{Application of the Laplacian Salp Swarm Algorithm to Continuous Wind Farm Layout Optimization}",
              r"\title{A Hybrid Laplacian Salp Swarm and Variable Neighbourhood Search Algorithm for Continuous Wind Farm Layout Optimization}")

t = between(t, r"\begin{abstract}", r"\end{abstract}", r"""\begin{abstract}
The Wind Farm Layout Optimization Problem (WFLOP) is a nonlinear constrained optimization problem in which turbine locations interact through wake-induced velocity deficits. This paper proposes LX-SSA-VNS, a two-phase hybrid for continuous WFLOP. The Laplacian Salp Swarm Algorithm (LX-SSA), previously introduced by Solanki and Deep~\cite{Solanki2023}, explores the layout space during the first half of the evaluation budget. Variable neighbourhood search (VNS) then intensifies the best layout found, using nested turbine-perturbation neighbourhoods and an adaptive local search. Both phases share a penalized objective with Jensen wake interactions, direction-dependent Weibull wind statistics and explicit farm-boundary and minimum-spacing constraints. The hybrid is compared with LX-SSA, SSA, particle swarm optimization (PSO), differential evolution (DE), stand-alone VNS and a gradient-based multistart SLSQP solver on 68 benchmark cases (two wind data sets, three circular farms, 2--15 turbines). Each comparison uses 30 seed-paired runs at an equal budget of 6,030 objective evaluations and is analyzed with Friedman, Holm-corrected Wilcoxon signed-rank and effect-size statistics, together with feasibility rates, convergence curves and an ablation of the budget split. LX-SSA-VNS attains the second-best average rank (2.38) and a feasibility rate of 98.1\%. It is significantly better than LX-SSA in 36 cases and never worse, and it ranks significantly better than SSA, PSO, DE, LX-SSA and SLSQP. Stand-alone VNS remains the best method (average rank 1.85): it is statistically equivalent to the hybrid in 52 cases and better in 16 dense cases, and the ablation attributes the gain of the hybrid to its VNS phase. On the Horns Rev~1 offshore site the same ordering appears, no method improves on the installed layout within the budget, and only the constraint-aware SLSQP finds feasible layouts for the full 80-turbine farm. The algorithm ordering is insensitive to the power-curve model but sensitive to the wake model, and absolute energy values are strongly model-dependent.
""")

t = t.replace(r"Laplacian Salp Swarm Algorithm (LX-SSA), Optimization, Renewable energy, Wake effect, Wind farm layout.",
              r"Hybrid metaheuristic, Laplacian Salp Swarm Algorithm (LX-SSA), Variable neighbourhood search, Wake effect, Wind farm layout optimization.")

t = between(t, r"Solanki and Deep~\cite{Solanki2023} introduced", r"Renewable energy systems are increasingly", r"""Solanki and Deep~\cite{Solanki2023} introduced the Laplacian Salp Swarm Algorithm (LX-SSA) in 2023 as an enhanced variant of SSA in which a Laplace-distribution-based position update is used to improve exploration and population diversity. The design, motivation and benchmark-function validation of LX-SSA belong to that earlier publication. Like other population-based methods, LX-SSA spends many evaluations on small improvements once its swarm has contracted around the food source. Trajectory-based methods such as variable neighbourhood search (VNS)~\cite{Mladenovic1997} have the opposite profile: they refine a single incumbent efficiently but depend on a good starting point. The present paper combines the two in a hybrid, LX-SSA-VNS, and evaluates it on continuous WFLOP, where each candidate solution represents turbine coordinates and must satisfy wake-dependent and geometric feasibility requirements.

The standard WFLOP ingredients---Jensen's wake model, a Weibull wind representation, and a benchmark turbine power model---are retained deliberately to permit comparison with the established WFLOP literature.

\textbf{The contributions of this paper are summarized as follows:}
\begin{itemize}
\item a two-phase hybrid, LX-SSA-VNS, in which the previously published LX-SSA~\cite{Solanki2023} explores the continuous layout space and a VNS with nested turbine-perturbation neighbourhoods and an adaptive local search intensifies the best layout (Section~\ref{sec:hybrid});
\item a continuous turbine-coordinate formulation with explicit farm-boundary and minimum-spacing constraints, shared by both phases, within the standard Jensen--Weibull WFLOP benchmark framework;
\item a controlled comparison of the hybrid with LX-SSA, SSA, PSO, DE, stand-alone VNS and a gradient-based multistart SLSQP solver on 68 benchmark cases at equal objective-evaluation budgets, with 30 seed-paired runs per method and case, full statistical analysis (Friedman and Holm-corrected Wilcoxon signed-rank tests, effect sizes), feasibility rates and convergence curves (Section~\ref{sec:results});
\item an ablation that separates the contributions of the two phases and examines the budget split (Section~\ref{sec:ablation});
\item an application to the Horns Rev~1 offshore site with its real turbine, wind climate and outline, compared against the installed layout (Section~\ref{sec:hornsrev-site}); and
\item a sensitivity analysis of the benchmark assumptions (power curve, wake model, minimum spacing and boundary handling; Section~\ref{sec:robustness}).
\end{itemize}

""")

t = between(t, r"Against this background, the purpose of the present study", r"\section{Suitability of LX-SSA for WFLOP}", r"""Against this background, the purpose of the present study is deliberately narrower than proposing a new WFLO paradigm. LX-SSA was already introduced and benchmarked as a continuous optimizer in 2023~\cite{Solanki2023}. Here it is combined with VNS, and the hybrid is examined in a standard continuous WFLO benchmark under geometric constraints and wake-coupled objectives. Jensen's model and the benchmark Weibull/power formulation are retained primarily for comparability with the established literature; they are not presented as state-of-the-art wake physics or as novel modeling contributions.


""")

t = between(t, r"\section{Suitability of LX-SSA for WFLOP}", r"\section{\textbf{Problem Modeling}}", r"""\section{Motivation for the Hybrid Design}

Continuous WFLOP contains $2N$ spatial decision variables for $N$ turbines, nonlinear wake coupling, and geometric feasibility constraints. Such features can produce many locally attractive layouts and make diversity preservation relevant to a population-based search.

In the 2023 LX-SSA, the standard SSA follower dynamics are supplemented by a Laplace-distribution-based candidate update and greedy retention~\cite{Solanki2023,Deep2007}. The heavy-tailed Laplace perturbation can generate both short and occasional longer moves, which provides a plausible mechanism for exploring spatially separated layout configurations. In a swarm, however, these moves are spread over the whole population, so the final refinement of the best layout, which in WFLOP means small coordinated moves of individual turbines along the spacing constraints, proceeds slowly.

VNS~\cite{Mladenovic1997}, which has been applied successfully to discrete WFLOP~\cite{Cazzaro2022}, addresses exactly this refinement: it perturbs one, two or three turbines at a time with increasing step sizes and keeps only improvements. Its weakness is its dependence on the starting layout. A sequential hybrid, in which LX-SSA supplies the starting layout for VNS, is therefore a natural combination. Whether it pays off at a fixed budget, compared with running either component alone, is an empirical question: the experiments answer it by comparing the hybrid with LX-SSA, with stand-alone VNS and with four further baselines at equal objective-evaluation budgets (Section~\ref{sec:results}). Pseudo-gradient, surrogate and data-driven WFLO methods remain literature context.

The remainder of the paper presents the benchmark model, summarizes the previously published LX-SSA to the extent needed for reproducibility, defines its WFLOP mapping and constraint handling, introduces the hybrid, reports the experimental results, and finally states the limitations of the evidence.

""")
open(MS, "w").write(t)
print("front matter updated")

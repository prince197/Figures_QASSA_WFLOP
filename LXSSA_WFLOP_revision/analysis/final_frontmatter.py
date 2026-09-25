"""Front matter (abstract, contributions, motivation) for the final manuscript
(hybrid = LX-SSA + original basic VNS; component ablation)."""
MS = "../LXSSA_WFLOP_reviewer_revised.tex"
t = open(MS).read()


def between(t, start, end, new):
    i = t.index(start); j = t.index(end, i)
    return t[:i] + new + t[j:]


t = between(t, r"\begin{abstract}", r"\end{abstract}", r"""\begin{abstract}
The Wind Farm Layout Optimization Problem (WFLOP) is a nonlinear constrained optimization problem in which turbine locations interact through wake-induced velocity deficits. This paper proposes LX-SSA-VNS, a two-phase hybrid for continuous WFLOP. The Laplacian Salp Swarm Algorithm (LX-SSA), previously introduced by Solanki and Deep~\cite{Solanki2023}, explores the layout space during the first half of the evaluation budget. The original basic variable neighbourhood search (VNS) then intensifies the best layout found, alternating shaking of the whole layout in nested neighbourhoods with a complete local search. Both phases share a penalized objective with Jensen wake interactions, direction-dependent Weibull wind statistics and explicit farm-boundary and minimum-spacing constraints. The hybrid is compared with LX-SSA, SSA, particle swarm optimization (PSO), differential evolution (DE), stand-alone VNS and a gradient-based multistart SLSQP solver on 68 benchmark cases (two wind data sets, three circular farms, 2--15 turbines). Each comparison uses 30 seed-paired runs at an equal budget of 6,030 objective evaluations and is analyzed with Friedman, Holm-corrected Wilcoxon signed-rank and effect-size statistics, together with feasibility rates and convergence curves. LX-SSA-VNS attains the best average rank (1.89), ranks significantly better than all six other methods, including VNS, and returns a feasible layout in 99.9\% of the runs; it is significantly better than LX-SSA in 39 cases and never worse. A component ablation of SSA, LX-SSA, VNS, SSA-VNS and LX-SSA-VNS shows that the gain comes from combining a swarm exploration phase with the VNS phase, which removes on average 39\% of the remaining wake loss, whereas the Laplace step itself does not contribute: an SSA first phase performs at least as well. On the Horns Rev~1 offshore site the hybrid has the highest mean energy yield in a 16-turbine block, no method improves on the installed layout within the budget, and only the constraint-aware SLSQP finds feasible layouts for the full 80-turbine farm. The algorithm ordering is insensitive to the power-curve model but sensitive to the wake model, and absolute energy values are strongly model-dependent.
""")

R = [
    (r"""\item a two-phase hybrid, LX-SSA-VNS, in which the previously published LX-SSA~\cite{Solanki2023} explores the continuous layout space and a VNS with nested turbine-perturbation neighbourhoods and an adaptive local search intensifies the best layout (Section~\ref{sec:hybrid});""",
     r"""\item a two-phase hybrid, LX-SSA-VNS, in which the previously published LX-SSA~\cite{Solanki2023} explores the continuous layout space and the original basic VNS~\cite{Mladenovic1997,Mladenovic2008} intensifies the best layout (Section~\ref{sec:hybrid});"""),
    (r"""\item a controlled comparison of the hybrid with LX-SSA, SSA, PSO, DE, the original basic VNS, a modified stand-alone VNS and a gradient-based multistart SLSQP solver on 68 benchmark cases""",
     r"""\item a controlled comparison of the hybrid with LX-SSA, SSA, PSO, DE, stand-alone VNS and a gradient-based multistart SLSQP solver on 68 benchmark cases"""),
    (r"""\item a comparison of the original and the modified VNS (Section~\ref{sec:bvns}) and an ablation that separates the contributions of the two phases and examines the budget split (Section~\ref{sec:ablation});""",
     r"""\item a component ablation of SSA, LX-SSA, VNS, SSA-VNS and LX-SSA-VNS that separates the contributions of the swarm phase, the Laplace step and the VNS phase, and examines the budget split (Section~\ref{sec:ablation});"""),
    (r"""VNS~\cite{Mladenovic1997}, which has been applied successfully to discrete WFLOP~\cite{Cazzaro2022}, addresses exactly this refinement: it perturbs one, two or three turbines at a time with increasing step sizes and keeps only improvements.""",
     r"""VNS~\cite{Mladenovic1997,Hansen2001}, which has been applied successfully to discrete WFLOP~\cite{Cazzaro2022}, addresses exactly this refinement: it alternates random perturbations of increasing size with a complete local search and keeps only improvements."""),
    (r"""the experiments answer it by comparing the hybrid with LX-SSA, with stand-alone VNS in its original and modified forms and with four further baselines at equal objective-evaluation budgets (Section~\ref{sec:results}).""",
     r"""the experiments answer it by comparing the hybrid with LX-SSA, with stand-alone VNS and with four further baselines at equal objective-evaluation budgets, and by an ablation that also includes an SSA-VNS variant (Section~\ref{sec:results})."""),
]
for a, b in R:
    if a in t:
        t = t.replace(a, b)
    else:
        assert b in t, a[:60]
open(MS, "w").write(t)
print("front matter updated (final)")

"""Rebuild the manuscript around the proposed hybrid LX-SSA-VNS (eight methods, incl. the original basic VNS).
Numbers were checked against hybrid_summary.json, hybrid_case_tests.csv, hybrid_robust_summary.json
and the fresh_*.csv run files."""
import pandas as pd

MS = "../LXSSA_WFLOP_reviewer_revised.tex"


def blocks(path):
    txt = open(path).read()
    return {b.split("\\label{")[1].split("}")[0]: ("\\begin{table" + b).rstrip()
            for b in txt.split("\\begin{table") if "\\label{" in b}


def fig(name, label, caption, wide=True, width="\\textwidth"):
    env = "figure*" if wide else "figure"
    return (f"\\begin{{{env}}}[!t]\n\\centering\n\\includegraphics[width={width}]{{figures_hybrid/{name}.pdf}}\n"
            f"\\caption{{{caption}}}\n\\label{{{label}}}\n\\end{{{env}}}\n")


T = blocks("hybrid_tables.tex")
T["tab:wtl"] = T["tab:wtl"].replace("{2.5pt}", "{1.6pt}")
BV = blocks("bvns_tables.tex")
RB = blocks("hybrid_robust_table.tex")
CAP = blocks("robustness_tables.tex")
SP = blocks("authors_tables.tex")


def setup_text():
    s = open("sec_setup.tex").read()
    reps = [
        ("\\subsection{Methods Compared}\nSix methods are compared, all with the same objective~\\eqref{eq:penalty}, the same bounds $[-r,r]^{2N}$ and the same budget:\n\\begin{itemize}\n",
         "\\subsection{Methods Compared}\nThe proposed hybrid and seven reference methods are compared, all with the same objective~\\eqref{eq:penalty}, the same bounds $[-r,r]^{2N}$ and the same budget:\n\\begin{itemize}\n\\item \\emph{LX-SSA-VNS} (proposed): Algorithm~\\ref{alg:hybrid} with $\\rho=0.5$, i.e., 50 LX-SSA iterations (3,030 calls) followed by 3,000 calls of VNS.\n"),
        ("SSA, LX-SSA, PSO and DE use the original implementation; VNS and MS-SLSQP were implemented for this study.",
         "SSA, LX-SSA, PSO and DE use the original implementation. VNS, BVNS, MS-SLSQP and the hybrid were implemented for this study, and the hybrid calls the original LX-SSA code in its first phase. Comparing LX-SSA-VNS with LX-SSA and with stand-alone VNS isolates the contribution of each phase, and comparing VNS with BVNS shows the effect of the modifications made to the original VNS."),
        ("so runs with the same seed start from the \\emph{same} initial population for all six methods.",
         "so runs with the same seed start from the \\emph{same} initial population for all eight methods."),
        ("For LX-SSA this corresponds to 100 iterations, because each follower evaluates two candidates before the population is re-evaluated ($30+100\\times60$ calls); SSA, PSO and DE reach the same count in 200 iterations ($30+200\\times30$).",
         "For LX-SSA this corresponds to 100 iterations, because each follower evaluates two candidates before the population is re-evaluated ($30+100\\times60$ calls); SSA, PSO and DE reach the same count in 200 iterations ($30+200\\times30$), and LX-SSA-VNS in $3{,}030+3{,}000$ calls."),
        ("Within each test case, the 30 seed-paired runs are compared by a run-level Friedman test and by two-sided Wilcoxon signed-rank tests of LX-SSA against each of the other five methods. The five $p$ values of each case are adjusted by Holm's procedure~\\cite{Holm1979},",
         "Within each test case, the 30 seed-paired runs are compared by a run-level Friedman test and by two-sided Wilcoxon signed-rank tests of LX-SSA-VNS against each of the other seven methods. The seven $p$ values of each case are adjusted by Holm's procedure~\\cite{Holm1979},"),
        ("\\item \\emph{VNS}~\\cite{Mladenovic1997,Cazzaro2022}: basic variable neighbourhood search adapted to continuous coordinates. It starts from the best member of the initial population. Shaking moves $k\\in\\{1,2,3\\}$ randomly chosen turbines by Gaussian steps with standard deviation $0.05r$, $0.15r$ and $0.40r$. Each shake is followed by a ten-call first-improvement local search of single-turbine Gaussian moves with an adaptive step (enlarged by 1.5 after a success, reduced by 0.8 after a failure). An improvement is accepted and resets $k=1$; otherwise $k$ increases cyclically.",
         "\\item \\emph{VNS} (modified; this study): a variant of basic VNS~\\cite{Mladenovic1997,Cazzaro2022} tailored to turbine coordinates. It starts from the best member of the initial population. Shaking moves $k\\in\\{1,2,3\\}$ randomly chosen turbines by Gaussian steps with standard deviation $0.05r$, $0.15r$ and $0.40r$. Instead of a complete local search, each shake is followed by only ten first-improvement trial moves of single turbines with an adaptive step (enlarged by 1.5 after a success, reduced by 0.8 after a failure). An improvement is accepted and resets $k=1$; otherwise $k$ increases cyclically."),
        ("minimizes the wake loss with the spacing and boundary constraints imposed explicitly and forward-difference gradients (1-m step). It starts from the members of the initial population in order of $F_p$ and restarts until the budget is exhausted.",
         "minimizes the wake loss with the spacing and boundary constraints imposed explicitly and forward-difference gradients (1-m step). It starts from the members of the initial population in order of $F_p$ and restarts until the budget is exhausted.\n\\item \\emph{BVNS} (original basic VNS): the basic VNS of Mladenovi\\'c and Hansen~\\cite{Mladenovic1997,Hansen2001} in the continuous form of Mladenovi\\'c et al.~\\cite{Mladenovic2008}, with the exterior penalty $F_p$ for the constraints. It uses $k_{\\max}=5$ nested $\\ell_\\infty$ neighbourhoods of the whole layout, $\\rho_{k-1}<\\lVert y-x\\rVert_\\infty\\le\\rho_k$ with $\\rho_k=0.1kr$; shaking draws a point uniformly from the $k$-th neighbourhood, so all $2N$ coordinates move. It is followed by a complete best-improvement local search to a local minimum: a derivative-free compass search evaluates all $4N$ moves $\\pm h$ along the coordinates, moves to the best improving one, and halves $h$ (from $0.05r$) when none improves, until $h<10^{-3}r$. (Gradient descent, used for smooth problems in~\\cite{Mladenovic2008}, is not applicable because the Jensen top-hat wake makes $F_p$ discontinuous.) BVNS first applies this local search to the best member of the initial population; neighbourhood change is as in VNS, and the returned layout is the best point evaluated."),
        ("Re-running the original SSA, LX-SSA, PSO and DE code at the archived settings reproduces the archived distributions",
         "Re-running the original SSA, LX-SSA, PSO and DE code at the archived settings reproduces the archived distributions"),
    ]
    for a, b in reps:
        assert s.count(a) == 1, a[:60]
        s = s.replace(a, b)
    return s


results = r"""\section{Numerical Results}
\label{sec:results}
\label{sec:equalbudget}
This section reports 16,320 runs: 68 test cases, eight methods and 30 seed-paired runs, all at 6,030 objective calls. It is followed by a comparison of the original and modified VNS, an ablation of the hybrid and the Horns Rev~1 site case. Detailed per-case results are given in Tables~\ref{tab:res-1-500}--\ref{tab:res-2-1000}.

\subsection{Overall Comparison}
\label{sec:overall}
""" + T["tab:friedman68"] + "\n" + fig("avg_ranks", "fig:avgranks",
    "Average rank of the eight methods over the 68 benchmark cases (ranks of the mean feasible objective within each case; 1 = best).",
    wide=False, width="0.95\\columnwidth") + T["tab:wtl"] + r"""
The case-level Friedman test over the 68 cases rejects equal performance of the eight methods ($\chi^2_F=212.5$, 7 d.f., $p=2.6\times10^{-42}$; Iman--Davenport $F_F=54.0$). Table~\ref{tab:friedman68} and Fig.~\ref{fig:avgranks} give the average ranks (in parentheses: cases in which the method alone has the highest mean):
\begin{itemize}
\item VNS (modified) 2.09 (36 cases);
\item \textbf{LX-SSA-VNS 2.74} (5 cases);
\item BVNS (original VNS) 3.62 (10 cases);
\item SSA 4.79, PSO 4.84, MS-SLSQP 5.38, LX-SSA 5.68 and DE 6.86.
\end{itemize}
In the Holm-adjusted post hoc tests, LX-SSA-VNS ranks significantly better than SSA, PSO, MS-SLSQP, LX-SSA and DE ($p_{\rm Holm}\le2.9\times10^{-6}$). Its differences from VNS ($p_{\rm Holm}=0.123$) and from BVNS ($p_{\rm Holm}=0.068$) are not significant at the case level, whereas VNS ranks significantly better than BVNS ($p_{\rm Holm}=5.1\times10^{-4}$).

The run-level pairwise tests (Table~\ref{tab:wtl}) show where the hybrid gains and loses:
\begin{itemize}
\item \emph{LX-SSA-VNS vs.\ LX-SSA.} The hybrid is significantly better in 35 cases and never worse ($r_{rb}$ from $0.50$ to $0.98$, median $0.77$). The 33 cases without a significant difference are mostly small cases ($N\le7$), in which LX-SSA is already close to the best layouts, plus four cases with $N=8$ or 10, including both 500-m ten-turbine cases, in which many runs of both methods are infeasible.
\item \emph{LX-SSA-VNS vs.\ SSA and PSO.} The hybrid is significantly better in 27 and 32 cases, respectively, and never worse.
\item \emph{LX-SSA-VNS vs.\ DE.} The hybrid is better in 52 cases and worse in only one three-turbine case.
\item \emph{LX-SSA-VNS vs.\ MS-SLSQP.} The hybrid is better in 34 cases and worse in five, all dense 500-m cases (Data Set~I, $N=9,10$; Data Set~II, $N=8$--10), where MS-SLSQP's explicit constraint handling pays off.
\item \emph{LX-SSA-VNS vs.\ VNS.} The two are statistically indistinguishable in 53 cases. VNS is significantly better in 15 cases, all with $N\ge7$ (median $r_{rb}=-0.62$), and the hybrid is never significantly better than VNS.
\item \emph{LX-SSA-VNS vs.\ BVNS.} The hybrid is significantly better in 12 cases, all with $N=5$ to 11, and worse in seven: five cases with $N=3$ or 4, and the two 500-m ten-turbine cases, in which BVNS is feasible in every run.
\end{itemize}
Hybridization therefore lifts LX-SSA from the seventh to the second average rank among the eight methods, but it does not surpass the modified VNS at this budget.

\subsection{Solution Quality by Farm Size and Turbine Count}
""" + fig("wakeloss_vs_n", "fig:wakeloss",
    "Mean wake loss (percentage of the wake-free objective) of the feasible runs versus the prescribed number of turbines; a point is omitted when a method has no feasible run.") + r"""
Figure~\ref{fig:wakeloss} and Tables~\ref{tab:res-1-500}--\ref{tab:res-2-1000} show that the wake loss grows with $N$, falls with the farm radius, and is higher for the broader wind rose of Data Set~II. For $N\le3$ all methods stay within 0.3\% of the wake-free objective. With growing $N$ the curve of the hybrid separates from that of LX-SSA and moves close to the VNS and MS-SLSQP curves, as in these examples:
\begin{itemize}
\item Data Set~I, 1000~m/15: LX-SSA-VNS 2.08\%, VNS 1.82\%, BVNS 1.99\%, MS-SLSQP 2.20\%, SSA 2.80\%, LX-SSA 3.32\%;
\item Data Set~II, 750~m/12: LX-SSA-VNS 5.95\%, VNS 5.76\%, BVNS 6.17\%, MS-SLSQP 5.60\%, SSA 6.79\%, LX-SSA 7.34\%.
\end{itemize}
In every largest case the hybrid reduces the mean wake loss of LX-SSA by 1.2 to 2.0 percentage points, while VNS or MS-SLSQP retains the lowest value.

""" + "\n".join(T[f"tab:res-{d}-{r}"] for d in (1, 2) for r in (500, 750, 1000)) + r"""

\subsection{Feasibility}
""" + fig("feasibility_vs_n", "fig:feasibility",
    "Percentage of the 30 runs that end with a feasible layout (all turbines inside the farm and at least 308~m apart).") + r"""
Figure~\ref{fig:feasibility} shows that LX-SSA-VNS returns a feasible layout in 98.1\% of all runs, compared with 96.9\% for LX-SSA, 97.6\% for SSA, 86.0\% for PSO, 71.4\% for DE, 99.6\% for VNS and 99.9\% for MS-SLSQP. BVNS is the only method that is feasible in every run. The only cases in which the hybrid misses feasibility are the densest 500-m cases: 26 of 30 runs are feasible for $N=9$ and 15 of 30 for $N=10$, compared with 21 and 13 for LX-SSA, 29 and 27 for VNS, and 30 and 30 for BVNS. In these cases the LX-SSA phase frequently ends with an infeasible layout, and the VNS phase, which works on a single incumbent, cannot always repair it in the remaining budget. The feasibility counts of all penalty-based methods are identical for the two wind data sets, because while every candidate is infeasible the penalty term dominates $F_p$.

\subsection{Convergence Behaviour}
""" + fig("convergence_mid", "fig:conv-mid",
    "Convergence for moderate turbine counts: median (line) and interquartile range (band) of the best feasible wake loss found so far over the 30 runs, versus the number of objective calls (logarithmic scale). A curve starts once at least half of the runs have found a feasible layout; the LX-SSA-VNS curve switches from the LX-SSA to the VNS phase at 3,030 calls.") + \
    fig("convergence_max", "fig:conv-max",
    "Convergence for the largest turbine count of each farm (as Fig.~\\ref{fig:conv-mid}). Methods for which fewer than half of the runs find a feasible layout have no curve.") + r"""
Figures~\ref{fig:conv-mid} and~\ref{fig:conv-max} show the effect of the hybridization directly. Until 3,030 calls the hybrid runs a shortened LX-SSA and is therefore slightly behind the full LX-SSA run. After the switch its curve drops steeply, as VNS continues to improve the incumbent that LX-SSA has found, whereas stand-alone LX-SSA and SSA stagnate after about 3,000 calls. Over all feasible runs, the VNS phase removes on average 36.4\% (median 31.1\%) of the wake loss present at the end of the LX-SSA phase, and it turns 82 runs that were still infeasible at the switch into feasible ones. Stand-alone VNS reaches low wake losses earlier because it spends its whole budget on intensification, and MS-SLSQP is the fastest method to reach feasibility. BVNS reaches feasibility in the dense 500-m cases earlier than VNS, but afterwards it improves in fewer and larger steps, because every sweep of its complete local search costs $4N$ calls.

\subsection{Distribution of Results}
""" + fig("boxplots_mid", "fig:box",
    "Distribution of the final wake loss of the feasible runs for moderate turbine counts (30 runs per method).") + r"""
Figure~\ref{fig:box} shows the run-to-run distributions for the moderate cases. LX-SSA-VNS has the lowest median wake loss in both 500-m cases and the second-lowest, after VNS, in the other four, with a smaller interquartile range than LX-SSA in five of the six cases. DE has the highest median in all six cases.

\subsection{Original versus Modified VNS}
\label{sec:bvns}
""" + BV["tab:bvns"] + r"""
The VNS used as a baseline and inside the hybrid departs from the original algorithm in two ways: it perturbs only one to three turbines, and it replaces the complete local search by ten trial moves. Table~\ref{tab:bvns} shows how the original basic VNS (BVNS), with shaking of the whole layout and a complete best-improvement local search, compares with every other method.
\begin{itemize}
\item \emph{BVNS vs.\ the population-based methods and MS-SLSQP.} BVNS is significantly better than LX-SSA in 32 cases, SSA in 24, PSO in 30, DE in 53 and MS-SLSQP in 24, and worse in at most four cases against any of them.
\item \emph{BVNS vs.\ the modified VNS.} BVNS is better in eight cases, all with $N=3$ or 4, and worse in 27 cases, all with $N\ge6$ (median $r_{rb}=-0.78$). For small layouts a complete descent over all coordinate moves is affordable and pays off; for larger layouts every sweep costs $4N$ calls, so few VNS iterations fit into the budget, and the cheaper turbine-wise moves of the modified VNS make better use of it.
\item \emph{BVNS vs.\ LX-SSA-VNS.} BVNS is better in eight cases and worse in 11. These counts differ slightly from Table~\ref{tab:wtl} because the Holm adjustment is applied to a different family of comparisons.
\item \emph{Feasibility.} BVNS is the only method that returns a feasible layout in all 2,040 runs, including the dense 500-m cases.
\end{itemize}
The modifications therefore make VNS stronger on the larger layouts that matter most in practice, which supports using the modified VNS as the intensification phase of the hybrid; the original BVNS nevertheless remains a strong baseline.

\subsection{Ablation of the Hybrid}
\label{sec:ablation}
""" + T["tab:split"] + r"""
Two ablations clarify the role of each phase:
\begin{itemize}
\item \emph{Phase contribution.} Comparing LX-SSA-VNS with LX-SSA (same first phase) and with VNS (same second phase) shows that nearly all of the improvement over LX-SSA comes from the VNS phase: the hybrid beats LX-SSA in 35 cases and never loses. The LX-SSA phase, however, does not provide a better starting point than the best member of a random initial population, which stand-alone VNS uses: stand-alone VNS remains better in 15 dense cases.
\item \emph{Budget split.} Table~\ref{tab:split} varies $\rho$ over 25\%, 50\% and 75\% on twelve representative cases. Giving more budget to LX-SSA is harmful: the 50\% split has a higher mean than the 75\% split in 11 of 12 cases, significantly so in four. Giving more budget to VNS helps slightly: the 25\% split has the higher mean in 9 of 12 cases, but the difference is significant in only one.
\end{itemize}
Both ablations point the same way: at this budget the intensification phase is the valuable part of the hybrid. A larger total budget, or a switch triggered by stagnation of the swarm instead of a fixed split, would be needed for the exploration phase to pay off.

\subsection{Optimized Layouts}
""" + fig("layouts_max", "fig:layouts",
    "Best layouts for the largest turbine count of each farm: best LX-SSA-VNS run (filled markers) and best run of any method (open squares, method in the legend). Circle: farm boundary. Coordinates are listed in Appendix~\\ref{app:coordinates}.") + r"""
Figure~\ref{fig:layouts} compares, for the largest turbine count of each farm, the best LX-SSA-VNS layout with the best layout found by any method, which comes from BVNS in three cases, from MS-SLSQP in two and from VNS in one. Most turbines are placed close to the boundary, which maximizes their mutual distances, while the remaining turbines occupy the interior at positions staggered with respect to the dominant wind directions.

\subsection{Computational Cost}
\label{sec:runtime}
At an equal number of objective calls all methods need almost the same time: 0.37~ms per call for BVNS, 0.39~ms for LX-SSA-VNS, 0.40~ms for LX-SSA and SSA, 0.41~ms for VNS and PSO, 0.47~ms for DE and 0.49~ms for MS-SLSQP. The objective evaluation dominates the cost, so neither the Laplace step nor the VNS bookkeeping adds measurable overhead. A 6,030-call run of LX-SSA-VNS takes 1.7--3.5~s on average, depending on $N$.

\subsection{Horns Rev 1 Site Case}
\label{sec:hornsrev-site}
""" + T["tab:hr-site"] + "\n" + fig("hr_convergence", "fig:hr-conv",
    "Horns Rev~1 site case: median best feasible AEP versus objective calls (16-turbine block: 30 runs; 80-turbine farm: 10 runs); dotted line: installed layout. Only MS-SLSQP finds feasible layouts for the 80-turbine farm.") + \
    fig("hr_layouts", "fig:hr-layouts",
    "Horns Rev~1: installed layout (crosses) and best optimized layout (open markers) inside the installed farm outline; local coordinates in metres.") + r"""
To bring the validation closer to a real project, the Horns Rev~1 offshore farm was modeled with its real turbine, wind climate and site outline, all taken from the PyWake distribution~\cite{PyWake}:
\begin{itemize}
\item the installed layout of 80 Vestas V80 turbines (rotor diameter 80~m, 2~MW) and the V80 power and thrust-coefficient curves;
\item the measured 12-sector wind climate;
\item the parallelogram spanned by the installed turbines as the site outline.
\end{itemize}
The farm model keeps the Jensen top-hat wake (with $k=0.04$, a value commonly used offshore) and root-sum-square superposition. It uses the tabulated power curve with a 25~m/s cut-out, a speed-dependent thrust coefficient at the free-stream speed, 5$^\circ$ direction bins and 1~m/s speed bins. For the installed layout it gives 664.6~GWh/yr and a wake loss of 10.95\%, against 662.5~GWh/yr and 10.96\% from PyWake's Jensen model with the same settings, a difference below 0.3\%. The 16-turbine block that PyWake uses as an example (30 seeds) and the complete 80-turbine farm (10 seeds) were optimized with all eight methods at 6,030 calls and the $4D$ (320~m) minimum spacing.

Table~\ref{tab:hr-site} and Figs.~\ref{fig:hr-conv} and~\ref{fig:hr-layouts} lead to three observations.
\begin{itemize}
\item \emph{The benchmark ordering carries over to the site.} In the 16-turbine block VNS attains the highest mean AEP (137.97~GWh/yr), followed by BVNS (137.67~GWh/yr, feasible in all 30 runs) and LX-SSA-VNS (137.59~GWh/yr). The hybrid is significantly better than LX-SSA ($p_{\rm Holm}=0.003$) and MS-SLSQP ($p_{\rm Holm}=0.040$); its differences from SSA and VNS ($p_{\rm Holm}=0.059$ each) and from BVNS ($p_{\rm Holm}=0.191$) are not significant. PSO and DE never find a feasible layout inside the parallelogram.
\item \emph{No method improves on the installed layout within the budget.} The installed $7D$ grid reaches 139.51~GWh/yr.
\item \emph{The full farm exceeds what the penalty-based methods can handle at this budget.} With 160 variables and 3,160 spacing constraints, none of the runs of LX-SSA-VNS, LX-SSA, SSA, PSO, DE, VNS or BVNS returns a feasible layout. MS-SLSQP, which handles the constraints explicitly, is feasible in all ten runs, but its mean AEP (653.9~GWh/yr) remains 1.6\% below the installed layout.
\end{itemize}
Realistic farm sizes therefore require larger budgets, feasible initialization or explicit constraint handling, and a regular installed layout is a strong reference.

"""

sensitivity = r"""\section{Sensitivity to the Modelling Assumptions}
\label{sec:robustness}
\subsection{Power Curve and Wake Model}
\label{sec:powercurve}
\label{sec:gaussian}
""" + RB["tab:robust-hybrid"] + r"""
To test how far the comparison depends on the benchmark power curve and wake model, every feasible final layout was re-evaluated, without re-optimization, with three alternative models:
\begin{itemize}
\item a cubic power curve~\cite{Carrillo2013}, $f_{\rm c}(s)=P_{\text{rated}}(s^3-s_{\text{cut-in}}^3)/(s_{\text{rated}}^3-s_{\text{cut-in}}^3)$ between cut-in and rated speed;
\item the same curve with a 25~m/s cut-out;
\item the Gaussian wake model of Bastankhah and Port\'e-Agel~\cite{Bastankhah2014},
\begin{align*}
\delta_{ij}&=\Bigl(1-\sqrt{1-C_T/(8\sigma^2/D^2)}\Bigr)\exp\!\bigl(-r_\perp^2/(2\sigma^2)\bigr),\\
\sigma/D&=k^*x/D+0.2\sqrt{\beta},
\end{align*}
with $\beta=\tfrac12(1+\sqrt{1-C_T})/\sqrt{1-C_T}$ and $k^*=0.04$, evaluated at the rotor center with the same superposition and Weibull scaling as the benchmark.
\end{itemize}

Table~\ref{tab:robust-hybrid} shows two kinds of sensitivity.
\begin{itemize}
\item \emph{Absolute energy is strongly model-dependent.} The linear benchmark curve overstates expected power relative to the cubic curve by 17.7\% for Data Set~I and 36.9\% for Data Set~II (21.5\% and 37.1\% with the cut-out). The absolute energy values of this paper are therefore benchmark values, not turbine-specific AEP predictions.
\item \emph{The algorithm comparison is robust to the power curve but less so to the wake model.} Under the cubic curve the within-case orderings are almost unchanged (mean Kendall $\bar\tau=0.96$), and VNS (2.13), LX-SSA-VNS (2.78) and BVNS (3.64) remain the three best methods. Under the Gaussian model the orderings change more ($\bar\tau=0.57$; same best method in 46\% of the cases), because layouts optimized for the sharp Jensen wake cone are penalized differently by the smooth Gaussian deficit. VNS (2.75) and LX-SSA-VNS (3.25) nevertheless keep the two best average ranks, whereas BVNS falls to 4.62, behind PSO (3.81) and MS-SLSQP (4.15).
\end{itemize}

\subsection{Minimum Spacing}
\label{sec:spacing}
""" + CAP["tab:capacity"] + "\n" + SP["tab:spacing-authors"] + r"""
The $4D$ spacing is a legacy benchmark value. Its influence was examined in three ways.

\emph{Geometric capacity.} Table~\ref{tab:capacity} gives the largest turbine count for which a multi-start packing search constructed a feasible layout; these are constructive lower bounds, and the 500-m values coincide with the known optimal circle-packing results. At $5D$ the 500-m farm holds at most 8 turbines and at $6D$ only 7, so the tested 500-m cases with $N\ge9$ (at $5D$) or $N\ge8$ (at $6D$) would be infeasible.

\emph{Final layouts.} The optimized layouts are frequently active at the $4D$ constraint. Of the 11,247 feasible final layouts of LX-SSA, SSA, PSO, DE, VNS and MS-SLSQP, 48.7\% also satisfy $5D$ and 35.7\% satisfy $6D$. For $N=6$--10 these fractions fall to 29.6\% and 11.3\%, and for $N=11$--15 to 3.1\% and 0\%. Optimized layouts therefore cannot be transferred to a site with larger spacing requirements without re-optimization.

\emph{Re-optimization.} Six representative cases were re-optimized with LX-SSA and SSA at $5D$ and $6D$ (Table~\ref{tab:spacing-authors}). In the 750-m eight-turbine cases the mean LX-SSA objective falls from $4D$ to $6D$ by 0.95\% (Data Set~I) and 0.22\% (Data Set~II), and the wake loss rises from 0.88\% to 1.85\% and from 3.18\% to 3.47\%. The spacing rule does not reverse the LX-SSA--SSA comparison.

\subsection{Boundary Handling}
\label{sec:boundary}
All methods clip coordinates to the bounding square and penalize points outside the circle (Section~\ref{sec:constraints}). An otherwise identical SSA that instead projects every turbine radially onto the circle attains higher mean objectives than the original SSA in six representative cases at 3,030 calls, by 26 to 649 objective units (Mann--Whitney $p<10^{-4}$ in every case). The boundary-handling rule therefore affects results noticeably and is held fixed for all methods.

"""

complexity = r"""\section{Computational Complexity}
Let $N_t$ be the number of turbines, $N_\theta$ the number of wind-direction bins and $N_s$ the number of wind-speed bins. A direct wake evaluation contains the pairwise turbine interactions for every direction and speed bin, so one objective evaluation costs
\[
C_{\mathrm{obj}}=\mathcal{O}\!\left(N_t^2N_\theta N_s\right),
\]
and a run with a budget of $B$ objective calls costs $\mathcal{O}(B\,N_t^2N_\theta N_s)$ for every method considered here, including both phases of LX-SSA-VNS. Box repair is $\mathcal{O}(N_t)$ per candidate and the spacing check $\mathcal{O}(N_t^2)$; the Laplace draw of LX-SSA and the shaking and local-search moves of VNS are $\mathcal{O}(N_t)$. All of these are dominated by the wake evaluation, which is consistent with the almost identical measured cost per call (Section~\ref{sec:runtime}).

"""

limitations = r"""\section{Scope, Limitations, and Validity of the Evidence}
The results should be interpreted within the following boundaries:
\begin{itemize}
\item \textbf{Relative performance of the hybrid:} LX-SSA-VNS improves substantially on LX-SSA, SSA, PSO, DE and MS-SLSQP, but it does not outperform the modified VNS at the tested budget, and at the case level its advantage over the original BVNS is not significant. The ablation indicates that its advantage over LX-SSA comes from the VNS phase.
\item \textbf{Budget and farm size:} all methods were compared at 6,030 objective calls; larger budgets could change the ordering, in particular the value of the exploration phase. For an 80-turbine farm the penalty-based methods, including the hybrid, did not find feasible layouts at this budget.
\item \textbf{Algorithm settings:} the LX-SSA parameters ($\phi=0$, $\chi=1$), the VNS and BVNS settings and the split $\rho=0.5$ are fixed defaults and were not tuned. Pseudo-gradient and surrogate-assisted methods are discussed but not benchmarked.
\item \textbf{Wake model and power curve:} Jensen's model and the piecewise-linear power curve are benchmark choices. Absolute energy values depend strongly on them, and the within-case orderings change under a Gaussian wake model (Section~\ref{sec:gaussian}); layouts were not re-optimized with the alternative models.
\item \textbf{Spacing and boundary handling:} $4D$ is a legacy constraint that governs the feasible turbine density, and the boundary-handling rule also affects absolute results.
\item \textbf{Scope:} $N$ is prescribed in every run. Terrain, heterogeneous turbines, electrical collector systems, grid, acoustic and environmental constraints, and economic objectives are not modeled.
\end{itemize}

"""

conclusion = r"""\section{Conclusion}
This paper proposes LX-SSA-VNS, a two-phase hybrid in which the previously published Laplacian Salp Swarm Algorithm (LX-SSA)~\cite{Solanki2023} explores the continuous layout space and variable neighbourhood search then intensifies the best layout found. The hybrid is applied to the continuous WFLOP with a Jensen--Weibull benchmark model and explicit constraint handling. It is evaluated against LX-SSA, SSA, PSO, DE, the original basic VNS (BVNS), a modified VNS and a gradient-based multistart SLSQP solver on 68 benchmark cases (16,320 seed-paired runs at equal budgets) and on the Horns Rev~1 site.

The results support five conclusions.
\begin{itemize}
\item \emph{Hybridization substantially improves LX-SSA.} LX-SSA-VNS is significantly better than LX-SSA in 35 of 68 cases and never worse. It lifts LX-SSA from the seventh to the second average rank among eight methods (5.68 to 2.74) and achieves a feasibility rate of 98.1\%.
\item \emph{The hybrid outperforms the other population-based and gradient-based references.} It is significantly better than SSA, PSO, DE and MS-SLSQP in most cases, and at the case level it ranks significantly better than all of them.
\item \emph{The modified VNS remains the best method at this budget.} The hybrid is statistically indistinguishable from it in 53 cases but worse in 15 dense cases. The ablation shows that the gain comes from the intensification phase, and that assigning more budget to LX-SSA is harmful.
\item \emph{The original basic VNS is a strong baseline.} BVNS ranks third (3.62), is feasible in every run, and is better than the hybrid in seven cases (small or very dense layouts) and worse in 12. The modified VNS beats it in 27 cases with $N\ge6$, while BVNS is better in eight cases with $N\le4$, so the cheaper turbine-wise search pays off for larger layouts.
\item \emph{A similar ordering holds at a real site.} In a 16-turbine block of Horns Rev~1 the hybrid is better than LX-SSA and MS-SLSQP and statistically indistinguishable from SSA, VNS and BVNS. No method improves on the installed layout within the budget, and for the full 80-turbine farm only the constraint-aware MS-SLSQP finds feasible layouts.
\item \emph{The comparison depends on the wake model more than on the power curve.} The algorithm ordering is insensitive to the power curve but changes considerably under a Gaussian wake model, and absolute energy values depend strongly on both.
\end{itemize}

Future work should use a stagnation-triggered switch between the two phases, feasibility-preserving initialization for dense farms, larger budgets for realistic farm sizes, and re-optimization under higher-fidelity wake and power models.

"""


def appendix():
    B = pd.read_csv("hybrid_best_layouts_maxN.csv")
    B = B.sort_values("Objective").groupby(["Dataset", "Radius", "Turbines"]).tail(1).sort_values(["Dataset", "Radius"])
    lab = {"LXVNS": "LX-SSA-VNS", "LXSSA": "LX-SSA", "SSA": "SSA", "PSO": "PSO", "DE": "DE", "VNS": "VNS", "SLSQP": "MS-SLSQP", "BVNS": "BVNS"}
    rows = []
    for _, r in B.iterrows():
        pts = [tuple(map(float, q.split())) for q in r.Coordinates.split(";")]
        cells = [f"({x:.1f}, {y:.1f})" for x, y in pts]
        lines = ["; ".join(cells[k:k + 4]) for k in range(0, len(cells), 4)]
        rows.append(f"{r.Dataset} & {r.Radius} & {r.Turbines} & {lab[r.Algorithm]} & {r.Objective:.1f} & "
                    + " \\\\\n & & & & & ".join(lines) + " \\\\")
    return r"""\appendices
\section{Coordinates of the Best Layouts}
\label{app:coordinates}
Table~\ref{tab:coords} lists the turbine coordinates (m, farm center at the origin) of the best layout found by any method for the largest turbine count of each farm (Fig.~\ref{fig:layouts}). The final coordinates and convergence curves of all 16,320 runs and the Horns Rev~1 layouts are provided as machine-readable files with the revision.

\begin{table*}[!t]
\centering
\caption{Coordinates of the best layout for the largest turbine count of each farm (method and benchmark objective given).}
\label{tab:coords}
\scriptsize\setlength{\tabcolsep}{3pt}
\begin{tabular}{ccccc l}
\toprule
DS & $r$ (m) & $N$ & Method & Objective & Coordinates $(\rho_i,\sigma_i)$ (m) \\
\midrule
""" + "\n\\midrule\n".join(rows) + r"""
\bottomrule
\end{tabular}
\end{table*}

"""


def main():
    t = open(MS).read()
    start = r"\section{Hybrid LX-SSA-VNS Algorithm}"
    i = t.index(start) if start in t else t.index(r"\section{Experimental Setup}")
    j = t.index(r"\appendices")
    t = t[:i] + open("sec_hybrid.tex").read() + setup_text() + results + sensitivity + complexity + limitations + conclusion + t[j:]
    i = t.index(r"\appendices")
    j = t.index("\\vspace{1em}\n\\noindent\\textbf{Compliance with Ethical Standards}")
    t = t[:i] + appendix() + t[j:]
    open(MS, "w").write(t)
    print("manuscript rebuilt (hybrid)")


if __name__ == "__main__":
    main()

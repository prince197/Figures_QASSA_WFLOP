"""Replace all earlier numerical results in the manuscript by the fresh study.

Inputs: sec_setup.tex, fresh_tables.tex, fresh_robust_table.tex, robustness_tables.tex
(capacity table), authors_tables.tex (spacing table), fresh_best_layouts_maxN.csv, figures in
../figures_fresh. Every number quoted below was checked against fresh_*.csv / fresh_summary.json.
"""
import pandas as pd

MS = "../LXSSA_WFLOP_reviewer_revised.tex"


def blocks(path):
    txt = open(path).read()
    return {b.split("\\label{")[1].split("}")[0]: ("\\begin{table" + b).rstrip()
            for b in txt.split("\\begin{table") if "\\label{" in b}


def fig(name, label, caption, wide=True, width="\\textwidth"):
    env = "figure*" if wide else "figure"
    return (f"\\begin{{{env}}}[!t]\n\\centering\n\\includegraphics[width={width}]{{figures_fresh/{name}.pdf}}\n"
            f"\\caption{{{caption}}}\n\\label{{{label}}}\n\\end{{{env}}}\n")


T = blocks("fresh_tables.tex")
RB = blocks("fresh_robust_table.tex")
CAP = blocks("robustness_tables.tex")
SP = blocks("authors_tables.tex")

results = r"""\section{Numerical Results}
\label{sec:results}
\label{sec:equalbudget}
This section reports the results of 12,240 runs: 68 test cases, six methods and 30 seed-paired runs, all at 6,030 objective calls. It closes with the Horns Rev~1 site case. Detailed per-case results are given in Tables~\ref{tab:res-1-500}--\ref{tab:res-2-1000}.

\subsection{Overall Comparison}
\label{sec:overall}
""" + T["tab:friedman68"] + "\n" + fig("avg_ranks", "fig:avgranks",
    "Average rank of the six methods over the 68 benchmark cases (ranks of the mean feasible objective within each case; 1 = best).",
    wide=False, width="0.95\\columnwidth") + T["tab:wtl"] + r"""
The case-level Friedman test over the 68 test cases rejects the hypothesis that the six methods perform equally ($\chi^2_F=135.3$, 5 d.f., $p=1.8\times10^{-27}$; Iman--Davenport $F_F=44.3$). VNS is clearly the strongest method (Table~\ref{tab:friedman68}, Fig.~\ref{fig:avgranks}). It has an average rank of 1.60 and alone attains the highest mean objective in 44 of the 68 cases. It is followed by SSA (3.15), PSO (3.33), MS-SLSQP (3.87), LX-SSA (3.97) and DE (5.08). The Holm-adjusted post hoc tests find VNS significantly better than every other method ($p_{\rm Holm}\le1.5\times10^{-6}$). The average rank of LX-SSA is significantly worse than that of SSA ($p_{\rm Holm}=0.031$) and significantly better than that of DE ($p_{\rm Holm}=0.002$); its differences from PSO and MS-SLSQP are not significant. LX-SSA never alone attains the highest mean.

The run-level pairwise tests (Table~\ref{tab:wtl}) give a consistent picture:
\begin{itemize}
\item \emph{LX-SSA vs.\ SSA (mechanism-level ablation).} LX-SSA is never significantly better. It is statistically indistinguishable in 59 cases and significantly worse in nine ($r_{rb}$ between $-0.48$ and $-0.65$), which are mostly dense cases: Data Set~I at 750~m/12 and 1000~m/14--15, and Data Set~II at 500~m/7 and 9, 750~m/7, 11 and 12, and 1000~m/14. At equal cost, the Laplace follower candidate therefore brings no benefit over the SSA follower update on this benchmark. The likely reason is that it spends two objective calls per follower and iteration, which halves the number of iterations.
\item \emph{LX-SSA vs.\ VNS.} VNS is significantly better in 41 cases, with large effects (median $r_{rb}=-0.87$). The 27 cases without a significant difference all have $N\le7$ and are mostly the small cases, in which all methods come close to the wake-free optimum.
\item \emph{LX-SSA vs.\ PSO} (18 better, 39 not different, 11 worse). LX-SSA wins mainly in the dense cases, where PSO often fails to find a feasible layout. It loses in small and medium cases ($N=3$--7), where PSO converges to better layouts.
\item \emph{LX-SSA vs.\ MS-SLSQP} (18/31/19). LX-SSA is better at intermediate turbine counts, mostly for the 1000-m farm, and MS-SLSQP is better in the densest cases ($N\ge7$ at 500~m, $N\ge10$ at 750~m and $N\ge13$ at 1000~m).
\item \emph{LX-SSA vs.\ DE} (46/18/4). LX-SSA is better in most cases; the four losses are all three-turbine cases.
\end{itemize}

\subsection{Solution Quality by Farm Size and Turbine Count}
""" + fig("wakeloss_vs_n", "fig:wakeloss",
    "Mean wake loss (percentage of the wake-free objective) of the feasible runs versus the prescribed number of turbines; a point is omitted when a method has no feasible run.") + r"""
Figure~\ref{fig:wakeloss} and Tables~\ref{tab:res-1-500}--\ref{tab:res-2-1000} show how the solution quality depends on the farm radius and the turbine count. For $N\le3$ all methods stay within 0.3\% of the wake-free objective. The wake loss then grows with $N$ and falls with the farm radius, and it is higher for Data Set~II, whose wind comes from a broader range of directions, than for the strongly directional Data Set~I. The differences between the methods grow with $N$. In the largest cases the lowest mean wake loss is obtained by VNS or MS-SLSQP; for example:
\begin{itemize}
\item Data Set~I, 1000~m/15: VNS 1.82\%, MS-SLSQP 2.20\%, SSA 2.80\%, LX-SSA 3.32\%;
\item Data Set~I, 500~m/10: MS-SLSQP 3.56\%, VNS 4.54\%, SSA 5.58\%, LX-SSA 6.76\%;
\item Data Set~II, 1000~m/15: VNS 5.19\%, MS-SLSQP 5.36\%, SSA 6.36\%, LX-SSA 6.75\%.
\end{itemize}
In every largest case, LX-SSA has a higher mean wake loss than SSA.

""" + "\n".join(T[f"tab:res-{d}-{r}"] for d in (1, 2) for r in (500, 750, 1000)) + r"""

\subsection{Feasibility}
""" + fig("feasibility_vs_n", "fig:feasibility",
    "Percentage of the 30 runs that end with a feasible layout (all turbines inside the farm and at least 308~m apart).") + r"""
Figure~\ref{fig:feasibility} reports how often each method ends with a feasible layout. Up to $N=6$ (500~m), $N=8$ (750~m) and $N=10$ (1000~m) all methods are always feasible. Beyond these counts the methods separate clearly:
\begin{itemize}
\item \emph{DE and PSO} lose feasibility first. DE finds no feasible layout at 500~m with $N=9$, at 750~m with $N=12$ or at 1000~m with $N\ge14$. PSO succeeds in 2 of 30 runs at 500~m/9 and in none at 500~m/10.
\item \emph{LX-SSA and SSA} remain feasible much longer. They drop to 13 and 12 feasible runs at 500~m/10 and 21 and 25 at 500~m/9, but stay at 28--30 in all 750-m and 1000-m cases.
\item \emph{VNS and MS-SLSQP} are the most reliable, with at least 27 feasible runs in every case.
\end{itemize}
For LX-SSA, SSA, PSO, DE and VNS the feasibility counts are identical for the two wind data sets. The constraints do not depend on the wind data, and while every candidate is infeasible the penalty term dominates $F_p$, so the search follows the same path until a feasible layout is found. MS-SLSQP, which uses the wind-dependent wake loss from the start, differs slightly between the data sets.

\subsection{Convergence Behaviour}
""" + fig("convergence_mid", "fig:conv-mid",
    "Convergence for moderate turbine counts: median (line) and interquartile range (band) over the 30 runs of the best feasible wake loss found so far, versus the number of objective calls (logarithmic scale). A curve starts once at least half of the runs have found a feasible layout.") + \
    fig("convergence_max", "fig:conv-max",
    "Convergence for the largest turbine count of each farm (as Fig.~\\ref{fig:conv-mid}). Methods for which fewer than half of the runs find a feasible layout within the budget have no curve (PSO and DE, and LX-SSA and SSA at 500~m/10).") + r"""
Figures~\ref{fig:conv-mid} and~\ref{fig:conv-max} show the median best feasible wake loss as a function of the number of objective calls.
\begin{itemize}
\item \emph{LX-SSA and SSA} converge early and then stagnate. In the six moderate cases they reach within 0.1\% of their final objective after a median of 2,900--3,800 calls, so the second half of the budget brings little improvement.
\item \emph{PSO and DE} converge more slowly: PSO needs 4,300--5,100 calls and DE 3,600--5,600.
\item \emph{VNS} needs 2,900--4,600 calls and still improves late in the budget, which is consistent with its local-search refinement.
\item \emph{MS-SLSQP} finds feasible layouts earliest. In the largest cases its median first feasible layout appears after 150--690 calls, compared with 1,500--2,900 for VNS and 2,800--3,700 for LX-SSA and SSA.
\end{itemize}
In the largest cases the median LX-SSA curve lies above the SSA curve for most of the budget.

\subsection{Distribution of Results}
""" + fig("boxplots_mid", "fig:box",
    "Distribution of the final wake loss of the feasible runs for moderate turbine counts (30 runs per method).") + r"""
Figure~\ref{fig:box} shows the distribution of the final wake loss for the moderate cases. VNS has the lowest median in five of the six cases (PSO has it in Data Set~I at 500~m/6), and DE has the highest median in all six. The LX-SSA and SSA distributions overlap strongly: LX-SSA has the lower median in three cases and SSA in the other three, in line with the mostly non-significant paired tests between them.

\subsection{Optimized Layouts}
""" + fig("layouts_max", "fig:layouts",
    "Best layouts for the largest turbine count of each farm: best LX-SSA run (filled circles) and best run of any method (open squares, method in the legend). Circle: farm boundary. Coordinates are listed in Appendix~\\ref{app:coordinates}.") + r"""
Figure~\ref{fig:layouts} compares, for the largest turbine count of each farm, the best LX-SSA layout with the best layout found by any method; the best overall layout comes from VNS or MS-SLSQP in all six cases. Most turbines are placed close to the boundary, which maximizes their mutual distances, while the remaining turbines occupy the interior at positions staggered with respect to the dominant wind directions.

\subsection{Computational Cost}
\label{sec:runtime}
At an equal number of objective calls the methods need almost the same time: 0.40~ms per call for LX-SSA and SSA, 0.41~ms for VNS and PSO, 0.47~ms for DE and 0.49~ms for MS-SLSQP, whose quadratic subproblems add a small overhead. A 6,030-call run of LX-SSA takes 1.7--3.6~s on average, depending on $N$; the slowest single run of any method (MS-SLSQP) took 6.3~s. The objective evaluation dominates the cost, so the Laplace step of LX-SSA adds no measurable overhead. The complete study of 12,240 runs took 2.2~h on four cores.

\subsection{Horns Rev 1 Site Case}
\label{sec:hornsrev-site}
""" + T["tab:hr-site"] + "\n" + fig("hr_convergence", "fig:hr-conv",
    "Horns Rev~1 site case: median best feasible AEP versus objective calls (16-turbine block: 30 runs; 80-turbine farm: 10 runs); dotted line: installed layout. Only MS-SLSQP finds feasible layouts for the 80-turbine farm.") + \
    fig("hr_layouts", "fig:hr-layouts",
    "Horns Rev~1: installed layout (crosses) and best optimized layout (open markers) inside the installed farm outline; local coordinates in metres.") + r"""
To bring the validation closer to a real project, the Horns Rev~1 offshore farm was modeled with its real turbine, wind climate and site outline, all taken from the PyWake distribution~\cite{PyWake}:
\begin{itemize}
\item the installed layout of 80 Vestas V80 turbines (rotor diameter 80~m, 2~MW) and the V80 power and thrust-coefficient curves;
\item the measured 12-sector wind climate (Weibull parameters and frequency per sector);
\item the parallelogram spanned by the installed turbines as the site outline.
\end{itemize}
The farm model keeps the Jensen top-hat wake (with $k=0.04$, a value commonly used offshore) and root-sum-square superposition. It uses the tabulated power curve with a 25~m/s cut-out, a speed-dependent thrust coefficient evaluated at the free-stream speed, 5$^\circ$ direction bins and 1~m/s speed bins. As an independent check, this model gives the installed layout an AEP of 664.6~GWh/yr with a wake loss of 10.95\%, compared with 662.5~GWh/yr and 10.96\% from PyWake's Jensen model with the same settings, a difference below 0.3\%. Two problems were solved with the same six methods, 6,030 calls per run and the $4D$ (320~m) minimum spacing: the 16-turbine block that PyWake uses as an example (30 seeds) and the complete 80-turbine farm (10 seeds, because one run takes about 100~s).

Table~\ref{tab:hr-site} and Figs.~\ref{fig:hr-conv} and~\ref{fig:hr-layouts} lead to three observations.
\begin{itemize}
\item \emph{No method improves on the installed layout within the budget.} In the 16-turbine block VNS reaches 137.97~GWh/yr on average and LX-SSA 136.46~GWh/yr, compared with 139.51~GWh/yr for the installed $7D$ grid.
\item \emph{The algorithm ordering matches the benchmark results.} VNS is significantly better than LX-SSA ($p_{\rm Holm}=1.7\times10^{-4}$); LX-SSA does not differ significantly from SSA ($p_{\rm Holm}=0.26$) or MS-SLSQP ($p_{\rm Holm}=0.085$); and PSO and DE never find a feasible layout inside the parallelogram. The convergence curves show the same behaviour as for the benchmark: VNS improves steadily, whereas LX-SSA and SSA stagnate after about 4,000 calls.
\item \emph{The full farm exceeds what the penalty-based methods can handle at this budget.} With 160 variables and 3,160 spacing constraints, none of the 50 runs of LX-SSA, SSA, PSO, DE and VNS returns a feasible layout. MS-SLSQP, which handles the constraints explicitly, is feasible in all ten runs, but its mean AEP (653.9~GWh/yr) remains 1.6\% below that of the installed layout.
\end{itemize}
The benchmark-scale conclusions therefore do not transfer directly to realistic farm sizes: such farms require larger budgets, feasible initialization or explicit constraint handling, and a regular installed layout is a strong reference that none of the tested methods beats within 6,030 calls.

"""

sensitivity = r"""\section{Sensitivity to the Modelling Assumptions}
\label{sec:robustness}
\subsection{Power Curve and Wake Model}
\label{sec:powercurve}
\label{sec:gaussian}
""" + RB["tab:robust-fresh"] + r"""
The benchmark uses a piecewise-linear power curve and Jensen's wake model. To test how far the comparison depends on these choices, every feasible final layout of the study was re-evaluated with three alternative models; the layouts were not re-optimized:
\begin{itemize}
\item a cubic power curve~\cite{Carrillo2013},
\[
f_{\rm c}(s)=P_{\text{rated}}\,\frac{s^3-s_{\text{cut-in}}^3}{s_{\text{rated}}^3-s_{\text{cut-in}}^3}\quad(s_{\text{cut-in}}\le s\le s_{\text{rated}}),
\]
with the same cut-in speed, rated speed and rated power;
\item the same cubic curve with a 25~m/s cut-out;
\item the Gaussian wake model of Bastankhah and Port\'e-Agel~\cite{Bastankhah2014},
\begin{align*}
\delta_{ij}&=\Bigl(1-\sqrt{1-C_T/(8\sigma^2/D^2)}\Bigr)\exp\!\bigl(-r_\perp^2/(2\sigma^2)\bigr),\\
\sigma/D&=k^*x/D+0.2\sqrt{\beta},
\end{align*}
with $\beta=\tfrac12(1+\sqrt{1-C_T})/\sqrt{1-C_T}$ and $k^*=0.04$, evaluated at the rotor center with the same superposition and Weibull scaling as the benchmark.
\end{itemize}

Table~\ref{tab:robust-fresh} shows two different kinds of sensitivity.
\begin{itemize}
\item \emph{Absolute energy is strongly model-dependent.} The linear benchmark curve overstates expected power relative to the cubic curve by 17.7\% for Data Set~I and 36.9\% for Data Set~II (21.5\% and 37.1\% with the cut-out). The absolute energy values of this paper are therefore benchmark values, not turbine-specific AEP predictions.
\item \emph{The algorithm comparison is robust to the power curve but not to the wake model.} Under the cubic curve the ordering of the methods within a case is almost unchanged (mean Kendall $\bar\tau=0.97$; same best method in 96\% of the cases). Under the Gaussian wake model the within-case orderings change substantially ($\bar\tau=0.58$; same best method in only 49\% of the cases), because layouts optimized for the sharp Jensen wake cone are penalized differently by the smooth Gaussian deficit. Nevertheless VNS keeps the best average rank (2.22), and LX-SSA remains fifth (4.19). Whether the ordering persists when layouts are \emph{optimized} with a higher-fidelity wake model remains open.
\end{itemize}

\subsection{Minimum Spacing}
\label{sec:spacing}
""" + CAP["tab:capacity"] + "\n" + SP["tab:spacing-authors"] + r"""
The $4D$ spacing is a legacy benchmark value. Its influence was examined in three ways.

\emph{Geometric capacity.} Table~\ref{tab:capacity} reports the largest turbine count for which a multi-start packing search constructed a layout satisfying both the circular boundary and the minimum distance. These are constructive lower bounds, and the 500-m values coincide with the known optimal circle-packing results. Under $4D$, all tested turbine counts are well below capacity. At $5D$ the 500-m farm holds at most 8 turbines and at $6D$ only 7, so the tested 500-m cases with $N\ge9$ (at $5D$) or $N\ge8$ (at $6D$) would be infeasible. All tested 750-m and 1000-m cases remain geometrically feasible up to $6D$.

\emph{Final layouts.} The optimized layouts are frequently active at the $4D$ constraint. Of the 11,247 feasible final layouts of the study, 48.7\% also satisfy $5D$ and 35.7\% satisfy $6D$. For $N=6$--10 these fractions fall to 29.6\% and 11.3\%, and for $N=11$--15 to 3.1\% and 0\%. Optimized layouts therefore cannot be transferred to a site with larger spacing requirements without re-optimization.

\emph{Re-optimization.} Six representative cases were re-optimized with LX-SSA and SSA at $5D$ and $6D$ (6,030 calls, 30 seed-paired runs). Table~\ref{tab:spacing-authors} shows that larger spacing has little effect in the four-turbine 500-m cases. In the 750-m eight-turbine cases, however, the mean LX-SSA objective falls from $4D$ to $6D$ by 0.95\% (Data Set~I) and 0.22\% (Data Set~II), and the wake loss rises from 0.88\% to 1.85\% and from 3.18\% to 3.47\%. The spacing rule does not reverse the LX-SSA--SSA comparison: the only significant differences (at $6D$) favor SSA.

\subsection{Boundary Handling}
\label{sec:boundary}
The implementation clips coordinates to the bounding square and penalizes points outside the circle (Section~\ref{sec:constraints}). An otherwise identical SSA that instead projects every turbine radially onto the circle, and ranks candidates by the feasibility-first rules, attains higher mean objectives than the original SSA in six representative cases at 3,030 calls. The gains range from 26 to 649 objective units (Mann--Whitney $p<10^{-4}$ in every case). The boundary-handling rule therefore affects results by more than the differences between the salp-swarm algorithms, and it is held fixed for all methods in this study.

"""

complexity = r"""\section{Computational Complexity}
Let $N_t$ be the number of turbines, $N_p$ the population size, $T$ the number of iterations, $N_\theta$ the number of wind-direction bins and $N_s$ the number of wind-speed bins. A direct wake evaluation contains the pairwise turbine interactions for every direction and speed bin, so the dominant cost of one objective evaluation is
\[
C_{\mathrm{obj}}=\mathcal{O}\!\left(N_t^2N_\theta N_s\right),
\]
and a population-based run costs $\mathcal{O}(TN_pN_t^2N_\theta N_s)$, or equivalently $\mathcal{O}(B\,N_t^2N_\theta N_s)$ for a budget of $B$ objective calls. The boundary repair is $\mathcal{O}(N_t)$ per candidate and the spacing check $\mathcal{O}(N_t^2)$; both are dominated by the wake evaluation. The Laplace draw and candidate update of LX-SSA are linear in the decision dimension. LX-SSA does, however, spend two objective calls per follower and iteration, so at a fixed budget it performs half as many iterations as SSA. The measured costs in Section~\ref{sec:runtime} confirm that time per call is almost identical across methods and grows with $N_t$.

"""

limitations = r"""\section{Scope, Limitations, and Validity of the Evidence}
The results should be interpreted within the following boundaries:
\begin{itemize}
\item \textbf{Wake-model fidelity:} Jensen's wake model is retained for benchmark comparability and computational efficiency. Gaussian, CFD-based and data-driven wake models represent wake physics more accurately~\cite{Tao2020,Yang2023,Wang2024}. Re-evaluation with a Gaussian model changes the within-case orderings substantially (Section~\ref{sec:gaussian}); layouts were not re-optimized with it.
\item \textbf{Power curve:} the piecewise-linear benchmark curve overstates expected power by about 18\% (Data Set~I) and 37\% (Data Set~II) relative to a cubic curve, so absolute energy values are benchmark values, not AEP predictions. The algorithm ordering is insensitive to this choice.
\item \textbf{Spacing and boundary handling:} $4D$ is a legacy benchmark constraint; larger spacing reduces the constructible turbine count and changes the optimized layouts (Section~\ref{sec:spacing}). The boundary-handling rule also affects absolute results and is held fixed across methods.
\item \textbf{Budget and farm size:} all methods were compared at 6,030 objective calls. Larger budgets could change the ordering. At realistic farm sizes (80 turbines) the penalty-based methods did not find feasible layouts at this budget.
\item \textbf{Algorithm settings:} the LX-SSA parameters ($\phi=0$, $\chi=1$) and the settings of the other methods are fixed defaults and were not tuned. Pseudo-gradient and surrogate-assisted methods are discussed but not benchmarked.
\item \textbf{Fixed turbine count and engineering scope:} $N$ is prescribed in every run. Terrain, heterogeneous turbines, electrical collector systems, grid, acoustic and environmental constraints, and economic objectives are not modeled.
\end{itemize}

"""

conclusion = r"""\section{Conclusion}
This work applies the previously published Laplacian Salp Swarm Algorithm (LX-SSA)~\cite{Solanki2023} to continuous WFLOP. LX-SSA was developed in 2023 and is not proposed, redesigned or claimed as a new algorithm here. The contribution is an application study: the continuous coordinate encoding, the wake-coupled objective with explicit constraint handling, and a controlled, statistically analyzed comparison with five alternative methods on the standard Jensen--Weibull benchmark and on the Horns Rev~1 site.

In 68 benchmark cases and 12,240 seed-paired runs at equal budgets, the comparison gives five main findings:
\begin{itemize}
\item VNS is the best method, with an average rank of 1.60 and the highest mean in 44 cases.
\item SSA, PSO and MS-SLSQP follow, and LX-SSA ranks fifth (3.97), ahead of DE.
\item In the mechanism-level comparison with SSA, LX-SSA is never significantly better and is significantly worse in nine cases, so the Laplace follower candidate does not provide an equal-budget improvement on this benchmark.
\item LX-SSA is clearly better than DE and more reliable than PSO and DE in reaching feasible layouts for dense farms, but less reliable than VNS and MS-SLSQP.
\item At Horns Rev~1, with its real turbine, wind climate and outline, the same ordering appears. No method improves on the installed layout within the budget, and only the constraint-aware MS-SLSQP finds feasible layouts for the full 80-turbine farm.
\end{itemize}
Under the tested settings, LX-SSA is therefore not the preferred optimizer for this problem.

The algorithm ordering is insensitive to the power-curve model but changes considerably under a Gaussian wake model, and the absolute energy values depend strongly on both models. Further work should examine whether tuned Laplace parameters, a cheaper selection step or feasibility-preserving initialization improve LX-SSA, re-optimize layouts under higher-fidelity wake and power models, and study realistic farm sizes with correspondingly larger budgets and explicit constraint handling.

"""


def appendix():
    B = pd.read_csv("fresh_best_layouts_maxN.csv")
    B = B.sort_values("Objective").groupby(["Dataset", "Radius", "Turbines"]).tail(1).sort_values(["Dataset", "Radius"])
    lab = {"LXSSA": "LX-SSA", "SSA": "SSA", "PSO": "PSO", "DE": "DE", "VNS": "VNS", "SLSQP": "MS-SLSQP"}
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
Table~\ref{tab:coords} lists the turbine coordinates (m, farm center at the origin) of the best layout found by any method for the largest turbine count of each farm (Fig.~\ref{fig:layouts}). The final coordinates of all 12,240 runs, their convergence curves and the Horns Rev~1 layouts are provided as machine-readable files with the revision.

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
    setup = open("sec_setup.tex").read()
    i = t.index(r"\section{\textbf{Experimental Analysis}}")
    j = t.index(r"\section{Scope, Limitations, and Validity of the Evidence}")
    t = t[:i] + setup + results + sensitivity + complexity + t[j:]
    i = t.index(r"\section{Scope, Limitations, and Validity of the Evidence}")
    j = t.index(r"\section{Conclusion}")
    t = t[:i] + limitations + t[j:]
    i = t.index(r"\section{Conclusion}")
    j = t.index(r"\appendices")
    t = t[:i] + conclusion + t[j:]
    i = t.index(r"\appendices")
    j = t.index("\\vspace{1em}\n\\noindent\\textbf{Compliance with Ethical Standards}")
    t = t[:i] + appendix() + t[j:]
    open(MS, "w").write(t)
    print("manuscript rebuilt")


if __name__ == "__main__":
    main()

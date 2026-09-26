"""Final manuscript: hybrid LX-SSA-VNS = LX-SSA + original basic VNS, compared with LX-SSA, SSA, PSO, DE,
the original basic VNS and MS-SLSQP, plus the component ablation (SSA, LX-SSA, VNS, SSA-VNS, LX-SSA-VNS).
Rewrites everything from the hybrid-algorithm section to the appendix. Numbers were checked against
final_summary.json, final_case_tests.csv, final_ablation_tests.csv, final_robust_summary.json and the
run files (fresh_grid.csv, fresh_vgrid.csv, fresh_bgrid.csv, fresh_bsplit.csv, Horns Rev files)."""
import pandas as pd

MS = "../LXSSA_WFLOP_reviewer_revised.tex"


def blocks(path):
    txt = open(path).read()
    return {b.split("\\label{")[1].split("}")[0]: ("\\begin{table" + b).rstrip()
            for b in txt.split("\\begin{table") if "\\label{" in b}


def fig(name, label, caption, wide=True, width="\\textwidth"):
    env = "figure*" if wide else "figure"
    return (f"\\begin{{{env}}}[!t]\n\\centering\n\\includegraphics[width={width}]{{figures_final/{name}.pdf}}\n"
            f"\\caption{{{caption}}}\n\\label{{{label}}}\n\\end{{{env}}}\n")


T = blocks("final_tables.tex")
T["tab:wtl"] = T["tab:wtl"].replace("{2.5pt}", "{2pt}")
RB = blocks("final_robust_table.tex")
CAP = blocks("robustness_tables.tex")
SP = blocks("authors_tables.tex")


def setup_text():
    s = open("sec_setup.tex").read()
    old_vns = "\\item \\emph{VNS}" + s.split("\\item \\emph{VNS}")[1].split("\n")[0]
    reps = [
        ("\\subsection{Methods Compared}\nSix methods are compared, all with the same objective~\\eqref{eq:penalty}, the same bounds $[-r,r]^{2N}$ and the same budget:\n\\begin{itemize}\n",
         "\\subsection{Methods Compared}\nThe proposed hybrid and six reference methods are compared, all with the same objective~\\eqref{eq:penalty}, the same bounds $[-r,r]^{2N}$ and the same budget:\n\\begin{itemize}\n\\item \\emph{LX-SSA-VNS} (proposed): Algorithm~\\ref{alg:hybrid} with $\\rho=0.5$, i.e., 50 LX-SSA iterations (3,030 calls) followed by 3,000 calls of basic VNS.\n"),
        (old_vns,
         "\\item \\emph{VNS}: the original basic VNS of Mladenovi\\'c and Hansen~\\cite{Mladenovic1997,Hansen2001} in the continuous form of Mladenovi\\'c et al.~\\cite{Mladenovic2008}, i.e., exactly Phase~2 of Algorithm~\\ref{alg:hybrid} ($\\ell_\\infty$ shaking of the whole layout with $k_{\\max}=5$ and $\\rho_k=0.1kr$; complete best-improvement compass local search), started from the best member of the initial population instead of the LX-SSA food source."),
        ("SSA, LX-SSA, PSO and DE use the original implementation; VNS and MS-SLSQP were implemented for this study.",
         "SSA, LX-SSA, PSO and DE use the original implementation. VNS, MS-SLSQP and the hybrid were implemented for this study, and the hybrid calls the original LX-SSA code in its first phase. For the ablation (Section~\\ref{sec:ablation}) a fifth variant, \\emph{SSA-VNS}, replaces LX-SSA by SSA in Phase~1 (100 SSA iterations, 3,030 calls) and is otherwise identical to LX-SSA-VNS."),
        ("so runs with the same seed start from the \\emph{same} initial population for all six methods.",
         "so runs with the same seed start from the \\emph{same} initial population for all methods."),
        ("For LX-SSA this corresponds to 100 iterations, because each follower evaluates two candidates before the population is re-evaluated ($30+100\\times60$ calls); SSA, PSO and DE reach the same count in 200 iterations ($30+200\\times30$).",
         "For LX-SSA this corresponds to 100 iterations, because each follower evaluates two candidates before the population is re-evaluated ($30+100\\times60$ calls); SSA, PSO and DE reach the same count in 200 iterations ($30+200\\times30$), and LX-SSA-VNS in $3{,}030+3{,}000$ calls."),
        ("Within each test case, the 30 seed-paired runs are compared by a run-level Friedman test and by two-sided Wilcoxon signed-rank tests of LX-SSA against each of the other five methods. The five $p$ values of each case are adjusted by Holm's procedure~\\cite{Holm1979},",
         "Within each test case, the 30 seed-paired runs are compared by a run-level Friedman test and by two-sided Wilcoxon signed-rank tests of LX-SSA-VNS against each of the other six methods. The six $p$ values of each case are adjusted by Holm's procedure~\\cite{Holm1979},"),
    ]
    for a, b in reps:
        assert s.count(a) == 1, a[:60]
        s = s.replace(a, b)
    return s


results = r"""\section{Numerical Results}
\label{sec:results}
\label{sec:equalbudget}
This section reports 14,280 runs: 68 test cases, seven methods and 30 seed-paired runs, all at 6,030 objective calls. It is followed by the component ablation of the hybrid (a further 2,040 runs of SSA-VNS and 720 runs with other budget splits) and the Horns Rev~1 site case. Detailed per-case results are given in Tables~\ref{tab:res-1-500}--\ref{tab:res-2-1000}.

\subsection{Overall Comparison}
\label{sec:overall}
""" + T["tab:friedman68"] + "\n" + fig("avg_ranks", "fig:avgranks",
    "Average rank of the seven methods over the 68 benchmark cases (ranks of the mean feasible objective within each case; 1 = best).",
    wide=False, width="0.95\\columnwidth") + T["tab:wtl"] + r"""
The case-level Friedman test over the 68 cases rejects equal performance of the seven methods ($\chi^2_F=168.4$, 6 d.f., $p=9.9\times10^{-34}$; Iman--Davenport $F_F=47.1$). Table~\ref{tab:friedman68} and Fig.~\ref{fig:avgranks} give the average ranks (in parentheses: cases in which the method alone has the highest mean):
\begin{itemize}
\item \textbf{LX-SSA-VNS 1.89} (31 cases);
\item VNS 2.85 (13 cases);
\item SSA 3.85, PSO 4.14, MS-SLSQP 4.44 (12 cases), LX-SSA 4.76 and DE 6.07.
\end{itemize}
In the Holm-adjusted post hoc tests, LX-SSA-VNS ranks significantly better than every other method, including VNS ($p_{\rm Holm}=0.010$); against the remaining five methods $p_{\rm Holm}\le2.3\times10^{-7}$.

The run-level pairwise tests (Table~\ref{tab:wtl}) show where the hybrid gains:
\begin{itemize}
\item \emph{LX-SSA-VNS vs.\ LX-SSA.} The hybrid is significantly better in 39 cases and never worse ($r_{rb}$ from $0.53$ to $1.00$, median $0.74$). Of the 29 cases without a significant difference, 22 have $N\le7$, where LX-SSA already comes close to the best layouts.
\item \emph{LX-SSA-VNS vs.\ SSA, PSO and DE.} The hybrid is significantly better in 32, 36 and 52 cases, respectively, and never worse.
\item \emph{LX-SSA-VNS vs.\ MS-SLSQP.} The hybrid is better in 36 cases and worse in four dense cases (Data Set~I, 500~m, $N=8$--10; Data Set~II, 750~m, $N=12$), where MS-SLSQP's explicit constraint handling pays off.
\item \emph{LX-SSA-VNS vs.\ VNS.} The two are statistically indistinguishable in 60 cases. The hybrid is significantly better in six cases with $N=5$ to 9 ($r_{rb}$ from $0.51$ to $0.63$) and worse in the two three-turbine 500-m cases.
\end{itemize}
Hybridization therefore lifts LX-SSA from the sixth to the first average rank. The margin over VNS, which forms the second phase of the hybrid, is small at the level of single cases but consistent across the benchmark.

\subsection{Solution Quality by Farm Size and Turbine Count}
""" + fig("wakeloss_vs_n", "fig:wakeloss",
    "Mean wake loss (percentage of the wake-free objective) of the feasible runs versus the prescribed number of turbines; a point is omitted when a method has no feasible run.") + r"""
Figure~\ref{fig:wakeloss} and Tables~\ref{tab:res-1-500}--\ref{tab:res-2-1000} show that the wake loss grows with $N$, falls with the farm radius, and is higher for the broader wind rose of Data Set~II. For $N\le3$ all methods stay within 0.3\% of the wake-free objective. With growing $N$ the curve of the hybrid separates from that of LX-SSA and joins the VNS and MS-SLSQP curves, as in these examples:
\begin{itemize}
\item Data Set~I, 1000~m/15: LX-SSA-VNS 1.78\%, VNS 1.99\%, MS-SLSQP 2.20\%, SSA 2.80\%, LX-SSA 3.32\%;
\item Data Set~II, 750~m/12: LX-SSA-VNS 5.95\%, VNS 6.17\%, MS-SLSQP 5.60\%, SSA 6.79\%, LX-SSA 7.34\%.
\end{itemize}
In the largest case of every farm the hybrid reduces the mean wake loss of LX-SSA by 1.4 to 1.9 percentage points. It attains the lowest mean wake loss in both 1000-m fifteen-turbine cases, while MS-SLSQP has the lowest value in the other four largest cases.

""" + "\n".join(T[f"tab:res-{d}-{r}"] for d in (1, 2) for r in (500, 750, 1000)) + r"""

\subsection{Feasibility}
""" + fig("feasibility_vs_n", "fig:feasibility",
    "Percentage of the 30 runs that end with a feasible layout (all turbines inside the farm and at least 308~m apart).") + r"""
Figure~\ref{fig:feasibility} shows that LX-SSA-VNS returns a feasible layout in 99.9\% of all runs: only one run in each of the two 500-m ten-turbine cases ends infeasible. The other methods reach 96.9\% (LX-SSA), 97.6\% (SSA), 86.0\% (PSO), 71.4\% (DE), 100\% (VNS) and 99.9\% (MS-SLSQP). In the densest cases (500~m, $N=10$) LX-SSA is feasible in only 13 of 30 runs; the complete local search of the VNS phase repairs almost all of the layouts it hands over. The feasibility counts of the penalty-based methods are identical for the two wind data sets, because while every candidate is infeasible the penalty term dominates $F_p$.

\subsection{Convergence Behaviour}
""" + fig("convergence_mid", "fig:conv-mid",
    "Convergence for moderate turbine counts: median (line) and interquartile range (band) of the best feasible wake loss found so far over the 30 runs, versus the number of objective calls (logarithmic scale). A curve starts once at least half of the runs have found a feasible layout; the LX-SSA-VNS curve switches from the LX-SSA to the VNS phase at 3,030 calls.") + \
    fig("convergence_max", "fig:conv-max",
    "Convergence for the largest turbine count of each farm (as Fig.~\\ref{fig:conv-mid}). Methods for which fewer than half of the runs find a feasible layout have no curve.") + r"""
Figures~\ref{fig:conv-mid} and~\ref{fig:conv-max} show the effect of the hybridization directly. Until 3,030 calls the hybrid runs a shortened LX-SSA and is therefore behind the full LX-SSA run. After the switch its curve drops steeply, whereas stand-alone LX-SSA and SSA stagnate after about 3,000 calls. Over all feasible runs, the VNS phase removes on average 39.0\% (median 31.9\%) of the wake loss present at the end of the LX-SSA phase, and it turns 118 runs that were still infeasible at the switch into feasible ones. Stand-alone VNS reaches low wake losses earlier because it spends its whole budget on intensification, but the hybrid catches up after the switch and ends with a lower mean wake loss in five of the six largest cases. MS-SLSQP is the fastest method to reach feasibility.

\subsection{Distribution of Results}
""" + fig("boxplots_mid", "fig:box",
    "Distribution of the final wake loss of the feasible runs for moderate turbine counts (30 runs per method).") + r"""
Figure~\ref{fig:box} shows the run-to-run distributions for the moderate cases. LX-SSA-VNS has the lowest median wake loss in all three Data Set~II cases and the second-lowest in the 500-m Data Set~I case. Its interquartile range is smaller than that of LX-SSA in five of the six cases. DE has the highest median in five cases.

\subsection{Component Ablation}
\label{sec:ablation}
""" + T["tab:ablation"] + "\n" + fig("ablation_convergence", "fig:ablation",
    "Component ablation, largest turbine count of each farm: median best feasible wake loss over 30 runs versus objective calls for SSA, LX-SSA, VNS and LX-SSA-VNS (dotted vertical line: switch from LX-SSA to VNS in the hybrid at 3,030 calls).") + \
    T["tab:split"] + r"""
The ablation separates the three ingredients of the hybrid: the swarm phase, the Laplace step inside that swarm phase, and the VNS phase. It compares SSA, LX-SSA, VNS, SSA-VNS and LX-SSA-VNS on all 68 cases with the same seeds and budget (Table~\ref{tab:ablation}, Fig.~\ref{fig:ablation}). The Friedman test over the five variants is highly significant ($\chi^2_F=156.5$, 4 d.f., $p=8.1\times10^{-33}$).
\begin{itemize}
\item \emph{VNS phase.} Adding the VNS phase is the decisive step. LX-SSA-VNS beats LX-SSA in 39 cases and SSA-VNS beats SSA in 40 cases, with no losses; the mean wake loss falls by 0.43 and 0.33 percentage points.
\item \emph{Swarm phase.} Starting VNS from the food source of a swarm instead of the best initial point helps moderately. Both two-phase variants have a better average rank than VNS (2.25 and 1.62 vs.\ 2.90); LX-SSA-VNS beats VNS in five cases and loses two, and SSA-VNS beats it in ten cases and loses one.
\item \emph{Laplace step.} The Laplace step of LX-SSA does not help. Without VNS, LX-SSA is worse than SSA in seven cases and never better. Inside the hybrid, LX-SSA-VNS and SSA-VNS are indistinguishable in 66 cases, SSA-VNS is better in two, and SSA-VNS has the better average rank (1.62 vs.\ 2.25).
\item \emph{Budget split.} Varying the share of LX-SSA between 25\%, 50\% and 75\% of the budget changes little (Table~\ref{tab:split}). The 50\% split has the highest mean in eight of 12 cases, but only one difference is significant (Data Set~I, 1000~m/15: 50\% better than 25\%, $p_{\rm Holm}=0.024$).
\end{itemize}
Taken together, the improvement of LX-SSA-VNS over its components comes from combining a swarm exploration phase with a complete VNS intensification phase. The specific Laplace mechanism of LX-SSA is not the source of the gain: a standard SSA first phase performs at least as well. The counts against VNS differ slightly from Table~\ref{tab:wtl} because the Holm adjustment is applied to a different family of comparisons.

\subsection{Optimized Layouts}
""" + fig("layouts_max", "fig:layouts",
    "Best layouts for the largest turbine count of each farm: best LX-SSA-VNS run (filled blue circles) and best run of all methods (open black squares; the method is named above each panel). Where the two coincide, the circles sit inside the squares. Circle: farm boundary. Coordinates are listed in Appendix~\\ref{app:coordinates}.") + r"""
Figure~\ref{fig:layouts} compares, for the largest turbine count of each farm, the best LX-SSA-VNS layout with the best layout found by any method. The best overall layout comes from LX-SSA-VNS in four cases, from MS-SLSQP in one (Data Set~I, 500~m) and from VNS in one (Data Set~II, 750~m). Most turbines are placed close to the boundary, which maximizes their mutual distances, while the remaining turbines occupy the interior at positions staggered with respect to the dominant wind directions.

\subsection{Computational Cost}
\label{sec:runtime}
At an equal number of objective calls all methods need almost the same time: 0.37~ms per call for VNS, 0.39~ms for LX-SSA-VNS, 0.40~ms for LX-SSA and SSA, 0.41~ms for PSO, 0.47~ms for DE and 0.49~ms for MS-SLSQP. The objective evaluation dominates the cost, so neither the Laplace step nor the VNS bookkeeping adds measurable overhead. A 6,030-call run of LX-SSA-VNS takes 1.7--3.5~s on average, depending on $N$.

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
The farm model keeps the Jensen top-hat wake (with $k=0.04$, a value commonly used offshore) and root-sum-square superposition. It uses the tabulated power curve with a 25~m/s cut-out, a speed-dependent thrust coefficient at the free-stream speed, 5$^\circ$ direction bins and 1~m/s speed bins. For the installed layout it gives 664.6~GWh/yr and a wake loss of 10.95\%, against 662.5~GWh/yr and 10.96\% from PyWake's Jensen model with the same settings, a difference below 0.3\%. The 16-turbine block that PyWake uses as an example (30 seeds) and the complete 80-turbine farm (10 seeds) were optimized with all seven methods at 6,030 calls and the $4D$ (320~m) minimum spacing.

Table~\ref{tab:hr-site} and Figs.~\ref{fig:hr-conv} and~\ref{fig:hr-layouts} lead to three observations.
\begin{itemize}
\item \emph{The benchmark ordering carries over to the site.} In the 16-turbine block LX-SSA-VNS attains the highest mean AEP (137.77~GWh/yr, feasible in all 30 runs). It is significantly better than LX-SSA ($p_{\rm Holm}=5.7\times10^{-6}$), SSA ($p_{\rm Holm}=3.1\times10^{-4}$) and MS-SLSQP ($p_{\rm Holm}=2.4\times10^{-4}$); its difference from VNS (137.67~GWh/yr) is not significant ($p_{\rm Holm}=0.685$). PSO and DE never find a feasible layout inside the parallelogram.
\item \emph{No method improves on the installed layout within the budget.} The installed $7D$ grid reaches 139.51~GWh/yr.
\item \emph{The full farm exceeds what the penalty-based methods can handle at this budget.} With 160 variables and 3,160 spacing constraints, none of the runs of LX-SSA-VNS, LX-SSA, SSA, PSO, DE or VNS returns a feasible layout. MS-SLSQP, which handles the constraints explicitly, is feasible in all ten runs, but its mean AEP (653.9~GWh/yr) remains 1.6\% below the installed layout.
\end{itemize}
Realistic farm sizes therefore require larger budgets, feasible initialization or explicit constraint handling, and a regular installed layout is a strong reference.

"""

sensitivity = r"""\section{Sensitivity to the Modelling Assumptions}
\label{sec:robustness}
\subsection{Power Curve and Wake Model}
\label{sec:powercurve}
\label{sec:gaussian}
""" + RB["tab:robust-final"] + r"""
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

Table~\ref{tab:robust-final} shows two kinds of sensitivity.
\begin{itemize}
\item \emph{Absolute energy is strongly model-dependent.} The linear benchmark curve overstates expected power relative to the cubic curve by 17.7\% for Data Set~I and 36.9\% for Data Set~II (21.5\% and 37.1\% with the cut-out). The absolute energy values of this paper are therefore benchmark values, not turbine-specific AEP predictions.
\item \emph{The algorithm comparison is robust to the power curve but less so to the wake model.} Under the cubic curve the within-case orderings are almost unchanged (mean Kendall $\bar\tau=0.97$), and LX-SSA-VNS (1.95) and VNS (2.85) remain the two best methods. Under the Gaussian model the orderings change more ($\bar\tau=0.54$; same best method in 43\% of the cases), because layouts optimized for the sharp Jensen wake cone are penalized differently by the smooth Gaussian deficit. LX-SSA-VNS keeps the best average rank (2.82), followed by PSO (3.16), MS-SLSQP (3.43) and VNS (3.76).
\end{itemize}

\subsection{Minimum Spacing}
\label{sec:spacing}
""" + CAP["tab:capacity"] + "\n" + SP["tab:spacing-authors"] + r"""
The $4D$ spacing is a legacy benchmark value. Its influence was examined in three ways.

\emph{Geometric capacity.} Table~\ref{tab:capacity} gives the largest turbine count for which a multi-start packing search constructed a feasible layout; these are constructive lower bounds, and the 500-m values coincide with the known optimal circle-packing results. At $5D$ the 500-m farm holds at most 8 turbines and at $6D$ only 7, so the tested 500-m cases with $N\ge9$ (at $5D$) or $N\ge8$ (at $6D$) would be infeasible.

\emph{Final layouts.} The optimized layouts are frequently active at the $4D$ constraint. Of the 13,293 feasible final layouts of the seven methods, 48.6\% also satisfy $5D$ and 35.8\% satisfy $6D$. For $N=6$--10 these fractions fall to 30.1\% and 11.6\%, and for $N=11$--15 to 3.0\% and 0.1\%. Optimized layouts therefore cannot be transferred to a site with larger spacing requirements without re-optimization.

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
and a run with a budget of $B$ objective calls costs $\mathcal{O}(B\,N_t^2N_\theta N_s)$ for every method considered here, including both phases of LX-SSA-VNS. Box repair is $\mathcal{O}(N_t)$ per candidate and the spacing check $\mathcal{O}(N_t^2)$; the Laplace draw of LX-SSA, the shaking step of VNS and each compass move are $\mathcal{O}(N_t)$. All of these are dominated by the wake evaluation, which is consistent with the almost identical measured cost per call (Section~\ref{sec:runtime}).

"""

limitations = r"""\section{Scope, Limitations, and Validity of the Evidence}
The results should be interpreted within the following boundaries:
\begin{itemize}
\item \textbf{Source of the improvement:} the ablation shows that the gain of LX-SSA-VNS comes from the combination of a swarm phase with a complete VNS phase. The Laplace step of LX-SSA does not contribute: an SSA-VNS variant is statistically equivalent in 66 of 68 cases and better in two, and its average rank is slightly better.
\item \textbf{Margin over VNS:} the hybrid is significantly better than stand-alone VNS at the case level, but at the level of single cases the two differ significantly in only eight of 68 cases.
\item \textbf{Budget and farm size:} all methods were compared at 6,030 objective calls; larger budgets could change the ordering. For an 80-turbine farm the penalty-based methods, including the hybrid, did not find feasible layouts at this budget.
\item \textbf{Algorithm settings:} the LX-SSA parameters ($\phi=0$, $\chi=1$), the VNS settings and the split $\rho=0.5$ are fixed defaults and were not tuned. Pseudo-gradient and surrogate-assisted methods are discussed but not benchmarked.
\item \textbf{Wake model and power curve:} Jensen's model and the piecewise-linear power curve are benchmark choices. Absolute energy values depend strongly on them, and the within-case orderings change under a Gaussian wake model (Section~\ref{sec:gaussian}); layouts were not re-optimized with the alternative models.
\item \textbf{Spacing and boundary handling:} $4D$ is a legacy constraint that governs the feasible turbine density, and the boundary-handling rule also affects absolute results.
\item \textbf{Scope:} $N$ is prescribed in every run. Terrain, heterogeneous turbines, electrical collector systems, grid, acoustic and environmental constraints, and economic objectives are not modeled.
\end{itemize}

"""

conclusion = r"""\section{Conclusion}
This paper proposes LX-SSA-VNS, a two-phase hybrid in which the previously published Laplacian Salp Swarm Algorithm (LX-SSA)~\cite{Solanki2023} explores the continuous layout space and the original basic variable neighbourhood search then intensifies the best layout found. The hybrid is applied to the continuous WFLOP with a Jensen--Weibull benchmark model and explicit constraint handling. It is evaluated against LX-SSA, SSA, PSO, DE, VNS and a gradient-based multistart SLSQP solver on 68 benchmark cases (14,280 seed-paired runs at equal budgets), in a component ablation, and on the Horns Rev~1 site.

The results support five conclusions.
\begin{itemize}
\item \emph{The hybrid is the best of the tested methods.} LX-SSA-VNS has the best average rank (1.89), is significantly better than every other method at the case level (VNS: $p_{\rm Holm}=0.010$), and returns a feasible layout in 99.9\% of the runs.
\item \emph{Hybridization substantially improves LX-SSA.} The hybrid is significantly better than LX-SSA in 39 of 68 cases and never worse, and it lifts LX-SSA from the sixth to the first average rank.
\item \emph{The gain comes from the two-phase design, not from the Laplace step.} The VNS phase removes on average 39\% of the wake loss left by the swarm phase, and starting VNS from a swarm solution is better than starting it from the best initial point. Replacing LX-SSA by standard SSA in the first phase, however, performs at least as well.
\item \emph{A similar ordering holds at a real site.} In a 16-turbine block of Horns Rev~1 the hybrid has the highest mean AEP; it is significantly better than LX-SSA, SSA and MS-SLSQP and statistically indistinguishable from VNS. No method improves on the installed layout within the budget, and for the full 80-turbine farm only the constraint-aware MS-SLSQP finds feasible layouts.
\item \emph{The comparison depends on the wake model more than on the power curve.} The algorithm ordering is insensitive to the power curve but changes considerably under a Gaussian wake model, although the hybrid keeps the best average rank; absolute energy values depend strongly on both.
\end{itemize}

Future work should use a stagnation-triggered switch between the two phases, examine which properties of the swarm phase make it a good starting point for VNS, and extend the approach with feasibility-preserving initialization, larger budgets for realistic farm sizes and higher-fidelity wake and power models.

"""


def appendix():
    B = pd.read_csv("final_best_layouts_maxN.csv")
    B = B.sort_values("Objective").groupby(["Dataset", "Radius", "Turbines"]).tail(1).sort_values(["Dataset", "Radius"])
    lab = {"LXBV": "LX-SSA-VNS", "LXSSA": "LX-SSA", "SSA": "SSA", "PSO": "PSO", "DE": "DE", "BVNS": "VNS",
           "SLSQP": "MS-SLSQP"}
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
Table~\ref{tab:coords} lists the turbine coordinates (m, farm center at the origin) of the best layout found by any method for the largest turbine count of each farm (Fig.~\ref{fig:layouts}). The final coordinates and convergence curves of all runs and the Horns Rev~1 layouts are provided as machine-readable files with the revision.

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
    i = t.index(r"\section{Hybrid LX-SSA-VNS Algorithm}")
    j = t.index(r"\appendices")
    t = t[:i] + open("sec_hybrid.tex").read() + setup_text() + results + sensitivity + complexity + limitations + conclusion + t[j:]
    i = t.index(r"\appendices")
    j = t.index("\\vspace{1em}\n\\noindent\\textbf{Compliance with Ethical Standards}")
    t = t[:i] + appendix() + t[j:]
    open(MS, "w").write(t)
    print("manuscript rebuilt (final)")


if __name__ == "__main__":
    main()

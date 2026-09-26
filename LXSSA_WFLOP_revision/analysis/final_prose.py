"""Final prose pass: rewrites the running text of the manuscript in the authors' own style (first person
plural, paragraphs rather than bullet lists, plain wording) and switches Figs. 1-3 to the redrawn
vector figures. Tables, algorithms, equations and all numbers are kept exactly; floats are moved
verbatim by label. Run once on the output of build_final_manuscript.py + final_frontmatter.py."""
import re

MS = "../LXSSA_WFLOP_reviewer_revised.tex"
t = open(MS).read()


def span(start, end):
    i = t.index(start); j = t.index(end, i + len(start))
    return i, j


def replace(start, end, new):
    global t
    i, j = span(start, end)
    t = t[:i] + new + t[j:]


def flt(label):
    """Float environment (table/figure/algorithm) that contains \\label{label}, verbatim."""
    k = t.index("\\label{" + label + "}")
    best = None
    for env in ("table*", "table", "figure*", "figure", "algorithm"):
        b = t.rfind("\\begin{" + env + "}", 0, k)
        if b >= 0 and (best is None or b > best[0]):
            best = (b, env)
    b, env = best
    e = t.index("\\end{" + env + "}", k) + len("\\end{" + env + "}")
    return t[b:e] + "\n"


def maths(start, end):
    """Display-math blocks of a span, in order (optionally wrapped in {\\small ...})."""
    i, j = span(start, end)
    pat = re.compile(r"(\{\\small\n)?\\begin\{(align\*?|multline|equation)\}.*?\\end\{\2\}(\n\})?", re.S)
    return [m.group(0) for m in pat.finditer(t[i:j])]


# ------------------------------------------------------------------ figures 1-3
FIG1 = r"""\begin{figure}[!t]
\centering
\includegraphics[width=\columnwidth]{figures_final/fig_wind_farm.pdf}
\caption{Circular wind farm of radius $r$. The wind direction $\theta$ is measured counter-clockwise from east. The points show the best LX-SSA-VNS layout for 12 turbines in the 750-m farm (Wind Data Set~I); the shaded discs have radius $2D$, so two discs touch when their turbines are exactly $4D$ apart.}
\label{fig:wind farm}
\end{figure}
"""
FIG2 = r"""\begin{figure}[!t]
\centering
\includegraphics[width=\columnwidth]{figures_final/fig_wake_model.pdf}
\caption{Jensen wake behind a rotor of radius $R$ (schematic, not to scale). The wake widens linearly with the spreading constant $K$, so that its width at the downstream distance $d$ is $2(R+Kd)$, and the wind speed inside it drops from $s_{\text{up}}$ to $s_{\text{down}}$.}
\label{fig:wake-model}
\end{figure}
"""
FIG3 = r"""\begin{figure}[!t]
\centering
\includegraphics[width=\columnwidth]{figures_final/fig_half_cone.pdf}
\caption{Wake half cone of turbine $j$ (schematic, $K$ exaggerated). The virtual vertex $A$ lies a distance $R/K$ upstream of the rotor, the half-cone angle is $\alpha=\arctan K$, and a downstream turbine $i$ lies in the wake of $j$ when $\beta_{ij}<\alpha$.}
\label{fig:half-cone}
\end{figure}
"""

# ------------------------------------------------------------------ abstract
replace(r"\begin{abstract}", r"\end{abstract}", r"""\begin{abstract}
The wind farm layout optimization problem (WFLOP) asks where a given number of turbines should be placed so that as little energy as possible is lost in the wakes, while every turbine stays inside the site and keeps a minimum distance from the others. In this work, we propose LX-SSA-VNS, a two-phase hybrid algorithm for the continuous WFLOP. In the first phase, the Laplacian Salp Swarm Algorithm (LX-SSA), introduced by Solanki and Deep~\cite{Solanki2023}, explores the layout space for half of the evaluation budget. In the second phase, the basic variable neighborhood search (VNS) refines the best layout found by the swarm by alternating shaking of the whole layout with a complete local search. Wake losses are computed with the Jensen model and direction-dependent Weibull wind data, and the boundary and spacing constraints are handled through a penalized objective. We compare LX-SSA-VNS with LX-SSA, SSA, particle swarm optimization (PSO), differential evolution (DE), VNS and a multistart SLSQP solver on 68 benchmark cases (two wind data sets, three circular farms, 2--15 turbines), using 30 seed-paired runs of 6,030 objective evaluations per method and case, and analyze the results with Friedman, Holm-corrected Wilcoxon and effect-size statistics. LX-SSA-VNS obtains the best average rank (1.89), is significantly better than each of the other six methods and finds a feasible layout in 99.9\% of the runs. An ablation with SSA, LX-SSA, VNS and SSA-VNS shows that the improvement comes from pairing a swarm phase with VNS, which removes on average 39\% of the wake loss left by the swarm, whereas the Laplace step itself adds nothing: an SSA first phase does at least as well. On the Horns Rev~1 offshore site, the hybrid gives the highest mean energy yield for a 16-turbine block, although no method improves on the installed layout within the budget and only the constraint-aware SLSQP solver finds feasible layouts for the full 80-turbine farm. Finally, we show that the ranking of the algorithms is stable under a different power curve but changes under a Gaussian wake model.
""")

# ------------------------------------------------------------------ introduction (incl. former Section II)
replace(r"\section{Introduction}", r"\section{\textbf{Problem Modeling}}", r"""\section{Introduction}

Wind energy has become one of the main sources of low-carbon electricity. Wind turbines convert the kinetic energy of moving air into electrical energy without direct greenhouse gas emissions, and the installed capacity has grown quickly over the last two decades. Since most of the cost of a wind project is fixed at the planning stage, the decisions taken before construction have a lasting effect on the energy a farm will deliver, and one of the most important of these decisions is where the turbines are placed.

The Wind Farm Layout Optimization Problem (WFLOP) deals with exactly this question. An upstream turbine extracts momentum from the flow and leaves behind a wake of slower and more turbulent air, so that the turbines downstream produce less power. The wake widens and recovers with distance, which suggests spreading the turbines as far apart as possible. The site, however, is limited, and the turbines must keep a minimum distance from each other. Placing more turbines in the same area raises the installed capacity but also strengthens the wake interaction, so that the energy produced per turbine falls. The objective is nonlinear, the wakes couple the positions of all turbines, and the feasible region is bounded by many nonconvex spacing constraints. For these reasons, WFLOP has become a popular test bed for metaheuristic search.

Mosetti et al.~\cite{Mosetti1994} and Grady et al.~\cite{Grady2005} were among the first to address WFLOP, using a binary-coded genetic algorithm (GA) on a grid to maximize the energy output at minimum installation cost. Ozturk and Norman~\cite{Ozturk2004} extended this work with a different objective and a heuristic that places turbines in continuous space, taking wind speed, wind direction and wake interactions into account. Lackner and Elkinton~\cite{Lackner2007} proposed an analytical model for offshore farms in which the annual power output depends directly on the turbine positions and wakes. Castro et al.~\cite{Mora2007} designed layouts with a GA under spacing and wind-speed constraints, and Huang~\cite{Huang2007} used a distributed GA to maximize the annual profit of large farms from their energy output, capital cost and operating cost. Elkinton et al.~\cite{Elkinton2008} compared GA, particle swarm optimization (PSO), differential evolution (DE) and simulated annealing (SA) for offshore layouts. SA was also used by Bilbao and Alba~\cite{Bilbao2009} to maximize the annual profit and by Rivas et al.~\cite{Rivas2009}, whose version combines three local-search moves for large offshore farms.

Emami and Noghreh~\cite{Emami2010} added the construction cost to the objective and solved the resulting problem with a GA, while Sisbot et al.~\cite{Sisbot2010} used a multi-objective GA to trade energy production against installation cost. Kusiak and Song~\cite{Kusiak2010} introduced the wind-distribution-based placement model that we also use in this paper. They treated the wake loss as a function of the turbine locations and the wind direction and solved a bi-criteria version of the problem with a multi-objective evolution strategy. Eroglu et al.~\cite{Eroglu2012,Eroglu2013} applied ant colony optimization (ACO) and particle filtering to the same model, Samorani~\cite{Samorani2013} studied the trade-off between power output and wake effect, and Bansal and Farswan~\cite{Bansal2017} used biogeography-based optimization (BBO) and proposed an upper limit on the number of turbines for a given farm. ACO has also been used to design the collector cable system of offshore farms~\cite{Srikakulapu2018}.

Later studies refined both the algorithms and the models. Ju and Liu~\cite{Ju2019} proposed adaptive and self-informed variants of the GA, and Jin et al.~\cite{Jin2019} combined a wake-aware power-loss cost model with an adaptive PSO. Hou et al.~\cite{Hou2019} reviewed the layout and electrical-system design of wind farms, including the routing of cables. Dhiman et al.~\cite{Dhiman2020} used lidar simulations to locate the wakes and steer the yaw of downstream turbines, Nagpal et al.~\cite{Nagpal2021} refined the layouts of a distributed GA with a continuous local method, and PSO has been applied to reduce wake losses and increase the power output~\cite{Asaah2021}. Bai et al.~\cite{Bai2022} integrated Monte Carlo tree search into an adaptive GA and applied it to a planned wind farm in New Jersey, and Cazzaro and Pisinger~\cite{Cazzaro2022} developed a variable neighborhood search (VNS) heuristic that handles spacing constraints and construction costs. Beyond classical metaheuristics, pseudo-gradient methods have been developed for large design spaces~\cite{Quaeghebeur2021}, three-dimensional Gaussian wake models have been coupled with layout optimization~\cite{Tao2020}, and machine-learning~\cite{Yang2023}, CFD-based Kriging~\cite{Wang2024} and data-driven or physics-informed~\cite{Yang2024,Li2025} surrogates have been introduced to reduce the cost of higher-fidelity evaluations. No single optimizer or wake model dominates this literature; each study balances model fidelity, computing cost and scalability in its own way.

Swarm methods such as the Salp Swarm Algorithm (SSA)~\cite{Mirjalili2017} are attractive for WFLOP because they are simple and cheap per iteration. Solanki and Deep~\cite{Solanki2023} introduced the Laplacian Salp Swarm Algorithm (LX-SSA), in which each follower salp also tries a Laplace-distributed move relative to the food source and keeps the better of its two positions. The heavy tail of the Laplace distribution produces mostly short moves and occasionally long ones, which helps the swarm to visit separated regions of the layout space. In WFLOP, however, we observed that a swarm spends most of the second half of its budget on small improvements once it has contracted around the best layout, because the final refinement of a layout requires coordinated moves of single turbines along the spacing constraints. Trajectory methods such as VNS~\cite{Mladenovic1997,Hansen2001} behave the other way round: they refine one layout very efficiently by alternating random perturbations of growing size with a local search, but the result depends on where the search starts.

This complementary behavior motivates the present work. We propose LX-SSA-VNS, a two-phase hybrid in which LX-SSA explores the layout space for the first half of the evaluation budget and the basic VNS then refines the best layout found by the swarm. Both phases work on the same continuous turbine-coordinate representation and the same penalized objective. We keep the standard ingredients of the WFLOP benchmark, namely the Jensen wake model, a Weibull description of the wind and the Kusiak--Song power model, so that our results can be compared with the earlier literature, and we examine the effect of these modeling choices separately in Section~\ref{sec:robustness}. Our contributions are summarized as follows:
\begin{enumerate}
\item We propose LX-SSA-VNS, a two-phase hybrid in which the previously published LX-SSA~\cite{Solanki2023} explores the layout space and the basic VNS~\cite{Mladenovic1997,Mladenovic2008} refines the best layout (Section~\ref{sec:hybrid}).
\item We formulate the continuous WFLOP with explicit boundary and minimum-spacing constraints, handled through a penalized objective that both phases of the hybrid share (Section~\ref{sec:model}).
\item We compare the hybrid with LX-SSA, SSA, PSO, DE, VNS and a gradient-based multistart SLSQP solver on 68 benchmark cases at equal numbers of objective evaluations, with 30 seed-paired runs per method and case, and support the comparison with Friedman, Holm-corrected Wilcoxon and effect-size statistics, feasibility rates and convergence curves (Section~\ref{sec:results}).
\item We carry out a component ablation with SSA, LX-SSA, VNS, SSA-VNS and LX-SSA-VNS that separates the effects of the swarm phase, the Laplace step and the VNS phase, and we examine how the budget should be split between the two phases (Section~\ref{sec:ablation}).
\item We apply all methods to the Horns Rev~1 offshore farm with its real turbine, wind climate and site outline (Section~\ref{sec:hornsrev-site}), and we test how the results depend on the power curve, the wake model, the minimum spacing and the boundary handling (Section~\ref{sec:robustness}).
\end{enumerate}

The remainder of this paper is organized as follows. Section~\ref{sec:model} describes the WFLOP model and the constraint handling. Section~\ref{sec:lxssa} summarizes LX-SSA, and Section~\ref{sec:lxssa-wflop} explains how it is applied to WFLOP. Section~\ref{sec:hybrid} introduces the hybrid LX-SSA-VNS. Section~\ref{sec:setup} describes the experimental setup, and Section~\ref{sec:results} reports the numerical results, including the ablation and the Horns Rev~1 case. Section~\ref{sec:robustness} examines the sensitivity of the results to the modeling assumptions, Section~\ref{sec:complexity} discusses the computational complexity, and Section~\ref{sec:limits} states the limitations of the study. Finally, Section~\ref{sec:conclusion} concludes the paper.

""")

# ------------------------------------------------------------------ problem modeling
replace(r"\section{\textbf{Problem Modeling}}", r"\subsection{\textbf{The Wake Effect Model}}", r"""\section{Problem Modeling}
\label{sec:model}

We use the following definitions and assumptions for the design of a wind farm.

\begin{enumerate}
    \item \textbf{Prescribed turbine count:} The nameplate capacity of a farm of identical turbines grows linearly with the number of turbines, but the expected farm power does not, because the wake losses depend on the layout. In every optimization run the number of turbines $N$ is fixed in advance, and only the $2N$ turbine coordinates are optimized. Results for several values of $N$ therefore form a sweep over the turbine count; $N$ and the turbine positions are not optimized jointly.

    \item \textbf{Turbine location:} Each turbine is located by its planar coordinates $(\rho,\sigma)$, measured from the center of the farm.

    \item \textbf{Identical turbines:} All turbines have the same power curve, rated capacity, rotor size and hub height.

    \item \textbf{Wind speed:} For a given site, height and wind direction, the wind speed $s$ follows a Weibull distribution with probability density function
    \[
        p_s (s;\zeta,\psi) = \frac{\zeta}{\psi} \left( \frac{s}{\psi} \right)^{\zeta-1} \exp\left(-\left( \frac{s}{\psi} \right)^\zeta \right),
    \]
    where $\zeta$ is the shape parameter and $\psi$ the scale parameter.

    \item \textbf{Wind direction:} The wind speed and both Weibull parameters depend on the wind direction $\theta$, i.e., $s = s(\theta)$, $\zeta = \zeta(\theta)$ and $\psi = \psi(\theta)$ with $0^\circ \leq \theta \leq 360^\circ$. Fig.~\ref{fig:wind farm} shows the convention used throughout the paper: $\theta=0^\circ$ points east, $90^\circ$ north, $180^\circ$ west and $270^\circ$ south.

    \item \textbf{Minimum spacing:} Any two turbines must be at least $4D=8R$ apart, where $D=2R$ is the rotor diameter:
    \[
        (\rho_i - \rho_j)^2 + (\sigma_i - \sigma_j)^2 \geq 64R^2 .
    \]
    We use this value because it is the one adopted in the benchmark studies we compare with~\cite{Kusiak2010,Bansal2017}. Real projects often require larger and direction-dependent distances, and the spacing limits how many turbines fit into a farm; Section~\ref{sec:spacing} studies the effect of $5D$ and $6D$ spacing.

    \item \textbf{Farm shape:} The farm is circular, as in the benchmark studies~\cite{Kusiak2010,Bansal2017}. A circle has no preferred orientation, so the boundary itself favors no wind direction. We make no claim that a circular site uses land better than other shapes.

    \item \textbf{Farm boundary:} Every turbine $T_i (\rho_i, \sigma_i)$ must lie inside the farm of radius $r$,
    \[
        \rho_i^2 + \sigma_i^2 \leq r^2.
    \]
    We study farms with $r=500$, 750 and 1000~m.

    \item \textbf{Search space:} The search space is the circular farm area shown in Fig.~\ref{fig:wind farm}. Because it is continuous, the turbines are not restricted to the nodes of a grid.

    \item \textbf{Model:} The WFLOP model consists of a wake model, which describes the speed reduction behind each turbine, and a power model. We use the Jensen wake model~\cite{Jensen1983, Katic1986} and the power model of Kusiak and Song~\cite{Kusiak2010}.

    \item \textbf{Objective:} The aim is to maximize the expected power of the farm, i.e., to minimize the wake losses, subject to the spacing and boundary constraints of assumptions~6 and~8.
\end{enumerate}

""" + FIG1 + "\n")

# ------------------------------------------------------------------ wake model
replace(r"\subsection{\textbf{The Wake Effect Model}}", r"\subsection{The Power Model}", r"""\subsection{Wake Effect Model}

A turbine extracts momentum from the wind, so the air behind it moves more slowly and the downstream turbines produce less power. Accounting for these wake losses is therefore at the core of layout design. We use Jensen's wake model~\cite{Jensen1983,Katic1986}, which is cheap to evaluate and has been used in most of the benchmark studies we compare with. Gaussian and three-dimensional wake models describe the velocity field more realistically and have also been used for layout optimization~\cite{Tao2020,Cao2022}; for this reason, Section~\ref{sec:gaussian} re-evaluates all optimized layouts with a Gaussian wake model.

In the Jensen model, the wake behind a rotor expands linearly with the downstream distance, and the wind speed inside it drops from the upstream value $s_{\text{up}}$ to $s_{\text{down}}$ (Fig.~\ref{fig:wake-model}). Let $K$ be the wake spreading constant, $R$ the rotor radius and $d_{ij}$ the distance between turbines $i$ and $j$ along the wind direction. The velocity deficit induced at turbine $j$ by the wake of turbine $i$ is
\begin{align*}
    \delta s_{ij} = 1-(s_{\text{down}}/s_{\text{up}})=\frac{1 - \sqrt{1 - C_T}}{\left(1 + \frac{K d_{ij}}{R}\right)^2},
\end{align*}
where $C_T$ is the thrust coefficient of the turbine.

""" + FIG2 + r"""
When a turbine lies in the wakes of several upstream turbines, the individual deficits are combined by the root-sum-square rule. For a turbine $i$ that is affected by the wakes of the other turbines,
\begin{align}\label{eq:veldeficit2}
    \delta s_i = \sqrt{ \sum_{\substack{j=1 \\ j \ne i}}^{N} (\delta s_{ij})^2}.
\end{align}

Each rotor faces the incoming wind direction $\theta$. The wake of a turbine at $(\rho,\sigma)$ can then be idealized as a half cone whose virtual vertex $A$ lies a distance $R/K$ upstream of the rotor (Fig.~\ref{fig:half-cone}). The half-cone angle $\alpha$ ($0 \leq \alpha \leq \pi/2$), i.e., the angle between the axis and the side of the cone, is $\alpha = \arctan(K)$. A turbine lies in the wake of another if it falls inside this cone, which the following two lemmas make precise.

""" + FIG3 + r"""
\textbf{Lemma 1.} For a given wind direction \( \theta \), let \( \beta_{ij} \) (\( 0 \leq \beta_{ij} \leq \pi \)) be the angle at the vertex \( A \) of the wake of turbine \( j \) between the wake axis and the line from \( A \) to turbine \( i \). Then
\[
\beta_{ij} = \cos^{-1}\!\left(\frac{u_{ij}}{\sqrt{u_{ij}^2+v_{ij}^2}}\right),
\]
with
\begin{align*}
u_{ij}&=(\rho_i - \rho_j)\cos\theta + (\sigma_i - \sigma_j)\sin\theta + \tfrac{R}{K},\\
v_{ij}&=-(\rho_i - \rho_j)\sin\theta + (\sigma_i - \sigma_j)\cos\theta ,
\end{align*}
where \( R/K \) is the distance from the rotor center to the vertex \( A \).

\textbf{Lemma 2.} If turbine \( i \) lies in the wake of turbine \( j \), the distance between the two along the wind direction \( \theta \) is
\begin{align*}
d_{ij} = \left| (\rho_i - \rho_j)\cos\theta + (\sigma_i - \sigma_j)\sin\theta \right|.
\end{align*}

Restricting the sum in \eqref{eq:veldeficit2} to the turbines \( j \) whose wake contains turbine \( i \), i.e., to those with \( \beta_{ij} < \alpha \), gives
\begin{align}\label{eq5}
\delta s_i = \sqrt{ \sum_{\substack{j = 1 \\ j \ne i,\, \beta_{ij} < \alpha}}^{N} \left( \delta s_{ij} \right)^2 }.
\end{align}
Equation~\eqref{eq5} shows that \( \delta s_i \) depends on the wind direction \( \theta \) and on the position \( (\rho_i, \sigma_i) \) of turbine \( i \).

The wake deficit enters the energy calculation through the Weibull scale parameter. For a direction bin $\theta$, let the free-stream speed be
\[
S_\theta \sim \mathrm{Weibull}\!\left(\zeta(\theta),\psi(\theta)\right).
\]
In the Jensen benchmark, the combined deficit $\delta s_i(\theta)$ acts as a fixed multiplicative reduction within the bin, so the local inflow at turbine $i$ is
\[
S_{i,\theta}=\left[1-\delta s_i(\theta)\right]S_\theta .
\]
Since $aS\sim\mathrm{Weibull}(\zeta,a\psi)$ for any constant $a>0$ when $S\sim\mathrm{Weibull}(\zeta,\psi)$, the shape parameter is unchanged and the scale parameter becomes
\begin{align}
\psi_i(\theta) &= \psi(\theta)\left[1-\delta s_i(\theta)\right],
\quad i=1,2,\ldots,N. \label{eq:weibull-scale}
\end{align}
This relation holds under the assumption that the wake acts as a deterministic speed multiplier in each direction bin. More detailed wake and turbulence models would change the whole speed distribution and not only its scale.

""")

# ------------------------------------------------------------------ power model
M = maths(r"\subsection{The Power Model}", r"\subsection{Constraint Handling Mechanism}")
assert len(M) == 5, len(M)
replace(r"\subsection{The Power Model}", r"\subsection{Constraint Handling Mechanism}", r"""\subsection{Power Model}
\label{sec:powermodel}

We use the simplified power curve of the Kusiak--Song benchmark~\cite{Kusiak2010},
""" + M[0] + r"""
where $s_{\text{cut-in}}$ and $s_{\text{rated}}$ are the cut-in and rated speeds, and $\lambda'$ and $\eta$ define the linear part of the curve. The benchmark does not include a cut-out speed, and we use the model exactly as written. A linear curve is only a rough approximation of a real turbine, so in Section~\ref{sec:powercurve} we re-evaluate all layouts with a cubic power curve, with and without a 25~m/s cut-out.

For a turbine at \( (\rho, \sigma) \) and wind direction \( \theta \), the expected power output is
""" + M[1] + r"""
The layout problem is then to maximize the total expected power of the farm subject to the spacing and boundary constraints of assumptions~6 and~8:
""" + M[2] + r"""
where $E(P_i)$ is the expected power output of the $i$th turbine. Substituting \eqref{eq7} and the Weibull density into \eqref{eq8} gives
""" + M[3] + r"""

Following Kusiak and Song~\cite{Kusiak2010}, we evaluate these integrals with a Riemann sum~\cite{Stroock1999}. The wind direction is divided into bins $\theta_0 = 0^\circ, \theta_1, \ldots, \theta_{N_\theta+1} = 360^\circ$, and the speed range between $s_{\text{cut-in}}$ and $s_{\text{rated}}$ into bins $s_0 = s_{\text{cut-in}}, s_1, \ldots, s_{N_s+1} = s_{\text{rated}}$, all of equal width. The expected power of the $i$th turbine is then approximated by
""" + M[4] + r"""
where $N_s$ and $N_\theta$ are the numbers of speed and direction intervals and $\omega_{l-1}$ is the probability that the wind blows from the $(l-1)$th direction bin.

\paragraph{Benchmark objective and energy}
Like the earlier studies, we keep the bin-width factor $(\theta_l-\theta_{l-1})$ in the sum. The objective values in our tables are therefore benchmark values and not electrical power bounded by $N P_{\text{rated}}$. Because all 24 direction bins are $15^\circ$ wide and $\sum_l\omega_{l-1}=1$, the objective equals 15 times the expected farm power in kW, and the annual energy production follows as $\mathrm{AEP}\,[\mathrm{MWh}] = (\text{objective}/15)\times 8.76$. For example, the wake-free objective of a single turbine under Wind Data Set~I is $14045.74$, which corresponds to an expected power of $936.38$~kW, a capacity factor of $0.624$ and an AEP of $8202.7$~MWh. We confirmed this conversion with the independent evaluator described in Section~\ref{sec:validation}.


""")

# ------------------------------------------------------------------ constraint handling
M = maths(r"\subsection{Constraint Handling Mechanism}", r"\section{\textbf{Previously Developed")
assert len(M) == 3, len(M)
replace(r"\subsection{Constraint Handling Mechanism}", r"\section{\textbf{Previously Developed", r"""\subsection{Constraint Handling}
\label{sec:constraints}

Every layout must satisfy the circular boundary and the minimum spacing. In all our runs, both constraints are handled by two simple operations, a box repair of the coordinates and a large additive penalty in the minimized objective; no other repair operator is used.

\subsubsection{Box Repair}
After every position update, each coordinate is clipped to the square that encloses the farm,
""" + M[0] + r"""
This keeps the decision variables bounded, but it does not guarantee $\rho_i^2+\sigma_i^2\le r^2$; the circular boundary is enforced by the penalty below.

\subsubsection{Constraint Violations}
The spacing and boundary constraints are $d_{ij}=\sqrt{(\rho_i-\rho_j)^2+(\sigma_i-\sigma_j)^2}\ge d_{\min}=8R$ for all $i\ne j$ and $\rho_i^2+\sigma_i^2\le r^2$ for all $i$. We measure their violations by
""" + M[1] + r"""
and call a layout feasible if $g^{\rm s}_{ij}\le0$ and $g^{\rm b}_i\le0$ for all $i$ and $j$.

\subsubsection{Penalized Objective}
All penalty-based optimizers minimize the penalized wake loss
""" + M[2] + r"""
where $N E_{\rm ideal}$ is the wake-free objective, so that the first term is the wake loss reported in our tables. Any violation larger than about $10^{-7}$ (m or m$^2$) makes the penalty larger than the largest possible wake loss. Comparing $F_p$ values therefore gives a feasibility-first selection: a feasible layout is always preferred to an infeasible one, two feasible layouts are compared by their wake loss, and two infeasible layouts by their total penalty. Box repair and $F_p$ are used in every update, comparison, sort and greedy selection of Algorithms~\ref{alg:lx-ssa} and~\ref{alg:lx-ssa-wflop}. Section~\ref{sec:boundary} shows that a radial projection onto the circle, used instead of box clipping, changes the results noticeably, so we keep the same boundary rule for all methods.

""")

# ------------------------------------------------------------------ LX-SSA preliminaries
M = maths(r"\section{\textbf{Previously Developed", r"\section{\textbf{Application of LX-SSA")
assert len(M) == 4, len(M)
ALG1 = flt("alg:lx-ssa")
replace(r"\section{\textbf{Previously Developed", r"\section{\textbf{Application of LX-SSA", r"""\section{Preliminaries: Laplacian Salp Swarm Algorithm}
\label{sec:lxssa}

LX-SSA was introduced by Solanki and Deep~\cite{Solanki2023}, and we summarize here only the operators needed to follow and reproduce its use for WFLOP. As in SSA~\cite{Mirjalili2017}, the sorted population is split into a leader half and a follower half, and the best solution found so far is kept as the food source $H$. LX-SSA keeps these SSA dynamics and adds a Laplace-based candidate with greedy selection to the follower update~\cite{Solanki2023,Deep2007}.

A leader salp $i\le N_p/2$ updates its $j$th coordinate around the food source as
""" + M[0] + r"""
where $r_2,r_3\sim U(0,1)$ and
\[
r_1=2\exp\!\left[-\left(\frac{4l}{L}\right)^2\right],
\]
with iteration $l$ and maximum iteration $L$. A follower salp $i>N_p/2$ first forms the standard SSA position from itself and the preceding salp,
""" + M[1] + r"""

In addition, LX-SSA forms a second follower candidate
""" + M[2] + r"""
where $x_j^i$ is the current coordinate of the follower and
""" + M[3] + r"""
with $z\sim U(0,1)$. This is the inverse-distribution sample of a Laplace variable with location $\phi$ and scale $\chi$ (mean $\phi$, variance $2\chi^2$). The location $\phi$ centers the step, and the scale $\chi>0$ controls its spread: a small $\chi$ keeps the moves short, while a large $\chi$ makes long exploratory moves more likely. We use $\phi=0$ and $\chi=1$ in all experiments, so that $\gamma$ is symmetric about zero and $P(|\gamma|>1)=e^{-1}\approx0.37$. Most candidates therefore move towards $H$ ($\gamma>0$) or away from it ($\gamma<0$) by a fraction of their current distance, while the heavy tail occasionally produces a move beyond $H$ or far away from it. These values were not tuned for WFLOP. The follower keeps the better of its two candidates. Since each follower evaluates both candidates before the population is re-evaluated, one LX-SSA iteration costs $N_p/2\times2+N_p=2N_p$ objective calls, twice the cost of an SSA iteration, and all comparisons in Section~\ref{sec:results} are therefore made at equal numbers of calls.

Algorithm~\ref{alg:lx-ssa} summarizes the search. Whether one candidate is better than another is decided by the objective and, for WFLOP, by the feasibility rules of Section~\ref{sec:constraints}. Because SSA and LX-SSA differ only in the Laplace-based follower candidate, a comparison between the two isolates the effect of this mechanism.

""" + ALG1 + "\n")

# ------------------------------------------------------------------ application of LX-SSA to WFLOP
ALG2 = flt("alg:lx-ssa-wflop")
i = t.index(r"\begin{table*}[htbp]" + "\n" + r"\centering" + "\n" + r"\caption{Wind Data Sets I and II")
j = t.index(r"\end{table*}", i) + len(r"\end{table*}")
WIND = t[i:j] + "\n"
replace(r"\section{\textbf{Application of LX-SSA", r"\section{Hybrid LX-SSA-VNS Algorithm}", r"""\section{Application of LX-SSA to WFLOP}
\label{sec:lxssa-wflop}

We apply the published LX-SSA~\cite{Solanki2023} to WFLOP without changing the algorithm. For a prescribed number of turbines $N$, candidate $j$ is the vector of all turbine coordinates,
\[
\mathbf{x}_j=(\rho_j^1,\sigma_j^1,\rho_j^2,\sigma_j^2,\ldots,\rho_j^N,\sigma_j^N)\in\mathbb{R}^{2N},
\]
so the optimizer searches only the turbine positions, and $N$ stays fixed during a run.

Algorithm~\ref{alg:lx-ssa-wflop} lists the steps of the implementation. The candidate coordinates are updated first and box-repaired; the wake deficits and Weibull scale factors are then recomputed for every direction bin, the constraint violations are measured, and the penalized objective $F_p$ is evaluated. All comparisons use $F_p$, which gives the feasibility-first selection of Section~\ref{sec:constraints}. The order of these steps matters because the energy of each turbine depends on the wake field of the whole layout.

""" + ALG2 + "\n" + WIND + "\n\n")

# ------------------------------------------------------------------ hybrid section
ALG3 = flt("alg:hybrid")
replace(r"\section{Hybrid LX-SSA-VNS Algorithm}", r"\section{Experimental Setup}", r"""\section{Hybrid LX-SSA-VNS Algorithm}
\label{sec:hybrid}
Swarm methods such as LX-SSA spread their search over many regions of the layout space, but once the swarm has contracted around the food source, most of the remaining evaluations bring only small gains. VNS~\cite{Mladenovic1997,Hansen2001} behaves the other way round: it refines a single layout systematically, alternating random perturbations of increasing size with a complete local search, but its result depends on the starting point. We combine the two in a two-phase hybrid, LX-SSA-VNS, in which LX-SSA~\cite{Solanki2023} explores the layout space and the basic VNS refines the best layout it finds.

\subsection{Phase 1: Exploration by LX-SSA}
For a budget of $B$ objective calls and a split fraction $\rho\in(0,1)$, LX-SSA (Algorithm~\ref{alg:lx-ssa-wflop}) runs with population $N_p$ for
\[
T_1=\operatorname{round}\!\left(\frac{\rho B-N_p}{2N_p}\right)
\]
iterations, since one LX-SSA iteration costs $2N_p$ calls. Its parameter $r_1$ decreases over these $T_1$ iterations, so the swarm itself moves from exploration to exploitation within Phase~1. The food source $H$ at the end of Phase~1 becomes the starting point of Phase~2.

\subsection{Phase 2: Refinement by Basic VNS}
Phase~2 is the basic VNS of Mladenovi\'c and Hansen~\cite{Mladenovic1997,Hansen2001} in the continuous form of Mladenovi\'c et al.~\cite{Mladenovic2008}, applied to the penalized objective $F_p$ in~\eqref{eq:penalty}. It uses $k_{\max}=5$ nested $\ell_\infty$ neighborhoods of the whole layout,
\[
\mathcal N_k(x)=\{y:\rho_{k-1}<\lVert y-x\rVert_\infty\le\rho_k\},\qquad \rho_k=0.1\,k\,r .
\]
In the shaking step, a point is drawn uniformly from $\mathcal N_k(x)$, so that all $2N$ coordinates move, and box repair is applied. The shaken layout is then improved by a complete best-improvement local search. Since the Jensen top-hat wake makes $F_p$ discontinuous, we use a derivative-free compass search: it evaluates all $4N$ moves $\pm h$ along the coordinate axes, moves to the best improving one, and halves $h$ (starting from $0.05r$) when no move improves, until $h<10^{-3}r$. If the resulting local minimum is better than the incumbent, it is accepted and $k$ is reset to 1; otherwise $k$ is increased by one, and it returns to 1 after $k_{\max}$.

Phase~2 first applies the local search to $H$ and then repeats shaking, local search and neighborhood change until the budget $B$ is used up. It returns the best layout evaluated. Box repair and the penalized objective are used exactly as in LX-SSA, so the feasibility-first rules of Section~\ref{sec:constraints} hold in both phases. Algorithm~\ref{alg:hybrid} gives the complete procedure.

""" + ALG3 + r"""
\subsection{Parameters and Cost}
The hybrid uses $N_p=30$, $\phi=0$, $\chi=1$ and the VNS settings above, and it adds a single parameter, the split $\rho$. With the default $\rho=0.5$, LX-SSA receives 3,030 of the 6,030 calls ($T_1=50$) and VNS the remaining 3,000; the effect of $\rho$ is studied in Section~\ref{sec:ablation}. The cost of both phases is dominated by the objective evaluations, so for a budget of $B$ calls the hybrid has the same $\mathcal O(B\,N^2N_\theta N_s)$ cost as the other methods. Because the random-number generator is seeded once and Phase~1 draws the initial population first, the hybrid starts from the same initial population as all other methods for a given seed, which keeps our comparisons paired.

""")

# ------------------------------------------------------------------ experimental setup
replace(r"\section{Experimental Setup}", r"\section{Numerical Results}", r"""\section{Experimental Setup}
\label{sec:setup}

\subsection{Test Problems}
We use two wind data sets (Table~\ref{tab:wind_data_I}). In Wind Data Set~I, the Weibull parameters are the same in all directions ($\zeta=2$, $\psi=13$) and the wind is strongly directional: $\omega_0=\omega_{23}=0$, $\omega_5=0.2$, $\omega_6=0.6$, and $\omega_l=0.01$ for all other intervals. Wind Data Set~II keeps $\zeta=2$, but the scale parameter varies with direction, from $\psi=2.6$ (255--270$^\circ$) to $\psi=10$ (180--195$^\circ$), and the wind is spread over more directions. The wind direction is divided into $N_\theta=24$ bins of $15^\circ$, and the speed range from $s_{\text{cut-in}}$ to $s_{\text{rated}}$ into $N_s=21$ intervals of 0.5~m/s.

Each data set is combined with three circular farms of radius $r=500$, 750 and 1000~m, and the turbine count is prescribed for each run: $N=2,\ldots,10$ for 500~m, $N=2,\ldots,12$ for 750~m and $N=2,\ldots,15$ for 1000~m. This gives $2\times(9+11+14)=68$ test cases. The turbine and wake parameters are:
\begin{itemize}
    \item rotor radius $R=38.5$~m and minimum spacing $4D=8R=308$~m;
    \item cut-in speed $s_{\text{cut-in}}=3.5$~m/s, rated speed $s_{\text{rated}}=14$~m/s and rated power $P_{\text{rated}}=1500$~kW;
    \item linear power-curve parameters $\lambda'=140.86$ and $\eta=-500$;
    \item thrust coefficient $C_T=0.8$ and wake spreading constant $K=0.075$;
    \item constraint penalty $\lambda=10^{10}$ in~\eqref{eq:penalty}.
\end{itemize}

\subsection{Methods Compared}
We compare the proposed hybrid with six reference methods. All of them are evaluated with the objective~\eqref{eq:penalty}, the bounds $[-r,r]^{2N}$ and the same budget:
\begin{itemize}
\item \emph{LX-SSA-VNS} (proposed): Algorithm~\ref{alg:hybrid} with $\rho=0.5$, i.e., 50 LX-SSA iterations (3,030 calls) followed by 3,000 calls of basic VNS.
\item \emph{LX-SSA}~\cite{Solanki2023}: Algorithm~\ref{alg:lx-ssa-wflop} with $\phi=0$, $\chi=1$, population 30 and 100 iterations.
\item \emph{SSA}~\cite{Mirjalili2017}: the same leader and follower dynamics without the Laplace candidate; population 30 and 200 iterations.
\item \emph{PSO}: inertia weight $w=0.7$ and acceleration coefficients $c_1=c_2=2$, with unbounded velocities and positions clipped to the bounds; population 30 and 200 iterations.
\item \emph{DE}: DE/rand/1/bin with $F=0.5$ and $CR=0.9$; population 30 and 200 iterations.
\item \emph{VNS}: the basic VNS of Mladenovi\'c and Hansen~\cite{Mladenovic1997,Hansen2001} in the continuous form of Mladenovi\'c et al.~\cite{Mladenovic2008}. It is exactly Phase~2 of Algorithm~\ref{alg:hybrid} ($\ell_\infty$ shaking of the whole layout with $k_{\max}=5$ and $\rho_k=0.1kr$, and a complete best-improvement compass search), started from the best member of the initial population instead of the LX-SSA food source.
\item \emph{Multistart SLSQP (MS-SLSQP)}, a gradient-based method: the sequential least-squares quadratic programming solver of SciPy~\cite{Kraft1988} minimizes the wake loss with the spacing and boundary constraints imposed explicitly and forward-difference gradients (1-m step). It starts from the members of the initial population in the order of their $F_p$ values and restarts until the budget is used up.
\end{itemize}
SSA, LX-SSA, PSO and DE use our original implementation. We implemented VNS, MS-SLSQP and the hybrid for this study; the hybrid calls the original LX-SSA code in its first phase. For the ablation in Section~\ref{sec:ablation}, we also use SSA-VNS, which replaces LX-SSA by SSA in Phase~1 (100 SSA iterations, 3,030 calls) and is otherwise identical to LX-SSA-VNS.

\subsection{Budget, Seeds and Pairing}
Every run has a budget of exactly 6,030 objective calls, counted by a wrapper around the objective function. The count includes the initial population and, for MS-SLSQP, the calls used for the finite-difference gradients. For LX-SSA this corresponds to 100 iterations, because each follower evaluates two candidates before the population is re-evaluated ($30+100\times60$ calls); SSA, PSO and DE reach the same count in 200 iterations ($30+200\times30$), and LX-SSA-VNS in $3{,}030+3{,}000$ calls. Each method is run with the seeds 1--30. Every optimizer seeds the random-number generator and draws its initial population first, so runs with the same seed start from the \emph{same} initial population for all methods, and the runs of a test case are paired by seed.

\subsection{Performance Measures and Statistical Analysis}
For each run we record the final objective and the wake loss (as a percentage of the wake-free objective), whether the final layout is feasible (spacing and boundary satisfied within $10^{-6}$~m), the convergence curve, i.e., the best feasible objective after every 30 calls, and the wall-clock time. Means and standard deviations are computed over the feasible runs.

Within each test case, the 30 seed-paired runs are compared by a run-level Friedman test and by two-sided Wilcoxon signed-rank tests of LX-SSA-VNS against each of the other six methods. The six $p$ values of each case are adjusted with Holm's procedure~\cite{Holm1979}, and the effect sizes are reported as matched-pairs rank-biserial correlations. For the ranking, a feasible run beats an infeasible one, and infeasible runs are compared by their spacing shortfall. Across the 68 cases, the methods are ranked within each case by their mean feasible objective and compared by a case-level Friedman test with the Iman--Davenport correction and Holm-adjusted post hoc tests on the average ranks, following Dem\v{s}ar~\cite{Demsar2006}.

\subsection{Implementation Verification}
\label{sec:validation}
We vectorized the objective function and checked it against the original implementation; the relative difference is below $10^{-6}$ on random feasible and infeasible layouts. It also reproduces the objective values of 720 archived final layouts from an earlier run of the original code, with a maximum absolute deviation of $8.7\times10^{-11}$. Re-running the original SSA, LX-SSA, PSO and DE code with the archived settings reproduces the archived result distributions: 22 of 24 Mann--Whitney comparisons give $p\ge0.05$, and the means differ by at most 0.21\%. This also confirms the Laplace parameters $\phi=0$ and $\chi=1$.

\subsection{Computing Environment}
All runs were carried out in a Linux container with an Intel Xeon processor (2.80~GHz, 4 virtual cores), four runs in parallel, using Python~3.11.15, NumPy~2.4.6 and SciPy~1.17.1. Each run was timed individually. The scripts, the per-run results (including the final coordinates and convergence curves) and the code that generates the figures are provided with the revision.

""")

# ------------------------------------------------------------------ numerical results
F = {k: flt(k) for k in ["tab:friedman68", "fig:avgranks", "tab:wtl", "fig:wakeloss", "tab:res-1-500", "tab:res-1-750",
                         "tab:res-1-1000", "tab:res-2-500", "tab:res-2-750", "tab:res-2-1000", "fig:feasibility",
                         "fig:conv-mid", "fig:conv-max", "fig:box", "tab:ablation", "fig:ablation", "tab:split",
                         "fig:layouts", "tab:hr-site", "fig:hr-conv", "fig:hr-layouts"]}
replace(r"\section{Numerical Results}", r"\section{Sensitivity to the Modelling Assumptions}", r"""\section{Numerical Results}
\label{sec:results}
\label{sec:equalbudget}
In this section, we report 14,280 runs (68 test cases, seven methods and 30 seed-paired runs, all with 6,030 objective calls). We then present the component ablation of the hybrid, which adds 2,040 runs of SSA-VNS and 720 runs with other budget splits, and the Horns Rev~1 site case. The detailed results of every case are given in Tables~\ref{tab:res-1-500}--\ref{tab:res-2-1000}.

\subsection{Overall Comparison}
\label{sec:overall}
""" + F["tab:friedman68"] + F["fig:avgranks"] + F["tab:wtl"] + r"""
Table~\ref{tab:friedman68} and Fig.~\ref{fig:avgranks} summarize the comparison over all cases. The Friedman test over the 68 cases rejects the hypothesis that the seven methods perform equally ($\chi^2_F=168.4$, 6 d.f., $p=9.9\times10^{-34}$; Iman--Davenport $F_F=47.1$). LX-SSA-VNS has the best average rank (1.89) and alone attains the highest mean objective in 31 cases. VNS follows with 2.85 (13 cases), and then come SSA (3.85), PSO (4.14), MS-SLSQP (4.44, 12 cases), LX-SSA (4.76) and DE (6.07). In the Holm-adjusted post hoc tests, LX-SSA-VNS ranks significantly better than every other method, including VNS ($p_{\rm Holm}=0.010$); for the remaining five methods, $p_{\rm Holm}\le2.3\times10^{-7}$.

The run-level tests in Table~\ref{tab:wtl} show where these differences come from. Against LX-SSA, the hybrid is significantly better in 39 cases and never worse, with rank-biserial correlations between 0.53 and 1.00 (median 0.74). Of the 29 cases without a significant difference, 22 have $N\le7$, where LX-SSA already comes close to the best layouts. Against SSA, PSO and DE, the hybrid is significantly better in 32, 36 and 52 cases, respectively, and never worse. MS-SLSQP is beaten in 36 cases but wins four dense cases (Data Set~I, 500~m, $N=8$--10, and Data Set~II, 750~m, $N=12$), where its explicit handling of the constraints pays off. The comparison with VNS is the closest one: the two methods are statistically indistinguishable in 60 cases, the hybrid is better in six cases with $N=5$ to 9 ($r_{rb}$ from 0.51 to 0.63), and VNS is better in the two three-turbine cases of the 500-m farm. Hybridization thus moves LX-SSA from the sixth to the first average rank. Its margin over VNS, which forms its own second phase, is small in individual cases but consistent over the whole benchmark.

\subsection{Solution Quality by Farm Size and Turbine Count}
""" + F["fig:wakeloss"] + r"""
Fig.~\ref{fig:wakeloss} and Tables~\ref{tab:res-1-500}--\ref{tab:res-2-1000} show how the wake loss depends on the case. It grows with $N$, falls as the farm becomes larger, and is higher for the broader wind rose of Data Set~II. For $N\le3$, all methods stay within 0.3\% of the wake-free objective. As $N$ grows, the curve of the hybrid separates from that of LX-SSA and joins the curves of VNS and MS-SLSQP. For example, for 15 turbines in the 1000-m farm (Data Set~I), the mean wake loss is 1.78\% for LX-SSA-VNS, 1.99\% for VNS, 2.20\% for MS-SLSQP, 2.80\% for SSA and 3.32\% for LX-SSA; for 12 turbines in the 750-m farm (Data Set~II), the corresponding values are 5.95\%, 6.17\%, 5.60\%, 6.79\% and 7.34\%. In the largest case of every farm, the hybrid reduces the mean wake loss of LX-SSA by 1.4 to 1.9 percentage points. It has the lowest mean wake loss in both 15-turbine cases of the 1000-m farm, while MS-SLSQP has the lowest value in the other four largest cases.

""" + "".join(F[k] for k in ["tab:res-1-500", "tab:res-1-750", "tab:res-1-1000", "tab:res-2-500", "tab:res-2-750",
                             "tab:res-2-1000"]) + r"""
\subsection{Feasibility}
""" + F["fig:feasibility"] + r"""
Fig.~\ref{fig:feasibility} shows the share of runs that end with a feasible layout. LX-SSA-VNS is feasible in 99.9\% of all runs; only one run in each of the two 500-m ten-turbine cases fails. The other methods reach 96.9\% (LX-SSA), 97.6\% (SSA), 86.0\% (PSO), 71.4\% (DE), 100\% (VNS) and 99.9\% (MS-SLSQP). In the densest cases (500~m, $N=10$), LX-SSA is feasible in only 13 of 30 runs, and the complete local search of the VNS phase repairs almost all of the layouts it receives. The feasibility counts of the penalty-based methods are the same for both wind data sets, because as long as every candidate is infeasible, the penalty dominates $F_p$.

\subsection{Convergence Behavior}
""" + F["fig:conv-mid"] + F["fig:conv-max"] + r"""
Figs.~\ref{fig:conv-mid} and~\ref{fig:conv-max} show the effect of the hybridization directly. Up to 3,030 calls, the hybrid runs a shortened LX-SSA and is therefore behind the full LX-SSA run. After the switch, its curve drops steeply, whereas LX-SSA and SSA stagnate after about 3,000 calls. Over all feasible runs, the VNS phase removes on average 39.0\% (median 31.9\%) of the wake loss left at the end of the LX-SSA phase, and it turns 118 runs that were still infeasible at the switch into feasible ones. VNS alone reaches low wake losses earlier because it spends its whole budget on refinement, but the hybrid catches up after the switch and ends with a lower mean wake loss in five of the six largest cases. MS-SLSQP is the fastest method to find a feasible layout.

\subsection{Distribution of Results}
""" + F["fig:box"] + r"""
Fig.~\ref{fig:box} shows the spread of the final wake loss over the 30 runs for moderate turbine counts. LX-SSA-VNS has the lowest median in all three cases of Data Set~II and the second-lowest in the 500-m case of Data Set~I, and its interquartile range is smaller than that of LX-SSA in five of the six cases. DE has the highest median in five cases.

\subsection{Component Ablation}
\label{sec:ablation}
""" + F["tab:ablation"] + F["fig:ablation"] + F["tab:split"] + r"""
To understand where the improvement of the hybrid comes from, we carry out a component ablation. The hybrid has three ingredients: a swarm phase, the Laplace step inside that swarm phase, and the VNS phase. We therefore compare SSA, LX-SSA, VNS, SSA-VNS and LX-SSA-VNS on all 68 cases with the same seeds and budget (Table~\ref{tab:ablation} and Fig.~\ref{fig:ablation}). The Friedman test over the five variants is highly significant ($\chi^2_F=156.5$, 4 d.f., $p=8.1\times10^{-33}$).

Adding the VNS phase is the most important step. LX-SSA-VNS beats LX-SSA in 39 cases and SSA-VNS beats SSA in 40 cases, without a single loss, and the mean wake loss falls by 0.43 and 0.33 percentage points, respectively. Starting VNS from the food source of a swarm, rather than from the best initial point, helps moderately: both two-phase variants have a better average rank than VNS (2.25 and 1.62 against 2.90), LX-SSA-VNS beats VNS in five cases and loses two, and SSA-VNS beats it in ten cases and loses one.

The Laplace step of LX-SSA, on the other hand, does not help. Without VNS, LX-SSA is worse than SSA in seven cases and never better. Inside the hybrid, LX-SSA-VNS and SSA-VNS are statistically indistinguishable in 66 cases, SSA-VNS is better in the remaining two, and SSA-VNS has the better average rank (1.62 against 2.25). The split of the budget between the two phases matters little (Table~\ref{tab:split}). The 50\% split gives the highest mean in eight of the 12 cases, but only one difference is significant (Data Set~I, 1000~m, $N=15$, where the 50\% split is better than the 25\% split, $p_{\rm Holm}=0.024$).

Taken together, the ablation shows that the gain of LX-SSA-VNS comes from combining a swarm exploration phase with a complete VNS refinement phase, and not from the particular Laplace mechanism of LX-SSA, since a standard SSA first phase performs at least as well. The counts against VNS differ slightly from Table~\ref{tab:wtl} because the Holm correction is applied to a different family of comparisons.

\subsection{Optimized Layouts}
""" + F["fig:layouts"] + r"""
Fig.~\ref{fig:layouts} compares, for the largest turbine count of each farm, the best LX-SSA-VNS layout with the best layout found by any method. The best layout overall comes from LX-SSA-VNS in four cases, from MS-SLSQP in one (Data Set~I, 500~m) and from VNS in one (Data Set~II, 750~m). In all of them, most turbines lie close to the boundary, which maximizes their mutual distances, while the remaining turbines sit in the interior at positions staggered with respect to the dominant wind directions.

\subsection{Computational Cost}
\label{sec:runtime}
At equal numbers of objective calls, all methods need almost the same time: 0.37~ms per call for VNS, 0.39~ms for LX-SSA-VNS, 0.40~ms for LX-SSA and SSA, 0.41~ms for PSO, 0.47~ms for DE and 0.49~ms for MS-SLSQP. The objective evaluation dominates the cost, so neither the Laplace step nor the VNS bookkeeping adds a measurable overhead. One 6,030-call run of LX-SSA-VNS takes 1.7--3.5~s on average, depending on $N$.

\subsection{Horns Rev 1 Site Case}
\label{sec:hornsrev-site}
""" + F["tab:hr-site"] + F["fig:hr-conv"] + F["fig:hr-layouts"] + r"""
To test the methods closer to a real project, we modeled the Horns Rev~1 offshore farm with its real turbine, wind climate and site outline, all taken from the PyWake distribution~\cite{PyWake}. The model contains the installed layout of 80 Vestas V80 turbines (rotor diameter 80~m, 2~MW) with the V80 power and thrust curves, the measured 12-sector wind climate, and the parallelogram spanned by the installed turbines as the site boundary. We keep the Jensen top-hat wake, with $k=0.04$ as commonly used offshore, and root-sum-square superposition, and we use the tabulated power curve with a 25~m/s cut-out, a speed-dependent thrust coefficient at the free-stream speed, 5$^\circ$ direction bins and 1~m/s speed bins. For the installed layout, our model gives 664.6~GWh/yr and a wake loss of 10.95\%, against 662.5~GWh/yr and 10.96\% from PyWake's Jensen model with the same settings, a difference below 0.3\%. We optimized the 16-turbine block that PyWake uses as an example (30 seeds) and the complete 80-turbine farm (10 seeds) with all seven methods, 6,030 calls and the $4D$ (320~m) minimum spacing.

The results are given in Table~\ref{tab:hr-site} and Figs.~\ref{fig:hr-conv} and~\ref{fig:hr-layouts}. In the 16-turbine block, LX-SSA-VNS attains the highest mean AEP (137.77~GWh/yr) and is feasible in all 30 runs. It is significantly better than LX-SSA ($p_{\rm Holm}=5.7\times10^{-6}$), SSA ($p_{\rm Holm}=3.1\times10^{-4}$) and MS-SLSQP ($p_{\rm Holm}=2.4\times10^{-4}$), while its difference from VNS (137.67~GWh/yr) is not significant ($p_{\rm Holm}=0.685$). PSO and DE never find a feasible layout inside the parallelogram, and none of the methods improves on the installed $7D$ grid, which reaches 139.51~GWh/yr. The full farm is a much harder problem. With 160 variables and 3,160 spacing constraints, no run of LX-SSA-VNS, LX-SSA, SSA, PSO, DE or VNS ends with a feasible layout. MS-SLSQP, which handles the constraints explicitly, is feasible in all ten runs, but its mean AEP (653.9~GWh/yr) remains 1.6\% below that of the installed layout. Realistic farm sizes therefore call for larger budgets, feasible initialization or explicit constraint handling, and a regular installed layout is a strong reference.

""")

# ------------------------------------------------------------------ sensitivity, complexity, limitations, conclusion
G = {k: flt(k) for k in ["tab:robust-final", "tab:capacity", "tab:spacing-authors"]}
replace(r"\section{Sensitivity to the Modelling Assumptions}", r"\appendices", r"""\section{Sensitivity to the Modeling Assumptions}
\label{sec:robustness}
\subsection{Power Curve and Wake Model}
\label{sec:powercurve}
\label{sec:gaussian}
""" + G["tab:robust-final"] + r"""
To see how far our conclusions depend on the benchmark power curve and wake model, we re-evaluated every feasible final layout, without re-optimization, with three alternative models:
\begin{enumerate}
\item a cubic power curve~\cite{Carrillo2013}, $f_{\rm c}(s)=P_{\text{rated}}(s^3-s_{\text{cut-in}}^3)/(s_{\text{rated}}^3-s_{\text{cut-in}}^3)$ between the cut-in and rated speeds;
\item the same curve with a 25~m/s cut-out;
\item the Gaussian wake model of Bastankhah and Port\'e-Agel~\cite{Bastankhah2014},
\begin{align*}
\delta_{ij}&=\Bigl(1-\sqrt{1-C_T/(8\sigma^2/D^2)}\Bigr)\exp\!\bigl(-r_\perp^2/(2\sigma^2)\bigr),\\
\sigma/D&=k^*x/D+0.2\sqrt{\beta},
\end{align*}
with $\beta=\tfrac12(1+\sqrt{1-C_T})/\sqrt{1-C_T}$ and $k^*=0.04$, evaluated at the rotor center with the same superposition and Weibull scaling as the benchmark.
\end{enumerate}

Table~\ref{tab:robust-final} shows two kinds of sensitivity. First, the absolute energy depends strongly on the model: the linear benchmark curve overstates the expected power relative to the cubic curve by 17.7\% for Data Set~I and 36.9\% for Data Set~II (21.5\% and 37.1\% with the cut-out). The energy values in this paper should therefore be read as benchmark values and not as AEP predictions for a specific turbine. Second, the ranking of the algorithms is stable under the cubic curve (mean Kendall $\bar\tau=0.97$), with LX-SSA-VNS (1.95) and VNS (2.85) still the two best methods, but it changes more under the Gaussian model ($\bar\tau=0.54$; the best method is the same in 43\% of the cases), because layouts optimized for the sharp Jensen wake cone are judged differently by a smooth Gaussian deficit. LX-SSA-VNS nevertheless keeps the best average rank under the Gaussian model (2.82), followed by PSO (3.16), MS-SLSQP (3.43) and VNS (3.76).

\subsection{Minimum Spacing}
\label{sec:spacing}
""" + G["tab:capacity"] + G["tab:spacing-authors"] + r"""
The $4D$ spacing is inherited from the benchmark literature, and we examined its influence in three ways. First, Table~\ref{tab:capacity} gives the largest turbine count for which a multi-start packing search constructed a feasible layout. These are constructive lower bounds, and the 500-m values agree with the known optimal circle-packing results. With $5D$ spacing, the 500-m farm holds at most 8 turbines and with $6D$ only 7, so the 500-m cases with $N\ge9$ (at $5D$) or $N\ge8$ (at $6D$) would be infeasible.

Second, the optimized layouts often lie on the $4D$ constraint. Of the 13,293 feasible final layouts of the seven methods, 48.6\% also satisfy $5D$ and 35.8\% satisfy $6D$; for $N=6$--10 these shares fall to 30.1\% and 11.6\%, and for $N=11$--15 to 3.0\% and 0.1\%. Layouts optimized for $4D$ therefore cannot be transferred to a site with larger spacing requirements without re-optimization.

Third, we re-optimized six representative cases with LX-SSA and SSA at $5D$ and $6D$ (Table~\ref{tab:spacing-authors}). In the 750-m eight-turbine cases, the mean LX-SSA objective falls from $4D$ to $6D$ by 0.95\% (Data Set~I) and 0.22\% (Data Set~II), and the wake loss rises from 0.88\% to 1.85\% and from 3.18\% to 3.47\%. The spacing rule does not reverse the comparison between LX-SSA and SSA.

\subsection{Boundary Handling}
\label{sec:boundary}
All methods clip the coordinates to the bounding square and penalize turbines outside the circle (Section~\ref{sec:constraints}). An otherwise identical SSA that instead projects every turbine radially onto the circle attains higher mean objectives than the original SSA in six representative cases at 3,030 calls, by 26 to 649 objective units (Mann--Whitney $p<10^{-4}$ in every case). The boundary rule therefore has a noticeable effect, and we keep it fixed for all methods.

\section{Computational Complexity}
\label{sec:complexity}
Let $N_t$ be the number of turbines, $N_\theta$ the number of wind-direction bins and $N_s$ the number of wind-speed bins. A direct wake evaluation covers the pairwise turbine interactions for every direction and speed bin, so one objective evaluation costs
\[
C_{\mathrm{obj}}=\mathcal{O}\!\left(N_t^2N_\theta N_s\right),
\]
and a run with a budget of $B$ objective calls costs $\mathcal{O}(B\,N_t^2N_\theta N_s)$ for every method considered here, including both phases of LX-SSA-VNS. Box repair is $\mathcal{O}(N_t)$ per candidate and the spacing check $\mathcal{O}(N_t^2)$, and the Laplace draw of LX-SSA, the shaking step of VNS and each compass move are $\mathcal{O}(N_t)$. All of these are dominated by the wake evaluation, which agrees with the almost identical measured cost per call (Section~\ref{sec:runtime}).

\section{Limitations}
\label{sec:limits}
Our results should be read within the following limits. The ablation shows that the gain of LX-SSA-VNS comes from the combination of a swarm phase with a complete VNS phase, and not from the Laplace step of LX-SSA: SSA-VNS is statistically equivalent in 66 of the 68 cases, better in two, and has a slightly better average rank. The hybrid is significantly better than VNS over the whole benchmark, but in single cases the two differ significantly in only eight of the 68 cases.

All methods were compared at 6,030 objective calls, and larger budgets could change the ranking; for an 80-turbine farm, none of the penalty-based methods, including the hybrid, found a feasible layout within this budget. The LX-SSA parameters ($\phi=0$, $\chi=1$), the VNS settings and the split $\rho=0.5$ are fixed defaults and were not tuned, and pseudo-gradient and surrogate-based methods are discussed but not benchmarked.

The Jensen model and the linear power curve are benchmark choices. The absolute energy depends strongly on them, the ranking of the methods changes under a Gaussian wake model, and we did not re-optimize the layouts with the alternative models. The $4D$ spacing and the boundary rule also affect the absolute results. Finally, $N$ is prescribed in every run, and terrain, mixed turbine types, collector cables, grid connection, noise and environmental constraints, and economic objectives are outside the scope of this study.

\section{Conclusion}
\label{sec:conclusion}
This work presents LX-SSA-VNS, a two-phase hybrid for the continuous wind farm layout problem in which the Laplacian Salp Swarm Algorithm~\cite{Solanki2023} explores the layout space and the basic variable neighborhood search then refines the best layout found. We evaluated the hybrid on the Jensen--Weibull benchmark with explicit constraint handling against LX-SSA, SSA, PSO, DE, VNS and a multistart SLSQP solver on 68 cases (14,280 seed-paired runs at equal budgets), in a component ablation, and on the Horns Rev~1 offshore farm.

LX-SSA-VNS obtained the best average rank (1.89), was significantly better than every other method over the whole benchmark, including VNS ($p_{\rm Holm}=0.010$), and returned a feasible layout in 99.9\% of the runs. Compared with LX-SSA alone, it was significantly better in 39 of the 68 cases and never worse. The ablation showed that this gain comes from the two-phase design: the VNS phase removes on average 39\% of the wake loss left by the swarm, and starting VNS from a swarm solution is better than starting it from the best initial point. The Laplace step, however, does not contribute, since an SSA first phase performs at least as well. At the Horns Rev~1 site, the hybrid gave the highest mean energy yield for a 16-turbine block, although no method improved on the installed layout within the budget and only the constraint-aware SLSQP solver found feasible layouts for the full 80-turbine farm. The ranking of the methods was stable under a cubic power curve but changed considerably under a Gaussian wake model, where the hybrid nevertheless kept the best average rank.

Despite these results, the hybrid shares the limits of penalty-based search on large farms, and its evaluation is tied to a benchmark wake model. In future work, we will study a switch between the two phases that is triggered by the stagnation of the swarm, investigate which properties of the swarm phase make it a good starting point for VNS, and extend the approach with feasibility-preserving initialization, larger budgets for realistic farm sizes, and higher-fidelity wake and power models.

""")

# ------------------------------------------------------------------ appendix, spelling
t = t.replace(r"""Table~\ref{tab:coords} lists the turbine coordinates (m, farm center at the origin) of the best layout found by any method for the largest turbine count of each farm (Fig.~\ref{fig:layouts}). The final coordinates and convergence curves of all runs and the Horns Rev~1 layouts are provided as machine-readable files with the revision.""",
              r"""Table~\ref{tab:coords} lists the turbine coordinates (in m, with the farm center at the origin) of the best layout found by any method for the largest turbine count of each farm (Fig.~\ref{fig:layouts}). The final coordinates and convergence curves of all runs and the Horns Rev~1 layouts are provided as machine-readable files.""")
t = t.replace("neighbourhood", "neighborhood").replace("Neighbourhood", "Neighborhood")
open(MS, "w").write(t)
print("prose rewritten")

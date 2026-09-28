# R4 — Wind-energy engineering review (wake modeling, AEP, farm design practice, IEA37)

Scope: optA/*.tex (as compiled into MPCE_PSO_VNS.pdf), MPCE_PSO_VNS_supplement.tex, analysis/wflop_model.py, analysis/hornsrev_model.py, analysis/iea37_model.py, and the Horns Rev runs in analysis/mpce_*.csv. The PDFs could not be rendered in this container (no poppler), so the review was done on the LaTeX sources and generated tables/macros that feed the PDF.

## Overall assessment

As a benchmarking paper the revision is careful. It explains the Kusiak–Song objective more honestly than most of the WFLOP literature: the factor 15, the AEP conversion, the negative power just above cut-in, the jump at rated power and the missing cut-out are all stated, and the arithmetic checks out (E_ideal = 14045.74 → 936.38 kW → 8202.7 MWh, reproduced from wflop_model.py). The IEA37 model matches the official calculator. The direction conventions (θ = 270° − θ_met, "towards", CCW from east) are consistent across text, supplement and code.

From a wind-engineering point of view, three things weaken it:
1. The Horns Rev 1 model has a direction-binning bug. The 12 sectors get unequal weights (±17 %) and the total probability is 1.003. Once corrected, the headline result "mean PSO-VNS AEP at 30,030 calls exceeds the installed layout" no longer holds, and the best-run gain shrinks from 0.63 % to 0.48 %.
2. The claim that optimized layouts "exceed the installed Horns Rev 1 block" is stated in the abstract, the contributions and the conclusion without the fairness caveats: 4D vs 7D spacing, an isolated corner block, no cables, loads or costs, and a gain smaller than the model uncertainty.
3. Several known limitations of the benchmark physics are not disclosed: the upstream-wake artefact of the Kusiak–Song cone test, the gaps between the coarse direction bins, constant C_T above rated, the unrealistic wind climate of Data Set I, and turbulence and loads at 4D. The practitioner framing (power-system planners) is also not supported by any result.

None of this invalidates the algorithmic conclusions. It does change what can be claimed about real sites. Recommendation: **major revision** on the wind-energy side.

## Issues

1. **[major] analysis/hornsrev_model.py (drives Table hr-site, 08_beyond.tex)**
   - Quote: `_SEC = (np.round(WD / 30.0).astype(int)) % 12  # nearest sector` … `F_WD = SEC_F[_SEC] / 6.0`
   - Problem: `np.round` rounds half to even, so the 5° bins at 15°, 45°, 75°, … are assigned alternately. Sectors 0°, 60°, …, 300° receive 7 bins and 30°, 90°, …, 330° receive 5 (checked: `np.bincount(_SEC)` = [7 5 7 5 …]). Each sector frequency is therefore scaled by 7/6 or 5/6 and ΣF_WD = 1.0031.
     - The row-aligned sectors 90° and 270° are among those under-weighted.
     - I re-evaluated the stored layouts with 6 bins per sector. Installed: 139.51 → 139.27 GWh/yr. Best PSO-VNS run at 6,030 calls: +0.63 % → +0.48 % over the installed layout. **Mean over the 30 PSO-VNS runs at 30,030 calls: +0.08 % → −0.06 %**, i.e. below the installed layout.
     - The PyWake agreement "within 0.3 %" is about the same size as this +0.31 % probability excess.
   - Fix: replace the rounding by `_SEC = (((WD + 15.0) // 30.0).astype(int)) % 12`, or interpolate A, k and f between sectors as PyWake's `Hornsrev1Site` does. Re-run the HR16 study and regenerate Table hr-site, the NHR* macros and CHECK-FINAL [C19]/[C30]. State the PyWake settings used in the 0.3 % check (wake model, rotor averaging, direction resolution, C_T evaluated at local or free-stream speed).

2. **[major] optA/01_front.tex (abstract); also 02_intro.tex item 4 and 10_limits_concl.tex (conclusion)**
   - Quote: "is feasible in every run and exceeds the energy yield of the installed Horns Rev~1 block"
   - Problem: The claim is unqualified. Only 6 of 30 runs at 6,030 calls exceed the installed layout, the mean does not, and the best gain (0.6 %, 0.5 % after the fix in issue 1) is smaller than the modeling uncertainty. The comparison is also structurally unfair, because the optimizer may use 4D spacing while the installed grid is 7D (see issue 3).
   - Fix (abstract): "…is feasible in every run and, under the same Jensen model, its best Horns Rev~1 block layout (with $4D$ instead of the installed $7D$ spacing and without cable or load constraints) reaches an AEP within \NHR…\% of the installed layout; this difference is smaller than the modeling uncertainty and does not imply a better real design."
   - Apply the same wording in the contributions and the conclusion.

3. **[major] optA/08_beyond.tex (Horns Rev 1)**
   - Quote: "the installed layout, a regular $7D$ grid, reaches \NHRInstalled~GWh/yr"
   - Problem: The comparison needs its engineering context, and none is given. Specifically:
     - The optimized layouts may use 4D spacing (best run: 4.15D minimum), whereas the installed layout keeps 7D.
     - The 16-turbine block is the NW corner of an 80-turbine farm and is evaluated without the wakes of the other 64 turbines.
     - The parallelogram boundary is the hull of the installed turbines, not the lease area.
     - The installed regular grid was also chosen for cabling, foundations, navigation, installation and loads, none of which is modeled.
   - Fix: add a control run with $\ell_{\min}=7D$ (the installed spacing), or state explicitly that none was done. Add after the last sentence: "The optimized layouts use the benchmark spacing $4D$ (minimum \dots D in the best run), whereas the installed grid keeps $7D$; the block is evaluated without the wakes of the other 64 turbines, and array cabling, foundations, fatigue loads and navigation corridors, which also shaped the installed grid, are not modeled. The AEP differences of a few tenths of a percent are therefore a test of the optimizers, not evidence that the installed layout could be improved."

4. **[major] optA/03_model.tex (Wake Effect Model)**
   - Quote: "turbine $i$ lies in the wake if $\beta_{ij}<\alpha$. Its downstream distance is $d_{ij}=\lvert u_{ij}-R/K\rvert$"
   - Problem: Because the cone vertex lies R/K = 513 m upstream of the rotor and d_ij uses an absolute value, a turbine that lies 4D–6.7D *upstream* of j, within a few metres of the axis, is counted as waked by j. I checked this with `wflop_model.jensen_deficits`: a turbine 350 m upstream gets Δ = 0.195, the same deficit as the one 350 m downstream. The artefact is inherited from Kusiak–Song but is non-physical. The Gaussian re-evaluation (x > 0 only) does not have it, which may partly explain the lower τ there.
   - Fix: add after the deficit equation: "As in~\cite{Kusiak2010}, the cone test also admits a turbine $i$ that lies upstream of $j$ by less than $R/K$ and within $R-K\lvert u_{ij}-R/K\rvert$ of its axis; this non-physical upstream deficit is kept for comparability with the published results. The model is a top-hat wake with a hub-centre in-wake test (no partial-wake overlap), so the deficit changes discontinuously at the cone edge."

5. **[major] optA/03_model.tex / 05_setup.tex (direction discretization)**
   - Quote: "evaluated, as in~\cite{Kusiak2010}, by a Riemann sum over $N_\theta$ direction bins"
   - Problem: The top-hat cone half-angle is only arctan 0.075 = 4.3°, while the bins are 15° apart. Beyond roughly 9D, the wake of one turbine no longer covers the angular gap between two bin centres, so a downstream turbine can sit "between" directions and receive no deficit at all. In large farms, layouts optimized on coarse bins are known to exploit this; IEA37 CS1 has the same issue with 22.5° bins. Real wind is continuous, and wind-direction uncertainty smears wakes.
   - Fix: add a sentence to the Limitations and, ideally, a robustness check that re-evaluates the final layouts with 1°–5° bins (interpolated rose). Text: "With 15° (benchmark) and 22.5° (IEA37) bins and narrow top-hat or Gaussian wakes, optimized layouts can place turbines in the angular gaps between bin centres; their energy is therefore optimistic for a continuous wind rose, and a re-evaluation with finer direction bins is advisable before any layout is used."

6. **[major] optA/02_intro.tex / 10_limits_concl.tex (practical relevance)**
   - Quote: "so the layout sets the wake loss and thereby the energy yield on which capacity planning and the design of the collector and grid connection rely"
   - Problem: The introduction motivates the work with power-system planning, but no result gives net AEP, capacity factor, variability or collector design. The benchmark farms have 2–15 turbines (36 at most), well below utility offshore farms of 50–150+ turbines, where gradient-based tools with analytic derivatives (e.g., TOPFARM, FLORIS/FLOWFarm) dominate practice. Planners reading the recommendations may over-generalize.
   - Fix: replace with "The layout sets the wake loss and thereby the gross energy yield of the farm~\cite{…}; this paper addresses only this wake-loss part of layout design, for farms of up to 36 turbines." Add to Section discussion: "Our results do not show how PSO-VNS scales to farms of 50–150 turbines, where gradient-based optimization with analytic derivatives is standard practice."

7. **[major] optA/10_limits_concl.tex (Limitations)**
   - Quote: "Finally, terrain, mixed turbine types, cables, grid connection, noise, environmental constraints and economic objectives such as the levelized cost of energy are outside our scope."
   - Problem: The list omits the effects that matter most at 4D spacing:
     - wake-added turbulence and fatigue loads (IEC 61400-1 effective turbulence);
     - wake steering and yaw control, which change the optimal layout (control co-design);
     - blockage and deep-array or wake-recovery effects in large offshore farms;
     - availability and electrical losses (the reported "AEP" is gross);
     - wind-direction and wind-climate uncertainty.
   - Fix: replace the sentence with "Finally, wake-added turbulence and the resulting fatigue loads (critical at $4D$), wake steering and other control strategies, blockage and deep-array effects, terrain, mixed turbine types, array cables, grid connection, availability and electrical losses, noise, environmental constraints and economic objectives such as the levelized cost of energy are outside our scope; all reported AEP values are gross, wake-only values."

8. **[minor] optA/03_model.tex (Power Model)**
   - Quote: "We use it as published: it is slightly negative just above the cut-in speed, jumps from 1472 to 1500~kW at the rated speed and has no cut-out speed."
   - Problem: This is good, but two more limitations of the benchmark physics should be stated. First, $C_T=0.8$ is constant at all speeds, whereas real C_T falls to 0.1–0.3 above rated; this overstates the deficits above rated, and the Weibull-scale shortcut is exact only because C_T is constant. Second, Data Set I (ψ = 13 m/s, ζ = 2 in all bins) gives a capacity factor of 62 % (936/1500 kW), which is far above realistic sites. Data Set II gives 32.5 %.
   - Fix: append "The thrust coefficient is constant ($C_T=0.8$) at all speeds, which overstates the deficits above rated speed and makes the Weibull-scale form~\eqref{eq:weibull-scale} exact; Data Set~I implies a capacity factor of 62\% and Data Set~II of 33\%, so the benchmark is a synthetic test case, not a representative site."

9. **[minor] optA/03_model.tex**
   - Quote: "$\mathrm{AEP}\,[\mathrm{MWh}] = (\text{objective}/15)\times 8.76$"
   - Problem: The conversion is correct, but it gives gross AEP (100 % availability, no electrical losses) under a linear power curve that the paper itself shows overstates energy by 18–37 % (Table robust-final).
   - Fix: "…gives the gross (wake-only) AEP of the benchmark, $\mathrm{AEP}\,[\mathrm{MWh}] = (\text{objective}/15)\times 8.76$, which is not an energy estimate for a real turbine (Section~\ref{sec:robustness})."

10. **[minor] optA/05_setup.tex (Horns Rev model description)**
    - Quote: "V80 turbines (2~MW, $D=80$~m, 25-m/s cut-out), the measured 12-sector wind climate, the parallelogram through the corner turbines as boundary, $\ell_{\min}=4D$, and the Jensen wake with $K=0.04$"
    - Problem: Several modeling choices affect the result but are not stated:
      - the 12 sectors are resampled to 5° bins (see issue 1);
      - C_T is taken at the free-stream speed, not at the upstream turbine's local speed as in PyWake's NOJ;
      - the in-wake test is hub-centre with no rotor-area overlap;
      - speeds use 1 m/s bins from 3 to 25 m/s;
      - K = 0.04 is not justified. The typical offshore range is 0.04–0.05, and PyWake's NOJ default is 0.1.
    - Fix: append "The 12 sectors are resampled to 5° bins with the parameters of the enclosing sector, speeds use 1-m/s bins, the thrust coefficient is taken at the free-stream speed and the in-wake test uses the hub centre; $K=0.04$ is the usual offshore value~\cite{…}." Add a sensitivity check with K = 0.05, or a PyWake Gaussian model, for the installed-vs-optimized comparison.

11. **[minor] optA/08_beyond.tex**
    - Quote: "At 30,030 evaluations, the mean AEP of PSO-VNS (\NHRMeanThirtyKPSOVNS~GWh/yr, loss \NHRLossThirtyKPSOVNS\%) exceeds the installed layout"
    - Problem: +0.08 % (0.11 GWh/yr) is reported as exceeding the installed layout, with no uncertainty estimate. The difference is within the binning error of issue 1 (it reverses to −0.06 % after the fix) and within the 0.3 % PyWake discrepancy.
    - Fix: after re-running, report the difference with a bootstrap CI and write "is within \dots\% of the installed layout" unless the CI excludes zero.

12. **[minor] optA/08_beyond.tex, analysis/mpce_tab_iea37.tex (IEA37)**
    - Quote: "Participant~4 in both scenarios (SNOPT with wake expansion continuation \TBD{verify against Baker et al.\ (2019)})"
    - Problem: The optimizer attribution is unverified (TBD in two places). The 64-turbine scenario of Case Study 1 is silently omitted, and it is the scenario where optimizer differences were largest in Baker et al. The reader is not told that the published CS1 layouts are also known to exploit the 16-direction discretization.
    - Fix: verify the attribution against Baker et al. (2019), Table 4, and the repository `iea37-cs1-results`, then remove the TBDs. Add: "The 64-turbine scenario was not run. As the published layouts, ours are optimized for 16 discrete directions and a single wind speed; re-evaluation with finer direction bins would lower all AEP values."

13. **[minor] optA/05_setup.tex**
    - Quote: "16 directions at 9.8~m/s and a simplified Bastankhah Gaussian wake"
    - Problem: This is correct, but the other features that make CS1 an optimization test and not a realistic site should be stated: one speed equal to rated, constant C_T = 8/9, fixed k* = 0.0324555, no turbulence dependence, 2D minimum spacing.
    - Fix: "…16 directions with a single free-stream speed equal to the rated speed (9.8~m/s), constant $C_T=8/9$ and a simplified Bastankhah Gaussian wake with fixed expansion $k^*=0.0325$; with $\ell_{\min}=2D$ the case study is an optimizer test rather than a realistic site."

14. **[minor] optA/09_robust.tex / supplement Section S-gauss**
    - Quote: "the Gaussian wake model of~\cite{Bastankhah2014} with growth rate $K_{\rm G}=0.04$"
    - Problem: The Jensen and Gaussian models differ in both profile shape and expansion rate (K = 0.075 vs K_G = 0.04, not calibrated to each other), so the lower τ cannot be attributed to "a smooth deficit" alone. The Gaussian model is also not valid in the near wake below about 2–3D, which does not matter at ℓ_min = 4D but should be said.
    - Fix: replace "since layouts tuned to the sharp Jensen wake cone are judged differently by a smooth deficit" with "since the Gaussian model differs from the Jensen cone both in its smooth lateral profile and in its (uncalibrated) expansion rate, layouts tuned to the Jensen cone are judged differently". Optionally add a re-evaluation with K_G chosen to match the Jensen wake width at 7D.

15. **[minor] optA/10_limits_concl.tex (Discussion, practitioners)**
    - Quote: "PSO-VNS (preferably with $\omega=0.75$, Section~\ref{sec:split}) or a constriction PSO is a sound first choice"
    - Problem: $\omega=0.75$ was tested on only 12 cases. The recommendation also reads as general advice for industrial layout design, yet it rests on a synthetic Jensen benchmark with at most 36 turbines, no cable, load or cost terms, and 4D spacing.
    - Fix: "For benchmark-type problems of this size (up to 36 turbines) with random starts and a penalized objective, PSO-VNS (with $\omega=0.75$ in our twelve-case test) or a constriction PSO is a sound first choice. For engineering layout design, the objective should include at least a site-calibrated wake model with turbulence, the site spacing, and cable and cost terms, and the choice of optimizer should be re-checked on that objective."

16. **[minor] optA/09_robust.tex (spacing)**
    - Quote: "so they cannot be transferred to sites with a larger spacing without re-optimization"
    - Problem: Good caveat, but the engineering consequence is missing. Offshore minimum spacings are typically 5–7D in the cross-wind direction and larger in the prevailing direction, and 4D spacing raises wake-added turbulence and fatigue loads.
    - Fix: append "Offshore practice typically uses at least 5–7D, partly because of wake-added turbulence and fatigue loads at close spacing, which the benchmark does not model."

17. **[minor] optA/01_front.tex (Nomenclature) and optA/03_model.tex**
    - Quote: "$\theta$ Wind direction (blowing towards, CCW from east)."
    - Problem: The convention is consistent and correct (checked θ = 270° − θ_met against code for Horns Rev and IEA37). However, the Kusiak–Song data sets are not in meteorological convention, and Table S-wind gives only "towards" bins. Readers comparing with other codes need the dominant direction in both conventions.
    - Fix: add to the caption of Table S-wind: "The dominant bin of Data Set~I (90°–105° towards, i.e., wind from 165°–180° meteorological) …". Check and state the equivalent for Data Set II.

18. **[minor] optA/08_beyond.tex (layout figure caption)**
    - Quote: "Right: Horns Rev~1 16-turbine block, installed layout and best PSO-VNS layout at 6,030 calls."
    - Problem: The figure does not show what drives the gain, namely turbines moved off the row-aligned directions, turbines on the boundary, and the minimum spacing, nor the wind rose.
    - Fix: overlay the 12-sector wind rose and the 4D and 7D spacing circles, and extend the caption: "…; circles mark $4D$ (allowed) and $7D$ (installed) spacing; inset: Horns Rev~1 wind rose."

## Verification notes

- The Horns Rev re-evaluation in issue 1 used the stored PSO-VNS layouts in analysis/mpce_b*.csv and a corrected sector map. The rest of the model was unchanged.
- The upstream-wake check in issue 4 used `wflop_model.jensen_deficits` with θ = 0 and pairs 350 m and 600 m apart.
- The ideal values (14045.74 for Data Set I and 7315.38 for Data Set II) and ΣDS1_W = 1, ΣDS2_W = 0.9999 were reproduced.
- PyWake is not installed in the container, so the 0.3 % agreement claim could not be checked.

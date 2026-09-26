# wind-layout-yield

The agent has three years of met-mast data from a site and twenty years of hourly data from a reanalysis grid node and an airport station. From these it must build the 2004-2023 hub-height wind climate, optimise a 24-turbine layout inside a 1 km-radius circular lease under a Jensen top-hat wake model, and report the long-term net yield that layout actually delivers.

## Difficulty

The data is **synthetic**, from a disclosed, seeded generator (`authoring/provenance/generate.py`, seed 2026). The generator first writes the true hourly hub-height wind for all of 2004-2023. Separate measurement models then produce the mast, reanalysis and station records from it.

The task is a pre-construction energy-yield assessment: several dependent resource-assessment decisions feed a hard layout-optimisation problem. The numbers below are errors on the shipped instance, with the 10-seed calibration range in brackets.

1. **The long-term correction method (MCP).** The mast covers only 2021-2023, a campaign about 10% windier than the long term. Both references correlate only moderately with the site hour by hour (r ≈ 0.84 and 0.81), as real reanalysis and airport data do.

   Regression-type MCP (linear regression on one or both references, with or without resampled residuals) and conditional analog/nearest-neighbour methods suffer regression dilution. Because the predictor carries noise, the fitted relationship pulls the long-term estimate back toward the anomalous concurrent period. The result is +3.5% (+1.0 to +4.9%) too high. Leave-one-year-out checks inside a three-year campaign do not reveal this bias, because the held-out years scatter on both sides of the training mean.

   The expert choice is a variance-preserving mapping, such as a variance ratio on the reference or on a combined reference predictor. Such a mapping transfers the full long-term shift in the reference to the site.

   Skipping MCP entirely (mast record used as the climate) is +18% too high.
2. **Reference homogeneity.** The reanalysis node has an input-stream change in mid-2012 and reads about 6% low before it. The step shows clearly in the reanalysis/station ratio. Correcting against the raw reanalysis gives −2.2% (−1.9 to −2.9%).
3. **Tower shadow.** The two 60 m cups sit on booms pointing 232° and 52° (see `site.json`). Each is waked by the lattice tower in the opposite sector, and the prevailing wind is from about 232°. Averaging the cups gives −4.3% (−4.2 to −5.0%).
4. **Icing.** About 5% of mast hours are iced unheated cups: all cups flat near zero under a frozen vane. Keeping them gives −8.7% (−8 to −10%).
5. **Shear.** The exponent must be measured from the same-boom 40 m / 60 m A pair, and it varies from about 0.11 by day to 0.27 by night. Assuming 1/7 gives −2.0% (−1.8 to −3.0%).
6. **Reference veer.** The reanalysis direction is rotated about 14° from the site vane. Optimising on the uncorrected long-term rose puts the layout 0.90% (0.6 to 1.4%) below the best-known one.
7. **Optimisation and honest yield.** 24 turbines in a 1 km radius at 4D spacing is a nonconvex, discontinuous problem. A hand-made ring layout falls 2.3-2.5% short. Coarse-sector (30°) models over-predict their own yield by 13-16%.

Items 2-7 are careful data work. Item 1 needs the conceptual point that, with a noisy reference and an anomalous concurrent period, a regression predicts the conditional mean for the concurrent climate, not the long-term climate.

## Reference solution

`solution/solve.py` uses NumPy and pandas and runs in about 90 s on 2 CPUs. The choices it makes:

* **Mast QC.** Drops logger gaps. Flags an iced episode when all three cups are below 1 m/s, the vane is unchanged, the run lasts at least 3 h, and the event touches air below 1 °C.
* **Tower shadow.** At 60 m, uses the cup whose boom points into the wind.
* **Shear.** Takes a power-law exponent for each hour of the day from the 40 m / 60 m A pair and extrapolates the merged 60 m speed to 80 m.
* **Homogenisation.** Finds a single break in the monthly reanalysis/station ratio by least-squares change-point search and rescales the reanalysis before it.
* **MCP.** Estimates veer as the circular mean of the vane minus the reanalysis direction. In each 30° sector of the veer-corrected reanalysis direction, builds a least-squares combination of reanalysis and station speed as the predictor. Maps that predictor to hub height by **variance ratio** (matching concurrent mean and standard deviation), giving an hourly 80 m series for 2004-2023.
* **Layout.** Uses a 1° wind rose with empirical power tables. Runs 6 random-search starts with incremental evaluation, then a local polish and a simulated-annealing polish.
* **Yield.** Evaluates the final layout hour by hour over the long-term series.

Other valid strategies exist, and the verifier requires none of these choices. The independent implementation (`authoring/evidence/indep.py`) takes a different route at every step:

* rolling flat-line icing test
* 60 m cup chosen from the sector-wise A/B ratio
* per-record shear
* CUSUM break detection
* variance ratio on an equal-weight, mean-normalised reanalysis+station predictor in 45° sectors
* simulated annealing on empirical 1° tables

It passes all 30 calibration seeds. Quantile mapping on a combined predictor also mostly passes (9/10, +1.8% median).

## Verification

The **ground truth** is the clean hourly hub-height speed and direction for 2004-2023 (`tests/truth.npz`). The generator writes it before any measurement model, and it never comes from a solver. The verifier re-evaluates the submitted layout for all 175,320 hours with the wake model stated in `environment/data/README.md`.

`P_best` (`tests/benchmark.json`, 12,710.22 kW) is the best true mean power found on the shipped instance by any of four searches:

* a truth-informed search: 6 random-search starts plus simulated annealing, run directly on the sealed truth
* the reference solver with seed 20260925
* the reference solver with seed 11 and 12 starts
* the independent solver

Gates, all required. In the table, "reference" and "independent" give the worst case over 30 regenerated seeds (401-430).

| gate | threshold | reference worst | independent worst | nearest wrong route (median) | why here |
|---|---|---|---|---|---|
| feasibility | r ≤ 1000 m and spacing ≥ 308 m, 0.1 m tolerance | pass | pass | n/a | stated hard constraint |
| layout quality: 1 − P_true/P_best | ≤ 0.75% | 0.30% | 0.31% | no veer correction: 0.83% (0/10) | 2.4× above the correct-route spread; below veer-less, station-only, coarse-sector and ring routes |
| yield accuracy: \|P_reported/P_true − 1\| | ≤ 2.0% | +1.05% | −1.17% | regression+residual MCP: +3.08% (1/10); analog MCP: +3.11% (1/10); no homogenisation: −2.42% (2/10) | 1.7-1.9× above the correct routes' worst case; below every major resource error |

Calibration used 30 regenerated seeds: reference, independent solver and truth-informed benchmark on each, with the full ablation ladder on 10 of them (`authoring/evidence/calibration/`). The shipped instance was also run through the whole ladder (`authoring/evidence/shipped_layouts/`). The verifier never executes agent code; it parses one CSV and one JSON file.

## Ablation ladder

Medians over 10 seeds with the pass count at the gates above. The last column is the shipped instance.

| route | shortfall % | yield error % | passes (10 seeds) | shipped: shortfall / yield error |
|---|---|---|---|---|
| full reference (VR on combined predictor) | 0.09 | −0.08 | 30/30 | 0.29 / +0.39 |
| independent implementation (VR, equal-weight predictor) | 0.20 | −0.26 | 30/30 | 0.38 / +0.34 |
| quantile mapping on combined predictor | 0.19 | +1.81 | 9/10 | 0.27 / +1.74 |
| linear regression on both references + resampled residuals | 0.20 | +3.08 | 1/10 | 0.29 / +3.46 |
| nearest-neighbour analog on both references | 0.23 | +3.11 | 1/10 | 0.42 / +3.48 |
| MCP against the station only | 1.53 | +2.12 | 0/10 | 1.76 / +5.26 |
| no reference homogenisation | 0.25 | −2.42 | 2/10 | 0.28 / −2.17 |
| mast record as long-term (no MCP) | 0.20 | +13.64 | 0/10 | 0.18 / +18.12 |
| average of the two 60 m cups | 0.23 | −4.56 | 0/10 | 0.31 / −4.34 |
| 60 m cup A only | 0.26 | −2.17 | 4/10 | 0.21 / −1.85 |
| no icing QC | 0.23 | −9.01 | 0/10 | 0.27 / −8.73 |
| shear assumed 1/7 | 0.24 | −2.42 | 2/10 | 0.30 / −2.00 |
| no veer correction of reference directions | 0.83 | +0.98 | 0/10 | 0.90 / +1.20 |
| 30° sectors, yield from the model | 3.53 | +14.09 | 0/10 | 3.32 / +14.00 |
| weak optimiser (1 start, 3k moves) | 0.69 | −0.15 | 7/10 | — |
| concentric-ring layout | 2.29 | −0.20 | 0/10 | 2.47 / +0.20 |

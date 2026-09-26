# wind-layout-yield

The agent has three years of met-mast data from a site and twenty years of hourly data from a reanalysis grid node and an airport station. From these it must build the 2004-2023 hub-height wind climate, optimise a 24-turbine layout inside a 1 km-radius circular lease under a Jensen top-hat wake model, and report the long-term net yield that layout actually delivers.

## Difficulty

The data is **synthetic**, from a disclosed, seeded generator (`authoring/provenance/generate.py`, seed 2026). The generator first writes the true hourly hub-height wind for all of 2004-2023. Separate measurement models then produce the mast, reanalysis and station records from it.

The task is a pre-construction energy-yield assessment: several dependent resource-assessment decisions feed a hard layout-optimisation problem. Each decision below is individually standard practice for a wind-resource engineer, and a generalist tends to miss at least one. The numbers are errors on the shipped instance, with the 10-seed calibration range in brackets.

1. **Long-term correction (MCP).** The mast covers only 2021-2023, which were windy years with more flow from the weaker north-east lobe. Using the mast record as the long-term climate over-predicts by +10.0% (+1.6 to +12.5%).
2. **Reference homogeneity.** The reanalysis node has an input-stream change in mid-2012 and reads about 6% low before it. The station is noisier but homogeneous, so the step shows up clearly in the reanalysis/station ratio. MCP against the raw reanalysis under-predicts by −5.2% (−4.8 to −5.4%). Nothing in the data or instruction flags the step; the expert habit is to check the reference's consistency before using it.
3. **Tower shadow.** The two 60 m cups sit on booms pointing 232° and 52°. The site metadata gives the boom azimuths, and the A/B ratio against direction shows the shadow directly. Each cup is waked by the lattice tower when the wind comes from its own boom's opposite side, and the prevailing wind is from about 232°. Averaging the cups under-predicts by −5.3% (−4.8 to −5.3%). Using cup A alone gives −2.7% (−1.8 to −2.9%).
4. **Icing.** About 5% of mast hours are iced unheated cups: all cups flat near zero under a frozen vane. Keeping them biases the MCP and gives −8.6% (−7 to −11.5%).
5. **Shear.** The exponent must be measured from the same-boom 40 m / 60 m A pair, and it varies from about 0.11 by day to 0.27 by night. Assuming 1/7 gives −2.9% (−2.7 to −3.1%).
6. **Reference veer.** The reanalysis direction is rotated about 14° from the site vane. Optimising on the uncorrected long-term rose puts the layout 1.04% (0.8 to 1.4%) below the best-known one.
7. **Optimisation and honest yield.** 24 turbines in a 1 km radius at 4D spacing is a nonconvex, discontinuous problem. A weak optimiser or a hand-made ring layout falls 0.8-2.5% short. Coarse-sector (30°) models over-predict their own yield by 14-16%.

## Reference solution

`solution/solve.py` uses NumPy and pandas and runs in about 95 s on 2 CPUs. The choices it makes:

* **Mast QC.** Drops logger gaps. Flags an iced episode when all three cups are below 1 m/s, the vane is unchanged, the run lasts at least 3 h, and the event touches air below 1 °C.
* **Tower shadow.** At 60 m, uses the cup whose boom points into the wind.
* **Shear.** Takes a power-law exponent for each hour of the day from the 40 m / 60 m A pair (same boom, so the shadow cancels in the ratio) and extrapolates the merged 60 m speed to 80 m.
* **Homogenisation.** Finds a single break in the monthly reanalysis/station mean-speed ratio by least-squares change-point search, then rescales the reanalysis before the break.
* **MCP.** Estimates veer as the circular mean of the vane minus the reanalysis direction. Applies a variance-ratio regression per 30° sector of the veer-corrected reanalysis direction to produce an hourly 80 m series for 2004-2023.
* **Layout.** Uses a 1° wind rose with empirical power tables. Runs 6 random-search starts with incremental evaluation, then a local polish and a simulated-annealing polish.
* **Yield.** Evaluates the final layout hour by hour over the long-term series.

Other valid strategies exist, and the verifier requires none of these choices. The independent implementation (`authoring/evidence/indep.py`) takes a different route at every step:

* rolling flat-line icing test
* 60 m cup chosen from the sector-wise A/B ratio
* per-record shear
* CUSUM break detection
* per-sector quantile-mapping MCP
* simulated annealing on empirical 1° tables

It passes all 30 calibration seeds. OLS-per-sector MCP and the max of the two cups are also acceptable routes (OLS passes 10/10).

## Verification

The **ground truth** is the clean hourly hub-height speed and direction for 2004-2023 (`tests/truth.npz`). The generator writes it before any measurement model, and it never comes from a solver. The verifier re-evaluates the submitted layout for all 175,320 hours with the wake model stated in `environment/data/README.md`.

`P_best` (`tests/benchmark.json`, 12,515.97 kW) is the best true mean power found on the shipped instance by any of four searches:

* a truth-informed search: 6 random-search starts plus simulated annealing, run directly on the sealed truth
* the reference solver with seed 20260925
* the reference solver with seed 11 and 12 starts
* the independent solver

Gates, all required. In the table, "reference" and "independent" give the worst case over 30 regenerated seeds (301-330).

| gate | threshold | reference worst | independent worst | nearest wrong route (median) | why here |
|---|---|---|---|---|---|
| feasibility | r ≤ 1000 m and spacing ≥ 308 m, 0.1 m tolerance | pass | pass | n/a | stated hard constraint |
| layout quality: 1 − P_true/P_best | ≤ 0.75% | 0.31% | 0.45% | no veer correction: 0.97% (0/10) | ≥ 1.6× above the correct-route spread; below veer-less, station-only, coarse-sector and ring routes |
| yield accuracy: \|P_reported/P_true − 1\| | ≤ 2.5% | 0.79% | 1.68% | 1/7 shear: −2.79%; no homogenisation: −5.05%; averaged cups: −4.97% | 1.5× above the independent route's worst case, which is driven by step-factor sampling noise; below every major resource error |

Calibration used 30 regenerated seeds: reference, independent solver and truth-informed benchmark on each, with the full ablation ladder on 10 of them (`authoring/evidence/calibration/`). The shipped instance was also run through the whole ladder (`authoring/evidence/shipped_layouts/`). The verifier never executes agent code; it parses one CSV and one JSON file.

## Ablation ladder

Medians over 10 seeds with the pass count at the gates above. The last column is the shipped instance.

| route | shortfall % | yield error % | passes (10 seeds) | shipped: shortfall / yield error |
|---|---|---|---|---|
| full reference | 0.10 | −0.32 | 30/30 | 0.08 / −0.42 |
| independent implementation | 0.20 | +0.49 | 30/30 | 0.22 / +0.17 |
| OLS MCP instead of variance ratio | 0.20 | −1.29 | 10/10 (acceptable) | 0.16 / −1.20 |
| no reference homogenisation | 0.23 | −5.05 | 0/10 | 0.39 / −5.21 |
| mast record as long-term (no MCP) | 0.24 | +7.22 | 1/10 | 0.18 / +9.99 |
| average of the two 60 m cups | 0.23 | −4.97 | 0/10 | 0.18 / −5.26 |
| 60 m cup A only | 0.24 | −2.54 | 4/10 | 0.12 / −2.71 |
| no icing QC | 0.24 | −8.60 | 0/10 | 0.36 / −8.56 |
| shear assumed 1/7 | 0.23 | −2.79 | 2/10 | 0.17 / −2.88 |
| no veer correction of reference directions | 0.97 | +0.66 | 0/10 | 1.04 / +0.63 |
| MCP against the station only | 1.56 | +2.09 | 0/10 | — |
| 30° sectors, yield from the model | 4.30 | +15.02 | 0/10 | 3.51 / +13.86 |
| weak optimiser (1 start, 3k moves) | 0.79 | −0.41 | 4/10 | — |
| concentric-ring layout | 2.33 | −0.53 | 0/10 | 2.46 / −0.66 |

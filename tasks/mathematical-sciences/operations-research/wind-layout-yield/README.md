# wind-layout-yield

Starting from a raw two-level met-mast record, the agent builds a hub-height wind climate, optimises a 24-turbine layout inside a 1 km-radius circular lease under a Jensen top-hat wake model, and reports a net-yield estimate that the layout actually delivers.

## Difficulty

The data is **synthetic**, from a disclosed, seeded generator (`authoring/provenance/generate.py`, seed 2026). The generator first produces the clean hourly hub-height wind, which is the truth. It then applies a measurement model: power-law shear with an hourly-varying exponent, cup and vane noise, unheated-cup icing episodes, and logger gaps.

An expert has to get four dependent decisions right. A careful generalist gets at least one of them wrong.

1. **Hub-height resource from a two-level mast.** The cups sit at 40 and 60 m and the hub is at 80 m. The shear exponent (about 0.20) has to be measured from the cup pair. Using the 60 m speed directly under-predicts yield by about 10%. Assuming the textbook 1/7 exponent under-predicts it by 2-3.4%.
2. **Icing QC.** About 5% of hours are iced-cup episodes: both cups flat-lined near zero with a frozen vane, starting in sub-zero air, and sometimes lasting after the air has warmed. They look like calms. Keeping them biases yield by −4.4 to −5.4%. They are identifiable from the data: a frozen vane, flat cups at both heights, and the temperature column.
3. **Directional resolution of the optimisation objective.** This is the core judgment call. The textbook WFLOP recipe, as in the Kusiak-Song benchmark, bins the rose into 12-36 sectors and evaluates each at its centre line. With a top-hat wake and a fairly narrow prevailing lobe, that lets the optimiser park turbines just off each sector's centre-line wake. The model then reports near-ideal efficiency (0.91-0.97), but against the real hourly directions the layout is 1.0-1.5% (10° sectors) or 3.7-5.7% (30° sectors) worse than a properly optimised one. The 30° layouts are worse than a random layout. This is the same artefact that lets published small-N WFLOP layouts hit "ideal" power. Nothing in the data or instruction says which resolution to use. The expert has to notice that the wake cone subtends only a few degrees at farm spacings and that the grading runs on hourly directions.
4. **Honest yield of the submitted layout.** Reporting the optimiser's own objective from a coarse-sector model over-predicts by 4.6% (10° sectors) to 15-17% (30° sectors). Evaluating power at each sector's mean wind speed under-predicts by about 20%. The yield must be recomputed for the final layout at full directional resolution, over the speed distribution rather than at mean speeds.

On top of this, 24 turbines in a 1 km radius at 4D spacing is a nonconvex, discontinuous problem. A poor optimiser (the weak_optimiser rung) or a hand-made ring layout loses 0.7-2.3% even with a correct climate. Several of the steps are individually textbook, but only the combination is graded, and the directional-resolution choice has to be reasoned out, not looked up.

## Reference solution

`solution/solve.py` uses NumPy only and runs in about 60 s on 2 CPUs. The choices it makes:

* Drops logger gaps. Flags iced episodes: both cups below 1 m/s, vane unchanged, run of at least 3 h, and the event touches air below 1 °C.
* Takes the shear exponent from the ratio of mean 60 m and 40 m speeds over hours with both cups above 3 m/s, then extrapolates the 60 m cup to 80 m.
* Builds a 1° wind rose, using vane directions directly. For each bin it tabulates the empirical mean turbine power as a function of the wake speed factor (1 − deficit), so the model integrates the full speed distribution without fitting one.
* Runs 6 random-search starts of 60k moves each, then a local polish. Moves stay feasible, and each is evaluated incrementally at O(bins × N) cost.
* Reports `net_mean_power_kw` by evaluating the final layout hour by hour on the cleaned hub-height series.

Other valid strategies exist, and the verifier requires none of these choices. The independent implementation (`authoring/evidence/indep.py`) takes a different route at every step:

* A rolling flat-line icing test.
* Sector-median shear applied per record.
* A parametric climate: Weibull MLE per 10° sector, with a von Mises KDE (3° bandwidth) for direction on a 1° grid.
* Simulated annealing from a ring seed.
* A yield taken from the parametric model.

It passes on every calibration seed.

## Verification

The **ground truth** is the clean hourly hub-height speed and direction written by the generator before the measurement model is applied (`tests/truth.npz`). It never comes from a solver. The verifier re-evaluates the submitted layout hour by hour over all 87,660 hours with the wake model stated in `environment/data/README.md`.

`P_best` (`tests/benchmark.json`, 12,838.06 kW) is the best true mean power found on the shipped instance by any of: a truth-informed random search that optimises directly on the sealed truth (4 starts), the reference solver (seed 20260925, and a 12-start run with seed 11), and the independent solver. The shipped oracle comes 0.03% short of it.

Gates, all required. In the table, "reference" and "independent" give the worst case over 30 regenerated seeds (201-230).

| gate | threshold | reference worst | independent worst | nearest wrong route (median) | why here |
|---|---|---|---|---|---|
| feasibility | r ≤ 1000 m and spacing ≥ 308 m, 0.1 m tolerance | pass | pass | n/a | stated hard constraint |
| layout quality: 1 − P_true/P_best | ≤ 0.75% | 0.20% | 0.36% | 10° sectors: 1.09% (0/10 pass) | ≥ 2× above correct-route spread; below every coarse-sector route |
| yield accuracy: \|P_reported/P_true − 1\| | ≤ 2.0% | 0.94% | 1.24% | 1/7 shear: −2.62% (0/10); 10° model: +4.56% | ≥ 1.6× above correct-route error, which comes from sampling of the removed hours plus about +0.4% optimiser's-curse bias |

Calibration used 30 regenerated seeds (reference, independent solver and truth-informed benchmark on each), with the full ablation ladder on 10 of them (`authoring/evidence/calibration/`). An earlier 5-year version of the record put the reference yield error within 1.3× of the gate. That is why the record is 10 years with shorter logger outages; both summaries are kept.

The shipped instance was checked end to end in containers. The oracle scores reward 1 in 3 of 3 harbor runs and nop scores 0 in 3 of 3 (`authoring/evidence/harbor/`). Random, ring, 30°-sector, 10°-sector and no-QC outputs each fail the intended gate. The verifier never executes agent code; it parses one CSV and one JSON file.

## Ablation ladder

Medians over 10 seeds, with the pass count at the gates above. The last two columns are the values on the shipped instance where they were measured.

| route | shortfall % | yield error % | passes | shipped: shortfall / yield error |
|---|---|---|---|---|
| full reference | 0.00 | +0.38 | 30/30 | 0.03 / −0.09 |
| independent implementation | 0.08 | +0.55 | 30/30 | 0.32 / +0.09 |
| skip icing QC | 0.02 | −4.55 | 0/10 | 0.03 / −5.21 |
| skip shear (60 m as hub) | 0.00 | −10.34 | 0/10 | — |
| assume shear 1/7 | 0.06 | −2.62 | 0/10 | — |
| 30° sectors, yield from the model | 3.72 | +15.26 | 0/10 | 5.62 / +16.45 |
| 30° sectors, yield re-evaluated hourly | 3.72 | +0.18 | 0/10 | — |
| 10° sectors, yield from the model | 1.09 | +4.56 | 0/10 | — |
| 10° sectors, yield re-evaluated hourly | 1.09 | +0.13 | 0/10 | 1.53 / −0.26 |
| 5° sectors, yield re-evaluated hourly | 0.72 | +0.20 | 7/10 (borderline by design) | — |
| yield at sector-mean wind speeds | 0.02 | −20.30 | 0/10 | — |
| weak optimiser (1 start, 3k moves) | 0.69 | +0.27 | 5/10 | — |
| concentric-ring layout, honest yield | 2.17 | +0.12 | 0/10 | 2.30 / −0.40 |
| naive: random layout, "no wake loss" yield | 5.26 | +19.63 | 0/10 | 5.36 / +19.12 |

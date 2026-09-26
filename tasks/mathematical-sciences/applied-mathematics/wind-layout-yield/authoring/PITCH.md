# Pitch (template section 0)

- **Workflow.** A pre-construction energy-yield assessment. Turn a short, imperfect met-mast campaign plus long-term reference data into a long-term hub-height wind climate, optimise a turbine layout under a wake model, and decide the net yield it will deliver.
- **Who does this.** A wind-resource / layout engineer at a developer or consultancy before financial close. In research, anyone benchmarking WFLOP metaheuristics on realistic resource data.
- **The hard part.**
  - Long-term correction: the mast years are anomalous, and the reanalysis has an inhomogeneity that only shows up against the independent station.
  - Instrument handling: tower shadow on opposite-boom cups, iced cups, and diurnal shear.
  - Direction veer between the reference and the site.
  - A hard layout optimisation on top of all this.
- **What a wrong answer looks like.** MCP straight onto the raw reanalysis with averaged cups: about -5% each, and they compound. Or the mast record used as the climate: +6 to +12%. Or a layout optimised on an uncorrected, veered rose: about 1% of true power lost.
- **Failure mode.** Silent systematic errors in the resource chain that a careful coder does not see without domain checks (reference homogeneity, tower shadow, long-term representativeness), combined with a nonconvex optimisation. This differs from contact-sync-3d, which tested robust outlier rejection in SE(3) synchronisation.
- **Data.** Synthetic, from a seeded generator. The ground truth is the clean 2004-2023 hub-height wind, written before any measurement model and never taken from a solver.

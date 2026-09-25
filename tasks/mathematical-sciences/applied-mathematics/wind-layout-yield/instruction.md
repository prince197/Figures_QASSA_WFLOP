# Wind Farm Layout and Yield Estimation

## Task

Design a layout for 24 turbines within the circular onshore lease area and estimate its long-term electrical output.

The only wind data available are ten years of hourly measurements from a met mast beside the site.

## Input files

- `/app/data/README.md` — columns, units, coordinate frame, and the exact wake model used for grading.
- `/app/data/mast_timeseries.csv` — hourly mast measurements.
- `/app/data/power_curve.csv` — turbine power curve.
- `/app/data/site.json` — turbine specifications, wake-model parameters, and lease limits.

## Files to submit

1. `/app/output/layout.csv` with the header `turbine_id,x_m,y_m` and exactly 24 turbine rows. Give coordinates in metres using the frame specified in the README.
2. `/app/output/yield.json` containing a JSON object with these two numeric fields:
   - `net_mean_power_kw`: estimated long-term mean electrical output of the submitted layout after wake losses, in kW.
   - `gross_mean_power_kw`: estimated long-term mean electrical output with wakes ignored, assuming 24 free-standing turbines, in kW.

Treat the full ten years in the mast record as the long-term wind climate. Assume the turbines operate in every hour of that period, including hours when the mast was down or its readings were unreliable.

## Evaluation

The submitted layout is evaluated for every hour in the ten-year period using that hour's true hub-height wind speed and direction and the wake model specified in the README.

- All 24 turbines must be inside the lease circle. Every pair of turbines must meet the minimum spacing requirement. There is a 0.1 m tolerance for both limits. Any violation fails the task.
- The true net mean power of your layout is compared with the best layout known for this site.
- Your `net_mean_power_kw` is compared with the true net mean power of the layout you submit.

`gross_mean_power_kw` is used for diagnostics only.

## Additional instructions

- JSON numbers may be integers or floats. Extra keys in `yield.json` are ignored.

You have 9000 seconds to complete this task. Do not cheat by using online solutions or hints specific to this task.

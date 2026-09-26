# Wind Farm Layout and Yield Estimation

## Task

Design a layout for 24 turbines within the circular onshore lease area and estimate its long-term electrical output.

On-site wind data cover only three years of hourly measurements from a met mast beside the site. Twenty years of hourly data are available from a reanalysis grid node and from an airport station.

## Input files

- `/app/data/README.md` — columns, units, coordinate frame, and the exact wake model used for grading.
- `/app/data/mast_timeseries.csv` — hourly mast measurements, 2021–2023.
- `/app/data/reference_timeseries.csv` — hourly reanalysis node data, 2004–2023.
- `/app/data/station_timeseries.csv` — hourly airport station data, 2004–2023.
- `/app/data/power_curve.csv` — turbine power curve.
- `/app/data/site.json` — turbine specifications, wake-model parameters, lease limits, and instrument metadata.

## Files to submit

1. `/app/output/layout.csv` with the header `turbine_id,x_m,y_m` and exactly 24 turbine rows. Give coordinates in metres using the frame specified in the README.
2. `/app/output/yield.json` containing a JSON object with these two numeric fields:
   - `net_mean_power_kw`: estimated long-term mean electrical output of the submitted layout after wake losses, in kW.
   - `gross_mean_power_kw`: estimated long-term mean electrical output with wakes ignored, assuming 24 free-standing turbines, in kW.

The long-term wind climate is the period 2004–2023. Assume the turbines operate in every hour of that period.

## Evaluation

The submitted layout is evaluated for every hour of 2004–2023 using that hour's true hub-height wind speed and direction and the wake model specified in the README.

- All 24 turbines must be inside the lease circle. Every pair of turbines must meet the minimum spacing requirement. There is a 0.1 m tolerance for both limits. Any violation fails the task.
- The true net mean power of your layout is compared with the best layout known for this site.
- Your `net_mean_power_kw` is compared with the true net mean power of the layout you submit.

`gross_mean_power_kw` is used for diagnostics only.

## Additional instructions

- JSON numbers may be integers or floats. Extra keys in `yield.json` are ignored.

You have 9000 seconds to complete this task. Do not cheat by using online solutions or hints specific to this task.

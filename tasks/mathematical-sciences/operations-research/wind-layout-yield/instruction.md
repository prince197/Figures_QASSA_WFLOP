## Context

We are laying out a 24-turbine onshore farm inside a circular lease area. The only wind data is ten years of hourly records from a met mast next to the site. The mast has cups at 40 m and 60 m and a vane at 58 m, and the turbines have an 80 m hub. A previous consultant sent us a layout claiming almost no wake loss, and nobody here believed the number. We need our own layout and a yield estimate we can stand behind.

## Inputs

- `/root/data/mast_timeseries.csv`: the hourly mast record.
- `/root/data/power_curve.csv`: the turbine power curve.
- `/root/data/site.json`: turbine, wake-model and lease parameters.
- `/root/data/README.md`: column meanings, units, the coordinate frame, and the exact wake model used for grading. Read it first.

## Deliverable

1. `/root/output/layout.csv` with header `turbine_id,x_m,y_m` and exactly 24 rows. Coordinates are in metres in the frame given in the data README.
2. `/root/output/yield.json` containing an object with two numbers:
   - `net_mean_power_kw`: your estimate of the farm's long-term mean electrical output for this layout, after wake losses, in kW.
   - `gross_mean_power_kw`: the same with wakes ignored (24 free-standing turbines).

Take the ten years covered by the mast record as the long-term climate. The turbines run every hour of that period, including hours when the mast was down or reading badly.

## What is graded

Your layout is evaluated hour by hour over the full ten-year period, using the true hub-height wind speed and direction of each hour and the wake model in the data README.

- **Layout feasibility (hard constraint).** Every turbine must lie inside the lease circle, and every pair must be at least the minimum spacing apart. A tolerance of 0.1 m is allowed on both. Any violation fails the task.
- **Layout quality.** The true net mean power of your layout is compared with the best layout we know of for this site.
- **Yield accuracy.** `net_mean_power_kw` is compared with the true net mean power of the layout you submitted.

`gross_mean_power_kw` is reported for diagnostics only.

## Notes

- Units are kW for power and metres for distances. Directions follow the meteorological convention: the direction the wind comes from, measured clockwise from north.
- Numbers may be written as integers or floats. Extra keys in `yield.json` are ignored.

You have 7200 seconds to complete this task. Do not cheat by using online solutions or hints specific to this task.

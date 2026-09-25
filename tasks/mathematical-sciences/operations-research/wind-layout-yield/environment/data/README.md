# Site data

## mast_timeseries.csv
Hourly means from the site met mast, 2019-01-01 00:00 to 2023-12-31 (UTC, hour-beginning).

| column | unit | meaning |
|---|---|---|
| timestamp_utc | ISO 8601 | start of the averaging hour |
| ws_40m | m/s | cup anemometer at 40 m |
| ws_60m | m/s | cup anemometer at 60 m |
| wd_58m | deg | wind vane at 58 m, direction the wind blows **from**, clockwise from true north, 1-degree resolution |
| temp_2m | degC | air temperature at 2 m |
| rh_2m | % | relative humidity at 2 m |

Empty fields mean the logger recorded nothing that hour. The cups are unheated and read 0.00 below their starting threshold (about 0.3 m/s). The vane sits close enough to hub height that no veer correction is needed. The mast stands in flat, open terrain next to the planned farm, so the hub-height wind is taken to be the same at every turbine position.

## power_curve.csv
Power curve of the planned turbine at hub-height free-stream speed, valid at site air density. Interpolate linearly between table points. Output is zero above 25 m/s.

## site.json
Farm and turbine parameters: turbine count, rotor diameter, hub height, thrust coefficient, wake decay constant, boundary radius, minimum spacing, and the mast instrument heights.

## Coordinates
x points east and y points north, in metres. The origin is the centre of the circular lease area. Each turbine centre must satisfy x^2 + y^2 <= R^2, where R = `boundary_radius_m`. Every pair of turbines must be at least `min_spacing_m` apart (4 rotor diameters, centre to centre).

## Wake model (the one used for grading)
The farm uses the Jensen/Park top-hat model:

* For free-stream speed U from direction theta, take turbine j as the source and turbine i as the target. Let x be the distance from j to i measured along the wind's direction of travel, and y the distance across it.
* i is waked by j only if x > 0 and |y| <= r0 + k*x, where r0 = D/2 and k = `wake_decay_constant`. The test is on i's hub centre, and there is no partial-overlap weighting.
* The fractional velocity deficit from one source is d_ij = (1 - sqrt(1 - Ct)) * (r0 / (r0 + k*x))^2, with Ct = `thrust_coefficient` held constant at all wind speeds.
* Deficits combine as a root sum of squares. Turbine i sees U_i = U * (1 - sqrt(sum_j d_ij^2)), and its power is the power curve evaluated at U_i.
* Farm power in an hour is the sum over turbines. Availability is 100%, and there are no electrical or other losses.

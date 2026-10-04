# rev3 precision: prespecification of the full-precision reruns (item C4; precision parts of R2 / A8)

Written before any rerun was launched, after the audit of the stored coordinates (rev3_precision_audit.py).

**Question.** Can the strict feasibility labels (tolerance 1e-6 m) and the objectives of the stored run records be
verified independently, given that the stored coordinates are rounded to 3 decimals (benchmark, IEA37) or 2 decimals
(Horns Rev 1, Lillgrund)?

**These are not new optimization runs.** No new seeds, methods, budgets or comparisons. Each rerun re-executes one
stored record (same driver, same task tuple, same seed, same budget, same initialization) through the unchanged
drivers; only the final layout is additionally written with 17 significant digits (rev3_precision_rerun.py captures
it; no driver is edited). No statistic of the paper is recomputed from rerun data; reruns only decide whether a stored
label/objective is verified.

**Records.**
1. Determinism sample: one record per study x site kind x method (160 records), the cheapest class (ii)/(iii) record
   where one exists, else the cheapest seed-1 record.
2. All records of class (ii) (label contradicted by the rounded coordinates beyond the rounding bound) and class (iii)
   (undecidable within the rounding bound) of rev3_precision_audit_records.csv: 0 and 5,714 records (34.7 CPU-hours by
   the stored Seconds). Order: the studies named by the reviewer first (rev2_grad, rev2_lg16, rev2_lg16b,
   rev2_laplace, rev2_spacing), then the remaining studies from the cheapest; whatever exceeds about 6 CPU-hours
   locally (2 worker processes) goes into rev3_precision_plan.json for cloud workers.

**Outcomes per rerun (prespecified).**
- bit-identical: Objective, MinSpacing, Calls, Feasible, rounded Coordinates and Curve strings equal to the stored
  record exactly;
- reproduced to float noise: rounded Coordinates, Calls, Curve and label equal, |Objective difference| <= 1e-9 relative
  (reported separately, never merged with bit-identical);
- diverged: anything else (reported with the record key).
- strict label from the full-precision layout with the drivers' own formula (min spacing >= s_min - 1e-6 and maximum
  distance outside the circle / parallelogram <= 1e-6 m): confirmed or changed (every change listed);
- objective replayed from the full-precision layout with the record's evaluator; difference to the stored objective.

**Decision rule for "verified at 1e-6 m".** A stored label counts as independently verified if (a) the rounded
coordinates certify it with slack beyond the rounding bound (class i), or (b) a bit-identical rerun gives the same
strict label from its full-precision layout. Labels from reruns that are reproduced only to float noise are reported
as a separate category with the margin of their slacks (|slack| compared with 1e-9 m).

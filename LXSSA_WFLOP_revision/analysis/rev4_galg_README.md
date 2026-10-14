# GA on the Lillgrund 4x4 block (rev4_galg) — reproduction

1. Run from this `analysis/` directory (Python 3.11, numpy/scipy/pandas, 4 cores); no other file is modified.
2. `python3 rev4_galg.py galg 0 1 --procs=4`   -> `rev4_galg_s0of1.csv`  (GA, seeds 1..30, budget 6,030)
3. `python3 rev4_galg.py galgb 0 1 --procs=4`  -> `rev4_galgb_s0of1.csv` (GA, seeds 1..30, budget 30,030)
4. Each run is the unchanged `rev2_site.run_lg` (the Lillgrund runner of the other nine methods, `rev2_lg16*.csv`:
   same model, l_min = 3 D, box, Init = random, feasibility test, Tracker, columns) executed via `rev2_ga._with_ga`,
   i.e. with `rev2_ga.run_method_ga` / `rev2_ga.RCGA` (unchanged; as for `rev2_ga.py gahr`) in place of run_method.
5. GA: N_p = 30, binary tournament, SBX eta_c = 15 p_c = 0.9, polynomial mutation eta_m = 20 p_m = 1/n, clipping, elitism.
6. Seeding: np.random.seed(seed) then the common init_pop draw, so the initial population is seed-paired with the
   other methods (checked for seeds 1, 2: first 30 evaluations identical for all ten methods).
7. `python3 rev4_galg_summary.py` -> prints the table and writes `rev4_galg_summary.json` (all ten methods,
   PSO-VNS vs GA exact Wilcoxon signed-rank test over seeds where both are feasible).
8. Coordinates are written with 17 significant digits (record_io.encode_coordinates, as current rev2_ga.py);
   the archived rev2_lg16*.csv predate that change and store 2 decimals. Timings: see `rev4_galg_run.log`.

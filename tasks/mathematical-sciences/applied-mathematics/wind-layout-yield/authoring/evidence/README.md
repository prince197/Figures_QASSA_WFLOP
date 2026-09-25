# Evidence

| file | what it is |
|---|---|
| `indep.py` | Independent implementation: rolling flat-line icing test, sector-median shear, Weibull-per-10-deg + von Mises KDE climate, simulated annealing, parametric yield. |
| `bench.py`, `opt.py` | Truth-informed benchmark: random search run directly on the sealed truth at 1 deg. Sets a lower bound for P_best. |
| `ablate.py` | Ablation routes: the reference pipeline with one decision skipped or swapped. |
| `truth_eval.py`, `score.py` | Author-side truth evaluator (same model as `tests/test_yield.py`) and per-run scorer. |
| `calib.py`, `summarize.py` | Multi-seed calibration driver and gate summary. |
| `calibration/` | Raw per-seed results and summaries. v1 is the 5-year record, which was rejected because the reference yield error came within 1.3x of the gate. v2 is the 10-year record as shipped, seeds 201-230. |
| `shipped_layouts/` | Layouts and yields on the shipped instance (seed 2026). The P_best of 12,838.06 kW comes from `ref2` (reference solver, seed 11, 12 starts). |
| `harbor/` | `harbor run` results: oracle 5/5 at reward 1 and nop 5/5 at reward 0. The two 22-xx runs use the final `/app` paths. |
| `verifier_local_run.sh` | Runs the built verifier image against an output directory. |
| `anti_cheat.md` | Laziest passing attempts and why they fail. |

To reproduce the calibration, put `../provenance/generate.py` and `../../solution/solve.py` next to these scripts, then run:

`python calib.py 201,202,...,230 10 results.jsonl && python summarize.py results.jsonl 0.75 2.0`

The run takes about 1.5 h on 4 cores.

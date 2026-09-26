# Evidence

| file | what it is |
|---|---|
| `indep.py` | Independent implementation: rolling flat-line icing test, sector-wise A/B-ratio cup selection, per-record shear, CUSUM break detection, variance-ratio MCP on an equal-weight reanalysis+station predictor (45-deg sectors), simulated annealing on empirical 1-deg tables. |
| `resource.py` | Long-term resource pipeline with switchable right and wrong routes, used by the ablations. |
| `ablate.py` | Ablation routes: the reference pipeline with one decision skipped or swapped. |
| `bench.py` | Truth-informed benchmark: the reference optimiser run directly on the sealed truth. Contributes to P_best. |
| `truth_eval.py`, `score.py` | Author-side truth evaluator (same model as `tests/test_yield.py`) and per-run scorer. |
| `calib.py`, `summarize.py` | Multi-seed calibration driver and gate summary. |
| `calibration/` | Per-seed results (seeds 401-430, ablations on 401-410) and the summary at the shipped gates (0.75% / 2.0%). |
| `shipped_layouts/` | Layouts and yields on the shipped instance (seed 2026) for the reference, a 12-start reference run (`ref2`, which sets P_best = 12,710.22 kW), the independent solver, the truth-informed benchmark and every ablation route; `summary.json` has the scores. |
| `harbor/` | `harbor run` on the shipped bundle: oracle 2/2 at reward 1 and nop 2/2 at reward 0. |
| `anti_cheat.md` | Laziest passing attempts and why they fail. |

To reproduce, put `../provenance/generate.py` and `../../solution/solve.py` next to these scripts, then run:

`python calib.py 401,...,430 10 results.jsonl && python summarize.py results.jsonl 0.75 2.0`

The run takes about 2.5 h on 4 cores.

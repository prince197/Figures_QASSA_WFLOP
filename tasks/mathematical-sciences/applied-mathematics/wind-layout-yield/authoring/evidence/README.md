# Evidence (v3)

| file | what it is |
|---|---|
| `indep.py` | Independent implementation: rolling flat-line icing test, sector-wise A/B-ratio cup selection, per-record shear, CUSUM break detection, per-sector quantile-mapping MCP, simulated annealing on empirical 1-deg tables. |
| `resource.py` | Long-term resource pipeline with switchable right and wrong routes, used by the ablations. |
| `ablate.py` | Ablation routes: the reference pipeline with one decision skipped or swapped. |
| `bench.py` | Truth-informed benchmark: the reference optimiser run directly on the sealed truth. Contributes to P_best. |
| `truth_eval.py`, `score.py` | Author-side truth evaluator (same model as `tests/test_yield.py`) and per-run scorer. |
| `calib.py`, `summarize.py` | Multi-seed calibration driver and gate summary. |
| `calibration/` | v3 per-seed results (seeds 301-330, ablations on 301-310) and the summary at the shipped gates (0.75% / 2.5%). The independent-solver rows for seeds 301-304 and 306 are from its final (quantile-mapping) version; the first run of those seeds used an analog MCP that carried about +1.2% bias, because the concurrent period is anomalous. |
| `shipped_layouts/` | Layouts and yields on the shipped instance (seed 2026) for the reference, a 12-start reference run (`ref2`, which sets P_best = 12,515.97 kW), the independent solver, the truth-informed benchmark, and every ablation route. |
| `harbor/` | `harbor run` on v3: oracle 2/2 at reward 1 and nop 2/2 at reward 0. |
| `anti_cheat.md` | Laziest passing attempts and why they fail. |
| `v2_superseded/` | The earlier long-mast version. Frontier agents solved it 3/3, so its data and gates were replaced; kept for the record. |

To reproduce, put `../provenance/generate.py` and `../../solution/solve.py` next to these scripts, then run:

`python calib.py 301,...,330 10 results.jsonl && python summarize.py results.jsonl 0.75 2.5`

The run takes about 2.5 h on 4 cores.

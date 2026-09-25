# Analysis supporting the revised manuscript

All scripts run with Python 3.11, NumPy, SciPy and pandas, from inside this folder.
Input: `../selected_30_run_data.csv` (720 stored follow-up runs with coordinates).

| Script | Manuscript item | Output |
|---|---|---|
| `wflop_model.py` | Benchmark evaluator (Jensen/Weibull, Kusiak–Song discretization), cubic power curve, Gaussian wake | library |
| `validate_evaluator.py` | Sec. VII-A: reproduces all 720 recorded objectives (max abs. error 8.7e-11) | stdout |
| `reevaluate_layouts.py` | Sec. VII-B/C: cubic power curve, Gaussian wake, 5D/6D check of stored layouts | `layout_reevaluation.csv`, `robustness_summary.csv` |
| `followup_statistics.py` | Sec. VI-D: case-level Friedman + Holm, bootstrap CIs, A12, batch runtimes, best layouts | `friedman_case_level.csv`, `bootstrap_ci.csv`, `batch_runtime_per_run.csv`, `best_layouts_selected_cases.csv` |
| `ssa_reference.py` | Sec. VII-D: reference SSA at 4D/5D/6D and 6,030-evaluation budget (~10 min) | `ssa_reference_runs.csv` |
| `packing_capacity.py` | Sec. VII-D, Table XIII: constructible turbine count per radius/spacing | `packing_capacity.csv` |
| `make_robustness_tables.py` | Generates the LaTeX tables of Sec. VII from the CSVs | `robustness_tables.tex` |

`robustness_text.tex` is the source text of Section VII; it is already inlined into the main manuscript.

Important caveats:
- LX-SSA is **not** re-implemented, because its Laplace parameters (phi, chi) are not in the archive.
- The reference SSA scores higher than both the recorded SSA and the recorded LX-SSA at equal budget. It is used only for comparisons within its own implementation (spacing, budget).
- `ssa_reference_runs.csv` also contains a 4D greedy-retention variant (`Greedy=True`), which the manuscript mentions in one sentence.

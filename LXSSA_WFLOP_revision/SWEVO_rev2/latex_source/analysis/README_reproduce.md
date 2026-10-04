# Available and blocked reproduction

Run commands from the archive root unless stated otherwise.

## Available

`python3 experiments/audit_archive.py` reads the eight supplied CSV families (11,640 records) and regenerates `validation/archive_audit.json`, raw-data SHA-256 checksums, Lillgrund success/conditional-mean summaries, paired all-run outcomes and two LaTeX tables. It uses the supplied standalone Lillgrund model and requires only NumPy/Pandas beyond the Python standard library. It does not call an optimizer or import the missing original modules.

`python3 experiments/check_revision.py` verifies the new float serialization, exact sign-test examples, a Holm example, Wilson boundary symmetry, PSO second-moment spectral radii, the zero-coefficient proposition counterexample and all 140 CSV checksums.

`sh latex_source/analysis/build_swevo.sh` compiles the documents and refreshes main/supplement labels; it does not regenerate research results. See `START_HERE.md` for TeX dependencies and the canonical source workflow.

## Blocked

The original `fresh_*.csv` / `mpce_*.csv` run records, summaries, original models, optimizers and most analysis/check scripts are absent. Fifteen directly identified project modules are missing: authors_objective, authors_optimizers, extra_baselines, feasible_init, hornsrev_model, hybrid_lxssa_bvns, iea37_experiments, iea37_model, init_hook, mpce_experiments, mpce_inference_extra, mpce_results, original_vns, rs_vns and wflop_model. This list is not an assertion that supplying only these fifteen files will restore every historical pipeline.

`rev2_analysis.py` depends on missing original analysis modules and records. The optimization drivers and launch scripts likewise require missing models and optimizers. Their former launch commands are retained as historical source, not presented as working full-reproduction instructions. Do not substitute newly written optimizers and label their outputs as reproductions of the original study.

## Raw results and precision

All supplied CSV contents are unchanged. Legacy coordinate strings are rounded to 0.001 m or 0.01 m; their original objective values and feasibility labels refer to predecessor positions whose full precision is not preserved. The audit reports the incompatibility with a 10^-6 m feasibility tolerance. Updated writers call `record_io.encode_coordinates` before serialization by delegated wrappers, preserving 17 significant digits prospectively. Search arithmetic, random draws and budgets were not intentionally changed. Optimizer behaviour could not be rerun because its original dependencies are absent.

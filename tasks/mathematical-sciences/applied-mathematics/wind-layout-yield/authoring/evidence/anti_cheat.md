# Anti-cheat notes (never executed by the pipeline)

What an agent can see: the noisy, iced, tower-shadowed 2021-2023 mast record; the inhomogeneous, veered reanalysis; the station; the power curve; and site.json.
What stays sealed in tests/: the clean hourly 2004-2023 hub-height speed and direction (truth.npz), plus the benchmark power.

| Attempt | Why it fails |
|---|---|
| Nop (no files) | `test_outputs_parse` cannot open /app/output/layout.csv, so the reward is 0. |
| Hard-code a literature or ring/grid layout | This lease, turbine count, spacing and wind rose appear nowhere else. A ring layout falls 2.3-2.5% short, and the gate is 0.75%. |
| Report the mast-period yield | +14 to +22% error (mast_only rung). |
| Report a regression-MCP yield | About +3% error (ols2_resid_mcp rung); without homogenisation about -2.4%. |
| Game the feasibility tolerance (turbines at r = R + 0.1 m, spacing 307.9 m) | Allowed by design. The benefit is about 1e-5 of farm power. |
| Write extra keys, strings like "12345", NaN or Infinity | Keys are normalised and numbers parsed. Non-finite values raise, and the test fails. |
| Look for truth in the environment image | There is none. environment/ holds only post-measurement-model products. The truth and benchmark live in the verifier image. |
| Recover the truth by inverting the generator | The generator is not shipped, and the noise draws cannot be inverted without the seed. |
| Write into /logs/verifier from the agent container | The verifier runs in a separate container and reads only the two declared artifacts. |

The verifier never executes agent code. It parses a CSV and a JSON file only.

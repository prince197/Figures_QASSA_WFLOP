# Anti-cheat notes (never executed by the pipeline)

What an agent can see: the noisy, iced, gappy mast record, the power curve and site.json.
What stays sealed in tests/: the clean hourly hub-height speed and true direction (truth.npz), plus the benchmark power.

| Attempt | Why it fails |
|---|---|
| Nop (no files) | `test_outputs_parse` cannot open /app/output/layout.csv, so the reward is 0. |
| Hard-code a layout from memory or the literature (Kusiak-Song or a grid) | This lease, turbine count, spacing and wind rose appear nowhere else. A generic ring or grid layout falls 2.2-2.3% short (the ring_layout rung), and the gate is 0.75%. |
| Random feasible layout, report 24 x free-stream power ("no wake loss") | About 5.3% shortfall and a yield error of about +19.6%, so both gates fail. |
| Match yield.json to its own layout without optimising | The yield gate passes, but the layout gate fails because the shortfall is 1.1-5.7% (sector routes) or 2.2-5.4% (ring or random). |
| Write the model's own coarse-sector prediction | The yield error is +4.6% to +17% (sector10_model / sector30_model rungs). |
| Game the feasibility tolerance (turbines at r = R + 0.1 m, spacing 307.9 m) | Allowed by design. The benefit is tiny, since 0.1 m changes power by about 1e-5. |
| Write extra keys, strings like "12345", NaN or Infinity | Keys are normalised and numbers parsed. Non-finite values raise, and the test fails. |
| Look for truth in the environment image | There is none. environment/ holds only generator outputs after the measurement model. The truth and benchmark live in the verifier image, which the agent never sees. |
| Recover the truth by inverting the generator | The generator is not shipped. Even with it, the hourly cup and vane noise draws cannot be inverted without the seed. |
| Write into /logs/verifier from the agent container | The verifier runs in a separate container and reads only the two declared artifacts. |

The verifier never executes agent code. It parses a CSV and a JSON file only.

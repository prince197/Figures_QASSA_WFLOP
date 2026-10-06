# Published values for the Kusiak–Song circular-farm WFLOP benchmark

Status: 2026-09-28. Companion file: `literature_values.csv`. Our values come from `fresh_bgrid.csv`, Algorithm `SSABV` (SSA-VNS), 30 seeds, 6,030 objective calls per run.

## 1. What could be verified, and what could not

**Could not access any full text.** Every publisher and repository host I tried was blocked by this environment's egress proxy: sciencedirect.com, researchgate.net, semanticscholar.org (web and API), osti.gov, arxiv.org, mdpi.com, link.springer.com, onlinelibrary.wiley.com, tandfonline.com, academia.edu, scitepress.org, core.ac.uk, scholar.archive.org and api.unpaywall.org. Web-search snippets gave only bibliographic facts and abstract-level statements. **No published objective value in this folder has been checked against the original paper.**

Verified from search records (bibliographic and abstract level only):

| Source | Verified facts |
|---|---|
| Kusiak & Song 2010, Renew. Energy 35(3):685–694, doi:10.1016/j.renene.2009.08.019 | Wind-distribution-based placement model with wake loss that depends on turbine location and wind direction. The problem is recast as a bi-criteria problem (maximize expected energy, minimize constraint violation) and solved with a multi-objective evolution strategy on a circular farm. |
| Eroglu & Seckiner 2012, Renew. Energy 44:53–62, doi:10.1016/j.renene.2011.12.013 | ACO applied to the Kusiak–Song energy-output model. The abstract says it beats the Kusiak–Song evolution strategy on expected energy and wake loss. |
| Eroglu & Seckiner 2013, Renew. Energy 58:95–107, doi:10.1016/j.renene.2013.02.019 | Particle filtering on the same model. Reported to beat both EA and ACO. |
| Bansal & Farswan 2017, Renew. Energy 107:386–402, doi:10.1016/j.renene.2017.01.064 | BBO on circular farms of radius 500, 750 and 1000 m with two wind data sets. It finds (a) layouts for a given N and (b) the largest number of turbines for a given farm. |

Other papers I located. None gave numbers I could read, so they appear in the CSV as "NA":
- Rehman, Ali & Khan, *Wind Farm Layout Design Using Cuckoo Search Algorithms*, Applied Artificial Intelligence 30(10), doi:10.1080/08839514.2017.1279043. It cites Kusiak–Song as the benchmark, but I could not check which cases it reports.
- Eroglu & Seckiner, *Wind farm layout optimization using ant colony and particle filtering approaches* (ResearchGate 318361195).
- Wang, Yuan, Cholette et al. 2018, J. Wind Eng. Ind. Aerodyn. 180:148–155. It compares discretization with Monte Carlo evaluation in the Kusiak–Song circular-farm setting.
- Kumar et al. 2024, J. Electr. Comput. Eng., doi:10.1155/2024/9406519 (a review).
- Not the same benchmark, so excluded: Bansal, Farswan & Nagar 2018 (EAAI 71:45–59, 2000 m × 2000 m square, non-uniform turbines); Pouraltafi-Kheljan et al. 2018 (JCP, 3850 m square, economic objective); Aggarwal et al. (IETE J. Res. 69(5), large grid farms); binary GSA papers (Mosetti grid).

## 2. Where the numbers in `literature_values.csv` come from

The EA/ACO/PF/BBO numbers are **transcriptions of unknown provenance** taken from this repository:
- **Primary:** the authors' archived manuscript, git commit `f5bb1a8` (manuscript source of that commit), tables `tab:500mI`, `tab:500mII`, `tab:750mI`, `tab:750mII`, `tab:1000mI` and `tab:1000mII`. They are labelled "Objective / Wake Loss (Mean)" with no SD, and the archived text calls them "literature benchmarks". The table does not cite the paper or table each number came from.
- **Secondary:** `/home/user/Figures_QASSA_WFLOP/index.html`, uploaded 2026-04-09. It agrees with the primary source except for BBO at 1000 m, Data Set II, N = 7–15, where the values differ. Both versions are in the CSV, and the conflict is flagged.

Assumed attribution: EA = Kusiak & Song 2010, ACO = Eroglu & Seckiner 2012, PF = Eroglu & Seckiner 2013, BBO = Bansal & Farswan 2017. EA/ACO/PF are given only for r = 500 m (EA: N = 2–6; ACO/PF: N = 2–8). BBO is given for 500 m (N = 2–8), 750 m (N = 2–12) and 1000 m (N = 2–15). For N = 9–10 at 500 m and N ≥ 7 for EA, the archive says "Infeasible", which most likely means the case was not reported. It is possible that the EA/ACO/PF values were copied from Bansal & Farswan's comparison table rather than from the original papers, but this is not confirmed.

Defects in the transcription:
1. BBO, Data Set II, 500 m, N = 3: 42137.21 is the Data Set I value and exceeds the ideal 21947.06. This is a copy error that the archived manuscript already flagged.
2. BBO, Data Set I, 500 m, N = 7: value 98280.5 with wake loss 239.71 does not add up to the ideal 98320.16; it is off by 200. Either the value should be 98080.45 or the wake loss should be 39.66.
3. BBO, Data Set I, 750 m, N = 10: off by 0.37, which is minor.
4. BBO, Data Set II, 1000 m, N ≥ 7: `index.html` and `f5bb1a8` disagree.
5. It is not known whether each value is the best of several runs, the mean, or a single run. The runs, budgets and feasibility handling of the source papers are also unverified.

## 3. Scaling

- **Data Set I:** same scale as our objective. The archived ideal for N = 2 is 28091.47, the same as ours (14045.74 per turbine), so these values include the 15° bin-width factor and need no conversion. Expected farm power in kW = objective / 15.
- **Data Set II:** the archived ideal is 14631.37 for N = 2, i.e. 7315.69 per turbine; ours is 7315.38 per turbine. The difference is +0.31 per turbine (+0.0042 %), so the source's Data Set II table probably differs slightly from ours, e.g. by rounding of the direction probabilities (ours sum to 0.9999). Compare Data Set II through the wake loss (ideal − value) or the efficiency (value / ideal), not the raw values. The offset is at most about 5 objective units at N = 15.

## 4. Provisional comparison (only if the transcriptions are verified)

These cases share the objective scale. Literature values are the unverified transcriptions. "Ours" is SSA-VNS over 30 seeds: mean of all runs, best feasible run, and the wake loss of that best run. The last column is the BBO wake loss as transcribed. At 500 m, N = 10 (not in this table), 26 of 30 runs are feasible for Data Set I and 25 of 30 for Data Set II.

| DS | r (m) | N | EA (Kusiak–Song) | ACO (Eroglu 2012) | PF (Eroglu 2013) | BBO (Bansal–Farswan) | SSA-VNS mean | SSA-VNS best | SSA-VNS best wake loss | BBO wake loss |
|---|---|---|---|---|---|---|---|---|---|---|
| I | 500 | 2 | 28083.42 | 28091.47 | 28091.47 | 28091.47 | 28091.47 | 28091.47 | 0.0 | 0.0 |
| I | 500 | 3 | 42101.06 | 42130.87 | 42128.32 | 42137.21 | 42128.48 | 42137.21 | -0.0 | 0.0 |
| I | 500 | 4 | 56057.77 | 56128.00 | 56135.28 | 56171.10 | 56101.80 | 56157.44 | 25.5 | 11.9 |
| I | 500 | 5 | 69922.97 | 70086.29 | 70085.66 | 70207.20 | 69967.18 | 70092.84 | 135.8 | 21.5 |
| I | 500 | 6 | 83758.79 | 84006.37 | 84007.84 | 84167.40 | 83780.79 | 83966.11 | 308.3 | 107.0 |
| I | 500 | 7 | - | 97822.66 | 97869.73 | 98280.50 | 96529.89 | 97559.43 | 760.7 | 239.7 |
| I | 500 | 8 | - | 111414.82 | 111498.83 | 111894.00 | 109316.36 | 111107.59 | 1258.3 | 472.0 |
| I | 750 | 2 | - | - | - | 28091.47 | 28091.47 | 28091.47 | 0.0 | 0.0 |
| I | 750 | 3 | - | - | - | 42137.21 | 42137.02 | 42137.21 | -0.0 | 0.0 |
| I | 750 | 4 | - | - | - | 56182.95 | 56162.33 | 56182.95 | -0.0 | 0.0 |
| I | 750 | 5 | - | - | - | 70228.69 | 70150.93 | 70228.69 | 0.0 | 0.0 |
| I | 750 | 6 | - | - | - | 84261.58 | 84095.59 | 84165.79 | 108.6 | 12.8 |
| I | 750 | 7 | - | - | - | 98298.66 | 97910.84 | 98090.86 | 229.3 | 21.5 |
| I | 750 | 8 | - | - | - | 112319.27 | 111683.88 | 112014.55 | 351.3 | 46.6 |
| I | 750 | 9 | - | - | - | 126313.25 | 125207.20 | 125862.75 | 548.9 | 98.4 |
| I | 750 | 10 | - | - | - | 140303.54 | 138304.27 | 139414.83 | 1042.5 | 153.5 |
| I | 750 | 11 | - | - | - | 154111.64 | 151286.31 | 153427.83 | 1075.3 | 391.5 |
| I | 750 | 12 | - | - | - | 167967.46 | 164279.23 | 166020.80 | 2528.1 | 581.4 |
| I | 1000 | 2 | - | - | - | 28091.47 | 28091.47 | 28091.47 | 0.0 | 0.0 |
| I | 1000 | 3 | - | - | - | 42137.21 | 42137.21 | 42137.21 | -0.0 | 0.0 |
| I | 1000 | 4 | - | - | - | 56182.95 | 56179.95 | 56182.95 | -0.0 | 0.0 |
| I | 1000 | 5 | - | - | - | 70228.69 | 70202.06 | 70228.69 | 0.0 | 0.0 |
| I | 1000 | 6 | - | - | - | 84274.42 | 84175.90 | 84230.26 | 44.2 | 0.0 |
| I | 1000 | 7 | - | - | - | 98320.16 | 98128.36 | 98230.62 | 89.5 | 0.0 |
| I | 1000 | 8 | - | - | - | 112338.35 | 112008.26 | 112166.86 | 199.0 | 27.6 |
| I | 1000 | 9 | - | - | - | 126333.03 | 125841.98 | 126106.91 | 304.7 | 78.6 |
| I | 1000 | 10 | - | - | - | 140354.84 | 139561.66 | 139990.76 | 466.6 | 102.5 |
| I | 1000 | 11 | - | - | - | 154379.54 | 153184.64 | 153727.01 | 776.1 | 123.6 |
| I | 1000 | 12 | - | - | - | 168355.41 | 166530.11 | 167512.75 | 1036.1 | 193.4 |
| I | 1000 | 13 | - | - | - | 182304.01 | 179973.54 | 181264.71 | 1329.9 | 290.6 |
| I | 1000 | 14 | - | - | - | 196264.54 | 193260.78 | 194722.33 | 1918.0 | 375.8 |
| I | 1000 | 15 | - | - | - | 210237.57 | 207119.29 | 208599.65 | 2086.4 | 448.5 |
| II | 500 | 2 | 14631.21 | 14631.37 | 14631.37 | 14631.37 | 14630.76 | 14630.76 | 0.0 | 0.0 |
| II | 500 | 3 | 21925.16 | 21928.07 | 21915.78 | 42137.21 | 21907.09 | 21946.14 | -0.0 | 0.0 |
| II | 500 | 4 | 29113.71 | 29174.20 | 29182.38 | 29232.54 | 29067.53 | 29193.87 | 67.6 | 30.2 |
| II | 500 | 5 | 36316.23 | 36256.10 | 36284.07 | 36507.35 | 35958.00 | 36206.26 | 370.6 | 71.1 |
| II | 500 | 6 | 43195.84 | 43125.19 | 43181.70 | 43795.59 | 42638.75 | 42899.29 | 993.0 | 98.5 |
| II | 500 | 7 | - | 49763.76 | 49819.71 | 50993.02 | 49008.59 | 49409.45 | 1798.2 | 216.8 |
| II | 500 | 8 | - | 56316.15 | 56498.03 | 58128.58 | 55274.87 | 55746.88 | 2776.1 | 396.9 |
| II | 750 | 2 | - | - | - | 14631.37 | 14630.76 | 14630.76 | 0.0 | 0.0 |
| II | 750 | 3 | - | - | - | 21947.06 | 21946.14 | 21946.14 | -0.0 | 0.0 |
| II | 750 | 4 | - | - | - | 29262.75 | 29210.47 | 29261.51 | 0.0 | 0.0 |
| II | 750 | 5 | - | - | - | 36578.44 | 36417.06 | 36526.37 | 50.5 | 0.0 |
| II | 750 | 6 | - | - | - | 43852.33 | 43500.11 | 43752.56 | 139.7 | 41.8 |
| II | 750 | 7 | - | - | - | 51154.82 | 50370.05 | 50674.29 | 533.4 | 55.0 |
| II | 750 | 8 | - | - | - | 58428.48 | 57033.92 | 57451.27 | 1071.8 | 97.0 |
| II | 750 | 9 | - | - | - | 65742.68 | 63632.10 | 64191.03 | 1647.4 | 98.5 |
| II | 750 | 10 | - | - | - | 72942.28 | 70000.50 | 70628.00 | 2525.8 | 214.6 |
| II | 750 | 11 | - | - | - | 80090.17 | 76264.96 | 77477.27 | 2991.9 | 382.4 |
| II | 750 | 12 | - | - | - | 87276.02 | 82500.77 | 83703.67 | 4080.9 | 512.2 |
| II | 1000 | 2 | - | - | - | 14631.37 | 14630.76 | 14630.76 | 0.0 | 0.0 |
| II | 1000 | 3 | - | - | - | 21947.06 | 21946.14 | 21946.14 | -0.0 | 0.0 |
| II | 1000 | 4 | - | - | - | 29262.75 | 29247.61 | 29261.51 | 0.0 | 0.0 |
| II | 1000 | 5 | - | - | - | 36578.44 | 36510.64 | 36557.18 | 19.7 | 0.0 |
| II | 1000 | 6 | - | - | - | 43894.12 | 43736.76 | 43846.17 | 46.1 | 0.0 |
| II | 1000 | 7 | - | - | - | 51182.14 | 50800.03 | 51067.13 | 140.5 | 27.7 |
| II | 1000 | 8 | - | - | - | 58475.62 | 57883.41 | 58213.88 | 309.1 | 49.9 |
| II | 1000 | 9 | - | - | - | 65762.46 | 64730.42 | 65033.59 | 804.8 | 78.8 |
| II | 1000 | 10 | - | - | - | 73046.09 | 71541.85 | 72076.66 | 1077.1 | 110.8 |
| II | 1000 | 11 | - | - | - | 80337.34 | 78133.13 | 78973.27 | 1495.9 | 135.2 |
| II | 1000 | 12 | - | - | - | 87546.56 | 84647.14 | 85677.87 | 2106.7 | 241.7 |
| II | 1000 | 13 | - | - | - | 94785.18 | 91179.39 | 92418.49 | 2681.4 | 318.8 |
| II | 1000 | 14 | - | - | - | 102023.84 | 97590.40 | 98469.48 | 3945.8 | 395.8 |
| II | 1000 | 15 | - | - | - | 109266.60 | 104035.26 | 104866.26 | 4864.4 | 468.7 |

Notes: in the row DS II / 500 m / N = 3, the BBO value 42137.21 is the copy error described in §2. For Data Set II, compare the wake-loss columns rather than the values.

## 5. What this means for the manuscript (important)

- **Data Set I:** SSA-VNS beats the transcribed EA (Kusiak–Song) and matches or slightly beats ACO/PF only on the best run at N ≤ 5. At N = 6–8 its best run is below ACO and PF. It is **below the transcribed BBO value in every case where BBO reports a non-zero wake loss**, often by a large margin: at 1000 m, N = 15, the best-run wake loss is 2086 versus 448 for BBO.
- **Data Set II, 500 m:** the SSA-VNS best run is even below Kusiak & Song's EA at N = 5–6, and below ACO and PF at N = 5–8. No method in our whole grid reaches these values either (best over all nine methods × 30 seeds, 500 m, N = 6: 43013.5 against EA's 43195.8). This points to a **model mismatch for Data Set II**, not just an algorithm gap. Candidates: differences in the Data Set II table (see the ideal offset above), the direction-bin convention (our bins are evaluated at 7.5°, 22.5°, …), or the penalty/feasibility handling (the Kusiak–Song ES traded power against constraint violation).
- The transcribed BBO values are also above the best value found by any of our nine methods in every non-trivial case. Our budget is small (6,030 calls), and the BBO budget is unknown.
- **Recommendation:** do not put these numbers in the paper until someone with library access has checked each one against the PDFs, recording the table and page, best versus mean, number of runs, budget, spacing and feasibility rule. If they are confirmed, present them as "reported values (different budgets, not re-run)", and compare by wake loss or efficiency. If they cannot be confirmed, cite these papers qualitatively only.
- The SSA/LX-SSA/QA-SSA values in the archived tables are the authors' own earlier runs, not published results: QA-SSA (Solanki & Deep 2024, OPSEARCH) is a function-optimization paper. They are far above the fresh re-runs (e.g. Data Set I, 500 m, N = 10: archived SSA 140416.92, best feasible fresh run of any method 138171.16), so they should not be used either.

import numpy as np, pandas as pd
from wflop_model import farm_objective, parse_coords
df = pd.read_csv("../selected_30_run_data.csv")
err, errw = [], []
for _, r in df.iterrows():
    obj, ideal = farm_objective(parse_coords(r.Coordinates), int(r.Dataset))
    err.append(obj - r.EnergyProduction); errw.append((ideal - obj) - r.WakeLoss)
err, errw = np.abs(err), np.abs(errw)
print("records", len(df), "max|dObj|", err.max(), "max|dWake|", errw.max(), "median", np.median(err))
for n, ds in [(2, 1), (3, 1), (4, 1), (2, 2), (3, 2), (4, 2)]:
    print(ds, n, farm_objective(np.array([[0, 0]] * 1 + [[1e5 * (i + 1), 3e5 * (i + 1)] for i in range(n - 1)]), ds)[1])

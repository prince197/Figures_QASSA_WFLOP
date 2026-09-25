import json
import os
import subprocess
import sys
from multiprocessing import Pool

PY = "/root/venv/bin/python"
ROUTES = ["no_icing_qc", "no_shear", "shear_one_seventh", "sector30_model", "sector30_series",
          "sector10_model", "sector10_series", "sector5_series", "sector_mean_speed",
          "weak_optimiser", "ring_layout", "random_layout"]


def sh(cmd, env=None):
    e = dict(os.environ, **(env or {}))
    r = subprocess.run(cmd, shell=True, capture_output=True, text=True, env=e)
    if r.returncode:
        raise RuntimeError(cmd + "\n" + r.stderr[-2000:])
    return r.stdout.strip().splitlines()[-1] if r.stdout.strip() else ""


def score(d, o):
    return json.loads(sh(f"{PY} score.py {d} {o}"))


def one(args):
    seed, do_ablate = args
    d = f"cal/s{seed}"
    if not os.path.exists(f"{d}/truth.npz"):
        sh(f"{PY} generate.py --seed {seed} --env-dir {d}/env --truth {d}/truth.npz")
    res = dict(seed=seed)
    if not os.path.exists(f"{d}/bench.json"):
        open(f"{d}/bench.json", "w").write(sh(f"{PY} bench.py {d} 4"))
    res["bench"] = json.load(open(f"{d}/bench.json"))["bench"]
    env = dict(WLY_DATA=f"{d}/env")
    sh(f"{PY} solve.py", dict(env, WLY_OUT=f"{d}/ref"))
    res["ref"] = score(d, f"{d}/ref")
    sh(f"{PY} indep.py", dict(env, WLY_OUT=f"{d}/ind"))
    res["ind"] = score(d, f"{d}/ind")
    if do_ablate:
        for r in ROUTES:
            sh(f"{PY} ablate.py {r} {d}/ab_{r}", env)
            res[r] = score(d, f"{d}/ab_{r}")
    return res


if __name__ == "__main__":
    seeds = [int(s) for s in sys.argv[1].split(",")]
    n_ab = int(sys.argv[2])
    jobs = [(s, i < n_ab) for i, s in enumerate(seeds)]
    with Pool(4) as p, open(sys.argv[3], "a") as fh:
        for r in p.imap_unordered(one, jobs):
            fh.write(json.dumps(r) + "\n")
            fh.flush()
            print(r["seed"], "done", flush=True)

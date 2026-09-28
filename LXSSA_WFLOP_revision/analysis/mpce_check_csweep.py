"""Check the statements drawn from the PSO coefficient sweep / bound-handling experiment (mpce_csweep.py) against
mpce_summary_csweep.json.

Usage (from analysis/):  python3 mpce_check_csweep.py [--summary mpce_summary_csweep.json]

Every sentence of the paper or supplement that rests on this experiment should carry a comment
% CHECK-CSWEEP [Snn]: <condition>. This script evaluates the conditions below and prints PASS / FAIL / PENDING
(data missing). A FAIL means the sentence next to the comment must be rewritten. Exit code 1 if any FAIL.
"""
import os, sys, json, argparse

HERE = os.path.dirname(os.path.abspath(__file__))
SWEEP = ["PSOW07C12", "PSOW07C14", "PSOW07C16", "PSOW07C17", "PSOW07C18", "PSOW07C19", "PSOW07C20"]


def g(d, *path):
    for p in path:
        if d is None:
            raise KeyError("/".join(map(str, path)))
        d = d[str(p)] if isinstance(d, dict) and str(p) in d else d[p]
    if d is None:
        raise KeyError("/".join(map(str, path)))
    return d


def st(s, a, k):
    return g(s, "settings", a, k)


CONDITIONS = []  # filled below (see CONDITIONS += ...)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--summary", default=os.path.join(HERE, "mpce_summary_csweep.json"))
    a = ap.parse_args()
    if not os.path.exists(a.summary):
        print(f"PENDING: {a.summary} missing (run mpce_csweep.py)"); sys.exit(0)
    s = json.load(open(a.summary))
    n = dict(PASS=0, FAIL=0, PENDING=0)
    for cid, text, cond in CONDITIONS:
        try:
            r = "PASS" if cond(s) else "FAIL"
        except (KeyError, TypeError, IndexError, ValueError):
            r = "PENDING"
        n[r] += 1
        print(f"{cid} {r:7s} {text}")
    print(f"\n{n['PASS']} PASS / {n['FAIL']} FAIL / {n['PENDING']} PENDING")
    sys.exit(1 if n["FAIL"] else 0)


if __name__ == "__main__":
    main()

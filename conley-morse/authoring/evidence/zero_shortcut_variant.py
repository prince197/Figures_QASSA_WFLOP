#!/usr/bin/env python3
"""Ablation: the reference solver with the invalid 'zero relation is absorbing'
shortcut re-inserted (the defect of the previous oracle generation).

Composing the zero relation {(0,0)} with an edge relation S gives {(0,w) : (0,w) in S},
which is nonzero whenever a later backward map has a kernel.  This variant stops
composing as soon as the running relation is zero, exactly as the retired oracle
generator did.  It must fail the regenerated sealed oracles.
"""
import importlib.util, json, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location('ref', ROOT / 'solution' / 'trajectory_solver.py')
ref = importlib.util.module_from_spec(spec); spec.loader.exec_module(ref)

_compose = ref._compose_linear_relations


def _absorbing_compose(R, du, dv, S, dw):
    if not R:
        return []
    return _compose(R, du, dv, S, dw)


def _between(dims, maps, left, right):
    R = [(1 << i) | ((1 << i) << dims[left]) for i in range(dims[left])]; rd = dims[left]
    for k in range(left, right):
        E = ref._edge_linear_relation(dims[k], dims[k + 1], maps[k][0], maps[k][1])
        R = _compose(R, dims[left], rd, E, dims[k + 1]); rd = dims[k + 1]
        if not R:
            break
    return R


_words = ref.holonomy_word_signature


def _holonomy(*a, **k):
    ref._compose_linear_relations = _absorbing_compose
    try:
        return _words(*a, **k)
    finally:
        ref._compose_linear_relations = _compose


ref._relation_between_nodes = _between
ref.holonomy_word_signature = _holonomy

if __name__ == '__main__':
    data = json.load(open(sys.argv[1]))
    with open(sys.argv[2], 'w') as f:
        json.dump(ref.solve_case(data), f, separators=(',', ':'))

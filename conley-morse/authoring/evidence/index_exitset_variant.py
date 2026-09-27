#!/usr/bin/env python3
"""Ablation: Conley index maps computed with the naive exit set E = P minus the
Morse boxes of the interval (instead of the boxes of P that can no longer reach
the interval).  The reference is otherwise unchanged.  For an interval with
connecting orbits the naive pair is not forward invariant, so the selector image
of E-chains is not confined to E; such probes are reported with every divisor
list computed after discarding the offending cells, which is what a solver that
ignores the exit-set definition would output."""
import importlib.util, json, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location('ref', ROOT / 'solution' / 'trajectory_solver.py')
ref = importlib.util.module_from_spec(spec); spec.loader.exec_module(ref)

_pair = ref.interval_index_pair


def _naive_pair(frame, ana, ids):
    P, _ = _pair(frame, ana, ids)
    B = set()
    for m in ids:
        B.update(ana['recurrent'][m])
    return P, P - B


def _record(frame, ana, ids, nx):
    P, E = _naive_pair(frame, ana, ids)
    H = ref.relative_homology_data_boxes(P, E, nx)
    try:
        mats = ref.index_map_matrices(P, E, ref.box_rectangles(frame, nx), nx, H)
    except AssertionError:
        mats = [[0] * len(H[d]['H']) for d in range(3)]
    return [{'dimension': d, 'space_dimension_F2': len(mats[d]), 'rank_F2': ref.map_rank(mats[d]),
             'elementary_divisors_F2': ref.kr_elementary_divisors(mats[d], len(mats[d]))} for d in range(3)]


ref.conley_index_record = _record

if __name__ == '__main__':
    data = json.load(open(sys.argv[1]))
    with open(sys.argv[2], 'w') as f:
        json.dump(ref.solve_case(data), f, separators=(',', ':'))

from pathlib import Path
from verifier_core import evaluate, write_metrics

def test_submission():
    ok,msg,diag=evaluate(Path('/app/trajectory_solver.py'))
    write_metrics(diag)
    assert ok,msg

from __future__ import annotations

import json
import math
import os
import shutil
import signal
import stat
import subprocess
import tempfile
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
HIDDEN = sorted((HERE / "heldout").glob("heldout_*.json"))
ORACLES = {p.name: json.loads((HERE / "oracles" / p.name).read_text()) for p in HIDDEN}
ARTIFACT = Path("/app/trajectory_solver.py")
CASE_TIMEOUT_S = 300
REQUIRED_CASES = 4


class DuplicateKey(ValueError):
    pass


def no_dupes(pairs):
    out = {}
    for k, v in pairs:
        if k in out:
            raise DuplicateKey(k)
        out[k] = v
    return out


def _normalize_json_scalar(x):
    if x is None or isinstance(x, (str, bool, int)):
        return x
    # The project template asks the verifier to accept harmless numeric spelling
    # differences such as 100 versus 100.0. Non-integral floats are still invalid.
    if isinstance(x, float):
        if not math.isfinite(x) or not x.is_integer():
            raise ValueError("JSON contains a non-integral numeric value")
        return int(x)
    if isinstance(x, list):
        return [_normalize_json_scalar(v) for v in x]
    if isinstance(x, dict):
        return {k: _normalize_json_scalar(v) for k, v in x.items()}
    raise ValueError("unsupported JSON value")


def load_json(path: Path):
    # Open once with O_NOFOLLOW and validate that exact descriptor. This avoids a
    # symlink-swap race between checking the path and reading its contents.
    flags = os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0)
    try:
        fd = os.open(path, flags)
    except FileNotFoundError:
        raise ValueError("missing output")
    except OSError as e:
        raise ValueError(f"cannot open output safely: {e}")
    try:
        st = os.fstat(fd)
        if not stat.S_ISREG(st.st_mode):
            raise ValueError("output must be a regular file")
        with os.fdopen(fd, "r", encoding="utf-8") as f:
            fd = -1
            raw = json.load(
                f,
                object_pairs_hook=no_dupes,
                parse_constant=lambda s: (_ for _ in ()).throw(
                    ValueError("non-finite JSON token")
                ),
            )
    finally:
        if fd >= 0:
            os.close(fd)
    return _normalize_json_scalar(raw)


def _fraction_equal(got, expected):
    """Return exact-item accuracy for two sequences, including length mistakes."""
    if not isinstance(got, list):
        return 0.0
    n = max(len(got), len(expected))
    if n == 0:
        return 1.0
    good = sum(1 for a, b in zip(got, expected) if a == b)
    return good / n


def _records_by_id(records, key):
    if not isinstance(records, list):
        return None
    out = {}
    for r in records:
        if not isinstance(r, dict) or key not in r:
            return None
        rid = r[key]
        if not isinstance(rid, int) or rid in out:
            return None
        out[rid] = r
    return out


def _record_accuracy(got, expected, key):
    gm = _records_by_id(got, key)
    em = _records_by_id(expected, key)
    if gm is None or em is None:
        return 0.0
    ids = set(gm) | set(em)
    if not ids:
        return 1.0
    return sum(gm.get(i) == em.get(i) for i in ids) / len(ids)


def compare_case(got, expected):
    """
    Compare only the scientific contract, not JSON formatting or redundant echoes.

    Compare the exact scientific records while ignoring JSON formatting and key order.
    Every graded quantity is a discrete invariant over F2 (dimensions, ranks,
    zigzag barcodes, and Kronecker invariants of continuation relations), so hidden
    cases require exact agreement on each record.
    """
    diag = {"query_metrics": []}
    if not isinstance(got, dict):
        return False, {"reason": "top-level output is not an object"}
    if got.get("case_id") != expected.get("case_id"):
        return False, {"reason": "case_id mismatch"}
    gq = got.get("queries")
    eq = expected.get("queries")
    if not isinstance(gq, list) or len(gq) != len(eq):
        return False, {"reason": "query count mismatch"}

    case_ok = True
    for gi, ei in zip(gq, eq):
        if not isinstance(gi, dict) or gi.get("query_id") != ei.get("query_id"):
            case_ok = False
            diag["query_metrics"].append({"query_id": ei.get("query_id"), "reason": "query_id mismatch"})
            continue

        sel_acc = _fraction_equal(
            gi.get("selected_morse_ids_by_frame"), ei["selected_morse_ids_by_frame"]
        )
        conley_dim_acc = _fraction_equal(
            gi.get("conley_node_dimensions_F2"), ei["conley_node_dimensions_F2"]
        )
        conley_rank_acc = _record_accuracy(
            gi.get("conley_generalized_rank_queries"),
            ei["conley_generalized_rank_queries"],
            "window_id",
        )
        conley_barcode_acc = _fraction_equal(
            gi.get("conley_zigzag_barcodes_F2"), ei["conley_zigzag_barcodes_F2"]
        )
        conley_loop_acc = _record_accuracy(
            gi.get("conley_loop_signatures"), ei["conley_loop_signatures"], "loop_id"
        )
        conley_holonomy_acc = _record_accuracy(
            gi.get("conley_holonomy_word_signatures"), ei["conley_holonomy_word_signatures"], "probe_id"
        )
        graph_dim_acc = _fraction_equal(
            gi.get("morse_graph_node_dimensions_F2"), ei["morse_graph_node_dimensions_F2"]
        )
        graph_rank_acc = _record_accuracy(
            gi.get("morse_graph_generalized_rank_queries"),
            ei["morse_graph_generalized_rank_queries"],
            "window_id",
        )
        graph_barcode_acc = _fraction_equal(
            gi.get("morse_graph_zigzag_barcodes_F2"), ei["morse_graph_zigzag_barcodes_F2"]
        )
        graph_loop_acc = _record_accuracy(
            gi.get("morse_graph_loop_signatures"),
            ei["morse_graph_loop_signatures"],
            "loop_id",
        )
        graph_holonomy_acc = _record_accuracy(
            gi.get("morse_graph_holonomy_word_signatures"),
            ei["morse_graph_holonomy_word_signatures"],
            "probe_id",
        )

        # Every reported invariant is discrete and exactly determined by the public
        # contract. Hidden cases therefore require every scientific record exactly.
        q_ok = (
            sel_acc == 1.0
            and conley_dim_acc == 1.0
            and graph_dim_acc == 1.0
            and conley_rank_acc == 1.0
            and conley_barcode_acc == 1.0
            and conley_loop_acc == 1.0
            and conley_holonomy_acc == 1.0
            and graph_rank_acc == 1.0
            and graph_barcode_acc == 1.0
            and graph_loop_acc == 1.0
            and graph_holonomy_acc == 1.0
        )
        case_ok &= q_ok
        diag["query_metrics"].append(
            {
                "query_id": ei["query_id"],
                "passed": q_ok,
                "selection_accuracy": round(sel_acc, 6),
                "conley_dimension_accuracy": round(conley_dim_acc, 6),
                "conley_rank_record_accuracy": round(conley_rank_acc, 6),
                "conley_barcode_degree_accuracy": round(conley_barcode_acc, 6),
                "conley_loop_record_accuracy": round(conley_loop_acc, 6),
                "conley_holonomy_record_accuracy": round(conley_holonomy_acc, 6),
                "graph_dimension_accuracy": round(graph_dim_acc, 6),
                "graph_rank_record_accuracy": round(graph_rank_acc, 6),
                "graph_barcode_degree_accuracy": round(graph_barcode_acc, 6),
                "graph_loop_record_accuracy": round(graph_loop_acc, 6),
                "graph_holonomy_record_accuracy": round(graph_holonomy_acc, 6),
            }
        )
    return case_ok, diag


def _syntax_check(path: Path):
    # Compile in memory. py_compile tries to write __pycache__ beside the artifact,
    # which is not guaranteed to be writable in a separate verifier mount.
    src = path.read_text(encoding="utf-8")
    compile(src, str(path), "exec")


def _kill_unprivileged_leftovers():
    # The verifier container has no service running as uid 65534. Kill any child
    # that tried to detach from the main process group so it cannot persist into
    # the next hidden case.
    proc = Path("/proc")
    for entry in proc.iterdir():
        if not entry.name.isdigit():
            continue
        try:
            status = (entry / "status").read_text(errors="ignore")
            uid_line = next((x for x in status.splitlines() if x.startswith("Uid:")), "")
            fields = uid_line.split()
            if len(fields) >= 2 and int(fields[1]) == 65534:
                os.kill(int(entry.name), signal.SIGKILL)
        except (FileNotFoundError, ProcessLookupError, PermissionError, ValueError, StopIteration):
            pass


def run_one(artifact: Path, hidden: Path, timeout_s: int = CASE_TIMEOUT_S):
    run_dir = Path(tempfile.mkdtemp(prefix="cm_hidden_"))
    try:
        os.chmod(run_dir, 0o755)
        inp = run_dir / "input.json"
        out = run_dir / "output.json"
        solver = run_dir / "solver.py"
        shutil.copyfile(hidden, inp)
        shutil.copyfile(artifact, solver)
        os.chmod(inp, 0o644)
        os.chmod(solver, 0o555)
        # The child runs as nobody. Give it write permission only to this sandbox.
        os.chmod(run_dir, 0o777)

        env = {
            "PATH": "/usr/local/bin:/usr/bin:/bin",
            "HOME": str(run_dir),
            "PYTHONNOUSERSITE": "1",
            "PYTHONDONTWRITEBYTECODE": "1",
            "PYTHONSAFEPATH": "1",
            "PYTHONPATH": "",
        }
        cmd = ["python3", str(solver), str(inp), str(out)]
        t0 = time.monotonic()
        # Capture child output in anonymous regular files rather than pipes.
        # A detached grandchild can inherit a pipe and keep communicate() blocked
        # after the solver process exits.  Regular-file redirection makes waiting
        # depend only on the solver process, while still keeping diagnostics hidden.
        with tempfile.TemporaryFile(mode="w+b") as stdout_file, tempfile.TemporaryFile(mode="w+b") as stderr_file:
            proc = subprocess.Popen(
                cmd,
                cwd=run_dir,
                stdout=stdout_file,
                stderr=stderr_file,
                env=env,
                start_new_session=True,
                user=65534,
                group=65534,
                extra_groups=(),
            )
            try:
                proc.wait(timeout=timeout_s)
            except subprocess.TimeoutExpired:
                try:
                    os.killpg(proc.pid, signal.SIGKILL)
                except ProcessLookupError:
                    pass
                try:
                    proc.wait(timeout=5)
                except subprocess.TimeoutExpired:
                    proc.kill()
                    proc.wait()
                _kill_unprivileged_leftovers()
                stderr_file.flush()
                stderr_file.seek(0)
                se = stderr_file.read().decode("utf-8", errors="replace")
                return False, timeout_s, "timeout", se[-2000:], {}

            # Remove any process that deliberately detached before reading the
            # capture files or advancing to the next hidden case.
            _kill_unprivileged_leftovers()
            stderr_file.flush()
            stderr_file.seek(0)
            se = stderr_file.read().decode("utf-8", errors="replace")

        elapsed = time.monotonic() - t0
        if proc.returncode != 0:
            return False, elapsed, f"exit {proc.returncode}", se[-2000:], {}
        try:
            got = load_json(out)
        except Exception as e:
            return False, elapsed, f"invalid output: {type(e).__name__}: {e}", se[-2000:], {}
        ok, compare_diag = compare_case(got, ORACLES[hidden.name])
        return ok, elapsed, "ok" if ok else "scientific result mismatch", se[-2000:], compare_diag
    finally:
        _kill_unprivileged_leftovers()
        shutil.rmtree(run_dir, ignore_errors=True)


def evaluate(artifact: Path = ARTIFACT):
    diagnostics = {
        "artifact": str(artifact),
        "cases": [],
        "passed_cases": 0,
        "required_cases": REQUIRED_CASES,
    }
    try:
        if not artifact.exists() or artifact.is_symlink():
            return False, "trajectory_solver.py must be a regular file", diagnostics
        st = artifact.stat()
        if not stat.S_ISREG(st.st_mode):
            return False, "trajectory_solver.py must be a regular file", diagnostics
        _syntax_check(artifact)

        passed = 0
        for hidden in HIDDEN:
            ok, elapsed, msg, stderr, compare_diag = run_one(artifact, hidden)
            diagnostics["cases"].append(
                {
                    "case": hidden.stem,
                    "passed": ok,
                    "runtime_sec": round(elapsed, 3),
                    "message": msg,
                    "stderr_tail": stderr,
                    "scientific_metrics": compare_diag,
                }
            )
            passed += int(ok)
        diagnostics["passed_cases"] = passed
        ok = passed >= REQUIRED_CASES
        return ok, f"{passed}/{len(HIDDEN)} held-out cases passed; need at least {REQUIRED_CASES}", diagnostics
    except Exception as e:
        return False, f"{type(e).__name__}: {e}", diagnostics


def write_metrics(diag):
    p = Path("/logs/verifier/metrics.json")
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(diag, indent=2, sort_keys=True) + "\n", encoding="utf-8")

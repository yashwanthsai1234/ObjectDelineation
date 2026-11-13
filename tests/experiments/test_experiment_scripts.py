"""Every script has a command line: `--help` prints usage and exits 0 without touching data (B4-7): the command
modules of objdel.cli and every script under tools/.

Scripts run as modules in an isolated interpreter, as docs/contributing.md prescribes for anything that reads data/
(`python -I -m tools.<folder>.<script>`; -I drops the current directory from sys.path, so this also checks that
the editable install makes `tools` importable). One `-I` subprocess runs every script through `runpy`, so the
heavy imports (torch, TRL) are paid once. A script whose optional extra is not installed (the CI `python` job has
only `[dev]`) is skipped, not failed.

References: https://docs.python.org/3/using/cmdline.html#cmdoption-I,
https://docs.python.org/3/library/runpy.html#runpy.run_module, https://docs.python.org/3/library/argparse.html
"""
import json
import subprocess
import sys

import pytest

from tests.support.paths import ROOT

SCRIPTS = [".".join(p.relative_to(root).with_suffix("").parts)
           for root, pattern in ((ROOT, "tools/**/*.py"), (ROOT / "src", "objdel/cli/*/*.py"))
           for p in sorted(root.glob(pattern)) if p.stem != "__init__"]
OPTIONAL = {"torch", "transformers", "trl", "datasets", "peft", "accelerate", "pydicom", "pylidc"}  # extras

RUNNER = """
import contextlib, io, json, runpy, sys
path, scripts, out = sys.argv[1], sys.argv[2:], {}
for s in scripts:
    buf, sys.argv = io.StringIO(), [s, "--help"]
    try:
        with contextlib.redirect_stdout(buf):
            runpy.run_module(s, run_name="__main__", alter_sys=True)
        out[s] = {"code": "returned without exiting", "stdout": buf.getvalue()}
    except SystemExit as e:
        out[s] = {"code": e.code or 0, "stdout": buf.getvalue()}
    except ModuleNotFoundError as e:
        out[s] = {"missing": e.name}
    except Exception as e:
        out[s] = {"code": repr(e), "stdout": buf.getvalue()}
open(path, "w").write(json.dumps(out))
"""


@pytest.fixture(scope="module")
def helps(tmp_path_factory):
    out = tmp_path_factory.mktemp("help") / "help.json"
    r = subprocess.run([sys.executable, "-I", "-c", RUNNER, str(out), *SCRIPTS], cwd=ROOT, capture_output=True,
                       text=True, timeout=600)
    assert r.returncode == 0 and out.is_file(), r.stderr[-3000:]
    return json.loads(out.read_text())


def test_scripts_found():
    assert {f"objdel.cli.{s}" for s in ("eval.score", "data.prepare", "data.store", "eval.run")} | {
        "tools.checker.eval", "tools.data.verify_target_hashes"} <= set(SCRIPTS)


@pytest.mark.parametrize("script", SCRIPTS)
def test_script_help(helps, script):
    h = helps[script]
    if "missing" in h:
        if h["missing"].split(".")[0] in OPTIONAL:
            pytest.skip(f"optional dependency {h['missing']} not installed")
        pytest.fail(f"{script}: missing module {h['missing']}")
    assert h["code"] == 0, h
    assert h["stdout"].startswith("usage:")

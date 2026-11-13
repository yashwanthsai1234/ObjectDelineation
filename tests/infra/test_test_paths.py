"""E8 step 3 / E-G15: every `tests/...py` path (and `::test` name) written in live code and config names a file (and a
function) that exists, so moving a test file cannot leave a stale path behind."""

import ast
import re
import subprocess
from pathlib import Path

from tests.support.paths import ROOT

PATH_RE = re.compile(r"tests/[\w/]+\.py(?:::\w+)?")
SOURCES = ("pyproject.toml", "Makefile", ".github", "ops", "tools", "specs", "src")


def live_files() -> list[Path]:
    """Tracked files only: git-ignored local state (vast/code_files.txt from sync_up, run outputs) is not live code."""
    names = subprocess.run(["git", "ls-files", "-z", *SOURCES], cwd=ROOT, capture_output=True, text=True,
                           check=True).stdout.split("\0")
    return [ROOT / n for n in names if n and Path(n).suffix in (".py", ".sh", ".yaml", ".yml", ".toml", ".txt", "")]


def missing(paths: dict[str, list[str]]) -> list[str]:
    bad = []
    for ref, where in sorted(paths.items()):
        file, _, func = ref.partition("::")
        f = ROOT / file
        if not f.is_file():
            bad.append(f"{ref} (in {where[0]})")
        elif func and func not in {n.name for n in ast.walk(ast.parse(f.read_text()))
                                   if isinstance(n, ast.FunctionDef)}:
            bad.append(f"{ref}: no function {func} (in {where[0]})")
    return bad


def test_every_test_path_in_live_code_exists():
    refs: dict[str, list[str]] = {}
    for f in live_files():
        for m in PATH_RE.findall(f.read_text(errors="ignore")):
            refs.setdefault(m, []).append(str(f.relative_to(ROOT)))
    assert refs and missing(refs) == []


def test_a_stale_path_is_reported(tmp_path):
    assert missing({"tests/no_such_file.py": ["x"]}) == ["tests/no_such_file.py (in x)"]
    assert missing({"tests/infra/test_test_paths.py::nope": ["x"]})[0].endswith("no function nope (in x)")

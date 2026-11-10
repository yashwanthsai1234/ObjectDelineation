"""E0: the import-linter contracts and the ruff TID251 ban bind, shown by planting each violation in a scratch copy
of the tree (src/ plus pyproject.toml) and running the real tool."""
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
# The editable install redirects `objdel` to the real tree; the scratch copy must win, so the runner drops that finder.
RUNNER = """
import sys
sys.meta_path[:] = [f for f in sys.meta_path if "skbc" not in repr(f) and "ditable" not in repr(f)]
sys.path.insert(0, "src")
from importlinter.cli import lint_imports_command
sys.argv = ["lint-imports", "--no-cache", "--no-logo"]
lint_imports_command()
"""


def lint_imports(tree: Path) -> subprocess.CompletedProcess:
    return subprocess.run([sys.executable, "-c", RUNNER], cwd=tree, capture_output=True, text=True)


@pytest.fixture(scope="module")
def scratch(tmp_path_factory) -> Path:
    tree = tmp_path_factory.mktemp("tree")
    shutil.copytree(ROOT / "src", tree / "src", ignore=shutil.ignore_patterns("__pycache__", "*.so"))
    shutil.copy(ROOT / "pyproject.toml", tree)
    return tree


def plant(tree: Path, rel: str, line: str):
    """Append `line` to a module of the scratch tree; the returned function restores the module."""
    path = tree / "src/objdel" / rel
    before = path.read_text()
    path.write_text(before + "\n" + line + "\n")
    return lambda: path.write_text(before)


def test_the_real_tree_keeps_every_contract():
    r = lint_imports(ROOT)
    assert r.returncode == 0, r.stdout[-3000:]


@pytest.mark.parametrize(
    ("rel", "line", "broken"),
    [
        ("data/__init__.py", "import torch", "never imports the ML stack"),  # (a)
        ("data/__init__.py", "def f():\n    import torch", "never imports the ML stack"),  # (d) function level
        ("geom/__init__.py", "import objdel.train", "objdel layers (current packages) BROKEN"),  # (b)
        # (e): the ignored edge's target from a module not on the list
        ("data/lidc.py", "import objdel.train.dataset", "objdel layers (current packages) BROKEN"),
    ],
)
def test_a_planted_violation_breaks_a_contract(scratch, rel, line, broken):
    undo = plant(scratch, rel, line)
    try:
        r = lint_imports(scratch)
    finally:
        undo()
    assert r.returncode == 1 and broken in r.stdout, r.stdout[-3000:]


def test_ruff_bans_the_private_trl_module_outside_the_allowed_files():
    code = "from trl.trainer import grpo_trainer\n"
    ruff = [sys.executable, "-m", "ruff", "check", "--no-cache", "--select", "TID251", "--stdin-filename"]
    bad = subprocess.run([*ruff, "src/objdel/rl_new.py", "-"], cwd=ROOT, input=code, capture_output=True, text=True)
    ok = subprocess.run([*ruff, "src/objdel/train/fast.py", "-"], cwd=ROOT, input=code, capture_output=True, text=True)
    assert bad.returncode == 1 and "TID251" in bad.stdout, bad.stdout
    assert ok.returncode == 0, ok.stdout


SCRATCH_PYPROJECT = """
[tool.importlinter]
root_package = "pkg"

[[tool.importlinter.contracts]]
name = "optional and independent"
type = "layers"
containers = ["pkg"]
layers = ["top", "(x) | (y) | (absent)", "base"]
"""


def test_optional_and_independent_layer_syntax_is_accepted(tmp_path):
    """(c): 2.16 accepts `(x) | (y)` (optional, independent) in one layer; a missing optional module is no error,
    and the siblings may not import each other. So the target contract can say `(rl) | (sft)` before E5."""
    pkg = tmp_path / "src/pkg"
    pkg.mkdir(parents=True)
    for name in ("__init__", "top", "x", "y", "base"):
        (pkg / f"{name}.py").write_text("import pkg.base\n" if name in ("x", "y") else "")
    (tmp_path / "pyproject.toml").write_text(SCRATCH_PYPROJECT)
    ok = lint_imports(tmp_path)
    (pkg / "x.py").write_text("import pkg.y\n")
    bad = lint_imports(tmp_path)
    assert ok.returncode == 0, ok.stdout[-3000:]
    assert bad.returncode == 1 and "BROKEN" in bad.stdout, bad.stdout[-3000:]

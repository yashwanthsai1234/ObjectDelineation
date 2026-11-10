"""E0/E8.4: the import-linter contracts and the ruff TID251 ban bind, shown by planting each violation in a scratch
copy of the tree (src/ plus pyproject.toml) and running the real tool (I1); the contract config itself is the
exhaustive E8 target form (U1)."""
import shutil
import subprocess
import sys
import tomllib
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


LAYERS = ["infra", "rl | sft", "policy", "data", "verify", "geom", "core"]
EXTENSIONS = ["_gpucheck"]  # compiled modules import-linter cannot read (exhaustive_ignores)
PURE_CORE = ["objdel.core", "objdel.geom", "objdel.verify", "objdel.data"]
LAYER_NAMES = {n.strip() for layer in LAYERS for n in layer.split("|")}


def contracts(pyproject: Path = ROOT / "pyproject.toml") -> dict:
    return {c["id"]: c for c in tomllib.loads(pyproject.read_text())["tool"]["importlinter"]["contracts"]}


def children(src: Path = ROOT / "src/objdel") -> set:
    """The packages and modules (extension modules included) directly under src/objdel, by import name."""
    out = {p.name for p in src.iterdir() if p.is_dir() and (p / "__init__.py").exists()}
    out |= {p.name.split(".")[0] for p in src.iterdir() if p.is_file() and p.suffix in (".py", ".so", ".pyd")}
    return out - {"__init__"}


def test_layers_contract_is_the_exhaustive_target_form():
    """U1 (E8.4): one exhaustive `layers` contract over the eight layers, no ignore_imports, guidance on every
    contract, and pure-core over the four ML-free packages. A loosened config disables the architecture silently."""
    cs = contracts()
    layers = cs["layers"]
    assert (layers["type"], layers["containers"], layers["layers"]) == ("layers", ["objdel"], LAYERS)
    assert layers["exhaustive"] is True and layers["exhaustive_ignores"] == EXTENSIONS
    assert not any("(" in layer for layer in layers["layers"]), "every layer is non-optional"
    assert [c["id"] for c in cs.values() if c.get("ignore_imports")] == []
    assert [c["id"] for c in cs.values() if not c.get("broken_contract_guidance")] == []
    assert cs["pure-core"]["type"] == "forbidden" and cs["pure-core"]["source_modules"] == PURE_CORE


def test_src_objdel_children_are_exactly_the_layers():
    """U1 precondition: the children of src/objdel are the eight layers plus the ignored extensions; an extra child
    is named here, before lint-imports reports it as an unlisted module."""
    kids = children()
    assert sorted(kids - LAYER_NAMES - set(EXTENSIONS)) == [] and sorted(LAYER_NAMES - kids) == []


def test_the_real_tree_keeps_every_contract():
    r = lint_imports(ROOT)
    assert r.returncode == 0, r.stdout[-3000:]


LAYERS_BROKEN = "objdel layers BROKEN"
PURE_BROKEN = "never imports the ML stack BROKEN"


@pytest.mark.parametrize(
    ("rel", "line", "broken"),
    [
        ("data/__init__.py", "import torch", PURE_BROKEN),  # (a)
        ("data/__init__.py", "def f():\n    import torch", PURE_BROKEN),  # (d) function level
        ("core/__init__.py", "import torch", PURE_BROKEN),  # E8.4: core joins pure-core
        ("geom/__init__.py", "import objdel.rl", LAYERS_BROKEN),  # (b)
        ("core/__init__.py", "import objdel.rl", LAYERS_BROKEN),
        ("verify/__init__.py", "import objdel.data", LAYERS_BROKEN),  # verify below data
        ("verify/__init__.py", "import objdel.policy", LAYERS_BROKEN),
        ("geom/__init__.py", "import objdel.data", LAYERS_BROKEN),
        # E5 gate: rl and sft are independent, and policy sits below both
        ("sft/__init__.py", "import objdel.rl", LAYERS_BROKEN),
        ("rl/__init__.py", "import objdel.sft", LAYERS_BROKEN),
        ("policy/__init__.py", "import objdel.rl", LAYERS_BROKEN),
        ("core/__init__.py", "import objdel.infra", LAYERS_BROKEN),  # E8: infra is the top layer
    ],
)
def test_a_planted_violation_breaks_a_contract(scratch, rel, line, broken):
    """I1: the output names the broken contract and prints its guidance."""
    undo = plant(scratch, rel, line)
    try:
        r = lint_imports(scratch)
    finally:
        undo()
    assert r.returncode == 1 and broken in r.stdout, r.stdout[-3000:]
    guidance = contracts(scratch / "pyproject.toml")["layers" if broken == LAYERS_BROKEN else "pure-core"]
    assert guidance["broken_contract_guidance"] in " ".join(r.stdout.split()), r.stdout[-3000:]


def test_an_unlisted_package_breaks_the_exhaustive_contract(scratch):
    """I1 (E8.4): a new package under objdel that is in no layer fails `exhaustive = true`."""
    new = scratch / "src/objdel/newpkg"
    new.mkdir()
    (new / "__init__.py").write_text("")
    try:
        r = lint_imports(scratch)
    finally:
        shutil.rmtree(new)
    assert r.returncode == 1 and LAYERS_BROKEN in r.stdout and "newpkg" in r.stdout, r.stdout[-3000:]


def test_ruff_bans_the_private_trl_module_outside_the_allowed_files():
    code = "from trl.trainer import grpo_trainer\n"
    ruff = [sys.executable, "-m", "ruff", "check", "--no-cache", "--select", "TID251", "--stdin-filename"]
    bad = subprocess.run([*ruff, "src/objdel/rl_new.py", "-"], cwd=ROOT, input=code, capture_output=True, text=True)
    ok = subprocess.run([*ruff, "src/objdel/rl/trl_adapter.py", "-"], cwd=ROOT, input=code, capture_output=True,
                        text=True)
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

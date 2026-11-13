"""U34 (E2.5): objdel.core.provenance finds the repo root by a search, not a fixed parents[n], and the code hash covers
the source tree (a hash of 0 files would let the launch gate miss code edits)."""
import ast
import importlib.util
import shutil
from pathlib import Path

import pytest

from objdel.core import provenance

SRC = Path(provenance.__file__).resolve().parents[1]  # src/objdel (this test's own path search is allowed)
E0_FILES = 73  # code_hashes()["files"] at be52b1c [plan E2.5 step 5]
# __file__-relative roots not yet on repo_root(), each with the task that removes it (none left)
PENDING: dict[str, str] = {}


def _copy_module(dest: Path):
    """provenance.py copied to dest/src/objdel/core and imported from there (its own __file__ drives the search)."""
    core = dest / "src" / "objdel" / "core"
    core.mkdir(parents=True)
    shutil.copy(provenance.__file__, core / "provenance.py")
    spec = importlib.util.spec_from_file_location("provenance_copy", core / "provenance.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_repo_root_holds_pyproject_and_src(monkeypatch):
    monkeypatch.delenv("OBJDEL_REPO", raising=False)
    root = provenance.repo_root()
    assert (root / "pyproject.toml").is_file() and (root / "src" / "objdel").is_dir()
    assert root == SRC.parents[1]


def test_env_wins_and_needs_pyproject(tmp_path, monkeypatch):
    monkeypatch.setenv("OBJDEL_REPO", str(tmp_path))
    with pytest.raises(FileNotFoundError, match="OBJDEL_REPO"):
        provenance.repo_root()
    (tmp_path / "pyproject.toml").write_text("")
    assert provenance.repo_root() == tmp_path


def test_search_finds_a_tree_without_git_from_any_cwd(tmp_path, monkeypatch):
    monkeypatch.delenv("OBJDEL_REPO", raising=False)
    mod = _copy_module(tmp_path)
    (tmp_path / "pyproject.toml").write_text("")
    assert not (tmp_path / ".git").exists()
    assert mod.repo_root() == tmp_path
    monkeypatch.chdir(tmp_path / "src")
    assert mod.repo_root() == tmp_path


def test_no_tree_and_no_env_is_a_clear_error(tmp_path, monkeypatch):
    monkeypatch.delenv("OBJDEL_REPO", raising=False)
    mod = _copy_module(tmp_path)  # no pyproject.toml anywhere above
    with pytest.raises(FileNotFoundError, match="pyproject.toml"):
        mod.repo_root()


def test_code_hash_covers_the_tree(monkeypatch):
    monkeypatch.delenv("OBJDEL_REPO", raising=False)
    h = provenance.code_hashes()
    assert h["files"] >= E0_FILES
    assert len(h["objdel_and_script_sha256"]) == 64


def test_an_edited_source_file_changes_the_hash(tmp_path, monkeypatch):
    root = provenance.repo_root()
    for p in sorted({p for g in provenance.CODE for p in root.glob(g)}):
        (tmp_path / p.relative_to(root)).parent.mkdir(parents=True, exist_ok=True)
        shutil.copy(p, tmp_path / p.relative_to(root))
    monkeypatch.setenv("OBJDEL_REPO", str(tmp_path))
    before = provenance.code_hashes()
    assert before["git_head"] is None and before["files"] >= E0_FILES  # a box tree: no .git
    with open(tmp_path / "src/objdel/core/provenance.py", "a") as f:
        f.write("\n")
    assert provenance.code_hashes()["objdel_and_script_sha256"] != before["objdel_and_script_sha256"]


def _file_roots(path: Path) -> list[int]:
    """Lines of __file__-relative paths: `parents[`, `.parent /`, `with_name` in an expression that names __file__."""
    hits = []
    for node in ast.walk(ast.parse(path.read_text())):
        if not isinstance(node, ast.expr) or not any(isinstance(n, ast.Name) and n.id == "__file__"
                                                     for n in ast.walk(node)):
            continue
        if isinstance(node, ast.Subscript) and isinstance(node.value, ast.Attribute) and node.value.attr == "parents":
            hits.append(node.lineno)
        elif isinstance(node, ast.Attribute) and node.attr == "with_name":
            hits.append(node.lineno)
        elif isinstance(node, ast.BinOp) and isinstance(node.op, ast.Div) and isinstance(node.left, ast.Attribute) \
                and node.left.attr == "parent":
            hits.append(node.lineno)
    return hits


def test_no_file_relative_roots_in_src():
    found = {p.relative_to(SRC).as_posix(): _file_roots(p) for p in SRC.rglob("*.py")}
    found = {k: v for k, v in found.items() if v and k != "core/provenance.py"}
    assert sorted(found) == sorted(PENDING), found


def test_file_root_guard_sees_a_planted_violation(tmp_path):
    bad = tmp_path / "bad.py"
    bad.write_text("from pathlib import Path\nA = Path(__file__).resolve().parents[2]\n"
                   "B = Path(__file__).parent / 'x'\nC = Path(__file__).with_name('y')\n")
    assert _file_roots(bad) == [2, 3, 4]


def test_reference_dir_data_root_and_run_dir_follow_the_environment(tmp_path, monkeypatch):
    (tmp_path / "pyproject.toml").write_text("")
    monkeypatch.setenv("OBJDEL_REPO", str(tmp_path))
    assert provenance.reference_dir() == tmp_path / "reference"
    for var, get, default in (("OBJDEL_DATA_ROOT", provenance.data_root, "data"),
                              ("OBJDEL_RUN_DIR", provenance.run_dir, "runs")):
        monkeypatch.delenv(var, raising=False)
        assert get() == Path(default) and get("other") == Path("other")
        monkeypatch.setenv(var, "")  # an empty value counts as unset
        assert get() == Path(default)
        monkeypatch.setenv(var, "/x/y")
        assert get() == Path("/x/y") and get("other") == Path("/x/y")


def test_val_halves_are_read_from_the_repo_root_the_environment_names(tmp_path, monkeypatch):
    from objdel.data import splits

    (tmp_path / "pyproject.toml").write_text("")
    (tmp_path / "reference").mkdir()
    (tmp_path / "reference" / "val_halves.json").write_text('{"halves": {"LIDC-IDRI-0001": "shape"}}')
    monkeypatch.setenv("OBJDEL_REPO", str(tmp_path))
    splits._val_halves.cache_clear()
    try:
        assert splits.val_half("LIDC-IDRI-0001") == "shape"
    finally:
        splits._val_halves.cache_clear()

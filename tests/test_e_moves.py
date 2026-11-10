"""Phase E move tool (experiments/phase2/tools/e_moves.py): path and module names, the import rewrite on a planted
tree, the long-line wrap, and idempotence (a second run changes nothing)."""

import subprocess
import sys
from pathlib import Path

from experiments.phase2.tools import e_moves as E

ROOT = Path(__file__).resolve().parents[1]


def test_module_and_path_names_round_trip():
    assert E.mod_of("src/objdel/train/fast.py") == "objdel.train.fast"
    assert E.mod_of("src/objdel/reward/__init__.py") == "objdel.reward"
    assert E.mod_of("tests/oracles/cuda_parity.py") == "tests.oracles.cuda_parity"
    assert E.path_of("objdel.geom.rings") == "src/objdel/geom/rings.py"


def test_every_split_name_has_one_home():
    table = E.json.loads(E.TABLE.read_text())
    for src, targets in table["splits"].items():
        names = [n for ns in targets.values() for n in ns]
        assert len(names) == len(set(names)), src


def test_wrap_keeps_code_and_shortens_prose():
    long_doc = '"""' + "word " * 30 + '"""\n'
    text = (f"def f():\n    {long_doc}    x = g(aaaa, {'b' * 50}, {'c' * 50})  # {'note ' * 4}\n"
            f"    return h(aaaa, {'d' * 60}, {'e' * 60})\n")
    new, left = E.wrap_prose(text)
    assert left == [] and all(len(x) <= 120 for x in new.splitlines())
    assert sorted(new.split()) == sorted(text.split())  # the same words, only moved (the comment goes above)
    assert new.splitlines()[3].lstrip().startswith("# note")
    compile(new, "x", "exec")


def test_planted_tree_moves_once(tmp_path):
    """A tiny tree with one moved module, one split module and a caller: the run moves, splits and rewrites; a second
    run prints nothing to do."""
    root = tmp_path
    subprocess.run(["git", "init", "-q", str(root)], check=True)
    (root / "src/objdel/geometry").mkdir(parents=True)
    (root / "src/objdel/geometry/__init__.py").write_text('"""g."""\n')
    (root / "src/objdel/geometry/text.py").write_text("SCHEMAS = {}\n")
    (root / "src/objdel/train").mkdir()
    (root / "src/objdel/train/__init__.py").write_text('"""t."""\n')
    (root / "src/objdel/train/rewards.py").write_text(
        "def stop_ids():\n    return 1\n\n\ndef _score():\n    return stop_ids()\n")
    (root / "tests").mkdir()
    (root / "tests/test_x.py").write_text(
        "from objdel.geometry.text import SCHEMAS\nfrom objdel.train import rewards\n\n\n"
        "def test_x():\n    assert rewards.stop_ids() == 1 and rewards._score() == 1 and SCHEMAS == {}\n")
    subprocess.run(["git", "-C", str(root), "add", "-A"], check=True)
    table = {"files": {"src/objdel/geometry/text.py": "src/objdel/geom/text.py"},
             "splits": {"src/objdel/train/rewards.py": {"objdel.policy.models": ["stop_ids"],
                                                        "objdel.rl.reward": ["_score"]}},
             "drop_reexports": [], "delete_packages": ["src/objdel/geometry", "src/objdel/train"],
             "packages": {"src/objdel/geom/__init__.py": '"""g."""\n', "src/objdel/policy/__init__.py": '"""p."""\n',
                          "src/objdel/rl/__init__.py": '"""r."""\n'},
             "engine": {"path": "src/objdel/geom/engine.py", "names": []}, "gpucheck_keep": [],
             "scan_roots": ["src", "tests"], "text_globs": ["*.py"], "text_extra": [], "docstring_to": {}}
    tab = root / "table.json"
    tab.write_text(E.json.dumps(table))
    code = (f"import sys; sys.path.insert(0, {str(ROOT)!r}); from pathlib import Path; "
            "from experiments.phase2.tools import e_moves as E; "
            f"E.TABLE = Path({str(tab)!r}); sys.exit(E.main(sys.argv[1:]))")
    run = [sys.executable, "-I", "-c", code, "--root", str(root)]
    first = subprocess.run(run, capture_output=True, text=True, cwd=root)
    assert first.returncode == 0, first.stdout + first.stderr
    assert not (root / "src/objdel/geometry").exists() and not (root / "src/objdel/train").exists()
    assert "def stop_ids" in (root / "src/objdel/policy/models.py").read_text()
    assert "def _score" in (root / "src/objdel/rl/reward.py").read_text()
    test = (root / "tests/test_x.py").read_text()
    assert "objdel.geom.text" in test and "objdel.train" not in test and "objdel.geometry" not in test
    again = subprocess.run(run, capture_output=True, text=True, cwd=root)
    assert again.returncode == 0 and "nothing to do" in again.stdout, again.stdout + again.stderr


def test_stages_move_a_split_module_in_parts(tmp_path):
    """--stage E3b moves only the names whose new module is in verify; the rest stays in the old module, which imports
    the moved name it uses. --stage E5 then moves the rest, and the caller ends on the new modules."""
    root = tmp_path
    subprocess.run(["git", "init", "-q", str(root)], check=True)
    (root / "src/objdel/train").mkdir(parents=True)
    (root / "src/objdel/train/__init__.py").write_text('"""t."""\n')
    (root / "src/objdel/train/rewards.py").write_text(
        '"""Rewards."""\n\n\ndef terms():\n    return 1\n\n\ndef _score():\n    return terms()\n')
    (root / "tests").mkdir()
    (root / "tests/test_x.py").write_text(
        "from objdel.train import rewards as R\n\n\ndef test_x():\n    assert R.terms() == 1 and R._score() == 1\n")
    subprocess.run(["git", "-C", str(root), "add", "-A"], check=True)
    table = {"files": {}, "splits": {"src/objdel/train/rewards.py": {"objdel.verify.score": ["terms"],
                                                                     "objdel.rl.reward": ["_score"]}},
             "drop_reexports": [], "delete_packages": ["src/objdel/train"],
             "packages": {"src/objdel/verify/__init__.py": '"""v."""\n', "src/objdel/rl/__init__.py": '"""r."""\n'},
             "engine": {"path": "src/objdel/geom/engine.py", "names": []}, "gpucheck_keep": [],
             "scan_roots": ["src", "tests"], "text_globs": ["*.py"], "text_extra": [], "docstring_to": {}}
    tab = root / "table.json"
    tab.write_text(E.json.dumps(table))
    code = (f"import sys; sys.path.insert(0, {str(ROOT)!r}); from pathlib import Path; "
            "from experiments.phase2.tools import e_moves as E; "
            f"E.TABLE = Path({str(tab)!r}); sys.exit(E.main(sys.argv[1:]))")
    for stage in ("E3b", "E5"):
        run = subprocess.run([sys.executable, "-I", "-c", code, "--root", str(root), "--stage", stage],
                             capture_output=True, text=True, cwd=root)
        assert run.returncode == 0, run.stdout + run.stderr
        if stage == "E3b":
            old = (root / "src/objdel/train/rewards.py").read_text()
            assert "def _score" in old and "def terms" not in old and "from objdel.verify.score import terms" in old
            assert not (root / "src/objdel/rl").exists()
            assert "from objdel.verify import score" in (root / "tests/test_x.py").read_text()
    assert not (root / "src/objdel/train").exists()
    assert "def _score" in (root / "src/objdel/rl/reward.py").read_text()
    test = (root / "tests/test_x.py").read_text()
    assert "objdel.train" not in test and "verify_score.terms()" in test and "R._score()" in test
    paths = [str(root / "src"), str(root / "tests")]
    check = f"import sys; sys.path[:0] = {paths!r}; import test_x; test_x.test_x()"
    test_run = subprocess.run([sys.executable, "-I", "-S", "-c", check], capture_output=True, text=True)
    assert test_run.returncode == 0, test_run.stderr

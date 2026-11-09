"""CLEAN phase hygiene tests (joined into tests/test_repo_hygiene.py at the merge)."""

import ast
import json
import re
import subprocess
import sys
import tomllib
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
REWARDCORE = {f"rewardcore/{n}.py" for n in ("advantages", "batch", "fixtures", "reward", "topology")}


def tracked(*paths: str) -> list[str]:  # tracked files, as the working tree has them
    out = subprocess.run(["git", "ls-files", "-z", *paths], cwd=ROOT, capture_output=True, text=True, check=True)
    return [f for f in out.stdout.split("\0") if f]


def test_vendor_copies_deleted_and_nothing_names_them():
    """CLEAN3: only the 5 vendor files that differ from their twins may remain (CLEAN5 moved them with the lab)."""
    left = {f.split("/vendor/", 1)[1] for f in tracked("archive/phase1_lab/t_stack/vendor") if "/vendor/" in f}
    assert left <= REWARDCORE
    names = [f for f in tracked("src", "tests", "experiments/phase2") if f != "tests/test_repo_hygiene_clean.py"
             and "t_stack/vendor" in (ROOT / f).read_text(errors="ignore")]  # fmt: skip
    assert names == []


# --- CLEAN5: the archive moves ---------------------------------------------------------------------------------------

RETIRED = ("experiments/phase0", "experiments/phase1/lab", "experiments/phase2/lab")
GONE_SCRIPTS = ("census_sizerule", "eval_demo", "sample_eval", "smoke_grpo")
CLEAN5_ARCHIVES = ("phase0", "phase1_lab", "phase1_scripts", "phase2_lab")
HEADINGS = ("## What it was", "## What replaced it", "## Retiring commit", "## How to restore")
TAG = "pre-clean-2026-11-08"
EXEMPT_DOCS = ("docs/changelog.md", "docs/research/", "docs/implementation2.0/")  # frozen records and the plan itself
# Phase 0 records moved to docs/research/records/phase0/ (owner question 3 default); the plan keeps them in docs/.
STALE = ("experiments/phase0", "experiments/phase1/lab", "experiments/phase2/lab")


def _exempt(f: str) -> bool:
    name = f.rsplit("/", 1)[-1]
    return (f.startswith(EXEMPT_DOCS) or "/impl/" in f or f.startswith("experiments/phase2/perf/")
            or (f.startswith("experiments/phase") and f.count("/") == 2 and name.endswith(".md")))  # fmt: skip


def stale_paths(files: list[str], root: Path) -> list[str]:
    """Live text files that still name a retired tree."""
    out = []
    for f in files:
        if _exempt(f) or f.endswith((".png", ".npz", ".pdf")):
            continue
        text = (root / f).read_text(errors="ignore")
        out += [f"{f}: {s}" for s in STALE if s in text]
    return out


def import_leaks(files: list[Path], src: bool) -> list[str]:
    """tests/: imports of archive or lab modules, lab or archive on sys.path, lab paths in code.
    src/: imports of experiments or tests."""
    out = []
    for f in files:
        tree = ast.parse(f.read_text())
        scopes = (ast.Module, ast.ClassDef, ast.FunctionDef, ast.AsyncFunctionDef)
        docs = {id(n.body[0].value) for n in ast.walk(tree) if isinstance(n, scopes) and ast.get_docstring(n)}
        for n in ast.walk(tree):
            mods = (
                [a.name for a in n.names]
                if isinstance(n, ast.Import)
                else [n.module or ""]
                if isinstance(n, ast.ImportFrom)
                else []
            )
            for m in mods:
                if (src and m.split(".")[0] in ("experiments", "tests")) or (
                    not src and (m.split(".")[0] == "archive" or "lab" in m.split("."))):  # fmt: skip
                    out.append(f"{f.name}: import {m}")
            if (
                isinstance(n, ast.Call)
                and not src
                and ast.unparse(n.func).endswith(("sys.path.insert", "sys.path.append"))
            ):
                if "lab" in ast.unparse(n) or "archive" in ast.unparse(n):
                    out.append(f"{f.name}: sys.path {ast.unparse(n)}")
            if isinstance(n, ast.Constant) and isinstance(n.value, str) and id(n) not in docs and not src:
                out += [f"{f.name}: {s}" for s in ("experiments/phase1/lab", "experiments/phase2/lab") if s in n.value]
    return out


def test_retired_trees_hold_nothing():
    """CT-37: no tracked file, no folder on disk, the live phase 1 modules still import, the retired ones do not."""
    for p in RETIRED:
        assert tracked(p) == [] and not (ROOT / p).exists(), p
    probe = ("import importlib, sys\n"
             "mods = 'build_data build_targets fetch_lidc_pixels prepare_rl eval_run reward_gates bench.xarch_hashes'\n"
             "for m in mods.split():\n"
             "    importlib.import_module('experiments.phase1.' + m)\n"
             "assert not [n for n in sys.modules if n.endswith('.lab') or '.lab.' in n]\n"
             "for m in ('experiments.phase0', 'experiments.phase1.lab', 'experiments.phase2.lab'):\n"
             "    try:\n        importlib.import_module(m)\n    except ModuleNotFoundError:\n        continue\n"
             "    raise SystemExit(m + ' imports')\n")  # fmt: skip
    r = subprocess.run([sys.executable, "-I", "-c", probe], cwd=ROOT, capture_output=True, text=True)
    assert r.returncode == 0, r.stderr[-800:]


def test_phase1_scripts_archived_and_their_tests_listed():
    """CT-83 and CT-40: the four one-off scripts left the live list; each removed test is named in the README."""
    live = {p.stem for p in (ROOT / "experiments/phase1").glob("*.py")}
    assert not live & set(GONE_SCRIPTS) and {"eval_run", "prepare_rl", "build_data"} <= live
    arch = ROOT / "archive/phase1_scripts"
    assert {p.stem for p in arch.glob("*.py")} == set(GONE_SCRIPTS)
    readme = (arch / "README.md").read_text()
    names = [n.name for f in (arch / "tests").glob("test_*.py") for n in ast.parse(f.read_text()).body
             if isinstance(n, ast.FunctionDef) and n.name.startswith("test_")]  # fmt: skip
    assert len(names) == 3 and all(f"::{n}`" in readme for n in names), names
    for f in ("tests/test_experiment_scripts.py", "tests/test_eval.py"):
        live_names = {n.name for n in ast.walk(ast.parse((ROOT / f).read_text())) if isinstance(n, ast.FunctionDef)}
        assert not live_names & set(names), f
    wiring = (ROOT / "tests/test_breakit_wiring.py").read_text()
    assert "census_sizerule" not in wiring


def test_archive_readmes_follow_the_contract():
    """CT-38: four headings in order, a retiring hash that is a prefix of the tag's commit, restore names the tag."""
    full = subprocess.run(
        ["git", "rev-parse", TAG], cwd=ROOT, capture_output=True, text=True, check=True
    ).stdout.strip()
    for name in CLEAN5_ARCHIVES:
        text = (ROOT / "archive" / name / "README.md").read_text()
        pos = [text.find("\n" + h + "\n") for h in HEADINGS]
        assert all(p >= 0 for p in pos) and pos == sorted(pos), name
        retiring = text.split(HEADINGS[2])[1].split("\n## ")[0]
        assert any(full.startswith(h) for h in re.findall(r"\b[0-9a-f]{7,40}\b", retiring)), name
        assert TAG in text.split(HEADINGS[3])[1], name
    assert not (ROOT / "archive/__init__.py").exists()
    assert "run it from the tag" in (ROOT / "archive/phase1_lab/README.md").read_text()
    assert "stay on the owner's disk only" in (ROOT / "archive/phase1_lab/README.md").read_text()
    assert "stay on the owner's disk only" in (ROOT / "docs/research/records/phase0/README.md").read_text()


def test_phase0_scripts_and_records_moved():
    """CT-87: 8 python files in the archive, records beside the research notes, the oracle parser stays in tests."""
    assert len(list((ROOT / "archive/phase0").glob("*.py"))) == 8
    rec = ROOT / "docs/research/records/phase0"
    assert (rec / "results").is_dir() and len(list(rec.glob("*.md"))) >= 5
    assert "class _Lab" in (ROOT / "tests/test_geometry_oracle.py").read_text()
    for f in tracked("src"):
        for m in re.findall(r"archive/[\w./-]+\.py", (ROOT / f).read_text()):
            assert (ROOT / m).exists(), (f, m)


def test_archive_is_not_collected_or_linted():
    """CT-39: pytest only walks tests/, ruff only skips archive (key extend-exclude keeps the default excludes)."""
    cfg = tomllib.loads((ROOT / "pyproject.toml").read_text())
    assert cfg["tool"]["pytest"]["ini_options"]["testpaths"] == ["tests"]
    assert cfg["tool"]["ruff"]["extend-exclude"] == ["archive"]


def test_fixtures_moved_and_digest_inline():
    """CT-41: the phase 2 lab fixtures live in tests/fixtures; nothing under tests/ or src/ reads a lab path."""
    fx = ROOT / "tests/fixtures"
    assert (fx / "vast_snapshot.json").is_file() and list(fx.glob("topo_probe_v3e_seed*.json"))
    for f in fx.glob("*.json"):
        json.loads(f.read_text(), parse_constant=lambda c, f=f: pytest.fail(f"{f.name}: {c}"))
    digest = re.search(r"sha256:[0-9a-f]{64}", (ROOT / "experiments/phase2/docker/Dockerfile.phase2").read_text())[0]
    assert digest in (ROOT / "tests/test_vast.py").read_text()


def test_tests_and_src_import_no_lab_or_archive():
    """CT-14, with a negative control for each kind of violation."""
    tests_py = [p for p in (ROOT / "tests").rglob("*.py") if p.name != Path(__file__).name]
    assert import_leaks(tests_py, src=False) == []
    assert import_leaks(sorted((ROOT / "src/objdel").rglob("*.py")), src=True) == []


def test_import_scanner_reports_each_kind(tmp_path):
    bad = tmp_path / "bad.py"
    bad.write_text('"""doc: archive/x and experiments/phase1/lab are fine here."""\n'
                   "import sys\nimport archive.x as a\nfrom experiments.phase2.lab.B2 import score\n"
                   'sys.path.insert(0, str(ROOT / "archive" / "x"))\nP = "experiments/phase1/lab/t"\n')  # fmt: skip
    assert len(import_leaks([bad], src=False)) == 4
    bad.write_text("def f():\n    import experiments.phase1.x\n")
    assert import_leaks([bad], src=True) == ["bad.py: import experiments.phase1.x"]
    bad.write_text('"""archive/x"""\nimport os\n')
    assert import_leaks([bad], src=False) == []


def test_no_live_file_names_a_retired_tree(tmp_path):
    """CT-15."""
    live = [f for f in tracked() if f.endswith((".py", ".md", ".yaml", ".yml", ".toml", ".txt", "Makefile", ".sh"))
            and not f.startswith("archive/") and f != "tests/test_repo_hygiene_clean.py"]  # fmt: skip
    assert stale_paths(live, ROOT) == []
    (tmp_path / "docs/research").mkdir(parents=True)
    (tmp_path / "a.md").write_text("see experiments/phase0/x")
    (tmp_path / "docs/research/n.md").write_text("see experiments/phase0/x")
    assert stale_paths(["a.md", "docs/research/n.md"], tmp_path) == ["a.md: experiments/phase0"]

"""Static and repository checks of implementation 2.0 (Implementation.md A.9.4) and of the CLEAN phase."""

import ast
import gzip
import json
import re
import subprocess
import sys
import tomllib
from pathlib import Path

import pytest

from objdel.sft import train as sft_train

ROOT = Path(__file__).resolve().parents[2]
NEEDS_ARCHIVE = pytest.mark.skipif(not (ROOT / "archive").exists(), reason="archive/ is local-only (git-ignored)")
LOCAL_DOCS = pytest.mark.skipif(not (ROOT / "docs").is_dir(), reason="docs/ is local-only (git-ignored, 2026-11-11)")
GOLDENS = ROOT / "tests" / "goldens"


def test_r5_test_golden_is_committed_and_holds_only_its_three_keys():
    """S09 (A1 part): the gate CI tests need the R5 test items in the repo, since runs/ is git-ignored."""
    names = ["tests/goldens/r5_test_items.jsonl.gz", "tests/goldens/r5_test_items.meta.json"]
    tracked = subprocess.run(["git", "ls-files", *names], cwd=ROOT, capture_output=True, text=True, check=True)
    assert tracked.stdout.split() == names
    with gzip.open(GOLDENS / "r5_test_items.jsonl.gz", "rt") as f:
        rows = [json.loads(line) for line in f]
    assert len(rows) == 2200 and all(list(r) == ["item_id", "group_id", "iou_tc"] for r in rows)
    assert all(re.fullmatch(r"LIDC-IDRI-\d{4}", r["group_id"]) for r in rows)
    assert set(json.loads((GOLDENS / "r5_test_items.meta.json").read_text())) == {"numpy", "machine"}


CLIS = ("final_eval", "ceiling", "outline_lengths", "gate_power")  # A6 adds nsys_baseline with its module


def test_clis_run_isolated_from_a_planted_module_in_the_cwd(tmp_path):
    """S02: python -I -m <cli> --help never imports a json.py planted in the working directory."""
    import sys

    (tmp_path / "json.py").write_text("raise SystemExit('planted json.py ran')\n")
    for name in CLIS:
        r = subprocess.run([sys.executable, "-I", "-m", f"experiments.phase2.{name}", "--help"], cwd=tmp_path,
                           capture_output=True, text=True)
        assert r.returncode == 0 and "usage" in r.stdout, (name, r.stderr)


@NEEDS_ARCHIVE
def test_the_old_ceiling_is_archived_and_not_imported():
    """S05: rater_ceiling moved to archive/phase2_lab/A4 (replaced by objdel.verify.ceiling); nothing imports it."""
    new = ROOT / "archive/phase2_lab/A4"  # the old tree is gone (test_retired_trees_hold_nothing)
    assert (new / "rater_ceiling.py").is_file() and (new / "rater_ceiling.json").is_file()
    assert "replaced by objdel.verify.ceiling (Phase A2)" in (ROOT / "archive/README.md").read_text()
    for top in ("src", "experiments", "tests"):
        for p in (ROOT / top).rglob("*.py"):
            assert not re.search(r"^\s*(from|import)\s.*rater_ceiling", p.read_text(), re.M), p


def test_metric_definitions_are_in_the_code():
    """S11 (A5): the metrics module docstring states GED's squared form and K, HM-IoU8, and the official boundary
    IoU band; the band widths are BIOU_D."""
    import importlib

    metrics = importlib.import_module("objdel.verify.metrics")
    doc = metrics.__doc__
    assert re.search(r"D(2|²)_GED[^\n]*K = 8", doc) and re.search(r"HM-IoU8[^\n]*not comparable", doc)
    assert re.search(r"[Bb]oundary IoU[^\n]*mask_to_boundary", doc)
    assert metrics.BIOU_D == {"nodule": 2, "disc": 9, "cup": 16}


# --- CLEAN phase ----------------------------------------------------------------------------------------------------
REWARDCORE = {f"rewardcore/{n}.py" for n in ("advantages", "batch", "fixtures", "reward", "topology")}


def tracked(*paths: str) -> list[str]:  # tracked files, as the working tree has them
    out = subprocess.run(["git", "ls-files", "-z", *paths], cwd=ROOT, capture_output=True, text=True, check=True)
    return [f for f in out.stdout.split("\0") if f]


@NEEDS_ARCHIVE
def test_vendor_copies_deleted_and_nothing_names_them():
    """CLEAN3: only the 5 vendor files that differ from their twins may remain (CLEAN5 moved them with the lab)."""
    left = {f.split("/vendor/", 1)[1] for f in tracked("archive/phase1_lab/t_stack/vendor") if "/vendor/" in f}
    assert left <= REWARDCORE
    names = [f for f in tracked("src", "tests", "experiments/phase2") if f != "tests/experiments/test_repo_hygiene.py"
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


@NEEDS_ARCHIVE
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
    for f in ("tests/experiments/test_experiment_scripts.py", "tests/verify/test_eval.py"):
        live_names = {n.name for n in ast.walk(ast.parse((ROOT / f).read_text())) if isinstance(n, ast.FunctionDef)}
        assert not live_names & set(names), f
    wiring = (ROOT / "tests/experiments/test_breakit_wiring.py").read_text()
    assert "census_sizerule" not in wiring


@NEEDS_ARCHIVE
@LOCAL_DOCS
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


@NEEDS_ARCHIVE
@LOCAL_DOCS
def test_phase0_scripts_and_records_moved():
    """CT-87: 8 python files in the archive, records beside the research notes, the oracle parser stays in tests.
    Owner's layout (a255988): README.md, RESULTS.md and RESULTS_0b.md are tracked; the two plans (PLAN.md,
    PLAN_0b.md) stay on the owner's machine, git-ignored, so a worktree may not hold them."""
    assert len(list((ROOT / "archive/phase0").glob("*.py"))) == 8
    rec = ROOT / "docs/research/records/phase0"  # local only since 2026-11-11, like the rest of docs/
    md = {p.name for p in rec.glob("*.md")}
    assert (rec / "results").is_dir() and {"README.md", "RESULTS.md", "RESULTS_0b.md"} <= md
    assert not tracked("docs")
    plans = [f"docs/research/records/phase0/{n}" for n in ("PLAN.md", "PLAN_0b.md")]
    ignored = subprocess.run(["git", "check-ignore", *plans], cwd=ROOT, capture_output=True, text=True).stdout.split()
    assert ignored == plans
    assert "class _Lab" in (ROOT / "tests/geom/test_geometry_oracle.py").read_text()
    for f in tracked("src"):
        for m in re.findall(r"archive/[\w./-]+\.py", (ROOT / f).read_text()):
            assert (ROOT / m).exists(), (f, m)


@NEEDS_ARCHIVE
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
    assert digest in (ROOT / "tests/infra/test_vast.py").read_text()


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
    """CT-15. The files that implement this rule name the trees on purpose: this test and the final sweep with its
    test (the sweep's own SELF set; importing it would make this test import experiments)."""
    checkers = {"tests/experiments/test_repo_hygiene.py", "experiments/phase2/tools/final_sweep.py",
                "tests/experiments/test_final_sweep.py"}
    live = [f for f in tracked() if f.endswith((".py", ".md", ".yaml", ".yml", ".toml", ".txt", "Makefile", ".sh"))
            and not f.startswith("archive/") and f not in checkers]  # fmt: skip
    assert stale_paths(live, ROOT) == []
    (tmp_path / "docs/research").mkdir(parents=True)
    (tmp_path / "a.md").write_text("see experiments/phase0/x")
    (tmp_path / "docs/research/n.md").write_text("see experiments/phase0/x")
    assert stale_paths(["a.md", "docs/research/n.md"], tmp_path) == ["a.md: experiments/phase0"]


# --- CLEAN9: small debts ---------------------------------------------------------------------------------------------


def test_make_preflight_default_config_loads_as_a_run_spec():
    """U19 after E6.8: the default is the live RL spec, and load_spec accepts it."""
    cfg = re.search(r"^CONFIG \?= (\S+)", (ROOT / "Makefile").read_text(), re.M).group(1)
    assert cfg == "specs/b1_rows.yaml"
    from objdel.core.spec import load_spec

    assert load_spec(ROOT / cfg).bundle_view == "views/lidc-bundle-v3"


def test_train_sft_names_its_loader_load_items():
    from experiments.phase2 import train_sft

    assert callable(sft_train.load_items) and not hasattr(train_sft, "preload")


def test_topology_any_docstring_does_not_claim_the_reward_uses_it():
    from objdel.geom import volume

    assert "reward compares" not in volume.__doc__ and "phase G" in volume.__doc__


def test_one_stated_rl_learning_rate():
    from objdel.policy import models

    assert "prereg.G_PROBE" in models.__doc__ and "1e-6 full, 1e-5 LoRA" not in models.__doc__


# --- CLEAN12: dead-code gate and dependency direction ----------------------------------------------------------------


def test_no_dead_code():
    """Fails (no skip) when vulture is missing: a skip would let the gate pass silently."""
    try:
        import vulture  # noqa: F401
    except ImportError:
        pytest.fail("vulture is not installed: install the dev extra (uv pip install -e '.[dev]')")
    r = subprocess.run([sys.executable, "-m", "vulture"], cwd=ROOT, capture_output=True, text=True)
    assert (r.returncode, r.stdout) == (0, "")


def test_vulture_is_pinned_and_configured():
    cfg = tomllib.loads((ROOT / "pyproject.toml").read_text())
    assert "vulture==2.16" in cfg["project"]["optional-dependencies"]["dev"]
    assert cfg["tool"]["vulture"]["paths"] == ["src", "experiments", "tests"]
    assert not (ROOT / "tests" / "vulture_whitelist.py").exists()
    assert re.search(r"^dead:.*\n\t\$\(PY\) -m vulture$", (ROOT / "Makefile").read_text(), re.M)


# --- CLEAN13: the Vast CLI is pinned outside the project environment ------------------------------------------------


def test_vastai_pinned_in_runbook_and_not_a_project_dependency():
    if (ROOT / "docs").is_dir():  # the runbook is local only (docs/, 2026-11-11)
        assert "uv tool install vastai==1.11.0" in (ROOT / "docs/runbooks/gpu-run.md").read_text()
    deps = tomllib.loads((ROOT / "pyproject.toml").read_text())["project"]
    assert not any("vastai" in d for d in deps["dependencies"] + sum(deps["optional-dependencies"].values(), []))


# --- CLEAN14: the library keep-list ADR ------------------------------------------------------------------------------

ADR_DIR = ROOT / "docs/architecture/decisions"
ADR_HEADINGS = ["Context", "Decision", "Options considered", "Trade-off analysis", "Consequences", "Action items"]


def _norm(name: str) -> str:  # PEP 503
    return re.sub(r"[-_.]+", "-", name).lower()


def _adr_sections(text: str) -> dict[str, str]:
    parts = re.split(r"^### (KEEP|DEFERRED|NOT ADOPTED)\s*$", text, flags=re.M)
    return {parts[i]: re.split(r"^#{1,2} ", parts[i + 1], flags=re.M)[0] for i in range(1, len(parts), 2)}


def _first_names(section: str) -> list[str]:  # first backticked token of each bullet
    return [_norm(m) for m in re.findall(r"^- `([^`]+)`", section, flags=re.M)]


@LOCAL_DOCS
def test_library_keep_list_adr_matches_pyproject():
    (path,) = ADR_DIR.glob("0021-*.md")
    text = path.read_text()
    assert text.startswith("# ADR-0021:")
    for h in ["Status", "Date", "Deciders"]:
        assert f"**{h}:**" in text
    for h in ADR_HEADINGS:
        assert re.search(rf"^## {h}$", text, flags=re.M), h
    sec = _adr_sections(text)
    assert set(sec) == {"KEEP", "DEFERRED", "NOT ADOPTED"}
    names = {k: _first_names(v) for k, v in sec.items()}
    flat = sum(names.values(), [])
    assert len(flat) == len(set(flat)) and all(names.values())  # one section per library
    deps = tomllib.loads((ROOT / "pyproject.toml").read_text())["project"]
    every = deps["dependencies"] + sum(deps["optional-dependencies"].values(), [])
    declared = {_norm(re.split(r"[<>=!~\[ ;]", d)[0]) for d in every}
    assert not declared & set(names["NOT ADOPTED"])
    rows = re.split(r"^- ", sec["DEFERRED"], flags=re.M)[1:]
    assert all("Trigger:" in b or re.search(r"phase [A-Z]", b) for b in rows)
    assert path.name in (ADR_DIR / "README.md").read_text()


WAVE = 2  # CLEAN11: D3 (the device `trace`) has landed
SRC = ROOT / "src" / "objdel"


def _functions(path: Path) -> set[str]:
    return {n.name for n in ast.parse(path.read_text()).body if isinstance(n, ast.FunctionDef)}


@pytest.mark.skipif(WAVE < 2, reason="D3 not landed")
def test_polygonizer_seam_is_device_only():
    """CLEAN11 steps 1-4: the numpy tracers live only as oracles; the rings module runs on check_masks; certify
    and RIGA reach it through `trace` and `certify_many`."""
    oracles = ROOT / "tests" / "oracles"
    names = set()
    for o in ("polygonize_oracle", "bitquad_oracle"):
        assert subprocess.run(["git", "ls-files", "--error-unmatch", f"tests/oracles/{o}.py"], cwd=ROOT,
                              capture_output=True).returncode == 0
        names |= _functions(oracles / f"{o}.py")
    for gone in ("checker/polygonize.py", "checker/bitquad.py", "geom/polygonize.py", "geom/bitquad.py"):
        assert not (SRC / gone).exists(), gone
    allowed = {"polygonize", "trace"}
    for f in SRC.rglob("*.py"):
        text = f.read_text()
        assert "tests.oracles" not in text and not re.search(r"^\s*(from|import) oracles", text, re.M), f
        assert not ((_functions(f) & names) - allowed), f
    rings = ast.parse((SRC / "geom" / "rings.py").read_text())
    assert any(getattr(n, "attr", getattr(n, "id", None)) == "check_masks" for n in ast.walk(rings))
    pkg = (SRC / "geom" / "__init__.py").read_text()
    assert "from objdel.geom.rings import polygonize" in pkg
    cert = ast.parse((SRC / "geom" / "certify.py").read_text())
    assert any(isinstance(n, ast.ImportFrom) and any(a.name == "trace" for a in n.names) for n in ast.walk(cert))
    assert any(isinstance(n, ast.Call) and getattr(n.func, "id", None) == "trace" for n in ast.walk(cert))
    assert "certify_many(" in (SRC / "data" / "riga.py").read_text()


def test_retired_models_absent():
    """CLEAN12 and CLEAN8: no retired model name in src/ (case-insensitive)."""
    retired = re.compile("smol" + "vlm|(?<![a-z])l" + "fm|liquid" + "ai", re.I)  # split: the CT-07 scan reads this file
    assert [str(f.relative_to(ROOT)) for f in (ROOT / "src").rglob("*.py") if retired.search(f.read_text())] == []


# --- A.9 S01 and S03, D8: dependencies, markers, test targets ---------------------------------------------------------

PRE_A_DEPENDENCIES = {"numpy", "scipy", "pillow", "opencv-python-headless", "pyarrow", "numba"}
E6_DEPENDENCIES = {"pydantic", "pydantic-settings", "pyyaml"}  # the typed run spec: the only runtime change since A


def test_runtime_dependencies_are_the_pre_a_list_plus_the_run_spec():
    """S01: Phase A adds no dependency; the later additions are the run spec's (E6, plan section 2)."""
    deps = tomllib.loads((ROOT / "pyproject.toml").read_text())["project"]["dependencies"]
    names = {re.split(r"[<>=!~ ]", d)[0] for d in deps}
    assert names == PRE_A_DEPENDENCIES | E6_DEPENDENCIES


def test_every_marker_is_registered_strictly_and_cuda_is_not_one():
    """S03: --strict-markers is on, every marker a test uses is registered, and no test uses `cuda`."""
    cfg = tomllib.loads((ROOT / "pyproject.toml").read_text())["tool"]["pytest"]["ini_options"]
    assert "--strict-markers" in cfg["addopts"]
    registered = {m.split(":")[0] for m in cfg["markers"]}
    builtin = {"parametrize", "skip", "skipif", "xfail", "usefixtures", "filterwarnings"}
    used = {m for f in (ROOT / "tests").rglob("*.py") for m in re.findall(r"pytest\.mark\.(\w+)", f.read_text())}
    assert used - builtin <= registered and "cuda" not in used


def test_makefile_has_the_data_and_order_targets_and_ci_runs_the_order_variants():
    """S03 and D8: `test-data` runs the data tests under -I, `test-order` runs the variants under -I, both are
    phony, and CI runs the variants after pytest. `make test` keeps the data tests (owner, 2026-11-10)."""
    make = (ROOT / "Makefile").read_text()
    assert re.search(r"^test-data:.*\n\t\$\(PY\) -I -m pytest -q -m data$", make, re.M)
    assert re.search(r"^test-order: \$\(EXT_STAMP\).*\n\t\$\(PY\) -I -m tests\.variants$", make, re.M)
    assert all(re.search(rf"^\.PHONY:.* {t}\b", make, re.M) for t in ("test-data", "test-order"))
    assert all("-I" in line for line in make.splitlines() if "pytest" in line and "-m data" in line)
    assert re.search(r"^test:.*\n\t\$\(PY\) -m pytest -q -m \"not slow\" -n auto$", make, re.M)
    steps = [s.strip() for s in (ROOT / ".github/workflows/tests.yml").read_text().split("- run: ")]
    python = [s.split("\n")[0].split("  #")[0] for s in steps]
    assert python.index(".venv/bin/python -I -m tests.variants") > python.index(
        '.venv/bin/python -m pytest -q -m "not slow"')

"""final_sweep: each check on a planted tree (positive and negative control), and the CLI."""

import json
import subprocess
import sys
from pathlib import Path

import pytest

from experiments.phase2.tools import final_sweep as fs

ROOT = Path(__file__).resolve().parents[1]


def make(tmp_path: Path, files: dict[str, str]) -> fs.Ctx:
    """A tree of text files; every file counts as public (no git needed)."""
    for rel, text in files.items():
        p = tmp_path / rel
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(text)
    return fs.Ctx(tmp_path, files=list(files))


def test_retired_trees_flags_live_files_and_ignores_frozen_ones(tmp_path):
    ctx = make(tmp_path, {"README.md": "see experiments/phase2/" + "lab/x",
                          "docs/research/n.md": "experiments/" + "phase0/y",
                          "src/a.py": "x = 1\n", "experiments/" + "phase0/old.py": "x = 1\n"})
    found = fs.retired_trees(ctx)
    assert any(f.startswith("README.md:1:") for f in found)
    assert not any(f.startswith("docs/research/") for f in found)
    assert "experiments/" + "phase0/old.py: file in a retired tree" in found


def test_retired_models_are_split_by_wave(tmp_path):
    ctx = make(tmp_path, {"src/objdel/m.py": "K = 'smol" + "vlm-256m'\n", "tests/test_x.py": "M = 'lfm2.5-vl-450m'\n",
                          "docs/a.md": "Smol" + "VLM is retired\n", "experiments/phase2/g.py": "# LiquidAI card\n",
                          "experiments/phase2/" + "lab/z.py": "smol" + "vlm\n"})
    assert [f.split(":")[0] for f in fs.retired_256m(ctx)] == ["src/objdel/m.py"]  # docs/ is out of scope
    assert sorted(f.split(":")[0] for f in fs.retired_lfm(ctx)) == ["experiments/phase2/g.py", "tests/test_x.py"]
    assert fs.retired_lfm(make(tmp_path / "clean", {"src/objdel/m.py": "x = 1\n"})) == []


def test_controller_check_reads_src_and_phase2_only(tmp_path):
    ctx = make(tmp_path, {"src/objdel/c.py": "class Controller: ...\n", "experiments/phase2/r.py": "alarms_for(x)\n",
                          "tests/t.py": "Controller\n", "docs/x.md": "Controller\n"})
    assert sorted(f.split(":")[0] for f in fs.controller(ctx)) == ["experiments/phase2/r.py", "src/objdel/c.py"]


def test_pool_check_finds_names_flags_and_stage_fields(tmp_path):
    ctx = make(tmp_path, {
        "experiments/phase2/ceiling.py": "pool = serve.fork_pool(2)\n",
        "experiments/phase2/select_ckpt.py": "ap.add_argument('--procs', type=int)\n",
        "experiments/phase2/build_eval_views.py": "ap.add_argument('--procs', type=int)\n",  # keeps its own flag
        "src/objdel/s.py": "cpus = 4\nclass SelectStage:\n    procs: int = 1\n",
        "src/objdel/ok.py": "'''quota and cpus in prose'''\nx = 1\n"})
    found = fs.pool_leftovers(ctx)
    text = "\n".join(found)
    for want in ("ceiling.py:1: pool name", "select_ckpt.py:1: --procs flag", "s.py:1: pool identifier `cpus`",
                 "s.py:3: SelectStage.procs"):
        assert want in text, (want, text)
    assert "build_eval_views" not in text and "ok.py" not in text


def test_doc_links_separates_ok_local_only_and_missing(tmp_path):
    (tmp_path / "plan").mkdir()
    (tmp_path / "plan/secret.md").write_text("internal")  # on disk but not in the public list
    ctx = make(tmp_path, {"docs/index.md": "[a](guide.md#top) [b](../plan/secret.md) [c](nope.md) [d](https://x.org) "
                                            "[e](<my%20file.md>) ![i](img.png)\n",
                          "docs/guide.md": "x", "docs/my file.md": "x"})
    found = fs.doc_links(ctx)
    assert sorted(f.split(": ", 1)[1] for f in found) == ["local-only: ../plan/secret.md", "missing: img.png",
                                                          "missing: nope.md"]
    assert fs.doc_links(make(tmp_path / "ok", {"README.md": "[x](docs/a.md)", "docs/a.md": "x"})) == []


def test_doc_mentions_reports_backticked_paths_that_are_not_public(tmp_path):
    ctx = make(tmp_path, {"docs/a.md": "`src/objdel/real.py:12` `src/objdel/gone.py` `docs/plan/` and `tests/<x>.py`\n",
                          "src/objdel/real.py": "x = 1\n"})
    (tmp_path / "docs/plan").mkdir()
    (tmp_path / "docs/plan/p.md").write_text("x")
    assert sorted(fs.doc_mentions(ctx)) == ["docs/a.md:1: local-only: docs/plan/",
                                            "docs/a.md:1: missing: src/objdel/gone.py"]


def test_unimported_distinguishes_unused_test_only_and_used(tmp_path):
    ctx = make(tmp_path, {
        "src/objdel/__init__.py": "", "src/objdel/used.py": "x = 1\n", "src/objdel/orphan.py": "x = 1\n",
        "src/objdel/surface.py": "x = 1\n", "src/objdel/main.py": "from objdel import used\nfrom . import named\n",
        "src/objdel/named.py": "x = 1\n", "tests/test_s.py": "from objdel.surface import x\n",
        "Makefile": "run:\n\tpython -m experiments.phase2.job\n", "experiments/phase2/job.py": "x = 1\n"})
    found = fs.unimported(ctx)
    text = "\n".join(found)
    assert "src/objdel/orphan.py: nothing imports or names it" in text
    assert "src/objdel/surface.py: only tests use it" in text
    assert "src/objdel/main.py: nothing imports" in text  # nothing imports the importer itself
    for ok in ("used.py", "named.py", "job.py"):
        assert ok not in text, text


def test_docstring_gate_needs_the_pyproject_selection_and_clean_src(tmp_path):
    base = '[tool.ruff.lint]\nselect = ["E"]\n'
    ctx = make(tmp_path, {"pyproject.toml": base, "src/objdel/__init__.py": "'''pkg.'''\n",
                          "src/objdel/m.py": "'''mod.'''\ndef f():\n    return 1\n"})
    found = fs.docstrings(ctx)
    assert any("select lacks" in f for f in found) and any("D103" in f and "m.py:2" in f for f in found)
    ctx = make(tmp_path / "ok", {"pyproject.toml": '[tool.ruff.lint]\nselect = ["E", "D100", "D101", "D102", "D103", '
                                 '"D104", "D106"]\n', "src/objdel/__init__.py": "'''pkg.'''\n",
                                 "src/objdel/m.py": "'''mod.'''\ndef f():\n    '''Doc.'''\n"})
    assert fs.docstrings(ctx) == []


def test_vulture_gate_lists_dead_code(tmp_path):
    ctx = make(tmp_path, {"pyproject.toml": '[tool.vulture]\npaths = ["src"]\nmin_confidence = 60\n',
                          "src/m.py": "def dead_fn():\n    return 1\n"})
    assert any("dead_fn" in f for f in fs.vulture_gate(ctx))
    assert fs.vulture_gate(make(tmp_path / "ok", {"pyproject.toml": '[tool.vulture]\npaths = ["src"]\n',
                                                  "src/m.py": "print(1)\n"})) == []


def test_compare_flags_a_changed_fingerprint_and_leftover_trees():
    before = {"fingerprint": "a" * 64}
    assert fs.compare(before, {"fingerprint": "a" * 64, "lab": 0, "phase0": 0, "claude": 0, "archive": 0}) == []
    found = fs.compare(before, {"fingerprint": "b" * 64, "lab": 3, "archive": 1})
    assert len(found) == 3 and found[0].startswith("fingerprint changed")
    assert fs.compare(before, {"fingerprint_error": "boom"})[0].startswith("no fingerprint now")


def test_run_fails_only_on_a_hard_finding_and_ignores_info(tmp_path, monkeypatch):
    ctx = make(tmp_path, {})
    monkeypatch.setattr(fs, "CHECKS", [("hard_ok", True, lambda c: ["info: fine"]), ("soft", False, lambda c: ["x"]),
                                        ("hard_bad", True, lambda c: ["bad"])])
    res = {r["name"]: r["failed"] for r in fs.run(ctx)}
    assert res == {"hard_ok": False, "soft": False, "hard_bad": True}
    assert [r["name"] for r in fs.run(ctx, only={"soft"})] == ["soft"]
    assert [r["name"] for r in fs.run(ctx, skip={"soft"})] == ["hard_ok", "hard_bad"]


def test_wave_below_3_makes_the_late_checks_advisory(tmp_path, monkeypatch):
    ctx = make(tmp_path, {})
    bad = lambda c: ["bad"]  # noqa: E731
    monkeypatch.setattr(fs, "CHECKS", [("retired_lfm", True, bad), ("controller", True, bad)])
    assert {r["name"]: r["failed"] for r in fs.run(ctx, wave=2)} == {"retired_lfm": False, "controller": True}
    assert {r["name"]: r["failed"] for r in fs.run(ctx)} == {"retired_lfm": True, "controller": True}


def test_main_exit_code_json_and_unknown_check(monkeypatch, capsys):
    monkeypatch.setattr(fs, "CHECKS", [("a", True, lambda c: []), ("b", True, lambda c: ["bad"])])
    monkeypatch.setattr(fs, "Ctx", lambda *a, **k: None)
    assert fs.main(["--only", "a"]) == 0 and "SWEEP PASSED" in capsys.readouterr().out
    assert fs.main(["--only", "b", "--json"]) == 1
    assert json.loads(capsys.readouterr().out)[0]["failed"] is True
    with pytest.raises(SystemExit):
        fs.main(["--only", "nope"])


def test_every_check_is_registered_once_with_a_docstring():
    names = [n for n, _, _ in fs.CHECKS]
    assert len(names) == len(set(names))
    assert all(fn.__doc__ for _, _, fn in fs.CHECKS)


def test_cli_runs_isolated_on_the_real_tree_for_the_cheap_checks():
    """`python -I -m` works from a foreign cwd; the cheap checks run and print their verdict lines."""
    r = subprocess.run([sys.executable, "-I", "-m", "experiments.phase2.tools.final_sweep", "--only",
                        "retired_trees,controller,doc_links", "--limit", "2"], cwd=ROOT, capture_output=True, text=True)
    assert r.returncode in (0, 1), r.stderr
    for name in ("retired_trees", "controller", "doc_links"):
        assert f"] {name}:" in r.stdout, r.stdout
    assert ("SWEEP PASSED" in r.stdout) == (r.returncode == 0)


def test_ignored_but_absent_file_is_local_only_not_missing(tmp_path):
    """A fresh worktree lacks git-ignored files: git's ignore rules, not existence, decide `local-only`."""
    subprocess.run(["git", "init", "-q", str(tmp_path)], check=True)
    (tmp_path / ".gitignore").write_text("/plan/\n")
    (tmp_path / "docs").mkdir()
    (tmp_path / "docs/a.md").write_text("[p](../plan/p.md) `docs/research/notes/n.md` [q](q.md)\n")
    ctx = fs.Ctx(tmp_path, files=["docs/a.md"])
    assert fs.doc_links(ctx) == ["docs/a.md:1: local-only: ../plan/p.md", "docs/a.md:1: missing: q.md"]

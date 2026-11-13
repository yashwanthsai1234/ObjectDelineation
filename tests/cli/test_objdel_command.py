"""The `objdel` command (objdel.cli.main): group and command lookup, help, the box group's Vast commands, and the
re-exec under `python -I`, so a module planted in the working directory never runs."""
import subprocess
import sys
from pathlib import Path

import pytest

from objdel.cli import main as M
from objdel.core.spec import CLI, MODULE


def test_every_stage_module_is_a_command():
    for module in MODULE.values():
        group, name = module.removeprefix(f"{CLI}.").split(".")
        assert M.resolve([group, name, "--x", "1"]) == ["-m", module, "--x", "1"]


def test_a_dash_stands_for_an_underscore_and_box_forwards_to_vast():
    assert M.resolve(["data", "eval-views", "--help"]) == ["-m", f"{CLI}.data.eval_views", "--help"]
    assert M.resolve(["box", "preflight", "--config", "c"]) == ["-m", f"{CLI}.box.preflight", "--config", "c"]
    assert M.resolve(["box", "sync-up", "r1", "--code-only"]) == ["-m", "objdel.infra.vast", "sync_up", "r1",
                                                                  "--code-only"]


@pytest.mark.parametrize(("argv", "code"), [([], 2), (["--help"], 0), (["train"], 2), (["train", "-h"], 0),
                                            (["nope"], 2), (["train", "nope"], 2), (["box", "nope"], 2)])
def test_help_and_unknown_names(argv, code, capsys):
    assert M.main(argv) == code
    out = capsys.readouterr()
    assert (out.err if code else out.out).startswith(("usage: objdel", "objdel: no "))


def test_groups_and_commands_are_found_on_disk():
    assert set(M.groups()) == {"box", "data", "eval", "hub", "train"}
    assert {"rl", "sft", "select", "probe"} <= set(M.commands("train"))
    assert {"preflight", "offers", "launch", "sync-up", "teardown", "sweep-orphans"} <= set(M.commands("box"))
    assert all(M.commands("data").values())  # every command has a one-line summary


def test_summary_is_the_first_paragraph_on_one_line():
    assert M.summary("One\nparagraph.\n\nNext.") == "One paragraph."
    assert M.summary("word " * 40, width=20) == "word word word ..."


def test_console_script_runs_the_module_isolated(tmp_path):
    script = Path(sys.executable).with_name("objdel")
    assert script.is_file(), "reinstall: uv pip install -e . (adds the objdel console script)"
    (tmp_path / "json.py").write_text("raise SystemExit('planted json.py ran')\n")
    r = subprocess.run([str(script), "eval", "ceiling", "--help"], cwd=tmp_path, capture_output=True, text=True,
                       timeout=120)
    assert r.returncode == 0 and r.stdout.startswith("usage:"), r.stderr[-2000:]

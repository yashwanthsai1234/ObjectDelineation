"""Static and repository checks of implementation 2.0 (Implementation.md A.9.4; CLEAN12 extends this file)."""

import gzip
import json
import re
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
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


def test_the_old_ceiling_is_archived_and_not_imported():
    """S05: rater_ceiling moved to archive/phase2_lab/A4 (replaced by objdel.eval.ceiling); nothing imports it."""
    old, new = ROOT / "experiments/phase2/lab/A4", ROOT / "archive/phase2_lab/A4"
    assert (new / "rater_ceiling.py").is_file() and (new / "rater_ceiling.json").is_file()
    assert not (old / "rater_ceiling.py").exists()
    assert "replaced by objdel.eval.ceiling (Phase A2)" in (ROOT / "archive/README.md").read_text()
    for top in ("src", "experiments", "tests"):
        for p in (ROOT / top).rglob("*.py"):
            assert not re.search(r"^\s*(from|import)\s.*rater_ceiling", p.read_text(), re.M), p


def test_metric_definitions_are_in_the_code():
    """S11 (A5): the metrics module docstring states GED's squared form and K, HM-IoU8, and the official boundary
    IoU band; the band widths are BIOU_D."""
    import importlib

    metrics = importlib.import_module("objdel.eval.metrics")  # the package exports the function of that name
    doc = metrics.__doc__
    assert re.search(r"D(2|²)_GED[^\n]*K = 8", doc) and re.search(r"HM-IoU8[^\n]*not comparable", doc)
    assert re.search(r"[Bb]oundary IoU[^\n]*mask_to_boundary", doc)
    assert metrics.BIOU_D == {"nodule": 2, "disc": 9, "cup": 16}

"""CLEAN0: the baseline records and the local safety tag exist (plan CLEAN0 steps 1-6).

The records live in `runs/clean/` (git-ignored), so this test runs only on the owner's machine, where `runs/` exists.
"""

import json
import subprocess
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
CLEAN = REPO / "runs" / "clean"
TAG = "pre-clean-2026-11-08"

pytestmark = pytest.mark.skipif(not (REPO / "runs").is_dir(), reason="no runs/ folder (not the owner's machine)")


def test_before_json_has_the_gate_counts():
    before = json.loads((CLEAN / "before.json").read_text())
    for key in (
        "fingerprint",
        "code_files",
        "claude",
        "archive",
        "tracked",
        "lab",
        "phase0",
        "vulture_lines",
        "checks",
    ):
        assert key in before, key
    assert before["tracked"] > 0 and isinstance(before["checks"], dict)


def test_baseline_files_are_written():
    assert "passed" in (CLEAN / "baseline.txt").read_text()
    for name in ("freeze.txt", "vulture_before.txt"):
        assert (CLEAN / name).stat().st_size > 0, name


def test_safety_tag_is_a_local_lightweight_tag():
    kinds = subprocess.run(
        ["git", "-C", str(REPO), "cat-file", "-t", TAG], capture_output=True, text=True, check=False
    ).stdout.strip()
    assert kinds == "commit"  # a lightweight tag points straight at a commit, an annotated tag would say "tag"

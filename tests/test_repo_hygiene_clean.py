"""CLEAN phase hygiene tests (joined into tests/test_repo_hygiene.py at the merge)."""

import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
REWARDCORE = {f"rewardcore/{n}.py" for n in ("advantages", "batch", "fixtures", "reward", "topology")}


def tracked(*paths: str) -> list[str]:
    out = subprocess.run(["git", "ls-files", "-z", *paths], cwd=ROOT, capture_output=True, text=True, check=True)
    return [f for f in out.stdout.split("\0") if f]


def test_vendor_copies_deleted_and_nothing_names_them():
    """CLEAN3: only the 5 vendor files that differ from their twins may remain (CLEAN5 moves them with the lab)."""
    left = {f.split("/vendor/", 1)[1] for f in tracked("experiments/phase1/lab/t_stack/vendor") if "/vendor/" in f}
    assert left <= REWARDCORE
    names = [f for f in tracked("src", "tests", "experiments/phase2") if f != "tests/test_repo_hygiene_clean.py"
             and "t_stack/vendor" in (ROOT / f).read_text(errors="ignore")]  # fmt: skip
    assert names == []

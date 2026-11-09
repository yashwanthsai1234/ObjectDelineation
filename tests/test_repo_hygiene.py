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

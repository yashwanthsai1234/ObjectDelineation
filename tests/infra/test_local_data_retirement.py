"""CLEAN10 / CT-52: the three superseded local data views were moved (not copied) to data/_retired/2026-11-08/.

data/ is local and git-ignored, so this test runs only on the owner's machine after the move.
"""

import hashlib
import re
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
DATA = REPO / "data"
RETIRED = DATA / "_retired/2026-11-08"
MOVED = ["views/lidc-bundle", "views/lidc-bundle-v2", "eval_views_archive/probe-precision-r1"]
V3_ROWS_SHA256 = "fa62897dc38c28ad971e4671a487e3879d3d3d77a4b7ceae3ae91bc20c1d269c"  # recorded before the move

pytestmark = pytest.mark.skipif(not (RETIRED / "LIST.txt").is_file(), reason="data/ is local: nothing retired here")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_the_views_moved_and_none_is_left_behind():
    for rel in MOVED:
        assert not (DATA / rel).exists(), rel
        assert (RETIRED / Path(rel).name).is_dir(), rel


def test_the_live_views_stay_and_keep_their_hash():
    assert sha256(DATA / "views/lidc-bundle-v3/rows.jsonl") == V3_ROWS_SHA256  # lidc-bundle is a prefix of -v3
    assert (DATA / "eval_views/lidc-val-rotated").is_dir()
    for name in ("shifted", "rescaled", "rotated"):
        assert (DATA / f"eval_views/lidc-val-{name}").is_dir(), name


def test_list_has_one_line_per_folder_with_true_size_and_hash():
    lines = [ln.split("\t") for ln in (RETIRED / "LIST.txt").read_text().splitlines()]
    assert [p for p, _, _ in lines] == [f"data/{m}" for m in MOVED]
    for path, size, digest in lines:
        folder = RETIRED / Path(path).name
        assert int(size) == sum(f.stat().st_size for f in folder.rglob("*") if f.is_file())
        rows = folder / "rows.jsonl"
        assert digest == (sha256(rows) if rows.is_file() else "-")


def test_no_live_spec_or_library_file_names_a_retired_view():
    gone = re.compile(r"lidc-bundle(?!-v3)|probe-precision-r1")
    live = [*REPO.glob("specs/*.yaml"), *REPO.glob("src/**/*.py")]
    assert live and [p.name for p in live if gone.search(p.read_text())] == []

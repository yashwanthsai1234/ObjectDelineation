"""E0: `tests/` has no `__init__.py`, so pytest's default import mode needs unique test basenames."""
from collections import defaultdict
from pathlib import Path

TESTS = Path(__file__).resolve().parents[1]


def test_test_basenames_are_unique():
    seen = defaultdict(list)
    for p in TESTS.rglob("test_*.py"):
        seen[p.name].append(p.relative_to(TESTS).as_posix())
    assert {k: v for k, v in seen.items() if len(v) > 1} == {}

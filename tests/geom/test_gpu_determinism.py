"""Determinism on the GPU (Implementation.md D8): each call of the sanitizer corpus, run 20 times, returns
byte-identical outputs. racecheck cannot see global-memory races, so this test covers them. eval_sets joins with D7."""
import pytest

from tests.geom.test_gpu_sanitize import CORPUS

pytestmark = pytest.mark.gpu


@pytest.mark.parametrize("name", CORPUS)
def test_twenty_runs_are_byte_identical(name):
    first = {k: v.tobytes() for k, v in CORPUS[name]().items()}
    for _ in range(19):
        assert {k: v.tobytes() for k, v in CORPUS[name]().items()} == first

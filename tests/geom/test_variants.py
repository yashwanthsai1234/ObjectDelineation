"""I10 (D-G16): the CPU test variant mutfg4 (K1_MUTATE_FG4, tests/variants.py) labels and traces foreground
4-connected, so its rings and Betti numbers disagree with the true (8, 4) topology on diagonal contacts, and certify,
whose b0, b1 come from a second labeller, must catch it. The mutant must exist when OBJDEL_REQUIRE_VARIANTS=1
(make test-order); it is skipped otherwise. The variant is a non-editable install, so rebuild it after a Python change
(python -I -m tests.variants)."""
import os

import numpy as np
import pytest

from objdel.geom import engine as gc
from tests import variants

DIAG = np.eye(2, dtype=np.uint8)               # two pixels touching at a corner: one 8-connected part
NOTCH = np.ones((5, 5), np.uint8)
NOTCH[0, 0] = NOTCH[1, 1] = 0                  # a hole whose background touches the border at a corner only
CODE = f"""
import numpy as np
from objdel import _gpucheck as g
for m in ({DIAG.tolist()}, {NOTCH.tolist()}):
    m = np.array(m, np.uint8)
    out = g.check_masks(m.ravel(), np.array([m.shape], np.int32), np.array([0, m.size], np.int64))
    print(*out["betti"][0], *out["ring_area"])
"""
CERTIFY = f"""
from objdel.geom.certify import certify
for m in ({DIAG.tolist()}, {NOTCH.tolist()}):
    r = certify(m)
    print(r.certified, r.b0, r.b1, any(f.startswith("rings:") for f in r.failures))
"""


@pytest.fixture
def mutant():
    if not variants.built("mutfg4"):
        if os.environ.get("OBJDEL_REQUIRE_VARIANTS") == "1":
            pytest.fail("variant mutfg4 not built (python -I -m tests.variants)")
        pytest.skip("variant mutfg4 not built (python -I -m tests.variants)")


def _rows(masks):
    rows = []
    for m in masks:
        out = gc.check_masks(m.ravel(), np.array([m.shape], np.int32), np.array([0, m.size], np.int64))
        rows.append([*out["betti"][0].tolist(), *out["ring_area"].tolist()])
    return rows


def test_a_four_adjacent_ring_labeller_breaks_the_topology(mutant):
    assert _rows([DIAG, NOTCH]) == [[1, 0, 2], [1, 1, 24, -1]]  # this build: (8, 4)
    got = [[int(v) for v in line.split()] for line in variants.run_code("mutfg4", CODE, {"K1_MUTATE_FG4"}).splitlines()]
    assert got == [[2, 0, 1, 1], [1, 0, 23]]  # the mutant: two parts, and no hole


def test_certify_catches_a_four_adjacent_ring_labeller(mutant):
    from objdel.geom.certify import certify
    assert [certify(m).certified for m in (DIAG, NOTCH)] == [True, True]  # this build
    got = variants.run_code("mutfg4", CERTIFY, {"K1_MUTATE_FG4"}).split()
    assert got == ["False", "1", "0", "True", "False", "1", "1", "True"]  # the pixel labeller: one part, NOTCH's hole

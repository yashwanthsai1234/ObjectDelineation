import importlib
import json

import numpy as np
import pytest
from scipy.stats import ks_2samp

from objdel.data import coverage as C
from objdel.data.windows import Window


def target(*rings, holes=()):
    return json.dumps({"nodule": [{"outer": r, "holes": list(holes)} for r in rings]})


BIG = [[200, 200], [800, 200], [800, 800], [200, 800]]
SMALL = [[480, 480], [520, 480], [520, 520], [480, 520]]


def test_props_bucket_extent_position_topology_grid():
    p = C.props(target(BIG), 100, "id")
    assert (p["ext"], p["pos"], p["topo"], p["verts"], p["grid10"]) == (">=60%", "centre", "simple", "<8", True)
    q = C.props(target(SMALL, [[900, 900], [950, 900], [950, 950]], holes=[[[1, 1], [2, 2], [3, 1]]]), 50, "r90")
    assert (q["ext"], q["topo"], q["grid10"], q["W"], q["op"]) == ("30-60%", "parts>1", False, 50, "r90")
    edge = C.props(target([[900, 900], [950, 900], [950, 950]]), 100, "id")
    assert (edge["ext"], edge["pos"]) == ("<10%", "edge")
    assert C.props(json.dumps({"nodule": []}), 100, "id") is None


def test_gaps_flag_exam_categories_training_lacks():
    train = C.distribution([C.props(target(BIG), 100, "id")] * 99 + [C.props(target(SMALL), 100, "id")])
    exam = C.distribution([C.props(target(BIG), 100, "id")] * 8 + [C.props(target(SMALL), 100, "id")] * 2)
    assert C.gaps(train, exam) == []  # small is 1 % of training: not below the 1 % floor
    train = C.distribution([C.props(target(BIG), 100, "id")] * 200 + [C.props(target(SMALL), 100, "id")])
    got = C.gaps(train, exam)
    assert [(g["key"], g["category"]) for g in got] == [("ext", "<10%")]
    assert got[0]["exam"] == 0.2 and got[0]["train"] < C.TRAIN_MAX


def test_gaps_ignore_rare_exam_categories():
    train = C.distribution([C.props(target(BIG), 100, "id")] * 10)
    exam = C.distribution([C.props(target(BIG), 100, "id")] * 96 + [C.props(target(SMALL), 100, "id")] * 4)
    assert C.gaps(train, exam) == []  # 4 % < 5 %


# ------------------------------------------------------------------------------- module move (U25), cell_coverage (U18)
def test_coverage_lives_in_the_data_layer_only():
    with pytest.raises(ModuleNotFoundError):
        importlib.import_module("objdel.eval.coverage")
    importlib.import_module("objdel.data.audit")


def _rows(counts: dict, kind="shifted", offset=0.1):
    """Windows whose cells are the keys of counts: size bin b -> extent 5, 15, 30, 50 px; offset bin from offset."""
    ext = (5, 15, 30, 50)
    offs = (0.05, 0.17, 0.25, 0.35)
    return [Window(f"lidc/s/n0/z{i}", "p", "lidc/s/n0", kind, 0, 0, 0, 100, 100,
                   offs[c[1]] if offset is None else offset, ext[c[0]], False)
            for c, n in counts.items() for i in range(n)]


def test_cell_coverage_literals():
    cells = [(a, b) for a in range(4) for b in range(4)]
    uniform = np.full((4, 4), 1 / 16)
    offs = list(np.linspace(0, 0.45, 50))
    full = _rows({c: 10 for c in cells}, offset=None)
    got = C.cell_coverage(full, uniform, offs)
    assert got["worst_ratio"] == pytest.approx(1.0) and got["tv"] == pytest.approx(0.0)
    assert got["ks"] == pytest.approx(ks_2samp([w.offset for w in full], offs).statistic)
    one_gone = _rows({c: 10 for c in cells[1:]}, offset=None)
    got = C.cell_coverage(one_gone, uniform, offs)
    assert got["worst_ratio"] == 0.0 and got["tv"] == pytest.approx(0.0625)
    t = np.zeros((4, 4))
    t[0, 0] = t[1, 1] = 0.5
    rows = _rows({(0, 0): 60, (1, 1): 40}, offset=None)
    got = C.cell_coverage(rows, t, [w.offset for w in rows])
    assert got["worst_ratio"] == pytest.approx(0.8) and got["tv"] == pytest.approx(0.1) and got["ks"] == 0.0
    rare = t * 0.995
    rare[3, 3] = 0.005  # a cell under 1 % of the target with no rows: left out of the worst ratio
    assert C.cell_coverage(rows, rare, [0.1])["worst_ratio"] == pytest.approx(0.4 / 0.4975)
    centred = _rows({(2, 2): 100}, kind="centred")
    assert C.cell_coverage(rows + centred, t, [0.1]) == C.cell_coverage(rows, t, [0.1])


def test_cell_coverage_of_centred_rows_only_fails_without_nan():
    got = C.cell_coverage(_rows({(0, 0): 5}, kind="centred"), np.full((4, 4), 1 / 16), [0.1])
    assert got["worst_ratio"] == 0.0 and got["tv"] == 1.0 and got["ks"] == 1.0


def test_cell_coverage_reads_the_apparent_size(monkeypatch):  # U59
    t = np.zeros((4, 4))
    t[2, 0] = 1.0
    rows = _rows({(2, 0): 10})
    assert C.cell_coverage(rows, t, [0.1])["tv"] == 0.0
    monkeypatch.setattr(Window, "apparent_px", property(lambda w: 5.0))
    assert C.cell_coverage(rows, t, [0.1])["tv"] == 1.0

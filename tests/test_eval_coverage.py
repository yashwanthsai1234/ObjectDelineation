import json

from objdel.eval import coverage as C


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

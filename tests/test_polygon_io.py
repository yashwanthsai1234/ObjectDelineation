from objdel.checker import betti, certify
from objdel.checker.polygon_io import as_points, extract_json, rasterize, ring_problems


def test_polygon_io():
    reply = 'Sure!\n```json\n{"disc": [[10, 10], [90, 10], [90, 90], [10, 90]], "cup": [0, 0, 5]}\n```'
    d = extract_json(reply)
    assert ring_problems(as_points(d["disc"])) == []
    assert as_points(d["cup"]) is None and ring_problems(None) == ["unparseable"]
    assert "self_crossing" in ring_problems(as_points([[0, 0], [10, 10], [10, 0], [0, 10]]))
    assert ring_problems(as_points([[0, 0], [10, 0], [10, 10], [5, -5], [0, 10]])) == ["self_crossing"]
    assert ring_problems(as_points([[0, 0], [1, 1]])) == ["fewer_than_3_points"]
    assert ring_problems(as_points([[0, 0], [10, 0], [10, 10], [0, 10], [0, 0]])) == []
    m = rasterize(as_points(d["disc"]), 100, 100)
    assert certify(m).certified and betti(m) == (1, 0) and m[50, 50] and not m[5, 5]
    assert extract_json("no json here") is None
    stair = as_points([350, 350, 350, 450, 450, 450, 450, 550, 350, 550, 350, 650])
    assert ring_problems(stair) == ["self_touching"]
    assert ring_problems(as_points([[0, 0], [1200, 0], [0, 10]]), bounds=(0, 1000)) == ["out_of_range"]

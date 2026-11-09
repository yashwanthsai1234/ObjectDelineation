"""K2: validate labels == oracle.ring_problems; fine classes == Python reference; targeted cases."""
import collections
import random

import numpy as np
import pytest

import oracle as O
import polycore as pc
from conftest import env_int
from helpers import counts_as_dicts, k2_ring, mutate, ref_validate, star

N_K2 = env_int("POLYCORE_K2_N", 10000)
D = {n: i for i, n in enumerate(pc.DEFECTS)}


def test_k2_10k_rings_labels_match_oracle():
    st = collections.Counter()
    bad = []
    for i in range(N_K2):
        fam, ring = k2_ring(i)
        res = pc.validate([ring])
        mine = pc.ring_labels(res["ring"][0])
        ref = O.ring_problems(ring)
        st[fam] += 1
        for lab in ref:
            st["label:" + lab] += 1
        if not ref:
            st["label:valid"] += 1
        if mine != ref:
            bad.append((i, ring, mine, ref))
        brute = pc.validate([ring], method="brute")
        assert (brute["ring"] == res["ring"]).all() and (brute["total"] == res["total"]).all(), i
        assert res["area2"][0] == O.twice_area(ring) or res["ring"][0][D["fewer_than_3_points"]]
    print(dict(st))
    f = N_K2 / 10000  # label coverage (10k: ~4.5k crossing, 3k touching, 0.3k zero area, 0.7k fewer, 3.6k valid)
    assert st["label:self_crossing"] > 2000 * f and st["label:self_touching"] > 1500 * f
    assert st["label:zero_area"] > 100 * f and st["label:fewer_than_3_points"] > 300 * f and st["label:valid"] > 2000 * f
    assert not bad, bad[:5]


def test_fine_classes_match_reference():
    """Exact per-class pair counts (pinch / overlap / spike / crossing) vs the slow Python reference."""
    for i in range(0, 3000):
        _, ring = k2_ring(i, seed=31)
        ref, _, _ = ref_validate([ring])
        assert counts_as_dicts(pc.validate([ring])) == ref, (i, ring)


def _multi(rng):
    cx, cy, R = rng.randint(250, 750), rng.randint(250, 750), rng.randint(40, 250)
    rings, roles = [], []
    for _ in range(rng.randint(1, 2)):
        o = len(rings)
        ox, oy = cx + rng.randint(-200, 200), cy + rng.randint(-200, 200)
        rings.append(star(rng, rng.randint(3, 16), ox, oy, R))
        roles.append(-1)
        for _ in range(rng.randint(0, 3)):
            rings.append(star(rng, rng.randint(3, 10), ox + rng.randint(-R, R) // 2, oy + rng.randint(-R, R) // 2,
                              R * rng.uniform(.05, .6))[::-1])
            roles.append(o)
    if rng.random() < .3:
        k = rng.randrange(len(rings))
        rings[k] = rings[k][::-1]
    if rng.random() < .3:
        k = rng.randrange(len(rings))
        rings[k] = mutate(rng, rings[k])
    if rng.random() < .2:  # coarse snapping: many touches
        rings = [[(x // 50 * 50, y // 50 * 50) for x, y in r] for r in rings]
    return rings, roles


def test_multi_ring_matches_reference():
    rng = random.Random(5)
    seen = collections.Counter()
    for i in range(1500):
        rings, roles = _multi(rng)
        res = pc.validate(rings, roles)
        ref, pcx, ptx = ref_validate(rings, roles)
        assert counts_as_dicts(res) == ref, (i, rings, roles)
        assert (res["pair_cross"] == pcx).all() and (res["pair_touch"] == ptx).all()
        assert (pc.validate(rings, roles, method="brute")["ring"] == res["ring"]).all()
        for k, v in zip(pc.DEFECTS, res["total"]):
            seen[k] += int(v > 0)
    print(dict(seen))
    for k in ("ring_ring_crossing", "ring_ring_touching", "hole_outside_outer", "nested_holes",
              "nested_outers", "wrong_winding", "self_crossing"):
        assert seen[k] > 10, k


def test_geos_consistency():
    """Anything we call fully valid is GEOS-valid; GEOS-valid single rings have no self defects."""
    shapely = pytest.importorskip("shapely")
    rng = random.Random(11)
    n_valid = 0
    for i in range(1500):
        rings, roles = _multi(rng)
        res = pc.validate(rings, roles)
        if res["total"].sum() == 0:
            n_valid += 1
            outers = [k for k, r in enumerate(roles) if r == -1]
            polys = [shapely.Polygon(rings[o], [rings[h] for h, r in enumerate(roles) if r == o]) for o in outers]
            g = polys[0] if len(polys) == 1 else shapely.MultiPolygon(polys)
            assert shapely.is_valid(g), (rings, roles, shapely.is_valid_reason(g))
        for k, ring in enumerate(rings):
            clean = [p for j, p in enumerate(ring) if j == 0 or p != ring[j - 1]]
            if len(set(clean)) < 3:
                continue
            if shapely.is_valid(shapely.Polygon(ring)):
                row = res["ring"][k]
                assert sum(row[D[c]] for c in ("self_crossing", "vertex_pinch", "collinear_overlap", "spike", "zero_area")) == 0, ring
    assert n_valid > 100


# ------------------------------------------------------------------ targeted cases
def V(rings, roles=None):
    res = pc.validate(rings, roles)
    return [dict((k, int(v)) for k, v in zip(pc.DEFECTS, row) if v) for row in res["ring"]], res


SQ = [(100, 100), (900, 100), (900, 900), (100, 900)]          # clockwise on screen (area > 0)
HOLE = [(400, 400), (400, 600), (600, 600), (600, 400)]        # counter-clockwise


def test_valid_square_and_annulus():
    assert V([SQ])[0] == [{}]
    assert V([SQ, HOLE], [-1, 0])[0] == [{}, {}]


def test_bowtie():
    c, _ = V([[(0, 0), (1000, 1000), (1000, 0), (0, 1000)]])
    assert c[0] == {"self_crossing": 1, "zero_area": 1}
    c, _ = V([[(0, 0), (1000, 1000), (1000, 0), (0, 800)]])
    assert c[0] == {"self_crossing": 1, "wrong_winding": 1}  # net shoelace < 0: the larger lobe is ccw
    c, _ = V([[(0, 800), (1000, 0), (1000, 1000), (0, 0)]])  # reversed: net area > 0
    assert c[0] == {"self_crossing": 1}


def test_figure8_vertex_pinch():
    fig8 = [(0, 0), (10, 0), (5, 5), (10, 10), (0, 10), (5, 5)]
    c, _ = V([fig8])
    assert set(c[0]) == {"vertex_pinch"} and O.ring_problems(fig8) == ["self_touching"]
    # vertex on a non-adjacent edge interior (T-junction)
    tj = [(0, 0), (10, 0), (10, 10), (5, 0), (0, 10)]
    c, _ = V([tj])
    assert c[0].get("vertex_pinch", 0) > 0 and "self_crossing" not in c[0]


def test_spike_and_overlap():
    spike = [(100, 100), (500, 100), (700, 100), (500, 100), (900, 100), (900, 900), (100, 900)]
    c, _ = V([spike])
    assert c[0].get("spike", 0) >= 1
    out_back = [(100, 100), (900, 100), (900, 900), (950, 950), (900, 900), (100, 900)]
    assert V([out_back])[0][0].get("spike") == 1
    # non-adjacent collinear overlap: stair that runs back along itself
    ov = [(0, 0), (6, 0), (6, 4), (3, 4), (3, 0), (9, 0), (9, 8), (0, 8)]  # edges 0 and 4 share [3, 6] x {0}
    c, _ = V([ov])
    assert c[0].get("collinear_overlap", 0) >= 1, c


def test_hole_touching_outer_at_one_point():
    hole = [(100, 100), (300, 400), (400, 300)]
    hole = hole if O.twice_area(hole) < 0 else hole[::-1]
    c, res = V([SQ, hole], [-1, 0])
    # both hole edges at the shared corner touch both square edges there: 4 edge pairs
    assert c[0] == {"ring_ring_touching": 4} and c[1] == {"ring_ring_touching": 4}
    assert res["pair_touch"][0, 1] == 4 and res["total"][D["ring_ring_touching"]] == 4


def test_hole_outside_and_crossing():
    far = [(950, 950), (950, 990), (990, 990), (990, 950)]
    c, _ = V([SQ, far], [-1, 0])
    assert c[1] == {"hole_outside_outer": 1}
    crossing = [(800, 400), (800, 600), (980, 600), (980, 400)]
    c, res = V([SQ, crossing], [-1, 0])
    assert c[1]["ring_ring_crossing"] == 2 and res["pair_cross"][0, 1] == 2
    assert c[1].get("hole_outside_outer") == 1


def test_nested_holes_and_island():
    h1 = [(200, 200), (200, 800), (800, 800), (800, 200)]
    h2 = [(300, 300), (300, 700), (700, 700), (700, 300)]
    c, _ = V([SQ, h1, h2], [-1, 0, 0])
    assert c[1] == {"nested_holes": 1} and c[2] == {"nested_holes": 1}
    island = [(350, 350), (650, 350), (650, 650), (350, 650)]  # outer inside the hole: valid
    c, _ = V([SQ, h1, island], [-1, 0, -1])
    assert c == [{}, {}, {}]
    c, res = V([SQ, island])  # two outers, one inside the other's filled region
    assert c[0] == {"nested_outers": 1} and res["total"][D["nested_outers"]] == 1


def test_wrong_winding_and_out_of_range():
    c, _ = V([SQ[::-1], HOLE[::-1]], [-1, 0])
    assert c[0] == {"wrong_winding": 1} and c[1] == {"wrong_winding": 1}
    c, _ = V([[(-5, 0), (1000, 0), (1001, 1000), (0, 1000)]])
    assert c[0] == {"out_of_range": 2}


def test_degenerate_rings():
    c, _ = V([[], [(5, 5)], [(1, 1), (2, 2)], [(1, 1), (2, 2), (1, 1), (2, 2)]])
    assert c == [{"fewer_than_3_points": 1}] * 4  # degenerate rings take part in no ring-ring test
    assert V([[(0, 0), (5, 0), (10, 0)]])[0] == [{"zero_area": 1, "spike": 2}]
    assert V([[(0, 0), (10, 0), (10, 10), (0, 10), (0, 0)]])[0] == [{}]  # closing duplicate allowed


def test_signed_area2():
    assert pc.signed_area2(SQ) == O.twice_area(SQ) > 0
    assert pc.signed_area2(HOLE) < 0
    big = [(-2**31, -2**31), (2**31 - 1, -2**31), (2**31 - 1, 2**31 - 1)]
    assert pc.signed_area2(big) == O.twice_area(big)
    assert pc.signed_area2([]) == 0


def test_roles_are_checked():
    with pytest.raises(ValueError):
        pc.validate([SQ, HOLE], [-1, 5])
    with pytest.raises(ValueError):
        pc.validate([SQ, HOLE, HOLE], [-1, 0, 1])  # parent is a hole
    with pytest.raises(ValueError):
        pc.validate([SQ, HOLE], [-1])
    np.testing.assert_array_equal(pc.validate([SQ, HOLE], np.array([-1, 0]))["total"], 0)

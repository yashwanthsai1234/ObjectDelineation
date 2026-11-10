import random
import numpy as np
import oracle as O

def test_square_full():
    sq = [(0, 0), (1000, 0), (1000, 1000), (0, 1000)]
    assert O.rasterize([sq], 10, 7).all()
    assert O.rasterize([sq], 10, 7, rule="nonzero").all()

def test_annulus_rules():
    out, inn = O.annulus(hole_ccw=True)
    m_eo = O.rasterize([out, inn], 40, 40); m_nz = O.rasterize([out, inn], 40, 40, rule="nonzero")
    assert not m_eo[20, 20] and not m_nz[20, 20] and m_eo.sum() == m_nz.sum()
    out, inn = O.annulus(hole_ccw=False)  # same orientation: nonzero fills the hole
    assert not O.rasterize([out, inn], 40, 40)[20, 20] and O.rasterize([out, inn], 40, 40, rule="nonzero")[20, 20]

def test_problems():
    assert "self_crossing" in O.ring_problems([(0, 0), (1000, 1000), (1000, 0), (0, 1000)])  # symmetric bow-tie: also zero net area
    assert O.ring_problems([(0, 0), (1000, 1000), (1000, 0), (0, 800)]) == ["self_crossing"]
    stair = [(350,350),(350,450),(450,450),(450,550),(350,550),(350,650)]
    assert "self_touching" in O.ring_problems(stair)
    assert O.ring_problems([(0,0),(10,0),(10,10),(0,10),(0,0)]) == []
    assert O.ring_problems([(0,0),(1,1)]) == ["fewer_than_3_points"]
    assert O.ring_problems([(0,0),(5,0),(10,0)]) == ["fewer_than_3_points"] or "zero_area" in O.ring_problems([(0,0),(5,0),(10,0)])
    # vertex touching a non-adjacent edge (pinch)
    assert O.ring_problems([(0,0),(10,0),(5,5),(10,10),(0,10),(5,5)]) != []

def test_matches_float_crossing_on_generic_points():
    rng = random.Random(0)
    for _ in range(40):
        ring = O.random_polygon_soup(rng, 9)
        W, H = 37, 23
        m = O.rasterize([ring], W, H)
        xs = (np.arange(W) + 0.5) * 1000 / W
        ys = (np.arange(H) + 0.5) * 1000 / H
        ref = np.zeros((H, W), bool)
        for r, py in enumerate(ys):
            for c, px in enumerate(xs):
                inside = False
                for i in range(len(ring)):
                    (x0, y0), (x1, y1) = ring[i], ring[(i + 1) % len(ring)]
                    if (y0 <= py < y1) or (y1 <= py < y0):
                        if x0 + (py - y0) * (x1 - x0) / (y1 - y0) > px:
                            inside = not inside
                ref[r, c] = inside
        assert (ref == m).all()

def test_simple_generator_valid():
    rng = random.Random(1)
    for _ in range(50):
        assert O.ring_problems(O.random_simple_polygon(rng, 10)) == []

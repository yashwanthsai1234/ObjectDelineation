"""The verifier service API (E7), CPU build: one scoring path with the trainer's reward."""
import dataclasses
import json
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import pytest
from conftest import KEY, SIDE, build_root
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

from objdel.geom import engine as _gpucheck
from objdel.geom.text import MAX_CHARS
from objdel.verify.api import API_VERSION, ScoreItem, ScoreResult, Stamp, Verifier, VerifyError
from objdel.verify.reward.advantages import logged_terms
from objdel.verify.reward.cache import ExpertCache
from objdel.verify.reward.score import score_batch
from objdel.verify.reward.spec import REWARD_V4, fingerprint


def test_score_has_the_bits_of_the_shared_scoring_function(view_root, answers):
    v = Verifier("reward_v4", view_root)
    items = [ScoreItem(f"i{i}", KEY, a, False) for i, a in enumerate(answers)]
    items += [ScoreItem("cut", KEY, answers[0], True), ScoreItem("bad", KEY, "not json", False)]
    got = v.score(items)
    cache = ExpertCache.load(view_root / "caches" / KEY / "nodule.npz")
    T = score_batch([(i.text, i.truncated, ("nodule",), {"nodule": cache}, SIDE, SIDE) for i in items])
    terms = logged_terms(T)
    assert [r.id for r in got] == ["i0", "i1", "i2", "cut", "bad"]
    for j, r in enumerate(got[:3]):
        assert r.error is None and r.defects == {}
        assert {k: x.hex() for k, x in r.terms.items()} == {k: float(t[j]).hex() for k, t in terms.items()}
        want = sum(w * terms[k][j] for w, k in zip(REWARD_V4.weights, ("iou_tc", "f5", "legal", "g_area"),
                                                   strict=True))
        assert r.reward.hex() == float(want).hex()
    for r, code in zip(got[3:], T["error"][3:], strict=True):  # a bad outline is a result, not an exception
        assert r.reward is None and r.error == code != ""


def test_unknown_key_refuses_the_whole_batch(view_root, answers):
    v = Verifier("reward_v4", view_root)
    with pytest.raises(VerifyError) as e:
        v.score([ScoreItem("a", KEY, answers[0], False), ScoreItem("b", "lidc/nope", answers[0], False)])
    assert e.value.code == "unknown_key"


def test_unknown_reward_is_a_spec_mismatch(view_root):
    with pytest.raises(VerifyError) as e:
        Verifier("reward_v3", view_root)
    assert e.value.code == "spec_mismatch"


def test_stamp_and_spec_hold_the_fingerprint(view_root):
    v = Verifier("reward_v4", view_root)
    s = v.stamp
    assert (s.api_version, s.reward_name, s.reward_fingerprint, s.checker_build) == (
        API_VERSION, "reward_v4", fingerprint(), _gpucheck.build())
    assert s.code_version
    assert v.spec()["stamp"] == dataclasses.asdict(s) and v.spec()["reward"]["group"] == REWARD_V4.group


# --- E-G12 properties of the verifier (P4, P5, P6, U32, I25, I26, U45) -----------------------------------------------

def box(x0, y0, x1, y1, holes=()):
    """Outline text of one nodule part: the clockwise grid box (x0, y0)-(x1, y1) with hole rings `holes`."""
    ring = [[x0, y0], [x1, y0], [x1, y1], [x0, y1]]
    return json.dumps({"nodule": [{"outer": ring, "holes": [list(h) for h in holes]}]}, separators=(",", ":"))


def square(x0, y0, x1, y1, side=SIDE):
    m = np.zeros((side, side), bool)
    m[y0:y1, x0:x1] = True
    return m


def bits(r):
    """A ScoreResult with every float as hex: equal bits, not equal values."""
    return (r.id, None if r.reward is None else r.reward.hex(), {k: v.hex() for k, v in r.terms.items()},
            dict(r.defects), r.error)


BROKEN = ["", "null", "[]", "{", "not json", "[[1,2]]", '{"nodule":[{"outer":[[1,1],[2,1],[2,99999999999999999999]],'
          '"holes":[]}]}', '{"nodule":[{"outer":[[1.5,1],[2000.0,1],[2,2]],"holes":[]}]}',
          '{"nodule":[[[[1,1],[2,1],[2,2]]]]}', "\x00\xff\x9c garbage", "\ud800", "[" * 5_000_000,
          "[" * 100_000 + "]" * 100_000]


def corpus(answers):
    rng = np.random.default_rng(0)
    boxes = [box(*(int(v) for v in (*xy, *(xy + wh)))) for xy, wh in zip(rng.integers(0, 600, (12, 2)),
                                                                          rng.integers(20, 400, (12, 2)), strict=True)]
    return [*answers, *boxes, *BROKEN[:11]]


@settings(max_examples=25, deadline=None, suppress_health_check=[HealthCheck.function_scoped_fixture])
@given(st.data())
def test_score_of_a_batch_equals_score_of_each_item_alone(view_root, answers, data):
    """P4: duplicates allowed; a permutation permutes the results; an empty batch is []."""
    texts = corpus(answers)
    v = Verifier("reward_v4", view_root)
    picks = data.draw(st.lists(st.integers(0, len(texts) - 1), min_size=1, max_size=20))
    items = [ScoreItem(f"i{j}", KEY, texts[p], p % 5 == 4) for j, p in enumerate(picks)]
    got = [bits(r) for r in v.score(items)]
    assert got == [bits(v.score([i])[0]) for i in items]
    perm = data.draw(st.permutations(range(len(items))))
    assert [bits(r) for r in v.score([items[k] for k in perm])] == [got[k] for k in perm]
    assert v.score([]) == []


@pytest.fixture
def square_roots(tmp_path, clean_store):
    """Two roots whose three raters all drew the square (10, 10)-(30, 30); the second shifted by (7, 4) pixels."""
    return (build_root(tmp_path / "a", {KEY: [square(10, 10, 30, 30)] * 3}),
            build_root(tmp_path / "b", {KEY: [square(17, 14, 37, 34)] * 3}))


def test_verifier_math_symmetries(square_roots):
    """P5. One pixel is 20 grid units (SIDE 50 on the 0-1000 grid)."""
    a, b = (Verifier("reward_v4", r) for r in square_roots)
    same = a.score([ScoreItem("s", KEY, box(200, 200, 600, 600), False)])[0]
    assert same.terms["iou_tc"] == 1.0 and same.reward == 1.0 and same.defects == {}
    # shifting polygon and reference together by whole pixels changes nothing
    for x0, y0, x1, y1 in [(200, 200, 600, 600), (240, 160, 500, 640), (100, 300, 400, 420)]:
        pa = a.score([ScoreItem("p", KEY, box(x0, y0, x1, y1), False)])[0]
        pb = b.score([ScoreItem("p", KEY, box(x0 + 140, y0 + 80, x1 + 140, y1 + 80), False)])[0]
        assert bits(pa) == bits(pb)
    # the reward stays inside REWARD_V4's range [0, sum of weights] over 1,000 random boxes
    rng = np.random.default_rng(1)
    xy, wh = rng.integers(0, 900, (1000, 2)), rng.integers(1, 100, (1000, 2)) * 20
    got = a.score([ScoreItem(str(i), KEY, box(*map(int, (*p, *np.minimum(p + s, 1000)))), False)
                   for i, (p, s) in enumerate(zip(xy, wh, strict=True))])
    assert all(r.error is None and 0.0 <= r.reward <= sum(REWARD_V4.weights) for r in got)
    # shrinking a box that covers the reference never raises IoU
    shrink = a.score([ScoreItem(str(k), KEY, box(200 + 20 * k, 200 + 20 * k, 600 - 20 * k, 600 - 20 * k), False)
                      for k in range(9)])
    iou = [r.terms["iou_tc"] for r in shrink]
    assert iou[0] == 1.0 and iou == sorted(iou, reverse=True) and iou[-1] < iou[0]
    # swapping hole and outer gives a defect
    outer, hole = (200, 200, 600, 600), [[300, 300], [300, 500], [500, 500], [500, 300]]
    ok = a.score([ScoreItem("h", KEY, box(*outer, holes=[hole]), False)])[0]
    swapped = json.dumps({"nodule": [{"outer": hole[::-1], "holes": [[[200, 200], [200, 600], [600, 600],
                                                                      [600, 200]]]}]}, separators=(",", ":"))
    bad = a.score([ScoreItem("w", KEY, swapped, False)])[0]
    assert ok.defects == {} and bad.defects and bad.terms["legal"] == 0.0


def test_verifier_is_thread_safe(view_root, answers):
    """P6: 8 threads, 200 concurrent calls on different batches; each equals its single-thread result bitwise."""
    texts = corpus(answers)
    v = Verifier("reward_v4", view_root)
    rng = np.random.default_rng(2)
    batches = [[ScoreItem(f"b{b}i{j}", KEY, texts[p], False) for j, p in enumerate(rng.integers(0, len(texts), n))]
               for b, n in enumerate(rng.integers(1, 12, 200))]
    alone = [[bits(r) for r in v.score(b)] for b in batches]
    with ThreadPoolExecutor(8) as ex:
        together = list(ex.map(lambda b: [bits(r) for r in v.score(b)], batches))
    assert together == alone


def test_a_text_that_does_not_parse_is_a_result_not_an_exception(view_root, answers):
    """U32: each broken text is a ScoreResult with an error code, reward None and terms {}; the neighbour in the same
    batch scores as it does alone."""
    v = Verifier("reward_v4", view_root)
    good = ScoreItem("good", KEY, answers[0], False)
    alone = bits(v.score([good])[0])
    for text in BROKEN:
        bad, neighbour = v.score([ScoreItem("bad", KEY, text, False), good])
        assert bad.error and bad.reward is None and bad.terms == {} and bad.defects == {}, text[:40]
        assert bits(neighbour) == alone


@pytest.fixture
def big_root(tmp_path, clean_store):
    """A 512 x 512 image whose raters drew a centred square."""
    return build_root(tmp_path, {KEY: [square(128, 128, 384, 384, side=512)] * 3})


def test_verifier_borders_and_sizes(big_root):
    """I25: the whole image, 10,000 vertices, MAX_CHARS +/- 1 character, coordinates 0 and 1000 (valid) and 1001
    (the parser's "range" code); no crash."""
    v = Verifier("reward_v4", big_root)
    whole = box(0, 0, 1000, 1000)
    ring = [[x, 0] for x in range(1000)] + [[1000, y] for y in range(1000)]
    ring += [[x, 1000] for x in range(1000, 0, -1)] + [[0, y] for y in range(1000, 0, -1)]
    many = json.dumps({"nodule": [{"outer": (ring * 3)[:10_000], "holes": []}]}, separators=(",", ":"))
    pad = MAX_CHARS - len(whole)
    texts = {"whole": whole, "many": many, "under_cap": whole + " " * (pad - 1), "at_cap": whole + " " * pad,
             "over_cap": whole + " " * (pad + 1), "1001": box(0, 0, 1001, 1000)}
    got = {r.id: r for r in v.score([ScoreItem(k, KEY, t, False) for k, t in texts.items()])}
    assert got["whole"].error is None and got["whole"].defects == {} and got["whole"].terms["iou_tc"] == 0.25
    assert bits(got["at_cap"])[1:] == bits(got["under_cap"])[1:] == bits(got["whole"])[1:]
    assert (got["over_cap"].error, got["many"].error, got["1001"].error) == ("too_long", "too_long", "range")


def test_verify_api_is_stdlib_only():
    """I26: the library API loads neither pydantic, torch nor fastapi; its dataclasses are frozen."""
    code = "import sys, objdel.verify.api\nprint(sorted({'pydantic', 'torch', 'fastapi'} & set(sys.modules)))"
    out = subprocess.run([sys.executable, "-I", "-c", code], capture_output=True, text=True, check=True).stdout
    assert out.strip() == "[]"
    for obj, field in ((ScoreItem("a", KEY, "", False), "text"), (ScoreResult("a", 1.0, {}, {}, None), "reward"),
                       (Stamp("1.0", "reward_v4", "f", "cpu", "0"), "checker_build")):
        with pytest.raises(dataclasses.FrozenInstanceError):
            setattr(obj, field, None)


def test_a_device_error_becomes_a_verify_error(view_root, answers, monkeypatch):
    """U45: with the engine's checker raising RuntimeError("CUDA error"), score raises VerifyError("device_error"),
    never a bare RuntimeError; the next call works once the engine recovers."""
    v = Verifier("reward_v4", view_root)
    items = [ScoreItem("a", KEY, answers[0], False)]
    want = [bits(r) for r in v.score(items)]
    real = _gpucheck.score_polygons

    def broken(*a, **k):
        raise RuntimeError("CUDA error")

    monkeypatch.setattr(_gpucheck, "score_polygons", broken)
    with pytest.raises(VerifyError) as e:
        v.score(items)
    assert e.value.code == "device_error" and "CUDA error" in str(e.value)
    monkeypatch.setattr(_gpucheck, "score_polygons", real)
    assert [bits(r) for r in v.score(items)] == want

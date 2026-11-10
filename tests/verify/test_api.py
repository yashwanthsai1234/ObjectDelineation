"""The verifier service API (E7), CPU build: one scoring path with the trainer's reward."""
import dataclasses

import pytest
from conftest import KEY, SIDE

from objdel import _gpucheck
from objdel.reward import REWARD_V4, ExpertCache, fingerprint
from objdel.train import rewards
from objdel.verify.api import API_VERSION, ScoreItem, Verifier, VerifyError


def test_score_has_the_bits_of_the_shared_scoring_function(view_root, answers):
    v = Verifier("reward_v4", view_root)
    items = [ScoreItem(f"i{i}", KEY, a, False) for i, a in enumerate(answers)]
    items += [ScoreItem("cut", KEY, answers[0], True), ScoreItem("bad", KEY, "not json", False)]
    got = v.score(items)
    cache = ExpertCache.load(view_root / "caches" / KEY / "nodule.npz")
    T = rewards.score_jobs([(i.text, i.truncated, ("nodule",), {"nodule": cache}, SIDE, SIDE) for i in items])
    terms = rewards.logged_terms(T)
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

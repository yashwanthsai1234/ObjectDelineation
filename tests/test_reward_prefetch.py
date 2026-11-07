"""TRLReward.prefetch: rewards start scoring (parse, one checker call, pool) as soon as vLLM returns the answers, and
overlap TRL's old-logprob and reference forwards on the GPU; the reward call takes the result only for the same
answers (objdel.train.fast.FastGRPO._generate)."""

import threading

from objdel import _gpucheck
from objdel.train import rewards as R
from tests.test_reward_one_checker_call import CACHES, EOS, EXPECTED, PAD, batch


def reward():
    rw = R.TRLReward(CACHES.__getitem__, stop_ids={EOS, PAD})
    return rw, rw.funcs()[0]


def test_prefetched_rewards_equal_the_direct_path(monkeypatch):
    calls = []
    real = _gpucheck.check_polygons
    monkeypatch.setattr(_gpucheck, "check_polygons", lambda *a: calls.append(threading.current_thread()) or real(*a))
    rw, funcs = reward()
    kw = batch()
    rw.prefetch(kw["completions"], kw["completion_ids"], kw["key"], kw["W"], kw["H"], kw["names"])
    assert {f.__name__: f(**kw) for f in funcs} == EXPECTED
    assert len(calls) == 1 and calls[0] is not threading.main_thread()  # scored once, off the main thread


def test_a_prefetch_of_other_answers_is_not_used(monkeypatch):
    rw, funcs = reward()
    other = batch(truncate=())  # same texts, other truncation: other rewards
    rw.prefetch(other["completions"], other["completion_ids"], other["key"], other["W"], other["H"], other["names"])
    assert {f.__name__: f(**batch()) for f in funcs} == EXPECTED

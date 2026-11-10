"""nvtx.count_syncs (plan C0 step 1): torch's sync-debug warnings counted per block, every occurrence, mode restored."""

import warnings

import pytest

torch = pytest.importorskip("torch")
from objdel.policy import marks as nvtx  # noqa: E402


@pytest.fixture
def fake_cuda(monkeypatch):
    mode = {"now": 0, "set": []}
    monkeypatch.setattr(torch.cuda, "is_available", lambda: True)
    monkeypatch.setattr(torch.cuda, "get_sync_debug_mode", lambda: mode["now"])

    def set_mode(m):
        mode["set"].append(m)
        mode["now"] = m

    monkeypatch.setattr(torch.cuda, "set_sync_debug_mode", set_mode)
    return mode


def sync():
    warnings.warn("called a synchronizing CUDA operation", UserWarning, stacklevel=1)


def test_every_sync_of_a_loop_counts_and_other_warnings_do_not(fake_cuda):
    with nvtx.count_syncs() as n:
        for _ in range(1000):
            sync()  # one code line: Python's default filter would record it once
        warnings.warn("something else", UserWarning, stacklevel=1)
    assert n.syncs == 1000 and fake_cuda["set"] == ["warn", 0] and fake_cuda["now"] == 0


def test_the_previous_mode_comes_back_after_an_exception(fake_cuda):
    fake_cuda["now"] = 1
    with pytest.raises(KeyError), nvtx.count_syncs():
        raise KeyError
    assert fake_cuda["now"] == 1


def test_without_cuda_the_count_is_zero_and_nothing_is_set(monkeypatch):
    monkeypatch.setattr(torch.cuda, "is_available", lambda: False)
    with nvtx.count_syncs() as n:
        sync()
    assert n.syncs == 0


def test_fast_grpo_logs_host_syncs_per_phase_and_the_vllm_range_apart(fake_cuda, monkeypatch):
    from collections import defaultdict

    pytest.importorskip("trl")
    from objdel.rl import trl_adapter as fast

    monkeypatch.setattr(torch.cuda, "synchronize", lambda: None)
    t = object.__new__(fast.FastGRPO)
    t.count_syncs, t._metrics = True, {"train": defaultdict(list)}
    t.model = torch.nn.Linear(1, 1).train()
    with t.timed("batch"):
        sync()
        with t.timed("generate"):
            sync(), sync()
    m = t._metrics["train"]
    assert m["perf/host_syncs/batch"] == [1] and m["perf/host_syncs_vllm"] == [2] and "time/batch_s" in m

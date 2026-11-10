"""C5, C6, C9 flags: train_rl --fast-sync / --group-requests / --no-ckpt-layers (all off by default: the B1 control),
their make_trainer wiring, and callbacks.NoCkptLayers."""

from types import SimpleNamespace

import pytest

torch = pytest.importorskip("torch")
from test_nvtx import built  # noqa: E402, F401 (fixture)

BASE = ["--model", "m", "--ckpt", "c", "--train", "t", "--out", "o", "--steps", "20"]


def test_train_rl_flags_default_off():
    from experiments.phase2 import train_rl

    a = train_rl.parse(BASE)
    assert (a.fast_sync, a.group_requests, a.no_ckpt_layers, a.lora_sync) == (False, False, 0, None)
    a = train_rl.parse([*BASE, "--fast-sync", "--group-requests", "--no-ckpt-layers", "4"])
    assert (a.fast_sync, a.group_requests, a.no_ckpt_layers) == (True, True, 4)


def test_no_ckpt_layers_is_a_range_checked_count_with_no_clamp(capsys):
    """U-CF1 / U-CO3: 0..28 (Qwen3-VL-2B's text decoder layers); -1 and 29 exit 2."""
    from experiments.phase2 import train_rl

    assert train_rl.parse([*BASE, "--no-ckpt-layers", "28"]).no_ckpt_layers == 28
    for bad in ("-1", "29"):
        with pytest.raises(SystemExit) as e:
            train_rl.parse([*BASE, "--no-ckpt-layers", bad])
        assert e.value.code == 2
        assert "--no-ckpt-layers" in capsys.readouterr().err


def test_make_trainer_adds_no_ckpt_layers_only_when_asked(built):  # noqa: F811
    from objdel.train.callbacks import NoCkptLayers

    assert built() == []
    (cb,) = built(no_ckpt_layers=8)
    assert isinstance(cb, NoCkptLayers) and cb.k == 8


def test_make_trainer_fast_sync_refuses_a_trainer_without_colocate_vllm(built):  # noqa: F811
    built(fast_sync=False)
    with pytest.raises(ValueError, match="--fast-sync"):
        built(fast_sync=True)  # the stub trainer has no vllm_generation


def test_make_trainer_lora_sync_refuses_a_trainer_without_colocate_vllm(built, monkeypatch):  # noqa: F811
    from pathlib import Path

    import trl.generation.vllm_generation as vgen

    monkeypatch.setattr(vgen, "LLM", lambda **kw: None, raising=False)  # vLLM is not installed on the Mac

    assert train_rl_lora_dir("--lora-sync", "d") == Path("d")
    with pytest.raises(ValueError, match="--lora-sync"):
        built(lora_sync=Path("d"))


def train_rl_lora_dir(*extra):
    from experiments.phase2 import train_rl

    return train_rl.parse([*BASE, *extra]).lora_sync


def test_no_ckpt_layers_turns_off_the_last_k_after_trainer_enables_all():
    from test_train_fast import tiny_model

    from objdel.train.callbacks import NoCkptLayers

    m = tiny_model()
    m.gradient_checkpointing_enable()  # Trainer.train does this just before on_train_begin
    layers = m.model.language_model.layers
    NoCkptLayers(1).on_train_begin(None, None, None, model=m)
    assert [x.gradient_checkpointing for x in layers] == [True] * (len(layers) - 1) + [False]
    with pytest.raises(ValueError, match="--no-ckpt-layers"):
        NoCkptLayers(len(layers) + 1).on_train_begin(None, None, None, model=SimpleNamespace(
            named_modules=m.named_modules))  # fmt: skip


def test_no_ckpt_layers_zero_changes_nothing():
    from test_train_fast import tiny_model

    from objdel.train.callbacks import NoCkptLayers

    m = tiny_model()
    m.gradient_checkpointing_enable()
    NoCkptLayers(0).on_train_begin(None, None, None, model=m)
    assert all(x.gradient_checkpointing for x in m.model.language_model.layers)


def test_no_ckpt_layers_changes_exactly_k_layers_and_neither_logprobs_nor_grads():
    """U-CO3 (the tiny model has 3 layers, so K = 2 stands in for 4): exactly K layers keep their activations, and the
    logprobs and gradients are torch.equal to K = 0 (no dropout: the recompute repeats the forward bit for bit)."""
    from test_train_fast import tiny_model

    from objdel.train.callbacks import NoCkptLayers

    ids = torch.randint(0, 99, (2, 12), generator=torch.Generator().manual_seed(3))
    out = {}
    for k in (0, 2):
        m = tiny_model().train()
        m.gradient_checkpointing_enable()
        NoCkptLayers(k).on_train_begin(None, None, None, model=m)
        assert sum(not x.gradient_checkpointing for x in m.model.language_model.layers) == k
        lp = m(input_ids=ids).logits.float().log_softmax(-1).gather(-1, ids[..., None]).squeeze(-1)
        lp.sum().backward()
        out[k] = (lp.detach(), {n: p.grad.clone() for n, p in m.named_parameters() if p.grad is not None})
    assert torch.equal(out[0][0], out[2][0]) and out[0][1].keys() == out[2][1].keys() and out[0][1]
    assert all(torch.equal(g, out[2][1][n]) for n, g in out[0][1].items())

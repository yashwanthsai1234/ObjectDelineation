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
    assert (a.fast_sync, a.group_requests, a.no_ckpt_layers) == (False, False, None)
    a = train_rl.parse([*BASE, "--fast-sync", "--group-requests", "--no-ckpt-layers", "4"])
    assert (a.fast_sync, a.group_requests, a.no_ckpt_layers) == (True, True, 4)


def test_make_trainer_adds_no_ckpt_layers_only_when_asked(built):  # noqa: F811
    from objdel.train.callbacks import NoCkptLayers

    assert built() == []
    (cb,) = built(no_ckpt_layers=8)
    assert isinstance(cb, NoCkptLayers) and cb.k == 8


def test_make_trainer_fast_sync_refuses_a_trainer_without_colocate_vllm(built):  # noqa: F811
    built(fast_sync=False)
    with pytest.raises(ValueError, match="--fast-sync"):
        built(fast_sync=True)  # the stub trainer has no vllm_generation


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

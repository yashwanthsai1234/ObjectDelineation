"""objdel.train.sample.sample: K answers per row with their completion token counts (U46)."""

from types import SimpleNamespace

import pytest

pytest.importorskip("trl")
from objdel.train import sample  # noqa: E402

EOS, PAD = 1, 0


def fake_trainer(lengths):
    """A trainer whose `_generate` returns completion ids in TRL's tuple order: answer j of a call has
    lengths[j % len(lengths)] ids, ending in EOS when that length is odd (so truncation varies)."""
    calls = []

    def _generate(prompts):
        calls.append(len(prompts))
        ids = []
        for j in range(len(prompts)):
            n = lengths[j % len(lengths)]
            ids.append([5] * (n - 1) + [EOS if n % 2 else 7] if n else [])
        return None, ids, None, [f"{p}#{j}" for j, p in enumerate(prompts)]

    model = SimpleNamespace(training=True, eval=lambda: None, train=lambda mode: None,
                            generation_config=SimpleNamespace(eos_token_id=EOS))
    tok = SimpleNamespace(eos_token_id=EOS, pad_token_id=PAD)
    return SimpleNamespace(model=model, processing_class=tok, _generate=_generate), calls


def test_tokens_are_the_completion_id_lengths():
    lengths = [5, 1536, 7, 0, 1, 2, 9, 4]
    trainer, calls = fake_trainer(lengths)
    rows = [{"prompt": f"p{i}", "key": f"k{i}"} for i in range(40)]
    out = sample.sample(trainer, rows, 8, items=16)
    assert calls == [128, 128, 64]
    assert [r["item_id"] for r in out] == [f"k{i}" for i in range(40)]
    for i, r in enumerate(out):
        assert r["tokens"] == lengths  # ids, not characters of the decoded text
        assert len(r["samples"]) == len(r["truncated"]) == len(r["tokens"]) == 8
        assert r["samples"][0] == f"p{i}#{(i % 16) * 8}"  # answers stay in item and answer order
        # the existing rule, frozen: truncated unless the last id is EOS or pad (an empty answer has no end)
        assert r["truncated"] == [False, True, False, True, False, True, False, True]
    assert all(type(t) is int for r in out for t in r["tokens"])


def test_no_rows_no_generate_call():
    trainer, calls = fake_trainer([3])
    assert sample.sample(trainer, [], 8) == [] and calls == []

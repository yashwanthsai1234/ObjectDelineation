"""Stubs for the data-quality CLI tests (audit_data, rl_rows): the ceiling, the tokenizer and the val half."""
import json
import sys
import types

from objdel.core.hashing import file_sha256
from objdel.data import quality as Q


class Tok:
    """A tokenizer stub: the first text is n tokens long, every other text 3."""

    def __init__(self, n):
        self.n = n

    def __call__(self, texts, add_special_tokens=False):
        assert not add_special_tokens
        return {"input_ids": [[0] * (self.n if i == 0 else 3) for i, _ in enumerate(texts)]}


def stub_ceiling(monkeypatch, values):
    """objdel.verify.ceiling with a loo_caches whose left-out doctor r scores values[r % len(values)] (score_mask
    patched to return the stub cache's value)."""
    from objdel.verify.reward import score

    def loo_caches(cache):
        return [types.SimpleNamespace(index=r, cache=values[r % len(values)]) for r in range(cache.k)]
    monkeypatch.setitem(sys.modules, "objdel.verify.ceiling", types.SimpleNamespace(loo_caches=loo_caches))
    monkeypatch.setattr(score, "score_mask", lambda mask, cache: types.SimpleNamespace(iou_tc=cache))


def quality_stubs(world, tmp_path, monkeypatch):
    """Stub the ceiling, the tokenizer (Tok(50)) and splits.val_half, and write the three inputs both CLIs read next
    to `tmp_path`: train_exclusions.json, frozen_views.json (the world's views) and z.json (the world's store).
    Returns SimpleNamespace(ex, frozen, z) of those paths."""
    from objdel.data import splits
    from objdel.policy import models

    stub_ceiling(monkeypatch, [0.9, 0.8])
    monkeypatch.setattr(models, "tokenizer", lambda model: (Tok(50), 1536))
    monkeypatch.setattr(splits, "val_half", lambda g: "shape" if int(g[-1]) % 2 else "select")  # fixture patients
    ex = tmp_path / "train_exclusions.json"
    ex.write_text("{}")
    frozen = tmp_path / "frozen_views.json"
    frozen.write_text(json.dumps({v: {"content_sha256": Q.view_content_sha256(world.views / v),
                                      "target_set_sha256": Q.target_set_sha256(world.views / v)} for v in Q.FROZEN}))
    z = tmp_path / "z.json"
    z.write_text(json.dumps({"train_sha256": file_sha256(world.store / "lidc/train.parquet"), "series": {},
                             "step2": []}))
    return types.SimpleNamespace(ex=ex, frozen=frozen, z=z)

"""B2.2 fp16 Mac check (U53): a stub model with a vision block, two decoder layers and logits."""
import torch
from torch import nn

from tools.studies import fp16_check


class Block(nn.Module):
    def __init__(self, scale):
        super().__init__()
        self.scale = scale

    def forward(self, x):
        return x * self.scale


class Stub(nn.Module):
    def __init__(self, vision=1.0, layer=1.0, logits=1.0):
        super().__init__()
        self.visual = nn.Module()
        self.visual.blocks = nn.ModuleList([Block(vision)])
        self.language_model = nn.Module()
        self.language_model.layers = nn.ModuleList([Block(1.0), Block(layer)])
        self.logits_scale = logits

    def forward(self, x):
        x = self.visual.blocks[0](x)
        for blk in self.language_model.layers:
            x = blk(x)
        return x * self.logits_scale


def run(model):
    return fp16_check.check(model, [(torch.full((1, 4), 2.0),)], lambda m, b: m(*b))


def test_clean_model_passes_and_reports_one_maximum_per_block():
    r = run(Stub())
    assert r["passed"] and r["nonfinite_logits"] == 0
    assert set(r["max_abs"]) == {"visual.blocks.0", "language_model.layers.0", "language_model.layers.1"}
    assert r["max_abs"]["language_model.layers.1"] == 2.0


def test_a_layer_above_half_the_fp16_range_fails():
    r = run(Stub(layer=20000.0))
    assert not r["passed"] and r["over"] == ["language_model.layers.1"]


def test_inf_in_the_logits_fails_and_is_counted():
    r = run(Stub(logits=float("inf")))
    assert not r["passed"] and r["nonfinite_logits"] == 4


def test_main_exit_code_follows_passed():
    assert fp16_check.exit_code({"passed": True}) == 0 and fp16_check.exit_code({"passed": False}) == 1

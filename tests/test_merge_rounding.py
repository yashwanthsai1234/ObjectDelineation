"""U64 (B2.0): the merge-rounding table on a planted case. Hand values: W = 1 in bf16 has ulp 2^-7, so 1 + 2^-10
rounds back to 1 (erased) and 1 + 2^-6 is exact; the error vector is (-2^-10, 0), so the relative error is
2^-10 / sqrt(2^-20 + 2^-12) = 1 / sqrt(257). fp16 (ulp 2^-10 at 1) keeps both."""

import json
import math

import pytest

torch = pytest.importorskip("torch")
from safetensors.torch import save_file  # noqa: E402

from experiments.phase2 import merge_rounding as mr  # noqa: E402

W = torch.ones(2, 2, dtype=torch.bfloat16)
D = torch.tensor([[2.0**-10, 0.0], [0.0, 2.0**-6]])


def test_planted_case_bf16():
    row = mr.table_row([(W, D)], torch.bfloat16)
    assert row["erased_share"] == 0.5  # 1 of the 2 non-zero delta elements; the zeros are not counted
    assert row["rel_err"] == pytest.approx(1 / math.sqrt(257), rel=1e-6)
    assert row["median_d_over_ulp"] == pytest.approx((0.125 + 2.0) / 2)  # |d| / 2^-7 over the non-zero elements


def test_fp16_keeps_the_planted_delta():
    row = mr.table_row([(W, D)], torch.float16)
    assert row["erased_share"] == 0.0 and row["rel_err"] == 0.0 < mr.table_row([(W, D)], torch.bfloat16)["rel_err"]


def test_zero_delta_gives_zeros_not_nan():
    row = mr.table_row([(W, torch.zeros(2, 2))], torch.bfloat16)
    assert row == {"erased_share": 0.0, "rel_err": 0.0, "median_d_over_ulp": 0.0}


def test_erased_uses_the_exact_rounding_not_half_an_ulp():
    """A tie (|d| = ulp/2) rounds to even and is erased at W = 1. Below W = 1 the spacing is 2^-8, half the ulp
    above, so d = -0.6 x 2^-8 moves W to 1 - 2^-8 although |d| < ulp(1)/2 = 2^-8 (math-stats note 8)."""
    w = torch.tensor([1.0, 1.0], dtype=torch.bfloat16)
    assert mr.table_row([(w, torch.tensor([2.0**-8, 0.0]))], torch.bfloat16)["erased_share"] == 1.0
    w = torch.tensor([1.0], dtype=torch.bfloat16)
    assert mr.table_row([(w, torch.tensor([-0.6 * 2.0**-8]))], torch.bfloat16)["erased_share"] == 0.0


def test_main_names_step_layers_and_modules(tmp_path, capsys):
    """A two-layer model folder and two adapter checkpoints in the train_rl layout (r 2, alpha 4: scale 2)."""
    mods = ["self_attn.q_proj", "mlp.down_proj"]
    g = torch.Generator().manual_seed(0)
    model = {f"model.language_model.layers.{i}.{m}.weight": torch.randn(4, 4, generator=g).to(torch.bfloat16)
             for i in range(2) for m in mods}  # fmt: skip
    (tmp_path / "final").mkdir()
    save_file(model, tmp_path / "final/model.safetensors")
    for step in (10, 20):
        ck = tmp_path / f"trainer/checkpoint-{step}"
        ck.mkdir(parents=True)
        (ck / "adapter_config.json").write_text(json.dumps({"r": 2, "lora_alpha": 4}))
        save_file({f"base_model.model.model.language_model.layers.{i}.{m}.lora_{x}.weight":
                   torch.randn(2, 4, generator=g) * 1e-3 if x == "A" else torch.randn(4, 2, generator=g) * 1e-3
                   for i in range(2) for m in mods for x in "AB"}, ck / "adapter_model.safetensors")  # fmt: skip
    out = tmp_path / "table.json"
    assert mr.main(["--run", str(tmp_path), "--steps", "10", "20", "--out", str(out)]) == 0
    doc = json.loads(out.read_text())
    assert [r["step"] for r in doc["rows"]] == [10, 20]
    assert doc["layers"] == 2 and doc["modules"] == ["down_proj", "q_proj"]
    for r in doc["rows"]:
        assert set(r) == {"step", "bf16", "fp16"} and 0 <= r["bf16"]["erased_share"] <= 1
    printed = capsys.readouterr().out
    assert "2 layers x 2 modules" in printed and "\n  10 " in printed and "\n  20 " in printed

"""B2.3b Mac part (G04's CPU half): the four trainer-side ways of scoring the same answers on a tiny random Qwen3-VL
(tests/experiments/test_trl_adapter_equivalence.tiny_model) with a LoRA adapter; the real 2B run is the `data` test at
the end."""

import hashlib
import math
import os
from pathlib import Path

import pytest

torch = pytest.importorskip("torch")
pytest.importorskip("peft")
from peft import LoraConfig, get_peft_model  # noqa: E402

from experiments.phase2 import b2_attribution as B  # noqa: E402
from tests.experiments.test_trl_adapter_equivalence import GRPOTrainer, group_batch, tiny_model  # noqa: E402

WAYS = ("a", "b", "c", "d")


def peft_model(dtype=torch.float32, delta=True):
    """tiny_model with base weights in `dtype` and an fp32 LoRA (r 2) on q_proj / v_proj; delta: random lora_B (PEFT
    starts it at 0)."""
    m = get_peft_model(tiny_model().to(dtype), LoraConfig(r=2, lora_alpha=4, target_modules=["q_proj", "v_proj"]))
    if delta:
        g = torch.Generator().manual_seed(1)
        for n, p in m.named_parameters():
            if "lora_B" in n:
                p.data = torch.randn(p.shape, generator=g) * 0.05
    return m.eval()


def batches():
    return [group_batch(groups=1, G=4, lengths=[3, 7, 1, 5], seed=s) for s in (3, 4)]


def test_zero_lora_delta_merged_equals_lora():
    lp = B.score_ways(peft_model(delta=False), batches(), lambda: tiny_model())
    for x, y in zip(lp["c"], lp["b"], strict=True):
        torch.testing.assert_close(x, y, atol=1e-6, rtol=0)


def test_ways_are_finite_repeatable_and_leave_the_model_as_it_was():
    m = peft_model(torch.bfloat16)
    before = {n: p.detach().clone() for n, p in m.named_parameters()}
    load = lambda: tiny_model().to(torch.bfloat16)  # noqa: E731
    one, two = B.score_ways(m, batches(), load), B.score_ways(m, batches(), load)
    assert set(one) == set(WAYS)
    for w in WAYS:
        for x, y in zip(one[w], two[w], strict=True):
            assert x.dtype == torch.float32 and torch.isfinite(x).all()
            torch.testing.assert_close(x, y, atol=1e-4, rtol=0)
    after = dict(m.named_parameters())
    assert all(after[n].dtype == p.dtype and torch.equal(after[n], p) for n, p in before.items())


def test_d_is_the_plain_fp32_forward_without_autocast():
    """Oracle: TRL's own forward on an fp32 copy of the same weights, no autocast."""
    m = peft_model(torch.bfloat16)
    got = B.score_ways(m, batches(), lambda: tiny_model().to(torch.bfloat16))["d"]
    ref = peft_model(torch.bfloat16)
    ref.float()
    t = B.forward_obj(GRPOTrainer, bf16=False)
    for (ids, mask, T, kw), x in zip(batches(), got, strict=True):
        with torch.no_grad():
            want = t._get_per_token_logps_and_entropies(ref, ids, mask, T, batch_size=len(ids), **kw)[0]
        torch.testing.assert_close(x, want.float(), atol=1e-6, rtol=0)


def test_b_and_c_differ_from_d_under_bf16():
    lp = B.score_ways(peft_model(torch.bfloat16), batches(), lambda: tiny_model().to(torch.bfloat16))
    mask = [m[:, -T:].bool() for _, m, T, _ in batches()]
    for w in ("b", "c"):
        assert any((x - d)[k].abs().max() > 1e-4 for x, d, k in zip(lp[w], lp["d"], mask, strict=True)), w


def test_report_signed_and_absolute_gaps_by_64_token_bins():
    """Two rows of 70 answer tokens, the second row's last 6 masked; way x is d + 0.1 on tokens 0-63 and d - 0.3 on
    64-69: signed mean (128 x 0.1 - 6 x 0.3) / 134, mean |gap| (128 x 0.1 + 6 x 0.3) / 134."""
    mask = torch.ones(2, 70, dtype=torch.long)
    mask[1, 64:] = 0
    d = torch.randn(2, 70, generator=torch.Generator().manual_seed(0))
    x = d.clone()
    x[:, :64] += 0.1
    x[:, 64:] -= 0.3
    r = B.report({"x": [x], "d": [d]}, [mask], pairs=(("x", "d"),))["x-d"]
    assert r["signed"] == pytest.approx((128 * 0.1 - 6 * 0.3) / 134, abs=1e-6)
    assert r["abs"] == pytest.approx((128 * 0.1 + 6 * 0.3) / 134, abs=1e-6)
    assert [b["n"] for b in r["bins"]] == [128, 6]
    assert r["bins"][0]["signed"] == pytest.approx(0.1, abs=1e-6)
    assert r["bins"][1]["abs"] == pytest.approx(0.3, abs=1e-6)
    assert all(math.isfinite(v) for b in r["bins"] for v in (b["signed"], b["abs"]))


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


@pytest.mark.data
def test_g04_mac_part_on_the_full_run(tmp_path):
    """G04 CPU half on the real checkpoint (1 prompt x 2 answers): finite gaps for every pair, adapter unchanged."""
    root = Path(os.environ.get("OBJDEL_RUNS", Path(__file__).resolve().parents[2] / "runs"))
    adapter = root / "rl/train_rl/trainer/checkpoint-100"
    if not (adapter / "adapter_model.safetensors").exists():
        pytest.skip(f"no adapter at {adapter}")
    h = sha(adapter / "adapter_model.safetensors")
    out = tmp_path / "attr.json"
    assert B.main(["--adapters", str(adapter), "--prompts", "1", "--k", "2", "--out", str(out)]) == 0
    doc = __import__("json").loads(out.read_text())
    pairs = doc["adapters"][str(adapter)]
    assert set(pairs) == {f"{a}-{b}" for a, b in B.PAIRS}
    assert all(math.isfinite(p["signed"]) and math.isfinite(p["abs"]) and p["bins"] for p in pairs.values())
    assert sha(adapter / "adapter_model.safetensors") == h


def test_default_adapters_are_the_plan_checkpoints_100_and_770():
    import inspect

    src = inspect.getsource(B.main)
    assert 'RUN / "checkpoint-100", RUN / "checkpoint-770"' in src and "checkpoint-750" not in src

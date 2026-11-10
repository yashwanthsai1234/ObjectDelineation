"""FastGRPO._prepare_inputs (E5.3): the update micro-batches it cuts and when it falls back to TRL's own method.

I8: on the tiny offline model (CPU, seed 0, dr_grpo) the micro-batch index lists and the logged losses equal
tests/goldens/phase_e/microbatches.json for shared_prefix off and on, and the old swap path (TRL's _prepare_inputs
with its module global swapped for the order, removed in E5.3) gives the same lists and losses.
I9: under each fallback condition FastGRPO's output equals GRPOTrainer's with the same seed; the dr_grpo boundaries
do not fall back; TRL's module global shuffle_sequence_dict is never left changed."""
import functools
import json
import platform
import sys
from pathlib import Path

import pytest

torch = pytest.importorskip("torch")
pytest.importorskip("trl")
import trl.trainer.grpo_trainer as G  # noqa: E402
from trl import GRPOTrainer  # noqa: E402

from objdel.rl import forward as rl_forward  # noqa: E402
from objdel.rl import trl_adapter as rl_trl_adapter  # noqa: E402
from tests.experiments.test_trl_adapter_equivalence import (  # noqa: E402
    Plain,
    c10_trainer,
    generation_batch,
    order_spy,
    tiny_loop_vllm_fake,
    train_logs,
    trainer,
)

GOLDEN = Path(__file__).resolve().parents[1] / "goldens/phase_e/microbatches.json"
PLATFORM = f"{sys.platform}-{platform.machine()}"
SPG = 4  # micro-batches per generation batch: 4 prompts x 4 answers, 4 rows per micro-batch
ORIGINAL_SHUFFLE = G.shuffle_sequence_dict


# I8 ------------------------------------------------------------------------------------------------------------------

def index_order(seen: list, order):
    """`order` (by_length or by_group_length) that also records each generation batch's micro-batches as row-index
    lists (the rows in the order _generate_and_score_completions returned them; split_tensor_dict cuts consecutive
    equal chunks)."""
    def spy(batch, *a, **k):
        n = len(batch["completion_ids"])
        out = order({**batch, "row": torch.arange(n)}, *a, **k)
        rows, m = out.pop("row").tolist(), n // SPG
        seen.append([rows[i * m:(i + 1) * m] for i in range(SPG)])
        return out

    return spy


def microbatch_run(monkeypatch, tmp_path, shared_prefix: bool, swap: bool = False) -> dict:
    """Two steps of the tiny loop with 4 micro-batches per step; swap: the pre-E5.3 path instead of the override."""
    seen = []
    monkeypatch.setattr(rl_trl_adapter, "by_length", index_order(seen, rl_forward.by_length))
    monkeypatch.setattr(rl_trl_adapter, "by_group_length", index_order(seen, rl_forward.by_group_length))
    if swap:
        def old_path(self, generation_batch):
            order = (functools.partial(rl_trl_adapter.by_group_length, G=self.num_generations) if self.shared_prefix
                     else rl_trl_adapter.by_length)  # fmt: skip
            with monkeypatch.context() as m:
                m.setattr(G, "shuffle_sequence_dict", order)
                return GRPOTrainer._prepare_inputs(self, generation_batch)

        monkeypatch.setattr(rl_trl_adapter.FastGRPO, "_prepare_inputs", old_path)
    t = tiny_loop_vllm_fake(monkeypatch, tmp_path, steps=2, shared_prefix=shared_prefix, num_generations=4,
                            per_device_train_batch_size=4, gradient_accumulation_steps=SPG)
    assert G.shuffle_sequence_dict is ORIGINAL_SHUFFLE
    return {"index_lists": seen, "loss": [float(h["loss"]).hex() for h in train_logs(t)]}


def record() -> None:
    """Write microbatches.json (index lists for every platform, losses for this one).

    See tests/goldens/phase_e/README.
    """
    import tempfile

    from _pytest.monkeypatch import MonkeyPatch

    old = json.loads(GOLDEN.read_text()) if GOLDEN.exists() else {"index_lists": {}, "loss": {}}
    for sp in (False, True):
        with MonkeyPatch.context() as mp, tempfile.TemporaryDirectory() as tmp:
            got = microbatch_run(mp, Path(tmp), sp)
        old["index_lists"][f"shared_prefix={sp}"] = got["index_lists"]
        old["loss"].setdefault(PLATFORM, {})[f"shared_prefix={sp}"] = got["loss"]
    GOLDEN.write_text(json.dumps(old, indent=1, sort_keys=True) + "\n")


@pytest.mark.parametrize("shared_prefix", [False, True])
def test_i8_micro_batches_and_loss_equal_the_golden(monkeypatch, tmp_path, shared_prefix):
    golden, key = json.loads(GOLDEN.read_text()), f"shared_prefix={shared_prefix}"
    got = microbatch_run(monkeypatch, tmp_path, shared_prefix)
    assert got["index_lists"] == golden["index_lists"][key]
    if shared_prefix:  # whole groups of 4 answers (one prompt) per micro-batch
        assert all(sorted(b) == list(range(b[0] // 4 * 4, b[0] // 4 * 4 + 4)) for g in got["index_lists"] for b in g)
    if PLATFORM not in golden["loss"]:
        pytest.skip(f"no loss golden recorded for {PLATFORM}")
    assert got["loss"] == golden["loss"][PLATFORM][key]


@pytest.mark.parametrize("shared_prefix", [False, True])
def test_i8_old_swap_path_gives_the_same_micro_batches_and_loss(monkeypatch, tmp_path, shared_prefix):
    new = microbatch_run(monkeypatch, tmp_path / "new", shared_prefix)
    assert microbatch_run(monkeypatch, tmp_path / "old", shared_prefix, swap=True) == new


# I9 ------------------------------------------------------------------------------------------------------------------

def twins(monkeypatch, case: str | None, gas: int = 4):
    """A FastGRPO and a plain GRPOTrainer with the same attributes (c10_trainer's), 8 rows per generation batch in 4
    micro-batches, with one fallback condition set (None: none)."""
    batch = generation_batch()
    monkeypatch.setattr(GRPOTrainer, "_generate_and_score_completions", lambda self, b: dict(batch))
    fast = c10_trainer(steps_per_generation=SPG)
    fast.args.gradient_accumulation_steps, fast.args.per_device_train_batch_size = gas, 2
    fast.shared_prefix, fast.num_generations, fast._step, fast._buffered_inputs = False, 4, 0, None
    if case == "eval":
        fast.model.training = False
    elif case == "loss_type":
        fast.loss_type = "bnpo"
    elif case == "entropy":
        fast.top_entropy_quantile = 0.2
    elif case == "accum":
        fast.args.gradient_accumulation_steps = 6  # 6 % 4: the micro-batches span two optimizer steps
    attrs = {k: v for k, v in vars(fast).items() if k not in ("accelerator", "_metrics", "state", "model")}
    plain = trainer(Plain, **attrs)
    plain.model.training = fast.model.training
    return fast, plain


def outputs(t) -> list[dict]:
    """t._prepare_inputs for every micro-batch of one generation batch, seed 0 (TRL's shuffle draws torch.randperm)."""
    torch.manual_seed(0)
    out = []
    for step in range(SPG):
        t._step = step
        out.append(t._prepare_inputs({}))
    return out


def assert_same(a: dict, b: dict):
    assert a.keys() == b.keys()
    for k in a:
        assert torch.equal(a[k], b[k]) if isinstance(a[k], torch.Tensor) else a[k] == b[k], k


@pytest.mark.parametrize("case", ["eval", "loss_type", "entropy", "accum"])
def test_i9_fallback_equals_trl(monkeypatch, case):
    fast, plain = twins(monkeypatch, case)
    for a, b in zip(outputs(fast), outputs(plain), strict=True):
        assert_same(a, b)
    assert G.shuffle_sequence_dict is ORIGINAL_SHUFFLE


@pytest.mark.parametrize("gas", [4, 8])  # gas == spg and gas == 2 * spg: every micro-batch in one optimizer step
def test_i9_dr_grpo_boundaries_do_not_fall_back(monkeypatch, gas):
    t = c10_trainer(steps_per_generation=SPG)
    t.args.gradient_accumulation_steps = gas
    assert order_spy(monkeypatch, t)() == "by_length"


def test_i9_without_a_condition_the_rows_are_a_reordering_of_trls(monkeypatch):
    fast, plain = twins(monkeypatch, None)
    a, b = outputs(fast), outputs(plain)
    ids = [torch.cat([x["completion_ids"] for x in o]) for o in (a, b)]
    assert not torch.equal(*ids)  # by_length, not TRL's seeded shuffle
    assert sorted(map(tuple, ids[0].tolist())) == sorted(map(tuple, ids[1].tolist()))


def test_i9_trls_global_is_untouched_after_a_fast_build_and_a_plain_trainer_uses_it(monkeypatch, tmp_path):
    tiny_loop_vllm_fake(monkeypatch, tmp_path, steps=1)  # builds and trains a real FastGRPO
    assert G.shuffle_sequence_dict is ORIGINAL_SHUFFLE
    _, plain = twins(monkeypatch, None)
    calls = []
    monkeypatch.setattr(G, "shuffle_sequence_dict", lambda b: calls.append(1) or ORIGINAL_SHUFFLE(b))
    outputs(plain)
    assert calls == [1]  # one shuffle per generation batch, through TRL's module global


@pytest.mark.parametrize("where", ["generate", "order"])
def test_i9_a_callee_raising_mid_prepare_inputs_leaves_trls_global(monkeypatch, where):
    fast, _ = twins(monkeypatch, None)

    def boom(*a, **k):
        raise RuntimeError(where)

    if where == "generate":
        monkeypatch.setattr(GRPOTrainer, "_generate_and_score_completions", boom)
    else:
        monkeypatch.setattr(rl_trl_adapter, "by_length", boom)
    with pytest.raises(RuntimeError, match=where):
        fast._prepare_inputs({})
    assert G.shuffle_sequence_dict is ORIGINAL_SHUFFLE

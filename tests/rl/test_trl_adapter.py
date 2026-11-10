"""objdel.rl.trl_adapter: the start-up check of the private TRL names and silent couplings FastGRPO relies on (A7).

Every case patches the real TRL classes with monkeypatch (never the library files). A passing check alone cannot catch a
typo in a checked substring; the negative cases do.
"""

import inspect
import math
import os
import re
import subprocess
import sys
import textwrap
from pathlib import Path
from types import SimpleNamespace

import pytest

torch = pytest.importorskip("torch")
trl = pytest.importorskip("trl")
from trl import GRPOTrainer  # noqa: E402
from trl.generation.vllm_generation import VLLMGeneration  # noqa: E402

from objdel.policy import caches as policy_caches  # noqa: E402
from objdel.rl import forward as rl_forward  # noqa: E402
from objdel.rl import trl_adapter as T  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]


def test_the_environment_has_the_pinned_trl():
    assert trl.__version__ == T.PIN == "1.17.1"


def test_check_passes_on_the_pin_without_vllm_and_quickly():
    code = ("import sys, time; sys.modules['vllm'] = None\n"
            "import trl.trainer.grpo_trainer\n"
            "from objdel.rl import trl_adapter as T\n"
            "for _ in range(2):\n"
            "    t = time.perf_counter(); r = T.check(); dt = time.perf_counter() - t\n"
            "    assert r is None and dt < 1.0, (r, dt)\n"
            "assert sys.modules['vllm'] is None\n")  # fmt: skip
    env = {k: v for k, v in os.environ.items() if k != "OBJDEL_TRL_UNPINNED"}
    p = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, env=env)
    assert p.returncode == 0, p.stderr[-2000:]


def test_all_items_is_the_one_count_of_the_contract():
    parts = (T.METHODS, T.GLOBALS, T.VLLM_METHODS, T.PRIVATE_ATTRS, T.ATTRS, T.KEYS, T.BEHAVIOURS, T.LORA_POINTS,
             T.SOURCE_SHA256, T.TRAINER_ATTRS)
    keys = [(kind, name) for kind, name, _ in T.ALL_ITEMS]
    assert len(keys) == sum(map(len, parts)) and len(set(keys)) == len(keys)
    # Phase A 29 (6 + 1 + 1 + 2 + 13 + 4 + 2), plus C5's _llm_weights_sleeping, plus BR1-rl's 3 reward-kwarg behaviours,
    # plus B2.1's set_is_ratio reads (attribute off_policy_mask_threshold, key advantages), plus B2.2 F step 3's 4
    # option-F vLLM names (add_lora, remove_lora, lora_index_to_id, lora_a_stacked / lora_b_stacked; lora-seqmask only),
    # plus E5's source hash of each of the 6 overridden methods and the 2 trainer attributes of _prepare_inputs, minus
    # the global shuffle_sequence_dict (E5 step 3: no swap any more; the _prepare_inputs hash guards the copy)
    assert set(T.SOURCE_SHA256) == set(T.METHODS) and set(T.TRAINER_ATTRS) == {"_step", "_buffered_inputs"}
    assert T.GLOBALS == {} and len(T.ALL_ITEMS) == 46
    doc = ROOT / "docs/reference/trl-coupling.md"  # A9 writes it: one table row per item (header and rule excluded)
    assert sum(x.startswith("|") for x in doc.read_text().splitlines()) - 2 == len(T.ALL_ITEMS)


def test_the_coupling_page_has_one_row_per_item_with_the_case_that_breaks_it():
    """trl-coupling.md: each row's item, kind and place are the code's; its failing test is a real BREAKS case."""
    rows = [[c.strip().strip("`") for c in ln.strip("|\n").split("|")] for ln in
            (ROOT / "docs/reference/trl-coupling.md").read_text().splitlines() if ln.startswith("| `")]  # fmt: skip
    code = {(name, kind): place for kind, name, place in T.ALL_ITEMS}
    assert {(r[0], r[1]) for r in rows} == set(code) and len(rows) == len(code)
    cases = {b[0] for b in BREAKS}
    live = {"option-F vLLM name": {"test_check_lora_names_the_missing_option_f_vllm_name"},  # one case per name
            "trainer attribute": {"test_check_trainer_names_the_missing_attribute"}}
    for name, kind, _, place, tests in rows:
        assert place == code[name, kind].replace(", ", "`, `"), name
        assert all(t.strip("`") in live.get(kind, cases) for t in tests.split(", ")), name


def test_no_live_python_file_says_trl_1_16():
    stale = [p for d in ("src", "experiments") for p in (ROOT / d).rglob("*.py")
             if "lab" not in p.relative_to(ROOT).parts and "TRL " + "1.16" in p.read_text()]  # fmt: skip
    assert stale == []


def drop(monkeypatch, owner, needle):
    """inspect.getsource(owner) without `needle` (every other object's source unchanged)."""
    real = inspect.getsource
    monkeypatch.setattr(inspect, "getsource", lambda o: real(o).replace(needle, "") if o is owner else real(o))


def setsig(monkeypatch, name, fn):
    monkeypatch.setattr(GRPOTrainer, name, fn)


def edit(monkeypatch, owner):
    """inspect.getsource(owner) with one more line: a body that TRL changed (every other object's source unchanged)."""
    real = inspect.getsource
    monkeypatch.setattr(inspect, "getsource", lambda o: real(o) + "        pass\n" if o is owner else real(o))


METHOD_NAMES = ("_generate_and_score_completions", "_prepare_inputs", "_compute_loss",
                "_get_per_token_logps_and_entropies", "_generate", "_calculate_rewards")
GASC, LOSS, CALC = (GRPOTrainer._generate_and_score_completions, GRPOTrainer._compute_loss,
                    GRPOTrainer._calculate_rewards)
ATTR_NAMES = ("use_vllm", "loss_type", "num_iterations", "num_generations", "top_entropy_quantile",
              "vllm_importance_sampling_correction", "vllm_importance_sampling_mode",
              "vllm_importance_sampling_clip_min", "vllm_importance_sampling_clip_max", "vllm_generation",
              "off_policy_mask_threshold")  # fmt: skip
BREAKS = [
    ("delete _prepare_inputs", lambda m: m.delattr(GRPOTrainer, "_prepare_inputs"), "_prepare_inputs"),
    ("_compute_loss extra param", lambda m: setsig(m, "_compute_loss", lambda self, model, inputs, extra: 0),
     "_compute_loss"),
    ("_compute_loss renamed param", lambda m: setsig(m, "_compute_loss", lambda self, model, batch: 0),
     "_compute_loss"),
    ("_compute_loss reordered", lambda m: setsig(m, "_compute_loss", lambda self, inputs, model: 0), "_compute_loss"),
    ("logps without sampling_mask", lambda m: setsig(m, "_get_per_token_logps_and_entropies",
     lambda self, model, input_ids, attention_mask, logits_to_keep, batch_size=None, compute_entropy=False,
     compute_aux_loss=False, pixel_values=None, image_grid_thw=None, num_images=None, pixel_attention_mask=None,
     spatial_shapes=None, num_tiles=None, image_sizes=None, token_type_ids=None, mm_token_type_ids=None,
     image_position_ids=None: 0), "_get_per_token_logps_and_entropies"),
    ("logps two names swapped", lambda m: setsig(m, "_get_per_token_logps_and_entropies",
     lambda self, model, attention_mask, input_ids, logits_to_keep, batch_size=None, compute_entropy=False,
     compute_aux_loss=False, pixel_values=None, image_grid_thw=None, num_images=None, pixel_attention_mask=None,
     spatial_shapes=None, num_tiles=None, image_sizes=None, token_type_ids=None, mm_token_type_ids=None,
     image_position_ids=None, sampling_mask=None: 0), "_get_per_token_logps_and_entropies"),
    ("delete _iter_named_params", lambda m: m.delattr(VLLMGeneration, "_iter_named_params"), "_iter_named_params"),
    ("__init__ without self._metrics", lambda m: drop(m, GRPOTrainer.__init__, "self._metrics"), "_metrics"),
    ("VLLMGeneration without self._dist", lambda m: drop(m, VLLMGeneration, "self._dist"), "_dist"),
    ("VLLMGeneration without sync_weights", lambda m: drop(m, VLLMGeneration, "def sync_weights"), "sync_weights"),
    ("VLLMGeneration without self.llm", lambda m: drop(m, VLLMGeneration, "self.llm"), "llm"),
    ("VLLMGeneration without self.model", lambda m: drop(m, VLLMGeneration, "self.model"), "model"),
    *[(f"no self.{n} =", lambda m, n=n: drop(m, GRPOTrainer, f"self.{n} ="), n) for n in ATTR_NAMES],
    *[(f"no {k}", lambda m, k=k: drop(m, GASC, f'"{k}"'), k)
      for k in ("sampling_per_token_logps", "completion_mask", "tool_mask", "importance_sampling_ratio", "advantages")],
    ("_compute_loss source without the ratio key", lambda m: drop(m, LOSS, '"importance_sampling_ratio"'),
     "importance_sampling_ratio"),
    ("no C10 gate", lambda m: drop(m, GASC, "self.use_vllm and self.vllm_importance_sampling_correction"),
     "old-logprob pass"),
    ("no ratio block", lambda m: drop(m, GASC, "per_token_logps_diff = (old_per_token_logps"), "IS ratio block"),
    ("_generate_and_score_completions extra param",
     lambda m: setsig(m, "_generate_and_score_completions", lambda self, inputs, extra: 0),
     "_generate_and_score_completions"),
    ("_generate renamed param", lambda m: setsig(m, "_generate", lambda self, prompt: 0), "_generate"),
    ("_calculate_rewards without completion_ids_list",
     lambda m: setsig(m, "_calculate_rewards", lambda self, inputs, prompts, completions: 0), "_calculate_rewards"),
    ("VLLMGeneration without self._llm_weights_sleeping",
     lambda m: drop(m, VLLMGeneration, "self._llm_weights_sleeping ="), "_llm_weights_sleeping"),
    *[(f"no {n} to reward functions", lambda m, t=t: drop(m, CALC, t), f"{n} to reward functions")
      for n, t in (("log_metric", 'reward_kwargs["log_metric"]'), ("log_extra", 'reward_kwargs["log_extra"]'),
                   ("completion_ids", "completion_ids=completion_ids_list"))],
    *[(f"changed source of {n}", lambda m, n=n: edit(m, getattr(GRPOTrainer, n)), n) for n in METHOD_NAMES],
]  # fmt: skip


@pytest.mark.parametrize("breaker, needle", [b[1:] for b in BREAKS], ids=[b[0] for b in BREAKS])
def test_each_broken_item_raises_naming_it(monkeypatch, breaker, needle):
    monkeypatch.delenv("OBJDEL_TRL_UNPINNED", raising=False)
    breaker(monkeypatch)
    with pytest.raises(T.TRLContractError, match=rf"\b{re.escape(needle)}\b.*Pin is 1\.17\.1"):
        T.check()


@pytest.mark.parametrize("first, later, needle", [("delete _prepare_inputs", "no tool_mask", "_prepare_inputs"),
                                                  ("no self.use_vllm =", "changed source of _prepare_inputs",
                                                   "use_vllm")])  # fmt: skip
def test_the_first_broken_item_in_check_order_is_named(monkeypatch, first, later, needle):
    breaks = {name: fn for name, fn, _ in BREAKS}
    breaks[later](monkeypatch)
    breaks[first](monkeypatch)
    with pytest.raises(T.TRLContractError, match=needle):
        T.check()


@pytest.mark.parametrize("env", [None, "0", "true"])
@pytest.mark.parametrize("version", ["1.18.0", "1.17.1.dev0", "1.17.1+cu128"])  # strict equality with the pin
def test_another_version_raises_unless_opted_in(monkeypatch, env, version):
    monkeypatch.delenv("OBJDEL_TRL_UNPINNED", raising=False)
    if env is not None:
        monkeypatch.setenv("OBJDEL_TRL_UNPINNED", env)
    monkeypatch.setattr(trl, "__version__", version)
    with pytest.raises(T.TRLContractError, match=rf"{re.escape(version)}.*1\.17\.1"):
        T.check()


def test_the_opt_in_waives_only_the_version_with_one_warning(monkeypatch):
    monkeypatch.setenv("OBJDEL_TRL_UNPINNED", "1")
    monkeypatch.setattr(trl, "__version__", "1.18.0")
    with pytest.warns(UserWarning, match="1.18.0") as w:
        assert T.check() is None
    assert len(w) == 1
    monkeypatch.delattr(GRPOTrainer, "_prepare_inputs")  # a broken item is never waived, and wins over the version
    with pytest.raises(T.TRLContractError, match="_prepare_inputs"):
        T.check()


def test_the_version_is_checked_last(monkeypatch):
    monkeypatch.delenv("OBJDEL_TRL_UNPINNED", raising=False)
    monkeypatch.setattr(trl, "__version__", "1.18.0")
    drop(monkeypatch, GASC, '"tool_mask"')
    with pytest.raises(T.TRLContractError, match="tool_mask") as e:
        T.check()
    assert "1.18.0" not in str(e.value)


def trl_ratio_block() -> str:
    """TRL's own importance-ratio code: from the anchor line to the `else:` before the unknown-mode error."""
    src = inspect.getsource(GASC)
    for anchor in (T.BEHAVIOURS["IS ratio block"][1], "Unknown vLLM importance sampling level"):
        assert src.count(anchor) == 1, f"TRL ratio block anchor moved: {anchor}"
    lines = src.splitlines()
    i0 = next(i for i, x in enumerate(lines) if T.BEHAVIOURS["IS ratio block"][1] in x)
    i1 = next(i for i, x in enumerate(lines) if "Unknown vLLM importance sampling level" in x)
    while lines[i1].strip() != "else:":
        i1 -= 1
    return textwrap.dedent("\n".join(lines[i0:i1]))


@pytest.mark.parametrize("lo", [None, 0.1])
@pytest.mark.parametrize("mode", ["sequence_mask", "sequence_truncate", "token_mask", "token_truncate"])
def test_trls_ratio_block_equals_fast_is_ratio(mode, lo):
    code, g = compile(trl_ratio_block(), "trl_ratio_block", "exec"), torch.Generator().manual_seed(0)
    stub = SimpleNamespace(vllm_importance_sampling_mode=mode, vllm_importance_sampling_clip_min=lo,
                           vllm_importance_sampling_clip_max=3.0)  # fmt: skip
    for _ in range(200):
        old = torch.randn(4, 12, generator=g) * 5
        samp = old + (torch.rand(4, 12, generator=g) * 2 - 1) * 15 * torch.rand(4, 1, generator=g) ** 4
        samp[1, 3] = float("nan")  # a token vLLM could not score
        mask = (torch.rand(4, 12, generator=g) < 0.9).float()
        mask[0], mask[2, 8:] = 0.0, 0.0  # an all-masked row, a masked tail
        ns = dict(old_per_token_logps=old, sampling_per_token_logps=samp, mask=mask, torch=torch, math=math, self=stub)
        exec(code, ns)
        assert torch.equal(ns["vllm_importance_sampling_ratio"], rl_forward.is_ratio(old, samp, mask, mode, lo, 3.0))


def test_fast_grpo_checks_the_contract_before_anything_else(monkeypatch):
    calls = []

    def init(self, *a, **k):
        calls.append("GRPOTrainer.__init__")
        raise RuntimeError("stop after the order is recorded")

    monkeypatch.setattr(GRPOTrainer, "__init__", init)
    monkeypatch.setattr(T, "check", lambda: calls.append("check"))
    with pytest.raises(RuntimeError, match="stop after"):
        T.FastGRPO()
    assert calls == ["check", "GRPOTrainer.__init__"]

    def broken():
        raise T.TRLContractError("boom")

    calls.clear()
    monkeypatch.setattr(T, "check", broken)
    with pytest.raises(T.TRLContractError, match="boom"):
        T.FastGRPO()
    assert calls == []


def test_fast_grpo_refuses_another_trl_version_without_the_opt_in(monkeypatch):
    monkeypatch.delenv("OBJDEL_TRL_UNPINNED", raising=False)
    monkeypatch.setattr(trl, "__version__", "1.18.0")
    with pytest.raises(T.TRLContractError, match="1.18.0"):
        T.FastGRPO()


def lora_engine(**drop):
    """A fake vg.llm carrying the four option-F vLLM names; `drop` sets one of them to None or a wrong value."""
    stacked = SimpleNamespace(lora_a_stacked=(), lora_b_stacked=())
    model = SimpleNamespace(named_modules=lambda: iter([("", SimpleNamespace()), ("layers.0.qkv_proj", stacked)]))
    if "stacked" in drop:
        model = SimpleNamespace(named_modules=lambda: iter([("layers.0.qkv_proj", SimpleNamespace(lora_a_stacked=()))]))
    manager = SimpleNamespace(lora_index_to_id=drop.get("index", [None, 1]))
    runner = SimpleNamespace(model=model, lora_manager=SimpleNamespace(_adapter_manager=manager))
    engine = SimpleNamespace(add_lora=drop.get("add", lambda r: True), remove_lora=drop.get("remove", lambda i: True),
                             model_executor=SimpleNamespace(driver_worker=SimpleNamespace(model_runner=runner)))
    return SimpleNamespace(llm_engine=engine)


@pytest.mark.parametrize("drop,item", [({"add": None}, "add_lora"), ({"remove": 3}, "remove_lora"),
                                       ({"index": None}, "lora_index_to_id"),
                                       ({"stacked": True}, "lora_a_stacked / lora_b_stacked")])
def test_check_lora_names_the_missing_option_f_vllm_name(drop, item):
    T.check_lora(lora_engine())
    with pytest.raises(T.TRLContractError, match=f"option-F vLLM name {re.escape(item)}:"):
        T.check_lora(lora_engine(**drop))


def test_fast_grpo_checks_the_option_f_names_only_under_lora_seqmask(monkeypatch):
    import contextlib

    seen = []
    monkeypatch.setattr(T, "check_lora", lambda llm: seen.append(llm))
    monkeypatch.setattr(GRPOTrainer, "__init__", lambda self, *a, **k: None)
    monkeypatch.setattr(T, "check", lambda: None)
    # no vLLM LLM class on the Mac
    monkeypatch.setattr(T, "vllm_kwargs", lambda **kw: contextlib.nullcontext())
    for name in ("check_engine_dtype", "drop_lora_input_cast", "sync_without_merge"):
        monkeypatch.setattr(T.FastGRPO, name, lambda self: None)
    monkeypatch.setattr(policy_caches.VisionCache, "install", classmethod(lambda cls, m: None))
    monkeypatch.setattr(policy_caches.RopeCache, "install", classmethod(lambda cls, m: None))
    monkeypatch.setattr(policy_caches.ImageMemo, "install", classmethod(lambda cls, p: None))
    for lora in (False, True):
        t = object.__new__(T.FastGRPO)
        t.args, t.vllm_generation, t.model, t.processing_class = SimpleNamespace(vllm_lora=lora), \
            SimpleNamespace(llm="engine"), None, None
        vars(t).update(_step=0, _buffered_inputs=None)  # TRL's __init__ sets both (check_trainer)
        T.FastGRPO.__init__(t, args=t.args)
    assert seen == ["engine"]


def test_a_changed_method_body_names_the_method_its_override_and_the_hash(monkeypatch):
    monkeypatch.delenv("OBJDEL_TRL_UNPINNED", raising=False)
    edit(monkeypatch, GRPOTrainer._prepare_inputs)
    with pytest.raises(T.TRLContractError, match=r"method source _prepare_inputs: expected '[0-9a-f]{64}'.*"
                                                 r"FastGRPO\._prepare_inputs\. Pin is 1\.17\.1"):
        T.check()


@pytest.mark.parametrize("item", ["_step", "_buffered_inputs"])
def test_check_trainer_names_the_missing_attribute(item):
    """At start-up FastGRPO checks the built trainer for the attributes the copied _prepare_inputs body reads."""
    T.check_trainer(SimpleNamespace(_step=0, _buffered_inputs=None))
    with pytest.raises(T.TRLContractError, match=rf"trainer attribute {item}:.*Pin is 1\.17\.1"):
        T.check_trainer(SimpleNamespace(**{n: 0 for n in ("_step", "_buffered_inputs") if n != item}))


def test_a_built_grpo_trainer_has_the_trainer_attributes(tmp_path):
    """The attribute check on a tiny trainer built here (TRL sets both in GRPOTrainer.__init__)."""
    from tests.rl.test_trl_adapter_equivalence import tiny_loop_vllm_fake

    with pytest.MonkeyPatch.context() as mp:
        T.check_trainer(tiny_loop_vllm_fake(mp, tmp_path, steps=1))


def text_model(*layer_types):
    text = SimpleNamespace(layer_types=list(layer_types))
    return SimpleNamespace(config=SimpleNamespace(get_text_config=lambda: text))


def test_shared_prefix_refuses_linear_attention_layers():
    """packed_chunk passes a 4-D mask, and linear-attention layers get none (M2 rec 5): a hybrid model with
    shared_prefix raises; full attention only, or shared_prefix off, passes."""
    with pytest.raises(NotImplementedError, match="linear_attention"):
        T.check_attention(text_model("linear_attention", "full_attention"), shared_prefix=True)
    T.check_attention(text_model("full_attention", "full_attention"), shared_prefix=True)
    T.check_attention(text_model("linear_attention", "full_attention"), shared_prefix=False)
    T.check_attention(SimpleNamespace(config=SimpleNamespace(get_text_config=lambda: SimpleNamespace())), True)
    fwd = SimpleNamespace(shared_prefix=True)
    with pytest.raises(NotImplementedError, match="linear_attention"):
        T.FastForward._get_per_token_logps_and_entropies(fwd, text_model("linear_attention"), None, None, 1)

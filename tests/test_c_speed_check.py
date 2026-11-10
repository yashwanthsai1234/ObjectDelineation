"""c_speed_check --sync: the pass rules of G-SY1..G-SY3, and the whole check on a stateful fake engine (the box run
itself needs the real engine)."""

import re
from types import SimpleNamespace

import pytest

torch = pytest.importorskip("torch")
pytest.importorskip("trl")
from test_train_fast import _lora_bf16  # noqa: E402
from test_vllm_sync import FakeEngine, lora_request  # noqa: E402, F401  (the fixture)

from experiments.phase2 import c_speed_check as C  # noqa: E402
from objdel.train import fast, vllm_sync  # noqa: E402

MB = 1 << 20


OK = {"seconds": [0.05] * 20, "mem": [100 * MB] * 50, "equal": {"a": True, "b": True}, "rises": [10 * MB] * 3}


def test_sync_passes_a_fast_flat_series():
    v = C.sync_verdict(**OK)
    assert v["passed"] and v["median_s"] == 0.05 and v["peak_rise_b"] == 10 * MB


@pytest.mark.parametrize(
    "kw",
    [
        {"seconds": [0.3] * 20},  # G-SY2: median over 0.2 s
        {"mem": [100 * MB] * 49 + [102 * MB]},  # G-SY3: a leak of 2 MB
        {"equal": {"a": True, "b": False}},  # G-SY1: one push differs from TRL's
        {"rises": [10 * MB, 100 * MB]},  # G-SY1: max_memory_allocated rose by 100 MB or more
        {"seconds": [0.05] * 19 + [0.5]},  # G-SY3: seconds not stable
    ],
)
def test_sync_fails_each_rule(kw):
    assert not C.sync_verdict(**OK | kw)["passed"]


def test_hash_params_sees_one_changed_value():
    a = {"w": torch.ones(4, dtype=torch.bfloat16)}
    b = {"w": torch.ones(4, dtype=torch.bfloat16)}
    assert C.hash_params(a.items()) == C.hash_params(b.items())
    b["w"][2] = 2
    assert C.hash_params(a.items()) != C.hash_params(b.items())


def test_parse_takes_sync_and_the_train_rl_flags():
    a = C.parse(["--sync", "--model", "qwen3vl2b", "--ckpt", "c", "--out", "o", "--train", "t", "--steps", "1"])
    assert a.sync and a.fast_sync is False


# the whole check on a stateful fake engine -------------------------------------------------------------------------

LAYER = re.compile(r"layers\.(\d+)\.(self_attn|mlp)\.(\w+)\.lora_([AB])\.weight$")
VLLM_OF = {"q_proj": ("qkv_proj", 0), "k_proj": ("qkv_proj", 1), "v_proj": ("qkv_proj", 2), "o_proj": ("o_proj", 0),
           "gate_proj": ("gate_up_proj", 0), "up_proj": ("gate_up_proj", 1), "down_proj": ("down_proj", 0)}  # fmt: skip


class StateEngine(FakeEngine):
    """Keeps the weights `load_weights` writes; level-2 sleep scrambles them and empties the LoRA slots; `generate`
    is a function of the weights; `add_lora` packs the adapter folder the way vLLM does (bf16, scaling in lora_b)."""

    def __init__(self, log, model, fold=2.0):
        super().__init__(log)
        self.params = {n: torch.zeros_like(t) for n, t in fast.lora_merged_params(model)}
        self.asleep, self.fold, self.adds, self.modules = False, fold, 0, {}
        runner = self.llm_engine.model_executor.driver_worker.model_runner
        runner.model.named_parameters = lambda: self.params.items()
        runner.lora_manager = SimpleNamespace(_adapter_manager=SimpleNamespace(
            modules=self.modules, lora_index_to_id=[1]))  # fmt: skip
        self.llm_engine.add_lora = self.add_lora

    def load_weights(self, weights):
        assert not self.asleep, "load_weights into slept memory"
        for n, t in weights:
            self.params[n].copy_(t)
        super().load_weights([])

    def wake_up(self, tags=None):
        self.asleep = self.asleep and "weights" not in tags
        super().wake_up(tags)

    def sleep(self, level=1):
        self.asleep = True
        self.modules.clear()
        for p in self.params.values():
            p.fill_(7.0)

    def reset_prefix_cache(self):
        super().reset_prefix_cache()

    def add_lora(self, req):
        from safetensors.torch import load_file

        self.adds += 1
        self.modules.clear()
        sd = load_file(str(req.lora_path) + "/adapter_model.safetensors")
        for k, t in sd.items():
            i, block, proj, ab = LAYER.search(k).groups()
            name, part = VLLM_OF[proj]
            mod = self.modules.setdefault(f"language_model.model.layers.{i}.{block}.{name}",
                                          SimpleNamespace(lora_a_stacked={}, lora_b_stacked={}))  # fmt: skip
            (mod.lora_a_stacked if ab == "A" else mod.lora_b_stacked)[part] = t if ab == "A" else t * self.fold
        for mod in self.modules.values():  # tuples of (slots, 1, ...) per slice, as vLLM stacks them
            mod.lora_a_stacked = tuple(mod.lora_a_stacked[j][None, None] for j in sorted(mod.lora_a_stacked))
            mod.lora_b_stacked = tuple(mod.lora_b_stacked[j][None, None] for j in sorted(mod.lora_b_stacked))

    def generate(self, prompts, sampling_params=None, **_):
        assert not self.asleep
        digest = int(sum(float(p.detach().float().sum()) for p in self.params.values()) * 1000)
        return [SimpleNamespace(outputs=[SimpleNamespace(token_ids=[digest % 1000 + i, i])])
                for i, _ in enumerate(prompts)]  # fmt: skip


class FakeVG:
    """TRL's VLLMGeneration as far as the check reads it: the per-parameter sync_weights loop (class level)."""

    def __init__(self, model, engine, sleep_mode):
        self.llm, self.model, self.enable_sleep_mode = engine, model, sleep_mode
        self._llm_weights_sleeping = sleep_mode
        if sleep_mode:
            engine.sleep(2)

    def sync_weights(self):
        if self.enable_sleep_mode:
            self.llm.wake_up(tags=["weights"])
            self._llm_weights_sleeping = False
        for n, t in fast.lora_merged_params(self.model):
            self.llm.load_weights([(n, t)])
        self.llm.reset_prefix_cache()


@pytest.fixture
def box(monkeypatch, lora_request):  # noqa: F811
    """A tiny LoRA model, a stateful fake engine and fake CUDA memory counters; `box(sleep_mode)` -> (trainer, log)."""
    import sys

    monkeypatch.setattr("trl.generation.vllm_generation.empty_cache", lambda: None)
    sys.modules["vllm"].SamplingParams = lambda **kw: SimpleNamespace(**kw)
    clock = iter(x * 0.01 for x in range(10**6))  # every timed push takes 0.01 s: the host's noise does not decide
    monkeypatch.setattr("time.perf_counter", lambda: next(clock))
    state = {"rise": 0}
    monkeypatch.setattr(torch.cuda, "memory_allocated", lambda *a: 100 * MB)
    monkeypatch.setattr(torch.cuda, "reset_peak_memory_stats", lambda *a: None)
    monkeypatch.setattr(torch.cuda, "max_memory_allocated", lambda *a: 100 * MB + state["rise"])

    def make(sleep_mode, fold=2.0):
        m = _lora_bf16()
        eng = StateEngine([], m, fold)
        vg = FakeVG(m, eng, sleep_mode)
        tok = SimpleNamespace(encode=lambda text: [ord(c) % 50 for c in text])
        return SimpleNamespace(vllm_generation=vg, processing_class=SimpleNamespace(tokenizer=tok), state=state)

    return make


def install_merged(t):
    vg = t.vllm_generation
    vg.sync_weights = lambda: vllm_sync.push_weights(vg, vg.model)


def install_f(t, folder):
    vg = t.vllm_generation
    vg.sync_weights = lambda: vllm_sync.sync_lora(vg, vg.model, folder)


@pytest.mark.parametrize("sleep_mode", [False, True])
def test_merged_push_check_passes_awake_and_after_sleep_wake_cycles(box, sleep_mode):
    t = box(sleep_mode)
    install_merged(t)
    out = C.run_sync(t, lora=False)
    assert out["passed"], out
    assert out["touched_tensors"] == sum(1 for m in t.vllm_generation.model.modules() if hasattr(m, "lora_A"))
    assert out["touched_s"] is not None and isinstance(out["readd_delta_push"], bool)
    assert 16 == len(C.prompts(t))


def test_a_push_that_changes_nothing_fails_g_sy1(box):
    t = box(False)
    t.vllm_generation.sync_weights = lambda: vllm_sync.PushStats(0, 0, 0.0)
    out = C.run_sync(t, lora=False)
    assert not out["G-SY1"] and not out["passed"]


def test_a_push_that_drops_the_delta_fails_g_sy1(box):
    t = box(False)
    vg = t.vllm_generation
    vg.sync_weights = lambda: vllm_sync.push_weights(vg, vg.model, values="base")
    assert not C.run_sync(t, lora=False)["G-SY1"]


def test_a_peak_memory_rise_of_100_mb_fails_g_sy1(box):
    t = box(False)
    install_merged(t)
    t.state["rise"] = 100 * MB
    out = C.run_sync(t, lora=False)
    assert not out["G-SY1"] and out["peak_rise_b"] == 100 * MB


@pytest.mark.parametrize("sleep_mode", [False, True])
def test_lora_check_passes_packed_slices_and_one_add_lora_per_sync(box, tmp_path, sleep_mode):
    t = box(sleep_mode)
    install_f(t, tmp_path / "lora")
    out = C.run_sync(t, lora=True)
    assert out["passed"], out
    assert out["add_lora_per_sync"] == [1] * C.PERTURBED
    assert out["trl_loop_s"] is None  # TRL's loop would write W + d into the base the adapter rides on


def test_lora_check_fails_a_lora_b_slice_without_the_folded_scaling(box, tmp_path):
    t = box(False, fold=1.0)
    install_f(t, tmp_path / "lora")
    assert not C.run_sync(t, lora=True)["G-SY1"]


def test_lora_check_fails_two_add_lora_calls_in_one_sync(box, tmp_path):
    t = box(False)
    vg = t.vllm_generation
    vg.sync_weights = lambda: [vllm_sync.sync_lora(vg, vg.model, tmp_path / "lora") for _ in range(2)][-1]
    out = C.run_sync(t, lora=True)
    assert out["add_lora_per_sync"] == [2] * C.PERTURBED and not out["G-SY1"]


def test_lora_check_fails_a_sync_that_skips_the_base_push_after_sleep(box, tmp_path):
    t = box(True)
    vg = t.vllm_generation

    def sync():
        vg.enable_sleep_mode = False  # sync_lora then pushes no base W into the engine the sleep scrambled
        try:
            return vllm_sync.sync_lora(vg, vg.model, tmp_path / "lora")
        finally:
            vg.enable_sleep_mode = True

    vg.sync_weights = sync
    assert not C.run_sync(t, lora=True)["G-SY1"]


def test_packed_mismatches_names_each_bad_tensor_and_never_passes_vacuously():
    m = _lora_bf16()
    from peft import get_peft_model_state_dict

    slices = {}
    for k, t in get_peft_model_state_dict(m).items():
        i, block, proj, ab = LAYER.search(k).groups()
        name, part = VLLM_OF[proj]
        a, b = slices.get((int(i), name, part), (None, None))
        t = t.to(torch.bfloat16) if ab == "A" else 2 * t.to(torch.bfloat16)
        slices[(int(i), name, part)] = (t, b) if ab == "A" else (a, t)
    assert C.packed_mismatches(m, slices) == []
    key = next(iter(slices))
    slices[key] = (slices[key][0], slices[key][1] + 1)
    assert len(C.packed_mismatches(m, slices)) == 1
    del slices[key]
    assert len(C.packed_mismatches(m, slices)) == 2  # a missing slice: both its tensors
    assert C.packed_mismatches(m, {})  # nothing read from the engine: not a pass

"""c_speed_check: --sync's pass rules of G-SY1..G-SY3 and the whole check on a stateful fake engine (the box run
itself needs the real engine); --pixels's five batches and pixel check (U-CB1, G-EQ1)."""

import re
from types import SimpleNamespace

import numpy as np
import PIL.Image
import pytest

torch = pytest.importorskip("torch")
pytest.importorskip("trl")
from datasets import Dataset, Image  # noqa: E402

from objdel.core.hashing import pixels_sha256  # noqa: E402
from objdel.rl import rows as dataset  # noqa: E402
from objdel.rl import trl_adapter as fast  # noqa: E402
from objdel.rl import vllm_sync  # noqa: E402
from tests.experiments.test_trl_adapter_equivalence import _lora_bf16  # noqa: E402
from tests.rl.test_vllm_sync import FakeEngine, lora_request  # noqa: E402, F401  (the fixture)
from tools.studies import c_speed_check as C  # noqa: E402

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
    flag, a = C.parse(["--sync", "--model", "qwen3vl2b", "--ckpt", "c", "--out", "o", "--train", "t", "--steps", "1"])
    assert flag.sync and not flag.pixels and a.fast_sync is False


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


# --- C1: the five batches and the pixel check of `--pixels` (U-CB1, G-EQ1) ----------------------------------------

G = 8


def bundle(n=170, side=32):
    """A fixture rows set: n train rows with one random RGB image each; (OpenOnce over them, key -> pixel SHA-256)."""
    imgs = [PIL.Image.fromarray(np.random.default_rng(i).integers(0, 256, (side, side, 3), dtype=np.uint8))
            for i in range(n)]
    rows = Dataset.from_list([{"key": f"item{i}", "image": im} for i, im in enumerate(imgs)])
    rows = rows.cast_column("image", Image())
    return dataset.OpenOnce(rows), {f"item{i}": pixels_sha256(np.asarray(im)) for i, im in enumerate(imgs)}


def fake_generate(seen):
    """The trainer's job without an engine: a record in the --record-batches format for the 256 rows it is given."""
    def generate(inputs, step):
        seen.append((step, inputs))
        n = len(inputs)
        ids = {k: torch.zeros(n, 4, dtype=torch.long) for k in ("prompt_ids", "completion_ids")}
        mask = {k: torch.ones(n, 4, dtype=torch.long) for k in ("prompt_mask", "completion_mask")}
        first = {id(r["image"]): i for i, r in reversed(list(enumerate(inputs)))}
        meta = {"image_id": torch.tensor([sorted(first.values()).index(first[id(r["image"])]) for r in inputs]),
                "group": torch.arange(n) // G}
        return {**ids, **mask, "image_grid_thw": torch.tensor([[1, 2, 2]] * n), "num_images": torch.ones(n).long(),
                "meta": meta}
    return generate


def test_the_draw_is_seeded_distinct_train_rows_in_five_batches_of_32():
    a = C.draw_batches(8320)
    assert a.shape == (5, 32) and len(set(a.ravel())) == 160 and a.max() < 8320
    assert (C.draw_batches(8320) == a).all() and (C.draw_batches(8320, seed=1) != a).any()
    with pytest.raises(ValueError, match="160"):
        C.draw_batches(159)


@pytest.mark.parametrize("path", ["data/eval_views/lidc-val-shifted", "/d/lidc-val-edge", "x/val", "eval_test"])
def test_a_val_view_path_raises(path):
    with pytest.raises(ValueError, match="val"):
        C.refuse_val(path)
    from pathlib import Path

    for ok in ("/workspace/data/lidc/rows/lidc-train-b1", "data/prepared", Path("runs/interval/x")):
        C.refuse_val(ok)


def test_main_refuses_a_val_root_before_building_anything():
    with pytest.raises(ValueError, match="val"):
        C.main(["--pixels", "--model", "m", "--ckpt", "c", "--out", "o", "--steps", "20",
                "--train", "data/eval_views/lidc-val-shifted"])


def test_the_advantages_are_seeded_normal_values_centred_per_group_of_8():
    a = C.centred_advantages(5)
    assert a.shape == (5, 256) and torch.allclose(a.view(5, 32, G).sum(-1), torch.zeros(5, 32), atol=1e-5)
    assert torch.equal(a, C.centred_advantages(5)) and (a[0] != a[1]).any()
    raw = torch.randn(5 * 256, generator=torch.Generator().manual_seed(0)).view(5, 32, G)
    assert torch.allclose(a.view(5, 32, G), raw - raw.mean(-1, keepdim=True))  # seed 0, as the plan says


def test_the_batch_maker_writes_exactly_five_files_that_weights_only_loads(tmp_path):
    """U-CB1: b0..b4 in the --record-batches format plus advantages, image keys and SHA-256 values."""
    rows, sha = bundle()
    seen = []
    paths = C.make_batches(fake_generate(seen), rows, len(rows), tmp_path / "batches")
    assert sorted(p.name for p in (tmp_path / "batches").iterdir()) == [f"b{i}.pt" for i in range(5)]
    assert [p.name for p in paths] == [f"b{i}.pt" for i in range(5)]
    assert [s for s, _ in seen] == [0, 1, 2, 3, 4] and all(len(x) == 256 for _, x in seen)
    for _, inputs in seen:  # groups of 8 identical rows share one image object, 32 distinct images a batch
        assert all(inputs[i]["image"] is inputs[i // G * G]["image"] for i in range(256))
        assert len({id(r["image"]) for r in inputs}) == 32
    chosen = [r["key"] for _, x in seen for r in x[::G]]
    assert len(set(chosen)) == 160  # five batches of 32 distinct rows, no row twice
    for i, p in enumerate(paths):
        rec = torch.load(p, weights_only=True)
        assert {"prompt_ids", "completion_ids", "meta", "advantages", "image_key", "image_sha256"} <= set(rec)
        assert rec["image_key"] == chosen[i * 32:(i + 1) * 32]
        assert rec["image_sha256"] == [sha[k] for k in rec["image_key"]]  # the fixture bundle's own hashes
        assert rec["advantages"].shape == (256,)
        assert torch.allclose(rec["advantages"].view(32, G).sum(-1), torch.zeros(32), atol=1e-5)
        assert torch.equal(rec["advantages"], C.centred_advantages(5)[i])


def test_the_batch_maker_leaves_a_non_empty_folder_alone(tmp_path):
    rows, _ = bundle()
    (tmp_path / "b").mkdir()
    (tmp_path / "b/old").write_text("x")
    with pytest.raises(FileExistsError):
        C.make_batches(fake_generate([]), rows, len(rows), tmp_path / "b")


def test_pixel_compare_counts_mismatches_and_names_the_first():
    """G-EQ1's comparison on CPU: the torch backend equals the HF processor on every view; a planted processor that
    changes one image shows up as one mismatch at its index."""
    from transformers import Qwen2VLImageProcessor

    ip = Qwen2VLImageProcessor(patch_size=16, temporal_patch_size=2, merge_size=2)
    imgs = [PIL.Image.fromarray(np.random.default_rng(i).integers(0, 256, (64, 64, 3), dtype=np.uint8))
            for i in range(70)]
    ok = C.compare_pixels(imgs, ip, "cpu", "torch", chunk=32)
    assert (ok["views"], ok["mismatches"], ok["first"]) == (70, 0, None) and ok["median_ms"] >= 0

    class Off(Qwen2VLImageProcessor):
        def __call__(self, images, *a, **k):
            out = super().__call__(images, *a, **k)
            if len(images) > 40 or len(images) == 6:  # chunk 32 -> 32, 32, 6: only the last chunk is touched
                out["pixel_values"][0, 0] += 1
            return out

    bad = Off(patch_size=16, temporal_patch_size=2, merge_size=2)
    r = C.compare_pixels(imgs, bad, "cpu", "torch", chunk=32)
    assert (r["mismatches"], r["first"]) == (1, 64)  # the chunk starting at view 64: its first image differs


# --- C3: --encoder (G-EN0, G-EN1, G-IM2); the CUDA memory measurement itself runs on the box ------------------------


def test_floors_read_only_the_answer_tokens():
    a = torch.tensor([[0.0, 1.0, 5.0], [2.0, 2.0, 2.0]])
    b = torch.tensor([[0.0, 1.5, 0.0], [2.0, 1.0, 2.0]])
    mask = torch.tensor([[1, 1, 0], [1, 1, 1]])  # the 5.0 difference is off the answer
    fmax, p99 = C.floors(a, b, mask)
    assert fmax == 1.0 and p99 == pytest.approx(np.quantile([0, 0.5, 0, 1, 0], 0.99))
    assert C.floors(a, b, torch.zeros_like(mask)) == (0.0, 0.0)


@pytest.mark.parametrize("diff, floor, ok", [
    ((0.0, 0.0), (0.0, 0.0), True), ((1e-9, 0.0), (0.0, 0.0), False),  # floor 0: equality only
    ((0.2, 0.02), (0.1, 0.01), True), ((0.21, 0.02), (0.1, 0.01), False), ((0.2, 0.021), (0.1, 0.01), False)])
def test_within_floor_is_two_times_both_floors(diff, floor, ok):
    assert C.within_floor(diff, floor) is ok


def test_feature_diff_sees_one_changed_deepstack_value():
    def entry(seed):
        g = torch.Generator().manual_seed(seed)
        return torch.randn(4, 8, generator=g), [torch.randn(4, 8, generator=g)]

    a = [entry(0), entry(1)]
    b = [(e.clone(), [x.clone() for x in lv]) for e, lv in a]
    assert C.feature_diff(a, b) == (True, 0.0)
    b[1][1][0][2, 3] += 0.25
    eq, mx = C.feature_diff(a, b)
    assert not eq and mx == pytest.approx(0.25)


def test_batch_tensors_type_only_real_prompt_image_tokens():
    from tests.experiments.test_trl_adapter_equivalence import IMG

    rec = {"prompt_ids": torch.tensor([[IMG, 5, IMG, IMG]]), "prompt_mask": torch.tensor([[0, 1, 1, 1]]),
           "completion_ids": torch.tensor([[IMG, 7]]), "completion_mask": torch.tensor([[1, 0]]),
           "meta": {f: torch.tensor([1]) for f in ("prompt_len", "completion_len", "image_id", "group", "layout",
                                                    "rope_bad", "pixel_rows")}}  # fmt: skip
    ids, mask, T, mm, meta = C.batch_tensors(rec, IMG)
    assert ids.tolist() == [[IMG, 5, IMG, IMG, IMG, 7]] and mask.tolist() == [[0, 1, 1, 1, 1, 0]] and T == 2
    assert mm.tolist() == [[0, 0, 1, 1, 0, 0]]  # the masked pad position and the answer's image id are text
    assert meta.image_id.dtype == np.int32 and meta.pixel_rows.tolist() == [1]


def _tiny_check(shared_prefix):
    """The tiny Qwen3-VL (LoRA, frozen tower) behind a trainer with VisionCache and RopeCache, on 4 groups of 4."""
    from peft import get_peft_model

    from objdel.policy import caches
    from objdel.policy import models as configs
    from tests.experiments.test_encoder_once import BUCKETS, distinct, four_groups
    from tests.experiments.test_trl_adapter_equivalence import Fast, meta_for, tiny_model, trainer

    m = get_peft_model(tiny_model(), configs.lora("qwen3-vl-2b"))
    g = torch.Generator().manual_seed(1)
    for n, p in m.named_parameters():
        if "lora_" in n:
            p.data = torch.randn(p.shape, generator=g) * 0.05
    t = trainer(Fast, buckets=BUCKETS, vision=caches.VisionCache.install(m), rope=caches.RopeCache.install(m),
                shared_prefix=shared_prefix)  # fmt: skip
    ids, mask, T, kw = four_groups()
    adv = C.centred_advantages(1)[0, :16]
    return t, m, (ids, mask, T, kw, meta_for(mask, T, 4), distinct(kw), adv)


@pytest.mark.parametrize("shared_prefix", [False, True])
def test_grouping_check_passes_on_the_tiny_model_and_writes_what_encode_check_reads(shared_prefix, tmp_path):
    """G-EN0 + G-EN1 on CPU: the encoder-once path equals the pixel path fed the one-call features (reference pass,
    every micro-batch, LoRA grads), one tower call; the record's bits_equal is a JSON boolean train_rl accepts."""
    import json

    from objdel.policy.callbacks import ExpectStep1
    from objdel.rl.run import save_json

    t, m, args = _tiny_check(shared_prefix)
    r = C.grouping_check(t, m, *args, chunk=4)
    assert r["deterministic"] and r["G-EN1"]["passed"], r
    assert r["G-EN1"]["tower_calls"] == 1 and t.encoder is not None and t.vision is None
    assert type(r["bits_equal"]) is bool and r["max_abs"] >= 0 and r["floor_max"] >= r["floor_p99"] >= 0
    assert r["bits_equal"] == (r["max_abs"] == 0)
    save_json(tmp_path / "encoder.json", r, ("passed",))
    from objdel.policy.callbacks import STEP1_CLOSE, STEP1_EXACT

    log = {"log": [{"step": 1} | dict.fromkeys(STEP1_EXACT + STEP1_CLOSE, 1.0)]}
    (tmp_path / "rl.json").write_text(json.dumps(log))
    got = ExpectStep1.from_files(tmp_path / "rl.json", tmp_path / "encoder.json", tmp_path)
    assert got.bits_equal is r["bits_equal"]


def test_grouping_check_fails_an_encoder_path_that_gives_rows_the_wrong_image(monkeypatch):
    from objdel.policy import caches

    t, m, args = _tiny_check(False)
    real = caches.EncoderCache.outputs
    monkeypatch.setattr(caches.EncoderCache, "outputs", lambda self, ids: real(self, [(int(i) + 1) % 4 for i in ids]))
    r = C.grouping_check(t, m, *args, chunk=4)
    assert not r["G-EN1"]["equal"] and not r["G-EN1"]["passed"]


def test_grouping_check_fails_a_second_tower_call(monkeypatch):
    from objdel.policy import caches

    t, m, args = _tiny_check(False)
    real = caches.EncoderCache.encode
    monkeypatch.setattr(caches.EncoderCache, "encode", lambda self, pv, grid: [real(self, pv, grid) for _ in "ab"][0])
    r = C.grouping_check(t, m, *args, chunk=4)
    assert r["G-EN1"]["tower_calls"] == 2 and r["G-EN1"]["equal"] and not r["G-EN1"]["passed"]


MEM = {"memo_b": 200 * MB, "distinct_b": 200 * MB, "call_b": 300 * MB, "out_b": 201 * MB, "full": 0,
       "peak_b": 20 << 30}  # fmt: skip


def test_memory_verdict_passes_the_plan_budget():
    assert C.memory_verdict(**MEM)["passed"]


@pytest.mark.parametrize("kw", [
    {"memo_b": 221 * MB},  # the memo holds more than the distinct images + 10 %
    {"call_b": 400 * MB},  # the pass-2 call raised max_memory_allocated by 400 MB
    {"call_b": 302 * MB, "out_b": 201 * MB},  # 101 MB beyond its own output
    {"full": 1},  # a 256-row patch tensor lives
    {"peak_b": (30 << 30) + 1}])  # fmt: skip
def test_memory_verdict_fails_each_rule(kw):
    assert not C.memory_verdict(**MEM | kw)["passed"]


def test_encoder_needs_the_pixels_batch(tmp_path):
    with pytest.raises(FileNotFoundError, match="--pixels"):
        C.encoder_batch(None, None, 170, tmp_path)


def test_parse_takes_encoder_and_main_refuses_a_val_root():
    flag, _ = C.parse(["--encoder", "--model", "qwen3vl2b", "--ckpt", "c", "--out", "o", "--train", "t",
                        "--steps", "1"])
    assert flag.encoder and not flag.pixels
    with pytest.raises(ValueError, match="val"):
        C.main(["--encoder", "--model", "m", "--ckpt", "c", "--out", "o", "--steps", "20",
                "--train", "data/eval_views/lidc-val-shifted"])


def test_main_runs_encoder_and_writes_encoder_json(tmp_path, monkeypatch):
    """main's --encoder writes RUN/c_check/encoder.json (the file --encode-check names) and exits on `passed`."""
    import json

    from objdel.rl import build

    monkeypatch.setattr(dataset, "grpo_rows", lambda *a: [])
    monkeypatch.setattr(dataset, "OpenOnce", lambda rows: rows)
    monkeypatch.setattr(build, "build_trainer", lambda *a: SimpleNamespace(trainer=None))
    for passed, code in ((True, 0), (False, 1)):
        monkeypatch.setattr(C, "run_encoder", lambda *a, p=passed: {"bits_equal": True, "passed": p})
        assert C.main(["--encoder", "--model", "qwen3-vl-2b", "--ckpt", "c", "--out", str(tmp_path), "--steps", "20",
                       "--train", "t"]) == code  # fmt: skip
        assert json.loads((tmp_path / "c_check" / "encoder.json").read_text())["passed"] is passed

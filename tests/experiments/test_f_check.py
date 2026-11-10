"""B2.2 F step 5, B2.4 rule 0 (U76, Mac twin of G05): the F check record, its pass rules and the `train_rl --f-check`
flag, with fake probes; then the engine half with a fake vLLM engine (no vLLM on the Mac): the comparator's order,
prompts, timing and restore, the packed-layout equality, the after-sync probe, way (b) scoring and train_rl's wiring."""
import json
from types import SimpleNamespace

import pytest

from experiments.phase2 import train_rl
from objdel.rl import f_check as fc

BASE = ["--model", "qwen3-vl-2b", "--ckpt", "c", "--train", "t", "--out", "o", "--steps", "20"]


def lora_fixed(**change):
    return {"comparator": fc.comparator(0.0176, 0.0096, 0.0096), "base_hash_ok": True, "lora_over_base": 1.05,
            "lora_equal": True} | change


def run(tmp_path, preset="lora-seqmask", steps=15, nan_at=None, loads_per_step=1, **kw):
    """Drive the callback like TRL: one on_log and one on_step_end per optimizer step."""
    loads = fc.count_calls(lambda: None)
    probes = ({"engine_dtype_ok": lambda: kw.pop("dtype_ok", True)} if preset == "fp16-seqmask" else {})
    cb = fc.FCheck(tmp_path, preset, probes, loads=loads if preset == "lora-seqmask" else None,
                   fixed=kw.pop("fixed", lora_fixed() if preset == "lora-seqmask" else {}))
    for step in range(1, steps + 1):
        for _ in range(loads_per_step):
            loads()
        cb.on_log(None, SimpleNamespace(global_step=step), None, {"loss": float("nan") if step == nan_at else 0.1})
        cb.on_step_end(None, SimpleNamespace(global_step=step), None)
    return json.loads((tmp_path / "f_check.json").read_text())


def test_all_checks_pass_writes_the_record_after_step_15_and_the_run_goes_on(tmp_path):
    rec = run(tmp_path, steps=20)
    assert rec["passed"] is True and rec["preset"] == "lora-seqmask" and rec["nan"] == 0 and rec["load_count"] == 1
    fc.check_record(rec)


@pytest.mark.parametrize("fixed,flag", [
    ({"lora_equal": False}, "lora_equal"), ({"base_hash_ok": False}, "base_hash_ok"),
    ({"lora_over_base": 1.151}, "lora_over_base"),
    ({"comparator": fc.comparator(0.0176, 0.0140, 0.0096)}, "comparator")])
def test_a_failed_check_writes_passed_false_then_raises(tmp_path, fixed, flag):
    with pytest.raises(RuntimeError, match="F check failed"):
        run(tmp_path, fixed=lora_fixed(**fixed))
    rec = json.loads((tmp_path / "f_check.json").read_text())
    assert rec["passed"] is False and not (rec[flag]["ok"] if flag == "comparator" else rec[flag] is True)


def test_nan_at_step_9_raises_there_and_two_loads_in_one_step_fail(tmp_path):
    with pytest.raises(RuntimeError):
        run(tmp_path, nan_at=9)
    assert json.loads((tmp_path / "f_check.json").read_text())["nan"] == 1
    with pytest.raises(RuntimeError):
        run(tmp_path / "two", loads_per_step=2)
    assert json.loads((tmp_path / "two/f_check.json").read_text())["load_count"] == 2


def test_fp16_record_is_nan_and_the_engine_dtype(tmp_path):
    rec = run(tmp_path / "ok", "fp16-seqmask")
    assert rec == {"preset": "fp16-seqmask", "nan": 0, "engine_dtype_ok": True, "passed": True}
    with pytest.raises(RuntimeError):
        run(tmp_path / "bad", "fp16-seqmask", dtype_ok=False)
    assert json.loads((tmp_path / "bad/f_check.json").read_text())["passed"] is False


@pytest.mark.parametrize("f0,merged,lora,ok,informative", [
    (0.0096, 0.0176, 0.0096, True, True),
    (0.0096, 0.0176, 0.0140, False, True),  # ratio 0.55
    (0.0096, 0.0176, 0.0130, True, True),  # ratio 0.425
    (0.0096, 0.0115, 0.0115, True, False)])  # merged - f0 = 0.0019: cannot tell F from no F
def test_comparator_subtracts_the_floor(f0, merged, lora, ok, informative):
    c = fc.comparator(merged, lora, f0)
    assert (c["ok"], c["informative"]) == (ok, informative) and (c["ratio"] is None) == (not informative)


def test_timing_ratio_is_the_median_and_1_15_passes():
    assert fc.lora_over_base([1.15] * 6, [1.0] * 6) == 1.15
    assert fc.passes({"preset": "lora-seqmask", "nan": 0, "load_count": 1, "base_hash_ok": True, "lora_equal": True,
                      "comparator": {"ok": True}, "lora_over_base": 1.15})
    assert fc.lora_over_base([1.0, 1.0, 9.0], [1.0, 1.0, 1.0]) == 1.0


def test_check_record_never_trusts_passed_alone():
    good = {"preset": "lora-seqmask", "nan": 0, "load_count": 1, "base_hash_ok": True, "lora_equal": True,
            "comparator": {"ok": True}, "lora_over_base": 1.0, "passed": True}
    fc.check_record(good)
    with pytest.raises(AssertionError):
        fc.check_record(good | {"lora_equal": False})  # one flag false, `passed` true
    with pytest.raises(AssertionError):
        fc.check_record(good | {"lora_equal": False, "passed": False})  # consistent, but a failure


def test_f_check_flag_is_on_by_default_under_fp16_and_an_argparse_error_under_bf16(capsys):
    assert train_rl.parse([*BASE, "--numerics", "fp16-seqmask"]).f_check is True
    assert train_rl.parse([*BASE, "--numerics", "fp16-seqmask", "--no-f-check"]).f_check is False
    assert train_rl.parse(BASE).f_check is False
    with pytest.raises(SystemExit):
        train_rl.parse([*BASE, "--f-check"])
    assert "--f-check" in capsys.readouterr().err


# U76, the engine half (fakes; no vLLM on the Mac): the comparator's order, prompts and restore, the packed-layout
# equality and the after-sync probe.
torch = pytest.importorskip("torch")
pytest.importorskip("trl")  # before the vllm stubs: trl's import probes vllm's spec
from objdel.rl import run as rl_run  # noqa: E402
from objdel.rl import trl_adapter as fast  # noqa: E402
from objdel.rl import vllm_sync  # noqa: E402

IMG, START, END = 100, 102, 103  # tiny_model's image pad and vision start / end ids


@pytest.fixture
def lora_request(monkeypatch):
    """vllm.lora.request.LoRARequest stand-in."""
    import sys
    import types

    class LoRARequest(SimpleNamespace):
        def __init__(self, lora_name, lora_int_id, lora_path, load_inplace=False):
            super().__init__(lora_name=lora_name, lora_int_id=lora_int_id, lora_path=lora_path,
                             load_inplace=load_inplace)  # fmt: skip

    mods = {k: types.ModuleType(k) for k in ("vllm", "vllm.lora", "vllm.lora.request")}
    mods["vllm.lora.request"].LoRARequest = LoRARequest
    for k, m in mods.items():
        monkeypatch.setitem(sys.modules, k, m)
    return LoRARequest


class Engine:
    """vg.llm: stores pushed tensors by name (load_weights), logs every call, generates tagged answers; vLLM's image
    expansion multiplies each pad by `expand`."""

    def __init__(self, log, base: dict, expand: int = 2):
        self.log, self.weights, self.expand = log, {n: t.clone() for n, t in base.items()}, expand
        model = SimpleNamespace(load_weights=self.load_weights, named_parameters=lambda: iter(self.weights.items()))
        runner = SimpleNamespace(model=model)
        self.llm_engine = SimpleNamespace(model_executor=SimpleNamespace(driver_worker=SimpleNamespace(
            model_runner=runner)), add_lora=lambda r: log.append(("add_lora", r.lora_name, r.lora_int_id)),
            remove_lora=lambda i: log.append(("remove_lora", i)))  # fmt: skip

    def load_weights(self, pairs):
        for n, t in pairs:
            self.weights[n] = t.clone()
        self.log.append(("load_weights",))

    def wake_up(self, tags=None):
        self.log.append(("wake_up", tuple(tags)))

    def sleep(self, level):
        self.log.append(("sleep", level))

    def reset_prefix_cache(self):
        self.log.append(("reset_prefix_cache",))

    def generate(self, rows, sampling_params=None, lora_request=None, **_):
        seed, fixed = sampling_params
        lora = lora_request and lora_request.lora_name
        self.log.append(("generate", lora, seed, fixed, [r["prompt_token_ids"] for r in rows]))
        n = fc.FIXED_LEN if fixed else 3
        return [SimpleNamespace(prompt_token_ids=[x for t in r["prompt_token_ids"] for x in [t] * (
            self.expand if t == IMG else 1)], lora=lora, outputs=[SimpleNamespace(token_ids=[1] * n)]) for r in rows]





def comparator_run(tmp_path, monkeypatch, sleep_mode: bool, expand: int = 2, gap=None, skip=None):
    from tests.experiments.test_trl_adapter_equivalence import _lora_bf16

    m = _lora_bf16()
    m.save_pretrained(tmp_path / "b1")  # the B1 adapter: this adapter's values
    for n, p in m.named_parameters():  # the run's own fresh adapter differs from B1
        if "lora_B" in n:
            p.data.zero_()
    log = []
    monkeypatch.setattr("trl.generation.vllm_generation.empty_cache", lambda: log.append(("empty_cache",)))
    base = dict(fast.lora_merged_params(m, with_delta=False))
    vg = SimpleNamespace(llm=Engine(log, base, expand), model=m, _llm_weights_sleeping=sleep_mode,
                         enable_sleep_mode=sleep_mode, mode="colocate")  # fmt: skip
    prompts = [{"prompt_token_ids": [START, IMG, IMG, END], "multi_modal_data": {"image": "a"}}] * 4
    gap = gap or (lambda outs, adapter: 0.0096 if not adapter else 0.0100 if outs[0].lora == "b1" else 0.0176)
    clock = iter(range(10_000)).__next__  # every timed generate takes 1 tick
    fixed = fc.run_comparator(vg, m, prompts, IMG, tmp_path / "b1", tmp_path / "tmp", gap,
                              sampling=lambda seed, fixed=False: (seed, fixed), clock=clock, skip=skip)
    return fixed, log, vg, m, base


@pytest.mark.parametrize("sleep_mode", [True, False])
def test_comparator_order_restore_and_sleep_last(tmp_path, monkeypatch, lora_request, sleep_mode):
    """B2.2 F 4b.4/4b.5: (ii) on slot 2 first, the floor arm, the A B A B timing (1 warm-up per arm, seeds 0-2),
    the merged push last; the restore pushes the base, hashes it equal, removes slot 2, drops "b1", resets the
    prefix cache and, in mode B only, sleeps at level 2 last with TRL's flag set."""
    fixed, log, vg, m, base = comparator_run(tmp_path, monkeypatch, sleep_mode)
    gens = [(e[1], e[2], e[3]) for e in log if e[0] == "generate"]
    timing = [("b1", 0, True), (None, 0, True)] + [(lo, s, True) for s in fc.SEEDS for _ in range(2)
                                                   for lo in ("b1", None)]  # fmt: skip
    assert gens == [("b1", 0, False), (None, 0, False), *timing, (None, 0, False)]
    names = [e[0] for e in log]
    assert names.index("add_lora") < names.index("generate")
    assert log[names.index("add_lora")][1:] == ("b1", fc.B1_SLOT)
    last_gen = len(names) - 1 - names[::-1].index("generate")
    assert names[last_gen - 2:last_gen] == ["load_weights", "reset_prefix_cache"]  # the merged push, then arm (i)
    assert "load_weights" not in names[names.index("generate"):last_gen - 2]  # base-W arms all ran before it
    tail = names[last_gen + 1:]
    want = ["load_weights", "reset_prefix_cache", "remove_lora", "reset_prefix_cache"]
    want += ["sleep"] if sleep_mode else []
    assert tail == want and vg._llm_weights_sleeping is sleep_mode
    if sleep_mode:
        assert log[-1] == ("sleep", 2) and ("wake_up", ("kv_cache",)) in log
        assert names.index("wake_up") < names.index("add_lora")  # weights woken by sync_lora before the load
    else:
        assert "wake_up" not in names and "sleep" not in names
    assert fixed["base_hash_ok"] is True and fixed["prompt_expansion"] is True and fixed["fixed_len_ok"] is True
    assert fixed["lora_over_base"] == 1.0 and fixed["comparator"]["ok"] and fixed["comparator"]["f0"] == 0.0096
    assert list(m.peft_config) == ["default"] and m.active_adapter == "default"
    assert all(torch.equal(vg.llm.weights[n], t) for n, t in base.items())


def test_comparator_sends_one_pad_per_image_and_fails_on_another_expansion(tmp_path, monkeypatch, lora_request):
    _, log, *_ = comparator_run(tmp_path, monkeypatch, True)
    assert all(ids == [START, IMG, END] for e in log if e[0] == "generate" for ids in e[4])
    fixed, log, *_ = comparator_run(tmp_path / "x", monkeypatch, True, expand=3)
    assert fixed["prompt_expansion"] is False and "3 pads" in fixed["error"]
    rec = {"preset": "lora-seqmask", "nan": 0} | fixed
    assert fc.passes(rec) is False  # no KeyError on the missing values


def test_comparator_base_hash_catches_a_wrong_restore(tmp_path, monkeypatch, lora_request):
    push = vllm_sync.push_weights
    monkeypatch.setattr(vllm_sync, "push_weights", lambda vg, model, values="merged": push(vg, model, "merged"))
    fixed, *_ = comparator_run(tmp_path, monkeypatch, False)
    assert fixed["base_hash_ok"] is False


def served_modules(m, idx: int, scale_b: float = 2.0, adapter="default"):
    """vLLM's LoRA-wrapped modules of a PEFT model: packed qkv_proj / gate_up_proj, slot `idx` of 2."""
    from peft.tuners.lora import LoraLayer

    groups = {}
    for name, mod in m.named_modules():
        if isinstance(mod, LoraLayer):
            parent, leaf = name.rsplit(".", 1)
            packed, k = fc.PACKED.get(leaf, (leaf, 0))
            groups.setdefault(f"vllm.{parent.split('language_model.')[1]}.{packed}", {})[k] = mod
    out = {}
    for name, ks in groups.items():
        a, b = [], []
        for k in sorted(ks):
            A, B = ks[k].lora_A[adapter].weight, ks[k].lora_B[adapter].weight
            sa = torch.zeros(2, 1, *A.shape, dtype=torch.bfloat16)
            sb = torch.zeros(2, 1, *B.shape, dtype=torch.bfloat16)
            sa[idx, 0], sb[idx, 0] = A.to(torch.bfloat16), B.to(torch.bfloat16) * scale_b
            a.append(sa)
            b.append(sb)
        out[name] = SimpleNamespace(lora_a_stacked=tuple(a), lora_b_stacked=tuple(b))
    return out


def vg_with(mods, index_to_id):
    model = SimpleNamespace(named_modules=lambda: iter(mods.items()))
    runner = SimpleNamespace(model=model, lora_manager=SimpleNamespace(
        _adapter_manager=SimpleNamespace(lora_index_to_id=index_to_id)))
    return SimpleNamespace(llm=SimpleNamespace(llm_engine=SimpleNamespace(model_executor=SimpleNamespace(
        driver_worker=SimpleNamespace(model_runner=runner)))))


def test_lora_equal_reads_the_packed_slices_of_the_adapters_slot():
    from tests.experiments.test_trl_adapter_equivalence import _lora_bf16

    m = _lora_bf16()
    assert fc.lora_equal(vg_with(served_modules(m, 1), [None, 1]), m)
    assert not fc.lora_equal(vg_with(served_modules(m, 0), [None, 1]), m)  # the values sit in another slot
    assert not fc.lora_equal(vg_with(served_modules(m, 1, scale_b=1.0), [None, 1]), m)  # alpha / r not folded
    mods = served_modules(m, 1)
    mods.pop(next(k for k in mods if k.endswith("gate_up_proj")))
    assert not fc.lora_equal(vg_with(mods, [None, 1]), m)  # a PEFT module without a vLLM twin
    assert not fc.lora_equal(vg_with({}, [None, 1]), m)


def test_after_syncs_checks_right_after_syncs_1_and_2():
    seen = []
    vg = SimpleNamespace(sync_weights=lambda: seen.append("sync"))
    probe = fc.AfterSyncs(vg, lambda: seen.append("check") or True, 2)
    assert not probe.ok()
    for _ in range(4):
        vg.sync_weights()
    assert seen == ["sync", "check", "sync", "check", "sync", "sync"] and probe.ok()
    bad = fc.AfterSyncs(SimpleNamespace(sync_weights=lambda: None), lambda: False, 2)
    bad()
    bad()
    assert not bad.ok()


def test_lora_seqmask_f_check_is_on_by_default_and_needs_the_b1_adapter(capsys):
    with pytest.raises(SystemExit):
        train_rl.parse([*BASE, "--numerics", "lora-seqmask"])
    assert "--b1-adapter" in capsys.readouterr().err
    a = train_rl.parse([*BASE, "--numerics", "lora-seqmask", "--b1-adapter", "b1"])
    assert a.f_check is True and str(a.b1_adapter) == "b1"
    assert train_rl.parse([*BASE, "--numerics", "lora-seqmask", "--no-f-check"]).f_check is False


@pytest.mark.parametrize("sleep_mode", [True, False])
def test_train_rl_runs_the_full_lora_check(tmp_path, monkeypatch, sleep_mode):
    """(d): under lora-seqmask the callback carries every field (comparator values, the load count on add_lora
    wrapped after the comparator, lora_equal after syncs 1 and 2 (mode B) or 1-3 (mode A)); the record passes."""
    calls = {}

    def comparator(vg, model, prompts, pad, b1, tmp, gap, **kw):
        calls.update(prompts=prompts, pad=pad, b1=b1, gap=gap(["o"], True), skip=kw["skip"])
        vg.llm.llm_engine.add_lora("b1")  # the comparator's own load: not counted
        return lora_fixed() | {"prompt_expansion": True}

    monkeypatch.setattr(fc, "run_comparator", comparator)
    monkeypatch.setattr(fc, "lora_equal", lambda vg, model: True)
    engine = SimpleNamespace(add_lora=lambda req: None)
    vg = SimpleNamespace(llm=SimpleNamespace(llm_engine=engine), model=SimpleNamespace(
        config=SimpleNamespace(image_token_id=IMG)), enable_sleep_mode=sleep_mode, sync_weights=lambda: None)
    trainer = SimpleNamespace(vllm_generation=vg, num_generations=8, args=SimpleNamespace(generation_batch_size=16),
                              lora_dir=tmp_path / "lora", comparator_prompts=lambda rows, G: [r for r in rows
                                                                                               for _ in range(G)],
                              way_b_gap=lambda prompts, outs, adapter: 0.5)  # fmt: skip
    a = train_rl.parse([*BASE, "--numerics", "lora-seqmask", "--b1-adapter", "b1", "--out", str(tmp_path)])
    cb = rl_run.f_check_callback(a, trainer, [{"k": 0}, {"k": 1}, {"k": 2}])
    assert calls["prompts"] == [{"k": 0}] * 8 + [{"k": 1}] * 8 and calls["pad"] == IMG and calls["gap"] == 0.5
    assert calls["skip"] == "c/hf_revision.json missing: --ckpt is not a fetch_hf pull of rl-full-3ep"
    for step in range(1, 16):
        vg.sync_weights()  # TRL's sync at each step: sync_lora adds the adapter once
        engine.add_lora("policy")
        cb.on_log(None, SimpleNamespace(global_step=step), None, {"loss": 0.1})
        cb.on_step_end(None, SimpleNamespace(global_step=step), None)
    rec = json.loads((tmp_path / "f_check.json").read_text())
    assert set(fc.LORA_KEYS) <= rec.keys() and rec["load_count"] == 1 and rec["lora_equal"] is True
    fc.check_record(rec)


def test_train_rl_expansion_error_writes_the_record_and_raises(tmp_path, monkeypatch):
    monkeypatch.setattr(fc, "run_comparator", lambda *a, **k: {"prompt_expansion": False, "error": "3 pads"})
    vg = SimpleNamespace(llm=SimpleNamespace(llm_engine=SimpleNamespace(add_lora=None)),
                         model=SimpleNamespace(config=SimpleNamespace(image_token_id=IMG)), enable_sleep_mode=True)
    trainer = SimpleNamespace(vllm_generation=vg, num_generations=8, args=SimpleNamespace(generation_batch_size=8),
                              lora_dir=tmp_path / "lora", comparator_prompts=lambda rows, G: rows, way_b_gap=None)
    a = train_rl.parse([*BASE, "--numerics", "lora-seqmask", "--b1-adapter", "b1", "--out", str(tmp_path)])
    with pytest.raises(RuntimeError, match="F check failed"):
        rl_run.f_check_callback(a, trainer, [{"k": 0}])
    rec = json.loads((tmp_path / "f_check.json").read_text())
    assert rec["prompt_expansion"] is False and rec["passed"] is False


def test_way_b_gap_scores_each_prompt_group_with_the_adapter_on_or_off(monkeypatch):
    """Mean |vLLM logp - trainer logp| over the answer tokens only; one forward per prompt group; the image prepared
    once per row of the group; adapter False runs under disable_adapter."""
    import contextlib

    from trl import GRPOTrainer

    seen = []

    def forward(self, model, ids, mask, T, batch_size=None, **kw):
        seen.append((ids.shape, T, batch_size, model.off, kw.get("num_images"), int(kw["mm_token_type_ids"].sum())))
        return torch.full((batch_size, T), -1.0), None

    monkeypatch.setattr(GRPOTrainer, "_get_per_token_logps_and_entropies", forward)

    class Model:
        off = False
        config = SimpleNamespace(image_token_id=IMG)

        @contextlib.contextmanager
        def disable_adapter(self):
            Model.off = True
            yield
            Model.off = False

    t = SimpleNamespace(model=Model(), _tokenizer=SimpleNamespace(pad_token_id=0),
                        accelerator=SimpleNamespace(device="cpu"),
                        processing_class=SimpleNamespace(image_processor=lambda images, return_tensors: {
                            "pixel_values": torch.zeros(len(images), 4), "image_grid_thw": torch.ones(len(images), 3)}))
    img = object()
    p1 = {"prompt_token_ids": [START, IMG, IMG, END], "multi_modal_data": {"image": img}}
    p2 = {"prompt_token_ids": [START, END]}
    lp = lambda *v: SimpleNamespace(token_ids=list(range(1, len(v) + 1)),  # noqa: E731
                                    logprobs=[{i + 1: SimpleNamespace(logprob=x)} for i, x in enumerate(v)])
    outs = [SimpleNamespace(outputs=[o]) for o in (lp(-1.0, -1.5), lp(-0.5), lp(-1.0), lp(-2.0, -1.0, -1.0))]
    gap = fast.FastGRPO.way_b_gap(t, [p1, p1, p2, p2], outs, False)
    assert gap == pytest.approx((0 + 0.5 + 0.5 + 0 + 1.0 + 0 + 0) / 7)  # padded positions never count
    assert seen == [((2, 6), 2, 2, True, [1, 1], 4), ((2, 5), 3, 2, True, None, 0)]
    fast.FastGRPO.way_b_gap(t, [p1], outs[:1], True)
    assert seen[-1][3] is False


def test_comparator_skip_reads_the_ckpt_revision_marker(tmp_path):
    """4b.1: the output-level check runs only on the --ckpt that hf_pull fetched as rl-full-3ep."""
    assert "missing" in fc.comparator_skip(tmp_path)
    (tmp_path / "hf_revision.json").write_text(json.dumps({"repo": "r", "revision": "r1-qwen3vl2b"}))
    assert fc.comparator_skip(tmp_path) == "--ckpt is HF revision 'r1-qwen3vl2b', not rl-full-3ep"
    (tmp_path / "hf_revision.json").write_text(json.dumps({"repo": "r", "revision": "rl-full-3ep"}))
    assert fc.comparator_skip(tmp_path) is None


def test_on_another_model_the_output_arms_are_skipped_and_the_reason_recorded(tmp_path, monkeypatch, lora_request):
    """Only the timing (2 warm-ups + 12 timed) and the restore run; no merged push; the record still passes on the
    other checks, with the reason in f_check.json."""
    fixed, log, vg, m, base = comparator_run(tmp_path, monkeypatch, True, skip="--ckpt is HF revision 'x'",
                                             gap=lambda *a: pytest.fail("no scoring on another model"))
    gens = [e for e in log if e[0] == "generate"]
    assert len(gens) == 14 and all(e[3] for e in gens)  # every generate is a fixed-length timing one
    names = [e[0] for e in log]
    assert names.count("load_weights") == 2  # sync_lora's base push (mode B) and the restore: no merged push
    assert fixed["comparator"] == {"skipped": "--ckpt is HF revision 'x'", "informative": False, "ok": True}
    assert fixed["base_hash_ok"] is True and log[-1] == ("sleep", 2)
    cb = fc.FCheck(tmp_path / "out", "lora-seqmask", {"lora_equal": lambda: True}, loads=fc.count_calls(lambda: None),
                   fixed=fixed)
    cb.loads()
    cb.finish()
    rec = json.loads((tmp_path / "out/f_check.json").read_text())
    assert rec["comparator"]["skipped"] == "--ckpt is HF revision 'x'" and rec["passed"] is True

"""C5: vllm_sync with a fake in-process vLLM engine (no vLLM, no GPU): one load_weights call per push, the same
names and values as fast.lora_merged_params, no trainer parameter written, TRL's wake order, the option-F sync."""

import json
import sys
import types
from types import SimpleNamespace

import pytest

torch = pytest.importorskip("torch")
pytest.importorskip("trl")
from objdel.rl import trl_adapter as fast  # noqa: E402
from objdel.rl import vllm_sync  # noqa: E402
from tests.rl.test_trl_adapter_equivalence import _lora_bf16  # noqa: E402


class FakeEngine:
    """vg.llm with the in-process path TRL uses; every call is appended to `log`."""

    def __init__(self, log):
        self.log, self.loads = log, []
        model = SimpleNamespace(load_weights=self.load_weights)
        executor = SimpleNamespace(driver_worker=SimpleNamespace(model_runner=SimpleNamespace(model=model)))
        self.llm_engine = SimpleNamespace(model_executor=executor, add_lora=lambda req: log.append(("add_lora", req)))

    def load_weights(self, weights):
        pairs = [(n, t.clone()) for n, t in weights]  # consumes the generator, as vLLM does
        self.loads.append(pairs)
        self.log.append(("load_weights", len(pairs)))

    def wake_up(self, tags=None):
        self.log.append(("wake_up", tuple(tags)))

    def reset_prefix_cache(self):
        self.log.append(("reset_prefix_cache",))


def fake_vg(model, monkeypatch, sleeping: bool, sleep_mode: bool | None = None):
    from trl.generation import vllm_generation

    log = []
    monkeypatch.setattr(vllm_generation, "empty_cache", lambda: log.append(("empty_cache",)))
    vg = SimpleNamespace(llm=FakeEngine(log), model=model, _llm_weights_sleeping=sleeping, mode="colocate",
                         enable_sleep_mode=sleeping if sleep_mode is None else sleep_mode,
                         _dist=SimpleNamespace(is_fsdp=False))  # fmt: skip
    return vg, log


@pytest.fixture
def lora_request(monkeypatch):
    """vllm.lora.request.LoRARequest stand-in (vLLM is not installed on the Mac)."""

    class LoRARequest(SimpleNamespace):
        def __init__(self, lora_name, lora_int_id, lora_path, load_inplace=False):
            super().__init__(lora_name=lora_name, lora_int_id=lora_int_id, lora_path=lora_path,
                             load_inplace=load_inplace)  # fmt: skip

    mods = {"vllm": types.ModuleType("vllm"), "vllm.lora": types.ModuleType("vllm.lora"),
            "vllm.lora.request": types.ModuleType("vllm.lora.request")}  # fmt: skip
    mods["vllm.lora.request"].LoRARequest = LoRARequest
    for k, m in mods.items():
        monkeypatch.setitem(sys.modules, k, m)
    return LoRARequest


def test_merged_push_is_one_call_of_every_merged_weight_and_writes_no_parameter(monkeypatch):
    m = _lora_bf16()
    before = {n: p.detach().clone() for n, p in m.named_parameters()}
    want = {n: t.clone() for n, t in fast.lora_merged_params(m)}
    vg, log = fake_vg(m, monkeypatch, sleeping=False)
    stats = vllm_sync.push_weights(vg, m)
    assert log == [("load_weights", len(want)), ("reset_prefix_cache",)]
    got = dict(vg.llm.loads[0])
    assert got.keys() == want.keys() and all(torch.equal(got[n], want[n]) for n in want)
    assert stats.n_tensors == len(want) and stats.bytes == sum(t.numel() * 2 for t in want.values())
    assert all(torch.equal(p, before[n]) for n, p in m.named_parameters())


def test_a_sleeping_engine_wakes_its_weights_before_the_push_in_trls_order(monkeypatch):
    m = _lora_bf16()
    vg, log = fake_vg(m, monkeypatch, sleeping=True)
    vllm_sync.push_weights(vg, m)
    assert [e[0] for e in log] == ["empty_cache", "wake_up", "load_weights", "reset_prefix_cache"]
    assert log[1] == ("wake_up", ("weights",)) and vg._llm_weights_sleeping is False
    vllm_sync.push_weights(vg, m)  # awake now: no second wake
    assert [e[0] for e in log[4:]] == ["load_weights", "reset_prefix_cache"]


def test_base_values_are_the_parameters_without_the_delta(monkeypatch):
    m = _lora_bf16()
    vg, _ = fake_vg(m, monkeypatch, sleeping=False)
    vllm_sync.push_weights(vg, m, values="base")
    got, merged = dict(vg.llm.loads[0]), dict(fast.lora_merged_params(m))
    params = {n.removeprefix("base_model.model.").replace(".base_layer", ""): p for n, p in m.named_parameters()}
    assert got.keys() == merged.keys() and all(torch.equal(got[n], params[n].data) for n in got)
    assert any(not torch.equal(got[n], merged[n]) for n in got)  # the adapter moves some weights: not vacuous
    with pytest.raises(ValueError, match="values"):
        vllm_sync.push_weights(vg, m, values="delta")


def test_install_replaces_the_sync_on_one_colocate_process_only(monkeypatch):
    m = _lora_bf16()
    vg, log = fake_vg(m, monkeypatch, sleeping=False)
    t = SimpleNamespace(vllm_generation=vg, accelerator=SimpleNamespace(num_processes=1))
    assert vllm_sync.install(t)
    vg.sync_weights()
    assert [e[0] for e in log] == ["load_weights", "reset_prefix_cache"]
    for bad in (SimpleNamespace(vllm_generation=None, accelerator=t.accelerator),
                SimpleNamespace(vllm_generation=vg, accelerator=SimpleNamespace(num_processes=2))):
        assert not vllm_sync.install(bad)
    vg.mode = "server"
    assert not vllm_sync.install(t)


@pytest.mark.parametrize("sleep_mode", [True, False])
def test_option_f_sync_order_and_adapter_folder(monkeypatch, tmp_path, lora_request, sleep_mode):
    """U-SY8: mode B wakes, pushes the base W once (no delta), writes the bf16 PEFT folder, adds the adapter once
    with load_inplace=True and resets the prefix cache once; mode A pushes no base."""
    from safetensors.torch import load_file

    m = _lora_bf16()
    vg, log = fake_vg(m, monkeypatch, sleeping=sleep_mode)
    vllm_sync.sync_lora(vg, m, tmp_path / "lora")
    names = [e[0] for e in log]
    if sleep_mode:
        assert names == ["empty_cache", "wake_up", "load_weights", "add_lora", "reset_prefix_cache"]
        params = {n.removeprefix("base_model.model.").replace(".base_layer", ""): p for n, p in m.named_parameters()}
        assert all(torch.equal(t, params[n].data) for n, t in vg.llm.loads[0])
    else:
        assert names == ["add_lora", "reset_prefix_cache"]
    req = log[names.index("add_lora")][1]
    assert (req.lora_name, req.lora_int_id, req.load_inplace) == ("policy", 1, True)
    sd = load_file(str(tmp_path / "lora/adapter_model.safetensors"))
    assert sd and all(t.dtype == torch.bfloat16 for t in sd.values()) and not any(".default." in k for k in sd)
    assert json.loads((tmp_path / "lora/adapter_config.json").read_text())["r"] == m.peft_config["default"].r


def test_install_lora_seqmask_routes_to_the_f_sync(monkeypatch, tmp_path, lora_request):
    m = _lora_bf16()
    vg, log = fake_vg(m, monkeypatch, sleeping=False)
    vg.llm.generate = lambda *a, **kw: None
    t = SimpleNamespace(vllm_generation=vg, accelerator=SimpleNamespace(num_processes=1))
    with pytest.raises(ValueError, match="adapter_dir"):
        vllm_sync.install(t, "lora-seqmask")
    assert vllm_sync.install(t, "lora-seqmask", tmp_path)
    vg.sync_weights()
    assert [e[0] for e in log] == ["add_lora", "reset_prefix_cache"]


def test_install_leaves_generate_alone(monkeypatch, tmp_path, lora_request):
    """install changes the sync only; under option F the LoRARequest rides on FastGRPO.lora_dir (patch_vllm_prompts)."""
    m = _lora_bf16()
    vg, _ = fake_vg(m, monkeypatch, sleeping=False)
    gen = vg.llm.generate = lambda prompts, **kw: kw
    t = SimpleNamespace(vllm_generation=vg, accelerator=SimpleNamespace(num_processes=1))
    vllm_sync.install(t)
    assert vg.llm.generate is gen
    vllm_sync.install(t, "lora-seqmask", tmp_path)
    assert vg.llm.generate is gen

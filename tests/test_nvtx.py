"""A6: objdel.train.nvtx (inert import, cached CUDA query, balanced ranges), the FastGRPO ranges and counters,
callbacks.NvtxSteps (step ranges and the capture window), the make_trainer / train_rl wiring, the `checker` range on
the reward path, and S12 (`reward` never imports `train`). CUDA and NVTX are fakes (monkeypatch); no GPU."""

import importlib
import os
import subprocess
import sys
import time

import numpy as np
import pytest
import torch

from objdel.train import nvtx


@pytest.fixture
def fresh(monkeypatch):
    """reload(env, cuda) -> the nvtx module re-imported with OBJDEL_NVTX = env (None: unset) and a counted
    torch.cuda.is_available returning `cuda`; the module is re-imported with the test's starting env afterwards."""
    def reload(env, cuda=True):
        calls = []
        if env is None:
            monkeypatch.delenv("OBJDEL_NVTX", raising=False)
        else:
            monkeypatch.setenv("OBJDEL_NVTX", env)
        monkeypatch.setattr(torch.cuda, "is_available", lambda: calls.append(1) or cuda)
        return importlib.reload(nvtx), calls

    yield reload
    monkeypatch.undo()
    importlib.reload(nvtx)


class Rec:
    def __init__(self):
        self.ev = []

    def range_push(self, n):
        self.ev.append(("push", n))

    def range_pop(self):
        self.ev.append(("pop",))

    def start(self):
        self.ev.append(("start",))

    def stop(self):
        self.ev.append(("stop",))


@pytest.fixture
def rec(monkeypatch):
    """NVTX enabled, CUDA faked, every torch.cuda.nvtx / profiler call recorded in one event list."""
    r = Rec()
    monkeypatch.setattr(nvtx, "enabled", lambda: True)
    monkeypatch.setattr(torch.cuda, "is_available", lambda: True)
    for mod, names in ((torch.cuda.nvtx, ("range_push", "range_pop")), (torch.cuda.profiler, ("start", "stop"))):
        for n in names:
            monkeypatch.setattr(mod, n, getattr(r, n), raising=False)
    return r


@pytest.fixture
def no_nvtx_calls(monkeypatch):
    """Every torch.cuda.nvtx and profiler function raises: a disabled path must call none."""
    def boom(*a, **k):
        raise AssertionError("NVTX / profiler call while disabled")

    for mod, names in ((torch.cuda.nvtx, ("range_push", "range_pop")), (torch.cuda.profiler, ("start", "stop"))):
        for n in names:
            monkeypatch.setattr(mod, n, boom, raising=False)


# U25 ------------------------------------------------------------------------------------------------------------


def test_import_is_inert_it_imports_no_torch():
    """objdel.reward is imported before serve.fork_pool: an import-time CUDA query would poison the fork."""
    code = "import sys, objdel.train.nvtx as n; assert 'torch' not in sys.modules, 'torch imported'; print(n._ON)"
    out = subprocess.run([sys.executable, "-I", "-c", code], capture_output=True, text=True,
                         env=os.environ | {"OBJDEL_NVTX": "1"}, cwd=os.path.dirname(__file__))  # fmt: skip
    assert out.returncode == 0, out.stderr
    assert out.stdout.strip() == "True"


@pytest.mark.parametrize("env", [None, "0", "true", ""])
def test_enabled_is_false_without_the_env_and_never_asks_cuda(fresh, env):
    n, calls = fresh(env)
    assert n.enabled() is False and calls == []


def test_enabled_needs_cuda_and_asks_once(fresh, monkeypatch):
    n, calls = fresh("1", cuda=False)
    assert n.enabled() is False
    n, calls = fresh("1", cuda=True)
    assert all(n.enabled() for _ in range(1000)) and len(calls) == 1
    monkeypatch.setenv("OBJDEL_NVTX", "0")  # read once, at import
    assert n.enabled() is True and len(calls) == 1


# U26 ------------------------------------------------------------------------------------------------------------


def test_range_pushes_and_pops_also_when_the_body_raises(rec):
    with nvtx.range("a"):
        with nvtx.range("b"):
            pass
    with pytest.raises(ValueError), nvtx.range("c"):
        raise ValueError
    nvtx.push("d")
    nvtx.pop()
    assert rec.ev == [("push", "a"), ("push", "b"), ("pop",), ("pop",), ("push", "c"), ("pop",), ("push", "d"),
                      ("pop",)]  # fmt: skip


@pytest.mark.parametrize("env,cuda", [(None, True), ("1", False)])
def test_disabled_range_push_pop_call_no_nvtx(fresh, no_nvtx_calls, env, cuda):
    n, _ = fresh(env, cuda)
    with n.range("x"):
        n.push("y")
        n.pop()


# Q03 (cost when off) ---------------------------------------------------------------------------------------------


def test_nvtx_is_free_when_off(fresh):
    slack = max(1.0, os.getloadavg()[0] / (os.cpu_count() or 1))  # load-tolerant: 0.52 s seen at load 40 on 8 cores
    n, _ = fresh(None)
    t = time.perf_counter()
    for _ in range(100_000):
        with n.range("x"):
            pass
    assert time.perf_counter() - t < 0.5 * slack
    n, _ = fresh("1", cuda=True)
    n.enabled()
    t = time.perf_counter()
    for _ in range(100_000):
        n.enabled()
    assert (time.perf_counter() - t) / 100_000 < 1e-6 * slack


# U27, U28: NvtxSteps ---------------------------------------------------------------------------------------------


def drive(cb, n, first=0, train_begin=True):
    """The Trainer's callback order for n optimizer steps after `first` (a resume starts later)."""
    import types

    st = types.SimpleNamespace(global_step=first)
    if train_begin:
        cb.on_train_begin(None, st, None)
    for i in range(first + 1, first + n + 1):
        cb.on_pre_optimizer_step(None, st, None)
        cb.on_optimizer_step(None, st, None)
        st.global_step = i
        cb.on_step_end(None, st, None)
    cb.on_train_end(None, st, None)
    return st


@pytest.mark.parametrize("n", [0, 1, 2, 10, 30])
def test_pushes_and_pops_balance_and_nest(rec, n):
    from objdel.train.callbacks import NvtxSteps

    drive(NvtxSteps(None), n)
    stack = []
    for e in rec.ev:
        if e[0] == "push":
            stack.append(e[1])
        elif e[0] == "pop":
            stack.pop()
        assert len(stack) <= 2 and (not stack or stack[0] == "step")
    assert stack == [] and not any(e[0] in ("start", "stop") for e in rec.ev)
    assert rec.ev.count(("push", "step")) == n + 1 and rec.ev.count(("push", "optimizer")) == n
    if n:
        assert rec.ev[:4] == [("push", "step"), ("push", "optimizer"), ("pop",), ("pop",)]


def test_unpaired_hooks_pop_nothing(rec):
    import types

    from objdel.train.callbacks import NvtxSteps

    cb, st = NvtxSteps(None), types.SimpleNamespace(global_step=3)
    cb.on_optimizer_step(None, st, None)  # no pre-step (a skipped optimizer step)
    cb.on_train_end(None, st, None)  # no train begin (an earlier exception)
    assert rec.ev == []


def test_profiler_window_starts_after_pop_before_push_and_stops_once(rec):
    from objdel.train.callbacks import NvtxSteps

    drive(NvtxSteps((4, 6)), 8)
    names = [e[0] for e in rec.ev]
    assert names.count("start") == 1 and names.count("stop") == 1
    i, j = names.index("start"), names.index("stop")
    assert names[i - 1] == "pop" and names[i + 1] == "push" and j > i and names[j - 1] == "pop"
    assert names[j + 1] == "push"  # training goes on after the window
    # the capture starts at the end of global_step 3: pushes of "step" before it = train_begin, after step 1, after 2
    assert len([e for e in rec.ev[:i] if e == ("push", "step")]) == 3
    # it holds exactly the step ranges of steps 4..6
    assert len([e for e in rec.ev[i:j] if e == ("push", "step")]) == 3


def test_profiler_window_edges(rec, monkeypatch):
    import types

    from objdel.train.callbacks import NvtxSteps

    cb, st = NvtxSteps((4, 6)), types.SimpleNamespace(global_step=0)
    cb.on_train_begin(None, st, None)
    for _ in range(2):  # the same global_step twice: one start
        st.global_step = 3
        cb.on_step_end(None, st, None)
    assert [e for e in rec.ev if e[0] == "start"] == [("start",)]
    rec.ev.clear()
    drive(NvtxSteps((4, 6)), 2, first=5)  # resumed inside the window (BR1-rl): train begin starts, step 6 stops
    assert [e[0] for e in rec.ev if e[0] in ("start", "stop")] == ["start", "stop"]
    rec.ev.clear()
    drive(NvtxSteps((4, 6)), 2, first=6)  # resumed past the window: no start, no stop
    assert not any(e[0] in ("start", "stop") for e in rec.ev)
    rec.ev.clear()
    drive(NvtxSteps((4, 6)), 4)  # the run ends at step 4, inside the window: on_train_end stops the capture
    names = [e[0] for e in rec.ev]
    assert names.count("start") == 1 and names[-2:] == ["pop", "stop"]
    for bad in ((1, 3), (0, 2), (5, 4)):
        with pytest.raises(ValueError):
            NvtxSteps(bad)


def test_no_cuda_no_profiler_and_no_ranges(monkeypatch, no_nvtx_calls):
    from objdel.train.callbacks import NvtxSteps

    monkeypatch.setattr(nvtx, "enabled", lambda: False)
    monkeypatch.setattr(torch.cuda, "is_available", lambda: False)
    drive(NvtxSteps((4, 6)), 8)
    drive(NvtxSteps(None), 3)


# U30, U56: FastGRPO.timed ----------------------------------------------------------------------------------------

NAMES = {"step", "loader", "micro", "batch", "generate", "sync", "wake", "sleep", "rewards", "checker", "logps",
         "forward", "optimizer", "log", "vision_cache", "rope_cache", "image_memo"}


def fake_trainer(**attrs):
    from collections import defaultdict
    from types import SimpleNamespace

    from objdel.train import fast

    t = object.__new__(fast.FastGRPO)
    t.model, t._metrics = SimpleNamespace(training=True), {"train": defaultdict(list), "eval": defaultdict(list)}
    t.off_policy_mask_threshold = None  # TRL's attribute
    t.accelerator = SimpleNamespace(is_main_process=False)  # log writes no seq_gaps file
    for k, v in attrs.items():
        setattr(t, k, v)
    return t


class FakeCuda:
    """torch.cuda with a sync counter and a CUDA event recorder (each event's time is set at record())."""

    def __init__(self, monkeypatch, clock):
        self.syncs, self.events = 0, []
        cuda = self

        class Event:
            def __init__(self, enable_timing=False):
                assert enable_timing
                self.t = None
                cuda.events.append(self)

            def record(self, stream=None):
                self.t = clock()

            def synchronize(self):
                pass

            def elapsed_time(self, end):
                return (end.t - self.t) * 1e3  # ms

        monkeypatch.setattr(torch.cuda, "is_available", lambda: True)
        monkeypatch.setattr(torch.cuda, "synchronize", self.sync)
        monkeypatch.setattr(torch.cuda, "Event", Event)

    def sync(self, device=None):
        self.syncs += 1


class Clock:
    def __init__(self, *ticks):
        self.ticks = list(ticks)

    def __call__(self):
        return self.ticks.pop(0)


def capture_log(monkeypatch) -> list:
    """GRPOTrainer.log replaced by a recorder of the train metrics it would log."""
    from objdel.train import fast

    logged = []

    def log(self, logs, start_time=None):
        logged.append(dict(self._metrics["train"]))

    monkeypatch.setattr(fast.GRPOTrainer, "log", log)
    return logged


def test_timed_opens_the_range_and_records_the_phase_also_when_it_raises(rec, monkeypatch):
    monkeypatch.setattr(torch.cuda, "synchronize", lambda *a: None)
    t = fake_trainer()
    with t.timed("generate"):
        pass
    with pytest.raises(ValueError), t.timed("generate"):
        raise ValueError
    assert rec.ev == [("push", "generate"), ("pop",)] * 2
    got = t._metrics["train"]["time/generate_s"]
    assert len(got) == 2 and all(isinstance(x, float) and x >= 0 for x in got)


def test_timed_syncs_once_by_default_and_uses_events_when_sync_is_off(monkeypatch):
    from objdel.train import fast

    clock = Clock(0.0, 2.0, 10.0, 13.5)
    cuda = FakeCuda(monkeypatch, clock)
    monkeypatch.setattr(fast.time, "perf_counter", clock)
    monkeypatch.setattr(nvtx, "SYNC", True)
    t = fake_trainer()
    with t.timed("rewards"):
        pass
    assert cuda.syncs == 1 and cuda.events == [] and t._metrics["train"]["time/rewards_s"] == [2.0]
    monkeypatch.setattr(nvtx, "SYNC", False)
    t = fake_trainer()
    with t.timed("generate"):
        pass
    assert cuda.syncs == 1 and len(cuda.events) == 2  # one event pair, no sync
    assert "time/generate_s" not in t._metrics["train"]  # read later, in log's host read
    logged = capture_log(monkeypatch)
    monkeypatch.setattr(fast.time, "perf_counter", lambda: 0.0)
    t.log({})
    assert cuda.syncs == 1 and logged[0]["time/generate_s"] == [pytest.approx(3.5)]


def test_timed_without_cuda_uses_the_host_clock_even_with_sync_off(monkeypatch):
    from objdel.train import fast

    monkeypatch.setattr(torch.cuda, "is_available", lambda: False)
    monkeypatch.setattr(torch.cuda, "Event", lambda **k: pytest.fail("event without CUDA"))
    monkeypatch.setattr(fast.time, "perf_counter", Clock(1.0, 4.0))
    monkeypatch.setattr(nvtx, "SYNC", False)
    t = fake_trainer()
    with t.timed("sync"):
        pass
    assert t._metrics["train"]["time/sync_s"] == [3.0]


def range_names():
    """Every string literal passed to nvtx.range / nvtx.push / .timed under src/objdel/train (AST scan), and the
    phase names of FastGRPO.time_vllm_phases's `hooks` table (its wrappers pass them to `timed` as `_phase`)."""
    import ast
    from pathlib import Path

    import objdel.train

    names, bad = set(), []
    for f in Path(objdel.train.__file__).parent.glob("*.py"):
        for node in ast.walk(ast.parse(f.read_text())):
            call = isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute)
            if call and (node.func.attr == "timed" or node.func.attr in ("range", "push")
                         and getattr(node.func.value, "id", "") == "nvtx"):  # fmt: skip
                arg = node.args[0] if node.args else None
                if isinstance(arg, ast.Constant) and isinstance(arg.value, str):
                    names.add(arg.value)
                elif f.name != "nvtx.py" and getattr(arg, "id", "") not in ("phase", "_phase"):  # timed, its wrappers
                    bad.append(f"{f.name}:{node.lineno}")
            if isinstance(node, ast.Assign) and getattr(node.targets[0], "id", "") == "hooks":
                names |= {e.elts[2].value for e in node.value.elts}
    return names, bad


def test_every_fixed_range_name_is_used_and_no_other():
    names, bad = range_names()
    assert names == NAMES, (names - NAMES, NAMES - names)
    assert bad == [], "range names must be literals"


# U31, U32: the IS ratio split and its share counters --------------------------------------------------------------

MODES = ("sequence_mask", "sequence_truncate", "token_mask", "token_truncate")


@pytest.mark.parametrize("mode", MODES)
@pytest.mark.parametrize("bounds", [(None, 3.0), (0.1, 3.0)])  # TRL: *_truncate needs a bound; clip_max default 3
def test_raw_ratio_then_bounded_is_is_ratio(mode, bounds):
    from objdel.train import fast

    g = torch.Generator().manual_seed(0)
    for _ in range(50):
        B, T = 6, 40
        logps = -torch.rand(B, T, generator=g) * 5
        sampling = logps - (torch.rand(B, T, generator=g) * 40 - 20)  # gaps in [-20, 20]: exp overflows, underflows
        sampling[0, 3] = -torch.inf
        mask = torch.arange(T)[None] < torch.randint(0, T + 1, (B, 1), generator=g)
        mask[1], mask[2] = False, True  # an all-masked and a fully valid row
        lo, hi = bounds
        want = fast.is_ratio(logps, sampling, mask, mode, lo, hi)
        got = fast.bounded(fast.raw_ratio(logps, sampling, mask, mode), mode, lo, hi)
        assert torch.equal(got, want)


RAW = torch.tensor([0.0, 0.05, 0.5, 1.0, 2.9, 3.0, 3.1, 100.0])


@pytest.mark.parametrize("mode,lo,want", [("sequence_mask", 0.1, [4, 2, 2]), ("sequence_truncate", None, [1, 2, 1]),
                                          ("sequence_truncate", 0.1, [0, 2, 2])])  # fmt: skip
def test_share_counts_equal_numpy_counts(mode, lo, want, monkeypatch):
    from objdel.train import fast

    raw = RAW[:, None]
    r = fast.bounded(raw, mode, lo, 3.0)
    for name in ("item", "tolist", "__float__"):  # no host read in the computation
        monkeypatch.setattr(torch.Tensor, name, lambda *a: pytest.fail("host read"))
    got = fast.clip_counts(raw, r, torch.ones_like(raw, dtype=torch.bool), lo, 3.0)
    monkeypatch.undo()
    assert isinstance(got, torch.Tensor) and got.tolist() == want
    n = RAW.numpy()
    lo_ = -np.inf if lo is None else lo
    assert want == [int((r.numpy() == 0).sum()), int((n > 3.0).sum()), int(((n < lo_) | (n == 0)).sum())]


def test_token_counts_leave_padding_out():
    from objdel.train import fast

    raw = torch.tensor([[0.5, 4.0, 0.0], [1.0, 0.0, 0.0]])  # padding (mask 0) underflows to 0
    keep = torch.tensor([[True, True, True], [True, False, False]])
    r = fast.bounded(raw, "token_mask", None, 3.0)
    assert fast.clip_counts(raw, r, keep, None, 3.0).tolist() == [2, 1, 1]
    assert fast.clip_counts(raw, r, torch.zeros_like(keep), None, None).tolist() == [0, 0, 0]


def test_log_writes_the_shares_from_one_read(monkeypatch):
    """The update micro-batches' counts stay on the device; log turns them into shares of the counted entries."""
    logged = capture_log(monkeypatch)
    t = fake_trainer(vllm_importance_sampling_mode="sequence_mask", vllm_importance_sampling_clip_min=0.1,
                     vllm_importance_sampling_clip_max=3.0)  # fmt: skip
    gaps = torch.log(RAW[:, None]).nan_to_num(neginf=-200.0)  # sequence gaps whose exp is RAW (0: underflow)
    for rows in (slice(0, 5), slice(5, 8)):  # two micro-batches
        inputs = {"sampling_per_token_logps": -gaps[rows], "completion_mask": torch.ones_like(gaps[rows]),
                  "advantages": torch.ones(len(gaps[rows]))}
        t.set_is_ratio(torch.zeros_like(gaps[rows]), inputs)
    t.log({})
    m = logged[0]
    shares = {k: m[f"sampling/importance_sampling_ratio/{k}_share"] for k in ("zero", "hi", "lo")}
    assert shares == {"zero": [4 / 8], "hi": [2 / 8], "lo": [2 / 8]}
    assert m["sampling/importance_sampling_ratio/min"] == [0.0]


# U33, Q03: time/step_s and no new host read in log ----------------------------------------------------------------


def test_step_time_is_the_gap_between_logs_without_a_sync(monkeypatch):
    from objdel.train import fast

    logged = capture_log(monkeypatch)
    monkeypatch.setattr(fast.time, "perf_counter", Clock(10.0, 40.0, 75.0))
    monkeypatch.setattr(torch.cuda, "synchronize", lambda *a: pytest.fail("sync in log"))
    t = fake_trainer()
    for _ in range(3):
        t.log({})
        t._metrics["train"].clear()  # as TRL's log does
    assert [x.get("time/step_s") for x in logged] == [None, [30.0], [35.0]]


def test_log_adds_no_host_read(monkeypatch):
    """Q03: with the counters, the step time and the NVTX range, log reads the device once (the stats), as before."""
    from objdel.train import fast

    monkeypatch.setattr(fast.GRPOTrainer, "log", lambda self, logs, start=None: None)
    reads = []
    for name in ("item", "cpu", "tolist"):
        orig = getattr(torch.Tensor, name)
        monkeypatch.setattr(torch.Tensor, name, lambda self, *a, _o=orig, _n=name: reads.append(_n) or _o(self, *a))
    monkeypatch.setattr(torch.cuda, "synchronize", lambda *a: reads.append("sync"))
    t = fake_trainer(vllm_importance_sampling_mode="token_mask", vllm_importance_sampling_clip_min=None,
                     vllm_importance_sampling_clip_max=3.0)  # fmt: skip
    inputs = {"sampling_per_token_logps": torch.zeros(2, 3), "completion_mask": torch.ones(2, 3),
              "advantages": torch.ones(2)}
    t.set_is_ratio(torch.zeros(2, 3), inputs)
    assert reads == []
    t.log({})
    assert reads == ["tolist"]


# I11: the Trainer-loop ranges on a real FastGRPO -----------------------------------------------------------------


@pytest.fixture(scope="module")
def tiny_grpo(tmp_path_factory):
    """A FastGRPO on a 1-layer Llama (CPU, HF generation), groups of 8 answers to one prompt."""
    from datasets import Dataset
    from tokenizers import Tokenizer, models, pre_tokenizers
    from transformers import LlamaConfig, LlamaForCausalLM, PreTrainedTokenizerFast
    from trl import GRPOConfig

    from objdel.train import fast

    vocab = {w: i for i, w in enumerate(["<pad>", "<eos>", "prompt", *map(str, range(40))])}
    tk = Tokenizer(models.WordLevel(vocab, unk_token="<pad>"))
    tk.pre_tokenizer = pre_tokenizers.WhitespaceSplit()
    tok = PreTrainedTokenizerFast(tokenizer_object=tk, pad_token="<pad>", eos_token="<eos>", padding_side="left")
    torch.manual_seed(0)
    model = LlamaForCausalLM(LlamaConfig(vocab_size=len(vocab), hidden_size=16, intermediate_size=32,
                                         num_hidden_layers=1, num_attention_heads=2, num_key_value_heads=2,
                                         pad_token_id=0, eos_token_id=1, bos_token_id=1))  # fmt: skip
    args = GRPOConfig(output_dir=str(tmp_path_factory.mktemp("grpo")), max_steps=1, per_device_train_batch_size=8,
                      num_generations=8, max_completion_length=4, use_cpu=True, bf16=False, report_to="none",
                      logging_steps=1, save_strategy="no", seed=0, remove_unused_columns=False)  # fmt: skip
    t = fast.FastGRPO(model=model, reward_funcs=lambda completions, **kw: [sum(map(ord, c)) % 7.0 for c in completions],
                      args=args, train_dataset=Dataset.from_list([{"prompt": "prompt 1"}] * 2), processing_class=tok)
    t.current_gradient_accumulation_steps = 1  # set by the Trainer loop
    return t


def test_trainer_loop_ranges_do_not_change_results(tiny_grpo, monkeypatch):
    from trl import GRPOTrainer

    t = tiny_grpo

    def step():
        torch.manual_seed(1)
        t.model.zero_grad()
        batches, n = t.get_batch_samples(iter(GRPOTrainer.get_train_dataloader(t)), 1, torch.device("cpu"))
        loss = t.training_step(t.model, batches[0], n)
        return batches, loss, [p.grad.clone() for p in t.model.parameters() if p.grad is not None]

    monkeypatch.setattr(nvtx, "enabled", lambda: False)
    b0, l0, g0 = step()
    r = Rec()
    monkeypatch.setattr(nvtx, "enabled", lambda: True)
    monkeypatch.setattr(torch.cuda.nvtx, "range_push", r.range_push)
    monkeypatch.setattr(torch.cuda.nvtx, "range_pop", r.range_pop)
    b1, l1, g1 = step()
    assert torch.equal(l0, l1) and len(g0) == len(g1) > 0
    assert all(torch.equal(a, b) for a, b in zip(g0, g1, strict=True))
    assert b0 == b1
    pushed = [e[1] for e in r.ev if e[0] == "push"]
    assert pushed.count("loader") == 1 and pushed.count("micro") == 1
    assert {"batch", "generate", "rewards", "forward", "logps"} <= set(pushed)
    assert len(pushed) == r.ev.count(("pop",))
    r.ev.clear()
    t.log({"loss": 0.0})
    assert r.ev == [("push", "log"), ("pop",)]


# U57: the checker range on the reward path ----------------------------------------------------------------------


def test_checker_range_opens_on_the_thread_that_scores(monkeypatch):
    import threading

    import numpy as np

    from objdel.train import rewards as R

    seen = []
    monkeypatch.setattr(nvtx, "range", lambda name: seen.append((name, threading.current_thread().name))
                        or __import__("contextlib").nullcontext())
    T = {"format_ok": np.ones(8), "legal": np.ones(8), "iou_tc": np.linspace(0, 1, 8), "f5": np.ones(8),
         "g_area": np.ones(8), "topo_any": np.ones(8), "topo_freq": np.ones(8), "sig_id": np.zeros(8),
         "loc": np.zeros(8), "error": [""] * 8, "defects": [{}] * 8}  # fmt: skip
    monkeypatch.setattr(R, "score_batch", lambda jobs: T)
    rw = R.TRLReward(lambda k: {}, stop_ids={1})
    args = ([[{"content": "x"}]] * 8, [[1]] * 8, ["k"] * 8, [5] * 8, [5] * 8, [["nodule"]] * 8)
    direct = rw._terms(*args, None, None)
    assert seen == [("checker", threading.current_thread().name)]
    rw._last = (None, None)
    rw.prefetch(*args)
    ahead = rw._terms(*args, None, None)
    assert seen[1][0] == "checker" and seen[1][1].startswith("reward-prefetch") and len(seen) == 2
    assert ahead == direct


# U29: make_trainer and train_rl ----------------------------------------------------------------------------------


@pytest.fixture
def built(monkeypatch, tmp_path):
    """make(**a) -> the callbacks g_probe.make_trainer hands the trainer class (model, processor, config, rewards
    and the trainer itself are stubs)."""
    from types import SimpleNamespace

    from experiments.phase2 import g_probe
    from objdel.train import configs, fast, rewards, serve

    got = []

    class Trainer:
        def __init__(self, **kw):
            got.append(kw["callbacks"])
            self.model = None

    stub = SimpleNamespace(from_pretrained=lambda *a, **k: SimpleNamespace())
    monkeypatch.setattr("transformers.AutoModelForImageTextToText", stub)  # by path: imports above swap the module
    monkeypatch.setattr("transformers.AutoProcessor", stub)
    monkeypatch.setattr(configs, "grpo_config", lambda *a, **k: SimpleNamespace())
    monkeypatch.setattr(rewards, "TRLReward", lambda *a, **k: SimpleNamespace(funcs=lambda: ([], [])))
    monkeypatch.setattr(rewards, "stop_ids", lambda *a: set())
    monkeypatch.setattr(serve, "loader", lambda *roots: None)
    monkeypatch.setattr(fast, "FastGRPO", Trainer)

    def make(**over):
        a = SimpleNamespace(model="qwen3-vl-2b", out=tmp_path, train=tmp_path, grpo="{}")
        for k, v in over.items():
            setattr(a, k, v)
        g_probe.make_trainer(a, tmp_path, None, 3)
        return got.pop()

    return make


def test_make_trainer_adds_nvtx_steps_only_with_a_window(built):
    from objdel.train.callbacks import NvtxSteps, TraceStep

    assert built() == [] and built(profile_window=None, trace=None) == []
    (cb,) = built(profile_window=[4, 6])
    assert isinstance(cb, NvtxSteps) and cb.window == (4, 6)
    cbs = built(profile_window=[4, 6], trace="T")
    assert [type(c) for c in cbs] == [NvtxSteps, TraceStep] and cbs[1].wait == 2
    (cb,) = built(trace="T")
    assert isinstance(cb, TraceStep) and cb.wait == 1
    assert built(profile_window=[2, 5], trace="T")[1].wait == 0
    with pytest.raises(ValueError):
        built(profile_window=[1, 3], trace="T")


def test_train_rl_flags():
    from experiments.phase2 import train_rl

    base = ["--model", "m", "--ckpt", "c", "--train", "t", "--out", "o", "--steps", "20"]
    a = train_rl.parse(base)
    assert a.profile_window is None and a.trace is None
    a = train_rl.parse([*base, "--profile-window", "4", "6", "--trace", "D"])
    assert a.profile_window == [4, 6] and a.trace == "D"
    for bad in (["1", "2"], ["6", "4"]):
        with pytest.raises(SystemExit) as e:
            train_rl.parse([*base, "--profile-window", *bad])
        assert e.value.code == 2


# S12: objdel.reward never imports objdel.train -------------------------------------------------------------------


def train_imports(folder):
    """Imports of objdel.train anywhere in the modules under `folder` (function bodies too), as file:line."""
    import ast
    from pathlib import Path

    hits = []
    for f in sorted(Path(folder).rglob("*.py")):
        pkg = ["objdel", *f.relative_to(Path(folder).parent).parent.parts]  # the module's package, e.g. objdel.reward
        for node in ast.walk(ast.parse(f.read_text())):
            if isinstance(node, ast.Import):
                mods = [n.name for n in node.names]
            elif isinstance(node, ast.ImportFrom):
                base = pkg[:len(pkg) - node.level + 1] if node.level else []
                mod = ".".join([*base, *([node.module] if node.module else [])])
                mods = [mod] + [f"{mod}.{n.name}" for n in node.names]
            else:
                continue
            if any(m == "objdel.train" or m.startswith("objdel.train.") for m in mods):
                hits.append(f"{f.name}:{node.lineno}")
    return hits


def test_reward_never_imports_train(tmp_path):
    import shutil
    from pathlib import Path

    import objdel.reward

    src = Path(objdel.reward.__file__).parent
    assert train_imports(src) == []
    bad = tmp_path / "reward"
    shutil.copytree(src, bad, ignore=shutil.ignore_patterns("__pycache__"))
    score = bad / "score.py"
    score.write_text(score.read_text() + "\n\ndef _planted():\n    from objdel.train import nvtx\n    return nvtx\n")
    assert [h.split(":")[0] for h in train_imports(bad)] == ["score.py"]
    score.write_text(score.read_text().replace("from objdel.train import nvtx", "from ..train import nvtx"))
    assert [h.split(":")[0] for h in train_imports(bad)] == ["score.py"]

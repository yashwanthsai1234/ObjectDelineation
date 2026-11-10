"""profile: the RL arms (g_probe timing runs), their summary, the traced step's GPU busy share, nvidia-smi means;
the shared torch.profiler callback (train_sft `profile:`, g_probe --trace; a CPU trace; refused outside timing runs)."""

import gzip
import json

import numpy as np
import pytest

from experiments.phase2 import profile as PS
from objdel.sft import train as sft_train


def test_rl_arms_are_short_g_probe_runs_of_256_answers(tmp_path):
    cmd = PS.arm_cmd("mb32", ckpt="CK", prepared="PREP", probe="VIEW", out=tmp_path)
    assert cmd[:3] == ["-I", "-m", "experiments.phase2.g_probe"] and "--ckpt" in cmd and "CK" in cmd
    grpo = json.loads(cmd[cmd.index("--grpo") + 1])
    assert grpo["max_steps"] == PS.STEPS and grpo["logging_steps"] == 1
    assert grpo["per_device_train_batch_size"] * grpo["gradient_accumulation_steps"] == 256
    assert cmd[cmd.index("--items") + 1] == str(PS.ITEMS) and cmd[cmd.index("--out") + 1] == str(tmp_path / "mb32")
    assert "--plain" in PS.arm_cmd("trl", "CK", "PREP", "VIEW", tmp_path)
    assert PS.arm_cmd("graphs", "CK", "PREP", "VIEW", tmp_path)[-2:] == ["--compile", "reduce-overhead"]
    assert "--trace" not in PS.arm_cmd("fast", "CK", "PREP", "VIEW", tmp_path)  # the profiler would slow it
    assert PS.arm_cmd("trace", "CK", "PREP", "VIEW", tmp_path)[-2:] == ["--trace", str(tmp_path / "trace" / "trace")]
    assert list(PS.ARMS)[-1] == "trace" and "trace" not in PS.TIMED
    for arm in PS.ARMS.values():
        g = {"per_device_train_batch_size": 16, "gradient_accumulation_steps": 16} | arm.get("grpo", {})
        assert g["per_device_train_batch_size"] * g["gradient_accumulation_steps"] == 256


def test_summary_reads_the_phase_times_from_the_g_probe_log():
    log = [{"step_time": 60.0, "time/batch_s": 40.0, "time/generate_s": 20.0, "time/rewards_s": 5.0,
            "completions/mean_length": 400.0, "sampling/importance_sampling_ratio/min": 0.5},
           {"step_time": 30.0, "time/batch_s": 18.0, "time/generate_s": 10.0, "time/rewards_s": 2.0,
            "completions/mean_length": 300.0, "sampling/importance_sampling_ratio/min": 0.0},
           {"step_time": 34.0, "time/batch_s": 20.0, "time/generate_s": 12.0, "time/rewards_s": 4.0,
            "completions/mean_length": 500.0, "sampling/importance_sampling_ratio/min": 0.7},
           {"train_runtime": 1.0}]  # fmt: skip
    s = PS.summary(log)
    assert s["s_per_step"] == 32.0  # step 1 (vLLM and allocator warm-up) left out
    assert s["split_s"] == {"generate": 11.0, "rewards": 3.0, "logprobs": 5.0, "update": 13.0}
    assert s["mean_length"] == 400.0 and s["is_ratio_min"] == 0.0
    plain = PS.summary([{k: v for k, v in h.items() if not k.startswith("time/")} for h in log])  # TRL as shipped
    assert plain["s_per_step"] == 32.0 and "split_s" not in plain


def test_an_arm_that_crashes_is_recorded_and_the_next_arm_runs(tmp_path, monkeypatch):
    calls = []
    smi = type("P", (), {"terminate": lambda s: None, "wait": lambda s: None})()
    monkeypatch.setattr(PS.subprocess, "Popen", lambda *a, **k: smi)
    failed = type("R", (), {"returncode": 1})()
    monkeypatch.setattr(PS.subprocess, "run", lambda cmd, **k: calls.append(cmd) or failed)
    argv = ["--ckpt", "C", "--prepared", "P", "--probe", "V", "--out", str(tmp_path), "--arms", "trl", "fast"]
    assert PS.main(argv) == 1
    report = json.loads((tmp_path / "profile.json").read_text())
    assert [report[a]["exit"] for a in ("trl", "fast")] == [1, 1] and len(calls) == 2


def test_busy_is_the_union_of_gpu_intervals_inside_the_traced_step():
    ev = [{"ph": "X", "name": "ProfilerStep#5", "cat": "user_annotation", "ts": 100, "dur": 1000},
          {"ph": "X", "name": "k1", "cat": "kernel", "ts": 50, "dur": 150},      # 100-200 inside
          {"ph": "X", "name": "k2", "cat": "kernel", "ts": 150, "dur": 100},     # overlaps k1: +50
          {"ph": "X", "name": "copy", "cat": "gpu_memcpy", "ts": 500, "dur": 100},
          {"ph": "X", "name": "k1", "cat": "kernel", "ts": 1050, "dur": 200},    # clipped at 1100: +50
          {"ph": "X", "name": "cudaLaunchKernel", "cat": "cuda_runtime", "ts": 120, "dur": 2000},
          {"ph": "X", "name": "aten::_local_scalar_dense", "cat": "cpu_op", "ts": 300, "dur": 5}]
    b = PS.busy({"traceEvents": ev})
    assert b["step_ms"] == 1.0 and b["gpu_busy"] == pytest.approx(0.3) and b["kernels"] == 3
    assert b["cuda_api_ms"] == 2.0 and list(b["top_kernels_ms"]) == ["k1", "k2"]
    assert b["syncs"]["aten::_local_scalar_dense"] == 1 and b["syncs"]["cudaStreamSynchronize"] == 0
    with pytest.raises(ValueError, match="ProfilerStep"):
        PS.busy({"traceEvents": ev[1:]})


def test_smi_mean_skips_unparsable_rows(tmp_path):
    p = tmp_path / "smi.csv"
    p.write_text("2026/11/06 10:00:00.000, 60, 200.5, 6900\n2026/11/06 10:00:00.500, [N/A], 0, 0\n"
                 "2026/11/06 10:00:01.000, 80, 250.5, 7000\n")
    assert PS.smi_mean(p) == {"util_mean": 70.0, "power_w_mean": 225.5, "mem_max_mib": 7000.0}
    p.write_text("")
    assert PS.smi_mean(p) == {}


@pytest.mark.slow
def test_train_sft_profile_writes_a_trace_of_one_step(tmp_path, synthetic_items):
    pytest.importorskip("trl")
    from test_p2_v1_train import _tiny_args

    prof = tmp_path / "prof"
    args = _tiny_args(sft_train, tmp_path, synthetic_items, branches=[], steps_per_epoch=7, profile=str(prof))
    assert sft_train.train(args) == 0
    with gzip.open(prof / "trace.json.gz", "rt") as f:
        b = PS.busy(json.load(f))
    assert b["step_ms"] > 0 and b["kernels"] == 0 and b["gpu_busy"] == 0  # CPU run: no CUDA activity
    assert "Self CPU" in (prof / "ops.txt").read_text()


def test_train_sft_refuses_profile_with_branches(tmp_path):
    pytest.importorskip("trl")
    from objdel.core.spec import SftStage

    with pytest.raises(ValueError, match="profile is for timing runs"):
        sft_train.train(SftStage(model="m", view="v", output_dir=str(tmp_path / "o"), vision=False,
                          profile=str(tmp_path / "p")))


def test_p1_chain_checks_the_card_and_the_checker_before_rl():
    from experiments.phase2.vast import policy

    steps = policy.chain(policy.ROOT / "experiments/phase2/configs/p1_rl.yaml", "p1")
    names = [" ".join(s) for s in steps]
    assert "env_check.sh /workspace/runs/p1/ckpt" in names[1] and "run_box.sh" in names[2]
    assert names[4].startswith("-m experiments.phase2.profile --ckpt /workspace/runs/p1/ckpt")
    assert [s[0] for s in steps[1:]] == ["bash", "bash", "~", "-m"]  # env and checker failures stop the chain
    assert not any("hf_sample" in n or "hf_rate" in n for n in names)  # vLLM only (owner 2026-11-06)


def test_memcheck_has_its_sanitizer_and_fails_loudly():
    """P1's memcheck step found no compute-sanitizer and passed silently. K1 (RTX 3090, 2026-11-06): the toolkit
    extra's 13.0.85 wheel lacks the TreeLauncher executables ("Target application terminated before first
    instrumented API call"); the 13.4.92 wheel ran memcheck and initcheck on 40 GPU tests, 0 errors."""
    from experiments.phase2.vast import policy

    box = (policy.ROOT / "experiments/cuda_checker/run_box.sh").read_text()
    assert '"cuda-toolkit[nvcc,cccl,nvtx]==13.0.3" "nvidia-cuda-sanitizer-api==13.4.92"' in box
    assert "sanitizer]" not in box
    mem = " ".join(policy.chain(policy.ROOT / "experiments/phase2/configs/p1_rl.yaml", "p1")[3])
    assert "--error-exitcode 1" in mem and "no compute-sanitizer" in mem
    # torch's caching allocator hides out-of-bounds accesses from memcheck; arena reuse is initcheck's case
    assert "PYTORCH_NO_CUDA_MEMORY_CACHING=1" in mem and "for t in memcheck initcheck" in mem
    assert "timeout 1200 $S --tool $t" in mem  # memcheck instruments torch's kernels too: each pass is bounded
    # R2 (2026-11-07): the step sanitized only test_cuda_checker_gpu.py (6 tests); K1's 40 are run_box.sh's group,
    # and the device scorer's kernels run only in test_gpu_scoring.py
    # the default (non --gpu-check) group: the last pytest call before "|| tests=1" (D8's --gpu-check calls come first)
    group = box.split("|| tests=1", 1)[0].rsplit("pytest -q -p no:cacheprovider", 1)[1]
    files = [*group.replace("\\", " ").split(), "tests/test_gpu_scoring.py"]
    assert len(files) == 5 and all(f in mem.split() for f in files), (files, mem)


def test_r2_chain_checks_the_card_and_checker_then_tunes_e2_with_the_preregistered_g_probe():
    """R2 (PLAN §5.1 G-probe, 20 GRPO steps at the pre-registered LR on R1's e2): P1's card and checker steps
    verbatim (so the two chains cannot drift), then g_probe on the fast trainer; vLLM only."""
    from experiments.phase2.vast import policy

    cfg = policy.ROOT / "experiments/phase2/configs"
    p1, r2 = (policy.chain(cfg / f, "r2") for f in ("p1_rl.yaml", "r2_qwen3vl2b.yaml"))
    assert r2[:4] == p1[:4]
    gp = r2[4]
    assert gp[:2] == ["-m", "experiments.phase2.g_probe"] and len(r2) == 5
    args = dict(zip(gp[2::2], gp[3::2], strict=True))
    data = "/workspace/data/lidc"
    assert args == {"--model": "qwen3-vl-2b", "--ckpt": "/workspace/runs/r2/ckpt", "--train": f"{data}/prepared",
                    "--probe": f"{data}/eval_views/lidc-val-shifted", "--out": "/workspace/runs/r2/g_probe"}
    assert not any(f in gp for f in ("--plain", "--grpo", "--items", "--select"))  # the pre-registered G-probe as is


def test_full_rl_chain_checks_the_card_and_checker_then_tunes_e2_resumably():
    """The full RL run: P1's card, checker and sanitizer steps verbatim, then train_rl on R1's e2, saving every 50
    steps and always with --resume, so rerunning the chain on a new box continues from the last checkpoint."""
    from experiments.phase2.vast import policy

    cfg = policy.ROOT / "experiments/phase2/configs"
    p1, rl = (policy.chain(cfg / f, "rl") for f in ("p1_rl.yaml", "rl_full.yaml"))
    assert rl[:4] == p1[:4] and len(rl) == 5
    tr = rl[4]
    assert tr[:2] == ["-m", "experiments.phase2.train_rl"] and tr[-1] == "--resume"
    args = dict(zip(tr[2:-1:2], tr[3:-1:2], strict=True))
    data = "/workspace/data/lidc"
    assert args == {"--model": "qwen3-vl-2b", "--ckpt": "/workspace/runs/rl/ckpt", "--train": f"{data}/prepared",
                    "--out": "/workspace/runs/rl/train_rl", "--steps": "770", "--save-steps": "50"}


def test_checker_gpu_chain_builds_then_sanitizes_before_the_long_checks():
    """A checker-only rental (no model, no vLLM, no data bundle): build for the card, memcheck + initcheck first (the
    step P1 never ran), then P1's equivalence and timing step if the cap leaves time."""
    from experiments.phase2.vast import policy

    cfg = policy.ROOT / "experiments/phase2/configs"
    p1, k = (policy.chain(cfg / f, "k1") for f in ("p1_rl.yaml", "checker_gpu.yaml"))
    assert k[0] == p1[0]  # box facts
    assert k[1] == ["bash", "-c", "mkdir -p /workspace/runs/k1/cuda_checker && "
                    "bash experiments/cuda_checker/run_box.sh --build-only"]  # fmt: skip
    assert k[2] == p1[3] and k[3] == ["~", *p1[2]] and len(k) == 4


def test_speed_chain_compares_vllm_asleep_and_awake_on_e2():
    """R4: the full run's card facts and env check, then two short train_rl arms on e2 with today's code: vLLM asleep
    between steps (the full run's setting), then kept awake. The awake arm may run out of memory: it must not stop the
    chain, and neither arm saves checkpoints."""
    import json

    from experiments.phase2.vast import policy

    cfg = policy.ROOT / "experiments/phase2/configs"
    full, r4 = (policy.chain(cfg / f, "r4") for f in ("rl_full.yaml", "r4_speed.yaml"))
    assert r4[:2] == full[:2] and len(r4) == 4
    asleep, awake = r4[2], r4[3]
    assert asleep[:2] == ["-m", "experiments.phase2.train_rl"]
    assert awake[:3] == ["~", "-m", "experiments.phase2.train_rl"]
    a, b = (dict(zip(s[i::2], s[i + 1 :: 2], strict=True)) for s, i in ((asleep, 2), (awake, 3)))
    assert a["--ckpt"] == b["--ckpt"] == "/workspace/runs/r4/ckpt" and a["--steps"] == b["--steps"] == "15"
    assert int(a["--save-steps"]) > 15 and int(b["--save-steps"]) > 15 and a["--out"] != b["--out"]
    assert "--grpo" not in a and json.loads(b["--grpo"]) == {"vllm_enable_sleep_mode": False}


def test_shared_prefix_chain_compares_the_packed_forward_on_e2():
    """R5: R4's box facts and env check, then two short train_rl arms on e2: today's forward, then --shared-prefix.
    The packed arm must not stop the chain, and neither arm saves checkpoints."""
    from experiments.phase2.vast import policy

    cfg = policy.ROOT / "experiments/phase2/configs"
    r4, r5 = (policy.chain(cfg / f, "r5") for f in ("r4_speed.yaml", "r5_shared.yaml"))
    assert r5[:2] == r4[:2] and len(r5) == 4
    out = r5[2].index("--out") + 1
    assert r5[2][:out] + r5[2][out + 1:] == r4[2][:out] + r4[2][out + 1:]  # R4's base arm, own folder
    base, packed = r5[2], r5[3]
    assert packed[:3] == ["~", "-m", "experiments.phase2.train_rl"] and packed[-1] == "--shared-prefix"
    assert "--shared-prefix" not in base
    a, b = (dict(zip(s[i::2], s[i + 1 :: 2], strict=False)) for s, i in ((base, 2), (packed[:-1], 3)))
    assert a["--ckpt"] == b["--ckpt"] == "/workspace/runs/r5/ckpt" and a["--steps"] == b["--steps"] == "15"
    assert int(a["--save-steps"]) > 15 and int(b["--save-steps"]) > 15 and a["--out"] != b["--out"]


def test_covered_hand_cases():
    """A6: the interval-union length clipped to [t0, t1], shared by busy (torch.profiler) and nsys_baseline."""
    c = PS.covered
    assert c([(0, 2), (1, 3), (5, 6)], 0, 10) == 4.0
    assert c([(0, 2), (1, 3), (5, 6)], 1.5, 5.5) == 2.0
    assert c([(0, 10), (2, 3)], 0, 10) == 10 and c([(0, 1), (1, 2)], 0, 10) == 2.0
    assert c([(-5, 2), (8, 20)], 0, 10) == 4
    for case in ([], [(3, 3)], [(20, 30)], [(-9, -1)]):
        assert c(case, 0, 10) == 0.0
    assert c([(0, 2)], 5, 5) == 0.0 and c([(0, 2)], 5, 1) == 0.0  # t1 <= t0: never negative
    assert c([(5, 6), (1, 3), (0, 2)], 0, 10) == c([(0, 2), (1, 3), (5, 6)], 0, 10)
    assert c(((a, a + 1) for a in (0, 4)), 0, 10) == 2.0  # an iterator, consumed once
    assert c([(3, 1)], 0, 10) == 0.0 and c([(3, 1), (0, 2)], 0, 10) == 2.0  # a reversed interval counts nothing (pin)


def test_covered_equals_the_rasterised_union():
    rng = np.random.default_rng(0)
    for _ in range(500):
        n = int(rng.integers(0, 8))
        a = rng.integers(-5, 25, n)
        iv = [(int(x), int(x + d)) for x, d in zip(a, rng.integers(0, 10, n), strict=True)]
        t0 = int(rng.integers(-2, 10))
        t1 = t0 + int(rng.integers(0, 15))
        grid = np.zeros(40, bool)
        for x, y in iv:
            grid[max(x, t0) + 10:max(min(y, t1), t0) + 10] = True
        got = PS.covered(iv, t0, t1)
        assert got == grid.sum() and 0 <= got <= t1 - t0
        half = len(iv) // 2
        A, B = PS.covered(iv[:half], t0, t1), PS.covered(iv[half:], t0, t1)
        assert max(A, B) <= got <= A + B
        assert PS.covered(iv + iv, t0, t1) == got and PS.covered(iv, t0, t1 + 3) >= got


def test_busy_golden_is_unchanged():
    """profile.busy before A6 moved its interval loop into covered (frozen outputs of 5 hand inputs)."""
    def ev(name, cat, ts, dur):
        return {"ph": "X", "name": name, "cat": cat, "ts": ts, "dur": dur}

    def S(ts, dur):
        return ev("ProfilerStep#3", "user_annotation", ts, dur)

    cases = [
        [S(0, 100), ev("k", "kernel", 10, 20), ev("k", "kernel", 20, 30)],
        [S(0, 100), ev("a", "kernel", -10, 30), ev("m", "gpu_memcpy", 90, 30),
         ev("cudaStreamSynchronize", "cuda_runtime", 5, 2)],
        [S(50, 50), S(100, 100), ev("k", "kernel", 60, 10), ev("s", "gpu_memset", 150, 5),
         ev("aten::nonzero", "cpu_op", 70, 1)],
        [S(0, 1000), ev("k1", "kernel", 0, 1000), ev("k2", "kernel", 100, 10), ev("k2", "kernel", 500, 10)],
        [S(0, 10), ev("x", "cpu_op", 0, 10)],
    ]  # fmt: skip
    zero = dict.fromkeys(PS.SYNCS, 0)
    want = [
        {"step_ms": 0.1, "gpu_busy": 0.4, "kernels": 2, "cuda_api_ms": 0.0, "syncs": zero,
         "top_kernels_ms": {"k": 0.05}},
        {"step_ms": 0.1, "gpu_busy": 0.3, "kernels": 1, "cuda_api_ms": 0.002,
         "syncs": zero | {"cudaStreamSynchronize": 1}, "top_kernels_ms": {"a": 0.03}},
        {"step_ms": 0.15, "gpu_busy": 0.1, "kernels": 1, "cuda_api_ms": 0.0, "syncs": zero | {"aten::nonzero": 1},
         "top_kernels_ms": {"k": 0.01}},
        {"step_ms": 1.0, "gpu_busy": 1.0, "kernels": 3, "cuda_api_ms": 0.0, "syncs": zero,
         "top_kernels_ms": {"k1": 1.0, "k2": 0.02}},
        {"step_ms": 0.01, "gpu_busy": 0.0, "kernels": 0, "cuda_api_ms": 0.0, "syncs": zero, "top_kernels_ms": {}},
    ]  # fmt: skip
    assert [PS.busy({"traceEvents": c}) for c in cases] == want

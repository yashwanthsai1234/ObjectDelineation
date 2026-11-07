"""profile: the RL arms (g_probe timing runs), their summary, the traced step's GPU busy share, nvidia-smi means;
the shared torch.profiler callback (train_sft `profile:`, g_probe --trace; a CPU trace; refused outside timing runs)."""

import gzip
import json

import pytest

from experiments.phase2 import profile as PS


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

    from experiments.phase2 import train_sft as TS

    prof = tmp_path / "prof"
    assert TS.train(_tiny_args(TS, tmp_path, synthetic_items, branches=[], steps_per_epoch=7, profile=str(prof))) == 0
    with gzip.open(prof / "trace.json.gz", "rt") as f:
        b = PS.busy(json.load(f))
    assert b["step_ms"] > 0 and b["kernels"] == 0 and b["gpu_busy"] == 0  # CPU run: no CUDA activity
    assert "Self CPU" in (prof / "ops.txt").read_text()


def test_train_sft_refuses_profile_with_branches(tmp_path):
    pytest.importorskip("trl")
    from experiments.phase2 import train_sft as TS

    with pytest.raises(ValueError, match="profile is for timing runs"):
        TS.train(TS.RunArgs("m", "v", str(tmp_path / "o"), vision=False, profile=str(tmp_path / "p")))


def test_p1_chain_checks_the_card_and_the_checker_before_rl():
    from experiments.phase2.vast import policy

    steps = policy.chain(policy.ROOT / "experiments/phase2/configs/p1_rl.yaml", "p1")
    names = [" ".join(s) for s in steps]
    assert "env_check.sh /workspace/runs/p1/ckpt" in names[1] and "run_box.sh" in names[2]
    assert names[4].startswith("-m experiments.phase2.profile --ckpt /workspace/runs/p1/ckpt")
    assert [s[0] for s in steps[1:]] == ["bash", "bash", "~", "-m"]  # env and checker failures stop the chain
    assert not any("hf_sample" in n or "hf_rate" in n for n in names)  # vLLM only (owner 2026-11-06)


def test_memcheck_has_its_sanitizer_and_fails_loudly():
    """P1's memcheck step found no compute-sanitizer (the toolkit install lacked the extra) and passed silently;
    run_box.sh now installs the toolkit's own sanitizer pin and the step reports a missing binary or any error."""
    from experiments.phase2.vast import policy

    box = (policy.ROOT / "experiments/cuda_checker/run_box.sh").read_text()
    assert '"cuda-toolkit[nvcc,cccl,sanitizer]==13.0.3"' in box
    mem = " ".join(policy.chain(policy.ROOT / "experiments/phase2/configs/p1_rl.yaml", "p1")[3])
    assert "--error-exitcode 1" in mem and "no compute-sanitizer" in mem


def test_r2_chain_checks_the_card_and_checker_then_tunes_e2_with_the_preregistered_g_probe():
    """R2 (PLAN §5.1 G-probe, 20 GRPO steps at the pre-registered LR on R1's e2): P1's card and checker steps
    verbatim (so the two chains cannot drift), then g_probe on the fast trainer; vLLM only."""
    from experiments.phase2.vast import policy

    cfg = policy.ROOT / "experiments/phase2/configs"
    p1, r2 = (policy.chain(cfg / f, "r2") for f in ("p1_rl.yaml", "r2_qwen3vl2b.yaml"))
    assert r2[:4] == p1[:4]
    gp = r2[4]
    assert gp[:2] == ["-m", "experiments.phase2.g_probe"] and len(r2) == 5
    args = dict(zip(gp[2::2], gp[3::2]))
    assert args == {"--model": "qwen3-vl-2b", "--ckpt": "/workspace/runs/r2/ckpt",
                    "--train": "/workspace/data/lidc/prepared", "--probe": "/workspace/data/lidc/eval_views/lidc-val-shifted",
                    "--out": "/workspace/runs/r2/g_probe"}  # fmt: skip
    assert not any(f in gp for f in ("--plain", "--grpo", "--items", "--select"))  # the pre-registered G-probe as is

"""profile_sft: arm configs (R1 as timing runs), the traced step's GPU busy share, nvidia-smi means, and train_sft's
profile option (a CPU trace; refused outside timing runs)."""

import gzip
import json

import pytest

from experiments.phase2 import profile_sft as PS


def test_arm_configs_are_r1_timing_runs(tmp_path):
    pytest.importorskip("yaml")
    base, arm = PS.arm_config("base", tmp_path), PS.arm_config("no_ckpt_mb8", tmp_path)
    assert base["branches"] == [] and base["vision"] == "off" and base["steps_per_epoch"] == PS.STEPS
    assert base["model"] == "qwen3-vl-2b" and base["profile"] == str(tmp_path / "base")
    assert "sampler" not in base and base["sft"]["per_device_train_batch_size"] == 4
    assert arm["sft"]["gradient_checkpointing"] is False and arm["sft"]["per_device_train_batch_size"] == 8
    assert arm["sft"]["gradient_accumulation_steps"] == 4  # the effective batch stays 32
    for a in PS.ARMS.values():
        assert a.get("per_device_train_batch_size", 4) * a.get("gradient_accumulation_steps", 8) == 32


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
    assert PS.smi_mean(p) == {"util_mean": 70.0, "power_w_mean": 225.5}
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

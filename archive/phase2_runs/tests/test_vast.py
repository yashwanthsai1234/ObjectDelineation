"""Archived by CLEAN8 from tests/test_vast.py (retired scripts and plans; nothing runs these)."""


R0 = ROOT / "experiments/phase2/configs/r0.yaml"


def test_r0_is_one_chain_of_the_plan_steps(E):
    """PLAN §4.7 as one run.sh chain for Qwen3-VL-2B: image check, x86 hashes, idle K3, sampler parity and HF sampling
    rate (a failure is recorded, not blocking), the G-vision reference, the 30-step timing arm, then the CUDA checker
    on the GPU (last, not blocking: `make cuda` reinstalls objdel)."""
    steps = P.chain(R0, "v1")
    assert steps[0][:3] == ["~", "bash", "-c"] and "nvidia-smi -q" in steps[0][3] and "cpu.max" in steps[0][3]  # step 2
    assert steps[1] == ["bash", "experiments/phase2/docker/image_check.sh", "/workspace/runs/v1/image_check.json"]
    mods = [after(s, "-m") for s in steps if "-m" in s]
    assert mods[:2] == ["experiments.phase1.bench.xarch_hashes", "experiments.phase2.k3_harness"]
    assert mods.count("experiments.phase2.r0_parity") == 1 and mods.count("experiments.phase2.train_sft") == 2
    k3 = next(s for s in steps if "experiments.phase2.k3_harness" in s)
    assert after(k3, "--root") == "/workspace/data/lidc/prepared"
    cfgs = [after(s, "--config") for s in steps if "--config" in s]
    assert all((ROOT / c).exists() for c in cfgs) and len(cfgs) == 2
    refs = [after(s, "--config") for s in steps if "--vision_reference_only" in s]
    assert refs == ["experiments/phase2/configs/r1_qwen3vl2b.yaml"]
    assert cfgs[1] == "experiments/phase2/configs/r0_qwen3vl_arm.yaml"
    E.launched()
    line = E.tmux(R0.relative_to(ROOT))
    assert "bash experiments/phase2/docker/image_check.sh /workspace/runs/r1/image_check.json" in line
    assert "/opt/venv/bin/python -I -m experiments.phase1.bench.xarch_hashes" in line
    assert "DATA=/workspace/data/lidc/checker_eval bash experiments/cuda_checker/run_box.sh" in steps[-1][3]
    assert steps[-1][0] == "~" and "/workspace/runs/v1/cuda_checker" in steps[-1][3]
    assert line.count("|| echo") == 4  # box record, parity, HF rate and the CUDA checker go on failing; the rest stop

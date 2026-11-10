"""E0: move-invariance goldens. A file move or a split must not change one bit of what these compute. Recorded once on
the E0 commit (not on an older one: phases A-D change rewards, targets and the trainer on purpose). The recording
script is the `record` function here; run it as in tests/goldens/phase_e/README.

Platform-independent values are in invariance.json (goldens 1, 2, 5, 6). Floats that depend on the compiler or the CPU
kernels (goldens 3, 4) are in invariance.<sys.platform>-<machine>.json, and a platform without a file skips."""
import hashlib
import json
import platform
import sys
from dataclasses import astuple
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
GOLDENS = ROOT / "tests/goldens/phase_e"
PLATFORM = f"{sys.platform}-{platform.machine()}"
# E6.2: every plan that is live on the E0 commit and goes through policy.chain.
PLANS = ("b1_rows", "rl_full", "r5_shared", "r1_qwen3vl2b")
# Fields of GRPOConfig.to_dict() that hold a path, a time or a host name: not part of the run's behaviour.
VOLATILE = {"output_dir", "logging_dir", "run_name", "hub_model_id", "hub_token", "push_to_hub_model_id",
            "push_to_hub_organization", "push_to_hub_token", "local_rank", "_n_gpu", "world_size", "process_index"}


def reference_value() -> str:
    """Golden 1: the value part of reward.spec.fingerprint(), without library versions."""
    from objdel.geom.topology import min_area
    from objdel.verify.reward.spec import REWARD_V4

    value = (astuple(REWARD_V4), [min_area(a) for a in (0, 999, 1001, 10**4, 10**6)])
    return hashlib.sha256(repr(value).encode()).hexdigest()


def target_hash() -> str:
    """Golden 2: target_set_hash of the 300 fixture targets of tests/data/test_targets_golden.py."""
    from objdel.data.targets.build import build_target, target_set_hash
    from tests.data.test_targets_golden import fixtures

    rows = []
    for i, (masks, names, ma) in enumerate(fixtures()):
        try:
            rows.append((str(i), "-", build_target(masks, names, min_area=ma).text))
        except ValueError as e:
            rows.append((str(i), "-", f"ValueError: {e}"))
    return target_set_hash(rows)


def _box_text(rng, names) -> str:
    def ring():
        x0, y0 = (int(v) for v in rng.integers(0, 700, 2))
        x1, y1 = x0 + int(rng.integers(50, 300)), y0 + int(rng.integers(50, 300))
        return {"outer": [[x0, y0], [x1, y0], [x1, y1], [x0, y1]], "holes": []}

    return json.dumps({n: [ring()] for n in names}, separators=(",", ":"))


def score_corpus() -> dict[str, list[str]]:
    """Golden 3: score_batch on 64 polygon texts (box rings, some broken) against fixture caches, each float as hex."""
    from objdel.verify.reward.cache import ExpertCache
    from objdel.verify.reward.score import NUMERIC, score_batch

    def stack(boxes, W=50):
        m = np.zeros((len(boxes), W, W), bool)
        for i, (a, b, c, d) in enumerate(boxes):
            m[i, b:d, a:c] = True
        return m

    caches = {
        "a": {"nodule": ExpertCache.build(stack([(10, 10, 30, 30), (11, 11, 31, 31), (10, 10, 32, 32)]))},
        "r": {"disc": ExpertCache.build(stack([(5, 5, 45, 45), (6, 5, 45, 44), (5, 6, 44, 45)])),
              "cup": ExpertCache.build(stack([(18, 18, 32, 32), (19, 18, 31, 32), (18, 19, 32, 31)]))},
    }
    names = {"a": ("nodule",), "r": ("disc", "cup")}
    rng, jobs = np.random.default_rng(20261109), []
    for i in range(64):
        k = "a" if i % 2 == 0 else "r"
        text = {6: "garbage", 7: '{"nodule":[]}'}.get(i % 8) or _box_text(rng, names[k])
        jobs.append((text, bool(i % 16 == 15), names[k], caches[k], 50, 50))
    out = score_batch(jobs)
    return {k: [float(v).hex() for v in out[k]] for k in NUMERIC} | {"error": out["error"]}


def update_norms(monkeypatch, tmp_path) -> dict[str, str]:
    """Golden 4: one FastGRPO update on the tiny offline model of tests/experiments/test_trl_adapter_equivalence.py
    (CPU, seed 0): the loss, the gradient L2 norm and the L2 norm of the parameter change, as float hex."""
    import torch

    from objdel.rl import trl_adapter as fast
    from tests.experiments.test_trl_adapter_equivalence import tiny_loop_vllm_fake

    before, train = {}, fast.FastGRPO.train

    def snapshot_then_train(self, *a, **k):
        before.update({n: p.detach().clone() for n, p in self.model.named_parameters()})
        return train(self, *a, **k)

    monkeypatch.setattr(fast.FastGRPO, "train", snapshot_then_train)
    t = tiny_loop_vllm_fake(monkeypatch, tmp_path, steps=1)
    log = next(h for h in t.state.log_history if "loss" in h)
    change = torch.sqrt(sum(((p.detach() - before[n]) ** 2).sum() for n, p in t.model.named_parameters()))
    return {"loss": float(log["loss"]).hex(), "grad_norm": float(log["grad_norm"]).hex(),
            "param_change_norm": float(change).hex()}


def chain_argv() -> dict[str, dict]:
    """Golden 5: C.chain() of each live plan as an ordered list of argv lists (stage index = position, a
    leading "~" kept), and the G-vision references the plan reads."""
    from objdel.infra.vast import chain as C

    cfgs = ROOT / "specs"  # E6.2: the run specs replaced the plan files
    return {p: {"argv": C.chain(cfgs / f"{p}.yaml", "x"),
                "vision_refs": C.vision_refs(cfgs / f"{p}.yaml", "x")} for p in PLANS}


class _Captured(Exception):
    pass


def grpo_config_dict(argv: list[str], **over) -> dict:
    """The GRPOConfig that build_trainer resolves for train_rl's stage of `argv` (objdel.rl.run.checked of the
    RlStage), as to_dict() without VOLATILE keys."""
    from objdel.core.spec import RlStage, cli_stage
    from objdel.rl.build import build_trainer
    from objdel.rl.run import checked

    budget = ["--session-start", "0", "--rate-usd-per-h", "0.5"] if "--budget-usd" in argv else []  # as run.sh
    a, seen = checked(cli_stage(RlStage, [*argv, *budget])), []

    def stop(cfg):
        seen.append(cfg)
        raise _Captured

    with pytest.raises(_Captured):
        build_trainer(a.model, Path(a.ckpt), None, Path(a.out), [], a.steps or 20, a, check=stop, **over)
    return {k: v for k, v in seen[0].to_dict().items() if k not in VOLATILE}


def configs_golden() -> dict[str, dict]:
    """Golden 6: the resolved GRPOConfig of (a) the train_rl stage of b1_rows.yaml, (b) a G-probe with default
    arguments; both on qwen3-vl-2b. (a) adds what train_rl.main adds for a schedule root."""
    from objdel.infra.vast import chain as C

    steps = C.chain(ROOT / "specs/b1_rows.yaml", "x")
    stage = next(s for s in steps if "experiments.phase2.train_rl" in s)
    argv = stage[stage.index("experiments.phase2.train_rl") + 1:]
    a_over = {"save_strategy": "steps", "save_steps": 50, "shuffle_dataset": False}
    probe = ["--model", "qwen3-vl-2b", "--ckpt", "c", "--train", "t", "--out", "o", "--steps", "20"]
    return {"train_rl_b1_rows": grpo_config_dict(argv, **a_over), "g_probe_default": grpo_config_dict(probe)}


def jsonable(d):
    return json.loads(json.dumps(d, default=str, sort_keys=True))


def record() -> None:
    """Write the golden files. Run on the E0 commit only (see tests/goldens/phase_e/README)."""
    from _pytest.monkeypatch import MonkeyPatch

    GOLDENS.mkdir(parents=True, exist_ok=True)
    import cv2
    import scipy

    from objdel.verify.reward.spec import fingerprint

    flat = {"reference_value": reference_value(), "target_set_hash": target_hash(), "chain_argv": chain_argv(),
            "dropped_flags": golden().get("dropped_flags", {}), "renamed_paths": golden().get("renamed_paths", {}),
            "grpo_config": jsonable(configs_golden()),
            "fingerprint": {"value": fingerprint(), "versions": [np.__version__, cv2.__version__, scipy.__version__]}}
    (GOLDENS / "invariance.json").write_text(json.dumps(flat, indent=1, sort_keys=True) + "\n")
    import tempfile

    with MonkeyPatch.context() as mp, tempfile.TemporaryDirectory() as tmp:
        per = {"score_batch": score_corpus(), "update_norms": update_norms(mp, Path(tmp))}
    (GOLDENS / f"invariance.{PLATFORM}.json").write_text(json.dumps(per, indent=1, sort_keys=True) + "\n")


def golden(name: str = "invariance.json") -> dict:
    return json.loads((GOLDENS / name).read_text())


def later_edits(argvs: list[list[str]]) -> list[list[str]]:
    """Golden 5's argv lists with the edits of later commits, applied like a module rename: the test files a later
    commit moved get their new path (golden()["renamed_paths"]: {old: new}; E5 splits test_train_fast.py), and the
    flags a later commit deleted go, each with its value (golden()["dropped_flags"]: {module: [flag, ...]}; E4 step 3
    deletes final_eval's --train)."""
    dropped, renamed, out = golden()["dropped_flags"], golden()["renamed_paths"], []
    for argv in argvs:
        a = list(argv)
        for old, new in renamed.items():
            a = [x.replace(old, new) for x in a]
        for flag in next((f for m, f in dropped.items() if m in a), ()):
            while flag in a:
                del a[a.index(flag):a.index(flag) + 2]
        out.append(a)
    return out


def test_golden_1_reward_value_part_of_the_fingerprint():
    assert reference_value() == golden()["reference_value"]


def test_golden_1_full_fingerprint_when_the_libraries_match():
    import cv2
    import scipy

    from objdel.verify.reward.spec import fingerprint

    rec = golden().get("fingerprint")
    if rec is None or rec["versions"] != [np.__version__, cv2.__version__, scipy.__version__]:
        pytest.skip("numpy, OpenCV or scipy differ from the recorded versions")
    assert fingerprint() == rec["value"]


def test_golden_2_target_set_hash():
    assert target_hash() == golden()["target_set_hash"]


def test_golden_3_score_batch_bits():
    f = GOLDENS / f"invariance.{PLATFORM}.json"
    if not f.exists():
        pytest.skip(f"no golden recorded for {PLATFORM}")
    assert score_corpus() == golden(f.name)["score_batch"]


def test_golden_4_one_update(monkeypatch, tmp_path):
    f = GOLDENS / f"invariance.{PLATFORM}.json"
    if not f.exists():
        pytest.skip(f"no golden recorded for {PLATFORM}")
    pytest.importorskip("trl")
    assert update_norms(monkeypatch, tmp_path) == golden(f.name)["update_norms"]


def test_golden_5_chain_argv_of_the_live_plans():
    got, want = json.loads(json.dumps(jsonable(chain_argv()))), json.loads(json.dumps(golden()["chain_argv"]))
    for g in want.values():
        g["argv"] = later_edits(g["argv"])
    for g in (got, want):  # r1 train: `--config <yaml>` before E6, SftStage flags after (test_chain compares by value)
        del g["r1_qwen3vl2b"]["argv"][1]
    assert got == want


def test_golden_6_resolved_grpo_config():
    pytest.importorskip("trl")
    assert jsonable(configs_golden()) == golden()["grpo_config"]


def test_golden_6_sees_the_learning_rate():
    assert golden()["grpo_config"]["train_rl_b1_rows"]["learning_rate"] == 5e-5
    assert golden()["grpo_config"]["g_probe_default"]["learning_rate"] == 5e-5

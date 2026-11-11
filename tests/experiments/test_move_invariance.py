"""E0: move-invariance goldens. A file move or a split must not change one bit of what these compute. Recorded once on
the E0 commit (not on an older one: phases A-D change rewards, targets and the trainer on purpose). The recording
script is the `record` function here; run it as in tests/goldens/phase_e/README.

Platform-independent values are in invariance.json (goldens 1, 2, 5, 6). Floats that depend on the compiler or the CPU
kernels (goldens 3, 4) are in invariance.<sys.platform>-<machine>.json, and a platform without a file skips."""
import hashlib
import json
import platform
import re
import sys
from dataclasses import astuple
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
GOLDENS = ROOT / "tests/goldens/phase_e"
PLATFORM = f"{sys.platform}-{platform.machine()}"
# E6.2: every plan that is live and goes through policy.chain; the goldens also hold plans that left the repo
PLANS = ("b1_rows", "r1_qwen3vl2b")
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


def _ring(outer, holes=()) -> dict:
    return {"outer": outer, "holes": list(holes)}


def _sq(x0, y0, x1, y1) -> list:
    """A clockwise (on screen) axis-aligned square ring; [::-1] gives a hole's winding."""
    return [[x0, y0], [x1, y0], [x1, y1], [x0, y1]]


def _nod(*rings) -> str:
    return json.dumps({"nodule": list(rings)}, separators=(",", ":"))


_CIRCLE = [[round(500 + 300 * np.cos(-2 * np.pi * k / 10_000)), round(500 + 300 * np.sin(-2 * np.pi * k / 10_000))]
           for k in range(10_000)]
# I12 (golden 3 must not be vacuous): the hard cases of the plan, each a named category. They take the first slots
# of the corpus whose text is a random box (i % 8 not 6 or 7), in this order; the other slots keep their bits.
HARD = {
    "empty_outline": ("a", _nod()),
    "single_pixel_box": ("a", _nod(_ring(_sq(400, 400, 420, 420)))),  # 20 grid units = 1 pixel of 50
    "touches_all_borders": ("a", _nod(_ring(_sq(0, 0, 1000, 1000)))),
    "collinear_zero_area": ("a", _nod(_ring([[100, 100], [500, 100], [900, 100]]))),
    "bow_tie": ("a", _nod(_ring([[200, 200], [600, 600], [600, 200], [200, 600]]))),
    "duplicate_vertices": ("a", _nod(_ring([[200, 200], [200, 200], [600, 200], [600, 600], [200, 600]]))),
    "two_distinct_points": ("a", _nod(_ring([[200, 200], [200, 200], [600, 600]]))),
    "vertex_pinch": ("a", _nod(_ring([[200, 200], [400, 400], [600, 200], [600, 600], [400, 400], [200, 600]]))),
    "spike": ("a", _nod(_ring([[200, 200], [600, 200], [600, 400], [800, 400], [600, 400], [600, 600],
                               [200, 600]]))),
    "collinear_overlap": ("a", _nod(_ring([[200, 200], [600, 200], [600, 600], [500, 600], [500, 200], [300, 200],
                                           [300, 600], [200, 600]]))),
    "hole": ("a", _nod(_ring(_sq(200, 200, 600, 600), [_sq(300, 300, 400, 400)[::-1]]))),
    "hole_touching_outer": ("a", _nod(_ring(_sq(200, 200, 600, 600), [_sq(200, 300, 400, 400)[::-1]]))),
    "hole_outside_outer": ("a", _nod(_ring(_sq(200, 200, 600, 600), [_sq(700, 700, 800, 800)[::-1]]))),
    "nested_holes": ("a", _nod(_ring(_sq(100, 100, 900, 900), [_sq(200, 200, 800, 800)[::-1],
                                                                _sq(300, 300, 700, 700)[::-1]]))),
    "two_nodules": ("a", _nod(_ring(_sq(100, 100, 300, 300)), _ring(_sq(500, 500, 700, 700)))),
    "nested_outers": ("a", _nod(_ring(_sq(100, 100, 900, 900)), _ring(_sq(300, 300, 700, 700)))),
    "rings_crossing": ("a", _nod(_ring(_sq(100, 100, 500, 500)), _ring(_sq(300, 300, 700, 700)))),
    "rings_touching": ("a", _nod(_ring(_sq(100, 100, 500, 500)), _ring(_sq(500, 100, 900, 500)))),
    "wrong_winding": ("a", _nod(_ring(_sq(200, 200, 600, 600)[::-1]))),
    "coordinate_1001": ("a", _nod(_ring(_sq(200, 200, 1001, 600)))),
    "coordinate_minus_1": ("a", _nod(_ring(_sq(-1, 200, 600, 600)))),
    "coordinate_float": ("a", _nod(_ring([[200.5, 200], [600, 200], [600, 600], [200, 600]]))),
    "coordinate_string": ("a", _nod(_ring([["200", 200], [600, 200], [600, 600], [200, 600]]))),
    "truncated_json": ("a", _nod(_ring(_sq(200, 200, 600, 600)))[:-7]),
    "non_json": ("a", "the nodule is in the centre"),
    "text_1_5_mb": ("a", _nod(_ring(_sq(200, 200, 600, 600))) + " " * 1_500_000),
    "vertices_10000": ("a", _nod(_ring(_CIRCLE))),
    "cup_outside_disc": ("r", json.dumps({"disc": [_ring(_sq(100, 100, 500, 500))],
                                          "cup": [_ring(_sq(600, 600, 800, 800))]}, separators=(",", ":"))),
}


def corpus_jobs() -> tuple[list[str], list[tuple]]:
    """Golden 3's 64 score_batch jobs against fixture caches and the category of each: "box" (random box rings),
    "garbage", "keys" (a nodule text on the disc/cup schema) or a HARD name."""
    from objdel.verify.reward.cache import ExpertCache

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
    rng, jobs, cats, hard = np.random.default_rng(20261109), [], [], iter(HARD.items())
    for i in range(64):
        k = "a" if i % 2 == 0 else "r"
        cat = {6: "garbage", 7: "keys"}.get(i % 8, "box")
        text = {6: "garbage", 7: '{"nodule":[]}'}.get(i % 8) or _box_text(rng, names[k])  # rng drawn as recorded
        if cat == "box" and (h := next(hard, None)):
            cat, (k, text) = h
        jobs.append((text, bool(i % 16 == 15), names[k], caches[k], 50, 50))
        cats.append(cat)
    return cats, jobs


def score_corpus() -> dict[str, list[str]]:
    """Golden 3: score_batch on the 64 texts of corpus_jobs, each float as hex."""
    from objdel.verify.reward.score import NUMERIC, score_batch

    out = score_batch(corpus_jobs()[1])
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
    from objdel.core.spec import MODULE
    from objdel.infra.vast import chain as C

    steps = C.chain(ROOT / "specs/b1_rows.yaml", "x")
    stage = next(s for s in steps if MODULE["rl"] in s)
    argv = stage[stage.index(MODULE["rl"]) + 1:]
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


# (old pin, new pin) of the rows root swap in later_edits: HF revision, then manifest SHA-256 (owner, 2026-11-10)
ROWS_780 = (("34213dceb068128f1a1abf9c67aa65dd59bc1220", "d2c8edb8e15534363ebbced8f839cecace49c495"),
            ("472738141f47bb0f1bdb1e98c679076cb73b979a4b8e256d7e2051c63e352fb4",
             "a0f2c1d212b8c685e2b85aaac1aa9fa3de9c004ab310af13412d5aa630d4a6b9"))
# {old: new} of the modules and box scripts that left experiments/ (owner decision D4, 2026-11-11: the command modules
# into objdel.cli, run as `python -I -m objdel.cli.<group>.<command>`; the tool scripts into tools/, the box scripts
# into ops/)
MODULE_MOVES = {"experiments.phase2.train_rl": "objdel.cli.train.rl",
                "experiments.phase2.train_sft": "objdel.cli.train.sft",
                "experiments.phase2.select_ckpt": "objdel.cli.train.select",
                "experiments.phase2.g_probe": "objdel.cli.train.probe",
                "experiments.phase2.final_eval": "objdel.cli.eval.run",
                "experiments.phase2.c_speed_check": "tools.studies.c_speed_check",
                "experiments/phase2/docker/nsys_capture.sh": "tools/profiling/nsys_capture.sh",
                "experiments/phase2/docker/env_check.sh": "ops/docker/env_check.sh",
                "experiments/cuda_checker/run_box.sh": "ops/box/run_box.sh"}


def later_edits(argvs: list[list[str]]) -> list[list[str]]:
    """Golden 5's argv lists with the edits of later commits, applied like a module rename: the test files a later
    commit moved get their new path (golden()["renamed_paths"]: {old: new}; E5 splits test_train_fast.py), and the
    flags a later commit deleted go, each with its value (golden()["dropped_flags"]: {module: [flag, ...]}; E4 step 3
    deletes final_eval's --train), and the moved modules and scripts get their new name last (MODULE_MOVES)."""
    dropped, renamed, out = golden()["dropped_flags"], golden()["renamed_paths"], []
    for argv in argvs:
        a = list(argv)
        for old, new in renamed.items():
            a = [x.replace(old, new) for x in a]
        for flag in next((f for m, f in dropped.items() if m in a), ()):
            while flag in a:
                del a[a.index(flag):a.index(flag) + 2]
        # E2.7 (owner, 2026-11-10): a speed arm saves once at its last step, not "never" (--save-steps 1000), and
        # c0_base got the 30 min nsys limit of the other arms
        if "--save-steps" in a and "--steps" in a and a[a.index("--save-steps") + 1] == "1000":
            a[a.index("--save-steps") + 1] = a[a.index("--steps") + 1]
        a = [re.sub(r"--steps (\d+)(\s+)--save-steps 1000", r"--steps \1\2--save-steps \1", x) for x in a]
        a = [re.sub(r"^(bash experiments/phase2/docker/nsys_capture\.sh \S+ c0_base )", r"NSYS_TIMEOUT=1800 \1", x)
             for x in a]
        # the 780-step rows root (owner, 2026-11-10): every rl stage pins rows/lidc-train-b1-780 at its HF revision and
        # manifest SHA-256 (the 260-step root is no longer run); the full RL run (b1-rows) trains 780 steps with
        # BudgetStop $8.00 of the owner's $10
        for was, now in ROWS_780:
            a = [x.replace(was, now) for x in a]
        a = [re.sub(r"lidc-train-b1(?![-\w])", "lidc-train-b1-780", x) for x in a]
        if "/workspace/runs/x/b1-rows" in a and "--steps" in a and "--budget-usd" in a:
            a[a.index("--steps") + 1] = "780"
            a[a.index("--budget-usd") + 1] = "8.0"
        for was, now in MODULE_MOVES.items():
            a = [x.replace(was, now) for x in a]
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
    """Every float equal by hex, the key set equal. A platform without a file skips, except darwin-arm64 (the
    binding gate), where a missing file fails."""
    f = GOLDENS / f"invariance.{PLATFORM}.json"
    if not f.exists() and PLATFORM != "darwin-arm64":
        pytest.skip(f"no golden recorded for {PLATFORM}")
    want, got = golden(f.name)["score_batch"], score_corpus()
    assert set(got) == set(want)
    assert got == want


def test_golden_3_corpus_is_not_vacuous():
    """I12: every hard category is in the corpus, and the union of the checker defects covers every class the text
    can reach. out_of_range cannot be reached: the parser refuses a coordinate outside 0-1000 first (error "range")."""
    from objdel.geom import engine
    from objdel.verify.reward.score import score_batch

    cats, jobs = corpus_jobs()
    assert len(jobs) == 64 and set(HARD) <= set(cats) and {"box", "garbage", "keys"} <= set(cats)
    out = score_batch(jobs)
    seen = {d for defects in out["defects"] for d, n in defects.items() if n}
    assert seen == set(engine.DEFECTS) - {"out_of_range"}
    errors = dict(zip(cats, out["error"], strict=True))
    assert (errors["coordinate_1001"], errors["coordinate_minus_1"]) == ("range", "range")
    assert errors["text_1_5_mb"] == errors["vertices_10000"] == "too_long"
    assert errors["truncated_json"] == "syntax" and errors["coordinate_float"] == "not_int"


def test_golden_4_one_update(monkeypatch, tmp_path):
    f = GOLDENS / f"invariance.{PLATFORM}.json"
    if not f.exists():
        pytest.skip(f"no golden recorded for {PLATFORM}")
    pytest.importorskip("trl")
    assert update_norms(monkeypatch, tmp_path) == golden(f.name)["update_norms"]


def test_golden_5_chain_argv_of_the_live_plans():
    got, want = json.loads(json.dumps(jsonable(chain_argv()))), json.loads(json.dumps(golden()["chain_argv"]))
    want = {p: want[p] for p in PLANS}  # the golden also holds the retired plans rl_full and r5_shared
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

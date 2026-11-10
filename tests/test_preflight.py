"""The preflight driver on a run spec (E6.8): stamp, staleness (refusal), fail-fast by stage, the inputs each check
reads. The checks themselves run on data/."""

import ast
import hashlib
import importlib.util
import json
import types
from pathlib import Path

import pytest
import yaml

from objdel.infra import preflight as PF

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = b'{"schedule": {}}\n'
CMD = {"kind": "cmd", "name": "box", "argv": ["bash", "-c", "true"]}


def tiny_data(tmp: Path) -> Path:
    """A data root with a training view, prepared and the exam views' and checker eval set's manifests: enough to
    hash."""
    data = tmp / "data"
    for d in ("views/v", "checker_eval", "prepared", *(f"eval_views/{v}" for v in PF.EXAMS)):
        (data / d).mkdir(parents=True)
        (data / d / "manifest.json").write_text(json.dumps({"view_hash": "h", "dir": d}))
    return data


def tiny_spec(tmp: Path, stages=(CMD,), gpu="RTX_5090", cap_hours=3.0, **top) -> Path:
    spec = {"schema_version": 1, "name": "pf_run", "compute": {"gpu": gpu, "cap_hours": cap_hours},
            "bundle_view": "views/v", "stages": list(stages)} | top
    path = tmp / "pf_run.yaml"
    path.write_text(yaml.safe_dump(spec, sort_keys=False))
    return path


@pytest.fixture
def cfg(tmp_path, monkeypatch):
    monkeypatch.setenv("OBJDEL_RUN_DIR", str(tmp_path / "runs"))
    monkeypatch.setenv("OBJDEL_DATA_ROOT", str(tiny_data(tmp_path)))
    return tiny_spec(tmp_path)


def run(cfg, monkeypatch, checks):
    monkeypatch.setattr(PF, "CHECKS", checks)
    code = PF.run(cfg)
    return code, json.loads(PF.stamp_path(cfg).read_text())


def test_a_passing_run_writes_a_stamp_the_gate_accepts_until_something_changes(cfg, monkeypatch):
    code, stamp = run(cfg, monkeypatch, [("DATA a", lambda s: []), ("CODE b", lambda s: [])])
    assert code == 0 and stamp["pass"] and set(stamp["hashes"]) == {"code", "data", "config"}
    assert stamp["config"] == str(cfg.resolve()) and stamp["hashes"]["config"] == hashlib.sha256(
        cfg.read_bytes()).hexdigest()  # the stamp holds the spec path and its sha256
    assert PF.refusal(cfg) is None
    text = cfg.read_text()
    cfg.write_text(text + " \n")
    assert "stale (config" in PF.refusal(cfg)
    (cfg.parent / "data/eval_views/probe-topology/manifest.json").write_text("{}")
    assert "stale (config, data" in PF.refusal(cfg)


def test_the_stamp_follows_the_spec_bytes_and_path(cfg, monkeypatch, tmp_path):
    """I15: reverted bytes pass, an mtime-only touch passes, the same bytes at another path are refused, a truncated
    stamp and a missing spec are refusals, never a crash."""
    run(cfg, monkeypatch, [("DATA a", lambda s: [])])
    text = cfg.read_text()
    cfg.write_text(text + "\n")
    assert "stale (config" in PF.refusal(cfg)
    cfg.write_text(text)
    cfg.touch()
    assert PF.refusal(cfg) is None
    other = tmp_path / "elsewhere" / cfg.name
    other.parent.mkdir()
    other.write_text(text)
    assert "another spec path" in PF.refusal(other)
    stamp = PF.stamp_path(cfg)
    stamp.write_text(stamp.read_text()[:20])
    assert "cannot check" in PF.refusal(cfg)
    run(cfg, monkeypatch, [("DATA a", lambda s: [])])
    cfg.unlink()
    assert "cannot check" in PF.refusal(cfg)


def test_a_changed_checker_eval_set_makes_the_stamp_stale(cfg, monkeypatch):
    run(cfg, monkeypatch, [("DATA a", lambda s: [])])
    assert PF.refusal(cfg) is None
    (cfg.parent / "data/checker_eval/manifest.json").write_text("{}")
    assert "stale (data" in PF.refusal(cfg)


def test_a_failure_finishes_its_stage_and_skips_later_ones(cfg, monkeypatch):
    ran = []
    checks = [("DATA a", lambda s: ["bad a"]), ("DATA b", lambda s: ran.append("b") or []),
              ("CHECKER c", lambda s: ran.append("c") or [])]  # fmt: skip
    code, stamp = run(cfg, monkeypatch, checks)
    assert code == 1 and not stamp["pass"] and ran == ["b"] and list(stamp["checks"]) == ["DATA a", "DATA b"]
    assert "FAILED" in PF.refusal(cfg)


def test_no_stamp_is_a_refusal(cfg):
    assert "no preflight stamp" in PF.refusal(cfg)


def test_a_crashing_check_fails_the_preflight(cfg, monkeypatch):
    def boom(s):
        raise KeyError("qwen3-vl-2b")
    code, stamp = run(cfg, monkeypatch, [("CODE boom", boom)])
    assert code == 1 and not stamp["pass"] and "check crashed: KeyError" in stamp["checks"]["CODE boom"]["problems"][0]


def test_every_check_resolves_its_lazy_imports():
    """E6.6: no check imports a module that is gone (projection imported the archived project_r1)."""
    tree = ast.parse(Path(PF.__file__).read_text())
    mods = {n.module for f in tree.body if isinstance(f, ast.FunctionDef) for n in ast.walk(f)
            if isinstance(n, ast.ImportFrom) and n.module}
    assert mods and [m for m in sorted(mods) if importlib.util.find_spec(m) is None] == []
    assert not hasattr(PF, "projection") and not hasattr(PF, "ceiling") and not hasattr(PF, "HF_RATE")
    assert "CODE projected GPU hours and dollars" not in dict(PF.CHECKS)


# ------------------------------------------------------------------------------- which input each check reads (I15)
def test_local_maps_data_templates_on_the_mac(monkeypatch):
    monkeypatch.delenv("OBJDEL_DATA_ROOT", raising=False)
    assert PF.local("{data}/rows/x") == Path("data/views/x")
    assert PF.local("{data}/eval_views/v") == Path("data/eval_views/v")
    assert PF.local("{rd}/ckpt") is None
    for bad in ("/abs", "{data}/../x", "{x}/y", "{data}/a/{rd}"):
        with pytest.raises(ValueError):
            PF.local(bad)
    monkeypatch.setenv("OBJDEL_DATA_ROOT", "/d")
    assert PF.local("{data}/rows/x") == Path("/d/views/x") and PF.local("{data}/prepared") == Path("/d/prepared")


def test_local_rows_root_maps_only_the_literal_box_path():
    assert PF.local_rows_root("{data}/rows/lidc-train-b1") == Path("data/views/lidc-train-b1")
    assert PF.local_rows_root("{data}/rows/x", "/d") == Path("/d/views/x")
    for bad in ("{data}/rows/../x", "{data}/views/lidc-train-b1", "/workspace/data/rows/x", "{data}/rows/a/b",
                "{data}/rows/.."):
        with pytest.raises(ValueError):
            PF.local_rows_root(bad)


def b1_rows_copy(tmp_path, sha) -> Path:
    """specs/b1_rows.yaml with the rows pin of the fixture root."""
    raw = yaml.safe_load((ROOT / "specs/b1_rows.yaml").read_text())
    raw["stages"][3]["rows_manifest_sha256"] = sha
    path = tmp_path / "b1_rows.yaml"
    path.write_text(yaml.safe_dump(raw, sort_keys=False))
    return path


@pytest.fixture
def rows_cfg(tmp_path, monkeypatch):
    """specs/b1_rows.yaml over fixture data: its train_rl stage names {data}/rows/lidc-train-b1, built locally at
    <data root>/views/lidc-train-b1 (one image and the manifest); bundle_view views/lidc-bundle-v3."""
    data = tiny_data(tmp_path)
    monkeypatch.setenv("OBJDEL_RUN_DIR", str(tmp_path / "runs"))
    monkeypatch.setenv("OBJDEL_DATA_ROOT", str(data))
    root = data / "views" / "lidc-train-b1"
    (root / "images").mkdir(parents=True)
    (root / "images" / "a.png").write_bytes(b"png")
    (root / "manifest.json").write_bytes(MANIFEST)
    (data / "views/lidc-bundle-v3").mkdir()
    (data / "views/lidc-bundle-v3/manifest.json").write_text(json.dumps({"view_hash": "h"}))
    (data / "eval_views/lidc-val-edge").mkdir()
    return types.SimpleNamespace(cfg=b1_rows_copy(tmp_path, hashlib.sha256(MANIFEST).hexdigest()), root=root,
                                 data=data)


def test_b1_rows_reads_the_bundle_view_and_one_rows_root(rows_cfg):
    s = PF.inputs(rows_cfg.cfg)
    assert s.views == (rows_cfg.data / "views/lidc-bundle-v3",)
    assert [(st.name, PF.local(st.train)) for st in s.rows] == [("train", rows_cfg.root)]
    assert PF.given_inputs(s) == []
    assert PF.spec_inputs(s.spec) == ["{data}/rows/lidc-train-b1", "{data}/eval_views/lidc-val-shifted",
                                      "{data}/eval_views/lidc-val-edge"]


def test_specs_b1_rows_passes_on_fixture_data_and_refuses_an_edit(rows_cfg, monkeypatch):
    from test_data_quality import clean_report
    gates_report(monkeypatch, clean_report())
    data_checks = [c for c in PF.CHECKS if c[0] in ("DATA spec inputs", "DATA ROWS root")]
    code, stamp = run(rows_cfg.cfg, monkeypatch, data_checks)
    assert code == 0 and PF.refusal(rows_cfg.cfg) is None
    rows_cfg.cfg.write_text(rows_cfg.cfg.read_text().replace("steps: 260", "steps: 261"))
    assert "stale (config" in PF.refusal(rows_cfg.cfg)


def test_a_missing_data_input_fails_spec_inputs(rows_cfg):
    (rows_cfg.data / "eval_views/lidc-val-edge").rmdir()
    assert PF.given_inputs(PF.inputs(rows_cfg.cfg)) == [
        f"{{data}}/eval_views/lidc-val-edge: no folder {rows_cfg.data / 'eval_views/lidc-val-edge'}"]


def test_a_cmd_and_eval_spec_has_no_training_view_and_no_rows_root(cfg, tmp_path):
    ev = {"kind": "eval", "name": "eval", "ckpt": "{rd}/ckpt",
          "view": ["{data}/eval_views/lidc-val-shifted"], "out": ["{rd}/eval"]}
    spec = tiny_spec(tmp_path, stages=(CMD, ev), bundle_view=None)
    s = PF.inputs(spec)
    assert s.views == () and s.rows == () and PF.given_inputs(s) == []
    assert [f(s) for f in (PF.patients, PF.ct_window, PF.targets, PF.coverage, PF.rows_root)] == [[]] * 5


def test_a_planted_shared_patient_fails_disjointness(cfg, tmp_path):
    data = tmp_path / "data"
    rows = [{"item_id": "lidc/s1/a", "group_id": "P1"}]
    for split, recs in (("train", rows), ("val", [{"item_id": "lidc/s2/b", "group_id": "P2"}])):
        (data / "prepared/lidc").mkdir(parents=True, exist_ok=True)
        (data / f"prepared/lidc/{split}.jsonl").write_text("".join(json.dumps(r) + "\n" for r in recs))
    (data / "views/v/rows.jsonl").write_text(json.dumps({"item_id": "lidc/s1/a"}) + "\n")
    for v in PF.EXAMS:
        (data / f"eval_views/{v}/lidc").mkdir()
        row = {"item_id": "lidc/s3/c", "group_id": "P3"}
        (data / f"eval_views/{v}/lidc/val.jsonl").write_text(json.dumps(row) + "\n")
    assert PF.patients(PF.inputs(cfg)) == []
    (data / f"eval_views/{PF.EXAMS[0]}/lidc/val.jsonl").write_text(
        json.dumps({"item_id": "lidc/s9/c", "group_id": "P1"}) + "\n")
    assert PF.patients(PF.inputs(cfg)) == [f"{PF.EXAMS[0]}: 1 patients also in training, e.g. ['P1']"]


# ------------------------------------------------------------------------------- ROWS root and revision (B1.3, I19)
def gates_report(monkeypatch, report):
    from objdel.data import quality as Q
    from objdel.policy import models
    seen = []
    monkeypatch.setattr(models, "tokenizer", lambda model: (None, 1536))
    monkeypatch.setattr(Q, "audit", lambda root, *a, **k: seen.append(root) or report)
    return seen


def test_rows_root_check_reads_the_stage_and_refuses_a_hard_gate_failure(rows_cfg, monkeypatch):
    from test_data_quality import clean_report
    s = PF.inputs(rows_cfg.cfg)
    seen = gates_report(monkeypatch, clean_report())
    assert PF.rows_root(s) == [] and seen == [rows_cfg.root]
    report = clean_report()
    report["windows"]["clear"] -= 1  # a crowded policy frame
    report["agreement"]["share_lt_0.3"] = 1.0  # DG14b, a report: never refuses
    gates_report(monkeypatch, report)
    assert [p.split(":")[0] for p in PF.rows_root(s)] == ["DG5"]


def test_rows_root_check_refuses_a_missing_or_changed_manifest(rows_cfg, monkeypatch, tmp_path):
    from test_data_quality import clean_report
    s = PF.inputs(rows_cfg.cfg)
    gates_report(monkeypatch, clean_report())
    (rows_cfg.root / "manifest.json").write_bytes(MANIFEST + b" ")
    assert "SHA-256" in PF.rows_root(s)[0]
    (rows_cfg.root / "manifest.json").unlink()
    assert "no manifest" in PF.rows_root(s)[0]
    (rows_cfg.root / "manifest.json").write_bytes(MANIFEST)
    assert PF.rows_root(s) == []
    assert "is not stage train's 000" in PF.rows_root(PF.inputs(b1_rows_copy(tmp_path, "0" * 64)))[0]


def test_an_edited_rows_manifest_changes_the_data_hash(rows_cfg):
    before = PF.hashes(rows_cfg.cfg)["data"]
    (rows_cfg.root / "manifest.json").write_bytes(MANIFEST.replace(b"{}", b"{ }"))
    assert PF.hashes(rows_cfg.cfg)["data"] != before


class FakeHub:
    def __init__(self, files, manifest):
        self.files, self.manifest, self.calls = files, manifest, []

    def api(self, token=None):
        return types.SimpleNamespace(list_repo_files=lambda repo, revision: self.calls.append(revision) or self.files)

    def download(self, tmp):
        def get(repo, filename, revision, token=None):
            self.calls.append((filename, revision))
            p = tmp / "hub_manifest.json"
            p.write_bytes(self.manifest)
            return str(p)
        return get


@pytest.mark.parametrize("files,manifest,ok", [
    (["rows/lidc-train-b1/images/a.png", "rows/lidc-train-b1/manifest.json", "README.md"], MANIFEST, True),
    (["rows/lidc-train-b1/manifest.json"], MANIFEST, False),  # an interrupted upload: an image missing
    (["rows/lidc-train-b1/images/a.png", "rows/lidc-train-b1/manifest.json"], b"{}", False)])
def test_rows_revision_compares_the_pinned_hf_revision(rows_cfg, monkeypatch, tmp_path, files, manifest, ok):
    import huggingface_hub
    hub = FakeHub(files, manifest)
    monkeypatch.setattr(huggingface_hub, "HfApi", hub.api)
    monkeypatch.setattr(huggingface_hub, "hf_hub_download", hub.download(tmp_path))
    assert (PF.rows_revision(PF.inputs(rows_cfg.cfg)) == []) is ok
    assert hub.calls[0] == "34213dceb068128f1a1abf9c67aa65dd59bc1220"  # the revision specs/b1_rows.yaml pins


def test_the_trl_contract_is_the_first_code_check(monkeypatch):
    from objdel.rl import trl_adapter as T

    name, fn = next(c for c in PF.CHECKS if c[0].startswith("CODE"))
    assert name == "CODE TRL private names" and fn(None) == []

    def broken():
        raise T.TRLContractError("boom: _compute_loss")

    monkeypatch.setattr(T, "check", broken)
    assert fn(None) == ["boom: _compute_loss"]

    def missing():
        raise ImportError("no trl")

    monkeypatch.setattr(T, "check", missing)
    with pytest.raises(ImportError):  # any other error propagates: main() records it as a crashed (failed) check
        fn(None)


@pytest.mark.parametrize("out, code, ok", [("1.11.0\n", 0, True), ("vastai 1.11.0\n", 0, True), ("1.10.0", 0, False),
                                           ("1.11.1", 0, False), ("1.11.01", 0, False), ("11.11.0", 0, False),
                                           ("1.11", 0, False), ("", 0, False), ("1.11.0", 3, False)])
def test_vastai_version_check_wants_the_pinned_token(tmp_path, monkeypatch, out, code, ok):
    """CT-46: the check passes on exactly 1.11.0 and fails with a message naming both versions (no traceback)."""
    bin_ = tmp_path / "bin"
    bin_.mkdir()
    fake = bin_ / "vastai"
    fake.write_text(f"#!/bin/sh\nprintf '%s' '{out}'\nexit {code}\n")
    fake.chmod(0o755)
    monkeypatch.setenv("PATH", str(bin_))
    bad = PF.vastai_problems()
    assert PF.VASTAI_VERSION == "1.11.0"
    assert (bad == []) is ok
    assert ok or ("1.11.0" in bad[0] and "vastai" in bad[0])


def test_vastai_version_check_fails_without_the_binary(tmp_path, monkeypatch):
    monkeypatch.setenv("PATH", str(tmp_path))
    bad = PF.vastai_problems()
    assert len(bad) == 1 and "1.11.0" in bad[0] and "uv tool install vastai==1.11.0" in bad[0]
    assert any(n == "CODE vastai version" for n, _ in PF.CHECKS)


# --- C0: the three V2 specs (U-CF4) --------------------------------------------------------------------------------

C_PLANS = ("c_speed.yaml", "c_speed_b.yaml", "c_speed_c6.yaml")
SPECS = ROOT / "specs"
BOX_ROWS = "/workspace/data/lidc/rows/lidc-train-b1"
G_GR1, G_RC0 = "tests/test_c_gpu.py::test_g_gr1_requests_independent", "tests/test_c_gpu.py::test_g_rc0_per_row_caps"


def plan_steps(name):
    from test_b1_plans import steps

    return steps(name)


def train_rl_steps(name):
    return [(s, PF.train_rl_argv(s)) for s in plan_steps(name) if PF.train_rl_argv(s) is not None]


def b1_pins():
    """(rows_revision, rows_manifest_sha256) of specs/b1_rows.yaml's rl stage: the pins every C arm must carry."""
    from objdel.core.spec import load_spec

    st = next(s for s in load_spec(SPECS / "b1_rows.yaml").stages if s.kind == "rl" and not s.merge_last)
    return st.rows_revision, st.rows_manifest_sha256


def check_train_rl(argv, ref):
    """U-CF4 on one `train_rl` command: it parses with train_rl's own parser, trains the B1 rows root with the same
    rows pins as specs/b1_rows.yaml, in whole 20-step blocks, and never names --train-views."""
    from experiments.phase2 import train_rl

    a = train_rl.parse(argv)  # SystemExit(2) on an unknown or missing flag
    assert "--train-views" not in argv
    assert str(a.train) == BOX_ROWS, a.train
    assert (a.rows_revision, a.rows_manifest_sha256) == ref and a.steps % 20 == 0, (a.rows_revision, a.steps)
    return a


def test_train_rl_argv_finds_a_direct_step_and_one_inside_a_bash_string():
    direct = ["~", "-m", "experiments.phase2.train_rl", "--model", "m"]
    wrapped = ["bash", "-c", "NSYS_TIMEOUT=1800 bash nsys_capture.sh /rd c 16 17 -- /opt/venv/bin/python -I -m "
                             "experiments.phase2.train_rl --model m --grpo '{\"a\": false}'"]
    assert PF.train_rl_argv(direct) == ["--model", "m"]
    assert PF.train_rl_argv(wrapped) == ["--model", "m", "--grpo", '{"a": false}']
    assert PF.train_rl_argv(["-m", "experiments.phase2.final_eval", "--ckpt", "c"]) is None
    assert PF.train_rl_argv(["bash", "-c", "echo experiments.phase2.train_rl"]) is None  # named, not run with -m
    assert PF.train_rl_argv(["bash", "env_check.sh", "/rd/ckpt"]) is None


def test_c_speed_is_the_launch_config_of_v2():
    """E6.6: c_launch.yaml (view, rows_plan) is archived; c_speed.yaml carries the B1 bundle view and an RTX 5090."""
    from objdel.core.spec import load_spec

    spec, b1 = load_spec(SPECS / "c_speed.yaml"), load_spec(SPECS / "b1_rows.yaml")
    assert spec.bundle_view == b1.bundle_view and spec.compute.gpu == b1.compute.gpu
    assert not [s for s in spec.stages if s.kind in ("sft", "rl")]  # no model kept: nothing to push


def test_the_rows_check_reads_the_pins_of_the_nsys_wrapped_arms(monkeypatch, tmp_path):
    """c_speed.yaml has no rl stage: preflight's rows checks read the pins inside the arms' bash strings, as the rl
    stage each arm's words parse to (rl_stages)."""
    monkeypatch.setenv("OBJDEL_DATA_ROOT", str(tmp_path))
    s = PF.inputs(SPECS / "c_speed.yaml")
    assert [st.name for st in s.rows] == ["c0_base", "c_fast"]
    assert {(st.rows_revision, st.rows_manifest_sha256) for st in s.rows} == {b1_pins()}
    assert {PF.local(st.train) for st in s.rows} == {tmp_path / "views/lidc-train-b1"}
    fast = s.rows[1]
    assert (fast.checks_from_step, fast.expect_step1, fast.record_batches) == (
        18, "{rd}/c0_base/rl.json", "{rd}/c_check/arm_batches")
    assert [st.name for st in PF.inputs(SPECS / "c_speed_c6.yaml").rows] == ["c6"]


def test_an_arm_whose_train_rl_words_do_not_parse_fails_the_inputs(tmp_path):
    bad = {"kind": "cmd", "name": "arm", "argv": ["bash", "-c", "bash nsys_capture.sh {rd} a 16 17 -- python -I -m "
                                                  "experiments.phase2.train_rl --model m --no-such-flag"]}
    with pytest.raises(ValueError, match="stage arm"):
        PF.inputs(tiny_spec(tmp_path, stages=[bad]))


@pytest.mark.parametrize("plan", C_PLANS)
def test_u_cf4_every_module_step_parses_with_its_module(plan):
    pytest.importorskip("trl")
    from test_b1_plans import module_steps, parse_step

    found = module_steps(plan)
    assert found or plan == "c_speed_b.yaml"  # the mode B plan is one nsys-wrapped train_rl step: the next test
    for s in found:
        parse_step(s)


@pytest.mark.parametrize("plan", C_PLANS)
def test_u_cf4_every_train_rl_step_trains_the_b1_rows_in_whole_blocks(plan):
    pytest.importorskip("trl")
    ref = b1_pins()
    found = train_rl_steps(plan)
    assert found
    for _, argv in found:
        check_train_rl(argv, ref)


def test_u_cf4_a_planted_train_rl_step_fails_each_rule():
    pytest.importorskip("trl")
    ref = b1_pins()
    good = train_rl_steps("c_speed.yaml")[0][1]
    check_train_rl(good, ref)
    with pytest.raises(SystemExit) as e:
        check_train_rl([*good, "--train-views", "x"], ref)
    assert e.value.code == 2
    bad_root = [BOX_ROWS.replace("rows/lidc-train-b1", "prepared") if w == BOX_ROWS else w for w in good]
    with pytest.raises(AssertionError):
        check_train_rl(bad_root, ref)
    with pytest.raises(AssertionError):
        check_train_rl([("17" if good[i - 1] == "--steps" else w) for i, w in enumerate(good)], ref)
    with pytest.raises(AssertionError):
        check_train_rl(good, ("0" * 40, ref[1]))


def flag(argv, name):
    return argv[argv.index(name) + 1]


def test_u_cf4_the_arms_carry_the_check_flags_of_c0_step_3():
    pytest.importorskip("trl")
    rd = "/workspace/runs/r"
    (base_step, base), (fast_step, fast) = train_rl_steps("c_speed.yaml")
    (b_step, mode_b), = train_rl_steps("c_speed_b.yaml")
    assert flag(base, "--out") == f"{rd}/c0_base" and flag(base, "--steps") == "20" and "--expect-step1" not in base
    assert "bash experiments/phase2/docker/nsys_capture.sh " + f"{rd} c0_base 16 17 -- " in base_step[2]
    assert "NSYS_TIMEOUT=1800" in fast_step[2] and f"{rd} c_fast 16 17 -- " in fast_step[2]
    for arm, folder in ((fast, "arm_batches"), (mode_b, "arm_batches_b")):
        assert flag(arm, "--expect-step1") == f"{rd}/c0_base/rl.json"
        assert flag(arm, "--encode-check") == f"{rd}/c_check/encoder.json"
        assert flag(arm, "--record-batches") == f"{rd}/c_check/{folder}"
        assert (flag(arm, "--steps"), flag(arm, "--checks-from-step")) == ("40", "18")
        assert flag(arm, "--check-rope-every") == flag(arm, "--check-meta-every") == "1"
        assert {"--fast-sync", "--host-meta", "--prompt-memo", "--shared-prefix"} <= set(arm)
        assert flag(arm, "--gpu-pixels") == "auto"
    assert json.loads(flag(fast, "--grpo")) == {"vllm_enable_sleep_mode": False}  # mode A: vLLM awake
    assert json.loads(flag(mode_b, "--grpo")) == {"vllm_enable_sleep_mode": True}  # mode B: sleep
    assert f"{rd} c_fast_b 16 17 -- " in b_step[2] and flag(mode_b, "--out") == f"{rd}/c_fast_b"


def test_u_cf4_the_c6_plan_runs_g_gr1_first_then_the_group_requests_arm():
    pytest.importorskip("trl")
    first, arm = plan_steps("c_speed_c6.yaml")[:2]
    assert first[-1] == G_GR1 and first[:5] == ["-m", "pytest", "-q", "-p", "no:cacheprovider"]
    assert "--group-requests" in PF.train_rl_argv(arm) and flag(PF.train_rl_argv(arm), "--steps") == "20"
    assert not [a for p in C_PLANS for _, a in train_rl_steps(p) if "--no-ckpt-layers" in a]  # no C9 arm


def test_u_cf4_c_speed_starts_with_the_cuda_build_and_checks_the_box_first():
    steps = plan_steps("c_speed.yaml")
    assert steps[0] == plan_steps("b1_rows.yaml")[0]  # step 0 is b1_rows.yaml's CUDA build (front 6.1 step 4d)
    check = next(s for s in steps if "experiments.phase2.c_speed_check" in s)
    assert {"--pixels", "--sync"} <= set(check) and flag(check, "--out") == "/workspace/runs/r"
    assert steps.index(check) < min(steps.index(s) for s, _ in train_rl_steps("c_speed.yaml"))  # item 2 before 3
    assert not [s for s in plan_steps("c_speed_b.yaml") + plan_steps("c_speed_c6.yaml")[1:]
                if "run_box.sh" in " ".join(s)]  # no build step: they run after c_speed.yaml


@pytest.mark.parametrize("node", [G_GR1, G_RC0])
def test_u_cf4_the_two_vllm_box_tests_collect_exactly_one_item(node):
    import subprocess
    import sys

    r = subprocess.run([sys.executable, "-m", "pytest", "--collect-only", "-q", "-p", "no:cacheprovider", node],
                       cwd=ROOT, capture_output=True, text=True, timeout=300)
    assert r.returncode == 0, r.stdout[-1500:]
    assert [ln for ln in r.stdout.splitlines() if "::" in ln] == [node]

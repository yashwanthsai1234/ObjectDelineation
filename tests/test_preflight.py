"""The preflight driver: stamp, staleness (refusal) and fail-fast by stage. The checks themselves run on data/."""

import hashlib
import json
import types
from pathlib import Path

import pytest
import yaml

from experiments.phase2 import preflight as PF


@pytest.fixture
def cfg(tmp_path, monkeypatch):
    monkeypatch.setenv("OBJDEL_RUN_DIR", str(tmp_path / "runs"))
    data = tmp_path / "data"
    for d in ("views/v", "checker_eval", *(f"eval_views/{v}" for v in PF.EXAMS)):
        (data / d).mkdir(parents=True)
        (data / d / "manifest.json").write_text(json.dumps({"view_hash": "h"}))
    c = tmp_path / "c.yaml"
    c.write_text(f"model: lfm2.5-vl-450m\nview: views/v\noutput_dir: x\ndata_root: {data}\n")
    return c


def run(cfg, monkeypatch, checks):
    monkeypatch.setattr(PF, "CHECKS", checks)
    code = PF.main(["--config", str(cfg)])
    return code, json.loads(PF.stamp_path(cfg).read_text())


def test_a_passing_run_writes_a_stamp_the_gate_accepts_until_something_changes(cfg, monkeypatch):
    code, stamp = run(cfg, monkeypatch, [("DATA a", lambda a, c: []), ("CODE b", lambda a, c: [])])
    assert code == 0 and stamp["pass"] and set(stamp["hashes"]) == {"code", "data", "config"}
    assert PF.refusal(cfg) is None
    cfg.write_text(cfg.read_text() + "seed: 2\n")
    assert "stale (config" in PF.refusal(cfg)
    (cfg.parent / "data/eval_views/probe-topology/manifest.json").write_text("{}")
    assert "stale (config, data" in PF.refusal(cfg)


def test_a_changed_checker_eval_set_makes_the_stamp_stale(cfg, monkeypatch):
    run(cfg, monkeypatch, [("DATA a", lambda a, c: [])])
    assert PF.refusal(cfg) is None
    (cfg.parent / "data/checker_eval/manifest.json").write_text("{}")
    assert "stale (data" in PF.refusal(cfg)


def test_a_failure_finishes_its_stage_and_skips_later_ones(cfg, monkeypatch):
    ran = []
    checks = [("DATA a", lambda a, c: ["bad a"]), ("DATA b", lambda a, c: ran.append("b") or []),
              ("CHECKER c", lambda a, c: ran.append("c") or [])]  # fmt: skip
    code, stamp = run(cfg, monkeypatch, checks)
    assert code == 1 and not stamp["pass"] and ran == ["b"] and list(stamp["checks"]) == ["DATA a", "DATA b"]
    assert "FAILED" in PF.refusal(cfg)


def test_no_stamp_is_a_refusal(cfg):
    assert "no preflight stamp" in PF.refusal(cfg)


def test_a_crashing_check_fails_the_preflight(cfg, monkeypatch):
    def boom(a, c):
        raise KeyError("qwen3-vl-2b")
    code, stamp = run(cfg, monkeypatch, [("CODE boom", boom)])
    assert code == 1 and not stamp["pass"] and "check crashed: KeyError" in stamp["checks"]["CODE boom"]["problems"][0]


# ------------------------------------------------------------------------------- ROWS root and revision (B1.3, I19)
MANIFEST = b'{"schedule": {}}\n'


@pytest.fixture
def rows_cfg(cfg, tmp_path):
    """The launch config with rows_plan: a plan whose train_rl stage names {data}/rows/lidc-train-b1, built locally at
    <data_root>/views/lidc-train-b1 (one image and the manifest)."""
    root = tmp_path / "data" / "views" / "lidc-train-b1"
    (root / "images").mkdir(parents=True)
    (root / "images" / "a.png").write_bytes(b"png")
    (root / "manifest.json").write_bytes(MANIFEST)
    plan = tmp_path / "b1_rows.yaml"
    step = ["-m", "experiments.phase2.train_rl", "--train", "{data}/rows/lidc-train-b1", "--rows-revision", "a" * 40,
            "--rows-manifest-sha256", hashlib.sha256(MANIFEST).hexdigest()]
    plan.write_text(yaml.safe_dump({"steps": [step]}))
    cfg.write_text(cfg.read_text() + f"rows_plan: {plan}\n")
    return types.SimpleNamespace(cfg=cfg, plan=plan, root=root)


def test_a_plan_file_never_loads_as_a_run_config(rows_cfg):
    with pytest.raises(TypeError, match="steps"):
        PF.load(rows_cfg.plan)
    a, _ = PF.load(rows_cfg.cfg)
    assert a.rows_plan == str(rows_cfg.plan)


def test_local_rows_root_maps_only_the_literal_box_path():
    assert PF.local_rows_root("{data}/rows/lidc-train-b1") == Path("data/views/lidc-train-b1")
    assert PF.local_rows_root("{data}/rows/x", "/d") == Path("/d/views/x")
    for bad in ("{data}/rows/../x", "{data}/views/lidc-train-b1", "/workspace/data/rows/x", "{data}/rows/a/b",
                "{data}/rows/.."):
        with pytest.raises(ValueError):
            PF.local_rows_root(bad)


def gates_report(monkeypatch, report):
    from experiments.phase2 import audit_data
    seen = []
    monkeypatch.setattr(audit_data, "audit", lambda root, *a, **k: seen.append(root) or report)
    return seen


def test_rows_root_check_reads_the_plan_and_refuses_a_hard_gate_failure(rows_cfg, monkeypatch):
    from test_data_quality import clean_report
    a, _ = PF.load(rows_cfg.cfg)
    seen = gates_report(monkeypatch, clean_report())
    assert PF.rows_root(a, rows_cfg.cfg) == [] and seen == [rows_cfg.root]
    report = clean_report()
    report["windows"]["clear"] -= 1  # a crowded policy frame
    report["agreement"]["share_lt_0.3"] = 1.0  # DG14b, a report: never refuses
    gates_report(monkeypatch, report)
    assert [p.split(":")[0] for p in PF.rows_root(a, rows_cfg.cfg)] == ["DG5"]


def test_rows_root_check_refuses_a_missing_or_changed_manifest(rows_cfg, monkeypatch):
    from test_data_quality import clean_report
    a, _ = PF.load(rows_cfg.cfg)
    gates_report(monkeypatch, clean_report())
    (rows_cfg.root / "manifest.json").write_bytes(MANIFEST + b" ")
    assert "SHA-256" in PF.rows_root(a, rows_cfg.cfg)[0]
    (rows_cfg.root / "manifest.json").unlink()
    assert "no manifest" in PF.rows_root(a, rows_cfg.cfg)[0]


def test_no_rows_plan_no_rows_checks(cfg):
    a, _ = PF.load(cfg)
    assert PF.rows_root(a, cfg) == [] and PF.rows_revision(a, cfg) == []


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
    a, _ = PF.load(rows_cfg.cfg)
    assert (PF.rows_revision(a, rows_cfg.cfg) == []) is ok
    assert hub.calls[0] == "a" * 40


def test_the_trl_contract_is_the_first_code_check(monkeypatch):
    from objdel.train import trl_adapter as T

    name, fn = next(c for c in PF.CHECKS if c[0].startswith("CODE"))
    assert name == "CODE TRL private names" and fn(None, None) == []

    def broken():
        raise T.TRLContractError("boom: _compute_loss")

    monkeypatch.setattr(T, "check", broken)
    assert fn(None, None) == ["boom: _compute_loss"]

    def missing():
        raise ImportError("no trl")

    monkeypatch.setattr(T, "check", missing)
    with pytest.raises(ImportError):  # any other error propagates: main() records it as a crashed (failed) check
        fn(None, None)

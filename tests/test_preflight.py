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

from experiments.phase2 import preflight as PF

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
    code = PF.main(["--config", str(cfg)])
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
    assert PF.spec_inputs(s.spec) == ["{data}/rows/lidc-train-b1", "{data}/prepared",
                                      "{data}/eval_views/lidc-val-shifted", "{data}/eval_views/lidc-val-edge"]


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
    ev = {"kind": "eval", "name": "eval", "ckpt": "{rd}/ckpt", "train": "{data}/prepared",
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
    from experiments.phase2 import audit_data
    seen = []
    monkeypatch.setattr(audit_data, "audit", lambda root, *a, **k: seen.append(root) or report)
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
    from objdel.train import trl_adapter as T

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

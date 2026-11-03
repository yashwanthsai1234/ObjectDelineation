"""The preflight driver: stamp, staleness (refusal) and fail-fast by stage. The checks themselves run on data/."""

import json

import pytest

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

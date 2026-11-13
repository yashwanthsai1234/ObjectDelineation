"""Break-it BK2 (Phase 2 Stage 1, wave 2): the scoring and decision path, JSONL bridge -> serve -> metrics / summarize
-> probes -> signal -> prereg -> gates.phase2_entry -> select_ckpt / eval_run. Findings are recorded in
docs/archive/experiments/phase2/impl/BK2.md;
a real defect is an xfail(strict=True) "BK2-<n>", every other test is a guard that passes today.
Real-data tests are marked `data` (data/eval_views)."""

import json
import shutil
from pathlib import Path

import numpy as np
import probe_models as M
import pytest
from oracles.eval_oracle import pairwise_iou

from objdel.cli.eval import baselines
from objdel.cli.train import prereg as P
from objdel.core.hashing import file_sha256
from objdel.data import probes
from objdel.policy import evaluate
from objdel.policy import sampler as V
from objdel.sft import select as sft_select
from objdel.verify import gates
from objdel.verify import stats as S
from objdel.verify import store as serve
from objdel.verify.metrics import _ged, metrics, min_mode_retention, summarize, tv_js
from objdel.verify.reward.cache import ExpertCache
from tests.support.shapes import disc as disc_hw
from tests.support.targets import rater_texts

ROOT = Path(__file__).resolve().parents[2]
VIEWS = ROOT / "data" / "eval_views"
SIDE = 50
STOP = {1}


def disk(cy, cx, r, side=SIDE):
    return disc_hw(side, side, cy, cx, r)


def item():
    raters = {"nodule": np.stack([disk(25, 25, 10), disk(25, 26, 11), disk(24, 25, 9)])}
    return raters, rater_texts(raters)


def m_of(samples, truncated=None, tokens=None):
    raters, _ = item()
    return metrics(samples, raters, SIDE, SIDE, ("nodule",), truncated, tokens=tokens)


# ================================================================ JSONL bridge (sample.import_completions)


def _write(path, recs):
    path.write_text("".join(json.dumps(r, ensure_ascii=False) + "\n" for r in recs))
    return path


def _rec(i, k, text="x", ids=(5, 1), finish="stop"):
    return {"id": i, "k": k, "text": text, "token_ids": list(ids), "finish_reason": finish}


@pytest.fixture
def sample_mod():
    pytest.importorskip("trl")
    from objdel.policy import sampler

    return sampler


def test_bridge_unicode_empty_and_length_round_trip(sample_mod, tmp_path):
    recs = [_rec("lidc/a", 0, "ü ∅ — 💥", (5, 1)), _rec("lidc/a", 1, "", ()),
            _rec("lidc/a", 2, "abc", (5, 6), "length")]
    rows = sample_mod.import_completions(_write(tmp_path / "c.jsonl", recs), STOP, n=3)
    assert rows == [{"item_id": "lidc/a", "samples": ["ü ∅ — 💥", "", "abc"], "truncated": [False, True, True],
                     "tokens": [2, 0, 2]}]


def test_bridge_gap_in_k_is_refused(sample_mod, tmp_path):
    with pytest.raises(ValueError):
        sample_mod.import_completions(_write(tmp_path / "c.jsonl", [_rec("a", 0), _rec("a", 2)]), STOP)


def test_bridge_duplicate_id_k_is_refused(sample_mod, tmp_path):
    recs = [_rec("a", 0, "first"), _rec("a", 1), _rec("a", 0, "second")]
    with pytest.raises(ValueError):
        sample_mod.import_completions(_write(tmp_path / "c.jsonl", recs), STOP)


def test_bridge_short_item_is_refused(sample_mod, tmp_path):
    recs = [_rec("a", k) for k in range(8)] + [_rec("b", k) for k in range(7)]  # file cut at a line boundary
    with pytest.raises(ValueError):
        sample_mod.import_completions(_write(tmp_path / "c.jsonl", recs), STOP)


def test_bridge_truncated_last_line_is_loud(sample_mod, tmp_path):
    p = tmp_path / "c.jsonl"
    p.write_text(json.dumps(_rec("a", 0)) + "\n" + json.dumps(_rec("a", 1))[:20])
    with pytest.raises(json.JSONDecodeError):
        sample_mod.import_completions(p, STOP)


def test_sampler_outputs_need_one_result_per_prompt():
    class Out:
        def __init__(self, n):
            self.outputs = [type("C", (), {"text": "t", "token_ids": [1], "finish_reason": "stop", "logprobs": None})()
                            for _ in range(n)]
    lines = [{"id": "a"}, {"id": "b"}]
    with pytest.raises(ValueError):
        list(V.outputs(lines, [Out(8)]))  # vLLM dropped a request: strict zip
    assert [r["k"] for r in V.outputs(lines[:1], [Out(3)])] == [0, 1, 2]


# ================================================================ serve: preload, CacheRef, fork pool


def _mini_root(out: Path, n=3) -> Path:
    raters, _ = item()
    rows = []
    for i in range(n):
        key = f"lidc/1.2.{i}/n000/z000"
        (out / "caches" / key).mkdir(parents=True)
        d = ExpertCache.build(np.roll(raters["nodule"], i, axis=2)).save(out / "caches" / key / "nodule.npz")
        rows.append({"item_id": key, "rl": True, "caches": {"nodule": d}, "W": SIDE, "H": SIDE})
    p = out / "lidc" / "val.jsonl"
    p.parent.mkdir(parents=True)
    p.write_text("".join(json.dumps(r) + "\n" for r in rows))
    (out / "manifest.json").write_text(json.dumps({"files": {"lidc/val.jsonl": {
        "sha256": file_sha256(p)}}}))
    return out


pytestmark = pytest.mark.usefixtures("clean_store")  # conftest: an empty serve.STORE per test


def test_cacheref_for_a_key_not_preloaded(tmp_path):
    root = _mini_root(tmp_path / "r")
    with pytest.raises(KeyError):
        serve.loader(root)("lidc/1.2.0/n000/z000")
    serve.preload(root, ("lidc",))
    with pytest.raises(KeyError):
        serve.loader(root)("lidc/nope")
    with pytest.raises(KeyError):
        serve.CacheRef("lidc/1.2.0/n000/z000", str(tmp_path / "other"))
    shutil.copytree(root, tmp_path / "r2")
    serve.preload(tmp_path / "r2", ("lidc",))
    assert len(serve.loader(tmp_path / "r2")("lidc/1.2.0/n000/z000")) == 1


def test_preload_refuses_tampered_jsonl_stray_jsonl_and_cache(tmp_path):
    root = _mini_root(tmp_path / "r")
    p = root / "lidc" / "val.jsonl"
    good = p.read_bytes()
    p.write_bytes(good.replace(b'"W": 50', b'"W": 51', 1))
    with pytest.raises(ValueError, match="SHA-256"):
        serve.preload(root, ("lidc",))
    p.write_bytes(good)
    (root / "lidc" / "extra.jsonl").write_text("{}\n")
    with pytest.raises(ValueError, match="SHA-256"):
        serve.preload(root, ("lidc",))
    (root / "lidc" / "extra.jsonl").unlink()
    npz = next((root / "caches").rglob("*.npz"))
    shutil.copy(next(x for x in (root / "caches").rglob("*.npz") if x != npz), npz)  # another item's cache
    with pytest.raises(ValueError):
        serve.preload(root, ("lidc",))


# ================================================================ metrics / summarize


def test_all_failed_group_is_finite_and_json_strict():
    m = m_of([""] * 8, [True] * 8, tokens=[0] * 8)
    assert m["samples"]["legal"] == [0.0] * 8 and m["groups"]["signal"] == [False]
    assert m["objects"]["nodule"]["hm_iou"] == 0.0 and m["pass_at_k"]["8"] == 0.0
    s = summarize([{"item_id": "lidc/a", "dataset": "lidc", **m}])
    json.dumps(s, allow_nan=False)
    assert s["datasets"]["lidc"]["iou_tc_pass"] == 0.0 and s["datasets"]["lidc"]["tokens"]["mean"] == 0.0


def test_identical_group_has_no_signal_and_full_duplicates():
    _, texts = item()
    m = m_of([texts[0]] * 8)
    assert m["groups"]["dup"] == [7 / 8] and m["groups"]["signal"] == [False]
    assert max(abs(a) for a in m["samples"]["advantage"]) == 0.0


def test_k_not_8_is_one_group_and_gates_refuse_it():
    _, texts = item()
    m = m_of((texts * 4)[:12])
    assert m["group"] == 12 and set(m["pass_at_k"]) == {"1", "2", "4", "8"}
    s = summarize([{"item_id": "lidc/a", "dataset": "lidc", **m}])
    rep = gates.phase2_entry({"views": {"shifted": {"summary": s}}}, ["lidc"])
    assert rep["gates"]["G-signal"]["values"]["sigma_share/shifted/lidc"] is None
    assert rep["passed"]["G-signal"] is False


def test_select_ckpt_eligibility_needs_groups_of_8():
    summary = {"format": 1.0, "legal": 1.0, "dup_rate": 0.0, "signal_share": 1.0, "group": 12}
    assert gates.eligible({"datasets": {"lidc": summary}}, {"pass": True}) is False  # select_ckpt's D12 call
    assert gates.eligible({"datasets": {"lidc": summary | {"group": 8}}}, {"pass": True}) is True


def test_select_ckpt_eligibility_needs_a_vision_result():
    summary = {"format": 1.0, "legal": 1.0, "dup_rate": 0.0, "signal_share": 1.0, "group": 8}
    for vision in ({}, sft_select.vision_check(Path("b"), None)):  # unmeasured; --no-vision
        assert gates.eligible({"datasets": {"lidc": summary}}, vision) is False


def test_summarize_mixed_tokens_and_mixed_k():
    _, texts = item()
    a = {"item_id": "lidc/a", "dataset": "lidc", **m_of(texts * 3 + texts[:1] * 7, tokens=[3] * 16)}
    b = {"item_id": "lidc/b", "dataset": "lidc", **m_of(texts * 2 + texts[:2])}
    s = summarize([a, b])["datasets"]["lidc"]
    assert s["tokens"] is None and s["group"] == 8 and set(s["pass_at_k"]) == {"1", "2", "4", "8"}


def test_ged_is_never_negative():
    m = np.random.default_rng(1).random((4, 12, 12)) > 0.5
    s = np.concatenate([m, m])  # 8 samples = each of the 4 raters twice: GED is exactly 0
    assert _ged(pairwise_iou(np.concatenate([s, m])), 8) >= 0.0  # actually -1.1e-16


def test_divergences_and_retention_stay_in_range():
    rng = np.random.default_rng(1)
    for _ in range(3000):
        k = int(rng.integers(2, 7))
        a = int(rng.integers(1, k + 1))
        f = {0: a / k} | ({1: (k - a) / k} if a < k else {})
        c = rng.multinomial(8, [0.4, 0.4, 0.2])
        p = {i: c[i] / 8 for i in range(3) if c[i]}
        tv, js = tv_js(p, f)
        assert 0.0 <= tv <= 1.0 and 0.0 <= js <= 1.0
        assert 0.0 <= min_mode_retention(p, f) <= 1.0


# ================================================================ probes


def _probe_rows(root, n=None, model=M.model_exact):
    recs = probes.records(root)[:n]
    rng = np.random.default_rng(0)
    return [{"item_id": r["item_id"], "samples": [model(r, rng) for _ in range(8)]} for r in recs]


@pytest.mark.data
def test_topology_probe_exact_model_passes_on_the_full_probe():
    root = VIEWS / "probe-topology"
    if not (root / "lidc" / "val.jsonl").is_file():
        pytest.skip("no probe-topology view")
    res = probes.score_topology(_probe_rows(root), root, B=200)
    assert res["pass"] and res["controls_correct"] == 1.0 and res["hole_recall"]["n_bins"] == [120, 120, 120]


@pytest.mark.data
def test_topology_probe_on_select_ckpt_items_subset_can_pass(monkeypatch, tmp_path):
    root = VIEWS / "probe-topology"
    if not (root / "lidc" / "val.jsonl").is_file():
        pytest.skip("no probe-topology view")
    from types import SimpleNamespace

    exact, seen = {r["item_id"]: r for r in _probe_rows(root)}, {}

    def sample_view(a, branch, rows, out, k):  # the sampler, replaced by an exact model
        seen[out.name] = len(rows)
        return {r["key"]: exact.get(r["key"], {"item_id": r["key"], "samples": []}) for r in rows}

    monkeypatch.setattr(sft_select, "sample_view", sample_view)
    monkeypatch.setattr(sft_select, "score_view", lambda *a: {"summary": {}, "panel": {}, "iou": {}})
    monkeypatch.setattr(probes, "score_precision", lambda rows, root: {"rows": len(rows)})
    topo = probes.score_topology
    monkeypatch.setattr(probes, "score_topology", lambda rows, root: topo(rows, root, B=200))
    a = SimpleNamespace(run=tmp_path, items=256, k=8, shifted=VIEWS / "lidc-val-shifted", rescaled=None,
                        rotated=None, prepared=None, probes=VIEWS, dataset="lidc")  # --items 256 (the HF path)
    rep = sft_select.judge(a, {"name": "e1", "vision": {}, "panel": {}}, tmp_path, baselines)
    assert seen["probe-topology"] == 680 and rep["probes"]["topology"]["pass"]  # every probe record, a perfect model


@pytest.mark.data
def test_select_ckpt_probe_rows_are_the_probe_items():
    pytest.importorskip("datasets")
    root = VIEWS / "probe-topology"
    if not (root / "lidc" / "val.jsonl").is_file():
        pytest.skip("no probe-topology view")
    rows = evaluate.eval_rows(root, "lidc", rl_only=False)  # select_ckpt.judge's probe rows (probe records: rl false)
    assert len(rows) == len(probes.records(root)) == 680


# ================================================================ signal: R-B


def _branch(name, depth, sd, iou=0.6, n=40, rstd=None):
    ids = [f"lidc/1.2.{i}/n000/z000" for i in range(n)]
    per = {"sd_geom": [float(sd)] * n, "r_std": [float(sd if rstd is None else rstd)] * n}
    return {"name": name, "depth": depth, "eligible": True,
            "panel": {"ids": ids, "per_item": per, "iou_tc": [[iou] * 8] * n}}


def test_rb_single_and_no_eligible_branch():
    one = _branch("e1", 1, 0.1)
    assert S.select_rb([one]) == ("e1", {"eligible": ["e1"], "pairs": [], "selected_r_std": "e1", "disagree": False})
    name, why = S.select_rb([one | {"eligible": False}])
    assert name is None and why["eligible"] == []


def test_rb_veto_fires_at_plus_16_percent_and_iou_ci_touching_zero():
    deep, shallow = _branch("e2", 2, 0.10), _branch("e1", 1, 0.116)
    name, why = S.select_rb([deep, shallow])
    assert name == "e1" and why["pairs"][0]["iou_median"] == [0.0, 0.0, 0.0]  # hi = 0: not a significant loss


def test_rb_veto_does_not_fire_when_ci_lower_bound_is_zero(monkeypatch):
    real = S.paired

    def fake(a, b, key, **kw):
        return [0.05, 0.0, 0.1] if key == "sd_geom" else real(a, b, key, **kw)
    monkeypatch.setattr(S, "paired", fake)
    name, why = S.select_rb([_branch("e2", 2, 0.10), _branch("e1", 1, 0.2)])
    assert name == "e2" and why["pairs"][0]["sd_geom"]["fires"] is False


def test_rb_veto_fires_at_exactly_15_percent():
    name, why = S.select_rb([_branch("e2", 2, 0.30), _branch("e1", 1, 0.345)])
    assert why["pairs"][0]["sd_geom"]["fires"] and name == "e1"


def test_rb_refuses_duplicate_depths():
    with pytest.raises(ValueError):
        S.select_rb([_branch("a", 1, 0.1), _branch("b", 1, 0.2)])


# ================================================================ prereg / C2


def test_c2_boundaries_exact():
    assert P.c2_seed_clause([0.5, 0.5, 0.3]) and not P.c2_seed_clause([0.5, 0.5, 0.2999])
    assert not P.c2_seed_clause([0.5, 0.49, 0.49])
    assert P.c2_seed_clause([0.5, 0.5, 0.5, 0.31]) and not P.c2_seed_clause([0.5, 0.4, 0.4, 0.9])
    assert not P.c2_seed_clause([0.9, 0.9, 0.9, 0.29])
    with pytest.raises(ValueError):
        P.c2_seed_clause([0.9, 0.9])
    ok = dict.fromkeys("12345", True)
    assert P.c2_outcome(0.30, ok | {"1": False}) == ("partial", "R-x8")
    assert P.c2_outcome(0.2999, ok) == ("fail", "R-ctrl")
    assert P.c2_outcome(0.50, ok) == ("pass", None)
    assert P.c2_outcome(0.50, ok | {"4": False}) == ("partial", None)
    assert P.c2_clauses(0.5, 0.30, [0.5] * 3, [0.3] * 3, 0.05, 0.90)["1"] is False  # lo must be > 0.30
    assert all(P.c2_clauses(0.5, 0.3001, [0.5] * 3, [0.3] * 3, 0.05, 0.90).values())


def test_fourth_seed_eligibility_boundaries():
    other = {k: True for k in "1345"}
    assert P.fourth_seed_eligible([0.5, 0.49, 0.6], 0.55, other)
    assert not P.fourth_seed_eligible([0.5, 0.5, 0.5], 0.5, other)  # 3 of 3: no single seed decides
    assert not P.fourth_seed_eligible([0.5, 0.49, 0.29], 0.45, other)  # a seed below the floor fails anyway
    assert not P.fourth_seed_eligible([0.5, 0.49, 0.6, 0.7], 0.55, other)


# ================================================================ gates <-> prereg <-> select_ckpt


def test_gate_names_agree_with_prereg():
    rep = gates.phase2_entry({})
    assert set(gates.REQUIREMENTS) == set(P.C1_GATES) == set(rep["gates"])
    assert set(rep["passed"]) == set(P.C1_GATES) | {"G-probe"}
    assert set(P.REPORT_ONLY) <= set(rep["reported"])
    assert P.RB["eligible"]["parse"] == gates.FORMAT and P.RB["eligible"]["legal"] == gates.LEGAL
    assert P.RB["eligible"]["dup_max"] == gates.DUPLICATES and P.RB["eligible"]["sigma_share"] == gates.SIGNAL


def test_empty_report_passes_no_c1_gate():
    rep = gates.phase2_entry({"vision": {"pass": True}})  # e.g. eval_run --view centred, or a missing shifted view
    assert rep["pass"] is False
    assert P.c1_stop([("qwen", rep["passed"])])[0] is None


def _good_report(win_share=0.0):
    ds = {"format": 1.0, "legal": 1.0, "iou_tc_pass": 0.9, "group": 8, "signal_share": 0.9, "dup_rate": 0.0,
          "objects": {"nodule": {"iou_tc_median": 0.8}}, "truncated": 0.0, "tokens": None}
    view = {"summary": {"datasets": {"lidc": ds}}}
    prior = {"lidc": {"nodule": {"diff": 0.1, "ci": [0.05, 0.15]}}}
    centred = {"lidc": {"nodule": {"win_share": win_share}}}
    return {"views": {"shifted": view, "rescaled": view},
            "prior": {"shifted": prior, "rescaled": prior, "centred": centred},
            "probes": {"precision": {"pass": True}, "topology": {"pass": True, "hole_recall": {"pooled": [0.1]}}},
            "vision": {"pass": True}, "g_probe": {"verdict": "INCONCLUSIVE"}}


def test_prior_a_stays_report_only_and_probe_warning():
    rep = gates.phase2_entry(_good_report(win_share=0.0), ["lidc"])
    assert rep["pass"] and rep["reported"]["G-prior-a"] == {"lidc/nodule": 0.0}
    assert P.c1_stop([("qwen", rep["passed"])]) == ("qwen", ["G-probe warning"])


def test_prior_a_threshold_is_applied_when_set():
    rep = gates.phase2_entry(_good_report(win_share=0.9), ["lidc"], prior_a_threshold=0.6)
    assert P.c1_stop([("qwen", rep["passed"])], prior_a_threshold=0.6)[0] == "qwen"
    rep = gates.phase2_entry(_good_report(win_share=0.5), ["lidc"], prior_a_threshold=0.6)
    assert not rep["pass"] and P.c1_stop([("qwen", rep["passed"])], prior_a_threshold=0.6)[0] is None


def test_nan_and_missing_values_fail_their_gates():
    r = _good_report()
    ds = r["views"]["shifted"]["summary"]["datasets"]["lidc"]
    ds |= {"format": float("nan"), "dup_rate": float("nan")}
    r["prior"]["rescaled"]["lidc"]["nodule"]["ci"] = [float("nan"), 1.0]
    rep = gates.phase2_entry(r, ["lidc"])
    assert not rep["passed"]["G-format"] and not rep["passed"]["G-signal"] and not rep["passed"]["G-prior-b"]


def test_lidc_dropped_reads_no_lidc_value():
    r = _good_report()
    riga = {"format": 1.0, "legal": 1.0, "iou_tc_pass": 0.9, "group": 8, "signal_share": 0.9, "dup_rate": 0.0,
            "objects": {"disc": {"iou_tc_median": 0.9}, "cup": {"iou_tc_median": 0.8}}}
    for v in ("shifted", "rescaled"):
        r["views"][v]["summary"]["datasets"] = {"riga": riga, "lidc": {"format": 0.0}}
        r["prior"][v] = {"riga": {n: {"ci": [0.01, 0.1]} for n in ("disc", "cup")}}
    r["probes"] = {"precision": {"pass": False}, "topology": {"pass": False}}
    rep = gates.phase2_entry(r, ["riga"])
    assert rep["pass"] and rep["datasets"] == {"riga": True}
    assert not any("lidc" in k for g in rep["gates"].values() for k in g["values"])

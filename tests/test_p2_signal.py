"""Phase 2 signal panel, pre-registered rules and the C2 analysis (PLAN §7 M11 and the M8 stats parts)."""

import itertools
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest
from scipy import stats

from experiments.phase2 import c2_analyze as C2
from experiments.phase2 import prereg as P
from objdel.eval import signal as S

ROOT = Path(__file__).resolve().parents[1]
FIX = json.loads((Path(__file__).parent / "p2_signal_fixture.json").read_text())  # B2 shifted LFM, 10 items
SHALLOW, DEEP = "lfm2.5-vl-450m_s0", "lfm2.5-vl-450m_s0_n150"


def rand_groups(rng, n=40, k=8):
    out = []
    for i in range(n):
        legal = (rng.random(k) < rng.random()).astype(float)
        R = rng.random(k) * legal
        out.append(
            {
                "item_id": f"lidc/1.2.{i // 3}/n000/z{i:03d}",
                "R": R.tolist(),
                "legal": legal.tolist(),
                "iou_tc": (R * 1.5).clip(0, 1).tolist(),
                "advantage": (R - R.mean()).tolist(),
                "sigma": float(R.std() > 0.05),
                "dup": 0.0,
                "distinct": 1.0,
            }
        )
    return out


# ---------------------------------------------------------------- variance split


def test_split_identity_random_groups():
    rng = np.random.default_rng(0)
    for _ in range(200):
        k = int(rng.integers(2, 17))
        legal = rng.random(k) < rng.random()
        R = np.where(legal, rng.random(k), 0.0)
        s = S.split(R[None], legal[None])
        assert abs(s["var"][0] - s["var_legality"][0] - s["var_geom"][0]) < 1e-12
        assert s["var"][0] == pytest.approx(R.var(), abs=1e-15)


def test_sd_geom_edge_cases():
    legal = np.array([[1, 1, 1, 0, 0, 1, 1, 0]], float)
    R = np.where(legal > 0, 0.4, 0.0)  # every legal sample equal
    assert S.split(R, legal)["sd_geom"][0] == 0.0
    R = np.random.default_rng(1).random((5, 8))  # all legal: sd_geom = r_std
    s = S.split(R, np.ones_like(R))
    np.testing.assert_allclose(s["sd_geom"], R.std(1), atol=1e-15)
    assert np.all(s["var_legality"] == 0)
    s = S.split(np.zeros((1, 8)), np.zeros((1, 8)))  # nothing legal
    assert s["sd_geom"][0] == 0 and s["var"][0] == 0 and s["geom_share"][0] == 0


# ---------------------------------------------------------------- metrics = B2


@pytest.mark.parametrize("run", [SHALLOW, DEEP])
def test_per_item_equals_b2(run):
    got = S.per_item(S.rows(FIX["runs"][run]))
    for k, x in FIX["expected"][run]["per_item"].items():
        np.testing.assert_allclose(got[k], x, atol=1e-12, err_msg=k)


@pytest.mark.parametrize("run", [SHALLOW, DEEP])
def test_panel_cis_equal_b2(run):
    p = S.panel(S.rows(FIX["runs"][run]))
    for k, v in FIX["expected"][run]["ci"].items():
        assert [round(x, 4) for x in p["summary"][k]] == v, k
    assert p["items"] == 10 and p["nodules"] == 10


def test_paired_equals_b2():
    a, b = (S.panel(S.rows(FIX["runs"][r])) for r in (DEEP, SHALLOW))
    assert [round(x, 4) for x in S.paired(a, b, "sd_geom")] == FIX["expected"]["paired_sd_geom_150_minus_36"]


def test_live_b2_functions_and_group_scoring():
    """B2's score / analyze / posthoc imported live (they pull torch through lab B1) on fresh random groups, and
    S.group = B2's score_view on synthetic items scored through objdel.eval.metrics."""
    pytest.importorskip("torch")
    from experiments.phase2.lab.B2 import analyze as A
    from experiments.phase2.lab.B2 import posthoc as PH
    from experiments.phase2.lab.B2.score import score_view
    from objdel.reward import ExpertCache
    from objdel.targets import build_target

    g = rand_groups(np.random.default_rng(3))
    v = {k: [x[k] for x in g] for k in g[0]}
    got = S.per_item(g)
    for k, x in (A.group_metrics(v) | PH.parts(v)).items():
        np.testing.assert_allclose(got[k], x, atol=1e-12, err_msg=k)

    yy, xx = np.ogrid[:50, :50]
    disk = lambda cy, cx, r: (yy - cy) ** 2 + (xx - cx) ** 2 <= r * r  # noqa: E731
    raters = {"nodule": np.stack([disk(25, 25, 12), disk(25, 27, 14), disk(24, 25, 11)])}
    ma = {"nodule": ExpertCache.build(raters["nodule"]).min_area}
    texts = [build_target({"nodule": disk(25, 25, r)}, ("nodule",), min_area=ma).text for r in (9, 11, 12, 13, 15)]
    samples = texts + ["garbage", texts[0], "{"]
    caches = {"nodule": ExpertCache.build(raters["nodule"])}
    item = {"item_id": "lidc/x/n000/z001", "caches": caches, "W": 50, "H": 50, "names": ("nodule",)}
    raw = score_view([{"item_id": item["item_id"], "samples": samples, "truncated": [False] * 8}], [item])
    from objdel.eval import metrics

    m = metrics(samples, {"nodule": raters["nodule"]}, 50, 50, ("nodule",), [False] * 8, caches)
    mine = S.group(item["item_id"], samples, m)
    for k in ("R", "legal", "iou_tc", "advantage", "sigma", "dup", "distinct"):
        assert mine[k] == raw[k][0], k


# ---------------------------------------------------------------- R-B


def branch(name, depth, sd, rstd, iou, eligible=True, n=60):
    """A synthetic branch panel: per-item sd_geom / r_std and IoU samples around the given levels (fixed noise)."""
    rng = np.random.default_rng(int(depth * 10))
    ids = [f"lidc/1.2.{i // 2}/n000/z{i:03d}" for i in range(n)]
    per = {
        "sd_geom": (sd + 0.002 * rng.standard_normal(n)).tolist(),
        "r_std": (rstd + 0.002 * rng.standard_normal(n)).tolist(),
    }
    return {
        "name": name,
        "depth": depth,
        "eligible": eligible,
        "panel": {"ids": ids, "per_item": per, "iou_tc": (iou + 0.01 * rng.standard_normal((n, 8))).tolist()},
    }


def test_rb_deepest_eligible_and_ties_stay_deep():
    bs = [branch("e0.5", 0.5, 0.10, 0.10, 0.6), branch("e1", 1, 0.10, 0.10, 0.6), branch("e2", 2, 0.10, 0.10, 0.6)]
    name, why = S.select_rb(bs)
    assert name == "e2" and not why["disagree"]
    bs[2]["eligible"] = False
    assert S.select_rb(bs)[0] == "e1"
    assert S.select_rb([dict(b, eligible=False) for b in bs])[0] is None


def test_rb_veto_needs_ci_and_margin_and_no_iou_loss():
    deep = branch("e2", 2, 0.10, 0.10, 0.6)
    assert S.select_rb([branch("e1", 1, 0.13, 0.13, 0.6), deep])[0] == "e1"  # +30 %, CI > 0, IoU equal
    assert S.select_rb([branch("e1", 1, 0.11, 0.11, 0.6), deep])[0] == "e2"  # CI > 0 but +10 % < 15 %
    assert S.select_rb([branch("e1", 1, 0.13, 0.13, 0.5), deep])[0] == "e2"  # significant IoU loss
    spike = branch("e1", 1, 0.10, 0.10, 0.6)  # +40 % in the mean, all from one item: the CI includes 0
    sd = np.array(deep["panel"]["per_item"]["sd_geom"])
    spike["panel"]["per_item"]["sd_geom"] = (sd + np.r_[2.4, np.zeros(59)]).tolist()
    name, why = S.select_rb([spike, deep])
    assert name == "e2" and why["pairs"][0]["sd_geom"]["rel"] > 0.15 and why["pairs"][0]["sd_geom"]["ci"][1] <= 0


def test_rb_pairwise_downward_and_disagreement_flag():
    bs = [branch("e0.5", 0.5, 0.20, 0.20, 0.6), branch("e1", 1, 0.13, 0.13, 0.6), branch("e2", 2, 0.10, 0.10, 0.6)]
    name, why = S.select_rb(bs)
    assert name == "e0.5" and [p["shallow"] for p in why["pairs"]] == ["e1", "e0.5"]
    bs = [branch("e1", 1, 0.13, 0.10, 0.6), branch("e2", 2, 0.10, 0.10, 0.6)]  # sd_geom fires, r_std not
    name, why = S.select_rb(bs)
    assert name == "e1" and why["selected_r_std"] == "e2" and why["disagree"]


# ---------------------------------------------------------------- Spearman


def test_spearman_equals_scipy():
    rng = np.random.default_rng(5)
    x = rng.random(128)
    y = x + rng.random(128)
    y[:5] = y[0]  # ties
    rho, lo, hi = S.spearman_ci(x, y, [f"lidc/a/n{i // 2:03d}/z1" for i in range(128)])
    assert rho == pytest.approx(stats.spearmanr(x, y).statistic, abs=1e-12)
    assert lo < rho < hi


# ---------------------------------------------------------------- c2_analyze


def test_c2_stratified_bootstrap_equals_power_c2():
    from experiments.phase2.lab.B3 import power_c2

    rng = np.random.default_rng(11)
    bins = [rng.random(120) for _ in range(3)]
    ref = power_c2.pooled(bins, np.random.default_rng(0))
    assert C2.pooled_ci(bins, seed=0) == pytest.approx(ref, abs=1e-12)  # rounded to 12 decimals (BK2-4)
    assert C2.WEIGHTS == tuple(power_c2.W)


def fake_seed(rng, level, legal=0.95, fns=0.01):
    return {
        "bins": [np.clip(level + 0.1 * rng.standard_normal(120), 0, 1).tolist() for _ in range(3)],
        "false_nonsimple": fns,
        "legal": legal,
    }


def test_c2_analyze_cli_and_outcomes(tmp_path, monkeypatch, capsys):
    rng = np.random.default_rng(2)
    paths = []
    for i, level in enumerate((0.7, 0.68, 0.72)):
        p = tmp_path / f"s{i}.json"
        p.write_text(json.dumps({"probes": {"topology": fake_seed(rng, level) | {"ids": list(range(120))}}}))
        paths.append(str(p))
    out = tmp_path / "c2.json"
    monkeypatch.setattr(sys, "argv", ["c2_analyze", "--seeds", *paths, "--out", str(out)])
    C2.main()
    monkeypatch.setattr(sys, "argv", ["prereg"])
    P.main()
    assert '"pass"' in capsys.readouterr().out
    r = json.loads(out.read_text())
    assert r["outcome"] == "pass" and r["conditional"] is None and r["clauses"]["2"] and not r["seed_discordant"]
    r = C2.analyze([fake_seed(rng, v) for v in (0.7, 0.45, 0.4)])  # 1 of 3 seeds >= 0.50
    assert r["outcome"] == "partial" and r["conditional"] == "R-x8" and r["fourth_seed_eligible"]
    r = C2.analyze([fake_seed(rng, v) for v in (0.2, 0.25, 0.2)])
    assert r["outcome"] == "fail" and r["conditional"] == "R-ctrl"
    r4 = C2.analyze([fake_seed(rng, v) for v in (0.7, 0.6, 0.45, 0.65)])  # 3 of 4 >= 0.50
    assert r4["clauses"]["2"] and r4["outcome"] == "pass"


def test_scripts_help():
    for mod in ("experiments.phase2.c2_analyze", "experiments.phase2.prereg"):
        r = subprocess.run([sys.executable, "-m", mod, "--help"], capture_output=True, text=True, cwd=ROOT)
        assert r.returncode == 0 and "usage" in r.stdout.lower()


# ---------------------------------------------------------------- prereg


def test_c2_outcome_table_exhaustive():
    seen = set()
    for pooled, lo_ok, seeds_ok, bins_ok, fns_ok, legal_ok in itertools.product(
        (0.1, 0.29, 0.30, 0.4, 0.5, 0.7), *[(True, False)] * 5
    ):
        cl = {"1": pooled >= 0.5 and lo_ok, "2": seeds_ok, "3": bins_ok, "4": fns_ok, "5": legal_ok}
        out, cond = P.c2_outcome(pooled, cl)
        assert out in ("pass", "partial", "fail")
        assert cond in (None, "R-x8", "R-ctrl")
        if out == "pass":
            assert all(cl.values()) and cond is None
        elif out == "fail":
            assert pooled < 0.30 and cond == "R-ctrl"
        else:
            assert pooled >= 0.30 and not all(cl.values())
            assert cond == ("R-x8" if fns_ok else None)
        seen.add((out, cond))
    assert seen == {("pass", None), ("partial", "R-x8"), ("partial", None), ("fail", "R-ctrl")}


def test_c2_seed_clause_and_fourth_seed():
    assert P.c2_seed_clause([0.6, 0.55, 0.4]) and not P.c2_seed_clause([0.6, 0.55, 0.29])
    assert not P.c2_seed_clause([0.6, 0.45, 0.4])
    assert P.c2_seed_clause([0.6, 0.55, 0.51, 0.4]) and not P.c2_seed_clause([0.6, 0.55, 0.4, 0.45])
    with pytest.raises(ValueError):
        P.c2_seed_clause([0.6, 0.6])
    ok = {"1": True, "3": True, "4": True, "5": True}
    assert P.fourth_seed_eligible([0.6, 0.55, 0.4], 0.55, ok)  # decided by one seed's side of 0.50
    assert not P.fourth_seed_eligible([0.6, 0.55, 0.7], 0.6, ok)  # 3 of 3: robust
    assert not P.fourth_seed_eligible([0.6, 0.55, 0.4], 0.2, ok)  # fail regardless of clause (2)
    assert not P.fourth_seed_eligible([0.6, 0.55, 0.4, 0.5], 0.55, ok)  # never a fifth seed


def test_g_probe_verdicts():
    null = (-0.01, 0.01)
    assert P.g_probe_verdict((0.005, 0.03), null, 0.02, 0.9) == "PASS"
    assert P.g_probe_verdict((0.005, 0.03), null, 0.06, 0.9) == "FAIL"  # legal falls > 5 pts
    assert P.g_probe_verdict((0.005, 0.03), null, 0.0, 0.4) == "FAIL"  # sigma-share < 50 %
    assert P.g_probe_verdict((-0.05, -0.02), null, 0.0, 0.9) == "FAIL"  # entirely below the null CI
    assert P.g_probe_verdict((-0.02, 0.01), null, 0.0, 0.9) == "INCONCLUSIVE"
    rng = np.random.default_rng(0)
    r0a, r0b = rng.random(128), rng.random(128)
    ids = [f"lidc/a/n{i:03d}/z1" for i in range(128)]
    res = P.g_probe(r0a + 0.3, r0a, r0b, ids, legal0=0.95, legal20=0.94, sigma20=0.8)
    assert res["verdict"] == "PASS" and res["h0"] > 0


def test_c1_stop_rule_and_prior_a_slot():
    assert P.G_PRIOR_A_THRESHOLD is None and "G-prior-a" not in P.c1_gates()
    gates = dict.fromkeys(P.c1_gates(), True)
    rungs = [
        ("lfm", gates | {"G-geometry": False}),
        ("qwen3-vl-2b", gates | {"G-topology-b": False, "G-probe": False}),
        ("larger", gates),
    ]
    assert P.c1_stop(rungs) == ("qwen3-vl-2b", ["G-probe warning"])
    assert P.c1_stop(rungs[:1]) == (None, [])
    assert "G-prior-a" in P.c1_gates(prior_a_threshold=0.6)

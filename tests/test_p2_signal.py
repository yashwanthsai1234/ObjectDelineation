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


# ---------------------------------------------------------------- A1: cluster draws, patient key, CR1 interval

GOLD = Path(__file__).parent / "goldens" / "r5_test_items.jsonl.gz"
GOLD_META = json.loads((Path(__file__).parent / "goldens" / "r5_test_items.meta.json").read_text())


def gold_rows():
    import gzip

    with gzip.open(GOLD, "rt") as f:
        return [json.loads(line) for line in f]


def legacy_draws(ids, boot, seed):  # verbatim copy of the body before A1, frozen (P01)
    nod = [S.cluster(i) for i in ids]
    groups = [np.array([j for j, x in enumerate(nod) if x == u]) for u in dict.fromkeys(nod)]
    rng = np.random.default_rng(seed)
    return [np.concatenate([groups[g] for g in rng.integers(0, len(groups), len(groups))]) for _ in range(boot)]


def same_draws(new, old):
    return len(new) == len(old) and all(np.array_equal(a, b) and a.dtype == b.dtype
                                        for a, b in zip(new, old, strict=True))


def test_draws_default_key_is_bit_identical_to_legacy():
    ids = [r["item_id"] for r in gold_rows()]
    for boot, seed in [(50, 0), (200, 3)]:
        assert same_draws(S.draws(ids, boot, seed), legacy_draws(ids, boot, seed))
    toy = ["lidc/s/n000/z1", "lidc/s/n000/z2", "x", "lidc/s/n001/z1", "lidc/s/n000/z3", "y", "lidc/t/n000/z1"]
    assert same_draws(S.draws(toy, 50, 0), legacy_draws(toy, 50, 0))
    rng = np.random.default_rng(0)
    for _ in range(50):
        ids = [f"lidc/s{u}/n{u % 3:03d}/z{z}" for u in range(int(rng.integers(1, 30)))
               for z in range(int(rng.integers(1, 7)))]  # fmt: skip
        ids += [f"other{j}" for j in range(int(rng.integers(0, 4)))]
        rng.shuffle(ids)
        for boot, seed in [(1, 0), (7, 1), (100, 12345)]:
            assert same_draws(S.draws(ids, boot, seed), legacy_draws(ids, boot, seed))


def test_one_item_per_cluster_equals_integers_rows():
    for n in (1, 2, 10, 17, 401):
        want = np.random.default_rng(5).integers(0, n, (100, n))
        assert np.array_equal(np.stack(S.draws([f"x{i}" for i in range(n)], 100, 5)), want)
        assert np.array_equal(np.stack(S.draws([f"lidc/s{i}/n000/z1" for i in range(n)], 100, 5)), want)


def test_draws_with_a_key_keep_whole_clusters_in_first_appearance_order():
    ids = ["a1", "a2", "b1", "c1", "c2", "c3"]
    blocks = [[0, 1], [2], [3, 4, 5]]
    want = np.concatenate([blocks[i] for i in np.random.default_rng(0).integers(0, 3, 3)])
    assert np.array_equal(S.draws(ids, 1, 0, key=lambda s: s[0])[0], want)
    # interleaved ids: a drawn block is the cluster's indices in ascending order
    for d in S.draws(["a", "b", "a", "b"], 20, 1, key=lambda s: s):
        assert all(list(d[j:j + 2]) in ([0, 2], [1, 3]) for j in range(0, len(d), 2))
    # random clustered ids: every draw is a concatenation of whole clusters
    rng = np.random.default_rng(2)
    sizes = rng.integers(1, 10, 12)
    ids = [f"c{c}/{j}" for c, n in enumerate(sizes) for j in range(n)]
    rng.shuffle(ids)
    members = {}
    for j, i in enumerate(ids):
        members.setdefault(i.split("/")[0], []).append(j)
    for d in S.draws(ids, 30, 0, key=lambda s: s.split("/")[0]):
        assert 0 <= d.min() and d.max() < len(ids)
        rest = list(d)
        while rest:
            block = members[ids[rest[0]].split("/")[0]]
            assert rest[: len(block)] == block
            rest = rest[len(block):]
    # patient key: two nodules of one patient form one cluster; no data folder is read
    pids = ["lidc/s1/n000/z1", "lidc/s1/n000/z2", "lidc/s1/n001/z5", "lidc/s2/n000/z9"]
    groups = dict(zip(pids, ["P1", "P1", "P1", "P2"], strict=True))
    seen = set()
    for d in S.draws(pids, 50, 0, key=groups.__getitem__):
        rest = list(d)
        while rest:
            block = [0, 1, 2] if rest[0] == 0 else [3]
            assert rest[: len(block)] == block
            rest = rest[len(block):]
            seen.add(tuple(block))
    assert seen == {(0, 1, 2), (3,)}
    for d in S.draws(pids, 5, 0, key=lambda s: "one"):
        assert np.array_equal(d, np.arange(4))


def test_draws_edges():
    with pytest.raises(ValueError, match="no ids"):
        S.draws([], 3)
    assert all(np.array_equal(d, [0, 1]) for d in S.draws(["lidc/s/n0/z1", "lidc/s/n0/z2"], 5))
    assert [d.tolist() for d in S.draws(["a"], 4)] == [[0]] * 4
    assert S.draws(["a", "b"], 0) == []


def test_cluster_parsing():
    assert S.cluster("lidc/1.3.6.1/n000/z185") == "lidc/1.3.6.1/n000"
    assert S.cluster("riga/x/img1") == "riga/x/img1"
    assert S.cluster("lidc/") == "lidc"  # pin


def test_draws_are_deterministic_and_leave_the_global_rng_alone():
    ids = [f"lidc/s{i % 5}/n000/z{i}" for i in range(20)]
    state = np.random.get_state()
    a, b, c = S.draws(ids, 10, 0), S.draws(ids, 10, 0), S.draws(ids, 10, 1)
    after = np.random.get_state()
    assert same_draws(a, b) and not same_draws(a, c)
    assert state[0] == after[0] and np.array_equal(state[1], after[1]) and state[2:] == after[2:]


def test_cluster_ci_wider_than_item_ci_under_within_cluster_correlation():
    # math-stats §11 item 4: >= 5 clusters, equal values inside each cluster, unequal sizes
    sizes, vals = [1, 2, 3, 4, 5, 6], [0.1, 0.9, 0.3, 0.7, 0.2, 0.8]
    ids = [f"lidc/s{c}/n000/z{j}" for c, n in enumerate(sizes) for j in range(n)]
    x = np.repeat(vals, sizes)
    idx = np.random.default_rng(0).integers(0, len(x), (2000, len(x)))
    item = np.percentile(x[idx].mean(1), [2.5, 97.5])
    nod = S.ci(x, S.draws(ids, 2000, 0))
    assert nod[2] - nod[1] > item[1] - item[0]


def test_patient_ci_wider_than_nodule_ci_when_one_patient_shares_values():
    # U53 with 3 nodules per patient, not 2: with 4 nodules, P(all drawn from one patient) = 1/16 > 2.5 %, so both
    # CIs are [0.2, 0.8] (math-stats §11 item 4, the few-cluster case)
    ids = [f"lidc/{p}/n00{n}/z{s}" for p in ("p1", "p2") for n in (0, 1, 2) for s in (0, 1)]
    x = np.array([0.2] * 6 + [0.8] * 6)
    groups = {i: i.split("/")[1] for i in ids}
    Dp = S.draws(ids, 2000, 0, key=groups.__getitem__)
    pat, nod = S.ci(x, Dp), S.ci(x, S.draws(ids, 2000, 0))
    assert pat[2] - pat[1] > nod[2] - nod[1]
    means = np.array([x[d].mean() for d in Dp])
    assert pat[1:] == [float(np.quantile(means, 0.025)), float(np.quantile(means, 0.975))]


def test_cluster_t_equal_sizes_is_the_t_interval_of_the_patient_means():
    pm = np.array([0.1, 0.2, 0.3, 0.4, 0.5])
    x = np.repeat(pm, 2) + np.tile([-0.05, 0.05], 5)
    g = np.repeat(["a", "b", "c", "d", "e"], 2).tolist()
    m, lo, hi = S.cluster_t(x, g, 0.95)
    h = stats.t.ppf(0.975, 4) * pm.std(ddof=1) / np.sqrt(5)
    assert abs(m - pm.mean()) < 1e-12 and abs(lo - (pm.mean() - h)) < 1e-12 and abs(hi - (pm.mean() + h)) < 1e-12


def test_cluster_t_is_the_cr1_interval_centred_on_the_item_mean():
    x = np.array([0.2, 0.4, 0.6, 0.9, 0.9, 0.8, 0.7, 1.0])
    g = ["a", "b", "b", "c", "c", "c", "c", "c"]
    m, lo, hi = S.cluster_t(x, g, 0.95)
    assert round(m, 12) == round(x.mean(), 12) == 0.6875  # not the patient mean 0.52
    Sp, Np, P = np.array([0.2, 1.0, 4.3]), np.array([1.0, 2.0, 5.0]), 3  # per-patient sums and counts, by hand
    se = np.sqrt(P / (P - 1) * ((Sp - m * Np) ** 2).sum() / Np.sum() ** 2)
    assert abs((hi - lo) / 2 - stats.t.ppf(0.975, P - 1) * se) < 1e-12
    assert abs((hi + lo) / 2 - m) < 1e-12
    with pytest.raises(ValueError, match="2 patients"):
        S.cluster_t(x, ["a"] * 8)


def test_coverage_and_width_order_of_the_three_cis():
    """P03. Nodule effect N(0, 1) shared by 6 slices of 200 nodules: the nodule CI covers the true mean 0, the item
    CI does not. The draws depend on the ids only, so one set serves every replication."""
    rng = np.random.default_rng(7)
    ids = [f"lidc/s{c}/n000/z{j}" for c in range(200) for j in range(6)]
    Dn = np.stack(S.draws(ids, 400, 0))
    Di = np.random.default_rng(0).integers(0, len(ids), (400, len(ids)))
    cover_n = cover_i = 0
    widths = []
    for _ in range(300):
        x = np.repeat(rng.normal(0, 1, 200), 6) + rng.normal(0, 0.1, 1200)
        n_lo, n_hi = np.percentile(x[Dn].mean(1), [2.5, 97.5])
        i_lo, i_hi = np.percentile(x[Di].mean(1), [2.5, 97.5])
        cover_n += n_lo <= 0 <= n_hi
        cover_i += i_lo <= 0 <= i_hi
        widths.append((n_hi - n_lo) / (i_hi - i_lo))
    assert 0.90 <= cover_n / 300 <= 0.99 and cover_i / 300 < 0.80 and np.mean(widths) > 1.8
    ids = [f"lidc/p{p}x/n{n:03d}/z{s}" for p in range(40) for n in range(3) for s in range(4)]
    groups = {i: i.split("/")[1] for i in ids}
    for seed in range(5):
        r = np.random.default_rng(seed)
        x = np.repeat(r.normal(0, 1, 40), 12) + np.repeat(r.normal(0, 0.5, 120), 4) + r.normal(0, 0.1, 480)
        w = {}
        for name, D in [("patient", S.draws(ids, 400, seed, key=groups.__getitem__)),
                        ("nodule", S.draws(ids, 400, seed)),
                        ("item", list(np.random.default_rng(seed).integers(0, 480, (400, 480))))]:
            _, lo, hi = S.ci(x, D)
            w[name] = hi - lo
        assert w["patient"] > w["nodule"] > w["item"]


def test_draws_speed_on_the_golden_ids():
    import time

    ids = [r["item_id"] for r in gold_rows()]
    t = time.perf_counter()
    legacy_draws(ids, 2000, 0)
    old = time.perf_counter() - t
    t = time.perf_counter()
    S.draws(ids, 2000, 0)
    new = time.perf_counter() - t
    assert new < 2.0 and new <= 1.5 * old + 0.05


def test_gate_a_g1_exact_numbers_on_r5_test_items():
    """I01: the item and nodule CIs of R5's test IoU_tc are reproduced, and the patient CI is the gate CI; `==` on the
    numpy version and machine that computed them, 1e-12 elsewhere."""
    import platform

    from experiments.phase2 import final_eval as FE

    rows = gold_rows()
    assert len(rows) == 2200 and len({r["item_id"] for r in rows}) == 2200
    assert len({S.cluster(r["item_id"]) for r in rows}) == 401 and len({r["group_id"] for r in rows}) == 158
    assert all(0 <= r["iou_tc"] <= 1 for r in rows)
    items = [{"item_id": r["item_id"], **{m: r["iou_tc"] for m in FE.METRICS}} for r in rows]
    exact = GOLD_META == {"numpy": np.__version__, "machine": platform.machine()}

    def same(got, want):
        return got == want if exact else all(abs(a - b) <= 1e-12 for a, b in zip(got, want, strict=True))

    s = FE.summarize(items, draws=2000, seed=0)
    assert same(s["iou_tc"]["ci95_item"], [0.6334010210209712, 0.6561989221793691])
    assert same(s["iou_tc"]["ci95_nodule"], [0.6279017988891661, 0.6620867640136667])
    assert s["miss_rate"]["mean"] == 0.10818181818181818
    assert (s["items"], s["nodules"]) == (2200, 401)
    p = FE.summarize(items, groups={r["item_id"]: r["group_id"] for r in rows}, draws=2000, seed=0)
    assert same(p["iou_tc"]["ci95_patient"], [0.6252288200704239, 0.6644106494958667])
    assert p["patients"] == 158 and p["gate_ci"] == "ci95_patient"
    assert p["iou_tc"]["ci95_nodule"] == s["iou_tc"]["ci95_nodule"]


def test_summary_speed_memory_and_two_draw_calls(monkeypatch):
    """Q01 (A1 part): the golden's 2,200 ids with every metric: one nodule and one patient draw set, reused."""
    import time
    import tracemalloc

    from experiments.phase2 import final_eval as FE

    rows = gold_rows()
    items = [{"item_id": r["item_id"], **{m: r["iou_tc"] for m in FE.METRICS}} for r in rows]
    groups = {r["item_id"]: r["group_id"] for r in rows}
    calls = []
    real = S.draws
    monkeypatch.setattr(S, "draws", lambda *a, **k: calls.append(1) or real(*a, **k))
    tracemalloc.start()
    t = time.perf_counter()
    FE.summarize(items, groups)
    seconds = time.perf_counter() - t
    peak = tracemalloc.get_traced_memory()[1]
    tracemalloc.stop()
    assert len(calls) == 2 and seconds < 20 and peak < 600e6

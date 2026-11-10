"""M4b probe tests (PLAN §7 M4, probe parts): objdel.data.probes and the experiments/phase2 probes / baselines CLIs.

Synthetic tests build small probe directories on synthetic lung patches; the `data` tests use the real store."""

import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import probe_models as M
import pytest

from objdel.data import probes as P

ROOT = Path(__file__).resolve().parents[1]
FIXTURES = ROOT / "tests" / "fixtures"


def _patches(n=6, side=125, seed=0):
    rng = np.random.default_rng(seed)
    groups = [f"LIDC-IDRI-{i:04d}" for i in range(n)]
    return [(g, rng.normal(-850, 40, (side, side)).astype(np.int16)) for g in groups]


def _diams(n=4000, seed=0):  # log-normal fitted to the LIDC train p10 / p50 (6.9 / 12.2 px)
    return np.exp(np.random.default_rng(seed).normal(np.log(12.2), 0.444, n))


# hole / two-part records with the real p50s (hole_census.json): part diam 24 px, wall 1.24 px, offset 0.46
HOLES = [
    {"part_diam": d, "area": a, "wall": w, "offset": o}
    for d, a, w, o in [
        (24, 6, 1.2, 0.45),
        (22, 8, 0.41, 0.5),
        (26, 5, 1.0, 0.4),
        (25, 15, 1.5, 0.4),
        (28, 22, 2.0, 0.3),
        (24, 12, 1.24, 0.46),
        (34, 40, 2.5, 0.3),
        (40, 55, 3.0, 0.25),
        (30, 35, 1.5, 0.35),
    ]
]
MULTIS = [{"diam": d, "part_ratio": r, "gap": g} for d, r, g in [(18, 0.57, 6.8), (22, 0.4, 3.0), (16, 0.8, 12)]]


def _rows(root, model, k=1, seed=7):
    rng = np.random.default_rng(seed)
    return [{"item_id": r["item_id"], "samples": [model(r, rng) for _ in range(k)]} for r in P.records(root)]


@pytest.fixture(scope="module")
def precision_dir(tmp_path_factory):
    out = tmp_path_factory.mktemp("probe-precision")
    P.build_precision(out, _patches(), _diams(), n=500, seed=2026)
    return out


@pytest.fixture(scope="module")
def topology_dir(tmp_path_factory):
    out = tmp_path_factory.mktemp("probe-topology")
    P.build_topology(out, _patches(), HOLES, MULTIS, _diams(), seed=3, per_bin=8, n_two=8, n_simple=8)
    return out


def test_layout_is_prepared_layout(precision_dir):
    recs = P.records(precision_dir)
    assert len(recs) == 500 and (precision_dir / "manifest.json").is_file()
    r = recs[0]
    assert r["item_id"].startswith("lidc/") and r["view"] == {"offset": [0, 0], "side": r["W"]}
    img = precision_dir / "images" / (r["item_id"] + ".png")
    assert img.is_file() and M.image_sha(img) == r["image_sha256"]
    man = json.loads((precision_dir / "manifest.json").read_text())
    assert man["sets"]["lidc/val"]["items"] == 500 and "sha256" in man["files"]["lidc/val.jsonl"]


def test_exact_on_ground_truth(precision_dir):
    for r in P.records(precision_dir)[:200]:
        gt = P.gt_mask(r)
        s = P.surface(gt, gt)
        assert s["assd"] == 0 and s["hd95"] == 0 and s["nsd1"] == 1
        a = P.score_answer(r, r["targets"][0])
        assert a["legal"] and a["iou"] == 1 and a["assd"] == 0 and a["hd95"] == 0 and a["centroid"] == 0


def test_surface_equals_brute_force(precision_dir):
    rng = np.random.default_rng(0)
    for r in P.records(precision_dir)[:60]:
        gt = P.gt_mask(r)
        pred = P.gt_mask({**r, "targets": [M.model_noise(1.0)(r, rng)]})
        if not pred.any():
            continue
        a, b = np.argwhere(P.boundary(gt)), np.argwhere(P.boundary(pred))
        D = np.sqrt(((a[:, None, :] - b[None, :, :]) ** 2).sum(-1))
        dg, dp = D.min(1), D.min(0)
        s = P.surface(gt, pred)
        assert abs(s["assd"] - np.concatenate([dg, dp]).mean()) < 1e-9
        assert abs(s["hd95"] - max(np.percentile(dg, 95), np.percentile(dp, 95))) < 1e-9


def test_calibration_disc_iou(precision_dir):
    iou = [P.score_answer(r, P.model_centroid_disc(r, None)).get("iou", 0.0) for r in P.records(precision_dir)]
    assert abs(np.median(iou) - 0.806) <= 0.03


@pytest.mark.parametrize(
    ("name", "should_pass"),
    [
        ("exact", True),
        ("sigma_0.25", True),
        ("sigma_0.5", True),
        ("sigma_1", False),
        ("sigma_2", False),
        ("centroid_disc", False),
        ("centred_prior", False),
        ("shift_1px", False),
    ],
)
def test_precision_discrimination(precision_dir, name, should_pass):
    models = {
        "exact": M.model_exact,
        "sigma_0.25": M.model_noise(0.25),
        "sigma_0.5": M.model_noise(0.5),
        "sigma_1": M.model_noise(1.0),
        "sigma_2": M.model_noise(2.0),
        "centroid_disc": P.model_centroid_disc,
        "centred_prior": M.model_centred_prior({50: 8.0, 100: 15.0, 125: 20.0}),
        "shift_1px": M.model_shift,
    }
    res = P.score_precision(_rows(precision_dir, models[name]), precision_dir)
    assert res["n_items"] >= 200
    assert res["pass"] is should_pass, res
    if name == "centroid_disc":  # passes every geometry check; only the paired win rule fails it (B3 P2)
        assert res["win_share"] == 0 and res["checks"]["legal"] and res["checks"]["iou"]


def test_one_px_shift_detected(precision_dir):
    for r in P.records(precision_dir)[:100]:
        a = P.score_answer(r, M.model_shift(r, np.random.default_rng(0)))
        assert a["assd"] > 0 and abs(a["centroid"] - 1) < 1e-9


def test_topology_probe_composition(topology_dir):
    recs = P.records(topology_dir)
    fam = [r["probe"]["family"] for r in recs]
    assert fam.count("hole") == 24 and fam.count("two") == 8 and fam.count("simple") == 8
    holes = [r for r in recs if r["probe"]["family"] == "hole"]
    assert [sum(r["probe"]["bin"] == b for r in holes) for b in range(3)] == [8, 8, 8]
    for r in holes:  # every hole survives the size rule: signature (1,) raw and at min_area(|GT|)
        gt = P.gt_mask(r)
        assert P.signature(gt) == (1,) == P.signature(gt, P.min_area(int(gt.sum())))
        assert P.bin_of(r["probe"]["hole_area"]) == r["probe"]["bin"]


@pytest.mark.parametrize(
    ("name", "recall", "controls"), [("exact", 1.0, 1.0), ("a0", 0.0, 1.0), ("centroid_disc", 0.0, 1.0)]
)
def test_topology_sanity(topology_dir, name, recall, controls):
    model = {"exact": M.model_exact, "a0": M.model_a0, "centroid_disc": P.model_centroid_disc}[name]
    res = P.score_topology(_rows(topology_dir, model, k=2), topology_dir)
    assert res["hole_recall"]["bins"] == [recall] * 3 and res["hole_recall"]["pooled"][0] == pytest.approx(recall)
    assert res["controls_correct"] == controls and res["false_nonsimple"] == 0 and res["pass"] is True
    assert res["part_recall"] == {"exact": 1.0, "a0": 1.0, "centroid_disc": 0.0}[name]


def test_topology_noise_between(topology_dir):
    res = P.score_topology(_rows(topology_dir, M.model_noise(0.5), k=8), topology_dir)
    lo, hi = res["hole_recall"]["pooled"][1:]
    assert 0 < res["hole_recall"]["pooled"][0] < 1 and lo <= res["hole_recall"]["pooled"][0] <= hi
    assert res["hole_recall"]["any_of_k"][0] >= res["hole_recall"]["pooled"][0]


def test_stratified_bootstrap_equals_power_c2():
    ref = json.loads((ROOT / "tests" / "goldens" / "p2_lab_refs.json").read_text())["probes_pooled"]
    rng = np.random.default_rng(1)
    per_bin = [rng.random(120), rng.random(120) * 0.5, (rng.random(120) > 0.3).astype(float)]
    got = P.stratified_ci(per_bin, P.BIN_WEIGHTS, np.random.default_rng(5))
    assert np.allclose(got, ref, rtol=0, atol=1e-12)


@pytest.mark.parametrize("script", ["probes", "baselines"])
def test_cli_help(script):
    out = subprocess.run(
        [sys.executable, "-I", "-m", f"experiments.phase2.{script}", "--help"],
        cwd=ROOT,
        capture_output=True,
        text=True,
        timeout=120,
    )
    assert out.returncode == 0 and "usage" in out.stdout


def test_baselines_prior_on_fixture(tmp_path, disc_mask):
    from experiments.phase2 import baselines
    from objdel.verify.reward.cache import ExpertCache

    recs = []
    for i, (cx, r) in enumerate([(25, 8), (25, 10), (20, 6)]):
        iid = f"lidc/1.2.{i}/n000/z001"
        st = np.stack([disc_mask(50, 50, 25, cx, r), disc_mask(50, 50, 25, cx, r + 1)])
        c = ExpertCache.build(st)
        (tmp_path / "caches" / iid).mkdir(parents=True)
        recs.append(
            {
                "item_id": iid,
                "group_id": f"P{i}",
                "W": 50,
                "H": 50,
                "rl": True,
                "caches": {"nodule": c.save(tmp_path / "caches" / iid / "nodule.npz")},
                "targets": [M.mask_text(st[0]), M.mask_text(st[1])],
                "raters": ["1", "2"],
            }
        )
    (tmp_path / "lidc").mkdir()
    (tmp_path / "lidc" / "val.jsonl").write_text("".join(json.dumps(r) + "\n" for r in recs))
    prior = baselines.fit_prior(tmp_path, "lidc", "val")
    assert set(prior) == {50} and prior[50]["nodule"] == pytest.approx(
        np.median([2 * np.sqrt(P.gt_mask({"W": 50, "targets": [t]}).sum() / np.pi) for r in recs for t in r["targets"]])
    )
    rows = baselines.score_prior(tmp_path, "lidc", "val", prior)
    assert [r["item_id"] for r in rows] == [r["item_id"] for r in recs]
    assert rows[0]["iou_tc"]["nodule"] > 0.5 and rows[2]["iou_tc"]["nodule"] < rows[0]["iou_tc"]["nodule"]
    other = [{"item_id": r["item_id"], "iou_tc": {"nodule": 1.0}} for r in rows]
    w = baselines.win_share(other, rows, "nodule")
    assert w["share"] == 1.0 and w["n"] == 3 and 0 < w["sign_p"] <= 1


@pytest.mark.data
def test_real_v3e_fidelity_reproduces_b3(store_dir):
    stats = P.lidc_stats(store_dir)
    assert len(stats["holes"]) == 285 and len(stats["multis"]) == 182
    shares = np.bincount([P.bin_of(h["area"]) for h in stats["holes"]], minlength=3) / 285
    assert np.allclose(shares, P.BIN_WEIGHTS, atol=1e-4)
    for seed, ok in ((7, False), (8, True)):
        ref = json.loads((FIXTURES / f"topo_probe_v3e_seed{seed}.json").read_text())["fidelity"]
        got = P.fidelity(stats["holes"], stats["multis"], seed)
        assert got["ks_p"] == ref["ks_p"] and got["gen_bin_shares"] == ref["gen_bin_shares"]
        assert got["pass"] is ok
    # the recorded P5.2 residual: seed 7 fails on hole area (p 0.004), seed 8 passes


@pytest.mark.data
def test_real_built_probe_dirs_patient_disjoint():
    root = ROOT / "data" / "eval_views" / "probe-precision"
    if not (root / "lidc" / "val.jsonl").is_file():
        pytest.skip("probe-precision not built")
    from objdel.data import store

    val = {it.group_id for it in store.load("lidc", "val", ROOT / "data" / "store")}
    assert {r["group_id"] for r in P.records(root)} <= val  # probe-test: val patients only, never SFT-trained
    if (root / "lidc" / "train.jsonl").is_file():
        assert not val & {r["group_id"] for r in P.records(root, split="train")}
    topo = ROOT / "data" / "eval_views" / "probe-topology"
    if (topo / "lidc" / "val.jsonl").is_file():
        assert {r["group_id"] for r in P.records(topo)} <= val


def test_probe_dir_accepted_by_preload(precision_dir):
    serve = pytest.importorskip("objdel.verify.store")
    assert serve.preload(precision_dir, ["lidc"]) == 0  # single-rater items: no RL caches, records match manifest


def test_probe_outputs_through_phase2_entry(precision_dir, topology_dir):
    from objdel.verify import gates

    good = P.score_precision(_rows(precision_dir, M.model_noise(0.5)), precision_dir)
    disc = _rows(precision_dir, P.model_centroid_disc)
    topo = P.score_topology(_rows(topology_dir, M.model_a0), topology_dir)
    for prec, ok in (
        (good, True),
        (P.score_precision(disc, precision_dir), False),
        (P.score_precision(disc, precision_dir, capacity=_rows(precision_dir, M.model_exact)), True),
    ):
        out = gates.phase2_entry({"probes": {"precision": prec, "topology": topo}}, ("lidc",))
        assert out["gates"]["G-precision"]["pass"] is ok and out["gates"]["G-topology-a"]["pass"] is True
        assert out["reported"]["G-topology-b"] == topo["hole_recall"]
    assert any("capacity" in w for w in out["warnings"])

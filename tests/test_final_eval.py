"""experiments.phase2.final_eval: per-item K-sample means and the split summary (schema 2: mean, item, nodule and
patient bootstrap intervals, median), and --resummarize of a finished run."""

import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from experiments.phase2 import final_eval as FE
from objdel.eval import signal as S

ROOT = Path(__file__).resolve().parents[1]
CI_KEYS = {"mean", "ci95_nodule", "ci95_patient", "ci95_item", "higher_is_better"}


def result(item, legal, iou, f5, p, r, d, sample_iou=None):
    K = len(legal)
    return {"item_id": item, "samples": {"legal": legal, "iou_tc": iou, "f5": f5},
            "objects": {"nodule": {"precision_vote": p, "recall_vote": r, "dice_vote": d, "nsd": d, "boundary_iou": d,
                                   "boundary_iou_vote": d, "soft_dice": .5, "hm_iou": .5, "ged": .5, "vote_area": 10,
                                   "sample_iou": np.eye(K).tolist() if sample_iou is None else sample_iou}}}


def items_of(ids, values):
    return [{"item_id": i, **{m: float(v) for m in FE.METRICS}} for i, v in zip(ids, values, strict=True)]


def layout(patients, nodules, slices):
    """Ids lidc/p<p>/n<n>/z<s> in that order and their {item_id: patient} map."""
    ids = [f"lidc/p{p}/n{n:03d}/z{s}" for p in range(patients) for n in range(nodules) for s in range(slices)]
    return ids, {i: i.split("/")[1] for i in ids}


def walk(x):
    if isinstance(x, dict):
        for k, v in x.items():
            yield k, v
            yield from walk(v)
    elif isinstance(x, list):
        for v in x:
            yield from walk(v)


def test_per_item_selection_keys():
    """U07: K = 4; medoid scores .45, .5, .15 over the legal answers 0-2 -> answer 1; the illegal answer 3 is never
    the medoid although its row has IoU 1 with itself."""
    M = [[1, .8, .1, 0], [.8, 1, .2, 0], [.1, .2, 1, 0], [0, 0, 0, 1]]
    it = FE.per_item(result("lidc/a", [1, 1, 1, 0], [.9, .8, .1, 0], [0] * 4, [0] * 4, [0] * 4, [0] * 4, M))
    assert {k: it[k] for k in FE.SELECTION if k != "sel_vote"} == pytest.approx({
        "sel_mean_of_k": .45, "sel_legal_mean": 0.6, "sel_medoid": .8, "sel_medoid_legal": 1.0, "sel_oracle": .9},
        rel=0, abs=1e-15)
    assert (it["sel_medoid"], it["sel_medoid_legal"], it["sel_oracle"]) == (.8, 1.0, .9)
    assert "sel_vote" not in it  # absent, not None, without its input (9.8 rule 2)
    it = FE.per_item(result("lidc/a", [0] * 4, [0] * 4, [0] * 4, [0] * 4, [0] * 4, [0] * 4, M))
    assert (it["sel_legal_mean"], it["sel_medoid_legal"], it["sel_medoid"]) == (0.0, 0.0, 0.0)
    r = result("lidc/a", [1, 1, 1, 0], [.9, .8, .1, 0], [0] * 4, [0] * 4, [0] * 4, [0] * 4, M)
    r["objects"]["nodule"]["pixel_vote_iou_tc"] = 0.7
    assert FE.per_item(r)["sel_vote"] == 0.7
    assert all(type(v) is float for k, v in FE.per_item(r).items() if k not in ("item_id", "vote_area"))


def test_summary_schema_2():
    """U08 (A1 part): 20 items of 8 patients."""
    ids = [f"lidc/s{i % 8}/n{i % 3:03d}/z{i}" for i in range(20)]
    groups = {i: f"P{i.split('/')[1]}" for i in ids}
    v = np.random.default_rng(0).random(20)
    s = FE.summarize(items_of(ids, v), groups)
    assert (s["schema"], s["boot"], s["ci_level"], s["gate_ci"]) == (2, 2000, 0.95, "ci95_patient")
    assert (s["items"], s["nodules"], s["patients"]) == (20, len({S.cluster(i) for i in ids}), 8)
    for m, (_, higher) in FE.METRICS.items():
        assert set(s[m]) == CI_KEYS | ({"median"} if m == "iou_tc" else set()), m
        assert s[m]["higher_is_better"] is higher
    assert s["iou_tc"]["median"] == float(np.median(v))
    assert set(s["miss_rate"]) == CI_KEYS - {"higher_is_better"}
    assert "ci95" not in {k for k, _ in walk(s)}
    json.dumps(s, allow_nan=False)
    s = FE.summarize(items_of(ids, v))
    assert s["patients"] is None and all(s[m]["ci95_patient"] is None for m in [*FE.METRICS, "miss_rate"])


def test_summary_cis_recomputed_independently():
    """U09: 30 items, 10 nodules, 4 patients, 300 draws."""
    ids = [f"lidc/p{n % 4}/n{n:03d}/z{s}" for n in range(10) for s in range(3)]
    groups = {i: i.split("/")[1] for i in ids}
    v = np.random.default_rng(3).random(30)
    s = FE.summarize(items_of(ids, v), groups, draws=300, seed=0)["f5"]
    idx = np.random.default_rng(0).integers(0, 30, (300, 30))
    assert s["mean"] == v.mean()
    assert s["ci95_item"] == np.percentile(v[idx].mean(1), [2.5, 97.5]).tolist()
    for key, D in [("ci95_nodule", S.draws(ids, 300, 0)),
                   ("ci95_patient", S.draws(ids, 300, 0, key=groups.__getitem__))]:
        means = np.array([v[d].mean() for d in D])
        assert s[key] == [float(np.quantile(means, 0.025)), float(np.quantile(means, 0.975))], key
    # one item per cluster: the nodule and item CIs are equal, exactly, for every metric
    s = FE.summarize(items_of([f"riga/x/{i}" for i in range(60)], np.random.default_rng(4).random(60)))
    for m in FE.METRICS:
        assert s[m]["ci95_nodule"] == s[m]["ci95_item"], m
        assert s[m]["ci95_item"][0] <= s[m]["mean"] <= s[m]["ci95_item"][1], m
    s = FE.summarize(items_of(ids, [0.5] * 30), groups)
    assert all(s["iou_tc"][k] == [0.5, 0.5] for k in ("ci95_item", "ci95_nodule", "ci95_patient"))


def test_summary_of_one_item_has_a_point_interval():
    s = FE.summarize(items_of(["lidc/s/n000/z1"], [0.7]), {"lidc/s/n000/z1": "P"})
    assert s["iou_tc"]["ci95_item"] == s["iou_tc"]["ci95_nodule"] == s["iou_tc"]["ci95_patient"] == [0.7, 0.7]


def test_miss_rate_is_the_share_of_items_below_0_1():
    s = FE.summarize(items_of([f"x{i}" for i in range(4)], [0.05, 0.5, 0.0999999, 0.1]))
    assert s["miss_rate"]["mean"] == 0.5


def test_bad_input_raises_and_names_it():
    """U12: a NaN metric names the item and the metric; a missing group names the item; a metric that some items
    hold and others do not names the first item without it."""
    ids = [f"x{i}" for i in range(4)]
    items = items_of(ids, [0.2, 0.4, 0.6, 0.8])
    items[2]["f5"] = float("nan")
    with pytest.raises(ValueError, match="x2.*f5"):
        FE.summarize(items)
    with pytest.raises(KeyError, match="x3"):
        FE.summarize(items_of(ids, [0.2] * 4), {"x0": "a", "x1": "a", "x2": "b", "extra": "c"})
    items = items_of(ids, [0.2] * 4)
    del items[1]["dice_vote"]
    with pytest.raises(ValueError, match="x1.*dice_vote"):
        FE.summarize(items)


def test_summary_metrics_absent_from_every_item_are_left_out():
    items = [{k: v for k, v in it.items() if k != "dice_vote"} for it in items_of(["a", "b"], [0.3, 0.6])]
    assert "dice_vote" not in FE.summarize(items)


def test_summary_reruns_and_item_order():
    """I23: reruns are identical; shuffled items keep the means (1e-12) but change the CIs (draws follow id order)."""
    ids, groups = layout(6, 3, 4)
    v = np.random.default_rng(5).random(len(ids))
    items = items_of(ids, v)
    a, b = FE.summarize(items, groups), FE.summarize(items, groups)
    assert a == b
    perm = np.random.default_rng(1).permutation(len(ids))
    c = FE.summarize([items[j] for j in perm], groups)
    assert abs(c["iou_tc"]["mean"] - a["iou_tc"]["mean"]) < 1e-12
    assert c["iou_tc"]["ci95_nodule"] != a["iou_tc"]["ci95_nodule"]


def test_outputs_are_json_native():
    """I24 (A1 part): only float, int, bool, str, list, dict and None."""
    ids, groups = layout(3, 2, 2)
    s = FE.summarize(items_of(ids, np.linspace(0, 1, len(ids))), groups)
    it = FE.per_item(result("lidc/a", [1, 0], [0.8, 0.0], [0.9, 0.0], [0.6, 0.0], [1.0, 0.0], [0.75, 0.0]))
    for doc in (s, it):
        json.dumps(doc, allow_nan=False)
        assert all(type(v) in (float, int, bool, str, list, dict, type(None)) for _, v in walk(doc))


def test_group_id_stays_out_of_the_prompt_rows(tmp_path, synthetic_items):
    """U50: TRL forwards extra prompt columns to the reward functions, so group_id comes from the records."""
    from objdel.data import store
    from objdel.train import dataset as D

    store.write(synthetic_items(), tmp_path / "store")
    D.prepare(tmp_path / "store", tmp_path / "out", splits=("train",), procs=1)
    rows = D.prompts(tmp_path / "out", "lidc", "train", rl_only=False)
    assert rows and all("group_id" not in r for r in rows)
    recs = D.records(tmp_path / "out", "lidc", "train")
    assert FE.view_groups(tmp_path / "out", "train") == {r["item_id"]: r["group_id"] for r in recs}
    assert all(r["group_id"] for r in recs)


def finished_run(tmp_path, with_group=True):
    """A schema-1 run folder (20 items of 8 patients, six metrics, as R5's) and a fake LIDC store of its split."""
    import pyarrow as pa
    import pyarrow.parquet as pq

    ids = [f"lidc/s{i % 8}/n{i % 3:03d}/z{i}" for i in range(20)]
    run = tmp_path / "run"
    run.mkdir()
    v = np.random.default_rng(2).random(20)
    items = [{"item_id": i, "legal": 1.0, "iou_tc": float(x), "f5": float(x)} for i, x in zip(ids, v, strict=True)]
    (run / "items.jsonl").write_text("".join(json.dumps(it) + "\n" for it in items))
    (run / "final_eval.json").write_text(json.dumps({"ckpt": "/workspace/runs/r5/final_model", "view": "/w/v",
                                                     "split": "test", "K": 8, "seed": 11, "items": 20,
                                                     "iou_tc": {"mean": 0.5, "ci95": [0.4, 0.6]}}))
    store = tmp_path / "store"
    (store / "lidc").mkdir(parents=True)
    keep = ids if with_group else ids[1:]
    pq.write_table(pa.table({"item_id": keep, "group_id": [f"LIDC-IDRI-{i.split('/')[1][1:]:>04}" for i in keep]}),
                   store / "lidc" / "test.parquet")
    return run, store, items


def tree(path):
    return {p: p.read_bytes() for p in sorted(path.rglob("*")) if p.is_file()}


def test_resummarize_writes_schema_2_beside_the_run(tmp_path, monkeypatch):
    """U08b: no model is loaded, DIR is untouched, metrics_absent lists the METRICS keys the items lack."""
    from experiments.phase2 import g_probe

    monkeypatch.setattr(g_probe, "make_trainer", lambda *a, **k: (_ for _ in ()).throw(AssertionError("model")))
    run, store, items = finished_run(tmp_path)
    before = tree(run)
    FE.main(["--resummarize", str(run), "--store", str(store), "--out", str(tmp_path / "out")])
    assert tree(run) == before
    doc = json.loads((tmp_path / "out" / "final_eval.json").read_text())
    assert (doc["schema"], doc["gate_ci"], doc["patients"]) == (2, "ci95_patient", 8)
    assert doc["resummarized_from"] == str(run)
    assert doc["metrics_absent"] == TWELVE[3:]
    assert {k: doc[k] for k in ("ckpt", "view", "split", "K", "seed")} == {
        "ckpt": "/workspace/runs/r5/final_model", "view": "/w/v", "split": "test", "K": 8, "seed": 11}
    groups = {i: f"LIDC-IDRI-{i.split('/')[1][1:]:>04}" for i in (it["item_id"] for it in items)}
    assert doc["iou_tc"] == FE.summarize(items, groups)["iou_tc"]


def test_resummarize_exits_2_on_an_item_without_a_group(tmp_path, capsys):
    run, store, _ = finished_run(tmp_path, with_group=False)
    with pytest.raises(SystemExit) as e:
        FE.main(["--resummarize", str(run), "--store", str(store), "--out", str(tmp_path / "out")])
    assert e.value.code == 2 and "lidc/s0/n000/z0" in capsys.readouterr().err
    assert not (tmp_path / "out").exists()
    with pytest.raises(SystemExit) as e:
        FE.main(["--resummarize", str(run), "--store", str(store), "--out", str(run)])
    assert e.value.code == 2


def test_resummarize_is_repeatable_across_processes(tmp_path):
    """E02 (A1 part): fresh processes with other hash seeds and thread counts write the same bytes."""
    run, store, _ = finished_run(tmp_path)
    outs = []
    for seed, threads in (("1", "1"), ("2", "4")):
        out = tmp_path / f"out{seed}"
        env = os.environ | {"PYTHONHASHSEED": seed, "OMP_NUM_THREADS": threads}
        subprocess.run([sys.executable, "-I", "-m", "experiments.phase2.final_eval", "--resummarize", str(run),
                        "--store", str(store), "--out", str(out)], cwd=ROOT, env=env, check=True)
        outs.append((out / "final_eval.json").read_bytes())
    assert outs[0] == outs[1]


def test_cli_needs_a_checkpoint_view_and_training_data(capsys):
    with pytest.raises(SystemExit):
        FE.main(["--view", "v", "--train", "t", "--out", "o"])
    assert "--ckpt" in capsys.readouterr().err


def sel_items(ids, mean, medoid, vote=None):
    return [{"item_id": i, "iou_tc": a, "sel_mean_of_k": a, "sel_medoid": b,
             **({} if vote is None else {"sel_vote": vote[j]})}
            for j, (i, a, b) in enumerate(zip(ids, mean, medoid, strict=True))]


def test_paired_differences():
    """U11: mean and both CIs of the per-item differences under the summary's draws; the unit is not swapped."""
    s = FE.summarize(sel_items([f"x{i}" for i in range(4)], [.5, .5, .5, .5], [.6, .6, .45, .5]))
    p = s["paired"]["medoid_minus_mean"]
    assert set(p) == {"mean", "ci95_patient", "ci95_nodule"} and abs(p["mean"] - 0.0375) < 1e-12
    assert s["paired"]["vote_minus_mean"] is None
    ids, groups = layout(3, 2, 2)
    diff = np.array([.1, .1, 0, 0, .05, .05, -.05, -.05, .2, .2, .1, .1])
    s = FE.summarize(sel_items(ids, [.5] * 12, .5 + diff), groups)["paired"]["medoid_minus_mean"]
    want = {k: S.ci(.5 + diff - .5, D)[1:] for k, D in (
        ("ci95_patient", S.draws(ids, 2000, 0, key=groups.__getitem__)), ("ci95_nodule", S.draws(ids, 2000, 0)))}
    for k, w in want.items():
        assert np.allclose(s[k], w, rtol=0, atol=1e-12) and s[k][0] <= s["mean"] <= s[k][1], k
    assert abs(s["ci95_nodule"][1] - 0.1333) < 1e-4 and abs(s["ci95_patient"][1] - 0.15) < 1e-12
    v = np.random.default_rng(0).random(12)
    s = FE.summarize(sel_items(ids, v, v + 0.05, v), groups)["paired"]
    assert abs(s["medoid_minus_mean"]["mean"] - 0.05) < 1e-12
    assert all(abs(x - 0.05) < 1e-12 for k in ("ci95_patient", "ci95_nodule") for x in s["medoid_minus_mean"][k])
    assert set(s["vote_minus_mean"]) == {"mean", "ci95_patient", "ci95_nodule"}
    items = sel_items(ids, v, v, v)
    del items[3]["sel_vote"]
    with pytest.raises(ValueError, match=ids[3]):
        FE.summarize(items, groups)


def test_paired_ci_is_narrower_than_the_unpaired_one_on_correlated_arms():
    ids, groups = layout(6, 3, 2)
    base = np.random.default_rng(1).random(len(ids))
    s = FE.summarize(sel_items(ids, base, base + 0.02), groups)
    paired = np.diff(s["paired"]["medoid_minus_mean"]["ci95_patient"])[0]
    unpaired = np.diff(s["sel_medoid"]["ci95_patient"])[0] + np.diff(s["sel_mean_of_k"]["ci95_patient"])[0]
    assert paired < unpaired


def summary(mean, patient, nodule=(0.0, 0.1), vote=None):
    pair = {"mean": mean, "ci95_patient": None if patient is None else list(patient), "ci95_nodule": list(nodule)}
    return {"gate_ci": "ci95_patient", "paired": {"medoid_minus_mean": pair, "vote_minus_mean": vote}}


def test_gate_a_rule():
    """U51: the val patient CI decides; the nodule CI and the test summary never do."""
    assert FE.gate_a(summary(0.015, (0.002, 0.03)), None)["change_default"] is True
    assert FE.gate_a(summary(0.015, (0.002, 0.03)), None)["method"] == "medoid"
    assert FE.gate_a(summary(0.0149, (0.002, 0.03)), None)["change_default"] is False
    assert FE.gate_a(summary(0.02, (0.0, 0.04)), None)["change_default"] is False
    g = FE.gate_a(summary(0.02, (-0.003, 0.043), nodule=(0.004, 0.036)), None)
    assert (g["method"], g["change_default"]) == ("mean_of_k", False)
    with pytest.raises(ValueError, match="ci95_patient"):
        FE.gate_a(summary(0.02, None, nodule=(0.004, 0.036)), None)
    vote = {"mean": 0.05, "ci95_patient": [0.03, 0.07], "ci95_nodule": [0.03, 0.07]}
    g = FE.gate_a(summary(0.0, (-0.01, 0.01), vote=vote), None)
    assert g["change_default"] is False and g["method"] == "mean_of_k" and "Phase D" in g["reason"]
    val, test = summary(0.02, (0.005, 0.035)), summary(-0.01, (-0.02, 0.0))
    a, b = FE.gate_a(val, None), FE.gate_a(val, test)
    assert a["method"] == b["method"] == "medoid" and a["test_same_sign"] is None and b["test_same_sign"] is False
    assert FE.gate_a(val, summary(0.01, (0.0, 0.02)))["test_same_sign"] is True


def test_selection_key_invariants():
    """P08: 200 random items, K = 8."""
    rng = np.random.default_rng(8)
    for _ in range(200):
        legal = (rng.random(8) < 0.7).astype(float)
        iou = (rng.random(8) * legal).tolist()
        a = rng.random((8, 8))
        M = (np.triu(a) + np.triu(a, 1).T).tolist()
        it = FE.per_item(result("lidc/a", legal.tolist(), iou, iou, iou, iou, iou, M))
        assert it["sel_oracle"] >= it["sel_medoid"] >= 0 and it["sel_oracle"] >= it["sel_legal_mean"]
        assert it["sel_oracle"] >= it["sel_mean_of_k"] and it["sel_legal_mean"] >= it["sel_mean_of_k"]
        assert it["sel_medoid_legal"] in (0.0, 1.0) and all(0 <= it[k] <= 1 for k in FE.SELECTION if k in it)


TWELVE = ["legal", "iou_tc", "f5", "precision_vote", "recall_vote", "dice_vote", "nsd", "boundary_iou",
          "boundary_iou_vote", "soft_dice", "hm_iou", "ged"]


def test_metrics_table():
    """U05."""
    assert list(FE.METRICS) == TWELVE and not hasattr(FE, "TERMS")
    assert {m for m, (_, h) in FE.METRICS.items() if not h} == {"ged"}
    assert {s for s, _ in FE.METRICS.values()} == {"samples", "objects", "item"}


def k4_result(legal=(1, 1, 0, 1), iou=(.8, .6, 0, .4), obj="nodule"):
    K = len(legal)
    o = {"precision_vote": [.5, .7, 0, .8], "recall_vote": [.9, .5, 0, .6], "dice_vote": [.6, .6, 0, .6],
         "nsd": [.3, .2, 0, .3], "boundary_iou": [.5, .5, 0, .5], "boundary_iou_vote": [.2, .2, 0, .2],
         "soft_dice": .33, "hm_iou": .44, "ged": .55, "vote_area": 40, "sample_iou": np.eye(K).tolist()}
    o = {k: v[:K] if isinstance(v, list) and k != "sample_iou" else v for k, v in o.items()}
    return {"item_id": "lidc/a", "samples": {"legal": list(legal), "iou_tc": list(iou), "f5": [.9, .7, 0, .5][:K]},
            "objects": {obj: o}}


def test_per_item_sources_averages_and_miss():
    """U06."""
    it = FE.per_item(k4_result())
    want = {"legal": .75, "iou_tc": .45, "f5": .525, "precision_vote": .5, "recall_vote": .5, "dice_vote": .45,
            "nsd": .2, "boundary_iou": .375, "boundary_iou_vote": .15, "soft_dice": .33, "hm_iou": .44, "ged": .55}
    assert {m: it[m] for m in FE.METRICS} == pytest.approx(want, rel=0, abs=1e-15)
    assert (it["soft_dice"], it["hm_iou"], it["ged"], it["vote_area"]) == (.33, .44, .55, 40)
    assert all(type(it[m]) is float for m in FE.METRICS) and it["miss"] == 0.0
    for x, miss in ((0.0999999, 1.0), (0.1, 0.0), (0.1000001, 0.0)):
        assert FE.per_item(k4_result(iou=(x,) * 4))["miss"] == miss
    # literal pin: numpy's mean of eight 0.1 is 0.1 exactly (pairwise sum), so K = 8 answers at 0.1 are found
    assert np.mean([0.1] * 8) == 0.1 and FE.per_item(k4_result(iou=(0.1,) * 8, legal=(1,) * 8))["miss"] == 0.0
    assert FE.per_item(k4_result(obj="cup"), "cup")["ged"] == .55
    with pytest.raises(KeyError):
        FE.per_item(k4_result(obj="cup"))
    for K in (1, 8):
        it = FE.per_item(k4_result(legal=(0,) * K, iou=(0,) * K))
        assert (it["legal"], it["sel_legal_mean"], it["miss"]) == (0.0, 0.0, 1.0)


def test_miss_rate_and_iou_found():
    """U10: iou_found is the mean over found items, its CIs drawn over the found subset's ids."""
    s = FE.summarize(items_of([f"x{i}" for i in range(4)], [0.05, 0.5, 0.7, 0.0]))
    assert s["miss_rate"]["mean"] == .5 and abs(s["iou_found"]["mean"] - .6) < 1e-12
    ids = [f"lidc/p{p}/n{n}/s{z}" for p in range(3) for n in range(2) for z in range(2)]
    groups = {i: i.split("/")[1] for i in ids}
    v = np.array([.05, .5, .7, 0, .3, .6, .8, .4, .9, .2, .65, .55])
    f = FE.summarize(items_of(ids, v), groups)["iou_found"]
    found = [i for i, x in zip(ids, v, strict=True) if x >= 0.1]
    x = v[v >= 0.1]
    want_p = S.ci(x, S.draws(found, 2000, 0, key=groups.__getitem__))[1:]
    want_n = S.ci(x, S.draws(found, 2000, 0))[1:]
    assert np.allclose(f["ci95_patient"], want_p, rtol=0, atol=1e-12)
    assert np.allclose(f["ci95_nodule"], want_n, rtol=0, atol=1e-12) and f["ci95_patient"] != f["ci95_nodule"]
    assert np.allclose(f["ci95_nodule"], [0.5, 0.62], atol=1e-12) and np.allclose(f["ci95_patient"], [0.525, 0.6])
    assert FE.summarize(items_of(ids, v + 0.1))["miss_rate"]["mean"] == 0.0
    s = FE.summarize(items_of(ids, v * 0), groups)
    assert s["iou_found"] is None
    json.dumps(s, allow_nan=False)


def test_summary_has_every_metric_and_selection_key():
    """U08 (A4 / A5 part): every METRICS and selection key the items hold, the paired block, iou_found."""
    ids = [f"lidc/s{i % 8}/n{i % 3:03d}/z{i}" for i in range(20)]
    groups = {i: f"P{i.split('/')[1]}" for i in ids}
    rng = np.random.default_rng(0)
    items = [FE.per_item(k4_result(iou=tuple(rng.random(4)))) | {"item_id": i} for i in ids]
    s = FE.summarize(items, groups)
    for m in [*FE.METRICS, *FE.SELECTION[:-1]]:
        assert set(s[m]) == CI_KEYS | ({"median"} if m == "iou_tc" else set()), m
        assert s[m]["higher_is_better"] is (m != "ged"), m
    assert {"miss_rate", "iou_found", "paired"} <= set(s) and "sel_vote" not in s
    json.dumps(s, allow_nan=False)


def test_size_terciles():
    """U52: terciles of the vote-mask area with their means; ties stay in one tercile; fewer than 3 items raise."""
    items = [{"item_id": f"x{a}", "vote_area": a, "iou_tc": a / 100, "boundary_iou": a / 200}
             for a in range(10, 100, 10)]
    t = FE.size_terciles(items)
    assert t["bounds"] == [30, 60] and [g["items"] for g in t["terciles"]] == [3, 3, 3]
    assert [g["iou_tc"] for g in t["terciles"]] == pytest.approx([.2, .5, .8])
    assert [g["boundary_iou"] for g in t["terciles"]] == pytest.approx([.1, .25, .4])
    tied = [it | {"vote_area": 30} if it["vote_area"] == 40 else it for it in items]
    assert [g["items"] for g in FE.size_terciles(tied)["terciles"]] == [4, 2, 3]
    with pytest.raises(ValueError, match="3 items"):
        FE.size_terciles(items[:2])

"""experiments.phase2.final_eval: per-item K-sample means and the split summary (schema 2: mean, item, nodule and
patient bootstrap intervals, median), and --resummarize of a finished run."""

import hashlib
import json
import os
import re
import subprocess
import sys
import time
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from experiments.phase2 import final_eval as CLI
from objdel.policy import evaluate as EV
from objdel.policy import final_eval as FE
from objdel.verify import stats as S

ROOT = Path(__file__).resolve().parents[2]
CI_KEYS = {"mean", "ci95_nodule", "ci95_patient", "ci95_item", "higher_is_better"}


def result(item, legal, iou, f5, p, r, d, sample_iou=None):
    K = len(legal)
    return {"item_id": item, "samples": {"legal": legal, "iou_tc": iou, "f5": f5},
            "objects": {"nodule": {"precision_vote": p, "recall_vote": r, "dice_vote": d, "nsd": d, "boundary_iou": d,
                                   "boundary_iou_vote": d, "soft_dice": .5, "hm_iou": .5, "ged": .5, "vote_area": 10,
                                   "sample_iou": np.eye(K).tolist() if sample_iou is None else sample_iou}}}


def items_of(ids, values):
    return [{"item_id": i, **{m: float(v) for m in EV.METRICS}} for i, v in zip(ids, values, strict=True)]


def layout(patients, nodules, slices):
    """Ids lidc/p<p>/n<n>/z<s> in that order and their {item_id: patient} map."""
    ids = [f"lidc/p{p}/n{n:03d}/z{s}" for p in range(patients) for n in range(nodules) for s in range(slices)]
    return ids, {i: i.split("/")[1] for i in ids}


def git(repo, *args):
    who = {"GIT_AUTHOR_NAME": "t", "GIT_AUTHOR_EMAIL": "t@t", "GIT_COMMITTER_NAME": "t", "GIT_COMMITTER_EMAIL": "t@t"}
    return subprocess.run(["git", *args], cwd=repo, capture_output=True, text=True, check=True,
                          env=os.environ | who).stdout


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
    it = EV.per_item(result("lidc/a", [1, 1, 1, 0], [.9, .8, .1, 0], [0] * 4, [0] * 4, [0] * 4, [0] * 4, M))
    assert {k: it[k] for k in EV.SELECTION if k != "sel_vote"} == pytest.approx({
        "sel_mean_of_k": .45, "sel_legal_mean": 0.6, "sel_medoid": .8, "sel_medoid_legal": 1.0, "sel_oracle": .9},
        rel=0, abs=1e-15)
    assert (it["sel_medoid"], it["sel_medoid_legal"], it["sel_oracle"]) == (.8, 1.0, .9)
    assert "sel_vote" not in it  # absent, not None, without its input (9.8 rule 2)
    it = EV.per_item(result("lidc/a", [0] * 4, [0] * 4, [0] * 4, [0] * 4, [0] * 4, [0] * 4, M))
    assert (it["sel_legal_mean"], it["sel_medoid_legal"], it["sel_medoid"]) == (0.0, 0.0, 0.0)
    r = result("lidc/a", [1, 1, 1, 0], [.9, .8, .1, 0], [0] * 4, [0] * 4, [0] * 4, [0] * 4, M)
    r["objects"]["nodule"]["pixel_vote_iou_tc"] = 0.7
    assert EV.per_item(r)["sel_vote"] == 0.7
    assert all(type(v) is float for k, v in EV.per_item(r).items() if k not in ("item_id", "vote_area"))


def test_summary_schema_2():
    """U08 (A1 part): 20 items of 8 patients."""
    ids = [f"lidc/s{i % 8}/n{i % 3:03d}/z{i}" for i in range(20)]
    groups = {i: f"P{i.split('/')[1]}" for i in ids}
    v = np.random.default_rng(0).random(20)
    s = EV.summarize(items_of(ids, v), groups)
    assert (s["schema"], s["boot"], s["ci_level"], s["gate_ci"]) == (2, 2000, 0.95, "ci95_patient")
    assert (s["items"], s["nodules"], s["patients"]) == (20, len({S.cluster(i) for i in ids}), 8)
    for m, (_, higher) in EV.METRICS.items():
        assert set(s[m]) == CI_KEYS | ({"median"} if m == "iou_tc" else set()), m
        assert s[m]["higher_is_better"] is higher
    assert s["iou_tc"]["median"] == float(np.median(v))
    assert set(s["miss_rate"]) == CI_KEYS - {"higher_is_better"}
    assert "ci95" not in {k for k, _ in walk(s)}
    json.dumps(s, allow_nan=False)
    s = EV.summarize(items_of(ids, v))
    assert s["patients"] is None and all(s[m]["ci95_patient"] is None for m in [*EV.METRICS, "miss_rate"])


def test_summary_cis_recomputed_independently():
    """U09: 30 items, 10 nodules, 4 patients, 300 draws."""
    ids = [f"lidc/p{n % 4}/n{n:03d}/z{s}" for n in range(10) for s in range(3)]
    groups = {i: i.split("/")[1] for i in ids}
    v = np.random.default_rng(3).random(30)
    s = EV.summarize(items_of(ids, v), groups, draws=300, seed=0)["f5"]
    idx = np.random.default_rng(0).integers(0, 30, (300, 30))
    assert s["mean"] == v.mean()
    assert s["ci95_item"] == np.percentile(v[idx].mean(1), [2.5, 97.5]).tolist()
    for key, D in [("ci95_nodule", S.draws(ids, 300, 0)),
                   ("ci95_patient", S.draws(ids, 300, 0, key=groups.__getitem__))]:
        means = np.array([v[d].mean() for d in D])
        assert s[key] == [float(np.quantile(means, 0.025)), float(np.quantile(means, 0.975))], key
    # one item per cluster: the nodule and item CIs are equal, exactly, for every metric
    s = EV.summarize(items_of([f"riga/x/{i}" for i in range(60)], np.random.default_rng(4).random(60)))
    for m in EV.METRICS:
        assert s[m]["ci95_nodule"] == s[m]["ci95_item"], m
        assert s[m]["ci95_item"][0] <= s[m]["mean"] <= s[m]["ci95_item"][1], m
    s = EV.summarize(items_of(ids, [0.5] * 30), groups)
    assert all(s["iou_tc"][k] == [0.5, 0.5] for k in ("ci95_item", "ci95_nodule", "ci95_patient"))


def test_summary_of_one_item_has_a_point_interval():
    s = EV.summarize(items_of(["lidc/s/n000/z1"], [0.7]), {"lidc/s/n000/z1": "P"})
    assert s["iou_tc"]["ci95_item"] == s["iou_tc"]["ci95_nodule"] == s["iou_tc"]["ci95_patient"] == [0.7, 0.7]


def test_miss_rate_is_the_share_of_items_below_0_1():
    s = EV.summarize(items_of([f"x{i}" for i in range(4)], [0.05, 0.5, 0.0999999, 0.1]))
    assert s["miss_rate"]["mean"] == 0.5


def test_bad_input_raises_and_names_it():
    """U12: a NaN metric names the item and the metric; a missing group names the item; a metric that some items
    hold and others do not names the first item without it."""
    ids = [f"x{i}" for i in range(4)]
    items = items_of(ids, [0.2, 0.4, 0.6, 0.8])
    items[2]["f5"] = float("nan")
    with pytest.raises(ValueError, match="x2.*f5"):
        EV.summarize(items)
    with pytest.raises(KeyError, match="x3"):
        EV.summarize(items_of(ids, [0.2] * 4), {"x0": "a", "x1": "a", "x2": "b", "extra": "c"})
    items = items_of(ids, [0.2] * 4)
    del items[1]["dice_vote"]
    with pytest.raises(ValueError, match="x1.*dice_vote"):
        EV.summarize(items)


def test_summary_metrics_absent_from_every_item_are_left_out():
    items = [{k: v for k, v in it.items() if k != "dice_vote"} for it in items_of(["a", "b"], [0.3, 0.6])]
    assert "dice_vote" not in EV.summarize(items)


def test_summary_reruns_and_item_order():
    """I23: reruns are identical; shuffled items keep the means (1e-12) but change the CIs (draws follow id order)."""
    ids, groups = layout(6, 3, 4)
    v = np.random.default_rng(5).random(len(ids))
    items = items_of(ids, v)
    a, b = EV.summarize(items, groups), EV.summarize(items, groups)
    assert a == b
    perm = np.random.default_rng(1).permutation(len(ids))
    c = EV.summarize([items[j] for j in perm], groups)
    assert abs(c["iou_tc"]["mean"] - a["iou_tc"]["mean"]) < 1e-12
    assert c["iou_tc"]["ci95_nodule"] != a["iou_tc"]["ci95_nodule"]


def test_outputs_are_json_native():
    """I24 (A1 part): only float, int, bool, str, list, dict and None."""
    ids, groups = layout(3, 2, 2)
    s = EV.summarize(items_of(ids, np.linspace(0, 1, len(ids))), groups)
    it = EV.per_item(result("lidc/a", [1, 0], [0.8, 0.0], [0.9, 0.0], [0.6, 0.0], [1.0, 0.0], [0.75, 0.0]))
    for doc in (s, it):
        json.dumps(doc, allow_nan=False)
        assert all(type(v) in (float, int, bool, str, list, dict, type(None)) for _, v in walk(doc))


def test_group_id_stays_out_of_the_prompt_rows(tmp_path, synthetic_items):
    """U50: TRL forwards extra prompt columns to the reward functions, so group_id comes from the records."""
    from objdel.data import prepare as data_prepare
    from objdel.data import store
    from objdel.policy import prompt as policy_prompt
    from objdel.policy import rows as policy_rows

    store.write(synthetic_items(), tmp_path / "store")
    data_prepare.prepare(tmp_path / "store", tmp_path / "out", splits=("train",), procs=1,
        prompts=policy_prompt.PROMPTS)
    rows = policy_rows.prompts(tmp_path / "out", "lidc", "train", rl_only=False)
    assert rows and all("group_id" not in r for r in rows)
    recs = policy_rows.records(tmp_path / "out", "lidc", "train")
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
    from objdel.policy import sampler as sample

    monkeypatch.setattr(sample, "VllmSampler", lambda *a, **k: (_ for _ in ()).throw(AssertionError("model")))
    run, store, items = finished_run(tmp_path)
    before = tree(run)
    CLI.main(["--resummarize", str(run), "--store", str(store), "--out", str(tmp_path / "out")])
    assert tree(run) == before
    doc = json.loads((tmp_path / "out" / "final_eval.json").read_text())
    assert (doc["schema"], doc["gate_ci"], doc["patients"]) == (2, "ci95_patient", 8)
    assert doc["resummarized_from"] == str(run)
    assert doc["metrics_absent"] == TWELVE[3:]
    assert {k: doc[k] for k in ("ckpt", "view", "split", "K", "seed")} == {
        "ckpt": "/workspace/runs/r5/final_model", "view": "/w/v", "split": "test", "K": 8, "seed": 11}
    groups = {i: f"LIDC-IDRI-{i.split('/')[1][1:]:>04}" for i in (it["item_id"] for it in items)}
    assert doc["iou_tc"] == EV.summarize(items, groups)["iou_tc"]


def test_resummarize_exits_2_on_an_item_without_a_group(tmp_path, capsys):
    run, store, _ = finished_run(tmp_path, with_group=False)
    with pytest.raises(SystemExit) as e:
        CLI.main(["--resummarize", str(run), "--store", str(store), "--out", str(tmp_path / "out")])
    assert e.value.code == 2 and "lidc/s0/n000/z0" in capsys.readouterr().err
    assert not (tmp_path / "out").exists()
    with pytest.raises(SystemExit) as e:
        CLI.main(["--resummarize", str(run), "--store", str(store), "--out", str(run)])
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


def sel_items(ids, mean, medoid, vote=None):
    return [{"item_id": i, "iou_tc": a, "sel_mean_of_k": a, "sel_medoid": b,
             **({} if vote is None else {"sel_vote": vote[j]})}
            for j, (i, a, b) in enumerate(zip(ids, mean, medoid, strict=True))]


def test_paired_differences():
    """U11: mean and both CIs of the per-item differences under the summary's draws; the unit is not swapped."""
    s = EV.summarize(sel_items([f"x{i}" for i in range(4)], [.5, .5, .5, .5], [.6, .6, .45, .5]))
    p = s["paired"]["medoid_minus_mean"]
    assert set(p) == {"mean", "ci95_patient", "ci95_nodule"} and abs(p["mean"] - 0.0375) < 1e-12
    assert s["paired"]["vote_minus_mean"] is None
    ids, groups = layout(3, 2, 2)
    diff = np.array([.1, .1, 0, 0, .05, .05, -.05, -.05, .2, .2, .1, .1])
    s = EV.summarize(sel_items(ids, [.5] * 12, .5 + diff), groups)["paired"]["medoid_minus_mean"]
    want = {k: S.ci(.5 + diff - .5, D)[1:] for k, D in (
        ("ci95_patient", S.draws(ids, 2000, 0, key=groups.__getitem__)), ("ci95_nodule", S.draws(ids, 2000, 0)))}
    for k, w in want.items():
        assert np.allclose(s[k], w, rtol=0, atol=1e-12) and s[k][0] <= s["mean"] <= s[k][1], k
    assert abs(s["ci95_nodule"][1] - 0.1333) < 1e-4 and abs(s["ci95_patient"][1] - 0.15) < 1e-12
    v = np.random.default_rng(0).random(12)
    s = EV.summarize(sel_items(ids, v, v + 0.05, v), groups)["paired"]
    assert abs(s["medoid_minus_mean"]["mean"] - 0.05) < 1e-12
    assert all(abs(x - 0.05) < 1e-12 for k in ("ci95_patient", "ci95_nodule") for x in s["medoid_minus_mean"][k])
    assert set(s["vote_minus_mean"]) == {"mean", "ci95_patient", "ci95_nodule"}
    items = sel_items(ids, v, v, v)
    del items[3]["sel_vote"]
    with pytest.raises(ValueError, match=ids[3]):
        EV.summarize(items, groups)


def test_paired_ci_is_narrower_than_the_unpaired_one_on_correlated_arms():
    ids, groups = layout(6, 3, 2)
    base = np.random.default_rng(1).random(len(ids))
    s = EV.summarize(sel_items(ids, base, base + 0.02), groups)
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
        it = EV.per_item(result("lidc/a", legal.tolist(), iou, iou, iou, iou, iou, M))
        assert it["sel_oracle"] >= it["sel_medoid"] >= 0 and it["sel_oracle"] >= it["sel_legal_mean"]
        assert it["sel_oracle"] >= it["sel_mean_of_k"] and it["sel_legal_mean"] >= it["sel_mean_of_k"]
        assert it["sel_medoid_legal"] in (0.0, 1.0) and all(0 <= it[k] <= 1 for k in EV.SELECTION if k in it)


TWELVE = ["legal", "iou_tc", "f5", "precision_vote", "recall_vote", "dice_vote", "nsd", "boundary_iou",
          "boundary_iou_vote", "soft_dice", "hm_iou", "ged"]


def test_metrics_table():
    """U05."""
    assert list(EV.METRICS) == TWELVE and not hasattr(EV, "TERMS")
    assert {m for m, (_, h) in EV.METRICS.items() if not h} == {"ged"}
    assert {s for s, _ in EV.METRICS.values()} == {"samples", "objects", "item"}


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
    it = EV.per_item(k4_result())
    want = {"legal": .75, "iou_tc": .45, "f5": .525, "precision_vote": .5, "recall_vote": .5, "dice_vote": .45,
            "nsd": .2, "boundary_iou": .375, "boundary_iou_vote": .15, "soft_dice": .33, "hm_iou": .44, "ged": .55}
    assert {m: it[m] for m in EV.METRICS} == pytest.approx(want, rel=0, abs=1e-15)
    assert (it["soft_dice"], it["hm_iou"], it["ged"], it["vote_area"]) == (.33, .44, .55, 40)
    assert all(type(it[m]) is float for m in EV.METRICS) and it["miss"] == 0.0
    for x, miss in ((0.0999999, 1.0), (0.1, 0.0), (0.1000001, 0.0)):
        assert EV.per_item(k4_result(iou=(x,) * 4))["miss"] == miss
    # literal pin: numpy's mean of eight 0.1 is 0.1 exactly (pairwise sum), so K = 8 answers at 0.1 are found
    assert np.mean([0.1] * 8) == 0.1 and EV.per_item(k4_result(iou=(0.1,) * 8, legal=(1,) * 8))["miss"] == 0.0
    assert EV.per_item(k4_result(obj="cup"), "cup")["ged"] == .55
    with pytest.raises(KeyError):
        EV.per_item(k4_result(obj="cup"))
    for K in (1, 8):
        it = EV.per_item(k4_result(legal=(0,) * K, iou=(0,) * K))
        assert (it["legal"], it["sel_legal_mean"], it["miss"]) == (0.0, 0.0, 1.0)


def test_miss_rate_and_iou_found():
    """U10: iou_found is the mean over found items, its CIs drawn over the found subset's ids."""
    s = EV.summarize(items_of([f"x{i}" for i in range(4)], [0.05, 0.5, 0.7, 0.0]))
    assert s["miss_rate"]["mean"] == .5 and abs(s["iou_found"]["mean"] - .6) < 1e-12
    ids = [f"lidc/p{p}/n{n}/s{z}" for p in range(3) for n in range(2) for z in range(2)]
    groups = {i: i.split("/")[1] for i in ids}
    v = np.array([.05, .5, .7, 0, .3, .6, .8, .4, .9, .2, .65, .55])
    f = EV.summarize(items_of(ids, v), groups)["iou_found"]
    found = [i for i, x in zip(ids, v, strict=True) if x >= 0.1]
    x = v[v >= 0.1]
    want_p = S.ci(x, S.draws(found, 2000, 0, key=groups.__getitem__))[1:]
    want_n = S.ci(x, S.draws(found, 2000, 0))[1:]
    assert np.allclose(f["ci95_patient"], want_p, rtol=0, atol=1e-12)
    assert np.allclose(f["ci95_nodule"], want_n, rtol=0, atol=1e-12) and f["ci95_patient"] != f["ci95_nodule"]
    assert np.allclose(f["ci95_nodule"], [0.5, 0.62], atol=1e-12) and np.allclose(f["ci95_patient"], [0.525, 0.6])
    assert EV.summarize(items_of(ids, v + 0.1))["miss_rate"]["mean"] == 0.0
    s = EV.summarize(items_of(ids, v * 0), groups)
    assert s["iou_found"] is None
    json.dumps(s, allow_nan=False)


def test_summary_has_every_metric_and_selection_key():
    """U08 (A4 / A5 part): every METRICS and selection key the items hold, the paired block, iou_found."""
    ids = [f"lidc/s{i % 8}/n{i % 3:03d}/z{i}" for i in range(20)]
    groups = {i: f"P{i.split('/')[1]}" for i in ids}
    rng = np.random.default_rng(0)
    items = [EV.per_item(k4_result(iou=tuple(rng.random(4)))) | {"item_id": i} for i in ids]
    s = EV.summarize(items, groups)
    for m in [*EV.METRICS, *EV.SELECTION[:-1]]:
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


# A3: saved answers, write order, rescore, code identity, view hashes, noise floor, view identity ----------------------

BOWTIE = '{"nodule":[{"outer":[[100,100],[900,900],[900,100],[100,900]],"holes":[]}]}'  # self-intersecting: illegal


@pytest.fixture(scope="module")
def ev_root(tmp_path_factory, synthetic_items):
    """A prepared val root of 6 LIDC items in 3 patients (test_ceiling's CLI view)."""
    from tests.experiments.test_ceiling import _view

    return _view(tmp_path_factory.mktemp("fe"), synthetic_items)


@pytest.fixture
def ev(ev_root, clean_store, tmp_path):
    """A fresh copy of the view per test (tests change its manifest)."""
    import shutil

    return Path(shutil.copytree(ev_root, tmp_path / "lidc-val-fix"))


def answer_texts(rec, K=8):
    """K texts: the record's targets (legal), one self-intersecting, one unparseable, one empty list; answer 7 is a
    legal target flagged truncated (fake_sampler)."""
    good = [t for t in rec["targets"] if t]
    return ((good * K)[:K - 4] + [BOWTIE, "garbage", "[]", good[0]])[:K] if K >= 4 else (good * K)[:K]


def fake_sampler(view, K=8, calls=None, loads=None):
    """A load_sampler stand-in: a Sampler answering answer_texts (answer 7 flagged truncated), tokens len + 1."""
    from objdel.policy import rows as policy_rows
    from objdel.policy.sampler import Sample

    recs = {r["item_id"]: r for r in policy_rows.records(view, "lidc", "val")}

    class Fake:
        def sample(self, rows, params):
            if calls is not None:
                calls.append([r["key"] for r in rows])
            k, out = params["n"], []
            for r in rows:
                texts = answer_texts(recs[r["key"]], k)
                out.append(Sample(r["key"], tuple(texts), tuple(tuple(range(len(t) + 1)) for t in texts),
                                  tuple(i == 7 for i in range(k))))
            return out

    def load_sampler(a):
        if loads is not None:
            loads.append(1)
        return Fake(), {"n": a.k}

    return load_sampler


def run_fe(ev, out, monkeypatch, *extra, scorer=None, K=8, loads=None):
    monkeypatch.setattr(FE, "load_sampler", fake_sampler(ev, K, loads=loads))
    if scorer:
        monkeypatch.setattr(FE, "score_items", scorer)  # --rescore
        monkeypatch.setattr(EV, "score_items", scorer)  # evaluate
    return CLI.main(["--ckpt", "c", "--view", str(ev), "--split", "val", "--out", str(out),
                    "-k", str(K), *extra])


def jsonl(path):
    return [json.loads(x) for x in Path(path).read_text().splitlines()]


RUN_SHA_KEYS = {"ckpt", "view", "split", "K", "seed", "view_manifest_sha256", "view_content_sha256", "code_sha",
                "checker_build"}


def test_answers_are_on_disk_before_scoring_fails(ev, tmp_path, monkeypatch):
    """I05: the paid sampling survives a scoring crash."""
    def boom(*a, **k):
        raise RuntimeError("boom")

    out = tmp_path / "out"
    with pytest.raises(RuntimeError, match="boom"):
        run_fe(ev, out, monkeypatch, scorer=boom)
    lines = jsonl(out / "answers.jsonl")
    assert len(lines) == 6 and all(set(r) == {"item_id", "group_id", "W", "H", "samples", "truncated", "tokens"}
                                   for r in lines)
    assert all(len(r[k]) == 8 for r in lines for k in ("samples", "truncated", "tokens"))
    fe = json.loads((out / "final_eval.json").read_text())
    assert fe["partial"] is True and RUN_SHA_KEYS <= set(fe) and set(fe["seconds"]) == {"load_and_sample"}
    assert not (out / "metrics.jsonl").exists() and not (out / "items.jsonl").exists()


def test_final_eval_json_is_never_half_written(ev, tmp_path, monkeypatch):
    """I06: a crash while writing leaves final_eval.json absent or valid JSON (temp file, then os.replace)."""
    real, n = Path.write_text, [0]

    def flaky(self, *a, **k):
        n[0] += 1
        if n[0] == 2:
            real(self, '{"half": ')
            raise OSError("disk full")
        return real(self, *a, **k)

    monkeypatch.setattr(Path, "write_text", flaky)
    out = tmp_path / "out"
    with pytest.raises(OSError):
        run_fe(ev, out, monkeypatch)
    monkeypatch.setattr(Path, "write_text", real)
    f = out / "final_eval.json"
    assert not f.exists() or json.loads(f.read_text())


def test_full_run_replaces_the_partial(ev, tmp_path, monkeypatch):
    """I07, I10, I21: four files; items.jsonl is per_item of metrics.jsonl; answers keep types and order; the
    truncated legal answer scores 0; tokens reach metrics."""
    from objdel.geom import engine as _gpucheck
    from objdel.policy import rows as policy_rows

    out = tmp_path / "out"
    run_fe(ev, out, monkeypatch)
    assert sorted(p.name for p in out.iterdir()) == ["answers.jsonl", "final_eval.json", "items.jsonl",
                                                    "metrics.jsonl"]
    fe = json.loads((out / "final_eval.json").read_text())
    assert fe["partial"] is False and fe["schema"] == 2 and RUN_SHA_KEYS <= set(fe)
    assert (fe["K"], fe["seed"], fe["split"], fe["items"], fe["patients"]) == (8, 11, "val", 6, 3)
    assert fe["checker_build"] == _gpucheck.build() and set(fe["noise_floor"]) == {"item_sd", "mean_sd"}
    assert fe["view_manifest_sha256"] == hashlib.sha256((ev / "manifest.json").read_bytes()).hexdigest()
    assert fe["view_content_sha256"] == FE.Q.view_content_sha256(ev)
    ans, met, items = jsonl(out / "answers.jsonl"), jsonl(out / "metrics.jsonl"), jsonl(out / "items.jsonl")
    ids = [r["item_id"] for r in policy_rows.records(ev, "lidc", "val") if r["rl"]]
    assert [r["item_id"] for r in ans] == [r["item_id"] for r in met] == [r["item_id"] for r in items] == ids
    assert items == [EV.per_item(m) for m in met]
    for r, m in zip(ans, met, strict=True):
        assert all(type(t) is int for t in r["tokens"]) and all(type(t) is bool for t in r["truncated"])
        assert type(r["W"]) is int and type(r["H"]) is int and all(type(s) is str for s in r["samples"])
        assert m["samples"]["tokens"] == r["tokens"]  # I21: score_items passes tokens to metrics
        s = m["samples"]
        assert s["legal"][0] == 1 and s["legal"][4:] == [0, 0, 0, 0]  # bowtie, garbage, [] and truncated: illegal
        assert s["iou_tc"][7] == 0 and s["iou_tc"][0] > 0  # I10: the truncated copy of a legal answer scores 0
    assert jsonl(out / "answers.jsonl") == ans


def rescore(src, ev, out, split="val"):
    return CLI.main(["--rescore", str(src), "--view", str(ev), "--split", split, "--out", str(out)])


def strip_seconds(doc):
    return {k: v for k, v in doc.items() if k not in ("seconds", "rescored_from", "engine")}


def test_rescore_is_bit_identical(ev, tmp_path, monkeypatch):
    """I08, I23: per-sample legal, iou_tc, f5, g_area equal (==); items and metrics lines equal; answers order kept."""
    src, out = tmp_path / "src", tmp_path / "re"
    run_fe(ev, src, monkeypatch)
    assert rescore(src, ev, out) in (None, 0)
    a, b = jsonl(src / "metrics.jsonl"), jsonl(out / "metrics.jsonl")
    for x, y in zip(a, b, strict=True):
        for k in ("legal", "iou_tc", "f5", "g_area"):
            assert x["samples"][k] == y["samples"][k] and np.array_equal(x["samples"][k], y["samples"][k]), k
    assert [json.dumps(x, sort_keys=True) for x in a] == [json.dumps(y, sort_keys=True) for y in b]
    assert (src / "items.jsonl").read_bytes() == (out / "items.jsonl").read_bytes()
    assert (src / "answers.jsonl").read_bytes() == (out / "answers.jsonl").read_bytes()
    fa, fb = (json.loads((d / "final_eval.json").read_text()) for d in (src, out))
    assert strip_seconds(fa) == strip_seconds(fb) and fb["rescored_from"] == str(src)


def test_rescore_loads_no_model(ev, tmp_path, monkeypatch):
    """I08: a rescore in a fresh process imports no torch, vllm, transformers, trl or g_probe."""
    src = tmp_path / "src"
    run_fe(ev, src, monkeypatch)
    code = ("import sys\n"
            "from experiments.phase2 import final_eval as CLI\n"
            f"CLI.main(['--rescore', {str(src)!r}, '--view', {str(ev)!r}, '--split', 'val', '--out', "
            f"{str(tmp_path / 're')!r}])\n"
            "print(sorted(m for m in ('torch', 'vllm', 'transformers', 'trl', 'experiments.phase2.g_probe')"
            " if m in sys.modules))\n")
    p = subprocess.run([sys.executable, "-I", "-c", code], cwd=ROOT, check=True, capture_output=True, text=True)
    assert p.stdout.splitlines()[-1] == "[]" and (tmp_path / "re" / "final_eval.json").is_file()


def test_rescore_guards_the_view_and_the_source(ev, tmp_path, monkeypatch, capsys):
    """I09: a one-byte manifest change exits 2 with no output; a partial source is accepted and the result is not
    partial; a source without final_eval.json or without view_manifest_sha256 exits non-zero with a message."""
    src, real = tmp_path / "src", FE.score_items
    with pytest.raises(RuntimeError):
        run_fe(ev, src, monkeypatch, scorer=lambda *a, **k: (_ for _ in ()).throw(RuntimeError("x")))
    monkeypatch.setattr(FE, "score_items", real)
    assert rescore(src, ev, tmp_path / "ok") in (None, 0)  # a partial source
    assert json.loads((tmp_path / "ok" / "final_eval.json").read_text())["partial"] is False
    m = ev / "manifest.json"
    raw = m.read_bytes()
    for changed in (raw + b" ", raw.rstrip(b"\n") if raw.endswith(b"\n") else raw + b"\n",
                    json.dumps(json.loads(raw)).encode()):
        m.write_bytes(changed)
        with pytest.raises(SystemExit) as e:
            rescore(src, ev, tmp_path / "bad")
        assert e.value.code == 2 and not (tmp_path / "bad").exists()
    m.write_bytes(raw)
    for broken in ("missing", "no_sha"):
        d = tmp_path / broken
        d.mkdir()
        (d / "answers.jsonl").write_bytes((src / "answers.jsonl").read_bytes())
        if broken == "no_sha":
            (d / "final_eval.json").write_text(json.dumps({"partial": True}))
        with pytest.raises(SystemExit) as e:
            rescore(d, ev, tmp_path / f"o_{broken}")
        assert e.value.code not in (0, None) and "final_eval.json" in str(e.value.code) + capsys.readouterr().err
        assert not (tmp_path / f"o_{broken}").exists()


def test_cli_modes_need_exactly_one_source(capsys):
    """A3: exactly one of --ckpt and --rescore; --train is gone (E4: no trainer); exit 2."""
    for argv in (["--view", "v", "--out", "o"],
                 ["--ckpt", "c", "--rescore", "r", "--view", "v", "--out", "o"],
                 ["--ckpt", "c", "--train", "t", "--view", "v", "--out", "o"]):
        with pytest.raises(SystemExit) as e:
            CLI.main(argv)
        assert e.value.code == 2
    assert "--ckpt" in capsys.readouterr().err


def test_hash_view_needs_its_out_and_rescore_parses_without_a_ckpt(tmp_path):
    """U53: --hash-view without --out exits 2; the A3 rescore form is an EvalStage with ckpt None."""
    from objdel.core.spec import EvalStage, cli_parse

    with pytest.raises(SystemExit) as e:
        CLI.main(["--hash-view", str(tmp_path)])
    assert e.value.code == 2
    st = cli_parse(EvalStage, ["--rescore", "d", "--view", "v", "--split", "val", "--out", "o"])
    assert (st.ckpt, st.rescore, st.view, st.split, st.out) == (None, "d", ("v",), ("val",), ("o",))


def test_k_follows_the_flag(ev, tmp_path, monkeypatch):
    """I17 (A3 part): K = 2 gives lists of 2; K = 1 is refused before any sampling."""
    run_fe(ev, tmp_path / "k2", monkeypatch, K=2)
    assert all(len(r["samples"]["iou_tc"]) == 2 for r in jsonl(tmp_path / "k2" / "metrics.jsonl"))
    with pytest.raises(SystemExit):
        run_fe(ev, tmp_path / "k1", monkeypatch, K=1)


def test_limit_rows_is_an_even_deterministic_slice():
    """--limit (the GPU smoke): None or a limit past the end keeps every row; otherwise at most N rows, evenly spaced
    from the first row, the same on every call."""
    rows = [{"key": f"i{i}"} for i in range(100)]
    assert FE.limit_rows(rows, None) is rows and FE.limit_rows(rows, 100) is rows and FE.limit_rows(rows, 500) is rows
    got = FE.limit_rows(rows, 32)
    assert len(got) == 25 and got[0] is rows[0] and [r["key"] for r in got[:3]] == ["i0", "i4", "i8"]  # stride 4
    assert got == FE.limit_rows(rows, 32)
    assert [r["key"] for r in FE.limit_rows(rows, 1)] == ["i0"] and len(FE.limit_rows(rows, 99)) == 50


def test_limit_flag_samples_and_scores_the_slice_only(ev, tmp_path, monkeypatch, capsys):
    calls = []
    monkeypatch.setattr(FE, "load_sampler", fake_sampler(ev, calls=calls))
    CLI.main(["--ckpt", "c", "--view", str(ev), "--split", "val", "--out", str(tmp_path / "o"),
             "--limit", "2"])
    assert len(calls) == 1 and len(calls[0]) == 2
    doc = json.loads((tmp_path / "o" / "final_eval.json").read_text())
    assert doc["limit"] == 2 and doc["items"] == 2 and len(jsonl(tmp_path / "o" / "answers.jsonl")) == 2
    with pytest.raises(SystemExit) as e:
        CLI.main(["--ckpt", "c", "--view", str(ev), "--out", str(tmp_path / "x"), "--limit", "0"])
    assert e.value.code == 2 and "--limit" in capsys.readouterr().err


def test_several_views_in_one_process(ev, tmp_path, monkeypatch):
    """U59: one sampler and one preload per root; one folder per view, each equal to a single-view call; flag
    counts that do not match exit 2 with no output."""
    import shutil

    from objdel.verify import store as serve

    v2 = Path(shutil.copytree(ev, tmp_path / "lidc-val-two"))
    loads, roots, calls = [], [], []
    real = serve.preload
    monkeypatch.setattr(serve, "preload", lambda root, *a, **k: roots.append(Path(root).resolve()) or real(root, *a,
                                                                                                          **k))
    monkeypatch.setattr(FE, "load_sampler", fake_sampler(ev, calls=calls, loads=loads))
    base = ["--ckpt", "c"]
    CLI.main([*base, "--view", str(ev), "--split", "val", "--view", str(v2), "--split", "val",
             "--out", str(tmp_path / "O")])
    assert loads == [1] and sorted(roots) == sorted({ev.resolve(), v2.resolve()}) and len(calls) == 2
    for v in (ev, v2):
        d = tmp_path / "O" / v.name
        assert {"answers.jsonl", "items.jsonl", "final_eval.json"} <= {p.name for p in d.iterdir()}
        assert "load_and_sample" in json.loads((d / "final_eval.json").read_text())["seconds"]
        serve.STORE.clear()
        CLI.main([*base, "--view", str(v), "--split", "val", "--out", str(tmp_path / "single" / v.name)])
        assert (d / "items.jsonl").read_bytes() == (tmp_path / "single" / v.name / "items.jsonl").read_bytes()
    serve.STORE.clear()
    CLI.main([*base, "--view", str(ev), "--split", "val", "--view", str(v2), "--split", "val",
             "--out", str(tmp_path / "O1"), "--out", str(tmp_path / "O2")])
    assert (tmp_path / "O1" / "items.jsonl").is_file() and (tmp_path / "O2" / "items.jsonl").is_file()
    for bad in (["--view", str(ev), "--split", "val", "--view", str(v2), "--out", str(tmp_path / "X")],
                ["--view", str(ev), "--view", str(v2), "--split", "val", "--split", "val",
                 "--out", str(tmp_path / "X"), "--out", str(tmp_path / "Y"), "--out", str(tmp_path / "Z")]):
        with pytest.raises(SystemExit) as e:
            CLI.main([*base, *bad])
        assert e.value.code == 2 and not any((tmp_path / n).exists() for n in "XYZ")


def test_code_sha_precedence(tmp_path, monkeypatch):
    """U47: CODE_SHA (stripped) wins; an empty one counts as absent; then objdel.core.io.git_sha with a .git (a clean
    tree: HEAD); else None, also when git is missing."""
    (tmp_path / ".git").mkdir()
    (tmp_path / "CODE_SHA").write_text("abc123-dirty-0123456789ab\n")
    assert FE.code_sha(tmp_path, required=False) == "abc123-dirty-0123456789ab"
    seen = []

    def git(cmd, **k):
        seen.append(cmd)
        return SimpleNamespace(stdout=b"f" * 40 + b"\n" if cmd[1] == "rev-parse" else b"", returncode=0)

    monkeypatch.setattr("objdel.core.io.subprocess.run", git)  # code_sha calls git through objdel.core.io
    (tmp_path / "CODE_SHA").write_text("")
    assert FE.code_sha(tmp_path, required=False) == "f" * 40
    assert seen == [["git", "rev-parse", "HEAD"], ["git", "status", "--porcelain"]]
    (tmp_path / "CODE_SHA").unlink()
    assert FE.code_sha(tmp_path, required=False) == "f" * 40

    def missing(cmd, **k):
        raise FileNotFoundError("git")

    monkeypatch.setattr("objdel.core.io.subprocess.run", missing)
    assert FE.code_sha(tmp_path, required=False) is None
    (tmp_path / ".git").rmdir()
    assert FE.code_sha(tmp_path, required=False) is None and json.dumps(FE.code_sha(tmp_path, required=False)) == "null"


def test_code_sha_of_a_dirty_tree_is_the_same_locally_and_on_the_box(tmp_path):
    """X2-3: final_eval's code_sha of a dirty tree on the Mac (git, no CODE_SHA) carries the -dirty-<sha12> suffix and
    equals the CODE_SHA sync_up writes for the box (git_sha), which the box (no .git) reads back."""
    from objdel.core.io import git_sha

    repo, box = tmp_path / "repo", tmp_path / "box"
    repo.mkdir()
    box.mkdir()
    git(repo, "init", "-q")
    (repo / "a.py").write_text("a = 1\n")
    git(repo, "add", "a.py")
    git(repo, "commit", "-qm", "one")
    (repo / "a.py").write_text("a = 2\n")
    local = FE.code_sha(repo, required=False)
    assert local.startswith(git(repo, "rev-parse", "HEAD").strip() + "-dirty-") and len(local) == 40 + 7 + 12
    (box / "CODE_SHA").write_text(git_sha(repo) + "\n")
    assert FE.code_sha(box, required=False) == local


def test_code_sha_of_a_real_repo(tmp_path):
    """U47 (b): a tmp repo with one commit."""
    git(tmp_path, "init", "-q")
    (tmp_path / "a.txt").write_text("a")
    git(tmp_path, "add", "a.txt")
    git(tmp_path, "commit", "-qm", "one")
    assert FE.code_sha(tmp_path, required=False) == git(tmp_path, "rev-parse", "HEAD").strip()


def test_view_manifest_sha_is_over_the_raw_bytes(tmp_path):
    """U48: CRLF, unicode and a reformatted manifest each hash as their bytes."""
    from objdel.geom import engine as _gpucheck

    got = set()
    for raw in (b'{"a": 1}\n', b'{"a": 1}\r\n', '{"a": "é"}\n'.encode(), b'{"a":1}'):
        (tmp_path / "manifest.json").write_bytes(raw)
        assert FE.view_manifest_sha256(tmp_path) == hashlib.sha256(raw).hexdigest()
        got.add(FE.view_manifest_sha256(tmp_path))
    assert len(got) == 4 and _gpucheck.build() in ("cuda", "cpu")


def content_view(tmp_path, recs, sets=None):
    v = tmp_path / "view"
    (v / "lidc").mkdir(parents=True, exist_ok=True)
    (v / "lidc" / "test.jsonl").write_text("".join(json.dumps(r) + "\n" for r in recs))
    (v / "manifest.json").write_text(json.dumps({"versions": {"python": "3.12"}, "reward_fingerprint": "f",
                                                 "sets": sets or {"lidc/test": {"target_set_sha256": "t" * 64}}}))
    return v


def test_view_content_hash(tmp_path, monkeypatch, capsys):
    """U54: SHA-256 of the newline-joined (item_id, offset, side, image_sha256) tuples; manifest versions and the
    reward fingerprint do not enter it; an offset, an image SHA or the record order does. --hash-view writes its
    record and exits 0 without a model."""
    from objdel.policy import sampler as sample

    recs = [{"item_id": f"lidc/p/n000/z{i}", "view": {"offset": [i, 2 * i], "side": 50}, "image_sha256": f"{i:064x}",
             "W": 50, "H": 50} for i in range(3)]
    v = content_view(tmp_path, recs)
    want = hashlib.sha256("\n".join(json.dumps([r["item_id"], r["view"]["offset"], r["view"]["side"],
                                                r["image_sha256"]]) for r in recs).encode()).hexdigest()
    assert FE.Q.view_content_sha256(v) == want
    m = json.loads((v / "manifest.json").read_text())
    (v / "manifest.json").write_text(json.dumps(m | {"versions": {"python": "3.13"}, "reward_fingerprint": "g"}))
    assert FE.Q.view_content_sha256(v) == want
    for change in ("offset", "image", "order"):
        r2 = [dict(r, view=dict(r["view"])) for r in recs]
        if change == "offset":
            r2[1]["view"]["offset"] = [2, 2]
        elif change == "image":
            r2[2]["image_sha256"] = "f" * 64
        else:
            r2[0], r2[1] = r2[1], r2[0]
        assert FE.Q.view_content_sha256(content_view(tmp_path / change, r2)) != want, change
    monkeypatch.setattr(sample, "VllmSampler", lambda *a, **k: (_ for _ in ()).throw(AssertionError("model")))
    out = tmp_path / "pre" / "pre_rental.json"
    assert CLI.main(["--hash-view", str(v), "--out", str(out)]) in (None, 0)
    doc = json.loads(out.read_text())
    assert set(doc) == {"view", "items", "content_sha256", "target_set_sha256", "generator_commit", "render_files",
                        "r5_commit", "note"}
    assert (doc["items"], doc["content_sha256"], doc["target_set_sha256"], doc["r5_commit"]) == (3, want, "t" * 64,
                                                                                                None)
    assert set(doc["render_files"]) == set(CLI.render_files()) and "experiments/phase2/build_eval_views.py" in doc[
        "render_files"]
    assert doc["generator_commit"] == doc["render_files"]["experiments/phase2/build_eval_views.py"]
    assert all(c is None or len(c) == 40 for c in doc["render_files"].values()) and "not recorded" in doc["note"]


def halves(first, second, K=8):
    return [first] * (K // 2) + [second] * (K // 2)


def test_noise_floor(tmp_path):
    """U55: equal halves give 0; the hand case gives sd(d) / sqrt 2 and the patient-draw SD of mean(d) / sqrt 2."""
    ids = [f"lidc/p{i // 2}/n000/z{i}" for i in range(4)]
    groups = {i: i.split("/")[1] for i in ids}
    same = [{"item_id": i, "samples": {"iou_tc": [.3, .5, .7, .1] * 2}} for i in ids]
    assert FE.noise_floor(same, groups) == {"item_sd": 0.0, "mean_sd": 0.0}
    pairs = [(.6, .5), (.4, .4), (.7, .5), (.2, .3)]
    res = [{"item_id": i, "samples": {"iou_tc": halves(*p)}} for i, p in zip(ids, pairs, strict=True)]
    d = np.array([np.mean(halves(*p)[:4]) - np.mean(halves(*p)[4:]) for p in pairs])
    assert abs(d.mean() - 0.05) < 1e-12
    D = S.draws(ids, 2000, 0, key=groups.__getitem__)
    nf = FE.noise_floor(res, groups)
    assert nf["item_sd"] == np.std(d, ddof=1) / np.sqrt(2)
    assert nf["mean_sd"] == np.std([d[x].mean() for x in D], ddof=1) / np.sqrt(2)


def identity_fixture(seed=0):
    """40 items in 8 patients, K = 8 answers around a per-item level."""
    rng = np.random.default_rng(seed)
    ids = [f"lidc/p{i % 8}/n{i:03d}/z1" for i in range(40)]
    level = rng.uniform(.1, .9, 40)
    res = [{"item_id": i, "samples": {"iou_tc": np.clip(x + rng.normal(0, .08, 8), 0, 1).tolist()}}
           for i, x in zip(ids, level, strict=True)]
    return ids, {i: i.split("/")[1] for i in ids}, res


def test_view_identity():
    """U58: r_sh = 2 r_h / (1 + r_h) from the half means; a same-view R5 column passes; the same values shuffled
    across items fail, and the message names the correlation and the bound."""
    ids, groups, res = identity_fixture()
    iou = np.array([r["samples"]["iou_tc"] for r in res])
    r_h = np.corrcoef(iou[:, :4].mean(1), iou[:, 4:].mean(1))[0, 1]
    old = {i: float(m + n) for i, m, n in zip(ids, iou.mean(1), np.random.default_rng(1).normal(0, .03, 40),
                                              strict=True)}
    v = FE.view_identity(res, old, groups)
    assert abs(v["r_sh"] - 2 * r_h / (1 + r_h)) < 1e-12 and v["pass"] is True
    assert abs(v["r_obs"] - np.corrcoef(iou.mean(1), [old[i] for i in ids])[0, 1]) < 1e-12
    shuffled = dict(zip(ids, np.random.default_rng(2).permutation(list(old.values())).tolist(), strict=True))
    v = FE.view_identity(res, shuffled, groups)
    assert v["pass"] is False and f"{v['r_obs']:.4f}" in v["message"] and f"{v['r_sh']:.4f}" in v["message"]
    with pytest.raises(ValueError, match=ids[0]):
        FE.view_identity(res, {k: x for k, x in old.items() if k != ids[0]}, groups)


def test_tost():
    """A-G3: the 90 % patient CI of the mean paired difference inside +-bound passes; a wider CI is 'not shown'."""
    ids, groups, _ = identity_fixture()
    rng = np.random.default_rng(3)
    old = {i: float(x) for i, x in zip(ids, rng.uniform(.2, .8, 40), strict=True)}
    near = {i: x + float(rng.normal(0, .005)) for i, x in old.items()}
    t = FE.tost(near, old, groups)
    D = S.draws(ids, 2000, 0, key=groups.__getitem__)
    d = np.array([near[i] - old[i] for i in ids])
    assert t["ci90_patient"] == np.quantile([d[x].mean() for x in D], [.05, .95]).tolist()
    assert t["pass"] is True and (t["items"], t["patients"]) == (40, 8)
    far = {i: x + float(rng.normal(0, .1)) for i, x in old.items()}
    assert FE.tost(far, old, groups)["pass"] is False
    with pytest.raises(ValueError, match=ids[3]):
        FE.tost({k: x for k, x in near.items() if k != ids[3]}, old, groups)


def test_code_sha_writer_marks_a_dirty_tree(tmp_path):
    """U49: vast sync_up's CODE_SHA text (objdel.core.io.git_sha): HEAD when clean; HEAD-dirty-<sha256(D)[:12]> with D =
    git diff HEAD plus the sorted untracked, not-ignored paths; staged-only changes count; ignored files do not; git
    failing raises."""
    from objdel.core.io import git_sha

    git(tmp_path, "init", "-q")
    (tmp_path / ".gitignore").write_text("ignored.txt\n")
    (tmp_path / "a.py").write_text("a = 1\n")
    git(tmp_path, "add", ".")
    git(tmp_path, "commit", "-qm", "one")
    head = git(tmp_path, "rev-parse", "HEAD").strip()
    assert git_sha(tmp_path) == head and len(head) == 40
    (tmp_path / "ignored.txt").write_text("x")
    assert git_sha(tmp_path) == head

    def suffix():
        d = subprocess.run(["git", "diff", "HEAD"], cwd=tmp_path, capture_output=True, check=True).stdout
        return hashlib.sha256(d).hexdigest()[:12]

    (tmp_path / "a.py").write_text("a = 2\n")
    first = git_sha(tmp_path)
    assert first == f"{head}-dirty-{suffix()}" and git_sha(tmp_path) == first
    (tmp_path / "a.py").write_text("a = 2\n#\n")
    assert git_sha(tmp_path) not in (first, head)
    (tmp_path / "a.py").write_text("a = 1\n")
    assert git_sha(tmp_path) == head  # modified, then reverted: clean
    (tmp_path / "a.py").write_text("a = 3\n")
    git(tmp_path, "add", "a.py")
    assert git_sha(tmp_path) == f"{head}-dirty-{suffix()}"  # staged only
    git(tmp_path, "reset", "-q", "--hard")
    (tmp_path / "new.py").write_text("b = 1\n")
    assert git_sha(tmp_path) == f"{head}-dirty-{hashlib.sha256(b'new.py\n').hexdigest()[:12]}"
    with pytest.raises((subprocess.CalledProcessError, OSError)):  # never an empty or garbage CODE_SHA
        git_sha(tmp_path / "missing")


SPECS = ROOT / "specs"


def test_a3_specs_are_valid_plans(monkeypatch):
    """I20 (E6.2: the plans are run specs): a3_measure.yaml holds one final_eval stage on lidc-val-edge then
    lidc-val-shifted (each split val and its own out), after the build-if-missing stage; a3_test.yaml one final_eval
    stage on lidc-test-shifted; no nsys stage; CH.chain fills every placeholder and the runner adds -I to the -m
    stages."""
    import time

    from objdel.core.spec import load_spec
    from objdel.infra.vast import __main__ as P
    from objdel.infra.vast import chain as CH

    sent, evals = [], {}
    run = {"deadline": time.time() + 1e5, "instance": 1, "dph": 0.5, "started": 0.0}
    monkeypatch.setattr(P, "get_run", lambda r: run)
    monkeypatch.setattr(P, "remote", lambda run, cmd, *a, **k: sent.append(cmd))
    monkeypatch.setattr(P, "ledger", lambda **k: None)
    for name in ("a3_measure", "a3_test"):
        spec = load_spec(SPECS / f"{name}.yaml")
        (evals[name],) = [s for s in spec.stages if s.kind == "eval"]
        chained = CH.chain(SPECS / f"{name}.yaml", "r")
        assert not any("{" in x or "}" in x or "nsys_capture.sh" in x for s in chained for x in s), name
        P.cmd_run(SimpleNamespace(run_id="r", module=f"specs/{name}.yaml", args=[], cap_hours=1))
        assert "/opt/venv/bin/python -I -m experiments.phase2.final_eval" in sent[-1]
    m, t = evals["a3_measure"], evals["a3_test"]
    assert m.view == ("{data}/eval_views/lidc-val-edge", "{data}/eval_views/lidc-val-shifted")
    assert m.split == ("val", "val") and m.out == ("{rd}/eval_val_edge", "{rd}/eval_val")
    hard = [s for s in load_spec(SPECS / "a3_measure.yaml").stages if not s.soft]
    assert "build() == \"cuda\"' || (bash experiments/cuda_checker/run_box.sh --build-only" in hard[0].argv[2]
    assert hard.index(m) == 1
    assert (t.view, t.split) == (("{data}/eval_views/lidc-test-shifted",), ("test",))


def test_answers_reader_memory(tmp_path):
    """Q04 (memory): 2,200 items x 8 answers of about 1,000 characters (about 18 MB) read under 200 MB peak."""
    import tracemalloc

    line = json.dumps({"item_id": "lidc/x", "group_id": "g", "W": 50, "H": 50, "samples": ["9" * 1000] * 8,
                       "truncated": [False] * 8, "tokens": [400] * 8}) + "\n"
    f = tmp_path / "answers.jsonl"
    f.write_text(line * 2200)
    tracemalloc.start()
    rows = FE.read_jsonl(f)
    peak = tracemalloc.get_traced_memory()[1]
    tracemalloc.stop()
    assert len(rows) == 2200 and peak < 200e6


@pytest.mark.slow
def test_rescore_of_300_items_is_fast(tmp_path, clean_store, monkeypatch, synthetic_items):
    """Q04 (time): the rescore of 300 synthetic items takes under 60 s."""
    from objdel.data import prepare as data_prepare
    from objdel.data import store
    from objdel.data.store import Item
    from objdel.policy import prompt as policy_prompt

    y, x = np.ogrid[:100, :100]
    rng = np.random.default_rng(0)
    items = [Item(f"lidc/1.2.{i // 6}/n000/z{i:03d}", "lidc", f"LIDC-IDRI-{i // 6:04d}", "val",
                  rng.integers(-1000, 400, (100, 100), dtype=np.int16),
                  {"nodule": np.stack([(y - 50) ** 2 + (x - 50 - d) ** 2 <= (8 + i % 7) ** 2 for d in (0, 1)])},
                  ["1", "2"], {}) for i in range(300)]
    riga = [it for it in synthetic_items(("val",)) if it.dataset == "riga"]  # prepare builds every dataset
    store.write(riga + sorted(items, key=lambda it: it.item_id), tmp_path / "store")
    data_prepare.prepare(tmp_path / "store", tmp_path / "view", splits=("val",), procs=4, prompts=policy_prompt.PROMPTS)
    view = tmp_path / "view"
    run_fe(view, tmp_path / "src", monkeypatch)
    t = time.time()
    rescore(tmp_path / "src", view, tmp_path / "re")
    assert time.time() - t < 60 and len(jsonl(tmp_path / "re" / "items.jsonl")) == 300


def test_rescore_is_repeatable_across_processes(ev, tmp_path, monkeypatch):
    """E02 (A3 part): two fresh processes with other hash seeds and thread counts write the same final_eval.json,
    apart from seconds, code_sha and the engine record's eval_sets wall time."""
    src = tmp_path / "src"
    run_fe(ev, src, monkeypatch)
    docs = []
    for seed, threads in (("1", "1"), ("2", "4")):
        out = tmp_path / f"re{seed}"
        env = os.environ | {"PYTHONHASHSEED": seed, "OMP_NUM_THREADS": threads}
        subprocess.run([sys.executable, "-I", "-m", "experiments.phase2.final_eval", "--rescore", str(src), "--view",
                        str(ev), "--split", "val", "--out", str(out)], cwd=ROOT, env=env, check=True,
                       capture_output=True)
        doc = json.loads((out / "final_eval.json").read_text())
        doc["engine"].pop("eval_sets_s")  # a wall time, as seconds
        docs.append({k: v for k, v in doc.items() if k not in ("seconds", "code_sha")})
    assert docs[0] == docs[1]


# --------------------------------------------------------------------------------- frozen views (B1.6, E05, I28)
def _view(root, name="lidc-test-shifted", split="test", n=3, target_set="t" * 64):
    """A view folder with the fields the content hash reads (item_id, view.offset, view.side, image_sha256)."""
    d = root / name
    (d / "lidc").mkdir(parents=True)
    recs = [{"item_id": f"lidc/1.2.{i}/n000/z001", "view": {"offset": [i, -i], "side": 50, "op": "id"},
             "image_sha256": f"{i:064x}", "targets": ["x"]} for i in range(n)]
    (d / "lidc" / f"{split}.jsonl").write_text("".join(json.dumps(r) + "\n" for r in recs))
    m = {"view": name, "versions": {"python": "3.12.0"}, "reward_fingerprint": "f0",
         "sets": {f"lidc/{split}": {"items": n, "target_set_sha256": target_set}}}
    (d / "manifest.json").write_text(json.dumps(m))
    return d


def _content(d):  # independent of objdel.data.quality: the plan's newline-joined tuples
    recs = [json.loads(x) for x in next((d / "lidc").glob("*.jsonl")).read_text().splitlines()]
    return hashlib.sha256("\n".join(json.dumps([r["item_id"], r["view"]["offset"], r["view"]["side"],
                                                 r["image_sha256"]]) for r in recs).encode()).hexdigest()


@pytest.fixture
def frozen(tmp_path, monkeypatch):
    """(view dir, frozen_views.json path): a test view and its entry, final_eval pointed at the file."""
    d = _view(tmp_path / "views")
    f = tmp_path / "frozen_views.json"
    f.write_text(json.dumps({d.name: {"content_sha256": _content(d), "target_set_sha256": "t" * 64}}))
    monkeypatch.setattr(FE, "FROZEN_VIEWS", f)
    return d, f


def _edit(path, fn):
    path.write_text(fn(path.read_text()))


def _swap_first_two(text):
    a, b, *rest = text.splitlines(keepends=True)
    return "".join([b, a, *rest])


@pytest.mark.parametrize("mutate, key", [
    (lambda d: _edit(d / "lidc/test.jsonl", lambda t: t.replace('"offset": [1, -1]', '"offset": [2, -1]')),
     "content_sha256"),
    (lambda d: _edit(d / "lidc/test.jsonl", lambda t: t.replace(f'{1:064x}', f'{9:064x}')), "content_sha256"),
    (lambda d: _edit(d / "lidc/test.jsonl", _swap_first_two), "content_sha256"),
    (lambda d: _edit(d / "manifest.json", lambda t: t.replace("t" * 64, "u" * 64)), "target_set_sha256"),
])
def test_final_eval_refuses_a_changed_frozen_view(frozen, mutate, key, tmp_path):
    """E05: one offset +1, one image sha, two records reordered, the target set: each stops before any model load."""
    d, _ = frozen
    mutate(d)
    with pytest.raises(SystemExit) as e:
        CLI.main(["--ckpt", "c", "--view", str(d), "--split", "test", "--out", str(tmp_path / "o")])
    assert key in str(e.value.code) and d.name in str(e.value.code)
    assert not (tmp_path / "o").exists()


def test_frozen_check_ignores_manifest_versions_and_reward_fingerprint(frozen):
    d, _ = frozen
    _edit(d / "manifest.json", lambda t: t.replace("3.12.0", "3.13.1").replace('"f0"', '"f1"'))
    assert FE.check_view(d, "test") is None


def test_frozen_entry_is_a_content_hash_not_the_manifest_sha(frozen):
    """A frozen entry that holds view_manifest_sha256 (the manifest file's SHA-256) never passes."""
    d, f = frozen
    f.write_text(json.dumps({d.name: {"content_sha256": hashlib.sha256((d / "manifest.json").read_bytes()).hexdigest(),
                                      "target_set_sha256": "t" * 64}}))
    with pytest.raises(SystemExit, match="content_sha256"):
        FE.check_view(d, "test")


def test_frozen_view_without_its_record_is_refused(frozen, tmp_path, monkeypatch):
    """No frozen_views.json (a box copy without docs/reference): a clear refusal, not a KeyError trace."""
    d, _ = frozen
    monkeypatch.setattr(FE, "FROZEN_VIEWS", tmp_path / "missing.json")
    assert FE.frozen_entry(d.name) is None
    with pytest.raises(SystemExit, match="frozen_views.json"):
        FE.check_view(d, "test")


def test_a_view_outside_the_frozen_list_is_not_hash_checked(tmp_path, monkeypatch):
    monkeypatch.setattr(FE, "FROZEN_VIEWS", tmp_path / "missing.json")
    assert FE.check_view(_view(tmp_path, "lidc-val-rescaled", "val"), "val") is None


def test_the_test_view_needs_the_explicit_split_flag(frozen, tmp_path):
    """I28: --split defaults to val, so a test view is read only with --split test."""
    d, _ = frozen
    with pytest.raises(SystemExit, match="--split test"):
        CLI.main(["--ckpt", "c", "--view", str(d), "--out", str(tmp_path / "o")])


def test_no_phase_b_script_but_final_eval_scores_the_test_view():
    """I28 (grep): the scoring code (experiments/phase2/*.py; objdel.verify, policy, rl, sft; data audit, probes and
    prepare) names the test view only in final_eval (freeze_views hashes it, never scores it), and every config line
    that names it is a final_eval step."""
    lib = [*(f for d in ("verify", "policy", "rl", "sft") for f in (ROOT / "src/objdel" / d).rglob("*.py")),
           *(ROOT / "src/objdel/data" / f for f in ("audit.py", "probes.py", "prepare.py"))]
    hits = [p.relative_to(ROOT).as_posix() for p in [*(ROOT / "experiments/phase2").glob("*.py"), *lib]
            if re.search(r"lidc-test|--split\W+test", p.read_text())]
    assert sorted(hits) == ["experiments/phase2/final_eval.py", "experiments/phase2/freeze_views.py",
                            "src/objdel/policy/final_eval.py"]
    from objdel.infra.vast import chain as C

    for y in [*(ROOT / "experiments/phase2/configs").glob("*.yaml"), *(ROOT / "specs").glob("*.yaml")]:
        for step in C.chain(y, "r"):
            if "lidc-test" in " ".join(step):
                assert "experiments.phase2.final_eval" in step, (y.name, step)


def adapter_dir(tmp_path, *, config=True):
    d = tmp_path / "adapter"
    d.mkdir()
    (d / "adapter_model.safetensors").write_bytes(b"weights")
    if config:
        (d / "adapter_config.json").write_text("{}")
    return d


def test_adapter_run_records_the_folder_and_its_hash(ev, tmp_path, monkeypatch):
    """I48: with --adapter final_eval.json holds the folder and the SHA-256 of its weights; without, no such key."""
    import hashlib

    d = adapter_dir(tmp_path)
    run_fe(ev, tmp_path / "o", monkeypatch, "--adapter", str(d))
    doc = json.loads((tmp_path / "o" / "final_eval.json").read_text())
    assert doc["adapter"] == str(d) and doc["adapter_sha256"] == hashlib.sha256(b"weights").hexdigest()
    run_fe(ev, tmp_path / "o2", monkeypatch)
    assert "adapter" not in json.loads((tmp_path / "o2" / "final_eval.json").read_text())


def test_adapter_folder_without_its_config_exits_2_naming_the_file(ev, tmp_path, monkeypatch, capsys):
    with pytest.raises(SystemExit) as e:
        run_fe(ev, tmp_path / "o", monkeypatch, "--adapter", str(adapter_dir(tmp_path, config=False)))
    assert e.value.code == 2 and "adapter_config.json" in capsys.readouterr().err


def test_load_sampler_serves_the_adapter_unmerged_with_the_checkpoints_stop_ids(tmp_path, monkeypatch):
    """E4 step 3: one VllmSampler of --ckpt (no trainer), serving --adapter unmerged (option F) when given; params_for
    at K and seed with the checkpoint's stop ids."""
    from types import SimpleNamespace as NS

    from objdel.policy import sampler as sample

    made = []

    class Fake:
        processor = "P"

        def __init__(self, ckpt, model, adapter=None):
            made.append((ckpt, model, adapter))

    monkeypatch.setattr(sample, "VllmSampler", Fake)
    monkeypatch.setattr(sample, "stop_ids", lambda path, proc: {9, 2} if (path, proc) == (Path("c"), "P") else None)
    base = dict(model="qwen3-vl-2b", ckpt="c", k=8, seed=11)
    s, params = FE.load_sampler(NS(**base, adapter=None))
    assert made == [(Path("c"), "qwen3-vl-2b", None)]
    assert params == sample.params_for("qwen3-vl-2b", n=8, seed=11, stop_token_ids=[2, 9])
    FE.load_sampler(NS(**base, adapter=str(tmp_path)))
    assert made[1] == (Path("c"), "qwen3-vl-2b", tmp_path)

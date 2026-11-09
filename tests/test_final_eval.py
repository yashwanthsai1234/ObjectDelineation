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


def result(item, legal, iou, f5, p, r, d):
    return {"item_id": item, "samples": {"legal": legal, "iou_tc": iou, "f5": f5},
            "objects": {"nodule": {"precision_vote": p, "recall_vote": r, "dice_vote": d}}}


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


def test_per_item_takes_the_k_sample_means():
    it = FE.per_item(result("lidc/a", [1, 0], [0.8, 0.0], [0.9, 0.0], [0.6, 0.0], [1.0, 0.0], [0.75, 0.0]))
    assert it == {"item_id": "lidc/a", "legal": 0.5, "iou_tc": 0.4, "f5": 0.45, "precision_vote": 0.3,
                  "recall_vote": 0.5, "dice_vote": 0.375}


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
    assert doc["metrics_absent"] == ["precision_vote", "recall_vote", "dice_vote"]
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

"""DOC6: docs/results/headline.json is copied from final_eval files, never typed."""

import json

import pytest

from experiments.phase2 import headline as H

KEYS = ("iou_tc", "f5", "precision_vote", "recall_vote", "dice_vote")


def fake(split, view, items, base=0.5, schema=1):
    d = {"ckpt": "/c", "view": f"/d/{view}", "split": split, "K": 8, "seed": 11, "items": items,
         "legal": {"mean": 0.99412345, "ci95": [0.98, 1.0], "ci95_patient": [0.97, 1.0]}}  # fmt: skip
    for i, k in enumerate(KEYS):
        d[k] = {"mean": base + i / 100 + 1e-9, "ci95": [0.1, 0.9], "ci95_patient": [0.2, 0.8]}
    if schema == 2:
        d |= {"schema": 2, "gate_ci": "ci95_patient"}
    return d


def write(tmp_path, name, doc):
    p = tmp_path / name
    p.write_text(json.dumps(doc))
    return str(p)


def run(tmp_path, test, val, *extra):
    out = tmp_path / "h.json"
    code = H.main(["--test", test, "--val", val, "--model-ref", "r/m@rev", "--out", str(out), *extra])
    return code, out


def pair(tmp_path, schema=1):
    return (write(tmp_path, "t.json", fake("test", "lidc-test-shifted", 2200, schema=schema)),
            write(tmp_path, "v.json", fake("val", "lidc-val-shifted", 1084, 0.4, schema=schema)))  # fmt: skip


def test_writes_schema_with_unrounded_values(tmp_path):
    t, v = pair(tmp_path)
    code, out = run(tmp_path, t, v)
    h = json.loads(out.read_text())
    assert code == 0 and h["model_ref"] == "r/m@rev" and h["ci_unit"] == "item" and h["K"] == 8
    assert h["test"]["view"] == "lidc-test-shifted" and h["test"]["items"] == 2200 and h["val"]["items"] == 1084
    assert h["test"]["metrics"]["iou_tc"] == {"mean": 0.5 + 1e-9, "ci95": [0.1, 0.9]}
    assert set(h["test"]["metrics"]) == {*KEYS, "legal"}


def test_a_val_file_in_the_test_slot_is_refused(tmp_path):
    t, v = pair(tmp_path)
    code, out = run(tmp_path, v, t)
    assert code != 0 and not out.exists()


def test_a_file_without_view_is_refused(tmp_path):
    doc = fake("test", "x", 1)
    del doc["view"]
    _, v = pair(tmp_path)
    code, out = run(tmp_path, write(tmp_path, "n.json", doc), v)
    assert code != 0 and not out.exists()


def test_schema_2_pair_copies_the_patient_ci(tmp_path):
    t, v = pair(tmp_path, schema=2)
    code, out = run(tmp_path, t, v)
    h = json.loads(out.read_text())
    assert code == 0 and h["ci_unit"] == "patient"
    assert h["test"]["metrics"]["f5"]["ci95"] == [0.2, 0.8]


def test_mixed_schemas_are_refused(tmp_path):
    t, _ = pair(tmp_path)
    v = write(tmp_path, "v2.json", fake("val", "lidc-val-shifted", 9, schema=2))
    code, out = run(tmp_path, t, v)
    assert code != 0 and not out.exists()


def test_a_null_gate_ci_is_refused(tmp_path):
    doc = fake("test", "lidc-test-shifted", 5, schema=2)
    doc["f5"]["ci95_patient"] = None
    _, v = pair(tmp_path, schema=2)
    code, out = run(tmp_path, write(tmp_path, "n.json", doc), v)
    assert code != 0 and not out.exists()


def test_run_and_paths_together_or_neither_exit_2(tmp_path):
    t, v = pair(tmp_path)
    assert H.main(["--run", "r5", "--test", t, "--val", v, "--model-ref", "m"]) == 2
    assert H.main(["--model-ref", "m"]) == 2


def test_run_short_form_reads_the_r5_layout(tmp_path, monkeypatch):
    for s, doc in (("test", fake("test", "lidc-test-shifted", 3)), ("val", fake("val", "lidc-val-shifted", 4))):
        (tmp_path / "runs/r5/eval" / f"eval_{s}").mkdir(parents=True)
        (tmp_path / "runs/r5/eval" / f"eval_{s}/final_eval.json").write_text(json.dumps(doc))
    monkeypatch.chdir(tmp_path)
    assert H.main(["--run", "r5", "--model-ref", "m", "--out", "h.json"]) == 0
    assert json.loads((tmp_path / "h.json").read_text())["run"] == "r5"


def test_ci_unit_other_than_item_adds_a_string():
    h = {"ci_unit": "patient", "test": {"view": "lidc-test-shifted", "metrics": {
        k: {"mean": 0.5} for k in (*KEYS, "legal")}}}  # fmt: skip
    from docs_check import headline_strings

    assert "patient-level" in headline_strings(h)
    assert "patient-level" not in headline_strings(h | {"ci_unit": "item"})

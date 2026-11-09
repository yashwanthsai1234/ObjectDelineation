"""experiments.phase2.final_eval: per-item K-sample means and the split summary (mean, 95 % bootstrap, median)."""

import hashlib
import json
import re
from pathlib import Path

import pytest

from experiments.phase2 import final_eval as FE

ROOT = Path(__file__).resolve().parents[1]


def result(item, legal, iou, f5, p, r, d):
    return {"item_id": item, "samples": {"legal": legal, "iou_tc": iou, "f5": f5},
            "objects": {"nodule": {"precision_vote": p, "recall_vote": r, "dice_vote": d}}}


def test_per_item_takes_the_k_sample_means():
    it = FE.per_item(result("lidc/a", [1, 0], [0.8, 0.0], [0.9, 0.0], [0.6, 0.0], [1.0, 0.0], [0.75, 0.0]))
    assert it == {"item_id": "lidc/a", "legal": 0.5, "iou_tc": 0.4, "f5": 0.45, "precision_vote": 0.3,
                  "recall_vote": 0.5, "dice_vote": 0.375}


def test_summary_is_the_mean_over_items_with_an_interval_around_it():
    items = [{"item_id": str(i), "legal": 1.0, "iou_tc": v, "f5": v, "precision_vote": v, "recall_vote": v,
              "dice_vote": v} for i, v in enumerate([0.2, 0.4, 0.6, 0.8, 1.0])]
    s = FE.summarize(items)
    assert s["items"] == 5 and s["iou_tc"]["mean"] == pytest.approx(0.6) and s["iou_tc"]["median"] == 0.6
    lo, hi = s["iou_tc"]["ci95"]
    assert 0.2 <= lo < 0.6 < hi <= 1.0
    assert s["legal"] == {"mean": 1.0, "ci95": [1.0, 1.0]}


def test_cli_needs_a_checkpoint_view_and_training_data(capsys):
    with pytest.raises(SystemExit):
        FE.main(["--view", "v", "--train", "t", "--out", "o"])
    assert "--ckpt" in capsys.readouterr().err


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
        FE.main(["--ckpt", "c", "--view", str(d), "--split", "test", "--train", "t", "--out", str(tmp_path / "o")])
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
        FE.main(["--ckpt", "c", "--view", str(d), "--train", "t", "--out", str(tmp_path / "o")])


def test_no_phase_b_script_but_final_eval_scores_the_test_view():
    """I28 (grep): the scoring code (experiments/phase2/*.py, objdel.eval, objdel.train) names the test view only in
    final_eval, and every config line that names it is a final_eval step."""
    hits = [p.relative_to(ROOT).as_posix() for p in [*(ROOT / "experiments/phase2").glob("*.py"),
                                                       *(ROOT / "src/objdel/eval").glob("*.py"),
                                                       *(ROOT / "src/objdel/train").glob("*.py")]
            if re.search(r"lidc-test|--split\W+test", p.read_text())]
    assert hits == ["experiments/phase2/final_eval.py"]
    for y in (ROOT / "experiments/phase2/configs").glob("*.yaml"):
        for line in y.read_text().splitlines():
            if "lidc-test" in line and not line.lstrip().startswith("#"):
                assert "experiments.phase2.final_eval" in line, (y.name, line)

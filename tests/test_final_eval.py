"""experiments.phase2.final_eval: per-item K-sample means and the split summary (mean, 95 % bootstrap, median)."""

import pytest

from experiments.phase2 import final_eval as FE


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

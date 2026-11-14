"""tools.model.outline_lengths (A8): doctor outline lengths per root and the model answers above each limit."""

import hashlib
import json
import subprocess
from pathlib import Path

import pytest

from objdel.policy import rows as policy_rows
from objdel.verify import store as verify_store
from tests.support.paths import ROOT
from tools.model import outline_lengths as OL


def test_count_takes_the_lengths_as_given_and_skips_none_targets():
    """T11.1: count applies tokens_of as is; None targets are skipped."""
    recs = [{"item_id": "A", "targets": ["abcd", None, "ab"]}, {"item_id": "B", "targets": ["abcdefghij"]}]
    c = OL.count(recs, len)
    assert (c["n_records"], c["n_targets"], c["max"], c["p50"]) == (2, 3, 10, 4.0)
    assert c["mean"] == pytest.approx(16 / 3)
    assert c["above"] == {256: 0, 512: 0, 800: 0, 1000: 0, 1200: 0, 1536: 0}
    plus = OL.count(recs, lambda s: len(s) + 1)
    assert (plus["max"], plus["mean"]) == (11, pytest.approx(19 / 3))


def test_count_above_each_limit_is_strict():
    recs = [{"item_id": "A", "targets": ["x" * n for n in (256, 257, 1200, 1201, 1537)]}]
    assert OL.count(recs, len)["above"] == {256: 4, 512: 3, 800: 3, 1000: 3, 1200: 2, 1536: 1}


def test_count_voting_targets_only():
    """T11.2: only voting indices that hold a target; voting=None counts every target as voting."""
    recs = [{"item_id": "A", "targets": ["x" * 300, "x" * 900, None, "x"]},
            {"item_id": "B", "targets": ["x" * 1300, "x", "x" * 600, "x"]}]
    c = OL.count(recs, len, voting=lambda r: {0, 2})
    assert (c["n_targets"], c["n_voting_targets"]) == (7, 3)
    assert c["above_voting"] == {256: 3, 512: 2, 800: 1, 1000: 1, 1200: 1, 1536: 0}
    assert c["above"][800] == 2
    every = OL.count(recs, len)
    assert every["n_voting_targets"] == every["n_targets"] and every["above_voting"] == every["above"]


def test_longest_keeps_the_20_largest_in_decreasing_order_ties_in_record_order():
    """T11.3."""
    recs = [{"item_id": f"i{j}", "targets": ["x" * (j % 15 + 1), "y" * (j % 15 + 1)]} for j in range(15)]
    top = OL.count(recs, len)["longest"]
    assert len(top) == 20 and top[0] == {"item_id": "i14", "rater": 0, "tokens": 15}
    assert top[1] == {"item_id": "i14", "rater": 1, "tokens": 15}
    assert [t["tokens"] for t in top] == sorted((t["tokens"] for t in top), reverse=True)
    assert len(OL.count(recs[:3], len)["longest"]) == 6


def test_count_of_nothing_has_zero_counts_and_null_statistics():
    """T11.4: no records, or only None targets."""
    for recs in ([], [{"item_id": "A", "targets": [None, None]}]):
        c = OL.count(recs, len)
        assert c["n_targets"] == 0 and c["longest"] == [] and c["above"][256] == 0
        assert all(c[k] is None for k in ("mean", "p50", "p90", "p99", "p99_9", "max"))
    assert OL.count([], len)["n_records"] == 0


def metric(item, legal, iou):
    return {"item_id": item, "samples": {"legal": legal, "iou_tc": iou}}


def test_answers_hand_fixture():
    """T11.5: shares over all 6 answers; legal share, mean IoU_tc and legal-and-good count over those above L."""
    ans = [{"item_id": "i1", "tokens": [100, 300, 600, 1600]}, {"item_id": "i2", "tokens": [50, 1537]}]
    met = [metric("i2", [1, 0], [0.9, 0.0]), metric("i1", [1, 1, 0, 0], [0.2, 0.5, 0.0, 0.0])]
    a = OL.answers(ans, met)
    assert a["n_answers"] == 6
    assert {L: a[L]["share_above"] for L in OL.LIMITS} == {256: 4 / 6, 512: 3 / 6, 800: 2 / 6, 1000: 2 / 6,
                                                           1200: 2 / 6, 1536: 2 / 6}
    assert a[256] == {"share_above": 4 / 6, "legal_share": 0.25, "mean_iou_tc": 0.125, "n_legal_good": 1}
    assert a[1536] == {"share_above": 2 / 6, "legal_share": 0.0, "mean_iou_tc": 0.0, "n_legal_good": 0}
    none = OL.answers([{"item_id": "i2", "tokens": [50, 60]}], [metric("i2", [1, 1], [0.9, 0.5])])
    assert none[256] == {"share_above": 0.0, "legal_share": None, "mean_iou_tc": None, "n_legal_good": 0}


def test_answers_refuse_a_missing_or_repeated_item():
    """T11.5 pins: an answer with no metric row, or a repeated item, raises."""
    with pytest.raises(KeyError, match="i2"):
        OL.answers([{"item_id": "i2", "tokens": [1, 2]}], [metric("i1", [1, 1], [0, 0])])
    with pytest.raises(ValueError, match="i1"):
        OL.answers([{"item_id": "i1", "tokens": [1, 2]}] * 2, [metric("i1", [1, 1], [0, 0])])


def test_dataset_names_are_public_with_no_alias():
    """T11.6: target_tokens and voting replace _tokens and _voting; no source still calls the old names."""
    assert callable(policy_rows.target_tokens) and callable(verify_store.voting)
    assert not any(hasattr(m, n) for m in (policy_rows, verify_store) for n in ("_tokens", "_voting"))
    # [^A-Za-z0-9_] rather than \b: macOS git grep -E does not read \b, so a local pass hid a hit CI found
    hits = subprocess.run(["git", "grep", "-nE", r"(^|[^A-Za-z0-9_])_(tokens|voting)\(", "--", "src"],
                          cwd=ROOT, capture_output=True, text=True).stdout
    assert hits == ""


def test_cli_tokens_add_the_end_token_and_no_special_tokens():
    """T11.7: len(input_ids without special tokens) + 1."""
    calls = []

    def tok(text, add_special_tokens):
        calls.append(add_special_tokens)
        return {"input_ids": list(text)}

    assert OL.tokens_of(tok)("abcde") == 6 and calls == [False]


def write_root(root: Path, rows: dict, manifest: bool):
    for (ds, split), recs in rows.items():
        (root / ds).mkdir(parents=True, exist_ok=True)
        (root / ds / f"{split}.jsonl").write_text("".join(json.dumps(r) + "\n" for r in recs))
    if manifest:
        (root / "manifest.json").write_text('{"view": 1}\n')


def test_report_keys_by_root_dataset_split_with_each_view_manifest_sha(tmp_path):
    """T11.8 (Mac half): rows per existing (dataset, split) file; the view roots carry their manifest SHA-256;
    model_answers per answers file, joined with the metrics.jsonl beside it."""
    rec = {"item_id": "lidc/a/n000/z1", "raters": ["1", "2"], "caches": {}, "targets": ["ab", "abc"]}
    write_root(tmp_path / "prep", {("lidc", "train"): [rec], ("lidc", "val"): [rec]}, manifest=False)
    write_root(tmp_path / "shift", {("lidc", "train"): [rec, rec]}, manifest=True)
    run = tmp_path / "eval_val"
    run.mkdir()
    (run / "answers.jsonl").write_text(json.dumps({"item_id": "x", "tokens": [10, 300]}) + "\n")
    (run / "metrics.jsonl").write_text(json.dumps(metric("x", [1, 1], [0.6, 0.7])) + "\n")
    doc = OL.report(tmp_path / "prep", [tmp_path / "shift"], ["lidc", "riga"], ["train", "val"], len,
                    [run / "answers.jsonl"])
    prep, shift = doc["roots"][str(tmp_path / "prep")], doc["roots"][str(tmp_path / "shift")]
    assert prep["manifest_sha256"] is None and set(prep["lidc"]) == {"train", "val"} and "riga" not in prep
    assert prep["lidc"]["train"]["max"] == 3 and shift["lidc"]["train"]["n_records"] == 2
    assert shift["manifest_sha256"] == hashlib.sha256(b'{"view": 1}\n').hexdigest()
    assert doc["model_answers"][str(run / "answers.jsonl")][256]["n_legal_good"] == 1

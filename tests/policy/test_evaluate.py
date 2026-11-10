"""E4 / I22: evaluation without a trainer (objdel.policy.evaluate.evaluate). The trainer layer (TRL's GRPOTrainer and
FastGRPO: __init__ and _generate) is patched to raise; the sampler is a stub and the scoring is real (score_items and
per_item on the preloaded caches of a tiny view root built in tmp_path)."""
import ast
import hashlib
import json
from pathlib import Path

import numpy as np
import pytest

from objdel.policy import evaluate as EV
from objdel.policy.sampler import Sample
from objdel.verify import store as serve
from objdel.verify.api import VerifyError
from objdel.verify.reward.cache import ExpertCache
from objdel.verify.reward.score import score_batch

ROOT = Path(__file__).resolve().parents[2]
SIDE, CELL = 50, 20  # 50 x 50 pixels; the 0-1000 grid has 20 units per pixel
KEY = "lidc/1.2.0/n000/z000"
ROWS = [{"key": KEY, "W": SIDE, "H": SIDE, "names": ["nodule"]}]


def square(x0, y0, side=20) -> np.ndarray:
    m = np.zeros((SIDE, SIDE), bool)
    m[y0:y0 + side, x0:x0 + side] = True
    return m


def text(x0, y0, side=20) -> str:
    """The outline text of the pixel square [x0, x0 + side) x [y0, y0 + side), clockwise on screen."""
    a, b, c, d = x0 * CELL, y0 * CELL, (x0 + side) * CELL, (y0 + side) * CELL
    return json.dumps({"nodule": [{"outer": [[a, b], [c, b], [c, d], [a, d]], "holes": []}]}, separators=(",", ":"))


ON, SHIFTED = text(10, 10), text(20, 20)  # the raters' square, and one moved by 10 pixels in x and y


@pytest.fixture
def view_root(tmp_path, clean_store) -> Path:
    """Three raters drawing the same 20 x 20 square at pixel (10, 10), preloaded like final_eval does."""
    cache = ExpertCache.build(np.stack([square(10, 10)] * 3))
    (tmp_path / "caches" / KEY).mkdir(parents=True)
    digest = cache.save(tmp_path / "caches" / KEY / "nodule.npz")
    records = tmp_path / "lidc" / "val.jsonl"
    records.parent.mkdir()
    records.write_text(json.dumps({"item_id": KEY, "rl": True, "caches": {"nodule": digest}}) + "\n")
    (tmp_path / "manifest.json").write_text(json.dumps({"files": {"lidc/val.jsonl": {
        "sha256": hashlib.sha256(records.read_bytes()).hexdigest()}}}))
    serve.preload(tmp_path, ["lidc"])
    return tmp_path.resolve()


@pytest.fixture(autouse=True)
def no_trainer(monkeypatch):
    """Any use of the trainer layer fails the test."""
    pytest.importorskip("trl")
    from trl import GRPOTrainer

    from objdel.rl import trl_adapter

    def boom(*a, **k):
        raise AssertionError("evaluate touched the trainer")

    for cls in (GRPOTrainer, trl_adapter.FastGRPO):
        monkeypatch.setattr(cls, "__init__", boom)
        monkeypatch.setattr(cls, "_generate", boom)


class StubSampler:
    def __init__(self, samples):
        self.samples, self.calls = samples, []

    def sample(self, rows, params):
        self.calls.append((rows, params))
        return self.samples


def sample(texts, truncated=None, key=KEY) -> Sample:
    return Sample(key, tuple(texts), tuple((1,) * (i + 1) for i in range(len(texts))),
                  tuple(truncated or (False,) * len(texts)))


def test_evaluate_scores_a_stub_samplers_answers_without_a_trainer(view_root):
    """The sampler is called once with (rows, params); the answers are kept before scoring; IoU of two 20 x 20 squares
    overlapping in 10 x 10 is 100 / 700."""
    s, kept, params = StubSampler([sample([SHIFTED, ON])]), [], {"n": 2, "seed": 3}
    got = EV.evaluate(s, ROWS, params, view_root, keep=kept.append)
    assert len(s.calls) == 1 and s.calls[0][0] is ROWS and s.calls[0][1] is params
    assert kept == [got["answers"]] == [[{"item_id": KEY, "samples": [SHIFTED, ON], "truncated": [False, False],
                                          "tokens": [1, 2]}]]
    assert got["results"][0]["samples"]["iou_tc"] == pytest.approx([100 / 700, 1.0], abs=1e-12)
    item = got["items"][0]
    assert item["item_id"] == KEY and item["legal"] == 1.0
    assert item["iou_tc"] == pytest.approx((100 / 700 + 1) / 2, abs=1e-12)
    assert item["sel_oracle"] == 1.0 and item["sel_mean_of_k"] == pytest.approx((100 / 700 + 1) / 2, abs=1e-12)


def test_truncated_and_unparseable_answers_score_as_the_trainer_path_scores_them(view_root):
    """A truncated answer and text that does not parse are results (every term 0, illegal), the same values as the
    trainer's score_batch gives them; nothing crashes."""
    texts, cut = [SHIFTED, SHIFTED, "not json", ""], [False, True, False, False]
    got = EV.evaluate(StubSampler([sample(texts, cut)]), ROWS, {"n": 4}, view_root)
    cache = serve.CacheRef(KEY, str(view_root))
    T = score_batch([(t, c, ("nodule",), cache, SIDE, SIDE) for t, c in zip(texts, cut, strict=True)])
    samples = got["results"][0]["samples"]
    for k in ("iou_tc", "legal", "f5"):
        assert samples[k] == pytest.approx(T[k].tolist(), abs=1e-12), k
    assert samples["iou_tc"][1:] == [0.0, 0.0, 0.0] and samples["legal"] == [1.0, 0.0, 0.0, 0.0]
    assert got["items"][0]["legal"] == 0.25 and np.isfinite(list(got["items"][0].values())[1:]).all()


def test_all_empty_text_is_scored_not_a_crash(view_root):
    got = EV.evaluate(StubSampler([sample(["", ""])]), ROWS, {"n": 2}, view_root)
    item = got["items"][0]
    assert item["iou_tc"] == item["legal"] == item["sel_oracle"] == 0.0


def test_an_unknown_view_key_is_a_verify_error(view_root):
    rows = [{**ROWS[0], "key": "lidc/1.2.0/n999/z000"}]
    with pytest.raises(VerifyError) as e:
        EV.evaluate(StubSampler([sample([ON, ON], key=rows[0]["key"])]), rows, {"n": 2}, view_root)
    assert e.value.code == "unknown_key"


@pytest.mark.parametrize("samples", [
    [sample([ON, ON]), sample([ON, ON], key="lidc/1.2.0/n001/z000")],  # an item_id that is not in the rows
    [sample([ON, ON], key="lidc/1.2.0/n001/z000")],  # the row's sample under another item_id
    [],  # fewer samples than rows
])
def test_samples_that_do_not_match_the_rows_raise(view_root, samples):
    """No silent cut: the samples are one per row, in row order."""
    with pytest.raises(ValueError, match="one Sample per row"):
        EV.evaluate(StubSampler(samples), ROWS, {"n": 2}, view_root)


def test_zero_rows_give_an_empty_result(view_root):
    assert EV.evaluate(StubSampler([]), [], {"n": 2}, view_root) == {"answers": [], "results": [], "items": []}


def test_one_answer_per_row_is_refused(view_root):
    """n = 1: the item metrics need at least 2 answers (metrics' rule); refused, not a wrong number."""
    with pytest.raises(ValueError, match=">= 2 samples"):
        EV.evaluate(StubSampler([sample([ON])]), ROWS, {"n": 1}, view_root)


def test_final_eval_has_no_trainer_path():
    """AST: objdel.policy.final_eval imports no build_trainer and names no _generate."""
    tree = ast.parse((ROOT / "src/objdel/policy/final_eval.py").read_text())
    imported = {a.asname or a.name for n in ast.walk(tree) if isinstance(n, ast.Import | ast.ImportFrom)
                for a in n.names}
    names = {n.attr for n in ast.walk(tree) if isinstance(n, ast.Attribute)}
    names |= {n.id for n in ast.walk(tree) if isinstance(n, ast.Name)}
    names |= {n.name for n in ast.walk(tree) if isinstance(n, ast.FunctionDef | ast.AsyncFunctionDef)}
    assert "build_trainer" not in imported and "build_trainer" not in names and "_generate" not in names
    assert not any("objdel.rl" in (n.module or "") for n in ast.walk(tree) if isinstance(n, ast.ImportFrom))

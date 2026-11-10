"""E4: evaluation without a trainer (objdel.train.evaluate.evaluate; policy/evaluate.py after the E5 move): a sampler's
answers are kept before scoring and scored with score_items and per_item."""
from pathlib import Path

from objdel.train import evaluate as EV
from objdel.train.sample import Sample


class FakeSampler:
    def __init__(self):
        self.calls = []

    def sample(self, rows, params):
        self.calls.append(([r["key"] for r in rows], dict(params)))
        return [Sample(r["key"], ("a", "b"), ((1, 7), (1, 1, 1)), (False, True)) for r in rows]


def test_evaluate_keeps_the_answers_then_scores_them(monkeypatch):
    rows, order = [{"key": "lidc/p/n1/z0"}, {"key": "lidc/p/n2/z0"}], []

    def score_items(rs, done, root):
        order.append("score")
        assert [r["key"] for r in rs] == [r["key"] for r in rows] and root == Path("v")
        return [{"item_id": k, "n": len(done[k]["samples"])} for k in done], []

    monkeypatch.setattr(EV, "score_items", score_items)
    monkeypatch.setattr(EV, "per_item", lambda r: {"item_id": r["item_id"], "twice": 2 * r["n"]})
    s = FakeSampler()
    got = EV.evaluate(s, rows, {"n": 2, "seed": 3}, Path("v"), keep=lambda answers: order.append(("keep", answers)))
    answers = [{"item_id": r["key"], "samples": ["a", "b"], "truncated": [False, True], "tokens": [2, 3]}
               for r in rows]
    assert s.calls == [([r["key"] for r in rows], {"n": 2, "seed": 3})]
    assert order == [("keep", answers), "score"]
    assert got["answers"] == answers
    assert got["results"] == [{"item_id": r["key"], "n": 2} for r in rows]
    assert got["items"] == [{"item_id": r["key"], "twice": 4} for r in rows]

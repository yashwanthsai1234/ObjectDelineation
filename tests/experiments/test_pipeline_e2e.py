"""E01 (Implementation.md A.9.6a): the Mac pipeline end to end on 12 synthetic items (3 patients, 2 nodules each, 2
slices): final_eval (full path), then ceiling --answers on its answers.jsonl, then outline_lengths with a stub
tokenizer. The item counts agree across every file and all JSON dumps with allow_nan=False."""

import json

import numpy as np
import pytest

from objdel.cli.eval import ceiling as CEILING
from objdel.cli.eval import run as CLI
from objdel.policy import final_eval as FE
from objdel.policy import rows as policy_rows
from objdel.policy.sampler import Sample
from tests.experiments.test_ceiling import text_of
from tools.model import outline_lengths as OL

S, K = 100, 8
EMPTY_POLYGON = '{"nodule":[{"outer":[],"holes":[]}]}'


def _disc(cy, cx, r):
    y, x = np.ogrid[:S, :S]
    return (y - cy) ** 2 + (x - cx) ** 2 <= r * r


def _masks(i):
    """The two doctors' masks of item i: nearby discs whose centre moves with i."""
    cy, cx = 35 + 4 * (i % 4), 35 + 5 * (i // 4)
    return [_disc(cy, cx, 12), _disc(cy + 1, cx, 11)]


@pytest.fixture(scope="module")
def view(tmp_path_factory, synthetic_items):
    from objdel.data import prepare as data_prepare
    from objdel.data import store
    from objdel.data.store import Item
    from objdel.policy import prompt as policy_prompt

    tmp, rng = tmp_path_factory.mktemp("e01"), np.random.default_rng(0)
    items = [Item(f"lidc/1.2.{p}/n00{n}/z00{s}", "lidc", f"LIDC-IDRI-000{p}", "val",
                  rng.integers(-1000, 400, (S, S), dtype=np.int16), {"nodule": np.stack(_masks(i))}, ["1", "2"], {})
             for i, (p, n, s) in enumerate((p, n, s) for p in range(3) for n in range(2) for s in range(2))]
    store.write([it for it in synthetic_items(("val",)) if it.dataset == "riga"] + items, tmp / "store")
    data_prepare.prepare(tmp / "store", tmp / "view", splits=("val",), procs=1, prompts=policy_prompt.PROMPTS)
    return tmp / "view"


def _texts(rec, i):
    """The K answers: truth, truth shifted 3 px, garbage, empty polygon, truth, truth, small blob, truth."""
    truth = rec["targets"][0]
    shifted = text_of(np.roll(_masks(i)[0], 3, axis=1))
    return [truth, shifted, "garbage", EMPTY_POLYGON, truth, truth, text_of(_disc(20, 20, 2)), truth]


def _fake_sampler(view):
    recs = {r["item_id"]: (i, r) for i, r in enumerate(policy_rows.records(view, "lidc", "val"))}

    class Fake:
        def sample(self, rows, params):
            out = []
            for r in rows:
                texts = _texts(recs[r["key"]][1], recs[r["key"]][0])
                out.append(Sample(r["key"], tuple(texts), tuple(tuple(range(len(t) + 1)) for t in texts),
                                  (False,) * K))
            return out

    return lambda a: (Fake(), {"n": a.k})


def jsonl(path):
    return [json.loads(x) for x in path.read_text().splitlines()]


def test_mini_pipeline_counts_and_rules_agree_across_files(view, tmp_path, clean_store, monkeypatch):
    monkeypatch.setattr(FE, "load_sampler", _fake_sampler(view))
    run = tmp_path / "eval"
    CLI.main(["--ckpt", "c", "--view", str(view), "--split", "val", "--out", str(run), "-k", str(K)])
    doc = json.loads((run / "final_eval.json").read_text())
    assert (doc["items"], doc["nodules"], doc["patients"]) == (12, 6, 3) and doc["gate_ci"] == "ci95_patient"
    assert doc["legal"]["mean"] == 6 / 8  # pinned: the checker refuses garbage and an empty polygon, not a small blob

    items, metrics = jsonl(run / "items.jsonl"), jsonl(run / "metrics.jsonl")
    assert all(m["samples"]["legal"] == [1, 1, 0, 0, 1, 1, 1, 1] for m in metrics)
    assert len(items) == len(metrics) == len(jsonl(run / "answers.jsonl")) == 12
    assert all(it["sel_oracle"] >= it["sel_mean_of_k"] for it in items)
    truth_iou = {m["item_id"]: m["samples"]["iou_tc"][0] for m in metrics}  # answer 0 is a doctor's outline
    assert all(it["sel_medoid"] == truth_iou[it["item_id"]] for it in items)  # the medoid picks a truth copy

    out = tmp_path / "ceiling"
    CEILING.main(["--view", str(view), "--split", "val", "--answers", str(run / "answers.jsonl"), "--out", str(out)])
    ceiling = json.loads((out / "ceiling.json").read_text())
    rows = jsonl(out / "ceiling_items.jsonl")
    assert ceiling["all"]["n"] == len(rows) == 12 and ceiling["all"]["n_patients"] == 3
    assert all(np.isfinite(r["model"]) and r["diff"] == r["model"] - r["doctor"] for r in rows)

    lengths = OL.report(view, [], ["lidc"], ["val"], len, [run / "answers.jsonl"])
    root, answers = lengths["roots"][str(view)]["lidc"]["val"], lengths["model_answers"][str(run / "answers.jsonl")]
    assert root["n_records"] == 12 and answers["n_answers"] == 12 * K

    for x in (doc, ceiling, lengths, *items, *metrics, *rows):
        json.dumps(x, allow_nan=False)  # no NaN or inf anywhere

"""The old train/dataset.py names: objdel.policy.rows, rl.rows, sft.rows; data.prepare; verify.store.cache_loader.

Unit tests run on a tiny synthetic store; the R1 round trip runs on the prepared train set (data/prepared)."""
import hashlib
import json
import shutil
from pathlib import Path

import numpy as np
import pytest

pytest.importorskip("pyarrow")
pytest.importorskip("datasets")

from objdel.data import prepare as data_prepare  # noqa: E402
from objdel.data import store  # noqa: E402
from objdel.data.store import Item  # noqa: E402
from objdel.geom.text import SCHEMAS  # noqa: E402
from objdel.policy import prompt as policy_prompt  # noqa: E402
from objdel.policy import rows as policy_rows  # noqa: E402
from objdel.rl import rows as rl_rows  # noqa: E402
from objdel.sft import rows as sft_rows  # noqa: E402
from objdel.verify import store as verify_store  # noqa: E402
from objdel.verify.reward.score import rater_answers, score_completion  # noqa: E402

PREPARED = Path(__file__).resolve().parents[1] / "data" / "prepared"


@pytest.fixture(scope="module")
def prepared(tmp_path_factory, synthetic_items):
    root = tmp_path_factory.mktemp("prep")
    store.write(synthetic_items(), root / "store")
    manifest = data_prepare.prepare(root / "store", root / "out", splits=("train",), procs=1,
        prompts=policy_prompt.PROMPTS)
    return root, manifest


def test_prepare_counts_skips_and_layout(prepared):
    root, m = prepared
    riga, lidc = m["sets"]["riga/train"], m["sets"]["lidc/train"]
    assert (riga["items"], riga["targets"], riga["rl_items"], riga["disagreement_items"]) == (1, 3, 1, 0)
    assert (lidc["items"], lidc["targets"], lidc["rl_items"], lidc["disagreement_items"]) == (2, 3, 1, 1)
    assert lidc["skipped_no_image"] == 1 and riga["skipped_no_image"] == 0
    assert lidc["target_errors"] == [] and riga["target_errors"] == []
    out = root / "out"
    assert (out / "images" / "riga/MESSIDOR/image001.png").is_file()
    assert (out / "caches" / "riga/MESSIDOR/image001" / "cup.npz").is_file()
    assert not (out / "caches" / "lidc/1.2.1/n000/z003").exists()  # one rater: not RL-eligible, no cache
    assert json.loads((out / "manifest.json").read_text()) == m
    for rel, entry in m["files"].items():
        assert hashlib.sha256((out / rel).read_bytes()).hexdigest() == entry["sha256"]


def test_records_and_cache_loader(prepared):
    root, _ = prepared
    out = root / "out"
    rec = {r["item_id"]: r for r in policy_rows.records(out, "lidc", "train")}
    assert set(rec) == {"lidc/1.2.1/n000/z001", "lidc/1.2.1/n000/z003"}
    r = rec["lidc/1.2.1/n000/z001"]
    assert r["W"] == r["H"] and 1000 % r["W"] == 0 and r["disagreement"] and r["rl"]
    caches = verify_store.cache_loader(out)("lidc/1.2.1/n000/z001")
    assert set(caches) == {"nodule"} and caches["nodule"].k == 2
    assert caches["nodule"].digest == r["caches"]["nodule"]
    assert caches["nodule"].raters.shape == (2, r["H"], r["W"])


def test_cache_loader_rejects_a_swapped_cache(prepared, tmp_path):
    root, _ = prepared
    out = tmp_path / "out"
    shutil.copytree(root / "out", out)
    c = out / "caches" / "riga/MESSIDOR/image001"
    (c / "disc.npz").write_bytes((c / "cup.npz").read_bytes())  # a valid cache file, but the wrong one
    with pytest.raises(ValueError, match="riga/MESSIDOR/image001"):
        verify_store.cache_loader(out)("riga/MESSIDOR/image001")


def test_sft_rows_schema(prepared):
    root, _ = prepared
    rows = sft_rows.sft_rows(root / "out", "lidc", "train")
    assert set(rows.column_names) == {"messages", "image", "key", "rater"}
    assert len(rows) == 3  # 2 raters of z001 + the single rater of z003
    texts = {(r["key"], r["rater"]): r for r in rows}
    row = texts[("lidc/1.2.1/n000/z001", "2")]
    user, assistant = row["messages"]
    assert user["role"] == "user" and assistant["role"] == "assistant"
    assert [c["type"] for c in user["content"]] == ["image",
        "text"] and user["content"][1]["text"] == policy_prompt.PROMPTS["lidc"]
    assert assistant["content"] == [{"type": "text", "text": assistant["content"][0]["text"]}]
    assert row["image"].mode == "RGB" and row["image"].size == (512, 512)  # LIDC model image (frame contract)
    assert '"holes":[[' in assistant["content"][0]["text"]  # rater 2 drew the ring


def test_sft_rows_drop_answers_beyond_the_token_budget(prepared):
    root, _ = prepared
    lengths = sorted(len(r["messages"][1]["content"][0]["text"]) for r in sft_rows.sft_rows(root / "out", "riga",
        "train"))
    tok = lambda text, add_special_tokens=False: {"input_ids": list(text)}  # one token per character
    kept = sft_rows.sft_rows(root / "out", "riga", "train", tokenizer=tok, max_tokens=lengths[1] + 1)  # +1: end of turn
    assert len(kept) == sum(n <= lengths[1] for n in lengths) < len(lengths)


def test_sft_rows_teach_only_voting_raters(tmp_path, disc_mask, synthetic_items):
    """A rater whose whole mask is below the size filter abstains from the topology vote (B2-4): its own answer can
    never be topo-any, so SFT does not teach it; RL keeps the item."""
    rng = np.random.default_rng(1)
    speck = np.zeros((100, 100), bool)
    speck[50:52, 50:52] = True  # 4 px < min_area 5
    item = Item("lidc/1.2.9/n000/z001", "lidc", "LIDC-IDRI-0009", "train",
                rng.integers(-1000, 400, (100, 100), dtype=np.int16),
                {"nodule": np.stack([disc_mask(100, 100, 50, 50, 12), disc_mask(100, 100, 50, 50, 11), speck])},
                ["1", "2", "3"], {})
    riga = [it for it in synthetic_items() if it.dataset == "riga"]  # prepare reads both datasets
    store.write([item, *riga], tmp_path / "store")
    data_prepare.prepare(tmp_path / "store", tmp_path / "out", splits=("train",), procs=1,
        prompts=policy_prompt.PROMPTS)
    (rec,) = policy_rows.records(tmp_path / "out", "lidc", "train")
    assert rec["targets"][2] is not None  # certifiable: dropped for abstaining, not for lack of a target
    assert verify_store.cache_loader(tmp_path / "out")(rec["item_id"])["nodule"].voters == (0, 1)
    assert sorted(sft_rows.sft_rows(tmp_path / "out", "lidc", "train")["rater"]) == ["1", "2"]
    assert rl_rows.grpo_rows(tmp_path / "out", "lidc", "train")["key"] == [rec["item_id"]]


def test_grpo_rows_schema_and_oversampling(prepared):
    root, _ = prepared
    out = root / "out"
    rows = rl_rows.grpo_rows(out, "lidc", "train")
    assert set(rows.column_names) == {"prompt", "image", "key", "W", "H", "names"}
    assert rows["key"] == ["lidc/1.2.1/n000/z001"]  # the one-rater item is not RL-eligible
    (row,) = rows
    assert [c["type"] for c in row["prompt"][0]["content"]] == ["image", "text"] and row["names"] == ["nodule"]
    assert rl_rows.grpo_rows(out, "lidc", "train", oversample=3)["key"] == ["lidc/1.2.1/n000/z001"] * 3
    assert rl_rows.grpo_rows(out, "riga", "train", oversample=3)["key"] == ["riga/MESSIDOR/image001"]  # no disagreement


def rows_hash(rows):
    return hashlib.sha256(json.dumps([{k: v for k, v in r.items() if k != "image"} for r in rows],
                                     sort_keys=True).encode()).hexdigest()


def test_prepare_takes_the_prompt_table(prepared, tmp_path):
    """U39 (E3c): prepare's `prompts` is required (dependency injection: data names no prompt module) and the manifest
    stores the table it is given, so a different table gives a different manifest."""
    import inspect

    root, m = prepared
    assert inspect.signature(data_prepare.prepare).parameters["prompts"].default is inspect.Parameter.empty
    assert m["prompts"] == policy_prompt.PROMPTS
    table = {ds: f"outline the {ds} structures" for ds in policy_prompt.PROMPTS}
    again = data_prepare.prepare(root / "store", tmp_path, splits=("train",), procs=1, prompts=table)
    assert again["prompts"] == table and again["files"] == m["files"]


def test_prepare_and_rows_are_deterministic(prepared, tmp_path):
    root, m = prepared
    # 1 vs 2 workers
    again = data_prepare.prepare(root / "store", tmp_path, splits=("train",), procs=2, prompts=policy_prompt.PROMPTS)
    assert again["files"] == m["files"] and again["sets"] == m["sets"]
    for ds in ("riga", "lidc"):
        assert rows_hash(sft_rows.sft_rows(root / "out", ds, "train")) == rows_hash(sft_rows.sft_rows(tmp_path, ds,
            "train"))
        assert rows_hash(rl_rows.grpo_rows(root / "out", ds, "train", 4)) == rows_hash(rl_rows.grpo_rows(tmp_path, ds,
            "train", 4))


def test_prepare_workers_survive_a_parent_that_used_cuda(prepared, tmp_path, monkeypatch):
    """GPU run (box_a): a worker forked from a parent that has used CUDA fails in the engine ("CUDA: initialization
    error"). Simulated on the Mac: Linux's fork start method, and an engine that fails in every process forked from
    this one. prepare's workers are fresh processes (spawn), so they build the same manifest."""
    import multiprocessing
    import os

    from objdel.geom import engine

    root, m = prepared
    parent, real = os.getpid(), engine.check_masks

    def check_masks(*a):
        if os.getpid() != parent:
            raise RuntimeError("CUDA: initialization error")
        return real(*a)

    monkeypatch.setattr(engine, "check_masks", check_masks)
    start = multiprocessing.get_start_method()
    multiprocessing.set_start_method("fork", force=True)
    try:
        again = data_prepare.prepare(root / "store", tmp_path, splits=("train",), procs=1,
            prompts=policy_prompt.PROMPTS)
    finally:
        multiprocessing.set_start_method(start, force=True)
    assert again["files"] == m["files"] and again["sets"] == m["sets"]


def test_every_process_pool_spawns():
    """The same rule in every builder: a process pool is made by get_context("spawn") (build_rl_rows, for one, runs the
    engine in its parent before it starts the pool); never the platform default (fork on Linux)."""
    import re

    root = Path(__file__).resolve().parents[1]
    files = sorted([*(root / "src").rglob("*.py"), *(root / "experiments").rglob("*.py")])
    bad = [f"{p.relative_to(root)}:{n}" for p in files for n, line in enumerate(p.read_text().splitlines(), 1)
           if re.search(r"\bPool\(", line) and 'get_context("spawn").Pool(' not in line and "`" not in line]
    assert bad == []


def round_trip(root, dataset, split):
    """Every SFT target of every RL item (voting raters, B2-4), scored against its own item's caches (R1); returns
    (scored, abstaining raters' targets, failures)."""
    load = verify_store.cache_loader(root)
    bad, n, abstained = [], 0, 0
    for r in policy_rows.records(root, dataset, split):
        if not r["rl"]:
            continue
        caches = load(r["item_id"])
        voting = {i for rs in rater_answers(caches, SCHEMAS[dataset]).values() for i in rs}
        for i, (rater, text) in enumerate(zip(r["raters"], r["targets"], strict=True)):
            if text is None:
                continue
            if i not in voting:
                abstained += 1
                continue
            s = score_completion(text, False, SCHEMAS[dataset], caches, r["W"], r["H"])
            n += 1
            if not (s["format_ok"] == 1 and s["legal"] == 1 and s["iou_tc"] >= 0.99 and s["topo_any"] == 1):
                bad.append((r["item_id"], rater, s))
    return n, abstained, bad


def test_targets_score_as_perfect_answers(prepared):
    root, _ = prepared
    for ds in ("riga", "lidc"):
        n, _, bad = round_trip(root / "out", ds, "train")
        assert n and not bad


@pytest.mark.data
@pytest.mark.parametrize("dataset", ["riga", "lidc"])
def test_r1_round_trip_on_the_prepared_train_set(dataset):
    if not (PREPARED / "manifest.json").is_file():
        pytest.skip(f"no prepared set at {PREPARED} (run experiments/phase1/prepare_rl.py)")
    m = json.loads((PREPARED / "manifest.json").read_text())
    rel = f"{dataset}/train.jsonl"
    assert hashlib.sha256((PREPARED / rel).read_bytes()).hexdigest() == m["files"][rel]["sha256"]
    n, abstained, bad = round_trip(PREPARED, dataset, "train")
    assert n + abstained == m["sets"][f"{dataset}/train"]["rl_targets"] and not bad[:5]
    rl = {r["item_id"] for r in policy_rows.records(PREPARED, dataset, "train") if r["rl"]}
    assert sum(k in rl for k in sft_rows.sft_rows(PREPARED, dataset, "train")["key"]) == n  # SFT teaches exactly these

"""experiments.phase2.build_rl_rows (Implementation.md B1.4) on the test_data_quality fixture store: 30 synthetic train
items, the val and test views, one 13-step block of 2 rows (shares 0.5 / 0.5 / 0, seed 4, draws 1, weight 1), so the
off-centre rows are the train shifted law and every hard gate passes on a clean build."""
import json
import types

import cv2
import pytest
from test_data_quality import PLAN, SHARES, Tok, stub_ceiling, world  # noqa: F401 (world: the module fixture)

from objdel.core.hashing import pixels_sha256
from objdel.data import quality as Q
from objdel.data import windows as Wn


@pytest.fixture
def run(world, tmp_path, monkeypatch):  # noqa: F811
    from experiments.phase2 import build_rl_rows as R
    from objdel.data import splits
    from objdel.policy import models
    stub_ceiling(monkeypatch, [0.9, 0.8])
    monkeypatch.setattr(models, "tokenizer", lambda model: (Tok(50), 1536))
    monkeypatch.setattr(splits, "val_half", lambda g: "shape" if int(g[-1]) % 2 else "select")
    ex = tmp_path / "train_exclusions.json"
    ex.write_text("{}")
    frozen = tmp_path / "frozen_views.json"
    frozen.write_text(json.dumps({v: {"content_sha256": Q.view_content_sha256(world.views / v),
                                      "target_set_sha256": Q.target_set_sha256(world.views / v)} for v in Q.FROZEN}))
    z = tmp_path / "z.json"
    z.write_text(json.dumps({"train_sha256": Q.file_sha256(world.store / "lidc/train.parquet"), "series": {},
                             "step2": []}))
    out = tmp_path / "rows"

    def go(*extra, steps=PLAN["steps"]):
        return R.main(["--store", str(world.store), "--out", str(out), "--steps", str(steps),
                       "--per-step", str(PLAN["per_step"]), "--block-steps", str(PLAN["block_steps"]),
                       "--draws-options", "1", "--shares", *map(str, SHARES.values()), "--seed", str(PLAN["seed"]),
                       "--disagreement-weight", "1.0", "--procs", "1", "--exclude", str(ex),
                       "--eval-views", str(world.views), "--review", str(tmp_path / "none.json"), "--frozen",
                       str(frozen), "--series-z", str(z), "--audit-out", str(tmp_path / "png"), *extra])
    return types.SimpleNamespace(go=go, out=out, z=z, R=R)


@pytest.fixture
def built(run):
    assert run.go() == 0
    return run


def test_manifest_and_rows_follow_the_schedule(built, world):  # noqa: F811
    m = json.loads((built.out / "manifest.json").read_text())
    s = m["schedule"]
    recs = Q.records(built.out)
    want = [w.window_id for w in world.sched]  # the same pool and schedule arguments as the fixture's
    assert [r["item_id"] for r in recs] == want
    assert s["ids_sha256"] == built.R._sha("\n".join(want).encode())
    assert {k: s[k] for k in ("steps", "planned_steps", "block_steps", "per_step", "seed", "draws")} == {
        "steps": 13, "planned_steps": 13, "block_steps": 13, "per_step": 2, "seed": 4, "draws": 1}
    for k in ("shares", "target_law", "target_sha256", "dropped", "disagreement_weight", "disagreement_max_reuse",
              "exclusions_sha256", "excluded_ids", "gates", "files_count", "files_sha256"):
        assert k in s
    assert all(g["ok"] for g in s["gates"] if g["hard"])
    assert (s["files_count"], s["files_sha256"]) == built.R.files_digest(built.out)
    assert m["files"]["lidc/train.jsonl"]["sha256"] == built.R._sha((built.out / "lidc/train.jsonl").read_bytes())


def test_each_window_has_its_own_png_and_caches_and_they_load(built):
    from objdel.verify.reward.cache import ExpertCache
    recs = Q.records(built.out)
    uniq = {r["item_id"]: r for r in recs}
    pngs = {built.out / "images" / f"{i}.png" for i in uniq}
    assert len(pngs) == len(uniq) and all(p.is_file() for p in pngs)
    for i, r in uniq.items():
        im = cv2.imread(str(built.out / "images" / f"{i}.png"), cv2.IMREAD_COLOR)
        assert pixels_sha256(cv2.cvtColor(im, cv2.COLOR_BGR2RGB)) == r["image_sha256"]
        assert r["view"]["kind"] in Wn.KINDS and r["view"]["op"] == "id" and "~" in i
        for name, digest in r["caches"].items():
            ExpertCache.load(built.out / "caches" / i / f"{name}.npz", digest)
    assert sum(r["rl"] for r in uniq.values()) == len(uniq)  # every window is of an RL item


def test_cache_loader_and_grpo_rows_read_the_root_in_schedule_order(built, world):  # noqa: F811
    from objdel.rl import rows as D
    from objdel.verify import store as verify_store
    rows = D.grpo_rows(built.out, "lidc", "train", oversample=1)
    assert list(rows["key"]) == [w.window_id for w in world.sched]
    load = verify_store.cache_loader(built.out)
    assert load(rows["key"][0])


def test_no_manifest_when_a_hard_gate_fails(run, world):  # noqa: F811
    rec_item = world.sched[0].item_id
    s = rec_item.split("/")[1]
    run.z.write_text(json.dumps({"train_sha256": Q.file_sha256(world.store / "lidc/train.parquet"), "step2": [s],
                                 "series": {s: {"split": "train", "items": 1, "item_ids": [rec_item]}}}))
    assert run.go() == 1
    assert (run.out / "lidc/train.jsonl").is_file() and not (run.out / "manifest.json").exists()


def test_steps_not_a_multiple_of_block_steps_is_refused(run):
    with pytest.raises(SystemExit, match="multiple of block_steps"):
        run.R.main(["--out", str(run.out), "--steps", "30", "--exclude", "{}"])
    assert not run.out.exists()


def test_exclude_is_required(run):
    with pytest.raises(SystemExit):
        run.R.parse_args(["--out", str(run.out), "--steps", "20"])


def test_choose_prefers_no_reuse_then_smallest_passing():
    R = __import__("experiments.phase2.build_rl_rows", fromlist=["choose"])
    t = {1: {"ok": False, "max_reuse": None}, 2: {"ok": True, "max_reuse": 3}, 4: {"ok": True, "max_reuse": 1},
         8: {"ok": True, "max_reuse": 1}}
    assert R.choose(t) == 4
    assert R.choose({2: {"ok": True, "max_reuse": 3}, 4: {"ok": True, "max_reuse": 2}}) == 2
    assert R.choose({1: {"ok": False, "max_reuse": None}}) is None


def test_hash_gate_rebuilds_a_rows_root(built, world):  # noqa: F811
    """verify_target_hashes follows source_item_id: a rows root keys records by window id, several per store item."""
    from experiments.phase2 import verify_target_hashes as V
    assert V.main([str(built.out / "manifest.json"), "--store", str(world.store), "--procs", "1"]) == 0
    m = json.loads((built.out / "manifest.json").read_text())
    m["sets"]["lidc/train"]["target_set_sha256"] = "0" * 64
    (built.out / "manifest.json").write_text(json.dumps(m))
    assert V.main([str(built.out / "manifest.json"), "--store", str(world.store), "--procs", "1"]) == 1

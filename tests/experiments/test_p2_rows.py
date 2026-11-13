"""objdel.cli.data.rl_rows (Implementation.md B1.4) on the test_data_quality fixture store: 30 synthetic train
items, the val and test views, one 13-step block of 2 rows (shares 0.5 / 0.5 / 0, seed 4, draws 1, weight 1), so the
off-centre rows are the train shifted law and every hard gate passes on a clean build."""
import io
import json
import re
import shutil
import sys
import types
from dataclasses import replace
from pathlib import Path

import cv2
import pytest
from test_data_quality import PLAN, SHARES, all_items, world, write_store  # noqa: F401

from objdel.core.hashing import bytes_sha256, pixels_sha256
from objdel.data import quality as Q
from objdel.data import windows as Wn
from tests.support.quality import quality_stubs


@pytest.fixture
def run(world, tmp_path, monkeypatch):  # noqa: F811
    from objdel.cli.data import rl_rows as R
    stubs = quality_stubs(world, tmp_path, monkeypatch)
    ex, frozen, z = stubs.ex, stubs.frozen, stubs.z
    out = tmp_path / "rows"

    def go(*extra, steps=PLAN["steps"], store=world.store):
        return R.main(["--store", str(store), "--out", str(out), "--steps", str(steps),
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
    assert s["ids_sha256"] == bytes_sha256("\n".join(want).encode())
    assert {k: s[k] for k in ("steps", "planned_steps", "block_steps", "per_step", "seed", "draws")} == {
        "steps": 13, "planned_steps": 13, "block_steps": 13, "per_step": 2, "seed": 4, "draws": 1}
    for k in ("shares", "target_law", "target_sha256", "dropped", "disagreement_weight", "disagreement_max_reuse",
              "exclusions_sha256", "excluded_ids", "gates", "files_count", "files_sha256"):
        assert k in s
    assert all(g["ok"] for g in s["gates"] if g["hard"])
    assert (s["files_count"], s["files_sha256"]) == built.R.files_digest(built.out)
    assert m["files"]["lidc/train.jsonl"]["sha256"] == bytes_sha256((built.out / "lidc/train.jsonl").read_bytes())


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
    R = __import__("objdel.cli.data.rl_rows", fromlist=["choose"])
    t = {1: {"ok": False, "max_reuse": None}, 2: {"ok": True, "max_reuse": 3}, 4: {"ok": True, "max_reuse": 1},
         8: {"ok": True, "max_reuse": 1}}
    assert R.choose(t) == 4
    assert R.choose({2: {"ok": True, "max_reuse": 3}, 4: {"ok": True, "max_reuse": 2}}) == 2
    assert R.choose({1: {"ok": False, "max_reuse": None}}) is None


def test_hash_gate_rebuilds_a_rows_root(built, world):  # noqa: F811
    """verify_target_hashes follows source_item_id: a rows root keys records by window id, several per store item."""
    from tools.data import verify_target_hashes as V
    assert V.main([str(built.out / "manifest.json"), "--store", str(world.store), "--procs", "1"]) == 0
    m = json.loads((built.out / "manifest.json").read_text())
    m["sets"]["lidc/train"]["target_set_sha256"] = "0" * 64
    (built.out / "manifest.json").write_text(json.dumps(m))
    assert V.main([str(built.out / "manifest.json"), "--store", str(world.store), "--procs", "1"]) == 1


def test_e02_cross_split_patient_fails_dg1_then_a_clean_rerun_writes_the_manifest(run, tmp_path, capsys):
    """E02: a val item of a train patient fails DG1 (exit 1, no manifest, stderr names DG1); the re-run into the same
    --out on the clean store writes the manifest and keeps nothing of the failed run."""
    train, other = all_items()
    other[0] = replace(other[0], group_id=train[0].group_id)  # val item 0 belongs to train patient 0
    assert run.go(store=write_store(tmp_path / "bad_store", train + other)) == 1
    assert re.search(r"\bDG1\b", capsys.readouterr().err) and not (run.out / "manifest.json").exists()
    (run.out / "images" / "stale.png").write_bytes(b"")
    assert run.go() == 0 and (run.out / "manifest.json").is_file()
    assert not (run.out / "images" / "stale.png").exists()


def test_e02_a_finished_root_or_a_foreign_folder_is_never_touched(built, tmp_path):
    with pytest.raises(SystemExit, match="manifest.json"):
        built.go()
    foreign = tmp_path / "foreign"
    (foreign / "notes").mkdir(parents=True)
    with pytest.raises(SystemExit, match="notes"):
        built.R.main(["--out", str(foreign), "--steps", "20", "--exclude", "{}"])
    assert (foreign / "notes").is_dir()


def test_i20_e02_a_pool_too_small_for_every_draws_value_names_the_last_kind_and_cell(run, capsys):
    """I20, E02: every draws value fails (the pool is too small for --steps): non-zero exit, no manifest, and the
    message names the last value's failing kind and cell."""
    assert run.go("--draws-options", "1", "2", steps=13 * 40) != 0
    err = capsys.readouterr().err
    assert re.search(r"draws 2: no (free )?(centred|shifted|rescaled) window in cell \(?\d", err), err
    assert not (run.out / "manifest.json").exists()


class Hub:
    """The HfApi calls hf.push_rows makes, recorded; nothing leaves the machine."""

    def __init__(self, private=True):
        self.private, self.calls, self.remote, self.blobs = private, [], [], {}

    def repo_info(self, repo, repo_type):
        return types.SimpleNamespace(private=self.private)

    def create_branch(self, repo, branch, repo_type, exist_ok):
        self.calls.append(("create_branch", repo, branch))

    def upload_file(self, *, path_or_fileobj, path_in_repo, repo_id, revision, repo_type, commit_message):
        self.calls.append(("upload_file", repo_id, path_in_repo, revision))
        self.remote.append(path_in_repo)
        self.blobs[path_in_repo] = Path(path_or_fileobj).read_bytes()  # the temp tar is gone after push_rows
        return types.SimpleNamespace(oid=("b" if path_in_repo.endswith("manifest.json") else "a") * 40)


TOKEN = "hf_secretTokenValue123"


def push(monkeypatch, root, hub):
    from objdel.cli.hub import push as hf_push
    from objdel.infra import hf
    monkeypatch.setattr(hf, "api", lambda token: hub if token == TOKEN else None)
    monkeypatch.setattr(sys, "stdin", io.StringIO(TOKEN + "\n"))
    return hf_push.main(["--rows", str(root), "--repo", hf.REPO, "--revision", "main"])


def test_i21_rows_upload_goes_to_the_one_place_with_the_manifest_last(built, tmp_path, monkeypatch, capsys):
    """I21, B1.4 step 5: the root without its manifest as one tar (hf.rows_tar: the free repo holds <= 20,000 files),
    then the manifest; the printed rows_revision is the last commit (it holds both) and rows_manifest_sha256 the
    manifest's SHA-256. The tar holds every other file of the root, byte for byte, and is the same bytes again."""
    import tarfile

    from objdel.core.hashing import file_sha256
    from objdel.infra import hf
    root = shutil.copytree(built.out, tmp_path / "lidc-train-b1")
    hub = Hub()
    assert push(monkeypatch, root, hub) == 0
    assert [c[0] for c in hub.calls] == ["create_branch", "upload_file", "upload_file"]
    assert hub.remote == ["rows/lidc-train-b1/root.tar", "rows/lidc-train-b1/manifest.json"]
    assert all(c[1] == hf.REPO and c[-1] == "main" for c in hub.calls)
    blob = hub.blobs["rows/lidc-train-b1/root.tar"]
    with tarfile.open(fileobj=io.BytesIO(blob)) as tf:
        packed = {m.name: tf.extractfile(m).read() for m in tf.getmembers()}
    local = {p.relative_to(root).as_posix(): p.read_bytes() for p in root.rglob("*") if p.is_file()}
    assert packed == {k: v for k, v in local.items() if k != "manifest.json"} and len(packed) > 1
    again = io.BytesIO()
    hf.rows_tar(root, again)
    assert again.getvalue() == blob
    out, err = capsys.readouterr()
    assert f"rows_revision: {'b' * 40}" in out and f"rows_manifest_sha256: {file_sha256(root / 'manifest.json')}" in out
    assert TOKEN not in out + err


def test_rows_tar_and_files_digest_keep_the_posix_string_order(tmp_path):
    """Pinned hashes: a tree with `a-b/x` and `a/x` (string order and Path order differ) packs and digests in the
    order the pre-files_under code gave: sorted by the relative posix string, manifest.json left out."""
    import tarfile

    from objdel.cli.data import rl_rows
    from objdel.infra import hf
    root = tmp_path / "root"
    for rel in ("a/x", "a-b/x", "manifest.json", "lidc/train.jsonl"):
        (root / rel).parent.mkdir(parents=True, exist_ok=True)
        (root / rel).write_text(rel)
    old = sorted(p.relative_to(root).as_posix() for p in root.rglob("*") if p.is_file())
    old = [r for r in old if r != "manifest.json"]
    assert old == ["a-b/x", "a/x", "lidc/train.jsonl"]
    assert rl_rows.files_digest(root) == (3, bytes_sha256("\n".join(old).encode()))
    buf = io.BytesIO()
    hf.rows_tar(root, buf)
    with tarfile.open(fileobj=io.BytesIO(buf.getvalue())) as tf:
        assert tf.getnames() == old


def test_i21_nothing_is_uploaded_after_a_hard_gate_failure(run, world, monkeypatch, capsys):  # noqa: F811
    item = world.sched[0].item_id
    s = item.split("/")[1]
    run.z.write_text(json.dumps({"train_sha256": Q.file_sha256(world.store / "lidc/train.parquet"), "step2": [s],
                                 "series": {s: {"split": "train", "items": 1, "item_ids": [item]}}}))
    assert run.go() == 1
    hub = Hub()
    with pytest.raises(SystemExit, match="manifest.json"):
        push(monkeypatch, run.out, hub)
    assert hub.calls == [] and TOKEN not in "".join(capsys.readouterr())


def test_i21_nothing_is_uploaded_to_a_public_repo(built, monkeypatch, capsys):
    hub = Hub(private=False)
    with pytest.raises(SystemExit, match="not private"):
        push(monkeypatch, built.out, hub)
    assert hub.calls == [] and TOKEN not in "".join(capsys.readouterr())

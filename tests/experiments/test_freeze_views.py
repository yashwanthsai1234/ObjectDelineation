"""B1.6: reference/frozen_views.json (experiments.phase2.freeze_views): one entry per frozen view with its content
hash, target-set hash, build arguments, generator commit, render files and known label issues (I05, I38)."""
import hashlib
import json
import subprocess
from pathlib import Path

import pytest

from experiments.phase2 import freeze_views as FV
from objdel.core.io import code_sha, git_sha
from objdel.data import quality as Q
from objdel.data import store

ROOT = Path(__file__).resolve().parents[2]
RECORD = ROOT / "reference/frozen_views.json"
KEYS = {"view", "items", "content_sha256", "target_set_sha256", "seed", "kinds", "store_manifest_sha256",
        "generator_commit", "render_files", "known_label_issues", "fallback", "padded", "r5_code_commit"}


@pytest.fixture(scope="module")
def view(tmp_path_factory, synthetic_items):
    pytest.importorskip("datasets")
    from experiments.phase2 import build_eval_views as BE
    root = tmp_path_factory.mktemp("frozen")
    store.write(synthetic_items(("train", "val")), root / "store")
    BE.build(root / "store", root / "out", ("lidc",), "val", procs=1, kinds=("shifted", "edge"))
    return root / "out" / "lidc-val-edge"


def _content(d):  # the plan's rule, written out again: newline-joined (item_id, offset, side, image_sha256)
    recs = [json.loads(x) for x in (d / "lidc/val.jsonl").read_text().splitlines()]
    return hashlib.sha256("\n".join(json.dumps([r["item_id"], r["view"]["offset"], r["view"]["side"],
                                                 r["image_sha256"]]) for r in recs).encode()).hexdigest()


def test_entry_holds_the_content_hash_and_the_build_record(view):
    m = json.loads((view / "manifest.json").read_text())
    ids = [json.loads(x)["item_id"] for x in (view / "lidc/val.jsonl").read_text().splitlines()]
    e = FV.entry(view, known=[ids[0], "lidc/9.9/n000/z001"])
    assert set(e) == KEYS and e["view"] == "lidc-val-edge" and e["items"] == len(ids)
    assert e["content_sha256"] == _content(view) != hashlib.sha256((view / "manifest.json").read_bytes()).hexdigest()
    assert e["target_set_sha256"] == m["sets"]["lidc/val"]["target_set_sha256"]
    assert (e["seed"], e["kinds"], e["store_manifest_sha256"]) == (0, ["shifted", "edge"], m["store_manifest_sha256"])
    assert e["known_label_issues"] == [ids[0]]  # only the view's own items
    assert set(e["render_files"]) == set(FV.RENDER_FILES)
    assert e["generator_commit"] == git_sha(ROOT)  # X2-3: HEAD, -dirty-<sha12> while the worktree holds edits
    assert e["r5_code_commit"] == code_sha(ROOT / "runs/r5", required=False)


def test_render_file_drift_is_reported():
    """I38: each listed render module carries its last commit; a changed commit is reported by path."""
    files = FV.render_files()
    assert set(files) == set(FV.RENDER_FILES) and all(len(c) == 40 for c in files.values())
    assert FV.drift({"render_files": files}) == []
    planted = files | {"src/objdel/data/targets/frame.py": "0" * 40}
    assert FV.drift({"render_files": planted}) == ["src/objdel/data/targets/frame.py"]


def test_freeze_refuses_uncommitted_render_files(tmp_path):
    """A commit id beside an edited file would not name the code that built the view."""
    for f in FV.RENDER_FILES:
        (tmp_path / f).parent.mkdir(parents=True, exist_ok=True)
        (tmp_path / f).write_text("x = 1\n")
    git = ["git", "-c", "user.name=t", "-c", "user.email=t@t"]
    subprocess.run([*git, "init", "-q"], cwd=tmp_path, check=True)
    subprocess.run([*git, "add", "."], cwd=tmp_path, check=True)
    subprocess.run([*git, "commit", "-qm", "c"], cwd=tmp_path, check=True)
    assert FV.dirty(tmp_path) == []
    (tmp_path / FV.RENDER_FILES[0]).write_text("x = 2\n")
    assert FV.dirty(tmp_path) == [FV.RENDER_FILES[0]]


def test_cli_merges_entries_into_the_record(view, tmp_path, monkeypatch):
    out = tmp_path / "frozen_views.json"
    out.write_text(json.dumps({"lidc-test-shifted": {"content_sha256": "a"}}))
    monkeypatch.setattr(FV, "dirty", lambda root=ROOT: [])  # the worktree may hold edits while tests run
    FV.main([str(view), str(view.parent / "lidc-val-shifted"), "--out", str(out), "--known", str(tmp_path / "none")])
    got = json.loads(out.read_text())
    assert sorted(got) == ["lidc-test-shifted", "lidc-val-edge", "lidc-val-shifted"]
    assert got["lidc-test-shifted"] == {"content_sha256": "a"} and got["lidc-val-edge"]["known_label_issues"] == []


def test_the_record_has_every_frozen_view():
    """I05 on the committed record: one entry per frozen view with every key, the test view's R5 identity note."""
    got = json.loads(RECORD.read_text())
    assert set(Q.FROZEN) <= set(got)
    for v in Q.FROZEN:
        assert KEYS <= set(got[v]) and got[v]["view"] == v
    assert "not measured" in got["lidc-test-shifted"]["r5_identity"]
    assert got["lidc-test-shifted"]["known_label_issues"]  # LIDC-IDRI-0659 (B1.7, owner default: keep)


@pytest.mark.data
def test_the_record_matches_the_built_views():
    """DG15 on the Mac views the record was frozen from (data/eval_views)."""
    views = ROOT / "data/eval_views"
    if not views.is_dir():
        pytest.skip(f"{views} not found")
    assert Q.frozen(views, json.loads(RECORD.read_text())) == {}


def _pre(view, tmp_path):
    """B1.0 step 1: final_eval --hash-view of the val shifted view (the pre-rebuild record)."""
    from experiments.phase2 import final_eval as CLI
    pre = tmp_path / "val_pre_rebuild.json"
    CLI.main(["--hash-view", str(view), "--out", str(pre)])
    return pre


def test_i36_r5_stamp_check_goes_into_the_pre_rebuild_record(view, tmp_path, monkeypatch):
    """I36: the R5 stamp's data hash against preflight.hashes recomputed, written beside the view hashes; a planted
    manifest change makes it "not equal"."""
    from experiments.phase2 import val_rebuild as VR
    from objdel.infra import preflight
    shifted = view.parent / "lidc-val-shifted"
    data = lambda cfg: {"data": Q.file_sha256(shifted / "manifest.json")}  # noqa: E731  the stamp covers the manifest
    monkeypatch.setattr(preflight, "hashes", data)
    stamp = tmp_path / "r1.json"
    stamp.write_text(json.dumps({"config": "experiments/phase2/configs/r1_qwen3vl2b.yaml",
                                 "hashes": {"data": data(None)["data"]}}))
    pre = _pre(shifted, tmp_path)
    VR.main(["stamp", "--pre", str(pre), "--stamp", str(stamp)])
    doc = json.loads(pre.read_text())
    assert doc["content_sha256"] == Q.view_content_sha256(shifted)
    assert doc["r5_stamp_check"] == {"stamp_file": str(stamp), "stamp_data_hash": data(None)["data"],
                                     "recomputed_data_hash": data(None)["data"], "equal": True}
    m = shifted / "manifest.json"
    keep = m.read_bytes()
    try:
        m.write_bytes(keep + b"\n")
        VR.main(["stamp", "--pre", str(pre), "--stamp", str(stamp)])
        assert json.loads(pre.read_text())["r5_stamp_check"]["equal"] is False
    finally:
        m.write_bytes(keep)


def _copy_view(view, to):
    import shutil
    shutil.copytree(view, to)
    return to


def test_i36_rebuild_replaces_the_folder_only_on_equal_hashes(view, tmp_path, monkeypatch):
    """I36 (B1.0 step 2): an equal rebuild replaces the folder; a planted difference keeps both folders, freezes the
    rebuilt view under the view's name, records both hashes and marks A3-val the only Gate B1 baseline, so the R5
    hash chain (rungates.hash_chain_ok) fails."""
    from experiments.phase2 import val_rebuild as VR
    from objdel.verify.rungates import hash_chain_ok
    monkeypatch.setattr(FV, "dirty", lambda root=ROOT: [])
    old = _copy_view(view.parent / "lidc-val-shifted", tmp_path / "views" / "lidc-val-shifted")
    pre = _pre(old, tmp_path)
    doc = json.loads(pre.read_text())
    pre.write_text(json.dumps(doc | {"r5_stamp_check": {"equal": True}}))
    frozen = tmp_path / "frozen_views.json"
    args = ["swap", "--pre", str(pre), "--frozen", str(frozen), "--known", str(tmp_path / "none")]

    new = _copy_view(old, tmp_path / "views" / "rebuilt")
    VR.main([*args, "--rebuilt", str(new)])
    assert not new.exists() and sorted(p.name for p in old.parent.iterdir()) == ["lidc-val-shifted"]
    e = json.loads(frozen.read_text())["lidc-val-shifted"]
    assert e["content_sha256"] == doc["content_sha256"] and e["r5_baseline"] == "hash chain"
    assert hash_chain_ok(json.loads(pre.read_text()), e)

    new = _copy_view(old, tmp_path / "views" / "rebuilt")
    recs = [json.loads(x) for x in (new / "lidc/val.jsonl").read_text().splitlines()]
    recs[0]["image_sha256"] = "f" * 64  # a planted rendering difference
    (new / "lidc/val.jsonl").write_text("".join(json.dumps(r) + "\n" for r in recs))
    VR.main([*args, "--rebuilt", str(new)])
    kept = old.parent / "lidc-val-shifted-pre-rebuild"
    assert Q.view_content_sha256(kept) == doc["content_sha256"] != Q.view_content_sha256(old)
    assert not new.exists()
    e = json.loads(frozen.read_text())["lidc-val-shifted"]
    assert e["content_sha256"] == Q.view_content_sha256(old)
    assert e["pre_rebuild"] == {"folder": str(kept), "content_sha256": doc["content_sha256"],
                                "target_set_sha256": doc["target_set_sha256"]}
    assert e["r5_baseline"] == "report line only: A3-val is the only Gate B1 baseline (B1.0 step 2)"
    assert not hash_chain_ok(json.loads(pre.read_text()), e)

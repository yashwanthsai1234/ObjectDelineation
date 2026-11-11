"""B1.0: the cache fingerprint covers only what stored cache values depend on (U01, I01); val halves (U02)."""
import dataclasses
import json
import pathlib
import subprocess
import sys
from collections import Counter

import cv2
import numpy as np
import pytest
import scipy

from objdel.data import splits
from objdel.data.splits import val_half
from objdel.verify.reward import spec as S
from objdel.verify.reward.cache import ExpertCache, StaleCacheError
from objdel.verify.reward.spec import fingerprint

ROOT = pathlib.Path(__file__).resolve().parents[2]


@pytest.mark.parametrize("field,val", [("weights", (1.0, 0.0, 0.0, 0.0)), ("floors", (0.1,) * 4), ("c_topo", 5.0),
                                       ("lam", 0.9), ("eps", 0.2), ("group", 4), ("drop_spread", 0.05),
                                       ("iou_gate", 0.6)])
def test_reward_weights_do_not_change_fingerprint(monkeypatch, field, val):
    before = fingerprint()
    monkeypatch.setattr(S, "REWARD_V4", dataclasses.replace(S.REWARD_V4, **{field: val}))
    assert fingerprint() == before


def test_cache_inputs_change_fingerprint(monkeypatch):
    """7 inputs, 7 distinct digests besides the base one: the three CacheSpec constants, the size rule, numpy,
    scipy and OpenCV versions."""
    base = fingerprint()
    assert len(base) == 64 and int(base, 16) >= 0
    seen = {base}
    changes = [(S, "CACHE_SPEC", dataclasses.replace(S.CACHE_SPEC, expert_score=0.9)),
               (S, "CACHE_SPEC", dataclasses.replace(S.CACHE_SPEC, tau_floor=2.0)),
               (S, "CACHE_SPEC", dataclasses.replace(S.CACHE_SPEC, tau_quantile=90.0)),
               (S, "min_area", lambda a: 7), (np, "__version__", "0.0.0"), (scipy, "__version__", "0.0.0"),
               (cv2, "__version__", "0.0.0")]
    for mod, name, val in changes:
        with monkeypatch.context() as m:
            m.setattr(mod, name, val)
            seen.add(fingerprint())
    assert len(seen) == 8 and fingerprint() == base


def disk(n, r):
    y, x = np.mgrid[:n, :n]
    return np.hypot(y - (n / 2 - 0.5), x - (n / 2 - 0.5)) <= r


OLD_FINGERPRINT = "b1eb9f51824fe9617341a34e0aaa470a4b75646e4fad035679c4ff772d96c81f"  # data/prepared, before B1.0


def test_stored_values_survive_the_new_fingerprint(tmp_path):
    """I01: tau and T5 of the three-disk fixture are bit-equal to the values built before B1.0; only the digest
    changes, and a cache saved under the old fingerprint is stale until rebuilt."""
    c = ExpertCache.build(np.stack([disk(60, 18), disk(60, 19), disk(60, 20)]))
    assert (c.tau.hex(), c.t5.hex()) == ("0x1.1e377a0000000p+1", "0x1.8412c21ddf4f2p+1")  # captured before B1.0
    assert fingerprint() != OLD_FINGERPRINT
    assert c.digest != "5076f66ee82bea6f9a06b15011539affae8e732181d3917b2ac09e51ca672024"  # the pre-B1.0 digest
    from objdel.verify.reward import cache as cache_mod
    old = cache_mod._digest(c.raters, c.tau, c.t5, OLD_FINGERPRINT)
    np.savez_compressed(tmp_path / "old.npz", raters=c.raters, tau=c.tau, t5=c.t5, fingerprint=OLD_FINGERPRINT,
                        sha256=old)
    with pytest.raises(StaleCacheError, match="old.npz"):
        ExpertCache.load(tmp_path / "old.npz", expected_digest=old)
    d = ExpertCache.load(tmp_path / "old.npz", expected_digest=old, rebuild=True)
    assert (d.tau, d.t5, d.digest) == (c.tau, c.t5, c.digest)


# ------------------------------------------------------------------ val halves (U02)

HALVES = ("shape", 0.5), ("select", 0.5)


def _halves(records, s):
    """Test-local B1.0 search step from the existing primitives (deal, lidc_stratum)."""
    items, dis = Counter(r["group_id"] for r in records), Counter(r["group_id"] for r in records if r["disagreement"])
    strata = {}
    for g in items:  # B1.0 Patient form: (items, disagreement items) in the val view
        strata.setdefault(splits.lidc_stratum(splits.Patient(items[g], dis[g])), []).append(g)
    out = {}
    for gids in strata.values():
        out |= splits.deal(gids, HALVES, seed=f"val_half:{s}")
    n, d = Counter(), Counter()
    for g, h in out.items():
        n[h] += items[g]
        d[h] += dis[g]
    ok = abs(n["shape"] - n["select"]) <= 0.04 * len(records) and abs(d["shape"] - d["select"]) <= 4
    return out, ok


@pytest.fixture
def val_records():
    """Synthetic val view: 40 patients, one with 60 items, so seed "val_half:0" breaks the item balance."""
    recs = []
    for p in range(40):
        n = 60 if p == 0 else 2 + p % 9
        recs += [{"group_id": f"LIDC-IDRI-{p:04d}", "disagreement": p % 7 == 0 and i < 3, "iou_tc": 0.5}
                 for i in range(n)]
    return recs


def test_val_half_table_is_the_first_balanced_seed(val_records):
    assert not _halves(val_records, 0)[1]  # the fixture is built so that s = 0 fails
    t = splits.val_half_table(val_records)
    assert t["s"] >= 1 and t["patient_form"] == "items,disagreement_items"
    assert not any(_halves(val_records, s)[1] for s in range(t["s"]))  # every smaller s fails (a) or (b)
    want, ok = _halves(val_records, t["s"])
    assert ok and t["halves"] == want
    n = Counter()
    for r in val_records:
        n[want[r["group_id"]], r["disagreement"]] += 1
    assert t["items"] == {h: n[h, False] + n[h, True] for h in ("shape", "select")}
    assert t["disagreement_items"] == {h: n[h, True] for h in ("shape", "select")}


def test_val_half_table_reads_no_score(val_records):
    planted = [r | {"iou_tc": float(i % 3)} for i, r in enumerate(val_records)]
    assert splits.val_half_table(planted) == splits.val_half_table(val_records)


def test_val_half_table_raises_when_no_seed_passes(val_records):
    with pytest.raises(ValueError, match="no seed"):
        splits.val_half_table(val_records, max_s=1)


def test_val_half_reads_the_committed_table():
    t = json.loads((ROOT / "reference/val_halves.json").read_text())
    assert {g: val_half(g) for g in t["halves"]} == t["halves"]
    assert set(t["halves"].values()) == {"shape", "select"}
    assert t["patient_form"] == "items,disagreement_items" and len(t["r5_selection_items"]) == 128
    total = sum(t["items"].values())
    assert abs(t["items"]["shape"] - t["items"]["select"]) <= 0.04 * total
    assert abs(t["disagreement_items"]["shape"] - t["disagreement_items"]["select"]) <= 4
    for bad in ("LIDC-IDRI-0001", "lidc/LIDC-IDRI-0026"):  # a train patient; a prefixed id
        with pytest.raises(KeyError):
            val_half(bad)


def test_val_half_table_matches_the_val_view():
    """The committed table is the seed search on today's val shifted view (skip without data)."""
    path = ROOT / "data/eval_views/lidc-val-shifted/lidc/val.jsonl"
    if not path.is_file():
        pytest.skip("no local val view")
    recs = [json.loads(x) for x in path.read_text().splitlines()]
    t = json.loads((ROOT / "reference/val_halves.json").read_text())
    table = splits.val_half_table(recs)
    assert table == {k: t[k] for k in table} and table["s"] == 78  # s = 78: Implementation.md B1.0 [measured]


def test_val_half_independent_of_hashseed():
    want = json.loads((ROOT / "reference/val_halves.json").read_text())["halves"]["LIDC-IDRI-0026"]
    code = "from objdel.data.splits import val_half;print(val_half('LIDC-IDRI-0026'))"  # a val id, no prefix
    outs = {subprocess.run([sys.executable, "-I", "-c", code], capture_output=True, text=True, cwd=ROOT,
                           env={"PYTHONHASHSEED": s}).stdout for s in ("0", "1", "2")}
    assert outs == {want + "\n"}  # a KeyError prints nothing and fails here


def test_val_halves_script_writes_the_table(tmp_path, val_records):
    """experiments.phase2.val_halves: the table plus the sorted held-out ids that chose R5; an id outside the view
    is refused."""
    from experiments.phase2 import val_halves as VH
    recs = [r | {"item_id": f"lidc/{r['group_id']}/z{i:03d}"} for i, r in enumerate(val_records)]
    view = tmp_path / "view.jsonl"
    view.write_text("".join(json.dumps(r) + "\n" for r in recs))
    held = tmp_path / "heldout.json"
    held.write_text(json.dumps({"per_item_iou": {recs[5]["item_id"]: 0.7, recs[2]["item_id"]: 0.1}}))
    VH.main(["--view", str(view), "--heldout", str(held), "--out", str(tmp_path / "t.json")])
    t = json.loads((tmp_path / "t.json").read_text())
    assert t == splits.val_half_table(val_records) | {"r5_selection_items": [recs[2]["item_id"], recs[5]["item_id"]]}
    held.write_text(json.dumps({"per_item_iou": {"lidc/elsewhere/z000": 0.5}}))
    with pytest.raises(SystemExit, match="not in the val view"):
        VH.main(["--view", str(view), "--heldout", str(held), "--out", str(tmp_path / "u.json")])


def test_r5_view_stamp_chain(tmp_path, monkeypatch):
    """I46: the launch stamp's data hash pins the val shifted view: unchanged folder -> equal to the planted stamp;
    one byte of lidc/val.jsonl changed (manifest's file sha updated with it) -> different."""
    import hashlib

    from objdel.infra import preflight as PF
    data = tmp_path / "data"
    for d in ("views/v", "checker_eval", *(f"eval_views/{v}" for v in PF.EXAMS)):
        (data / d).mkdir(parents=True)
        (data / d / "manifest.json").write_text(json.dumps({"view_hash": "h"}))
    val = data / "eval_views/lidc-val-shifted"
    (val / "lidc").mkdir()

    def write(text):
        (val / "lidc/val.jsonl").write_text(text)
        sha = hashlib.sha256(text.encode()).hexdigest()
        (val / "manifest.json").write_text(json.dumps({"files": {"lidc/val.jsonl": {"sha256": sha}}}))

    write('{"item_id": "a", "view": {"offset": [3, 4]}}\n')
    cfg = tmp_path / "c.yaml"
    cfg.write_text("schema_version: 1\nname: c\ncompute: {gpu: RTX_5090, cap_hours: 1.0}\nbundle_view: views/v\n"
                   "stages:\n- {kind: cmd, name: box, argv: [bash, -c, 'true']}\n")
    monkeypatch.setenv("OBJDEL_DATA_ROOT", str(data))
    stamp = PF.hashes(cfg)["data"]
    assert PF.hashes(cfg)["data"] == stamp
    write('{"item_id": "a", "view": {"offset": [3, 5]}}\n')
    assert PF.hashes(cfg)["data"] != stamp

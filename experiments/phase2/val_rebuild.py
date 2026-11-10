"""B1.0 val view identity (Implementation.md B1.0 steps 1-2): the R5 stamp check beside the pre-rebuild hash record,
then the folder swap after the val rebuild.

    .venv/bin/python -I -m experiments.phase2.final_eval --hash-view data/eval_views/lidc-val-shifted \
        --out runs/a3/val_pre_rebuild.json
    .venv/bin/python -I -m experiments.phase2.val_rebuild stamp --pre runs/a3/val_pre_rebuild.json \
        --stamp runs/preflight/r1_qwen3vl2b.json
    (rebuild val into a new folder with build_eval_views)
    .venv/bin/python -I -m experiments.phase2.val_rebuild swap --pre runs/a3/val_pre_rebuild.json --rebuilt NEW

stamp: r5_stamp_check = the R5 launch stamp's data hash, objdel.infra.preflight.hashes of the stamp's config
recomputed now (it covers the view's manifest.json, which pins val.jsonl and target_set_sha256), and whether they are
equal. Gate B1 reads it for the R5 fallback (objdel.verify.rungates.hash_chain_ok).
swap: the rebuilt view replaces the pre-rebuild folder only when its content hash and target-set hash equal the
record's. On a difference both folders stay (the old one as <name>-pre-rebuild, the rebuilt one under the view's
name), and A3-val is the only Gate B1 baseline: R5's items.jsonl is a report line. Either way the view at the old path
is frozen (freeze_views; DG15) with the pre-rebuild hashes and the baseline rule in its entry.
Reads untrusted data: run with python -I.
"""
import argparse
import shutil
from pathlib import Path

from experiments.phase2 import freeze_views as FV
from objdel.core.io import read_json, write_json
from objdel.data import quality as Q

ROOT = FV.ROOT
KEYS = ("content_sha256", "target_set_sha256")


def stamp(pre: Path, stamp_file: Path) -> None:
    from objdel.infra import preflight

    s = read_json(stamp_file)
    want, got = s["hashes"]["data"], preflight.hashes(ROOT / s["config"])["data"]
    write_json(pre, read_json(pre) | {"r5_stamp_check": {"stamp_file": str(stamp_file), "stamp_data_hash": want,
                                                         "recomputed_data_hash": got, "equal": want == got}})


def swap(pre: Path, rebuilt: Path, frozen: Path, known: Path) -> None:
    rec = read_json(pre)
    old = Path(rec["view"])
    same = (Q.view_content_sha256(rebuilt), Q.target_set_sha256(rebuilt)) == tuple(rec[k] for k in KEYS)
    kept = old if same else old.with_name(f"{old.name}-pre-rebuild")
    if same:
        shutil.rmtree(old)
    else:
        old.rename(kept)
    rebuilt.rename(old)
    FV.freeze([old], frozen, known, {
        "pre_rebuild": {"folder": str(kept)} | {k: rec[k] for k in KEYS},
        "r5_baseline": "hash chain" if same else "report line only: A3-val is the only Gate B1 baseline (B1.0 step 2)"})
    print(f"{rebuilt} -> {old}: " + ("hashes equal, folder replaced" if same else
                                     f"hashes differ, the pre-rebuild view stays at {kept}; A3-val is the only baseline"))


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("stamp")
    s.add_argument("--pre", type=Path, required=True, help="final_eval --hash-view record of the view")
    s.add_argument("--stamp", type=Path, required=True, help="the R5 launch stamp")
    w = sub.add_parser("swap")
    w.add_argument("--pre", type=Path, required=True)
    w.add_argument("--rebuilt", type=Path, required=True)
    w.add_argument("--frozen", type=Path, default=ROOT / "docs/reference/frozen_views.json")
    w.add_argument("--known", type=Path, default=ROOT / "runs/b1/series_z.json", help="as freeze_views --known")
    a = ap.parse_args(argv)
    if a.cmd == "stamp":
        stamp(a.pre, a.stamp)
    else:
        swap(a.pre, a.rebuilt, a.frozen, a.known)


if __name__ == "__main__":
    main()

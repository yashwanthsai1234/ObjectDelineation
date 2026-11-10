"""Build the ordered RL rows root of B1 (Implementation.md B1.4): one record per scheduled train window, in schedule
order, and a manifest written only when every hard B data gate passes.

    .venv/bin/python -I -m experiments.phase2.build_rl_rows --store data/store --target-law train-shifted
        --out data/views/lidc-train-b1 --steps 260 --per-step 32 [--draws N] --shares 0.15 0.75 0.10 --seed 0
        --procs 6 --exclude docs/reference/train_exclusions.json

1. `pool(train items, draws, seed, exclude)` with the `--exclude` file's keys (required; an empty list is `{}`). The
   target is `cell_shares` of the pool's draw-0 shifted windows (the train shifted law; train store only).
2. Draws table: for each value of `--draws-options` (default 1, 2, 4, 8) the schedule is simulated (the windows of a
   pool are the windows of the larger pool with draw < value, so one pool pass serves the table). A value passes when
   the schedule does not raise and DG7-DG9 hold on every block prefix; the rule takes the smallest passing value
   that reuses no disagreement window, else the smallest passing value. `--draws N` uses N instead (it must pass).
3. Each scheduled window is materialized once by `objdel.train.dataset.window_record` (item id = window_id: its own
   PNG and cache folder) in a spawn `multiprocessing.Pool(procs)`.
4. `lidc/train.jsonl` in schedule order (a reused disagreement window repeats its record); the gates of
   `experiments.phase2.audit_data.audit` run on the root with this plan; the manifest (build_eval_views' fields plus
   `schedule`, with `files_count` and `files_sha256` of the root without manifest.json) is written only when
   `quality.failures` is empty. Exit 1 and no manifest otherwise.
No upload (B1.4 step 5 is a separate, owner-run step). Reads untrusted data: run with python -I.
"""
import argparse
import hashlib
import json
import platform
import sys
import time
from collections import Counter
from dataclasses import asdict
from importlib.metadata import version
from itertools import batched
from multiprocessing import get_context
from pathlib import Path

import numpy as np

from objdel.data import quality as Q
from objdel.data import store, windows
from objdel.data.coverage import cell_coverage
from objdel.io import write_json

ROOT = Path(__file__).resolve().parents[2]
REF = ROOT / "docs" / "reference"
DRAWS_OPTIONS = (1, 2, 4, 8)


def _sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def target_law(pool) -> tuple[np.ndarray, list[float]]:
    """The train shifted law: cell shares and offsets of the pool's draw-0 shifted windows."""
    law = [w for w in pool.windows if w.kind == "shifted" and w.draw == 0]
    return windows.cell_shares(law), [w.offset for w in law]


def simulate(pool_windows, target, offsets, *, steps, per_step, shares, seed, disagreement_weight,
             block_steps) -> dict:
    """One row of the draws table: the schedule (or the error it raised), DG7-DG9 on every block prefix and the
    largest reuse of one disagreement window."""
    try:
        s = windows.schedule(pool_windows, target, steps, per_step=per_step, shares=shares, seed=seed,
                             disagreement_weight=disagreement_weight, block_steps=block_steps)
    except ValueError as e:
        return {"ok": False, "error": str(e), "max_reuse": None, "scheduled": None}
    nb = per_step * block_steps
    pre = [cell_coverage(s[:k], target, offsets) for k in range(nb, len(s) + nb, nb)]
    worst, tv, ks = (min(c["worst_ratio"] for c in pre), max(c["tv"] for c in pre), max(c["ks"] for c in pre))
    reuse = max(Counter(w.window_id for w in s if w.disagreement).values(), default=0)
    return {"ok": worst >= 0.8 and tv <= 0.10 and ks <= 0.06, "error": None, "worst": worst, "tv": tv, "ks": ks,
            "max_reuse": reuse, "scheduled": s}


def choose(table: dict[int, dict]) -> int | None:
    """Smallest passing value that reuses no disagreement window, else the smallest passing value."""
    ok = [d for d in sorted(table) if table[d]["ok"]]
    return next((d for d in ok if table[d]["max_reuse"] <= 1), ok[0] if ok else None)


def _record(job):
    from objdel.train import dataset as D  # lazy: needs the train extra
    out, item, window = job
    return D.window_record(out, item, window)


def materialize(store_root, out: Path, scheduled, procs: int) -> dict[str, dict]:
    """window_id -> record of every distinct scheduled window (images and caches written under out)."""
    by_item: dict[str, list] = {}
    for w in {w.window_id: w for w in scheduled}.values():
        by_item.setdefault(w.item_id, []).append(w)
    jobs = ((str(out), it, w) for it in store.load("lidc", "train", store_root) if it.item_id in by_item
            for w in by_item[it.item_id])
    recs = {}
    with get_context("spawn").Pool(procs) as pool:
        for chunk in batched(jobs, 24):  # bounds parent memory
            recs |= {r["item_id"]: r for r in pool.map(_record, chunk)}
    return recs


def files_digest(root: Path) -> tuple[int, str]:
    """(count, SHA-256 of the newline-joined sorted relative paths) of the root's files without manifest.json."""
    rel = sorted(p.relative_to(root).as_posix() for p in root.rglob("*") if p.is_file())
    rel = [r for r in rel if r != "manifest.json"]
    return len(rel), _sha("\n".join(rel).encode())


def build(a) -> int:
    from experiments.phase2 import audit_data
    from objdel.targets import target_set_hash
    from objdel.train import dataset as D
    from objdel.verify.reward.spec import fingerprint

    if a.steps <= 0 or a.steps % a.block_steps:
        raise SystemExit(f"--steps {a.steps} is not a positive multiple of block_steps {a.block_steps}")
    out = Path(a.out)
    if out.exists() and any(out.iterdir()):
        raise SystemExit(f"{out} exists and is not empty: remove it first (a root holds one build)")
    t0 = time.time()
    exclude_bytes = Path(a.exclude).read_bytes()
    excluded = frozenset(json.loads(exclude_bytes))
    shares = dict(zip(windows.KINDS, a.shares, strict=True))
    options = sorted(set(a.draws_options) | ({a.draws} if a.draws else set()))
    big = windows.pool(store.load("lidc", "train", a.store), max(options), a.seed, excluded)
    target, offsets = target_law(big)
    kw = {"steps": a.steps, "per_step": a.per_step, "shares": shares, "seed": a.seed,
          "disagreement_weight": a.disagreement_weight, "block_steps": a.block_steps}
    table = {d: simulate([w for w in big.windows if w.draw < d], target, offsets, **kw) for d in options}
    print("draws  ok     max_disagreement_reuse  worst   tv      ks      error")
    for d, r in table.items():
        print(f"{d:5d}  {str(r['ok']):5s}  {r['max_reuse']!s:22s}  " + (
            f"{r['worst']:.4f}  {r['tv']:.4f}  {r['ks']:.4f}" if r["error"] is None else f"-       -       -       "
                                                                                       f"{r['error']}"))
    rule = choose(table)
    draws = a.draws or rule
    print(f"draws rule choice {rule}, used {draws}")
    if draws is None or not table[draws]["ok"]:
        print(f"draws {draws}: the schedule raises or fails DG7-DG9 on a block prefix", file=sys.stderr)
        return 1
    pool = windows.pool(store.load("lidc", "train", a.store), draws, a.seed, excluded)  # exact dropped counts
    scheduled = table[draws]["scheduled"]
    assert {w.window_id for w in scheduled} <= {w.window_id for w in pool.windows}
    ids = [w.window_id for w in scheduled]

    recs = materialize(a.store, out, scheduled, a.procs)
    rows = [recs[i] for i in ids]
    text = "".join(json.dumps(r, sort_keys=True) + "\n" for r in rows)
    rel = "lidc/train.jsonl"
    (out / "lidc").mkdir(parents=True, exist_ok=True)
    (out / rel).write_text(text)
    t_mat = time.time() - t0

    plan = {"steps": a.steps, "planned_steps": a.steps, "block_steps": a.block_steps, "per_step": a.per_step,
            "seed": a.seed, "draws": draws, "shares": shares, "target_law": a.target_law,
            "target_sha256": _sha(np.ascontiguousarray(target, "<f8").tobytes()),
            "ids_sha256": _sha("\n".join(ids).encode()), "dropped": dict(pool.dropped),
            "disagreement_weight": a.disagreement_weight, "disagreement_max_reuse": table[draws]["max_reuse"],
            "exclusions_sha256": _sha(exclude_bytes), "excluded_ids": sorted(excluded),
            "draws_table": {str(d): {k: v for k, v in r.items() if k != "scheduled"} for d, r in table.items()}}
    report = audit_data.audit(out, a.store, a.eval_views, a.val_view, a.exclude, a.review, a.frozen, a.series_z,
                              a.model, a.audit_out, plan=plan)
    gs = Q.gates(report)
    for g in gs:
        print(f"{g.name:6s} {'report' if not g.hard else 'ok' if g.ok else 'FAIL':6s} {g.value:10.4f}  {g.limit}")
    if a.json:
        Path(a.json).write_text(json.dumps({"gates": [asdict(g) for g in gs], "report": report}, indent=1) + "\n")
    if bad := Q.failures(gs):
        print(f"GATES FAILED: {', '.join(g.name for g in bad)}; no manifest written ({out})", file=sys.stderr)
        return 1
    count, digest = files_digest(out)
    uniq = list(recs.values())
    triples = [(r["item_id"], k, t) for r in uniq for k, t in zip(r["raters"], r["targets"], strict=True) if t]
    m = {"view": out.name, "seed": a.seed,
         "store_manifest_sha256": _sha(Path(a.store, "manifest.json").read_bytes()),
         "reward_fingerprint": fingerprint(), "prompts": D.PROMPTS,
         "versions": {"python": platform.python_version(), **{p: version(p) for p in D._PINS}},
         "files": {rel: {"sha256": _sha(text.encode())}},
         "sets": {"lidc/train": {
             "items": len(rows), "windows": len(uniq), "skipped_no_image": 0, "targets": len(triples),
             "target_set_sha256": target_set_hash(triples),
             "target_errors": [e for r in uniq for e in r["errors"]],
             "rl_items": sum(r["rl"] for r in rows),
             "rl_targets": sum(t is not None for r in uniq if r["rl"] for t in r["targets"]),
             "disagreement_items": sum(r["disagreement"] for r in rows),
             "kept_side": 0, "fallback": 0, "padded": report["windows"]["padded"],
             "kinds": dict(Counter(r["view"]["kind"] for r in rows)),
             "side": {str(s): n for s, n in sorted(Counter(r["view"]["side"] for r in rows).items())}}},
         "schedule": plan | {"gates": [asdict(g) for g in gs], "files_count": count, "files_sha256": digest}}
    write_json(out / "manifest.json", m)
    print(json.dumps({"root": str(out), "rows": len(rows), "windows": len(uniq), "files_count": count,
                      "files_sha256": digest, "manifest_sha256": _sha((out / "manifest.json").read_bytes()), "draws": draws,
                      "materialize_s": round(t_mat), "total_s": round(time.time() - t0)}, indent=1))
    return 0


def parse_args(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--store", type=Path, default=ROOT / "data/store")
    ap.add_argument("--target-law", choices=["train-shifted"], default="train-shifted")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--steps", type=int, required=True)
    ap.add_argument("--per-step", type=int, default=32)
    ap.add_argument("--block-steps", type=int, default=20)
    ap.add_argument("--draws", type=int, default=None, help="default: the draws-table rule (module docstring)")
    ap.add_argument("--draws-options", type=int, nargs="+", default=list(DRAWS_OPTIONS))
    ap.add_argument("--shares", type=float, nargs=3, default=[windows.SHARES[k] for k in windows.KINDS],
                    metavar=("CENTRED", "SHIFTED", "RESCALED"))
    ap.add_argument("--disagreement-weight", type=float, default=4.0)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--procs", type=int, default=6)
    ap.add_argument("--exclude", type=Path, required=True, help="train exclusions JSON (no default; none = {})")
    ap.add_argument("--eval-views", type=Path, default=ROOT / "data/eval_views", help="val and test views")
    ap.add_argument("--val-view", type=Path, default=None)
    ap.add_argument("--review", type=Path, default=REF / "neardup_review.json")
    ap.add_argument("--frozen", type=Path, default=REF / "frozen_views.json")
    ap.add_argument("--series-z", type=Path, default=ROOT / "runs/b1/series_z.json")
    ap.add_argument("--model", default="qwen3-vl-2b", help="tokenizer and token limit of DG12")
    ap.add_argument("--audit-out", type=Path, default=ROOT / "runs/b1/audit", help="near-duplicate PNGs (DG3)")
    ap.add_argument("--json", type=Path, default=None, help="gate report")
    return ap.parse_args(argv)


def main(argv=None) -> int:
    return build(parse_args(argv))


if __name__ == "__main__":
    raise SystemExit(main())

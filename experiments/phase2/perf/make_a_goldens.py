"""Mac golden of A.9.6 G02b: the fair ceiling's `doctor` and `model` (objdel.eval.ceiling.item) on 50 LIDC val items,
on the CPU build of the checker. The box test (tests/test_a_gpu.py::test_g02b_builds) scores the same items and texts
on the CUDA build and compares with ==, since one venv holds one build.

    .venv/bin/python -I -m experiments.phase2.perf.make_a_goldens [--prepared data/prepared] \
        [--out tests/goldens/a_g02b.json]

Items: the first 50 RL records of <prepared>/lidc/val.jsonl, sorted by item_id, that hold a target text. Texts (K = 8):
the item's target texts cycled to 6, plus 2 adversarial texts drawn with seed 0 from ADVERSARIAL; truncated all
False. The golden stores the ids, the texts, their SHA-256 and the two values as float hex. Reads untrusted data: run
with python -I.
"""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from objdel.eval import ceiling

ADVERSARIAL = tuple(f'{{"nodule":[{{"outer":{o},"holes":{h}}}]}}' for o, h in (
    ("[[100,100],[500,500],[900,100],[900,900],[500,500],[100,900]]", "[]"),  # self-touching (vertex pinch)
    ("[[100,100],[900,100],[900,900],[100,900]]", "[[[300,300],[300,700],[700,700],[700,300]]]"),  # a hole
    ("[[0,0],[1000,0],[1000,1000],[0,1000]]", "[]"),  # the whole frame
    ("[[500,500],[520,500],[520,520],[500,520]]", "[]"),  # one pixel of a 50 px frame
    ("[[0,0],[1000,0],[0,1000]]", "[]"),  # grid corners
))


def texts_sha256(texts: dict[str, list[str]]) -> str:
    return hashlib.sha256(json.dumps(texts, sort_keys=True).encode()).hexdigest()


def inputs(prepared: Path, n: int = 50) -> dict[str, list[str]]:
    """{item_id: 8 texts} of the golden's items (module docstring)."""
    from objdel.train import dataset

    recs = sorted((r for r in dataset.records(prepared, "lidc", "val") if r["rl"] and any(r["targets"])),
                  key=lambda r: r["item_id"])[:n]
    rng = np.random.default_rng(0)
    return {r["item_id"]: ([t for t in r["targets"] if t] * 6)[:6] + [ADVERSARIAL[i] for i in rng.choice(5, 2)]
            for r in recs}


def values(prepared: Path, texts: dict[str, list[str]]) -> dict[str, list[str]]:
    """{item_id: [doctor, model] as float hex} on this venv's build."""
    from objdel.train import dataset

    caches, recs = dataset.cache_loader(prepared), {r["item_id"]: r for r in dataset.records(prepared, "lidc", "val")}
    out = {}
    for i, t in texts.items():
        row = ceiling.item(caches(i)["nodule"], t, [False] * len(t), recs[i]["W"], recs[i]["H"])
        out[i] = [row["doctor"].hex(), row["model"].hex()]
    return out


def main(argv=None) -> None:
    from objdel.geom import engine as _gpucheck

    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--prepared", type=Path, default=Path("data/prepared"))
    ap.add_argument("--out", type=Path, default=Path("tests/goldens/a_g02b.json"))
    a = ap.parse_args(argv)
    if _gpucheck.build() != "cpu":
        raise SystemExit("make_a_goldens: the golden is the CPU build's (run on the Mac)")
    texts = inputs(a.prepared)
    a.out.write_text(json.dumps({"texts_sha256": texts_sha256(texts), "texts": texts,
                                 "values": values(a.prepared, texts)}, indent=1) + "\n")


if __name__ == "__main__":
    main()

"""The synthetic items of D-R05 (Implementation.md D.9): 50 items of K = 8 completions against their raters, 35
nodule items (k = 2..4, frames 33..64 px, raters with and without a hole) and 15 disc + cup items (k = 3..6). The
completions mix rater copies, shifted copies, exact and re-spelled duplicates, an empty object list, a
self-intersecting (illegal) ring, unparseable and truncated texts, and one item of 8 identical answers.
tests/goldens/d_r05_metrics.json holds metrics() of every item, written before D7 moved its pixel work to the engine.
"""
import numpy as np

from objdel.data.targets.build import build_target, min_areas

BOWTIE = {"nodule": '{"nodule":[{"outer":[[300,300],[700,700],[700,300],[300,700]],"holes":[]}]}',
          "riga": '{"disc":[{"outer":[[300,300],[700,700],[700,300],[300,700]],"holes":[]}],"cup":[]}'}


def _disk(H, W, cy, cx, r, hole=0.0):
    y, x = np.mgrid[:H, :W]
    d = np.hypot(y - cy, x - cx)
    return (d <= r) & (d >= hole)


def _text(masks, names, raters):
    return build_target(masks, names, min_area=min_areas(raters)).text


def _completions(rng, raters, names, kind):
    k = len(raters[names[0]])
    copy = [_text({n: raters[n][r] for n in names}, names, raters) for r in range(k)]
    shift = [_text({n: np.roll(raters[n][r], int(rng.integers(1, 4)), int(rng.integers(0, 2))) for n in names},
                   names, raters) for r in range(min(k, 2))]
    pool = [*copy, *shift, copy[0], copy[0].replace(",", ", "), BOWTIE[kind], "not json", "",
            '{"nodule":[]}' if kind == "nodule" else '{"disc":[],"cup":[]}']
    texts = [pool[int(i)] for i in rng.permutation(len(pool))[:8]]
    truncated = [bool(rng.random() < 0.1) for _ in texts]
    return texts, truncated


def items() -> list[tuple]:
    """(texts, truncated, raters {name: bool (k, H, W)}, W, H, names) of the 50 items."""
    rng = np.random.default_rng(31)
    out = []
    for i in range(35):
        H, W = [(40, 48), (64, 64), (33, 57)][i % 3]
        k = 2 + i % 3
        cy, cx = H / 2 + rng.uniform(-3, 3), W / 2 + rng.uniform(-3, 3)
        st = np.stack([_disk(H, W, cy + rng.uniform(-2, 2), cx + rng.uniform(-2, 2), rng.uniform(5, 12),
                             hole=2.5 if rng.random() < 0.3 else 0.0) for _ in range(k)])
        raters = {"nodule": st}
        texts, tr = _completions(rng, raters, ("nodule",), "nodule")
        out.append((texts, tr, raters, W, H, ("nodule",)))
    for i in range(15):
        H, W = [(100, 100), (80, 96)][i % 2]
        k = 3 + i % 4
        cy, cx = H / 2 + rng.uniform(-4, 4), W / 2 + rng.uniform(-4, 4)
        disc = np.stack([_disk(H, W, cy, cx, rng.uniform(26, 32)) for _ in range(k)])
        cup = np.stack([_disk(H, W, cy + rng.uniform(-2, 2), cx + rng.uniform(-2, 2), rng.uniform(9, 15))
                        for _ in range(k)])
        raters = {"disc": disc, "cup": cup}
        texts, tr = _completions(rng, raters, ("disc", "cup"), "riga")
        out.append((texts, tr, raters, W, H, ("disc", "cup")))
    _, _, raters, W, H, names = out[0]
    same = _text({n: raters[n][0] for n in names}, names, raters)
    out[1] = ([same] * 8, [False] * 8, raters, W, H, names)  # all-identical legal answers
    return out


def val_items(root, ds: str, n: int = 20) -> list[tuple]:
    """(item_id, texts, truncated, W, H, names, caches) of the first n RL records of `ds`/val.jsonl under a prepared
    root: K = 8 answers cycled from the rater targets, with a self-intersecting and an unparseable one (D-R05's real
    items). Caches loaded against their record digests."""
    import json
    from pathlib import Path

    from objdel.core import datasets
    from objdel.verify.reward.cache import ExpertCache

    root, names = Path(root), datasets.get(ds).names
    out = []
    for line in (root / ds / "val.jsonl").read_text().splitlines():
        r = json.loads(line)
        if not r["rl"]:
            continue
        good = [t for t in r["targets"] if t is not None]
        texts = [good[i % len(good)] for i in range(6)] + [BOWTIE["nodule" if ds == "lidc" else "riga"], "not json"]
        caches = {m: ExpertCache.load(root / "caches" / r["item_id"] / f"{m}.npz", d) for m, d in r["caches"].items()}
        out.append((r["item_id"], texts, [False] * 7 + [True], r["W"], r["H"], names, caches))
        if len(out) == n:
            break
    return out

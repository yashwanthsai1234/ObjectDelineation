"""For near-duplicate pairs of the same photo: are the 6 raters' masks copies or independent re-annotations?"""
import json, numpy as np
from pathlib import Path
from PIL import Image
O = Path(__file__).parent
def iou(a, b):
    u = (a | b).sum(); return float((a & b).sum() / u) if u else 1.0
rows = []
for p in json.load(open(O / "dedup/pairs.json")):
    if p["kind"] != "same_photo":
        continue
    for e in range(1, 7):
        fa, fb = [O / "out/masks" / f"{x.replace('/image', '_image')}-{e}.png" for x in (p["a"], p["b"])]
        if not (fa.exists() and fb.exists()):
            continue
        A, B = np.asarray(Image.open(fa)), np.asarray(Image.open(fb))
        if A.shape != B.shape:
            rows.append({"a": p["a"], "b": p["b"], "e": e, "shape_mismatch": True}); continue
        rows.append({"a": p["a"], "b": p["b"], "e": e, "identical": bool(np.array_equal(A, B)),
                     "disc_iou": iou(A >= 100, B >= 100), "cup_iou": iou(A >= 200, B >= 200)})
ok = [r for r in rows if "identical" in r]
d = np.array([r["disc_iou"] for r in ok]); c = np.array([r["cup_iou"] for r in ok])
s = {"pairs": len({(r["a"], r["b"]) for r in rows}), "rater_pairs": len(ok), "identical": sum(r["identical"] for r in ok),
     "shape_mismatch": sum("shape_mismatch" in r for r in rows),
     "disc_iou": {"min": float(d.min()), "median": float(np.median(d)), "max": float(d.max())},
     "cup_iou": {"min": float(c.min()), "median": float(np.median(c)), "max": float(c.max())},
     "pairs_all_identical": sorted({(r["a"], r["b"]) for r in ok if all(x["identical"] for x in ok if (x["a"], x["b"]) == (r["a"], r["b"]))})}
(O / "dedup/dup_masks.json").write_text(json.dumps({"summary": s, "rows": rows}, indent=1))
print(json.dumps(s, indent=1))

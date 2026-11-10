"""Overlay extracted disc (green) / cup (magenta) on the rater image. python -I overlay.py out.png sub:idx:e ..."""
import sys
import numpy as np
from pathlib import Path
from PIL import Image
from scipy import ndimage
O = Path(__file__).parent
R = Path(__file__).resolve().parents[4] / "data/raw/riga"
tiles = []
for spec in sys.argv[2:]:
    sub, idx, e = spec.split(":")
    a = np.asarray(Image.open(next((R / sub).glob(f"image{idx}-{e}.*"))).convert("RGB")).copy()
    v = np.asarray(Image.open(O / "out/masks" / f"{sub.replace('/', '_')}_image{idx}-{e}.png"))
    d, c = v >= 100, v >= 200
    ys, xs = np.nonzero(d)
    cy, cx = int(ys.mean()), int(xs.mean()); r = int(max(np.ptp(ys), np.ptp(xs)) * 0.75)
    ov = a.copy()
    ov[d & ~ndimage.binary_erosion(d, iterations=2)] = (0, 255, 0)
    ov[c & ~ndimage.binary_erosion(c, iterations=2)] = (255, 0, 255)
    sl = (slice(max(0, cy - r), cy + r), slice(max(0, cx - r), cx + r))
    tiles.append(np.asarray(Image.fromarray(np.concatenate([a[sl], ov[sl]], 1)).resize((600, 300))))
Image.fromarray(np.concatenate(tiles, 0)).save(sys.argv[1])

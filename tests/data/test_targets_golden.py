"""D-R02 (gate D-G4): build_target's text on 300 synthetic fixtures equals the golden written by the code before D3
(commit 7eaa9d8: rings from the engine, components from cv2.connectedComponents). Fixtures: LIDC-like nodules (discs
with holes, specks, several components; min_area 5, the mask's own rule, or above every component), RIGA-like disc
and cup pairs made in code (no real RIGA item; the rule of the rater's own masks or of the union), and sparse noise.
D4 made min_area required: the 100 fixtures that had none now pass the rules above, and the code before D4 (bb8a3f7)
gives the same 300 rows. A ValueError is part of the behaviour, so the golden keeps its message. Regenerate (only on
purpose): python -I -m tests.data.test_targets_golden --write"""
import hashlib
import json
import sys
from pathlib import Path

import numpy as np

from objdel.data.targets.build import build_target, min_areas

GOLDEN = Path(__file__).resolve().parents[1] / "goldens" / "d_r02_targets.json"


def _disc(h, w, cy, cx, r):
    y, x = np.ogrid[:h, :w]
    return (y - cy) ** 2 + (x - cx) ** 2 <= r * r


def fixtures(seed=20261109):
    """300 (masks, names, min_area) build_target inputs."""
    rng = np.random.default_rng(seed)
    out = []
    for i in range(300):
        kind = i % 6
        h, w = (int(v) for v in rng.integers(12, 90, 2))
        if kind < 3:  # nodule: discs minus discs; kind 2 adds specks
            m = np.zeros((h, w), bool)
            for _ in range(int(rng.integers(1, 5))):
                m |= _disc(h, w, *rng.integers(0, (h, w)), int(rng.integers(1, 12)))
            for _ in range(int(rng.integers(0, 3))):
                m &= ~_disc(h, w, *rng.integers(0, (h, w)), int(rng.integers(1, 5)))
            if kind == 2 and i % 3 != 2:
                m |= rng.random((h, w)) < 0.004
            ma = [{"nodule": 5}, min_areas({"nodule": m[None]}), {"nodule": int(m.sum()) + 1}][i % 3]
            out.append(({"nodule": m}, ("nodule",), ma))
        elif kind < 5:  # disc and cup of k raters; kind 4 puts the cup near the disc edge
            k = int(rng.integers(1, 4))
            cy, cx, r = h // 2, w // 2, max(4, min(h, w) // 2 - 2)
            discs = np.stack([_disc(h, w, cy + int(rng.integers(-1, 2)), cx, r - j) for j in range(k)])
            cups = np.stack([_disc(h, w, cy, cx + int(rng.integers(-1, 2)), max(1, r // 2 - j)) for j in range(k)])
            if kind == 4:
                cups = np.stack([_disc(h, w, cy, cx + r // 3, max(1, r // 2)) for _ in range(k)])
            ma = min_areas({"disc": discs[:1], "cup": cups[:1]} if i % 2 else {"disc": discs, "cup": cups})
            out.append(({"disc": discs[0], "cup": cups[0]}, ("disc", "cup"), ma))
        else:  # sparse noise
            m = rng.random((h, w)) < rng.uniform(0.01, 0.08)
            out.append(({"nodule": m}, ("nodule",), [{"nodule": 0}, min_areas({"nodule": m[None]})][i % 2]))
    return out


def results():
    rows = []
    for masks, names, ma in fixtures():
        try:
            rows.append(hashlib.sha256(build_target(masks, names, min_area=ma).text.encode()).hexdigest())
        except ValueError as e:
            rows.append(f"ValueError: {e}")
    return rows


def test_target_texts_equal_the_pre_d3_golden(monkeypatch):
    """D3 groups holes by the engine's ring part, so the ring step needs no cv2 labelling."""
    import cv2

    def no_cv2(*a, **k):
        raise AssertionError("cv2.connectedComponents called")
    monkeypatch.setattr(cv2, "connectedComponents", no_cv2)
    want = json.loads(GOLDEN.read_text())["rows"]
    assert len(want) == 300 and sum(not r.startswith("ValueError") for r in want) >= 250
    got = results()
    assert [i for i, (a, b) in enumerate(zip(got, want, strict=True)) if a != b] == []


if __name__ == "__main__" and sys.argv[1:] == ["--write"]:
    GOLDEN.write_text(json.dumps({"source": "objdel.data.targets.build.build_target before D3 (commit 7eaa9d8)",
                                  "rows": results()}, indent=0) + "\n")

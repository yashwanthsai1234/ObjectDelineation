"""No CPU geometry at runtime (Implementation.md D-G11, S01 runtime half): with scipy's labelling and EDT and OpenCV's
three labelling calls raising, the target builder (LIDC and RIGA), the cache build with its device tables, the cache
load, `certify` and `metrics()` run and give what they give unpatched. `scipy.optimize.linear_sum_assignment` stays
allowed. The static half (no `tests.oracles` import in src, no numpy polygonizer, the seam's names) is
tests/experiments/test_rings.py::test_the_numpy_oracles_left_the_package."""
import cv2
import pytest
from scipy import ndimage

from objdel.geom.certify import certify
from objdel.verify.metrics import metrics
from objdel.verify.reward.cache import ExpertCache
from tests.verify.test_eval import nodule_raters, riga_raters, targets

CASES = (("lidc", nodule_raters, 50), ("riga", riga_raters, 100))
BANNED = ((ndimage, "label"), (ndimage, "distance_transform_edt"), (cv2, "connectedComponents"),
          (cv2, "connectedComponentsWithAlgorithm"), (cv2, "connectedComponentsWithStats"))


def run(tmp_path) -> dict:
    out = {}
    for name, make, side in CASES:
        raters = make()
        names = tuple(raters)
        texts = targets(raters, names)
        caches = {n: ExpertCache.build(raters[n]) for n in names}
        for c in caches.values():
            c.device()
        digest = caches[names[0]].save(tmp_path / f"{name}.npz")
        out[name] = {"texts": texts, "digest": digest, "loaded": ExpertCache.load(tmp_path / f"{name}.npz").digest,
                     "certify": repr(certify(raters[names[0]][0])),
                     "metrics": repr(metrics(texts, raters, side, side, names, caches=caches))}
    return out


def test_no_cpu_labelling_or_edt_runs_and_the_results_are_unchanged(tmp_path, monkeypatch):
    for d in ("unpatched", "patched"):
        (tmp_path / d).mkdir()
    want = run(tmp_path / "unpatched")

    def boom(*a, **k):
        raise AssertionError("CPU geometry at runtime")

    for mod, name in BANNED:
        monkeypatch.setattr(mod, name, boom)
    with pytest.raises(AssertionError, match="CPU geometry"):
        ndimage.label([1])  # the patch is live
    assert run(tmp_path / "patched") == want

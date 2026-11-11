"""objdel.core.datasets: the registry holds today's LIDC values verbatim (each replaced a table in another layer),
lookups fail loudly, and every table that read a dataset now reads the entry."""
import numpy as np
import pytest

from objdel.core import datasets
from objdel.core.datasets import FrameConfig, Structure

LIDC_PROMPT = ('Outline the lung nodule in this CT crop. Answer with JSON only: '
               '{"nodule":[{"outer":[[x,y],...],"holes":[[[x,y],...]]}]}, '
               'integers 0-1000 relative to the crop, outer rings clockwise, holes counter-clockwise.')


def test_lidc_entry_equals_the_tables_it_replaced():
    from objdel.data.targets import augment, frame

    d = datasets.get("lidc")
    assert d.names == ("nodule",)  # geom.text SCHEMAS
    assert d.structures == (Structure("nodule", 2.0, 2, 0.60),)  # metrics NSD_TAU, BIOU_D; gates MEDIAN_IOU_TC
    assert d.prompt == LIDC_PROMPT  # policy.prompt PROMPTS
    assert d.pixels == ("int16", 2)  # store IMAGE
    assert d.window == frame.LUNG_WINDOW == (-1000, 400)
    assert (d.window_jitter, d.noise_max) == (((-300, 60), (1400, 80)), 10.0)  # augment WINDOW_JITTER, NOISE_MAX
    assert d.frame == FrameConfig(max_side=500, min_side=50, margin_px=8, margin_pct=25, image_side=512)  # FRAMES
    assert d.ops == augment.OPS == ("id", "r90", "r180", "r270", "fh", "fv", "tr", "at")  # DATASET_OPS
    assert d.rotated_op == "r90"  # build_eval_views ROTATION


def test_the_library_registers_lidc_only():
    """RIGA is retired: the tests' conftest registers its facts as their second dataset, the library does not."""
    import importlib
    import subprocess
    import sys

    assert datasets.DATASETS["lidc"] is datasets.LIDC
    out = subprocess.run([sys.executable, "-c", "from objdel.core import datasets; print(sorted(datasets.DATASETS))"],
                         capture_output=True, text=True, check=True).stdout
    assert out.strip() == "['lidc']"
    assert importlib.import_module("objdel.policy.prompt").PROMPTS["lidc"] == LIDC_PROMPT


def test_every_entry_is_consistent():
    from objdel.data.targets import augment, frame

    for d in datasets.DATASETS.values():
        assert d.names and len(set(d.names)) == len(d.names), d.name
        assert set(d.ops) <= set(augment.OPS) and d.ops[0] == "id", d.name
        assert d.rotated_op is None or d.rotated_op in augment.OPS, d.name
        assert d.frame.max_side in frame.DIVISORS and np.dtype(d.pixels[0]) and d.pixels[1] in (2, 3), d.name
        assert (d.window is None) == (d.pixels[0] == "uint8"), d.name  # a non-uint8 image needs a window


def test_lookups_fail_loudly(monkeypatch):
    with pytest.raises(KeyError, match="unknown dataset 'polyp'.*registered"):
        datasets.get("polyp")
    with pytest.raises(KeyError, match="'cornea': 0 different"):
        datasets.structure("cornea")
    assert datasets.structure("nodule").nsd_tau == 2.0
    other = datasets.Dataset("ct2", (Structure("nodule", 3.0, 3),), "p", ("int16", 2), (-1000, 400))
    monkeypatch.setitem(datasets.DATASETS, "ct2", other)
    with pytest.raises(KeyError, match="'nodule': 2 different"):  # one name, two thresholds: ambiguous
        datasets.structure("nodule")

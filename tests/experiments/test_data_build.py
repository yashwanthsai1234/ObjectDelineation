"""experiments/phase1/build_data: the report counts what the store gets (review #4)."""
from collections import Counter

import numpy as np

from experiments.phase1 import build_data
from objdel.data import riga


def test_report_counts_the_riga_images_the_store_gets():
    """The F7 item (reference "maxrg6:prime_is_rater<e>_copy", Magrabia Female image48) keeps its image; only the
    BinRushed1 reconstruction ("maxrg6") has none. The report counted only "prime" (94 test images, store 95)."""
    disc = np.zeros((2, 8, 8), bool)
    disc[:, 2:6, 2:6] = True
    cup = disc & np.pad(np.ones((2, 2, 2), bool), ((0, 0), (3, 3), (3, 3)))
    masks = {e: (disc[0], cup[0], "A") for e in (1, 2)}
    ref = np.zeros((8, 8, 3), np.uint8)
    stats = {"riga_stages": Counter(), "riga_reference": Counter(), "riga_with_image": Counter()}
    got = [build_data.riga_item(stats, riga.Photo("Magrabia/MagrabiFemale", i, None, {}), (ref, kind, masks), "test",
                                f"g{i}", [])
           for i, kind in enumerate(["prime", "maxrg6", "maxrg6:prime_is_rater1_copy"])]
    assert [it.image is not None for it in got] == [True, False, True]
    assert stats["riga_with_image"] == {"test": 2}
    assert stats["riga_reference"] == {"prime": 1, "maxrg6": 1, "maxrg6:prime_is_rater1_copy": 1}
    assert stats["riga_stages"] == {"A": 6}

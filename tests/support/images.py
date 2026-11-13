"""Random images for the tests."""
import numpy as np


def random_image(h=64, w=64, seed=0, high=256):
    """A PIL RGB image of seeded uniform noise (rng = default_rng(seed), pixels in [0, high))."""
    from PIL import Image

    return Image.fromarray(np.random.default_rng(seed).integers(0, high, (h, w, 3), dtype=np.uint8))

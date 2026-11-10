"""Geometry layer: outline text, shapes, rings, topology and the CUDA engine.

Re-exports only `polygonize` (the polygonizer seam).
"""

from objdel.geom.rings import polygonize

__all__ = ["polygonize"]

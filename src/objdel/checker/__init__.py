"""Exact geometry of binary masks: bit-quad invariants, Betti numbers, polygons and certificates.

Depends only on numpy, scipy and pillow so it can run inside any training loop.
"""
from .bitquad import bitquads
from .certify import certify
from .polygonize import polygonize
from .topology import betti, betti3d, size_filter

__all__ = ["betti", "betti3d", "bitquads", "certify", "polygonize", "size_filter"]

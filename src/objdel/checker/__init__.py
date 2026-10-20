"""Exact geometry of binary masks (bit-quad invariants, per-instance topology, polygons, certificates; numpy and
OpenCV only) and of integer polygons (`polycore`, the C++ extension: rasterize, validate, signed area)."""
from .bitquad import bitquads
from .certify import certify
from .polygonize import polygonize
from .topology import InstanceTopology, as_mask, betti, instance_topology, min_area, signature, size_filter

__all__ = [
    "InstanceTopology", "as_mask", "betti", "bitquads", "certify", "instance_topology", "min_area", "polygonize",
    "signature", "size_filter",
]

"""Exact geometry of binary masks (bit-quad invariants, per-instance topology, polygons, certificates) and of voxel
volumes (`volume`: 3-D Betti numbers of stacked slices). Integer polygons (rasterize, validate) and the reward's
scoring are the CUDA checker's, objdel._gpucheck (cuda/; objdel.geometry.to_masks, objdel.reward.score)."""
from .bitquad import bitquads
from .certify import certify
from .rings import polygonize
from .topology import InstanceTopology, as_mask, betti, instance_topology, min_area, signature, size_filter
from .volume import betti3d

__all__ = [
    "InstanceTopology", "as_mask", "betti", "betti3d", "bitquads", "certify", "instance_topology", "min_area",
    "polygonize", "signature", "size_filter",
]

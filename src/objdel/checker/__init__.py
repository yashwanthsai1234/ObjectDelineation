"""Exact geometry of binary masks (rings, bit-quads and Betti numbers from the CUDA engine, per-instance topology,
certificates) and of voxel volumes (`volume`: 3-D Betti numbers of stacked slices). Integer polygons (rasterize,
validate) and the reward's scoring are the CUDA checker's, objdel._gpucheck (cuda/; objdel.geometry.to_masks,
objdel.reward.score)."""
from .certify import certify, certify_many
from .rings import MaskRings, Ring, polygonize, trace
from .topology import (
    InstanceTopology,
    as_mask,
    betti,
    instance_topology,
    min_area,
    signature,
    signatures,
    size_filter,
    topologies,
)
from .volume import betti3d

__all__ = [
    "InstanceTopology", "MaskRings", "Ring", "as_mask", "betti", "betti3d", "certify", "certify_many",
    "instance_topology", "min_area", "polygonize", "signature", "signatures", "size_filter", "topologies", "trace",
]

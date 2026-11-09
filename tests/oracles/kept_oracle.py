"""The interval loop that objdel.targets.build used when min_area was unknown (`_kept`, removed by D4: every caller now
passes the reward's rule). Reference for the ladder's topology test (tests only), on the scipy topology oracle.

kept(a, b, lo, hi): signature(a, t) == signature(b, t) for every rule t in [lo, hi]. Between consecutive rules where
either filtered mask can change (a raw component of t' - 1 px is removed at t'; a hole of the filtered mask of t' - 1
px is filled at t'), the opening and the set of filled holes are constant, so one check per interval suffices."""
import numpy as np

from .topology_oracle import hole_areas, instance_topology, signature, size_filter


def kept(a, b, lo: int, hi: int) -> bool:
    if lo == hi:  # the reward's rule is known
        return signature(a, lo) == signature(b, lo)
    raw = [area for m in (a, b) for area, _ in instance_topology(m)]
    t = lo
    while t <= hi:
        fa, fb = size_filter(a, t), size_filter(b, t)
        if sorted(h for _, h in instance_topology(fa)) != sorted(h for _, h in instance_topology(fb)):
            return False
        sizes = np.array(raw + hole_areas(fa) + hole_areas(fb), np.int64)
        sizes = sizes[sizes >= t]
        if not sizes.size:
            return True
        t = int(sizes.min()) + 1
    return True

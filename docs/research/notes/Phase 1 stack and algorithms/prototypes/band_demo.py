import numpy as np
from scipy import ndimage as ndi
def disk(n, r): y, x = np.mgrid[:n, :n]; return np.hypot(y-n/2+.5, x-n/2+.5) <= r
def boundary(m): return m & ~ndi.binary_erosion(m, np.ones((3,3)))  # inner 8-boundary pixels
def scores(pred, anns, tau=2.0):
    s = np.sum(anns, 0); U = s > 0; I = s == len(anns)
    miss = (I & ~pred).sum(); spill = (pred & ~U).sum()
    F0 = max(0, 1 - (miss+spill)/U.sum())                               # current
    perim = boundary(U).sum()
    F1 = max(0, 1 - (miss+spill)/(tau*perim))                          # perimeter-normalised (mean px excursion / tau)
    dU = ndi.distance_transform_edt(~U); dI = ndi.distance_transform_edt(I)  # dist outside U; depth inside I
    w = np.minimum(1, dU/tau)*(pred & ~U) + np.minimum(1, dI/tau)*(I & ~pred)
    F2 = max(0, 1 - w.sum()/perim)                                       # distance-weighted (errors within tau partly forgiven)
    band = U & ~I; dband = ndi.distance_transform_edt(~band)
    bp = boundary(pred); F3 = (dband[bp] <= tau).mean() if bp.any() else 0.  # NSD-like: pred boundary within tau of band
    kA = 0.1*np.sqrt(U.sum()); Ud = ndi.distance_transform_edt(~U) <= kA; Id = ndi.distance_transform_edt(I) > kA
    F4 = max(0, 1 - ((pred & ~Ud).sum() + (Id & ~pred).sum())/U.sum())   # tolerance = k*sqrt(area), area-normalised
    iou_maj = (pred & (s*2 > len(anns))).sum()/ (pred | (s*2 > len(anns))).sum()
    return F0, F1, F2, F3, F4, iou_maj
print(f"{'r':>4} {'|U|':>6} {'case':>8} |  F0cur  F1perim F2dist F3nsd  F4sqrtA  IoUmaj")
for r in (4, 6, 10, 20, 40, 80):
    n = int(2*r + 30); anns = [disk(n, r-1), disk(n, r), disk(n, r+1)]
    for name, pr in [("spill1", r+2), ("spill3", r+4), ("miss1", r-2), ("miss3", r-4)]:
        if pr <= 0: continue
        f = scores(disk(n, pr), anns)
        print(f"{r:>4} {np.sum(anns,0).astype(bool).sum():>6} {name:>8} | " + " ".join(f"{x:6.2f}" for x in f))

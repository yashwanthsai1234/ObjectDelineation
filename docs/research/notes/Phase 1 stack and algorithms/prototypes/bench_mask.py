"""Micro-benchmark of mask-side primitives on crop sizes 64..1024 (single thread)."""
import time, os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
from scipy import ndimage as ndi
from scipy.optimize import linear_sum_assignment
import cv2, cc3d, edt
from skimage.measure import label as sk_label, euler_number as sk_euler

cv2.setNumThreads(1)
rng = np.random.default_rng(0)
FG8 = np.ones((3, 3), bool)
BG4 = ndi.generate_binary_structure(2, 1)


def make_mask(n, kind):
    if kind == "blob":  # one smooth object w/ occasional holes (polygon-like)
        f = ndi.gaussian_filter(rng.standard_normal((n, n)), n / 12)
        yy, xx = np.mgrid[:n, :n]
        r = np.hypot(yy - n / 2, xx - n / 2) / (n / 2)
        return (f / f.std() * 0.6 - 2.2 * r + 1.0) > 0
    # noisy: many components/holes (stress)
    f = ndi.gaussian_filter(rng.standard_normal((n, n)), 2.0)
    return f > 0


def timeit(fn, reps):
    fn()
    t = time.perf_counter()
    for _ in range(reps):
        fn()
    return (time.perf_counter() - t) / reps * 1e6  # us


# ---------- bit-quad Euler (8-conn FG) in numpy
def chi8_bitquad(m):
    p = np.pad(m, 1).astype(np.uint8)
    code = p[:-1, :-1] | (p[:-1, 1:] << 1) | (p[1:, :-1] << 2) | (p[1:, 1:] << 3)
    h = np.bincount(code.ravel(), minlength=16)
    q1 = h[[1, 2, 4, 8]].sum(); q3 = h[[7, 11, 13, 14]].sum(); qd = h[[6, 9]].sum()
    return (q1 - q3 - 2 * qd) // 4


# ---------- per-component holes via first-raster-pixel trick
from holes_fix import holes_per_component, holes_per_component_cv2


def holes_bruteforce(m):
    lab, n = ndi.label(m, FG8)
    out = []
    for k in range(1, n + 1):
        c = np.pad(lab == k, 1)
        _, nb = ndi.label(~c, BG4)
        out.append(nb - 1)
    return np.array(out)


def same_partition(a, b):
    # label images equal up to permutation
    if ((a > 0) != (b > 0)).any():
        return False
    pa, pb = a[a > 0].astype(np.int64), b[b > 0].astype(np.int64)
    pairs = np.unique(np.stack([pa, pb]), axis=1)
    return len(np.unique(pairs[0])) == pairs.shape[1] == len(np.unique(pairs[1]))


print("== correctness (200 random masks per kind) ==")
bad = 0
for kind in ("blob", "noisy"):
    for i in range(200):
        n = int(rng.choice([32, 64, 100]))
        m = make_mask(n, kind)
        ref, nr = ndi.label(m, FG8)
        for name, l in [
            ("cv2_spaghetti", cv2.connectedComponentsWithAlgorithm(m.view(np.uint8), 8, cv2.CV_32S, cv2.CCL_SPAGHETTI)[1]),
            ("cv2_sauf", cv2.connectedComponentsWithAlgorithm(m.view(np.uint8), 8, cv2.CV_32S, cv2.CCL_SAUF)[1]),
            ("cc3d", cc3d.connected_components(m, connectivity=8)),
            ("skimage", sk_label(m, connectivity=2)),
        ]:
            if not same_partition(ref, l):
                bad += 1; print("MISMATCH", name, kind, n)
        _, nbg = ndi.label(~np.pad(m, 1), BG4)
        b1 = nbg - 1
        if chi8_bitquad(m) != nr - b1: bad += 1; print("bitquad mismatch")
        if sk_euler(m, connectivity=2) != nr - b1: bad += 1; print("skimage euler(conn=2) mismatch", sk_euler(m, 2), nr - b1)
        lab, nn, h = holes_per_component(m)
        _, _, h2 = holes_per_component_cv2(m)
        if not np.array_equal(np.sort(h2), np.sort(holes_bruteforce(m))): bad += 1; print('cv2 holes mismatch')
        if not np.array_equal(np.sort(h), np.sort(holes_bruteforce(m))) or h.sum() != b1:
            bad += 1; print("holes mismatch")
print("mismatches:", bad)

# skimage connectivity semantics on a diagonal pair and a 4-conn ring
diag = np.array([[1, 0], [0, 1]], bool)
ring = np.zeros((5, 5), bool); ring[1, 2] = ring[2, 1] = ring[2, 3] = ring[3, 2] = True  # diamond: 8-conn ring w/ 4-conn hole
print("diag: sk conn2", sk_euler(diag, 2), "conn1", sk_euler(diag, 1), "| diamond: sk conn2", sk_euler(ring, 2), "conn1", sk_euler(ring, 1), "bitquad8", chi8_bitquad(ring))

print("\n== timings, us per call, single thread ==")
hdr = f"{'n':>5} {'kind':>6} | {'scipy':>8} {'cv2_spag':>8} {'cv2_sauf':>8} {'cv2_bbdt':>8} {'cc3d':>8} {'skimage':>8} | {'bitquad':>8} {'skEuler':>8} {'holes/comp':>10} | {'sp_edt':>8} {'edt':>8} {'cv2_dt':>8}"
print(hdr)
for n in (64, 128, 256, 512, 1024):
    for kind in ("blob", "noisy"):
        m = make_mask(n, kind)
        u8 = m.view(np.uint8)
        reps = max(5, int(2e5 / n / n * 20))
        r = [
            timeit(lambda: ndi.label(m, FG8), reps),
            timeit(lambda: cv2.connectedComponentsWithAlgorithm(u8, 8, cv2.CV_32S, cv2.CCL_SPAGHETTI), reps),
            timeit(lambda: cv2.connectedComponentsWithAlgorithm(u8, 8, cv2.CV_32S, cv2.CCL_SAUF), reps),
            timeit(lambda: cv2.connectedComponentsWithAlgorithm(u8, 8, cv2.CV_32S, cv2.CCL_BBDT), reps),
            timeit(lambda: cc3d.connected_components(m, connectivity=8), reps),
            timeit(lambda: sk_label(m, connectivity=2), reps),
            timeit(lambda: chi8_bitquad(m), reps),
            timeit(lambda: sk_euler(m, 2), reps),
            timeit(lambda: holes_per_component_cv2(m), reps),
            timeit(lambda: ndi.distance_transform_edt(~m), reps),
            timeit(lambda: edt.edt(~m, black_border=False, parallel=1), reps),
            timeit(lambda: cv2.distanceTransform((~m).view(np.uint8), cv2.DIST_L2, cv2.DIST_MASK_PRECISE), reps),
        ]
        print(f"{n:>5} {kind:>6} | " + " ".join(f"{x:8.1f}" for x in r[:6]) + " | " + f"{r[6]:8.1f} {r[7]:8.1f} {r[8]:10.1f} | " + " ".join(f"{x:8.1f}" for x in r[9:]))

print("\n== Hungarian (linear_sum_assignment) us ==")
for k in (2, 6, 20, 100):
    C = rng.random((k, k))
    print(k, f"{timeit(lambda: linear_sum_assignment(C, maximize=True), 2000):.1f}")

print("\n== 3D cc3d 26-conn and 6-conn bg, 64^3 / 128^3 nodule-like volumes ==")
for n in (32, 64, 128):
    f = ndi.gaussian_filter(rng.standard_normal((n, n, n)), n / 10)
    zz, yy, xx = np.mgrid[:n, :n, :n]
    r = np.sqrt((zz - n / 2) ** 2 + (yy - n / 2) ** 2 + (xx - n / 2) ** 2) / (n / 2)
    v = (f / f.std() * 0.5 - 2 * r + 0.8) > 0
    reps = 20
    t26 = timeit(lambda: cc3d.connected_components(v, connectivity=26), reps)
    t6 = timeit(lambda: cc3d.connected_components(~np.pad(v, 1), connectivity=6), reps)
    tsp = timeit(lambda: ndi.label(v, np.ones((3, 3, 3), bool)), reps)
    te = timeit(lambda: sk_euler(v, 3), reps)
    print(n, f"cc3d26={t26:.0f}us cc3d6bg={t6:.0f}us scipy26={tsp:.0f}us skEuler3D={te:.0f}us vox={v.sum()}")

import numpy as np, cc3d, cv2

def holes_per_component(m):
    """FG 8-conn, BG 4-conn. Returns (fg labels, n, holes per component[1..n])."""
    mp = np.pad(m, 1)
    lab, n = cc3d.connected_components(mp, connectivity=8, return_N=True)
    bl = cc3d.connected_components(~mp, connectivity=4)
    outside = bl[0, 0]
    flat = bl.ravel()
    idx = np.flatnonzero((flat != 0) & (flat != outside))     # raster-ordered hole pixels
    holes = np.zeros(n + 1, np.int64)
    if idx.size:
        _, first = np.unique(flat[idx], return_index=True)    # topmost-leftmost pixel of each hole
        enc = lab.ravel()[idx[first] - mp.shape[1]]           # pixel above is FG of enclosing comp
        np.add.at(holes, enc, 1)
    return lab[1:-1, 1:-1], n, holes[1:]

def holes_per_component_cv2(m):
    mp = np.pad(m, 1).view(np.uint8)
    n, lab = cv2.connectedComponentsWithAlgorithm(mp, 8, cv2.CV_32S, cv2.CCL_SPAGHETTI)
    nb, bl = cv2.connectedComponentsWithAlgorithm(1 - mp, 4, cv2.CV_32S, cv2.CCL_SPAGHETTI)
    outside = bl[0, 0]
    flat = bl.ravel()
    idx = np.flatnonzero((flat != 0) & (flat != outside))
    holes = np.zeros(n, np.int64)
    if idx.size:
        _, first = np.unique(flat[idx], return_index=True)
        np.add.at(holes, lab.ravel()[idx[first] - mp.shape[1]], 1)
    return lab[1:-1, 1:-1], n - 1, holes[1:]

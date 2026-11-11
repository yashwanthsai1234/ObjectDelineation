"""LIDC-IDRI reader outlines -> per-slice multi-reader nodule masks.

Each CT series has (usually) four unblinded reading sessions. A nodule >= 3 mm is outlined slice by slice as pixel
contours on the 512 x 512 grid; ``inclusion=FALSE`` contours cut regions out; smaller nodules are a single point.
Rules (parity with pylidc on all 1,018 scans: 6,859 / 6,859 annotations matched, 99.97 % of reader-slice masks
exact):

- P1  Single-point contours inside a >= 3 mm nodule are dropped (pylidc's database drops them too); they rasterized
      as stray 1-px components.
- P2  One outline rule for every contour: a contour covers its filled polygon including the outline pixels.
      Inclusions add their polygon, then exclusions remove theirs (interior and outline), as pylidc applies them.
      LIDC readers drew the outline on the first pixels outside the nodule, so masks are about 1 px generous;
      this matters most for nodules under 50 px.
- P3  Slice index = nearest scan slice z position (pylidc database), not the spacing of the annotated z values,
      which is wrong for overlapping reconstructions (6 / 67 scans of a sample).
- P4  Patient ids from the pylidc database (the XML carries none); 1,018 series, 1,010 patients, 8 with two scans.

Reader annotations are clustered into nodules by shared voxels (z index, y, x); an item is one slice of one nodule
outlined by at least two readers, and holds one mask per reader who drew on that slice. Two more rules:

- P5  Copied outlines. Readers whose masks on a slice are identical count once (the lowest session number
      is kept; ``Slice.copied`` records the others). In the XML, 160 nodule pairs of two sessions in 144 session
      pairs share point-for-point identical contours, and 158 of 160 carry different characteristics, so these are
      outlines one reader adopted from another in the unblinded phase (Armato et al. 2007: readers saw the others'
      marks), not duplicated sessions. 392 of the 425 identical slice masks come from such identical contours, the
      other 33 are all >= 20 px; the reward treats raters as independent draws, so a copy would count twice.
- P6  Other nodules on the slice. ``Slice.others`` lists the box of every other nodule outlined on that
      slice, by any number of readers, so a consumer can refuse a crop that shows an annotated nodule absent from
      its ground truth (``crowded``).

References:
- LIDC-IDRI collection, known issues and citation (DOI 10.7937/K9/TCIA.2015.LO9QL9SX):
  https://www.cancerimagingarchive.net/collection/lidc-idri/
- Armato et al., "The LIDC/IDRI ...: a completed reference database of lung nodules on CT scans", Med. Phys. 38(2),
  915-931 (2011), https://doi.org/10.1118/1.3528204
- Outline on "those pixels that were just outside" the nodule: Biancardi, Reeves et al., Proc. SPIE 6915 (2008),
  https://www.via.cornell.edu/publications/ABAR08a.pdf
- Armato et al., "The Lung Image Database Consortium (LIDC): ensuring the integrity of expert-defined truth",
  Acad. Radiol. 14(12), 1455-1463 (2007), https://pmc.ncbi.nlm.nih.gov/articles/PMC2151472 (unblinded phase:
  readers review their own marks with the others' anonymised marks; no forced consensus)
- DICOM PS3.5 §9.1 UID encoding rules (digits and dots, no leading zeros, <= 64 characters), checked on every
  series uid read from the XML: https://dicom.nema.org/medical/dicom/current/output/chtml/part05/chapter_9.html
- pylidc (Hancock & Magnan 2016) Annotation.boolean_mask and Scan.slice_zvals, the parity reference:
  https://pylidc.github.io/annotation.html , https://github.com/notmatthancock/pylidc/blob/master/pylidc/Annotation.py
"""
import importlib.util
import re
import sqlite3
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw

from objdel.geom.box import boxes_meet

from .store import Item

NS = {"n": "http://www.nih.gov"}
GRID = 512
MARGIN = 8  # crop margin of Slice.stack around the readers' union
# DICOM PS3.5 §9.1: UID = <org root>.<suffix>, components of digits 0-9 ([0-9], not \d, which also matches other
# scripts' digits) without a leading zero unless the component is "0", separated by ".", at most 64 characters.
UID = re.compile(r"(0|[1-9][0-9]*)(\.(0|[1-9][0-9]*))+")
UID_MAX = 64


@dataclass(frozen=True)
class Roi:
    """One reader's outline on one slice: z position, inclusion flag and (x, y) pixel vertices."""
    z: float
    inclusion: bool
    points: tuple[tuple[int, int], ...]  # (x, y) pixel coordinates


@dataclass(frozen=True)
class Scan:
    """Series geometry: patient id, ascending slice z positions (mm), pixel spacing and slice thickness."""
    patient: str
    zvals: np.ndarray  # slice z positions, ascending (mm)
    pixel_spacing: float
    slice_thickness: float


@dataclass
class Slice:
    """One multi-reader slice of one nodule: masks cropped to the readers' union plus MARGIN."""

    cluster: int
    k: int
    readers: list[int]  # 1-based reading-session numbers, one per stack layer
    cluster_readers: list[int]
    dup_reader: bool  # some reader has two annotations merged into this nodule
    y0: int
    x0: int
    stack: np.ndarray  # bool (len(readers), h, w)
    copied: dict[int, int]  # reader -> the kept reader whose identical mask it repeats on this slice
    others: list[tuple[int, int, int, int]]  # (x0, y0, x1, y1) half-open box of each other nodule on the slice

    def full(self) -> np.ndarray:
        """The reader masks as a bool (readers, GRID, GRID) stack: the cropped stack placed at (y0, x0) on an all-False
        512 x 512 grid.
        """
        out = np.zeros((len(self.readers), GRID, GRID), bool)
        out[:, self.y0:self.y0 + self.stack.shape[1], self.x0:self.x0 + self.stack.shape[2]] = self.stack
        return out


def _text(el, tag):
    t = el.find(f"n:{tag}", NS)
    return t.text.strip() if t is not None and t.text else None


def uid(s) -> str:
    """s if it is a DICOM UID (PS3.5 §9.1), else ValueError. Series UIDs name files and directories and go into the
    TCIA URL, and the XML is untrusted: this keeps "../..", absolute paths and query strings out."""
    if not (isinstance(s, str) and len(s) <= UID_MAX and UID.fullmatch(s)):
        raise ValueError(f"not a DICOM UID (digits and dots, <= {UID_MAX} characters): {s!r}")
    return s


def _series(root, path) -> str:
    """The message's series uid, the one place a series uid enters from the XML (every reader goes through it)."""
    s = _text(root.find("n:ResponseHeader", NS), "SeriesInstanceUid")
    try:
        return uid(s)
    except ValueError as e:
        raise ValueError(f"{path}: SeriesInstanceUid {e}") from None


def _messages(xml_root):
    """(path, series uid, root) of every LIDC read message under xml_root, by path (radiograph reads skipped)."""
    for p in sorted(Path(xml_root).rglob("*.xml")):
        if (root := ET.parse(p).getroot()).tag.endswith("LidcReadMessage"):
            yield p, _series(root, p), root


def parse(path):
    """(series uid, sessions) or None for a non-CT message. sessions: per reading session, its nodules >= 3 mm
    (contoured and characterised), each a list of Roi; single-point contours of those nodules dropped."""
    root = ET.parse(path).getroot()
    if not root.tag.endswith("LidcReadMessage"):
        return None
    sessions = []
    for s in root.findall("n:readingSession", NS):
        nodules = []
        for nod in s.findall("n:unblindedReadNodule", NS):
            rois = [Roi(float(_text(roi, "imageZposition")), _text(roi, "inclusion") == "TRUE",
                        tuple((int(e.find("n:xCoord", NS).text), int(e.find("n:yCoord", NS).text))
                              for e in roi.findall("n:edgeMap", NS)))
                    for roi in nod.findall("n:roi", NS)]
            rois = [r for r in rois if len(r.points) > 1]
            if rois and nod.find("n:characteristics", NS) is not None:
                nodules.append(rois)
        sessions.append(nodules)
    return _series(root, path), sessions


def draw(points) -> np.ndarray:
    """The closed polygon through `points`, outline pixels included."""
    img = Image.new("1", (GRID, GRID), 0)
    ImageDraw.Draw(img).polygon(list(points), fill=1, outline=1)
    # Not np.asarray(img, dtype=bool): for mode "1" that keeps byte value 255 for True, which code that reads the
    # bytes (OpenCV labeling in objdel.geom.topology) misreads.
    return np.asarray(img, dtype=np.uint8) > 0


def slice_index(zvals, z) -> int:
    """Index of the scan slice nearest to z."""
    return int(np.abs(np.asarray(zvals) - z).argmin())


def rasterize(rois, zvals) -> dict[int, np.ndarray]:
    """{slice index: mask} of one reader's nodule: inclusions first, then exclusions on slices that have one."""
    masks = {}
    for r in rois:
        if r.inclusion:
            k = slice_index(zvals, r.z)
            masks[k] = masks[k] | draw(r.points) if k in masks else draw(r.points)
    for r in rois:
        if not r.inclusion and (k := slice_index(zvals, r.z)) in masks:
            masks[k] = masks[k] & ~draw(r.points)
    return {k: m for k, m in masks.items() if m.any()}


def clusters(sessions, zvals) -> list[tuple[dict[int, dict[int, np.ndarray]], bool]]:
    """Nodules across readers: annotations sharing a voxel are one nodule (union-find). Each nodule is
    ({0-based reader: {slice index: mask}}, dup_reader): a reader's annotations in one nodule are merged, and
    dup_reader says that happened. Ordered by first slice, then by first annotation."""
    anns = [(r, m) for r, nods in enumerate(sessions) for rois in nods if (m := rasterize(rois, zvals))]
    keys, owner = [], []
    for i, (_, m) in enumerate(anns):
        for k, mask in m.items():
            keys.append(k * GRID * GRID + np.flatnonzero(mask))
            owner.append(np.full(keys[-1].size, i))
    parent = list(range(len(anns)))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    if anns:
        keys, owner = np.concatenate(keys), np.concatenate(owner)
        order = np.argsort(keys, kind="stable")
        keys, owner = keys[order], owner[order]
        same = keys[1:] == keys[:-1]
        for a, b in zip(owner[:-1][same], owner[1:][same], strict=True):
            parent[find(int(a))] = find(int(b))
    groups = {}
    for i in range(len(anns)):
        groups.setdefault(find(i), []).append(i)
    out = []
    for members in sorted(groups.values(), key=lambda g: min(min(anns[i][1]) for i in g)):
        per_reader = {}
        for i in members:
            r, m = anns[i]
            pr = per_reader.setdefault(r, {})
            for k, mask in m.items():
                pr[k] = pr[k] | mask if k in pr else mask
        out.append((dict(sorted(per_reader.items())), len(members) > len(per_reader)))
    return out


def _box(mask) -> tuple[int, int, int, int]:
    ys, xs = np.nonzero(mask)
    return int(xs.min()), int(ys.min()), int(xs.max()) + 1, int(ys.max()) + 1


def slices(sessions, zvals) -> list[Slice]:
    """Every slice of every nodule on which at least two readers drew distinct masks, with the boxes of the
    other nodules outlined on that slice."""
    cl = clusters(sessions, zvals)
    boxes = {}  # slice index -> [(cluster, box)]
    for ci, (per_reader, _) in enumerate(cl):
        for k in sorted({k for m in per_reader.values() for k in m}):
            boxes.setdefault(k, []).append((ci, _box(np.logical_or.reduce([m[k] for m in per_reader.values()
                                                                           if k in m]))))
    out = []
    for ci, (per_reader, dup) in enumerate(cl):
        if len(per_reader) < 2:
            continue
        for k in sorted({k for m in per_reader.values() for k in m}):
            readers, copied = [], {}
            for r in (r for r, m in per_reader.items() if k in m):
                twin = next((q for q in readers if np.array_equal(per_reader[q][k], per_reader[r][k])), None)
                if twin is None:
                    readers.append(r)
                else:
                    copied[r + 1] = twin + 1
            if len(readers) < 2:
                continue
            stack = np.stack([per_reader[r][k] for r in readers])
            ys, xs = np.nonzero(stack.any(0))
            y0, y1 = max(0, ys.min() - MARGIN), min(GRID, ys.max() + MARGIN + 1)
            x0, x1 = max(0, xs.min() - MARGIN), min(GRID, xs.max() + MARGIN + 1)
            out.append(Slice(ci, k, [r + 1 for r in readers], [r + 1 for r in per_reader], dup, int(y0), int(x0),
                             stack[:, y0:y1, x0:x1], copied, [b for cj, b in boxes[k] if cj != ci]))
    return out


def crowded(s: Slice, box) -> bool:
    """Does the half-open window box = (x0, y0, x1, y1) overlap another nodule's box on the slice?"""
    return boxes_meet(box, s.others)


def series_files(xml_root) -> dict[str, Path]:
    """Series uid -> its XML. A series in several files (byte-identical duplicates, or a resubmitted correction):
    the correction file, else the first path in sorted order."""
    by = {}
    for p, series, _ in _messages(xml_root):
        by.setdefault(series, []).append(p)
    return {s: next((p for p in ps if "correction" in p.name), ps[0]) for s, ps in by.items()}


def sop_z(xml_root) -> dict[str, dict[str, float]]:
    """Series uid -> {imageSOP_UID: imageZposition rounded to 0.01 mm} over every ROI of every XML of the series
    (what the CT fetch needs to pick slices; experiments/phase1/fetch_lidc_pixels)."""
    out = {}
    for _, series, root in _messages(xml_root):
        m = out.setdefault(series, {})
        for roi in root.iter(f"{{{NS['n']}}}roi"):
            m[_text(roi, "imageSOP_UID")] = round(float(_text(roi, "imageZposition")), 2)
    return out


def pylidc_db() -> Path:
    """The SQLite database bundled with pylidc 0.2.3, located without importing pylidc (whose import needs
    pkg_resources and numpy < 2)."""
    return Path(importlib.util.find_spec("pylidc").submodule_search_locations[0]) / "pylidc.sqlite"


def scans(db=None) -> dict[str, Scan]:
    """Series uid -> Scan (patient id, slice z positions, spacing), read-only from the pylidc database."""
    con = sqlite3.connect(f"file:{db or pylidc_db()}?mode=ro", uri=True)
    try:
        z = {}
        for sid, val in con.execute("select scan_id, val from zvals"):
            z.setdefault(sid, []).append(val)
        return {series: Scan(patient, np.sort(np.asarray(z.get(sid, []), float)), float(ps), float(st))
                for sid, series, patient, ps, st in con.execute(
                    "select id, series_instance_uid, patient_id, pixel_spacing, slice_thickness from scans")}
    finally:
        con.close()


def ct_slices(npz_path, zvals) -> dict[int, np.ndarray]:
    """{slice index: int16 Hounsfield units} from a slice archive (keys "z<mm, 2 decimals>"); keys that are
    not within 0.01 mm of a scan slice are ignored."""
    out = {}
    with np.load(npz_path, allow_pickle=False) as z:
        for key in z.files:
            k = slice_index(zvals, float(key[1:]))
            if abs(zvals[k] - float(key[1:])) < 0.01:
                out[k] = z[key]
    return out


def item(series: str, scan: Scan, s: Slice, split: str, image, sessions: int, xml: str) -> Item:
    """The store `Item` of one multi-reader slice: the pixels, one "nodule" mask per reader on the 512 x 512 grid and
    metadata.
    """
    return Item(item_id=f"lidc/{series}/n{s.cluster:03d}/z{s.k:03d}", dataset="lidc", group_id=scan.patient,
                split=split, image=image, masks={"nodule": s.full()}, raters=[str(r) for r in s.readers],
                meta={"series": series, "patient": scan.patient, "cluster": s.cluster, "slice": s.k,
                      "z_mm": float(scan.zvals[s.k]), "readers_in_cluster": [str(r) for r in s.cluster_readers],
                      "dup_reader_in_cluster": s.dup_reader,
                      "sessions": sessions, "xml": xml, "pixel_spacing_mm": scan.pixel_spacing,
                      "slice_thickness_mm": scan.slice_thickness,
                      "copied_raters": {str(r): str(q) for r, q in s.copied.items()},
                      "other_nodules": [list(b) for b in s.others]})

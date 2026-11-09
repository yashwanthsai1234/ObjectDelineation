"""Storage-format prototype: a small Parquet dataset (Hugging Face `datasets`).

One row per item (RIGA image, or LIDC lesion slice crop). Columns: id, dataset,
subset, group_id (patient / source image), split, image (PNG bytes, lossless;
16-bit PNG of HU+1024 for CT), image_sha256 (decoded pixels), per-rater masks
as exact polygons (or RLE in the comparison variant), meta (JSON string).
Writes <out>/<variant>/data.parquet + manifest.json (SHA-256 of files, rows,
images and masks), then reloads and checks bit-exact round trips and speed.

    python -I build_store.py <out_dir> <riga_extract_out> <splits_dir> <series_patient.csv> [n_riga_per_subset] [n_lidc_series]
"""
import csv
import hashlib
import io
import json
import sys
import time
from collections import defaultdict
from pathlib import Path

REPO = Path("/Users/sreeyashwanthsaivenkatesh/Workspace/projects/CUDA/ObjectDelineation")
sys.path[:0] = [str(REPO), str(REPO / "src"), str(Path(__file__).parent)]

import numpy as np  # noqa: E402
from PIL import Image  # noqa: E402

from maskcodec import mask_sha, poly_decode, poly_encode, rle_decode, rle_encode  # noqa: E402

RAW = REPO / "data" / "raw"


def png_bytes(arr):
    buf = io.BytesIO()
    Image.fromarray(arr).save(buf, format="PNG", compress_level=6)
    return buf.getvalue()


def pix_sha(arr):
    return hashlib.sha256(np.ascontiguousarray(arr).tobytes() + f"{arr.shape}{arr.dtype}".encode()).hexdigest()


def riga_items(extract_out, splits_dir, per_subset):
    split = {r["image_id"]: r for r in csv.DictReader(open(Path(splits_dir) / "riga_images.csv"))}
    ann = defaultdict(dict)
    for r in csv.DictReader(open(Path(extract_out) / "annotations.csv")):
        ann[(r["subset"], int(r["image"]))][int(r["expert"])] = r
    by_sub = defaultdict(list)
    for r in csv.DictReader(open(Path(extract_out) / "images.csv")):
        by_sub[r["subset"]].append(r)
    for sub in sorted(by_sub):
        for r in sorted(by_sub[sub], key=lambda r: int(r["image"]))[:per_subset]:
            idx = int(r["image"])
            stem = f"{sub.replace('/', '_')}_image{idx}"
            img = np.asarray(Image.open(Path(extract_out) / "refs" / f"{stem}.png").convert("RGB"))
            raters, src = [], {}
            for e in range(1, 7):
                p = Path(extract_out) / "masks" / f"{stem}-{e}.png"
                if p.exists():
                    v = np.asarray(Image.open(p))
                    src[e] = {"disc": v >= 100, "cup": v >= 200}
                raters.append((e, src.get(e), ann[(sub, idx)].get(e, {}).get("stage", "")))
            iid = f"{sub}/image{idx}"
            s = split[iid]
            yield {"id": f"riga/{iid}", "dataset": "riga", "subset": sub, "group_id": s["group_id"], "split": s["split"],
                   "img": img, "raters": raters,
                   "meta": {"ref": r["ref"], "stage": {e: st for e, _, st in raters}}}


def lidc_items(splits_dir, series_patient, n_series, crop=128):
    from experiments.phase0.lidc_census import nodule_clusters
    sp = {r["series_uid"]: r for r in csv.DictReader(open(series_patient))}
    split = {r["patient_id"]: r["split"] for r in csv.DictReader(open(Path(splits_dir) / "lidc_patients.csv"))}
    slices_dir = REPO / "data" / "lidc_slices"
    series = sorted(p.stem for p in slices_dir.glob("*.npz"))[:n_series]
    for s in series:
        xml = RAW / "lidc_xml" / sp[s]["xml_file"]
        _, sessions, zi, clusters, _ = nodule_clusters(xml)
        z_of = {}
        for z, k in zi.items():
            z_of.setdefault(k, z)
        with np.load(slices_dir / f"{s}.npz") as npz:
            hu = {k: npz[k] for k in npz.files}
        for ci, cl in enumerate(clusters):
            z0, y0, x0, z1, y1, x1 = cl["box"]
            cy, cx = (y0 + y1) // 2, (x0 + x1) // 2
            oy, ox = int(np.clip(cy - crop // 2, 0, 512 - crop)), int(np.clip(cx - crop // 2, 0, 512 - crop))
            for k in sorted({k for pr in cl["per_reader"].values() for k in pr}):
                key = f"z{z_of[k]:.2f}"
                if key not in hu:
                    continue
                img = (hu[key][oy:oy + crop, ox:ox + crop].astype(np.int32) + 1024).clip(0, 65535).astype(np.uint16)
                raters = []
                for r in range(len(sessions)):
                    m = cl["per_reader"].get(r, {}).get(k)
                    m = m[oy:oy + crop, ox:ox + crop] if m is not None else None
                    raters.append((r + 1, {"nodule": m} if m is not None else None, "outlined" if r in cl["readers"] else "none"))
                pid = sp[s]["patient_id"]
                yield {"id": f"lidc/{s}/c{ci}/{key}", "dataset": "lidc", "subset": "LIDC-IDRI", "group_id": pid,
                       "split": split[pid], "img": img, "raters": raters,
                       "meta": {"series": s, "cluster": ci, "z": z_of[k], "crop_origin_yx": [oy, ox], "hu_offset": 1024,
                                "readers_outlined_cluster": cl["readers"]}}


def to_row(it, codec):
    enc = poly_encode if codec == "poly" else rle_encode
    raters, shas = [], {}
    for rid, ms, note in it["raters"]:
        structs = []
        if ms:
            for name, m in ms.items():
                e = enc(m)
                structs.append({"name": name, "h": e["h"], "w": e["w"],
                                "rings": e.get("rings", []), "counts": e.get("counts", [])})
                shas[f"{rid}/{name}"] = mask_sha(m)
        raters.append({"rater_id": rid, "present": ms is not None, "structures": structs})
    return {"id": it["id"], "dataset": it["dataset"], "subset": it["subset"], "group_id": it["group_id"],
            "split": it["split"], "image": {"bytes": png_bytes(it["img"]), "path": None},
            "image_sha256": pix_sha(it["img"]), "h": int(it["img"].shape[0]), "w": int(it["img"].shape[1]),
            "raters": raters, "meta": json.dumps({**it["meta"], "mask_sha256": shas}, sort_keys=True, default=str)}


def features():
    from datasets import Features, Image as HFImage, List, Value
    struct = {"name": Value("string"), "h": Value("int32"), "w": Value("int32"),
              "rings": List(List(Value("int32"))), "counts": List(Value("uint32"))}
    return Features({"id": Value("string"), "dataset": Value("string"), "subset": Value("string"),
                     "group_id": Value("string"), "split": Value("string"), "image": HFImage(),
                     "image_sha256": Value("string"), "h": Value("int32"), "w": Value("int32"),
                     "raters": List({"rater_id": Value("int8"), "present": Value("bool"), "structures": List(struct)}),
                     "meta": Value("string")})


def file_sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def verify(path, codec):
    from datasets import Dataset
    t0 = time.perf_counter()
    ds = Dataset.from_parquet(str(path), keep_in_memory=True)
    t_open = time.perf_counter() - t0
    dec = poly_decode if codec == "poly" else rle_decode
    n_masks = bad_masks = bad_imgs = 0
    per = defaultdict(lambda: {"rows": 0, "masks": 0, "fetch_s": 0.0, "img_s": 0.0, "mask_s": 0.0})
    t0 = time.perf_counter()
    it = iter(ds)
    while True:
        a = time.perf_counter()
        try:
            row = next(it)
        except StopIteration:
            break
        b = time.perf_counter()
        img = np.asarray(row["image"])
        c = time.perf_counter()
        st_ = per[row["dataset"]]
        st_["rows"] += 1
        st_["fetch_s"] += b - a
        st_["img_s"] += c - b
        bad_imgs += pix_sha(img) != row["image_sha256"]
        shas = json.loads(row["meta"])["mask_sha256"]
        for r in row["raters"]:
            for st in r["structures"]:
                a = time.perf_counter()
                m = dec({"h": st["h"], "w": st["w"], "rings": st["rings"], "counts": st["counts"]})
                st_["mask_s"] += time.perf_counter() - a
                st_["masks"] += 1
                n_masks += 1
                bad_masks += mask_sha(m) != shas[f"{r['rater_id']}/{st['name']}"]
    t_all = time.perf_counter() - t0
    timing = {k: {"rows": v["rows"], "fetch_ms_per_row": round(1000 * v["fetch_s"] / v["rows"], 2),
                  "image_decode_ms_per_row": round(1000 * v["img_s"] / v["rows"], 2),
                  "mask_decode_ms_per_mask": round(1000 * v["mask_s"] / max(1, v["masks"]), 3),
                  "rows_per_s_incl_masks": round(v["rows"] / (v["fetch_s"] + v["img_s"] + v["mask_s"]), 1)}
              for k, v in per.items()}
    return {"rows": len(ds), "masks": n_masks, "mask_mismatch": bad_masks, "image_mismatch": bad_imgs,
            "open_s": round(t_open, 3), "iterate_all_s": round(t_all, 3), "timing": timing}


def main(out, extract_out, splits_dir, series_patient, per_subset=8, n_series=12):
    from datasets import Dataset
    out = Path(out)
    items = list(riga_items(extract_out, splits_dir, int(per_subset))) + list(lidc_items(splits_dir, series_patient, int(n_series)))
    report = {}
    for codec in ("poly", "rle"):
        d = out / codec
        d.mkdir(parents=True, exist_ok=True)
        rows = [to_row(it, codec) for it in items]
        ds = Dataset.from_list(rows, features=features())
        pq = d / "data.parquet"
        ds.to_parquet(str(pq))
        manifest = {"files": {"data.parquet": file_sha(pq)}, "rows": len(rows),
                    "row_sha256": {r["id"]: hashlib.sha256((r["image_sha256"] + r["meta"]).encode()).hexdigest() for r in rows}}
        (d / "manifest.json").write_text(json.dumps(manifest, indent=1, sort_keys=True))
        size = pq.stat().st_size
        by_ds = {}
        for name in ("riga", "lidc"):
            rs = [r for r in rows if r["dataset"] == name]
            img_b = sum(len(r["image"]["bytes"]) for r in rs)
            mask_b = sum(len(json.dumps(r["raters"])) for r in rs)
            nm = sum(len(st["rings"]) * 0 + 1 for r in rs for ra in r["raters"] for st in ra["structures"])
            enc_b = sum(4 * (sum(len(x) for x in st["rings"]) + len(st["counts"])) for r in rs for ra in r["raters"] for st in ra["structures"])
            by_ds[name] = {"rows": len(rs), "png_bytes_per_row": round(img_b / max(1, len(rs))), "masks": nm,
                           "mask_raw_int32_bytes_per_mask": round(enc_b / max(1, nm)), "mask_json_bytes_per_row": round(mask_b / max(1, len(rs)))}
        report[codec] = {"parquet_bytes": size, "bytes_per_row": round(size / len(rows)), "by_dataset": by_ds, "verify": verify(pq, codec)}
    # Parquet column sizes for the poly variant (compressed), via pyarrow metadata.
    import pyarrow.parquet as papq
    md = papq.ParquetFile(out / "poly" / "data.parquet").metadata
    cols = defaultdict(int)
    for g in range(md.num_row_groups):
        for c in range(md.num_columns):
            col = md.row_group(g).column(c)
            cols[col.path_in_schema.split(".")[0]] += col.total_compressed_size
    report["poly_column_bytes"] = dict(cols)
    md = papq.ParquetFile(out / "rle" / "data.parquet").metadata
    cols = defaultdict(int)
    for g in range(md.num_row_groups):
        for c in range(md.num_columns):
            col = md.row_group(g).column(c)
            cols[col.path_in_schema.split(".")[0]] += col.total_compressed_size
    report["rle_column_bytes"] = dict(cols)
    import datasets
    import pyarrow
    report["versions"] = {"datasets": datasets.__version__, "pyarrow": pyarrow.__version__}
    (out / "report.json").write_text(json.dumps(report, indent=1))
    print(json.dumps(report, indent=1))


if __name__ == "__main__":
    main(*sys.argv[1:])

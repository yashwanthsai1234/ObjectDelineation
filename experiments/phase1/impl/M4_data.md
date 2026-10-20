# M4: data (RIGA extractor, LIDC parser, splits, store)

Plan: D13, §4.5, gates K10 (data part) and K11. Promoted from `lab/t_data` (riga, lidc_full, lidc_parity, splits,
storage). TRL chat-row export moved to M7 (EXECUTION.md §2).

**Code:**
- `src/objdel/data/__init__.py` (exports `Item`)
- `src/objdel/data/riga.py`: staged extractor F1–F7, certification, near-duplicate descriptor, `item()`
- `src/objdel/data/lidc.py`: XML parser with fixes P1–P4, clustering, slice records, pylidc DB reader, `item()`
- `src/objdel/data/splits.py`: seeded group-aware splits, near-duplicate search, leak checks
- `src/objdel/data/store.py`: `Item`, `write`, `load`, `verify`, mask codec
- `tests/test_data_riga.py`, `tests/test_data_lidc.py`, `tests/test_data_splits.py`, `tests/test_data_store.py`
- `experiments/phase1/build_data.py`

## Result

| Gate / check | Result |
|---|---|
| RIGA masks | 750 photographs × 6 raters = **4,500 / 4,500** disc + cup pairs, all certified (disc and cup b0 = 1, b1 = 0, cup ⊆ disc, rim certified). Stages: A 4,456, B 41, C 3. References: prime 699, F1 (BinRushed1) 50, F7 (Magrabia F image48) 1. |
| RIGA vs lab extractor | 4,500 / 4,500 masks bit-exact (sandbox, lab extractor re-run) |
| RIGA vs Phase 0 | stage-A masks equal `data/riga_true` (data test on the val split) |
| RIGA duplicates | 37 near-duplicate pairs (31 same photo, 6 same eye) → 37 groups, all within one split; Magrabia test = 95 files, **88 unique photos** |
| RIGA split | train 591 / val 64 / test 95 (plan: 590 / 65 / 95; the val carve-out is per subset within ±1 by construction, and the group names, hence the hash order, differ from the lab's) |
| LIDC items | 1,018 series, 1,010 patients, **10,921** multi-reader slices; **10,921 / 10,921 bit-exact** vs the lab fixture `lidc_full.npz` (crop boxes, reader ids, session counts equal) |
| LIDC vs pylidc | 25-scan sample: annotation sets equal, 848 / 848 annotation-slice masks exact vs a pylidc-rule oracle (lab, all scans: 99.97 % reader-slice exact) |
| LIDC split | patients 710 / 99 / 201; items 7,623 / 1,085 / 2,213; 5-px disagreement slices 319 / 46 / 91 (70.0 / 10.1 / 20.0 %) |
| Plan numbers on the lab split | train disagreement slices: **313** (5-px filter); 308 under the size rule; **251** under the size rule without clusters where one reader's two annotations merged (how the lab's reward fixture selected its 251 R4 items) |
| LIDC CT pixels | 4,084 / 10,921 items have pixels (214 series in `data/lidc_slices`): train 2,557, val 483, test 1,044 |
| **K11** | **PASS**: 0 groups across splits (RIGA groups, LIDC patients incl. two-scan patients), 0 cross-split near-duplicate pairs (max cross-split NCC 0.282; largest NCC below the 0.45 threshold 0.326), 0 identical CT slices across splits |
| **K10** (store) | **PASS**: three full builds (2 and 3 workers) give identical manifests; load → write reproduces the file hash (data test) |
| Build | `build_data` 505 s with 2 workers (468 s with 3); peak total RSS 1.85 GB (main + 2 workers) |
| Store | 3.3 GB: riga 2.36 GB (train 1,696 MB, val 177 MB, test 487 MB), lidc 1.09 GB (train 687 MB, val 121 MB, test 283 MB) |
| Tests | `test_data_*.py`: 58 pass (49 unit, 9 data-backed, none skipped with the store present), 76 s; `make lint` clean; full `make test`: 364 passed, 1 failure in `test_reward_gates.py` (M6, not M4) |

SHA-256 (manifest of the final build, `data/store/manifest.json`):

| File | Rows | SHA-256 |
|---|---|---|
| riga/train.parquet | 591 | 83fb4f0ff59b2c517aca7ce838c7f0f6c9614fe74e4a90cf938c2c3b6699c077 |
| riga/val.parquet | 64 | a39bbe6f4e6cee37178dadbfd32a2168ffe0bb4b433dffa183ad094c81c2571a |
| riga/test.parquet | 95 | bcbbd7931536b6558859592e23f5340dab36d82de60dc7977f7934aa0fbcb010 |
| lidc/train.parquet | 7,623 | 0fdd68015df7a89340e1838bc393166562555c21b2da999874809c83737f51ab |
| lidc/val.parquet | 1,085 | 2f1cb5298f68b80d5ccc3a76344f7cdb2c66b1e1c8a3fc37256350ba5d64ad6e |
| lidc/test.parquet | 2,213 | 5dbb592cdb39a895701b0455edf8b3b34f143d19dbf2eba32652603d02fd1611 |

**Superseded by** the I9 CT fetch and `fixes/F4.md` (break-it B3): final store manifest SHA-256 `37fe7665…`;
RIGA 592 / 63 / 95 (keypoint duplicate groups, B3-8; 45 groups, 51 pairs; BinRushed1 `image=None`, B3-1); LIDC
10,625 items = 7,341 / 1,084 / 2,200, patients 709 / 100 / 201 (copied readers merged, B3-3: −80 slices; crowded
frames dropped, B3-2: −216), every item with pixels (B3-4); lab-split disagreement 296 (5 px) / 291 (rule) / 242
(rule, no merged-reader clusters). Per-file hashes: `data/store/manifest.json`, `fixes/F4.md`.

## Contract (EXECUTION.md §3, as delivered)

```python
from objdel.data import Item, store, riga, lidc, splits
Item(item_id, dataset, group_id, split, image, masks, raters, meta)   # dataclass, mutable
store.write(items, root) -> manifest     # items ascending by item_id within each (dataset, split)
store.load(dataset, split, root) -> Iterator[Item]   # verifies the file's SHA-256 and row count first
store.verify(root) -> manifest           # ValueError on any mismatch
```

Changes from the §3 sketch, with the evidence above:
- `image` is `np.ndarray | None`: uint8 (H, W, 3) RGB for RIGA; **int16 (H, W) Hounsfield units** for LIDC (lossless;
  the window to uint8 is a model-input choice for M5); `None` for the 6,837 LIDC items without local CT pixels.
  *(Superseded: after the I9 fetch every LIDC item has pixels; `None` now marks the RIGA BinRushed1
  reconstructions, B3-1. dtype/layout checked on write and load, B3-7.)*
- Frames are the originals: RIGA full photograph (W × H from 1,440 × 960 to 2,743 × 1,936 px), LIDC 512 × 512. `masks` share the
  image frame: RIGA `{"disc", "cup"}` (6, H, W); LIDC `{"nodule"}` (k, 512, 512), k = 2–4 readers on that slice.
- `item_id`: `riga/<subset>/image<NNN>`, `lidc/<series uid>/n<cluster NNN>/z<slice NNN>`. `group_id`:
  `riga:<smallest image id of its duplicate group>`; the LIDC patient id.
- `raters`: RIGA "1"–"6"; LIDC reading-session numbers "1"–"4" (not comparable across series).
- `meta` (JSON): RIGA `subset, source, index, reference (prime | maxrg6 | maxrg6:prime_is_rater1_copy), stages
  {rater: A | B:<ref>/b<k> | C:<ref>/b<k>}, duplicates [image ids]`; LIDC `series, patient, cluster, slice, z_mm,
  readers_in_cluster, dup_reader_in_cluster, sessions, xml, pixel_spacing_mm, slice_thickness_mm` (+
  `copied_raters`, `other_nodules` since F4: B3-3, B3-2). Later API additions: `lidc.sop_z`, `lidc.crowded`,
  `riga.keypoints`, `splits.keypoint_duplicates` (EXECUTION.md §3).

## Research

**LIDC XML and pylidc.**
- TCIA collection page ([LIDC-IDRI](https://www.cancerimagingarchive.net/collection/lidc-idri/)): 8 patients have
  two time-point scans ("users may use either or both"); 044.xml / 191.xml had a ≥ 3 mm nodule without
  characteristics (corrected 2018); XMLs "do not store radiologist annotations in a manner that allows for a
  comparison of individual radiologist reads across cases" (so reader ids are per-series session numbers only);
  citation DOI 10.7937/K9/TCIA.2015.LO9QL9SX.
- Outline semantics: per the LIDC specification the outline is "made of those pixels that were just outside the
  largest region that could be linked to the nodule" (Biancardi, Reeves et al., SPIE 6915, 2008,
  [PDF](https://www.via.cornell.edu/publications/ABAR08a.pdf)). The plan fixed "fill including the outline", so
  masks are about 1 px generous; documented in `lidc.py`.
- pylidc `Annotation.boolean_mask` (read from the installed 0.2.3 source and
  [GitHub](https://github.com/notmatthancock/pylidc/blob/master/pylidc/Annotation.py)): inclusion contours are
  closed, tested with `matplotlib.path.Path.contains_points` on integer pixel coordinates and OR-ed; with
  `include_contour_points=False` their points are then cleared; exclusions are applied after all inclusions.
  `Scan.slice_zvals` is the sorted `zvals` table of the bundled SQLite DB
  ([docs](https://pylidc.github.io/annotation.html)). pylidc 0.2.3 cannot be imported here (needs
  `pkg_resources` and numpy < 2), so `lidc.scans()` reads `pylidc.sqlite` directly with `sqlite3`, located with
  `importlib.util.find_spec` (no import).
- The parity oracle in `tests/test_data_lidc.py` re-implements pylidc's rasterizer (matplotlib Path) under rule P2.

**RIGA.** Almazroa et al., SPIE 10579 (2018), [DOI](https://doi.org/10.1117/12.2293584); six ophthalmologists
drew on copies of each photograph; licence forbids redistribution
([Deep Blue](https://deepblue.lib.umich.edu/data/concern/data_sets/3b591905z)). No published extraction audit
was found; the lab's inventory and F1–F7 are the evidence. Literature split (train BinRushed + MESSIDOR 655, test
Magrabia 95): MRNet, TAB, MSE-Nets (see the data-engineering research note).

**Hysteresis.** scikit-image `apply_hysteresis_threshold` (0.26 source): strict `>` for both thresholds,
`ndi.label` with the default (4-connected) structure. `riga.hysteresis` re-implements it in 5 lines of scipy so
scikit-image is not a dependency; a Hypothesis test checks equality against scikit-image when installed.

**Leakage.** Patient-level LIDC splits are the norm because adjacent slices are near-duplicates (ProSeg,
[arXiv 2512.00748](https://arxiv.org/abs/2512.00748)); slice-mixed splits inflate CT segmentation scores
(CTSCAN, [arXiv 2604.15561](https://arxiv.org/abs/2604.15561)). For RIGA, global perceptual hashes are dominated by
the shared fundus layout (lab: thousands of unrelated pairs within pHash 10), so the near-duplicate test is the
lab's vessel-texture NCC (high-pass central crop, shift-tolerant), threshold 0.45 in the empty gap between
same-eye (≥ 0.49) and unrelated (≤ 0.33) pairs.

**Parquet.** [Apache Parquet configuration](https://parquet.apache.org/docs/file-format/configurations/)
recommends 512 MB–1 GB row groups for analytics scans; HF `datasets` instead uses 100 rows per row group for image
data to bound the bytes read per random access
([config.py](https://github.com/huggingface/datasets/blob/main/src/datasets/config.py)), and the Hub viewer caps
row groups near 300 MB. pyarrow `write_table` defaults: snappy, dictionary on, statistics on
([docs](https://arrow.apache.org/docs/python/generated/pyarrow.parquet.write_table.html)). Chosen: zstd level 6,
dictionary and statistics off (binary payload, no predicate push-down needed), row groups closed at 128 MB of
payload.

**Masks.** COCO RLE is column-major run lengths with LEB128 string compression
([pycocotools mask.py](https://github.com/cocodataset/cocoapi/blob/master/PythonAPI/pycocotools/mask.py));
`np.packbits` is one call each way.

## Sandbox comparisons (scratch `impl/M4/`)

**Image encoding inside Parquet** (12 RIGA photographs, 48 LIDC slices; `codec_cmp.py`):

| Encoding | RIGA MB / item | RIGA decode ms | LIDC MB / item | LIDC decode ms |
|---|---|---|---|---|
| PNG (Pillow), LIDC as 16-bit PNG of HU + 32768 | 2.96 | 53.5 | 0.269 | 3.0 |
| NPY + Parquet zstd (default level) | 5.30 | 2.0 | 0.315 | 0.08 |
| NPY + zstd 3 | 4.42 | | 0.289 | |
| NPY + zstd 9 | 3.78 | 2.1 | 0.273 | 0.09 |

On the real RIGA validation file: zstd 3 / 6 / 9 = 198 / 177 / 172 MB, written in 3.4 / 5.4 / 9.6 s. Chosen:
NPY + zstd 6. It is lossless for any dtype (int16 HU needs no offset), needs no image codec at load (decoded pixels
do not depend on a JPEG/PNG library version), decodes 25× faster than PNG, and the full RIGA store came out at
2.8 MB per photograph, close to PNG.

**Mask encoding** (8 RIGA items × 6 raters, full frame):

| Encoding | Parquet KB / item | decode ms |
|---|---|---|
| bit-packed (`np.packbits`) | 3.6 | 2.6–5.1 |
| row-major RLE (lab) | 1.9–2.1 | 6.2–7.7 |
| NPY bool | 5.6 | — (8× the memory: 45 MB per RIGA row) |

Chosen: bit-packed. Within 2 KB of RLE after zstd, no custom codec, 8× smaller in memory than bool.

**Full image vs ROI** (`roi_need.py`). The lab's model crops (3× the object extent, centre jittered over the
middle 60 %) need the window centre ± 2.4 × extent: a median 38 % (min 21 %) of a RIGA photograph, and a median
72 px (max 365 px) window on a 512 px CT slice. An ROI would save about 60 % of RIGA bytes and most LIDC bytes, but
it bakes M5's crop policy into the data and adds an offset to every consumer. The full store is 3.3 GB, so items
keep the full frame.

**Parity runs:**
- RIGA, production `extract` vs the lab extractor re-run (`riga_parity.py`): 4,500 / 4,500 masks bit-exact,
  all certified; stages A 4,456, B 41, C 3; references prime 699, F1 50, F7 1. 3:52 with 3 processes
  (lab: 4:46).
- LIDC, production `slices` vs the lab fixture `lidc_full.npz` (`lidc_parity.py`): 10,921 / 10,921 slices
  bit-exact with identical crop boxes, reader ids and session counts; XML choice and patient id equal
  `series_patient.csv` for all 1,018 series. 11.5 s.
- LIDC vs pylidc DB contours (oracle in `test_data_lidc.py`, 25 random scans): annotation sets equal on every scan,
  848 / 848 annotation-slice masks exact. (Lab, all 1,018 scans per reader-slice: 99.97 %.)

## Decisions

- **Item frame.** Image and masks in the original frame: RIGA full photograph (uint8 RGB, the reference the raters
  drew on, so BinRushed1 is the F1 reconstruction), LIDC full 512 × 512 slice (int16 HU). Crop and scale are M5's.
- **LIDC pixels.** Only `data/lidc_slices` (the Phase 0 fetch: 214 series, annotated slices only) has CT pixels;
  the XML-only rest of LIDC has none locally. All 10,921 items are stored; `image is None` for the 6,837 without
  pixels. Fetching the other 804 series is a coordinator decision (see open issues).
- **Disagreement for stratification.** The 5-px-filter per-instance signature (`objdel.checker.signature(m, 5)`,
  M3), as the lab split used. The plan's size rule and the merged-reader flag are reported, not used to split.
- **Memory.** RIGA extraction runs in bounded batches (2 × workers photographs; a result is about 56 MB), row
  groups close at 128 MB of payload, and the default is 2 workers so the build stays at 3 processes.
- **Reader ids.** LIDC raters are reading-session numbers within the series ("1"–"4"; TCIA: reads are not
  comparable across cases). RIGA raters are "1"–"6".
- **Removed from the lab path:** hull-IoU audit vs `RIGA_preprocessed`, inventory, thumbnails, overlay, compare,
  pHash/dHash (NCC decides), HF `datasets` on write (plain pyarrow; the files stay readable by `datasets`),
  polygon codec, Phase 0 z-index heuristic.

## Tests

| File | Unit (synthetic) | Data-backed (`@pytest.mark.data`) |
|---|---|---|
| `test_data_riga.py` | discover (F2, F1 naming), max-RG reference (F1), hysteresis = scikit-image (Hypothesis), stage A, stage B bridge (F4), stage C split (F5), blank copy (F6), prime-is-a-copy (F7), F1 end to end, certification, item format | store counts (750 / 4,500 / 95 / 655 / 88 unique); stage-A masks = Phase 0 `riga_true` |
| `test_data_lidc.py` | P1, non-CT messages, P2 inclusion and exclusion order, P3 z mapping, clustering, merged-reader flag, slice records, duplicate-XML choice, pylidc DB reader, CT slice lookup, item format | all 10,921 store items = lab fixture; pylidc parity sample (25 scans) |
| `test_data_splits.py` | the lab's 11 tests on synthetic populations (patient completeness, series follow patients, literature split + groups, no cross-split near-duplicates, same / different seed, strata within 1, disagreement share within 2 %, val carve-out, two negative controls), near-duplicate detection | K11 on the store; 313 / 251 on the lab split (296 / 242 since F4) |
| `test_data_store.py` | exact round trip (RGB, int16, no image), mask codec (Hypothesis), K10 same bytes, manifest, flipped byte detected by `verify` and `load`, missing file, 5 contract violations, sorted ids, one-row row groups | `verify`; load → write gives identical file hash |

## Open issues

- **LIDC CT pixels exist for 214 of 1,018 series** (the Phase 0 fetch): 4,084 / 10,921 items have an image.
  *(Resolved: I9, `experiments/phase1/fetch_lidc_pixels.py`; with B3-4 every stored LIDC item has pixels.)*
  Fetching the other 804 series from TCIA (about 100 GB of DICOM streamed one series at a time, annotated slices
  kept: about 6 GB) needs a coordinator decision; `experiments/phase0/lidc_fetch.py` already does it per series, and
  a rebuild picks the new slices up with no code change. Until then SFT/GRPO on LIDC can only use the items with
  `image is not None`; reward and caches (M6) need masks only.
- **M3 robustness (reported, not fixed: not my file; since fixed by `topology.as_mask`, F1 / I1).** `objdel.checker.topology` reads mask bytes: a bool array
  whose True bytes are 255 (what `np.asarray(pil_mode_1_image, dtype=bool)` returns) loses holes
  (`instance_topology` gave holes 0 instead of 1 on 10 LIDC slices). `lidc.draw` now returns canonical bools and a
  regression test covers it; masks from `store.load` are canonical anyway (packbits). Suggest M3 normalise with
  `np.asarray(mask).view(np.uint8) != 0` or equivalent.
- RIGA split is 591 / 64 / 95, one image off the plan's 590 / 65 (see Result). *(Final: 592 / 63 / 95, F4.)*
- The near-duplicate threshold (0.45) and descriptor are from the lab audit (checked by eye there); this run
  reproduces the gap (largest sub-threshold NCC 0.326, smallest duplicate ≥ 0.45, 37 pairs as in the lab).
- `pylidc` (for its bundled database) and `matplotlib` (parity-test oracle, skipped if absent) are used but are not
  declared in `pyproject.toml`; `skimage` is used only by an optional differential test. *(Since declared:
  `pylidc` in the `data` extra, `matplotlib` in `dev`.)*
- Masks follow the LIDC "outline outside the nodule" convention filled inclusively (plan rule P2): about 1 px
  generous; matters for nodules under 50 px.

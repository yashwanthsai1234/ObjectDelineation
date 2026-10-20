# Break-it B3: data, targets, prepared set

Scope: `src/objdel/data/` (riga, lidc, splits, store), `src/objdel/targets/` (frame, build), `src/objdel/train/dataset.py`
(prepare, records, sft_rows, grpo_rows, cache_loader), `experiments/phase1/{build_data,build_targets,prepare_rl,fetch_lidc_pixels}.py`.
Data: `data/store` (manifest `da5b3dc0…`, CT pixels for 10,917 / 10,921 LIDC items) and `data/prepared` (rebuilt
from that store; its manifest pins the same store SHA-256). Reproducers: `tests/test_breakit_data.py` (strict xfails
`B3-<n>`, plus the passing robustness checks). Scratch scripts: `scratchpad/breakit/B3/scripts/`.

## Findings

| ID | Severity | Module | One line |
|---|---|---|---|
| B3-1 | high | data.riga (F1 reference) → targets.frame → prepared images | BinRushed1 model inputs keep the raters' pen strokes where all six copies are inked: rater outlines show in the SFT/GRPO image |
| B3-2 | medium | targets.frame (LIDC crop) | A frame can contain a second annotated nodule that is not in its ground truth: 114 prepared items, 48 with the other nodule wholly inside |
| B3-3 | medium | data.lidc (readers as raters) | Two reading sessions often hold byte-identical contours (copied outlines). They are stored as two independent raters: 420 items, 149 with ≥ 100 px |
| B3-4 | low | experiments/phase1/fetch_lidc_pixels (phase0 `lidc_fetch`) | CT slices are matched by the XML's SOP UIDs. LIDC-IDRI-0365's XML cites UIDs absent from TCIA, so 4 train items have no pixels although their z positions exist |
| B3-5 | low | train.dataset (sft_rows / grpo_rows) | Prepared PNGs are never checked against `image_sha256`. A swapped or edited image reaches the rows silently |
| B3-6 | low | data.store.write | `item_id` is unique only per file. The same id in two splits is written (a leak the store cannot see; `cache_loader` keys collide) |
| B3-7 | low | data.store / targets.frame | The image dtype contract (uint8 RGB RIGA, int16 HU LIDC) is not enforced, and `frame` picks the lung window by dtype |
| B3-8 | medium | data.splits (near-duplicate search, K11) | Same-eye photographs across train/val: BinRushed4/image002 (train) is the eye of the val group BinRushed3/image004 + 026. A 2.5° / 0.96× / 83 px reframing defeats the ±12 px NCC search |

### B3-1 (high): BinRushed1 inputs show the rater outlines
- **Reproduction:** `test_b3_1_binrushed1_inputs_carry_no_rater_ink`. For every prepared RIGA train item, it counts
  near-black pixels (R < 80, G < 80) within 2 px of any rater's disc/cup edge in the model image
  (`data/prepared/images/…png`).
- **Measured** (`scripts/ink.py`, train + val):

  | subset | items | median dark px on rater edges | items ≥ 30 | control band (8–14 px away), items ≥ 30 |
  |---|---|---|---|---|
  | BinRushed1 (reference = max-R+G reconstruction, F1) | 50 | 12 | **18** (16 train, 2 val; max 285) | 0 |
  | BinRushed2 / 3 / 4 (prime) | 47 / 47 / 51 | 0 | 1 / 0 / 0 | 2 / 1 / 1 |
  | MESSIDOR (prime) | 460 | 0 | 1 | 5 |

  `scratchpad/breakit/B3/viz_riga.png` shows it on image042 and image030: black arcs on the disc boundary.
- **Expected vs actual:** expected the model input to be the photograph without annotations, as for every other subset.
  Actually, `riga.max_rg` picks, per pixel, the copy with the largest R+G. Where all six raters' strokes overlap,
  every copy is inked, so the "reference" stays pen-coloured. `riga.item` stores it as `image`, and `frame` crops it
  into the SFT/GRPO input. The source code says so itself (`drawing`: "Where every copy is inked the F1 reference is
  pen-coloured too").
- **Why it matters:** the input leaks the answer on about 3 % of RIGA train images and on 2 val images. The model can
  learn to trace ink. Test images (Magrabia, clean primes) have none, so train and test differ, and val scores on
  BinRushed1 are optimistic. The ten over-long targets of I5 come from the same images (BinRushed1 image018, 030,
  042, 045, 048).
- **Suggested fix:** keep the reconstruction for mask extraction, but build the *model* image without ink. Inpaint the
  pixels that are pen-like in the reference (`riga.pen_like(ref)`, dilated 1–2 px), e.g. with
  `cv2.inpaint(..., INPAINT_TELEA)`
  ([OpenCV inpainting](https://docs.opencv.org/4.x/df/d3d/tutorial_py_inpainting.html)). Or exclude BinRushed1
  from the image-conditioned sets and keep it for reward statistics only. Add the test above as the gate.

### B3-2 (medium): a second annotated nodule inside the crop
- **Reproduction:** `test_b3_2_frame_contains_no_other_annotated_nodule` (removed after F4; regression `test_data_lidc.py::test_store_lidc_crowded_break_it_item_is_left_out`). Item
  `lidc/…137763212752154081977261297097/n000/z092` (LIDC-IDRI-0942, train) has a frame box of
  (67, 293, 167, 393). Its neighbour on the same slice, `…/n001/z092`, is a separate 306 px nodule that lies wholly
  inside that box and outside the item's masks (`scratchpad/breakit/B3/viz_pair.png`: the green nodule is visible
  and not in the ground truth).
- **Scale** (`scripts/lidc_align.py`): 114 prepared items (94 train, 20 val) have another item's nodule pixels
  inside their frame; for 48 the whole other nodule is inside; 11 of them are disagreement items (×4 in GRPO). The
  count is a lower bound: nodules outlined by a single reader, or under 3 mm, are not items, so they are not counted.
- **Expected vs actual:** the prompt says "Outline the lung nodule in this CT crop", singular. The answer should be
  unambiguous given the image. Actually, outlining the other visible nodule is scored as spill (IoU_tc, F5, area
  envelope), and SFT teaches the model to ignore a nodule it can see.
- **Root cause:** `frame` crops around one cluster's union only. `lidc.slices` keeps clusters separate, which is
  correct, but nothing checks the crop for other annotations.
- **Suggested fix (pick one):** (a) shrink or shift the window to exclude other clusters' pixels when possible, else
  drop the item from SFT/RL; (b) mask other annotated nodules out of the image (fill with local lung HU); (c) make
  the item "all nodules in the crop" (all clusters' readers, per-object instances). Keep the test as the gate.

### B3-3 (medium): copied contours count as independent raters
- **Reproduction:** `test_b3_3_reading_sessions_are_independent` (removed after F4; regression `test_data_lidc.py::test_store_lidc_copied_reader_of_0811_is_one_rater`). In `tcia-lidc-xml/188/294.xml` (LIDC-IDRI-0811),
  sessions 1 and 3 share 42 contours with identical `edgeMap` point lists. Store item
  `lidc/…070084384692/n000/z243` holds raters 1 and 3 as two identical 1,211-px layers, with no flag in `meta`.
- **Scale** (`scripts/dupread.py`, `scripts/xmlpts.py`, store scan):
  - 144 reader pairs in 126 series share at least one identical contour of more than 4 points. In 11 pairs one
    session's every contour equals the other's.
  - In the store, 430 rater pairs in 420 items (train 275, val 38, test 107) have identical masks. In 397 items
    the mask is ≥ 20 px, in 149 it is ≥ 100 px. 158 / 166 of the cluster-level pairs are byte-identical in the XML,
    not just identical after rasterizing.
- **Expected vs actual:** the LIDC protocol lets a reader, in the unblinded phase, "leave unchanged, delete, switch in
  lesion category" *their own* marks, or add marks (Armato et al. 2007,
  [PMC2151472](https://pmc.ncbi.nlm.nih.gov/articles/PMC2151472)). Two independent pen outlines cannot match point
  for point over 1,211 px. Actually, the copies enter the reward as independent experts:
  - a zero-distance pair lowers the pooled inter-rater τ (P95 of rater-pair distances);
  - the T5 calibration (`reward/cache.py`) assumes i.i.d. raters;
  - the topology histogram f and the vote weigh the copied outline twice.
- **Suggested fix:** in `lidc.slices`, detect layers whose source contours are identical (point lists, or equal
  masks ≥ 20 px). Then either merge them into one rater, or keep both and add a `meta["copied_raters"]` flag that
  `ExpertCache` uses to deduplicate. Report the count in `build_report.json`.

### B3-4 (low): 4 LIDC items without pixels because SOP UIDs do not match
- **Reproduction:** `test_b3_4_every_lidc_item_has_pixels` (removed after F4; regression `test_data_lidc.py::test_store_lidc_every_item_has_pixels`). `data/lidc_slices/<0365 series>.npz` exists but is
  empty. The series' contour z values (−858 … −849 mm) lie on its pylidc slice grid. `build_report` gives train 7,619
  items with an image out of 7,623.
- **Evidence (TCIA NBIA API, `getSOPInstanceUIDs` / `getSeries`):** none of the 9 SOP UIDs in 0365's XML exist in
  either of the patient's two CT series. The same happens for LIDC-IDRI-0017: readers 2 and 4 cite 5 SOP UIDs that
  are not in its only CT series; readers 1 and 3 cite the real ones at the same z. Contours there still align by z,
  as pylidc does.
- **Root cause:** `experiments/phase0/lidc_fetch.fetch` keeps a DICOM only if `SOPInstanceUID` is in the XML's
  SOP→z map. The store, however, maps contours to slices by z (`lidc.slice_index`, P3).
- **Suggested fix:** in the fetcher, also key by `ImagePositionPatient[2]`. Keep every slice whose z matches a
  contour z within 0.01 mm, the same rule as `ct_slices`. Then re-fetch the 5 series whose json reports
  `slices < annotated_sops` (0013, 0017, 0365, 0566, 0659).

### B3-5 (low): prepared images are not verified
- **Reproduction:** `test_b3_5_tampered_prepared_image_is_detected`. After swapping two PNGs of a tiny prepared
  set, `grpo_rows` returns the wrong image for the key, with no error. The decoded pixels' SHA-256 differs from the
  record's `image_sha256`.
- **Expected vs actual:** the module docstring says the manifest "pins every image and cache through their hashes".
  Caches are checked (`cache_loader`; a swap raises, which `test_cache_loader_detects_swapped_caches` confirms), but
  nothing reads `image_sha256`.
- **Suggested fix:** add a `dataset.verify(root)` (decode every PNG, compare `_pixels_sha`), called by `prepare_rl` and
  before training; or check in `sft_rows` / `grpo_rows`.

### B3-6 (low): duplicate item ids across splits
- **Reproduction:** `test_b3_6_store_refuses_one_item_id_in_two_splits` (removed after F4; regression `test_data_store.py::test_write_rejects_one_item_id_in_two_splits`). `store.write` accepts one `item_id` in
  `train` and `test`.
- **Root cause:** `write` keeps `last[rel]` per file only. `cache_loader` builds `{item_id: caches}` from every split's
  records, so a second split silently overrides the first.
- **Fix:** keep a set of ids across all files in `write` and raise. The real store has none: the B3 scan found
  10,921 + 750 unique ids.

### B3-7 (low): image dtype contract not enforced
- **Reproduction:** `test_b3_7_store_enforces_the_image_dtype_contract` (removed after F4; regressions `test_data_store.py::test_write_enforces_image_dtype_and_layout`, `test_load_enforces_the_item_contract`). An LIDC item with a uint16 (HU + 1024)
  slice is written without error.
- **Why it matters:** `frame` decides "Hounsfield units" by `image.dtype != np.uint8`. A uint8 LIDC slice would skip
  the lung window. A uint16 offset slice is windowed on the wrong scale: −800 HU stored as 224 maps to 1.0 instead
  of 0.14. The real store is correct: all 10,917 LIDC images are int16 512 × 512, and all RIGA images are uint8 RGB.
- **Fix:** in `store._check`, require uint8 (H, W, 3) for RIGA and int16 (H, W) for LIDC. Dispatch the window on
  `item.dataset`.

### B3-8 (medium): a same-eye pair crosses train/val (K11)
- **Reproduction:** `test_b3_8_same_eye_photographs_share_a_split` (removed after F4; regressions in `test_data_splits.py`: `test_reframed_same_eye_found_by_keypoints`, keypoint K11).
- **Found by** an independent duplicate search (`scripts/orb_k11.py`): ORB on the CLAHE green channel at 800 px
  width, Lowe ratio 0.75, RANSAC homography inliers, mirrored images included, over all 100,049 cross-split pairs.
  - Calibration: the store's 37 duplicate pairs give ≥ 207 inliers (median 991).
  - Exactly two cross-split pairs reach that range, both 269 inliers: BinRushed3/image004 and image026 (one val
    group, a same-photo pair) against BinRushed4/image002 (train, its own group).
  - The next best cross-split pair has 35 inliers.
- **Confirmation:** on the full photographs, ORB + RANSAC finds a similarity transform: scale 0.962, rotation
  2.47°, shift (83, −9) px. The pixels differ (median |Δ| 22 after alignment), so this is a second photograph of
  the same eye, not a copy. The store's own descriptor (`riga.descriptor`, shift-tolerant NCC) gives 0.271 raw,
  below `NCC_DUP` = 0.45. After the alignment it gives **0.825**. An unrelated photograph under the same
  transform gives 0.03.
- **Expected vs actual:** M4 groups same-eye pairs (6 of the 37) so that they share a split; K11 requires 0
  cross-split near-duplicates. Actually, the descriptor only tolerates translation of ±12 thumbnail px. It has no
  rotation or scale tolerance, so a reframed same-eye photograph lands in another split. The leak is val ↔ train,
  so 2 val images (3 % of RIGA val) are scored against an eye the model trained on.
- **Suggested fix:** add a keypoint check (ORB/SIFT + RANSAC inlier count, threshold set between unrelated pairs
  ≤ 35 and duplicates ≥ 207) to `splits.near_duplicates`, or search the NCC over a small grid of rotations and
  scales. Merge the groups and re-split. Fundus pairs of one eye under rotation, scale and shift are the standard
  registration setting, evaluated with keypoint correspondences on vessels and crossings (FIRE,
  Hernandez-Matas et al. 2017; https://arxiv.org/pdf/2006.01030 uses it for keypoint matching).

## Deepened known issues
- **B2-4 (sub-filter raters) on the real store:**
  - 21 items have a rater whose size-filtered signature is `()` (train 14, val 1, test 6).
  - In 18 of them, that rater is the only source of disagreement (train 12, val 1, test 5). That is 12 / 311 =
    3.9 % of the prepared LIDC train disagreement items, each repeated ×4 by `grpo_rows`.
  - On 3 items every rater is below the filter (2–4 px nodules): `…331179/n001/z096` and `…625409/n000/z189`
    (train), `…085942/n005/z371` (test).
- **I9 / I23:** after the fetch, 10,917 / 10,921 items have pixels. The 4 missing ones are B3-4. The prepared set
  reproduces the counts: LIDC train 7,619 items, 311 disagreement items, so `grpo_rows(…, 4)` has 8,552 rows; val
  1,085 / 46.
- **I5:** with each model's tokenizer, `sft_rows` drops exactly the over-budget answers: SmolVLM and Qwen 10, LFM 8,
  all on BinRushed1 image018, 030, 042, 045, 048 and BinRushed3 image014. Every dropped item stays in `grpo_rows`.

## Attacks that held
Store, real data (`scripts/scan_store.py`, all 11,671 rows streamed in 2-row batches; peak RSS 2.3 GB):
1. Every item: masks (k, H, W) with k = len(raters), the image shares (H, W), no empty rater layer, rater ids unique and
   sorted, `dataset` / `split` fields equal the file, item ids unique across all files, row counts equal the manifest.
2. RIGA: 750 × 6 raters. Cup ⊆ disc for every rater. Disc and cup are each one component without holes, and the rim is
   non-empty. `stages` keys equal the raters. `group_id` = smallest id of the duplicate group. Images are uint8 RGB.
   Split counts are 591 / 64 / 95, equal to `build_report`. No union touches the image border.
3. LIDC: `group_id` = `meta.patient` = pylidc patient. Item-id series / cluster / slice equal `meta`. `z_mm` equals the
   pylidc slice z. Raters ⊆ `readers_in_cluster`. `pixel_spacing_mm` equals pylidc. Every image is int16 512 × 512.
4. LIDC z mapping: all 41,406 contour ROIs lie exactly on a pylidc slice z (max |Δz| 0.0 mm). No reader has two
   inclusion z values mapped to one slice. pylidc's z tables have no duplicate z. That is the known LIDC duplicate-z
   pitfall: series 0085, 0146, 0418, 0572 and 0979 have DICOM slices sharing one z
   ([LIDC-IDRI-Lite notes](https://huggingface.co/datasets/YongchengYAO/LIDC-IDRI-Lite)).
5. CT pixels ↔ z: I re-downloaded LIDC-IDRI-0120 from TCIA and matched by `ImagePositionPatient` z. All 27 stored
   slices are bit-identical to `pixel_array · slope + intercept`. That includes 12 out-of-range pixels per slice
   (±32,7xx HU at fixed positions), which are in the source DICOM itself, not an int16 wrap.
6. Statistical alignment over all 10,917 items: the mean HU inside the union is above a 2–6 px ring around it on
   95.2 % of items, against 34.9 % for the same masks shifted 15 px. The lowest-contrast series are ground-glass or
   juxta-pleural nodules, checked by eye (`viz_lowcontrast.png`).
7. LIDC XML: duplicate XML files of 6 series carry identical annotations. The one resubmitted correction (0101)
   differs and is the file chosen. One XML (0566) has 8 sessions, 4 of them empty; one (0700) has 3.
8. HU ranges: minimum values −3024 / −2048 / −1024 (scanner padding), 0 slices with max < 500, 0 with min > −900.
   The extremes are the 0120 source pixels above.

K11, independent methods:
9. RIGA, ORB + RANSAC over all 100,049 cross-split pairs. Mirrored images were tried as well. Calibrated on the
   37 store duplicate pairs (≥ 207 inliers, median 991). One group pair crosses splits, which is **B3-8**. Every other
   cross-split pair has ≤ 35 inliers, including every train/val ↔ Magrabia-test pair, so the test split holds.
10. LIDC: no StudyInstanceUID (XML header) is shared by two patients, and no series has two. Near-identical CT across
    patients was searched on 32 × 32 lung-window thumbnails over 10,367 unique slices: max cross-patient Pearson
    r = 0.961 (two different train patients), 0 pairs above 0.98, 0 cross-split.

Determinism and tampering:
11. LIDC val rebuilt from raw XML + CT archives into scratch, in one process: SHA-256 `36973066…` equals the store's
    (278,172,830 bytes).
12. `store.load` gives the same items with reader batch 1, 5 and 10⁶ (`test_store_load_is_independent_of_the_reader_batch_size`).
13. Truncated file, appended byte, files swapped between splits, flipped last byte: `verify` and `load` raise
    (`test_store_tampering_is_detected`).
14. Cache files swapped between two items: `cache_loader` raises (`test_cache_loader_detects_swapped_caches`).

Frames, targets, prepared set (`scripts/check_prepared.py`, every train/val item):
15. `frame(item)` recomputed in a fresh process equals the record (W, H, box, scale) for all 8,704 LIDC and 655 RIGA
    items. The PNG equals the fresh frame image bit for bit and its `_pixels_sha` equals `image_sha256`. RIGA PNG =
    `item.image[y0:y1, x0:x1]` exactly (image crop = mask crop). LIDC PNGs are 512 × 512 × 3.
16. No crop loses a rater pixel: frame mask sums equal item mask sums for every structure.
17. Caches: for every key, `cache_loader` loads, `ExpertCache.raters` equals the frame masks in rater order, and k
    equals len(raters). The recomputed disagreement flag equals the record.
18. Prepared splits equal store splits, keys are unique, no record exists without pixels. `sft_rows` has one image
    placeholder and one image per row. `grpo_rows` oversampling = rl + 3 × disagreement
    (`test_prepared_rows_match_records_and_oversampling`).
19. **R1 on the full prepared train + val sets** via `score_completion` against each item's own caches:
    - LIDC: 27,978 targets; format, legal, own-rater signature and topo_any 100 %; IoU_tc min 1.0; F5 min 1.0.
    - RIGA: 3,930 targets; format, legal, own-rater signature and topo_any 100 %; IoU_tc min 0.99005;
      F5 ≥ 0.95 100 % (min 0.988).
20. Target determinism across processes: every 25th item's targets (347 LIDC, 25 RIGA items, all raters) were rebuilt
    with `build_target` in the checker process. Identical to the prepared texts made by pool workers.
21. Tiny objects: 19 items have a rater mask ≤ 4 px; all their targets certify (included in 19). Frame side
    S = 50 / 100 / 125 (LIDC) and 200 / 250 / 500 (RIGA); scale 1 everywhere.
22. Scale > 1 with a window larger than the image (padded both sides): the frame mask is the block majority of the
    padded item mask, and the padding is black (`test_frame_scale_gt_1_maps_grid_to_item_pixels_at_the_border`).
23. LIDC 512 upscale: a bright square drawn on the mask spans 512 / S × the mask extent in the model image, within
    1 px, at S = 50, 100 and 200 (`test_lidc_upscaled_image_lines_up_with_the_mask`).

## Sources
- TCIA LIDC-IDRI collection: https://www.cancerimagingarchive.net/collection/lidc-idri/
- TCIA NBIA REST API (getSeries, getSOPInstanceUIDs, getImage): https://services.cancerimagingarchive.net/nbia-api/services/v1/
- Armato et al. 2007, LIDC reading protocol and QA: https://pmc.ncbi.nlm.nih.gov/articles/PMC2151472
- LIDC-IDRI-Lite (duplicate-z series, pylidc index pitfalls): https://huggingface.co/datasets/YongchengYAO/LIDC-IDRI-Lite
- OpenCV inpainting: https://docs.opencv.org/4.x/df/d3d/tutorial_py_inpainting.html
- FIRE fundus registration pairs / keypoint matching on fundus images: https://arxiv.org/pdf/2006.01030

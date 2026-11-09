# Data reference

Downloaded data is untrusted. It lives in `data/` (git-ignored). Python that reads it runs with `-I`.
No image or annotation of LIDC-IDRI or RIGA is in the repository.

## 1. Datasets

| | LIDC-IDRI | RIGA |
|---|---|---|
| Domain | Lung CT, nodules | Fundus photographs, optic disc and cup |
| Source | TCIA, DOI 10.7937/K9/TCIA.2015.LO9QL9SX. Reader XML from TCIA, CT from DICOM | Almazroa et al. 2018, University of Michigan Deep Blue (`riga.py` docstring) |
| Raters | Up to four reading sessions per nodule | Six ophthalmologists per photograph, each drew on their own copy |
| Structures | `nodule` | `disc`, `cup` |
| Model image | int16 Hounsfield units in the store. The frame applies the lung window (-1000, 400) | uint8 RGB |
| Series or photographs | 1,018 series, 1,010 patients, 8 patients with two scans (`lidc.py` docstring) | 750 photographs in 7 subsets (BinRushed 1 to 4, MESSIDOR, Magrabia female and male) |
| Store rows | train 7,341, val 1,084, test 2,200 (10,625 items) | train 592, val 63, test 95 (750 photographs) |
| Used for RL | Yes. `train_rl` preloads `lidc` only | No |
| Licence | CC BY 3.0 (TCIA). May go on any host | RIGA's licence forbids redistribution. Goes only to Secure Cloud hosts, with written permission from the depositor and the MESSIDOR partners (`experiments/phase1/RESULTS.md` I44) |

Store row counts are read from `data/store/manifest.json` (rows per Parquet file). RIGA train and test are from the same manifest. The val count is 750 - 592 - 95 = 63.
The LIDC items are slices that at least two readers outlined. Prepared-set counts are in section 4.

### Rater rules that change the data

| Rule | Effect | Source |
|---|---|---|
| P1 | Single-point contours inside a nodule of 3 mm or more are dropped | `data/lidc.py` |
| P2 | One outline rule: a contour covers its filled polygon including the outline pixels. Exclusions subtract | same |
| P3 | Slice index is the nearest scan slice position from the pylidc database | same |
| P4 | Patient ids from the pylidc database | same |
| P5 | Readers whose masks on a slice are identical count once. The lowest session is kept. These are outlines one reader adopted from another | same |
| P6 | Other outlined nodules on the slice are recorded. A window that shows one is refused | same |
| RIGA F1 to F7 | Staged extractor from pen drawings. BinRushed1 has no clean reference photograph: it keeps masks but has `image = None`, so it is skipped for training | `data/riga.py` |

Items without a model image are skipped by `prepare`, and the manifest counts them: RIGA train 45, RIGA val 5.

## 2. Store (`data/store`, `objdel.data.store`)

```text
data/store/
  manifest.json            SHA-256, rows and bytes per file, library versions that wrote it
  lidc/{train,val,test}.parquet
  riga/{train,val,test}.parquet
  build_report.json        counts, histograms, K11 split checks, timings
```

- One Parquet file per dataset and split. Rows are sorted by `item_id` and written with fixed settings (zstd level 6, no dictionary, row groups of about 128 MiB), so the same items give the same bytes.
- Image: the decoded pixel array as NPY bytes. Lossless for any dtype. No JPEG or PNG library in the loop.
- Masks: per structure, the bool stack `(k, H, W)` bit-packed along W with `np.packbits`, as NPY bytes.
- Item id: `lidc/<series uid>/n<cluster>/z<slice>` or `riga/<subset>/image<NNN>`.
- Group id: the patient for LIDC. `riga:<group>` for RIGA, where a group is a set of near-duplicate photographs.
- Checks on write and on load: item ids unique across the whole store, row dataset and split match its file, images are uint8 `(H, W, 3)` for RIGA and int16 `(H, W)` for LIDC or `None`.
- `load` verifies the file's SHA-256 first and reads one row at a time (a RIGA row can be 16 MB).

## 3. Splits (`objdel.data.splits`)

Splits are group-aware, seeded and deterministic. Within a stratum, groups sort by SHA-256 of seed and group id, and a largest-remainder schedule deals them. The result depends only on the seed and the population.

| Dataset | Unit | Rule |
|---|---|---|
| LIDC | Patient. Both scans of a two-scan patient stay together | 70 / 10 / 20 (`LIDC_FRACTIONS`). Patients without reader disagreement are dealt per stratum of multi-reader slice count. Patients with disagreement are dealt by weighted greedy on their disagreement-slice count |
| RIGA | Near-duplicate group | Literature split: test = Magrabia (95). Train = BinRushed plus MESSIDOR with a 10 % validation carve-out stratified by subset. A group that touches the test set goes to test |

Near-duplicate photographs (RIGA) form one group when either test fires: normalised cross-correlation of vessel-texture descriptors of at least 0.45 (`NCC_DUP`), or at least 22 RANSAC inliers of ORB matches after Lowe's ratio test 0.75 (`ORB_DUP`, `RATIO`).
`leaks(split_of, group_of)` and `crossing(pairs, split_of)` return violations. The store build and `preflight` check both.

## 4. Prepared set (`data/prepared`, `objdel.train.dataset.prepare`)

Built by `experiments.phase1.prepare_rl`. One item with pixels gives:

```text
data/prepared/
  images/<dataset>/<item>.png       the model frame, uint8 RGB (lossless PNG)
  caches/<dataset>/<item>/<object>.npz   ExpertCache.save (only items with >= 2 raters)
  <dataset>/<split>.jsonl           one record per item
  manifest.json                     SHA-256 of every jsonl, target-set hash per split, counts, reward fingerprint, versions
```

A record holds: `item_id`, `group_id`, `W`, `H`, `box`, `scale`, `raters`, `targets` (one canonical text per rater, or `None`), `rl`, `disagreement`, the SHA-256 of the image pixels and each cache digest.
Every row builder decodes each image it uses and compares its pixels with the recorded hash. A swapped PNG raises `ValueError`.

Counts in this checkout (`data/prepared/manifest.json`):

| Set | Items | Targets | RL items | Disagreement items | Skipped (no image) |
|---|---|---|---|---|---|
| lidc/train | 7,341 | 23,469 | 7,341 | 289 | 0 |
| lidc/val | 1,084 | 3,444 | 1,084 | 44 | 0 |
| riga/train | 547 | 3,282 | 547 | 0 | 45 |
| riga/val | 58 | 348 | 58 | 0 | 5 |

RL rows for LIDC: 7,341 items, with 289 disagreement items repeated 4 times: 7,341 + 3 x 289 = 8,208 rows. This is the 8,208 in `rl_full.yaml`.
SFT rows are one per voting rater with a target (23,469 for LIDC train). A target longer than the model's answer budget is dropped from SFT.

## 5. Frames and views

### The frame (`objdel.targets.frame`)

A frame is the square region the model sees. The policy frame centres on the union box of all raters:

```text
extent = larger side of the union box
need   = extent + 2 * max(8 px, ceil(25 % of extent))      (margin_px = 8, margin_pct = 25)
f      = ceil(need / 500)
S      = smallest divisor of 1000 that is >= max(ceil(need / f), 50)
```

S is one of 50, 100, 125, 200, 250, 500. One pixel is exactly 1000 / S grid units, so a mask to polygon to mask round trip has no rounding error.
When f > 1 masks reduce by f x f block majority. No LIDC or RIGA item needs f > 1 (the largest union is 309 pixels).
The LIDC image is upscaled from S to 512 pixels with Pillow's bicubic filter (bit-identical on every CPU). Masks, W and H stay at frame resolution.

### View kinds

| View | Definition |
|---|---|
| `centred` | The policy frame. Object offset from the frame centre has a median of 0.018 of the side in the RL training rows |
| `shifted` | The policy side S, one uniform offset per item. The union stays inside with 2 pixels of clearance. No other outlined nodule appears |
| `rescaled` | The side moves to a neighbouring divisor (LIDC: 50, 100, 125) plus a shift. Items with no admissible neighbour keep S and are counted as `kept_side` |
| `rotated` | LIDC only. The shifted view's windows with a quarter turn (`r90`) applied. Report only |

D4 operations: `id`, `r90`, `r180`, `r270`, `fh`, `fv`, `tr`, `at`. LIDC uses all eight. RIGA uses `id` and `fh` because fundus orientation is physiological.
Training window draws also jitter the CT window (level -300 +- 60 HU, width 1,400 +- 80 HU) and add Gaussian noise of up to 10 HU.
Targets and expert tables are rebuilt from the transformed masks. Vertices are never transformed.

### View directories

A view directory holds the rows of one training or exam set. The word "view" also names the window kind above.

| Directory | Built by | Content |
|---|---|---|
| `data/views/lidc-bundle`, `-v2`, `-v3` | `experiments.phase2.build_views` | SFT row specs per epoch (`rows.jsonl`), contexts, `manifest.json` with `view_hash` |
| `data/views/lidc-train-shifted`, `-rescaled` | `build_eval_views --split train` | RL train windows with their own expert tables (`train_rl --train-views`) |
| `data/eval_views/lidc-val-shifted`, `-rescaled`, `-rotated` | `build_eval_views` | Exam views of the validation split, fixed seed |
| `data/eval_views/lidc-test-shifted` | `build_eval_views --split test` | The test split. Never scored before the final evaluation |
| `data/eval_views/probe-precision`, `probe-topology` | `experiments.phase2.probes build` | Synthetic probes on real lung patches |
| `data/checker_eval` | `experiments/cuda_checker` | Packed polygons and masks for the checker eval |

Two `build_views` builds are byte-identical. Eval views use fixed seeds. Row seeds come from `row_rng(seed, epoch, item_id, rater_index, rep)`, so they do not depend on item order.

## 6. Rebuild order

All commands run from the repository root with `.venv/bin/python -I -m`.

| Step | Command | Output |
|---|---|---|
| 1 | `experiments.phase1.fetch_lidc_pixels XML_ROOT data/lidc_slices` | CT slices for annotated positions. Resumable. `--refetch-incomplete` |
| 2 | `experiments.phase1.build_data --out data/store` | The store, splits, K11 checks |
| 3 | `experiments.phase1.prepare_rl --store data/store --out data/prepared` | Frames, targets, expert tables |
| 4 | `experiments.phase2.build_views --out data/views/lidc-bundle-v3 --epochs 4` | SFT training view |
| 5 | `experiments.phase2.build_eval_views --out data/eval_views` | Exam views of the validation split |
| 6 | `experiments.phase2.probes build` | Precision and topology probes |
| 7 | `experiments.phase2.vast.bundle lidc --view lidc-bundle-v3` | `build/bundles/bundle-lidc.tar` and its SHA-256 |

Raw inputs: `data/raw/riga` and `data/raw/lidc_xml`. `build_data` also reads the database bundled with pylidc.
`build_views` defaults to 3 epochs. The config `r1_qwen3vl2b.yaml` uses 4 view epochs because the e3 branch decays to 3.45 epochs.

## 7. Bundle (`experiments.phase2.vast.bundle`)

One deterministic tar per dataset. Members are sorted, mtime 0, owner 0, mode 0644, with `SHA256SUMS` first. Two builds of one tree are byte-identical.

| Member | Content |
|---|---|
| `prepared/` | `manifest.json` (rewritten to list only this dataset's files), dataset records, images, caches |
| `views/<view>/` | The `--view` names, or every `<dataset>-*` |
| `eval_views/<dataset>-*/` | LIDC also `probe-*`. LIDC requires `lidc-val-shifted`, `lidc-val-rescaled`, `probe-precision`, `probe-topology` |
| `checker_eval/` | The checker eval set (LIDC bundle only) |

The bundle refuses symlinks and anything under `data/raw`. RIGA streams into `age -r <recipient>` so only `bundle-riga.tar.age` touches the disk.
`vast.bundle lidc --view lidc-bundle-v3 --view lidc-train-shifted --view lidc-train-rescaled` ships the RL train windows.

## 8. Items that never reach training

| Reason | Where enforced |
|---|---|
| The model frame overlaps another outlined nodule on the slice | `data/lidc.crowded`, `build_data` |
| Only one rater | `rl_eligible`. `ExpertCache.build` raises |
| No model image (RIGA BinRushed1) | `prepare` skips and counts |
| A patient in both training and an exam view | `preflight` DATA stage |

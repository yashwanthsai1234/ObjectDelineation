# t_data: data pipeline validation

Saved by the coordinator from the agent's final message; the agent could not write this file.
Code is in `riga/`, `lidc_parity/`, `splits/` and `storage/`. Tests: splits 11/11, parity 8/8, storage 7/7.

## Verdict

| Check | Result |
|---|---|
| RIGA extraction | 750 × 6 = 4,500 / 4,500 masks (repo: 4,172 over 698 images). Nothing is missing from the download. |
| Mask validity | All certified. Disc and cup each b0 = 1, b1 = 0, cup inside disc. 2 drawings have the cup touching the disc. |
| Repo parity | 4,172 / 4,172 bit-exact |
| Train/test near-duplicates | 0. Max vessel-texture correlation 0.21 vs threshold 0.45 (gap 0.33–0.49). Min pHash distance 14. |
| Within-split near-duplicates | 37 pairs (31 same photo, 6 same eye). Magrabia test has 95 files but 88 unique photos. |
| LIDC vs pylidc (67 scans) | Clusters, reader sets and annotations match exactly (582 / 582). Masks exact on 99.57 % of reader-slices once the outline rule is aligned. |
| Splits | Deterministic (same seed → byte-identical files), group-aware. LIDC: 709 / 100 / 201 patients, disagreement slices split 69.7 / 10.1 / 20.2 %. RIGA: train 590, val 65, test 95. |
| Storage | Parquet (HF datasets). Bit-exact round trip. RLE masks are 3× smaller than polygons; polygons decode faster. |
| **K11** | **PASS** |

## Why the RIGA data was lost, and the fixes (staged extractor F1–F7)

| Cause | Fix |
|---|---|
| BinRushed1 has no prime image | Rebuild one from the 6 rater copies by per-pixel max R+G (validated on BinRushed2–4: disc IoU median 0.9985) |
| Magrabia Female image48 prime is named `image48.tif` and equals rater 1's copy | F2 accepts the name; F7 detects the copy |
| 22 masks failed extraction | Hysteresis threshold, gap-bridging before the blob filter, splitting a cup that touches the disc |

- The 44 masks recovered by the later stages have boundaries 1–2 px ragged (from JPEG). Their topology is clean, and a `stage` column records which masks they are.
- Duplicate photos were annotated independently by the same rater. Disc IoU median 0.946, cup median 0.81, giving a free measure of within-rater variability.

## LIDC quirks

- **Outline pixels:** pylidc never adds them; we fill including the outline. Pick one rule and document it.
- **Single-point ROIs:** 15 sit inside large nodules. Drop them.
- **z-index:** wrong in 6 / 67 scans. This matters only for 3D.
- **Odd XML files:** 188/011.xml has 8 sessions. 6 series have duplicate XMLs.
- **Patient IDs:** absent from the XML. `series_patient.csv` (from the pylidc DB) maps 1,018 series to 1,010 patients; 8 patients have two scans.
- **pylidc clustering tolerance:** uses slice thickness, not pixel spacing as the research note said.

## Suggested repo changes

1. Adopt the staged RIGA extractor and regenerate all 4,500 masks.
2. Group duplicate pairs in the RIGA split. Report Magrabia results on 95 and on 88 images.
3. LIDC parser:
   - drop single-point ROIs;
   - document the outline rule;
   - use the true slice index;
   - take patient IDs from pylidc.
4. Adopt `build_splits.py` and its tests as the K11 gate.
5. Storage:
   - keep the exact masks as the source of truth;
   - keep the original JPEG bytes;
   - use disc-centred crops for training;
   - keep the RIGA store private.

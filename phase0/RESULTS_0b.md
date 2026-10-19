# Phase 0b results — do models make topology errors that IoU misses?

Run 2026-10-18 on an Apple M3, 16 GB, inference only. Plan: `phase0/PLAN_0b.md`.
Numbers: `phase0/results/summary_0b.json`; per-item scores: `phase0/results/medsam_*.csv`, `vlm_*/scores.csv`.

## Summary

1. **Models do make topology errors that IoU misses.** MedSAM, given the exact box:
   - RIGA discs: 2.3 % have IoU ≥ 0.8 with the wrong topology (19 % in one subset).
   - LIDC nodules: 7.7 % of decent outputs (IoU ≥ 0.5) have the wrong topology.
2. **Topology errors concentrate where experts disagree.** On LIDC, 11.4 % of outputs are wrong on those slices vs 4.9 % elsewhere.
3. **The set-valued reward matters, measured.** Where radiologists disagree, a single-reader ground truth penalizes 46 % of outputs, a majority rule 25 %, the set-valued rule 11 %. Over half (54 %) of majority-rule penalties land on an answer a radiologist drew.
4. **A small VLM cannot output polygons zero-shot.** Phase 1 must start with supervised fine-tuning on polygon output.
5. **RIGA's expert ambiguity is boundary-only** (raw outlines confirm it). RIGA remains the band testbed; LIDC is the topology testbed.

## Q5 — raw expert outlines (RIGA, 699 images × 6 experts)

| | |
|---|---|
| Drawings that are two clean nested closed curves | 4,172 of 4,194 (99.5 %) |
| Others | 20 with a pen gap, 2 with touching curves |
| Certificate failures | 0 |
| Topology disagreement between experts | 0 % (disc, cup, rim) |
| Ambiguous band, median | disc 17 %, cup 55 %, rim 42 %: same as the convex-hull mirror |

The convex-hull mirror did not hide topology disagreement. Extraction found two non-pen artifacts in the raw data
(background change, a solid patch); they are filtered out.

## Q2 — MedSAM, box-prompted

Boxes come from the annotators' union (+5 %, or +5 px), so these are lower bounds on real-use errors.
"Topology wrong" = the output's (β0, β1) matches no annotator, after the size filter (100 px RIGA, 5 px LIDC).

### RIGA (698 images, 6 experts)

| Structure | IoU median | Topology wrong | IoU ≥ 0.8 **and** topology wrong | Error types |
|---|---|---|---|---|
| Disc | 0.92 | 17 (2.4 %) | **16 (2.3 % of all)** | 13 extra holes, 5 extra parts |
| Cup | 0.67 | 5 (0.7 %) | 0 | 3 extra holes, 2 extra parts |
| Rim | 0.78 | 21 (3.0 %) | 5 (0.7 %) | 11 extra holes, 9 extra parts, 4 missing holes |

The wrong-topology discs have median IoU 0.90. On the images, MedSAM **punches holes where blood vessels cross the disc**
and leaves detached fragments. Errors concentrate in **Magrabia Male: 9 of 47 discs (19 %)**; MESSIDOR is 4 of 460.

### LIDC (214 CT series, 4,087 nodule slices with ≥ 2 readers)

MedSAM sees a zoomed window around each nodule (3× the box, at least 64 px). On full 512 × 512 slices it lost small
nodules (IoU 0.22 for nodules < 300 px vs 0.28 zoomed). Full-slice results for the first 24 series are kept in
`phase0/results/medsam_lidc_fullslice/`.

| Slices | n | MedSAM IoU median | Reader-vs-readers IoU median | Topology wrong | Of IoU ≥ 0.5: wrong | Of IoU ≥ reader level: wrong |
|---|---|---|---|---|---|---|
| All | 4,087 | 0.33 | 0.74 | 220 (5.4 %) | 70 / 908 (7.7 %) | 9 / 104 (8.7 %) |
| Readers agree on topology | 3,761 | 0.33 | 0.75 | 183 (4.9 %) | 63 / 822 (7.7 %) | 6 / 69 |
| Readers disagree | 326 | 0.36 | 0.66 | **37 (11.4 %)** | 7 / 86 (8.1 %) | 3 / 35 |

Error types (all slices): 115 extra holes, 85 extra parts, 28 missing parts, 3 missing holes.

MedSAM is far below reader agreement on nodules (0.33 vs 0.74). Most LIDC outputs are simply poor segmentations.
Among the outputs that are decent, about 8 % still have the wrong topology.

### What the set-valued topology reward changes (LIDC, 326 slices where readers disagree)

| Topology target | Outputs penalized |
|---|---|
| One random reader's answer (common practice) | 151 (46 %, expected) |
| Majority answer | 81 (25 %) |
| **Any reader's answer (set-valued)** | **37 (11 %)** |

44 of the 81 majority-rule penalties (54 %) hit an answer that at least one radiologist drew.

## Q6 — Qwen3-VL-2B, zero-shot (RIGA)

**Polygon prompt**, 30 images:

| | Disc | Cup |
|---|---|---|
| Valid simple polygon | **1 / 30** | **2 / 30** |
| Unparseable (reply ran past 500 tokens, JSON never closed) | 20 | 20 |
| Self-touching / out of range (> 1000) / degenerate | 9 | 8 |
| IoU median | 0.00 | 0.00 |

Replies are degenerate coordinate runs, e.g. a vertical line `410,380 → 410,480 → … → 410,2180`. 30 images gave only
16 distinct replies, independent of the image. Median 24 s per image.

**Box prompt** (native `bbox_2d` grounding, scored as the inscribed ellipse), 30 images:

| Structure | Parsed + valid | IoU median | Topology wrong |
|---|---|---|---|
| Disc | 30 / 30 | 0.81 | 0 |
| Cup | 30 / 30 | 0.35 | 0 |
| Rim | — | 0.47 | **14 (47 %)**: missing holes (14), missing parts (4) |

The model finds the disc but not the cup. It often gives the cup the same box as the disc, so the rim loses its hole.

## Decision gates

| Gate | Result | Decision |
|---|---|---|
| High-IoU-wrong-topology ≥ ~5–10 % for some model × dataset | LIDC: 7.7 % of IoU ≥ 0.5 outputs, 8.7 % at reader-level IoU. RIGA Magrabia Male: 19 %. Overall RIGA disc: 2.3 % | **Go to Phase 1.** The error exists and IoU does not see it. |
| Set-valued reward needed? | Majority/single-reader targets wrongly penalize 54 % / ~4× more outputs on ambiguous slices | **Keep it as the headline claim.** |
| VLM polygon output usable zero-shot? | No: 1–2 valid polygons of 30, IoU 0 | **Phase 1 = SFT on polygons first, then GRPO.** |
| Raw outlines hide disagreement? | No | RIGA = band testbed, using our own fills. |

## Caveats

- Box prompts from the annotators' union make MedSAM's job easier than real use; real error rates are higher.
- LIDC reads are unblinded, so reader disagreement is a conservative estimate.
- The LIDC zoom (3×, ≥ 64 px) was chosen after seeing full-slice IoU on 24 series; it is a baseline fix, not tuned further.
- The band score is harsh for tiny objects: MedSAM's spill outside all readers is large relative to a 100-px nodule
  (LIDC band median ≈ 0). The band term needs a scale-aware form for small objects in Phase 1.
- The reader-vs-readers IoU yardstick compares each reader to the majority of the others; for 2-reader slices that is one other reader.
- VLM tests use 30 images: indicative, not precise.

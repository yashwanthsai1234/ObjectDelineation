> Frozen record (written 2026-10-18). Versions, budgets and paths were true then. Current facts: docs/README.md.

# Phase 0 — topology and disagreement census

Phase 0 trains nothing. It runs an exact geometry checker over expert annotations
to find out whether the problem the RL reward targets is real, and how large it is.

## Questions

| # | Question | Needs | Status |
|---|---|---|---|
| Q1 | How often do experts disagree on **topology** (objects β0, holes β1, presence), not only on boundaries? | Multi-annotator masks | Phase 0a (now) |
| Q2 | How often does a current model output **high IoU but wrong topology**? | Source images + a model | Phase 0b (needs images) |
| Q3 | Does the checker agree with itself on real data (certificate)? | Any masks | Phase 0a (now) |
| Q4 | How wide is the ambiguity **band** (union minus core) relative to the object? | Multi-annotator masks | Phase 0a (now) |

## Data

| Dataset | What we have | Annotators | Licence | Notes |
|---|---|---|---|---|
| RIGA via Segmentation-Ready-RIGA (Zenodo 23158636) | 4,470 disc/cup masks, no fundus images | up to 6 per image | CC BY 4.0 | Masks are **filled convex hulls** of the drawn outlines (see `preprocess_RIGA.py` upstream). Disc and cup are convex by construction, so disc/cup topology disagreement here is a processing artifact. Valid uses: band width, cup-to-disc ratio spread, rim breaks (cup reaching the disc edge). |
| LIDC-IDRI radiologist XML (TCIA `LIDC-XML-only.zip`) | 1,319 reading files: per-slice nodule contours + exclusion contours | up to 4 per CT | CC BY 3.0 | Contours are pixel lists on a 512×512 grid, so masks can be rasterised without the 100+ GB of CT. Has real topology content: exclusion regions (holes), multipart outlines, presence/absence per reader. |

Raw RIGA outline images (Deep Blue) and CT/fundus images are behind a browser
check or too large for this disk; they are Phase 0b.

## Checker v0 (Python reference, now `src/objdel/checker/`)

Convention used everywhere: **8-connected foreground, 4-connected background**.

1. Bit-quad histogram (16 codes, one vectorised pass) → χ8, area, perimeter, corners (Gray 1971).
2. β0 = 8-connected foreground components; β1 = 4-connected background components not touching the border.
3. Polygonizer: boundary tracing with left turn at saddles → closed rectilinear rings, outer (shoelace > 0) and holes (< 0), canonical start, ring hash.
4. Certificate: β0 − β1 = χ8; rings outer/hole counts = β0/β1; Σ ring areas = area; Σ ring vertices = corners.
5. Noise filter: drop objects and fill holes smaller than `min_size` pixels; report raw and filtered.
6. Annotator set: core (all), union (any), majority vote, band score, topology histogram with frequency-weighted topology score, typed error report.

Tests (now `tests/`): brute-force random grids against plain flood fill, known shapes (ring, nested ring, saddles, border contact), annotator-set scoring reproducing the artifact demo.

## Procedure

1. Build checker v0 and pass all tests.
2. RIGA: for each image and expert, extract disc, cup, rim (disc minus cup).
   Record β0/β1 per structure, certificate, areas, vertical cup-to-disc ratio (vCDR).
3. RIGA per image: rim topology agreement across experts; band width for disc and cup; vCDR range; pairwise IoU.
4. LIDC: parse every reading; rasterise each nodule ROI per slice (inclusion fills, exclusion subtracts).
5. LIDC: group the readers' annotations into nodule clusters by voxel overlap.
6. LIDC per cluster: number of readers who outlined it, marked it as small (<3 mm, point only) or did not mark it;
   per-reader 3D β0 and per-slice β0/β1; exclusion (hole) use.
7. Run both censuses raw and with the noise filter. Record every certificate failure.
8. Write `experiments/phase0/RESULTS.md`: one table per dataset, the decision gates below, caveats.

## Decision gates (starting values, to be revised after seeing distributions)

| Result | Decision |
|---|---|
| Any certificate failure | Fix the checker. Nothing else proceeds. |
| Topology disagreement ≥ ~5 % of objects after filtering | Keep set-valued topology reward as a headline claim. |
| Topology disagreement < ~2 % | Topology target = majority; the band carries the ambiguity contribution. |
| Band width small (core ≈ union) | Band reward adds little; reconsider dataset. |
| (0b) High-IoU-wrong-topology ≥ ~5–10 % for a baseline | Strong case for Phase 1 training. |

## Phase 0b (after 0a)

- Fetch RIGA fundus images and raw outlines (manual browser download from Deep Blue) → true outline topology + model baseline (Q2).
- One SA-1B shard: hole and multipart frequency in large-scale masks.
- A small model baseline (e.g. SAM or a small VLM) on RIGA images.

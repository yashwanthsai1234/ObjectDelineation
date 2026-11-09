# Phase 0a results — topology and disagreement census

Run 2026-10-18 with checker v0 (now `src/objdel/checker/`, 19 tests passing). Plan: `experiments/phase0/PLAN.md`.
Raw outputs: `experiments/phase0/results/{riga,lidc}/` (`summary.json` + per-item CSV).

## Q3 — does the checker agree with itself?

| Data | Masks certified | Certificate failures |
|---|---|---|
| RIGA (disc, cup, rim × 2 filter settings) | 26,820 | **0** |
| LIDC (every outlined nodule slice) | 40,189 | **0** |

Every real mask passed all four checks: β0 − β1 = χ8, ring counts = Betti numbers, ring area = bit-quad area = pixel count, ring vertices = bit-quad corners.

## LIDC-IDRI — real topology disagreement (Q1)

1,018 CT series, 4 radiologists each (final unblinded reads). Contours rasterised on the 512×512 grid; nodules clustered across readers by shared voxels.

Clustering sanity check: 2,654 nodules ≥ 3 mm found, 907 outlined by all four readers. Armato et al. 2011 report 2,669 and 928.

**Presence (is the object there?)**, 2,652 nodules from 4-reader series:

| Readers who outlined it | 1 | 2 | 3 | 4 |
|---|---|---|---|---|
| Nodules | 775 | 489 | 481 | 907 |

Only 34 % are outlined by all four. In 1,114 nodules at least one reader marked nothing. In 631 the other readers marked it only as "< 3 mm", which is a size call, not an existence call.

**Shape topology**, 1,877 nodules outlined by ≥ 2 readers:

| Disagreement among the readers who outlined it | Raw | Filtered (objects/holes < 5 px removed) |
|---|---|---|
| Any topology disagreement (below, or a reader split it into two nodules) | 225 (12.0 %) | **114 (6.1 %)** |
| Hole on some slice: some readers yes, some no | 179 (9.5 %) | 65 (3.5 %) |
| Several parts on some slice: some readers yes, some no | 70 (3.7 %) | 57 (3.0 %) |
| Different number of 3D components | 32 (1.7 %) | 16 (0.9 %) |
| Split into two nodules by one reader | 17 (0.9 %) | 17 (0.9 %) |
| Slice-level (β0, β1) differs, slices drawn by ≥ 2 readers | 632 / 10,916 (5.8 %) | 317 / 10,898 (2.9 %) |

**Holes always come with disagreement.** After filtering, 65 nodules have a hole in at least one reader's outline, and in all 65 at least one other reader drew no hole. 38 of these were outlined by all four readers.

Raw holes are partly raster artifacts: 39 of the 190 raw hole nodules have no exclusion contour. After the 5 px filter, 6 of 65 remain without one. **Use the filtered column.**

## RIGA — boundary disagreement, the band (Q4)

745 fundus images, 6 ophthalmologists each, 4,470 masks.

| Structure | Topology disagreement | Ambiguous band (union − core) / union, median [p10–p90] | Mean pairwise IoU, median |
|---|---|---|---|
| Disc | 0 % | 17 % [12–24 %] | 0.92 |
| Cup | 0 % | **56 % [35–79 %]** | 0.70 |
| Rim (disc − cup) | 0 % | 41 % [33–52 %] | 0.79 |

Vertical cup-to-disc ratio across the six experts: median spread 0.16. **29 % of images spread ≥ 0.2.**

Zero topology disagreement is **by construction**. This mirror fills the convex hull of each outline, so every disc and cup is one convex blob and every rim is one ring. RIGA here measures boundary ambiguity only.

## Decision gates

| Gate | Result | Decision |
|---|---|---|
| Certificate failures | 0 of 67,009 | Pass. The checker is trusted as a reward. |
| Topology disagreement ≥ ~5 % | LIDC 6.1 % filtered (12 % raw); presence 66 % | **Keep the set-valued topology reward.** LIDC is its testbed. |
| Band width not small | RIGA cup band 56 %, vCDR spread ≥ 0.2 in 29 % | **Keep the band reward.** RIGA is its testbed. |
| Q2 high-IoU-wrong-topology | Not measured | Phase 0b. Needs images and a model. |

## What changed in the plan

- The first training targets split by claim. **RIGA tests the band** (cup boundary, vCDR). **LIDC tests set-valued topology** (holes, parts, presence). RIGA alone cannot test topology from this mirror.
- **For LIDC, topology disagreement is mostly about holes and presence, not extra parts.** The reward's hole term and its "object present" term matter most.
- Expect the set-valued term to touch ~6 % of nodules. A training set for that claim should over-sample them.

## Caveats

- LIDC raster convention: inclusion = filled polygon with outline; exclusion = interior only. A different convention shifts the raw numbers; the filtered numbers are less sensitive.
- LIDC reads are unblinded: readers saw each other's marks before finalising. Real independent disagreement is likely higher.
- The nodule cluster box allows a 3 px / ±1 slice margin when matching "< 3 mm" point marks. That's an approximation.
- No pixel spacing without the CT headers, so sizes are in pixels/voxels, not mm.

## Phase 0b — next

1. RIGA fundus images and raw outlines from Deep Blue (browser download: the site blocks scripts). Gives true outline topology plus Q2.
2. LIDC CT slices for the 1,877 multi-reader nodules only, via the TCIA API, not the full 124 GB.
3. A baseline model on those images → Q2 (high IoU, wrong topology).
4. One SA-1B shard → hole/multipart frequency at scale.

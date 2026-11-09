> Frozen record (written 2026-10-18). Versions, budgets and paths were true then. Current facts: docs/README.md.

# Phase 0b — do models make topology errors that IoU misses?

Phase 0a measured experts (`experiments/phase0/RESULTS.md`). Phase 0b measures models, still with no training.
Machine: Apple M3, 16 GB, inference only.

## Questions

| # | Question |
|---|---|
| Q2 | How often does a model output a mask with **good IoU but a topology no expert gave**? (the decision number) |
| Q5 | Did the convex-hull RIGA mirror hide real topology disagreement in the experts' raw outlines? |
| Q6 | Can a small VLM output usable polygons with no training? (parse rate, valid rings, IoU) |

## Data

| Data | Items | Source |
|---|---|---|
| RIGA originals + raw expert outlines | 699 images × 6 experts (MESSIDOR 460, Magrabia 94, BinRushed 2–4 145) | Deep Blue, downloaded by hand. BinRushed1 has no originals, so it is excluded. |
| LIDC CT slices | 107 series with a reader topology disagreement + 107 random control series from the 800 multi-reader series (seed 0) | TCIA NBIA API. One series at a time; keep only the nodule slices (as .npz), delete the DICOM. |

## Steps

1. **Raw outline census (Q5).** Expert outline = annotated image − original (upstream rules per subset).
   Find closed curves with the checker: enclosed background regions of the outline. A clean annotation has
   2 closed curves (disc, cup) nested, giving 2 enclosed regions. Record open curves, extra strokes, touching curves.
   Fill the curves properly (no convex hull) → true disc/cup/rim masks → rerun the 0a topology and band census.
2. **LIDC fetch.** Download the selected series, match XML contours to slices by SOP Instance UID, save lung-window
   slices plus a crop box around each nodule cluster.
3. **Baseline A: MedSAM (box prompt).** RIGA: disc box and cup box from the experts' union, +5 % margin.
   LIDC: box per slice from the readers' union, +5 px. Output masks, threshold 0.5.
4. **Baseline B: Qwen3-VL-2B-Instruct (zero-shot).** RIGA, a fixed random 100 images. Prompt for the disc and cup
   as polygons in 0–1000 coordinates. Parse, check validity, rasterise.
5. **Score with the checker** against each item's annotator set (`AnnotatorSet`, min_size 5 px for LIDC, 100 px for RIGA):
   IoU vs majority, band, topology frequency, typed errors. RIGA also checks the rim: disc − cup must be one ring with one hole.
6. **Q2 metric:** share of outputs with IoU ≥ 0.8 **and** a topology no expert gave (topology frequency = 0).
   Also report it at IoU ≥ the median inter-expert IoU for that structure.
7. Write `experiments/phase0/RESULTS_0b.md`.

## Box prompts leak the answer's location

Box prompts come from the experts' masks. That is standard in promptable-segmentation evaluation, and it makes the
baseline optimistic: real use would have noisier boxes. Topology errors under a perfect box are a lower bound.

## Decision gates

| Result | Decision |
|---|---|
| High-IoU-wrong-topology ≥ ~5–10 % for some model × dataset | Phase 1 has a real target. Go. |
| ≈ 0 everywhere | Choose a harder task (thin structures, multipart) before training. |
| VLM parse/valid rate low | Phase 1 starts with SFT on polygon output before GRPO. |
| Raw outlines show topology disagreement the mirror hid | Use our own fill for RIGA in Phase 1. |

## Out of scope

SA-1B shard (no licence link yet; 11 GB). GPU checker. Any training.

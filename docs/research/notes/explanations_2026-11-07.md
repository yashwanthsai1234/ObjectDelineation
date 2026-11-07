# Explanations given to the owner, 2026-11-07

Collected from the day's conversation, with the corrections made during it. Related notes: cpu_gpu_trace.md (full
technical trace), image_opening_cpu_vs_gpu.md (opens, CPU vs GPU, the fix's measurements), image_opens_explained.md
(the card/stack draft), gpu_scoring.md (design of scoring on the GPU), report_cpu_gpu.md (CPU vs GPU in
one RL step).

Words: **image** (one 512 × 512 LIDC lung slice, one PNG file), **outline** (what the model draws), **card** (one image
plus the question; one card gives one outline), **stack** (the 256 cards of one step), **pile** (16 outlines; the GPU
learns from 16 piles per step), **open** (read a PNG from the disk and unpack it into pixels).

## 1. G-probe

A short, pre-registered RL check: 20 GRPO steps on the SFT checkpoint e2 (R1's Qwen3-VL-2B branch, decayed at SFT epoch
2). It measures the composite reward on 128 held-out items twice before training (to see the noise) and once after.
It passes when the gain is clearly above that noise. It answers "does RL improve this model at all?" before paying
for a full run. It saves no tuned weights.

R2 result (RTX 5090, billed $0.550): PASS. Reward 0.694 / 0.679 → 0.769, gain +0.082 (CI 0.067-0.097, step-0 noise
0.018); legal outlines 94.5 % / 93.5 % → 98.2 %; 56.6 s per step.

## 2. Exam metrics

The model's report card on held-out images it never trains on. They grade the model; they do not teach it.

| Question | Metric |
|---|---|
| Readable answer, legal outline? | format rate, legal rate |
| Overlap with the doctors' outline? | IoU |
| Edges close to the doctors' edges? | F-score (band) |
| Right size? | area score |
| Right parts and holes? | topology |
| Does one of the 8 outlines pass? | pass rate, pass@k |
| Edges within a few pixels of a doctor's edge? | surface distance (NSD) |
| Agrees with the doctors' blended outline? | soft Dice |
| Are the 8 outlines varied like the 4 doctors? | spread match (GED), overlap between the model's own outlines |
| Copies itself? | duplicate rate |

The first five are the RL reward's own terms (now on the GPU). The last four compare outlines with each other and with
each doctor; they still copy pictures to the CPU, once per image. They run a few times per run (the G-probe's three
probe draws, checkpoint selection, an optional held-out audit), never per training step.

## 3. Where the 256 cards come from

1. Once, at the start, our code builds a table: one row per training image (question, image file name, size). No
   pixels.
2. Every step, the training library (TRL) shuffles the rows, takes 32, writes each row number 8 times (256 numbers),
   and writes that list 16 times (4,096 numbers).
3. The data loader on the CPU reads 256 numbers at a time; each row it reads opens that row's image.

Why 8: we asked for 8 outlines per image, and the model writes one outline per card. Why 32: 256 ÷ 8. Why 16: the
model learns from the 256 outlines in 16 piles of 16, and the library asks for a stack per pile.

## 4. The 16 stacks, step 1 included

At the start of every step, step 1 included, the library collects all 16 stacks before the model draws. The model draws
256 outlines from stack 1 in one go: the CPU hands all 256 cards to the drawing engine (vLLM) in a single call, and the
GPU moves all 256 outlines forward one word at a time together, not 256 turns. The 16 piles are cut from those 256
outlines. Stacks 2-16 are thrown away on the CPU and never reach the GPU.

Correction made during the day: the waste is not "32 images opened per pile, 2 kept". Every stack opens all 256 cards'
images, and 15 whole stacks go unused.

## 5. The fix, part 1 and part 2

- Part 1: one stack per step (3,840 opens saved, most of the wait).
- Part 2: inside the stack, each image opened once and shared by its 8 cards (224 opens saved, under 1 s).

Neither changes results; both are speed. Built 2026-11-07 as dataset.OpenOnce: 4,096 → 32 opens and 12.9 s → 0.14 s
per step on this Mac with the real loader code.

## 6. Opening on the CPU vs the GPU

PNG unpacking is a chain, so the CPU is the right place for our 512 × 512 images. At 100,000 × 100,000 the answer flips:
tiled formats, read only the needed tiles, unpack them in parallel on the GPU. Numbers in image_opening_cpu_vs_gpu.md.

## 7. Corrections made during the day

- The step-time split "14 s writing + about 40 s learning" was wrong: learning is about 14 s, and about 13-20 s was the
  CPU opening images (report_cpu_gpu.md §1).
- "About 5 ms per open" was an estimate from an unmeasured 20 s; the Mac measures 2.6 ms per open, and the loader
  measures 12.9 s per step before the fix.
- "Fix the images" was the wrong word: the images were always correct; the issue is speed only.
- One word per thing: "scan", "photo" and "image" were mixed; this note uses "image" only.

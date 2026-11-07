# What runs on the CPU and what runs on the GPU in one RL step

> **Update, later on 2026-11-07.** Two things below have changed since this report was written:
>
> 1. Scoring (IoU, edges, area, shape) now runs on the GPU, next to the checker (branch `gpu-scorer-open-once`;
>    design in gpu_scoring.md). The pictures no longer go to the CPU for the reward.
> 2. The image-reading waste is fixed (dataset.OpenOnce): 4,096 opens → 32 per step. On this Mac, measured with the
>    real loader code: 12.9 s → 0.14 s per step. The "~20 s" below was an estimate; 12.9 s (Mac) and 12.96 s (P1 box)
>    are the measurements. Details: image_opening_cpu_vs_gpu.md.

This report covers one training step of the G-probe run on 2026-11-07, on an RTX 5090. Each step, the model looks
at 32 lung scans and writes 8 outlines for each one, so 256 answers per step. It then learns from how good those
answers were.

Every number is labelled:

- **measured**: from that run's log, or from the earlier profiling run on a 5090;
- **estimated**: worked out from the code plus measurements;
- **not checked**: could not be verified.

---

## 1. The short version

```
 one step = about 57 seconds (measured)

 ┌───────────────────────────────────────────────────────────────────────────────┐
 │ CPU  reading the same 256 scan images 16 times      ~20 s   GPU does nothing  │  ← biggest waste
 ├───────────────────────────────────────────────────────────────────────────────┤
 │ GPU  writing the 256 answers                        ~12 s                     │
 ├───────────────────────────────────────────────────────────────────────────────┤
 │ CPU  preparing the images (twice) + 1.6 GB copy     ~2-4 s  GPU mostly idle   │
 ├───────────────────────────────────────────────────────────────────────────────┤
 │ GPU  "reference" read of all answers                ~4-6 s                    │
 ├───────────────────────────────────────────────────────────────────────────────┤
 │ GPU  learning from the answers                      ~14 s                     │
 ├───────────────────────────────────────────────────────────────────────────────┤
 │ CPU  scoring the answers (IoU, edges, area, shape)  ~0.02 s, hidden           │  ← must move to GPU
 └───────────────────────────────────────────────────────────────────────────────┘
```

Four things to know:

1. **Scoring is on the CPU today.** The GPU checker draws each outline into a picture and checks it for defects.
   Then the pictures are copied to the CPU, and the CPU works out every reward number: overlap with the expert
   outline, edge closeness, size, and the count of parts and holes. That breaks the rule that checking and scoring
   stay on the GPU. The rewrite now in progress fixes it.
2. **Scoring is not what makes a step slow.** It takes about 0.02 s and runs while other work is happening.
3. **The biggest slowdown is the CPU reading images.** Every step, the CPU opens and decodes the same 256 scan
   images 16 times, throws away 15 of the 16 copies, and the GPU waits about 20 s for it (estimated; the earlier
   profiling run measured 13 s on a faster CPU).
4. **I gave you a wrong breakdown earlier.** I said the step was "14 s writing + about 40 s learning". The learning
   part is about 14 s. About 20 s is the CPU reading images.

---

## 2. One step, in order

```
 ─────────────────────────────────────────────── time ─────────────────────────────────────────────────►

 A. CPU: fetch the scan images ×16 (only 1 copy used)  │████████████████████│              ~20 s  (est.)
                                                        GPU idle ↑
 B. CPU: prepare images for the model (1st time)       │██│                                  ~1-2.5 s
 C. GPU: copy current weights into the answer-writer   │█│                                   < 1 s (not checked)
 D. GPU: write 256 answers                             │████████████│                         ~12 s  (measured)
         ├─ starts scoring in the background ─┐
 E. CPU: prepare images again (2nd time)       │       │██│                                  ~1.5-3 s
 F. CPU→GPU: copy 1.6 GB of image numbers      │       │▌│  (the GPU waits for the copy)
 G. GPU: reference read of all 256 answers     │       │█████│                              ~4-6 s
 H. scores needed ◄── already done ────────────┘       │ wait 0.004 s (measured)
 I. GPU: turn scores into "better or worse than its group"
 J. GPU: learn, in 16 slices of 16 answers             │██████████████│                       ~14 s  (measured)
 K. GPU: update the weights (small; LoRA only)         │                                       < 1 s
 L. CPU: log numbers, start the next step
```

**Example.** Take one scan, nodule n001, slice 54. In this step:

- It is opened and decoded from disk 16 times (A).
- It is prepared for the model twice (B, E), and a third time inside the answer-writer (not checked).
- The model writes 8 outlines for it (D).
- The GPU draws each outline into a 50 × 50 picture and checks it.
- The pictures go to the CPU, which compares each one with the experts' outlines. Say outline 3 gets overlap 0.61,
  edge score 0.70, size score 0.95, and the right shape (one part, no holes).
- Those numbers go back to the GPU and become a single "outline 3 was better than its group of 8" number.
- Learning (J) then makes outline 3's words a little more likely next time.

---

## 3. Every part, CPU or GPU

| Part of the step | Where | Time per step | Notes |
|---|---|---|---|
| Fetch scan images from disk and decode them | **CPU** | ~20 s (est.) | Done 16 times per step; 15 copies thrown away. The GPU waits. |
| Prepare images for the model, 1st time | **CPU** | ~1-2.5 s | Result thrown away; the answer-writer prepares them itself. |
| Copy the latest weights into the answer-writer | GPU | under 1 s (not checked) | Done every step: ~4.3 GB copied inside the GPU. |
| Write the 256 answers | **GPU** (+ CPU bookkeeping) | ~12 s (measured) | Set by the longest answer: 13 of 20 steps had one hitting the 1,536-word limit. |
| Prepare images, 2nd time | **CPU** | ~1.5-3 s | The images the learning step reads. |
| Copy those images to the GPU | CPU → GPU | ~0.1-0.2 s | 1.6 GB each step, copied the slow way (the GPU waits). |
| Read the answers back as text | CPU | tiny | |
| Turn answers into outlines (read the text) | **CPU** | ~0.002 s | Must stay on the CPU (it is text). |
| Pack outlines for the checker | **CPU** | ~0.002 s | |
| Draw outlines into pictures + find defects (**the checker**) | **GPU** | milliseconds | One call for all 256 answers. |
| Copy the pictures back to the CPU | GPU → CPU | ~190 KB | **Should not happen.** Removed by the rewrite. |
| Overlap with expert outline (IoU) | **CPU** | ~0.01 s for all 4 scores | **Moving to GPU.** |
| Edge closeness score (F-score) | **CPU** | (incl. above) | **Moving to GPU.** |
| Size score (not too big, not too small) | **CPU** | (incl. above) | **Moving to GPU.** |
| Shape: parts and holes, compared with the experts | **CPU** | (incl. above) | **Moving to GPU.** |
| Hand the pictures to 16 CPU helper processes and collect results | **CPU** | a few ms | Only exists because scoring is on the CPU; goes away. |
| Combine scores into one number per answer, compared within its group of 8 | **CPU** | ~0.0002 s | **Moving to GPU.** |
| Copy the reward numbers to the GPU | CPU → GPU | 8 small copies | GPU waits briefly for each. Goes away once scores are made on the GPU. |
| Reference read (the model with training switched off, for the "don't drift" penalty) | **GPU** | ~4-6 s (est.) | |
| Learn: 16 slices, forward and backward | **GPU** | ~14 s (measured) | Uses rope positions computed once per step on the CPU (~0.04 s). |
| Update the weights | **GPU** | small | |
| Log numbers | CPU | small | |

---

## 4. Scoring today and after the rewrite

**Today:**

```
 answer text ──CPU──► outline ──► GPU checker: draw + find defects
                                        │
                                        │  pictures copied GPU → CPU   ✗
                                        ▼
                            CPU: overlap, edges, size, shape          ✗
                            CPU: combine into one number per answer   ✗
                                        │
                                        │  numbers copied CPU → GPU   ✗
                                        ▼
                                  GPU: learning
```

**After the rewrite:**

```
 answer text ──CPU──► outline ──► GPU checker: draw + find defects
                                              + overlap, edges, size, shape   ✓
                                              + combine per answer            ✓
                                        │
                                        │  only ~10 numbers per answer stay on the GPU
                                        ▼
                                  GPU: learning
```

One part stays on the CPU by nature: reading the answer text and turning it into an outline. It is text, it takes
about 2 milliseconds, and the checker cannot start without it. **Your call:** whether this counts as "checking".

The expert outlines must sit on the GPU for this. They are about 0.33 GB for the training set, uploaded once at the
start.

---

## 5. Everything that crosses between CPU and GPU in one step

```
             CPU ═══════════════════════════════════════ GPU
 images      1.6 GB  ─────────────────────────────────►   (copied the slow, waiting way)
 answers     ~9 MB   ─────────────────────────────────►
 positions   ~11 MB  ─────────────────────────────────►
             ~27 MB  ◄─────────────────────────────────   positions, small reads
 pictures    190 KB  ◄─────────────────────────────────   ✗ goes away with the rewrite
 scores      ~10 KB  ─────────────────────────────────►   ✗ goes away with the rewrite
```

Besides the big copies, the CPU makes the GPU stop and report back about **1,300 times per step** (estimated from the
code). Most are small "what is this number?" questions in the learning and reference reads. The answer-writer adds
about one per word it writes (not checked). Each one makes the GPU finish everything queued before the CPU can go
on.

---

## 6. Other problems found

1. **The probe's scoring workers start their own GPU sessions.** When G-probe scores its 128 test scans, the
   checker runs inside each of up to 16 CPU helper processes. Each one may open its own GPU session, and those
   sessions stay open, using GPU memory, through all 20 training steps. Not checked on the box. The rewrite removes
   the helper processes.
2. **The frozen model may shift slightly each step.** Every step, the trained adjustments are added into the
   model's weights to update the answer-writer, then subtracted again. In the reduced precision used (bf16), adding
   and subtracting does not give back exactly the same numbers. Over 20 steps, the "frozen" weights (and the
   reference used for the drift penalty) can drift by rounding. Size not measured.
3. **Each scan image is prepared three times per step** (twice by the trainer, once inside the answer-writer, not
   checked) and decoded 16 times.
4. **The vision part of the model runs 16 times per step on 2 images each,** instead of once on the 32 different
   images.
5. **The "copy images without waiting" setting has no effect.** It only applies to a part of the pipeline that
   carries no images.

---

## 7. Where speed could come from (for you to decide)

| Idea | Saves per step (est.) | Moves work to |
|---|---|---|
| Read and decode each scan image once per step, not 16 times | ~15-20 s | stays CPU, but the GPU stops waiting |
| Prepare images once, not three times; copy them the non-waiting way | ~2-4 s | — |
| Move scoring to the GPU (in progress) | ~0.02 s, but fixes the rule and removes the helper processes | GPU |
| Stop the long answers early (one 1,536-word answer holds up all 256) | up to ~6 s | — |
| Fewer GPU stop-and-report points in learning | unknown until measured | — |

---

## 8. How this was worked out

- I followed one step through the training code line by line: the trainer library, our own trainer changes, the
  reward code and the checker.
- I matched the parts against the run's log (step time, writing time, scoring wait) and against the earlier
  profiling run on a 5090.
- The answer-writer's (vLLM's) inner workings were not available on this Mac. Everything about them is marked "not
  checked".
- The ~20 s image-reading figure is estimated. It is wall time minus the trainer's own step time (22.6 s, measured),
  matched against the 13 s the earlier profiling run measured on a faster CPU. The run did not log it directly.

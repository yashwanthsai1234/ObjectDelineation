# Opening images: CPU or GPU, and the 16-stack waste (2026-11-07)

Words used here:

- **image**: one 512 × 512 lung slice, stored as one PNG file on the disk.
- **open**: read the PNG file from the disk and unpack it into pixels. The model can only use pixels.
- **card**: one image plus the question. The model writes one outline from one card.
- **stack**: the 256 cards of one step (32 images × 8 cards).
- **pile**: 16 outlines. The GPU learns from one pile, then the next: 16 piles per step.

## 1. What one open costs (measured on this Mac, 2026-11-07)

| | |
|---|---|
| PNG file on disk | 146 KB |
| After unpacking | 512 × 512 pixels × 3 colours = 786 KB |
| One open | 2.6 ms (one CPU core) |

## 2. The waste: 16 stacks, 1 used

At the start of every step, step 1 included, the training library asks for 16 stacks, one per pile, before the model
draws anything (transformers `get_batch_samples`, gradient accumulation 16). TRL's sampler makes every stack the same
256 cards. The model draws from stack 1 only, and all 16 piles are cut from stack 1's 256 outlines. Stacks 2-16 are
thrown away on the CPU and never reach the GPU. Each card opens its image again, even when 8 cards share one image.

```
 BEFORE                                         AFTER (dataset.OpenOnce)
 stack  1: 256 cards, 256 opens ──► model       stack 1: 32 opens, each image shared by its 8 cards ──► model
 stack  2: 256 cards, 256 opens ──► thrown away stacks 2-16: the same stack handed back, 0 opens
   ...
 stack 16: 256 cards, 256 opens ──► thrown away
 4,096 opens per step                           32 opens per step
```

Measured with the real loader code (TRL 1.17.1, transformers, Accelerate) on the prepared LIDC training rows, this Mac:

| | Opens per step | CPU time per step, GPU waiting |
|---|---|---|
| Before | 4,096 | 12.9 s (P1 on the RTX 5090 box: 12.96 s) |
| After | 32 | 0.14 s |

The model gets the same pixels in the same order (test: `tests/test_open_once.py`, every row of every fetch compared).
Accelerate reads one fetch ahead, so a step's 32 opens happen at the end of the previous step's fetch.

The fix is in two parts:

- part 1: one stack per step, not 16 (3,840 opens saved);
- part 2: inside the stack, one open per image, not 8 (224 opens saved).

Both parts are speed only: the images were never wrong.

## 3. Is the CPU or the GPU better at opening a PNG?

The CPU. PNG uses zip-style compression (DEFLATE), and unpacking it is a chain: each piece depends on the one before.
A GPU is fast only when thousands of threads work on separate pieces, so it gains nothing on a chain. torchvision 0.29
(installed here) unpacks JPEG on the GPU (`decode_jpeg(device=...)`), but PNG only on the CPU (`decode_png` has no
device argument). The next stage, preparing the image for the model (resize, cut into patches), also runs on the CPU
and needs the pixels there.

After the fix, opening is 32 × 2.6 ms ≈ 0.08 s of a ≈ 40-57 s step: not worth moving.

## 4. What if the image were 100,000 × 100,000?

Scaled from the measurement above:

| | 512 × 512 (ours) | 100,000 × 100,000 |
|---|---|---|
| Pixels | 262,144 | 10,000,000,000 |
| Memory after unpacking (3 colours) | 0.8 MB | 30 GB |
| One open, one core | 2.6 ms | about 99 s |
| Pillow's safety limit (2 × 89,478,485 pixels) | under | 56 × over: refuses to open |
| Model tokens (Qwen3-VL-2B: patch 16, merge 2, so one token per 32 × 32 pixels) | 256 | 9.8 million (limit 262,144) |

So at that size:

1. The model cannot take the whole image: shrink it or cut it into tiles first.
2. A single PNG is the wrong format: to reach the end you unpack everything before it. Huge images (pathology slides
   are about this size) use tiled formats (tiled TIFF, JPEG 2000, pyramid files), and you read only the tiles you need.
3. Then the GPU wins: each tile unpacks on its own, so thousands unpack at the same time (NVIDIA's GPU decoders for
   JPEG and JPEG 2000 tiles; NVIDIA cuCIM for slide images), and the tiles can stay on the GPU for the model and the
   checker. 3-D voxel volumes of that size would use the same idea with 3-D blocks.

## 5. Does the CPU keep up, so the GPU never waits?

After the fix, the CPU work before the model can draw is 32 opens plus the reads of the rows: 0.14 s per step on this
Mac. The GPU part of a step is tens of seconds, so the wait is under 1 % of a step. Hiding even that would need a
background thread that opens the next step's 32 images while the GPU works. The larger CPU work left is preparing each
image for the model 2-3 times per step (about 2-4 s, estimated, not measured): see report_cpu_gpu.md §7.

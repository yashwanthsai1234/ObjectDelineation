# Issues solved on 2026-11-07

Branch `gpu-scorer-open-once` unless noted. "Mac" numbers were measured in this checkout on real `data/prepared`
LIDC rows; GPU results are in §8.

| # | Issue | Fix | Evidence |
|---|---|---|---|
| 1 | Reward scoring ran on the CPU: the GPU checker drew the outlines, then the pictures went to the CPU for IoU, F-score, area and topology, and the numbers came back | Scoring on the GPU next to the checker; expert tables uploaded once and kept on the device; one call for all 256 outlines (design: gpu_scoring.md) | Same scores, bit for bit, as the old CPU scorer on recorded goldens and a numpy oracle (tests/test_gpu_scoring.py, tests/test_checker_oracle.py) |
| 2 | The CPU checker (polycore) still lived beside the CUDA checker | Archived in archive/cpu_checker; its recorded results stay as test goldens | tests/goldens/polycore.npz, tests/oracles/score_numpy.py |
| 3 | The sanitizer step (memcheck + initcheck) covered 1 test file (6 tests) | Covers the 5 checker and scorer test files in all 3 run configs | tests/test_profile.py::test_memcheck_has_its_sanitizer_and_fails_loudly |
| 4 | 16 stacks of 256 cards made per step, 1 used: 4,096 image opens for 32 images, GPU waiting | dataset.OpenOnce: one stack per step, each image opened once and shared by its 8 cards; FastGRPO's loader uses it | 4,096 → 32 opens, 12.9 s → 0.14 s per step (Mac, real loader path); identical rows (tests/test_open_once.py) |
| 5 | Each image prepared for the model twice by TRL (prompt ids, then pixel values), over all 256 rows | fast.ImageMemo: the image processor runs once per distinct image per step; both passes reuse it | 2.12 s → 0.78 s per step (Mac, real processor); prompt ids and every pixel tensor identical (tests/test_image_memo.py) |
| 6 | EV-6: G-probe's probe sets were sampled outside the trainer | Probe sets sample through the trainer's vLLM colocate (commit 0c83810, branch rl-ready) | experiments/phase2/perf/issues.md |
| 7 | Did RL improve the SFT model at all? | R2 G-probe on e2, RTX 5090, $0.550: PASS, reward 0.694 / 0.679 → 0.769, legal 94.5 % / 93.5 % → 98.2 % | runs/r2/g_probe.json |

## What each fix does not cover

- **Image preparing inside vLLM (pass 3).** vLLM prepares the images again for drawing. Whether its cache removes the
  8 repeats is not checked.
- **The 1.6 GB pixel copy to the GPU.** Pass 2 still builds the pixel block for all 256 rows and copies it to the GPU
  with a waiting copy (about 0.1-0.2 s, estimated).
- **Exam metrics** (G-probe probe draws, checkpoint selection) still copy pictures to the CPU. They run a few times per
  run, never per training step.
- **No tuned weights.** The G-probe saves no weights; a full RL run needs its save step.

## 8. GPU test

(filled in after the Vast test)

# Verdict on final_r12_qa.md (fresh verifier, 2026-11-09)

## Skills used

| Skill | What it changed |
|---|---|
| nvidia-skills:nvidia-skill-finder | No NVIDIA catalog skill applies; it made me check which build (CPU vs CUDA) each D-G3 mask test runs on. |
| engineering:code-review | Read `experiments/cuda_checker/compare.py:21, 165-169`, `src/objdel/checker/__init__.py:6`, `tests/test_cuda_checker.py:35-38`, `tests/test_cuda_checker_gpu.py:46-51`, `tests/test_train_fast.py:24-35`, `experiments/phase2/vast/policy.py:483-490`. |
| engineering:debug | Traced `mask_reference` -> `objdel.checker.polygonize` after B1.4a -> `rings.polygonize` -> `_gpucheck.check_masks` (item 1); traced a failing hard step through the `&&` chain (item 2). |
| engineering:system-design | Checked that G-RC0 serves C11 (a run after Gate B3) yet gates the B3 run (item 2). |
| engineering:architecture | Checked who owns the oracle move (B1.4a, D3, CLEAN11, 00_front 2.2 row) for the item 1 fix's source list. |
| engineering:documentation | Checked the review logs in review-log.md for earlier handling of both items (none: qa-2 of round 11 only placed G-RC0). |
| huggingface-skills:huggingface-best | No model choice in scope; no change. |

## Verdicts (blocking/major items)

- qa-1 | confirmed (major) | After B1.4a `checker/__init__.py` exports `rings.polygonize` (B.md B1.4a "import `polygonize` from `.rings`"), and `compare.py:21` `from objdel.checker import bitquads, polygonize` feeds `mask_reference` (`compare.py:167`), which `tests/test_cuda_checker_gpu.py:50` uses, so D-G3's GPU mask-ring check compares `check_masks` with itself until D3 (after D2) moves the callers; `test_rings.py` (the only independent ring oracle) is CPU-build only, so D-G3 on R-D2 cannot fail on ring fields (an untestable paid gate) and no $0 check flags it.
- qa-2 | confirmed (major) | B.md B3.9 step 2 is a hard step (`[-m, pytest, ...]`, no `~`; policy.py:487-490 joins non-`~` steps with `&&`) before step 3 and `train_rl`, while C.md C.9.4 G-RC0 says only "Tiny model, vLLM 0.32.0" and `tests/test_train_fast.py:24-35` builds that model in memory (vocab 128, no tokenizer/processor files) that `vllm.LLM` cannot load by path; the engineer must guess a load path, I22b only runs `--collect-only` (the `gpu` test skips on the Mac), and a G-RC0 failure on V3 stops the B3 run even though C11 applies only after Gate B3.

## Minor items (filed as minor; not promoted)

- qa-3 | minor as filed | Caught at $0 by I22b's `--collect-only`, as the item itself says.
- qa-4 | minor as filed | U-CF3/U-CF5 are $0 checks that point at the flag list.
- qa-5 | minor as filed | Arithmetic margin of about 0.02 min; budgeting detail, no planned check, but the item stays minor as filed (not a blocking/major item to verify).
- qa-6 | minor as filed | Test-coverage nice-to-have.
- qa-7 | minor as filed | E6 golden 5 fails on the Mac at $0.
- qa-8 | minor as filed | Crash-window edge case for `--merge-last`.
- qa-9 | minor as filed | Round 11 qa-4 carry-over.
- qa-10 | minor as filed | Round 11 qa-5 carry-over.
- qa-11 | minor as filed | Round 11 qa-6 carry-over; `launch_refusal` refuses at $0 (6.1 step 4b).
- qa-12 | minor as filed | Round 11 qa-7 carry-over.
- qa-13 | minor as filed | Stale header and doc size.

Totals: 2 confirmed (qa-1 major, qa-2 major), 0 refuted, 0 downgraded.

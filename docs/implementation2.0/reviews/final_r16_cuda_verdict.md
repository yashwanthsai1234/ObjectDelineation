# Verdict on final_r16_cuda.md (fresh verifier, 2026-11-09)

## Skills used

| Skill | What it changed |
|---|---|
| nvidia-skills:nvidia-skill-finder | Nothing installed (owner approval needed). Confirmed no NVIDIA catalog skill changes the sanitizer-timing question. |
| engineering:code-review | Re-read `src/objdel/reward/cache.py:167-197` to count the EDT calls in one `ExpertCache.build`. |
| engineering:debug | Traced G01 `b_sets` -> `b_tables` -> `ExpertCache.build` -> `rater_tau` + `f5_tolerance` -> `edt`, and D7's move of `edt` to `_gpucheck.edt`. |
| engineering:system-design | Checked the R-D2 `--gpu-check` step chain (step 4 `sanitize` before steps 5-6) to see the cost of one timeout. |
| engineering:architecture | Checked whether any planned $0 check measures sanitizer time (none: Mac has no compute-sanitizer; `--dry` is collect-only). |
| engineering:documentation | One line per item with the quote and proof. |
| huggingface-skills:huggingface-best | No model choice in scope. Nothing changed. |

## Verdicts

- cuda-1 (major): CONFIRMED (major). D.md heading "D8" says "Once B3.1 lands, its `tests/test_scorer_sums.py` joins this corpus; its GPU cases carry `gpu` and `sanitizer`", and the only GPU case in the B.md "B.9.8" sketch is `test_g01_golden`, which runs the whole `b_sets` (10,881 objects; 100.9 s of Mac host work [measured: B.md B.9.5 G01]); each `ExpertCache.build` calls `edt` 6 times in `rater_tau` and 6 times via `_band` in `f5_tolerance` [code: src/objdel/reward/cache.py:167-197], which D7 moves to `_gpucheck.edt` (D.md "D7", cache.py `edt` bullet), so about 130k+ instrumented device calls run under each of four tools inside `timeout 300`; no $0 check measures sanitizer time, no review log handles it (log "cuda 18" only added the file to the corpus), and the doc contradicts its own "small fixed corpus so each tool finishes inside its cap"; a timeout fails D-G5 on paid R-D2 time.
- cuda-2 (minor): kept minor (reviewer's own severity; R-D1 is optional and dropped on the time-line order; not checked further).
- cuda-3 (minor): kept minor (wording of the stated reason; no wrong result; the reviewer agrees the child context is legal).
- cuda-4 (minor): kept minor (time estimate for optional V3 item 1, whose stop rule is a skip).

Counts: confirmed 1 (major), refuted 0, downgraded 0.

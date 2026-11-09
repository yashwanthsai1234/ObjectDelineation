# Verdict on final_r15_cuda.md (fresh skeptical verifier, 2026-11-09)

## Skills used
- nvidia-skills:nvidia-skill-finder: review only; no NVIDIA catalog skill covers engine counters or gate wording; nothing installed.
- engineering:code-review: re-read metrics.py `_samples`, score.py `score_many`, gpucheck.cpp `resident`/`Registry` to check item 1's code claims.
- engineering:debug: traced the post-D7 eval path: `metrics()` -> `_samples` -> `score_many` -> `score_polygons` stays in-process, so the eval's `arena_bytes` includes the checker/scorer arenas.
- engineering:system-design, engineering:architecture: checked whether any other part of D (D.4, D5c, D.9.6, D.9.7) defines "table bytes" or records it; none does.
- engineering:documentation: verdict format, one line per item.
- huggingface-skills:huggingface-best: no model choice in scope; nothing changed.

## Verdicts
- cuda-1: CONFIRMED (major). D.md D.6 row D-G14, S13 `test_b04` and B04 all gate "`arena_bytes` <= the `eval_sets` scratch cap + table bytes", but no record holds table bytes (D7 step 5 and S12 write only `eval_sets_s`, `arena_bytes`, `max_memory_reserved`; D.4 defines only "the current arena bytes" and `live_device_bytes`; `resident()` returns a count [code: cuda/gpucheck.cpp:1778 `i64 resident()`]), and the eval process's `arena_bytes` also holds the `score_polygons` arenas because D7 keeps `score_many` [code: src/objdel/eval/metrics.py:253-256 `score_many(...)`; src/objdel/reward/score.py:113 `_gpucheck.score_polygons(...)`], so the engineer writing `test_b04` must guess both sides; it is read only after the paid V3 run and no earlier review log handles it.
- cuda-2: MINOR (reviewer-rated minor; not re-assessed for upgrade). Missing `gpu_memory_utilization`/`max_model_len` for the post-E4 in-process `VllmSampler`.
- cuda-3: MINOR (reviewer-rated minor; not re-assessed for upgrade). D8 `--gpu-check` does not say whether a sanitizer failure stops the benches.
- cuda-4: MINOR (reviewer-rated minor; not re-assessed for upgrade). `drop_tables`/`resident` bindings keep the GIL while waiting on `e.mu`.

Counts: confirmed 1 (major 1), refuted 0, downgraded 0.

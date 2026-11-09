# Verdict: final review round 7, CUDA (final_r7_cuda.md)

## Skills used

| Skill | What it changed |
|---|---|
| nvidia-skills:nvidia-skill-finder | No catalog skill is needed to check these items. No change. |
| engineering:code-review | I read the cited code lines myself: gpucheck.cpp:1597-1616, reward/cache.py:84-86 and 241-257, reward/score.py:39-45 and 112-115, Makefile:41-44, final_eval.py:23-29, eval/metrics.py:298-311. |
| engineering:debug | I traced the D-G8 old arm import chain: score.py imports .cache, which imports cv2 and scipy (item 1 holds). I traced the V4 TOST input: iou_tc comes from score_polygons terms, which GP1 checks at V4 (item 3 fails). |
| engineering:system-design | I checked the V4 step order (E.8 steps 0-3) against the gates that guard step 3. |
| engineering:architecture | No design choice is in dispute. No change. |
| engineering:documentation | Each verdict names the doc line and the code line. |
| huggingface-skills:huggingface-best | No model choice is in scope. No change. |

## Verdicts

- cuda-1: CONFIRMED major. D.md D.6 says D-G8 uses "Same setup as D-G7" (the old venv `/opt/old` gets only numpy plus a `--no-deps` editable install [code: Makefile:42]). The replay must load caches and call `cache.device()`, and `reward/score.py:45` imports `.cache`, which imports `cv2` and `scipy` at module level [code: reward/cache.py:84-86]. So `--engine old` on the replay fails with ModuleNotFoundError in paid R-D2 time, and the engineer must guess a workaround.
- cuda-2: CONFIRMED major. D.md D.4 defines `label_pairs` as "launches of the scorer's second label pair". D5 step 1 puts the skip inside the kernel, per crop row. `label(which)` always launches three kernels over all `crows` [code: gpucheck.cpp:1597-1601, 1614]. A launch count is the same whether the skip works or not, so I14 and the D5 step 1 stale-flag test cannot be checked as written.
- cuda-3: REFUTED. The claim says V4 has no GPU equality check before the sampler eval. But V4 step 2 (E.8) runs GP1 `cuda_cpu_bit_equal` (`tests/verify/test_verifier_gpu.py`), which compares reward, terms and defects of the box CUDA build with a Mac CPU golden by hex. The TOST (GP5) uses the mean `iou_tc`, and `final_eval` reads `iou_tc` from the `score_polygons` terms [code: final_eval.py:29; eval/metrics.py:311 `t.iou_tc for t in ot`]. So the doc already guards the E-G10 input. A missing R-D2 (sanitizers, `eval_sets`) does not change the E-G10 result.
- cuda-4: minor (kept as filed). BudgetStop minute arithmetic: 1.95 / 0.01133 = 172.1, not 172.3. No money or gate changes.
- cuda-5: minor (kept as filed). Stale 9,600-row count in 00_front 2.4 and D7b. Wording only.
- cuda-6: minor (kept as filed). D-G13 priced at dph_total, not r_s. The bound shifts by $0.05; no gate result changes.
- cuda-7: minor (kept as filed). "pool size 1 (as A3)" is stale after D7 step 6. Wording only.
- cuda-8: minor (kept as filed). The `contexts` counter wording. The reviewer notes that no gate result changes.

Counts: confirmed 2 (major 2), refuted 1, downgraded 0, minor kept 5.

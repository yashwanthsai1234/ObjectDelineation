# Verdict on final_r5_cuda.md (fresh verifier, 2026-11-08)

## Skills used

| Skill | What it changed |
|---|---|
| nvidia-skills:nvidia-skill-finder | Routing only. No catalog skill was installed. The CUDA facts I checked (pool allocator, grow-only arenas) are read straight from repo code. |
| engineering:code-review | I read each cited code line: policy.py:374-383 (code_files drops `docs/`), bundle.py:48-70 (train views only when named), gpucheck.cpp:131-158 (grow-only `Block`), reward/cache.py:241-257 (tables uploaded per new cache), CMakeLists.txt:15 (CUDA OFF by default), policy.py:399 (CPU install). |
| engineering:debug | I traced each failure to its first break: R-D1 step 1 on a box with no `docs/` tree; `c_speed_check` with no train views in the bundle; D-G15/B04 on a step with new tables or a larger batch; a V4 eval on the CPU build. |
| engineering:system-design | I checked that no other section ships the kernels, names a V2 bundle command, or gives V4 a plan file (grep of all phase sources). |
| engineering:architecture | One-rule check: C1 "train views in the V2 bundle" against A.md's only bundle command and front 6.4 rule 7; D-G15 against B04 and the D1 counter list. |
| engineering:documentation | Each confirmed item below cites the source file and heading plus the code line that proves it. |
| huggingface-skills:huggingface-best | No model claim is in scope. No change. |

Web: WebSearch and WebFetch loaded. No web fact was needed: every item is decided by repo code.

## Items

- cuda-1: CONFIRMED (blocking). D.md D1 step 0 has `r_d1.yaml` run the K1/K2 runbooks from `docs/implementation2.0/kernels/` with "No manual `scp`", but `code_files()` drops `^(data|docs)/` [policy.py:379] (CLEAN1's new regex keeps that drop), and the LIDC bundle holds no `k1_bundle` [bundle.py:52-70]; the runbooks still hold scp/rsync lines (K1 GPU_RUN.md:51, K2:33, K3:16), so V2 item 4 ($0.37) and item 5 cannot run as written.
- cuda-2: CONFIRMED (major). C.md C1 draws from "the train views in the V2 bundle", but `members()` drops `lidc-train-*` unless named [bundle.py:55-56], A.md's only bundle command is `bundle lidc --view lidc-bundle-v3` (rl_shifted.yaml:14 also names the train views; A.md dropped them), and no V2 item names a bundle command or a rows fetch; front 6.4 rule 7 does not list train views.
- cuda-3: CONFIRMED (major). D.6 D-G15 counts "`cudaMalloc`/`cudaFree`" calls that D5c removes (vacuous), while the B04 row and the D1 counter read "device allocs/releases delta == 0" in steps 2-6, which fails on correct code: `Block::get` grows on a larger call [gpucheck.cpp:136-150] (I08 itself says "a larger shape allocates"), and each new image's tables are uploaded from the same pool on first use [reward/cache.py:247-257]; the engineer must guess which counter the gate reads.
- cuda-4: CONFIRMED (major). Front 6.3 "V4" and E.8 name no plan file and no CUDA build step, unlike V1-V3 (front 6.1 step 4d covers only `b1_rows.yaml`/`a3_measure.yaml`); `sync_up` installs the CPU build [policy.py:399; CMakeLists.txt:15 `OBJDEL_CUDA ... OFF`], so the V4 sampler arm scores on the CPU checker and the GP1-GP5 tests skip under the A.9 hook (`build() == "cuda"`), which spends V4 money for nothing and breaks the owner CUDA rule.
- cuda-5: minor (as filed). Wording of the pool-size-1 sentence; behaviour is as intended.
- cuda-6: minor (as filed). 1-2 min [estimate] inside a 2.7-min margin; no proof that the cap is crossed.
- cuda-7: minor (as filed). Timing text disagreement; `b3_v5.yaml` step 0 guards the reward path.
- cuda-8: minor (as filed). Placement-table wording; no V1-V4 plan runs `select_ckpt`/`g_probe` on a box.
- cuda-9: minor (as filed). $2.49 vs $2.51; no budget rule changes.
- cuda-10: minor (as filed). Stale numbers in FR1, D7b, C.8.2; 6.3 holds the one price.
- cuda-11: minor (as filed). D-G7 tar path and Python pin; a careful engineer can choose them without changing the result.
- cuda-12: minor (as filed). Key list for "IDENTICAL"; the comparable keys are implied.

Totals: 4 confirmed (1 blocking, 3 major), 0 refuted, 0 downgraded; 8 minors unchanged.

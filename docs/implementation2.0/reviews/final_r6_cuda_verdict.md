# Verdict on final_r6_cuda.md (fresh verifier, 2026-11-08)

## Skills used

| Skill | What it changed |
|---|---|
| nvidia-skills:nvidia-skill-finder | No NVIDIA catalog skill covers fork safety, rental timing or plan wiring. No verdict depends on it. |
| engineering:code-review | Read `reward/cache.py:241-257` (`device()` per pid), `reward/score.py:108-114` (`score_many` calls `device()`), `train/serve.py:163-177` (`fork_pool` forks with `get_context("fork")`), `vast/policy.py:308-316, 470-490, 527-536`, `final_eval.py:67-77`, `K1/cross_numpy.py:16`, `K1/build_cpu.sh`. |
| engineering:debug | Traced each paid path: V2 item 1 -> C0 arm `--ckpt {rd}/ckpt` (A.md A6 step 8), V3 item 0 -> `b3_v5.yaml` rows root, V4 cap -> eval start, R-D1 step 1 -> `cross_numpy` import after D3. |
| engineering:system-design | Checked process topology of the A2 ceiling CLI against the D7 step 6 pid guard. |
| engineering:architecture | Checked rule agreement: D7 step 6 / CLEAN16 vs A2 step 2 / I19; front time line step 5 vs D1 step 4; front 6 "One source" vs E6 step 8b (V4 RunSpec is its own stamped launch config). |
| engineering:documentation | Each confirmed line names the doc text (source file + heading) and the code line or ledger row. |
| huggingface-skills:huggingface-best | No model choice in scope; no change. |

Web tools loaded (WebSearch, WebFetch); no verdict needed a web fact (all claims are repo-internal).

## Verdicts

- cuda-1 (BLOCKING): CONFIRMED blocking. D.md "D7" step 6 makes `ExpertCache.device()` raise in any forked child and converts only `final_eval`, `select_ckpt`, `g_probe`; A.md "A2" step 2 parallelises the ceiling with `serve.fork_pool(serve.pool_size())` and A.9.3 I19 runs it with pool sizes 1 and 2; the worker path is `score_many` -> `caches[n].device()` [code: reward/score.py:114] in a child made by `get_context("fork")` after the parent imported `SCORERS` [code: serve.py:170-173], so the pid differs and I19 fails once D7 lands; no phase converts the ceiling (grep "ceiling" with pool/fork in D, E, CLEAN: none), and CLEAN16 step 1 then deletes `fork_pool` the ceiling imports.
- cuda-2 (MAJOR): CONFIRMED major. 00_front.md 6.3 V2 item 1 lists "launch, sync, `env_check`, toolchain and CUDA build" with no `fetch_hf.sh`, yet the C0 arm is `train_rl ... --ckpt {rd}/ckpt` (A.md A6, line 847) and C1 `c_speed_check` loads the same `--ckpt`; V3 item 0 fetches `rl-full-3ep` and `b1-rows --as b1` but not `--rows <rows_revision> lidc-train-b1`, while B.md B3.9 uses "the same B1 rows root" and B1.5 says "A plan step on the box cannot download" [code: `cmd_fetch_hf` is the only token path, policy.py:527-536]; V4's sampler-arm `final_eval` also has no fetch in front 6.3 or E.8.
- cuda-3 (MAJOR): CONFIRMED major for V2 and V3; the V4 part is refuted (E.md E6 step 8b: "the file is a `RunSpec` ... `make preflight` stamps it ... front section 0's 'launch config' and 'plan file' name this one file"). Front section 6 claims to be "the only place that lists ... the launch configs that `make preflight` stamps" but names only `b1_launch.yaml`; `cmd_launch` refuses without a stamped `--config` [code: policy.py:315-316]; for V3 the guess changes the outcome, because `push_hf` takes the HF revision from the launch config's `output_dir` [code: policy.py:511-524], so reusing `b1_launch.yaml` would push B3 onto revision `b1-rows`, and a `r1_qwen3vl2b.yaml` guess drops the "ROWS root" stamp check.
- cuda-4 (MAJOR): CONFIRMED major. 00_front.md 6.3 "V4" prices "30 min overhead + 3 + 9 = about 42 min = $0.58", while the front's own inputs table gives launch to first `run` up to 38.2 min [measured: ledger.csv r5 launch 15:35:45Z, run 16:13:57Z]; the cap from $0.58 is floor_to_0.05(0.48/0.680) = 0.70 h = 42 min and `run.sh` ends jobs at the cap [code: policy.py:480], leaving about 3.8 min for build 3 + env_check 1 + tests about 2 + eval 7 on an r5-like host, so the E-G10 eval never runs.
- cuda-5 (MAJOR): CONFIRMED major. 00_front.md section 1 "Time line" step 5 lands "D2-D8 code" before V2 (flow diagram the same), but D.md "D1" step 0/4 writes R-D1 for a pre-D3 tree; D.md "D3" drops `bitquads` from exports and moves `bitquad.py` to `tests/oracles/` without listing `docs/implementation2.0/kernels/K1/cross_numpy.py`, whose line 16 is `from objdel.checker import bitquads, polygonize` [code], so `r_d1.yaml` step 1 raises ImportError and `&&` [code: policy.py:490] stops the K1, K2 and bench steps; the front's own 6.3 line "Phase D, which lands after the B runs" contradicts step 5 too.
- cuda-6 (minor): as filed, minor (wording on the forked pool-size-1 worker; result is correct since `fork_pool` precedes `make_trainer` [code: final_eval.py:70-73]).
- cuda-7 (minor): as filed, minor (BudgetStop rule wording vs the V1 value).

Counts: confirmed 5 (1 blocking, 4 major); refuted 0 whole items (cuda-3's V4 part refuted, the item stands for V2/V3); downgraded 0; minor as filed 2.

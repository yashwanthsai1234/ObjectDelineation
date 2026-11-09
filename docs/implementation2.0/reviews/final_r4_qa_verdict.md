# Verdict on final_r4_qa.md (fresh skeptical verifier, 2026-11-08)

## Skills used

| Skill | What it changed |
|---|---|
| nvidia-skills:nvidia-skill-finder | No catalog lookup needed (no new NVIDIA product claim); it pointed me to check the box toolchain (`run_box.sh` installs nvcc, the image has none). |
| engineering:code-review | Read `CMakeLists.txt:15`, `Makefile:41-44`, `policy.py:399, 484-490`, `bundle.py:26, 53-62`, `preflight.py:28, 59-95`, `train_rl.py:10, 58-60`, `env_check.sh`, the run plans in `experiments/phase2/configs/`. |
| engineering:debug | Traced each claimed failure: which build is on the box after `sync_up`, which exit code a step gets, whether `&&` stops the chain. |
| engineering:system-design | Checked the V1/V2/V3 order against the B3 preset rule (item 1). |
| engineering:architecture | ADR-0013 rule (b) in B vs front 2.5 / 6.3, C header, 99_back decision 16. |
| engineering:documentation | Checked whether each problem is already handled in another part of the doc (none was, except item 7's runnability). |
| huggingface-skills:huggingface-best | No model question; used only for the HF upload content of item 8 (merged `final/` vs adapter). |
| Web | Not needed: every claim was checked against repo code or the doc text. |

## Verdicts (blocking/major items)

- qa-1 | confirmed (blocking) | 00_front 2.5 ("so F is the expected B3 preset"), 00_front 6.3 row 2b ("F fails: items 3 and B3 use the preset that rule 0 falls back to"), C.md:31-37/1059 and 99_back decision 16 contradict B.md:375/1152/1716 and U51 ("the B3 config is `bf16-seqmask` whatever the rule returns"; "no F arm is needed in V2"); no precedence rule exists, so the V2 item 2b spend and the C-arm preset are undecidable.
- qa-2 | confirmed (blocking) | `sync_up` runs `uv pip install ... -e .` (policy.py:399) and CMake defaults `OBJDEL_CUDA` OFF (CMakeLists.txt:15); only `run_box.sh`/`make cuda` builds CUDA, and front 6.2 / B1.5 `b1_rows.yaml` steps have no build step (rl_full.yaml:25 had one), so V1 rewards are scored by the CPU build (owner rule) and A.9's hook (A.md:1238) skips G02/G08 so the smoke passes by construction; env_check.sh has no build check.
- qa-3 | confirmed (blocking) | `b1_rows.yaml` step 2 selects `-k "G02 or G08"` (B.md:709, 1947) but B's naming conventions (B.9.0) and sketches (`test_sums_cpu_gpu_bit_equal`, `test_is_ratio_geo_mask_values`) carry no ids, so the selection can be empty (pytest exit 5, chain stopped by `&&`, policy.py:490), and G08 (B.md:1937) includes G01 `sum_x`/`sum_y` parity from B3.1, which order rule 2 keeps off the `v1-b1` tag.
- qa-4 | confirmed (blocking) | I22 (B.md:1903) requires `--budget-usd 1.90` and an "R5 on val edge" step, while B1.5 (B.md:713-724) and front 6.4 rule 5 set $2.05 and move R5 val edge to `a3_measure.yaml`; the test and the task contradict, and `make test` gates the tag.
- qa-5 | confirmed (major) | B1.3 change 2 (B.md:649-652) checks "the local root" taken from the plan's `--train {data}/rows/lidc-train-b1` (B.md:712), but the local root is `data/views/lidc-train-b1` (B.md:668) and no text maps one to the other; I19 does not pin it.
- qa-6 | confirmed (major) | Front 6.4 rule 7 (00_front.md:782-784) says a view not in the rental's table is moved out, but `bundle.members` dies without every `LIDC_EVAL` view (bundle.py:26, 58-60) and `preflight.hashes` reads `lidc-val-rescaled`, `lidc-val-rotated` and both probes (preflight.py:28, 66), so the literal rule refuses V1 (at $0) and the engineer must guess which views stay.
- qa-7 | downgraded (minor) | G05 (B.md:1934) describes a self-contained test that builds its own engine ("`vllm_kwargs(dtype=...)` yields an engine"), so a separate pytest step can run it before `train_rl`; only the runner wording "inside its first 15 steps" (B.md:1950) is unclear, and no code proves a wrong outcome.
- qa-8 | confirmed (major) | B1.5 uploads only `final/` to `b1-rows` (B.md:733-735), and `train_rl` saves `final/` after `merge_and_unload()` (train_rl.py:58-59), but B2.3b's box part (B.md:1131-1133) reads "the B1 final adapter" through `fetch_hf.sh <run> b1-rows`, so way (c) "LoRA delta merged" has no delta and the merge-rounding reading is vacuous; front 6.2 item 2 repeats "final adapter".

## Minor items (not re-graded; already minor)

- qa-9 .. qa-15 | unchanged minor | Already minor in the review.

Totals: confirmed 7 (4 blocking, 3 major), refuted 0, downgraded 1.

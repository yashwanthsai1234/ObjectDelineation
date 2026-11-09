# Verdict on final_r8_qa.md (round 8, QA lead)

## Skills used

| Skill | What it changed |
|---|---|
| nvidia-skills:nvidia-skill-finder | No new NVIDIA product claim is in scope; it focused me on the box CUDA-build order (which plan builds `_gpucheck` with nvcc before each step). Items 1 and 5. |
| engineering:code-review | Read `policy.py:396-399` (sync installs `-e .`), `:432-435` (`steps:` plans run as written), `:484-490` (`job`, `&&` chain) and `CMakeLists.txt:15` (`OBJDEL_CUDA` OFF by default) to check each build claim. |
| engineering:debug | Traced the V2 and V3 launch lines in order and asked which build is in `/opt/venv` at each step. Confirmed items 1 and 5. |
| engineering:system-design | Matched each box gate (D-G14/15/17, G08b) to a producing plan step or code task. None exists (items 3, 4). |
| engineering:architecture | One owner per rule: A6 and C0 give different argvs for the same C0 baseline arm (item 2). |
| engineering:documentation | Checked that each box test names its inputs (G02b model texts: item 6). |
| huggingface-skills:huggingface-best | No model choice in scope; no change. |
| engineering:debug (second listing) | Same skill as above; no extra effect. |
| Web (WebSearch, WebFetch loaded) | No web fact needed: every item turns on source text and repo code. |

Sources read: phases/00_front.md 6.3 (V2 table, V2 launch lines, V3 table and lines); A.md A3 plan files (lines 479-505), A6 step 8 (lines 852-857), A.9.6 G02b (1885-1898), A2 `item` (377); C.md C0 files and steps (330-395), U-CF4 (1647-1652); D.md D1 counters (215-247), D7 steps 5-6 (685-709), D.6 (935-938), D.8 (983), D.9 B04/B06, S05; B.md B1.5 step 0 (757-762), B2.3 (1335), B3.9 (1670-1678), I22b (2202), B.9.5 G08b (2237, 2248). Repo: experiments/phase2/vast/policy.py, CMakeLists.txt.

## Verdicts (blocking and major items)

- qa-1 | confirmed (blocking) | 00_front 6.3 "V2 launch lines" runs `a3_measure.yaml` (item 1b) before `c_speed.yaml`, whose step 0 is the only V2 CUDA build, while A3's `a3_measure.yaml` step 2 asserts `g.build() == 'cuda'`; `sync_up` installs the default CPU build [code: policy.py:398-399; CMakeLists.txt:15] and the `&&` chain stops [code: policy.py:490], so item 1b always fails and the paid minutes give no A3 val answers.
- qa-2 | confirmed (major) | A "#### A6" step 8 gives the C0 baseline arm the full argv `--train {data}/prepared ... --steps 17` while C "### C0" files require the B1 rows root, both rows flags and `--steps 20`, and C0 step 3 points to "A6's script (A6 step 8)" with only `...`; U-CF4 parses only `-m` steps, and both arms are `bash` nsys-wrapped steps [code: policy.py:486], so the guard cannot catch the stale argv.
- qa-3 | confirmed (major) | D.6 D-G14/D-G15/D-G17 and D.9 B04/B06 require `gc.stats()` deltas, `max_memory_reserved`, `pinned_grows` in the B3 run JSON and a 5 s sampler writing `runs/d/topology.json`, but D7 step 5 says only "In V3, record ...", and a grep of every phase source finds no task that makes `train_rl`/`final_eval` log these or writes `topology.json`; B3.9/I22b list `b3_v5.yaml` with no such step, so the gates cannot pass or fail.
- qa-4 | confirmed (major) | B.9.5 places G08b at "the step after the CUDA build of the first plan after B3.1 merges", but every V3 plan has a fixed step list with no room for it (D S05: `d_gpu_test.yaml` "has exactly 2 steps" naming G02b/G03; A3 `a3_test.yaml` has 2 steps; B3.9/I22b `b3_v5.yaml` has no smoke step; `b2_fp16.yaml` has no build), so an engineer must guess where the pre-B3 `sum_x`/`sum_y` box parity check runs, and front 6.3 item 4 prices no smoke minutes.
- qa-5 | confirmed (major) | B "#### B2.3" lists `b2_fp16.yaml` with no step 0 and I22b asserts no build for it (it does for `b3_v5.yaml`); V3 item 1 is "Skippable" and A3 itself says "no earlier V3 item is sure to build it", so with item 1 skipped the paid 60-step fp16 arm scores rewards on the CPU build [code: policy.py:399; CMakeLists.txt:15], against the owner's "CUDA is a must" rule; kept at major because both items are optional paths.
- qa-6 | confirmed (major) | A.9.6 G02b scores `model` of `ceiling.item` for "50 val items", which needs K answer texts (A2: "With `texts=None` it returns the doctor side only"), but no part names which 50 items or which texts, and D8 says the saved val answers are "Not in the payload", so the natural input (A3's `answers.jsonl`) is missing on the V3 box and the box half fails; a missing input that changes the outcome.

Minor items 7-10: not re-rated (already minor; no action required in this verdict).

Totals: 6 confirmed (1 blocking, 5 major), 0 refuted, 0 downgraded.

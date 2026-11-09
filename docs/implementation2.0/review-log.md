# Review log of the Implementation plan

Each phase keeps the record of what each review found and what changed. The plan is `Implementation.md`.


# Front matter

## Review log (final round 1)

Sources: `reviews/final_r1_{arch,cuda,data,fresh-arch,qa,rl}.md`, items addressed to 00_front. Facts were checked
against code and files on 2026-11-08 before each change.

| Review item | Sev | What changed | Check made |
|---|---|---|---|
| data 3 | blocking | Section 0: miss = K-mean of IoU_tc below 0.1; new word "missed answer"; found image = not a miss | `final_eval.py:26-30` stores K-means; test 10.8 % (0.10818), val 12.6 % recomputed from `items.jsonl` |
| cuda 1, qa 1, fresh-arch 1, cuda 11 | blocking / major | G1 removed. V1 = B1 run first (data first), then A3 in the same session. V2 = C checks and arms, R-D1, K3 only if C1 fails both backends. R-D2 moves to V3 after D2-D8 merge. All minutes and caps recomputed | D.8 rental table (R-D1 vs R-D2); B.6.2 pairs Gate B1 with R5 `items.jsonl`; B.R item 11 allows B1 on bf16-seqmask |
| cuda 2, qa 2, rl 15, fresh-arch 9 | blocking / major | Gate for C adoption = C-G-EQ1..6 + C-G-SP1 (<= 20.0 s A / <= 20.7 s B). 18 s (C-G-C) measured inside B. Diagram and section 1 changed | C.md lines 707, 715, 819 |
| arch 1, qa 3 | blocking | One numbered time line in section 1; order rule "every A-D task path is pre-E" | C.8.3, D.7 R10 agree; E.7/E.8 must follow (E fixer) |
| arch 2 | blocking | NVTX marks move from `core/marks.py` to `policy/marks.py` in 2.2, 2.3 and section 4 (qa 20 asked for `core`; arch 2 wins because `core` falls under the no-torch contract) | A6 `nvtx.py` imports torch inside `enabled()` (A.md) |
| arch 3 | blocking | Seam = D3's `geom.rings.polygonize`; tree lists `rings.py volume.py certify.py`, no `polygonize.py`; numpy code in `tests/oracles/` | D3 file list (D.md:283-312) |
| arch 4, fresh-arch 7, qa 7 | blocking / major | One ADR table 0012-0021 in section 4 (arch 4's order, time-ordered) | last ADR on disk is 0011 |
| cuda 7 | major | A3 no longer last in a tight cap; it rides V1 with its own skip rule; Mac half of C-G-EQ2/3 after V1 is stated | A.8 minutes |
| qa 1 (nsys) | blocking part | Corrected in round 2: A6 capture = C0 baseline arm in V2 (option a). The round-1 text ("A3 keeps its own 6-step window") was wrong | Round 2 grep: A.md has no `--profile-window 4 5`; A3 has no nsys step |
| cuda 5, data 13, fresh-arch 6, arch 9 | major | View bank removed from 2.2, 2.3, section 4; deferred with a trigger and a fork-safe design note; 2.1 reworded | grep: no task builds it; C calls it optional |
| cuda 12, arch 9, qa 25, rl 23, fresh-arch 10 | major / minor | 2.4 shows vLLM's CPU image lane and "C1 patchify (HF/torch; K3a fallback)"; 2.5 drops the sampler column and K3b, shows 256 requests with n=8 as the C6 arm, and no old-logp pass | C.3.3 (201 MB), C6 default off, `fast.py:540-547` |
| arch 9 | major | 2.2/2.3 regenerated from E.3.1 plus fixes; sampler interface `sample(rows, params) -> list[Sample]`; placement table | E.md:139-178, 665 |
| arch 6, 7, 8; fresh-arch 8; qa 20 | major / minor | Placement table: coverage in `data`, advantages in `verify`, two stores with purposes, homes for A/C/D/CLEAN new modules | B1.3, E.3.1 |
| arch 10, fresh-arch 5, 6 | major | Section 4 has an "Owning task" column; deferred rows named; EDT row = K2 on device (cuCIM not adopted); Vast row = `vastai` CLI 1.11.0, no SDK; xdist deferred; `data_revision` replaced by `rows_revision` + manifest SHA | D.2 line 96; CLEAN.md:139; D.md:551 |
| fresh-arch 6 | major | `linear_attention` guard assigned to E5 instead of "guard now" with no task | grep: no task had it |
| cuda 13, arch 14, qa 5 | major | Section 0 path rule; section 5 full repo paths and runbook path rule; FR13 widened | GPU_RUN.md files still contain `/[p]rivate/tmp` (grep) |
| cuda 14, qa 15, fresh-arch 9 | major | One planning-input table: $0.641/h, measured $0.544-0.628/h, 27.5 s base, 23.65 s, 20.0 s only after C-G-SP1 | ledger.csv r3/r4/r5 launch dph 0.5437/0.5926/0.6281; r3 billed $0.440; gpu-run.md:136-137 |
| cuda 15, qa 15 | minor | One stop rule: 35 min to first sync | ledger: launch-to-sync 16.5 (rl), 18 (r3), 20 (r4), 31 (r5) min |
| cuda 16 | minor | `--cap-hours` computed from the picked offer's dph | - |
| cuda 17, qa 27 | minor | Toolchain `cuda-toolkit[nvcc,cccl,nvtx]==13.0.3`, noted as after D1 | `run_box.sh:24` has `[nvcc,cccl]` today |
| data 4, rl 3 | major / blocking | 4B gate on val, not test; suggested val numbers; test read once | R5 val 0.618 / 12.6 %, test 0.645 / 10.8 % |
| data 11, rl 20 | major / minor | Reason 1 marked "[to be confirmed by A2]"; free check reworded: bounds found images only | - |
| rl 14 | major | Size-limit paragraph: all 4B candidates exceed 4.0 B; owner decides | HF API: MedMO-4B 4,437,815,808 |
| rl 18, fresh-arch 11 | minor | EGM: Tables 4 and 9; Instruct 84.2/87.3, Thinking 83.6/87.2, after 89.6/91.0; gap 3.6 -> 1.4. HuatuoGPT-Vision cited to its size (M3), not to 2604.27720 | arXiv HTML of 2601.13633; abstract of 2604.27720 names MedGemma and Lingshu only |
| rl 19 | minor | "Vision tower identical" -> same blocks, merger 2,048 vs 2,560 | HF `config.json` of both models |
| qa 18, rl 22 | minor | One medical model id, `MBZUAI/MedMO-4B` (99_back decision 13 must match) | HF API |
| qa 6 | major | Id prefix rule in section 0; rentals V1-V3; risks FR1-FR14; B data gates B-DG1..16 | - |
| qa 12 | major | Gate A box moved before V2/V3 (after V1, on the Mac) | - |
| qa 14 | major | Section 0: image vs view defined apart; code-name rows (window, rater/voter, completion); "test", never "exam" | - |
| qa 8 | major | Tests row: one marker `gpu` (+ `sanitizer`), registered in A3/A6, `--strict-markers` | - |
| qa 17 | major | Front sentences split to one action each in sections 3 and 6 | - |
| qa 19 | minor | Section 6.1 step 3 commands fixed: `build_cpu.sh "$REPO/cuda"`, `run_cpu.sh`, `make test-cpu` | K1 `build_cpu.sh` needs `$1`; K2 has `run_cpu.sh`; K3 Makefile has `test-cpu` |
| qa 28 | minor | Tags `[web:]`, `[unverified]` added; dated `[measured]` form | - |
| arch 22 | minor | Tree: K1, K2 merged into `gpucheck.cpp`; K3 only if adopted, as `cuda/k3/` torch extension | `pyproject.toml` builds one extension today (per arch 22) |
| arch 11 | blocking (E) | Tree and 2.1 keep `experiments/phase2/` (CLAUDE.md rule) | CLAUDE.md |
| fresh-arch 14 | minor | Skills row and section 4 name `nv-segment-ct`, `nv-segment-ct-finetune`, `dicom-series-to-volume` (recommendation only) | live catalog JSON contains all three |
| data 25 | minor | Data versioning row: one place (model repo `rows/`), pinned `rows_revision` + manifest SHA, LIDC only | B1.4 step 5 (B.md:350) |
| cuda 24, qa 4 | minor | Section 5 notes the kernel `build/` folders stay out of git (DOC1 owns the `.gitignore` line) | - |

Not handled here (owned by another section): data 25's `train_rl` SHA assert (B), C-G-SP4 wording (C), B.8.1's 17 s
column and its session split (B), 99_back decisions 1, 3, 11, 12 and 13 (back matter must follow sections 3 and 6),
E move tables (E), copying the notes into `docs/implementation2.0/` (DOC1).

## Review log (final round 2)

Sources: `reviews/final_r2_{arch,cuda,data,fresh-arch,qa,rl}.md`, items addressed to 00_front. Facts were checked
against code and files on 2026-11-08 before each change.

| Review item | Sev | What changed | Check made |
|---|---|---|---|
| fresh-arch 1 | blocking | New FR15; section 6 input row "Code budget guard"; 6.1 step 0 (`HARD_STOP` = $16.90 + $3.00 = $19.90, owner approval, unit test), step 4b (`spend()` + cap <= `HARD_STOP`, `--owner-approved`); word "session spend" in section 0; every "$" stop rule in 6.2 uses session spend | `policy.py:58` `ALARM, HARD_STOP = 10.0, 17.5`; `policy.py:324-329` launch refusals; `policy.py:588-592` sweep destroys all at `HARD_STOP`; sweep plist `StartInterval` 3600; ledger sum 16.9025 over 133 rows; `PULL_S = 600` (`policy.py:73`) |
| fresh-arch 2, data 1, qa 2, rl 10 | blocking / major | V1 table rebuilt: B1 run; B1 model on val shifted and val edge (not skippable); R5 on val edge (not skippable); A3-val (skippable). Totals 195.1 min = $2.08 at 27.5 s, 175.9 min = $1.88 at 23.65 s; cap stays $2.20 (owner may raise to $2.30). Noise floor = A3 step 7 halves. B2.3b and A3's test half move to V3. 6.1 step 4c: the bundle holds `lidc-val-edge`. "A3 is not needed for Gate B1" replaced by the four-input list | val edge view = whole val (B1.6 `--split val --kinds ...,edge`); 415 s includes the model load (`seconds.load_and_sample` 414.7 s); `a3_measure.yaml` has val and test shifted only (A.md); B2.3b is 15 min [estimate] (B.md) |
| cuda 2, qa 1 | major / blocking | V1 has no nsys window. A6 capture = C0 baseline arm in V2 item 3 (steps 16-17); A-G6 box half open until then. Round-1 log row corrected | A.md line 395 "The plan has no nsys step"; A-G6 reads the C0 baseline arm; no `--profile-window 4 5` in A.md (grep) |
| cuda 3 | major | Section 6 input "start-up per new training process, about 6 min"; V2 item 2 = 12 min, item 3 = 37 min with the check arm continued in the fast arm's process; items 1-3 = 89 min = $0.95; "fits only near 23.65 s" replaced by "does not fit at either step time" | ledger rows 125-126 (r4 31.3 min run to teardown) and 133-134 (r5 57.3 min); r5 `timing_s.train` 465.8 s and 354.8 s; r4 asleep 451.1 s |
| fresh-arch 8, qa 16 | minor / major | Order rule 2: owner tag `v1-b1` = be52b1c + A3, A6, B1 (B1.5) + CLEAN wave 1, with the tests that must pass before the tag. 6.1 step 4 names `b1_rows.yaml` (launch) and `a3_measure.yaml`. "today's code" removed | `BudgetStop` absent from `src/` today (review grep); `tests/test_p2_train_rl.py` exists; B E08 is the dry run (B.md) |
| fresh-arch 6 | major | Section 4 row "RL trainer": keep TRL + FastGRPO behind the adapter; verl considered; trigger named | Web 2026-11-08: verl v0.7.1 release notes list LoRA adapter-only refit and Qwen3-VL fixes; `merge=False` semantics with vLLM rollout not confirmed, so tagged as not verified |
| fresh-arch 5 (front part) | major | ADR row 0022 for C10; note that B and C shift 0012-0014 up by one; DOC4 drift check named | C.md C10 asks for "the number from the front's ADR allocation table" (grep) |
| arch 10, qa 14 | major | V4 (optional, after E, owner decision) in section 1 diagram, time line step 7 and 6.3: about 39 min = $0.42 standalone, because launch time is real; E4 keeps the trainer path until V4 passes | E.8 lines on "piggyback" (about 9 min) and "standalone 15 min"; launch to first sync 17-31 min [measured] |
| qa 5 | major | "B-DG1..16" -> "B-DG1..B-DG18 (21 gates incl. DG1b, DG9b, DG14b)" in sections 0, 1, 2.4, 6.1, FR11; prefix rule: inside a phase the prefix may drop, cross-phase references keep it | B.md gate table: DG1-DG18 plus DG1b, DG9b, DG14b = 21 rows (grep) |
| qa 8 | major | Tests row lists every marker with its owner (`slow`, `data`, `gpu`, `sanitizer`, `docs_build`, `post_commit`, `box` if kept) | `pyproject.toml:68-71` has `slow`, `data`; DOC.md registers `docs_build`, `post_commit`; A registers `gpu`, `sanitizer` |
| qa 13 | major | 6.3: D1's Mac parts gate D2; R-D1 numbers gate only D-G7; D-G14/15/17 read V3 | D.md (per review); front rental order |
| qa 15 | major | Time line step 1: DOC1 and DOC2 run before any phase (S06 needs DOC2's "TRL 1.16" edit) | DOC.md lines 58-61, 297 (the 1.16 edits) |
| cuda 10 | major | Section 5: the runbook path edit is D1 step 0 (before V2), with an R-D1 `run.sh` config and a preflight check for `/[p]rivate/tmp`; V2 item 4 names it | grep: K1 GPU_RUN.md:15-18, K2:33, K3:15 still name the scratchpad |
| arch 9 | major | Section 3 step 1: `policy/models.py` after E, `configs.py:73-77` before E; CLEAN8 keeps `_LORA`, `_GRAD_CKPT` | `configs.py:76-77` define `_LORA`, `_GRAD_CKPT` |
| rl 11 | major | Section 3: step 0 SFT for the 4B (sft_config, LoRA r 16 / alpha 32, 50-step timing probe); the probe starts from that checkpoint; gate text shows SFT + RL cost | `configs.py:90-95` (`LoraConfig(r=16, lora_alpha=32)`), `configs.py:148` (`sft_config`) |
| data 14 | minor | Section 0 "test" and section 3: test read once after the owner closes the 4B gate; V3 item 5 holds the test eval only then | - |
| data 7 (front part) | major | Time line step 2: B1.7 runs before A3 step 1 hashes any eval view | - |
| data 12 (front part) | minor | Section 1 adds B's union-box numbers: median 0.262, 36.7 % (807 of 2,200) | recomputed from `geom.json` (`off_euc`, test shifted): 0.2617, 807 / 2,200 |
| qa 17 | minor | Section 0 "test": A3's R5 re-score on test is an equivalence check, not a selection | - |
| qa 18 | minor | 6.1 step 3 sets `REPO=$(git rev-parse --show-toplevel)` and runs `test_k1_cpu`, `_o1`, `_o2` with `3000 41` | `K1/build_cpu.sh` only compiles and prints "build_cpu: OK"; `K1/DESIGN.md` section 7 names the test runs |
| qa 26 | minor | Order rule 1: no exception, B2/B3 land on pre-E paths; E.8 must drop its "otherwise" clause | - |
| rl 17 | minor | "B2 probe" removed from section 1, 6.3 and FR4; FR4 = B2.0 table + B1-run counterfactual stats + B2.4 preset; V3 lists the optional fp16 arm | B.md B.R2 item 4 removed the 4-arm probe; B2.0 exists (grep) |
| arch 14 (front part) | minor | `check()` (A7, raises `TRLContractError`) in the seam table and FR5 | A.md A7 steps 2-4 |
| arch 19 | minor | Tree adds `experiments/phase1/`, `experiments/phase2/vast/`, `docs/plan/`, `docs/results/headline.json`, `mkdocs.yml`, `CITATION.cff` | DOC.md names each (grep counts 12-22); `experiments/phase1` exists |
| fresh-arch 9 | minor | Section 2.5: "no separate old-logp pass (`FastGRPO.skip_old_logps`, existing)"; "push merged weights W + BA to vLLM (all tensors, one call, C5)" | `fast.py:540` `skip_old_logps`; `fast.py:696-717` pushes `p.data + d.to(p.dtype)` for every parameter |
| qa 24 | minor | No front change: section 0 already says 23.65 s; G06 and B.8.1 are A's and B's | `runs/r5/shared_rl.json` `per_step` 23.653 |

Not handled here (owned by another section): A's "G1a" text, A3 header config and `a3_measure.yaml` reduction to val
(A); B.6.2 noise floor, B.8.1 B-S1/B-S2 split, B.8.4 A3 row, B-S2 test removal, ADR renumbering in B and C (B, C);
C.8.2 start-up row, check-arm merge and rate (C); D1 step 0, D.4, D.8 (D); E.8 "otherwise" clause, E4 trainer path
until V4, `box` marker (E); DOC U19 constant and DOC4 ADR drift check (DOC); CLEAN8 post-E names (CLEAN); 99_back
decisions 1, 2, 4 and the V4 row (back matter must follow sections 3 and 6).

## Review log (final round 3)

Sources: `reviews/final_r3_{arch,cuda,data,fresh-arch,qa,rl}.md`, items addressed to 00_front. Facts were checked
against code and files on 2026-11-08 before each change. Where two reviews proposed different fixes, one coherent
choice was made and the reason is in the row.

| Review item | Sev | What changed | Check made |
|---|---|---|---|
| arch 1, qa 1 | blocking | Words "launch config" and "plan file" in section 0. V1 launch config = new `b1_launch.yaml` (RunArgs: `model`, `view: views/lidc-train-b1`, `output_dir: b1-rows`, no `branches`); 6.1 step 4 stamps only it, with a B.9 test (`load(b1_launch)` passes, `load(b1_rows)` raises `TypeError`); order rule 2 and 6.2 follow. arch 1's form chosen over qa 1 (a)'s `rows_plan` key, because `RunArgs(**raw)` would reject a new key too | `preflight.py:49-56` (`RunArgs(**raw)`); `train_sft.py:77-101` (no `steps` field); `cmd_launch` calls `preflight.refusal` (`policy.py:315`); `rl_shifted.yaml` header launches with `r1_qwen3vl2b.yaml` and runs the plan with `run.sh` |
| arch 2 (front part) | major | 6.2 launch lines: `fetch_hf.sh <run> --rows <rows_revision> lidc-train-b1`; the code-change list names `--rows`, `train_rl --rows-revision/--rows-manifest-sha256` (B must follow) | `cmd_fetch_hf` pulls one revision into `{rd}/ckpt` and streams the token on ssh stdin (`policy.py:527-536`) |
| data 1 | blocking | 6.2 is the one V1 list: `a3_measure.yaml` holds one R5 process (val edge = item 4, val shifted = item 5); its test step moves to `a3_test.yaml` (V3 item 3). New 6.1 step 4d: the V1 plan files hold exactly the 6.2 views and no `--split test` (hygiene grep test) | Front 6.2 vs the A3 plan described in the review; arithmetic re-done |
| data 2, qa 3, qa 5 | major / blocking | One primary baseline: item 5 (A3's R5 val items, same frozen view, same box); R5's `items.jsonl` is a second line and the fallback only if item 5 has no output, with a "hash only" view-identity rule. Item 5 is no longer skippable (see fresh-arch 5), so the qa 5 deadlock is gone; order rule 3 and V2 item 1b cover a failed R5 process. qa 3's "R5 file primary" not taken: data 2, B.6.2 and A.4 already use A3's items, and item 5 now always runs | 6.2 arithmetic: worst case $1.90 + 20.8 min = $2.12 <= $2.20 |
| fresh-arch 5 | major | `final_eval --view` repeats (A3): one process per checkpoint. V1 items 3 (10.4 min), 4 (6.9) + 5 (3.5) in one R5 process. Totals 188.3 min = $2.01 (27.5 s), 169.1 min = $1.81 (23.65 s). New section 6 input row (205 s fixed + 0.193 s per item) | `final_eval.py:49` (one `--view`), `:67-73` (preload, pool, trainer per call); 414.7 s / 630.5 s from `final_eval.json`; fit 0.1934 s, 205.1 s |
| fresh-arch 2, qa 6 | major | `BudgetStop` clock = launch time (`train_rl --session-start`, from the active-run `started`; `run.sh` passes it and `dph`). 6.4 rules 2 and 5: the 90 % rule is applied by the main agent between `run.sh` calls; an item that fits is not stopped | `policy.py:323` (`started`, `dph`, `deadline`); `policy.py:424` ("~" only lets a step fail); `cmd_launch` deadline = cap + `PULL_S` 600 s (`policy.py:73, 321`) |
| qa 7 | major | B1 final adapter leaves the box by `push_hf.sh <run> b1_launch.yaml` from the Mac (revision `b1-rows`; `hf_push` gains a `train_rl` mode); V3 step 0 fetches it into its own folder. Step-100 adapter dropped: B2.3b reads only the B1 final adapter. A plan step on the box was not used: the token reaches the box only on ssh stdin | `hf_push.py:1-13` pushes `branches/<selected>` only; `cmd_push_hf` revision = config `output_dir` (`policy.py:511-524`); B2.3b box part names "the B1 final adapter" |
| qa 4 | blocking | Order rule 2: `v1-b1` = be52b1c + A1-A8 + B1.0-B1.7, B1.9 + B2.0 + B2.1 (logging only) + 6.1 step 0 + CLEAN wave 1; tests add B2.1's P13 and R01. "B1 waits for B2.1 logging only" | B.md: B2.1 owns `cf/*` keys (R01), P13 bit-equality; B1.8 is "Clean after gate B1" |
| fresh-arch 1 (front part) | blocking | 6.2: B-S1 copies this table with these item numbers and no nsys window; a DOC4 drift check compares names and minutes | B.8.1 still lists the nsys window and $2.35 (B's edit) |
| fresh-arch 3, rl 11 | major / minor | Section 2.5 draws the push by preset (F expected: adapter in place via `add_lora`, base W in mode B; else merged W + BA). ADR 0022 title: "merged or adapter-only, as ADR 0013 decides". V2 item 2b = F arm before the C arms; C arms run under its preset. V2 required items 102 min = $1.09; all items 150 min = $1.60. Section 4 RL-trainer trigger: "the F arm fails twice on a TRL or vLLM bump" | B.md B2.0 table 0.531 at step 100; F arm 13 min = $0.14 (C.8.2); 10-30 % decode cost (C5); `configs.py:193` sleep mode on; TRL `sleep(level=2)` (`vllm_generation.py:377, 772`) |
| fresh-arch 4 | major | New section 1.1 decision table (pass; far-miss only; IoU or DiD fail; guard or legal fail), approved before V1 (6.1 step 6); FR2 points to it; 4B gate timing stated | Gate B1 hard rules in B.md U49 |
| cuda 3 (front part) | major | Section 1: C-G-SP1 relative to the same-rental baseline (- 3.4 s A, - 2.7 s B = 23.43 - 20.0 and 23.43 - 20.7). Section 6 row: 20.0 s is a planning value only | C.md G-SP1 row: pass = 23.43 s minus the midpoint |
| cuda 5 (front part) | major | 6.3: `c_speed_check` reads 5 Mac-recorded batches (`runs/a3/batches/*.pt`, with advantages), shipped in the V2 sync, with a Mac load test | C.md G-EQ1/G-EQ4 use recorded batches; box batches come only from item 3's check arm |
| rl 9 | major | Section 3 step 0: SFT the 4B with `r1_qwen3vl2b.yaml`, only `model` and `output_dir` changed; e2 branch only (`branches: [2.0]`); priced from the 2B SFT ledger rows x the 50-step ratio | `r1_qwen3vl2b.yaml` (view `lidc-bundle-v3`, seed 1, effective batch 32, branches 0.5-3); `full_rl_2026-11-08.md:3` (RL started from e2) |
| data 12 | minor | 6.4 rule 7 names the LIDC eval views of each rental's table; 6.1 step 4c: no `lidc-test-*` folder in the V1 bundle | `bundle.py:53-62` ships every `eval_views/lidc-*`; `LIDC_EVAL` has no test view (`bundle.py:26`) |
| arch 12 | minor | `views` dropped from 2.2 and 2.3; tree notes that the data builders and `build_views.py` / `build_eval_views.py` are not thin yet; the `experiments/phase2/` line names the CLIs that become thin | grep: no task creates `data/views.py`; both build scripts exist in `experiments/phase2/` |
| arch 16, qa 16 (V4) | minor | 6.3: 39 min = $0.42 is the one V4 price; E.8, E decision 7 and 99_back 3b must use it | ledger launch-to-sync 17-31 min |
| qa 16 (B2.3b) | minor | B2.3b box part = about 10 min, B's own estimate (the task owner), in 6.2 notes and V3 step 2 | B.md B2.3b "about 10 min [estimate]" |
| fresh-arch 7 | minor | Section 4 "3D" rows point to ADR 0020's "3D extension" paragraph | E.md skills row: ADR 0020 gains that paragraph |
| qa 13 | minor | Section 0 "outline length": 1,181 without the end token, 1,182 with it; one limit form `L - soft_punish_cache >= 1,182 + 32` | B.md header and DG12; A.md A8 (1,182 with the end token) |
| qa 19 (front part) | minor | Section 8: "Review logs of the final rounds may follow part 9" | - |
| qa 20 | minor | 6.1 step 0: `policy.launch_refusal(offer, cap_hours)` holds `cmd_launch`'s expression; the unit test calls it; step 4b uses it | `policy.py:321-326` |
| qa 21 (front part) | minor | New front text uses one rule per sentence and "Why:" lines; older table cells were not all split | - |

Not handled here (owned by another section): B1.5 code (launch config, rows flags, `fetch_hf --rows`, `hf_push`
`train_rl` mode, session-clock `BudgetStop`, U45/U46), B.6.2 baseline and fallback text, B.8.1 copy of the V1 table,
U49 fallback case (B); A3 multi-view `final_eval`, `a3_measure.yaml` reduction and `a3_test.yaml` (A); C.6 relative
G-SP1, C1/C3 recorded batches, C5 push under F, V2 order with item 2b (C); E.8 V4 price (E); DOC4 drift check of B-S1
(DOC); 99_back decisions 1 (102 min = $1.09), 2 (no $2.30 option), 3b, 13, 17 and the 4B SFT recipe (back matter);
rl 8 (what B3 credits under F) and qa 16's B-S2 price ($1.98 vs $1.94), which no front number fixes.

## Review log (consistency pass, round 4 prep)

Sources: `reviews/final_r3_{arch,cuda,data,fresh-arch,qa,rl}.md` (V1 lists, preflight stamps, `BudgetStop`, base step
time, ledger-based skip rule) and the owner's outline-limit question. Every fact was checked in the repo on
2026-11-08 before the change.

| Item | What changed | Check made |
|---|---|---|
| One source for rentals (r3 data 1, qa 2, fresh-arch 1) | Section 6 opens with "One source": it is the only place for rental items, order, minutes, dollars, caps, skip rules and stamped launch configs. A, B, C, D, E and 99_back now point here. V3 has its own table (items 0-5, B-S2 = items 0 + 4). Owner-decision totals moved into 6.3 | grep of the assembled plan: no "$" table outside section 6 (the new DOC4 check 16 rule, dry run) |
| Session spend left out bandwidth (new, found in this pass) | New inputs rows: bandwidth (b = $0.10 once, $0.039/h) and session rate r_s = $0.680/h. All V1, V2, V3, V4 dollars recomputed: V1 $2.23 / $2.01; V2 required 103 min = $1.26; B-S2 $2.11 / $1.90; V4 $0.54. Offer limit added | `policy.accrued` adds 15 GB down + 6 GB/h up at the offer's $/GB (`policy.py:61, 190-193`); r5: $1.0709 = $0.943 + $0.128 (ledger rows 131-134); highest launch $/GB 0.0065 (row 128) |
| V1 cap and `BudgetStop` (r3 fresh-arch 2, qa 6) | V1 amount $2.40 (was $2.20), cap `floor_to_0.05((A - b) / r_s)` = 3.35 h; `BudgetStop` = $2.05 on session spend incl. bandwidth (run.sh exports started, rate, fixed); 6.4 rule 5 gives the general formula (amount - $0.35 - later items). Worst case $2.49, guard check $19.39 <= $19.90, + reserve $2.74 <= $3 | `run.sh` cap = min(cap, deadline - now - 600 s) (`policy.py:480`); deadline = cap + `PULL_S` (`policy.py:321`) |
| Preflight stamp (r3 arch 1, qa 1) | The launch config keeps `view: views/lidc-bundle-v3` and adds B1.3's `rows_plan` field. The old front text put the rows root in `view`, which `preflight.hashes` cannot read | `preflight.hashes` reads `view_hash` from `<view>/manifest.json`, `view_integrity`/`patients` read `<view>/rows.jsonl` (`preflight.py:59-70, 101-139`); B1.4 writes `lidc/train.jsonl` and no `view_hash` |
| Adapter push | `push_hf.sh` runs `hf_push` on the box before teardown (its minutes in item 1) | `cmd_push_hf` calls `remote(...)` with the token on stdin (`policy.py:511-524`) |
| Ledger-based skip rule (r3 fresh-arch 1, qa 17) | No V1 item is skippable; 6.4 rule 8: no rule reads the ledger total | ledger total $16.90 (133 rows), re-summed |
| Base step time | 27.5 s stays the base; the stop trips only above 28.5 s per step (was 29.6 s) | (172.3 - 30) x 60 / 300 |
| `HARD_STOP` | 6.1 step 0 says it stays 17.5 in code until the owner approves, and lists the one approval message: $19.90, V1 $2.40 with `--budget-usd 2.05`, the item list, decision table 1.1 | `policy.py:58` |
| V2 inputs (C.8.2 vs front) | Item 2 = 13 min: `c_speed_check` makes its own batches (C1); optional rows 3c, 3d added; item 1b = the whole R5 process again (10.4 min); order rule 3: Mac half of C-G-EQ2/3 before V3, not before V2 (C.6) | C.md C1, C.6, C.8.2 |
| Outline limit (owner question) | Section 1: a higher limit lengthens the wait (5.2 s at 1,800 vs 4.0 s at 1,536 per runaway step); limit from the count, not raised | R5 steps 3-15 re-fit: 5.606 s + 4.328 ms, r = 0.918, 5 of 13 steps at 1,536 (`runs/r5/shared_rl.json`) |
| 1,182 claim | Section 0 tag fixed: `[measured: docs/implementation2.0/notes/A1.md F9 ...]` (was "B header; A8", not a file) | Re-ran `A1_tokcount_probe.py` on `data/prepared` with the Qwen3-VL-2B tokenizer: LIDC train max 1,182, 0 above 1,200; RIGA train max 2,246, 1 above 1,536 |
| Data first | Section 1 "Priority (owner rule)": dataset quality first; B1 and its coverage gates B-DG7..B-DG10 before any paid RL run; 6.1 step 2 names them | B.6.1 gate table |
| Section 8 | Review logs live in `docs/implementation2.0/review-log.md` (assembler); part table drops the review-log row; part 8 row says "pointer to section 6" | `assemble.py` `LOG_RE` |

What the owner must approve before V1 (front 6.1 step 0): `HARD_STOP` 17.5 -> 19.90; V1 $2.40 with `--budget-usd
2.05`; the V1 item list; decision table 1.1.

## Review log (final round 4)

Cross-file consistency fixer, round 4. Sources: `reviews/final_r4_{cuda,rl,qa,arch,fresh-arch}.md` and their verdicts.

| Item | Severity | Change | Evidence |
|---|---|---|---|
| cuda-2, qa-2 | blocking | 6.2: `b1_rows.yaml` step 0 = CUDA build (`run_box.sh --build-only`, `build() == 'cuda'` assert, `tests/test_gpu_scoring.py`); `a3_measure.yaml` asserts the build first; "Why step 0" note; `final_eval` pool size 1 on a CUDA build. 6.3: V2 item 1 and `c_speed.yaml` step 0 build; V3 item 3 (`a3_test.yaml`) and item 4 (`b3_v5.yaml`) build first. 6.1 step 4d checks it | `policy.py:399, 432-433, 450-452`; `CMakeLists.txt:15`; `rl_full.yaml:25`; `runs/k1/checker_gpu.log:12-22` |
| cuda-3 | major | 6.2 item 1 = everything before step 1, 53 min (launch to first `run` up to 38.2 min measured + build 3 + `env_check` 1 + smoke 5 + start-up 6). B1 run 260 steps: 53 + 119.2 = 172.2 <= 172.3. Totals $2.29 / $2.10; money left $0.46-0.65; decision-1 totals and the $2.20 fallback (220 steps) recomputed. New input row "Launch to first `run`". Rows root 8,320 rows (section 0). V3 item 0 = 39 min, item 4 = 260 steps + build + start-up; B-S2 $2.11 / $1.92 | `ledger.csv` rows 117-133 (17.7 / 24.1 / 38.2 min); `runs/r5/r5_shared.log:39` |
| rl-1, qa-1 | blocking | Rule (b) everywhere: section 1 diagram has no F arm; 2.5 says F is the expected preset of the runs after Gate B3, B1/B3 run `bf16-seqmask`; 6.3 "One preset for the C arms" replaces "Numerics before speed"; V2 item 2b deleted; item 3 under `bf16-seqmask`; required items 90 min = $1.12, items 1-5 138 min = $1.66; V3 item 2b fp16 re-runs in the first run that uses it; FR4 and section 4 say "F check" | B.md B2.4 "Scope", B.8.4, U51 |
| rl-2, fresh-arch-2 | blocking | Section 3 step 0 is the only 4B SFT recipe: `sft_qwen3vl4b.yaml` with five changed lines incl. `stages: [train]` and the 4B `vision.reference`; the reference is written first in the same rental; chain and preflight tests; 99_back decision 17 points here | `policy.py:441-450, 456-475`; `train_sft.py:160-166`; `r1_qwen3vl2b.yaml:9`; `r0.yaml` |
| rl-5 | major | V3 item 2a runs only when the step-100 error is < 0.2; the F comparator is measured inside the F run (B2.2 F step 4) | B.md B2.2 |
| rl-6 | major | Section 1: C-G-C is measured in the first run after Gate B3; B3 keeps 1,536 and no runaway cap (B3.9) | `configs.py:73`; `grpo_trainer.py:3286` |
| qa-3, fresh-arch-1 | blocking | 6.2 smoke step uses node ids `test_g02_is_summary_no_sync`, `test_g08a_smoke` with `-m, pytest`; G08a has no B3.1 fields; 6.1 step 4d collect-only check (2 items) | `policy.py:484-490`; B.9.5 |
| qa-6 | major | 6.4 rule 7 allows the views of `bundle.LIDC_EVAL` and `preflight.EXAMS`; only `lidc-test-*` is moved out. 6.1 step 4c checks the bundle list | `bundle.py:26, 58-60`; `preflight.py:28` |
| qa-8 | major | 6.2: `push_hf.sh` uploads `final/` and `final_adapter/`; item 2 names both | `train_rl.py:58-59` |

## Review log (final round 5)

Cross-file consistency fixer, round 5. Sources: `reviews/final_r5_{cuda,rl,data,arch,qa,fresh-arch}.md` and their
`*_verdict.md` files. Only the CONFIRMED items assigned to this fixer are listed.

| Item | Severity | Change | Evidence |
|---|---|---|---|
| fresh-arch-1 | blocking | 6.3 heading now "V2, V3 and V4: only after a budget raise". Money-left paragraph: no V2 subset with a gate result fits $0.46-0.65 (item 1 $0.58; 1 + 2 $0.73; 1 + 1b $0.70); under (b') no rental follows V1 and the money left is a second reserve. Totals (b') = V1 only. Time line step 6, diagram (V2 and V4 "after a raise"), 6.1 step 6, FR1 follow. V2 item 1 re-based on V1's basis: 42 min = $0.58; C adoption 92 min = $1.14; items 1-5 140 min = $1.69; (a) adds $1.14 ($5.79) | `policy.py:480` (cap formula); 6.2 item 1 basis (38.2 + 3 + 1 min); arithmetic at r_s |
| cuda-4, qa-3 | major / blocking | 6.3 V4: plan file `e_v4.yaml`, step 0 = `b1_rows.yaml` step 0 (CUDA build), then `env_check.sh`, the E.9 box list, the sampler-arm `final_eval` (pool 1). Price 30 + 3 + 9 = 42 min = $0.58; (a) total "V4 adds $0.58" | `policy.py:399`; `CMakeLists.txt:15`; `policy.py:432-433` |
| arch-1 | blocking | 6.3 V4: E4 keeps only the old samplers until V4; `final_eval` uses the new sampler from E4 step 3; the trainer path is in git history (tag `pre-e4-eval`); the TOST trainer arm is already on disk | `final_eval.py:73-75`; E4 gate (trainer arm = `runs/r5/eval/eval_val/items.jsonl`) |
| data-1, qa-2 | blocking | 1.1 is the only copy of the table (B.4.1 and 99_back 2b point here). Row 4 = offset < 0.15 guard, edge guard, legal share or the frozen-view hash (DG15); "view identity" removed. New line: the A-G3 correlation is not a Gate B1 rule (B.6.2) | B.6.2 "View identity"; A-G3 row |
| qa-1 | blocking | Order rule 2: pre-tag list adds B's I22a and I47; the tag has no `--numerics`; without it `train_rl` runs bf16, `sequence_mask`, cap 3. 6.1 step 4d now cites I22a | `train_rl.py:23-36` (no `--numerics`); `configs.py:183`; TRL `grpo_config.py:939-940` and clip_max default 3.0 |
| cuda-1 | blocking | Section 5: new paragraph "How the drafts reach the box" (D1 step 0 keep rule in `code_files()`; box paths under `/workspace/src`; no scp/rsync; K1 raw masks built on the box from `checker_eval/masks.npz`). 6.1 step 4 (V2) runs D's S09 box-input test. 6.3 item 4 names `r_d1.yaml`. 6.4 rule 7: "the K1 mask bundle" replaced by the bundle's `checker_eval/` set | `policy.py:379`; `bundle.py:52-70`; K1 GPU_RUN.md:51, K2:33, K3:16 (grep) |
| cuda-2 | major | 6.3 "Own batches" and item 1: the V2 bundle command names `lidc-bundle-v3`, `lidc-train-shifted`, `lidc-train-rescaled`; 6.1 step 4 checks the V2 manifest lists both train views; 6.4 rule 7 adds the LIDC train views named in the rental's table | `bundle.py:55-56, 121` (`--view` repeatable, replaces the default list) |
| data-2 | blocking | 6.2 code-change list: the rows root lands at `{data}/rows/lidc-train-b1/` and the fetch checks its manifest SHA on the box (B1.5 owns the command) | `hf_pull.py:22`; `_snapshot_download.py:152` |
| qa-4 | major | 6.3 V3 item 0: `fetch_hf.sh <run> rl-full-3ep` to `{rd}/ckpt` and `fetch_hf.sh <run> b1-rows --as b1` to `{rd}/b1`; B2.3b reads `{rd}/b1/final_adapter`, `b3_v5.yaml` reads `--ckpt {rd}/ckpt`. 6.2 code-change list names `--as` | `policy.py:534-535` |

## Review log (final round 6)

Cross-file consistency fixer, round 6. Sources: `reviews/final_r6_{cuda,rl,data,arch,qa,fresh-arch}.md` and their
`*_verdict.md` files. Only the CONFIRMED items assigned to this fixer are listed.

| Item | Severity | Change | Evidence |
|---|---|---|---|
| cuda-2, cuda-3, qa-4 | major | 6.3: V2 and V3 each get a launch config (`c_launch.yaml`; `b3_launch.yaml` with `output_dir: b3-v5` and `rows_plan: b3_v5.yaml`) and launch lines. V2 fetches `rl-full-3ep` to `{rd}/ckpt`; V3 item 0 adds `fetch_hf.sh <run> --rows <rows_revision> lidc-train-b1`; V3 pushes with `push_hf.sh <run> b3_launch.yaml`. Section 0 "launch config" names all four files; 6.1 step 4 stamps each before its rental | `policy.py:315-316` (stamp needed), `:511-524` (push revision = `output_dir`), `:527-536` (only fetch path) |
| arch-4, qa-3 | major | 6.3 V4: launch config = `e_v4.yaml` (E6 step 8b); launch lines with `fetch_hf.sh <run> rl-full-3ep` and no `push_hf.sh` | `env_check.sh:9` needs MODEL_DIR; `final_eval_2026-11-08.md:3` (`rl-full-3ep`) |
| cuda-4 | major | V4 re-priced on the V2/V3 basis: 38.2 + 3 + 1 + 2 + 7 + 2.6 = about 54 min = $0.71; approved $0.75 (cap 0.95 h); worst case $0.86. Money-left line and (a) totals follow | ledger.csv r5 launch 15:35:45Z, first `run` 16:13:57Z; `policy.py:480` |
| qa-1 | blocking | 6.3 V4: the TOST runs on the Mac after `pull.sh` (E.9 GP5); E-G10's box half = GP4 on the box + GP5 on the Mac; "GP1-GP5" -> "GP1-GP4" | `.gitignore:23`; `policy.py:376-384` |
| cuda-5 | major | Time line step 5, the flow diagram, V2 item 4 and the R-D1 sentence under the V2 table point to the D1 step 4 tree rule: R-D1 runs only before D2 merges, so it is dropped on the time-line order | `K1/cross_numpy.py:16`; `policy.py:490` |
| data-3 | major | 6.1 step 4: with `rows_plan` set, preflight checks the pinned HF revision (file list and manifest SHA-256, B1.3 change 2b); pass column names the "ROWS revision" check | `hf_api.py:5779-5784` (multi-commit upload) |
| qa-5 | major | Order rule 2: R01 runs on B's fake-vLLM fixture `tiny_loop_vllm_fake` | `test_p2_scripts.py:392-393`; `fast.py:568` |
| arch-5 | major | Sections 0, 5, FR13 and the review logs write the path patterns with brackets; section 0 says why | check 14 scan on the sources: 0 hits |
| fresh-arch-1 | blocking | New "End state under (b')" after time line step 8: what lands on `main`, what waits on the owner branch `after-raise` (C0 V2 tools, C1-C11, B2.2-B2.4, B3), the end state, and the merge rule after a raise. Step 5, step 7, order rule 1 (one named exception) and the flow diagram follow | grep of src/ and experiments/: no C flag and no v5 code exist today |

## Review log (final round 7)

Cross-file consistency fixer, round 7. Sources: `reviews/final_r7_{rl,arch,qa,fresh-arch}.md` and their
`*_verdict.md` files. Only the CONFIRMED items assigned to this fixer are listed.

| Item | Severity | Change | Evidence |
|---|---|---|---|
| qa-2 | major | 6.3 V3 item 1 adds A's G02b and G03 as `d_gpu_test.yaml` step 2: 27-46 min, $0.31-0.52; the launch-line comment names it | D8; A.9.6 |
| fresh-arch-1 | blocking | 6.3 V2: every C arm and the `c_speed_check` batches use the B1 rows root ("One data path"); V2 bundle is the default; item 1 fetches the rows root; item 3 43 min $0.49, 3c 16 min $0.18, 3d 15 min $0.17; required 98 min $1.21, items 1-5 146 min $1.75; owner totals $1.21 ($5.86); `c_launch.yaml` gets `rows_plan`; 6.1 step 4 checks the rows revision for V2 and runs C's U-CF4; 6.2 stamp rule covers plans with several `train_rl` stages; 6.4 rule 7 drops the V2 train views | r_s $0.01133/min, b $0.10; C.8.2 |
| fresh-arch-4 | major | 1.1 rows: "no data effect" -> B1.8 runs (the effect failed, not the path); "Stop. Owner review" -> B1.8 runs after the written review closes, unless the review drops the rows path | B.md B1.5, B1.8 |

## Review log (final round 8)

Cross-file consistency fixer, round 8. Sources: `reviews/final_r8_{cuda,data,arch,qa}.md` and their `*_verdict.md`
files. Only the CONFIRMED items assigned to this fixer are listed. Each fact was checked in the repo first.

| Item | Severity | Change | Evidence |
|---|---|---|---|
| cuda-1 | blocking | 6.3 V3 item 1: every R-D2 output goes to `{rd}/d/` (D8 `$OUT`), so `pull.sh` brings it back; price unchanged | `run_box.sh:16`; `policy.py:489, 505-507` |
| cuda-2, qa-3 | major | 6.3 V3 note: D7 step 5 code writes the D-G14/15/17 records inside the B3 run and its val eval; `pull.sh` brings them back; D's S13 reads them | D.md D7 step 5 |
| data-1 | blocking | No front number used 0.734; nothing changed here | grep of 00_front.md |
| arch-1 | blocking | Section 0 and 6.3: V4 plan file is `specs/e_v4.yaml` (launch config, preflight, `launch.sh`, `run.sh`) | `policy.py:469-470`; E I29 |
| arch-3 | major | 2.2 placement row: the CPU fork pool is never moved to `verify/`; one rule (D7 step 6, E3b R-2, CLEAN16 check) | `select_ckpt.py:374`; `g_probe.py:139`; `test_p2_scripts.py:369, 384` |
| arch-4 | major | 2.2 and 2.3: `core/provenance.py`, `policy/callbacks.py`, `data/audit.py`, `infra/licenses/`; no `infra/paths.py`, no `rl/callbacks.py`; "E's move tables win" line | E.md E2.5, E5, E8 step 1, E3c |
| qa-1 | blocking | 6.2 `a3_measure.yaml` text: build-if-missing step; 6.1 step 4d: every rental's first plan builds or follows a plan that builds; 6.3 item 1b 13.4 min, $0.15 (3 min build) | `policy.py:398-399`; `CMakeLists.txt:15` |
| qa-2 | major | 6.1 step 4: U-CF4 also parses `train_rl` commands inside `bash -c` strings (nsys-wrapped arms) | `policy.py:486`; C0 step 3 |
| qa-4, qa-5 | major | 6.3 item 2b: build 3 + G08b fp16 part 5 = 41.5 / 37.7 min, $0.47; item 4: G08b smoke 5 = 143.6 / 126.9 min, $1.63; B-S2 182.6 min = $2.17 / 165.9 = $1.98 (C code: 150.1 = $1.80); (a) totals $4.71 / $4.33, $5.92 | r_s $0.01133/min, b $0.10; V1 smoke 5 min |

## Review log (final round 9)

Cross-file consistency fixer, round 9. Sources: `reviews/final_r9_{cuda,rl,data,arch,qa,fresh-arch}.md` and their
`*_verdict.md` files. Only the CONFIRMED items assigned to this fixer are listed. Each fact was checked in the repo first.

| Item | Severity | Change | Evidence |
|---|---|---|---|
| arch-1 | major | Order rule 1 and "End state under (b')": after E6 the `after-raise` merge follows E R-4 "After E6" (specs port, launch configs folded, post-E0 flags as fields); V2/V3 launch section notes the `specs/<name>.yaml` lines after E6 | E.md R-4, E6.8, I29; `policy.py:315` |
| qa-6 | major | 6.3 V3: after `pull.sh`, `OBJDEL_V3_RUN=<run> make test-data`; S13 skips without the variable | D.9.6 S13 |

## Review log (final round 10)

Cross-file consistency fixer, round 10. Sources: `reviews/final_r10_{cuda,data,qa,fresh-arch}.md` and their
`*_verdict.md` files. Only the CONFIRMED items assigned to this fixer are listed. Each fact was checked in the repo first.

| Item | Severity | Change | Evidence |
|---|---|---|---|
| cuda-1 | blocking | 6.3 "After E6": every plan path of the V2/V3 launch lines moves to `specs/` (V2 1b `a3_measure`, 3c `c_speed_b`, 3d `c_speed_c6`; V3 1 `d_gpu_test`, 2b `b2_fp16`, 3 `a3_test`); V2 item 4 dropped after E; I29 loads each on the Mac | `policy.py:471` (`cmd_run` calls `chain` after launch); E I29, R-4 |
| data-1 | major | 6.2 new bullet "The edge rule under the fallback": item 4 present -> edge rule pairs with it; item 4 missing -> Gate B1 "open: edge baseline missing", never a pass; V2 item 1b closes it under a raise; under (b') or a second failure -> 1.1 row 4. 1.1 row 4 and V2 item 1b stop rule updated | `runs/r5/eval/` holds `eval_val`, `eval_test` and their logs only [measured: ls, 2026-11-09] |
| qa-2 | major | V2 items 3c and 3d get their own plans (`c_speed_b.yaml`, `c_speed_c6.yaml`) and `run.sh` calls; a `pull.sh` between calls gives the G-MEM peak for 3c; item 3 stop rule acts in-process through `train_rl --expect-step1` (exit 3 stops the `&&` chain); branch file list and 6.1 step 4 U-CF4 name the three plans | `policy.py:432-433` (steps as written), `:487-494` (one `&&` job), `:603` (30-min idle destroy after a pulled EXIT) |
| qa-3 | major | C9 arms: no row in V2 or V3 and no plan file; about 14 min = $0.16 each if the owner adds rows (one home: C.8.2) | C.8.2 8 min + 6 min start-up; r_s $0.01133/min |
| qa-4, fresh-arch-2 | major | One B-S2 definition (V3 items 0 + 4); items 1, 2a, 2b, 3, 5 optional by name (launch lines say "only if approved"); new "V3 amount rule" (spend at the end of 260 B3 steps + $0.35 + item 5): items 0+4 $2.40 (BudgetStop $2.05), items 0,1,3,4 $2.84 expected, $3.08 approved (BudgetStop $2.73); decision 1 (a) totals add $5.38 with items 1 and 3 | arithmetic at r_s: 0.10 + 231.8 x 0.01133 = 2.726; 39+46+13.6+143.6 = 242.2 min = $2.84 |
| qa-5 | major | V2: after `pull.sh`, `OBJDEL_V2_RUN=<run> make test-data`; G-RP2/G-RM1 read `runs/<run>/c_check/batches` and `c_check/arm_batches` and skip without the variable | `policy.py:498-507` (`pull.sh` copies `/workspace/runs/<id>/` to `runs/<id>/`); C.9.9 item 11 (`--record-batches` refuses a non-empty folder) |

## Review log (final round 11)

Skills used (cross-file consistency fixer): nvidia-skills:nvidia-skill-finder (no NVIDIA catalog skill covers these items), engineering:code-review (read the TRL, policy.py, train_rl.py and binding lines cited below), engineering:debug (traced the sleep flag, the `&&` chain after a crash or OOM, and the step-1 stop input), engineering:system-design (one owner per state and per check: TRL's flag; each box test in one plan), engineering:architecture (seam before B1.4, gate pairing in one rental), engineering:documentation (one source for prices: front 6.3), huggingface-skills:huggingface-best (no model choice in scope), mattpocock-skills:codebase-design (`rings.py` as the deep module's first part). Web tools loaded; every fact here is in the repo.

| Item | Severity | Change | Evidence |
|---|---|---|---|
| qa-1 | major | 6.3 V2 item 3 stop rule: the fast arm reads item 2's `{rd}/c_check/encoder.json` through `--encode-check` to pick E-EQ6's `grad_norm` rule | R5 step 1 `grad_norm` 0.0032658 vs 0.0032640 [measured: runs/r5]; C.md C0 files, C3 step 0 |
| qa-2 | major | 6.3 V3 item 0: G-RC0 runs as `b3_v5.yaml` step 2 (minute stays in item 0); item 3d: step 0 is G-GR1 by node id, 16 min = $0.18 | `policy.py:432-433` (`steps:` run as written); 7 x 20/17 + 6 + 1 |
| qa-3 | major | 6.3 item 3c: the mode B arm has its own check steps 18-32 (`--steps 40`, `arm_batches_b`), 25 min = $0.28; the test-data paragraph names `arm_batches_b` | `policy.py:490` (`&&` join); 8 x 40/17 + 6 = 24.8 min x $0.01133 |
| fresh-arch-1 | major | 6.3 item 4 gains `b3_v5.yaml` step 3, the B1 re-eval (13.0 min): 156.6 / 139.9 min, $1.77. B-S2 195.6 min = $2.32 (178.9 = $2.13; C code 163.1 = $1.95), amount $2.55, `BudgetStop` $2.20. Items 0, 1, 3, 4: 255.2 = $2.99, amount $3.23, `BudgetStop` $2.88. Decision 1 (a) $4.86 / $4.48; $5.53 with items 1 and 3; $6.07 with V2. V4 paragraph: a B run after E pairs within V3 | 205 s + 2 x 1,084 x 0.193 s + 2.6 min = 13.0 min; 0.10 + 185.2 x 0.01133 = 2.198; 0.10 + 244.8 x 0.01133 = 2.874 |
| fresh-arch-2 | major | 6.2: plan step (3) soft, new (3b) `train_rl --merge-last`; launch lines load the 20-min pull job before `b1_rows.yaml` and unload it after teardown, with a note (bandwidth already in r_s; 30-min idle destroy); 1.1 row "B1 stopped early"; order rule 2 lists `--merge-last` and I49; FR16; 6.1 step 0 states FR16 to the owner | `com.objdel.vast.pull.plist:11` (1200 s); `policy.py:59-61, 484-490, 502, 603`; `train_rl.py:55-60`; `launchctl list | grep objdel` empty (2026-11-09) |
| fresh-arch-3 | major | 2.2 seam row: B1.4a puts `rings.polygonize` over today's `check_masks` before B1.4; ADR 0015 title drops the exception; order rule 2 lists B1.4a and its hash gate | K1/results/cross_numpy_m3.txt (check_masks IDENTICAL, 0 mismatches on 37,093 masks); gpucheck_py.cpp:140-158 |

## Review log (final round 12)

| Item | Severity | Change | Evidence |
|---|---|---|---|
| qa-1 | major | 2.2 seam row: B1.4a moves the numpy code to `tests/oracles/`, so every ring reference (D-G3 too) stays numpy | B.md B1.4a |
| qa-2 | major | 6.3 V3 item 0: G-RC0 priced here (1 min) but runs as the last, soft step of `b3_v5.yaml`; up to 2 min fits the $0.35 tail; item 4 and the amount rule renumber the B1 re-eval to B3.9 step 2; 6.1 step 4: C's U-RC6 passes (not skips) before V3, and before V2 with item 3d | `policy.py:484-487`; C.md C.9.4 |
| fresh-arch-1 | blocking | Section 5: B1.6 adds the keep rule, D1 step 0 widens it; 6.1 step 4c and order rule 2 name B's I50; 6.4 rule 7 notes the shipped record (hashes, LIDC ids, commits) | `policy.py:379`; B.md B1.6 |

## Review log (final round 14)

| Item | Severity | Change | Evidence |
|---|---|---|---|
| rl-2 | major | 6.3 V3 launch lines load the 6.2 pull job after sync and unload it after teardown; note after the block (reason in 6.2, idle rule between `run.sh` calls, crash covered by B3.9 steps 3 and 3b). V3 item 4 row names step 3b. FR16 covers V1 and V3 | `policy.py:484-490` (`&&` join, `~` soft step); `train_rl.py:55-60` (`final/` only after `trainer.train`); front 6.2 pull-job lines |
| fresh-arch-1 | blocking | 6.3 V4: the old samplers go to `archive/phase2_samplers/` in the E4 step 3 commit, not after V4; tag `pre-e4-eval` keeps them; restore from the tag if V4 fails. "End state under (b')" now holds | `select_ckpt.py:93-96` (only live caller); `parity.py:5,71,89`; `r0_parity.py:7,44`; `sample_vllm.py:47-49,62-63,84,97`; grep of `tests/` (6 files import or name the old samplers) |

## Review log (final round 15)

| Item | Severity | Change | Evidence |
|---|---|---|---|
| arch-1 | major | 6.3 V4 paragraph: "restores the old samplers from that tag" replaced by a pointer to E4 step 1's V4-failure rule (fix forward, no restore). Replaces the round-14 fresh-arch-1 wording | E.md E4 step 1; `select_ckpt.py:93-96` |
| data-2 | major | Order rule 2 (I50), section 5 "How the drafts reach the box" and 6.4 rule 7 now name `docs/reference/val_halves.json` beside `frozen_views.json` (LIDC patient ids, item ids, counts; no pixels) | `policy.py:379-380` (drop `^(data|docs)/`, no keep today); `docs/reference/` holds no JSON today [measured: ls]; `policy.py:398-399` (editable install at `/workspace/src`); B.md B1.6 |


# Phase A

## A.R. Review log

Parts 2, 3 and 4 (the final r1, r2 and r3 reviews) are at the end of this file, after A.9.

Part 1. Source: the first phase review of A (30 items; a session file, so this log restates each item).
B = blocking, M = major, m = minor.

| # | Sev | Item | What I did |
|---|---|---|---|
| 1 | B | A-G6 cannot be computed from `nvtx_sum` | Applied. `nsys_capture.sh` exports sqlite and the CSV reports `nvtx_pushpop_trace, nvtx_gpu_proj_trace, cuda_gpu_trace, cuda_api_trace`. New `nsys_baseline.py` (coverage by containment nesting, GPU busy by interval union, host syncs with `_v` suffix removed), CSV fixture test, gate stated on `nsys_baseline.json` (A6 steps 8-9, A-G6) |
| 2 | M | no `--trace` in `train_rl`; wiring is in `make_trainer`; `profile.py` already computes busy share | Applied. Verified [code: train_rl.py:23-36; g_probe.py:92-95; profile.py:77-99]. `make_trainer` wires `NvtxSteps` and `TraceStep` (wait from the window); `train_rl` gains both flags; the interval union moves to `profile.covered`, used by both paths. I kept it in `profile.py` rather than a new library module: both callers are experiment scripts (A6 steps 5, 9, 10) |
| 3 | M | `step` range misses loader and log; backward has no range | Applied. Verified [code: transformers trainer.py:1856,1877,1929-1931]. `step` = `on_step_end` to `on_step_end`; profiler start at the end of step s - 1; new ranges `micro` (`training_step`), `loader` (`get_batch_samples`), `log` (A6 steps 3-4, 6) |
| 4 | M | version-first order makes the signature checks dead | Applied. Contract first, version last, `OBJDEL_TRL_UNPINNED=1` opt-in with a warning, plus its test (A7 step 2, 5) |
| 5 | M | silent couplings not checked | Applied. 13 `ATTRS`, 4 `KEYS`, 2 behaviours; each verified by substring count this session. The ratio behaviour test `exec`s TRL's own block; I prototyped it and all 4 modes equal `fast.is_ratio` [measured]. `llm.generate/wake_up/sleep` left to `env_check.sh` (vLLM public API, absent on the Mac). Count updated to 29 |
| 6 | M | import-time `is_available` breaks fork safety | Applied. Verified `SCORERS` imports `objdel.reward` pre-fork [code: serve.py:52] and the guard [code: serve.py:167-169]. `nvtx.py` reads only the env var at import; lazy cached CUDA check; spy test (A6 step 1, 11) |
| 7 | M | R5 fixture is git-ignored | Applied. Verified [code: .gitignore:23]. `tests/goldens/r5_test_items.jsonl.gz` (31,314 bytes measured; `tests/goldens/` is the existing folder, not `tests/fixtures/`); patient test uses a synthetic map (A1) |
| 8 | M | answers written after scoring | Applied. Answers + partial `final_eval.json` written before `score_items`; `--rescore` accepts a partial file; crash test (A3) |
| 9 | M | nsys hang has no time limit | Applied. `timeout --signal=INT --kill-after=120 900`, step stays last, `pkill` cleanup. Not applied: making `pull.sh` copy evals before the step; the answers are already on disk before it, so this adds runner complexity for no gain |
| 10 | M | nsys install vague; Blackwell unchecked | Applied with a correction: the devtools repo the review names stops at 2021.4.1 [measured]. Pinned `nsight-systems-2026.5.2` from the CUDA Debian 13 repo with URL, size and SHA-256 [measured]; `dpkg-deb -x` CLI unpack; `nsys status -e` gate; `curl -sI` on the Mac first. Blackwell support marked [estimate] |
| 11 | M | `zero_share` mixes two failures | Applied. `raw_ratio` + `bounded` split; `zero/hi/lo` counts (lo as one OR, so no double count, unlike the review's sum); verified TRL defaults clip_min None, clip_max 3.0 |
| 12 | M | A-G3 sanity rule too weak | Applied. Paired by `item_id` vs R5, CI of the difference contains 0 and |mean diff| <= 0.01; correlation and old rule as secondary lines |
| 13 | M | A8 ignores model answers | Applied. `--answers` with `answers()` (share > L, legal share, mean IoU_tc, legal-and-good count) |
| 14 | M | A-G8 cannot close in A | Applied, and simpler than the review assumed: the train views already exist on the Mac (7,341 records each, built 2026-11-08 [measured]). A8 runs on all three roots now and records the manifest SHAs; A-G8 closes in A |
| 15 | m | `sel_medoid_legal` gate is a tautology | Agreed. Replaced by the test-direction check; method chosen on val only |
| 16 | m | legacy check cannot fail; pooling unstated | Agreed. Two pooling orders; fail stops A2 until the owner accepts |
| 17 | m | medoid parsed fallback is dead logic | Agreed. Verified [code: metrics.py:306-307]. `parsed` dropped; empty-majority note and test |
| 18 | m | preflight check name and signature | Agreed. Verified [code: preflight.py:244-250,266-280]. `"CODE TRL private names"`, `fn(a, cfg)` returning a list |
| 19 | m | `configs.processor` does not exist | Agreed. Verified [code: configs.py:70,134]. Tokenizer built as in `sample.export_prompts`, offline |
| 20 | m | preflight command incomplete | Agreed. Verified [code: Makefile:29-30]. Full commands in the plan header |
| 21 | m | `code_sha` always null on the box | Agreed. Verified rsync without `.git` [code: policy.py:374-401]. `cmd_sync_up` writes `CODE_SHA`; `final_eval` reads it first |
| 22 | m | docs quote item CIs; `ci95` meaning changes | Agreed. Verified [measured: grep]. Keys renamed `ci95_nodule/_patient/_item`, bare `ci95` dropped; both docs added to A9 |
| 23 | m | "CT series" CI is not what A1 computes | Agreed. Re-measured with `group_id`: patient CI [0.6252, 0.6644], 158 patients vs 160 series on test [measured] |
| 24 | m | boundary IoU d = 2 on small nodules; margin | Agreed. Kept d = 2, added size-tercile report, margin clipped to the frame with an edge test |
| 25 | m | doctor masks vs polygon text | Agreed (no design change). Added the rasterize-equals-mask test to A2 |
| 26 | m | test 2 can be exact | Agreed. Re-measured row-for-row equality (n = 7, 2,200); test asserts exact equality |
| 27 | m | host-read marks redundant | Agreed. Dropped `mark("host_read")`; fewer `fast.py` edits |
| 28 | m | set metrics read masks from a second rasterizer | Partly disagreed: `to_masks` already calls `_gpucheck.check_polygons` [code: shape.py:136-148], and equality is tested [code: tests/test_gpu_scoring.py:51-57]. I stated this in §3.1-3.2 and named the test. Phase D swaps only the arithmetic |
| 29 | m | profile command without `-I`; graph kernels hidden | Agreed. `python -I -m` in `nsys_capture.sh`; default `--cuda-graph-trace=graph` and the note for `profiling.md` |
| 30 | m | no check of nsys overhead | Agreed. `time/step_s` logged by FastGRPO; steps 4-5 against 3 and 6, flag at > 10 % |

---

## Review log, part 2 (final r1 reviews, 2026-11-08)

Sources: the six final r1 review files (`final_r1_arch`, `_cuda`, `_data`, `_fresh-arch`, `_qa`, `_rl`), items
addressed to A. Each fact was checked against the code before the change. `_rl` has no item for A.
Sev: B = blocking, M = major, m = minor.

| Review # | Sev | Item | What I did |
|---|---|---|---|
| arch 2 | B | A6's `checker` range in `reward/score.py` makes `reward` import `train`, breaking E0's contracts | Applied. Verified `score_many` calls `score_polygons` at score.py:115 and `TRLReward` calls `score_batch` at rewards.py:119,143. The range moves to a `_score` helper in `train/rewards.py` (both the prefetch thread and the direct path). `reward/score.py` is no longer edited. New S12 (AST scan, with a planted negative case) and U57. E destination `train/nvtx.py -> policy/marks.py` stated in A.8 |
| arch 6 | M | E's move tables miss modules A creates | Applied for A's modules: A.8 lists the rows E must add (`eval/ceiling.py`, `eval/select.py`, `train/nvtx.py`, `train/trl_adapter.py`). The table edits themselves are E's |
| arch 15 | m | "TRL 1.16" fix has three owners | Applied. DOC2 owns the wording; A7 no longer edits docstrings; S06 asserts that no live `.py` file under `src/` or `experiments/` says "TRL 1.16". Grep found a 5th file the notes missed: `experiments/phase2/train_sft.py:49` |
| arch 17 | m | duplicate creation of `tests/test_repo_hygiene.py` | Applied. A creates it; CLEAN12 extends it (A.8, A.9.0) |
| cuda 6 | M | two nsys installs, undefined `OBJDEL_NVTX_SYNC`, two capture flags | Applied. `nsys_capture.sh RD NAME S E -- CMD` is the one mechanism (pinned `.deb`, SHA-256, `dpkg-deb -x`); C0 calls it; a Dockerfile layer may unpack the same `.deb`, no second source. `OBJDEL_NVTX_SYNC` defined in A6 step 2 (1 = today's sync, 0 = CUDA events read in `log`), with U56. `--profile-window` is the one capture flag; `--profile-steps` is not added (S04 checks) |
| cuda 7 | M | A3 runs last in G1, but C's gates need its answers | Applied on the A side: A3 is the first G1a item after setup (A.4, A3 step 3); A3 step 6 states that the Mac half of G-EQ2/G-EQ3 runs on A3's answers after the rental. The front table edit is the front's |
| cuda 3 | (D/E) | device work in forked pool workers | A side only: the ceiling CLI runs on the Mac (CPU build); on a CUDA build it uses pool size 1 (A2 step 2) |
| data 2 | B | gates resample nodules; the patient is the independent unit | Applied. Measured 78 val and 158 test patients. A1 step 5: `ci95_patient` is the gate CI (`gate_ci` key), cluster count printed, cluster t interval below 50 patients (`signal.cluster_t`). A-G3, A-G4, A2 strata and the selection pairs use patient draws; nodule CI is a second line. New U53 (perfect within-patient correlation: patient CI wider). The B and E gate texts are those phases' edits |
| data 3 | B (front) | "miss" has three meanings | Applied on the A side. Verified `per_item` stores the K-mean [code: final_eval.py:24-30]. The header "Words" block defines miss (K-mean rule), found image and missed answer; A5 step 4 names the K-mean rule; A9 adds both words to the glossary |
| data 5 | M | A-G4 lets test veto the method | Applied. "Same sign as on val" deleted from A4 and A-G4; the test difference is a direction line, not a pass condition; U51 states that `test_same_sign` never feeds `change_default` |
| data 11 | M (front, B) | the model decision uses the unfair ceiling | Applied on the A side: A2 adds `found` and `miss` strata with patient CIs on val, and states that the found number bounds found images only. The front's "[to be confirmed by A2]" text is the front's edit |
| data 14 | m | A-G3 "CI contains 0" is not an equivalence test | Applied. TOST: the 90 % patient CI of the paired difference inside ±0.01; a wider CI is "not shown" (A3 step 5, A-G3, G01, diagram) |
| data 21 | m | the R5 test view was never hashed | Applied. Verified no manifest SHA in R5 `final_eval.json` and no local view. A3 step 1 writes the B1.6 content hash, `target_set_sha256` and `generator_commit` to `runs/a3/pre_rental.json` with the note "equal by code identity, not measured"; `final_eval.json` gains `view_content_sha256`; new `--hash-view` flag and U54 |
| data 22 | m | A2 strata use raters, the left-out set uses voters | Applied. Verified `voters()` [code: cache.py:108-111]. k = `len(cache.raters)`; averages over voting doctors only; non-voters stay in the stacks; `n_voters_lt_raters` per stratum; new test case in A2 step 3 |
| fresh-arch 13 | m | CIs cover data variance only | Applied, adapted to $0: A1 states "CIs exclude training-seed variance"; A3 step 7 computes a sampling-noise floor from the two halves of the saved K = 8 answers (scaled by 1/√2, marked [estimate]); a second val sampling (415 s) is costed but not planned. U55 |
| qa 1 | B | the G1 rental cannot run as planned; A3 has its own launch block; A-G6 window | Applied on the A side. A3's header has no `launch.sh`/`--cap-hours`; it names G1a; `teardown.sh && sweep_orphans.sh` added. The "All" gate reads "A3 runs inside the G1a cap". A-G6 takes option (a): the C0 baseline arm capture (17 steps, `--profile-window 16 17`), overhead compared with steps 13-15 (the arm ends at 17, so the review's step 18 does not exist). `a3_measure.yaml` drops the nsys step; I20 and the sketch updated. The G1a/G1b split itself is the front's |
| qa 5 | M | the doc depends on scratch files | Applied for A. No `S[P]/` path remains: the spec and session notes are cited at their DOC-promoted repo paths (header); kernel designs at `docs/implementation2.0/kernels/K*/DESIGN.md` (exist in the repo tree); review files named, not linked. The copies are DOC's task |
| qa 8 | M | different pytest markers | Applied for A. One set: `gpu` and `data`. Verified `data` is registered [code: pyproject.toml:68-71] and `gpu` is not; A registers `gpu` and `sanitizer` with `tests/test_a_gpu.py` and adds `--strict-markers` (S03 checks) |
| qa 10 | M | open behaviours with tests that accept either outcome | Applied. A.9.8 now has one rule per item (default "raise"); U11, U12, U21, U29, U41, U47, U49, I10, I16 assert only that rule. Rule 6 follows the code: `parse` makes a truncated answer illegal [code: text.py:106-107] |
| qa 12 | M | the B gates do not name the K-aggregation | Applied on the A side: A4 states that every B run gate uses `sel_mean_of_k` whatever Gate A picks, and that Gate A is decided on the Mac after G1a, before the B spend |
| qa 14 | M | synonyms | Partly applied: the header "Words" block fixes one word per thing and names the code names (rater, sample/completion, window). Code identifiers keep their names. A full synonym sweep of A's prose was not done |
| qa 15 | M | several prices and step times | Applied for A. Verified $0.440 is R3's billed total, not a rate [code: gpu-run.md:129]; A.8 cites the front's section 6 rate and gives A3's cost at $0.641/h. R5 step = 23.65 s [measured: runs/r5/shared_rl.json] in the header, A.1 and A6 |
| qa 16 | M | "-I for data/" cannot run in the tests | Applied. Verified `make test` has no `-I` [code: Makefile:19-20]. A changes `test` to deselect `data` and adds `test-data` with `-I`; S03 checks the Makefile and that data-reading tests carry `data` |
| qa 17 | M | long sentences | Partly applied: the new and changed text uses short sentences and one action per step. A full rewrite of A was not done |
| qa 20 | m | two homes for the NVTX module | Applied with a different choice: `policy/marks.py`, per arch item 2 (an ML layer, because `enabled()` imports torch). `core/` would break `pure-core`; `infra/` is dropped |
| qa 21 | m | two phases create the same test files | Applied: CLEAN12 extends `test_repo_hygiene.py`; E5 `git mv tests/test_trl_adapter.py tests/rl/` (owner runs); E-G5 ">= 5 negative cases (A-G7)" stated in A.8 |
| qa 23 | m | bit-exact CIs need a numpy pin | Applied. Verified `numpy>=2.0` in pyproject. New `tests/goldens/r5_test_items.meta.json` (numpy version, machine); `==` binds there, else `abs <= 1e-12` (A1, A-G1, I01) |
| qa 24 | m | A and C disagree on sampling repeatability | Applied. Verified R5 step 1 equal across two arms [code: shared_prefix_2026-11-08.md:61-64]. A.7 row and G07 cite it; A-G3 keeps the ±0.01 TOST for another host |
| qa 26 | m | A3 rebuilds the test view before B1.6 freezes it | Applied with data 21: A3's content hash is the value B1.6 freezes and must reproduce (A3 step 1) |
| qa 28 | m | undefined tags | Applied. Header defines `[web: ...]` and `[unverified]`; "measured: this session" tags became `[measured: ..., 2026-11-08]` |
| qa 2, 3, 4, 6, 7, 9, 11, 13, 18, 19, 22, 25 | - | front, B, C, D, E, DOC, CLEAN or 99_back only | Not A's text; no change here |
| cuda 17, qa 27 | m | front toolchain lacks `nvtx` | Not A's text (front, D) |
| data 26 | m | DOC data card line cites item CIs | Not A's text (DOC). A9 already relabels the item CI and adds the patient CI |

## Review log, part 3 (final r2 reviews, 2026-11-08)

Sources: `reviews/final_r2_{arch,cuda,data,fresh-arch,qa,rl}.md`, items addressed to A. Each fact was checked against
the code or a file before the change. Sev: B = blocking, M = major, m = minor. Parts 1 and 2 above are history: their
"G1a" and "A3 first in the rental" texts are superseded by this part.

| Review # | Sev | Item | What I did |
|---|---|---|---|
| fresh-arch 3, cuda 2, qa 1 | B / M / B | A still plans the removed rental G1a; the V1 nsys window exists in no A text | Applied. Every live "G1a" is now V1 (A3) or V2 (C0 capture): A.3.1 diagram (new V2 box), A.4, A3 heading/header/steps, A4 Gate A line, A6 steps 8-9, A-G3, A-G6, "All" row, A.7, rollback, A.8, A.9.6, G06, coverage table. The three reviews disagree on the capture. I took option (a) of cuda 2 and qa 1 (2 of 3 reviewers): A-G6 is split; the Mac half passes in A; the box half reads the C0 baseline arm of V2 item 3 and stays open until then. Reason: no A config defines a V1 window (grep finds no `--profile-window 4 5`), and fresh-arch 3's window would add minutes to a V1 that is already near its cap. The front table and 99_back open question 8 must say the same (not A's text) |
| fresh-arch 3, qa 16 | B / M | A3 preflights `r1_qwen3vl2b.yaml`; V1 runs two configs | Applied with a correction. `preflight.load` builds `train_sft.RunArgs(**raw)` [code: preflight.py:49-56], so `make preflight CONFIG=a3_measure.yaml` (or any `steps:` plan) cannot work. R5 launched with `r1_qwen3vl2b.yaml` and ran its plan with `run.sh` [code: r5_shared.yaml:9-13]. A3 now names "the V1 launch config (front 6.1 step 4)", not a plan file. The preflight runs on the tree that holds A3, A6 and B1.5. The stamp's data hash covers only `preflight.EXAMS` [code: preflight.py:28]; test and edge views are checked by `final_eval` against `frozen_views.json`. I20 and the sketch assert no plan file follows `preflight CONFIG=` |
| fresh-arch 2, qa 2, data 1, rl 10 | B / B / B / M | V1 does not produce the R5 edge baseline that Gate B1 needs | Applied on the A side. `a3_measure.yaml` has three evals in this order: val shifted, val edge (new), test. The two val steps are required; the test step is last and may move to V3 (A-G3 test part "open (V3)", A.9.8 rule 14). Verified the bundle ships every `data/eval_views/lidc-*` folder [code: bundle.py:53-54], so `lidc-val-edge` reaches the box once B1.6 builds it. A.8 minutes: val part 13.8 min, test 10.6 min, all 24.4 min [measured: R5 `seconds.total` 415 s and 634 s; edge time an estimate]. The V1 table, the cap and B-S1 are the front's and B's edits |
| qa 3, data 2 | B / M | A and B use different decision intervals; no task builds BCa | Applied in A1 step 5: the A1 rule is the one rule for every gate in A to E; no BCa. Subset rules use the patients with items in the subset and print that count. Paired rules below 50 patients run `cluster_t` on per-patient mean differences. If the owner wants BCa, A1 adds `signal.bca` with a test against `scipy.stats.bootstrap(method="BCa")`. A9 `metrics.md` says "no BCa". The B and E texts are theirs |
| data 3 | M | the Gate B1 pairs have no proof of the same view | Applied on the A side. A3 step 1 states that A3 reads the frozen val views (B1.0 hashes the local val view, built 2026-10-20, `target_set_sha256` 2da645bf... [measured: manifest.json], before any rebuild; B1.6 freezes the edge view). New hard view-identity rule in A3 step 5, A-G3 and G01: per-item Pearson correlation with R5 >= the lower 95 % patient-bootstrap end of the Spearman-Brown `r_sh`. New U58. A.4 and A.8 say Gate B1 can use A3's val items as the paired baseline. The B1.6 key rename is B's |
| data 5 | M | the noise floor is undefined | Applied. A3 step 7 defines `noise_floor.item_sd` and `noise_floor.mean_sd`, says both are K = 8 sampling noise inside every paired patient CI, and no gate uses them as a pass rule. A1 and A9 say training-seed variance is not measured. The A.7 seed row now says "limitation; a second-seed B1 run is an owner decision, not in the $3". U55 tests both keys |
| data 7 | M | the test hash is frozen before the contour-to-slice check | Applied. A3 step 0: B1.7 runs before step 1 hashes any view; test hits keep their items, go to `known_label_issues`, and add a sensitivity line; a later label fix needs the owner in writing. The header comments start with `check_series_z` (I20 checks). DG17's text is B's |
| qa 7 | M | `make cov` and `make test-all` still run data tests without `-I` | Applied. Verified `Makefile:19-26` and 39 `data` marker uses [measured: grep]. `test-all` and `cov` deselect `data`; new `test-all-data`; `test-data` is `data and not slow`. A runs `make cov` once and records the percentage; below 90 % it becomes D-G12's floor. S03 checks all `pytest` lines; the "All" gate records the number |
| qa 11 | M | G06's NVTX on/off 3 % check has no data source | Applied. Clause deleted from G06 (and so the "23.6 s" spelling, qa 24). NVTX cost = Q03 plus the 16-17 vs 13-15 overhead line. A fourth `OBJDEL_NVTX=0` arm is named as an owner option in V2 |
| qa 15 | M | S06 and S08 have no run condition | Applied. S06: green when A7 lands because DOC1 and DOC2 run first (front time line step 1); if DOC2 is late, A7 waits; never skipped. S08: `pytest.skip("A9 pages not written yet")` until `metrics.md` holds "fair ceiling" |
| arch 6 | M | E5 changes A7's contract; S07 asserts the literal 29 | Applied on the A side. A7 adds `trl_adapter.ALL_ITEMS`; S07 and the `trl-coupling.md` row check read its length; the 29 is checked once as the Phase A value. A7 notes that E5 must replace the `shuffle_sequence_dict` negative case. E-G5 and E5 are E's |
| arch 5 | M | E5's table uses A8's old names | A side only: A.8 states that E's move tables and `old_names.json` must use `target_tokens` and `voting`. The table is E's |
| cuda 12 | m | "the devtools apt repo stops at 2021.4.1" is false | Applied. Verified `devtools/repos/ubuntu2404/amd64/Packages`: `nsight-systems-cli-2026.6.1`, 232,877,734 bytes, SHA-256 `e7144097...b56b`, `Depends: libc6, libglib2.0-0`, HTTP 200 [measured: curl]. A.0 and A6 corrected. The pin stays the Debian 13 package (the image is Debian 13); the CLI-only package is the script's pinned second try (S04 checks its SHA). Part 1 item 10 is history |
| cuda 11 | m | C0 adds `--cuda-memory-usage` to its own nsys line | A side: A6 step 8 says callers add no other nsys flag and names `--cuda-memory-usage` as excluded; S04 checks the script has none. C0's command line is C's edit |
| data 16 | m | "code identity" names only the generator files | Applied. `--hash-view` also writes `render_files` (targets `augment`, `frame`, `build`, `simplify`; `train/dataset.py`; `build_eval_views.py`, each with its last commit) and `r5_commit: null` with "not recorded" (no commit under `runs/r5/` [measured: grep]). Verified `build_target` is in `targets/build.py:175` |
| qa 17 | m | A3 reads the test split again | Applied on the A side: A3 step 5 says the test re-score is an equivalence check, not a selection. The section 0 exception is the front's |
| qa 19 | m | 1,182 vs 1,181 + 1 | Applied. A counts the end token (A8 step 2; A1.md mean 128.2 vs B1.md 127 without it [measured]). A.1, A8 step 5, A-G8 and G05 say "with the end token"; A8 says DG12 must not add "+1" again |
| qa 20 | m | test-source paths point at `tests/` | Applied: `docs/implementation2.0/tests/A_tests_{1,2}.md` |
| qa 21 | m | "exam" still used | Applied at A.3.2 and A.7 ("test evaluation"). Code docstrings keep their words |
| qa 25 | m | long sentences | Partly applied: the new text uses short sentences and bullet lists. No full rewrite of older text |
| fresh-arch 1, 4-9; cuda 1, 3-10, 13-20; data 4, 6, 8-15, 17, 18; qa 4-6, 8-10, 12-14, 18, 22, 23, 26; arch 1-4, 7-20; rl 1-9, 11-18 | - | front, B, C, D, E, DOC, CLEAN or 99_back only | Not A's text; no change here |

## Review log, part 4 (final r3 reviews, 2026-11-08)

Sources: `reviews/final_r3_{arch,cuda,data,fresh-arch,qa,rl}.md`, items addressed to A. Each fact was checked against
the code, a file or a fresh measurement before the change. Sev: B = blocking, M = major, m = minor. Part 3's
"A3 has val shifted, val edge and test" text is superseded by this part.

| Review # | Sev | Item | What I did |
|---|---|---|---|
| data 1, qa 2 | B / B | A plans a different V1 than front 6.2: `a3_measure.yaml` runs val edge again and the test split; V1 then costs $2.27, above the $2.20 cap | Applied. Front 6.2 is the one list. `a3_measure.yaml` keeps one `final_eval` step (R5 on `lidc-val-shifted`, front item 5). R5 on `lidc-val-edge` is front item 4, a `b1_rows.yaml` step on the A3 code (B1.5 names its folder). The test step moves to new `a3_test.yaml` (V3 item 3). Changed: A3 title, Files, plan text, header, steps 3-5, A.3.1 diagram, A.4, A4 Gate A line, A8 CLI, A-G3, "All" row, A.7 row, A.8 (A3 in V1 6.9 min = $0.07; test half 10.6 min = $0.11), G01, A.9.8 rules 14-15, I20 and sketch 8 (no `--split test`, no val edge in `a3_measure.yaml`). Verified 415 s and 634 s [measured: R5 `seconds.total`] |
| arch 1 | B | V1 cannot get a preflight stamp: front 6.1 step 4 names two `steps:` plans; A3 points back at it | Applied on the A side. Re-verified `preflight.load` builds `RunArgs(**raw)` [code: preflight.py:49-56]. A3's header and step 3 name "the V1 `RunArgs` launch config of front 6.1 step 4; never a plan file", and state that the review proposes `b1_launch.yaml` (B1.5). The front and B edits are theirs |
| qa 5 | B | When V1 skips A3-val, Gate C-adopt waits on V2 for answers V2 makes | Applied on the A side with option (a): A3 step 6 says the C items may enter V2 after the C1-C5 Mac tests; the Mac half of G-EQ2/3 runs on A3's val answers after V1 or after V2 item 3b, and must pass before V3. A.8 C line and A.9.8 rule 15 match. C.6, C.8.2 and front 6.3 own the rule text |
| data 2 | M | Gate B1's primary baseline is named two ways; A-G3 failure has no consequence | Applied on the A side. A.4 no longer claims the baseline; it says front 6.2 and B.6.2 own it. A3 step 5 and A-G3 state what a val failure means: view identity drops only the R5 second line; a TOST failure leaves Gate B1 on A3's items with the A3-R5 difference reported. The item-5 skip rule change is the front's |
| data 4 | M | `cluster_t` on per-patient means estimates the patient-weighted mean, not the item mean of the point rules | Applied. Re-measured on R5 val: item-weighted 0.6182, patient-weighted 0.6472; 2-78 items per patient; top 5 patients hold 24.9 %; CR1 interval [0.5837, 0.6527] [measured: 2026-11-08]. A1 step 5 defines `cluster_t` as the CR1 interval with its formula; paired and subset rules apply it to items. Tests: equal `N_p` equals the patient-mean t interval; unequal `N_p` is centred on the item mean (U53, sketch 1). A-G1 and A9 `metrics.md` updated. B.6.2 and B.8.4 wording is B's |
| data 9 | m | The gate CI has no golden test | Applied. The golden gains `group_id` (public LIDC patient id); re-measured size 33,761 bytes (gzip 9, `mtime=0`) and the patient CI [0.6252288200704239, 0.6644106494958667] from it, equal to A.1 [measured: 2026-11-08]. A1 Files and step 6, A-G1, I01, S09 and sketch 1 pin `ci95_patient` |
| cuda 4 | M | `gpu_busy` omits vLLM decode run as CUDA graphs; `host_syncs` misses driver-API syncs | Applied. Verified TRL builds vLLM with no `enforce_eager` [code: vllm_generation.py:356-372]. A web check found the nsys SQLite graph table name differs by version, so A6 does not guess it: `nsys_capture.sh` writes `${NAME}_tables.txt`; `nsys_baseline.py` adds rows of `GRAPH_TABLES` (start: `CUPTI_ACTIVITY_KIND_GRAPH_TRACE`, marked [unverified]), drops zero-time rows, and sets `graph_rows: "absent"` plus `gpu_busy_is_lower_bound` when graph launches exist but no table does. Graph rows are an upper bound on graph segments (stated). Driver sync names added; counts split vLLM / trainer by range time. Fixture `gpu_rows.json` (SQLite built in the test, no binary); I12, I14, A-G6, G06, A.3.2, A9 updated. C0's `gpu_active.py` is C's |
| fresh-arch 5 | M | each val eval pays about 205 s of fixed load | Applied as code only. Verified the fit: (630.5 - 414.7) / (2,200 - 1,084) = 0.193 s per item, 205 s fixed [measured: R5 `seconds.load_and_sample`]. `final_eval --view` now appends (one trainer build, one folder per view, own `seconds`); U59 tests it. The engine seed is set once [code: vllm_generation.py:366], so a second view's answers differ from a fresh process (stated; A-G3 is a TOST). Whether front 6.2 merges items 4 and 5 is the front's choice; A.8 gives the merged price (about 10.4 min) |
| qa 10 | M | `tests/test_a_gpu.py` has no per-test markers; G06 has no runner | Applied. A.9.0 and A.9 intro: one marker line per test, no module-level `pytestmark`. A.9.6 gives each test its marker, run place and command. G02 split into G02a (`data`, Mac rescore) and G02b (`gpu`, box, both builds). G06 becomes `data` on the Mac after `pull.sh` of V2 item 3; its 3-step toy is dropped (no V2 step runs pytest) |
| qa 9 | M | `pytest -m gpu` collects every phase's box tests | Applied on the A side: A.9 intro says `pytest -m gpu` with no file list never runs on a box; each A box test names its node id and rental item (G02b and G03 in V3 item 1 as their own command; G07 optional in V3 item 3; G08 with D's CUDA path); a missing `runs/` input skips with a reason. D-G3 and D8 are D's |
| data 6 | M | A val label hit has no rule that fits the frozen view | Applied on the A side: A3 step 0 states the val rule (keep items and labels; `known_label_issues`; excluded in both arms of a paired rule; a fix only before B1.0 hashes). DG17 is B's |
| rl 2, rl 5 | B / M (B) | add `llm_engine.add_lora` and the OPSM threshold attribute to the A7 check | Applied as a rule, not as items. Verified `self.off_policy_mask_threshold` at grpo_trainer.py:882 (set) and :3190 (read). A7 now says the commit that starts to use a name adds it: B adds `off_policy_mask_threshold` to `ATTRS`; `add_lora` is vLLM, so it goes to `env_check.sh`, not `ALL_ITEMS`. A adds neither (both belong to code B may not adopt); S07 reads the length |
| data 7 | m | Val already chose R5 | Applied on the A side. Verified the note picks the final model on 128 val shifted items [code: full_rl_2026-11-08.md:31,49]; 51 of 78 patients per the review. A3 step 5, A-G3 and a new A.7 row: the bias favours R5; a sensitivity line without the 128 items. B.6.2 and DOC11 are theirs |
| data 8 | m | B1.7 already has a test hit; the owner decision is due before A3 step 1 | Applied on the A side: A3 step 0 names LIDC-IDRI-0659 (4 test items) and the four train hits, tagged as the review's measurement (not re-measured here), and the fix-or-keep choice for the owner before step 1. B1.7 "Today" and 99_back are theirs |
| qa 12 | m | three GPU skip rules | Applied. A.9 adds one root `tests/conftest.py` hook (the file exists): `gpu` tests skip unless `objdel._gpucheck.build() == "cuda"` [code: cuda/gpucheck_py.cpp:113]; no torch import. C, D and E point at it (their edits) |
| qa 13 | m | "outline length" has two definitions | Applied on the A side: A8 step 5 names the one config limit `L - soft_punish_cache >= 1,182 + 32` (the answer cap counts the end token) and says front section 0 holds the definition. B3.7, C7 and decision 13 are theirs |
| data 12 | m | stop rule 7 does not list the eval views | A side only: A3 says the views are LIDC-only and that front 6.4 rule 7 must name them (the front's edit) |
| qa 19 | m | review logs follow part 9 | No change in A. The fix is front section 8 ("Review logs of the final rounds may follow part 9"). The A.R pointer now names parts 2-4 |
| qa 21 | m | long sentences | Partly applied: all new text uses short sentences and one rule per bullet. No full rewrite of older text |
| fresh-arch 1-4, 6, 7; cuda 1-3, 5-13; data 3, 5, 10, 11; qa 1, 3, 4, 6-8, 11, 14-18, 20; arch 2-17; rl 1, 3, 4, 6-15 | - | front, B, C, D, E, DOC, CLEAN or 99_back only | Not A's text; no change here |

## Review log (consistency pass, round 4 prep)

| Item | What changed |
|---|---|
| One source (r3 data 1, qa 2) | A.4 "Rental", A3 header comments, A.7 row and A.8 point to front 6.2-6.3 and hold no minutes or dollars |
| V1 items 4-5 | `a3_measure.yaml` holds both R5 val views in one process (edge, then shifted), with one `--out` per view (`eval_val_edge`, `eval_val`); the old text put item 4 in `b1_rows.yaml`, which the front no longer does |
| Skip rule | Item 5 is not skippable; the fallback (no output) is V2 item 1b, not "V2 item 3b"; A3 step 6, Gate A, A8, A-G3, G01, edge case 15 follow |
| Preflight | A3 header names `make preflight CONFIG=experiments/phase2/configs/b1_launch.yaml`; the stale "front names plan files" paragraph is gone |
| Tests | I20 and the pytest sketch expect one `final_eval` step with two val views; U59 covers repeated `--out` |

## Review log (final round 4)

Cross-file consistency fixer, round 4. Sources: `reviews/final_r4_{cuda,rl,qa,arch,fresh-arch}.md` and their verdicts.

| Item | Severity | Change | Evidence |
|---|---|---|---|
| cuda-2, qa-2 | blocking | A.9 hook: on a box (`OBJDEL_RUN_DIR` set) a `gpu` test on a CPU build fails, not skips; Mac unit test. `a3_measure.yaml` asserts `build() == 'cuda'` before `final_eval` (I20 extended); `a3_test.yaml` builds first. A3: `final_eval` caps the pool at 1 on a CUDA build | `policy.py:484-492`; `final_eval.py:69`; `serve.py:132` |

## Review log (final round 5)

Cross-file consistency fixer, round 5. Sources: `reviews/final_r5_{cuda,rl,data,arch,qa,fresh-arch}.md` and their
`*_verdict.md` files. Only the CONFIRMED items assigned to this fixer are listed.

| Item | Severity | Change | Evidence |
|---|---|---|---|
| qa-6 | major | G02b: the CPU side is a Mac golden `tests/goldens/a_g02b.json` from `experiments/phase2/perf/make_a_goldens.py` (D8 pattern); the box test compares the CUDA build with `==`; a missing golden fails | `gpucheck_py.cpp:113`; `CMakeLists.txt:15`; `policy.py:376-383` (`tests/` ships) |
| data-1, qa-2 | blocking | A3 step 5 and the A-G3 row: a failed TOST or correlation does not block Gate B1; a `view_content_sha256` mismatch is a Gate B1 hard rule (B.6.2; front 1.1: stop) | B.6.2 |
| fresh-arch-1 | blocking | A.3.1 diagram and A-G6: "open until V2 runs" (V2 needs a budget raise, front 6.3), not "money left" | front 6.3 |

## Review log (final round 6)

Cross-file consistency fixer, round 6. Sources: `reviews/final_r6_{cuda,rl,data,arch,qa,fresh-arch}.md` and their
`*_verdict.md` files. Only the CONFIRMED items assigned to this fixer are listed.

| Item | Severity | Change | Evidence |
|---|---|---|---|
| cuda-1 | blocking | A2 step 2: the fork-pool path holds only until D7 step 6, whose commit switches the ceiling CLI to in-process scoring; the A2 pool output is kept as the golden. I19 has two forms (pool before D7, golden match after). A.9.9 row 3 notes it | `serve.py:170-173`; `reward/score.py:114`; D7 step 6 guard |

## Review log (final round 7)

Cross-file consistency fixer, round 7. Sources: `reviews/final_r7_{rl,arch,qa,fresh-arch}.md` and their
`*_verdict.md` files. Only the CONFIRMED items assigned to this fixer are listed.

| Item | Severity | Change | Evidence |
|---|---|---|---|
| arch-4 | major | A3 `--rescore`: `--ckpt` and `--train` become optional; exactly one of `--ckpt` and `--rescore`; the I08 form stays; E6.4 keeps both A3 forms | `final_eval.py:48,51` (`required=True` today) |
| qa-2 | major | A.9.6 G02b and G03 run as `d_gpu_test.yaml` step 2 in V3 item 1 (D8 holds the argv); not in `D_GPU_TESTS` because G08 in the same file skips | D8; `policy.py:479, 486` |

## Review log (final round 8)

Cross-file consistency fixer, round 8. Sources: `reviews/final_r8_{cuda,data,arch,qa}.md` and their `*_verdict.md`
files. Only the CONFIRMED items assigned to this fixer are listed. Each fact was checked in the repo first.

| Item | Severity | Change | Evidence |
|---|---|---|---|
| data-1 | blocking | A.1 item 2, A2 step 4, A-G2 row, G04: val legacy reference 0.723 (whole val, item pooling); test 0.732; new found-stratum line 0.734 (val, 947) / 0.740 (test, 1,962) on R5's K-mean rule | re-run of `notes/r8data/legacy.py` (val 0.7232 / 0.7419, test 0.7319 / 0.7486) and `legacy_found.py` (found val 0.7339, test 0.7403) on 2026-11-09 |
| qa-1 | blocking | `a3_measure.yaml` step 2 builds if missing (assert, else `run_box.sh --build-only`, assert again); text and I20 updated; V1 adds no minutes | `policy.py:398-399`; `CMakeLists.txt:15`; front 6.3 V2 order |
| qa-2 | major | A6 step 8 no longer gives the baseline argv; C0 step 3 owns it; "17 steps" -> `--steps 20` with 17 read (A.3 diagram, `timeout` line, G06) | C0 files (B1 rows root, whole 20-step blocks) |
| qa-6 | major | G02b names its inputs: first 50 val rows by `item_id` with a target text; 8 model texts per item (targets cycled to 6 + 2 adversarial, seed 0) stored in the golden with a SHA-256; reads no `runs/` | A2 `item(texts=None)` gives the doctor side only; D8 payload has no saved answers |

## Review log (final round 9)

Cross-file consistency fixer, round 9. Sources: `reviews/final_r9_{cuda,rl,data,arch,qa,fresh-arch}.md` and their
`*_verdict.md` files. Only the CONFIRMED items assigned to this fixer are listed. Each fact was checked in the repo first.

| Item | Severity | Change | Evidence |
|---|---|---|---|
| data-1 | major | A1 step 3: `headline.py` reads both schemas (DOC6); new A1 step 5b `final_eval --resummarize DIR --out OUT` (Mac, no model; joins `group_id` from the store; schema 2; never writes DIR) for R5; U08b; A9 adds the `headline.json` regeneration with `ci_unit: patient` | `final_eval.py:41` (bare `ci95`); `ls runs/r5/eval`; keys of `runs/r5/eval/eval_test/items.jsonl` (6 headline metrics) |
| fresh-arch-1 | blocking | `a3_measure.yaml` step 2 uses a subshell `( ... )`, not a bash group; why-text names the `str.format` cause; I20 runs `policy.chain` on both A3 plans and a planted `{ true; }` case | `policy.py:434-435`; `policy.chain` on the doc block passes [measured 2026-11-09, impl/r9fix/] |

## Review log (final round 11)

Fixer for A, round 11. Sources: `reviews/final_r11_data.md` item 1 and `reviews/final_r11_data_verdict.md`
(CONFIRMED, major). Each fact was checked in the repo first. Skills: engineering:code-review and engineering:debug
(re-read `signal.draws`/`signal.ci` and recomputed the fixture CIs before pinning them); engineering:system-design
and engineering:architecture (one owner per rule: the phase section is the rule, the spec follows);
engineering:documentation (the key list is inline, not a pointer); nvidia-skills:nvidia-skill-finder and
huggingface-skills:huggingface-best: no NVIDIA product or model choice in scope, nothing changed.

| Item | Severity | Change | Evidence |
|---|---|---|---|
| data-1 | major | Header: "Where `A.spec.md` and this phase section differ, this phase section is the rule", plus the list of differences `phases/A.spec.md` must take before DOC1 step 0 copies it (§3.2 key list; §5 item 7 = A-G3 TOST on the 90 % patient CI, view identity, frozen-view hash; §5 item 8 = Gate A on the patient CI; revision line). A1: the pointer to A.spec §3.2 is replaced by the schema-2 key list (`paired.*` = `{mean, ci95_patient, ci95_nodule}`, `iou_found`/`miss_rate` with subset CIs, run keys incl. `view_content_sha256`). A4: `paired` keys named; `gate_a` reads `[gate_ci]` and raises on a null patient CI (new 9.8 rule 16). A-G4 names `paired.*.ci95_patient`; the nodule CI never decides. U10 and U11 pin `ci95_patient` and `ci95_nodule` against `signal.draws(..., key=groups.__getitem__)` on a 12-item, 3-patient fixture where the two lists differ. U51 adds the nodule-lower > 0 / patient-lower < 0 case (no change) and the null case (`ValueError`); sketch 8 follows. I22 checks the three `paired` keys and `gate_ci`; E01 reads `["mean"]`. U11's old `0.0125` corrected to 0.0375 (the mean of `[.1,.1,-.05,0]`). Not done here: the edit of `phases/A.spec.md` itself (this fixer edits only `A.md`); the header rule makes A.md win until that file is updated | `src/objdel/eval/signal.py:102-107` (`draws` has no `key` today; nodule only); `signal.ci` returns `[mean, lo, hi]` (`signal.py:110-119`); `docs/implementation2.0/notes/b_final/cluster_ci.json` (nodule width 0.069, patient 0.085, 635 items, 38 patients); A.spec.md §3.2 lines 69-71 and §5 items 7-8 (lines 326-329) still say nodule CI; fixture CIs recomputed 2026-11-09 (U10: nodule [0.5, 0.62], patient [0.525, 0.6]; U11: upper ends 0.1333 vs 0.15) |

## Review log (final round 12)

| Item | Severity | Change | Evidence |
|---|---|---|---|
| fresh-arch-1 | blocking | A3: one line, the record reaches the box through B1.6's keep rule; B's I50 checks it at $0 | B.md B1.6; `policy.py:379` |


# Phase B

## B.R. Review log

Source: the Phase B review (docs/implementation2.0/notes/reviews/B_phase_review.md) (27 items). Every item was checked against code before the change.

| # | Severity | Item | What I did |
|---|---|---|---|
| 1 | blocking | `prompts()` drops single-rater rows, steps shift | Confirmed [dataset.py:221, cache.py:103-105]. `pool` skips non-RL items and counts them (B1.1 step 1). `train_rl` asserts row count and ordered-id hash (B1.5). Tests added. |
| 2 | blocking | Shifted draw needs `side=S0`; rescaled "fallback" | Confirmed [augment.py:85-94, build_eval_views.py:57-62]. Exact calls written (B1.1 step 4). Kept-side windows are dropped and counted; DG6 = 0 and 0. S0 >= 0.75 test added. |
| 3 | blocking | Frozen hash breaks at the fingerprint change and across machines | Confirmed [build_eval_views.py:89-92, spec.py:42-50]. Content hash instead of manifest SHA (B1.6). Fingerprint split moved to new task B1.0, before any build. |
| 4 | major | Freeze contents and place left open | One file, `docs/reference/frozen_views.json`, with named fields; `final_eval` recomputes and stops (B1.6). |
| 5 | major | Neyman tilt down-weights the target cells, fails DG7, confounds B3 | Agreed. B3.8 removed with the reasons; `tilt` argument removed; B3 uses B1's untilted rows. |
| 6 | major | Batch-shape check tests the wrong quantity | Confirmed [grpo_trainer.py:1372-1377]. Check now uses `generation_batch_size // num_generations`, `steps_per_generation == gradient_accumulation_steps`, `num_iterations == 1` (B1.5). |
| 7 | major | New audit/coverage modules duplicate existing ones | Confirmed (`eval/coverage.py` used by preflight; `eval/gates.py` uses tuples). `eval/coverage.py` moves to `data/coverage.py` and is extended; gate module renamed `data/quality.py` (B1.3). |
| 8 | major | Window key and folder rule | Confirmed [dataset.py:102,110]. `window_record` passes the window id; separator "~", also because `signal.py:99` and `preflight.py:137` parse "/" (B1.1, B1.4). Unique-path test added. |
| 9 | major | B1.7 has no DICOM on the Mac | Confirmed (`data/raw` = `lidc_xml`, `riga`). XML-only $0 step first via `lidc.sop_z`; DICOM header download only on a hit, into its own folder (B1.7). |
| 10 | major | DG7-DG9 circular; val shapes and judges | Agreed and measured. Val split by patient now (B1.0; baselines in section 1 from valsplit.json). DG7-DG9 labelled implementation checks; DG16 report against the train-split law; "today" values against val-shape (also item 26). |
| 11 | major | Numerics confound B1; B1 blocked on B2 | Agreed. No passing arm: B1 runs on bf16-seqmask. A non-control winner is noted in ADR-0013/0015, bounded by the probe (B2.4). |
| 12 | major | RLOO without std changes the step size | Agreed. `adv_scale` constant in RewardV5, set from the median |A| ratio in B3.7; recorded in ADR-0014; test T13. |
| 13 | major | "Half the signal" unsupported | Reworded (section 1 item 5). Added B2.3b attribution. Note: at step 1 the adapter is fresh (LoRA B = 0 [g_probe.py:85-95]), so merge rounding is tested on the step-20 adapter, not at step 1 as the review proposed. |
| 14 | major | v4 does give gradient to all-miss groups | Confirmed [advantages.py:95-107, spec.py:28-34]. Section 1 item 6 reworded. B3.7 metric and gate B3 changed to Spearman(A, -d/r), v5 > 0.3 and > v4; variance-share report added. |
| 15 | major | fp16 vision overflow, unbuildable fallback, unproven patch | Agreed. Mac fp16 vision check before rental; mixed-dtype fallback dropped; dtype assertion after `__init__` with a CPU test (B2.2). |
| 16 | major | Found-image gate has composition bias | Agreed. Gate B3 uses images found by both runs, paired CI lower bound > -0.01. |
| 17 | minor | CI rules incomplete and unpaired | Agreed. All gates paired on val-select; B3 miss gate has a CI rule; B2 uses paired per-step differences (section 6.2, B2.3). |
| 18 | minor | Length-fair geo mask missing | Agreed. `sequence_geo_mask` in `is_ratio` (FastGRPO path only, since TRL raises on unknown modes [grpo_trainer.py:2823]); 4th probe arm bf16-geomask; SkyRL bounds checked on the web. |
| 19 | minor | alpha_loc claim false once g_area counts | Agreed. Claim reworded to "equal f5 and g_area" (B3.2, B3.3). B3.7 measures miss-beats-hit; `area_gate` switch if > 1 %. |
| 20 | minor | Old path deleted before gate B1 | Agreed. New task B1.8 deletes the old path and the no-schedule branch after gate B1. |
| 21 | minor | No task builds the guard | Agreed. B3.5b LocGuard with a synthetic-log test. |
| 22 | minor | Substitution breaks DG11 | Agreed. No substitution: `ValueError` naming (kind, cell); cell counts by largest remainder (B1.2); DG11 checks both. |
| 23 | minor | Optimistic cost; test read twice | Agreed. 27.5 s/step base case; rows root uploaded to HF before the rental; BudgetStop callback; test eval once on the final candidate (section 8). |
| 24 | minor | numpy polygonizer vs owner rule | Agreed. Exception in ADR-0015 with end condition (K1 parity); owner decision 8.3 item 2. |
| 25 | minor | Location kernel may be too short | Agreed. B3.7 reports d / r and sets `loc_width` or the heavy-tailed kernel; spill-distance alternative checked [gpucheck.cpp:1430-1440; full-frame EDT at cache.py:130-135]. |
| 26 | minor | DG9 baseline uses another reference | Agreed and measured: DG7-DG9 "today" values now against val-shape (0.448 / 0.157 / 0.040). |
| 27 | minor | IoU >= 0.66 not linked to the miss gate | Agreed. 0.66 is a target; decomposition added (+0.019 from the far-miss cut on val-select [estimate]); the CI rule is the gate. |

---

## B.R2. Review log (final review round 1)

Sources: the six final-review files `final_r1_{arch,cuda,data,fresh-arch,qa,rl}.md` (promoted to
`docs/implementation2.0/notes/reviews/`). Only items addressed to section B are handled here; items for other
sections are listed at the end. Every fact was checked against code or a measurement before the change:
- the mismatch-by-step table, grad norm 0.0018-0.022 (median 0.0048) and the KL spike (398,284 at step 746, grad
  norm 0.0041, median KL 0.32) [measured again: runs/rl/train_rl/rl.json];
- cell-law distances (train law vs val-select TV 0.063 / KS 0.037; val-shape vs val-select 0.169 / 0.099; today's
  off-centre rows vs train law 0.759 / 0.108 / 0.109) [measured: docs/implementation2.0/notes/b_final/trainlaw.json];
- R5 val-select CI width by item / nodule / patient 0.045 / 0.069 / 0.085 [measured: b_final/cluster_ci.json];
- 289 of 7,341 train items are disagreement items [measured: data/prepared/lidc/train.jsonl];
- `augment.shift` is uniform over the admissible rectangle [code: augment.py:58-71]; TRL OPSM and its use by
  FastGRPO's loss [code: grpo_trainer.py:3066-3089, 3190-3203; fast.py:567-571]; `sequence_mask` masks only ratios
  above 3 with `clip_min` None [code: grpo_trainer.py:2803-2818]; PEFT needs no reference model at any beta [code:
  grpo_trainer.py:1029-1036]; `get_soft_overlong_punishment` [code: trl/rewards/other_rewards.py:83]; v4 drop rule
  [code: advantages.py:103]; `train_rl` merges the final adapter and `make_trainer` adds a fresh LoRA [code:
  train_rl.py:58-59; g_probe.py:83-95]; Qwen3-VL-2B has 28 text layers [measured: HF config.json]; R5 val IoU 0.61819
  [measured: runs/r5/eval/eval_val/final_eval.json]; R5 has no eval on val rescaled / rotated [measured: ls
  runs/r5/eval]; generator files last changed in 734b1ec [measured: git log, read-only].

| # | Source item | Severity | What I did |
|---|---|---|---|
| 1 | data 1: target is a 40-patient sample | blocking | Target = train shifted law (B1.4, ADR-0015 with the measured distances). Gate B1 on the whole val; halves only for B3 (B.3.2, B.6.2). Diagram, DG7-DG9 "Today" (0.759 / 0.108 / 0.109) and DG16 re-measured. `--target-law train-shifted`; manifest records `target_law`. |
| 2 | data 2: gates resample nodules | blocking | Patient-level (`group_id`) paired bootstrap, BCa for the decision, cluster count shown (B.6.2); U50 rewritten; B3.7 Spearman CI by patient; B.8.4 A1 row. Width numbers re-measured. |
| 3 | data 3: "miss" has three meanings | blocking | Header defines miss (K-mean rule), found image, missed answer. B3.2, B3.5, B3.7 use "missed answer"; gates use the K-mean (B.6.2). |
| 4 | rl 1: the probe measures the mismatch where it does not exist | blocking | 4-arm probe removed. Counterfactual stats in the B1 run (B2.1); B2.3b on step-100 and final adapters; optional fp16 arm 60 steps, steps 41-60 (B2.3); Gate B2 and B2.4 rewritten; B.1 item 5 table; B.8.1 re-costed. |
| 5 | rl 2: drop rule removes far-miss groups | blocking | Drop rule on the span of each unweighted term incl. `loc` (B3.4); T14 test and sketch; B3.7 rule 2 drop share v5 <= v4 (Gate B3). |
| 6 | fresh-arch 1 (B part): paid order | blocking | B-S1 = B1 run first on control numerics; B-S2 waits for a budget decision (B.4 order, B.8.3 item 1 default (b')). |
| 7 | arch 4, fresh-arch 7, qa 7: ADR numbers collide | blocking / major | B keeps 0012 / 0013 / 0014 (two of the three proposed tables keep them). Header: the front-matter allocation table wins; rename in the same commit if it differs. |
| 8 | arch 1 (B part): phase order | blocking | B.4: "Every task path is pre-E; E's move tables define post-E paths." |
| 9 | fresh-arch 3, rl 5: `adv_scale` reason wrong under AdamW; RLOO = 8/7 mean-centring | major | `adv_scale` removed from RewardV5, B3.4, B3.5, tests (U01, U26, U29, P07, P08, T13). KL balance = `beta = 1e-3 / r` in the B3 config (B3.7 rule 6). ADR-0014 rewritten. |
| 10 | rl 4: KL weight has no owner | major | Owner decision B.8.3 item 5: B1 beta 1e-3 (clean data test vs R5, so my default differs from the reviewer's 0.0 for B1); B3 `1e-3 / r` or 0.0. B2.4 bullet. |
| 11 | rl 6: own mask where TRL has OPSM | major | ADR-0013 option E (TIS + OPSM, delta = q99.5), preset `bf16-tis-opsm`; option D built only under the counterfactual trigger; geo-mode tests made conditional. |
| 12 | rl 7: log counterfactual effects | major | `is_summary` counterfactual fields, `modified_token_share`, ESS (B2.1); U54, U62; B.1 item 5 corrects the "18 % masked" reading. |
| 13 | rl 8: attribution cannot separate dtype from engine | major | B2.3b adds (d) fp32 trainer and (e) vLLM teacher-forced bf16 / fp32; signed gap and 64-token bins; readings; G04. |
| 14 | rl 9: rule compares two units | major | B2.4 rule 3: higher ESS, tie <= 0.01 -> lower `modified_token_share`, then bf16-tis; U51 and sketch. |
| 15 | rl 10 (B part): vLLM LoRA serving | major | ADR-0013 option F, triggered by B2.3b (c) >> (b); B.8.3 item 9. The C5 text is Phase C's. |
| 16 | rl 11: LocGuard stops healthy runs | major | t > 3 on the location slope and IoU-slope upper bound < 0, on off-centre means, in two disjoint windows (B3.5b); U47, U58. |
| 17 | rl 12: Spearman gate undefined | major | Per-group Spearman, mean over groups, constant-A groups count 0, patient CI, paired v5 - v4 (B3.7 rule 1, Gate B3); U48. |
| 18 | rl 13: runaways get no gradient | major | Option O (soft overlong penalty as a term of R) as owner decision B.8.3 item 7 with C7; B3.7 rule 9 reports the share. Default off. |
| 19 | rl 17: fp16 check covers only the vision tower | major | `fp16_check` runs the whole model (28 text layers) on 32 prompts with saved answers; drop rule 0.5 x 65,504 or any inf/NaN (B2.2); U53, G03. |
| 20 | fresh-arch 4: shift law fitted to one law | major | Guard view added, but as the `edge` law: a "uniform" view would repeat the training law (`augment.shift` is already uniform). Hard floor in Gates B1 and B3; R5 baseline from A3; cost in B.8.1; B1.4 sentence "by choice". |
| 21 | data 9: trains on the test generator, deletes out-of-law checks | major | `--kinds` default all; val rescaled and rotated kept as report lines in Gate B1; Limitations line. |
| 22 | data 6: trained prefix never checked | major | 20-step blocks with exact counts (B1.2); DG7-DG9 on every block prefix; `planned_steps`; `train_rl` refuses partial blocks; P14, U61, sketch. |
| 23 | data 7: leakage gates miss the rows root | major | `rows_leakage` and hard gate DG1b, run in `build_rl_rows` and `preflight`; U55. |
| 24 | data 8: no label-noise audit | major | `agreement` report DG14b (needs A2) and owner decision B.8.3 item 4, default keep all; I31. |
| 25 | data 10: no power check | major | New task B1.9 `gate_power.py`; "owner may accept a near-miss" removed (B.6.2, B.7); E10. |
| 26 | data 11 (B part): ceiling claim | major | B.1 item 7 marked "[to be confirmed by A2]", bounds found images only; recompute in B.8.4. |
| 27 | data 12: z check covers 4 series | major | B1.7 checks every series ((a) shared z, (b) z distance, (c) 0396 fix); hard gate DG17; I17, I35. |
| 28 | qa 11: DG3 has no unit | major | DG3 = cross-split pairs, each listed and reviewed, 0 unreviewed; within-train rate informational; U24, I30. Also data 16 (test->val hit added). |
| 29 | qa 12 (B part): K-aggregation not named | major | B.6.2: every run gate uses the per-item mean over K = 8; Gate A method is a second line; U49. |
| 30 | qa 10 (B part): open behaviours | major | B.9.12 decided: DG12 = section 6.1 form; no padded windows (DG18 = 0, U60); invalid inputs raise. Owner confirms (B.8.3 item 10). |
| 31 | qa 6: one id means different things | major | Data gates renamed DG1..DG18; sessions B-S1 / B-S2; prefix rule in the header. |
| 32 | qa 14: synonyms | major | Header word table; "exam" replaced by "test view"; `K` for draws replaced by `draws`. |
| 33 | qa 15, cuda 14, fresh-arch 9 (B part): rates and step times | major | B.8.1 cites the front section 6 table; $0.641/h kept as the planning rate (two reviewers); 17 s column deleted; 20.0 s column only after G-SP1. |
| 34 | qa 16: `-I` rule has no way to run in tests | major | B.9.0: `data` marker, `make test-data` with `-I`, hygiene test. |
| 35 | qa 8 (B part): markers | major | B.9.0: one marker `gpu`. |
| 36 | qa 5, arch 14, cuda 13: scratchpad paths | major | Every `S[P]/` and absolute path replaced by `docs/implementation2.0/{notes,specs,seeds,kernels}/`; copy step listed under DOC in B.8.4. |
| 37 | fresh-arch 6 (B part): resume and data_revision unowned | major | Resume task in B1.5 (I33); `rows_revision` + manifest SHA pinned in `b1_rows.yaml` (B1.4 step 5, I34). |
| 38 | fresh-arch 8, arch 7, arch 8: layer placement | major / minor | `data/coverage.py` is not moved again; `advantages_v5` stays in `verify/reward/advantages.py`, `rl/` only calls it (B1.3, B3.4, B.8.4). |
| 39 | data 25: rows storage and pinning | minor | One place (owner decision B.8.3 item 6), `rows_revision`, manifest SHA asserted by `train_rl`, LIDC only. |
| 40 | data 15: DG13 samples 150 windows | minor | DG13 re-renders every window; I04, F06. |
| 41 | data 17: disagreement weight moves the mix | minor | Expected share 14 % stated; DG14 reports it; owner decision B.8.3 item 3 (4 or 1). |
| 42 | data 18: coverage on off-centre rows only | minor | DG9b report on the full mix (today 0.507 / 0.266 / 0.331). |
| 43 | data 19: size axis ignores rescale | minor | Apparent size `extent_px x S0 / side` (B1.1, ADR-0015); U59. |
| 44 | data 20: test already shaped the design | minor | B.1 item 3 uses val and train-law numbers; Limitations line in B.6.2. |
| 45 | data 21, qa 26: test view never hashed; order with A3 | minor | B1.6: A3's pre-rental hash must equal the frozen hash; `generator_commit`; "equal by code identity, not measured"; I05. |
| 46 | data 23: t bound over dependent steps | minor | Paired prompt bootstrap over the 640 prompts of steps 41-60 is the decision number; t bound second line (B2.3, I23). |
| 47 | data 24, fresh-arch 12, rl 24: extra-steps confound | minor | `--ckpt rl-full-3ep` stated; near-centre control and R5 plateau slope reported with Gate B1; ADR-0015 confound paragraph. |
| 48 | data 27: unbalanced val halves | minor | `val_half` = stratified `splits.deal`, stored in `docs/reference/val_halves.json`; B1.0 re-measures the half rows (not re-measured in this pass; the table says so). |
| 49 | fresh-arch 13: CIs exclude seed variance | minor | Stated in B.6.2; A3 two-seed noise floor; point gain must exceed it. |
| 50 | rl 21: KL spike row implies harm | minor | Reworded in B2.4 and B.7. |
| 51 | qa 22: sketch can pass without testing | minor | `hypothesis.event("raised")`, a feasible-pool test, keyword arguments everywhere. |
| 52 | qa 17: long sentences | major | All new and rewritten text uses one action per sentence; the merged test tables were not rewritten sentence by sentence (left for DOC's metric). Partial. |
| 53 | fresh-arch 2: decode tail has no task | major | Not a B task (C11 proposed). B.8.4 states that the decode tail and the runaway cut are Phase C's. |

Not handled here (addressed to other sections): data 4, 5, 13, 14, 22, 26; rl 3, 14, 15, 16, 18, 19, 20, 22, 23;
fresh-arch 5, 10, 11, 14; arch 2, 3, 5, 6, 9-13, 15-22; cuda 1-12, 15-25; qa 1-4, 9, 13, 18-21, 23-25, 27, 28.

---

## B.R3. Review log (final review round 2)

Sources: the six final-review files `final_r2_{arch,cuda,data,fresh-arch,qa,rl}.md` (promoted to
`docs/implementation2.0/notes/reviews/`). Only items addressed to section B are handled here. Rows B.R2 2, 7 and 49
describe round-1 choices; this log supersedes them (A1 interval rule, renumbered ADRs, split-half noise line).
Every fact was checked against code or a fresh measurement before the change (B.0 revision 4 lists them):
- merge rounding re-run on `runs/rl/train_rl/trainer/checkpoint-{50,100,300,500,770}` (bf16 erased share 0.692,
  0.643, 0.557, 0.512, 0.464; relative error 0.605, 0.531, 0.426, 0.379, 0.339) [measured: merge_rounding.py];
- runaway shares by step range, `lr` 5e-5, no eval key in 771 log entries [measured: runs/rl/train_rl/rl.json];
- held-out IoU_tc 0.511 / 0.634 / 0.659 [measured: heldout.json files; full_rl_2026-11-08.md:37];
- 5 of 999 fetch sidecars with `slices != annotated_sops` or `matched_by_z > 0` [measured: data/lidc_slices];
- apparent-axis coverage and null quantiles [measured: final_r2_data_check.py on geom.json];
- adapter size 17,432,576 parameters [measured: checkpoint-100 tensor shapes];
- `load_inplace` in vLLM's LoRA docs [web: https://docs.vllm.ai/en/latest/features/lora/].

| # | Source item | Severity | What I did |
|---|---|---|---|
| 1 | fresh-arch 2, data 1, qa 2, rl 10 (B parts): V1 does not make the Gate B1 inputs | blocking | B.8.1: B-S1 = the V1 item list: B1 on val shifted and val edge, R5 on val edge (never skipped), A3's val part under a skip rule (nsys first, then R5-val). Totals 205.6 min = $2.20 at 27.5 s, 186.4 = $1.99 at 23.65 s; items 1-5 alone $2.01 / $1.81. B2.3b (a)-(d) moved to the Mac; A3's test re-sampling leaves V1. B.8.3 item 11 (cap $2.20 or $2.35). B.8.4 A3 row lists only what A3 makes. The front 6.2 table, its "A3 is not needed for Gate B1" sentence and 99_back decision 2 must copy this list (front's edit). |
| 2 | data 1, data 5, fresh-arch 2, qa 2: "two sampling seeds" noise floor does not exist | blocking / major | B.6.2 "Noise floor": `noise_floor.mean_sd` from A3 step 7 is reported beside each gain; it is K = 8 sampling noise already inside the paired CI; training-seed variance is not measured. B.7 row rewritten (second-seed run = owner decision, not in $3). B1.9 no longer reads A3. |
| 3 | qa 3, data 2: two interval rules for the same gates; no BCa task | blocking / major | B.6.2 uses the A1 gate CI: percentile patient interval at >= 50 patients (Gate B1), `signal.cluster_t` below 50 (Gate B3); subset rules use per-patient means of the patients in the subset, with the count printed. B.8.4 A1 row, U49, U50 changed. No BCa anywhere outside the round-1 log. |
| 4 | arch 3, fresh-arch 5, qa 4, rl 9: ADR numbers collide with the front table | blocking / major | All 45 references renumbered in one pass, highest first in effect (0014 -> 0015, 0013 -> 0014, 0012 -> 0013), including combined forms (`ADR-0013/0014`). Header states the mapping and that 0012 is DOC's. The DOC4 drift check is DOC's. |
| 5 | arch 1 (B part): two archive folders for old plans | blocking | B1.5, I22 and the B.8.4 CLEAN row use `archive/phase2_runs/configs/`; B1.5 creates its README (Was / Replaced by / Commit); the live RL plan is `b1_rows.yaml`. E's own references are E's edit. |
| 6 | rl 1: merge rounding is measured and dominant, but tested late and conditional | blocking | Measurement re-run (B.1 item 5 table). New task B2.0 ($0, Mac, all 28 layers x 7 modules, planted 2x2 test U64). B2.3b (a)-(d) on the full run's checkpoints at $0. ADR-0013 option F first-class (vLLM `enable_lora`, `LoRARequest(..., load_inplace=True)`, adapter-only `save_pretrained`), preset `lora-seqmask`, F plumbing in B2.2 (U70, G05), A7 patch points. B2.4 rule 0. C5's push under F is noted in B.8.4 (C's edit). |
| 7 | rl 2: drop rule lost `legal` | blocking | B3.4 `flat` includes `legal`, with the 6-legal / 2-illegal example; T15 in B3.6, B.9.0, Gate B3; U63; P08 updated. Cited line corrected to `advantages.py:98-103`. |
| 8 | rl 3: masking truncated rows biases toward truncation | blocking | B.1 item 8 (derivation and the measured runaway table). B3.5 and `b3_v5.yaml`: `mask_truncated_completions=False`; B1 keeps True. Clipped-ratio slope in every run report; Gate B3 report line. ADR-0014, B.7, B.8.3 item 7; option O dropped. C11 note in B.8.4 (C's edit). I39. |
| 9 | rl 4: preset rule cannot be computed and is degenerate | blocking | B2.1 logs `cf/ess_C`, `cf/ess_B`, `cf/ess_E@δ` at token level, reported, never maximised. B2.4 rule 2 default E, rule 3 vetoes, rule 0 first. Gate B2 ESS line replaced. U51, U54 and the sketch rewritten. |
| 10 | rl 5: teacher-forced vLLM scoring is not the IS-ratio quantity | major | B2.3b way (f) generation-time logps is the vLLM line; (e)-bf16 kept as "prefill vs decode"; (e)-fp32 only in a fresh process with `VLLM_ATTENTION_BACKEND` named, no retry. G04. |
| 11 | rl 6: no check of the all-miss gradient share | major | B3.7 rule 10 (S), Gate B3 pre-run rule S_v5 >= 0.5 S_v4, B.7 row, U65. |
| 12 | rl 7: option O penalises correct long outlines | major | Option O dropped (item 8). B3.7 rule 9 states the constraint `L - cache >= 1,181 + 32` if it is ever revived. |
| 13 | rl 8: `c_topo` has a lower bound only | major | B3.7 rule 8 reports the jump and the topology share of the A variance; Gate B3 share <= 0.5; owner choice recorded in ADR-0014; U66. |
| 14 | data 3, qa 6: paired baseline has no proof of the same images; wrong hash key | major | B1.0 pre-rebuild hash of `data/eval_views/lidc-val-shifted` and a third frozen entry; Gate B1 pairs with A3's R5 val items (R5 items second line, or fallback only with an equal pre-rebuild hash); view-identity rule with `r_sh`; B1.6 key `view_content_sha256`; DG15; I05, I36. CLEAN10 note in B.8.4. A-G3 is A's edit. |
| 15 | data 4: confound as large as the effect; cited series does not exist | major | ADR-0015 confound with the held-out points (+0.025 over 520 steps; up to +0.014 for 300 [estimate]); hard difference-in-differences rule in Gate B1 and in B1.9; slope bullet replaced by the two held-out points; U49, sketch. |
| 16 | data 6: B1.9 inputs come after the spend | major | B1.9 noise SD `sqrt(mu_i (1 - mu_i) / 8)` (exact upper bound) and a noise-0 run; rules fixed from the bound run; A3 `item_sd` only for the record after B-S1; E10. |
| 17 | data 7: test hash frozen before the z check; no rule for test hits | major | B.4 step 1 and B1.6: B1.7 runs before A3 step 1. B1.7 test-hit rule (`known_label_issues`, sensitivity line, owner decision after the hash). B.7 row; I17, I37. A3 step 1 and front section 1 need the same order (their edits). |
| 18 | qa 5: gate ids "B-DG1..16" vs 21 gates | major | Header: 21 gates (DG1..DG18 plus 1b, 9b, 14b); B-DG prefix for other phases. The front lines are the front's edit. |
| 19 | qa 9: tests that accept either outcome | major | U12, U14, U27, U32, U39, U46 assert one outcome; B.9.12 item 5 lists the new rules. |
| 20 | qa 12 (B part): geom.json cited but not copied | major | B.1 cites geom.json as "not kept in the repo: 9.1 MB" and names the files DOC1 copies. |
| 21 | rl 12, cuda 15: quantiles cannot use min/max merges | minor | B2.1 device buffer of per-sequence gaps, one `torch.cat` + `torch.quantile` in `log`; U67. |
| 22 | rl 13: `NUMERICS` sketch not valid Python | minor | Named parts `_SEQMASK`, `_TIS`, `_FP16` first; U38. |
| 23 | rl 14: fp16 arm confounded by the IS mode | minor | Preset `fp16-seqmask` (control IS mode); `b2_fp16.yaml`, Gate B2, U40, I27, G05. |
| 24 | rl 15: OPSM delta from B1 is stale | minor | B2.4 rule 4: delta 0.5 until step 20, then q99.5 of the B3 run's own first 20 steps, frozen; U68. |
| 25 | rl 16: LocGuard reads an unwritten log | minor | B3.5 logs `kind/<k>/loc`; B3.5b names the keys; U44. |
| 26 | rl 18: learning rate not stated | minor | `--lr 5e-5` in `b1_rows.yaml`, `b2_fp16.yaml`, `b3_v5.yaml`; ADR-0013 records it; CLEAN9 sets one value (B.8.4); I22, I39. |
| 27 | arch 20: `quality.py` seed imports `transformers` | minor | B1.3 import rule and import test (U69). |
| 28 | arch 4, arch 8 (B parts) | major | B.8.4 Layers row names the three new import edges for E0's ignores and keeps `rows_revision` + `rows_manifest_sha256` as the rows pin. E0 and E6 text is E's edit. |
| 29 | data 8: free sidecar evidence ignored | minor | B1.7 check (d); 5 of 999 sidecars hit [measured]; step 2 takes them; DG17 "Today"; I17. |
| 30 | data 9: "Today" values on the rejected size axis | minor | DG7/DG8/DG9b and B.1 item 3 on the apparent axis (0.557 / 0.253 / 0.109; 0.375 / 0.350 / 0.331), raw values kept for reference; U59 checks callers use `apparent_px`. |
| 31 | data 10: rescaled kind trade-off unstated | minor | ADR-0015 trade-off with the measured shares; DG14 reports item reuse; owner decision B.8.3 item 3b. |
| 32 | data 11: DG16 has no noise reference | minor | DG16 prints null quantiles; B.6.1 row shows TV p95 0.100 and KS p95 0.052. The val/train worst-ratio orientation differs between the two scripts, so only TV and KS (symmetric) carry the null. |
| 33 | data 12: two offset definitions | minor | B.1 item 1 uses geom.json (0.262, 36.7 %, 807 of 2,200); test row 158 patients; `view.offset_frac` in every view record; I06. |
| 34 | data 13: val halves reuse the split seed | minor | B1.0 seed `"val_half:0"`, `deal_weighted` for `dis=1`; U02 and the B1.0 test. |
| 35 | data 14: test can be read for two candidates | minor | B.4 step 6, B3.9, B.8.1: test is read once after the 4B gate; B-S2 ends with the val evals. |
| 36 | data 15: edge fallback and padding undefined | minor | B1.6 edge draw rule (64 directions, policy-frame fallback counted), padded share reported, Gate B1 edge rule on non-padded non-fallback items; I32. |
| 37 | data 16: code identity names only generator files | minor | `render_files` and the R5 code commit in `frozen_views.json`; I38. |
| 38 | data 17 (B1.7 part) | minor | B1.7 (c) names 139.xml and position 1420 with the TCIA source. The DOC11 rows are DOC's edit. |
| 39 | data 18: val-select can shape B1 | minor | ADR-0015, B.6.2 and B.7: changes after Gate B1 use val-shape and logs only. |
| 40 | qa 19: 1,181 vs 1,182 | minor | Header defines outline length (no end token): 1,181, or 1,182 with it. DG12 keeps `max + 1`. |
| 41 | qa 20: test-source paths | minor | `docs/implementation2.0/tests/B_tests_{1,2}.md`. |
| 42 | qa 24: one R5 step time | minor | B.8.1 fp16 row uses 23.65. |
| 43 | qa 25: long sentences | minor | New and rewritten text uses short sentences; old test-table cells were not split. Partial. |

Not handled here (addressed to other sections): fresh-arch 1, 3, 4, 6-9; arch 2, 5-7, 9-19; cuda 1-14, 16+;
data 17 (DOC rows); qa 1, 7, 8, 10, 11, 13-18, 21-23, 26; rl 11, 17. Cross-phase follow-ups that this revision
needs: front 6.2 and 99_back decision 2 copy the B-S1 list (item 1); A3 step 1 after B1.7, A3 step 7 keys
`item_sd` / `mean_sd`, and the A-G3 view-identity rule (items 2, 14, 17); C5 and C11 under F and unmasked rows
(items 6, 8); E0 ignores and E6 rows fields (item 28); DOC4 ADR drift check (item 4).

---

## B.R4. Review log (final review round 3)

Sources: `reviews/final_r3_{data,fresh-arch,rl,qa,arch,cuda}.md`, items addressed to B. Facts checked before any
change: front 6.2 V1 table (5 items; 195.1 min = $2.08 at 27.5 s, 175.9 = $1.88; item 5 skipped 188.2 = $2.01)
[code: 00_front.md]; R5 preflight stamp `runs/preflight/r1_qwen3vl2b.json` data hash equals today's recomputed hash,
stamped 6 min before the R5 launch [measured]; val-half seed search (rev-4 rule: 454/630 items, 22/22 disagreement
items; first passing s = 78: 533/551, 24/20) [measured: S[P]/impl/b_r4/valhalf.py on the store, `-I`]; R5 val IoU
0.618 item-weighted vs 0.647 patient-weighted [measured]; held-out set 128 items, 51 patients [measured];
sidecar hits by split and the duplicate-z val patients [measured]; `b1_pixels.py:39-41, 67` Pearson rule [code];
`RunArgs(**raw)` rejects a `steps:` plan [code: preflight.py, train_sft.py]; RL AdamW comes from GRPOConfig, only
train_sft.py:384-388 builds its own optimizer [code]; TRL sleep level 2, OPSM sign and dr_grpo divisor [code: .venv
trl grpo_trainer.py, vllm_generation.py]; `LoRARequest.load_inplace` defaults to False and forces a reload when
True [web: github.com/vllm-project/vllm main vllm/lora/request.py, fetched 2026-11-08].

| # | Review item | Severity | Where handled |
|---|---|---|---|
| 1 | data 1, qa 2, qa 16, fresh-arch 1: B-S1 is not the front V1 list | blocking | B.8.1 table is a copy of front 6.2 items 1-5 (minutes and $), no nsys line, no 20.0 s column for V1; DOC4 drift check. |
| 2 | data 2, qa 3: two primary baselines; view identity has no input when item 5 is skipped | blocking | B.6.2 "Paired baseline" rule fixed before any B1 number: A3-val items primary; R5 items with the B1.0 stamp hash chain ("hash only, correlation not run") when item 5 was skipped; broken chain waits for V2 item 3b. Item 5 skip rule on session spend + $0.08 <= cap (B.8.1). U49, I46. |
| 3 | data 3: val-half rule does not balance items | major | B1.0 seed search with criteria (a) items within 4 %, (b) disagreement items within 4; s = 78; `val_halves.json` keeps s and counts. U02. |
| 4 | data 4: `cluster_t` estimates another number | major | B.6.2 and B1.9: CR1 interval of the item-weighted mean difference, formula and citation; B.8.4 A1 row asks A1 step 5 for it. U50. |
| 5 | data 5: DiD cannot separate extra training at the ceiling | major | ADR-0015 and B.6.2: A2 headroom informative check, "B1 vs R5" wording; checkpoint-250 option (B.8.3 item 12, $0.07). |
| 6 | data 6: val label hits | major | B1.7, B.6.2, B.7 row, DG17: val hits stay in the frozen view, `known_label_issues`, excluded in both arms. I17, I45. |
| 7 | data 7: val already chose R5 | minor | B.6.2 Limitations and sensitivity line without the 128 `r5_selection_items`. |
| 8 | data 8: B1.7 test hit (LIDC-IDRI-0659) | minor | B1.7 table, B.7 row, B.8.3 item 13 (due before A3 step 1). |
| 9 | data 10: near-duplicate rule | minor | B1.3 and DG3: mean-subtracted unit-norm dot product, -1,000 HU pad. U74. |
| 10 | data 12 (B side): bundle ships test views | minor | B.8.1 bundle note (no `lidc-test-*` folder; `bundle.py:53-56`). |
| 11 | fresh-arch 2, qa 6: BudgetStop on the wrong clock | major | B1.5: `BudgetStop(rate, cap, session_start)`, `run.sh` exports; missing start is an argparse error; B.7 row; diagram. U45, U46, I42, sketch. |
| 12 | fresh-arch 3, rl 7, rl 8: F picked for B3; cross-run timing; B3 changes four things | major | ADR-0013 rule (b): B3 keeps `bf16-seqmask`, F only after Gate B3; F check in-process ABAB ratio in the first 15 steps; B3.9, B3.10 credit line, Gate B2/B3, B.8.1 rows, B.8.3 item 9. U51. |
| 13 | fresh-arch 4: no decision tree after Gate B1 | major | New B.4.1 table; B1.8 runs only on its two rows; B.7 row; B.8.3 item 14. I43. |
| 14 | fresh-arch 5: fixed model load per eval | major | B.8.1 note (about 205 s per process; multi-view `final_eval` saves about 6.8 min); B.8.4 A3 row. |
| 15 | fresh-arch 6: no held-out curve | minor | B1.5 optional `--heldout` (owner item 15, default off), B.6.2 report line. I44. |
| 16 | rl 1, rl 2, rl 3: F base push under sleep, `load_inplace`, layout check | blocking | ADR-0013 F text and B2.2 F steps 2-4: unmerged base push in mode B, `add_lora` once per sync, `load_inplace=False` on generate, `reset_prefix_cache`, packed-layout equality after a sleep/wake cycle. U70, G05. |
| 17 | rl 4: OPSM sign | major | B2.1 `seq_kl` definition (sampling minus trainer), buffer stores it. U54. |
| 18 | rl 5: live delta unwatched | major | B2.4 rule 4: `on_step_end` at step 20, `opsm_delta.json`, A7 attribute, 2 % guard; B2.1 masked-share logs. U68, U75. |
| 19 | rl 6: Adam eps with smaller gradients | major | ADR-0014, B2.4 KL bullet, B3.4, B3.7 rule 6, B3.9: `adam_epsilon = 1e-8 / r`; B.7 row. U71. |
| 20 | rl 10: "no IS mode can undo" | minor | B2 intro and B.7 row: systematic gap; IS modes trade bias for variance. |
| 21 | rl 13: LocGuard on independent noise only | minor | B3.5b Newey-West lag 5; AR(1) simulation test. U72. |
| 22 | rl 14: truncated-row token weight | minor | B3.5 `completions/truncated_grad_share`, 10 % owner switch. U73. |
| 23 | rl 15: erased-share denominator | minor | B2.0: non-zero delta elements only. U64. |
| 24 | qa 1, arch 1 (B side): no preflight stamp for plan files | blocking | B1.3 and B1.5: `b1_launch.yaml` (RunArgs + `rows_plan`), "ROWS root" check, manifest SHA in the data hash; B.8.4 infra row. I19, I22. |
| 25 | qa 4: `v1-b1` lacks B2.1 | blocking | B.4: tag holds B1.0-B1.9, B2.0, B2.1 (logging only). |
| 26 | qa 7: adapters V1 -> V3 | major | B1.5 `pull.sh` then `hf_push --rl` to `b1-rows`; B2.3b reads the B1 final adapter. I41. |
| 27 | qa 8: B2.4 timing and code edits in V3 | major | B2.4 header (rules by time), no code edit inside a rental, `--numerics`, missing-input rule. U51. |
| 28 | qa 9 (B side): `pytest -m gpu` collects every phase | major | B.9.5 rental/command table; "never without a file list"; `-I`. |
| 29 | arch 2: rows root and manifest SHA onto the box | major | `cmd_fetch_hf --rows`, `train_rl --rows-revision --rows-manifest-sha256`. I34, I40. |
| 30 | qa 13: two outline-length limits | minor | Header note; B3.7 rule 9 uses 1,182 + 32. |
| 31 | qa 17: "ledger" for session spend | minor | Header definition; B.6.2, B.7, B.8.1, B.8.3 item 11 say session spend. |
| 32 | qa 18: option letters reuse phase letters | minor | Preset names in rules, gates, B.7; letters A-F only in the ADR table. |
| 33 | qa 19: heading order | minor | B.9.12 moved before B.9.13. |
| 34 | cuda 11: "Phase D GPU CI" | minor | B3.1 and G07 name R-D2's `make gpu-check` sanitize step. |
| 35 | qa 21: long sentences | minor | New and rewritten text uses short sentences; old test cells not split. Partial. |

Not handled here (addressed to other sections only): arch 3-17; cuda 1-10, 12, 13; qa 5, 10-12, 14, 15, 20; rl 9,
11, 12; data 9, 11. Cross-phase follow-ups this revision needs: the front must adopt the item 5 rule (session spend
+ $0.08 <= cap), `b1_launch.yaml` with `rows_plan` as the V1 launch config, the launch line and ADR-0013 rule (b)
(fresh-arch 3: V2 item 2b is then not needed); A3: `a3_measure.yaml` holds only val-shifted in V1, a multi-view
`final_eval`, and CR1 `cluster_t` in A1 step 5; 99_back decision 1 quotes B-S2 at $1.94; C: G-EQ5 uses the B2.2
packed-layout check after a sleep/wake cycle and C5's mode-B base push.

## Review log (consistency pass, round 4 prep)

| Item | What changed |
|---|---|
| One source (r3 fresh-arch 1, qa 2) | B.8.1 holds no rental table; it maps B-S1 = V1 (front 6.2) and B-S2 = V3 items 0, 2a, 2b, 4 (front 6.3). B.8.3 items 1, 8, 9, 11, 12 point to the front |
| Ledger-based skip rule (r3 fresh-arch 1, qa 17) | Removed: no V1 item is skippable; Gate B1 report lines no longer say "when the ledger allows" |
| `b1_rows.yaml` | Steps: `env_check`, GPU smoke G02/G08, `train_rl --budget-usd 2.05`, one `final_eval` with two views; R5 on val edge moved to `a3_measure.yaml` (front item 4) |
| `BudgetStop` (r3 fresh-arch 2, qa 6) | Adds `fixed_usd`; `run.sh` exports `OBJDEL_FIXED_USD` and a rate with the bandwidth term, so it equals `policy.accrued` (U46, I42 updated) |
| Adapter push | `push_hf.sh` runs `hf_push --rl` on the box before teardown (`policy.py:511-524`), not on the Mac |
| Gate B1 fallback | "item 5 skipped / V2 item 3b" -> "R5 process gave no output / V2 item 1b" (B1.0, B.6.2, U49, I36) |
| B2.3b box part | 15 min, B's own estimate, priced as V3 item 2a |

## Review log (final round 4)

Cross-file consistency fixer, round 4. Sources: `reviews/final_r4_{cuda,rl,qa,arch,fresh-arch}.md` and their verdicts.

| Item | Severity | Change | Evidence |
|---|---|---|---|
| cuda-2, qa-2 | blocking | B1.5 `b1_rows.yaml` step 0 = CUDA build; B3.9 `b3_v5.yaml` step 0 the same; I22 asserts both | `policy.py:399, 432-433`; `CMakeLists.txt:15` |
| cuda-3 | major | B1 run `--steps 260` (front 6.2); B1.2 note: 8,320 rows; B1-run windows "steps 50-260"; Gate B1 confound 260 steps, +0.013, (260/520); B.8.3 item 1; B3.9 same steps as B1 | front 6.2 "Why 260 steps" |
| qa-3, fresh-arch-1 | blocking | G08 split: G08a (V1: build assert + G02 no-sync run; today's fields via `test_gpu_scoring.py` in step 0) and G08b (after B3.1: G01 sums, fp16 dtype check). Smoke step by node ids with `-m, pytest`. B.9.0 rule: box test names start with `test_<id>_`. Sketch header renamed to `tests/test_train_fast.py`. Coverage map updated | `policy.py:484-490`; `grep sum_x` = 0 hits |
| qa-8 | major | `train_rl` saves `final_adapter/` before the merge; `hf_push --rl` uploads `final/` and `final_adapter/`; B2.3b and G04 load `final_adapter/` on `rl-full-3ep`; I41 checks both | `train_rl.py:58-59` |
| rl-5 | major | B2.2 F step 4: the F run measures its comparator before step 1 (merged push vs `sync_lora` of the B1 final adapter, both scored by way (b)); pass <= 0.5x. B2.4 rule 0 names it; B2.3b box part only when rule 0 does not fire | B2.3b condition; B2.0 seed 0.531 |
| rl-6 | major | B3.9: B3 keeps `max_completion_length` 1,536 and `--runaway-cap 0`; L and the cap apply after Gate B3; I22 checks `b3_v5.yaml` | `configs.py:73`; `grpo_trainer.py:3286` |
| fresh-arch-3 | major | Held-out curve: `num_generations_eval=8`, `--heldout-items 64` (512 answers per eval kept) | `configs.py:177-181`; `rewards.py:132-135` |
| qa-4 (side effect) | blocking | I22 now lists the five B1.5 steps with `--budget-usd 2.05`, `--steps 260` and no R5 step | B1.5 |
| rl-3 | major | B2.1: `log` also puts each step's per-sequence `seq_kl`, advantage and valid-token count (256 x 3 float32) into the one host read and appends them to `<out>/seq_gaps.f32` (+ `seq_gaps.json`; 0.8 MB for 260 steps), which `pull.sh` brings to the Mac. B2.4 rule 4: the callback keeps its own steps 1-20 device list and a 20-step ring buffer for the guard, apart from the list that `log` clears. Rules 3 and 4 pricing read the B1 run's `seq_gaps.f32`. U68 asserts the q99.5 over all 20 x 256 values with `logging_steps=1`; new R01b | `g_probe.py:88` (`logging_steps` 1); `policy.py:80` (pull excludes) |
| data-1, qa-4 | blocking | I22 already listed the five B1.5 steps at $2.05 with no R5 step (row above). Added: `--model qwen3-vl-2b` and `--train {data}/rows/lidc-train-b1` in step 3, `--train {data}/prepared` in step 4, and the expected `--budget-usd` computed as `V1_APPROVED_USD - TAIL_USD` (2.40 - 0.35, front 6.4 rule 5) from two named constants, not a literal | front 6.2 plan files; front 6.4 rule 5 |
| data-2 | major | B1.9 noise model: both arms get noise with variance `c x mu (1 - mu) / 8`; c = 0.077 from one measured pair (labelled as one same-seed pair, not proof of independence), c = 1 on both arms printed as the stress line, plus noise 0; rules fixed from c = 0.077; "< 0.5: change a rule" reads the c = 0.077 line. "Exact upper bound" removed. B.6.2 "Noise floor" gets the 128-item retest (SD 0.0569, per-arm variance 0.00162 vs bound 0.0211). E10 updated, new E10b | `r4data/retest.py` re-run 2026-11-08: var diff 0.00324, two-arm bound 0.04216, ratio 0.0768; `trl/generation/vllm_generation.py:366` (seed = `process_index`) |
| data-3 | major | B1.2 step 3: weighted draw replaced by a per-block disagreement quota per (kind, cell) slot set, `q = w p / (1 + (w - 1) p)` (14.1 % at w = 4, R5's share), largest-remainder rounded; step 4: a disagreement window may repeat up to 4 times only when the cell runs out, else `ValueError` naming (kind, cell, part). B1.4 step 2 prefers `draws` with no reuse and records `disagreement_max_reuse`. DG11 (hard, exact) now includes the quota counts in every block, so the gate count stays 21; DG14 text "q per cell; 14.1 % overall". P06 and the scheduler sketch assert the quota per block | `dataset.py:81, 227` (4x repeat); 289 / 7,341 train items; 7,341 - 289 + 4 x 289 = 8,208; `r4data/disagree.py` (block drift at draws 1-8) |
| data-5 | major | B1.5 step 3 adds `--model qwen3-vl-2b`; step 4 adds `--train {data}/prepared` (the `a3_measure.yaml` root). New I47 parses every `-m` step of `b1_rows.yaml` and `a3_measure.yaml` with the module's own argparse parser on the Mac | `train_rl.py:23` (`--model` required); `final_eval.py:51` (`--train` required); `policy.py:486` (only `-I` is added) |
| qa-5 | major | B1.3 change 2: `preflight.local_rows_root` maps the plan's box path `{data}/rows/<name>` to the local `data/views/<name>`; any other form raises `ValueError` (a failed "ROWS root" check). I19 adds the mapping cases | B1.4 `--out data/views/lidc-train-b1`; I40 box path `/workspace/data/lidc/rows/lidc-train-b1`; `policy.py:433` (`{data}` = `/workspace/data/lidc`) |

## Review log (final round 5)

Cross-file consistency fixer, round 5. Sources: `reviews/final_r5_{cuda,rl,data,arch,qa,fresh-arch}.md` and their
`*_verdict.md` files. Only the CONFIRMED items assigned to this fixer are listed.

| Item | Severity | Change | Evidence |
|---|---|---|---|
| data-1, qa-2 | blocking | B.4.1: table removed; it points to front 1.1. B.6.2: the correlation is not a Gate B1 hard rule; new hard rule "frozen-view hash" (each B1 and A3 `view_content_sha256` equals `frozen_views.json`, DG15). U49: correlation removed from the hard-rule list; planted hash mismatch fails with "stop". B.9.8 sketch: `view_hash_ok` replaces `view_corr`; a low correlation still passes and drops the R5 line | front 1.1; A-G3 |
| data-2 | blocking | B1.5 "Rows-root download": no `--subdir`; `hf_pull --out /workspace/data/lidc --allow "rows/ROOT/**"` (new `--allow` flag), so the root lands at `{data}/rows/lidc-train-b1/`; the same ssh call checks `manifest.json` and its SHA-256. I40 asserts `local_dir` and the box check; new I40b fakes the real layout and resolves the root | `hf_pull.py:18-22`; `.venv huggingface_hub/_snapshot_download.py:152`; `policy.py:433` |
| qa-1 | blocking | B1.5 step 3: `--numerics bf16-seqmask` removed (the tag holds no B2.2); `rl.json` writes the literal. I22 split: I22a (V1 configs, on the tag, no `--numerics` key, `--shared-prefix`) and I22b (`b2_fp16.yaml`, `b3_v5.yaml`, in the commits that create them). B.9.0: a config test lands with its file. I47 runs on the tag tree. Coverage rows follow | `train_rl.py:23-36`; front order rule 2 |
| qa-4 | major | B1.5 "Checkpoint folder name": `fetch_hf.sh <run> REV [--as NAME]` to `{rd}/<NAME>` (default `ckpt`); V3 fetches `b1-rows --as b1`; B2.3b reads `{rd}/b1/final_adapter`. I40b covers `--as` and refuses `/` or `..` | `policy.py:534-535` |
| rl-1 | major | B2.2 F step 2.3: adapter folder in PEFT form, `adapter_model.safetensors` from `peft.get_peft_model_state_dict` (no `.default.` in keys) plus `adapter_config.json`; ADR-0013 F text points there; U70 checks the key pattern and the config | checkpoint-100 header: 392 keys, 0 with `.default.`, `r` 16, `lora_alpha` 32 [measured 2026-11-08]; `.venv peft/utils/save_and_load.py:349` |
| rl-2 | major | B2.4 rule 1: the "s/step at most 1.05 x the B1 run's" clause removed; s/step of both runs is a report line (`b2_summary.py`), with rule 0's reason (two rentals; one vLLM engine holds one dtype, so no in-process bar). U51 text and the B.9.8 sketch: rule 1 has only the half-logp-diff boundary; a new case with 2x the B1 s/step still picks `fp16-seqmask` | B2.4 rule 0 text; front 6.3 (fp16 arm = V3 item 2b, B1 = V1) |
| rl-3 | major | B2.1 new "Train mode only" bullet: `set_is_ratio` appends to `seq_gaps` and the rule-4 callback lists only when `self.model.training`; `log` computes the quantiles and writes `seq_gaps.f32` only in train mode; eval `is_stats` stay in the eval metrics. New R01c (one eval between two train steps leaves 2 x 256 rows and the quantiles unchanged); coverage row | `.venv trl/trainer/grpo_trainer.py:3416-3420, 3425`; `fast.py:524-538, 567-577, 609` |
| rl-5 | major | B2.2 new F step 5: `train_rl --f-check` (default on under `lora-seqmask` and `fp16-seqmask`) runs the G05 / G-EQ5 tensor checks, the comparator, the load count, the ABAB timing and NaN/inf in the training process, writes `<out>/f_check.json` and raises on a failure. B2.4 rule 0 points there. G05 rows: the box part is the node-id record check `test_g05_f_check_record` after `train_rl` (no `-k`). New U76 (fake engine); coverage rows | `configs.py:192` (`vllm_gpu_memory_utilization` 0.35); B1.5 step 2 (`policy.py:490`, empty selection exits 5) |
| data-3 | major | B1.9: joint P(pass) = Gate B1 rules (1) val IoU, (2) far-miss point + CI, (3) offset < 0.15, (4) difference in differences; edge view, legal share and frozen-view hash printed as "not simulated" (no R5 edge values before V1). New report-only Gate B3 power line on val-select (B3 arm cuts 25 / 35 / 50 % of the B1 arm's remaining misses, same c, `cluster_t`), never in the joint value or the owner's rule-change decision. B.7 row; new E10c; coverage row | `ls runs/r5/eval` = `eval_val`, `eval_test` [measured 2026-11-08]; B.6.2 Gate B1 and Gate B3 lists |

## Review log (final round 6)

Cross-file consistency fixer, round 6. Sources: `reviews/final_r6_{cuda,rl,data,arch,qa,fresh-arch}.md` and their
`*_verdict.md` files. Only the CONFIRMED items assigned to this fixer are listed.

| Item | Severity | Change | Evidence |
|---|---|---|---|
| rl-1 | blocking | B2.2 F step 2.2: a sync never pushes the base in mode A; the comparator is the one exception and restores the base before step 1 in both modes (new step 4b.5, with a host hash check) | `fast.py:716-717` (merged values overwrite the base) |
| rl-2 | major | New B2.2 F step 4b: base = the run's `--ckpt`; B1 loads as adapter `"b1"` in slot 2 (`sync_lora(adapter_name, slot)`); KV wake before each comparator generation; direct `llm.generate` with `None` for arm (i); order (ii) then (i); restore and `sleep(level=2)` in mode B. B2.4 rule 0 ABAB timing runs inside the comparator. U70 and U76 gain the comparator cases | `peft save_and_load.py:148-151`; `fast.py:703-707`; `vllm_generation.py:376-377, 718-719` |
| data-2 | major | B1.8: grep-based test list (CT-33, CT-35 included) instead of line numbers; later phases state their form with and without B1.8 | `train_rl.py:26, 47`; `serve.py:107-121` |
| data-3 | major | B1.4 step 5: one `upload_folder` call, `rows_revision` = `CommitInfo.oid`, `files_count` and `files_sha256` in the manifest. New B1.3 change 2b ("ROWS revision" Mac check). B1.5 box check compares the file count. I19 gains the fake-Hub cases | `hf_api.py:5779-5784`; 14,684 files in `data/views/lidc-train-shifted` [measured]; `dataset.py:176-180` |
| cuda-3, qa-4 | major | B3.9: `b3_launch.yaml` (output_dir `b3-v5`, `rows_plan: b3_v5.yaml`) and `--out {rd}/b3-v5`; I22b checks both | `policy.py:511-524` |
| qa-5 | major | B.9.0: fixture `tiny_loop_vllm_fake` (use_vllm True with a fake `VLLMGeneration`, IS correction on, `sequence_mask`); R01, R01b, R01c run on it and size by `rows_per_step`, not 256 | `test_p2_scripts.py:392-393`; `fast.py:568, 576-578`; `vllm_generation.py:111` |
| arch-5 | major | Review-log scratchpad path written with brackets | check 14 scan: 0 hits |
| fresh-arch-1 | blocking | B3.10 and B.8.4 point to the front "End state under (b')" | front section 1 |
| data-1 | major | B1.1: new bullet, `Window.disagreement` is per item (slice), copied from the source item's prepared field (`dataset._prepare` rule), equal across the item's kinds and draws, never spread over the nodule; a nodule-level flag would give 1,394 items (19.0 %) and a quota of about 48 % instead of 14.1 %. B1.2: p is the per-slice share. U11 rewritten (per item; planted two-slice nodule with one True and one False; 289 flagged on the real records) | `dataset.py:114` (`len(rater_answers(caches, names)) > 1`, per record); 1,247 nodules, 100 with a disagreement slice, 100 of those mixed, 1,394 items, 289 / 7,341 flagged [measured: `python -I` on data/prepared/lidc/train.jsonl, 2026-11-08] |
| qa-2 | blocking | B2.4 rule 0 split into a start choice (`start_preset`, reads only the B2.0 error: >= 0.2 starts the first post-B3 run under `lora-seqmask` with `--f-check`) and a restart choice (`restart_preset`, reads `f_check.json`: passed keeps it; failed restarts under `fp16-seqmask` if rule 1 passed, else the rules 2-3 preset). "Missing inputs": only the fp16 arm can be missing; the F check is never a pre-run input. ADR-0013 records both presets before the run. U51 text and the B.9.8 sketch: `decide_numerics` replaced by the two functions; `(0.531, no F, no fp16)` now gives `lora-seqmask`; the F cases move to `restart_preset`; veto cases call `start_preset` | `train_rl.py:23-36` (no `--numerics`, no `--f-check` today; B2.2 F step 5 adds them); front 6.3 ("the F check is the first 15 steps of the first run after Gate B3") unchanged |

## Review log (final round 7)

Cross-file consistency fixer, round 7. Sources: `reviews/final_r7_{rl,arch,qa,fresh-arch,data}.md` and their
`*_verdict.md` files. Only the CONFIRMED items assigned to this fixer are listed.

Skills used for rows rl-2, rl-3, data-1, qa-1: nvidia-skills:nvidia-skill-finder (no catalog skill covers TRL/vLLM
comparator rules or val halves; nothing installed); engineering:code-review (read `splits.py`, `build_data.py`,
`configs.sampling()`; checked the new sketches parse and run against a reference rule); engineering:debug
(re-measured the d = 0 floor, the drift medians and the seed-0 vs s = 78 halves before editing); engineering:system-design
(fixed what each comparator number measures: floor-subtracted gap, fixed-length timing); engineering:architecture
(one owner per rule: `passed` is set by `train_rl`, `restart_preset` only reads it; B1.0 owns the halves rule, U02 and
B.9.8 follow it); engineering:documentation (one field list for `f_check.json`, in B2.2 F step 5);
huggingface-skills:huggingface-best (no model choice in scope). Web: vLLM `SamplingParams` docs for `ignore_eos`,
`min_tokens`, `seed`.

| Item | Severity | Change | Evidence |
|---|---|---|---|
| rl-1 | blocking | B2.2 F step 4b.5: the restore order is push base, host hash against `p.data`, `remove_lora(2)`, adapters, `reset_prefix_cache()`, then `sleep(level=2)` last in mode B; the hash never reads a sleeping engine. U76 gains the mode-B sleep-order case (fake tensors raise when read asleep) | `.venv trl/generation/vllm_generation.py:378` ("Sleep level 2 discards the weights"), `:497-499` |
| fresh-arch-1 | blocking | B1.8 says no plan trains without a schedule after it; C's arms use the rows root. 6.1 step 4's HF revision check runs before V2 too | front 6.3; C0 files |
| fresh-arch-4 | major | B1.5 deletes the `--train-views` path (`grpo_rows_views`, `VIEW_SEP`, `_tagged`, the flag, `g_probe.py:90`) in the commit that archives `rl_shifted.yaml`; the no-schedule branch calls `grpo_rows(a.train, ...)`. B1.8 deletes only the no-schedule branch, runs on rows 1-3, and on row 4 after the owner's written review (unless it drops the rows path). U42, I43 and the archive row follow | `train_rl.py:26,47,51`; `dataset.py:240-250`; `serve.py:108-123`; grep: only `rl_shifted.yaml:26` names `--train-views` |
| rl-2 | major | B2.2 F step 4 (iii): a floor arm (0) generates the batch on base W with `lora_request=None`; the trainer scores it with the adapter disabled; f0 = mean \|f_base - b_base\|. Pass: (f_lora gap - f0) <= 0.5 x (f_merged gap - f0); merged gap - f0 < 0.002 records `informative: false` and the output-level check passes. `f_check.json` gains `comparator` {`f_merged`, `f_lora`, `f0`, `ratio`, `informative`, `ok`}. Comparator cost now about 5 min [estimate]. B2.4 rule 0, Gate B2 and U76 (pass 0.0096/0.0176/0.0096, fail at ratio 0.55, `informative: false` at 0.0019) follow | re-measured `sampling/sampling_logp_difference/mean`: step 1 = 0.009563 (runs/r5/base/rl.json), 0.009515 (runs/r5/shared_rl.json), 0.009556 (runs/rl/train_rl/rl.json); full-run medians 0.01764 (240-280), 0.02361 (101-300), 0.01853 (301-500) |
| rl-3 | major | B2.2 F step 4b.4 and B2.4 rule 0: both timed arms use `SamplingParams(max_tokens=512, min_tokens=512, ignore_eos=True, seed=s)`; every timed row must decode 512 tokens; one untimed warm-up per arm; A B A B for seeds 0, 1, 2; `lora_over_base` = median of the 6 ratios. Arms (0), (i), (ii) keep the normal sampling. Order: (ii), (0), timing, then (i) (its merged push changes the base). U76 gains the planted length-difference case (1,536 vs 800 tokens does not move the ratio) and a 511-token row that fails. `min_tokens` added beside `ignore_eos` so no stop id can end a row early | docs/implementation2.0/notes/B_gpu_audit.md:40 (5.61 s + 4.3 ms x longest answer, r = 0.92); vLLM `SamplingParams` docs (`ignore_eos`, `min_tokens`, `seed`) [web: docs.vllm.ai api sampling_params, 2026-11-08] |
| data-1 | major | B1.0 names the `Patient` form for the halves: (item count, disagreement items) per patient in the frozen val view, not the split's XML form; adds `splits.val_half_table(records, max_s)` (pure search) and stores `patient_form` in `val_halves.json`; states that seed 0 fails both rules. B.9.8: the seed-0 sketch is replaced by `test_val_half_table_is_the_first_balanced_seed` (first s that passes, every smaller s fails, table equals `deal(..., seed=f"val_half:{s}")` per stratum) and `test_val_half_reads_the_committed_table` (`KeyError` on `LIDC-IDRI-0001` and `lidc/LIDC-IDRI-0026`); the hash-seed test prints `val_half('LIDC-IDRI-0026')` and asserts the table value, so a `KeyError` fails. U02 and the fixture list (`tmp_val_records`, built so that s = 0 fails) follow | re-measured `python -I` on `data/eval_views/lidc-val-shifted/lidc/val.jsonl`: seed 0 gives 40 / 38 patients, 567 / 517 items, 14 / 30 disagreement items (fails (a) and (b)); `"val_half:<s>"` first passes at s = 78: 533 / 551, 24 / 20; strata 41 / 23 / 2 / 12; ids unprefixed, `LIDC-IDRI-0001` absent from val; `splits.py:53-101`; `build_data.py:123-127`. The new `_halves` sketch body re-run on the real val gives s = 78 |
| qa-1 | major | B2.4 rule 0: `restart_preset(cf, f_check, fp16, b1)` reads `f_check["preset"]` (start preset) and `f_check["passed"]` only; passed keeps the start preset; failed from `lora-seqmask` with rule 1 passed -> `fp16-seqmask`; failed otherwise (including start `fp16-seqmask`) -> the rules 2-3 preset; a restart never returns the failed preset. B2.2 F step 5 lists the restore base hash as a check, names every `f_check.json` field and defines `passed` as the AND of the preset's flags. U51 moves the per-check boundaries to U76 and adds the comparator-only, base-hash, trusted-flag and fp16-start cases; U76 asserts `passed` = AND. B.9.8: `F(...)` is a full `f_check.json` record; 7 restart cases plus a veto case; each case asserts a failed record never returns its own preset | Sketch cases run against a reference `restart_preset` written from the new rule text: all 7 cases and the veto case pass. Front 6 unchanged (the F run is after Gate B3, outside the funded plan) |

## Review log (final round 8)

Cross-file consistency fixer, round 8. Sources: `reviews/final_r8_{cuda,data,arch,qa}.md` and their `*_verdict.md`
files. Only the CONFIRMED items assigned to this fixer are listed. Each fact was checked in the repo first.

| Item | Severity | Change | Evidence |
|---|---|---|---|
| data-1 | blocking | ADR-0015 "Confound": near items compared with their own legacy ceiling 0.731 (206 items; whole val 0.723), not 0.734 | re-run of `notes/r8data/near_ceiling.py`, 2026-11-09 |
| cuda-2, qa-3 | major | B3.9: `b3_v5.yaml` needs no step or flag for D-G14/15/17; D7 step 5 code writes them | D.md D7 step 5 |
| qa-4 | major | B3.9: `b3_v5.yaml` step 1 = `test_g08b_smoke` (`-m, pytest`), before `train_rl`; input = first 1,000 crops of `{data}/checker_eval/masks.npz`; G08b split into `test_g08b_smoke` and `test_g08b_fp16_dtype`; runner row names V3 items and the `-m` argv; I22b asserts both and a 1-item collect | `policy.py:486`; `bundle.py:66-70` |
| qa-5 | major | B2.3: `b2_fp16.yaml` step 0 = `b1_rows.yaml` step 0 (CUDA build), step 1 = `test_g08b_fp16_dtype`; I22b asserts both | `policy.py:398-399`; `CMakeLists.txt:15`; V3 item 1 skippable |

## Review log (final round 9)

Cross-file consistency fixer, round 9. Sources: `reviews/final_r9_{cuda,rl,data,arch,qa,fresh-arch}.md` and their
`*_verdict.md` files. Only the CONFIRMED items assigned to this fixer are listed. Each fact was checked in the repo first.

| Item | Severity | Change | Evidence |
|---|---|---|---|
| rl-1 | major | B2.2 F step 4 defines "the 32 x 8 batch" (rows 0-31 of the `--train` root, 8 answers each, 256, seed 0, TRL's processor-expanded ids); 4b.3: comparator and wrapper share `fast.vllm_generate_checked` (one pad per image, then the expansion check); `prompt_expansion` flag; U70 and U76 updated | `fast.py:407-415, 618-622, 680-688`; `.venv trl/generation/vllm_generation.py:720-743`; `configs.py:10` |
| rl-2 | major | B2.4 "Eval of a run trained under `lora-seqmask`" (the one rule): `final_eval --ckpt <base> --adapter <out>/final_adapter` with `enable_lora=True, max_lora_rank=16` and `LoRARequest("eval", 1, dir)` through `vllm_generate_checked`; `final/` only for a merged report line; ADR-0013 and B1.5 point to it; new I48; coverage row | `train_rl.py:58-59`; `final_eval.py:48`; grep `enable_lora` 0 hits |
| qa-1 | major | B3.9 step 1 (G08b smoke) is a `bash -c` step with `OBJDEL_DATA_ROOT={data}`; the test reads through `box_input`, `BOX_INPUTS` in `tests/test_train_fast.py`; runner row and I22b assert it | `policy.py:379, 486, 492` |
| qa-2 | major | G01: one file `tests/test_scorer_sums.py`; box reference = `np.nonzero` sums + Mac golden `tests/goldens/b_g01.json` (`make_b_goldens.py`, GPU02 form); G08b uses its `windows_1000` entry; B.9.8 drops `cpu_build`/`cuda_build` and the two-build sketch; `test_g08b_fp16_dtype` reads `$OBJDEL_RUN_DIR/ckpt`; runner rows split G01/G06 | `cuda/gpucheck_py.cpp:113`; `Makefile:41-45`; A.9.6 G02b |
| qa-6 | major | B3.9 names the B3 val eval outs `{rd}/eval_b3_val_shifted` and `{rd}/eval_b3_val_edge` (B1 form); I22b checks them; D's S13 reads the first | B1.5 step 4 form |
| fresh-arch-1 | blocking | I47 gains a format check: `policy.chain` on every `steps:` plan (bash steps too) returns without error; planted `{ true; }` fails; today all 9 `steps:` plans pass | `policy.py:434-435`; measured 2026-11-09 |

## Review log (final round 10)

Cross-file consistency fixer, round 10. Sources: `reviews/final_r10_{cuda,data,qa,fresh-arch}.md` and their
`*_verdict.md` files. Only the CONFIRMED items assigned to this fixer are listed. Each fact was checked in the repo first.

| Item | Severity | Change | Evidence |
|---|---|---|---|
| data-1 | major | B.6.2 "Paired baseline" gains the "Edge baseline" bullet (front 6.2 rule); the guard-law line points to it; U49 adds the two cases (items 4 and 5 missing: "open: edge baseline missing", row 4 with `raise_planned=False`; item 4 present, item 5 missing: edge pairs with item 4) | `runs/r5/eval/` [measured: ls, 2026-11-09]; front 6.2 |
| qa-4, fresh-arch-2 | major | B.4 step 4 and B.8.1: B-S2 = V3 items 0 and 4; 2a, 2b, 5 are optional V3 rows (front 6.3 one definition). B3.9 states `--budget-usd` from the front "V3 amount rule" ($2.05 / $2.73); I22b checks it against named constants | front 6.3, 6.4 rule 5 |
| qa-1 | major | B.9.5 G01 gains the "Tables rule": every object m (random mask, crop, or a polygon's even-odd raster from `check_polygons`) is scored against `b_tables(m) = ExpertCache.build(stack of m rolled +1 x, +1 y, -1 x)` (three raters in m's own frame; m is not a rater, so `iou_tc` is not trivially 1, a change from the review's fix, whose stack held m and gave crop `iou_tc` = 1.0 everywhere); crops pass as stored (`bits` one byte per pixel, `hw`, `start`; the review called them bit-packed, they are not), one handle per crop; polygons scored once per grid G (7, 1000, 2^30). One copy of inputs and rule: `b_sets` / `b_tables` / `random_masks` in `make_b_goldens.py`, imported by the test; golden keys `random`, `polygons_<G>`, `windows_1000`, each with the input SHA-256 and whole-call field hashes. G07, G08b, P12 and the B3.1 tests bullet point to the rule; the B.9.8 fixture line names `score.py:67-77` (not `test_cuda_checker.py`, which never calls `score_masks`); the sketch is rewritten with the real 4-argument call | `cuda/gpucheck_py.cpp:72-87, 104-106, 213-224`; `cuda/gpucheck.cpp:1174-1185`; `reward/score.py:67-77`; `reward/cache.py:218-266`; `geometry/shape.py:88, 143-147`; Mac CPU build runs (`notes/r10qa1/`): 0 build errors on 203 masks, 1,000 crops, 10,000 polygons; sketch `b_sets` 100.9 s; extracted sketches compile |

## Review log (final round 11)

| Item | Severity | Change | Evidence |
|---|---|---|---|
| rl-1 | major | B2.2 F step 2.1 (asleep = TRL's flag, cleared after the wake), 4b.3 (wake only through `sync_lora`/`push_weights`), 4b.5 (set the flag True after `sleep(level=2)`); U70 (one `add_lora` per RL step through TRL's real path) and U76 (flag True after the restore; fake raises on `load_weights` while asleep). C5 "Sleep state" is the one rule | `.venv trl/generation/vllm_generation.py:376-379, 500-503, 602-603, 771-773`; `grpo_trainer.py:1203, 1838-1841` |
| qa-2 | major | B3.9 step 2 = G-RC0 by node id; I22b asserts it and its collect count | `policy.py:432-433`; B.9.5 node-id rule |
| fresh-arch-1 | major | B3.9 step 3 = B1 re-eval on `{rd}/b1/final` (`eval_b1v3_val_*`), before `train_rl`; Gate B3 pairs with it, V1's eval is a report line; I22b asserts the step; B3.9 `--budget-usd` points to front 6.3 without copying values | E4 step 3 (`VllmSampler` in `final_eval`); D7 eval data flow; B1.5 `hf_push --rl` uploads `final/` |
| fresh-arch-2 | major | B1.5 `b1_rows.yaml` step 3 soft, step 3b `train_rl --merge-last` (merge last checkpoint, `partial: true`); Resume note (pull job); B.7 crash row; I43 partial case; I47 parses soft steps; new I49; test map | `policy.py:484-490`; `train_rl.py:55-60` |
| fresh-arch-3 | major | New B1.4a: `rings.polygonize` over today's `check_masks`, callers via `checker/__init__.py` and `certify.py`, numpy file kept only as oracle until D3, hash gate on the Mac before B1.4. B.3.3 row, ADR-0015 bullet, B.8.3 item 2 (removed, number kept), K1 note and D dependency row changed; test map row | K1/results/cross_numpy_m3.txt; bench_cpu_m3.txt (4,833.9 ms for 37,093 masks); issues_solved_2026-11-07.md:37 |

## Review log (final round 12)

| Item | Severity | Change | Evidence |
|---|---|---|---|
| qa-1 | major | B1.4a moves `checker/polygonize.py` to `tests/oracles/polygonize_oracle.py` now (D3 keeps only the bit-quad move) and points the reference callers at it (`compare.py` :21/167, `eval.py`, `xarch_hashes.py`, `cross_numpy.py`, `test_checker.py`, `test_targets.py`); `test_rings.py` gains the oracle guard (AST import check, function not shadowed, `mask_reference` unchanged when the seam is patched) | `compare.py:21, 167`; `test_cuda_checker.py:38`; `test_cuda_checker_gpu.py:50`; `checker/__init__.py:6` |
| qa-2 | major | B3.9: G-RC0 is the last, soft step of `b3_v5.yaml` (`~`, after the B3 val eval); order 0 build, 1 smoke, 2 B1 re-eval, 3 `train_rl`, 4 B3 val eval, 5 G-RC0; I22b asserts the leading `~`; Gate B3 cites B3.9 step 2 | `policy.py:484-487` (soft step); `test_train_fast.py:24-35` (tiny model in memory) |
| fresh-arch-1 | blocking | B1.6: `final_eval.frozen_entry` reads `parents[2]/docs/reference/frozen_views.json`; B1.6 adds the one `code_files()` keep rule for that file (D1 step 0 widens it); new I50 `test_frozen_views_ship` (temp repo, real repo, child `python -I` on a `code_files()` copy finds `lidc-val-edge`); E05 and test map point to I50 | `policy.py:379` (drops `docs/`); `bundle.py:49-75`; `final_eval.py` (no ROOT today) |
| data-1 | major | One exclusion list `docs/reference/train_exclusions.json` (`{item_id: {reason, date}}`), written by B1.7 drops and DG3 "drop train item". B1.1 `pool(..., exclude=frozenset())` skips excluded items first and counts them under `dropped["excluded"]` (five keys; U09 fixed from "3 keys"). B1.4 `build_rl_rows --exclude` (required, `{}` when empty) passes the keys to `pool` and `rows_leakage`, so the target law is computed without them; the manifest `schedule` records `exclusions_sha256` and `excluded_ids` (the box `preflight` reads these, because the code sync drops `docs/`). B1.7: resolutions "dropped" (every store item of the series in the list; today 27 items in 4 series) or "mapped" (fix, re-fetch for a (d) hit, rebuild store with `build_data` and prepared root with `prepare_rl`, val/test parquet SHA-256 unchanged, re-run `check_series_z` lists no hit); the report records the train-store SHA-256 and a stale report fails DG17. DG1b gains "no `source_item_id` in the list" and the Mac SHA check. B.4: B1.7 runs before B1.4. Tests: U09, U24, U55, I17, I24 (rows root holds no window of a dropped series), I30, E04; gate map DG17 adds U09, I24; DOC row lists the file | `reward/cache.py:103-105`; `vast/policy.py:379` (drops `docs/`); `phase1/fetch_lidc_pixels.py:5, 115-118` (`--refetch-incomplete`, sidecar rewrite); `phase1/build_data.py:3`, `phase1/prepare_rl.py:3`; 27 train items [measured: `data/store/lidc/train.parquet`, `python -I`, 2026-11-09: 380786 -> 15, 228428 -> 5, 179144 -> 4, 295628 -> 3] |

## Review log (final round 14)

| Item | Severity | Change | Evidence |
|---|---|---|---|
| rl-2 | major | B3.9: step 3 (`train_rl`) soft with `--save-steps 50`; new hard step 3b `train_rl --merge-last --out {rd}/b3-v5`; order line 0, 1, 2, 3 (soft), 3b, 4, 5; pull job pointer to front 6.3. I22b asserts the leading `~`, `--save-steps 50` and step 3b, and the val eval follows step 3b. I49 names B3.9. New failure-table row "B3 training crashes before step 260" | `policy.py:484-490` (`&&` join, `~` soft step); `train_rl.py:55-60` (`final/` only after `trainer.train`); front 6.2 pull-job lines |
| rl-3 | major | B2.2: new `train_rl --opsm-start-delta X` (only with `bf16-tis-opsm`, argparse error otherwise or for X <= 0) sets `off_policy_mask_threshold` after the preset; the only start-delta path (`--grpo` holding the key is refused); the run plan writes ADR-0013's `start_delta`; `rl.json` records `opsm_start_delta`; preset comment says 0.5 is the default start delta. B2.4 rule 3 names its q99.9 `start_delta`; rule 4: "Until step 20, delta = `start_delta`" (0.5, or the rule-3 value), `start_preset` / `restart_preset` return it as a number (`None` for other presets), never the string "q99.9"; the "Reason" line now says it is about the frozen step-20 value. U38 keeps 0.5 as the preset default only; U45 adds the flag cases; U51 text and the B.9.8 sketch assert `start_delta == approx(np.quantile(SEQ_KL, 0.999))` at `opsm05` 0.021, 0.5 at 0.020, the same value from `restart_preset`, `None` for `lora-seqmask` (`raise_delta_to` removed); U68 plants 0.7 and asserts 0.7 for steps 1-20, 0.5 without the flag | `.venv trl/trainer/grpo_trainer.py:882` (threshold copied once at init), `:3190-3203` (read at each loss call); `.venv trl/trainer/grpo_config.py:348-352` ("e.g., 0.5"); `experiments/phase2/train_rl.py:23-36` (no `--numerics`, `--grpo` JSON at :36); `configs.py:195` (`| overrides` last). Front 6 unchanged (runs after Gate B3 are outside the funded plan) |

## Review log (final round 15)

| Item | Severity | Change | Evidence |
|---|---|---|---|
| data-2 | major | B1.6 keep regex is now `^docs/reference/(frozen_views|val_halves)\.json$`, so the B1.5 held-out curve finds `val_half` on the box. B1.0 states where `val_half` reads the table (`parents[3]`, then `repo_root()` from E2.5). B1.5 names the need. I50 gains (a)/(b) for both files and (d): a child `python -I` on the `code_files()` copy prints `val_half("LIDC-IDRI-0026")`; a copy without the file fails. Test map B1.5 row adds I50 | `policy.py:379-380` (drop `^(data|docs)/`, no keep today); `docs/reference/` holds no JSON today [measured: ls]; `policy.py:398-399` (editable install at `/workspace/src`) |

## Review log (final round 16)

| Item | Severity | Change | Evidence |
|---|---|---|---|
| cuda-1 | major | B.9.8 sketch: `b_sets` gains `n_masks=203, n_polys=None, n_crops=1000` (defaults = G01's full sets, so the golden inputs do not change); new `test_g07_sanitize_sums` (markers `gpu` and `sanitizer`) calls `b_sets(..., 32, 64, 32)` and checks only the sums; `test_g01_golden` keeps `gpu` only. B.9.5 G07 row defines the subset and why (about 141,000 engine calls for the full G01 versus about 3,300). B3.1 Tests bullet and the G07 runner row name the one test. D8 and D-G5 changed to match | `src/objdel/reward/cache.py:168-197` (12 `edt` calls per build); `experiments/cuda_checker/compare.py:249` (`fuzz_polygons`); B.9.5 G01 set sizes 203, 170, 9,678, 152, 1,000 [measured: notes/r10qa1/run_sketch.py] |


# Phase C

## C.R. Review log (revision 2, answers `docs/implementation2.0/reviews/C_phase_review.md`)

| # | Sev. | Item | What I did |
|---|---|---|---|
| 1 | blocking | Reference pass runs inside TRL before C3 swaps `pixel_values` | Accepted, verified [code: grpo_trainer.py:2827-2860]. New mechanism: `RowMeta` + `prefill` + `EncoderCache.encode` in `FastGRPO._generate`; `FastForward` drops `pixel_values` and gathers by `self.row_image_ids` (from `self.meta` in the reference pass, from `host_meta` in our `_compute_loss` override, because TRL passes only named kwargs [code: grpo_trainer.py:3151-3169]). C3 depends on C4 steps 1-2. Tower-call counter test added to G-EQ4. Also found the split/unsplit patch is unneeded [code: utils.py:1169-1170] and removed it. |
| 2 | major | Pass 2 still tokenizes 256 rows on the CPU | Accepted with a measurement: pass-2 text part 0.065 s, renders 0.003 s on the Mac [measured now]. Added `PromptMemo` (C2 step 5, adopt if >= 50 ms on the box) and a Mac equality test. Rejected the "key by message parts, rendering is the cost" half: rendering 256 prompts costs 0.003 s, so the rendered text stays the key. §3.3 rows fixed. |
| 3 | major | CUDA `image_grid_thw` causes 256 syncs | Accepted, verified [code: processing_qwen3_vl.py:78]. Grid is always CPU int64; Mac test and a box `count_syncs()`=0 check added. |
| 4 | major | Delta push saves ms, adds risk | Accepted, verified the parameter count (1.41 B touched, config measured). C5 is now one call, full push, every step; `lora_touched_names`, `full=`, `--full-push-every`, the sleeping-delta error and the delta clause of G-EQ5 are gone. Micro-benchmark with a >= 50 ms re-add rule; vLLM `weight_transfer` recorded as rejected. |
| 5 | major | vLLM's own CPU preprocessing, H2D and duplicate encode missing | Accepted. §3.1 and §3.3 rows added; C0 step 6a measures it; new C6b investigation gated on vLLM DeepStack embed support [unverified] and a >= 0.3 s measurement. |
| 6 | major | Sync counter undercounts and includes vLLM | Accepted. `simplefilter("always")`; counting suspended in `generate` with a separate vLLM count; G-SP4 restated (< 200 and <= 25 % of C0's trainer-side count); nsys validation. |
| 7 | major | G-SP5 unreachable as scoped | Accepted. G-SP5 is trainer-side; vLLM H2D reported; C6b gate has "vLLM H2D < 10 MB". |
| 8 | major | R4 awake evidence not reconciled; wake double count | Accepted. Verified 8.83 vs 9.22 s from the R4 table. Mode A is now 1.0-1.7 s; C0 splits the weights and KV wake timers; A/B logs residual decode and longest answer; new risk row. |
| 9 | major | Baseline arm distorted by nsys | Accepted. 17-step arms, means steps 3-15, `--profile-steps 16-17` on both arms; §8.2 updated. |
| 10 | major | GPU-active share undefined; fallback cannot produce it | Accepted. Definition and `experiments/phase2/gpu_active.py` (nsys CSV or Chrome trace) with a Mac test; nsys baked into the image (C0 step 0, exact Debian repo path [unverified], checked at image build). |
| 11 | major | "50 recorded batches" unreachable on the box | Accepted. Separate 15-step check arm (checks every step, `--record-batches`), Mac batches from Phase A answers with the real processor; gates now ">= 50, >= 15 from the box". Not put into the timed fast arm, so check work never inflates G-SP1. |
| 12 | major | G-SP1 counts C6 | Accepted. G-SP1 from C1-C5 only, after items 8 and 14: A <= 20.0 s (fail > 21.2), B <= 20.7 s (fail > 21.9); rounding to 0.1 s (revision 1's "round up to 0.5 s" was not applied consistently). |
| 13 | major | 1.61 GB pixel tensor still built on GPU then discarded | Accepted. `ImageMemo.per_image` returns the 201 MB 32-image tensor under `--encoder-once`; verified TRL's InternVL check reads only `ndim == 4` [code: grpo_trainer.py:2662-2665]. No split/unsplit patch needed (item 1). G-MEM measures both ways. |
| 14 | minor | G-SP2 mode B too loose | Accepted. G-SP2 gates `PushStats.seconds` <= 0.2 s and mode B `sync_s` <= push + measured weights wake + 0.1 s; §8.1 recomputed with the measured p before the A/B. |
| 15 | minor | Silent CPU fallback | Accepted. Raise unless `--allow-cpu-pixels`; self-check on 8 images (centred, shifted, rescaled, grey, digit). |
| 16 | minor | Copying `_compute_loss` vs a library | Partly accepted. The 5-line `append` monkeypatch is impossible: `.item()` runs inside TRL's local `global_masked_mean` before `append` [code: grpo_trainer.py:3381-3387]. Kept the >= 0.1 s rule; added an upstream issue draft for the owner to file and an ADR record. |
| 17 | minor | Goal says C1-C6 bit-equal | Accepted (§1). |
| 18 | minor | G-EQ6 labels wrong | Accepted, verified: step-1 `reward` -2.98e-8, IoU_tc 0.5467, legal 0.8984, length 117.89, grad norm 0.003264, logp diff 0.009515. Note: the R5 base arm's grad norm is 0.003266, so shared-prefix vs base differ at the 4th digit; both Phase C arms use `--shared-prefix`. |
| 19 | minor | C5 numbering and cross-reference | Accepted; C5 rewritten with steps 1-6 and the rejected options moved into the ADR paragraph. |
| 20 | minor | No model-change risk row | Accepted (§7), plus a header note on Qwen3-VL-specific tasks. |
| 21 | minor | Grey planes, pinned buffer per step | Partly accepted. Persistent pinned buffer reused after a CUDA event: accepted. Grey planes: only for mode-`L` sources (RIGA); LIDC PNGs are RGB [measured now], so their H2D stays 25.2 MB. nvImageCodec claim kept with its web tag. |
| 22 | minor | `expand_outputs` not needed | Accepted, verified [code: vllm_generation.py:745-751, :87-88, :618; grpo_trainer.py:1844]. `group_requests` returns `prompts[::G]`; engine limits taken from code (256 seqs, 4,096 batched tokens). |
| 23 | minor | Decode tail has no task | Accepted. C0 step 6b measures tail seconds; §8.3 hands it to Phase B/D; C7 text notes truncated answers are pure waste. |
| 24 | minor | C10 list vague; rollback text goes stale | Accepted. Deletions named with fast.py lines; §7 note "after C10, rollback = `git revert <sha>`". |
| extra | - | Not in the review | Box stack differs from the Mac (torch 2.13.0, transformers 5.19.0 vs 2.14.1, 5.20.1) [code: docker/locks/vllm_native.txt; measured: runs/r5/env_check.json]. Revision 1 said both used the Mac pins. Checked that `mm_encoder_outputs` and the fused-factor function exist in 5.19.0 [web: v5.19.0 sources]; all bit gates run on the box stack. |

---

## C.R3. Review log (revision 3, answers the final review round 1, `docs/implementation2.0/reviews/final_r1_*.md`)

Sources: `rl` = final_r1_rl.md, `cuda` = final_r1_cuda.md, `fa` = final_r1_fresh-arch.md, `arch` = final_r1_arch.md,
`qa` = final_r1_qa.md, `data` = final_r1_data.md. Only items addressed to section C, or front items with a C side, are
listed. Each fact was checked in code or a log before the text changed.

| # | Source, sev. | Item | What I did |
|---|---|---|---|
| 1 | qa-1 blocking; cuda-1 blocking (C side) | One rental table; C plans its own session; G1 runs C items before the Mac work | C.8.2 is now "GPU items inside the front's validation rental", not a rental. C items enter it only after every C1-C5 Mac test passes (C.6 protocol). Costs given without setup (31 / 66 min) and, for reference, with setup. Revision 4: the rental is V2 (front 6.3); the baseline arm's steps 16-17 ARE A6's capture (C.R4 items 1, 9); costs now include process start-up (C.R4 item 6). |
| 2 | cuda-2, qa-2 blocking; rl-15, fa-9 major (C side) | Front uses 18 s as the gate for B to use C code | Defined "Gate C-adopt" = G-EQ1..6 + G-SP1 (+ the Mac half of G-EQ2/3) in the header and C.6; G-C (18 s) is measured inside B. The front text itself is outside this file. |
| 3 | arch-1, qa-3 blocking (C side) | Phase order: E before or after C | Header rule "Paths are pre-E"; C.8.3 keeps "land C first, then E moves the modules unchanged" and adds the full destination list, plus "if E5 lands first, map through this table". |
| 4 | rl-4 major | KL weight 0 has no owner | C8 now follows the back-matter owner decision "beta for the B runs" for every B run (default proposed by rl: 0.0). G-C is gated only if beta = 0; else reported as a computed number (C8, C.6, C.8.3). Verified: with PEFT TRL keeps no reference copy [code: grpo_trainer.py:1030-1037]. |
| 5 | rl-10 major | C5 rejects vLLM LoRA serving on a false numerics claim | Verified the trainer never merges [code: fast.py:495-504, 696-717]. Rewrote the C5 reason: not adopted for speed; it is the numerics candidate for merge rounding; ADR-0013 option F (front numbering; revision 3 wrote 0012); optional 15-step decode arm in V3 if B picks F, not in V2 (C.R4 items 2, 3). |
| 6 | rl-13 major | Runaways get no gradient; a lower limit worsens it | C7 step 5, option O: TRL's `get_soft_overlong_punishment(L, 200)` added to R as one term before centring [code: .venv trl/rewards/other_rewards.py:83]. Correction to the review wording: the penalty is 0 up to L - 200 and -1 at L, so it reaches answers that end just below the limit; truncated rows stay masked. Decode-tail seconds logged per run. Test U-CO4. |
| 7 | rl-16, cuda-8 major | C3 bit-equality depends on the GEMM batch shape | Verified `VisionCache` encodes a varying subset per call [code: fast.py:269-291]. Added C3 step 0 (encode-grouping check), a same-grouping G-EQ4 baseline, a logprob bound for the grouping effect (revision 4: 2 x the measured noise floor, not 2e-3; C.R4 item 10), and a relative 1e-3 tolerance for `grad_norm` and logp difference in G-EQ6 (R5 base vs shared differ by 6e-4). The three reward numbers stay exact (step-1 sampling does not use the trainer's encode). C.1 goal states the exception. Tests G-EN0, E-EQ6. |
| 8 | fa-2 major | The decode tail (4-5 s) has no task | New task C11 "Runaway cut": cap_i = min(L, ceil(1.5 T_i) + 64) per row, vLLM `SamplingParams` list [web: vLLM `LLM.generate` docstring; unverified for 0.32.0, box test G-RC0]. Changes against the review text: sizing uses A3 VAL answers only (the test split is not used for a choice, data review); G-RC2 is statistical, because a per-row cap changes the batch and so is not bit-equal; files are pre-E (`train/vllm_requests.py`, `fast.py`), not `rl.trl_adapter`. C.8.1 row, C.8.3 row and the "tail goes to B and D" line replaced. Tests U-RC1..RC4, G-RC0, E-RC2. |
| 9 | cuda-6 major | Two nsys installs, an undefined env name, two window flags | C0 step 0 uses A6's pinned .deb with `dpkg-deb -x` in the last image layer; devtools repo dropped. C0 defines `OBJDEL_NVTX_SYNC` (1 = today's sync, 0 = events). `--profile-steps` deleted everywhere; C uses A6's `--profile-window 16 17`. Tests U-DK1, U-CF1, U-CL1 updated. |
| 10 | cuda-7 major | A3 answers feed the C gates | The Mac half of G-EQ2/G-EQ3 runs on A3's val answers; Gate C-adopt includes it (C2 tests, C.6, C.8.3). Revision 4: A3 runs in V1, so this Mac half runs after V1 and before V2 (C.R4 item 1). |
| 11 | cuda-9 major | Pinned staging reuse race | Verified torch's caching host allocator records an event per non-blocking copy and does not reuse the block before it [code: .venv torch/include/ATen/core/CachingHostAllocator.h:202-215, 505]. C1 step 1 and C4 step 4 allocate per step / per chunk through it; the false "each is a cudaHostAlloc" sentence is gone; a persistent buffer only with an event per use and a test. Tests U-PX15, U-RM6, G-RM2. |
| 12 | cuda-10 major | G-SP4 not reachable with the optional tasks off | Verified DeepStack boolean indexing [code: modeling_qwen3_vl.py:847-855] and TRL's 7 metric reads. G-SP4 is now derived from the adopted set and report-only; E-SP4 matches. |
| 13 | cuda-11 major (C side) | K3 runbook unconditional in G1 | C was already conditional; the C.8.2 K3 row now names the trigger (`c_speed_check --pixels` fails both "hf" and "torch"). R-D1 placement is the front's. |
| 14 | cuda-13, qa-5, arch-14 major | Paths to a session scratchpad | Absolute paths removed from the header; path rule added; every kernel path now `docs/implementation2.0/kernels/K3/` (checked: the folder holds `k3_torch.cpp`, `GPU_RUN.md`, `build/test_hf_cpu.log`). Notes this phase cites (revision 4 cites each at its DOC1 step 0 copy): `docs/implementation2.0/notes/B_gpu_audit.md`, `docs/implementation2.0/notes/A1.md`, `docs/implementation2.0/notes/C1.md`, `docs/implementation2.0/notes/c_rev/render_bench.py`, `docs/implementation2.0/notes/c_rev/pass2_text_bench.py`, `docs/implementation2.0/notes/c_rev/C.v1.md`, `docs/implementation2.0/specs/C.spec.md`, `docs/implementation2.0/tests/C_tests_1.md`, `docs/implementation2.0/tests/C_tests_2.md`, `docs/implementation2.0/reviews/C_phase_review.md`, `docs/implementation2.0/notes/doc/nv_catalog.json`. The scratchpad lines inside the K3 runbook are outside this file (D1 step 0). |
| 15 | cuda-14, qa-15, fa-9 major | Several rates and step times | C.8.2 now uses one rate, $0.641/h (front 6), as revision 4 (C.R4 item 4); revision 3 had a second $0.60 column. Ledger launch prices $0.544-0.628/h [measured: ledger.csv]. Header defines R5 = 23.65 s (all steps) and 23.43 s (steps 3-15, C gates only) [measured: runs/r5/shared_rl.json]. |
| 16 | cuda-5, data-13, fa-6 major | View bank has no task; preload breaks fork-before-CUDA | View bank DEFERRED (C.3, C.8.3) with trigger (C0 decode > 0.2 s per step) and the rules for any later task: B1 rows root as source, every image checked by `_pixels_sha` [code: dataset.py:86-87], pageable pre-fork memory [code: serve.py:167-169], 32 planes per step to pinned staging, GPU-resident only with measured free memory. Trackio and `events.jsonl` (fa-6) are not added to C0: no C gate reads them; they stay DEFERRED. |
| 17 | arch-6 major, arch-18 minor | E homes for C modules; C names E homes E does not define | C.8.3 lists `policy/caches.py`, `rl/forward.py`, `policy/gpu_pixels.py`, `rl/vllm_sync.py`, `rl/vllm_requests.py`, and the NVTX home E sets; C.0 and C.4 no longer name `rl/caches.py` / `rl/sync.py`. |
| 18 | arch-2 blocking (C side) | NVTX module and import contracts | `count_syncs()` imports torch inside the function; the module home is left to E (C.8.3). |
| 19 | qa-8 major | Marker `cuda` vs `gpu` | One marker `gpu` in C.6, C.9.0, C.9.4, C.9.5; the strict-markers sentence corrected (no `addopts` today [code: pyproject.toml:66-71]); C adds no marker. |
| 20 | qa-10 major | Open behaviours with either-outcome tests | C.9.9 rewritten as ten rules; U-PX5, U-PX9, U-IM3, U-IM8, U-RP2, U-RM2, U-RM3, U-SY3, U-MS1, P-RM2, I-MS1, I-CF1 assert one rule each. Rule 3 and 7 verified in TRL [code: grpo_trainer.py:2604-2608]. |
| 21 | qa-13 major, cuda-22 minor | "$17.50" sentence; C8 not in the owner list | Sentence deleted; C.8.2 cites the binding $3. C8 names the owner decision; adding the row is a 99_back edit. |
| 22 | qa-6 major (C side) | One id, several meanings | Header: other phases cite C gates as `C-G-...`. Renaming the front's rentals is the front's. |
| 23 | qa-14 major (C side) | Synonyms | Header "Words" rule (answer = TRL completion; view; image; test); "crop" replaced in U-PX13. |
| 24 | qa-16 major (C side) | `-I` for data tests | C.9.0: `data` tests run with `.venv/bin/python -I -m pytest -m data`. |
| 25 | qa-17 major | Long sentences | Partly: all revision-3 text uses short sentences; the older table cells were not rewritten in this round. |
| 26 | qa-28 minor | Undefined tags | Header defines `[measured now]` and `[unverified]`; the one `[measured/estimate]` tag split in two. |
| 27 | cuda-20 minor | nsys sync validation includes non-torch syncs | C0 step 1 compares per thread with `cuda_api_trace` (trainer thread only); G-MS1 updated. `gc.stats()` named by the review does not exist in today's engine, so it is not cited. |
| 28 | cuda-21 minor | CPU sampling may be blocked | C0 step 6a falls back to the NVTX `generate` start when `nsys status -e` reports no CPU sampling. |
| 29 | cuda-25 minor | No guard for 32-image `pixel_values` with a 256-row grid | C3 step 2 `EncoderMismatch` guard; U-EN7 runs TRL's full step with it on. |
| 30 | arch-22 minor | No build home for K3 | C1 step 7: `cuda/k3/`, separate torch extension `objdel._k3` (revision 4: built on the box against the image's torch with `--no-build-isolation`, not through an extra; C.R4 item 13), TID251 ban with one exception (ruff selects E, F, I, B, UP today [code: pyproject.toml:88-89]). |
| - | rl-23, fa-10, cuda-12, qa-25 minor/major | Front diagrams show default-off C work | No change here: C.3.1 already draws vLLM's CPU lane and marks C6 as an arm; the front diagrams are outside this file. |
| - | qa-24 minor | A and C disagree on vLLM repeatability | No change in C: G-EQ6 already rests on R5's measured identical step 1 and E-EQ6 runs the baseline twice; the A.7 text is A's. |
| extra | - | ADR folder | C10 and U-CL1 named `docs/decisions/`; the repo folder is `docs/architecture/decisions/` (last file 0011) [measured: ls]. Fixed; the number comes from the front's allocation table. |

---

## C.R4. Review log (revision 4, answers the final review round 2, `docs/implementation2.0/reviews/final_r2_*.md`)

Sources: `arch` = final_r2_arch.md, `cuda` = final_r2_cuda.md, `fa` = final_r2_fresh-arch.md, `qa` = final_r2_qa.md,
`rl` = final_r2_rl.md (`data` has no C item). Only items addressed to section C, or front items with a C side, are
listed. Each fact was checked in code, a log or a measurement before the text changed.

| # | Source, sev. | Item | What I did |
|---|---|---|---|
| 1 | qa-1 blocking (C side) | Three plans for the A3 slot and the A6 capture | C0 step 3: "If the front picks QA option (a)" deleted; the baseline arm's steps 16-17 ARE A6's capture (`c0_base`), option (a). A3 is V1 item 4; the Mac half of G-EQ2/G-EQ3 runs after V1, before V2 (C2 tests, C.6, C.8.3). The front, A and 99_back edits are theirs. |
| 2 | arch-3 blocking; fa-5, qa-4, rl-9 major | C uses the colliding ADR number 0012 | Every C use of the old number 0012 for the number match -> `ADR-0013` (C5, C.8.3, C.R3 row 5); a grep for `ADR-?0012` in C now hits 0 lines. C10's ADR is ADR-0022, and C10 names the front row to add ("0022 \| GPU pixel backend, K3, one-call weight push \| C10"); the front table edit is the front's. |
| 3 | rl-1 blocking, part (e) | Under option F the push is not W + d | C5: under F the push is the adapter only, about 17.4 M parameters [estimate: r 16 x (in + out) x 7 modules x 28 layers = 622,592 x 28, shapes from C5], and G-EQ5 compares vLLM's LoRA tensors with the trainer's. The F decode arm moves to V3 before the B3 run, only if B picks F. B2.0 and rule 0 are B's (B.md has not changed yet, so C names them as "proposed by RL round 2 item 1"). |
| 4 | fa-4 major; cuda-17 minor | Stale rental map and rate | Header "Rentals" paragraph (V1, V2, V3). C.4 order, C0 files, C.6 protocol, C.8.2 title: "validation rental (V2, front 6.3)". One rate, $0.641/h (front 6, `docs/runbooks/gpu-run.md:136-137`); the $0.60 column is gone; payback recomputed (-18 min = -$0.19 per 300 steps). C11 step 3: "first rental after C11 lands = V3, before the B3 run". C.R3 row 15 annotated. |
| 5 | rl-3 blocking (C side) | Masking truncated or capped rows is not neutral | Verified the bias source: TRL zeroes only the loss mask [code: grpo_trainer.py:2603-2611], an unparseable or truncated answer scores 0 [code: src/objdel/reward/score.py:24]; runaways 0.00059 -> 0.00373 and 1,536-token steps 13.8 % -> 31.5 % (steps 21-100 vs 501-770) [measured: runs/rl/train_rl/rl.json, recomputed 2026-11-08]. C7 "Facts" now states the -b·∇P(not truncated) bias; the B3 setting is B's. C11 step 4: capped rows are not masked; `--runaway-cap > 0` requires `mask_truncated_completions=False` (U-RC5, C.9.9 rule 13); the optional "mask only the first min(len, T_i) tokens" is named, not default. C.7: new risk row; the C11 row no longer says "masked". Option O default off (item 19). |
| 6 | cuda-3 major | Arm minutes omit process start-up | Verified from the ledger: r4 31.3 min `run` to `teardown` for 15.0 min of steps in 2 arms; r5 57.3 min for 13.7 min of steps + 17.5 min eval. C.8.2 adds 6 min start-up per process [estimate]. The check arm continues in the fast arm's process (steps 18-32, new flag `--checks-from-step`, U-CF1), which saves one start-up. Required C items 49 min = $0.52 without setup, 89 min = $0.95 with V2's 40 min setup; this does not fit the $0.59-0.79 left after V1, so C adoption needs the owner's budget decision (stated in C.8.2). The front 6.3 table and decision 1 are the front's. |
| 7 | qa-10 major | G-SP1 has a gap between pass and fail and moves during the rental | G-SP1: pass <= 20.0 s (A) / <= 20.7 s (B); anything above = not adopted; above 21.2 s / 21.9 s also = re-profile. "Recomputed before the fast arm" deleted in G-SP1, C5 step 5, C.8.1 and E-SP1; the recompute is reported beside the gate and never changes it. |
| 8 | qa-12 major; arch-15 minor | C still cites the scratchpad | Header `S[P]/` paragraph deleted; path rule points at the front section 0 rule. Every scratchpad path replaced by its DOC1 step 0 copy: notes -> `docs/implementation2.0/notes/`, spec -> `docs/implementation2.0/specs/C.spec.md`, tests -> `docs/implementation2.0/tests/`, reviews -> `docs/implementation2.0/reviews/`, catalog -> `docs/implementation2.0/notes/doc/nv_catalog.json`, `c_rev/` -> `docs/implementation2.0/notes/c_rev/`. 0 scratchpad paths left in C [measured: grep]. The DOC1 step 0 row for `impl/c_rev/` is DOC's edit. |
| 9 | cuda-11 minor | Two nsys command lines; `--cuda-memory-usage` | C0 step 3 now calls A6's `nsys_capture.sh {rd} <arm> 16 17 -- ...` with no other nsys flags; peak memory from `torch.cuda.max_memory_allocated` and `nvidia-smi`. C0 step 5 reads the CSVs A6 exports. The fast arm sets `NSYS_TIMEOUT=1800` (it runs 32 steps; A6's default is 1,200 s). |
| 10 | cuda-4 major | The 2e-3 bound is not calibrated | C3 step 0 also encodes in two groupings of today's path (A, B) and records `floor_max` and `floor_p99` of the logprob difference. G-EQ4's second check: batched vs chunked-A `p99 <= 2 x floor_p99` and `max <= 2 x floor_max`; `torch.equal` if the floor is 0. Fixed 2e-3 dropped from G-EQ4, C.7, C.9.0; G-EN0 updated. |
| 11 | cuda-5 major | CUDA divides by a CPU scalar as a reciprocal multiply | Verified: 126 of 256 uint8 values differ in float32 [measured: numpy, 2026-11-08]; Qwen3-VL-2B mean = std = 0.5 [measured: local HF cache]. C1 interface "divisor rule": `(3, 1, 1)` float32 tensors on `device`. New tests U-PX16 (type and device of the divisor) and U-PX17 (the risk is real); new C.7 risk row. |
| 12 | cuda-6 major | `bad` rows read once in `log` cannot be routed in their step | Verified TRL fills the answer's `mm_token_type_ids` with zeros [code: grpo_trainer.py:2685-2697]. New host function `rope_bad_rows` (prompt ids vs template ids; vision token ids in the answer) routes rows to HF in the same chunk with no device read. The device `bad` is read only on `--check-rope-every` steps, and a device-bad, host-good row raises there. U-RP7, C.9.9 rule 12, C.7 row. |
| 13 | cuda-20 minor | K3 cannot be built "with the `train` extra" | Verified the build backend is scikit-build-core + nanobind only [code: pyproject.toml:35-37]. C1 step 7: `cuda/k3/setup.py` with `CUDAExtension` and `--fmad=false -arch=sm_120`, installed with `--no-build-isolation` in the image against its torch 2.13.0+cu130, or `torch.utils.cpp_extension.load` at start-up; lazy import, the self-check decides. |
| 14 | cuda-14 minor | G-EQ5's full device copy may OOM in mode A | C5 step 5(b) and G-EQ5: per-tensor SHA-256 of the int16 view, hashed on the host one tensor at a time after each push; no device snapshot. G-SY1 asserts the device memory rise stays under 100 MB. |
| 15 | cuda-12 minor | "The devtools repo stops at 2021.4.1" is false | Verified: `devtools/repos/ubuntu2404/amd64/Packages` lists `nsight-systems-cli-2026.6.1` (232,877,734 B, SHA-256 `e7144097...cb56b` in full in C0 step 0) [measured: curl, 2026-11-08]. C.0 and C0 step 0 now state that fact; the pinned-.deb decision stays, and switching to the CLI-only file is left to A6 (A6 owns the pin). |
| 16 | qa-9 major (C side) | U-CF2 accepts either outcome | U-CF2: an existing non-empty directory raises `FileExistsError`; nothing is written. C.9.9 rule 11. |
| 17 | qa-19 minor | 1,182 vs 1,181 | Header "Outline length" rule: lengths count the end token; 1,182 with, 1,181 without [measured: A1.md F9 probe adds 1; B1.md:90]. Found a real error while checking: C11 said `dataset.target_tokens` includes the end token, but `_tokens` counts none [code: dataset.py:185-186, 209]. C11 now uses `target_tokens(...) + 1`; U-RC1 states it. |
| 18 | qa-22 minor | C.9.0 describes A's marker rules as proposals | C.9.0: "A registers `gpu` with `--strict-markers`; `data` tests run only through `make test-data` (A.9)". |
| 19 | rl-7 major | Option O at L = 1,200 penalises the largest correct outlines | C7 step 5: option O default off; if kept, `L - c >= max doctor outline + 32` (1,214 today), so off at L = 1,200 and `c <= 322` at L = 1,536; the run config check raises otherwise (U-CO4 updated). Correction to the review: 1,536 - 1,214 = 322, not 290 [arithmetic]. |
| 20 | arch-5 major (C side) | C11's `target_tokens` home after E5 | C11 rule and C.8.3 Phase E row: after E5, C11 imports `target_tokens` from `policy.rows`, never `sft/rows.py` (`rl \| sft` contract). The E5 move table is E's. |
| 21 | qa-5 major (C side) | Gate prefix rule inside a phase | Header: inside C the `C-` prefix is dropped; gates of other phases keep their prefix (`A-G6`, Gate B1). C's cross-phase references already use prefixes [measured: grep]. |
| 22 | qa-25 minor (partly) | Long sentences in gate cells | Split G-SP4 (was one 96-word cell) and G-SP1 into one rule per sentence. Older task text not rewritten in this round. |
| - | outside C | Front and DOC edits named above | Front 6.3 start-up rows and decision 1 wording (item 6), front ADR row 0022 (item 2), DOC1 step 0 row for `impl/c_rev/` (item 8), A's G1a/A3 text (item 1): not edited here; this file edits only C. |

---

## C.R5. Review log (revision 5, answers the final review round 3, `docs/implementation2.0/reviews/final_r3_*.md`)

Sources: `rl` = final_r3_rl.md, `cuda` = final_r3_cuda.md, `fa` = final_r3_fresh-arch.md, `qa` = final_r3_qa.md
(`arch` and `data` have no C item). Only items addressed to section C, or front items with a C side, are listed. Each
fact was checked in code, a log or a web source before the text changed.

| # | Source, sev. | Item | What I did |
|---|---|---|---|
| 1 | qa-5 blocking | Gate C-adopt deadlocks when V1 skips A3-val | Option (a). The Mac half of G-EQ2/G-EQ3 is no longer a V2 entry condition. It runs on A3's val answers after V1 (or after V2 when item 3b makes them) and must pass before V3. Header, C2 tests, C.6 protocol and Gate C-adopt, C.8.2, C.8.3 A3 row. V2 order puts item 3b first. |
| 2 | rl-1 blocking (C side) | Under F, level-2 sleep discards the base weights that an adapter-only push never restores | Verified `vllm_enable_sleep_mode: True` [code: configs.py:193], `sleep(level=2)` after every generate and the re-push on the next sync [code: vllm_generation.py:771-773, 500-528, 602-603]. C5 "Under option F" now has the 4-step F sync: wake; mode B base push `values="base"` (p.data, no delta), mode A none; one `add_lora(load_inplace=True)`; `reset_prefix_cache`. New `push_weights(..., values=)`, `lora_merged_params(with_delta=)`, `sync_lora` (one home for every sync; B2.2 calls it). G-EQ5 under F runs after a sleep/wake cycle. U-SY8, G-SY1, C.7 rows. The B2.2 and U70 text is B's. |
| 3 | rl-3 blocking (C side) | G-EQ5 under F ("LoRA tensors equal") can never pass | Verified alpha / r = 32 / 16 = 2 [code: configs.py:95]. G-EQ5, C5 step 5(c) and G-SY1: packed map (q/k/v -> `qkv_proj` 0/1/2, gate/up -> `gate_up_proj` 0/1), `lora_a == A.to(bf16)`, `lora_b == 2 * B.to(bf16)`, `torch.equal`; the adapter is written in bf16; plus B's output-level check (way (f) vs (b) gap at most half of B1's), recorded by C. The vLLM internals (`lora_dtype` auto, `optimize()` folding, packing) are taken from the RL review's vLLM 0.32.0 source reads; G-SY1 on the box is the check. |
| 4 | fa-3 major | C measures the merged push while B2.0 already points to F | Header "Numerics before speed"; C.6: the F check (V2 item 2b, B's arm) runs before the C arms, unless B ran it in V1; both arms run under the preset it picks (`lora-seqmask`, else `bf16-seqmask`); a later `fp16-seqmask` pick re-runs G-EQ1/4/5 in V3. C.8.1 adds F rows (C1-C5 under F: A 1.47-3.97 s, B 1.25-3.0 s) and C.8.2 re-prices V2: 103 min = $1.10 with item 2b (the review's 102 min plus the 1-min batch step of item 7 below). ADR-0022 title per the review (item 11). C.3.2 step 3 states the F push. Front 2.5, front 6.3 and decision 1 are the front's. |
| 5 | rl-7 major (C side) | The F arm has two prices and a cross-host timing rule | C.8.2 uses one price, 15 steps + 6 min start-up = 12.9 min = $0.14, as item 2b; the in-process LoRA/base generate-time ratio is named as the F check's timing rule (B owns it). The "V3 decode arm" sentences in C5 and C.8.2 are gone. Where B finally places the arm (V1 start, V2 item 2b, or B3's own start) is B's; C needs only that it runs before the C arms. |
| 6 | cuda-3 major | G-SP1 is a fixed number while host and profiler move the step | G-SP1 and E-SP1 are relative to the same-rental baseline arm: merged `base_mean - 3.4 s` (A) / `- 2.7 s` (B); under F `- 2.7 s` / `- 2.1 s`. Re-profile margins `- 2.25 / 1.55 s` (merged), `- 1.5 / 1.25 s` (F). Both means and their ratio are reported; "host not comparable to R5" outside 23.43 s ± 10 %. C.8.1's last sentence: 20.0 s / 20.7 s are planning values only (the front row already says so). New U-SP1 with values where the old and new rules disagree. |
| 7 | cuda-5 major | `c_speed_check` needs recorded batches that V2 makes only later | Changed against the review's fix: A3's Mac batches do not exist before V2 in the base case (V1 at 27.5 s skips A3-val, item 1), so `c_speed_check` makes its own 5 batches on the box first (5 x 32 seeded train views from the V2 bundle; one vLLM generate each; `--record-batches` format; seeded advantages centred per group of 8, because the gradient check compares two code paths on the same loss). They also count as box batches for G-EQ2/3. Cost +1 min (c_speed_check 13 min). C1 paragraph, C3 step 0, G-EQ1, G-EQ4, G-PX1, G-EN0, G-EN1, new U-CB1. |
| 8 | cuda-4 major (C side) | GPU-busy share misses vLLM decode in CUDA graphs; driver-API syncs not counted | Verified TRL's `LLM(...)` has no `enforce_eager` [code: vllm_generation.py:356-372; 0 hits in the file]. C0 step 4: busy intervals from A6's `$NAME.sqlite` (kernel, memcpy, memset and any `GRAPH` table with `start`/`end`, found by `sqlite_master`, because no current nsys schema page names the table [web search, 2026-11-08]); box check once with and without graph rows; "decode not visible" when absent (G-SP6 is report-only). C0 step 1 adds `cuStreamSynchronize`, `cuCtxSynchronize`, `cuEventSynchronize`, `cuMemcpyDtoH*` to the vLLM count. New U-GA3. A6 step 9 is A's. |
| 9 | qa-9 major (C side) | `pytest -m gpu` collects every phase's box tests | C.9.0 and C.6 code gates: C's pytest box tests are only `tests/test_c_gpu.py`, run by name in V2 item 2; real-model checks are `c_speed_check` modes; artifact-dependent tests skip with a reason; `pytest -m gpu` without a file list never runs on a box. C.9.4 opens with a table mapping every C box test to its rental item and command (G-RP2 and the G-RM1 equality part move to the Mac as `data` tests). D-G3/D8 are D's. |
| 10 | rl-11 minor (C side) | Stale option-F text in C5 and C10 | C5 cites B2.4 rule 0 and the B2.0 seed value; the "RL review round 2 proposes ... B.md has not changed" sentence is deleted. ADR-0022 title: "GPU pixel backend, K3, one-call weight push (merged or adapter-only per ADR-0013)"; body keeps "not adopted for speed". Front table row text named for the front. |
| 11 | rl-12 minor | C8 tells engineers to set beta 0 in "the B run configs" | Verified decision 15 (B1 1e-3, B3 `1e-3 / r`) [99_back]. C8 Files: beta 0.0 only in `b3_v5.yaml` and only if the owner picks it; `b1_rows.yaml` keeps 1e-3. New U-CO5 asserts it. |
| 12 | cuda-6, qa-14 minor | C.8.2 / C.8.3 stale numbers | Money left $0.74 (27.5 s, item 5 skipped) / $0.87 (23.65 s) [front 6.3, checked]; the "front's V2 table omits start-up rows" clause deleted (front 6.3 has them); C.8.3 A3 row "front 6.2 item 5 (or V2 item 3b)". |
| 13 | cuda-8 minor (C side) | K3 production build flags differ from the tested build | Verified [code: K3/test_torch_gpu.py:41-45; K3/GPU_RUN.md:64]. C1 step 7 copies the recipe: `TORCH_CUDA_ARCH_LIST=12.0`, cxx `-O3 -ffp-contract=off`, nvcc `-O3 --fmad=false -lineinfo` plus the three `-U__CUDA_NO_*` flags; `-arch=sm_120` removed. E.5 is E's. |
| 14 | qa-12 minor (C side) | C says A adds a CUDA skip hook that A does not add | Checked A.md: no root `gpu` skip hook (grep `conftest`). C.9.0 now relies on the one root hook the review asks A.9 to add (`objdel._gpucheck.build() == "cuda"`, no torch import) and states that A does not have it yet; C adds no skip line. |
| 15 | qa-13 minor (C side) | Two definitions of outline length | Header: outline length = tokens without the end token (1,181); "with the end token" (1,182) is written where a count includes it; the answer cap and C11 count it, so C7's rule stays `1,182 + 32 = 1,214`. C7 "Facts", C7 step 5 and U-CO4 reworded. Front section 0 row is the front's. |
| 16 | qa-21 minor (C side) | Long sentences in G-SP4 | G-SP4 cell split into one rule per sentence; the two removed sources each in their own sentence. New text in this revision uses short sentences. |
| 17 | fa-4 major (C side) | C/D Mac work does not depend on the Gate B1 row | C.8.3 "Phase C blocks" states it: every row of the front's after-V1 decision table keeps the C Mac work; only the arm preset reads a B result (B2.0, $0). The decision table is the front's. |
| - | rl-8, rl-2 (B side) | F in B1 vs B3; `load_inplace` on the request | B's items. C5 states the call-site rule that C's `sync_lora` follows (`add_lora` once per sync; requests with `load_inplace=False`) and adds `LLMEngine.add_lora` to the A7 list (C.4, C.9.9 rule 10, U-LM1). |
| - | outside C | Front and B edits named above | Front 6.3 item 2b and the 13-min `c_speed_check` row, front 2.5 push label, front ADR row 0022 title, front section 0 outline-length row, B2.2/U70/G05 F text, A.9 root skip hook: not edited here; this file edits only C. |

## Review log (consistency pass, round 4 prep)

| Item | What changed |
|---|---|
| One source | C.8.2 holds no V2 table; it names C's rows in front 6.3 and keeps only the phase-local basis (start-up, arm minutes, payback in minutes) |
| Skip rule / base case | "At 27.5 s V1 skips A3-val" removed (C1, C2 data, C.6 entry, C.8.3): item 5 is not skippable; fallback V2 item 1b |
| C11 vLLM check | Named as V3 item 0 (front 6.3) |
| Outline limit (owner question) | C7 step 1 points to front section 1: 1,800 lengthens the wait; L comes from A8's count |

## Review log (final round 4)

Cross-file consistency fixer, round 4. Sources: `reviews/final_r4_{cuda,rl,qa,arch,fresh-arch}.md` and their verdicts.

| Item | Severity | Change | Evidence |
|---|---|---|---|
| rl-1, qa-1 | blocking | Header, C5, C.6 "Numerics first", baseline arm, C.8.2, C.8.3 row: C arms run under `bf16-seqmask` (rule (b)); no F arm in V2; the F check is the first 15 steps of the first run after Gate B3, its step time a report line | B.md B2.4 "Scope", B.8.4 |
| rl-5 | major | C5 step 5(c): output-level comparator measured inside the F run; B2.2 F step 4 holds the rule | B.md B2.2 |
| rl-6 | major | C.1, C7 title and Files, C11 step 5, C.6 G-C and G-RC2, E-GC: L and the runaway cap apply from the first run after Gate B3; `b3_v5.yaml` never sets L | B.md B3.9 |
| qa-2 | blocking | `c_speed.yaml` step 0 = CUDA build (C0 file list, header, C.9.4); C.9.0 notes the A.9 fail-on-box rule | `policy.py:432-433` |

## Review log (final round 5)

Cross-file consistency fixer, round 5. Sources: `reviews/final_r5_{cuda,rl,data,arch,qa,fresh-arch}.md` and their
`*_verdict.md` files. Only the CONFIRMED items assigned to this fixer are listed.

| Item | Severity | Change | Evidence |
|---|---|---|---|
| rl-1 | major | C5 F sync step 3: the adapter folder uses B2.2's PEFT form (`adapter_model.safetensors` without `.default.` keys, `adapter_config.json`). U-SY8 and the C5 Mac line check the file name, the key pattern and the config | B2.2 F step 2.3; checkpoint-100 header |
| cuda-2 | major | C1: the V2 bundle is built with `--view lidc-bundle-v3 --view lidc-train-shifted --view lidc-train-rescaled` (front 6.3); front 6.1 step 4 checks the manifest | `bundle.py:55-56` |
| fresh-arch-1 | blocking | C.8.2 "Budget": V2 runs only after a budget raise; under (b') it does not run | front 6.3 |

## Review log (final round 6)

Cross-file consistency fixer, round 6. Sources: `reviews/final_r6_{cuda,rl,data,arch,qa,fresh-arch}.md` and their
`*_verdict.md` files. Only the CONFIRMED items assigned to this fixer are listed.

| Item | Severity | Change | Evidence |
|---|---|---|---|
| rl-1, rl-2 | blocking / major | C5 F step 2 and the interface: the comparator's base restore is the one mode-A exception; `sync_lora(..., adapter_name, slot)`; step 5(c), U-SY8 and G-EQ5 point to B2.2 F step 4b and hash the base against `p.data` after the restore in both modes | `fast.py:716-717` |
| cuda-2, cuda-3, qa-4 | major | C0 files: every `c_speed.yaml` arm reads `--ckpt {rd}/ckpt` (fetched by the V2 lines); new `c_launch.yaml` with a load test in `tests/test_preflight.py` | `policy.py:315-316, 527-536` |
| fresh-arch-1 | blocking | C10 and the C.8.3 Phase E row point to the front "End state under (b')" | front section 1 |

## Review log (final round 7)

Cross-file consistency fixer, round 7. Sources: `reviews/final_r7_{rl,arch,qa,fresh-arch}.md` and their
`*_verdict.md` files. Only the CONFIRMED items assigned to this fixer are listed.

| Item | Severity | Change | Evidence |
|---|---|---|---|
| rl-1 | blocking | C5 step 5(c) and G-EQ5: the base hash after the comparator restore runs right after the restore push, before mode B's `sleep(level=2)`; B2.2 F step 4b.5 holds the order | `.venv trl/generation/vllm_generation.py:378` |
| fresh-arch-1 | blocking | C0 files: every C arm trains on the B1 rows root (`--train {data}/rows/lidc-train-b1` plus both rows flags); whole 20-step blocks, so 17-step arms run `--steps 20` and the fast arm `--steps 40` (read windows unchanged). `c_launch.yaml` gets `rows_plan`; C1 batches come from the rows root; C2, C.6, C.8.2 (item 3 43 min), C.8.3, C9, G-PX1, the protocol and E-C9 follow. New U-CF4: every `-m` step of `c_speed.yaml` parses with its module's parser; `train_rl` steps name the rows root and a multiple of 20 | B.md B1.3 (`--steps % block_steps == 0`, `block_steps=20`); B1.5 deletes `--train-views`, B1.8 the no-schedule branch |

## Review log (final round 8)

Cross-file consistency fixer, round 8. Sources: `reviews/final_r8_{cuda,data,arch,qa}.md` and their `*_verdict.md`
files. Only the CONFIRMED items assigned to this fixer are listed. Each fact was checked in the repo first.

| Item | Severity | Change | Evidence |
|---|---|---|---|
| qa-2 | major | C0 step 3 holds the one full baseline argv (B1 rows root, both rows flags, `--steps 20`); U-CF4 (C0 files and C.9) also parses `train_rl` commands inside `bash -c` strings with `shlex.split`; planted nsys-wrapped arm with `--train {data}/prepared` or `--steps 17` fails | `policy.py:486` |

## Review log (final round 9)

Cross-file consistency fixer, round 9. Sources: `reviews/final_r9_{cuda,rl,data,arch,qa,fresh-arch}.md` and their
`*_verdict.md` files. Only the CONFIRMED items assigned to this fixer are listed. Each fact was checked in the repo first.

| Item | Severity | Change | Evidence |
|---|---|---|---|
| rl-1 | major | C5 step 5(c): comparator prompts as B2.2 4b.3 (the 32 x 8 batch through `fast.vllm_generate_checked`, never raw) | `fast.py:618-622, 680-688` |
| qa-3 | blocking | C.9.5 E-SP1: deleted "std of steps 3-15 at most 0.3 s in each arm, else noise, rerun". Now the per-arm median and SD and the Welch 95 % interval of `base_mean - fast_mean` are report lines only; they never change G-SP1 (C.6 stays the one decision rule) and no arm is rerun. Added the reference SDs and the SE estimate (0.79 s, so the 3.4 s margin is about 4 SE) | recomputed 2026-11-09: `runs/r5/base_rl.json` steps 3-15 mean 30.78 s, SD 1.89 s; `runs/r5/shared_rl.json` mean 23.43 s, SD 2.11 s |
| qa-4 | major | C.9.5: E-EQ6's control baseline and bisect by flag, E-REP, E-MEM's forced-1,536 run and E-SP4's B = 128 arm are marked "optional, not in V2; priced in front 6.3 only if the owner adds it". Without them a step-1 mismatch fails G-EQ6, so Gate C-adopt fails (front 6.3 item 3 stop rule). E-REP is a report line; it is removed from the C.9.7 rows for G-EQ6 and G-SP1. Front 6.3 is unchanged | front 6.3 V2 item 3: one baseline arm and one fast arm with the check steps, 43 min, $0.49, basis "12 + 8 + 7 + 4 + 2 x 6 start-up"; C.8.2 lists the same arms |

## Review log (final round 10)

Cross-file consistency fixer, round 10. Sources: `reviews/final_r10_{cuda,data,qa,fresh-arch}.md` and their
`*_verdict.md` files. Only the CONFIRMED items assigned to this fixer are listed. Each fact was checked in the repo first.

| Item | Severity | Change | Evidence |
|---|---|---|---|
| qa-2 | major | C0 files: `c_speed.yaml` holds items 2-3 only; new `c_speed_b.yaml` (3c) and `c_speed_c6.yaml` (3d), one `run.sh` call each; new `train_rl --expect-step1 FILE` (E-EQ6 comparison after step 1, exit 3); fast arm `--out {rd}/c_fast`; C.6, E-EQ6, U-CF4 updated; new U-CF5; flag kept in the C10 list | `policy.py:432-433, 487-494` |
| qa-3 | major | C9 has one home (C.8.2): no plan file, no rental row; C0 arm list, C.8.2, C.8.3 row, G-CF1 runner row and E-C9 say so; `c_mem.yaml` is written only if the owner adds C9 rows to front 6.3 | front 6.3 V2/V3 tables hold no C9 row |
| qa-5 | major | Check arm writes `--record-batches {rd}/c_check/arm_batches` (apart from C1's `{rd}/c_check/batches`); C.9.4 runner row, G-RP2 and G-RM1 read `runs/$OBJDEL_V2_RUN/...` and skip without it; fewer than 15 box batches fails | C.9.9 item 11; `policy.py:498-507` |

## Review log (final round 11)

| Item | Severity | Change | Evidence |
|---|---|---|---|
| rl-1 | major | C5 Interface "Sleep state": TRL's `vg._llm_weights_sleeping` is the one sleep state; push/sync wake only when True and clear it; outside sleeps set it True. F sync step 1 follows. U-SY3/U-SY8 assert the flag; new U-SY9 runs a real TRL `VLLMGeneration` around the fake `LLM` for 3 steps (1 `load_weights`, 1 `add_lora` per step); G-SP2 adds `time/sync_calls` = 1 | `vllm_generation.py:379, 503, 602-603, 772-773`; `grpo_trainer.py:1203, 1838-1841` |
| qa-1 | major | C3 step 0 and G-EN0 write `{rd}/c_check/encoder.json` (`bits_equal`, `max_abs`, `floor_max`, `floor_p99`); `--expect-step1` requires `--encode-check`; the comparison takes `bits_equal`; E-EQ6 names the key; U-CF5 adds the 0.003266/0.003264 cases; U-CF4 asserts the flag | R5 `grad_norm` 0.003265786 vs 0.003263979 [measured: runs/r5/base_rl.json, shared_rl.json]: 4 s.f. fails, 1e-3 passes |
| qa-2 | major | G-RC0 = `b3_v5.yaml` step 2, G-GR1 = `c_speed_c6.yaml` step 0, both by node id; the V2 item 2 pytest step deselects them; C11 step 3, test entries and U-CF4 follow | `policy.py:432-433`; B.9.5 rule |
| qa-3 | major | `c_speed_b.yaml`'s mode B arm carries check steps 18-32 (`--steps 40`, `arm_batches_b`, `--expect-step1`, `--encode-check`); C.6, C.9.4 runner row and G-RP2 read the adopted arm's folder | `policy.py:490`; C5 step 2 (about 28 of 32 GB) |

## Review log (final round 12)

| Item | Severity | Change | Evidence |
|---|---|---|---|
| qa-2 | major | C11 step 3, C.9.0, C.9.4 map and table: G-RC0 is the last, soft step of `b3_v5.yaml`; G-RC0 and G-GR1 load `{rd}/ckpt` (HF `rl-full-3ep`, fetched in V2 and V3 item 0) with `max_model_len=256`, `gpu_memory_utilization=0.3`, `enforce_eager=True`; G-RC0 uses `ignore_eos=True` (exact caps); new U-RC6 (Mac, marker `data`): config, tokenizer and processor of `rl-full-3ep` load | `test_train_fast.py:24-35`; `runs/rl/train_rl/final` file list [measured: ls]; `rl_shifted.yaml:20` |


# Phase D

## D.R. Review log (revision 2, from the D phase review, 2026-11-08)

| # | Severity | Item | What I did |
|---|---|---|---|
| 1 | blocking | Certificate legs share one CCL | Verified (K1 header: RingMark/RingInit come from the same Merge/Flatten; certify.py:44-46). D3 now takes b0/b1 from the pixel label chain (D4b `topology`), states which legs share code, and never compares rings with `check_masks`' own `betti`. D5b dropped. Added `K1_MUTATE_FG4` and gate D-G16. |
| 2 | blocking | RIGA rim has no device path | Verified (build.py:239). Added pair queries (`pairs`, `ref_pairs`, `pair_min_area`) to `certify_polygons`, a rim test on 50 RIGA items, and the rim in spec B11. |
| 3 | major | Duplicate rule used drawn masks | Verified (metrics.py:265-268 hashes raw masks plus flags). Added `inter_raw` to `eval_sets`; new duplicate rule; deleted the output-format.md step (hash is internal, :335/:347). |
| 4 | major | D-G7 compared two hosts | Accepted with one change: the old engine is built into a separate venv and run in alternating subprocesses (two statically linked engines in one process is a risk, and the module name is a literal in `NB_MODULE(_gpucheck, …)` [code: gpucheck_py.cpp:110]). Pass rule 1.05x / 0.9x, median and IQR. |
| 5 | major | D-G8 on a 15-step noisy mean | Verified (mean 0.004603, median 0.004509, range 0.003903-0.005383). D-G8 is now a same-rental replay benchmark plus a wait count from engine counters. |
| 6 | major | `metrics()` runs inside RL | Verified (callbacks.py:162 -> audit.py:154). `eval_sets` on the scorer lane with a scratch cap, NVTX range, gate D-G14 in the next B/C rental; §1 item 9 corrected. Superseded in revision 4: no audit callback runs in live RL (round 2 log). |
| 7 | major | Expert caches are not offline | Verified (cache.py:222-238, 294). Added D7b `add_tables_from_raters`; host construction stays engine-free (fork safety, D-R11); corrected §3.2 and §3.3. |
| 8 | major | cv2 topology stays; D-G11 grep is weak | Verified (topology.py:77 and callers). Added D4b (device `topology`, cv2 into `topology_oracle.py`). D-G11 is now a runtime monkeypatch test. |
| 9 | major | No pool-memory task | Verified (gpucheck.cpp:78-84, 135-150, 1741-1748, 1774). Added D5c (stream-ordered pool, batched uploads) and gate D-G15, measured by engine counters instead of an nsys trace. |
| 10 | major | `-n 0` without xdist | Verified (pyproject.toml:29). Removed `-n 0`; added `make gpu-check DRY=1`, run by preflight for `d_gpu_test.yaml`. |
| 11 | major | Sanitizer chain undefined and uncapped | Accepted: named tests in `tests/test_gpu_sanitize.py`, small corpus, `timeout 300` per tool, times in `runs/d/sanitize.json`. Partly disagreed: racecheck on K1-draft does not apply (no shared memory [K1 GPU_RUN.md:92]); K2-draft's runbook already runs it. Racecheck coverage is the engine's `raster_kernel` buffer (gpucheck.cpp:418). |
| 12 | major | NVTX include path | Accepted: CMake include from the `nvidia.nvtx` wheel with a clear error; `run_box.sh` installs `[nvcc,cccl,nvtx]`, and `docker_check.sh` inherits it (it runs `run_box.sh --build-only`); docker_check gate before R-D1. Stated how C++ ranges relate to `OBJDEL_NVTX`. Superseded in revision 4: the wheel has no `nvidia.nvtx` package; D1 checks `nvtx3/nvtx3.hpp` under the CUDA include folder (round 2 log). |
| 13 | major | nsys may be missing | Accepted: per-stage CUDA events (`timing=True`), D5 step 4 rule ">= 20 % of score_words"; nsys optional. |
| 14 | major | "One call per level for all names" breaks `_build` control flow | Verified (build.py:220-243). Now per name, all remaining rungs in one call; the rim retry calls again. |
| 15 | major | Budget ignores failed launches | Verified ($0.2385 destroy, 18/20/31 min setups). Added a $0.25 reserve, a 35-min stop rule (not 20: r5's good launch needed 31 min), all ledger rows count; ceiling $1.40 with owner question Q4. |
| 16 | minor | Disabled `gpu.yml` is dead code | Agreed: no workflow in Phase D. |
| 17 | minor | Chunk/scratch split untested | Agreed: test-only `chunk_pixels`/`max_scratch_bytes` keywords and a byte-equality test. |
| 18 | minor | Grid reason covers model output only | Agreed: D1 logs rung vertex counts; D6 keeps the grid parked only if >= 512-vertex rungs are < 5 % of build time. |
| 19 | minor | Caller lists incomplete | Verified and agreed: added `test_breakit_geometry.py:17`, `cold_start.py:29`, `xarch_hashes.py:45`, phase0 scripts; extended `test_experiment_scripts.py`. |
| 20 | minor | Wrong cuCIM rejection reason | Agreed: §2 now gives the real reasons (Mac oracle, dependency size). |
| 21 | minor | Soft-Dice tolerance hides errors | Agreed and checked (bool mean = exact sum / K): `soft_dice` exact; 1e-12 for `nsd` only. |
| 22 | minor | "tests/goldens fixtures" do not exist | Verified (`tests/goldens/` holds only `polycore.npz`; R5 `items.jsonl` holds metrics, not completions). Named the inputs: 256 `lidc/val` target texts on their caches, later Phase A saved answers; 1,000 LIDC caches for the EDT gate; SHA-256 recorded. |
| 23 | minor | Device bit-packing multiplies H2D bytes | Agreed: D5 step 2 is conditional on the replay benchmark and comes after D7b. |
| 24 | minor | `int(n)` in the CUB scan | Agreed: `INT_MAX` guard; `NumItemsT` check in the docker_check container ($0). |
| 25 | minor | 3D and probes missing from §3.3 | Agreed: rows for `volume.py`, `probes.py`; `cache.py` boundary moves to the device in D7b. |

## Review log (final review round 1, revision 3)

Sources: `reviews/final_r1_{arch,cuda,qa,fresh-arch,data,rl}.md`, items addressed to D. Every code fact below was
checked against the repo on 2026-11-08. `data` and `rl` have no item for Phase D (their "D1-D16" are B data gates).

| Review item | Severity | What I did |
|---|---|---|
| cuda 1, qa 1 | blocking | Verified: D's `make gpu-check` needs D2-D8 and the pre-D2 engine. D now has no rental of its own. R-D1 rides the first validation rental (before B1); R-D2 rides the first B/C rental after D8 (header, D1 step 4, D8 steps, §8). If D1 has not merged before the validation rental, the bench moves to R-D2. Costs recomputed without setup: R-D1 32-35 min, R-D2 25-44 min, D share $0.52-0.83. Superseded in revision 4: R-D1 = V2 item 4 (optional), R-D2 = V3 item 1, priced at $0.641/h (round 2 log). |
| cuda 3 | blocking | Verified: `final_eval.py:69`, `g_probe.py:157`, `select_ckpt.py:390` fork the pool; `score_items` maps `metrics()` over it; `POOL_CAP = 16` [code: serve.py:132]. Added D7 step 6 (one process for device eval; `score_items` drops `pool`; forked child's `device()` raises), gate D-G17, tests I18 and B06, changed I06. Pool deletion owner: CLEAN16 step 1 (arch 17), after D7 step 6, not D9. |
| arch 1, qa 3, cuda 23 | blocking | D already used pre-E paths. Removed the "or D rebases onto E" option from D-R10; §8 Phase E bullet states one order and that every D path is pre-E. |
| arch 3 | blocking | D3's seam is the chosen design. §8 now states the E3a rows (`checker/rings.py`, `certify.py`, `topology.py`, `volume.py` -> `geom/`), that E has no rows for the oracle files, and that CLEAN11 must not archive the live `certify.py`. |
| arch 4, qa 7, fresh-arch 7 | blocking/major | ADR renamed 0012 -> 0016 (D9, S04, GAP-3, sketch S7). arch 4 and qa 7 both give 0016 to D9; fresh-arch 7 gives 0019. I used the majority, 0016, and cite the front section 4 allocation table. |
| cuda 4 | major | Verified: engine "Never destroyed" [code: gpucheck.cpp:168-172]; today's memcheck has no `--leak-check` [code: checker_gpu.yaml:16-21]. Dropped `--leak-check full` (D8, Z01). Added `live_device_bytes` counter, a leak test (D5c, I19, `test_san_no_leak`). |
| cuda 13, arch 14, qa 5 | major | Removed the REPO/SP header lines and every `S[P]/` and `/[p]rivate/tmp` path. Kernel drafts now cite `docs/implementation2.0/kernels/K*`; notes cite `docs/implementation2.0/notes/`; the spec is `docs/implementation2.0/specs/D.spec.md`. D1 step 4 copies logs to `K*/results/` and `runs/d1/`; D1 step 5 fixes the GPU_RUN.md scratch paths. "SP-K1/K2/K3" renamed K1-draft/K2-draft/K3-draft. |
| cuda 14, qa 15 | major | D cost text cites the front section 6 rate table; D numbers use the measured $0.544-0.628/h and give the maximum at $0.641/h ($0.84). D-R8 changed. |
| qa 6 | major | Prefixed ids: gates D-G1..D-G17, risks D-R1..D-R13, test gaps GAP-n, GPU tests GPU01-GPU03. Research-note citations renamed note-D1 / note-D2 / note-T3 / note-CL2, so they no longer read as tasks D1, D2. |
| qa 8 | major | Verified: `pyproject.toml` has only `slow`, `data`, no `--strict-markers` [code: pyproject.toml:65-71]. D uses `gpu` only; registration and `--strict-markers` belong to the first task that uses `gpu` (Phase A); D8 adds them only if still missing. |
| qa 16 | major | Verified: `make test` has no `-I` [code: Makefile]. D8 adds `make test-data` (`-I -m data`) if no earlier phase did; `make test` deselects `data`; conventions in §9 updated. |
| cuda 15, qa 15 | minor/major | D-R13 cites the one front section 6.4 stop rule (35 min). |
| cuda 17, qa 27 | minor | No change in D: D already installs `cuda-toolkit[nvcc,cccl,nvtx]==13.0.3`. The fix is in the front. |
| cuda 18 | minor | Verified arenas are never cleared [code: gpucheck.cpp:135-150]. D5 step 1 zeroes `changed` before `CropFill`; I14 adds the stale-flag check. |
| cuda 19 | minor | Verified: `Drain` syncs both lanes directly [code: gpucheck.cpp:177-184]. Also found the review's count is low: `score_polygons` makes 4 `wait()` calls (1802, 1574, 1626, 1639) plus 2 `Drain` syncs. Added `drain_syncs`; D-G8, I08, B03 and sketch S4 now use the D1-measured baseline (today 4 and 2), not the literal 3. |
| GAP-5 (own) | minor | `label_pairs` is now in the D1 counter list, so I14 checks the skip directly. |
| fresh-arch 2 | major | No D change needed: the decode tail is a C task (C11). §8 says D has no decoding task. |
| fresh-arch 5, arch 10 | major | D.2 already keeps the K2 EDT and makes cuCIM an optional oracle. Added the revisit trigger (EDT > 10 % of eval wall time); the front row must match. |
| arch 6 | major | §8 lists the D-created modules that E's move tables need: `checker/rings.py`, `eval/device.py`. |
| arch 17 | minor | Single owner of the pool deletion is CLEAN16 step 1 (D7 step 6, §8). |
| qa 14 | major | Added a "Words" line (answer = completion, crop, rater = doctor). The code names (`completions`) stay where they are identifiers. |
| qa 28 | minor | Bare `[measured]` and `[code]` tags replaced with file tags; `[web: ...]` and `[unverified]` are listed in the header. |
| qa 17 | major | Partly done: all new text uses short sentences. A full split of older long sentences in D is not done in this revision. |

## Review log (final review round 2, revision 4)

Sources: `reviews/final_r2_{cuda,qa,fresh-arch,arch,data,rl}.md`, items addressed to D. Every code fact was checked
against the repo on 2026-11-08 (read-only). `arch`, `data` and `rl` have no item for D. Skills used this round:
nvidia-skill-finder, code-review, debug, system-design, architecture, documentation, huggingface-best (D.0 rows,
"Revision 4"). Web and package checks are listed in D.0.

| Review item | Severity | What I did |
|---|---|---|
| cuda 1 | blocking | Verified: the `nvidia-nvtx 13.0.85` wheel has `nvidia/cu13/include/nvtx3/{nvtx3.hpp,nvToolsExt.h}` and no `nvidia/nvtx` [measured: `unzip -l`]. Changed the review's fix: FindCUDAToolkit (CMake 3.26.0 and 3.31.0) creates `CUDA::nvtx3` for every CUDA >= 10.0 with no header check, so `if(NOT TARGET CUDA::nvtx3)` never fires [web: FindCUDAToolkit.cmake source]. D1 now uses `find_file(... nvtx3/nvtx3.hpp PATHS ${CUDAToolkit_INCLUDE_DIRS})` with a `FATAL_ERROR`, then links `CUDA::nvtx3`. S06 and the D.R item 12 cell updated. |
| cuda 7 | major | Verified: `cudaMallocAsync` uses the stream device's current pool; `cudaMallocFromPoolAsync` takes a pool [web: CUDA runtime API, memory pools]. D5c: `cudaMemPoolCreate` + `cudaMallocFromPoolAsync` + `cudaFreeAsync`; never `cudaDeviceSetMemPool`. Added counters `pool_reserved` and `default_pool_reserved`, box test `test_san_pool`, D-G15 rule, GAP-17. |
| cuda 8 | major | Verified: in a container `nvidia-smi` shows host PIDs [web: nvidia-docker issue 1460]. D-G17 and B06 now check (a) a constant row count, (b) workers log `torch.cuda.is_initialized()` and `contexts` = 0, (c) peak `memory.used` <= 30 GB. Added the `contexts` counter. |
| cuda 9 | major | Verified: no `cpp` path in the repo (`ls`); `make cuda` is in `Makefile:41-44`; `readme = "README.md"` [code: pyproject.toml:5]. D-G7 archive is now `cuda src CMakeLists.txt pyproject.toml Makefile README.md`, with the box commands. `make gpu-check DRY=1` builds the tar on the Mac CPU (S06, GAP-16). Added: the old venv runs the current `bench_checker.py` by path (imports only numpy and `objdel._gpucheck`). |
| cuda 10 | major (D part) | Verified: K1 GPU_RUN.md:18, K2:33, K3:15 still hold `/[p]rivate/tmp` paths [measured: grep]. Moved the path and price fix to a new D1 step 0 (Mac, $0, before V2). Added `r_d1.yaml` (runbook steps through `run.sh`) and a preflight check for `/[p]rivate/tmp` under `docs/implementation2.0/kernels/`. The front section 5 and 6.1 text is not mine to edit. |
| cuda 13 | minor | Agreed: `GC_RANGE` now uses `GC_CAT(gc_range_, __LINE__)`; each stage range sits in its own `{ }` block. S06 compiles two ranges in one scope (GAP-15). |
| cuda 16 | minor | Verified: `Block::get` grows on demand and `release` calls `cudaFreeHost` [code: gpucheck.cpp:84, 135-150]. D5c sizes pinned staging at the first call from a stated bound (`MAX_VERTICES`, `MAX_RINGS` [code: geometry/text.py:59]; int32 x, y [code: geometry/shape.py:120]). Added `pinned_grows`; D-G15 and B04 read it over the whole B3 run (GAP-18). |
| cuda 18 | minor (D part) | D8's sanitizer corpus now includes B3.1's `tests/test_scorer_sums.py` once B3.1 lands. The B text is not mine to edit. |
| cuda 19 | minor | Verified: `score_many` calls `cache.device()` per cache [code: reward/score.py:114]; `add_tables` locks, copies and waits [code: gpucheck.cpp:1737, 1743-1744]. D7b: `TRLReward.prefetch` and `_terms` build a step's cold caches in one `add_tables_from_raters` call (`device_many`). D-G8 and B03 add a cold-cache replay (GAP-19). |
| fresh-arch 7 | minor | Verified: `g_probe.make_trainer` attaches only `TraceStep` [code: experiments/phase2/g_probe.py:93]; CLEAN6 archives `Controller` and `alarms_for`. Removed the audit claim from D.0, D.1 item 9, D.3.1, D.3.2, D7 steps 5-6, D-R12, D.8 and B04. D-G14 now reads the reward path in the B3 run and `eval_sets` in `final_eval`; the audit leg returns only with CLEAN6's `--audit` branch. |
| qa 7 | major (D part) | Verified: `cov` and `test-all` collect `data` tests without `-I` [code: Makefile:22-26]. D8 now changes `cov`, `test-all` and adds `test-all-data` if no earlier phase did. D-G12: run `make cov` once with `data` deselected; if below 90 %, that value is the floor. |
| qa 13 | major | Verified against front 6.3 and back decision 3. Header, D.4, D1 step 4, D8 steps and D.8: R-D1 = V2 item 4 (optional), R-D2 = V3 item 1. D.4: D1's Mac parts gate D2; R-D1 numbers gate only speed claims; D2-D8 use the CPU oracle; if R-D1 is dropped, D-G3 in R-D2 covers D-G1/D-G2. All rows priced at $0.641/h: R-D1 $0.34-0.37, R-D2 $0.27-0.47, total $0.61-0.84 (79 x 0.641 / 60 = $0.84). D-G13 to D-G17, B04, B06 say "V3". |
| qa 21 | minor | "exam metrics" in D7 step 6 now reads "evaluation metrics". No other "exam" in D. |
| qa 25 | minor | Partly done: new and rewritten text uses one rule per sentence with "Why:" lines; split four older long sentences. Count by a rough script: 131 of 1,130 sentences over 25 words before, 130 of 1,267 after. Test-table cells are not split. |

## Review log (final review round 3, revision 5)

Sources: `reviews/final_r3_{cuda,qa,arch,data,fresh-arch,rl}.md`, items addressed to D. `arch`, `data` and `rl` have
no item for D. `fresh-arch` 7 names D.3.3 but is addressed to E (ADR 0020) and the front; no D change. Every code fact
was checked against the repo on 2026-11-08 (read-only). Skills used this round: nvidia-skill-finder, code-review,
debug, system-design, architecture, documentation, huggingface-best (D.0 rows, "Revision 5"). Web check: D.0.

| Review item | Severity | What I did |
|---|---|---|
| cuda 1 | major | Verified: `cuda/docker_check.sh:17` copies `cpp` under `bash -euo pipefail`; `ls cpp` fails. D1 file list: drop `cpp` from the copy list. D1 step 3 names the fix. D.1 item 10 states the failure. S06: every copied path has a `git ls-files` match. |
| cuda 2 | major | Verified: `CUDACXX` is exported only inside `run_box.sh` [code: run_box.sh:25-29]; `PY := .venv/bin/python` [code: Makefile:1]; the box venv is `/opt/venv` [code: run_box.sh:17]; the sanitizer sits under `nvidia/cu13/bin` [code: checker_gpu.yaml:16]. D8: new `run_box.sh --gpu-check` mode (toolkit, `make cuda PY=`, CUDA-build assert, pytest, sanitize by full path, old-engine build, benches) in one shell. The Makefile target only calls it; `d_gpu_test.yaml` calls it with `PY=/opt/venv/bin/python`. Added: `--dry` runs before the RIGA refusal check, because the Mac holds RIGA data [code: run_box.sh:18-21]. D-G7 wording and S02 updated. |
| cuda 7 | minor | Verified by web (issue 1460: "empty process panel"; the driver "does not support pid-namespaces"). D-G17 (a) and B06: record the row count first; if it is 0, count `/proc` processes that hold `/dev/nvidia[0-9]*` open (>= 1 before scoring, no rise); `runs/d/topology.json` records the method. |
| cuda 9 | minor | Rejected. The claim "lanes[1] is dead" rests on a grep for the literal `lanes[1]`. `check_polygons` indexes `e.lanes[k]`, k = 0, 1, and runs a batch of >= 64 shapes as two halves on both lanes [code: gpucheck.cpp:1060 `kSplitMin = 64`, 1671-1692]. Removing lane 1 would remove that overlap. Added a note to the engine-counter text: both lanes stay; the scorer entries use lane 0, so `Drain` syncs an idle lane 1 there. D-G8 keeps the D1-measured `drain_syncs` (today 2). |
| cuda 10 | minor | Verified: `label(1), label(0)` is gpucheck.cpp:1614; :1615 is `zero(c.aux, ...)`. D5 step 1 now says :1614, and that the zero at :1615 and `Holes` at :1616 always run. |
| cuda 12 | minor | Verified: `K3/GPU_RUN.md:7` says torch 2.14.1; the box has torch 2.13.0+cu130 [measured: runs/r5/env_check.json]; `vllm_native.txt:838` pins `torch==2.13.0`. D1 step 0 fixes the K3 line and adds a preflight grep for `torch 2.14.1`. Also found: `run_box.sh:5` says the same wrong version; D1 fixes that comment. S06 updated. |
| cuda 13 | minor | Verified: `Lane::out` holds defects (:1799), boxes (:1572), `ObjOut` + `sig_start` (:1622-1623) and signatures (:1638); `ObjOut` is 48 B [code: gpucheck.cpp:1199-1202]. (a) D5c states the `out` bound per object as the largest of these. Changed the review's signature term: a signature holds one i32 per component (holes per component) [code: gpucheck.cpp:1526-1545], not one per hole, and the grammar does not bound the component count of a rasterized object. So C_max = 64 [estimate] or twice the D1 replay maximum; a growth logs n and the component total. (b) The sanitize command uses `"$PY" -I -m pytest`; the corpus uses synthetic fixtures in `tests/`, never `data/`. |
| qa 9 | major | Verified: one `gpu` marker for all phases; D-G3 and D8 ran a bare `pytest -m gpu`. D8 defines `D_GPU_TESTS` (D's files, plus `tests/test_scorer_sums.py` after B3.1). Every box pytest call of D names it, sanitizer runs included. D-G3, D-G5, D.8 R-D2 row and S02 updated. The front section 4 rule ("never a bare `pytest -m gpu` on a box") and the A/B/C/E rental tags are other phases' text. |
| qa 12 | minor (D part) | D8 drops its own `skipif(_gpus() == 0)` and uses A.9's root `tests/conftest.py` hook (skip unless `build() == "cuda"`, no torch import). D8 adds the hook only if A has not. Kept D's "a CPU build fails on the box" rule: `run_box.sh --gpu-check` asserts the CUDA build and fails on any skipped test. "BOTH" tests now use a `host` / `gpu` parameter pair, because the hook skips a whole `gpu`-marked test. GPU01 and the §9 conventions updated. |
| qa 15 | minor | Verified: front 6.3 item 5 runs the test eval only if the owner closes the 4B gate; item 4 is "the B3 run and its val evals", and B-S2 prices "val shifted and edge". D-G14, D-G17, D7 step 5, D.8 row, B04 and B06 now read the B3 run's val eval (`final_eval` on `lidc-val-shifted`). |

## Review log (consistency pass, round 4 prep)

| Item | What changed |
|---|---|
| One source | D.8 table keeps the work minutes; its cost cells point to front 6.3 (V2 item 4, V3 item 1). D-G13 stays <= $0.85 on the billed `dph_total`, which the front offer limit caps at $0.641/h. D-R8 points to front section 6 |

## Review log (final round 4)

| Item | Severity | What changed | Verified against |
|---|---|---|---|
| cuda-4 | major | D1 file list: `d_gpu_test.yaml` no longer says "code-only sync"; it reads only the R-D2 payload. D8: new "R-D2 payload" list (LIDC bundle from V3 item 0: `checker_eval/` masks, polygons and manifest; `prepared/lidc/val.jsonl`; `prepared/caches/lidc/`; plus the tracked digest golden `tests/goldens/d_gpu02.json`, written on the Mac by `experiments/phase2/perf/make_d_goldens.py`). `BOX_INPUTS` per box test module, a `box_input()` helper that fails (not skips) on a missing file, and the rule "`gpu` or `data`, never both". Step 3 checks the payload before pytest; `d_gpu_test.yaml` passes `OBJDEL_DATA_ROOT={data}`; `--dry` runs the new S08. D1 step 2: the box replay uses target texts (saved answers do not ship). D-G9: box half dropped (host floats from `eval_sets` integers; GPU02 covers `eval_sets`), now Mac only. D-G10 box half = GPU02 cache rows (first 1,000 sorted caches; scipy on the box). GPU02: corpora, entries and golden path named. D-R01 and D-R08: Mac only, box half is GPU02. Conventions, coverage map, GAP-20 updated. Front section 6 unchanged: R-D2 stays V3 item 1, after item 0's sync | `sync_up` default `--bundle lidc` (policy.py:604); `--code-only` returns before the bundle (policy.py:401); code sync drops `data/` and `docs/`, keeps `tests/` (policy.py:379); bundle unpacked at `/workspace/data/lidc` = `{data}` (policy.py:406-417, 430-433); LIDC bundle holds `prepared/lidc`, `prepared/caches/lidc` (bundle.py:52) and `checker_eval/` (bundle.py:66-70); `checker_eval` = masks.npz (37,093) + polygons.npz (144,133) + manifest.json; 8,425 LIDC cache files, 33 MB; `tests/goldens/polycore.npz` 547 KB tracked; scipy 1.18.1 in vllm_native.txt:762; `checker_gpu.yaml` header ("stops on the missing checker_eval set") |

## Review log (final round 5)

Cross-file consistency fixer, round 5. Sources: `reviews/final_r5_{cuda,rl,data,arch,qa,fresh-arch}.md` and their
`*_verdict.md` files. Only the CONFIRMED items assigned to this fixer are listed.

| Item | Severity | Change | Evidence |
|---|---|---|---|
| cuda-1 | blocking | D1 Files and step 0: `code_files()` gains a keep rule for `docs/implementation2.0/kernels/K[123]/` (not `build/`, `results/`) with a unit test; runbooks use `/workspace/src/...` and the run folder, every scp/rsync line deleted; `r_d1.yaml` step 1 builds K1's raw masks on the box from `{data}/checker_eval/masks.npz`; new S09 box-input test (§9.6); preflight rejects `scp`/`rsync` in a runbook. Step 4 copies logs from the pulled run folder | `policy.py:379`; `bundle.py:52-70`; K1 GPU_RUN.md:51, 126; K2:33, 117; K3:16, 93-94 |
| fresh-arch-1 | blocking | Header, D1 step 4 and D.8: V2 runs only after a budget raise; R-D1 is dropped first if the raised amount is too small | front 6.3 |
| cuda-3 | major | D1 counter list: the device/host alloc counters are replaced by `sync_alloc`/`sync_free` (every `cudaMalloc`, `cudaMallocHost`, `cudaFree`, `cudaFreeHost` call, counted inside the `alloc`/`release` helpers) and `pool_allocs`/`pool_frees` (`cudaMallocFromPoolAsync`/`cudaFreeAsync`); the CPU build maps device-kind calls to the pool counters. D-G15 now states its intent (no synchronizing allocation after warm-up): box rule = steps 2-6 delta of `sync_alloc` and `sync_free` is 0 plus `pinned_grows == 0`; pool counters recorded, not gated. B04 reads the same counters. D5c: Mac test that a larger call raises `pool_allocs` and leaves `sync_alloc` at 0; the D7b batch upload stages in pageable host memory, never a pinned `Block`. I08 and sketch S4 use the new names. Front section 6 unchanged (D-G15 still read in V3, first 6 steps of the B3 run) | `alloc`/`release` are the only allocation helpers (gpucheck.cpp:78-84) and the only call sites are :146, :1741, :1746, :1772; `Block::get` grows on a larger call (:131-150; `Block<true> dev, score` at :154); each cache uploads its own tables at first `device()` (reward/cache.py:247-257); `add_tables` stages in a `std::vector` (gpucheck.cpp:1726, 1743) |

## Review log (final round 6)

Cross-file consistency fixer, round 6. Sources: `reviews/final_r6_{cuda,rl,data,arch,qa,fresh-arch}.md` and their
`*_verdict.md` files. Only the CONFIRMED items assigned to this fixer are listed.

| Item | Severity | Change | Evidence |
|---|---|---|---|
| cuda-1 | blocking | D7 step 6 converts the ceiling CLI to in-process scoring in the same commit; I18's spy covers it | `serve.py:170-173`; `reward/score.py:114` |
| cuda-5 | major | D1 step 4 tree rule: R-D1 runs only before D2 merges; on the time-line order it is dropped and D-G1/D-G2 go to D-G3 in R-D2. Header, D.4 and D.8 follow. D3 adds `K1/cross_numpy.py:16` to the oracle callers | `K1/cross_numpy.py:16`; `policy.py:490` |
| qa-4 | major | D1 step 0 S09: a `make test` test run at front 6.1 step 1, not "preflight when the plan is r_d1.yaml" | front 6.1 step 4 (plan files are never stamped) |
| arch-5 | major | D1 step 0, S06 and the review logs write the path pattern with brackets | check 14 scan: 0 hits |

## Review log (final round 7)

Cross-file consistency fixer, round 7. Sources: `reviews/final_r7_{rl,arch,qa,fresh-arch,cuda}.md` and their
`*_verdict.md` files. Only the CONFIRMED items assigned to this fixer are listed.

| Item | Severity | Change | Evidence |
|---|---|---|---|
| qa-2 | major | D8 `d_gpu_test.yaml` step 2: `bash -c` with a CUDA-build assert, then pytest on the G02b and G03 node ids; the D.8 R-D2 row notes it (A's share). S05 checks 2 steps, both `bash -c`, and collects the two node ids | `policy.py:479` (`cmd_run` dies on a bare pytest), `:486` (non-bash argv runs as `python -I`) |
| cuda-1 | major | D1 file list: `bench_checker.py` gains `--dump-replay OUT` (main venv; `pack_objects` arrays, `obj_cache`, the 13 `add_tables` arguments per cache, SHA-256; `--synthetic N` for tests) and `--replay FILE` (numpy and `_gpucheck` only). D1 step 2, D8 step 6, D-G7, D-G8 (both arms and the cold arms run on the dump; the host cache build stays out) and B03 follow. New S10 (import-closure subprocess, `make test` and `make gpu-check DRY=1`) and GAP-21. No rental change: the dump runs inside the D-G7/D-G8 bench minutes of the R-D2 row | `reward/score.py:45` imports `.cache`; `reward/cache.py:84-86` imports `cv2`, `scipy`; `Makefile:42` (`--no-deps`); `cache.py:247-253` (the 13 `add_tables` arguments); `score.py:112-115` |
| cuda-2 | major | `label_pairs` is now "crop rows relabelled by the second label pair" (D.4), not a launch count. D5 step 1: a per-object `i64 relabelled[n]` in the `score` arena, zeroed with `changed`, `atomic_add` in the second-pair `LabelInit{c, 1}` per row not skipped, copied by `ObjTerms` into a new `ObjOut` field, so the existing D2H copy before the wait at :1626 returns it (no new wait; D-G8 stays 4 + 2). I14 asserts the sum of (bh + 2) over objects with area > 0 and `changed[o]`, 0 when nothing is removed. D5c `out` bound uses `sizeof(ObjOut)` (48 B today, 56 B with this field, 72 B with B3's sums; signatures still dominate). GAP-5 updated | `gpucheck.cpp:1597-1601, 1614` (three launches over all `crows` every call); `:1583` (bh + 2 crop rows per object with area > 0); `:209-211` (`atomic_add(i64*)` in both builds); `:1199-1202` (`ObjOut`); `:1622-1626` (copy, wait) |

## Review log (final round 8)

Cross-file consistency fixer, round 8. Sources: `reviews/final_r8_{cuda,data,arch,qa}.md` and their `*_verdict.md`
files. Only the CONFIRMED items assigned to this fixer are listed. Each fact was checked in the repo first.

| Item | Severity | Change | Evidence |
|---|---|---|---|
| cuda-1 | blocking | D8 `--gpu-check`: `OUT="${OBJDEL_RUN_DIR:?}/d"` (`--dry`: `mktemp -d`); sanitizer JSON, JUnit XML, replay dump, `bench_<engine>.json`, `replay_<engine>.json` under `$OUT`; `--out` required; `--compare "$OUT"` prints the D-G7/D-G8 table and exits 1 on a fail; D1 step 2, D-G5, D-G7, D-G8, B03, Z01 and GPU02 paths; new S11 | `run_box.sh:16`; `policy.py:489, 505-507`; `.gitignore:23` |
| cuda-2, qa-3 | major | D7 step 5 names the code: `FastGRPO.log` writes `engine/*` (every `gc.stats()` key, `max_memory_reserved`) on steps 1-6 and `engine/pinned_grows` at the end into `rl.json`; `final_eval` writes `final_eval.json["engine"]` and runs the 5 s sampler thread around the `item_counts` loop, writing `$OBJDEL_RUN_DIR/d/topology.json`; D-G14/15/17, B04, B06, D.8 updated; new S12 (fakes) and S13 (Mac read after `pull.sh`) | `fast.py:604`; `train_rl.py:63` |
| arch-3 | major | D7 step 6 also deletes `--procs` from `select_ckpt` and `g_probe` and `"--procs", "1"` from the preflight dry run; one pool rule (D7 step 6, E3b R-2, CLEAN16 check); D.8 CLEAN line | `select_ckpt.py:374, 390`; `g_probe.py:139, 157`; `final_eval.py:69`; `test_p2_scripts.py:369, 384` |
| qa-6 | major | D8 `d_gpu_test.yaml` step 2 inputs: the golden holds G02b's item ids and texts, so no saved answers are needed | A.9.6 G02b |
| cuda-3 | major | D-G7: old tree = owner tag `pre-d2` (owner command `git tag pre-d2 main` at D2 start; fallback `<D2 merge>^1` or `<first D2 commit>^`); new Makefile target `old-engine-tar` writes `old_engine.tar` + `.sha256` at the REPO ROOT (fails without the tag); D8 step 0 checks both files before the toolkit install; step 5 = `sha256sum -c`, `tar -xf old_engine.tar -C /opt/old_src`, `uv venv --python "$PY" /opt/old`, numpy under the image lock, `make -C /opt/old_src cuda`; `--dry` checks the tag, the files, `code_files()` and a CPU build of the extracted tree in `mktemp -d` (skips with a message before the tag exists, so D1 step 3 still passes); D8 Steps: `make old-engine-tar` + `make gpu-check DRY=1` on the Mac before the V3 launch; D9 runbook deletes both files after teardown; S06 reworded, new S14, GAP-22, coverage map. No rental change: extract and venv sit inside the R-D2 row's "old + new about 6" build minutes (front 6.3 V3 item 1 unchanged). Residual, not changed here: preflight runs `make gpu-check DRY=1` only for `d_gpu_test.yaml`, which never follows `preflight CONFIG=` (front 6.2 item 4), so the Mac step is a D8 Steps rule, and step 0 is the box guard | `policy.py:377-380` (sync list and drop regex), `:398` (rsync to `/workspace/src`); `.gitignore:16, 23` (`/build/`, `/runs/`); `git check-ignore old_engine.tar` exit 1 [measured 2026-11-09]; `run_box.sh:15-16, 22-23`; `Makefile:13-17, 41-45`; `bundle.py:113`; `git archive HEAD` of the six paths = 563,200 B [measured] |

## Review log (final round 9)

Cross-file consistency fixer, round 9. Sources: `reviews/final_r9_{cuda,rl,data,arch,qa,fresh-arch}.md` and their
`*_verdict.md` files. Only the CONFIRMED items assigned to this fixer are listed. Each fact was checked in the repo first.

| Item | Severity | Change | Evidence |
|---|---|---|---|
| cuda-1 | blocking | D1 file list: `bench_checker.py` resolves every input under `ROOT = Path(os.environ.get("OBJDEL_DATA_ROOT") or "data")` (`checker_eval/masks.npz`, `prepared/lidc/val.jsonl`, `prepared/caches/lidc/`), has no hard-coded `data/` path, checks its inputs first and exits 1 naming the missing full path, writes the root, resolved paths and SHA-256 into the JSON, and passes its environment to the old-engine subprocess. D1 step 2: inputs written relative to the data root. D1 step 0: `r_d1.yaml`'s bench step sets `OBJDEL_DATA_ROOT={data}` and `--out {rd}/d1/bench.json`. D8 step 0: `: "${OBJDEL_DATA_ROOT:?}"` guard (not dry) before the toolkit install. D8 step 6: every bench call inherits `OBJDEL_DATA_ROOT={data}` from `d_gpu_test.yaml` step 1; step 3 already checks the three inputs. New S15 ($0, no marker, reads no `data/`, runs in `make test` and `--dry`), GAP-23, coverage map. No rental change: front section 6 untouched (R-D2 stays V3 item 1 after item 0's LIDC sync) | `run_box.sh:15-16` (`set -euo pipefail`, `cd` to repo root); `policy.py:492` (`cd /workspace/src` in `cmd_run`); `policy.py:379` (drop `data/`); `policy.py:408, 417, 433` (bundle unpacked at `/workspace/data/lidc` = `{data}`); `bundle.py:52, 66-70` (bundle holds `prepared/lidc`, `prepared/caches/lidc`, `checker_eval/`); `train_sft.py:82` (same `OBJDEL_DATA_ROOT` pattern) |
| cuda-2 | major | D1 file list: `bench_checker.py` keeps `from objdel import _gpucheck` after E (E3a TID251 exception); S10 gains a third subprocess under a fake pre-E tree (`objdel/__init__.py` + `_gpucheck` stub, no `objdel.geom`) | E.md E3a; Makefile:42 (`--no-deps` old venv) |
| qa-1 | major | D8 rules: `box_input()` also serves B's G08b smoke in `tests/test_train_fast.py` (a `bash -c` step with `OBJDEL_DATA_ROOT={data}`); S08 checks that module's `BOX_INPUTS` once B3.9 lands | `policy.py:379, 486, 492` |
| qa-2 | major | D8 `D_GPU_TESTS`: `tests/test_scorer_sums.py` is the one G01 file (np.nonzero sums + Mac golden `tests/goldens/b_g01.json`, never two builds in one venv) | `cuda/gpucheck_py.cpp:113`; `Makefile:41-45` |
| qa-6 | major | S13: run folder = `runs/$OBJDEL_V3_RUN`; without the variable both tests skip (A.9.6 rule), so `make test-data` stays green under (b'); with it a missing file fails; `test_b04` reads `eval_b3_val_shifted/final_eval.json` (B3.9 names it) | A.9.6 skip rule; front 6.3 V3 note |

## Review log (final round 10)

Cross-file consistency fixer, round 10. Sources: `reviews/final_r10_{cuda,data,qa,fresh-arch}.md` and their
`*_verdict.md` files. Only the CONFIRMED items assigned to this fixer are listed. Each fact was checked in the repo first.

| Item | Severity | Change | Evidence |
|---|---|---|---|
| cuda-1 | blocking | D8: `d_gpu_test.yaml` is `specs/d_gpu_test.yaml` after E6 (front 6.3 "After E6", E6 step 2); preflight sentence names both paths | `policy.py:471`; E I29 |
| cuda-3 | major | D8 `--dry`: the `code_files()` assert line is rewritten by E8 step 1 to `objdel.infra.vast.remote.code_files()`; E's U8 fails on a missed rewrite | E8 name table (`code_files` -> `infra/vast/remote.py`) |

## Review log (final round 11)

| Item | Severity | Change | Evidence |
|---|---|---|---|
| fresh-arch-3 | major | D3 now extends B1.4a's `rings.py` (adds `trace`, `MaskRings`, `ring.part`), exports only the new names, extends `test_rings.py`, and re-runs B1.4a's `verify_target_hashes.py`; step 3 time uses the measured `check_masks` 4.83 s | K1/results/bench_cpu_m3.txt |

## Review log (final round 12)

| Item | Severity | Change | Evidence |
|---|---|---|---|
| qa-1 | major | D2: the fuzz and 37,093-crop GPU tests compare with `tests.oracles.polygonize_oracle` through `compare.mask_reference`, never with `objdel.checker.polygonize`; D3 moves only `bitquad.py` and its callers; D-G3 row, D9, the E paragraph and the test conventions name B1.4a as the polygonizer oracle's mover | `compare.py:21, 167`; B.md B1.4a |
| fresh-arch-1 | blocking | D1 Files and step 0: widen B1.6's keep regex (`frozen_views.json` plus `K[123]/`); the unit test allows `docs/reference/frozen_views.json` | B.md B1.6; `policy.py:379` |
| cuda-1 | major | D2 step 1: `Lane` gains a third device arena `Block<true> scan` (CUB scratch only) and the K1-draft map is stated (`devA` -> `L.dev`, `devB` -> `L.score`, `tmp` -> `L.scan`, `in` -> `L.in`, `out` -> `L.out`; never a live block as scan scratch). D2 step 7: test-only CPU define `GC_POISON_SCAN` (CPU `scan` fills `tmp.get(4 * (n + 1))` with 0xFF); `make test-order` runs `-k masks` a fourth time under it, so an aliased map fails on the Mac at $0. D1 `pool_allocs` note and D5c pool path name `dev`, `score` and `scan`. Front section 6 unchanged | K1/k1poly.cpp:140-144 (CUB scratch; CPU build ignores `tmp`), :162 (`devA, devB, tmp`), :824-825, 871-872, 876, 884 (six scans while A and B live); gpucheck.cpp:139 (`get` returns the same base), :154 (`dev, score`) |

## Review log (final round 13)

| Item | Severity | Change | Evidence |
|---|---|---|---|
| cuda-1 | major | D5c: each lane's pinned `in` and `out` have ONE bound each (the `score_polygons` input and output bounds), set once inside the first engine call of the process, whatever the entry; `score_polygons`, `check_polygons`, `check_masks`, `score_masks` stage through them; `add_tables_from_raters` and `eval_sets` copy through pageable host memory and never touch them; the "`eval_sets`: its scratch cap" bullet is gone. D7b: the returned scalars go to a `std::vector` after the batch's one wait, never `Lane::out`. D.4 `pinned_grows` text matches. New I20 (fresh subprocess, D7b batch first, then `score_polygons`, then 5 cold rounds: `pinned_grows == 0`, sync deltas 0; reverse order too; a control growth adds 1), GAP-24, D-G15 Mac half and the coverage map name it. Front section 6 unchanged | every D2H result uses `L.out` today [code: gpucheck.cpp:991, 1572, 1623, 1638, 1681, 1799]; `Block::get` 4,096 B minimum, free + alloc on growth [code: gpucheck.cpp:139-149]; `Lane` = `Block<true> dev, score; Block<false> in, out` [code: gpucheck.cpp:154-155]; D2H to pageable "might be synchronous with respect to host" [web: docs.nvidia.com/cuda/cuda-runtime-api/api-sync-behavior.html] |
| qa-1 | major | D2 Files: `CMakeLists.txt` cache string `OBJDEL_TEST_HOOKS`, binding `test_hooks()`, helper `tests/variants.py`. D2 step 7: each variant is a NON-editable `uv pip install --target build/variants/<name>/site` and runs in a `-I -S` child with `sys.path` set inside it (as CLEAN1 step 5); load guard (`__file__` under the variant folder, `test_hooks()` equals its defines, exit 3 otherwise; normal suite asserts `()`); planted `GC_ALIAS_SCAN_FOR_TEST` variant `alias` must make the poison run fail on every `make test-order` (red first recorded in D2's note). D8: Makefile target `test-order` and the CI lint + pytest step. D-G3 (Mac half), D-G12, D-G16 name it; new S16, GAP-25; D.9.9 item 1, sketch S3 and the `mutant_dir` fixture drop `PYTHONPATH` | editable finder inserted at `sys.meta_path[0]` with `{'objdel._gpucheck': 'objdel/_gpucheck.abi3.so'}` [code: .venv/lib/python3.12/site-packages/_editable_skbc_objdel.py:742, 805]; probe 2026-11-09: `PYTHONPATH` and `-I` + `sys.path.insert` load the installed `.so`, `-I -S` + `sys.path` loads the variant folder and `-k masks` passes there (3 passed); non-editable `--target` build writes the full package, 3.8 s wall [measured: scratch probe]; Makefile has no `test-order` [code: Makefile:6] |

## Review log (final round 14)

| Item | Severity | Change | Evidence |
|---|---|---|---|
| qa-1 | blocking | D-R01 always-run copy is LIDC only: 50 LIDC items in `tests/goldens/d_r01_lidc.json` (run-length input masks plus pre-D text), written by `make_d_goldens.py` on the Mac; the test asserts every item id starts with `lidc/` and none contains `riga`; RIGA items are rebuilt only in the Mac `data` half and never written to `tests/`. D8 `make_d_goldens.py` text names the new golden and its `lidc/` refusal. D-R02: "RIGA disc/cup/rim" becomes synthetic disc/cup/rim shapes (as the conftest generator), no real RIGA item. D-G4 (Mac, LIDC and RIGA hashes) unchanged. Front section 6 unchanged (no box or rental change) | owner rule "RIGA-derived results stay local: no permission to share RIGA data (repo public, owner 2026-11-06)" [code: .gitignore:37]; code sync keeps `tests/` and drops only RIGA names matching `riga[^/]*\.(csv|jsonl?)$` [code: experiments/phase2/vast/policy.py:379]; `test_code_snapshot_has_no_riga_derived_results` checks file names only [code: tests/test_p2_breakit_ops.py:297-300]; today `tests/goldens/` holds only `polycore.npz`; synthetic RIGA items `riga/MESSIDOR/image00<s>` [code: tests/conftest.py:44-56]; real ids are `lidc/...` and `riga/...` [code: src/objdel/data/lidc.py:316, src/objdel/data/riga.py:332] |

## Review log (final round 15)

| Item | Severity | Change | Evidence |
|---|---|---|---|
| data-2 | major | D1 Files line and step 0: the widened `keep` regex starts with `^docs/reference/(frozen_views|val_halves)\.json$`; the step-0 unit test allows both JSON files | `policy.py:379-380` (drop `^(data|docs)/`, no keep today); `docs/reference/` holds no JSON today [measured: ls]; `policy.py:398-399` (editable install at `/workspace/src`); B.md B1.6 |
| cuda-1 | major | D.4 "Engine counters" defines `arena_bytes` (sum of `cap` of every `Block<true>` of both lanes: `dev`, `score`, `scan`), `table_bytes` (live table blocks in the registry; new counter) and `eval_scratch_peak` (largest scratch one `eval_sets` chunk asked for), and `live_device_bytes == arena_bytes + table_bytes` after D5c; `reset_stats()` zeroes event counters and `eval_scratch_peak` only. D7 Files: `eval_sets` takes its chunk scratch from lane 0's `score` arena and raises `eval_scratch_peak`. D7 step 5 and S12: `final_eval` calls `reset_stats()` before each view and records `eval_sets_s`, `max_scratch_bytes`, `arena_bytes`, `table_bytes`, `eval_scratch_peak`, `live_device_bytes`, `max_memory_reserved`. D-G14 eval half, S13 `test_b04` and B04 regate on `eval_scratch_peak <= max_scratch_bytes` and `live_device_bytes <= arena_bytes + table_bytes`; `arena_bytes`, `table_bytes` recorded, not gated. $0 checks of the new counters: P03 (`eval_scratch_peak` <= each `eval_sets` cap) and I19 (`live_device_bytes == arena_bytes + table_bytes`, `table_bytes` back to start). RL half unchanged (it logs every `gc.stats()` key, so the new keys appear without a code change). D.spec.md G14 still carries the old audit-era rule; it is not an assembled source. Front section 6 unchanged (no box or rental change) | `Lane` = `Block<true> dev, score; Block<false> in, out`, `Block::get` rounds `cap` up [code: cuda/gpucheck.cpp:135-156]; `i64 resident()` returns a count [code: cuda/gpucheck.cpp:1778]; `metrics()` keeps `score_many` [code: src/objdel/eval/metrics.py:250-258], which calls `_gpucheck.score_polygons` [code: src/objdel/reward/score.py:115]; `score_polygons` sizes `L.dev` [code: gpucheck.cpp:1799], `score_words` sizes `L.score` [code: gpucheck.cpp:1565, 1592] |

## Review log (final round 16)

| Item | Severity | Change | Evidence |
|---|---|---|---|
| cuda-1 | major | D8 sanitizer corpus: only `tests/test_scorer_sums.py::test_g07_sanitize_sums` (markers `gpu` and `sanitizer`, the G01 subset defined in B.9.5 G07) joins after B3.1. `test_g01_golden` carries `gpu` only and runs in `--gpu-check` step 3, never under a sanitizer. D8's "reads no `data/`" sentence names the one exception (32 crops of `checker_eval/masks.npz` through `box_input`). D-G5 names the one test. Front section 6 unchanged (no box or rental change) | `ExpertCache.build` with 3 raters calls `edt` 6 times in `rater_tau` and 6 times in `f5_tolerance` [code: src/objdel/reward/cache.py:168-197]; full G01 = 10,881 caches, 100.9 s Mac host work [measured: B.9.5 G01]; subset = 256 caches, about 3,300 engine calls [arithmetic] |


# Phase E

## E.R. Review log (revision 1, 2026-11-08)

Source: docs/implementation2.0/reviews/E_phase_review.md (32 items). I checked each claim in the code at be52b1c before acting.
"Done" means the plan text changed as the review asked; "Done, changed" means I did it a different way and say why.

| # | Severity | Item | What I did |
|---|---|---|---|
| 1 | blocking | `prompts`, `records`, `cache_loader` in `rl/rows.py` break `rl \| sft` and the layer order | Done, changed. Verified select_ckpt.py:184, final_eval.py:70-71. New `policy/rows.py` (prompts, records, _verified, _image, _check, _caches, _voting, cache_loader, DISAGREEMENT_OVERSAMPLE). `rl/rows.py` keeps prompt_rows, grpo_rows, grpo_rows_views, OpenOnce. Correction to the review: `_ViewRows` (and `_Render`) serve only `dataset.view_rows` [code: dataset.py:392], so they go to `sft/rows.py`, not rl. All 29 grep-found names placed; `tests/core/test_move_complete.py` (E-G16). Noted the two different `view_rows` functions; the select_ckpt one is renamed `eval_rows` in its own commit. |
| 2 | blocking | `verify/parity.py` from compare.py makes src depend on tests | Done. compare.py -> `tests/oracles/cuda_parity.py`; AST test bans `tests`/`oracles`/`test_*` imports and `sys.path` changes in src. 5 test users verified. |
| 3 | blocking | `make lint` cannot pass with the planned TID251 exceptions | Done. Verified test_train_fast.py:427,483 and the `_gpucheck` importers. Count correction: 13 files under tests/ in all (11 test files + 2 oracles), not 13 + 2. E3a rewrites all of them plus compare.py; per-file-ignores for test_train_fast.py and the adapter tests; "tests never import experiments" is an AST test; E0 pastes a 0-error `ruff check`. |
| 4 | blocking | E1.5 breaks default `make preflight`; no task makes preflight read a spec | Done. Verified preflight.py:49-56, rl_shifted.yaml:17, policy.py:315. E1.5 default -> `r1_qwen3vl2b.yaml`. New E6.8: preflight takes a RunSpec, stamp keyed on spec path + sha, `refusal(spec_path)`, tests; default becomes `specs/rl_shifted.yaml` there. |
| 5 | major | Deleting G_PROBE LR breaks the pre-registration rule and a test | Done. Verified prereg.py:1-4 and test_p2_breakit_runners.py:60. G_PROBE unchanged; `rl_lr(model)` in configs; `tests/experiments/test_prereg_lr.py`. |
| 6 | major | Golden 4 cannot see an LR change | Done. Golden 6 (resolved `GRPOConfig.to_dict()` for the train_rl stage and a G-probe); golden 4 adds the parameter-change L2 norm. |
| 7 | major | Bit-equal goldens fail across CI platforms | Done. Verified pyproject.toml:8-9 and the 3-OS matrix. Golden 1 = value hash without versions (full fingerprint only when versions match); goldens 3 and 4 per platform; Mac binding; README holds the recording script. |
| 8 | major | Hard-coded test paths; CI never runs lint-imports | Done. Verified pyproject.toml:58, tests.yml:31,54, preflight.py:29-34; found one more config (`r2_qwen3vl2b.yaml`). E0 adds the CI lint-imports step; E8 adds `tests/infra/test_test_paths.py` and the gate. |
| 9 | major | `__file__.parents[n]` roots break after the move | Done. Verified secrets.py:14, policy.py:46-47; found preflight.py:26 too. `infra.paths.repo_root()` + 3 tests. |
| 10 | major | SamplingParams drops 4 fields; carry-over list missing; two vLLM paths; venv | Done. Verified configs.py:108-114, sample_vllm.py:17,34-101, select_ckpt.py:55,88-95. `params_for` returns a dict; carry-over list in E4.1; child entry `python -I -m objdel.policy.sampler`; sample_vllm.py and hf_sample.py archived in the same commit; vllm_native only (ADR 0020, IMPORT_CHECK). Spec updated (SamplingParams dataclass removed). |
| 11 | major | "CI contains 0" is weak; E4 adds an eval that section 8 denied | Done. Deterministic gate (prompt ids, SamplingParams fields) binding; TOST ±0.01 with n = 1,084, K 8, 90 % CI. Section 8 now counts one extra eval: 415 s measured (better than the review's 9 min estimate) + ~2 min engine start [estimate], ~$0.09. |
| 12 | major | No SelectStage/ProbeStage; g_probe would need a fake RlStage | Done. SelectStage, ProbeStage added; `build_trainer` takes explicit typed arguments; `TrainerFlags` base shared by RlStage and ProbeStage; r1_qwen3vl2b added to golden 5. Added: its train stage is compared by value, because today's argv names the yaml file that is archived. |
| 13 | major | `ComputeSpec.regions` is a second source against the owner rule | Done. Field removed; `REGION` constant is the one source; test that a `regions` key fails. |
| 14 | major | Module names in strings (serve.SCORERS) escape the import rewrite | Done. Verified serve.py:52,170-171. SCORERS removed; callers pass `modules=`. `test_module_strings_resolve` in E0, made static (no imports) so it runs in the CI job without torch; it also resolves attribute paths such as `objdel.train.dataset.records`. R-1 now names string rewrites. |
| 15 | major | k3_harness hides a dependency on tests | Done, second option. `digest`/`make_jobs` need `H.G` and `H.jitter` from the oracle [code: reward_gates.py:237-249], so they cannot go to verify. k3_harness (only plan: r0.yaml) moves to `tests/harness/`; code_files guard now checks `tests/test_gpu_scoring.py`; `-X importtime` test on all live CLIs. |
| 16 | major | import-linter pin 2.11 is stale | Done. Verified on PyPI: 2.16 (2026-11-06). Pin 2.16; `broken_contract_guidance` on each contract; negative tests rerun on 2.16; section 2 corrected. |
| 17 | major | The E0 contract checks nothing | Done, changed. Interim contracts on current packages from a fresh AST graph: `["train", "eval", "reward \| targets", "data", "checker \| geometry"]` (the review's list left out `data`, which imports checker and is imported by targets) with the one ignore `eval.vision -> train.dataset`. `pure-core` names only existing packages (2.16 rejects missing local modules and empty configs). `(x) \| (y)` tested on a scratch package. |
| 18 | minor | `import *` fails ruff F403/F405 | Done. Explicit names + `__all__`; grep found 10 names in use (adds `score_masks`, `resident`, `add_tables`, `drop_tables`); a test checks coverage. |
| 19 | minor | The grep gate cannot measure "one place" | Done. Replaced by the prereg test + `test_g_probe_lr_not_read`. |
| 20 | minor | A second GRPOConfig to read one number | Done (`rl_lr`). |
| 21 | minor | Copied `_prepare_inputs` reads unlisted private state; eval branch | Done. Verified grpo_trainer.py:991,995,998,1582 and fast.py:553-556. `_step`, `_buffered_inputs`, `num_iterations` in REQUIRED and PRIVATE; eval branch and fallback guard kept; eval-mode equivalence. |
| 22 | minor | Private-name scan flags our own classes | Done. Scan only `self.` in TRL subclasses, `trainer.`, `super()`, TRL module aliases. |
| 23 | minor | Host-sync regex counts numpy calls | Done. Per-file count in torch-importing modules; `host_syncs.json` with a move map. Spec B11 updated. |
| 24 | minor | `save_steps` default differs; dict field breaks hash | Done. `save_steps: int = Field(50, gt=0)` + `<= steps` rule; `grpo: Mapping[str, JsonValue]`; stages never hashed. |
| 25 | minor | Exact pydantic pin; phase2 image has no pydantic | Done. Verified train.txt (0) and `--no-deps`. Range in base deps; `import pydantic` in IMPORT_CHECK and image_check.sh; vllm_native in ADR 0018; owner decision 5 (archive phase2 image). |
| 26 | minor | Use pydantic-settings instead of a hand-written CLI | Done, with a spike. Verified the lock (2.16.2) and the docs (CliApp.run accepts BaseModel and `cli_args`). Half-day spike must round-trip golden-5 argv; fallback is a generated argparse. Stage models are flat (TrainerFlags inheritance) so flags stay `--shared-prefix`, not `--flags.shared-prefix`. |
| 27 | minor | serve.py is more than the pool | Done. Split into `verify/store.py` (stays) and `verify/pool.py` (Phase D deletes). |
| 28 | minor | `git ls-files -co` ships `.claude/` to boxes | Done. Drop pattern adds `^\.claude/\|^archive/`; test. |
| 29 | minor | Skipped HTTP test lowers coverage below the 90 % gate | Done. Train job installs `.[dev,data,train,serve]` (E7.5). |
| 30 | minor | env_check line would not write env_check.json | Done. `@check("trl_private_names")` inside the heredoc; E-G13 reads that key. |
| 31 | minor | Polygonizer seam for Phase D | Done. `geom.polygonize(mask)` seam in E3a; section 3.3 and section 5 updated. |
| 32 | minor | infra code inside the wheel | Done, as an owner decision. Verified no workflow publishes the wheel (CI uploads an artifact only). Default: keep infra in the package; exclude from the wheel only if one is published, after checking that the editable install on the box still sees it (ADR 0017, owner decision 6). |

Own findings while verifying (not in the review):
- Golden 5 strings change at E8, because module paths move (`experiments.phase2.train_rl` -> `experiments.rl.train_rl`).
  Golden 5 is now compared after a recorded module rename map.
- The E3b text said the `eval <-> train` cycle ends in E3. It ends in E5, when `vision` leaves `eval`; the interim
  contract ignores that one edge by name until then.
- E0 contracts must be edited in each move commit (package names change); E3 and E5 now say so.


---

## Review log (final review round 1, 2026-11-08)

Source: `docs/implementation2.0/reviews/final_r1_{arch,fresh-arch,cuda,qa,data,rl}.md`, every item that names E. Each
claim was checked in the code at be52b1c or in the current A-D, CLEAN and DOC texts before the change. "Done" = the
text changed as the item asked; "Done, changed" = done another way, with the reason; "No E change" = the item
addresses E only by reference, and another phase owns the edit. The `rl` review has no E item. Rows 27 and 31 of the
revision-1 log (E.R) are superseded by arch items 3 and 17 below; the ADR numbers in E.R were renumbered with the rest
of the file.

| Review, item | Severity | Item | What I did |
|---|---|---|---|
| arch 1, QA 3 | blocking | Three phase orders; E said "E5 before C" | Done. Adopted the front order: A-D and CLEAN wave 1 first, E0 after the last A-D commit that edits `fast.py`, `dataset.py`, `reward/`, `checker/`. New "When E starts" paragraph (E.4); E.7 row 3, R-4 and the E.8 dependency bullets rewritten; B2/B3 land on pre-E paths if written before E, else on post-E paths with the move tables as the map (reconciles QA's "B2/B3 in the E tree"). Superseded in round 2 (QA item 26): B2/B3 land before E on pre-E paths, no post-E exception. |
| arch 2 | blocking | A6's NVTX module breaks the E0 contracts | Done. `train/nvtx.py` -> `policy/marks.py` (E5 row); `marks.py` removed from `core` (E.3.1 box and text); E0 states the A6 precondition (`checker` range in `train/rewards.py`) and stops if the E0 graph shows `reward -> train`; E0 negative test (d) (function-level `import torch`). grimp docs list only `TYPE_CHECKING` as skipped [web], so the function-level case is marked [estimate] until test (d) runs. |
| QA 20 | minor | NVTX module: use `core/marks.py` everywhere | Done, changed: one home, but `policy/marks.py`, not `core/` (arch item 2's reason: lazy torch import in a `pure-core` source). Two stores explained in E.3.1. |
| arch 3, CUDA 23 | blocking / minor | Three polygonizer seams; seam order | Done. D3's design only: E3a rows `checker/rings.py -> geom/rings.py`, `certify/topology/volume -> geom/`; rows for `polygonize.py` and `bitquad.py` deleted (already oracles after D3); E3a "change commit" and E.3.3 paragraph 2 replaced ("D3 lands first on today's paths; E3a moves `checker.rings` unchanged"); E.5 K1 bullet rewritten; E.8 D bullet. CLEAN11 text is CLEAN's (noted in E.8). |
| arch 4, fresh-arch 7, QA 7 | blocking / major | ADR numbers collide | Done, with arch item 4's table (the most complete; fresh-arch and QA propose different tables, so the front section 4 table is the one source and wins if it differs): E's ADRs are 0017 layers + import-linter + wheel, 0018 typed run spec + CliApp + vllm_native image, 0019 TRL adapter, 0020 verifier API + sampler on vllm_native. Every `ADR 0012/0013/0014/0015` in E renumbered (E2, E4, E6.0, E7, E8.1, E8.5, E.8, U43, E.R). |
| arch 5 | major | E1.3 and CLEAN2 differ (`gc.freeze()` lost) | Done. Verified `gc.freeze()` at serve.py:172. E1.3 = CLEAN2; E1 only checks `--procs` gone and one `gc.freeze()` after the last preload; U17 updated. |
| arch 6 | major | No move rows for A-D / CLEAN modules | Done. Rows for `eval/{ceiling,select,device}.py -> verify/`, `geometry/box.py -> geom/`, `checker/rings.py -> geom/`, `train/{gpu_pixels -> policy, vllm_sync/vllm_requests -> rl, nvtx -> policy/marks, trl_adapter -> rl}`. `old_names.json`, goldens, host-sync baseline and the package graph recorded on the E0 commit; "a module created after E0 and before E5 adds its row in the same commit" (E5); U10 and E-G16 list the A-D modules. Added a guard for `ceiling`/`select`: if one imports `data`, its row moves to `data/`. |
| arch 7, fresh-arch 8 | major | `eval/coverage.py` has two homes; stores unclear | Done. `verify/coverage.py` row deleted; E.3.1 data box has `coverage, windows, quality`; E3c lists B1.3's files as "stay in place"; `data/store.py` (source store) vs `verify/store.py` (doctor-store tables) stated once. |
| arch 8 | minor | RLOO advantages: verify or rl | Done. `advantages_v5` stays in `verify/reward/advantages.py`; `rl/reward.py` only calls it (E3b note, E.8). |
| arch 9 | major | Sampler seam signature differs (front) | Done for the E part: E.3.2 diagram now says `sample(rows, params) -> list[Sample]` as E4 defines. The front 2.2-2.5 text is the front's. |
| arch 10 | major | `data_revision` row has no field | Done (option 1). `RunSpec.data_revision: str | None` (40 hex when set), written to `run_spec.json`, checked by preflight against the rows-root manifest; test in E6.7 and U21. Superseded in round 2 (arch item 8): no `data_revision`; `RlStage.rows_revision` + `rows_manifest_sha256` (B1.4/B1.5 names). |
| arch 11 | blocking | E8.2 default breaks "experiments/<phase>/" | Done. Default is now "no move" (E8.2, owner decision 2, E.7 row, E.3.1/E.3.2 paths, shell shims path); layer folders only after the owner edits the CLAUDE.md rule; `LIVE_CLIS` constant feeds golden 5's rename map, E-G15, `test_cli_smoke`, `test_cli_imports`, U38. |
| arch 12 | major | Nobody archives the 12 old YAML plans; wrong CLEAN refs | Done. E6.6 archives them itself, names listed (16 at be52b1c [measured: ls] minus the 4 ported), README rule, CLEAN8's LFM configs in the same commit; finished A-D plans join them. E8.2 cites "CLEAN items 5 and 7"; E.8 and E8.2 say E6.6 archives `r0.yaml`. |
| arch 14, QA 5, CUDA 13 | major | Doc depends on the scratchpad | Done in E. Header REPO/SP lines deleted; "paths are relative to the repo root"; every scratchpad path became `docs/implementation2.0/{notes,specs,reviews,kernels,tests}/...` (arch item 14's layout; QA proposes `docs/research/notes/impl_2026-11-08/` for notes: the DOC copy task decides, and a later search-and-replace follows it). 0 absolute local paths and 0 scratchpad prefixes left in this file. |
| arch 15 | minor | "TRL 1.16" text has three owners | Done, changed. DOC2 owns the text. Verified DOC2's list does not include the `fast.py:672` comment "grpo_trainer.py 1844", so E1.4 keeps only that one line (-> `_generate_single_turn` at 1832 [measured: grep of the installed TRL]); U18 checks only it. |
| arch 16 | major | Test locations disagree; no import-mode rule | Done. `tests/<layer>/` from E0; `tests/phase_e/` removed (E.9 conventions, sketch 1 path, fixtures in root conftest); unique-basename rule and `tests/infra/test_unique_basenames.py` (U47); `tests/conftest.py` stays at the root so `from oracles import` resolves (13 uses [measured]). |
| arch 17 | minor | Duplicate ownership | Done. E0's `test_src_never_imports_tests` deleted (CLEAN12's test, extended in E0 if needed; U5, E-G8); E8.1 drop pattern -> CLEAN1; E1.5 -> CLEAN9.3; pool deletion -> CLEAN16 step 1 after D7 (E.3.1, E.3.3, E3b, E.5, E.8). |
| arch 18 | minor | C names E homes E does not define | No E change: E now defines `policy/caches.py`, `policy/gpu_pixels.py`, `rl/vllm_sync.py`, `rl/vllm_requests.py`, the names the item asks C to use. |
| arch 19 | minor | CLEAN1's bundle test hard-codes module names | Done. E8.1 adds `LIVE_CLIS` beside `IMPORT_CHECK` and rewrites `tests/test_vast.py::test_bundle_imports_standalone` to read both. |
| arch 20 | minor | `dataset.py` name count 31 vs 29 | Done. AST count 29 [measured, this session]; E.0 revision-1 row corrected. |
| arch 21 | minor | 99_back labels | No E change (99_back owns it). |
| arch 22 | minor | No build home for K3 | Done (round 2, CUDA item 20, replaced "with the `train` extra" by a box build against the image torch). E.5: `cuda/k3/`, extension `objdel._k3`, imported only by `policy/gpu_pixels.py`, TID251 ban with one exception, `_k3` in `exhaustive_ignores`; K3b has no task. E0 ban note. |
| fresh-arch 6 | major | Front parts with no task (`linear_attention` guard -> E5; `uv.lock`/pyright/`ruff format` -> "CLEAN9 or E0") | Done for the guard: E5.2 adds it in `FastForward._get_per_token_logps_and_entropies` per M2 rec 5 [notes: M2.md:109], test U46. Not taken for `uv.lock`, pyright, `ruff format`: the front's "owning task" column decides; E0 adds only import-linter. |
| CUDA 3 | blocking | Device eval in forked workers; nobody deletes the pool | Done for E: one owner (CLEAN16 step 1, after D7, per arch item 17); E3b archives the pool part if no caller is left on the E0 commit. The process topology fix (D7 step 6, gate D-G17) is D's. |
| QA 6 | major | One id means different things | Done for E: gates are E-G1 to E-G16 everywhere in this file (one search-and-replace; A-G7 references kept). Test ids U/P/I/PF/GP/X are local to E.9. |
| QA 9 | major | Test-count gates use the 2026-11-08 baseline | Done. E0 writes `runs/E/baseline.txt` on the E0 commit; R-3, E-G1, I30 read "passed >= E baseline - removed + added"; 1447 kept as a be52b1c fact only. |
| QA 10 | major | Open behaviours with either-way tests | Done for E. E.9.10 is now "Decided rules" (7 rules); U21 (`schema_version` strict, `stages: []` raises), U23 (format rule), U32 (`terms = {}`), U33 (`threads=0` and `a::b::c` raise), I27 (empty `items` -> 200 `[]`; extra field rejected), sketch 4 (`ValueError` naming the key) assert one outcome. The owner confirms them in 99_back. |
| QA 15 | major | Several prices | Done for E. E.8 uses the front section 6 planning rate ($0.641/h, the rate QA proposes; section 6 wins if it differs): 9 min = $0.10, 15 min = $0.16 [estimate]. |
| QA 21 | minor | Same test files created twice; E-G5 count | Done. E5.5: `git mv tests/test_trl_adapter.py tests/rl/` (A7's file) and add 2 cases; E-G5 ">= 7 negative cases (A7's 5, A-G7, plus E's)". |
| data 2 | blocking | Gates resample nodules, not patients | Done for E: E4 TOST, E-G10 and GP5 use the 90 % patient-level (`group_id`) paired bootstrap, BCa for the decision, cluster count printed. Round 2 (QA item 3): the decision interval is now A1's rule (percentile at 50 patients or more; cluster t below), not BCa. |

Own findings while verifying (not in the reviews):
- A7 already creates `train/trl_adapter.py` with `check()`, `TRLContractError`, version checked last and the
  `OBJDEL_TRL_UNPINNED` opt-in [A.md A7]. E's old text defined `check_private_names()` / `TrlDriftError` with the
  version first. E now moves A7's file unchanged (R-5: no renames in moves) and adds only `SOURCE_SHA256` and the
  `_step` / `_buffered_inputs` / `num_iterations` checks in a change commit; all names in E follow A7; sketch 3 and
  I5 unset the opt-in.
- `chain()` fills placeholders with `str.format` [code: policy.py:435], so an unknown `{x}` raises and `{{rd}}` gives
  `{rd}`. U23 said "`{x}` left as is", which would break golden 5. `to_argv` now states the rule, and that JSON built
  from fields (`--vision`, `--grpo`) is never formatted.
- `k3_harness.py` is in `experiments/phase2/` [measured: ls], not `experiments/phase1/` (U38 fixed).
- CLEAN8 step 1 asks for "`RlSpec.max_completion_length` (phase E2 field)"; E2 has no such field. The cap goes through
  `TrainerFlags.grpo={"max_completion_length": 32}` (E2 note); the CLEAN text may name it so.
- A-D add plans (`a3_measure`, `b1_rows`, `b2_fp16` (B's name; was `b2_probe`), `b3_v5`, `c_speed`, `d_gpu_test`) and CLI flags before E. E6.2
  ports the live ones (golden 5 covers them) and E2 turns the new flags into stage fields.

---

## Review log (final review round 2, 2026-11-08)

Source: `docs/implementation2.0/reviews/final_r2_{arch,fresh-arch,cuda,qa,data,rl}.md`, every item that names E. Each
claim was checked in the code at be52b1c (`vast/policy.py`, `train/dataset.py`, `configs/r5_shared.yaml`,
`pyproject.toml`, `.gitignore`, `final_eval.py`, `runs/r5/eval/eval_val/`) or in the current A, B, C, DOC and front
texts before the change. "Done" = the text changed as the item asked; "Done, changed" = done another way, with the
reason; "E part only" = another phase owns the rest. The `data` and `rl` reviews have no E item.

| Review, item | Severity | Item | What I did |
|---|---|---|---|
| arch 1 | blocking | E ports and pins `rl_shifted.yaml`, which B1.5 archives before E | Done. Verified B1.5 (B.md) archives it and adds `b1_rows.yaml`. `b1_rows` replaces `rl_shifted` in goldens 5 and 6(a), E.3.2, E1.5, E6.2, E6.7, E6.8 (`CONFIG ?= specs/b1_rows.yaml`), E-G9, E.7, U19, U20, I10, I13, I15, sketches 4 and 8. One archive folder, `archive/phase2_runs/configs/`; E6.6 merges B's `archive/phase2_configs/` into it if B's text still names that folder (B owns B1.5). |
| arch 2 | blocking | After E8, Vast state resolves to an empty folder | Done. Verified `home()` falls back to `HERE` [policy.py:46,84-85], `spend()` returns 0.0 without a ledger [198-199], `watchdog.sh` read via `HERE` [271], `.gitignore:25-35`. E8.1: `home()` defaults to `repo_root()/experiments/phase2/vast`; state, plists, shims and `watchdog.sh` stay there. Tests (a)-(c) in `tests/infra/test_paths.py` (U48), E-G15, E.7 risk row. The missing-ledger refusal is a behaviour change, so it gets its own change commit (R-1). |
| arch 4 | major | E0's interim contract cannot pass: B1 adds upward edges from `data` | Done. Verified the B1.1/B1.3 imports in B.md (`targets.frame`/`augment`, `reward.cache.rl_eligible`, `eval.ceiling.loo_caches`). Four named ignores with wildcards, each tagged with the move that removes it (E3b, E3c); E0 re-derives the exact edges on the E0 commit; each move commit drops its ignores; negative test (e). Wildcards and the unmatched-ignore error checked on the import-linter docs [web]. |
| arch 5 | major | E5 `dataset.py` table is stale (A8 renames, B1.4 `window_record`, `target_tokens` home) | Done. Verified A8 (A.md: `target_tokens`, `voting`, no aliases), B1.4 `window_record`, C11's use of `target_tokens` (C.md). `target_tokens` -> `policy/rows.py` (sft and rl import it); `voting` -> `verify/store.py` (arch 12); `window_record` -> `data/prepare.py`; name-count sentence rewritten. |
| arch 6 | major | E5 changes A7's contract but A's S07, one A7 negative case, the coupling page and E-G5 are not updated | Done. Verified A7 step 5 (5 negative cases, one removes `shuffle_sequence_dict`) and S07 (29 items). E5.3 (a)-(c): replace that case with the `_prepare_inputs` source case, S07 reads `len(trl_adapter.ALL_ITEMS)`, coupling page count; E5.5 and E-G5 count "A7's 4 + source case + 2 attribute cases". |
| arch 7 | major | No stage can express the soft-fail `"~"` step | Done. Verified `cmd_run`'s `job()` [policy.py:484-486] and `r5_shared.yaml`: steps 1 (bash) and 4 (train_rl) start with `"~"` (the review named only step 4). `StageBase` (kind, name, `soft`) in E2.1; every stage inherits it in E6.1; `to_argv` puts `"~"` first; E6.0 spike, E6.7, U21, U23, I13, E-G9. |
| arch 8 | major | Rows pin has two names (`data_revision` vs `rows_revision` + manifest SHA) | Done. Verified B1.4 step 5 and B1.5. `RunSpec.data_revision` deleted; `RlStage.rows_revision` (40 hex) + `rows_manifest_sha256` (64 hex), both or neither; `to_argv` emits B1.5's flags; `run_spec.json` records them; preflight refuses a manifest-SHA mismatch. E6.5, E6.7, E6.8, U21, I15, E-G9; round-1 row arch 10 marked superseded. |
| arch 9 | major | CLEAN8 deletes `_LORA`, which `rl_lr` reads | E part only. E1.1 and E.8 now state that CLEAN8 must keep `_LORA` (as `frozenset(MODELS)`) and `_GRAD_CKPT` in `policy/models.py`. CLEAN and the front own the rest. |
| arch 10, QA 14 | major | E-G10 box part and E-G13 need a rental after E; none exists | Done, changed. E.8 rewritten around an optional, owner-approved V4 (owner decision 7). Estimate corrected: the reviews' "15 min, $0.16" leaves out launch to sync (17-31 min [measured: front 6.3 row 1]), so V4 is about 26-40 min, $0.28-0.43 [estimate]. E4 is not done until the box part passes; `sample_vllm.py`, `hf_sample.py` and `train/sample.py`'s `sample()` go in a follow-up commit after V4; E-G10, E-G13, E-G14 and E.9.9 say "open until V4". The front owns the 6.3 row and 99_back decision 1. |
| arch 11 | major | E8.1 splits `vast/policy.py` with no name table; preflight/hf CLIs contradict E8.2 | Done, extended. Verified all top-level names [policy.py:46-627]. The review's table left out about 25 names (`CARDS`, `RUN_ID_RE`, `BAD_STATUS`, `MARGIN_S`, `OPS`, `BAD_MACHINES`, `PULL_EXCLUDE`, `ONSTART_MAX`, `INET_CEIL`, `die`, `is_raw`, `governance`, `secret`, `vast`, `sh`, `destroy`, `onstart_cmd`, `offer_cols`, `get_run`, `billed`, `alert`, `FLAG`, `CLI`, ...); the E8.1 table places every one. `bundle.py` is `experiments/phase2/vast/bundle.py` [measured: ls], not `phase2/bundle.py`. Thin-CLI rule for preflight/hf_push/hf_pull; `vast/__main__.py`; `old_names.json` covers the five sources (E-G16, U10); `POLICY_MODULE` note for DOC13. |
| arch 12 | minor | `cache_loader` belongs in `verify`, not `policy` | Done. Verified it uses only `reward` and `geometry` [dataset.py:68-71,189-200,290-307] and its torch-free callers (`eval_run.py:90`, `smoke_grpo.py:124`, tests). `cache_loader`, `_caches`, `voting` move to `verify/store.py` in E3b ("the one doctor-store loader module"); E.3.1 box updated. |
| arch 13 | minor | E3c names `policy.prompt` before E5 creates it | Done. `prepare_rl` (in `experiments/phase1/` [measured: ls]) passes `objdel.train.dataset.PROMPTS` until E5; the E5 move commit switches it. |
| arch 14 | minor | Stale names | E part done: `b2_probe` -> `b2_fp16` (E6.2, own-findings line); E.0 K3 row says "Phase C's K3a fallback (C1)". Front and CLEAN own their parts. |
| arch 16 | minor | DOC4's config glob misses `specs/` | E part only: E8.5 hand-off states that DOC4 must glob `specs/*.yaml` and `docs/reference/configs.md` must document each spec. DOC owns DOC4 and DOC13. |
| arch 17 | minor | Two owners for ADRs 0017-0020 | Done. E8.5 and E.8: E writes each ADR in the commit of the task the front table names; DOC4 checks hold them; DOC9's E row should list them (DOC owns that table). |
| arch 18 | minor | The TOST trainer arm already exists on disk | Done. Verified `runs/r5/eval/eval_val/items.jsonl` and `final_eval.json` (`items` 1084, `K` 8, `seed` 11, `ckpt` = the R5 final model) and the trainer path in `final_eval.py:73-75`. E4, E.8 and GP5: the box runs only the sampler arm. Note: the file's `ckpt` key is `/workspace/runs/r5/final_model`; the review calls it `rl-full-3ep`. |
| arch 20 | minor | `data/quality.py` seed imports `transformers` | E part only: E0 states the B1.3 precondition (no transformers or torch in `data/quality.py`) and stops if `pure-core` fails; E3c row repeats it. B owns B1.3. |
| CUDA 20 | minor | `objdel._k3` cannot be built "with the `train` extra" | Done. Verified the build backend and its `requires` [pyproject.toml:33-36; `requires` at line 35]. E.5: box build against the image torch (`cuda/k3/setup.py`, `CUDAExtension`, `--fmad=false`, `-arch=sm_120`, `pip install --no-build-isolation`, or `cpp_extension.load`); lazy import, C1 self-check decides. |
| QA 3 | blocking | A and B/E use different decision intervals | E part done. Verified A1 step 5 (percentile at >= 50 patients, `signal.cluster_t` below; no BCa in A1). E4, E-G10 and GP5 use the A1 rule; val has 78 patients, so the percentile interval decides. Round-1 row data 2 annotated. |
| QA 8 | major | Unregistered `box` marker fails under `--strict-markers` | Done (first option). Verified A registers `gpu` with `--strict-markers` and the front Tests row names one GPU marker. No `box` marker: GP3-GP5 carry `gpu` (GP5 also `slow`) and skip with "box only: V4 rental" when the build is not CUDA (E.9 conventions). |
| QA 21 | minor | "exam" still used | Done for E: E.3.3 and E.5 now say "test evaluation". 0 hits left. |
| QA 25 | minor | Long sentences | Partly: all new text in this revision uses short sentences and lists. Existing long cells were not rewritten. |
| QA 26 | minor | B2/B3 path rule differs from the front | Done. The "otherwise ... post-E paths" clause is deleted in E.4 ("When E starts") and E.8; B2/B3 land before E on pre-E paths (front rule 1). Round-1 row arch 1 annotated. |
| fresh-arch 6 | major | No recorded decision on the RL trainer framework | E part done: ADR 0019 (E8.5) records the options (TRL + adapter chosen, verl, own loop), the reason and the trigger. The front owns its section 4 row. |

Own findings while verifying (not in the reviews):
- `r5_shared.yaml` has two soft steps (1: the box-info bash step; 4: the shared-prefix arm). So `soft` sits on
  `StageBase`, not only on `RlStage`.
- `cmd_launch` already refuses an allowlist with no `image` [code: policy.py:332-333]. Only the missing-ledger refusal
  is new behaviour.
- `secret()` uses `die()` and `objdel.secrets`, so it goes to `remote.py` with `die`, not to `ledger.py`.
- With E4's archive moved after V4, `train/sample.py` is the only file left in `train/` until that commit. E-G14
  lists it as a named exception.

---

## Review log (final review round 3, 2026-11-08)

Source: `docs/implementation2.0/reviews/final_r3_{arch,cuda,fresh-arch,qa,data,rl}.md`, every item that names E. Each
claim was checked in the code at be52b1c or by a command in this session before the change. "Done" = the text
changed as the item asked. "Done, changed" = done another way, with the reason. "E part only" = another phase owns
the rest. The `data` and `rl` reviews have no E item (the letter "E" in `final_r3_rl.md` is a B3 numerics option).

| Review, item | Severity | Item | What I did |
|---|---|---|---|
| arch 3 | blocking | DOC accepts a red `make test` from the first E move until DOC13 | E part done. R-1 now says: a move commit rewrites every backtick repo path in the DOC4 `doc_files` set with the same map; the E6.6 switch commit adds each `specs/*.yaml` to `docs/reference/configs.md`; the E8.1 move commit sets `POLICY_MODULE = "objdel.infra.vast.ledger"` (E8.1 text, E8.5 hand-off). DOC owns the DOC.7 row and DOC13 step 3; CLEAN owns its archive-commit rule. |
| arch 4 | blocking | E2.5 moves `code_hashes`; the preflight code hash then covers 0 files | Done, changed. Verified `train_sft.py:70,225-233`, `preflight.py:60,69`, `policy.py:377,399`; the `src/objdel/**/*.py` glob from `src/` gives 0 files [measured]. New `core/provenance.py` (E2.5): `repo_root()`, `CODE`, `code_hashes()`. Change: the root search walks up from the module's own file path, not `objdel.__file__`, so `core` imports nothing from objdel. `infra/paths.py` is dropped. The file-count floor is the E0 count (73 at be52b1c [measured]), not the review's ">= 100". AST test widened to all of `src/objdel/` (U34, E-G15, E.7 row). |
| arch 5 | major | E3b rename keeps `verify` above `data` | Done. Exact `layers` list per commit (E0, E3a, E3b, E3c, E5, E2, E8) in an E0 table; negative test "planted `verify -> data` fails from E3b on"; E5 step 1 points at the table. |
| arch 6 | major | `__init__` re-exports have no move rule; 54 `from objdel.eval import` lines | Done. Verified `eval/__init__.py:11-15` and the counts (54 lines, 13 for `metrics` [measured]). R-1: `__init__` files define no logic; callers import from the defining module; only `geom/__init__.py` re-exports `polygonize`. E3a row changed (no merged exports). E3b deletes the `eval/__init__.py` imports; `test_move_complete` counts definitions only; U10 extended. |
| arch 7 | major | DOC5's docstring gate fails the moves into `src/` | Done. Re-ran `ruff check --isolated --select D100,...,D106` on the 12 files: 59 (58 D103 + 1 D101) [measured, ruff 0.17.1]. New R-6; named in E2.4/E2.6 (new numbering), E4.1 and E8.1; E.8 Mac time + about 1 h. |
| arch 8 | major | `to_argv` has no flag order; sketch 8 compares a dict | Done, changed. Re-ran the pydantic probe (2.11.7): mixin fields first [measured]. `FLAG_ORDER` per kind in E6.1; golden 5 is an ordered list of argv lists; sketch 8 and I13 use list equality; new U49. Change: the review said "skip defaults that today's argv omits"; that would drop `--save-steps 50` from `rl_shifted.yaml` (50 is the default [code: train_rl.py:31]). So `to_argv` emits the fields in `model_fields_set`. |
| arch 9 | major | `hf_push`'s licence path points to nothing after E8 | Done. Verified `hf_push.py:60,75`, the file's SHA-256 `cfc7749b...` [measured] and the box call [code: policy.py:522-524]. E8.1 moves `licenses/` to `src/objdel/infra/licenses/`, reads it with `importlib.resources`; U50; E-G15 AST test covers `with_name` and `.parent /`. |
| arch 10 | major | `EvalStage` lacks live `final_eval` flags; E2.1's rule is too narrow | Done, changed. Verified `final_eval.py:48-56` (`--model`, `--seed` 11, `--per-call` 64) and A3's `--rescore`/`--hash-view` (A.md). `EvalStage` gains `model`, `seed`, `per_call`, `rescore`; `hash_view` stays out. E2.1's rule covers `train_rl`, `g_probe`, `final_eval`, `select_ckpt`, `train_sft`; the E6.0 spike also parses `eval`, `select`, `probe`. Change: E6 runs after E4, so `EvalStage` never has `train`; E4 deletes `--train`, and golden 5 records a dropped-flag map for it. |
| arch 12 | minor | `data/views.py` has no task; "thin CLIs only" has no plan for the builders | E part done. Verified the 907 lines [measured: wc]. `views` removed from the E.3.1 box and the E.3.3 row; E8.2 lists the thin CLIs and keeps the data builders in place. The front owns 2.2 and 2.3. |
| arch 13 | minor | `truncated` in `policy/sampler.py` makes the reward import the sampler | Done. E4.1: `truncated` goes to `policy/models.py`; sampler and reward import it from there. |
| arch 14 | minor | `rl/__init__.py` is created after the move that needs it | Done. E5 step 1: the move commit creates `policy`, `rl` and `sft` `__init__.py` files, each with a one-line docstring and no imports (U40). |
| arch 15 | minor | R-1 names no rewrite tool | Done. R-1: one libcst codemod, `tests/goldens/phase_e/rewrite.py`, input = the move table JSON; libcst 1.9.0 is the latest [measured: PyPI JSON]; `rope` is the fallback. |
| arch 16 | minor | V4 has two prices | Done. E.8 and owner decision 7: about 39 min = $0.42 (front 6.3). 99_back is not E's file. |
| CUDA 8 | minor | K3 production flags differ from the tested build | E part done. Verified `kernels/K3/test_torch_gpu.py:12,41-45`. E.5 copies the recipe: `TORCH_CUDA_ARCH_LIST=12.0`, no `-arch`, the three `-U` macros, host `-ffp-contract=off`. C owns C1 step 7. |
| fresh-arch 7 | minor | Long-lived interfaces have no 3D note | Done. Verified `checker/volume.py` (`betti3d`, `scipy.ndimage`) and D.3.3 ("deferred to the voxel phase"). ADR 0020 (E8.5) gets a "3D extension" paragraph, (a)-(d), no code. The front owns its section 4 link. |
| QA 9 | major | `pytest -m gpu` collects every phase's box tests | E part done. E.9 conventions: GP1-GP5 name V4 and one file list; E never runs `pytest -m gpu` without a file list; tests that need a run artifact skip with a reason when it is absent. D owns D-G3 and D8. |
| QA 12 | minor | Three different GPU skip rules | E part done. E.9 points at A.9's root `tests/conftest.py` hook (`objdel._gpucheck.build() == "cuda"`, no torch). GP1 loses its own skip line. If the hook is missing on the E0 commit, E0 adds exactly that hook. A owns A.9. |
| QA 16 | minor | Same rental items have different prices | E part done (see arch 16). |
| QA 21 | minor | Long sentences remain | Partly. All new text in this revision uses short sentences, one rule each, with a "Why:" where needed. Older cells were not rewritten. |

Own findings while verifying (not in the reviews):
- The code-hash floor ">= 100" in arch item 4 is wrong for be52b1c: the `CODE` globs match 73 files [measured]. The
  test uses the E0 count.
- `preflight.py` imports `train_sft as TS` and calls `TS.code_hashes()` [code: preflight.py:69]. The E2 move commit
  must rewrite that call, or preflight breaks before E8 (E2.5 now says so).
- E2's step numbers shifted (provenance is the new step 5), so E.9 task ids now read E2.1-E2.7.

## Review log (consistency pass, round 4 prep)

| Item | What changed |
|---|---|
| One source | E.8 and owner decision 7 point to front 6.3 ("V4") for the V4 price (now $0.54 with the bandwidth estimate) instead of quoting $0.42 |

## Review log (final round 4)

Cross-file consistency fixer, round 4. Sources: `reviews/final_r4_{cuda,rl,qa,arch,fresh-arch}.md` and their verdicts.

| Item | Severity | Change | Evidence |
|---|---|---|---|
| arch-7 | major | E6 step 8: the E6.6 commit deletes `projection`, `HF_RATE`, `ceiling`, `--gpu` and the `CHECKS` row from `preflight.py`; test: every `CHECKS` function's lazy imports resolve and a fixture preflight has no "check crashed" | `preflight.py:216-217, 250, 269-276`; `policy.py:313-316` |
| arch-1 | blocking | One rule for emptied packages, so E8.4's `exhaustive_ignores = ["_gpucheck"]` holds. R-1: the move commit that empties an old package deletes its `__init__.py` and directory (`checker/`, `geometry/` E3a; `reward/` E3b; `targets/` E3c; `eval/`, `train/` E5). `train/sample.py` now moves whole in the E5 move commit (arch-4), so `train/` is gone after E5, as the E5 table header says. `eval/audit.py` (CLEAN6 keeps `hack_audit` and its helpers) gets a home: `data/audit.py` in E3c (E3b row, E3c table, box diagram). E8.4 states the precondition and the `_k3` case (E.5); U1 lists the children of `src/objdel/` and names an extra one. TID251 comment, E0 and E5.4 scan text, E-G14 and the E5 gate now say E5, not E4 | `sample.py:46-137`; `eval/audit.py:54-155`; CLEAN.md CLEAN6 steps 1-2; `ls src/objdel/{checker,geometry,reward,targets,eval,train}`; import-linter layers docs (exhaustive: an unlisted child fails) |
| arch-2 | major | E2 step 6: name tables for `select_ckpt.py`, `train_sft.py`, `g_probe.py`, `train_rl.py`, `final_eval.py`. `view_rows`, `_score`, `score_items`, `SUBSET_SEED` -> `policy/evaluate.py` in E2 (g_probe and final_eval import them; E4 now extends the module); checkpoint `stop_ids` -> `policy/sampler.py`; the bridge, `hf_check`, `score_view`, `prior_ci`, `score_probe`, `vision_check`, `exited` -> `sft/select.py`; `main`, `judge`, `topology_only`, `vision_config`, `dump` stay in the CLI, so every `prereg` and `baselines` call stays in `experiments/` (U5, E-G8); `prereg.py` and `baselines.py` do not move. `train_sft`: `CODE`, `code_hashes`, `sha256` -> `core/provenance.py`, the rest -> `sft/train.py`, `main` stays; both `ROOT`s -> `repo_root()`. U10 covers the five scripts; E2 step 7 and the coverage map list the new homes; E4's test names `sft.select.sampler_cmd` | `select_ckpt.py:50-63,179-220,271,354,397-403`; `train_sft.py:68-263,299`; `g_probe.py:34-35,175-194`; `final_eval.py:23-46,63`; `wc -l` (prereg 160, baselines 213) |
| arch-3 | major | `TraceStep` -> `policy/callbacks.py` (not `rl/callbacks.py`); `stop_ids`, `truncated`, `_text` -> `policy/models.py` in the E5 move commit (not `rl/reward.py`); `rl/reward.py` keeps `TRLReward`, `LOGGED`, `GATED`, `AGGREGATOR`. E4 step 1 no longer moves `truncated` (it is there from E5). E5 gate: no `policy.*` or `sft.*` import of `objdel.rl`; I1 adds a `policy -> rl` case. Box diagram updated | `train_sft.py:390`; `g_probe.py:91,93`; `select_ckpt.py:187-200`; `callbacks.py:210-237` (TraceStep needs only torch and transformers); `rewards.py:68-86` |
| arch-4 | major | Name table for `train/sample.py` in the E5 move commit: `sample` -> `rl/trl_adapter.py` (calls `trainer._generate`, so E-G5 allows only the adapter; name kept, R-5); `G`, `export`, `write`, `_path`, `export_prompts`, `import_completions` -> `policy/sampler.py`. The V4 follow-up commit no longer deletes `sample()`: it stays as the live G-probe draw and the follow-up archives only `sample_vllm.py` and `hf_sample.py`. E4 step 1 extends `policy/sampler.py` and reuses `sample.write` (same JSONL job as `sample_vllm.write`, a name clash otherwise). U10 covers `sample` | `sample.py:46-158`; `g_probe.py:44-48`; `final_eval.py:75`; `select_ckpt.py:212-220`; `sample_vllm.py:68-69`; tests `test_train_trl.py:413-471`, `test_p2_bridge.py` |

## Review log (final round 5)

Cross-file consistency fixer, round 5. Sources: `reviews/final_r5_{cuda,rl,data,arch,qa,fresh-arch}.md` and their
`*_verdict.md` files. Only the CONFIRMED items assigned to this fixer are listed.

| Item | Severity | Change | Evidence |
|---|---|---|---|
| arch-1 | blocking | E4 step 3: the owner tags `pre-e4-eval` before the switch; the trainer path lives in git history, also before V4 (the TOST trainer arm is on disk) | `final_eval.py:73-75`; E4 gate |
| cuda-4, qa-3 | major / blocking | E.8 V4: plan file `e_v4.yaml` with step 0 = `b1_rows.yaml` step 0 (CUDA build), then `env_check.sh`, the box test list, the sampler arm (pool 1); new U51 checks it. E.9 "One skip rule": on a box the hook fails on a CPU build (A.9) | `policy.py:399, 432-433`; `CMakeLists.txt:15` |
| qa-6 | major | GP1: CUDA build against the Mac golden `tests/goldens/e_gp1.json` (`experiments/phase2/perf/make_e_goldens.py`); "scored with `checker_build` cuda and cpu" removed; missing golden fails | `gpucheck_py.cpp:113`; `CMakeLists.txt:15` |
| fresh-arch-1 | blocking | E.8 V4 and owner decision 7: V4 needs a budget raise (not "inside the $3") | front 6.3 |
| arch-2 | major | E6 step 1: `RunSpec.bundle_view` (the SFT training view; replaces the RunArgs launch configs' `view`) with two spec rules (rl stage without sft stage needs it; not beside an sft stage). E6 step 8: Mac data root (`OBJDEL_DATA_ROOT` or `data`), one `preflight.local()` for `{data}` (B1.3's `{data}/rows/<x>` -> `views/<x>` branch kept; `{rd}` -> not a data input; other forms raise), a check -> input -> runs-when table (six RunArgs-era checks on T; new "spec inputs" check; B1.3 "ROWS root" per rl rows root), the data hash, deletion of `preflight.load` and the `rows_plan` reader; ported specs set `bundle_view: views/lidc-bundle-v3`. Tests in step 8, I15, U21. Not taken: the review's "run `patients`, `coverage`, `ct_window`, `targets`, `prompts` on the rl rows root", because those checks read `rows.jsonl` keys (`epoch`, `window`, `target`, `side`) and the B1 rows root has `lidc/train.jsonl` and no `view_hash`; its gates are `quality.failures` (B1.3), so the rows root gets the ROWS root check only, as before E | `preflight.py:49-70, 101-201, 244-251`; `train_sft.py:82-84`; `policy.py:433`; `r1_qwen3vl2b.yaml:4`; `rl_full.yaml:14-16`; `r5_shared.yaml:9-11`; front 6.2 "Why `view` stays the SFT bundle"; B1.3 changes 1-4 |
| arch-3 | major | E6 step 0: the spike round-trips the `r1_qwen3vl2b` train stage through `CliApp.run(SftStage)`. Step 1: `SftStage` = RunArgs fields minus `stages`, `sampler`, `rows_plan`; `vision` before-validator replaces `__post_init__`; "SFT argv" rule for `to_argv(SftStage)` (kebab flags in `FLAG_ORDER["sft"]`, lists and dicts as JSON, no temp YAML). Step 4: `main` -> `CliApp.run`; `resolve` and the `branches_from` rebind use `model_copy(update=)`; manifest uses `model_dump`; the preflight dry-run test passes `to_argv(stage)[2:]` instead of `--config`. Step 7 and I14: round trip; new frozen-stage test. U22 subtracts `{kind, name, soft}` and the three runner-only fields. E2.6 table row notes the E6.4 change | `train_sft.py:72, 77-111, 150-156, 257-260, 276, 297, 301`; `policy.py:441`; `tests/test_p2_scripts.py:46, 63, 354`; 6 test files name `RunArgs` (grep -l); pydantic-settings docs (CLI lists and dicts accept JSON) |
| arch-4 | major | New E6 step 8b: `push_target(spec, run_id)` (last rl or sft stage; rl -> `out` basename and `hf_push --rl`; sft -> `output_dir`; none -> die before token or ssh); `vision_refs(spec)` over `SftStage.vision` and `SelectStage.vision`; `--cap-hours` stays on every launch (front section 6) and dies above `compute.cap_hours` or on an offer card other than `compute.gpu`; `b1_launch.yaml` archived in the E6.6 commit (the `v1-b1` tag keeps it); `e_v4.yaml` is a `RunSpec` at the path front 6.3 names (E.8 V4 line, U51 updated). E8.1 table lists `push_target`. Tests in step 8b and I15 | `policy.py:303, 321, 437-440, 458-461, 471-476, 517`; `train_sft.py:153-156`; front 6.2 launch lines; front section 6 cap rule; B1.5 `hf_push --rl` |

## Review log (final round 6)

Cross-file consistency fixer, round 6. Sources: `reviews/final_r6_{cuda,rl,data,arch,qa,fresh-arch}.md` and their
`*_verdict.md` files. Only the CONFIRMED items assigned to this fixer are listed.

| Item | Severity | Change | Evidence |
|---|---|---|---|
| arch-4, qa-3 | major | E4 gate: checkpoint = HF `rl-full-3ep` (the `ckpt` key is a dead box path). E.8 V4 step 1 = `env_check.sh {rd}/ckpt {rd}/env_check.json`; step 3 = `EvalStage(ckpt="{rd}/ckpt", ..., out="{rd}/eval_v4_sampler", k=8, seed=11)`; no push line. U51, GP3 and GP4 name the paths | `env_check.sh:9`; `runs/r5/eval/eval_val/final_eval.json` ckpt; `final_eval_2026-11-08.md:3` |
| qa-1 | blocking | GP5 becomes the Mac test `test_e_g10_tost` (`data` marker, `OBJDEL_V4_RUN`), run after the V4 pull; E-G10 row, E4 gate, E.9 conventions and the closing paragraph follow | `.gitignore:23`; `policy.py:376-384, 490` |
| data-2 | major | E3b `verify/store.py` row, U33 and `RlStage.train_views` / FLAG_ORDER state the form with and without B1.8 | B1.8 |
| fresh-arch-1 | blocking | "When E starts", R-4 and the E.7 merge-conflict row: E0 waits for the last `main` commit of the front (b') list; the `after-raise` branch merged after E0 uses the move tables | front section 1 |
| arch-1 | blocking | E6.1 `EvalStage`: `view`, `split`, `out` become tuples (`view` min 1; `split` empty or one per view, empty = "val"; `out` one or one per view; no comma in an item), read by index as one group. New "Eval argv" rule: `FLAG_ORDER["eval"]` = `ckpt, views, model, k, seed, per_call, rescore`; at `views` `to_argv` writes `--view v[i] [--split s[i]] [--out o[i]]` per index, or one `--out` after the last group. Found while probing: a one-letter field is `-k`, not `--k` (pydantic-settings; `--k 8` exits 2), so `to_argv` writes `-k` and E6.4 updates A3's I17 argv; a comma splits a tuple item, so it raises. E4 step 3 keeps A3's view loop. E.8 V4 step 3 uses the tuple form and its argv. Tests: E6.7 (b1_rows eval argv = B1.5 step 4 minus `--train`, plus round trip), U21, U49, U51, I14 | B.md B1.5 step 4; A.md A3 (`action="append"`, `--out` per view); `final_eval.py:47-56` (one `--view` today, `--split` default "val"); pydantic-settings 2.16.2 probes (interleaved argv -> 2-item tuples; `-k`; `a,b` -> two items) |
| arch-2 | major | E2.1 `StageBase`: `kind: CliSuppress[str]`, `name: CliSuppress[str] = ""`, `soft: CliSuppress[bool] = False`; every stage sets `kind: CliSuppress[Literal[...]] = "<kind>"` (E2.1 `RlStage`, E6.1 all others). `RunSpec._rules`: every stage name non-empty. E6 intro: one parser form `cli_parse` = `CliApp.run(Stage, cli_args=argv, cli_settings_source=CliSettingsSource(Stage, cli_kebab_case=True, cli_implicit_flags=True))` (the keyword form raises `extra_forbidden`), and one equality rule `same_stage` = dumps equal without `name` and `soft`. Spike, E6.4 (`main`, `--spec`/`--stage`), E6.7, E2.6 row, I14 (`same_stage(S1, S2)`) use them; new E6.7 test: every non-cmd ported stage parses back. U21 "missing kind" still raises (discriminated union: `union_tag_not_found`) | final_r6_arch probes 1-3 re-run (exit 2 "required: --kind, --name"; keyword form `extra_forbidden`; CliSettingsSource form `True`); own probe: `CliSuppress` hides help only, `--kind cmd` on `EvalStage` fails validation, a spec stage without `kind` raises `union_tag_not_found` |

## Review log (final round 7)

Cross-file consistency fixer, round 7. Sources: `reviews/final_r7_{rl,arch,qa,fresh-arch}.md` and their
`*_verdict.md` files. Only the CONFIRMED items assigned to this fixer are listed.

| Item | Severity | Change | Evidence |
|---|---|---|---|
| arch-2 | major | E6.1 `SelectStage` gains `model`, `k`, `seed`, `procs`, `timeout_s`, `out`, `no_vision` (validator: exactly one of `vision`/`no_vision` unless `only`); `ProbeStage` gains `procs`. New "Select and probe argv" rule with `FLAG_ORDER["select"]`/`["probe"]`; E6.4 names the flag spellings that change (`--timeout-s`, `-k`); the preflight dry run keeps `--procs 1` until CLEAN16. New U52 (argparse dest names equal stage fields) | `select_ckpt.py:362-379, 383-384`; `g_probe.py:139`; `tests/test_p2_scripts.py:369, 384`; `policy.py:440-446`; grep: no plan or test writes `--k` or `--timeout_s` |
| arch-4 | major | E6.1 `EvalStage.ckpt` optional with the validator "exactly one of `ckpt` and `rescore`"; E6.4 "final_eval after E6": `--hash-view` is read by a two-flag argparse before `cli_parse`. New U53 | A.md A3 (`--hash-view`, `--rescore`), I08 argv; `final_eval.py:48-56` |
| fresh-arch-4 | major | B1.5 now deletes `VIEW_SEP`, `_tagged` and `--train-views` before E0: the E3b `verify/store.py` row, `RlStage` (no `train_views`), `FLAG_ORDER["rl"]` and U33 state the one form | B.md B1.5 |
| arch-3 | major | E.8 V4 step 2 is now `CmdStage(name="box_tests", argv=("-m", "pytest", "-q", "-p", "no:cacheprovider", "-m", "gpu", <the three E.9 files>))` with the reason (a leading `pytest` runs as `python -I pytest`, fails, and the `&&` chain stops before the sampler arm); E.9 "Each box test names ..." writes the same argv; E6 step 8b points to E.8's argv; U51 adds "every `CmdStage.argv` of `e_v4.yaml` starts with `bash` or `-m`; step 2 starts with `-m, pytest`" (as B I22a for V1) | `policy.py:486` (`def job`: non-`bash` argv -> `/opt/venv/bin/python -I ...`), `policy.py:490` (`" && ".join`), `policy.py:451-452` and `preflight.py:205` (`-m pytest -q -p no:cacheprovider` form) |

## Review log (final round 8)

Cross-file consistency fixer, round 8. Sources: `reviews/final_r8_{cuda,data,arch,qa}.md` and their `*_verdict.md`
files. Only the CONFIRMED items assigned to this fixer are listed. Each fact was checked in the repo first.

| Item | Severity | Change | Evidence |
|---|---|---|---|
| arch-1 | blocking | V4 plan file moves to `specs/e_v4.yaml` (E6 step 8b, E.8, U51); I29: `chain()` outside `specs/` raises before any `remote` call, and `chain(specs/e_v4.yaml)` gives the U51 argv on the Mac; E-G9 names it | `policy.py:469-470` (`chain()` runs after launch, sync, fetch) |
| arch-3 | major | E3b pool row, E.3.3 row, box diagram, E.5 K2 line, E.8 line, U33: the pool is archived by E3b (R-2), never moved to `verify/pool.py`; E6.1 drops `SelectStage.procs`, `ProbeStage.procs` and `procs` in `FLAG_ORDER`; U52 now expects `--procs` to exit 2 | grep: only `final_eval.py:69`, `g_probe.py:157`, `select_ckpt.py:390`, `train_rl.py:50` call `fork_pool` |
| arch-2 | major | U52 no longer asserts plain equality: the golden file holds `dests`, `dropped`, `added` per CLI; (a) `dests - dropped == fields - added`, (b) `dropped` is a subset of `dests` and disjoint from the fields, (c) `added == TrainerFlags fields - dests` for `train_rl`/`g_probe`, empty elsewhere. `dropped`: `{train, hash_view}` (`final_eval`), `{stages, sampler, rows_plan}` (`train_sft`). `added`: `{compile, plain}` (`train_rl`), `{lr}` (`g_probe`). New E2.1 paragraph "Mixin flags that a CLI gains"; E6 step 4 names the flags gained (`train_rl --compile/--plain`, `g_probe --lr`) and lost (`train_sft --stages/--sampler/--rows_plan`, `final_eval --train`) at E6.4; `FLAG_ORDER["rl"]` ends `compile, plain`, `FLAG_ORDER["probe"]` ends `lr`. Not taken: dropping the mixin fields, because `build_trainer` reads all six for both callers | `train_rl.py:23-36` (has `--lr`, no `--compile`/`--plain`/`--trace`); `g_probe.py:131-145` (no `--lr`); `train_sft.py:100-101, 260` (`stages`, `sampler` are `RunArgs` fields, `TrlParser([RunArgs])`); A.md A6 step 3 adds `--trace` to `train_rl`; grep: no repo test or plan passes `--stages`/`--sampler`/`--rows_plan` to `train_sft` (every `--sampler` hit is `select_ckpt`) |

## Review log (final round 9)

Cross-file consistency fixer, round 9. Sources: `reviews/final_r9_{cuda,rl,data,arch,qa,fresh-arch}.md` and their
`*_verdict.md` files. Only the CONFIRMED items assigned to this fixer are listed. Each fact was checked in the repo first.

| Item | Severity | Change | Evidence |
|---|---|---|---|
| cuda-2 | major | E3a: `experiments/phase2/perf/bench_checker.py` is the one `_gpucheck` importer outside `src/` that is not rewritten; own TID251 per-file-ignore in the E0 toml (repeats the `experiments/**` codes); U4 allowlist and U41 allowed files name it | `ls src/objdel` (no `geom`, 2026-11-09); Makefile:33 (`ruff check src tests experiments`); D-G7 (`pre-d2` old venv) |
| arch-1 | major | R-4 "After E6" (the one port rule): the `after-raise` merge after E6 ports `c_speed`/`b2_fp16`/`b3_v5` to `specs/`, folds `c_launch`/`b3_launch` into them and archives them, adds post-E0 flags as stage fields with `FLAG_ORDER` and golden-5 rows, re-points I22b/U-CF4, adds the specs to I29; E.8 dependency text drops "no post-E exception" and points to R-4 | `policy.py:315` (`preflight.refusal` before launch); E6.4, E6.8, E6.8b, I29 |
| arch-2 | major | E2.1 new rule "Session-spend flags after E6.4": B1.5's `--budget-usd`, `--session-start`, `--rate-usd-per-h`, `--fixed-usd` become `RlStage` fields with plain `None` defaults (no env `default_factory`, no env-reading validator, so `to_argv`, golden 5, `same_stage` and I16 do not change); the env defaults and the "budget needs a session start and a rate" check move to one box-runtime function `rl.run.session_budget(stage, environ)`, which `run_rl` calls first, before any model load and `build_trainer`; not a model validator, so the Mac `load_spec(specs/b1_rows.yaml)` and E-G9 pass; `train_rl.main` maps `BudgetFlagError` to exit 2, so B's U45 keeps its exit code. E2.3 names the first call; E6.4 "Env defaults that leave the parser"; E5 move table: `BudgetStop` and `LocGuard` go with `TraceStep` to `policy/callbacks.py`. New U54; coverage map E2.1, E2.3, E6.4, E-G9. Not taken: the reviewer's env `default_factory` on the stage, because the stage would then differ between the Mac spec and the box CLI (`same_stage`, I16, `run_spec.json`) | B.md B1.5 run guard and V1 step 3 (`b1_rows.yaml`: `--budget-usd 2.05`, no `--session-start`), B U45; `reviews/r9arch_probe/env_probe.py` (`session_start=None` with `OBJDEL_SESSION_START` set; pydantic-settings 2.16.2 and 2.10.1); grep of E.md before the fix: 0 hits for `BudgetStop`, `OBJDEL_SESSION_START` |
| arch-3 | major | E8 step 1 new bullet "Coverage gate": E0 records in `runs/E/baseline.txt` the combined coverage of `objdel` plus the files E8.1 moves; if it is below the E0 `--cov-fail-under` value, one change commit before the move adds fake-only unit tests in `tests/infra/` (spend guard `active`/`accrued`/`spend` vs `HARD_STOP`/`ALARM`, `cmd_launch` refusal paths before `vast create`, each preflight `CHECKS` function pass/fail); lowering the floor, coverage `omit` and `pragma: no cover` on moved code are not allowed. E-G1 adds `make cov` on every E commit, the E8.1 move commit included. New U55; coverage map E8.1, E-G1 | `.github/workflows/tests.yml:44` and `Makefile:25-26` (`--cov=objdel --cov-fail-under=90`); `reviews/r9arch_cov.txt` (combined 4,608 stmts / 522 missed = 88.7 %; policy.py 50 %, preflight.py 49 %); this fixer's run of the four infra test files with `--cov-report=term-missing` (133 passed): missed policy.py:185-200, 309-367, 387-419, 468-496 and preflight.py:101-198; `policy.py:309-333` (launch refusals) |

## Review log (final round 10)

Cross-file consistency fixer, round 10. Sources: `reviews/final_r10_{cuda,data,qa,fresh-arch,arch}.md` and their
`*_verdict.md` files. Only the CONFIRMED items assigned to this fixer are listed. Each fact was checked in the repo first.

| Item | Severity | Change | Evidence |
|---|---|---|---|
| cuda-1 | blocking | E6 step 2: `a3_test` and `r_d1` join the A-D plan list; rule "a plan that a front 6.3 rental not yet run names is live"; R-4 (i) ports `c_speed_b`, `c_speed_c6` too, (iii) names `--expect-step1`; I29 gains `RENTAL_PLANS` (each front 6.3 "After E6" spec loads, its `configs/` path raises) | `policy.py:469-471` |
| cuda-3 | major | E8 step 1 "Other callers": the E8.1 move rewrites `run_box.sh`'s dry-run `code_files()` line and re-points `tests/test_box_inputs.py`; U8 also scans tracked `*.sh`, `Makefile`, `*.yaml` (outside `archive/`) for `experiments.phase2.vast.policy`/`.bundle` | `git grep`: today only `vast/vast.sh:6` and the comment in `rl_shifted.yaml:14` hold such strings |
| arch-2 | major | E8 step 1 "Coverage gate" widened: the gate binds first on the E2 move commit, not E8.1. E0's coverage record adds `--cov` for `train_sft`, `select_ckpt`, `g_probe`, `train_rl`, `final_eval`, `sample_vllm`, `hf_sample` beside the four E8.1 files and writes two figures in `runs/E/baseline.txt`: (a) `objdel` plus the names the E2/E4 tables move (per-function lines), (b) (a) plus the E8.1 files. If (a) is below the E0 floor, one change commit before the E2 move (after E2's R-6 docstring commit) adds fake/tiny-model tests in `tests/sft/` and `tests/rl/` (`train_sft.train` with a stub `Trainer` and a 2-row view, which reaches the nested `Rolling.on_save`/`on_step_end`; `vision_reference`, `replace_dir`, `latest`, `write_json`, `preload`, `make_trainer`, `select_ckpt.score_view`/`prior_ci`/`score_items`); same rule before the E4 move if the sampler bodies drop it again. The no-`omit`/no-`pragma`/no-lower-floor rule covers all three moves. E0 bullet names the record; E2 and E4 gates add `make cov`; E-G1 names the E2, E4 and E8.1 move commits; U55 renamed `moved_code_coverage_floor`, scope "E2 move commit and every later E commit", grep over all of `src/objdel`; coverage map adds U55 to E2.6, E4.1, E4.2 | `Makefile:25-26` and `.github/workflows/tests.yml:44` (`--cov=objdel --cov-fail-under=90`); `reviews/r10arch_cov.txt` (TOTAL 4,748 / 516; experiments rows: train_sft 293/166, select_ckpt 269/35, g_probe 121/40, train_rl 45/24, final_eval 57/25, so `objdel` 3,866/220); `reviews/r10arch_probe/cov.json` re-read by this fixer (`train` 97/74, `make_trainer` 20/20, `train_rl.main` 38/23, `vision_reference` 15/15, `score_view` 6/6, `prior_ci` 8/8, `score_items` 4/4); `e2_move_estimate.txt` (423 moved / 232 missed, 89.46-89.61 %); `train_sft.py:353-363` (`Rolling` nested in `train`); `grep -rn "pragma: no cover" src/objdel` = 0; no `omit` in `pyproject.toml`; E.md task order line (E5 -> E2 -> E4 -> E6 -> E8) |

## Review log (final round 11)

| Item | Severity | Change | Evidence |
|---|---|---|---|
| fresh-arch-1 | major | E4 Gate: a B run follows E only on the post-E raise path, where Gate B3 pairs within V3 (B3.9 step 3); V4 still closes E-G10. E.8: V4 checks ride on V3 only if named, and the main agent prices them | front "End state under (b')"; E4 step 3 |
| fresh-arch-3 | major | Polygonizer seam paragraph and the D-on-pre-E line name B1.4a first | B.md B1.4a |

## Review log (final round 12)

| Item | Severity | Change | Evidence |
|---|---|---|---|
| qa-1 | major | Seam paragraph and the move-table note: B1.4a moved `polygonize.py`, D3 `bitquad.py` | B.md B1.4a |
| qa-2 | major | B3.9 step 3 -> step 2 (the B1 re-eval) in the post-E pairing lines | B.md B3.9 |

## Review log (final round 14)

| Item | Severity | Change | Evidence |
|---|---|---|---|
| fresh-arch-1 | blocking | E4 step 1: "Archive in the switch commit" (E4 step 3 commit; tests that import the old samplers move or are re-pointed; rollback `git checkout pre-e4-eval -- ...`; V4 still closes E-G10's box half). Carry-over list adds `--logprobs`. E4 step 3: same commit archives and rewrites the `parity.py` / `r0_parity.py` producer lines. E4 Gate, E-G10, E-G14 (no exception), U2, E.8, the dependency line and owner decision 7 follow | `select_ckpt.py:93-96` (only live caller); `parity.py:5,71,89`; `r0_parity.py:7,44`; `sample_vllm.py:47-49,62-63,84,97`; grep of `tests/` (6 files import or name the old samplers) |

## Review log (final round 15)

| Item | Severity | Change | Evidence |
|---|---|---|---|
| arch-1 | major | E4 step 1: the restore command is removed. New V4-failure rule: fix forward in `objdel.policy.sampler` (deterministic part: change commit + failing row in `tests/policy/test_sampler.py` + re-run on the next owner-approved rental; TOST alone: "not shown", owner decides); label "sampler not box-verified (E-G10 open)" until V4 passes; tag `pre-e4-eval` is a read-only reference. E.7 row: mitigation says archived in the E4 step 3 commit; rollback is fix forward, not "revert E4 only". Owner decision 7 note says the same. Now agrees with E-G14 and U2 (no V4 exception) | `select_ckpt.py:93-96` (the old samplers' only caller, switched in E4 step 1); E-G14; U2 |


# Phase DOC

## DOC.R. Review log

Source: `docs/implementation2.0/reviews/DOC_phase_review.md` (30 items). Each item was checked against the code before the change.

| # | Sev. | What I did |
|---|---|---|
| 1 | blocking | Fixed, simpler than proposed. Path-based `FROZEN_GLOBS` (DOC4). Instead of a per-file lock map, `allowed_versions` = union of pins in `pyproject.toml` and all locks, plus a unit guard. Probe on today's tree: exactly the 13 real stale hits, 0 false (the review's map would also have flagged ADR-0006:23 and `vLLM 0.09 s` in `tests/test_p2_serve.py:288`). Negative tests added. |
| 2 | blocking | Fixed: archive READMEs dropped from check 2 (kept in 11); negative "not reported" test. Verified the 3 retired paths in `archive/cpu_checker/README.md`. |
| 3 | major | Fixed: DOC1 moves `docs/research/research2.0/` to `docs/plan/`, renames the file, updates `docs/README.md:44-46` and `Research2.0.md:10`; `docs/plan/**` in scope; nav "Plan"; DOC-G1 updated. Added a note that other phase pages use the old path. |
| 4 | blocking | Fixed with a change: verified `runs/r5` and `runs/rl` have no `select/` or `manifest.json`, and the model is revision `rl-full-3ep`. Instead of building a new card from a run, `--card-only` downloads the live `README.md`, replaces only Evaluation and Limitations (`card_update`), uploads one file from the Mac. Full-push `card()` unchanged. |
| 5 | major | Fixed (a) with a change: the review's regex matches none of today's cost lines (`gpu-run.md:115-116` are table rows), so check 7 also reads Alarm / Hard stop table rows and must find at least one of each. Fixed (b): DOC10 rule allows action-item boxes; DOC2 only ticks ADR-0009 item 3. |
| 6 | major | Fixed: no cost column in ROADMAP; $3 budget line in ROADMAP and `gpu-run.md`; `HARD_STOP` vs $3 listed as an owner decision. |
| 7 | major | Fixed: check 13 `test_nav_complete`; DOC3 step 3 corrected. |
| 8 | major | Fixed: `validation.links.anchors: warn`. Re-checked with `-c` as asked: clean "No issues found" exit 0; seeded bad anchor "1 issue found" exit 1. |
| 9 | major | Fixed: schema adds `view`, `ckpt`, `seed`, `ci_unit`, `model_ref`. Changed `"image"` to `"item"`, the word `final_eval.py:10,34` uses. Refuses input without `view`. |
| 10 | major | Fixed: DOC1 step 1 scan (ran it: prints nothing over 34 files), RIGA owner OK, owner reads `git diff --cached --stat`; new gate DOC-G13. |
| 11 | major | Fixed with a correction: the table already exists (`overview.md:159-174`), so DOC9 adds Target/Phase columns and the polygonizer, `pixel_values` and rope rows. |
| 12 | major | Fixed: `locks/docs.txt` via `uv pip compile`; CI installs from it. Did not make check 7 read zensical (not one of its 4 packages). |
| 13 | major | Fixed: URLs to `blob/v1.17.1/`; check 7 matches `/blob/vX.Y.Z/` for trl and transformers. |
| 14 | major | Fixed: risk row; order rule; PR template box; `runs.md` records commit and code hash. |
| 15 | minor | Agreed: DOC5 after CLEAN's removal list, gate in the same PR. Reconciled with item 14 (also before B1 preflight or after teardown). |
| 16 | minor | Agreed: DOC5 snippet is now a diff that keeps `ignore` and isort. |
| 17 | minor | Agreed: "reaches $10" also replaced; DOC-G5 = B9 = `\$\d` regex. |
| 18 | minor | Agreed, simpler: docs CI job installs only the lock and builds; `test_docs.py` runs in the existing `python` job. |
| 19 | minor | Agreed: weekly non-blocking Material build in `links.yml`. |
| 20 | minor | Agreed: `actions/cache` restore/save; `.lycheecache` and `/.cache/` git-ignored. |
| 21 | minor | Agreed: check 1 limited to the files outside `docs/`; Zensical covers `docs/`. Kept the stdlib check rather than lychee `--offline` (one tool fewer). |
| 22 | minor | Agreed: section 4 cited; CAD row "Parked". |
| 23 | minor | Agreed: `secrets` in the list; page checked for names only. |
| 24 | minor | Agreed: full precision in JSON; round only in `headline_strings`. |
| 25 | minor | Agreed: `test` and `val` blocks; val strings match README.md:121 ("F 0.806", not "boundary F"). |
| 26 | minor | Agreed: `docs/results/runs.md` everywhere. |
| 27 | minor | Agreed: "32 in scope today, 36 after DOC1 (98 under docs/)" with commands. |
| 28 | minor | Agreed: `.PHONY` and `headline` recipe with empty-`RUN` failure. |
| 29 | minor | Agreed and measured now: 1,484 collected; DOC1 records it again before the first edit. |
| 30 | minor | Agreed: line in `ci.md`; reproduced the traceback in SP. |

---

## Review log (final round 1)

Sources: `docs/implementation2.0/reviews/final_r1_{arch,cuda,data,fresh-arch,qa,rl}.md`, every item addressed to DOC
(`final_r1_rl.md` has none). Each fact was checked against the repo or a command on 2026-11-08 before the change.

| Review item | Sev | What changed | Check made |
|---|---|---|---|
| qa 4, arch 13 | blocking / major | DOC1: `docs/implementation2.0` added to the owner's `git add`; `.gitignore` gets `docs/implementation2.0/kernels/*/build/`; the evidence that `K2/DESIGN.md` and `K3/DESIGN.md` cite inside `build/` moves to `results/` first (as `.txt`), then `build/` is deleted; the `GPU_RUN.md` files and `K1/results/build_cpu.txt` use `$REPO`; scan (a) covers `docs/implementation2.0` and local paths; step 3 and DOC-G1 use `git ls-files --others --exclude-standard -- docs ':(exclude)docs/research'` (the `find` count also counted ignored files); DOC.1 item 1 and DOC.3.3 counts re-measured | `find`: 385 outside research, 348 in `docs/implementation2.0`, 298 in `build/`, 49 trackable; 8 tracked; 367 untracked not-ignored; `du`: 3.5 MB + 10 MB; 8 `build/...` citations in `DESIGN.md`; scan (a) 15 lines |
| arch 13 | major | DOC3: `exclude_docs` adds `implementation/kernels/**` and the five copied-evidence folders; `nav` Plan lists `implementation/README.md` and `implementation/Implementation.md`; step 3 records the build result | Built a copy of today's `docs/` plus `src/` with that config: "No issues found", exit 0, 3.4 s (`docs/implementation2.0/notes/doc/r1build/build.txt`) |
| arch 14, qa 5, cuda 13 | major | DOC1 step 0: one copy table from the scratchpad to `docs/implementation2.0/{specs,notes,seeds,tests,reviews}/` and `notes/doc/`, with skips (venv, paper copies, 12 MB outputs, `.log` renamed); run at plan hand-off; scrub of prefixes. Every scratchpad path in this page replaced (23); DOC.5 kernel links are repo paths. New check 14, gate DOC-G14, test U28, risk row. The `A.md` alias `docs/research/notes/impl_2026-11-08/` is named as a conflict: DOC1 uses `docs/implementation2.0/notes/` (33 uses in the phases, research folder frozen) | grep of the phase pages for scratchpad paths; `ls`/`du` of the scratchpad folders (notes 999 files, 29 MB, of which `venv/` 919 files) |
| cuda 24 | minor | Covered by the qa 4 change (`.gitignore` line, move, delete) | `file`: arm64 Mach-O in `K1/build` |
| (new, from verification) | major | `PLAN_GLOBS = ["docs/implementation2.0/**"]`: plan files skip checks 2, 3, 7, 8 and keep checks 1, 13, 14. Without it check 2 fails as soon as the plan is tracked | stdlib probe of the check 2 rule: 285 of 1,240 backtick paths in `docs/implementation2.0` do not exist; `docs/plan/`: 0 of 4 |
| arch 4, qa 7 | blocking / major | DOC10 keeps ADR 0012 (front section 4 table; arch 4's order). Check 4 also reports duplicate ADR numbers, with a negative test; DOC-G11 updated | `ls docs/architecture/decisions`: 0001-0011 |
| qa 6 | major | Every DOC gate id now carries the prefix: DOC-G1..DOC-G14 (128 replacements, tests and coverage maps included) | grep: all `G<n>` ids in this page were DOC gates |
| qa 9 | major | DOC-G7: passed >= baseline at the start of the task group (DOC1, or DOC13 step 0) - archived tests + new tests; "no test removed" removed; baselines in `runs/doc/baseline.txt` and `baseline13.txt`; I16 and the coverage map follow | DOC removes no test; CLEAN archives tests before DOC13 |
| arch 15 | minor | DOC2 states it is the one owner of the "TRL 1.16" text; A's S06 only asserts absence; E does not edit it | A.md (S06 text) and E.md (E1 item 4) already say this |
| arch 11 | blocking (E) | DOC13: E8.2's default keeps `experiments/phase2/`, so `CLAUDE.md` needs no edit; under the alternative the owner edits the Project rule first, then the run-command example | E.md risk row: "Default: no move (E8.2)" |
| data 26, data 9, data 20 | minor / major / minor (DOC parts) | DOC11 "Known issues": 8 items with sources (P5 population, P2 generous masks, 8 two-scan patients, B1.7 z check, one synthetic shift law, test informed the diagnosis; reader consistency and the TCIA 0396/0510 notes marked `[unverified]`). Card "Limitations" and the reproducibility line ("patient-cluster bootstrap CIs (A1)") updated | `docs/reference/data.md:15,28,31,61`; `package.md:125`; `data_engineering.md:94`; `build_eval_views.py:57,60`; B.md B1.7 / DG17 |
| (consistency with A) | minor | `ci_unit` that phase A writes is `"patient"`, not `"nodule"` (DOC6, owner table, dependencies, U20) | A.md "Words": gate CI = patient-cluster bootstrap, key `group_id` |
| qa 14 | major | One word for each thing: the DOC6 view string is "shifted views", and DOC7 changes README:121 from "(shifted crops, 8 samples each)" to "(shifted views, 8 answers each)"; glossary rows follow front section 0 (miss by the K-mean rule, found image, missed answer, window = B1 code name); X5 expects 15 terms; "crop (engine)" stays (another thing) | README.md:121; docs/glossary.md:53 |
| qa 8 | major (DOC part) | DOC4 registers markers `docs_build` and `post_commit` (A turns on `--strict-markers`) | pyproject.toml:68-71 lists `slow`, `data` only |
| fresh-arch 6 | major (DOC12 part) | DOC12: pre-commit and the `no-ai-attribution` commit-msg hook are DEFERRED, as front section 4 decided, with its trigger | front section 4 "Types, format, hooks" row |
| qa 17 | major | DOC.6 non-blocking metric: re-count sentences over 25 words after the edits (baseline 803 of about 4,965). New text in this revision uses short sentences; the older long table cells were not all split | review count cited, not re-measured |
| qa 28 | minor | Header tags follow front section 0 (`[web: source]`, `[unverified]`, dated `[measured: ..., 2026-11-08]`); 15 vague `[measured: review/grep/ls]` tags dated | - |

Also updated: the box diagram check list, DOC.0 skills rows (revision 3), DOC.8 time (DOC1 2 h, total about 30 h),
DOC.7 two risk rows (scratchpad cleared before step 0; binaries or paths entering git).

Not changed here (owned by another section): the `A.md` notes alias (main agent or A fixer), removal of the
scratchpad paths inside the other phase pages (each phase fixer; DOC1 step 0 and check 14 catch what remains), the
kernel `GPU_RUN.md` edits themselves (DOC1 does them when it runs; front section 5 states the rule).

## Review log (final round 2)

Sources: `docs/implementation2.0/reviews/final_r2_{arch,cuda,data,fresh-arch,qa,rl}.md`, every item addressed to DOC
or naming a DOC part. `final_r2_cuda.md` has none. Each fact was checked against the repo, the scratchpad or the web
on 2026-11-08 before the change.

| Review item | Sev | What changed | Check made |
|---|---|---|---|
| fresh-arch 1 (DOC part) | blocking | U19 no longer asserts `HARD_STOP == 17.5`; U18, check 8 and the sketch read the constants through `POLICY_MODULE`. DOC7 drops "the agent does not change it" and points at front 6.1 step 0 ($19.90, owner approval) and 99_back decision 4. DOC2 adds the rule: the commit that edits `policy.py:58` also changes the runbook Hard stop row. DOC2 rewrites the Sweeper row with the constant names. DOC.1 item 3, DOC.8 owner list and 9.10 decision 5 updated | `policy.py:58` = `10.0, 17.5`; refusal at `policy.py:324-326`, sweep at `588-592`; ledger 133 rows, $16.90; `gpu-run.md:115-118` (Sweeper row names $17.50 and $10); front 6.1 step 0 text |
| fresh-arch 5, arch 3, qa 4, rl 9 (DOC part) | blocking / major | New check 15 `test_adr_numbers_match_front_table` with helpers `adr_table`, `adr_mentions`; four rules (numbers, one file name per number, file name on disk, index title = table title); negative tests with run-time fixture strings; U29; DOC-G11, DOC-G3, I01, M01, M02, coverage maps updated. The renumbering itself is in B and C (qa 4: "Do not touch DOC") | Regex probe of all phase pages after the round 2 fixes: mentioned numbers are on disk or in the table; only 0012 has a backtick file name. My first draft of the negative test text tripped the probe (a literal file name and an out-of-table number), so the fixtures are now described, not written |
| arch 11 (DOC part) | major | DOC13 step 3 sets `POLICY_MODULE = "objdel.infra.vast.ledger"`; check 8 names it | E.md E8.1 move table: `ALARM`, `HARD_STOP` in `infra/vast/ledger.py`; E.md says DOC owns the edit |
| arch 15, qa 12 (DOC part) | minor / major | DOC1 step 0 rows: `impl/c_rev/` (4 files, skip `tmp_tf519/`) -> `notes/c_rev/`; `b1_out/neardup.json` and `shacheck.json` copied, `geom.json` (9.1 MB) and `pix.json` (3.6 MB) listed "not kept". Added `impl/r2rl/` -> `notes/r2rl/` (found while checking: B cites `notes/r2rl/merge_rounding.py`, but the folder sits at `impl/r2rl/`, outside `notes/`) | `ls -la impl/c_rev impl/notes/b1_out impl/r2rl`; grep of the phase pages for `c_rev`, `b1_out`, `r2rl`. C's own `S[P]/` paths and E.0's `nv_catalog.json` path are fixed in C and E |
| arch 16 | minor | `config_files` reads `specs/*.yaml` plus `experiments/**/configs/*.yaml` and raises when empty; check 7 scans `specs/*.yaml`; pin 6 and U12 updated; DOC13 step 3: `configs.md` documents every `specs/*.yaml` after E6 | E.md: `specs/b1_rows.yaml`, `tests/core/test_spec.py` loads every `specs/*.yaml` |
| arch 17 | minor | DOC9 owner table names each phase's ADR numbers (B1 0015, B2 0013, B3 0014, C 0022, D 0016, E 0017-0020, CLEAN 0021); one line: each phase writes its own ADRs, DOC only 0012 | Front section 4 table rows 0012-0022 |
| arch 6 (DOC part, from (c)) | major | DOC9 E row updates `reference/trl-coupling.md` (new item count after E5) | arch 6 fix (c) |
| data 17 | minor | DOC11 known issues 4 and 5 now cite the TCIA page items 2, 6, 7, with no `[unverified]` tag. New item 5b: the 0101 XML replaced 2012-03-21; 044.xml and 191.xml lack ratings. Changed from the review: "ratings are not used here" is only half true. The parser keeps a nodule only when it has a `characteristics` element. The two flagged nodules have none in the local copies, so their outlines are dropped. Items 5 and 5b state that | WebFetch of the TCIA LIDC-IDRI page (8 items, quoted); `lidc.py:140-147`; grep and awk of `data/raw/lidc_xml/tcia-lidc-xml/187/{044,191}.xml`: comment at lines 2599 and 6911, no `<characteristics` before the closing tag at 3043 and 7266 |
| data 5 (DOC part) | major | The sentence "Training-seed variance is not measured; a gain is a single-run result" goes into the card Limitations and the reproducibility page. It is not a data issue, so it is not in the data card's Known issues | data 5 fix text |
| qa 15 (DOC part) | major | DOC.4 order and DOC.8 dependencies: DOC1 and DOC2 run before any phase starts (same as front section 1, time line step 1); DOC2's ownership sentence matches. S08's skip rule belongs to A | Front section 1 item 1: "DOC1 and DOC2 run before any phase starts" |
| qa 25 | minor | DOC.6 non-blocking metric records the round 2 count (993 of 5,203). The new text in this revision uses short sentences | Count cited from the review, not re-measured |
| (found while checking) | minor | DOC3 said "check 7 imports `experiments.phase2.vast.policy`"; it is check 8, through `POLICY_MODULE` | DOC4 check list |
| qa 21 | minor | No change: DOC has no "exam" (the line-10981 hit is CLEAN) | `grep -n -w -i exam DOC.md`: 0 hits |

Not changed here (owned by another section): the B and C ADR renumbering; C's `S[P]/` paths and header; E.0's
`nv_catalog.json` path; front 2.3 target tree (arch 19); A's S08 skip rule (qa 15); the `policy.py:58` edit itself
(front 6.1 step 0, main agent).

## Review log (final round 3)

Sources: `docs/implementation2.0/reviews/final_r3_{arch,cuda,data,fresh-arch,qa,rl}.md`, every item addressed to DOC
(`final_r3_cuda.md` and `final_r3_rl.md` have none). Each fact was checked against the repo, a command or the web on
2026-11-08 before the change.

| Review item | Sev | What changed | Check made |
|---|---|---|---|
| arch 3 | blocking | DOC no longer accepts a red `make test` between E and DOC13. Section 4 order: each E move and CLEAN archive commit fixes its own doc paths, `configs.md` rows (E6), `POLICY_MODULE` (E8.1), `FROZEN_GLOBS` and the lock path (E R-1). Check 8 says the E8.1 move commit sets `POLICY_MODULE`. DOC13 steps 1 and 3 only re-run and check (X9). DOC.7 risk row and DOC.8 dependencies rewritten. DOC13 time 3 h -> 2 h | E.md R-1 (lines 312-325) already holds the rule and the `POLICY_MODULE` edit; E8.1 text (E.md:1114-1115) agrees |
| arch 7 (DOC part) | major | DOC5 "Scope": DOC writes the 63 docstrings already in `src/objdel`; E writes the 59 for the 12 moved files in a change commit before each move (E R-6). Section 4 order and DOC.8 name R-6 | ruff 0.17.1 `--select D100,D101,D102,D103,D104,D106`: `src/objdel` 63 (25 D102, 23 D103, 15 D101); the 12 files 59 (58 D103, 1 D101). E.md R-6 exists (line 345) |
| arch 11 (DOC part) | major | Check 2 states it is the only doc-path test; CLEAN adds none. DOC.8 CLEAN line says so. The CLEAN12 deletion itself is CLEAN's edit | CLEAN.md:771 still defines `test_doc_paths_exist` (owned by CLEAN); 285 of 1,240 plan paths missing (DOC4 measurement) |
| arch 17 | minor | `FROZEN_GLOBS`: dropped `experiments/phase1/**` and `experiments/cuda_checker/**`; added `experiments/cuda_checker/results/**` (11 globs). The records stay frozen through the `lab`, `impl`, `PLAN*`, `RESULTS*`, `EXECUTION` globs. Ran check 7 on the newly checked files: 2 "TRL 1.16" and 2 `blob/v1.16.0` hits in `experiments/phase1/{sample_eval,smoke_grpo}.py`; DOC2 fixes them; check 7, DOC-G5 and I02 counts now 10 + 7. Check 2 has no newly checked Markdown (the folders' `.md` files match record globs). U01 fixtures updated | `ls experiments/phase1 experiments/cuda_checker`; `policy.py:451-452` runs `run_box.sh`; regex probe on the 18 newly checked `.py`, `.sh` and `Dockerfile` files (the 4 `.md` files there stay frozen); `run_box.sh:5` names torch 2.14.1 (legal) |
| (found while checking arch 17) | major | `is_frozen` uses a segment-aware glob, not `fnmatch`: `fnmatch` lets `*` match `/`, so U01's near miss `experiments/phase2/sub/PLAN.md` would be frozen and U01 would fail. New finding F10 | `fnmatch("experiments/phase2/sub/PLAN.md", "experiments/*/PLAN*.md")` = True, Python 3.12.11 |
| qa 11 | major | DOC-G13 and DOC1 step 2 exempt `docs/implementation2.0/Implementation.md` by name, with its own 2 MB limit; baseline cell records 1,480,096 bytes. Split into per-phase files not chosen: check 15 and check 16 read one file, and the site already builds it | `ls -l`: 1,480,096 bytes. Web search for a GitHub size limit on rendered Markdown found no documented number, so the doc makes no claim about GitHub rendering |
| fresh-arch 1 (DOC part) | blocking | New check 16 `test_v1_table_matches_b_s1` with helper `table_after`: rows, item text and the 27.5 s / 23.65 s minutes of the B-S1 table equal the front V1 table; extra B columns allowed; raises on a missing heading or column. Negative tests, U30, M01, DOC-G3 (16 checks), coverage maps, P04 updated. DOC4 time 5 h -> 5.5 h | Headings `### 6.2 V1` (00_front.md:555) and `### B.8.1 Money and GPU minutes` (B.md:1424); the B table today has 6 rows and an extra "at 20.0 s" column, so check 16 would report it until B copies the front table |
| data 11 | minor | DOC11 item 5: the parser places contours by z, but the CT fetch picks slices by SOP UID (`lidc.sop_z`); B1.7 (c) checks 139.xml and DOC11 writes the result. Item 8: "train or val only" | `lidc.py:140-147` (z, `characteristics`), `lidc.py:272-279` (`sop_z` docstring names `fetch_lidc_pixels`); B.md:851 (139.xml check), B.md:1336 ("train or val only") |
| data 7 (DOC part) | minor | DOC11 item 8 and the card Limitations say that 128 val shifted items chose R5 over checkpoint-250 and served the G-probe; Gate B1 is conservative. DOC cites B.6.2's sensitivity line only if B prints one; else "sensitivity not measured" | `full_rl_2026-11-08.md` lines 31 and 49; the 51-of-78 patient count is the reviewer's measurement (cited, not re-run); B.md has no `r5_selection_items` key today, so DOC does not claim one |
| qa 21 (DOC part) | minor | DOC.6 metric: round 3 count 1,165 of 7,034 with its rule written beside it; noted that this rule splits list items, so it is not comparable with rounds 1-2; split gate cells and task steps first ("Why:" lines). New text in this revision keeps to short sentences | `final_r3_qa.md` lines 24-25 and item 21 |

Not changed here (owned by another section): CLEAN12's second test and the CLEAN task rule for archive commits
(CLEAN); E R-1 and R-6 text (already in E); the B-S1 table itself and the B.6.2 sensitivity line (B); the front 6.4
rule 7 eval-view wording (data 12, front).

## Review log (consistency pass, round 4 prep)

| Item | What changed |
|---|---|
| Check 16 | `test_v1_table_matches_b_s1` becomes `test_rental_tables_only_in_front`: a table with a "$" header outside front section 6 is reported. B no longer copies the V1 table, so the old rule would raise. Negative tests and U30 follow. Dry run on the assembled plan: 0 reports |

## Review log (final round 6)

Cross-file consistency fixer, round 6. Sources: `reviews/final_r6_{cuda,rl,data,arch,qa,fresh-arch}.md` and their
`*_verdict.md` files. Only the CONFIRMED items assigned to this fixer are listed.

| Item | Severity | Change | Evidence |
|---|---|---|---|
| arch-5 | major | DOC1 step 0: bracket every literal mention of the path patterns in the phase pages and review logs, then re-run the check 14 scan on the built pages until it prints nothing | check 14 scan on the sources: 0 hits |

## Review log (final round 7)

Fixer for the DOC page, round 7. Sources: `reviews/final_r7_arch.md` item 5 and `reviews/final_r7_arch_verdict.md`.
Only the CONFIRMED item assigned to this fixer is listed.

| Item | Severity | Change | Evidence |
|---|---|---|---|
| arch-5 | major | DOC1 kernel bullet: the scrub now covers short scratchpad paths too, with the check 14 patterns as its measure (16 lines in 7 files) and one replacement per hit: `REPO=$(git rev-parse --show-toplevel)`, `K<n>=$REPO/docs/implementation2.0/kernels/K<n>`, work files to `$REPO/data/k1_bundle`, DESIGN.md folder lines to `docs/implementation2.0/kernels/K<n>/`, `nv_catalog.json` to `docs/implementation2.0/notes/doc/`, deep findings to `docs/implementation2.0/notes/`; the `S[P]=` line dropped. DOC1 step 0 re-scan now also runs on `docs/implementation2.0/kernels/**`, so check 14 is green at the DOC4 commit | `grep -rnE` with the check 14 patterns on `docs/implementation2.0/kernels`: 16 lines in 7 files (`K1/GPU_RUN.md:18,19,20,25,31,126`, `K1/DESIGN.md:3,19,242,304`, `K1/results/build_cpu.txt:1`, `K2/GPU_RUN.md:33`, `K2/DESIGN.md:3`, `K3/GPU_RUN.md:15`, `K3/DESIGN.md:3,9`); `git check-ignore` confirms `data/k1_bundle/` is ignored (`.gitignore:13`); the step 0 table copies `nv_catalog.json` to `notes/doc/` and the deep findings to `notes/` |

## Review log (final round 8)

Cross-file consistency fixer, round 8. Sources: `reviews/final_r8_{cuda,data,arch,qa}.md` and their `*_verdict.md`
files. Only the CONFIRMED items assigned to this fixer are listed. Each fact was checked in the repo first.

| Item | Severity | Change | Evidence |
|---|---|---|---|
| data-1 | blocking | DOC1 step 0 copies `impl/r8data/` to `docs/implementation2.0/notes/r8data/`; DOC1 Files corrects `docs/plan/Research2.0.md:196` "All images" ceiling 0.734 -> 0.723 (the "76 %" sentence stays: 0.080 of 0.105) | `Research2.0.md:194-196`; `notes/r8data/legacy.py`, `legacy_found.py` |

## Review log (final round 9)

Cross-file consistency fixer, round 9. Sources: `reviews/final_r9_{cuda,rl,data,arch,qa,fresh-arch}.md` and their
`*_verdict.md` files. Only the CONFIRMED items assigned to this fixer are listed. Each fact was checked in the repo first.

| Item | Severity | Change | Evidence |
|---|---|---|---|
| data-1 | major | DOC6 `headline.py`: `--test PATH --val PATH` or `--run ID` (R5 short form); reads schema 1 (`ci95`, `item`) and schema 2 (`gate_ci`, `patient`); refuses mixed pairs and null gate CIs; Makefile `RUN=` or `TEST=`/`VAL=`; I03/I04 schema cases; anchor "A.md, A1 step 5"; owner row A and dependency line name A1 step 5b / A3 | A.md A1 step 3 (no bare `ci95` in schema 2), A3 (`{rd}/eval_val`, `{rd}/eval_test`); `ls runs/r5/eval` |

## Review log (final round 11)

| Item | Severity | Change | Evidence |
|---|---|---|---|
| fresh-arch-3 | major | ADR 0015 is "scheduled RL rows" only; DOC9's polygonizer row states the path before and after B1.4a | B.md B1.4a |


# Phase CLEAN

## CLEAN.R. Review log (revision 2, from docs/implementation2.0/reviews/CLEAN_phase_review.md)

Each item was re-checked against the code or a command on 2026-11-08 before acting.

| # | Severity | Item | What I did |
|---|---|---|---|
| 1 | blocking | CLEAN4 golden would publish RIGA masks | Accepted. Golden stores no masks: LIDC per-case hashes and vertex counts, one RIGA total and one hash. `crops_oracle.py` is always ported; "size rule" deleted. CLEAN-G16 guard added. Stricter than the review: no per-case RIGA vertex sum. |
| 2 | blocking | `rewards.py:107` is live (`fast.py:633`) | Accepted (verified). CLEAN6 keeps the line, changes the comment, adds an assertion; CLEAN-G17. |
| 3 | blocking | CLEAN5 archives tests of live code | Accepted (verified at test_p2_signal.py:91-120,189-196; test_p2_probes.py:174-180). Recorder + `p2_lab_refs.json`; tests rewritten, not archived; CLEAN-G18. |
| 4 | blocking | vastai pin would downgrade Pillow | Accepted (verified: `pillow==12.2.0` pin, `.venv` Pillow 13.0.0). `uv tool install`, no `pyproject.toml` change, preflight version check. |
| 5 | major | CLEAN0 count wrong (spaces) | Accepted (verified 1,038 vs 802). `-z`, `-C ROOT`, `phase0` count, vulture count, preflight checks added. |
| 6 | major | CLEAN2 drops `gc.freeze()` | Accepted (verified serve.py:172). Line 50 becomes `gc.freeze()`; spy test; troubleshooting row. |
| 7 | major | CLEAN6 test list incomplete | Accepted (verified by grep). Node-id list from grep; wiring test kept without the Controller; CLEAN6 before CLEAN7. |
| 8 | major | CLEAN5 step 3 assertion and range | Accepted (verified line 52 and the 66-100 span). Delete by name. |
| 9 | major | Stale paths in docs and docstrings | Accepted (verified: 13 docs, 22 src, 9 tests files). Sweep step with exempt records; `test_doc_paths_exist`. |
| 10 | major | `git add -A` sweeps owner's docs work | Accepted (verified git status). Clean-tree precondition in the task rules; `git add -p`. |
| 11 | major | Archive README and tag wording | Accepted. Tag `pre-clean-2026-11-08` + push; README text and `git restore --source` command. |
| 12 | major | Whitelist fails lint; test skips | Accepted, different fix: no whitelist file; `[tool.vulture]` `ignore_names` in `pyproject.toml` (removes the F821/B018 problem); test fails, not skips. |
| 13 | major | CLEAN-G9 baseline scope | Accepted. Measured the exact scope: 56 lines (not 6); 4 after the ignores. CLEAN0 records it. |
| 14 | major | No $0 bundle import check | Accepted, corrected fix: `-I` ignores `PYTHONPATH` and the copy lacks the built `.so`, so the test copies `_gpucheck` and sets `sys.path` in code, then asserts every module came from the copy. Collect-only dropped (it would resolve to the editable install). In CLEAN1. |
| 15 | major | Hidden conflicts with B and C | Accepted. Wave 1a = CLEAN2, 6, 7 first; Dependencies and risk row updated. |
| 16 | major | CLEAN11 import list incomplete | Accepted, plus one more importer found: `experiments/phase1/bench/xarch_hashes.py:45`. |
| 17 | major | CLEAN1 treats the symptom | Accepted. `/.claude/worktrees/` into `.gitignore`. |
| 18 | minor | `boxes_meet` in a loader | Accepted. `objdel.geometry.box`; `max_examples=2000`. |
| 19 | minor | pyyaml wording | Accepted (moot after 4). Deleted. |
| 20 | minor | SkyPilot evidence | Accepted. Keep-list reason is governance; SkyPilot row added to section 2. |
| 21 | minor | CLEAN-G13 FAIL = FAIL | Accepted per-check comparison. Rejected the config fix: `rl_shifted.yaml` cannot load in preflight (`TypeError`, measured), so its comment is right; Makefile default and CLEAN-G13 use `r1_qwen3vl2b.yaml`. |
| 22 | minor | Archived `prod.py` will not run | Accepted. README says run it from the tag. |
| 23 | minor | Same test twice | Accepted. Kept in CLEAN1 with a real-repo assertion; removed from CLEAN12. |
| 24 | minor | CLEAN9 before CLEAN10; hash with `-I` | Accepted. Order rule and precondition; `-I` hashing. |
| 25 | minor | CI skips the CLEAN4 gate | Accepted. Mac gate with `-rs`, must show passed and no skip. |
| 26 | minor | `sample_vllm.py` missing from CLEAN7 | Accepted (verified line 6). |
| 27 | minor | GPU versus CPU | No action needed; added the `gc.freeze` row to section 3.4. |
| new | (own) | CLEAN9 step 5 would delete fake-function parameters | Replaced with the 4 measured real findings; parameters go to `ignore_names`. The `reward_gates.py:138` call stays (it draws from `rng`). |

---

## Review log (final round 1)

Sources: `docs/implementation2.0/reviews/final_r1_{arch,cuda,data,fresh-arch,qa,rl}.md`, items addressed to CLEAN
(`data` and `rl` have none). Each fact was checked against the code, the front matter and the D and E phase texts on
2026-11-08 before the change.

Skills used in revision 3 (all invoked with the Skill tool):

| Skill | What it changed in revision 3 |
|---|---|
| nvidia-skills:nvidia-skill-finder | No NVIDIA product choice in CLEAN. cuCIM wording now says "not adopted" (front section 4, D.2). No install. |
| engineering:code-review | Re-read CLEAN11 against D3: revision 2 would archive `certify.py`, which RIGA calls after D3. CLEAN11 is now check-only. |
| engineering:debug | Ran the revised bundle probe on today's tree: with `-I` only, 33 of 44 first-party modules came from the repo through the editable install; with `-I -S`, 0 of 44; the canary fails as it should. The test text changed to match. |
| engineering:system-design | Aligned the waves with the front time line: Wave 2 after D3, Wave 3 after E and D7 step 6, DOC13 before Waves 2-3. |
| engineering:architecture | ADR number 0021 from the front table; CLEAN14 mirrors front section 4 (KEEP, DEFERRED with trigger, NOT ADOPTED). |
| engineering:documentation | Repo paths replace the scratchpad paths; gate ids carry `CLEAN-`; undefined tag forms fixed. |
| huggingface-skills:huggingface-best | No model choice in CLEAN; no change. |
| Web tools (WebSearch, WebFetch) | Loaded. No new web fact was needed: every change rests on repo files or a local measurement. |

| Review item | Sev | What changed | Check made |
|---|---|---|---|
| arch 3 | blocking | CLEAN11 rewritten as a check-only task: D3 owns the move to `tests/oracles/{polygonize,bitquad}_oracle.py`; `certify.py` stays live; nothing archived; new permanent `test_polygonizer_seam_is_device_only` (guard `WAVE >= 2`, module paths through `RINGS` so it survives E3a). Section 2 K1 row, 3.1 diagram, 3.4 row, CLEAN-G14, CT-60, CT-61, 9.7 updated | D.md D3 file list (lines 308-345); D.8 "CLEAN (seam and pool)"; `data/riga.py:226` `certify_many` per D3 |
| arch 4, fresh-arch 7, qa 7 | blocking / major | CLEAN14 file is `0021-library-keep-list.md` (front ADR table), not "next free"; CT-45 checks 0021 and unique ADR numbers. fresh-arch 7's own table (0018) is superseded by the front table | `ls docs/architecture/decisions`: last is 0011; front section 4 table |
| arch 5 | major | CLEAN2 states that its text wins; E1.3 only checks (`--procs` gone, one `gc.freeze()` after the last preload) | `serve.py:172` `gc.freeze()` inside `fork_pool` [code]; E.md item 3 "E1.3 = CLEAN2" |
| arch 6 | major | CLEAN9 step 1 and 3.2 name E's destination `geometry/box.py -> geom/box.py` | front 2.3 placement table row "A, C, D, CLEAN new modules"; E.md line 184 |
| arch 10, fresh-arch 5 | major | cuCIM: CLEAN.0 row, section 2 row and CLEAN14 say K2 EDT on the device (D7), cuCIM not adopted (optional box oracle). Vast row unchanged (`vastai` CLI 1.11.0 as `uv tool`, no SDK), now matching front section 4 | front section 4 EDT and GPU orchestration rows; D.2 |
| fresh-arch 6 | major | CLEAN14 DEFERRED list carries the front "none" rows with triggers: Trackio, `events.jsonl`, `uv.lock`, pyright, `ruff format`, pre-commit hook (candidate DOC12), pytest-xdist, view bank. CLEAN takes no new task for them (front chose DEFERRED, not CLEAN9) | front section 4 "Owning task" column |
| arch 12 | major | CLEAN8 step 4: E6.6 owns the move of the 12 old plans (16 minus the 4 ported `rl_shifted`, `rl_full`, `r5_shared`, `r1_qwen3vl2b`); CLEAN8's LFM-coupled scripts go in that commit; `hf_sample.py` left to E4. The wrong "CLEAN item 11" cross-reference is in E.md (E's fixer) | `ls experiments/phase2/configs` = 16 files [measured]; E-G9 names the 4 ported specs; E.8 "E4 archives sample_vllm.py and hf_sample.py" |
| arch 17, qa 21 | minor | CLEAN12 extends A's `tests/test_repo_hygiene.py` (S01-S11 kept); `test_src_never_imports_experiments_or_tests` is the only test of that edge; CLEAN1 is the only owner of the bundle drop rule; CLEAN16 step 1 is the single owner of the pool deletion | A.md A.9.4 (S01-S11), A.md line 932; E.md lines 429-430 |
| arch 17 (pool timing) | minor | CLEAN16 step 1, CLEAN2 step 6, 3.4 row, CLEAN.8 and CT-63 now say "after D7 step 6" (not D5); E3b archive branch (E R-2) noted | D.md D7 step 6 and D.8 ("runs after D7 step 6, not after D5"); E.md line 254 |
| arch 19 | minor | CLEAN1 step 5, CT-30, CT-31 and sketch B read the module list from `P.IMPORT_CHECK` plus a new `P.LIVE_CLIS` constant; E8.2 must list the test among the strings it rewrites | `policy.py:69` `IMPORT_CHECK`, used only at `policy.py:400` [code] |
| cuda 13, arch 14, qa 5 | major | Every scratchpad path replaced by its repo path (`docs/implementation2.0/{kernels,notes,reviews,tests}/`); header no longer defines SP; the note that DOC1 makes these copies is in the tag paragraph | a grep for scratchpad paths and absolute home or temp paths in this file finds 0 after the edit; `docs/implementation2.0/kernels/K1/DESIGN.md` exists [measured: ls] |
| cuda 24 | minor | CLEAN.5 says DOC1 keeps the kernel `build/` folders out of git (front section 5); CLEAN deletes nothing there | `du -sh`: K1/build 3.5 MB, K2/build 10 MB, K3/build 24 KB [measured] |
| qa 6 | major | Every CLEAN gate id is `CLEAN-G1`..`CLEAN-G18` (every use in this file); CLEAN.spec constraint ids written as "CLEAN.spec constraint Cn" so they do not read as Phase C tasks; "Direction phase G3" and the old test lists' "G1..G16" left as written | front section 0 id prefix rule |
| qa 9 | major | CLEAN0 step 1 records `runs/clean/baseline.txt` at phase start and at each wave start; CLEAN-G2, CT-70 and CT-80 compare with that baseline; 1447 is a 2026-11-08 reference only; task rule added | A-D add tests before CLEAN (front time line) |
| qa 16 | major | Task rule: `data` tests run only through `make test-data` (`-I -m data`), `make test` deselects them; CLEAN adds no target (A or D8 does). CLEAN4 gate command now uses `-I -m data`; L4 header and 9.2 updated | `Makefile:19-20` has no `-I` today; `Makefile:30` preflight already uses `-I` [code]; A "All" gate and D8 own the target |
| qa 28 | minor | `[pypi: ...]` defined as a dated `[web: pypi.org JSON]` tag; `[measured: list 2, this session]` replaced by a dated file tag | front section 0 tag list |
| qa 17 | major | New and rewritten text uses short sentences. A full split of older CLEAN text was not done here | - |
| (own) | major | The revision 2 bundle test could never pass on the Mac: the editable install's `.pth` finder served 33 of 44 modules from the repo, and namespace packages (`tests`, `tests.oracles`, `experiments.phase1`, no `__init__.py`) have `__file__` None. CLEAN1 step 5 and sketch B now run the child with `-I -S`, pass the site-packages folders as an argument, and check namespace `__path__` entries. New risk row | probe run on today's tree with the CLEAN1 drop rule simulated: `-I` 33 foreign, `-I -S` 0 foreign of 44, canary detected [measured: 2026-11-08]; `ls tests/__init__.py` absent; `_objdel_repo.pth` holds the repo root |
| (own) | minor | CLEAN5: A2 may already have created `archive/phase2_lab/A4/`; a plain `mv` of the lab folder would nest it as `archive/phase2_lab/lab/`. Move child by child and check | A.md S05; `mv` semantics |
| (own) | minor | CLEAN5 notes that the Wave 1 lab move covers D's requests (broken `polycore` imports, lab scripts using `bitquad` or `build_target` without `min_area`) before D3/D4 | `lab/A6/prof_reward.py:18`, `lab/A1/a0_ceiling.py:9` import `polycore` [code]; D.8 |
| (own) | minor | CLEAN.8: Wave 1 before D3; Wave 2 needs D3 and D-G4; Wave 3 needs E2, E6 and D7 step 6; DOC13 runs before Waves 2-3 (front time line step 8), and Waves 2-3 touch only archive READMEs and lines that `test_doc_paths_exist` flags | front section 1 time line |

Not handled here (owned by another section): E6.6/E8.2 cross-references to "CLEAN item 11" and "CLEAN items 5, 7, 12"
(E), the E move-table rows for `geometry/box.py` and the D modules (E), E8.1's duplicate drop bullet (E), front 2.3
tree entries (front), DOC1's `.gitignore` line for kernel `build/` folders and the copies into `docs/implementation2.0/`
(DOC).

## Review log (final round 2)

Sources: `docs/implementation2.0/reviews/final_r2_{arch,cuda,data,fresh-arch,qa,rl}.md`, items addressed to CLEAN or
naming a CLEAN task. Each fact was checked on 2026-11-08 against the code at the working tree and against the B, E and
front phase texts before the change. `cuda` names CLEAN16 only in passing (D-G17, a D item): no CLEAN change.

Skills used in revision 4 (all invoked with the Skill tool):

| Skill | What it changed in revision 4 |
|---|---|
| nvidia-skills:nvidia-skill-finder | No NVIDIA product choice in these items. No catalog skill fits model-table cleanup or test rules. No install. |
| engineering:code-review | Checked each changed line against the code: `configs.py:35,70-105,189`, `g_probe.py:86,93`, `prereg.py:63-64`, `train_rl.py:41,53`, `train_sft.py:339-351,384`, `preflight.py:35,49-56`, `hf_push.py:12`, `policy.py:58,372-383`, `tests/test_p2_scripts.py:306,329`. Found that CT-88 asserted `hf_sample.py` absent, against CLEAN8 step 4 (fixed). |
| engineering:debug | Traced what breaks after CLEAN8 at run time: deleting `_LORA` breaks `rl_lr` and `test_prereg_lr` (E1.1). Traced the steps-plan `TypeError` forward to front 6.1 step 4 (`b1_rows.yaml`), flagged in CLEAN9 step 3. |
| engineering:system-design | Placed Wave 1 inside the owner tag `v1-b1` (front order rule 2) and kept the `HARD_STOP` edit with front 6.1 step 0. |
| engineering:architecture | CLEAN8 now uses E's seams (`TrainerFlags.grpo`, `policy/models.py`); one archive folder for old plans (`archive/phase2_runs/configs/`). |
| engineering:documentation | One name per thing: `b1_rows.yaml` is the live RL plan after B1.5; "exam" replaced; the four open test questions became rules. |
| huggingface-skills:huggingface-best | No leaderboard query. Checked "a new model = one row" after CLEAN8: it holds, because `_LORA = frozenset(MODELS)` (a 4B row joins it by itself). |
| Web tools (WebSearch, WebFetch) | Loaded. No new web fact was needed: every change rests on repo files or the phase texts. |

| Review item | Sev | What changed | Check made |
|---|---|---|---|
| arch 9 | major | CLEAN8 rewritten on post-E names. Step 1: the test cap is `flags.grpo={"max_completion_length": 32}` (`TrainerFlags.grpo`, passed as `flags` to `build_trainer`). Step 3: in `policy/models.py` keep `_LORA` as `frozenset(MODELS)` and keep `_GRAD_CKPT`; delete only the full-fine-tune branches (`param_groups` and its caller, the float32 branch of `dtype`, the full-checkpoint branch of `_load_from_checkpoint`). `rl_lr`, `test_prereg_lr` and `G_PROBE["lr_full"]` stay. CT-62 matches | E.md E1.1 (`rl_lr` reads `_LORA`), E2.1 (`TrainerFlags.grpo`; "there is no separate `max_completion_length` field"), E2.2 signature, E5 table (`policy/models.py`); `configs.py:76-77,83-105`; `train_sft.py:384` calls `param_groups` |
| arch 14 | minor | Stale `RlSpec.max_completion_length` and `build_trainer(spec, ...)` removed from CLEAN8 and CT-62 (the "CLEAN.9" place of the review is CT-62 in this file). Step 4 names post-E homes of `HF_RATE` and the `hf_push` licence text (`objdel.infra.preflight`, `objdel.infra.hf` after E8); step 2 finds test files by grep, because E8.3 may move them | grep of this file: 0 `RlSpec` outside "no `RlSpec`" notes; E.md E8.1 thin CLI rule |
| arch 1 | blocking (E, B) | CLEAN side only: task rule "Live RL plan"; CLEAN8 step 4 counts the 12 old plans as 16 minus `rl_shifted.yaml` (B1.5) minus 3 ported, and names one archive folder `archive/phase2_runs/configs/`; CLEAN9 step 3 says E6.8 moves the default to `specs/b1_rows.yaml`; CLEAN2 step 7, CT-84 and CLEAN10 step 1 cover `b1_rows.yaml`; CT-88 names the folder | B.md B1.5 lines 622-625 (archive to `archive/phase2_runs/configs/`); E.md E6.2, E6.6, E1.5, U19 |
| data 3 | major | CLEAN10 step 4 keeps every `data/eval_views/lidc-val-*` folder until B1.0 records the `lidc-val-shifted` hashes in `runs/a3/val_pre_rebuild.json`. CT-52 checks that these folders and the `target_set_sha256` are unchanged | `ls data/eval_views`: `lidc-val-rescaled`, `lidc-val-rotated`, `lidc-val-shifted`, `probe-precision`, `probe-topology` [measured, 2026-11-08]; B.md B1.0 "Keep the R5 val view provable" |
| qa 23 | minor | CLEAN.9.1: Q-T1..Q-T4 are now rules (deeper `lab/` ships; file filter in `code_files()`; exempt `data/`, `build/`, `.venv/`, `.claude/`; skip with a reason while `docs/` is untracked). CLEAN1 step 1 adds the filter `[f for f in dict.fromkeys(files) if (ROOT / f).is_file()]`; CT-26 now passes with it; table rows 9 and 10 updated. 99_back decision 32 must add CLEAN.9.1 (back's edit) | `policy.py:372-383` (`code_files` uses `ROOT`, no file filter today) |
| qa 21 | minor | "live exam harness" -> "live test-evaluation harness" (section 1). No other "exam" in this file | grep `\bexam\b` = 0 after the edit |
| rl 18 | minor | New CLEAN9 step 6: fix the `configs.py:35` docstring so it states the real RL rate (G_PROBE 5e-5 LoRA, 5e-6 full) and that the 1e-5 / 1e-6 default is test-only; the default itself is left to E1.1 (`rl_lr`). Files list updated | `configs.py:35,189`; `g_probe.py:86`; `prereg.py:63-64`; `train_rl.py:41,53` (uses `make_trainer`) |
| fresh-arch 4 | blocking (front) | Not a CLEAN task: the front chose 6.1 step 0 (main agent edits `policy.py:58`). New task rule: CLEAN does not touch `ALARM, HARD_STOP`; CLEAN1 and CLEAN13 commit apart from that edit | front 6.1 step 0 (line 499); `policy.py:58` `ALARM, HARD_STOP = 10.0, 17.5` |
| fresh-arch 7 | minor (D) | CLEAN6 step 5 states that no live run attaches the Controller, so D-G14 has no audit leg by default; D adds it back only on the `--audit` branch | `g_probe.py:93` builds `TraceStep` only; `grep Controller experiments/phase2/*.py` = 0 [measured] |
| fresh-arch 8 | minor (front) | New task rule: Wave 1 is part of owner tag `v1-b1`; every Wave 1 commit lands before the tag, so the tag's checks run on CLEAN2's and CLEAN6's edits | front order rule 2 (lines 163-167) |
| (own) | minor | CT-88 asserted `hf_sample.py` absent, but CLEAN8 step 4 leaves it to E4 (after E-G10's box check). Removed from CT-88 | E.md lines 836, 859 |
| (own) | minor | CLEAN9 step 3 flags for front and B: `make preflight CONFIG=.../b1_rows.yaml` (front 6.1 step 4) raises the measured `TypeError` if `b1_rows.yaml` is a steps-only plan; CLEAN does not change `preflight.load` | `preflight.py:49-56` (`RunArgs(**raw)`); B.md B1.3 wires `quality.failures` into preflight but does not say how a steps plan loads |

Not handled here (owned by another section): front section 3 step 1 path for `_LORA` (front, arch 9); B1.5/B.9 I22/B.4
archive folder name (B, arch 1: B.md line 1365 and I22 still say `archive/phase2_configs/`); E6 ports (E, arch 1);
99_back decision 32 row for CLEAN.9.1 (back, qa 23); D-G14 and D.1 item 9 audit text (D, fresh-arch 7); DOC U19
literal 17.5 (DOC, fresh-arch 4); how preflight loads `b1_rows.yaml` (B and front).

## Review log (final round 3)

Sources: `docs/implementation2.0/reviews/final_r3_{arch,cuda,data,fresh-arch,qa,rl}.md`, items addressed to CLEAN or
naming a CLEAN task. `cuda`, `data` and `rl` name no CLEAN task. Each fact was checked on 2026-11-08 against the
working tree and against the B, DOC, E and front phase texts before the change.

Skills used in revision 5 (all invoked with the Skill tool):

| Skill | What it changed in revision 5 |
|---|---|
| nvidia-skills:nvidia-skill-finder | No NVIDIA product choice in these items. No catalog skill fits doc tests or archive READMEs. No install. |
| engineering:code-review | Checked each changed line: `callbacks.py:1-11` (the Controller reads TRL's held-out eval), `rewards.py:41,97,139-140` (`recorded`), `hf_push.py:12,60`, `dataset.py:235,240`, the `grpo_rows` callers, `experiments/phase1/lab/README.md` (19 lines). Found that the archive README text fails DOC U23. |
| engineering:debug | Traced each CLEAN commit through `make test` with DOC4 live: CLEAN12's doc-path test fails on the plan files; CLEAN5's `extend-exclude` edit breaks DOC I14; the CLEAN5 lab move puts two files at `archive/phase1_lab/README.md`. |
| engineering:system-design | Placed the docs rule inside each CLEAN commit (E R-1), not in DOC13. Wave 1 stays inside the tag `v1-b1`. |
| engineering:architecture | One rule, one test: DOC4 check 2 owns doc paths; DOC check 11 owns the archive README parts; CT-38 keeps only the CLEAN parts. |
| engineering:documentation | Archive README "Retiring commit" names the tag's commit hash; the lab's own README becomes "What it was". New text uses short sentences (qa item 21). |
| huggingface-skills:huggingface-best | No model choice. Checked that CLEAN8 leaves the HF card licence file to E8.1 (E U50). |
| Web tools (WebSearch, WebFetch) | Loaded. No new web fact was needed. TRL's `num_generations_eval` was checked in the local TRL source instead [code: .venv trl/trainer/grpo_config.py:489]. |

| Review item | Sev | What changed | Check made |
|---|---|---|---|
| arch 3 (CLEAN part) | blocking | New task rule "Docs stay green on every CLEAN commit". Each archive, move or rename commit rewrites the backtick paths in DOC4 `doc_files`, keeps DOC config assertions true, leaves `FROZEN_GLOBS` and `POLICY_MODULE` unchanged (or adds an empty glob to `ALLOW_EMPTY_GLOBS`), adds new doc pages to the indexes, and passes DOC check 11. CLEAN5 step 4 says its sweep is that rewrite. CLEAN.8 names check 2 | DOC.md task order ("Each E move commit and each CLEAN archive commit fixes its own doc paths (E R-1)"); E.md R-1 "Docs follow the move"; DOC1 `doc_files`, `FROZEN_GLOBS`, U01, X9 |
| arch 11 | major | CLEAN12 step 4 deletes `test_doc_paths_exist`; DOC4 check 2 is the only doc-path test. Step 5 gate checks the test is absent. CT-15 is now the stale-path scan only. Q-T3 and Q-T4 retired. Coverage row for CLEAN12 drops CT-15. CLEAN-G9 never listed the test, so it needs no change | DOC1 `PLAN_GLOBS`: 285 of 1,240 plan paths missing [DOC.md]; DOC4 check 2 exemptions; DOC U06 tests the checker |
| fresh-arch 6 (CLEAN part) | minor | New CLEAN6 step 4a: the Controller is the only reader of TRL's held-out eval, no live run attaches it, and B1.5 `--heldout` (TRL eval, `num_generations_eval=4`) replaces it. CLEAN6 keeps `grpo_rows` and `grpo_rows_views`; B1.5 never reads the deleted `TRLReward.recorded` | `callbacks.py:1-11`; `g_probe.py:93` (TraceStep only); B.md B1.5 held-out bullet; `grpo_config` has no `eval_*` key [code: grep] |
| arch 9 (CLEAN part) | minor | CLEAN8 step 4: the LFM text is a docstring line; CLEAN8 does not touch the Apache file; E8.1 moves `licenses/`; E U50 must pass on the CLEAN8 commit | `hf_push.py:12` (LFM line), `hf_push.py:60` (`APACHE` path); E.md E8.1, U50 |
| arch 1 (CLEAN part) | minor | CLEAN9 step 3's round-2 flag marked resolved: front stamps `b1_launch.yaml` (RunArgs with `rows_plan`); plan files never follow `preflight CONFIG=`. The Makefile default stays `r1_qwen3vl2b.yaml` in Wave 1. CLEAN10 step 1 also greps `b1_launch.yaml` (its view is `lidc-bundle-v3`) | 00_front 6.1 step 4 and section 0; B.md B1.3 (preflight wiring), B1.5 (`b1_launch.yaml` keys) |
| qa 21 | minor | New and changed CLEAN text uses one rule per sentence. Older text is not split here | - |
| (own) | major | CLEAN5 step 7's "Retiring commit" text had no hash, so DOC check 11 (U23: "the word commit without a hash fails") would fail on the CLEAN5 commit. The README now names `git rev-parse --short=12 pre-clean-2026-11-08`. CT-38 checks that hash and the tag; the four-part check stays with DOC check 11. CLEAN-G12 updated | DOC.md U23; DOC8 check 11 |
| (own) | minor | `experiments/phase1/lab/README.md` moves to the path of the new `archive/phase1_lab/README.md`. CLEAN5 step 7 writes its folder table under "What it was". `experiments/phase0/README.md` goes to `docs/research/records/phase0/` | `git ls-files` finds those 2 READMEs in the moved trees [measured] |
| (own) | minor | CLEAN5 step 6 drops `experiments/phase1/lab` from ruff `extend-exclude`; DOC I14 asserts it is there. The CLEAN5 commit edits that assertion | DOC.md I14 (DOC.9.3) |
| (own) | minor | CLEAN14 commit adds the ADR index row and the `nav` entry, so DOC4 checks 4 and 13 pass | `docs/architecture/decisions/README.md` exists [measured]; DOC3 `nav` covers `architecture/**` |

Not handled here (owned by another section): B1.5 creates `archive/phase2_runs/configs/README.md` with the headings
"Was / Replaced by / Commit" and a `<commit>` placeholder; DOC check 11 (four parts, 7-40 hex) will report it (B, DOC).
`archive/cpu_checker/README.md` has "What replaced each use" and no "Retiring commit" heading; DOC check 11 reports it
unless DOC8 fixes it (DOC). DOC I14's `extend-exclude` text should say "after CLEAN5: `archive` only" (DOC). The DOC
side of arch item 3 is already in DOC.md (task order rule and the DOC.7 row) [checked: DOC.md].

## Review log (final round 4)

Cross-file consistency fixer, round 4. Sources: `reviews/final_r4_{cuda,rl,qa,arch,fresh-arch}.md` and their verdicts.

| Item | Severity | Change | Evidence |
|---|---|---|---|
| arch-7 | major | CLEAN8 step 4: the E6.6 commit also removes preflight's `projection` check; E6 step 8 holds the list and the test | `preflight.py:216-217` |
| fresh-arch-3 | major | CLEAN6 step 4a: `--heldout-items 64`, `num_generations_eval=8` | `configs.py:177-181`; `rewards.py:132-135` |
| cuda-1 | blocking | CT-72, CT-55 and the task rules no longer diff `cuda/`, `CMakeLists.txt` and the two GPU test files against the Wave 1 tag. New task rule (CLEAN.4): no CLEAN commit touches these paths (CLEAN.spec C1 read as a rule on CLEAN commits), and each CLEAN commit subject starts with `CLEAN<n>:`. CT-73 checks the staged files of each CLEAN commit. CT-72 and CT-55 run `git log -E --grep='^CLEAN[0-9]+:' --name-only pre-clean-2026-11-08..HEAD -- <paths>`, so D's edits are not counted. The `build/.gpucheck-stamp` age rule is dropped. CT-72's freeze check compares with the record from the current wave's start (`freeze.txt`, `freeze_w2.txt`, `freeze_w3.txt`, added in CLEAN0 steps 2 and 5), and CLEAN's own `pyproject.toml` hunks add only `vulture`, `shapely`, `pyyaml`. CT-85 compares Wave 3 values with `before_w3.json`. | D.md D1 (`cuda/gpucheck.cpp` NVTX, `CMakeLists.txt` `find_file`), D2 (engine symbols deleted), D8 (`skipif` at `tests/test_cuda_checker_gpu.py:27` [code], GPU cases in `tests/test_gpu_scoring.py`); `Makefile:3-4,15-17` (the stamp depends on `cuda/` and `CMakeLists.txt`); the `git log` filter was tried read-only on the repo |
| arch-8 | major | CLEAN0 step 4 and CLEAN-G13 pick the preflight config from the tree. While `experiments/phase2/configs/r1_qwen3vl2b.yaml` exists (before E6.8), they use it. After E6.8 (always in Wave 3), they use `specs/r1_qwen3vl2b.yaml` and the stamp path that preflight prints. The gate compares `pass` only on check names present in both maps, against the start record of the current wave, and lists names that are in one map only (E6.6 removes "CODE projected GPU hours and dollars"). CT-43 holds only until E6.8, when E U19 replaces it. CT-80 and CT-85 follow the same rule | E.md E6 step 6 (archives the 4 ported originals), step 8 (preflight takes a `RunSpec`, default `specs/b1_rows.yaml`, `projection` removed), E-G9/U20 (`specs/r1_qwen3vl2b.yaml`); `preflight.py:73-75` (`stamp_path` = `runs/preflight/<stem>.json`), `:244-251` (`CHECKS`), `:284` (path print) [code] |

## Review log (final round 5)

Cross-file consistency fixer, round 5. Sources: `reviews/final_r5_{cuda,rl,data,arch,qa,fresh-arch}.md` and their
`*_verdict.md` files. Only the CONFIRMED items assigned to this fixer are listed.

| Item | Severity | Change | Evidence |
|---|---|---|---|
| cuda-1 | blocking | CLEAN1 step 1: note that D1 step 0 later adds a keep line for the kernel drafts; the drop regex and CLEAN's `docs/a.md` test case are unchanged | D1 step 0 |

## Review log (final round 6)

Cross-file consistency fixer, round 6. Sources: `reviews/final_r6_{cuda,rl,data,arch,qa,fresh-arch}.md` and their
`*_verdict.md` files. Only the CONFIRMED items assigned to this fixer are listed.

| Item | Severity | Change | Evidence |
|---|---|---|---|
| cuda-1 | blocking | CLEAN16 step 1: the ceiling CLI is a pool user; D7 step 6 already moved it in-process, and this step checks that it imports no pool name | D7 step 6 |
| data-2 | major | CT-33 and CT-35 (and sketch F) state their form after B1.8 (no `--train-views`, one preload); B1.8 updates them in its commit | B1.8 |

## Review log (final round 7)

Cross-file consistency fixer, round 7. Sources: `reviews/final_r7_{rl,arch,qa,fresh-arch}.md` and their
`*_verdict.md` files. Only the CONFIRMED items assigned to this fixer are listed.

| Item | Severity | Change | Evidence |
|---|---|---|---|
| arch-2 | major | CLEAN16 step 1 deletes `SelectStage.procs` and `ProbeStage.procs` with the flag, drops `--procs 1` from the dry run, and adds `procs` to E's U52 dropped set | E6.1; `tests/test_p2_scripts.py:369, 384` |
| fresh-arch-4 | major | CT-33, CT-35 and sketch F key on B1.5 (which now deletes `--train-views`), not B1.8 | B.md B1.5 |

## Review log (final round 8)

Cross-file consistency fixer, round 8. Sources: `reviews/final_r8_{cuda,data,arch,qa}.md` and their `*_verdict.md`
files. Only the CONFIRMED items assigned to this fixer are listed. Each fact was checked in the repo first.

| Item | Severity | Change | Evidence |
|---|---|---|---|
| arch-3 | major | CLEAN16 step 1 becomes a check (one pool rule: D7 step 6 removes callers and `--procs`, E3b archives, CLEAN16 checks); CLEAN.3.4 row and CLEAN2 step 6 follow | `select_ckpt.py:374`; `g_probe.py:139`; `test_p2_scripts.py:369, 384` |

## Review log (final round 10)

Cross-file consistency fixer, round 10. Sources: `reviews/final_r10_arch.md` and `reviews/final_r10_arch_verdict.md`.
Only the CONFIRMED items assigned to this fixer are listed. Each fact was checked in the repo first.

| Item | Severity | Change | Evidence |
|---|---|---|---|
| arch-3 | major | CLEAN8 step 3 no longer writes `MODELS`/`MAX_COMPLETION` literals: it deletes only the `lfm2.5-vl-450m` key from every per-model table and keeps every other key and value (C7's L, front section 3's `qwen3-vl-4b` row). CT-62 drops `set(MODELS) == {"qwen3-vl-2b"}` and pins `MODELS`/`MAX_COMPLETION` to the parent-commit values minus the LFM key (CT-06 literal pattern). CT-05 uses subset checks in Wave 1 and "2B present, LFM absent" after Wave 3. CT-06 notes that the C7 commit updates its `max_completion_length` literal | `src/objdel/train/configs.py:70-77` at be52b1c (`MAX_COMPLETION = {"smolvlm-256m": 1536, "lfm2.5-vl-450m": 1024, "qwen3-vl-2b": 1536}`); 00_front.md section 3 4B step 1; C.md C7 "Files" |

## Review log (final round 11)

| Item | Severity | Change | Evidence |
|---|---|---|---|
| fresh-arch-3 | major | CLEAN11 "Owner of the move": B1.4a creates `rings.py` with `polygonize`; D3 adds `trace` and moves the oracles | B.md B1.4a; D.md D3 |

## Review log (final round 12)

| Item | Severity | Change | Evidence |
|---|---|---|---|
| qa-1 | major | CLEAN11 owner line and the two tables: B1.4a moves the numpy polygonizer, D3 the bit-quads | B.md B1.4a; D.md D3 |
| fresh-arch-1 | blocking | CLEAN1 step 1 keeps B1.6's keep line; `test_bundle_imports_standalone` allows `docs/` paths that the keep rule ships (`re` imported) | B.md B1.6, I50 |

## Review log (final round 13)

Fixer for CLEAN, round 13. Sources: `reviews/final_r13_arch.md`, `reviews/final_r13_fresh-arch.md` and their
`*_verdict.md` files. Both items are the same defect and get one fix. Each fact was checked in the repo first
(read-only git; probe repo `impl/r13clean_probe/`).

| Item | Severity | Change | Evidence |
|---|---|---|---|
| arch-1 | blocking | CLEAN5 new step 0: record the 31 ignored files of the moving trees in `runs/clean/ignored_before.txt`; rewrite `.gitignore:38-49` to the new anchored paths before any `mv` (phase 0 lines follow owner question 3); add CT-89 in the same commit; README note with no backtick path. Step 8 stages `.gitignore` and runs only after the new step 9 RIGA gate (mapped paths exist and are ignored, none stageable, CT-89 passes). New gate CLEAN-G19, new test CT-89 (`test_no_riga_derived_file_is_committable` in `tests/test_p2_breakit_ops.py`: whole-repo `ls-files -co --exclude-standard`, `RIGA_BLOCK` match incl. `t[24]_*.json`, `riga` name rule, temp-repo meta cases). CT-73 grep widened to `csv|jsonl?` plus the RIGA-block match. Step 4 note, 9.4, 9.5, counts (66 tests), 9.7, risk row, owner question 6 (two tracked RIGA tables, allowed by literal by default) | `.gitignore:37-49` (12 anchored lines); `git ls-files -o -i --exclude-standard experiments/phase0 experiments/phase1/lab experiments/phase2/lab` = 11 + 20 + 0 files (2026-11-09); `policy.py:379` drops `docs/`; temp repo: the new 12 lines ignore 31 of 31 mapped paths, `git add -A --dry-run` adds only `.gitignore` and the 2 tracked tables; git-scm.com/docs/gitignore (leading `/` anchors to the `.gitignore` folder) |
| fresh-arch-1 | blocking | Same fix as arch-1. Chosen form: rewrite the anchored lines (files stay beside their records), not a move into `data/`; the verifier named either form plus the mapped-path `check-ignore` gate as the fix. CT-89 adds `riga_block_hits` because a `riga`-name rule alone misses `t_stack/results/t[24]_*.json` (arch verdict note) | `final_r13_arch_verdict.md` note for the fixer; `final_r13_fresh-arch_verdict.md` |

## Review log (final round 14)

| Item | Severity | Change | Evidence |
|---|---|---|---|
| fresh-arch-1 | blocking | CLEAN8 step 4 note and CT-88: E4 archives `hf_sample.py` with `sample_vllm.py` in the E4 step 3 commit (before CLEAN8), not after E-G10's box check; E's U2 checks it | E.md E4 step 1 |

## Review log (final round 15)

| Item | Severity | Change | Evidence |
|---|---|---|---|
| data-2 | major | CLEAN1 step 1 keep-line note names `val_halves.json`; `test_bundle_imports_standalone` `kept_docs` widened to `^docs/reference/(frozen_views|val_halves)\.json$` (same regex as B1.6) | `policy.py:379-380` (drop `^(data|docs)/`, no keep today); `docs/reference/` holds no JSON today [measured: ls]; `policy.py:398-399` (editable install at `/workspace/src`); B.md B1.6 |


# Back matter

## Review log (final round 1)

Sources: `reviews/final_r1_{arch,cuda,data,fresh-arch,qa,rl}.md`, items addressed to 99_back or to "back matter".
Each fact was checked against code, the current phase files, or the saved catalog on 2026-11-08 before a change.

| Review item | Sev | What changed | Check made |
|---|---|---|---|
| fresh-arch 1 | blocking | Decision 1 default is (b'): B1 300 steps first and whole; B2/B3 wait. Decision 2 approves V1 (B1 run then A3), cap $2.20 | front 6.2 V1 table; B.8.1 (B-S1 + B-S2 $4.07 / $3.66); B.8.3 item 1 |
| arch 11 | blocking | Decision 21 (E-Q2) default "no": live CLIs stay in `experiments/phase<n>/`; layer folders only after a CLAUDE.md edit | E.8 owner decision 2 and E8 step 2 (current E.md) |
| cuda 1, qa 1 | blocking | Decision 3: "share R-D2 inside G1" removed. R-D1 = V2 item 4; R-D2 = first item of V3, after D2-D8 merge | front 6.3; D.8 |
| cuda 11 | major | R-D1 placed after the C items in V2, and dropped first if money is short | front 6.3 (R-D1 does not fit at 27.5 s) |
| qa 13 | major | Added D-Q1, D-Q3 (D-Q4 is withdrawn in D.7), D-R6, C8 beta (decision 15), A-G2 acceptance (12), decided test rules (32), the `HARD_STOP` vs $3 question (4). E rows relabelled E-Q1..E-Q6 | D.md:697-701; C.md C8; A.md A.8 and A-G2 row; `policy.py:58` |
| arch 21 | minor | Same relabel; item 19 ADR number fixed; view-bank item marked conditional | E.8 list |
| arch 4, fresh-arch 7, qa 7 | blocking / major | ADR numbers from the front section 4 table: polygonizer exception 0015, infra 0017, numerics/beta 0013 | front section 4 ADR table; last ADR on disk 0011 (`ls docs/architecture/decisions`) |
| rl 4 | major | Decision 15: beta for the B runs. Default follows B.8.3 item 5 (B1 1e-3, B3 `1e-3/r`); RL review's 0.0 shown as the alternative, with the C-G-C consequence | B.8.3 item 5; C.md C8 line 722; B.md B2.4 |
| rl 3, data 4 | blocking / major | Decision 17: 4B gate on val (0.68 / 7 % suggested), never on test | front section 3 gate; R5 val 0.618 / 12.6 % |
| rl 14 | major | New decision 19: "up to 4B" means size class or hard 4.0 B | front section 3 size-limit paragraph (HF API 4.438 B / 4.660 B) |
| qa 18, rl 22 | minor | One medical id `MBZUAI/MedMO-4B` (decision 18) | front section 3; HF API per both reviews |
| cuda 5, data 13 | major | Decision 14 (PNG archive) deferred with the view bank; check every image by `_pixels_sha`. Open questions 4-6 marked DEFERRED | front section 4 "Data format"; `dataset.py:86-87, 108` |
| data 25 | minor | Decision 10: one storage place (model repo `rows/`), pinned `rows_revision` + manifest SHA-256 | B.8.3 item 6; front section 4 "Data versioning" |
| qa 12 | major | Decision 11 states that B run gates keep `sel_mean_of_k` | A-G4 row in A.md; B.md line 982 |
| qa 10 | major | Decision 32 lists the decided rules of A.9.8, B.9.12, C.9.9, E.9.10 for owner sign-off | headings exist in each phase file |
| arch 14, qa 5, cuda 13 | major | Sources point at `docs/implementation2.0/{notes,specs,seeds,kernels}/` (DOC1 step 0), not the scratchpad | DOC.md copy table (line 210-211); B.md header |
| fresh-arch 14 | minor | Open question 16: `dicom-series-to-volume`, `nv-segment-ct` for the 3D direction (recommend only) | saved catalog `nv_catalog.json` lists both |
| cuda 14, qa 15 | major | Money figures use the front section 6 planning table ($0.641/h; 27.5 s base, 23.65 s optimistic) | front section 6 |
| qa 6 | major | Rentals V1-V3; prefixed ids (B-DG14b, C-G-C, D-Q1, E-Q1) | front section 0 id prefix rule |
| fresh-arch 5 | major | No change: this file has no cuCIM line (the cited line is in CLEAN) | grep of 99_back |
| - | minor | Open questions 8, 9, 11, 12, 13 re-pointed to V1/V2 and deferred rows; DOC owner list no longer repeats `HARD_STOP` | front sections 4, 6 |

## Review log (final round 2)

Sources: `reviews/final_r2_{arch,cuda,data,fresh-arch,qa,rl}.md`, items addressed to 99_back or to the back matter,
plus items whose owner decision lands here. Each fact was checked on 2026-11-08 against code, the ledger, or the
current phase files (00_front round-2 text, B, C, D, E, CLEAN, A) before a change. Old decision numbers are kept;
new decisions use sub-numbers.

| Review item | Sev | What changed | Check made |
|---|---|---|---|
| fresh-arch 1 | blocking | Decision 4 rewritten: "keep 17.5" removed; default `HARD_STOP` = $16.90 + $3.00 = $19.90, `ALARM` 10.0, `--owner-approved` on every launch, unit test, DOC U19 reads the constant; a raise in decision 1 raises `HARD_STOP` too | `policy.py:58` (10.0, 17.5), `:73` (`PULL_S` 600), `:324-329`, `:588-592`; ledger 133 rows, sum 16.9025 (`python3 -I` csv sum); front 6.1 step 0 |
| fresh-arch 2, qa 2, data 1, rl 10 | blocking / major | Decision 2 copies the new V1 list (5 items, edge evals not skippable, A3-val skippable, no nsys, no B2.3b, no A3 test half), totals $2.08 / $1.88, cap $2.20 or $2.30, "session spend" defined | front 6.2 table and notes; arithmetic 195.1 min and 175.9 min x $0.641/60 |
| qa 2 (decision 1 part), cuda 3 | blocking / major | Decision 1 re-priced: V1 $2.08 / $1.88 ($2.01 with item 5 skipped); money left $0.74 / $0.87; V2 items 1-3 $0.95 (6 min start-up per process) do not fit $3 at either step time; option (a) $4.31, $5.26 with V2 items 1-3 | front 6.2-6.3 tables; B.8.1 B-S2 row ($1.98 / $1.77); arithmetic |
| arch 10, qa 14 | major | New decision 3b (E-Q7): V4 after E, about 39 min = $0.42, inside $3 unless raised; box halves of E-G10/E-G13 open until V4. Label line now E-Q1..E-Q7 | front 6.3 V4 paragraph; E.8 owner decision 7 (26-40 min, $0.28-0.43) |
| qa 1, cuda 2 | blocking / major | Open question 8: the capture runs in the C0 baseline arm of V2 item 3 (A6 capture); V1 has no nsys window | front 6.2 header; A.md A-G6 row ("open until V2 item 3") |
| qa 13 (back part) | major | Decision 3: D1's Mac parts gate D2, R-D1 gates only D-G7, dropping R-D1 does not block D; R-D1 priced 33 min, $0.35 | front 6.3 text and V2 item 4 |
| arch 8 | major | Decision 10 names both pins, `rows_revision` (40 hex) and `rows_manifest_sha256` (64 hex), checked by `train_rl` and `infra.preflight`; no `data_revision` | B.md and E.md both use `rows_manifest_sha256` (grep: 5 and 15 hits) |
| rl 11 | major | Decision 17 and open questions 1-2: the 4B path starts with SFT step 0 (50-step timing probe); the owner sees SFT + RL cost; test read once after the gate closes | front section 3 step 0 and gate text; `configs.py:90-95, 148` per the front log |
| data 14 | minor | Decision 17: no B session reads test before the 4B gate closes | front section 0 "test", section 3; B task order item 6 |
| qa 23 | minor | Decision 32 adds CLEAN.9.1 with one rule each for Q-T1..Q-T4 (CLEAN must write them) | CLEAN.md lines 980-987 still say "not decided here" |
| data 10 | minor | New decision 7b: target law for the rescaled kind, default (i) | B.md B1 "Trade-off of the rescaled kind" (62.0 %, 14.9 %, 38.0 %, 924 / 1,074 of 7,341) |
| data 5 | major | New decision 8b: a second-seed B1 run is the only measure of seed variance, not in $3, default no | A.md A.7 row "CIs exclude training-seed variance" |
| rl 8 | major | New decision 15b (only if it binds): smaller `eps` or accept a topology share above 0.5; ADR 0014 records it | B.md B3.7 rule 8 (lines 1114-1115) |
| rl 7 | major | Decision 13: option O needs `L - cache >= 1,214`, so it is off at L = 1,200 | C.md C7 step 5 (line 790) |
| fresh-arch 6 | major | New decision 16b: keep TRL + FastGRPO behind the adapter (ADR 0019), verl considered, trigger named | front section 4 row "RL trainer" |
| fresh-arch 7 | minor | Decision 20: live runs attach no `Controller`; with CLEAN6's `--audit` branch, D-G14 adds the audit leg back | `g_probe.py:93` per D.md; D.md D-G14 row |
| fresh-arch 5, rl 9 | major | No change: decisions 9 and 15 already use the front ADR table (0015, 0013); B now uses the same numbers | front section 4 ADR table; B.md lines 21-22, 949 |
| arch 15 | minor | No change: the catalog path `docs/implementation2.0/notes/doc/nv_catalog.json` is already the one DOC1 uses | grep of this file |
| qa (skills row) | - | No change: model ids and sizes agree with section 3 (`MBZUAI/MedMO-4B`, 4.438 B) | front section 3 |

## Review log (final round 3)

Sources: `reviews/final_r3_{arch,cuda,data,fresh-arch,qa,rl}.md`, items addressed to 99_back, plus the back-matter
items that the front round-3 log hands over ("99_back decisions 1, 2, 3b, 13, 17 and the 4B SFT recipe"). Each fact
was checked on 2026-11-08 against code (`policy.py`, `r1_qwen3vl2b.yaml`, `configs.py`, the ledger) or the current
front, A, B and E text before a change. Where B and the front still differ, the front value is used and the owning
section is named.

| Review item | Sev | What changed | Check made |
|---|---|---|---|
| arch 1, qa 1 | blocking | Decision 2 names the launch config `b1_launch.yaml` (the only file `make preflight CONFIG=` stamps); `b1_rows.yaml` and `a3_measure.yaml` are plan files run by `run.sh` | `preflight.load` builds `RunArgs(**raw)` (front log: `preflight.py:49-56`); `cmd_launch` calls `preflight.refusal` [code: policy.py:315]; front 6.1 step 4, 6.2 "Configs" |
| data 1 | blocking | Decision 2 is the front 6.2 list: items 1-5, one R5 process for items 4-5, no `--split test` in V1 (A3's test half is `a3_test.yaml` in V3) | front 6.2 table and "Out of V1" note; front 6.1 step 4d |
| data 2, qa 3, qa 5 | major / blocking | Decision 2: item 5 is no longer skippable; the fallback (R5 `items.jsonl`, hash-only view rule, A3-val as V2 item 1b) applies only if item 5 gives no output | front 6.2 notes (worst case $1.90 + 20.8 min = $2.12 <= $2.20); front log row "data 2, qa 3, qa 5" |
| fresh-arch 5 | major | Decision 2 re-priced: one eval process per checkpoint; 188.3 min = $2.01 (27.5 s), 169.1 min = $1.81 (23.65 s); the $2.30 option dropped; `--cap-hours` 3.40 h at $0.641/h | front 6.2 table (205 s fixed + 0.193 s per item); 2.20 / 0.641 = 3.43, floored to 3.40 |
| fresh-arch 2, qa 6 | major | Decision 2: `BudgetStop` clock = launch time (`train_rl --session-start`); open question 13 re-stated (step time now decides the B1 step count, not a skip; about 29.6 s [estimate]) | front 6.2 code-change list and 6.4 rule 5; `policy.py:323` (`started`, `dph`) |
| qa 7, arch 2 | major | Decision 2 item 1 includes the rows-root fetch (`fetch_hf.sh <run> --rows ...`) and the push of the B1 final adapter to revision `b1-rows`; decision 3 notes V3 step 0 fetches it | front 6.2 launch lines; front 6.3 V3 step 0 |
| qa 4 | blocking | Decision 2 states the tag content (B2.1 logging only on `v1-b1`) | front order rule 2, as copied in the front round-3 log |
| fresh-arch 3, qa 16 (B-S2) | major / minor | Decision 1 re-priced: V1 $2.01 / $1.81; money left $0.74 (69 min) / $0.94 (88 min); C adoption = V2 items 1, 2, 2b, 3 = 102 min = $1.09; all V2 items $1.60; B-S2 $1.94 / $1.73 from B.8.1 (B owns it; the old $1.98 was a third number); option (a) $4.20, $5.29 with the C items; merged run $1.61-2.01 | front 6.3 V2 table and money-left line; B.8.1 B-S2 and merged rows; arithmetic 2.01 + 1.94 + 0.25 = 4.20 |
| fresh-arch 4 | major | New decision 2b: the after-V1 decision table (4 rows) is approved with the V1 item list, before V1 | front section 1.1 and 6.1 step 6; B.4.1 (B's row 4 names only the edge guard and legal; the front row is copied) |
| fresh-arch 4 | major | Decision 17: the 4B gate is evaluated after the last 2B run the budget allows (B1 if B-S2 is not funded), and at once on the "no data effect" row; "single final candidate" points to decision 2b | front section 1.1 last paragraph; front section 3 gate text |
| rl 9 | major | Decision 17 records the 4B SFT recipe: `r1_qwen3vl2b.yaml` with only the model and output folder changed, e2 branch only, priced by the 50-step probe ratio. New finding: the `vision.reference` line names the 2B file, so the 4B needs its own reference (front section 3 step 0 must add it). Open question 2 points to the recipe | `r1_qwen3vl2b.yaml:3-12` (view `lidc-bundle-v3`, seed 1, batch 4 x 8, branches 0.5-3, `vision.reference: vision_ref/qwen3-vl-2b.json`); `configs.py:90-95` (r 16 / alpha 32); `full_rl_2026-11-08.md:3` (RL from e2) |
| arch 16, qa 16 (V4) | minor | Decision 3b drops the E.8 range; 39 min = $0.42 is the one V4 price | front 6.3 V4 paragraph; E.md lines 1387 and 1445 now use $0.42 |
| qa 13 | minor | Decision 13 uses the one limit form of front section 0: `L - soft_punish_cache >= 1,182 + 32` = 1,214, with the outline length 1,181 without the end token and 1,182 with it | front section 0 row "outline length" |
| rl 11 | minor | Decision 16: F is picked by B2.4 rule 0 (B2.0 step-100 relative error >= 0.2; seed 0.531), its arm is V2 item 2b ($0.14). Decision 16b trigger: "the F arm fails twice on a TRL or vLLM bump" | front section 4 row "RL trainer"; front 6.3 item 2b; B.8.3 item 9 |
| data 8 | minor | New decision 8c: the B1.7 test hit (LIDC-IDRI-0659, 4 items), the four train hits and the two val patients on the duplicate-z list; fix or keep, due before A3 step 1; default keep | A.md A3 step 0 (names the hit and sends the choice here); B.md B1.7 table (split, patient, items) |
| qa 20 | minor | Decision 4: the unit test calls `policy.launch_refusal(offer, cap_hours)`, which holds `cmd_launch`'s expression; line ref corrected to `policy.py:321-326` | `policy.py:321-326` (deadline = cap + `PULL_S`, `committed + worst > HARD_STOP`); front 6.1 steps 0 and 4b |
| qa (skills row), data (skills row) | - | No change: model ids, sizes and "test read once" agree with section 3 | front section 3 |
| fresh-arch 1 | blocking | No change here: decision 2 already defines session spend as `spend()` minus the spend at launch; the B-S1 copy is B's (B.8.1) | front 6.2; B.8.1 still lists the old 6-item draft, B must follow |

## Review log (consistency pass, round 4 prep)

| Item | What changed |
|---|---|
| One source | Decisions 1, 2, 3, 3b, 16 and open question 13 point to front 6.2-6.3 instead of quoting minutes and dollars |
| Decision 2 | Approve V1 at $2.40 with `--budget-usd 2.05` (front 6.1 step 0 message); $2.20 would cut the B1 run to 260 steps |
| Decision 4 | `HARD_STOP` stays 17.5 in code until the owner approves; the refusal uses the V1 worst case $2.49 |

## Review log (final round 4)

Cross-file consistency fixer, round 4. Sources: `reviews/final_r4_{cuda,rl,qa,arch,fresh-arch}.md` and their verdicts.

| Item | Severity | Change | Evidence |
|---|---|---|---|
| rl-1, qa-1 | blocking | Decision 16: F is for the runs after Gate B3; its check is the first 15 steps of that run; no V2 item; B1, B3 and the C arms run `bf16-seqmask`. 16b says "F check" | B.md B2.4 "Scope" |
| rl-2, fresh-arch-2 | blocking | Decision 17 holds no recipe copy; it points to front section 3 step 0 | front section 3 |
| cuda-3 | major | Decision 1: 260 steps. Decision 2: $2.40 and 260 steps; $2.20 gives 220; 300 steps need `BudgetStop` >= $2.26 and a recomputed plan. Open question 13 uses the 53-min start and 27.5 s | front 6.2 "Why 260 steps" |

## Review log (final round 5)

Cross-file consistency fixer, round 5. Sources: `reviews/final_r5_{cuda,rl,data,arch,qa,fresh-arch}.md` and their
`*_verdict.md` files. Only the CONFIRMED items assigned to this fixer are listed.

| Item | Severity | Change | Evidence |
|---|---|---|---|
| fresh-arch-1 | blocking | Decision 1 (b'): no rental follows V1; V2, V3 and V4 each wait for a raise; money left is a second reserve. Decision 3: V2 needs a raise; R-D1 is dropped first if the raised amount is too small. Decision 3b: V4 needs a raise. Skills line updated | front 6.3 money-left paragraph and totals |
| arch-1 | blocking | Decision 3b: E4 keeps the old samplers until V4; the trainer eval path is in git history (tag `pre-e4-eval`, E4 step 3) | E4 step 3; `final_eval.py:73-75` |
| cuda-4, qa-3 | major / blocking | Decision 3b names `e_v4.yaml`, which starts with the CUDA build | front 6.3 V4; E.8 |
| data-1, qa-2 | blocking | Decision 2b points to front 1.1 (the only copy); its row 4 = offset < 0.15 guard, edge guard, legal share or frozen-view hash; the A-G3 correlation is not a Gate B1 rule | front 1.1; B.6.2 |
| rl-4 | major | Decision 15: the alternative is now "beta 0.0 for B3 only"; B1 keeps 1e-3 in every case (C8, U-CO5, B.8.3 item 5) | C8; U-CO5; B.8.3 item 5 |

## Review log (final round 6)

Cross-file consistency fixer, round 6. Sources: `reviews/final_r6_{cuda,rl,data,arch,qa,fresh-arch}.md` and their
`*_verdict.md` files. Only the CONFIRMED items assigned to this fixer are listed.

| Item | Severity | Change | Evidence |
|---|---|---|---|
| fresh-arch-1 | blocking | Decision 1 (b') points to the front "End state under (b')" | front section 1 |
| cuda-4, arch-4 | major | Decision 3b: V4 time described as launch to first `run` plus about 16 min; front 6.3 holds the price, launch config and lines | front 6.3 |
| cuda-5 | major | Decision 3: R-D1 runs only before D2 merges (D1 step 4) | D1 step 4 |

## Review log (final round 7)

Cross-file consistency fixer, round 7. Sources: `reviews/final_r7_{rl,arch,qa,fresh-arch}.md` and their
`*_verdict.md` files. Only the CONFIRMED items assigned to this fixer are listed.

| Item | Severity | Change | Evidence |
|---|---|---|---|
| fresh-arch-4 | major | Decision 2b rows 3 and 4 follow front 1.1 (B1.8 runs on row 3; after the owner review on row 4); `--train-views` is gone before V1 (B1.5) | front 1.1; B.md B1.5 |

## Review log (final round 8)

Cross-file consistency fixer, round 8. Sources: `reviews/final_r8_{cuda,data,arch,qa}.md` and their `*_verdict.md`
files. Only the CONFIRMED items assigned to this fixer are listed. Each fact was checked in the repo first.

| Item | Severity | Change | Evidence |
|---|---|---|---|
| data-1 | blocking | Decision 12: val reference 0.723 (all items) / 0.732 (test); found-stratum line 0.734 / 0.740 | A2 step 4 |
| arch-1 | blocking | Decision 3b names `specs/e_v4.yaml` | E6 step 8b; front 6.3 |

## Review log (final round 9)

Cross-file consistency fixer, round 9. Sources: `reviews/final_r9_{cuda,rl,data,arch,qa,fresh-arch}.md` and their
`*_verdict.md` files. Only the CONFIRMED items assigned to this fixer are listed. Each fact was checked in the repo first.

| Item | Severity | Change | Evidence |
|---|---|---|---|
| rl-2 | major | RL settings 16: an F-trained run is evaluated unmerged with `final_eval --adapter <out>/final_adapter` (B2.4 eval rule) | `train_rl.py:58-59`; `final_eval.py:48` |

## Review log (final round 10)

Cross-file consistency fixer, round 10. Sources: `reviews/final_r10_{cuda,data,qa,fresh-arch}.md` and their
`*_verdict.md` files. Only the CONFIRMED items assigned to this fixer are listed. Each fact was checked in the repo first.

| Item | Severity | Change | Evidence |
|---|---|---|---|
| qa-4, fresh-arch-2 | major | Decision 1 (a): B-S2 = V3 items 0 + 4; other V3 items only by name; the amount follows front 6.3 "V3 amount rule". Decision 16 renamed "Optional V3 arms" (not part of B-S2) | front 6.3 |

## Review log (final round 11)

| Item | Severity | Change | Evidence |
|---|---|---|---|
| fresh-arch-3 | major | Decision 9 removed (number kept): no polygonizer exception; B1.4a seam | K1/results/cross_numpy_m3.txt |
| fresh-arch-1 | major | Decision 3b: a B run after E pairs within V3 (B3.9 step 3); V4 still closes E-G10 | E.md E4 Gate |

## Review log (final round 12)

| Item | Severity | Change | Evidence |
|---|---|---|---|
| qa-2 | major | Decision 3b: B3.9 step 3 -> step 2 (the B1 re-eval) | B.md B3.9 |

## Review log (final round 14)

| Item | Severity | Change | Evidence |
|---|---|---|---|
| fresh-arch-1 | blocking | Decision 3b: the old samplers are archived in the E4 step 3 commit either way; tag `pre-e4-eval` keeps them | E.md E4 step 1 |

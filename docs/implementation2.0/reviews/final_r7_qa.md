# Final review round 7: QA lead and technical writer (ASD-STE100)

Target: `docs/implementation2.0/Implementation.md` (round 7 build: 1,577,738 bytes), built from `docs/implementation2.0/*.md` by
`assemble.py`. I read the round-6 fix rows of every phase in `review-log.md`, then the text that each fix changed:
B2.2 (F plumbing, step 4b, step 5), B2.4 (rules 0-5, "Missing inputs"), B.9.0 (fixture `tiny_loop_vllm_fake`),
B.9.1 U51, B.9.5, B.9.8 (decision sketch), B1.3-B1.5, front sections 1, 1.1, 6, 6.1, 6.2, 6.3, 6.4, E4 gate,
E6.1 (`EvalStage`, "Eval argv"), E.8 (V4), E.9 (conventions, U21, U49, U51, GP1-GP5), A1 steps 1-5, A.6, A.9.0,
A.9.6, D8 (`D_GPU_TESTS`, R-D2 payload, `d_gpu_test.yaml`) and CLEAN16. Items fixed in rounds 1-6 are not repeated.

No REPO file was edited by hand. No git write. No rental. Python that read `data/` ran with `-I`. One note: I ran
`python3 -I assemble.py` once by mistake (inside a combined shell line). The sources had not changed since the
23:07 build, so the rebuild wrote the same bytes: `Implementation.md` 1,577,738 B and `review-log.md` 300,925 B
before and after [measured: `ls -l`].

## Skills used

| Skill | What it changed in this review |
|---|---|
| nvidia-skills:nvidia-skill-finder | The round-6 text adds no new NVIDIA product claim, so I did no catalog lookup. It pointed me at the box toolchain path: who runs which pytest under the CUDA build (`run_box.sh --gpu-check`, `D_GPU_TESTS`). That led to item 2. |
| engineering:code-review | I read the code behind each claim I checked: `policy.py:420-500` (`chain`, `cmd_run` and its module regex), `policy.py:505-536` (`cmd_push_hf`, `cmd_fetch_hf`), `policy.py:374-417` (`code_files`, `cmd_sync_up`), `env_check.sh:1-30`, `final_eval.py:44-56`, `fast.py:545, 568`, `tests/test_p2_scripts.py:391-393`, `huggingface_hub` 2.2.0 `hf_api.py:5770-5790` and `_upload_pipeline.py:476, 650, 684-711`. Items 1, 2, 3. |
| engineering:debug | I traced the B2.4 restart path step by step: the run fails its F check, `train_rl` raises, and then `restart_preset` runs. The seam has no input for the start preset or for the comparator result (item 1). I also traced how an A box test reaches the box: `run.sh` accepts only a module or a yaml, so a bare pytest command cannot be issued (item 2). |
| engineering:system-design | I put B2.4 rule 0, B2.2 F steps 4, 4b and 5, U51 and the B.9.8 sketch side by side. The F check has six parts, but the restart seam reads three (item 1). |
| engineering:architecture | Checked each rental's launch config against `cmd_launch` and `cmd_push_hf` (V1-V4): all four now have one (round-6 fix holds). No new item. |
| engineering:documentation | Checked "one rule, one place" for the V4 work minutes (E.8 says 9 min; the front sums 10 min of the same steps, item 8). Checked naming: `ci95_patient` used for a 90 % interval (item 6). |
| huggingface-skills:huggingface-best | No model choice changed. I used it only to check the HF upload claim in B1.4 step 5: `upload_folder` commits in parts and returns the last `CommitInfo` [code: `.venv huggingface_hub/_upload_pipeline.py:650, 684-711`]. The claim holds. No item. |
| engineering:testing-strategy | Matched each test to a runner, an input file and a rental step. A2/A3's G02b and G03 have no runner (item 2). GP5 has no `group_id` source (item 3). B's G06 has no file (item 7). |
| mattpocock-skills:tdd | Looked for tests that cannot fail or that encode the wrong rule. The `restart_preset` cases pass for a function that ignores the comparator. Such a function then restarts into the preset that just failed (item 1). |
| Web | pydantic-settings docs: "Aliases of length one are converted into short options", and lists accept repeated flags and lazy comma lists [web: https://pydantic.dev/docs/validation/latest/concepts/pydantic_settings/]. This confirms the E6.1 `-k` and comma rules. My own probe agrees (below). |

Measured or read during this review (2026-11-08):
- pydantic-settings 2.16.2 probe (`(planning probe, not kept: r7qa)/p1.py`, the round-6 probe venv): interleaved `--view a --split val --out
  o1 --view b --split val --out o2 -k 8` parses to `view=('a','b')`, `split=('val','val')`, `out=('o1','o2')`,
  `k=8`. One `--out` after two views gives `out=('o',)`. `--k 8` exits 2 ("unrecognized arguments: --k 8"). The
  E6.1 "Eval argv" text is correct.
- TOST noise floor (`(planning probe, not kept: r7qa)/noise.py`, `python -I`): split the K = 16 samples of
  `runs/r1/r1-qwen3vl2b/select/e2/shifted-k16/items.jsonl` into two K = 8 halves. Join `group_id` from
  `data/store/lidc/val.parquet`. Result: n = 1,084, 78 patients, mean difference 0.0021, item SD 0.104. The 90 %
  patient-bootstrap CI is [-0.0045, 0.0083] (half-width 0.0064, B = 4,000, seed 0). So the E4 TOST bound of ±0.01 can
  pass with pure sampling noise on an SFT checkpoint, but the margin is small (item 4).
- `runs/r5/eval/eval_val/items.jsonl`: 1,084 rows. The keys are `item_id, legal, iou_tc, f5, precision_vote,
  recall_vote, dice_vote`, and there is no `group_id` [measured: `head`, `wc -l`].
- `cmd_run` accepts only `experiments.phase2.<script>` (or R0's `xarch_hashes`) or a `configs/*.yaml` path; any
  other first argument dies [code: experiments/phase2/vast/policy.py, `cmd_run`, the `re.fullmatch` branch and the
  `die("module must be ...")` line].

Severity follows the task rules. Count: **0 blocking, 2 major, 7 minor**.

---

## Major

### 1. [major] [B] The `restart_preset` seam (U51, B.9.8) has no input for the start preset or for the comparator result, so it can restart into the preset that just failed

Problem. B2.4 says what a failed F check restarts under. The U51 text and the B.9.8 sketch encode a function that
cannot follow that rule.
- B "#### B2.4", rule 0: "If the F check fails, `train_rl` raises after writing `f_check.json`, and the run restarts
  under `fp16-seqmask` (if rule 1 passed) or else the rules 2-3 preset (...). A failed check in a run that started
  under `fp16-seqmask` (error < 0.2, rule 1) restarts under the rules 2-3 preset."
- B "#### B2.4", rule 0, the F check: "NaN/inf = 0; the LoRA tensors equal the trainer's under vLLM's layout (B2.2 F
  step 4); the output-level check against the comparator that the run measures before step 1 (B2.2 F step 4); and the
  timing rule." B "#### B2.2", F step 5 adds "the LoRA load count per RL step is 1". Step 4b.5 adds the base-hash
  restore check ("A mismatch fails the F check before step 1"). Under `fp16-seqmask` the checks are "the engine dtype
  assertion (fp16 step 3) and NaN/inf = 0". Step 5 also says `f_check.json` holds "each check's value and pass flag,
  `passed`".
- B "### B.9.1", U51: "`restart_preset` (reads `f_check.json`): F passed (NaN 0, LoRA tensors equal under the packed
  map, in-process median LoRA / base generate-time ratio <= 1.15) -> lora-seqmask kept; ... F failed and fp16 passed
  rule 1 -> fp16-seqmask".
- B "### B.9.8", sketch: `F = lambda ok, ratio: dict(nan=0 if ok else 1, lora_equal=ok, lora_over_base=ratio)` and
  `restart_preset(CF(), f_check=f, fp16=fp, b1=B1)`. There is no start-preset argument, and there are no fields for
  the comparator, the load count or the base hash.

Result.
(a) A run that fails only the comparator check (the check that "F exists to fix", B2.2 F step 4) has `nan=0`, equal
LoRA tensors and a ratio <= 1.15. A `restart_preset` built to U51 and the sketch returns `lora-seqmask`. That is the
preset that just failed.
(b) A run that started under `fp16-seqmask` (error < 0.2, rule 1) and fails its F check also gets the wrong answer.
`restart_preset(F(False, ...), FP(True, 0.010, ...))` returns `fp16-seqmask`, because the function cannot see the
start preset. The rule text says "rules 2-3 preset".
In both cases a paid restart (about 6 min start-up plus 15 steps, B.8.1) runs the failed preset again. This is a
contradiction between two parts that both state the rule. I rate it major, not blocking: it affects only the first
run after Gate B3. That run is outside the funded plan and lives on the `after-raise` branch (front "End state under
(b')").

Exact fix.
- B2.4 rule 0 and U51: `restart_preset(cf, start, f_check, fp16, b1)`, where `start` is the preset the run started
  under (from `f_check.json["preset"]`). It reads `f_check["passed"]` and never recomputes the pass from some of the
  fields. Rule: `passed` -> keep `start`; failed and `start == "lora-seqmask"` and fp16 passed rule 1 ->
  `fp16-seqmask`; failed otherwise (including `start == "fp16-seqmask"`) -> the rules 2-3 preset (`bf16-tis-opsm`,
  or `bf16-seqmask` by the rule 3 veto).
- B.9.8 sketch: make `F(...)` a full `f_check.json` record: `preset`, `nan`, `lora_equal`, `comparator_ratio` (mean
  |f_lora - b| / mean |f_merged - b|), `load_count`, `base_hash_ok`, `lora_over_base`, `passed`. Add the cases:
  comparator 0.51 with every other check passing -> not kept (`bf16-tis-opsm` with no fp16 arm);
  `start="fp16-seqmask"` failed with an fp16 arm that passed -> `bf16-tis-opsm`; `base_hash_ok=False` -> not kept.
- U76 (the Mac twin of the F check) asserts that `passed` is the AND of every check that step 5 lists, so that
  `restart_preset` can trust the one flag.

### 2. [major] [A, D, 00_front] A2/A3's box tests G02b and G03 have no plan step, and the stated command cannot be issued through `run.sh`

Problem.
- A "### A.9.6", G02b: "Run place: box, V3 item 1 (R-D2), as its own command after D's `make gpu-check` (about 2 min
  [estimate]; the front prices it inside R-D2, or it waits for the next box). Command: `/opt/venv/bin/python -I -m
  pytest -q tests/test_a_gpu.py::test_g02b_builds`." G03: "Run place and command: as G02b, node id
  `tests/test_a_gpu.py::test_g03_device_bytes`" (pixels stay on the device, the owner's GPU rule).
- D "### D8": `d_gpu_test.yaml` "calls `OBJDEL_DATA_ROOT={data} PY=/opt/venv/bin/python bash
  experiments/cuda_checker/run_box.sh --gpu-check` as one step", and `--gpu-check` runs "`pytest -q $D_GPU_TESTS -m
  gpu`". The `D_GPU_TESTS` list does not include `tests/test_a_gpu.py`.
- 00_front "### 6.3", V3 item 1: "R-D2: `make gpu-check` (build, `pytest -m gpu` on D's list, 4 sanitizer tools ...)".
  The V3 launch lines have one call for item 1: `run.sh <run> experiments/phase2/configs/d_gpu_test.yaml`. So "the
  front prices it inside R-D2" is not true.
- Code: `cmd_run` takes `experiments.phase2.<script>` or a `configs/*.yaml` path and dies on anything else
  (`die("module must be experiments.phase2.<script> or experiments/phase2/configs/<name>.yaml")`)
  [code: experiments/phase2/vast/policy.py, `cmd_run`]. A bare `python -m pytest ...` cannot be sent with `run.sh`.

Result. An engineer on the paid V3 box cannot run G02b or G03 by any documented step. They must invent a plan step
and its minutes. Otherwise the only box check that A2/A3 eval pixels stay on the device (G03) and that both builds
agree on the ceiling items (G02b) never runs. You also cannot add the whole file to `D_GPU_TESTS`.
`tests/test_a_gpu.py` holds G07 and G08 under the same `gpu` marker, and G08 "Skipped ... while the boundary IoU is
CPU-only". Step 3 of `--gpu-check` fails when `skipped > 0`.

Exact fix.
- D8, `d_gpu_test.yaml`: add a second step after the `--gpu-check` step. In it, `/opt/venv/bin/python -I -m pytest -q
  -p no:cacheprovider tests/test_a_gpu.py::test_g02b_builds tests/test_a_gpu.py::test_g03_device_bytes` runs with
  `OBJDEL_DATA_ROOT={data}` (the CUDA build from step 1 is still in `/opt/venv`). Add `prepared/lidc/val.jsonl` (G02b's
  50 val items) to the R-D2 payload note; it is already in the payload.
- A.9.6 G02b/G03: "Run place: V3 item 1, `d_gpu_test.yaml` step 2 (node ids)". Delete "as its own command" and "or it
  waits for the next box".
- 00_front 6.3 V3 item 1: add "+ A's G02b and G03 (2 min)" to the item text, and change 25-44 to 27-46 min ($0.29-0.51).
- D's S08 / `--gpu-check --dry`: also run `pytest --collect-only` on the two node ids, so that a renamed test fails
  on the Mac.

---

## Minor

### 3. [minor] [E] GP5 does not say where `group_id` comes from

E "### E.9.5", GP5: "the 90 % patient-level (`group_id`) paired bootstrap CI". The trainer-arm file
`runs/r5/eval/eval_val/items.jsonl` has no `group_id` key [measured: `head`]. A1 step 2 reads groups from
`dataset.records(a.view, "lidc", a.split)`. The A3 `answers.jsonl` carries `group_id` (A.md A3 step 4,
`{item_id, group_id, W, H, samples[K], ...}`).
Fix: GP5 text: "groups = `{item_id: group_id}` from the pulled sampler arm's `runs/<v4 run>/eval_v4_sampler/answers.jsonl`.
Every trainer-arm `item_id` must be in it, or the test fails."

### 4. [minor] [E] The TOST bound ±0.01 has a small margin over pure sampling noise; print the floor beside the result

Measured above: two independent K = 8 halves of the same SFT checkpoint give a 90 % patient CI of [-0.0045, 0.0083]
(half-width 0.0064). A sampler path with no real difference has a fair chance of a CI edge near ±0.01.
Fix: E4 gate, "Statistical" bullet: GP5 also prints this noise floor. Compute it from A3's K = 8 val answers of R5
(two K = 4 halves, A3 step 7), so the owner reads "not shown" against a known floor. No gate change.

### 5. [minor] [E] U51 (`e_v4.yaml`) does not pin the seed and split of the sampler arm

E "### E.9.1", U51 checks "one `final_eval` on `lidc-val-shifted` with K 8 (`-k 8`) ...". The TOST needs "the same
checkpoint, view, K and seed" (E4 gate). E.8 step 3 writes `--split val ... --seed 11`.
Fix: add to U51: "argv holds `--split val` and `--seed 11` (R5's seed, `runs/r5/eval/eval_val/final_eval.json` `seed`
11)".

### 6. [minor] [E] `ci95_patient` names a 95 % key but the TOST reads a 90 % interval

E "### E4", Decision interval: "resampled by patient (key `group_id`, `ci95_patient` at the 90 % level)". A1 step 3
defines `ci95_patient` as the 95 % CI in `final_eval.json`. A reader may take the stored 95 % value. That makes the
test stricter, not wrong.
Fix: "the 90 % patient-bootstrap interval (`signal.ci` with patient draws, `level=0.90`); not the stored
`ci95_patient`".

### 7. [minor] [B] G06 has no test file, and the B.9.5 runner row names the wrong file for G01

B "### B.9.5", runner table: "G01, G06 | R-D2 (V3 item 1), with D's list | `make gpu-check` (D8 file list:
`tests/test_gpu_scoring.py`, `tests/test_cuda_checker_gpu.py`)". The B.9.8 sketch puts G01 in
`tests/test_scorer_sums.py`. D8 adds that file to `D_GPU_TESTS` "once B3.1 lands". No file is named for G06.
Fix: the row reads "`tests/test_scorer_sums.py` (G01, G06; in `D_GPU_TESTS` after B3.1)". The G06 row names the same
file.

### 8. [minor] [E, 00_front] V4 work minutes: E.8 says 9 min; the front sums 10 min for the same steps

E "## E.8": "E's work inside it is about 9 min of checks [estimate]; the front adds launch to first `run`, the 3-min
build and the eval-process imports". Front 6.3 V4: "`env_check` 1 + box tests 2 + sampler eval 7" = 10 min. E.8's
last bullet also says "about 9 min".
Fix: E.8 says "about 10 min (front 6.3: `env_check` 1, box tests 2, sampler eval 7)" in both places, or points to
front 6.3 without a number.

### 9. [minor] [E] GP3 and E-G13 say "before training starts", but V4 trains nothing

E "### E.9.5", GP3: "`IMPORT_CHECK` exits 0 before training starts". The E-G13 row says the same. V4 has no `rl`
stage (U51).
Fix: "before the box tests of step 2". The same change goes in the E-G13 row.

---

Checked and found correct this round (no item): the V1-V4 money arithmetic (V4: 38.2 + 3 + 1 + 2 + 7 + 2.6 = 53.8
min -> $0.71; cap 0.95 h; worst $0.86. V2 92 min -> $1.14. B-S2 177.6 / 160.9 min -> $2.11 / $1.92. (a) $4.65). The
`cmd_push_hf` revision equals `output_dir` for `b1_launch.yaml` (`b1-rows`) and `b3_launch.yaml` (`b3-v5`). The
`env_check.sh {rd}/ckpt {rd}/env_check.json` form appears in B1.5 step 1 and E.8 step 1. The `tiny_loop_vllm_fake`
premise is right (`fast.py:568`; `GRPO_SMOKE` `use_vllm: False`). The B1.4 `CommitInfo.oid` rule holds. The E6.1
tuple and `-k` rules hold (probe above). The B2.4 start choice and its sketch cases agree.

# Final review round 4: fresh outside architect

Reviewer: a fresh outside architect with no earlier opinion. Doc: `docs/implementation2.0/Implementation.md` (round 4,
14,787 lines), built from `docs/implementation2.0/*.md`. I read the front matter (sections 0-8) in full. I also read B1.5,
B.6.2, B.9.5, D3, A3 (plan files and header), the 99_back owner decisions and the open questions. I checked each
claim below against the repo code, the R5 run files and PyPI. I did not repeat items that `review-log.md` already
lists as fixed.

## Skills used

| Skill | What it changed in this review |
|---|---|
| `nvidia-skills:nvidia-skill-finder` | Used for the CUDA toolchain check. The front says the catalog was already read on 2026-11-08, so I did no new catalog lookup and recommend no install. It made me check the toolchain pins on PyPI: `cuda-toolkit` 13.0.3 and `nvidia-cuda-sanitizer-api` 13.4.92 both exist (item 10). |
| `engineering:code-review` | Every blocking and major item quotes a code line: `policy.py:450, 469-480, 511-536`, `configs.py:179-181`, `rewards.py:132-135`, `train_rl.py:58-60`, `bundle.py:53-62`. I grepped `sum_x` in `cuda/` and `src/` (no match). |
| `engineering:debug` | Traced what each V1 plan step does on the box. `policy.chain` joins the steps with `&&`, so one failing smoke test stops the B1 run (item 1). Traced the 4B SFT recipe through `chain` and `vision_refs` (item 2). Read the R5 log for the time to save the final model (47.62 s, item 5) and the R5 `EXIT` code (137, item 8). |
| `engineering:system-design` | Judged the system as a whole: data path, rental order, money guards, eval path. The design is sound. My changes are about seams between parts: plan file to test file, recipe to runner, bundle to gate. |
| `engineering:architecture` | Checked the decision tables (sections 3 and 4) against code facts. The TRL row has a missing alternative: TRL 1.17.1 ships `trl/experimental/async_grpo`, found in the local `.venv` (item 11). |
| `engineering:documentation` | Looked for places where two parts each claim to be the rule: front 3 vs 99_back 17 (item 2), A3 header vs front 6.1 step 4c (item 4), B.9.5 G01 row vs G08 row (item 1). |
| `huggingface-skills:huggingface-best` | No Hub leaderboard covers polygon text for CT. The model choice (keep Qwen3-VL-2B, 4B only behind a val gate) is right. The review effort went into the 4B recipe, which cannot run as written (item 2). |
| `system-design-skills:system-design` | Stress-tested the V1 time budget: what happens when a step fails or runs long. That is where items 1, 5 and 8 come from. |
| `engineering:tech-debt` | Checked that CLEAN and D leave no second copy of anything. D3's hash check misses the newest data artifacts (item 7). |

Web: PyPI JSON for `nvidia-cuda-sanitizer-api` and `cuda-toolkit` (fetched 2026-11-08). The NVIDIA CUDA
minor-version compatibility page: it says nothing about tools such as compute-sanitizer.

## Verdict

The architecture is the right one for this budget and this stage. It puts data first, runs one paid experiment
with a decision table approved in advance, uses one CUDA engine source for the CPU and GPU builds, uses deep seams
(`geom.polygonize`, `policy.sampler`, `rl.trl_adapter`), and keeps TRL behind one adapter. A top lab would not
restructure it. A top lab would do five things differently. None of them needs a new system:

1. Make the first paid step fail only for real reasons. Today the V1 smoke test checks code that the V1 tag does
   not hold (item 1).
2. Make each conditional path runnable as written. The 4B recipe cannot run on the existing runner (item 2).
3. Run the V1 work as one chain on the box, with the measured save time counted. Today there are two `run.sh`
   calls with a Mac poll between them, and the slack is 1-3 minutes (item 5).
4. Take the eval path out of the trainer earlier. Every eval process pays about 205 s to build a trainer and vLLM
   (front 6 inputs). E4's `VllmSampler` fixes this, but only after E and V4 (item 12).
5. Record the overlap option (TRL's experimental async GRPO) as a named revisit trigger. The GPU is idle about 40 %
   of each step (item 11).

## Items

### 1. BLOCKING. B (B1.5, B.9.5), 00_front (6.2 plan files, order rule 2): the V1 GPU smoke test G08 includes G01, which tests B3.1 code that the V1 tag does not hold. The `&&` chain then stops before the B1 run.

Doc lines:
- B.md, "B1.5 Train on the ordered rows", step 2 of `b1_rows.yaml`: "The GPU smoke `pytest -m gpu
  tests/test_train_fast.py -k "G02 or G08"` (< 5 min; its minutes are in front item 1)."
- B.md, "B.9.5 GPU tests", row G08: "Preflight smoke gate (< 5 min): G01 parity on 1,000 windows, G02 no-sync run
  ... | B2.1, B2.2, B3.1".
- B.md, "B.9.5 GPU tests", row G01: "CPU build and CUDA build give bit-equal results on every field ... including
  `sum_x`, `sum_y` and all old fields; old fields equal the goldens from the previous build. | B3.1". The command
  table in the same section puts G01 in "R-D2 (V3 item 1)".
- 00_front, "1. Summary and order of work", order rule 2: the tag `v1-b1` "holds no C, D, B2.2-B2.4, B3 or kernel
  code."

Code and measurement:
- `sum_x` and `sum_y` do not exist today: `grep -n sum_x cuda/gpucheck.cpp src/objdel -r` finds nothing. B3.1 adds
  them (`kPart` 5 -> 7, new `ObjOut` fields).
- `policy.chain` returns the `steps:` list as is [code: experiments/phase2/vast/policy.py:434-435]. `cmd_run` joins
  the steps with `&&` [code: policy.py:490, `jobs = " && ".join(...)`]. Only a step that starts with `~` may
  fail and let the chain go on.

What goes wrong: G08 runs on the V1 tree. Either its G01 part fails, because there is no `sum_x` and there is no
"previous build" golden for the B3.1 fields, or it must be skipped, which the doc does not say. A failing step 2
stops the chain before step 3 (`train_rl`) and step 4 (the B1 val evals). V1 then pays item 1 ($0.44) and the R5
evals, but makes no B1 run and no Gate B1 input. The doc also contradicts itself: G01's own row says G01 runs in
V3, not V1. The "Every new GPU path runs before money is spent" goal stays valid. The V1 tree adds no new scorer
path, so V1 only needs today's scorer parity.

Exact fix (B.md B.9.5, B1.5 step 2; front 6.2 item 1 text):
- Split G08 into two tests. G08a (V1): G02 no-sync run, plus parity of today's scorer fields, CPU build against
  CUDA build, on 1,000 rows-root windows (the existing `tests/test_gpu_scoring.py` check). It has no `sum_x`, no
  `sum_y` and no B3.1 goldens. G08b (the first rental after B3.1 merges, V3): G01 sums, plus the fp16 `vllm_dtype`
  check when that arm runs.
- `b1_rows.yaml` step 2 becomes `[-m, pytest, -m, gpu, tests/test_train_fast.py, -k, "G02 or G08a"]`. Write it
  with `-m pytest`: `cmd_run` turns any step that does not start with `bash` into `/opt/venv/bin/python -I <argv>`
  [code: policy.py:484-487]. A bare `pytest` first token would run as a script path and fail.
- Add a Mac test to `tests/test_repo_hygiene.py`: every test id named in `b1_rows.yaml` exists on the `v1-b1` tree
  and imports no B3, C or D symbol.

### 2. BLOCKING. 00_front (3, "If the gate fires", step 0), 99_back (decision 17): the 4B SFT recipe cannot run as written. With `branches: [2.0]`, the runner also runs select and probe on val, and the 2B G-vision reference makes `run.sh` refuse the chain. The two parts disagree about what changes.

Doc lines:
- 00_front, "3. Model decision", step 0: "Use the 2B recipe unchanged: `experiments/phase2/configs/r1_qwen3vl2b.yaml`
  with only `model: qwen3-vl-4b` and `output_dir` changed. ... Train only the e2 branch (`branches: [2.0]`), so no
  branch selection runs and val stays unread."
- 99_back, "Owner decisions needed" > "Model" > 17: "The line `vision.reference: vision_ref/qwen3-vl-2b.json` is
  model-specific ... so this line changes too. Front section 3 step 0 lists only `model` and `output_dir` and must
  add it." The front was not changed. Both parts claim to be the recipe.

Code:
- `stages = c.get("stages", ["train", "select", "probe"]) if c.get("branches") or c.get("branches_from") else
  ["train"]` [code: experiments/phase2/vast/policy.py:450]. With `branches: [2.0]` and no `stages:` key, the chain
  adds `select_ckpt` on `lidc-val-shifted`, `-rescaled` and `-rotated`, then `g_probe` on `lidc-val-shifted`
  [code: policy.py:441-448]. It also adds the CUDA build and the GPU scoring test first [code: policy.py:451-452].
  So val is read, and GPU minutes are spent on stages that the recipe says do not run.
- `cmd_run` refuses the chain on the box when the G-vision reference file is missing: "missing G-vision reference
  ... run r0.yaml first" [code: policy.py:471-475 with `vision_refs`, :456-461]. `train_sft` requires G-vision;
  `vision: off` is only for discarded timing runs [code: experiments/phase2/train_sft.py:18-20, 103-107]. The 2B
  reference is the wrong file for a 4B projector, and no step writes a 4B reference.
- `r1_qwen3vl2b.yaml:7-9`: `branches: [0.5, 1.0, 2.0, 3.0]`, `vision: {prepared: prepared, reference:
  vision_ref/qwen3-vl-2b.json}`.

What goes wrong: an engineer who follows front 3 exactly, after the 4B gate fires, launches a rental. On the box,
`run.sh` refuses the chain (no 4B reference). If the engineer copies the 2B reference to get past the guard,
G-vision compares 4B features against a 2B baseline, and select and probe then read val. Either way GPU money is
spent and the stated rule ("val stays unread") is broken.

Exact fix (00_front section 3 step 0; 99_back decision 17 points to it and holds no copy):
- New file `experiments/phase2/configs/sft_qwen3vl4b.yaml` = `r1_qwen3vl2b.yaml` with these changes:
  `model: qwen3-vl-4b`, `output_dir: sft-qwen3vl4b`, `branches: [2.0]`, `stages: [train]`, and `vision: {prepared:
  prepared, reference: vision_ref/qwen3-vl-4b.json}`.
- Plan step 0a, in the same rental and before the chain: `run.sh <run> experiments.phase2.train_sft --config
  experiments/phase2/configs/sft_qwen3vl4b.yaml --vision_reference_only true`. It writes the 4B reference on CPU in
  float32, as R0 did. Add its minutes to the 4B price shown to the owner.
- A Mac test: `policy.chain(sft_qwen3vl4b.yaml)` returns exactly one stage, the `train_sft` stage.

### 3. MAJOR. B (B1.5 "Held-out curve"), CLEAN (CLEAN6 step 4a): the held-out eval sets `num_generations_eval=4`, which the config guard and the reward refuse. The tag's pre-check then cannot pass.

Doc lines:
- B.md, "B1.5 Train on the ordered rows", held-out bullet: "`train_rl --heldout <view> --heldout-items 128
  --heldout-every 50` sets `GRPOConfig(eval_strategy="steps", eval_steps=50, num_generations_eval=4)` ... Cost: 6
  evals x 512 answers ... Before the tag `v1-b1`, the tiny TRL loop test (`tests/test_p2_train_rl.py`) runs one
  held-out eval with FastGRPO's patches on."
- B.md, "B.9.1 Unit tests", U38: "`num_generations` or `num_generations_eval` not equal to `reward.group` raises,
  8 and 8 pass."
- CLEAN.md, CLEAN6 step 4a repeats `num_generations_eval=4`.

Code:
- `grpo_config` raises when `num_generations_eval` is not in (None, 8) [code: src/objdel/train/configs.py:179-181].
- `TRLReward` raises unless each reward call holds "whole groups of 8 completions per prompt" [code:
  src/objdel/train/rewards.py:132-135].

What goes wrong: once the held-out flag is built, even with its default off, the pre-tag test that runs one
held-out eval fails at config time. The engineer must then choose between 4 and 8. That choice changes the guard,
the reward grouping and the cost.

Exact fix: use `num_generations_eval=8` (= `reward.group`; no guard change). Keep 512 answers per eval with 64
items: `--heldout-items 64` x 8. Cost: 6 x 512 answers, about 2 min, as stated. If 128 items are wanted, the cost
is 1,024 answers per eval, about 4 min, $0.05. Change CLEAN6 step 4a to match.

### 4. MAJOR. A (A3 plan-file header and bundle note), 00_front (6.1 step 4c, 6.4 rule 7): the A3 command list ships the test view to the V1 box, while front 6.1 step 4c requires that no test folder is in the V1 bundle.

Doc lines:
- A.md, "A3. Save every answer ...", header comments, in order: `build_eval_views --split test --kinds shifted
  --out data/eval_views`, then `final_eval --hash-view data/eval_views/lidc-test-shifted`, then `vast.bundle lidc
  --view lidc-bundle-v3`. The same section: "`lidc-test-shifted` reaches the box when it exists on the Mac at
  bundle time."
- 00_front, "6.1 Before any rental", step 4c: "The V1 bundle holds the frozen `lidc-val-edge` view (B1.6) and
  `lidc-val-shifted`, and no `lidc-test-*` folder ... so the test view is moved out before bundling (or a bundle
  `--view` filter is added)". Pass: "equal; no test folder".

Code: `bundle.py` adds every `data/eval_views/lidc-*` folder; `--view` filters only `views/`, not `eval_views/`
[code: experiments/phase2/vast/bundle.py:53-62].

What goes wrong: the A3 list, followed in order, builds the test view into `data/eval_views` and bundles right
after. The V1 bundle then holds `lidc-test-shifted`, step 4c fails, and the engineer must guess which fix to use
(move the folder, or add an `eval_views` filter). The gate stops the error before any money is spent, so this is
major, not blocking.

Exact fix: pick one mechanism. Add `vast.bundle --exclude-eval lidc-test-*`, or have `build_eval_views --split
test` write to `data/eval_views_test/`, which the bundle never globs. Then add that step to the A3 header between
the hash and the bundle commands. The 4c test checks the bundle's file list, not the Mac folder. Change the A3
sentence "reaches the box when it exists" to say that the test view never reaches the V1 box.

### 5. MAJOR. 00_front (6.2 notes "Why items 3-5 always fit"; launch lines): the eval-tail budget leaves out the measured final-model save and the Mac poll between the two `run.sh` calls. That leaves 1-3 min of slack before `timeout` kills item 5.

Doc line: 00_front, "6.2 V1", notes: "The eval tail is items 3-5 (20.8 min) plus about 5.2 min for the untimed
imports of the two eval processes: 26 min ... 172.3 + 26 = 198.3 <= 201. For any offer inside the offer limit,
$2.40 - $2.05 = $0.35 buys at least 30.9 min, more than the 26 min of the tail plus the 3 min that the rounding of
the cap can take." Launch lines: `run.sh <run> b1_rows.yaml # items 2-3; wait for EXIT` then `run.sh <run>
a3_measure.yaml`.

Code and measurement:
- After the last step, `train_rl` merges the adapter and writes `final/` [code: experiments/phase2/train_rl.py:58-60].
  R5 measured "Writing model shards ... 47.62s" [measured: runs/r5/r5_shared.log, last lines]. This time is in no
  item of the tail.
- The second `run.sh` gets `cap = min(cap_hours, deadline - now - MARGIN_S)` [code: policy.py:480]. Every second
  between `EXIT` of the first chain and the start of the second chain comes out of item 5. That covers the Mac
  poll for `EXIT`, ssh and tmux start. The doc gives no bound for this gap.

Arithmetic: slack at planning values = 201 - 198.3 = 2.7 min. For the worst offer inside the limit: 30.9 - 29 =
1.9 min. Minus 0.8 min for the save, 1.1-1.9 min is left for an unbounded Mac poll. When the gap is longer, the
timeout kills item 5, the primary Gate B1 baseline and the A3-val answers. The fallback exists (6.2, V2 item 1b),
but it costs V2 money and weakens Gate B1 to "hash only".

Exact fix:
- Run V1 as one chain. Append the `a3_measure.yaml` step (one `final_eval` process, two `--view`) to
  `b1_rows.yaml` as its last step. V1 then has one `run.sh` call and no Mac gap. `a3_measure.yaml` stays for V2
  item 1b.
- Mark the `train_rl` step soft (`"~", -m, experiments.phase2.train_rl, ...`). A nonzero exit after training
  (item 8) then still lets items 3-5 run on the last saved `final/` or checkpoint.
- In the tail note, add "+ 0.8 min final merge and save [measured: runs/r5/r5_shared.log]". The tail becomes 26.8
  min, and 26.8 + 3 = 29.8 <= 30.9 still holds with `--budget-usd 2.05`.

### 6. MINOR. 00_front (6.2 "five code changes", `push_hf.sh`), B (B1.5 "B1 adapters to V3"): the doc does not say how `cmd_push_hf` passes `--rl`.

B1.5 says "B1.5 adds `hf_push --rl <train_rl out>` ... `push_hf.sh <run> b1_launch.yaml` runs it". But
`cmd_push_hf` builds a fixed command, `hf_push --run .../{output_dir} --repo ... --revision {output_dir}`, with no
`--rl` [code: policy.py:520-521]. Fix: `hf_push` finds the run type itself (`rl.json` present means an RL output,
so it uploads `final/`), or `cmd_push_hf` adds `--rl` when the config has `rows_plan`. Add a unit test on the built
command.

### 7. MINOR. D (D3 step 3): the target-hash check misses the B1 rows root and the frozen views.

D3 step 3 lists `data/prepared`, `data/targets`, `data/eval_views/lidc-val-{shifted,rescaled,rotated}` and
`data/views/lidc-train-{shifted,rescaled}`. It leaves out `lidc-val-edge` and `lidc-test-shifted` (frozen in
`docs/reference/frozen_views.json`, B1.6), and the B1 rows root, whose manifest SHA-256 is pinned on HF. ADR 0015
also says that B1.4 switches to the CUDA polygonizer after K1 parity. Fix: add these three to the
`verify_target_hashes.py` list, plus a rebuild of the rows root that must reproduce `rows_manifest_sha256` byte for
byte.

### 8. MINOR. 00_front (6 inputs, "Start-up per new training process"): the R5 plan job ended with `EXIT` = 137.

[measured: runs/r5/EXIT = 137; ledger row 133 `run ... cap_s=1310` at 16:13:57Z; last HEARTBEAT 16:34:57Z (21.0
min)]. The last `train_rl` step was soft (`~`) and had printed its `rl.json`. So the job was killed (128 + 9) near
its cap, after the work was done. `full_rl_2026-11-08.md:26` says exit code 0 for the vLLM `cumem_allocator`
message, but that is for the full run. Fix: before V1, the 6.1 step 4 tiny-model box check (or V1 item 1's smoke)
records the `train_rl` exit time and code after `rl.json`. Item 5's soft `train_rl` step guards V1 either way.

### 9. MINOR. 00_front (section 0, "outline length"): the "+ 32" in `L - soft_punish_cache >= 1,182 + 32` has no definition.

Section 0 and 99_back 13 use the form, but neither says what 32 is (a margin, or the C11 cap?). Fix: name it in
section 0, for example "32 = the C11 runaway margin", or remove it.

### 10. MINOR. 00_front (6.4, toolchain line): the sanitizer and toolkit versions differ.

`cuda-toolkit[nvcc,cccl,nvtx]==13.0.3` with `nvidia-cuda-sanitizer-api==13.4.92`. Both exist on PyPI [web: PyPI
JSON, 2026-11-08]. CUDA 13.4 maps to driver branch R615, and the box rule asks only for R580 or newer. NVIDIA's
minor-version compatibility page does not cover tools such as compute-sanitizer. Fix: pin the sanitizer to the
toolkit's minor version (`nvidia-cuda-sanitizer-api==13.0.85`, on PyPI). Or require an R615 or newer driver for
R-D2 offers.

### 11. MINOR (top-lab difference). 00_front (section 4, "RL trainer" row): the overlap option is missing.

The GPU is idle about 40 % of the 23.4 s step (BRIEF section 2; C.1). TRL 1.17.1 already ships
`trl/experimental/async_grpo` (`async_rollout_worker.py`, `vllm_client.py`) [code: .venv/lib/python3.12/site-packages/trl/experimental/async_grpo/].
With one GPU and colocate it cannot overlap much, so "not now" is right. Fix: add it to the row's "Considered" list
with the trigger "a second GPU or a server-mode vLLM on the same card". It needs no new library.

### 12. MINOR (top-lab difference). 00_front (6 inputs, "Eval process"), E (E4): evals pay for a trainer build.

About 205 s of fixed cost per eval process goes to the store preload, the pool fork, and the trainer and vLLM
build [estimate in front 6]. V1 pays it twice, and V3 more. A top lab would make the eval a vLLM-only sampler plus
the GPU scorer before the first paid eval. E4's `VllmSampler` is that design, but it lands after E and needs V4.
Fix (optional): move E4's `VllmSampler` and `final_eval --sampler vllm` before V3, as a pre-E path. The CPU
equality test (E-G10 Mac half) plus one V3 item can stand in for V4.

## Counts

Blocking 2 (items 1, 2). Major 3 (items 3, 4, 5). Minor 7 (items 6-12).

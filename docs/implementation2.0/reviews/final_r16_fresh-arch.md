# Final review round 16: fresh outside architect

Reviewer role: a fresh architect with no prior opinion. Scope: the whole plan `docs/implementation2.0/Implementation.md`
(round 16), read through its sources in `phases/`. Focus: is this the best end-to-end architecture up to Phase E
with GPU everywhere, and what a top lab would do differently. I read the r15 fresh-arch review and verdict and the
round 13-15 log entries in `docs/implementation2.0/review-log.md` first. I do not re-raise r14 fresh-arch-3..-8 or
r15 fresh-arch-2/-3 (minor, open), nor the (b') end state (decided; log rows 185, 210, 2265).

## Skills used

| Skill | What it changed |
|---|---|
| nvidia-skills:nvidia-skill-finder | Loaded. No item depends on an NVIDIA SDK choice. The front already names the catalog skills for the 3D step. No catalog lookup was needed, and nothing was installed. |
| engineering:code-review | Read `experiments/phase2/train_rl.py` (whole), `final_eval.py` (whole), `vast/policy.py` (`accrued`/`spend` :190-199, `cmd_launch` :315-345, `chain` :422-452, `cmd_run` :464-495, `cmd_pull` :499-508, `onstart_cmd` :268-283), `serve.fork_pool`/`pool_size`, `reward/advantages.py`, `reward/spec.py`, and TRL 1.17.1 `RepeatSampler` and `vllm_importance_sampling_mode`. Every claim I checked held, except the two minor ones in items 1 and 2. |
| engineering:debug | Traced where V1's A3 answers land on the Mac (`cmd_pull` copies to `runs/<run_id>/`). Then traced who reads them (`runs/a3/eval_val/answers.jsonl` in A2, A8, C2, C11). This found item 1. |
| engineering:system-design | Re-checked the V1 session clock end to end: `started` is set in `cmd_launch`, before `vast create` (:323). `run.sh` caps each plan at deadline minus `MARGIN_S` (:480). `BudgetStop` uses the same rate and fixed terms as `accrued`. I found no gap. |
| engineering:architecture | Wrote the overall judgement below as an ADR-style review: seams, layers, data-first order, and the (b') end state. |
| engineering:documentation | Each item quotes the source file and heading, and the code line or measurement. Item 3 is a wording fix to section 2.4. |
| huggingface-skills:huggingface-best | No model choice is in question (Qwen3-VL-2B now, 4B gate on val). I checked the HF side of the rows root: largest folder is 562 entries, far under the Hub's 10k-per-folder advice [measured: `find data/views/lidc-train-shifted -type d`, 2026-11-09]. |
| system-design-skills:system-design | Ran the loop on "one paid experiment decides the plan". Requirement: attribute a gain to the data fix. Constraint: $3. Failure mode: extra steps, not data, explain the gain. The plan's difference-in-differences rule is the cheapest sound control. A matched-compute control arm is better, but it is not affordable. No defect. |
| engineering:tech-debt | Checked what stays on `main` after E and CLEAN under (b'): one reward (v4), one pixel path (CPU), one polygonizer (engine). I found no duplicate. CLEAN2 removes an idle pool (train_rl never sends work to it; `rewards.py` has no pool reference) [code: train_rl.py:50, 61; grep]. |

Web checks:
- `vastai` on PyPI: latest 1.11.0, uploaded 2026-11-04T19:12:57 [web: https://pypi.org/pypi/vastai/json, fetched
  2026-11-09]. This matches the section 4 pin. `vastai-sdk` is a deprecated wrapper that installs `vastai`
  [web: https://pypi.org/project/vastai-sdk/]. The front row is correct.
- TRL `vllm_importance_sampling_mode`: the default is `sequence_mask`, and mask modes zero the ratios above the cap
  [web: TRL PR diffs via search; code: .venv trl/trainer/grpo_config.py:939-940]. This matches BRIEF and the
  order rule 2 text.

Mechanical check: I extracted all 328 distinct `path:line` citations to repo files in the built doc and checked
that each file exists and that each line is in range
[script: scratchpad `impl/reviews/r16fa/refs.py`]. Five "missing" hits are not errors: `.venv torch/...` paths
written without the `.venv` prefix, an archived lab path written with `...`, and a test fixture string. One real
error remains (item 2).

## Overall judgement

The architecture is sound. A small top-lab team would build something close to it under a $3 cap:
- **Data first answers the measured cause.** The rows root is built and gated on the Mac (21 B data gates). It is
  pinned by HF revision and manifest SHA-256, and checked again by `train_rl` on the box. The paid run trains in
  a fixed order that the gates audited.
- **One CUDA engine source** builds for the CPU and for nvcc. The CPU build is the oracle. The parity test
  `tests/test_gpu_scoring.py` is step 0 of every paid plan. B1.4a puts the polygonizer seam in front of today's
  engine before the targets are built, so the owner rule "CUDA for the polygonizer" holds from V1. D2 later
  swaps K1 in behind the seam, under a byte-identical target-hash gate (D3 step 3).
- **The three seams** (`geom.polygonize`, `policy.sampler`, `rl.trl_adapter`) are at the right places. The layer
  contract becomes a test (import-linter, E0).
- **Money is guarded by code, not prose.** `launch_refusal` and `HARD_STOP` = $19.90 enforce the cap. A
  session-spend `BudgetStop` uses the same formula as `policy.accrued`, so the eval tail is always paid for.

What a top lab would do differently. None of these is a defect under the $3 cap:
1. **A matched-compute control arm.** A top lab would run R5 + 260 steps on centred-only rows next to B1. That
   splits "more steps" from "new data" directly. The difference-in-differences rule (B.6.2) is the cheapest stand-in.
   The plan already says it can be uninformative (the A2 headroom check) and offers the checkpoint-250 variant.
2. **Seed replication before a claim.** A Gate B1 pass is one training run. The plan says so (B.6.2 "Noise floor",
   decision 8b).
3. **Windows sampled on the device from source slices**, not stored as PNG views (ADR-0015 option (b)). This
   matters at the 3D voxel step. The view-bank trigger and ADR-0015 already name it.
4. **A short executable runbook beside the plan.** The plan is 17,141 lines for one $3 rental. A top lab would
   also keep a one-page V1 runbook (commands, stop rules, expected outputs) generated from front 6.2. That keeps
   the paid-time path short for the person at the terminal. Item 4 proposes this as a minor.

## Items

### 1. MINOR. The V1 run id is never fixed, and the A3 answers' Mac path `runs/a3/eval_val/` is made by no step. Today's `runs/` and the ledger already use `v1`, `v2`, `v3`

Sources: 00_front (6.2 "Launch lines"), A (A2 step 1, A3 steps 4-5, A8 step 1), C (C2, C11 step 1).

Problem.
- 00_front 6.2: "**Launch lines** (main agent, after the 6.1 step 0 approval; `<run>` = the V1 run id)". No section
  gives the id.
- `a3_measure.yaml` writes on the box to `"{rd}/eval_val"` and `"{rd}/eval_val_edge"` (A3 plan text). `pull.sh`
  copies `/workspace/runs/<run_id>/` to the Mac folder `runs/<run_id>/`:
  `dst = Path(os.environ.get("VAST_RUNS", ROOT / "runs")) / run_id` [code: experiments/phase2/vast/policy.py:503].
- The Mac readers use another folder. A8 step 1: "`[--answers runs/a3/eval_val/answers.jsonl]`". A2 step 1 has
  the same path. C (C2) reads "`runs/a3/eval_val/answers.jsonl`". C11 step 1 reads "A3's val answers
  (`runs/a3/eval_val/answers.jsonl` with `tokens[K]`)". No step copies or links `runs/<V1 run id>/eval_val` to
  `runs/a3/eval_val`. A3 step 4 says only "on each pulled eval folder".
- Collision risk if the engineer picks the obvious id. `ls runs` gives `k1 p1 preflight r1 r2 r3 r4 r5 rl v1 v2 v3`.
  `runs/v1/` already holds `EXIT`, `HEARTBEAT`, `box.txt` and logs from 2026-11-02 [measured: ls -la,
  2026-11-09]. `ledger.csv` has many `launch`/`run`/`teardown` rows with `run_id` `v1` (RTX 4090, LFM runs)
  [measured: ledger.csv]. `cmd_pull` uses rsync without `--delete`, so V1 output under id `v1` mixes with the
  old files. DOC6 test X2 lists `v1, v2, v3` as today's runs and checks "cost cells equal ledger.csv per run". The
  V1 rental cost would then add to the old `v1` row.

Why minor. Nothing is spent for nothing. Gate B1 reads paths that the plan gives on the box. The first Mac command
that reads `runs/a3/eval_val` fails with "file not found" and points at the cause. An engineer fixes it in a minute.

Fix (00_front, A, C):
1. 00_front 6.2 "Launch lines": "`<run>` = `b1v1`. The id must not appear in `ls runs` or in the `run_id`
   column of `ledger.csv`. 6.1 step 4b checks this." Give V2/V3/V4 the ids `cv2`, `b3v3`, `ev4` in 6.3 in the same
   way.
2. A3 step 4, first line: "After `pull.sh`: `mkdir -p runs/a3 && cp -R runs/b1v1/eval_val runs/b1v1/eval_val_edge
   runs/a3/` (and `runs/b3v3/eval_test` after V3)". A2, A8, C2 and C11 keep `runs/a3/...`.
3. Add one line to the A3 unit tests: a hygiene test reads the run ids from 6.2-6.3 and checks that none appears
   in `ledger.csv`.

### 2. MINOR. B.1 item 6 cites `advantages.py:95-107`, but the file has 105 lines. The "zeroed only when" rule also leaves out the topology-signature condition

Sources: B (B.1 item 6; also the B review-log row 14, which is never graded).

Problem. B.1 item 6: "A group is zeroed only when every term spans less than 0.01 [code:
src/objdel/reward/advantages.py:95-107; src/objdel/reward/spec.py:28-34]". `advantages.py` has 105 lines
[measured: wc -l]. `advantages_v4` is lines 95-105. Its zero rule is `flat = (np.ptp(R, axis=1).max(-1) <
REWARD_V4.drop_spread) & (T["sig_id"] == T["sig_id"][:, :1]).all(1)` [code: advantages.py:103]. So a group is
zeroed only when, in addition, all 8 answers share one topology signature. `spec.py:33` says the same ("and one
topology answer"). The plan's conclusion ("an all-miss group still gets an advantage") does not change.

Fix (B): "[code: src/objdel/reward/advantages.py:95-105; ...]". Write the rule as "only when every term spans less
than 0.01 and all 8 answers have the same topology signature".

### 3. MINOR. Section 2.4 shows the C pixel path as the end state "after B1, C and D", but under the default (b') C never lands on `main`

Sources: 00_front (2.4 heading and diagram; 1 "End state under (b')").

Problem. 2.4 is headed "GPU-everywhere data path (after B1, C and D)". Its box shows "25.2 MB uint8 H2D (pinned),
not 1.61 GB fp32 | C1 patchify ... vision encoder once per image (C3)". Section 1, "End state under (b')", says
"End state on `main` under (b'): reward v4 only; today's CPU pixel path only, with no C flag ... The GPU pixel path
is not reached until a raise funds V2". 99_back decision 1 gives (b') as the default. A reader of section 2 alone
believes the default end state is GPU everywhere for pixels.

Fix (00_front): retitle 2.4 to "GPU-everywhere data path (after B1 and D; the C rows need a raise that funds V2)".
Add one line under the diagram: "Under (b') (the default), the C1/C3 rows stay on the `after-raise` branch. `main`
builds `pixel_values` on the CPU as today (section 1, End state under (b'))."

### 4. MINOR. Nice-to-have: generate a one-page V1 runbook from front 6.2 so the paid path is short for the operator

Sources: 00_front (6.1, 6.2), DOC.

Problem. The paid path is spread over about 1,400 lines of 00_front plus B1.5 and A3 plan text. During V1 the main
agent needs only: 6.1 steps 0-6, the five launch lines, the stop rules per item, and the expected output files
(`b1-rows/final`, `final_adapter`, `eval_b1_val_*`, `eval_val`, `eval_val_edge`, `rl.json`). This is not a defect:
6.2 is the one source, and the planned checks cover it. But a long plan raises the chance of operator error while
the box is billing.

Fix (DOC): `assemble.py` (or DOC4) writes `docs/runbooks/v1.md` from the tagged blocks of front 6.1-6.2. It holds
only commands, stop rules and the expected-file list, plus the run id of item 1. A drift check compares it with
6.2 on each build.

## Counts

blocking 0, major 0, minor 4.

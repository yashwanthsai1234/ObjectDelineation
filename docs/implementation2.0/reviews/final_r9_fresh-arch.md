# Final review round 9: fresh outside architect

Reviewed: `docs/implementation2.0/Implementation.md` (built from `impl/phases/*.md`), round 9. Focus: whether this is
the best end-to-end architecture to Phase E with GPU everywhere, plus any blocking item anywhere. I did not repeat
items that `docs/implementation2.0/review-log.md` already lists as fixed.

## Skills used

| Skill | What it changed |
|---|---|
| nvidia-skills:nvidia-skill-finder | No item needs an NVIDIA catalog skill. The front already names `compute-sanitizer-diagnose-cuda-bugs` and the 3D skills. I did no catalog lookup and installed nothing. I checked the "no GPU PNG decoder" claim on the web instead: NVIDIA lists only nvJPEG, nvJPEG2000 and nvTIFF as GPU codecs in nvImageCodec, so the claim stands. |
| engineering:code-review | I read the code behind every V1 plan step: `policy.chain` (the `str(x).format(rd=rd, data=data)` line), `cmd_run`, `code_files`, `run_box.sh`, `Makefile` `cuda`, `CMakeLists.txt`, `train_rl.py`, `final_eval.py`, `bundle.py`, `preflight.py`. This produced item 1. |
| engineering:debug | I ran the A3 plan step through Python's `str.format` with the runner's two keys. It raises `KeyError(' bash experiments/cuda_checker/run_box')` (item 1). |
| engineering:system-design | I traced the V1 box path in order (launch, sync, fetch, `b1_rows.yaml`, `a3_measure.yaml`, pull, push, teardown), then the (b') end state on `main` through E and CLEAN (item 3). |
| engineering:architecture | I compared the adoption rules of C and D under (b') (item 6) and checked whether the Gate B1 design separates data from extra steps (item 7). |
| engineering:documentation | Each quote below is checked against its source file under `impl/phases/`, with the line number. |
| huggingface-skills:huggingface-best | No Hub leaderboard covers polygon-text segmentation. Section 3's choice (keep Qwen3-VL-2B, gated 4B arm) is sound for $3. I fetched arXiv 2601.13633 to check that the cited paper exists (EGM, Qwen3-VL grounding). No change. |
| system-design-skills:system-design | Stress test: what breaks first in V1 if one step fails, and what the (b') end state leaves on `main` (items 1, 3, 9). |
| engineering:tech-debt | Old code beside new code under the default budget (item 3); the size of the plan against one $3 experiment (item 8). |
| WebSearch / WebFetch | nvImageCodec codec list (NVIDIA product page); arXiv 2601.13633 abstract page. |

## Verdict in one paragraph

The architecture is sound for the money. Data comes first. One CUDA engine source serves the CPU oracle and the GPU.
Only scores leave the device. TRL stays behind one adapter. Every rental item has a price and a stop rule, and the
(b') end state is stated honestly. One new blocking defect remains: the round-8 `a3_measure.yaml` step cannot run
through the plan runner (item 1). Two major items remain: a stale decision-interval rule in B1.9 (item 2), and old
samplers left on `main` indefinitely under (b') (item 3). What a top lab would do differently is in items 6-8. All of
those are minor at this budget.

## Items

### 1. BLOCKING: the `a3_measure.yaml` build-if-missing step kills `run.sh` (str.format braces)

Sources: A (A3, "Plan file `a3_measure.yaml`"), 00_front (6.2 plan files and launch lines; 6.3 item 1b).

Doc, A.md line 495 (A3, plan file `a3_measure.yaml`, step 2):
`- [bash, -c, "/opt/venv/bin/python -I -c 'import objdel._gpucheck as g; assert g.build() == \"cuda\"' || { bash experiments/cuda_checker/run_box.sh --build-only && /opt/venv/bin/python -I -c 'import objdel._gpucheck as g; assert g.build() == \"cuda\"'; }"]`

Front 6.2 launch line: `run.sh <run> experiments/phase2/configs/a3_measure.yaml    # items 4-5`.

Code: `cmd_run` calls `chain(ROOT / a.module, a.run_id)` on the Mac, and `chain` fills every step with
`return [[str(x).format(rd=rd, data=data) for x in s] for s in c["steps"]]` [code:
experiments/phase2/vast/policy.py, `chain`, the `if "steps" in c:` branch]. The bash group `{ ... ; }` is read as a
format field.

Measurement (2026-11-09): `s.format(rd='R', data='D')` on that exact string raises
`KeyError(' bash experiments/cuda_checker/run_box')`. E's own rule confirms the behaviour: "an unknown placeholder
such as `{x}` raises" (E.md E6, `to_argv`). The review log notes it too: "`chain()` fills placeholders with
`str.format` [...] so an unknown `{x}` raises" (E review log, own findings). Step 2 came in round 8 (qa-1), after
that note.

Effect: an engineer who writes the file as the doc gives it gets one of two outcomes.
- If I20 calls the runner's own function ("asserted through the runner's own function"), I20 fails on the Mac. The
  doc's plan text and its own test then disagree, and the V1 tag cannot pass `make test` (6.1 step 1).
- If I20 checks only with a regex, V1 reaches the box. After `b1_rows.yaml` ends, `run.sh ... a3_measure.yaml`
  raises on the Mac before any remote call. Items 4-5 are "Not skippable", and they are the Gate B1 paired baseline
  (primary rule), Gate A, the Mac half of C-G-EQ2/3, B3.7 and C11. Gate B1 falls back to R5's old `items.jsonl`
  ("hash only, correlation not run"), and V2 item 1b (the same file) fails the same way.

I47 does not catch it: it parses only `-m` steps.

Fix:
- A.md A3 plan file: use a subshell, which has no braces:
  `- [bash, -c, "/opt/venv/bin/python -I -c 'import objdel._gpucheck as g; assert g.build() == \"cuda\"' || ( bash experiments/cuda_checker/run_box.sh --build-only && /opt/venv/bin/python -I -c 'import objdel._gpucheck as g; assert g.build() == \"cuda\"' )"]`
  (`{{ ... ; }}` also works, but then E6's `to_argv` port must keep the escape; the subshell needs none.)
- A.md I20: add "`policy.chain(Path('experiments/phase2/configs/a3_measure.yaml'), 'r')` and the same call on
  `a3_test.yaml` return without an exception".
- B.md I47 and C.md U-CF4 (or one new test in `tests/test_preflight.py`): call `policy.chain(p, 'r')` on every
  `steps:` plan that the front's launch lines name (`b1_rows`, `a3_measure`, `a3_test`, `c_speed`, `r_d1`,
  `b2_fp16`, `b3_v5`, `d_gpu_test`) and expect no exception. This is a $0 guard for every rental.
- Front 6.1 step 4d: add "every plan file of 6.2-6.3 passes `policy.chain` on the Mac" to the pass column.

### 2. MAJOR: B1.9 gives two different decision intervals for the Gate B1 simulation

Source: B (B1.9 "Gate power check"; B.6.2 common rules), A (A1 step 5).

Doc, B.md line 1027 (B1.9): "apply every Gate B1 hard rule that reads val-shifted item values, exactly as B.6.2
writes it (paired, the A1 decision interval, same B and seed)".
Doc, B.md line 1046 (B1.9, last bullets): "The decision interval of the simulation is the CR1 interval of the
item-weighted mean (B.6.2), the same as the gate's."
Doc, B.md line 1767 (B.6.2): "the percentile patient interval when the gate set has >= 50 patients (Gate B1, whole
val, 78 patients), and `signal.cluster_t` when it has < 50 (Gate B3, val-select, about 39 patients)".
Doc, A.md line 292 (A1 step 5): "The percentile interval is the decision value on the whole val split (78
patients)".

So Gate B1 decides with the percentile patient bootstrap, but B1.9's last bullet simulates it with CR1. The review
log shows the cause: round 3 wrote "B.6.2 and B1.9: CR1 interval" (B review log, round 3, item 4), and the later
percentile/CR1 split changed B.6.2 but not that bullet.

Why it changes the outcome: B1.9's joint P(pass) decides whether the owner is asked to change a rule before V1 ("If
joint P(pass) [...] < 0.5 for the 25 % effect in the c = 0.077 run, the owner decides before the rental whether to
change a rule"). On R5 val the two intervals differ: 2-78 items per patient, and the top 5 patients hold 24.9 % of
the items (A1 step 5, measured). The engineer must guess which interval to use, and the answer near 0.5 can flip.

Fix: B.md B1.9, replace the bullet with: "The decision interval of the simulation is the gate's own interval
(B.6.2, A1 step 5): for Gate B1 (whole val, 78 patients) the percentile patient bootstrap, same B and seed; for the
report-only Gate B3 line (val-select, < 50 patients) `signal.cluster_t` (CR1). The other interval is printed
beside each P(pass)." Add a B1.9 test assert: the Gate B1 rules in `gate_power.py` call the same function as the
gate (`gate_ci`), not `cluster_t`.

### 3. MAJOR: under the default (b'), old samplers stay on `main` beside the new one with no end

Sources: E (E4 step 1 "Archive after the box check"; E4 gate; E.8), 00_front (section 1 "End state under (b')"),
99_back (decision 3b).

Doc, 00_front line 188 ("End state under (b')", "the one rule"): "`main` then never holds a temporary flag, a second
reward or a second pixel path, and no old code sits beside new code (owner archive rule)." The same section lands
"E" on `main` under (b').
Doc, E.md line 993 (E4 step 1): "`experiments/phase2/sample_vllm.py` and `experiments/phase2/hf_sample.py` go to
`archive/phase2_samplers/` [...] This is a follow-up commit, made only after E-G10's box part passes (V4, E.8).
Until then they stay".
Doc, 00_front 6.3: "V4 (optional, after E, owner decision; needs a budget raise [...])"; 99_back decision 3b: "[no
default: owner decision]".

Under (b') no rental follows V1 (front 6.3 "Rule"), so V4 never runs and the two old samplers stay on `main` after
E, the DOC sweep and CLEAN, with no end condition. After E4 step 3 they have no caller: `final_eval` uses
`VllmSampler`, and `sft.select.sampler_cmd` builds `-m objdel.policy.sampler --kind vllm|hf` (E4 step 1). Today
their only callers are the select bridge and docstrings [code: experiments/phase2/select_ckpt.py:10-13, 55, 88-95;
grep of src, experiments and tests, 2026-11-09]. That is old code beside new code, against BRIEF section 3 ("Old
code must be ARCHIVED (not left beside new code) [...] Remove dead code") and against the front's "one rule". The
engineer must guess whether CLEAN16's final sweep archives them or leaves them. Round 7 (fresh-arch-4, CONFIRMED
major) graded the same pattern major: a kept path with no trigger and no end.

E also contradicts itself on whether E4 is done: E.md line 1028 says "**E4 is not done until E-G10's box part
passes.**", while E.8 (line 1726) says "The E phase gate passes on the Mac halves. E4 is done, and its follow-up
commit archives [...] only after V4 passes."

Fix (pick one and write it in E4 step 1, E.8, front "End state under (b')" and 99_back 3b):
- (preferred) Under (b'), archive both samplers in E4 step 3's commit. The fallback stays in git history under the
  tag `pre-e4-eval`, which E4 step 3 already asks the owner to set. The next rental after a raise runs E-G10's box
  part first. If it fails, the owner restores from the tag. Front "End state under (b')" then lists "E4 archives the
  old samplers at step 3; E-G10's box half stays open until a rental".
- or name the end: "CLEAN16 archives `sample_vllm.py` and `hf_sample.py` if V4 has not run by then", and add that
  to CLEAN16's step list and to its hygiene test.
Make E.8 match the E4 gate: "E4's Mac half is done; E-G10's box half stays open until V4."

### 4. MINOR: section 2.4 rows-root size

Source: 00_front (2.4 diagram, line 366). The diagram says "B1 rows root (9,600 rows)". Section 0 says "8,320 rows =
260 steps x 32", and B1.2 says "The 9,600-window sizes and times in B1 are the 300-step upper bound". Fix: write
"(8,320 rows; B1.2)".

### 5. MINOR: V1 `BudgetStop` margin at the base step time is 0.1 min

Source: 00_front (6.2 "Why 260 steps"). "260 steps end at minute 53 + 119.2 = 172.2", and the stop is at 172.3. Any
slowdown (FR3, the shifted-row decode tail; B2.1 logging) trims the last 20-step block. The doc handles the money
("the adapter is saved and items 3-5 still run"), but it does not say how a cut block affects the coverage gates
B-DG7..B-DG10, which audit all 13 blocks. Fix: one sentence in B1.2/B.6.2: "If `BudgetStop` ends the run inside a
block, the report gives the trained step count and the coverage of the trained rows (`coverage.report` on the
first N x 32 rows)". The alternative is 240 steps (12 whole blocks), which gives about 9 min of margin at 27.5 s.

### 6. MINOR (architecture): two adoption rules for GPU code under (b')

Sources: 00_front (section 1 "End state under (b')"; 6.3 V3 item 1), D (D.4, D-G7), C (C.6).

Under (b'), D2-D8 land on `main` with no GPU run at all ("D's box gates stay 'open until a rental'"; R-D1 is
dropped; R-D2 is V3 item 1 and "Skippable"). C1-C5, which also have CPU-side oracles (K3's CPU reference "0 of 64.5
M values differ"), wait on `after-raise`. So the CUDA kernels that score the reward reach production untested on a
GPU, while the GPU pixel path does not. A top lab would use one rule: code with a bit-exact CPU oracle lands, and
the first rental runs that code's box gate before any training step. Fix: in front 6.3 make V3 item 1 (R-D2) "not
skippable when D2-D8 merged after the last GPU run", or move `tests/test_gpu_scoring.py` plus D's GPU02 subset into
`b3_v5.yaml` step 0 (it already holds `test_gpu_scoring.py`), and state the rule once in section 1.

### 7. MINOR (science design, top-lab view): Gate B1 has no matched control arm

Sources: B (B.6.2 Gate B1 "extra-steps control"), 00_front (section 1.1). B1 continues R5 (`--ckpt {rd}/ckpt`,
`rl-full-3ep`) for 260 more steps, so data and extra steps change together. The difference-in-differences proxy
assumes extra centred steps lift near and far items alike. Near items already sit at or above the doctor ceiling
(0.76-0.77 against 0.71-0.74 [measured: docs/research/notes/rl_plateau_diagnosis_2026-11-08.md]), so the near gain
is near 0 whatever the cause. The doc already flags this ("the report calls this rule uninformative"). A top lab
would buy a matched control (R5 plus 260 centred steps) or a second seed before buying B3. Fix: in 99_back decision
1, option (a), list "control arm, about $1.35 at 27.5 s (V1 item 2 price)" beside B-S2, so the owner can choose
attribution over a reward change.

### 8. MINOR (tech debt, top-lab view): plan size against one paid experiment

Sources: all. The plan is 16,124 lines for one $2.40 experiment, and E, DOC and CLEAN make up most of the Mac work.
A top lab would freeze V1's scope, run it, and re-plan C, B3 and E from the result. Fix: none required. A
"minimum path to V1" list in section 1 (DOC1-2, A1-A8, B1.0-B1.9, B2.0-B2.1, CLEAN1-2, 6.1 step 0) would help the
engineer.

### 9. MINOR: a V1 `env_check` failure has no stop rule

Source: 00_front (6.2 item 1; 6.4). `b1_rows.yaml` step 1 runs `env_check.sh`, which "exit[s] 1 unless every check
passes" [code: experiments/phase2/docker/env_check.sh:3], including `bf16_matmul >= 0.6 of the dense BF16 peak`. A
slow host then stops the `&&` chain after item 1 ($0.70 spent), and 6.4 names no action. Fix: add to 6.4 rule 3: "a
failed `env_check` or GPU smoke in item 1: tear down, record the host in the ledger, and pick another offer only
with the owner's approval of a second launch".

## Counts

blocking 1, major 2, minor 6.

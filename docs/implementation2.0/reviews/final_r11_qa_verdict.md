# Verdict on final_r11_qa.md (fresh verifier, 2026-11-09)

## Skills used
| Skill | What it changed |
|---|---|
| nvidia-skills:nvidia-skill-finder | No NVIDIA product choice in scope; no catalog skill recommended. No change to verdicts. |
| engineering:code-review | Read `experiments/phase2/vast/policy.py:420-510` (`chain` returns `steps:` literally; `cmd_run` joins all steps with `&&`) against items 2 and 3. |
| engineering:debug | Traced the `c_speed.yaml` chain for a mode A crash before step 18 (item 3) and the inputs `--expect-step1` receives (item 1). |
| engineering:system-design | Checked V2 item 3/3c and V3 item 0 rows and launch lines against the plan files that carry the steps. |
| engineering:architecture | One owner per check: G-RC0 and G-GR1 rows name a rental but no plan file owns the step (item 2). |
| engineering:documentation | Checked that the doc does not already define a machine-readable G-EN0 result or a mode B check arm elsewhere (grep over all phase sources). |
| huggingface-skills:huggingface-best | No model choice in scope. No change. |
Web tools loaded (WebSearch, WebFetch); no external claim needed checking.

Measurement this round: `runs/r5/base_rl.json` step 1 `grad_norm` 0.003265785751864314, `runs/r5/shared_rl.json`
0.00326397898606956. Rounded to 4 s.f. (0.003266 vs 0.003264) they fail `math.isclose(rel_tol=5e-4)` and pass 1e-3.

## Verdicts (blocking/major items)
- qa-1, confirmed (major): C.md C0 files says `--expect-step1 FILE` "applies the E-EQ6 comparison, the same function that the gate uses", but E-EQ6 (C.9.5) and C.6 G-EQ6 switch `grad_norm`/logp tolerance on the G-EN0 result, which C.9.4 G-EN0 writes only "to the check log" (grep finds no file/key a program can read), so the in-process stop has no input for the rule and the strict branch would kill the fast arm on a 5.5e-4 move (measured R5 0.0032658 vs 0.0032640); U-CF5 plants only an `iou_tc` mismatch, so no $0 check catches it.
- qa-2, confirmed (major): C.md C.9.4 row "G-RC0 | V3, before the B3 run | pytest `-k rc0`" and C11 step 3 "rental V3 item 0", but the V3 launch lines (00_front 6.3) have no item 0 `run.sh` call, `b3_v5.yaml` (B.md B3.9/I22b) has steps CUDA build, G08b, `train_rl`, `final_eval` only, and `chain` runs `steps:` literally [code: policy.py:432-433]; `grep rc0` finds only C.md rows; G-GR1 ("V2 optional C6 arm | pytest `-k gr1`") likewise has no step in `c_speed_c6.yaml`; `-k` with no file list also breaks B.9.5 "`pytest -m gpu` without a file list is never run on a box".
- qa-3, confirmed (major): C.md C.6 "Check arm = steps 18-32 of the fast arm's process (`--steps 40 ... --record-batches {rd}/c_check/arm_batches`)" and 00_front 6.3 item 3c "`c_speed_b.yaml` ... (17 steps read, `--steps 20`)" give mode B no check steps, while G-EQ2/G-EQ3 need ">= 15 come from the box" and C.9.4 says "fewer than 15 box batches fails"; a mode A OOM (C5 step 2 plans ~28 of 32 GB) ends the single `&&` job [code: policy.py:487-494] before step 18, so after 3c's $0.18 Gate C-adopt cannot pass and no text says what 3c must add.

## Minor items (not re-graded upward)
- qa-4 .. qa-10: stay minor as filed (sanitizer marker caught by D-G5/G07 setup, exit-code mapping, totals caught by 6.1 step 4b at $0, wording, stale header/configs, doc size).

Totals: 3 confirmed (all major), 0 refuted, 0 downgraded.

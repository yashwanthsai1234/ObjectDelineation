# Verdict on final_r6_arch.md (fresh verifier)

## Skills used

| Skill | What it changed |
|---|---|
| nvidia-skills:nvidia-skill-finder | No item makes an NVIDIA product claim, so no catalog skill was needed. I installed nothing. |
| engineering:code-review | I opened every cited code line myself: train_sft.py:150-157 (resolve assigns to args) and :272 (bare `resolve(args)`), env_check.sh:9 (MODEL_DIR required), policy.py:527-536 (fetch_hf pulls to `{rd}/ckpt`), runs/r5/eval/eval_val/final_eval.json `ckpt` = `/workspace/runs/r5/final_model`, docs/research/notes/final_eval_2026-11-08.md:3 (`rl-full-3ep`). |
| engineering:debug | I re-ran the reviewer's five probes in its scratch venv (pydantic 2.14, pydantic-settings 2.16.2): probe.py exits 2 "required: --kind, --name"; probe2 raises extra_forbidden for the keyword form; probe3 prints True; probe4 shows the dropped return value; probe5 parses interleaved views into tuples. |
| engineering:system-design | I traced the V4 rental end to end (launch, fetch, env_check, eval, push, teardown) against front 6.3 and E.8. That set item 4's severity. |
| engineering:architecture | I checked whether any rule in E, A or B already makes EvalStage multi-view (item 1). None does: E.md has 0 hits for a tuple view or a per-view loop. |
| engineering:documentation | I ran the check 14 regex on the built Implementation.md, review-log.md and kernels/, and counted hits per source page (00_front 5, B 1, D 4). |
| huggingface-skills:huggingface-best | No model choice is involved, so nothing changed. |

## Verdicts

- arch-1: confirmed (blocking). E.md E6 step 1 defines `EvalStage ... view: str; split: Literal["val","test"]; out: str` (extra="forbid"), while E6 step 2 and step 7 / I13 require `to_argv` on the ported `b1_rows` to equal golden 5 string for string, and B.md B1.5 step 4 has two interleaved `--view ... --split val --out ...` groups (A.md A3: `action="append"`); a single-field stage cannot hold the second view, so two parts that both claim to be the rule contradict (probe5 shows the tuple form works).
- arch-2: confirmed (major). E.md E2 step 1 `StageBase` has `kind: str` and `name: str` with no default, E6 step 1 `SftStage.kind: Literal["sft"]` has no default, and the Flag order rule never writes `kind`/`name`; probe.py reproduces `CliApp.run(SftStage, cli_args=...)` exiting 2 with "the following arguments are required: --kind, --name", so the E6.0 spike, I14 (`S1 == S2`) and the preflight dry run cannot pass as written and the engineer must invent the CLI default and the comparison rule.
- arch-3: downgraded (minor). The missing rebind (train_sft.py:272 `resolve(args)` drops the return value) is real, but the doc already says `resolve(stage) -> SftStage` returns a new copy, the one-line fix is implied by that signature, and the reviewer itself notes the Mac CPU dry run (test_cpu_smoke_train_then_select) fails before any rental, so there is no wrong result and no GPU spend.
- arch-4: confirmed (major). 00_front.md 6.3 "V4" gives no launch, fetch or teardown lines and E.md E.8 step 3 gives no `ckpt` for the EvalStage, while E4's gate says "checkpoint = that file's `ckpt` key", which is `/workspace/runs/r5/final_model` on the dead R5 box; env_check.sh:9 requires MODEL_DIR and nothing in V4 runs `fetch_hf.sh <run> rl-full-3ep`, so on a paid box the engineer must guess the checkpoint source and path (a wrong revision breaks the TOST), and copying V1's `pull && push_hf && teardown` line dies at push by E6.8b design and skips teardown.
- arch-5: confirmed (major). DOC.md DOC4 check 14 and DOC-G13 require 0 matches of `/[p]rivate/tmp|/[U]sers/|S[P]/impl` in plan files and `docs/implementation2.0/**/*.md`, and DOC1 step 1 says "Step 0 removes them", but step 0 only replaces path prefixes; the built Implementation.md still matches at lines 67, 544, 889 (00_front), 8498, 9254 (D) and review-log.md at 5 lines (00_front, B, D logs) as literal mentions, so the gate and `make test` stay red with no stated step to fix them.
- arch-6: minor as filed (not contested).
- arch-7: minor as filed (not contested).
- arch-8: minor as filed (not contested).
- arch-9: minor as filed (not contested).
- arch-10: minor as filed (not contested; round 5 minors, already minor).

Counts: confirmed 4 (blocking 1, major 3), refuted 0, downgraded 1.

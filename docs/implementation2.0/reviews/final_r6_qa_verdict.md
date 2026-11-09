# Verdict on final_r6_qa.md (round 6, QA lead)

## Skills used

| Skill | What it changed |
|---|---|
| nvidia-skills:nvidia-skill-finder | No NVIDIA product claim to check; pointed at the box toolchain path (sync_up installs the CPU build, step 0 builds CUDA) when re-reading the V4 step list. |
| engineering:code-review | Read the code behind each claim: policy.py:373-384 (code_files), 440-492 (chain, cmd_run `&&` chain), 517-536 (push/fetch), 308-316 (stamp), .gitignore:23, env_check.sh:1-10, final_eval.py:48-56, fast.py:567-595, tests/test_p2_scripts.py:391-393, tests/test_p2_train_rl.py:1-35. |
| engineering:debug | Traced V4 on a fresh box step by step (step 1 needs MODEL_DIR; GP5 inputs at step 2) and the B2.4 call order (decide_numerics before the first post-B3 run). |
| engineering:system-design | Put B2.4 rule 0, "Missing inputs", U51 and the B.9.8 sketch side by side. |
| engineering:architecture | Checked the launch-config role (front section 0, 6, 6.1 step 4) against every rental V1-V4. |
| engineering:documentation | Checked whether any other part of the doc already handles each item (none did). |
| huggingface-skills:huggingface-best | No model choice involved; only the checkpoint fetch for V4 (fetch_hf from the Mac). |
| Web | Not needed: every claim was settled by repo code and the review's own measurement (all-skipped pytest exits 0, consistent with pytest ExitCode.OK; 5 only when nothing is collected). |

## Verdicts (blocking and major items)

- qa-1 (blocking): CONFIRMED blocking. E.8 V4 runs GP5 at step 2 before the sampler arm (step 3), the trainer-arm file is under git-ignored `/runs/` [.gitignore:23] so `sync_up` (`git ls-files --exclude-standard`, policy.py:376-384) never ships it, GP5 skips by the E.9 rule and the step passes, and no step or Mac test computes the TOST afterwards, so the $0.58 V4 rental cannot close E-G10/E4.
- qa-2 (blocking): CONFIRMED blocking. B2.4 rule 0 prose starts the run under `lora-seqmask` and runs the F check inside it, while B2.4 "Missing inputs", U51 ("Inputs absent ... -> bf16-tis-opsm") and the B.9.8 case `(0.531, None, None, "bf16-tis-opsm")` return `bf16-tis-opsm` in exactly the pre-run state; two parts both claim to be the rule and ADR-0013 records the wrong preset.
- qa-3 (major): CONFIRMED major. Front 6.3 V4 and E.8 V4 give no launch lines, no `fetch_hf.sh`, and no argument for `env_check.sh` (which dies without MODEL_DIR: `MODEL=${1:?usage...}`, env_check.sh:9) or `--ckpt` for `final_eval` (required, final_eval.py:48); the engineer must guess the checkpoint fetch, else step 1 fails and V4 money is lost.
- qa-4 (major): CONFIRMED major. Front section 6 claims to list every stamped launch config but names only `b1_launch.yaml` (6.1 step 4, 6.2); V2 and V3 (both before E6, so RunArgs-era) have no launch config or launch lines, `cmd_launch` refuses without a stamp of its `--config` (policy.py:313-316), and `push_hf` takes its revision from that config's `output_dir` (policy.py:517), so the engineer must invent the file and the B3 push target.
- qa-5 (major): CONFIRMED major. R01b/R01c cite "the tiny TRL loop", which today runs `GRPO_SMOKE` with `use_vllm: False` and 8 rows per step (test_p2_scripts.py:392-393), and FastGRPO feeds `set_is_ratio` (the only writer of the IS stats, `seq_gaps` and `cf/*`) only under `self.use_vllm and self.vllm_importance_sampling_correction` (fast.py:568); no fake-vLLM fixture is named in B, so the "3 x 256 x 3" assertions cannot pass as written and the engineer must invent the harness.

## Minor items (not re-graded; kept minor as the review filed them)

- qa-6 minor (fallback A3 hash file): kept minor.
- qa-7 minor (A3 step 5 "stop"): kept minor.
- qa-8 minor (99_back 3b old V4 time): kept minor.
- qa-9 minor (`e_v4.yaml` outside `specs/`): kept minor.
- qa-10 minor (D1 "changes" a file D8 creates): kept minor.
- qa-11 minor (S09 used twice): kept minor.
- qa-12 minor (V4 item "E checks" name): kept minor.

Totals: 5 confirmed (2 blocking, 3 major), 0 refuted, 0 downgraded.

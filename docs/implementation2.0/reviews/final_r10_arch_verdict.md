# Verdict on final_r10_arch.md (round 10, architecture)

## Skills used

| Skill | What it changed |
|---|---|
| nvidia-skills:nvidia-skill-finder | No item makes an NVIDIA product claim. No catalog lookup, no install. |
| engineering:code-review | Read the cited plan text (CLEAN12, CLEAN8, CT-44, CT-62, E2, E8 "Coverage gate", E-G1, U55, C7, front section 3) and repo code (`Makefile:18-26`, `.github/workflows/tests.yml:30-44`, `pyproject.toml` markers). |
| engineering:debug | Re-checked the reviewer's probes: `reviews/r10arch_probe/out.txt` (vulture names `_save` at 60 %) and `reviews/r10arch_probe/cov.json` per-function data (train_sft `train` 97/74 missed, `make_trainer` 20/20, `train_rl.main` 38/23). Confirmed `make cov` runs `-m "not slow"`, so the slow X1 run does not cover the moved code. |
| engineering:system-design | Walked the front time line: 4B gate closes before E; E; then CLEAN wave 3. So a 4B row or a C7 limit can exist on the tree when CLEAN8 runs. |
| engineering:architecture | Ownership of the per-model table: front section 3 step 1 adds the 4B row "before E" or "after E"; CLEAN8 step 3 and CT-62 overwrite it. |
| engineering:documentation | Checked whether E or CLEAN already handles items 1-3 elsewhere (grep of E.md for `vulture`, `cov`, `U55`). |
| huggingface-skills:huggingface-best | No model choice in scope. Nothing changed. |
| WebFetch/WebSearch | Loaded; not needed (vulture behaviour proven by the local probe output). |

## Verdicts

- arch-1 (major) - downgraded to minor: the probe is right (vulture flags `@model_validator` methods), but the failure is loud and self-explaining, no test pins the exact ignore list (CT-44 checks only comment lines, paths, confidence and a canary), and the natural fix (add the decorator to `ignore_decorators` with a comment, the mechanism CLEAN12 already uses) does not change any outcome.
- arch-2 (major) - confirmed: E8 "Coverage gate" records the E0 figure only for the E8.1 infra files, yet E2/E4 move ~423 statements with ~232 missed (measured `r10arch_cov.txt`/`cov.json`: `train_sft.train` 97/74, `make_trainer` 20/20, `train_rl.main` 38/23) into `objdel` under `make cov --cov=objdel --cov-fail-under=90 -m "not slow"` [Makefile:25-26], giving about 89.5 %; E-G1 binds on every E commit and the E8 rule's own E0 base omits this code, so following it exactly can still leave E8.1 red.
- arch-3 (major) - confirmed: CLEAN8 step 3 writes `MODELS = {"qwen3-vl-2b": ...}` and `MAX_COMPLETION = {"qwen3-vl-2b": 1536}` and CT-62 asserts `set(MODELS) == {"qwen3-vl-2b"}`, while front section 3 step 1 adds a `qwen3-vl-4b` row before or after E (the 4B gate closes before E, CLEAN wave 3 runs after E) and C7 lets the owner change `MAX_COMPLETION`; following CLEAN8 exactly deletes the 4B row or resets L.
- arch-4 (minor) - stays minor: filed as minor; the I42 re-point is obvious once the test goes red.
- arch-5 (minor) - stays minor: both orders work, as the reviewer says; wording only.
- arch-6 (minor) - stays minor: wording of R-4; the `extra="forbid"` error names the missing field.
- arch-7 (minor) - stays minor: citation durability, no outcome change.
- arch-8 (minor) - stays minor: stale byte count, still under the limit.
- arch-9 (minor) - stays minor: timing words; E6 step 8 already deletes `HF_RATE` in E6.6.

Counts: confirmed 2 (arch-2, arch-3, both major), downgraded 1 (arch-1), refuted 0, minor as filed 6.

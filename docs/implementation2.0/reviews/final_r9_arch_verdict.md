# Verdict on final review round 9: architecture

## Skills used

| Skill | What it changed |
|---|---|
| nvidia-skills:nvidia-skill-finder | No item makes an NVIDIA product claim. No catalog lookup, nothing installed. |
| engineering:code-review | Opened `policy.py:312-317, 470-476`, `train_sft.py:78-90`, `tests.yml:34-44`, `Makefile:20-30`, `pyproject.toml:73-77`, `preflight.py` function list. |
| engineering:debug | Re-ran the env probe (`r9arch_probe/env_probe.py`) with pydantic-settings 2.10.1: `budget_usd=2.05 session_start=None`, same as the reviewer's 2.16.2 result. |
| engineering:system-design | Walked the raise-after-E0 path: front "End state under (b')" -> order rule 1 -> E R-4 -> E6 steps 2, 8, 8b -> I29. |
| engineering:architecture | Checked the coverage gate against the E8.1 move (combined 88.7 % measured on be52b1c). |
| engineering:documentation | Each verdict names the source heading and the code line or measurement. |
| huggingface-skills:huggingface-best | No model choice is in scope. Nothing changed. |
| WebFetch | pydantic-settings docs: `CliApp.run` on a `BaseModel` uses BaseSettings defaults that do not mention env sources; the local probe is the proof. |

## Verdicts

- arch-1: confirmed (major). E.md "## E.8" dependencies says "B2 ... and B3 ... land before E on pre-E paths (front order rule 1; no post-E exception)" while 00_front.md "Order rules" 1 and E R-4 name the `after-raise` exception as "a mechanical change"; after E6.8 `chain()` reads only `specs/*.yaml` (I29) and preflight reads only a `RunSpec`, so `c_launch.yaml`/`b3_launch.yaml` (front section 0) are refused by `preflight.refusal` [code: policy.py:315], and no text says how to port the branch's plans, launch configs and post-E0 flags.
- arch-2: confirmed (major). B.md B1.5 makes `OBJDEL_SESSION_START` etc. the `train_rl` flag defaults and `b1_rows.yaml` omits `--session-start`; E6's `CliApp.run` on a `BaseModel` stage reads no environment [measured: env_probe.py, pydantic-settings 2.10.1 here, 2.16.2 by reviewer: `session_start=None`], and E.md has no `BudgetStop`/budget/env rule (grep: 0 hits), unlike the explicit SFT `default_factory` rule (E6 "SFT argv"; train_sft.py:81-85), so the engineer must choose between breaking E-G9 and silently losing the paid-run spend guard.
- arch-3: confirmed (major). E.md "### E8" step 1 moves `vast/policy.py`, `preflight.py`, `bundle.py`, `hf_push/pull` into `src/objdel/infra`, which `--cov=objdel --cov-fail-under=90` [code: tests.yml:44; Makefile:25-26] then measures; combined coverage is 4,608 stmts / 522 missed = 88.7 % [measured: r9arch_cov.txt, policy.py 50 %, preflight.py 49 %]; E6.6's deletion of `projection`/`ceiling`/`HF_RATE` (about 35 statements) does not close the 62-statement gap, and E only handles the gate for `verify/http.py` (E7 step 5).
- arch-4: confirmed as minor (not escalated). Stale "pool size 1 (as A3)" wording at E.md E.8 V4 step 3; no outcome change.
- arch-5: confirmed as minor. Refusal text names an archived plan [code: policy.py:473-474]; message wording only.
- arch-6: confirmed as minor. `LIVE_CLIS` owner and list completeness; the E8 gate test catches a missing name at $0.
- arch-7: confirmed as minor. `VllmSampler` constants unstated for `final_eval`; V4 is optional and the values exist in code.
- arch-8: confirmed as minor. Diagram signature drift only.
- arch-9: confirmed as minor. U52 set detail; the test itself reports the set.

Counts: confirmed blocking/major 3, refuted 0, downgraded 0.

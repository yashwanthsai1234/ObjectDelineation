# Verdict: final_r14_fresh-arch

## Skills used
- nvidia-skills:nvidia-skill-finder: loaded. No item depends on an NVIDIA product, so I did no catalog lookup.
- engineering:code-review: read `experiments/phase2/select_ckpt.py:93-95` (the only callers of the old samplers), `parity.py:5,71,89`, `r0_parity.py:7`, and the E4 step 1/3 text, E-G14 and U2.
- engineering:debug: traced both cases. Item 1: what stays on `main` after E under (b'). Item 2: what happens after a V1 crash before step 50 (step 3b exits non-zero, `a3_measure.yaml` still runs, then pull/push/teardown).
- engineering:system-design: checked item 2 against the (b') rule "Under (b') no rental follows V1" (00_front, "End state under (b')").
- engineering:architecture: compared the sampler rule with the trainer-path rule (tag `pre-e4-eval`). The tag already gives rollback, so the files on disk add nothing.
- engineering:documentation: every confirmed item quotes the source heading and the code line.
- huggingface-skills:huggingface-best: no model choice in question.
- Web tools loaded (WebSearch, WebFetch, Exa). Not used: both items are about repo text and code only.

## Verdicts
- fresh-arch-1: CONFIRMED (blocking). 00_front "End state under (b') (the one rule ...)" says "Under (b') no rental follows V1" and "no old code sits beside new code (owner archive rule)", and it lists E as landing on `main`. E4 step 1 says "Archive after the box check ... only after E-G10's box part passes (V4, E.8). Until then they stay". E-G14 lists `sample_vllm.py` and `hf_sample.py` as allowed exceptions until V4. 00_front 6.3 and 99_back decision 1 (b') both say "V4 waits for a budget raise". So under the default (b') the two scripts stay on `main` beside `policy/sampler.py` for good. After E4 step 1 rewrites `sampler_cmd`, they have no caller [code: select_ckpt.py:94-95 are their only live callers; parity.py and r0_parity.py only read sample_vllm output files or go through `select_ckpt.sample_view`]. That breaks the binding owner rules "Old code must be ARCHIVED" and "Remove dead code". No planned check catches it, because E-G14/U2 whitelist both files. The fix does not make things worse: rollback through the tag `pre-e4-eval` matches the plan's own trainer-path rule (E4 step 3). The new argument this round (no caller, so the in-tree fallback is unreachable anyway) answers the r9 downgrade's "box-tested fallback" reason.
- fresh-arch-2: DOWNGRADED (minor). After a crash before step 50, the engineer has nothing to guess. B1.5 step 3b says the chain stops. The launch lines still run `a3_measure.yaml` (items 4-5). The (b') rule says "Under (b') no rental follows V1" (00_front "End state under (b')"). Every rental also needs the 6.1 step 0 owner approval and `--owner-approved`. So the outcome is fixed: no B1 result, no retry spend without a new owner decision. The missing 1.1 row is a one-line clarity fix. `--save-steps 20` is a nice-to-have that shrinks the crash window.
- fresh-arch-3: minor as filed (no change).
- fresh-arch-4: minor as filed (no change).
- fresh-arch-5: minor as filed (no change).
- fresh-arch-6: minor as filed (no change).
- fresh-arch-7: minor as filed (no change).
- fresh-arch-8: minor as filed (no change).

Counts: confirmed 1, refuted 0, downgraded 1.

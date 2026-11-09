# Verdict on final_r14_rl.md (fresh skeptical verifier, 2026-11-09)

## Skills used
- nvidia-skills:nvidia-skill-finder: no NVIDIA catalog skill covers TRL GRPO / OPSM / Vast run plans; nothing changed.
- engineering:code-review: read REPO experiments/phase2/train_rl.py:42-68 (final/ saved only after trainer.train returns) and experiments/phase2/vast/policy.py:484-490 (`&&` join, `~` soft steps).
- engineering:debug: traced a crash at B3 step ~200 through b3_v5.yaml (no soft step, no --merge-last) and a LocGuard stop through Gate B3.
- engineering:system-design: compared the V1 (6.2) and V3 (6.3) launch lines side by side (pull job present in V1 only).
- engineering:architecture: checked B2.4 rule 3 against rule 4, the B2.2 preset literal, U38/U51/U68 for one source of truth.
- engineering:documentation: each verdict cites source file + heading and the code line.
- huggingface-skills:huggingface-best: no model-choice item in this review; nothing changed.
- Web tools loaded (WebSearch/WebFetch); no web fact was needed to judge these items.

## Verdicts
- rl-1, downgraded (minor): the Gate B3 after-run rules are held-out val-select rules (miss share = IoU < 0.1 must fall >= 25 % with CI upper bound < 0, both-found IoU and edge IoU CI bounds, legal >= 0.99; B.md "B.6.2 Run gates" Gate B3), so the reviewer's mechanism (a gamed policy lowers miss share "without raising IoU") is false by definition, a LocGuard-stopped model can pass only by being better on held-out data, and the missing "partial run" label is a consistency nicety, not a wrong result.
- rl-2, confirmed (major): B.md "B3.9 B3 run" fixes the order "0 CUDA build, 1 G08b smoke, 2 B1 re-eval, 3 `train_rl`, 4 B3 val eval, 5 G-RC0 (soft)" with a hard `train_rl` step and no `--merge-last`, and 00_front "6.3" V3 launch lines have no `launchctl load` of the pull job (V1 6.2 lines 722-726 do); with `jobs = " && ".join(...)` [code: policy.py:489] and `model.save_pretrained(a.out / "final")` only after `trainer.train` [code: train_rl.py:55-59], a crash or lost host late in the ~2 h B3 train leaves no B3 eval and no checkpoints off the box, and I22b asserts the current list, so no planned $0 check catches it (V3 money spent for nothing).
- rl-3, confirmed (major): B.md "B2.4 Decide" rule 3 says "raise the starting delta to the q99.9 of `seq_kl` ..." while rule 4 says "Until step 20, delta = 0.5 (TRL docstring)", the B2.2 preset is `"bf16-tis-opsm": {**_TIS, "off_policy_mask_threshold": 0.5}`, U38/U68 test 0.5 and U51 tests only the string `raise_delta_to == "q99.9"`; TRL copies the threshold once at init [code: .venv trl/trainer/grpo_trainer.py:882] and no flag carries a raised value, so two parts both claim to be the rule for steps 1-20 and every planned $0 test passes either reading.
- rl-4, minor (as filed): NaN-denominator gap is near zero today; not re-judged.
- rl-5, minor (as filed): per-request seed wording; not re-judged.
- rl-6, minor (as filed): seq_gaps duplicate blocks after a rare resume; not re-judged.
- rl-7, minor (as filed): review-log bookkeeping; review-log text is never blocking or major.

Counts: confirmed 2 (rl-2, rl-3), refuted 0, downgraded 1 (rl-1).

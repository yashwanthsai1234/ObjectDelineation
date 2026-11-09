# Verdict on final_r5_arch.md (fresh verifier)

## Skills used

| Skill | What it changed |
|---|---|
| nvidia-skills:nvidia-skill-finder | No item makes an NVIDIA product claim; no catalog skill was needed. |
| engineering:code-review | I opened every cited code line myself: preflight.py:49-56,101-251; policy.py:313-330,425-461,466-524; train_sft.py:74-105,149-157,257-276; select_ckpt.py:93-96; hf_sample.py / sample_vllm.py argparse; tests/test_p2_scripts.py:306-345. |
| engineering:debug | Traced post-E paths: spec -> preflight DATA checks (no data_root/view); SftStage(frozen) -> train_sft.resolve (assigns args.view); spec -> cmd_push_hf (KeyError output_dir) and vision_refs (top-level key, empty list). |
| engineering:system-design | Checked whether a planned post-E rental (V4) actually hits each seam; this set the severity of items 2, 4 and 5. |
| engineering:architecture | Checked the rule hierarchy (front order rules 1-3) to see if one part already wins over the other for item 1; none does. |
| engineering:documentation | Checked that each quoted doc line exists verbatim in the sources (00_front 6.3, 99_back 3b, E4 step 3, E6 steps 1/4/8, CLEAN8 step 1). |
| huggingface-skills:huggingface-best | No model choice involved; no change. |

## Verdicts

- arch-1: confirmed (blocking). 00_front.md 6.3 V4 paragraph ("Until V4 passes, E4 keeps the trainer eval path.") and 99_back.md 3b ("E4 keeps the trainer eval path and the old samplers") contradict E.md E4 step 3 ("`final_eval` ... no longer calls `build_trainer` ... The same change commit deletes its `--train` flag") and I22's AST assert; no order rule decides between them, and the trainer path is final_eval.py:73-75 today.
- arch-2: confirmed (major). E.md E6 step 8 says only "Data and view hashes come from the spec stages", but every DATA check reads `a.data_root`/`a.view` rows.jsonl (preflight.py:101-216, CHECKS :244-251) and no RlStage/EvalStage has those fields, so for any post-E spec (including V4's, which needs `make preflight` first) the engineer must invent which checks run on which root and how `{data}` resolves on the Mac.
- arch-3: confirmed (major). E6.1/E6.4 make `train_sft` build a frozen `SftStage` via CliApp.run, but `resolve()` assigns `args.view, args.output_dir` (train_sft.py:152) and `args.vision` (:155), `train()` assigns `args.branches` (:276) and `__post_init__` sets `self.vision` (:105); E6 never defines `to_argv(SftStage)` (today `--config <yaml>`, policy.py:443, and the yaml is archived), and U22 subtracts only `{kind,name}` while StageBase has `kind,name,soft` (E.md:805-808); the CPU dry run that preflight requires would crash, refusing every rental until fixed.
- arch-4: confirmed (major). E8.1 moves `cmd_push_hf` and `vision_refs` unchanged, yet `cmd_push_hf` reads top-level `output_dir` (policy.py:517) and `vision_refs` reads top-level `vision` (policy.py:460), so on a RunSpec push crashes with KeyError (skipping teardown in the front 6.2 `&&` line) and the G-vision guard (policy.py:472-475) is silently empty; `b1_launch.yaml` has no fate in E (0 hits) although E6.8 makes it unreadable by preflight; the cap-source part alone would be minor.
- arch-5: downgraded (minor). CLEAN8 step 1 already instructs "the same commit gives its stage the same override and a test", the sampler interface already accepts a params mapping with `max_tokens` (E.md I20 uses `max_tokens=1`), and select already samples in a child today (select_ckpt.py:93-96), so only the field name is left to choose; a missing cap shows up as a slow Mac test, not as a wrong result or GPU spend.
- arch-6: minor as filed (not contested; front 2.2/2.3 drift, E's move tables are the declared source by order rule 1).
- arch-7: minor as filed (not contested).
- arch-8: minor as filed (not contested).
- arch-9: minor as filed (not contested; one-word file name fix).
- arch-10: minor as filed (not contested).
- arch-11: minor as filed (not contested).

Counts: confirmed 4 (blocking 1, major 3), refuted 0, downgraded 1.

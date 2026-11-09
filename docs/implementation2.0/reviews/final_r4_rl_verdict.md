# Verdict on final_r4_rl.md (round 4, RL review)

## Skills used

| Skill | What it changed |
|---|---|
| nvidia-skills:nvidia-skill-finder | No catalog skill covers vLLM LoRA serving numerics; no verdict depends on it. |
| engineering:code-review | Read the repo lines each item cites: `r1_qwen3vl2b.yaml:9`, `train_sft.py:159-165, 292`, `policy.py:428-452`, `g_probe.py:88`. Confirmed items 2 and 3 from them. |
| engineering:debug | Traced vLLM 0.32.0 (local uv cache, `vllm-0.32.0+cpu`): `base_linear.py:221-233` and `config/compilation.py:663` (`cudagraph_specialize_lora = True`), `lora_shrink_op.py:329-333` (grid z = `num_active_loras`). This refutes item 4. |
| engineering:system-design | Checked which run each rule feeds (V1 B1, V2 C arms, V3 item 2a, B3, runs after Gate B3) for items 1, 5, 6. |
| engineering:architecture | ADR-0013 rule (b) in B2.4 "Scope" against the front 6.3 "Numerics before speed" rule (item 1). |
| engineering:documentation | Each verdict names the source file and heading and the proving line. |
| huggingface-skills:huggingface-best | Not used for any verdict (model list not in dispute). |

## Verdicts

- rl-1: CONFIRMED blocking. 00_front 2.5 ("so F is the expected B3 preset", line 384) and 6.3 ("Numerics before speed ... F arm (item 2b) runs before the C arms", lines 703-717; V2 item 2b 13 min / $0.15) and C header "Numerics before speed" (C.md:31-37), C.6 (C.md:1058), C.8.3 row (C.md:1230), 99_back decision 16 (line 154) contradict B2.4 "Scope" ("The B3 run keeps B1's preset `bf16-seqmask`", B.md:1152) and B.8.4 ("no F arm is needed in V2", B.md:1716); front 6.3 claims to win, so V2 pays $0.15 for an arm B says is not needed and measures C-G-SP1 on a push B3 does not run.
- rl-2: CONFIRMED blocking. Front section 3 step 0 says "only `model: qwen3-vl-4b` and `output_dir` changed" (00_front.md:439-442), but `r1_qwen3vl2b.yaml:9` keeps `reference: vision_ref/qwen3-vl-2b.json` and `train_sft.py:163-164` raises `ValueError` when the JSON's model differs, called from `preload` at `train_sft.py:292` before training; also `policy.py:450` runs `["train","select","probe"]` when `branches` is set, so select/probe read val (`policy.py:441-448`) and "val stays unread" is false in both 00_front and 99_back decision 17 (99_back.md:174-175).
- rl-3: CONFIRMED major. B2.1 says `log` "clears the list" (B.md:984-985) and RL logs every step (`g_probe.py:88` `"logging_steps": 1`), so the buffer holds one step, yet B2.4 rule 4 needs q99.5 over steps 1-20 "from the B2.1 buffer" and rules 3/4 pricing need per-sequence B1 values over steps 50-300 (B.md:1176-1189) that no task saves; the B1 run happens once in V1, so the input is lost and the engineer must guess.
- rl-4: REFUTED. The base side does not pay LoRA compute on decode: vLLM 0.32.0 defaults `cudagraph_specialize_lora = True`, documented as "separate cuda graphs for cases with and without active LoRA adapters ... remove this overhead" (`config/compilation.py:663-669`), and the shrink kernel grid z is `num_active_loras` (0 for a base-only batch, `lora_shrink_op.py:329-333`); the cold/warm prefix-cache point is symmetric under ABAB (one cold and one warm per side), so the ratio is not biased.
- rl-5: CONFIRMED major. B2.2 F step 4 pass rule compares against "the B1 run's (f) against (b)" (B.md:1090-1092; mirrored C.md:811), which only the B2.3b box part measures, and that part runs "only if B2.4 rule 0 does not already pick `lora-seqmask` or `fp16-seqmask`" (B.md:1129) / "only if B2.4 rule 0 did not pick F or A" (00_front.md:738); with the seed value 0.531 >= 0.2 rule 0 picks F, so the F check has no comparator.
- rl-6: CONFIRMED major. C7 is "run inside Phase B" and adds `max_completion_length: L` to "the B config" (C.md:867, 883), G-C is measured "inside Phase B" (C.md:1104; 00_front.md:96) and C11 step 5 arms "inside a Phase B run" (C.md:1020), but B3.9 lists every B3 change without L (B.md:1417-1418) and B1 runs at 1,536 (`configs.py:73`); an L applied in B3 rescales dr_grpo (grpo_trainer.py:3286) and turns 1,200-1,493-token answers into R = 0 rows, a confound B3.10 would credit to v5, and nothing tells the engineer which way to go.
- rl-7: minor (as filed; estimate only, not escalated).
- rl-8: minor (as filed).
- rl-9: minor (as filed; follows from rl-1).
- rl-10: minor (as filed; citation hygiene).

Counts: confirmed 5 (blocking 2, major 3), refuted 1, downgraded 0.

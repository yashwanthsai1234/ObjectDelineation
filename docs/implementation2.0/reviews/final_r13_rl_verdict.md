# Verdict on final_r13_rl.md (fresh verifier, 2026-11-09)

## Skills used

| Skill | What it changed |
|---|---|
| nvidia-skills:nvidia-skill-finder | No NVIDIA skill covers TRL IS weights; no effect on the verdict. |
| engineering:code-review | Read `src/objdel/train/fast.py:351-362` (`is_ratio`, sequence modes sum the gaps; mask mode zeroes only r > 3) and `:581-602` (`set_is_ratio`); TRL 1.17.1 `grpo_trainer.py:2603-2611`, `:2776-2824`, `:3268-3269`, `:3285-3286`. The code mechanism of item 1 is correct. |
| engineering:debug | Re-measured `runs/rl/train_rl/rl.json` (770 steps): median step-mean ratio 0.823, median mean length 102.5 (mean 105.5); min ratio = 0 in 98.05 % of steps; steps with a 1,536-token answer have median step-mean ratio 0.808 vs 0.826 without (one row in 256, so it says nothing per row); corr(step-mean ratio, mean length) = 0.013. No per-row weight of a truncated row exists in any log. |
| engineering:system-design | Checked what the gate decides: Gate B3 pass rules are outcome rules (miss share, paired IoU, edge IoU, legal share); clipped share is "report only"; the 10 % switch is "the owner may". |
| engineering:architecture | ADR-0014 credit sentence names the package that changed ("v5 + unmasked truncation on bf16-seqmask"), not a magnitude; it stays true whatever the weight. |
| engineering:documentation | One line per item with a one-sentence reason. |
| huggingface-skills:huggingface-best | No model-choice item in this review; not used further. |
| Web (WebSearch/WebFetch loaded) | Not needed: the item rests on local TRL code and local logs. |

## Verdicts

- rl-1 (MAJOR -> minor): downgraded. The code path is right (TRL `grpo_trainer.py:3269` `per_token_loss * inputs["importance_sampling_ratio"]`; `fast.py:356-357` sums the gaps per sequence), but the weight of a 1,536-token row is an extrapolation, not a measurement, and the quoted bound is invalid (the step mean 0.823 includes rows masked to 0 above 3.0, so Jensen gives no bound from log 0.823; runaway repetition text may also have far smaller per-token gaps than the 105-token mean). It changes neither the B3 run nor any Gate B3 pass rule (B.md B.6.2: clipped share is "report only"; B3.5: the 10 % switch is "the owner may"); the ADR-0014 sentence credits the package, not a mechanism size. The same point was raised as minor in rounds 4, 5, 7 and 11 (`final_r11_rl.md` item 2); no engineer must guess anything that changes an outcome.
- rl-2 (minor): not a blocking/major item; left as minor (the three round-12 RL minors are wording and report-unit fixes; `final_r12_rl.md` items 1-3 have no verdict file and change no result).

Counts: confirmed 0, refuted 0, downgraded 1.

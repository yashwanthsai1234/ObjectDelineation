# Verdict: final round 6, fresh-arch review

## Skills used

| Skill | What it changed |
|---|---|
| nvidia-skills:nvidia-skill-finder | No item depends on an NVIDIA tool. I did no catalog lookup and installed nothing. |
| engineering:code-review | I checked the code facts in the repo. `src/objdel/train/fast.py` and `src/objdel/reward/` exist. `advantages_v4` and `REWARD_V4` exist (reward/__init__.py:7, train/rewards.py:61). There is no v5 code yet. A grep of src/ and experiments/ finds none of `--gpu-pixels`, `--encoder-once`, `--host-meta`, `--fast-sync` or `--prompt-memo`. So only the plan decides what happens to these items. |
| engineering:debug | I traced the default (b') path. V1 runs. Then front step 5 lands the C1-C5 and B3 code. V2 and V3 do not run, so Gate C-adopt and Gate B3 are never read. So C10 and B3.10 never fire. E0's start rule waits for them. |
| engineering:system-design | I judged the end state the default leads to. No source defines it. |
| engineering:architecture | I checked whether a higher rule already settles this. Front order rule 1 says "no exception". E R-4 escape covers only a "fix". C.8.3's Phase E row says "If E5 lands before C, every C path maps through this table", which conflicts with order rule 1 and does not cover B3.10. None of these settles item 1. |
| engineering:documentation | I opened each quoted line in its source file and heading. Every quote matches word for word. |
| huggingface-skills:huggingface-best | No model-choice item. No change. |

The verdicts cover only the blocking and major items (item 1). The review rates items 2-8 as minor, and I did not re-grade them.

## Verdicts

- fresh-arch-1: CONFIRMED (blocking). The quotes match their sources word for word:
  - 00_front section 1, time line step 5: "C1-C5, B2/B3 and D2-D8 code land on the Mac, on pre-E paths".
  - Step 6: "Under the default (b') neither runs".
  - Step 7 and E.md "When E starts": "E0 starts after the last A-D commit that edits `fast.py`, `dataset.py`, `reward/` or `checker/`".
  - Order rule 1: "There is no exception".
  - C.md C10: "Files and actions (after each gate passes): delete from `src/objdel/train/fast.py` ...; remove the temporary flags ...".
  - C.md C.8.2: V2 "under the default (b') it does not run".
  - B.md B3.10 is gated on Gate B3, which needs V3.

  So under the default (b'), C10 and B3.10 never run. That leaves two failures. First, the E0 start condition can never be met, or E starts and these tasks break order rule 1. Second, reward v4 stays beside v5, and the CPU pixel path stays beside the GPU path behind temporary flags. No source plans their removal, which breaks the owner rule "Old code must be ARCHIVED (not left beside new code)". No source (00_front, 99_back decision 1, E, CLEAN) defines this end state: the grep for C10 and B3.10 finds no (b') handling. The fix must touch 00_front (end-state paragraph, order rule 1, step 7), E.md ("When E starts", R-4 or E.7 row), C.md (C10, C.8.3), B.md (B3.10, B.8.4) and 99_back (decision 1 option (b')).
- fresh-arch-2..8: not re-graded (the review already rates them minor).

Counts: confirmed 1, refuted 0, downgraded 0.

# Verdict on final_r7_data.md (fresh verifier, 2026-11-08)

## Skills used

| Skill | What it changed |
|---|---|
| nvidia-skills:nvidia-skill-finder | No NVIDIA product in the data scope. No catalog lookup, no install. No verdict depends on it. |
| engineering:code-review | Read `splits.deal`, `splits.lidc_stratum`, `splits.lidc_split` [code: src/objdel/data/splits.py:59-101] and how the split built `Patient` [code: experiments/phase1/build_data.py:123-127], against the B.9.8 sketch (data-1). |
| engineering:debug | Re-ran the halves with `python -I` on `data/eval_views/lidc-val-shifted/lidc/val.jsonl`: seed 0 gives 40 / 38 patients, 567 / 517 items, 14 / 30 disagreement items (fails (a) and (b)); seed `"val_half:78"` gives 533 / 551, 24 / 20, as B1.0 says. |
| engineering:system-design | Checked which text each test row reads: U02 seeds `"val_half:<s>"`; the sketch seeds 0; B1.0 says nothing about the `Patient` source. |
| engineering:architecture | Checked rule ownership: B1.0 owns the halves rule, U02 is its test spec, and the B.9.8 header tells the engineer to "keep the assertions". Two places claim the rule and disagree. |
| engineering:documentation | One line per item, with the doc line and the code line or measurement. |
| huggingface-skills:huggingface-best | No model choice in scope. No verdict uses it. |
| Web tools | Loaded (WebSearch, WebFetch). Not needed: every claim in data-1 can be checked against the repo code and data. |

## Proof for data-1

- B.md "B1.0", Val halves: "deal the val patients with `splits.deal(groups, (("shape", 0.5), ("select", 0.5)),
  seed=f"val_half:{s}")` inside each `splits.lidc_stratum` stratum ... Today's val gives s = 78: 40 / 38 patients,
  533 / 551 items, 24 / 20 disagreement items".
- B.md "B.9.1" U02: "seeds `"val_half:<s>"` from s = 0 ... an unknown `group_id` raises `KeyError`."
- B.md "B.9.8 Pytest sketches": header "keep the assertions, change the call sites if the real seams differ";
  `test_val_half_is_the_stratified_deal` has `want |= splits.deal(gids, (("shape", 0.5), ("select", 0.5)), seed=0)`
  and `assert {g: val_half(g) for g in want} == want`; `test_val_half_independent_of_hashseed` prints
  `val_half('lidc/LIDC-IDRI-0001')` and asserts only `len(outs) == 1`.
- [measured, 2026-11-08, `python -I`]: seed 0 per `lidc_stratum` (record-count `Patient`) gives 567 / 517 items
  (diff 50 > 43.4) and 14 / 30 disagreement items (diff 16 > 4). `val_half:78` gives 533 / 551 and 24 / 20.
- [measured]: group ids have no prefix (`"group_id": "LIDC-IDRI-0026"` in the val view); LIDC-IDRI-0001 has 8 train
  records and 0 val records. So the hash-seed test sees a `KeyError` (empty stdout) in all three runs and passes.
- [code: experiments/phase1/build_data.py:125-127]: the split's `Patient` = (XML slices, `disagree_5px` count);
  B1.0's s = 78 needs `Patient(items, disagreement items)` from the val view. B1.0 does not name the form.

## Verdicts

| id | verdict | reason |
|---|---|---|
| data-1 | confirmed (major) | The B.9.8 sketch asserts `val_half` equals the seed-0 deal, which contradicts B1.0 and U02 (`"val_half:<s>"`, s = 78) and fails B1.0's balance rules on the real val (567 / 517 items, 14 / 30 disagreement, re-measured); the hash-seed sketch passes vacuously on a `KeyError`; B1.0 omits the `Patient` form that its s = 78 needs. |
| data-2 | minor (as filed) | Filed minor; reading val pixels for a leakage gate changes no outcome, and preflight re-runs all gates. |
| data-3 | minor (as filed) | Filed minor; a test-text gap on padded windows, no wrong training result. |
| data-4 | minor (as filed) | Filed minor; I20 test coverage gap, the B1.4 rule text is correct. |
| data-5 | minor (as filed) | Filed minor; changes only Gate B1 report wording. |
| data-6 | minor (as filed) | Filed minor; affects one branch of an E-phase TOST report. |
| data-7 | minor (as filed) | Filed minor; duplicate materialization of reused windows, a test-wording issue. |
| data-8 | minor (as filed) | Carried round-6 minors; wording only. |

Counts: confirmed 1 (major), refuted 0, downgraded 0.

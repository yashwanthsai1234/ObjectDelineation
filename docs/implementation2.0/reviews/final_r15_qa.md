# Final review round 15: QA and technical writing (Tests sections, gates, testability, cross-section consistency)

## Skills used

| Skill | What it changed |
|---|---|
| nvidia-skills:nvidia-skill-finder | This round makes no choice of NVIDIA product, so I recommend no catalog skill. The skill made me re-read the box toolchain lines (front 6.4 tail, D.4 "Common toolchain"). They agree with each other. No finding. |
| engineering:code-review | I read the code that the round 14 changes cite. `experiments/phase2/train_rl.py:20-62`: `--save-steps` default 50; `final/` is written only after `trainer.train`. `experiments/phase2/vast/policy.py:420-456` (`chain`), `:470-510` (`job`, the `&&` join, the `~` soft step, `cmd_pull`), `:527-536` (`cmd_fetch_hf` writes `{rd}/ckpt`) and `:61, 73, 80` (`UP_GB_PER_H`, `PULL_S`, `PULL_EXCLUDE`). `experiments/phase2/vast/bundle.py:20-70` (the `eval_views` glob). `experiments/phase2/select_ckpt.py:93-96` and `hf_rate.py:1-30`. `src/objdel/data/lidc.py:316` and `riga.py:332` (item id prefixes). The TRL lines `.venv trl/trainer/grpo_trainer.py:882, 3188-3204` and `grpo_config.py:346-353`. |
| engineering:debug | I traced three paths step by step. (1) The A3 header commands on the Mac: build `lidc-test-shifted`, then bundle, then front 6.1 step 4c (item 1). (2) When `make_d_goldens.py` runs, compared with when D-R01 needs its "pre-D" text (item 2). (3) A B3 run under `bf16-seqmask`: the `rl.json` key `opsm_start_delta` (item 3). |
| engineering:system-design | One owner per fact. The V1 launch lines have one owner (front 6.2). The A3 header keeps a second, shorter copy (item 1). |
| engineering:architecture | Licence boundary: D-R01's new golden is LIDC only (CC BY 3.0). Real item ids start with `lidc/` and `riga/` [code: src/objdel/data/lidc.py:316; src/objdel/data/riga.py:332], so the round 14 prefix assert works. No finding. |
| engineering:documentation | Carry-over check: I re-read the round 13 and round 14 qa minors. I checked each one with grep. The open ones are listed at the end. |
| huggingface-skills:huggingface-best | This round makes no model choice. No change. |
| engineering:testing-strategy | For each round 14 test change (D-R01, D-R02, I22b, I49, U38, U45, U51, U68, E-G14), I asked three questions. Can it fail? Is its expected value independent of the code under test? Does a named gate run it? D-R01's always-run copy can pass by construction (item 2). |
| mattpocock-skills:tdd | A golden must come from a source that is independent of the code under test. D-R01's always-run copy does not name such a source (item 2). U45's "absent gives 0.5" checks a value that means nothing under a preset without OPSM (item 3). |

Web tools: loaded (WebSearch, WebFetch, Exa). No external claim in the round 14 changes needed a new web check. The
TRL claims were checked against the installed package (lines above). The LIDC licence check of round 14 still holds.

Measurements this round ($0, Mac, 2026-11-09):
- `docs/implementation2.0/Implementation.md` is 1,760,475 bytes and `review-log.md` is 380,747 bytes (`ls -l`). Both
  are under the DOC-G13 limits (2 MB for the plan, 1 MB for every other file).
- No match of the check 14 pattern (`/[p]rivate/tmp`, `/[U]sers/`) in `Implementation.md` or `review-log.md`
  (`grep -cE`: 0 and 0).
- `data/eval_views/` holds today `lidc-val-rescaled`, `lidc-val-rotated`, `lidc-val-shifted`, `probe-precision` and
  `probe-topology` (`ls`). There is no `lidc-test-*` folder yet. A3's header creates one (item 1).
- These old-sampler references are outside `archive/`: `select_ckpt.py:93-96` (the only live caller),
  `hf_rate.py` (docstring only; it calls through `select_ckpt`), `parity.py:5,71,89` and `r0_parity.py`. In
  `tests/`, the six files that E4 step 1 names, plus `tests/test_profile.py:113`. That file only asserts that the
  names are absent, so it keeps passing after the archive (`grep -rlE "sample_vllm|hf_sample"`). E4 step 1's list
  is complete for the tests that break.

Round 14 changes checked, no finding:
- B3.9 steps 3 (soft, `--save-steps 50`) and 3b (`--merge-last`):
  - `job()` runs `~` steps as `{ cmd || echo ...; }` and joins every step with `&&` [code: policy.py:484-490].
  - `train_rl` has `--save-steps` (default 50) and writes `final/` only after `trainer.train` [code:
    train_rl.py:30, 55-60].
  - I22b, I49 and the B.7 crash row agree with the order line "0, 1, 2, 3 (soft), 3b, 4, 5".
- V3 money:
  - Item 4: 3 + 5 + 13.0 + 6 + 119.17 + 10.4 = 156.6 min (139.9 at 23.65 s).
  - The B-S2 row: 195.6 min = $2.32.
  - `BudgetStop`: $0.10 + 185.2 x $0.01133 = $2.20; the amount is $2.55.
  - The four-item row: 255.2 min = $2.99, `BudgetStop` $2.88, amount $3.23.
  - All of these agree with front 6.4 rule 5.
- The OPSM start delta (B2.2, B2.4 rules 3-4, U38, U45, U51, U68, the B.9.8 sketch):
  - TRL copies `off_policy_mask_threshold` once at init [code: grpo_trainer.py:882] and reads the attribute at
    each loss call [code: grpo_trainer.py:3190-3203].
  - The docstring example is 0.5 [code: grpo_config.py:352].
- E4 step 1 "Archive in the switch commit", E-G14 and CLEAN8 step 4 agree. `hf_rate.py` runs through
  `select_ckpt.sample_view` [code: hf_rate.py:6, 31], so it keeps working on the new sampler until CLEAN8 archives it.
- B.6.1 has 18 numbered gates plus DG1b, DG9b and DG14b. That is 21, the count in front 6.1 step 2.
- DOC-G7 and E-G1 baselines agree: 1,447 passed + 6 skipped + 31 deselected = 1,484 collected.

## Items

### 1. minor: A3's header commands build `lidc-test-shifted` and then bundle, so the test view goes into the V1 bundle; front 6.1 step 4c then fails (A, 00_front)

- Doc (A.md, A3, "Header comments, in this order"). Three lines in this order:
  - `build_eval_views --split test --kinds shifted --out data/eval_views --procs 6`
  - `final_eval --hash-view data/eval_views/lidc-test-shifted --out runs/a3/pre_rental.json`
  - `.venv/bin/python -I -m experiments.phase2.vast.bundle lidc --view lidc-bundle-v3`
  No line moves the test view out between them.
- Code: the bundle takes every `data/eval_views/lidc-*` folder: `{p.name for p in (data / kind).glob(f"{ds}-*")}`
  [code: experiments/phase2/vast/bundle.py:54]. So `lidc-test-shifted` goes into the V1 bundle.
- Doc (00_front 6.1 step 4c): "The V1 bundle holds ... no `lidc-test-*` folder ... so the test view is moved out
  before bundling (or a bundle `--view` filter is added)". Front 6.4 rule 7 states the same rule.
- So the two parts disagree:
  - A3's header is the ordered command list. I20 asserts it "in that order".
  - That list has no move step, and 6.1 step 4c offers two ways and picks neither.
  - V3 item 3 needs the test view back in the V3 bundle, and no line puts it back.
- The A3 header also keeps a short copy of the V1 launch lines: `sync_up.sh <run> && fetch_hf.sh <run>
  rl-full-3ep`, then `run.sh ... a3_measure.yaml && pull.sh <run>`. That copy has no `--rows` fetch, no pull job,
  no `b1_rows.yaml` run and no `push_hf.sh`. Front 6.2 says it "holds the only copy".
- Severity: 6.1 step 4c is a $0 check before V1, so the gap is caught before any spend.
- Fix:
  - A3 header, after the `--hash-view` line: `mkdir -p data/held_views && mv data/eval_views/lidc-test-shifted
    data/held_views/   # 6.1 step 4c; V3 item 3 moves it back before its bundle`. The content hash does not depend
    on the folder path, so `pre_rental.json` stays valid.
  - 6.1 step 4c: drop "(or a bundle `--view` filter is added)", so there is one way.
  - 6.3, before the V3 `sync_up` line: "If item 3 is approved: `mv data/held_views/lidc-test-shifted
    data/eval_views/` and re-bundle."
  - A3 header: replace the `sync_up` / `run.sh` / `teardown` lines with
    `# V1 launch lines: front 6.2 (the only copy); a3_measure.yaml is its second run.sh call`.
  - I20: change "end with `teardown.sh` and `sweep_orphans.sh`" so that it matches the new last line.

### 2. minor: D-R01's always-run golden holds "the stored pre-D text", but the script that writes it runs on the post-D build and names no pre-D source, so the test can pass by construction (D)

- Doc (D.md D.9.4, D-R01): "Always-run copy (LIDC only): 50 fixed LIDC items ... in `tests/goldens/d_r01_lidc.json`
  (input masks as run-length lists plus the stored pre-D text), written on the Mac by `make_d_goldens.py` (D8);
  rebuilt text equals the stored text."
- Doc (D.md D8, payload item 3): `make_d_goldens.py` writes GPU02's golden from "the Mac CPU build". GPU02 covers
  the entries of D2-D7b. Then: "The same script writes D-R01's always-run golden".
- Doc (D.4 "Order"): "D8 grows with each task".
- So the script runs after D2-D7b exist. Nothing says that its D-R01 part reads text made before D3. Suppose the
  engineer builds the "stored" text by calling today's `build_target` (the post-D3 engine). Then the always-run test
  compares the post-D code with itself and cannot fail.
- D-R02 names its source ("the golden from the archived pre-D code"); D-R01 does not.
- Severity: D-R01's `data` half (every `target_set_sha256` of the 7 manifests, rebuilt on the Mac) and D-G4 still
  catch a text change. Only the always-run copy, which CI runs, loses its value.
- Fix (D-R01 and D8 payload item 3): "The stored text is the `targets` field of the 50 items in
  `prepared/lidc/val.jsonl` as built before D3. B1.4a's hash check proved that field byte-identical to the numpy
  polygonizer. `make_d_goldens.py` copies it and never calls `build_target`. The golden records the SHA-256 of the
  source `val.jsonl`. The test asserts that the golden's text equals no fresh output that the script itself wrote."
  A simpler form: "the golden is written once, in the commit before D3, and later commits never rewrite it."

### 3. minor: `rl.json` records `opsm_start_delta` 0.5 when the flag is absent, also under presets with OPSM off (B1, B3), so the record names a delta that was never used (B)

- Doc (B.md B2.2): "`rl.json` records `opsm_start_delta` (0.5 when the flag is absent)". U45: "absent gives 0.5".
- B3 runs `--numerics bf16-seqmask` (B3.9). Under that preset `off_policy_mask_threshold` is not set, because
  `_SEQMASK` has no such key (B2.2 `NUMERICS`). So the B3 `rl.json` says `opsm_start_delta: 0.5` while OPSM was off.
- U51 already returns `start_delta` `None` for every preset other than `bf16-tis-opsm`. So the decision function
  and the run record disagree.
- Fix:
  - B2.2: "`rl.json` records `opsm_start_delta`: the flag value, or 0.5 when the flag is absent, under
    `bf16-tis-opsm`; `null` under every other preset."
  - U45: add "`--numerics bf16-seqmask` without the flag gives `opsm_start_delta` `null`".

## Carry-overs still open (minor; checked with grep this round; not re-argued)

- r14 qa-2: A2, A8, C2 and C11 still read `runs/a3/eval_val/answers.jsonl`. No source has `OBJDEL_V1_RUN`
  (A.md:428, 1109; C.md:651, 1093).
- r14 qa-3: E.md still does not name `tests/variants.py` in E8's list of hard-coded test paths.
- r14 qa-4: D5c still sizes `in` by "names per answer" (D.md:685). There is no `NAMES_MAX` and no
  `pinned_bounds`.
- r14 qa-5: `GC_ALIAS_SCAN_FOR_TEST` still does not name the aliased `scan` call (D.md:469).
- r13 qa-3: the B1.3 DG1b row still says "(box: the manifest's `excluded_ids`)" (B.md:1934).
- r13 qa-4: C.9.4 still has no `limit_mm_per_prompt={"image": 0}` (no match in C.md).
- r13 qa-5: C.md:1968 still says "15.4 min < 30.9 min".
- r13 qa-7: the C10 keep list still omits `--encode-check` and `--checks-from-step` (C.md, C10 "Keep ...").
- r13 qa-8: E.md:645 and E.md:747 still say "moved by D3". The CLEAN-G14 row (CLEAN.md:1052) still has the D3
  wording.
- r13 qa-9: 00_front.md line 4 still says "Revision: final round 3 plus the round-4 consistency pass".

## Totals

Blocking 0, major 0, minor 3 (new), plus 10 carry-over minors.

# Final review round 8: fresh outside architect

## Skills used

| Skill | What it changed |
|---|---|
| nvidia-skills:nvidia-skill-finder | No item needs an NVIDIA catalog skill. Front section 4 already lists `nv-segment-ct` and `dicom-series-to-volume` for the 3D work. I did no install and no catalog lookup. Item 4 asks for one CUDA run on sm_120, which is a plan question, not a tool question. |
| engineering:code-review | I read the code that the V1 path uses: `train_rl.py:19-68`, `final_eval.py:1-80` (defaults `--k 8`, `--seed 11`), `bundle.py:26-75`, `policy.py:58, 188-193, 315-330, 469-492, 601-606`, `configs.py:170-195`, `dataset.py:85-260` and TRL 1.17.1 `grpo_config.py:939`. Every number the plan states for these lines matched. `run.sh` defaults to `--cap-hours` 24, so the per-call cap is the deadline minus `MARGIN_S`, as the plan says. |
| engineering:debug | I followed V1 from the Mac to the box, one command at a time. This found item 1: A3's command list builds the test view into `data/eval_views/` and then bundles, but `bundle.members` ships every `lidc-*` eval view. No task owns the step that front 6.1 step 4c needs. |
| engineering:system-design | I checked whether the default end state is coherent end to end: data -> V1 -> Gate B1 -> D -> E -> CLEAN. It is, with one gap. Under (b') the CUDA engine changes of D2-D8 never run on a GPU (item 4). |
| engineering:architecture | I compared the plan with what a top lab would do: budget use, the long-lived branch, on-the-fly views, and numerics before data. These are recommendations (items 4-7), each with its trade-off. None is a defect in the plan, so each is minor. |
| engineering:documentation | I checked the front's target tree and diagrams against the phases that own the paths. This found items 2, 3, 8 and 9. |
| engineering:tech-debt | I checked CLEAN wave 1 against D3. CLEAN5 archives `experiments/phase0/` before V1, but D3 and D's S03 still treat those scripts as live callers (item 3). |
| huggingface-skills:huggingface-best | No Hub leaderboard covers polygon-text segmentation. The model choice (Qwen3-VL-2B now, 4B gate on val) is sound. I have no model item. |
| system-design-skills:system-design | Trade-off method: for each recommendation, what it solves, what it makes worse, and when to change. Items 4-7 use this form. |
| WebFetch / WebSearch | PyPI JSON for import-linter: 2.16 was uploaded on 2026-11-06, as E.2 says. A search on vLLM LoRA RL numerics found no 2026 fix that would change B2's plan. Sources: https://pypi.org/pypi/import-linter/json ; https://swift.readthedocs.io/en/v4.3/_sources/Instruction/GRPO/AdvancedResearch/training_inference_mismatch.md.txt ; https://verl.readthedocs.io/en/v0.4.1/advance/ppo_lora.html |

## Verdict on the architecture (short)

For a $3 budget the shape is right. The plan does data first: one paid run, with a decision table approved before the result. The CUDA engine has one source (CPU build as the oracle). There is one TRL adapter with a start-up check, a typed run spec, and an exhaustive import contract. Under the default (b'), no temporary code ends up on `main`. A top lab would keep this core. The main differences are in where effort goes (items 4-7). I found no blocking item. I found one major item (item 1), on the only paid path.

## Items

### 1. major: no task owns keeping the test view out of the V1 bundle; A says the opposite (A, 00_front)

Problem.
- 00_front, "6.1 Before any rental", step 4c: "The V1 bundle holds the frozen `lidc-val-edge` view (B1.6) and `lidc-val-shifted`, and no `lidc-test-*` folder. ... so the test view is moved out before bundling (or a bundle `--view` filter is added)." This gives two options, but no command, no task and no place to move the view to.
- A, "A3. Save every answer ...", the header command list builds the test view and then bundles, with no step in between:
  `build_eval_views --split test --kinds shifted --out data/eval_views ...`, then `final_eval --hash-view data/eval_views/lidc-test-shifted ...`, then `experiments.phase2.vast.bundle lidc --view lidc-bundle-v3`. The A3 text then says: "`lidc-test-shifted` reaches the box when it exists on the Mac at bundle time."
- A, "A.3.1 Boxes and data flow": the arrow to the V1 box is labelled "bundle: lidc-test-shifted, lidc-val-edge, val".
- 00_front, order rule 2 (the `v1-b1` tag list) names no bundle filter. B.8.1 no longer has the bundle note that the B review log cites ("B.8.1 bundle note (no `lidc-test-*` folder ...)"). A grep of B.md section B.8.1 for "bundle" finds 0 hits.
- Code: `bundle.py:51` keeps every `data/eval_views/lidc-*` folder: `names = ... {p.name for p in (data / kind).glob(f"{ds}-*")}`. Only `views` has a "train only when named" rule (`bundle.py:52-53`). `eval_views` has none.

Effect. An engineer who follows A3 exactly fails 6.1 step 4c on the V1 tag. Two choices are then open, and they lead to different results:
- Move the folder. The plan names no target place. B1.6's DG15 check and V3's `a3_test.yaml` read `{data}/eval_views/lidc-test-shifted`, so the folder must also come back before V3.
- Write a filter. The filter must be on the tagged tree, and no task in the tag list owns it.

If the engineer picks neither, the test view goes to every rental. That is against 6.4 rule 7.

Fix.
- A3 owns one code change on the `v1-b1` tag. In `bundle.members`, an eval view named `{ds}-test-*` ships only when it is named. Copy the existing train-view rule (`bundle.py:52-53`). Add a `--eval-view NAME` flag (repeatable) to `bundle.py`.
- Add a unit test in `tests/test_vast.py`. Build a temp tree with `eval_views/lidc-test-shifted`. With no flag, the bundle has no `eval_views/lidc-test-*` member. With `--eval-view lidc-test-shifted`, the bundle has it.
- 00_front 6.1 step 4c: replace "moved out before bundling (or a bundle `--view` filter is added)" with "`bundle.py` ships a `lidc-test-*` view only when `--eval-view` names it (A3)".
- 6.4 rule 7: replace "only a `lidc-test-*` view is moved out before a rental whose table does not name it" with "a `lidc-test-*` view ships only in a rental whose table names it (V3 items 3 and 5: `bundle lidc --view lidc-bundle-v3 --eval-view lidc-test-shifted`)".
- Add the filter to order rule 2's tag list.
- A3 header and A.3.1: drop "lidc-test-shifted" from the V1 bundle arrow. Change the sentence "reaches the box when it exists on the Mac at bundle time" to "reaches the box only when the bundle names it (V3)".

### 2. minor: E points to the wrong front rule for its start condition (E)

E.md, "Order of execution", "When E starts": "Phase E runs after A-D and after CLEAN wave 1, as front section 1 order rule 2 says: E0 starts after the last commit on `main` ...". Front order rule 2 is the `v1-b1` tag list. The E start rule is front time line step 7 and "End state under (b')". Fix: "as front section 1 time line step 7 and 'End state under (b')' say".

### 3. minor: D3 and D's S03 treat archived phase-0 and lab scripts as live callers (D, CLEAN)

- D.md, D3 "Files": "`experiments/phase2/lab/B3/crop_trace.py:43` and `experiments/phase0/{lidc_census,riga_census,riga_outlines}.py` import `objdel.checker.polygonize`/`certify` ...; they must still import". D3 also extends the `--help` check to `experiments/phase0/*.py`.
- D.9.6 S03 runs `--help` for `experiments/phase0/*.py`.
- But CLEAN.md, CLEAN5 (Wave 1) moves "`experiments/phase0/*.py` (7 scripts + `__init__.py`) -> `archive/phase0/`" and "`experiments/phase2/lab/` -> `archive/phase2_lab/`". Wave 1 is on the `v1-b1` tag (front order rule 2), and D2-D8 land after V1 (front time line step 5).
- So at D3 these paths no longer exist. Archived code "is never imported" (front section 0, "archive"). The S03 glob is then empty, and the test passes without checking anything.

Fix. In D3, delete the clause about `crop_trace.py` and phase0 "must still import". Drop `experiments/phase0/*.py` from the D3 `--help` extension and from S03. Add one line: "CLEAN5 archives phase 0 and the lab trees before D3 (Wave 1), so D3 has no caller there."

### 4. minor (architecture recommendation): under (b') the D2-D8 engine port never runs on a GPU; V1 is the only window for R-D1

- 00_front 6.3, V2 item 4: R-D1 "Runs only if D2 has not merged on the V2 tree (D1 step 4); on the time-line order it is dropped".
- "End state under (b')": D1-D9 land "(D's box gates stay 'open until a rental')".
- Under (b') no rental follows V1. So the CUDA polygonizer (D2/D3) and the scorer changes (D5, D5c pool memory) land in production code with zero GPU runs, and R-D2 (V3 item 1, also "Skippable") never runs. The owner's rule names CUDA for exactly these parts.

Trade-off. V1 already pays the launch (about 53 min, $0.70). After item 5, the K1 runbook alone (12-15 min, about $0.17 at r_s) would give the first sm_120 exactness check. V1 is before D2 merges, so the D1 tree rule allows it. Cost: V1 then holds D1 step 0 (runbook paths, the `code_files()` keep rule, `r_d1.yaml`), which order rule 2 now excludes ("no ... kernel code"). Money: at 27.5 s, V1 plus reserve is $2.54 of $3, which leaves $0.46. That is enough only when the B1 run is near 23.65 s per step.

What a top lab would do. Add a skippable V1 item 6, "K1 GPU_RUN steps 1-5". It runs only if session spend after item 5 plus the item's cost is at most the approved $2.40. This needs the owner to approve it with the V1 item list.

When to change. If the owner raises the budget before D2 merges, keep R-D1 in V2 as planned.

### 5. minor (architecture recommendation): do not write C1-C11 code on `after-raise` before a raise is decided

- 00_front, "End state under (b')": "C1-C11 code and their temporary flags, B2.2-B2.4 ... and B3.1-B3.10" wait on the owner branch, which "never merges under (b')".
- E then splits `fast.py` symbol by symbol into `policy/caches.py`, `rl/forward.py` and `rl/trl_adapter.py` (E5 move map). R-4 calls the later rebase "a mechanical change". For a branch that edits the bodies of `FastGRPO` and `packed_chunk`, the rebase becomes a hand port per symbol.

What a top lab would do. Freeze C at design plus C0's measurement tools. Write C1-C11 only after the owner funds V2. Then write it on post-E paths, which removes the rebase. This keeps the B3 code (reward v5, $0 lab round), because it decides a paid run.

Trade-off. If a raise comes soon, C then starts later. If no raise comes, many days of Mac work and one rebase across E5 are saved.

### 6. minor (architecture recommendation): RL rows do not need polygon-text targets; the ADR-0015 exception is only for DG12

- B.md B.3.3: "**Owner exception needed:** the targets use the numpy polygonizer, against the rule 'CUDA for the polygonizer'."
- Code: GRPO reads only the prompt, the image and the doctor caches (`dataset.prompts` builds `{prompt, image, key, W, H, names}` [code: dataset.py:217-222]). `targets` is used by `sft_rows` and by DG12 only.

Option. Build the RL rows with `targets: null` until D3 lands. Then run DG12 on the Mac after D3, with the engine CPU build, as a report on the same windows. This removes the owner exception from V1's data path.

Trade-off. DG12 moves from "before V1" to "after D3". This needs the owner's OK, because DG12 is a hard data gate today. Keep the exception if the owner wants every data gate before V1.

### 7. minor (architecture recommendation): numerics before data at the next raise

B.1 seed value: merge relative error 0.531 at step 100. FR4: IS-ratio minimum is 0 in 98.1 % of 770 steps. B1 and B3 both keep `bf16-seqmask` with the merged push (00_front 2.5). This keeps Gate B1 fair, because R5 had the same numerics.

What a top lab would do. A run with masked sequences wastes rollouts. Most RL stacks that use vLLM rollouts fix the trainer-rollout mismatch before they tune data or reward. The next paid run after V1 should take preset F (adapter serving) as its first change, before B3, rather than after Gate B3. The plan already orders "F check in the first 15 steps of the first run after Gate B3". The suggestion is to move that check into the B3 run when a raise funds V3.

Trade-off. B3 then changes two things at once (reward and numerics), so the B1-vs-B3 split weakens. Keep the plan's order if attribution matters more than sample use.

### 8. minor: front target tree and diagrams drift from the owning phases (00_front)

- 00_front "2.3 Target tree": `core/ spec.py`. E.3.1 adds `provenance.py`.
- `policy/` has no `callbacks.py`. E5 creates `policy/callbacks.py` (`TraceStep`).
- `data/` has no `audit.py`. E.3.1 lists "audit (hack_audit)".
- `mkdocs.yml ... (DOC10)`: DOC3 creates `mkdocs.yml`. DOC10 is the ADR process.
- "2.4 GPU-everywhere data path": "B1 rows root (9,600 rows)". Section 0 "rows root" says "8,320 rows = 260 steps x 32". B1.2 calls 9,600 "the 300-step upper bound".

Fix.
- Add `provenance.py`, `policy/callbacks.py` and `data/audit.py` to 2.3.
- Change "(DOC10)" to "(DOC3)".
- Change "9,600 rows" to "8,320 rows".
- Better: say once in 2.3 that E's move tables win over this tree.

### 9. minor: the date line under the title is stale (00_front)

00_front, first lines: "Revision: final round 3 plus the round-4 consistency pass." The doc is now at round 8 (review-log.md has rounds 5-7). Fix: "Revision: final round 8. The review logs are in `docs/implementation2.0/review-log.md`."

## Counts

blocking 0, major 1 (item 1), minor 8 (items 2-9; items 4-7 are architecture recommendations).

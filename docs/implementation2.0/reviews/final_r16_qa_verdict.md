# Verdict on final round 16 QA review (final_r16_qa.md)

## Skills used
- nvidia-skills:nvidia-skill-finder: no NVIDIA product choice in these items; no catalog lookup; nothing changed.
- engineering:code-review: read today's `Makefile:19-26` (no `test-data` target; `test` is `-m "not slow"`) and the existing `slow`+`data` tests (`tests/test_cuda_checker.py:155-156` `lidc_full`, `tests/test_data_splits.py:230-231`).
- engineering:debug: traced (1) A9's `make headline TEST=... VAL=...` through DOC6, F3 and I05; (2) D-R08 (`slow`,`data`) through A.9's targets and the D-G3 coverage row; checked what else would catch a real-data polygonizer defect (D-R01).
- engineering:system-design: one-owner check on the `headline` and `test-data` definitions; both have stale restatements, but the owning text (DOC6; A.9, with D8's "only if no earlier phase did") is explicit.
- engineering:architecture: checked which outcome each wrong reading produces; neither path spends GPU money or silently passes a wrong result.
- engineering:documentation: verified each quoted line in the sources (A.md:328-330, A.md:1161, DOC.md:159, DOC.md:674-676, DOC.md:1172, DOC.md:1177, DOC.md:1207; A.md:1331-1332, D.md:899-903, D.md:1336, D.md:1407, CLEAN.md:269, B.md:2293).
- huggingface-skills:huggingface-best: no model choice in scope; nothing changed.

## Verdicts
- qa-1 | downgraded (minor) | The contradiction is real (F3 `case "$(RUN)" in ""...exit 2` and I05 "without RUN ... exits non-zero" against DOC6's recipe "`$(if $(RUN),--run $(RUN),--test $(TEST) --val $(VAL))`" and A9's command). But DOC6 states the TEST/VAL form outright as "Phase A's patient-CI headline". Every wrong reading also fails loudly at $0 on the Mac: A9's own named command exits 2 or hits argparse on the missing `--model-ref`, `RUN=<a3 run>` makes headline.py refuse the missing file, and I10 goes red once the patient file is tracked. Nothing here rents a GPU, so the caught-before-spend rule applies. Still worth fixing: F3/I05 should check RUN only when it is set, and MODEL_REF should be added to the A.md commands.
- qa-2 | downgraded (minor) | The stale restatements are real: B.9.0, D8 and CLEAN.4 say `make test-data` is `-m data`, while A.9 (A lands first, and D8 "adds each change only if no earlier phase did") says `-m "data and not slow"`. D-R08 does run only under `test-all-data`. But skipping D-R08 gives no wrong result. D-R01 (`data`, not `slow`, so `make test-data` runs it) rebuilds every target set of the 7 manifests and needs `target_set_sha256` byte-identical after D3, D4 and D4b. P01/P02 check the polygonizer against the oracle and cv2 on fuzz masks. So a polygonizer defect that would change training targets is still caught at $0 before R-D2. The D-G3 Mac half also names test files that "pass on the Mac CPU build", and running those files directly does run D-R08.
- qa-3 | minor (kept) | Confirmed: CLEAN.md:1056 greps `src/objdel/train/rewards.py`, which E5 moves, so the grep fails visibly at $0.
- qa-4 | minor (kept) | P03's 1<<12 `eval_sets` cap and its peak assertion are underspecified. The test is Mac-only and fails visibly at $0.
- qa-5 | minor (kept) | D.md:250 defines `==` while D-G14/S13/B04 use `<=`. I19 (Mac) already asserts `==`, so an undercounting counter is caught at $0.
- qa-6 | minor (kept) | 00_front.md:738 says 3.39 h; 2.30/0.680 = 3.382 [arithmetic]. The cap 3.35 h is unaffected.
- carry-overs (13) | minor (not re-argued) | The review lists them as minor and they stay minor.

Totals: confirmed blocking/major 0; refuted 0; downgraded 2.

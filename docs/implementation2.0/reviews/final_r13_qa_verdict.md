# Verdict on final round 13 QA review (final_r13_qa.md)

## Skills used

| Skill | What it changed |
|---|---|
| nvidia-skills:nvidia-skill-finder | No NVIDIA product choice in scope; no catalog skill recommended. Kept focus on how the CPU test variants of the CUDA engine load. |
| engineering:code-review | Read D.md D2 step 7, D.9.9 item 1, sketch S3, D-G3, D-G12, D-G16, D8 Makefile list; read the editable finder `.venv/lib/python3.12/site-packages/_editable_skbc_objdel.py` (meta_path insert at index 0, `known_wheel_files` checked first in `find_spec`, lines 506-515, 734-742). |
| engineering:debug | Re-ran the variant-load probe myself (scratch `reviews/r13v/mut/objdel/_gpucheck.py`, run from `/`): `PYTHONPATH=mut python`, `python -I` + `sys.path.insert(0, mut)`, and `PYTHONPATH=mut python -I` all import `.venv/lib/python3.12/site-packages/objdel/_gpucheck.abi3.so`. |
| engineering:system-design | Checked whether any planned $0 gate catches an aliased scan map before R-D2 (none does). |
| engineering:architecture | One owner for the variant-build mechanism: D2 step 7 / D8; noted the simpler alternative already used in CLEAN1 step 5 (`-S` plus `sys.path` set in the child). |
| engineering:documentation | Checked that `make test-order` is named only in D.md D2 step 7 and D.9.9 item 1, defined nowhere, gated nowhere. |
| huggingface-skills:huggingface-best | No model choice in scope; no change. |

Web tools loaded (WebSearch, WebFetch); no web check was needed (the item is local code behaviour, measured).

## Verdicts

Only item 1 is blocking/major in the review; items 2-9 are already rated minor by the reviewer and are not re-rated here.

- qa-1 (major): CONFIRMED major. D.md D2 step 7 makes `make test-order` run the `K1_ORDER` and `GC_POISON_SCAN` builds, and D.9.9 item 1 / sketch S3 load a variant by putting its directory first on `PYTHONPATH` with `-I`; but `-I` ignores `PYTHONPATH`, and even without `-I` (or with `sys.path.insert(0, ...)`) the editable finder at `sys.meta_path[0]` maps `objdel._gpucheck` to `site-packages/objdel/_gpucheck.abi3.so` (re-measured: all three probes load the installed `.so`), and no Makefile target, D8 change, gate (D-G3, D-G12, D-G16) or CI job defines or runs `make test-order` [code: Makefile has setup, build, test, test-all, cov, preflight, lint, fmt, wheels, cuda]. So the round-12 poison check passes whatever the scan map is, an aliased map passes every $0 check and shows only in paid R-D2 (D-G3 on the RTX 5090), and the engineer must guess the variant build and run steps. Caught-before-spend does not apply to the poison and order runs (they stay green); only I10 goes red. Fix note: the review's separate non-editable venv per variant works, but a simpler fix that matches CLEAN1 step 5 is to run each variant child with `-I -S` and set `sys.path` in the child code (no `.pth`, so no editable finder); either way add the `test_hooks()` load guard, a red-first planted alias run, a defined `test-order` target in D8, and name it in D-G3's Mac half and D-G16.
- qa-2..qa-9: minor as filed (not re-rated).

Counts: confirmed 1, refuted 0, downgraded 0.

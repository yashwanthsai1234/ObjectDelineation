# Verdict: final review round 13, CUDA (final_r13_cuda.md)

## Skills used
- nvidia-skills:nvidia-skill-finder: loaded. Review only, so I installed nothing. No catalog skill covers the custom engine's pinned staging. Nothing changed.
- engineering:code-review: checked each claim against `Block::get`, `Lane` and every `L.out.get` site in cuda/gpucheck.cpp.
- engineering:debug: traced the RL reward path (prefetch -> score_batch -> score_many -> score_polygons) to see which engine entry runs first in the process.
- engineering:system-design / engineering:architecture: checked the D5c rule (pinned bound per block vs per entry) against the D7b data flow.
- engineering:documentation: one line per item, with the doc quote and the code.
- huggingface-skills:huggingface-best: no model choice in scope. Nothing changed.

## Verdicts
- cuda-1 (major): CONFIRMED major. D.md D5c says "the engine sizes each pinned block at its first call from a stated bound" and lists bounds only per entry (`score_polygons` in/out, `eval_sets`). D7b ("It returns the handles plus the scalars the host needs ...") never says how those scalars reach the host. Today every D2H result goes through the lane's grow-only pinned `L.out` [code: gpucheck.cpp:991, 1572, 1623, 1638, 1681, 1799; `Block::get` min 4,096 B and free+alloc on growth :139-149]. In RL, D7b's batch runs before the step's first `score_polygons`. So an engineer who ports D7b in the house style sizes `L.out` small (a few KB). The 256-answer `score_polygons` (defects alone 26,624 B) then grows it once, `pinned_grows` becomes 1, and D-G15 ("`pinned_grows == 0` over the whole run") fails on correct code. This shows only after the paid V3 B3 run (S13 `test_b04` reads the pulled records). No $0 test runs `add_tables_from_raters` and then `score_polygons` in one fresh process with a `pinned_grows` check (I05, I08, I15/GAP-18, test_tables_device.py, D-G8 cold replay do not). So the gap is not caught before spend. Fix in D.md D5c, D7b and D.9.3, as the review states.
- cuda-2 (minor): minor as filed. The measurement is right: `seconds` val 414.687/415.378, test 630.456/633.567 [measured: runs/r5/eval/eval_*/final_eval.json]. Only the wording changes.
- cuda-3 (minor): minor as filed. No measurement shows a failure.
- cuda-4 (minor): minor as filed. The wording is wrong, but nothing breaks.
- cuda-5 (minor): minor as filed. `grep` prints `158:  Block<true> devA, devB, tmp;` (the doc cites :162).
- cuda-6 (minor): minor as filed. The thread placement costs a few ms [estimate].

Counts: confirmed blocking/major 1, refuted 0, downgraded 0.

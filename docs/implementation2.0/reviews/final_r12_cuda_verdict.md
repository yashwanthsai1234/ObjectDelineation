# Verdict: final review round 12, CUDA (final_r12_cuda.md)

## Skills used

| Skill | What it changed |
|---|---|
| nvidia-skills:nvidia-skill-finder | Verification only; no catalog install. Kept the check on the repo CUDA source (CUB scan scratch, Lane arenas). |
| engineering:code-review | Read K1 `scan`/`Engine` (k1poly.cpp:136-165, 824-884) against gpucheck `Block`/`Lane` (gpucheck.cpp:131-160) line by line. |
| engineering:debug | Traced why an aliased scratch block is silent on the Mac: the CPU `scan` ignores `tmp` (k1poly.cpp:144). |
| engineering:system-design | Checked how gpucheck uses `L.dev`/`L.score` today (gpucheck.cpp:989, 1565, 1592, 1804) to judge which mappings an engineer might pick. |
| engineering:architecture | Checked K1 DESIGN.md section 9 item 1: it says "Reuse gpucheck's `Lane` ... instead of K1's own `Engine`" and also does not place `tmp`. |
| engineering:documentation | One line per item with the quoted doc text and code line. |
| huggingface-skills:huggingface-best | No model choice in these items. Nothing changed. |

## Verdicts

- cuda-1: CONFIRMED (major). D.md "D2" step 1 says "Use gpucheck's `Lane` ... and its arenas in place of K1-draft's `Engine`", but K1's `Engine` holds three device blocks (`Block<true> devA, devB, tmp;`, k1poly.cpp:162) and `Lane` holds two (`Block<true> dev, score;`, gpucheck.cpp:157); no source file says where CUB's scratch `tmp` goes (scan calls at k1poly.cpp:824-825, 871-872, 876, 884 run while layouts A and B are live), `Block::get` returns the same base when big enough (gpucheck.cpp:139), and the CPU `scan` ignores `tmp` (`(void)st, (void)tmp;`, k1poly.cpp:144), so an aliasing mapping passes every planned $0 check and fails only in paid R-D2 (D-G3); the fix (add `Block<true> scan` to `Lane`, the explicit mapping, D1 counters and D5c lists) touches D.md only. The `GC_POISON_SCAN` test is optional; the mapping line alone closes the gap.
- cuda-2: minor as filed (wording of what `time/rewards_s` measures; D-G8 logs it for information only). No action at blocking/major level.
- cuda-3: minor as filed (stale evidence paths for kernel logs; no wrong result or spend). No action at blocking/major level.
- cuda-4: minor as filed (stale runbook line numbers; the S09 and preflight greps stay the real check). No action at blocking/major level.

Counts: confirmed 1 (major), refuted 0, downgraded 0 (items 2-4 were already minor).

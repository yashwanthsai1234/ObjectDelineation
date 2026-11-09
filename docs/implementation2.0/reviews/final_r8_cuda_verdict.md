# Verdict: final review round 8, CUDA (final_r8_cuda.md)

## Skills used

| Skill | What it changed |
|---|---|
| nvidia-skills:nvidia-skill-finder | No catalog skill is needed to verify these items. Nothing installed. |
| engineering:code-review | I read run_box.sh:15-16, policy.py code_files :374-383, cmd_run :464-496, cmd_pull :499-508, .gitignore, checker_gpu.yaml. These lines confirm items 1 and 3. |
| engineering:debug | I traced each R-D2 output path from the cwd of run_box.sh (/workspace/src) to the folder that pull copies (/workspace/runs/<run>). This is item 1. |
| engineering:system-design | I checked sync, run, pull and teardown against the gate evidence each rental must bring back (items 1, 2). |
| engineering:architecture | No design change. Items 4-9 stay minor. |
| engineering:documentation | Each confirmed item quotes the source heading and a code line. |
| huggingface-skills:huggingface-best | No model choice is in scope. No change. |

No web check was needed. Every item depends on repo code and doc text, and I checked both directly.

## Verdicts

- cuda-1 (blocking): CONFIRMED. D.md "D8" writes `runs/d/sanitize.json`, `--junitxml runs/d/gpu.xml` and `--dump-replay runs/d/replay_inputs.npz`, and the benches have no `--out`. `run_box.sh:16` runs `cd "$(dirname "$0")/../.."`, so these files land in /workspace/src/runs/d. `cmd_run` (policy.py:489) also cds to /workspace/src. `cmd_pull` (policy.py:505-507) copies only /workspace/runs/<run>/. D1's `d_gpu_test.yaml` spec and front 6.3 V3 item 1 add no copy step, unlike `checker_gpu.yaml`, which copies by hand. Result: the D-G5/D-G7/D-G8 evidence is lost at teardown, and the bench minutes are paid for nothing.
- cuda-2 (major): CONFIRMED. D.md D7 step 5, D-G14, D-G15 and D-G17 say to "record" `gc.stats()` keys, step deltas for steps 2-6, `pinned_grows` and an nvidia-smi sampler "before scoring starts". No file list in A, B, C, D or E names the code that writes them: a grep for `arena_bytes`, `pinned_grows`, `engine/` and `topology.json` finds them only in D's gate text. B.md I22b's `b3_v5.yaml` adds nothing. The D-G17 scoring point is inside `final_eval`, so a sampler outside the process cannot find it. The engineer must invent the recording, and the V3 gates cannot be read as written.
- cuda-3 (major): CONFIRMED. D-G7 says "The tar ships in the code sync" but gives no Mac path. D8 step 5 uses the placeholders `make -C <old dir>` / `<dir>` and has no `tar -x` step. `<pre-D2 sha>` has no rule. `--dry` builds the tar in a temp folder, so it never proves the shipped path. The code sync drops `runs/` and `build/` (.gitignore:16, 23; policy.py:377-380). Only a repo-root tar ships, because `git check-ignore old_engine.tar` returns nothing. The result is a missing step plus an undefined value, which fails `set -euo pipefail` (run_box.sh:15) in paid time.
- cuda-4 (minor): stays minor. It is a real layout slip, but I14 catches it on the Mac at $0.
- cuda-5 (minor): stays minor. It is a stale citation and gives the same bytes (D-G10).
- cuda-6 (minor): stays minor. No measured or documented driver floor exists.
- cuda-7 (minor): stays minor. The implicit API seams are implementation detail.
- cuda-8 (minor): stays minor. It is wording only ("per step" in a replay).
- cuda-9 (minor): stays minor. Gate (b) has no workers after D7 step 6. This is wording.

Counts: confirmed 3 (1 blocking, 2 major), refuted 0, downgraded 0.

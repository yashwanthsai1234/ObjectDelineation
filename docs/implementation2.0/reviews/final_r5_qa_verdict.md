# Verdict on final_r5_qa.md (fresh skeptical verifier, 2026-11-08)

## Skills used

| Skill | What it changed |
|---|---|
| nvidia-skills:nvidia-skill-finder | No new NVIDIA product claim, so no catalog lookup. It made me check the box toolchain: CMake builds `_gpucheck` for the CPU unless `OBJDEL_CUDA=ON` (CMakeLists.txt:15), which is the base of qa-3 and qa-6. |
| engineering:code-review | Read `train_rl.py:23-36` (flag list, no `--numerics`), `vast/policy.py:399, 432-433, 484-492, 527-536`, `hf_pull.py:15-22`, `gpucheck_py.cpp:113`, `tests/test_gpu_scoring.py` (compares the device scorer with the numpy host oracle, not with a second build). |
| engineering:debug | Traced each failure: I47 parse of step 3 on the `v1-b1` tree; what `snapshot_download(local_dir=...)` does with two revisions in one folder; which build is on the V4 box after `sync_up`. |
| engineering:system-design | Put front 1.1, 99_back 2b, B.4.1, B.6.2, A-G3 and U49 side by side for the view-identity row (qa-2). |
| engineering:architecture | Checked order rule 2 (tag holds no B2.2) against the owner of `--numerics` (B2.2) for qa-1. |
| engineering:documentation | Checked whether each problem is already handled in another part of the doc: only qa-5 is (B.9.8 seam list names `BudgetStop(clock=...)`). |
| huggingface-skills:huggingface-best | No model choice. Used only for qa-4: where an HF revision lands on the box. |
| Web | Not needed: every claim was checked against repo code or the source text. |

## Verdicts (blocking/major items)

- qa-1 | confirmed (blocking) | B1.5 step 3 (B.md:744) passes `--numerics bf16-seqmask`, a B2.2 flag (B.md:1162) that order rule 2 keeps off the tag (00_front.md:195 "holds no C, D, B2.2-B2.4, B3"), and `train_rl` has no such flag today (train_rl.py:23-36), so I47 fails on the tag or argparse exits 2 on the paid box; I22 also checks `b2_fp16.yaml`/`b3_v5.yaml` that do not exist on the tag.
- qa-2 | confirmed (blocking) | Front 1.1 row 4 (00_front.md:216) and 99_back 2b (99_back.md:67) say a view-identity failure means "Stop, R5", while B.6.2 (B.md:1598-1602), B.4.1 (B.md:434, no view identity), A-G3 (A.md "A val failure does not block Gate B1") and the second half of U49 say Gate B1 continues on A3's items; U49's first half lists the correlation as a hard rule, so the next paid step (B-S2) depends on which part is read.
- qa-3 | confirmed (blocking) | Front 6.3 V4 (00_front.md:794-798) and E.8 (E.md:1449-1455) have no CUDA build step, while `sync_up` installs the CPU build (policy.py:399, CMakeLists.txt:15) and a `steps:` plan runs only its steps (policy.py:432-433), so the A.9 hook fails GP1-GP5 on the box and the E-G10 sampler arm is scored by the CPU build (owner rule: CUDA checker/scorer), spending V4 ($0.54) for nothing.
- qa-4 | confirmed (major) | Front 6.3 V3 item 0 (00_front.md:779) says "B1.5 gives `fetch_hf.sh` an output folder name", but B1.5 (B.md:762-777) defines only `--rows`, and `cmd_fetch_hf` always pulls to `{rd}/ckpt` (policy.py:534-535, `snapshot_download(local_dir=...)` hf_pull.py:22), so the engineer must invent the option and the paths B2.3b and `b3_v5.yaml` read; without it the b1-rows merged model (uploaded in the form `--ckpt` reads) lands in the same folder as `rl-full-3ep`.
- qa-5 | downgraded (minor) | The B.9.8 seam list (B.md:2069) already names `BudgetStop(clock=...)` as an assumed seam and U46 (B.md:1937) requires exact fake-clock trip times, so an engineer following the doc adds the clock; the dead `pytest.skip` branch and the missing `clock=` in the B1.5 signature are wording, not a wrong outcome.
- qa-6 | confirmed (major) | G02b (A.md:1874-1879, "equal with `==` under both builds") and GP1 (E.md:1709, "scored with `checker_build="cuda"` and `"cpu"`") need two builds in one process, but the venv holds one `_gpucheck` whose `build()` returns one value (gpucheck_py.cpp:113; CMake switch at install, CMakeLists.txt:15); no golden or second venv is named, so the reference must be guessed and the gate is untestable as written.

## Minor items (not re-graded; already minor)

- qa-7 .. qa-14 | unchanged minor | Already minor in the review.

Totals: confirmed 5 (3 blocking, 2 major), refuted 0, downgraded 1.

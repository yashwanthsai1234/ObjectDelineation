# Final review round 9: CUDA and GPU performance (Implementation.md)

## Skills used

| Skill | What it changed |
|---|---|
| nvidia-skills:nvidia-skill-finder | No catalog skill was needed. I checked the CUDA facts against PyPI and the repo instead. Nothing was installed. |
| engineering:code-review | I read `run_box.sh:1-36`, `Makefile`, `CMakeLists.txt`, `pyproject.toml`, `policy.py:374-508` (`code_files`, `sync_up`, `chain`, `cmd_run`, `cmd_pull`), `gpucheck.cpp:60-215, 1570-1645, 1677-1829`, `reward/cache.py:80-265`, `K1/k1poly.cpp:125-150`, and `fast.py:604-616`. These reads give items 1 and 3-4. |
| engineering:debug | I followed each R-D2 bench input from the box working folder (`/workspace/src`) to the folder where the bundle unpacks (`/workspace/data/lidc`). This gives item 1. I followed the old-engine venv's imports through E3a's importer rewrite. This gives item 2. |
| engineering:system-design | I checked the V3 order (item 0 sync, item 1 R-D2, item 4 B3 run) and the time line (V3 before or after E0) against the bench's inputs and imports. |
| engineering:architecture | I checked the E3a layer move (`objdel.geom.engine` is the only importer) against D's old-engine arm, which must import a pre-E tree. This gives item 2. |
| engineering:documentation | Each item quotes the source file and heading and gives one exact fix. |
| huggingface-skills:huggingface-best | No model choice is in scope. Nothing changed. |
| huggingface-skills:huggingface-papers | No paper claim needed a check in my area. Nothing changed. |

Web checks: PyPI `cuda-toolkit 13.0.3` JSON. It confirms the `nvtx` extra = `nvidia-nvtx==13.0.85.*`, `nvcc` = `nvidia-cuda-nvcc==13.0.88.*` and `cccl` = `nvidia-cuda-cccl==13.0.85.*`. The D "Common toolchain" line is correct.

Repo checks that passed (no item):
- `gpucheck.cpp` waits at :1574, :1626, :1639 and :1802, and `Drain` at :177-184. This gives 4 + 2 syncs per `score_polygons` call.
- `alloc`/`release` are at :78-84. `Block` is at :131-150. `ObjOut` is at :1199-1202.
- `K1/k1poly.cpp:140-142` has the `int(n)` cast.
- `OBJDEL_RUN_DIR` is exported by `cmd_run` (policy.py:489).
- Sanitizer 13.4.92 ran with 0 errors on driver 580.178.04 (runs/r2 log :59, `runs/r2/cuda_checker/memcheck.txt`).
- RL logs every step (771 log entries for 770 steps in `runs/rl/train_rl/rl.json`), so D7 step 5's per-step `engine/*` keys exist for steps 1-6.

---

## 1. BLOCKING: the R-D2 benches read `data/...` paths that do not exist on the box (D)

**Source:** D.md, "D1. Measure first ..." step 2 and Files (`bench_checker.py`). D.md, "D8 ..." `--gpu-check` step 6. D.md, D.6 D-G7 and D-G8.

**Doc lines:**
- D1 step 2: "the 37,093 crops of `data/checker_eval/masks.npz`" and "the first 256 rows of `data/prepared/lidc/val.jsonl` ... scored on their caches from `data/prepared/caches/lidc`".
- D8 step 6: `"$PY" -I -m experiments.phase2.perf.bench_checker --dump-replay "$OUT/replay_inputs.npz"`, then `bench_checker.py --engine old|current --out "$OUT/bench_<engine>.json"`. These lines have no data-root flag. The D1 file list for `bench_checker.py` names only `--engine`, `--stages`, `--dump-replay`, `--replay`, `--out` and `--compare`.
- D8 names `OBJDEL_DATA_ROOT` only for the pytest payload check (step 3) and the `box_input()` helper of test modules. It never names it for `bench_checker.py`.

**Code proof:**
- `run_box.sh:16` runs `cd "$(dirname "$0")/../.."`, so the cwd is `/workspace/src`. `cmd_run` also runs `cd /workspace/src` (policy.py:489).
- `code_files()` drops `^(data|docs)/` (policy.py:379), so `/workspace/src/data` does not exist.
- The bundle unpacks at `/workspace/data/lidc` (policy.py:406-417, and `data = "/workspace/data/lidc"` in `chain`).
- So `--dump-replay` (it needs `prepared/lidc/val.jsonl` and `prepared/caches/lidc/`) and the 37,093-crop workload (it needs `checker_eval/masks.npz`) fail with a missing file. They fail in step 6, after the build, `pytest -m gpu` and the four sanitizer tools have run.
- `set -euo pipefail` (run_box.sh:15) then ends the mode. D-G7 and D-G8 get no numbers, and the old-engine build minutes buy nothing.
- No $0 check sees this. I16 runs the bench on the Mac, where `data/` exists. S08 covers only test modules' `BOX_INPUTS`. S10 uses `--synthetic`. S11 greps only output paths.

**Fix (D.md, D1 Files and step 2, D8 step 6, D.9.6):**
1. In the D1 file list, `bench_checker.py` resolves every named input under `ROOT = Path(os.environ.get("OBJDEL_DATA_ROOT") or "data")`. The inputs are `ROOT/checker_eval/masks.npz`, `ROOT/prepared/lidc/val.jsonl` and `ROOT/prepared/caches/lidc/`. It writes the resolved path and the SHA-256 of each input into the output JSON. The old-engine subprocess inherits the variable.
2. D1 step 2: write the inputs relative to the data root ("`checker_eval/masks.npz` under the data root, `data` on the Mac").
3. D8 step 6: one sentence. `OBJDEL_DATA_ROOT={data}` comes from `d_gpu_test.yaml` step 1, and every bench call in `run_box.sh` inherits it. Before the first bench, the mode checks that the three inputs exist under it (the step 3 check already lists them).
4. Add a test to I16 or S11 (no marker, reads no `data/`). With `OBJDEL_DATA_ROOT=tmp`, which holds a fake 8-mask `checker_eval/masks.npz`, 4 `val.jsonl` rows and 4 caches, `--dump-replay` and the crops workload read only from `tmp`. With the variable pointing to an empty folder, the bench exits 1 and names the missing file.
5. Make the same change for `r_d1.yaml`'s bench call (D1 step 4) if R-D1 ever runs.

## 2. MAJOR: E3a rewrites `bench_checker.py` to `objdel.geom.engine`, which the pre-D2 old venv does not have (D, E)

**Source:** D.md, "D1 ..." Files (`bench_checker.py`) and D.6 D-G7. E.md, "E3a" (the TID251 ban paragraph). 00_front.md, section 1, time line items 5-7 and "End state under (b')".

**Doc lines:**
- D1: "Its module-level imports are only the standard library, numpy and `objdel._gpucheck`."
- D-G7: "the old venv runs the current script by path", and "`bench_checker.py --engine old` runs with that venv's Python". The old venv is built from `git archive ... pre-d2 cuda src ...`.
- E3a: "Then add the TID251 ban on `objdel._gpucheck`. In the same commit, rewrite every direct importer to `objdel.geom.engine`. The list is re-grepped on the E0 commit (D adds importers ...)."
- Front time line: "Under a raise that funds V2 or V3 before E, E also waits for C10 and B3.10." Also: "If E0 has started, the owner rebases the branch onto post-E paths". So V3 (R-D2) can run on a post-E tree.

**Code proof:**
- Today `src/objdel/` holds `checker data eval geometry reward targets train` and no `geom` package [measured: `ls src/objdel`, 2026-11-09]. The `pre-d2` tag is pre-E (D-G7: the tag is set when D2 starts, and D lands before E).
- The old venv installs only numpy and the old tree with `--no-deps` (Makefile:42).
- On a post-E tree, E3a's re-grep finds `bench_checker.py` and rewrites it to `from objdel.geom import engine`. The `ruff` TID251 ban forces this, because `make lint` covers `experiments` (Makefile `lint`). Then `/opt/old/bin/python bench_checker.py --engine old` raises `ModuleNotFoundError: objdel.geom` in R-D2 step 6. Paid time follows, and D-G7 and D-G8 cannot be read.
- S10 does not catch this. Its import hook allows any `objdel` name and runs against the current tree.

**Fix:**
- In E.md E3a, add `experiments/phase2/perf/bench_checker.py` as a named exception: "keeps `from objdel import _gpucheck`. Reason: its `--engine old` arm runs in a venv of the pre-E `pre-d2` tree (D-G7). It gets one `per-file-ignores` TID251 entry."
- In D.md D1, add one sentence that says the same. S10 gets one more case: run the script under a fake package tree that holds only `objdel/__init__.py` and `objdel/_gpucheck` (a stub with the old entries), and assert that it imports and runs `--engine current --replay`.
- U4 (E.9) adds this file to its allowlist with that reason.

## 3. MINOR: D-G7 says the arms run in "alternating subprocesses", but D8 step 6 shows one call per engine (D)

**Source:** D.md D.6 D-G7 ("`--engine old` runs with that venv's Python, `--engine current` with the main one, in alternating subprocesses, 3 rounds each") and B01 ("alternating old/new subprocess"). D8 step 6 shows `bench_checker.py --engine old|current --out "$OUT/bench_<engine>.json"`, which reads as two separate calls, one per engine.

**Problem:** Two sequential calls do not alternate. Clock and thermal drift across the 6-minute bench then biases the old/new ratio. The doc does not say which process launches the other engine's subprocess.

**Fix:** In D8 step 6, write one orchestrating call: `"$PY" -I -m experiments.phase2.perf.bench_checker --engines old,current --rounds 3 --out "$OUT"`. In each round it starts one subprocess per engine, in alternating order (old venv Python for `old`). It writes `bench_old.json` and `bench_current.json`. Do the same for `--replay`.

## 4. MINOR: D.4 lists every counter as "new in D1", but the pool counters exist only after D5c (D)

**Source:** D.md D.4 "Engine counters (new in D1, ...)". This list holds `pool_allocs`, `pool_frees`, `pinned_grows`, `pool_reserved`, `default_pool_reserved` and `label_pairs`.

**Problem:**
- The pre-D2 old engine (D-G7) has D1 only. D5 and D5c come later.
- So the doc does not say whether the D1 binding must already return these keys (as 0) or adds them later.
- `bench_checker --compare` and S12 read the keys. A missing key in the old arm's `stats` would raise a `KeyError` in the box compare.

**Fix:** Add one sentence to D.4: "D1 adds every key to `gc.stats()`. Keys whose mechanism lands later (`label_pairs` in D5, the pool keys and `pinned_grows` in D5c) report 0 until then. `bench_checker --compare` reads only `waits` and `drain_syncs` from the old arm."

---

Counts: blocking 1, major 1, minor 2.

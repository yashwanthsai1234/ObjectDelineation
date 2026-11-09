# Review of Phase A (`docs/implementation2.0/A.md`, "measure fairly")

Reviewer: fresh critical pass, 2026-11-08. REPO = $REPO.
Tags: [measured: ...] = I ran it this session; [code: file:line] = I read it; [web: URL]; [estimate].

## Skills used

| Skill | What it changed in this review |
|---|---|
| nvidia-skills:nvidia-skill-finder | I fetched the live catalog again (`raw.githubusercontent.com/NVIDIA/skills/main/skills.sh.json`, saved as `(planning probe, not kept: nv_catalog_live.json)`). The only profiling skill is `deepstream-profile-pipeline`, which does not apply. A.md's "no Nsight/NVTX skill" claim holds. Nothing to install. |
| engineering:code-review | Drove the line-by-line check of every file and line A.md names. It found the missing `--trace` flag, `configs.processor` (missing), the dead signature checks in A7, and the import-time CUDA probe in A6. |
| engineering:debug | Reproduce-first. I re-ran the R5 bootstrap (item, nodule), the miss rate, the RNG equality of `signal.draws`, and the polygon-vs-mask exactness of doctor targets. The checks confirmed some claims and removed one false lead (items 25-26). |
| engineering:system-design | Checked the Mac/box data flow for failure points: answers written after scoring, an nsys hang with no time limit, and the step range that misses the data loader. |
| engineering:architecture | Judged the inline ADRs (pin + signature check, CPU oracle for set metrics, duplicated-rater stack for k = 2). Pin-first ordering makes the signature ADR a no-op (item 4). |
| engineering:documentation | Checked A9. The existing docs that quote item CIs are not on its list (item 22). |
| huggingface-skills:huggingface-best | Not applicable: Phase A chooses no model. I used it only to confirm that no HF leaderboard covers multi-rater metrics, so A.md is right to take its metric definitions from papers. |
| engineering:testing-strategy | Checked that each test can run in a clean checkout. The R5 fixture lives in a git-ignored folder (item 7). Two gates are tautologies (items 15, 16). |
| engineering:tech-debt | Checked the archive and removal items (lab/A4, a0_ceiling) and found where A.md adds new code next to old code: `profile.py` already computes GPU busy share and syncs (item 2). |

Web checks: Nsight Systems user guide (`--cuda-graph-trace` values) [web: docs.nvidia.com/nsight-systems/UserGuide];
NVTX GPU projection report `nvtx_gpu_proj_trace` [web: forums.developer.nvidia.com/t/nvtx-with-gpu-timing/268356];
the `nsight-systems-cli` apt package in NVIDIA's devtools repo [web: developer.download.nvidia.com/devtools/repo-deb/amd64/Packages].

## Verified correct (no action)

- Item CI [0.6334010210209712, 0.6561989221793691] and nodule CI [0.6279017988891661, 0.6620867640136667], 401 nodules,
  miss rate 0.10818, IoU on found images 0.72188 [measured: runs/r5/eval/eval_test/items.jsonl, `signal.draws/ci`].
- TRL 1.17.1 signatures of the 5 checked methods, `shuffle_sequence_dict`, and `VLLMGeneration._dist` all exist
  [measured: `inspect.signature`].
- LIDC `group_id` is the patient (`scan.patient`) [code: src/objdel/data/lidc.py:316]. The eval view records carry it.
  Val has 1,084 RL items: k = 4: 473, k = 3: 330, k = 2: 281. There are 78 patients, one series each [measured].
- CPU and CUDA checker builds are bit-identical (`-ffp-contract=off`, `--fmad=false`) [code: docs/research/notes/gpu_scoring.md:72].
  So the `--rescore` bit-equality gate is feasible.

## Findings

1. **blocking: A-G6 (gate) and A6 step 8: the gate cannot be computed from the listed reports.**
   Problem: `nvtx_sum` gives totals per range name over the whole capture. It has no nesting and no per-step split, so
   "top-level ranges cover >= 90 % of `step` on steps 4-5" cannot be read from it. No method is given for the
   "GPU idle share" either. NVTX ranges are CPU push/pop times. GPU work queued inside a range runs later on the GPU
   timeline.
   Fix: in `nsys_capture.sh`, also run `nsys export --type sqlite -o {rd}/nsys/step45.sqlite {rd}/nsys/step45.nsys-rep`
   and `nsys stats --report nvtx_pushpop_trace,nvtx_gpu_proj_trace,cuda_gpu_trace,cuda_api_trace --format csv`.
   Add `experiments/phase2/nsys_baseline.py` (Mac, offline, numpy only). Per `step` instance it computes (a) the
   coverage: the summed duration of depth-1 child ranges / the step duration (from `nvtx_pushpop_trace`); (b) GPU busy =
   the length of the union of the kernel + memcpy + memset intervals inside the step window / the step duration (from
   `cuda_gpu_trace`), idle = 1 - busy; (c) host syncs = the count of `cudaStreamSynchronize`, `cudaDeviceSynchronize`,
   `cudaEventSynchronize` and synchronous `cudaMemcpy*` (no `Async`) in `cuda_api_trace` inside the window. Write
   `nsys_baseline.json`. Give the gate as numbers from that file. Before the rental, test the script on a tiny CSV fixture.

2. **major: A6 step 9 (fallback) and the A6 file list: the fallback flag does not exist, and the callback cannot be wired as written.**
   Problem: `train_rl` has no `--trace` [code: experiments/phase2/train_rl.py:23-36]. `TraceStep` is wired only in
   `g_probe` (`getattr(a, "trace")`) and `train_sft` [code: g_probe.py:93,145]. `NvtxSteps` must go through
   `g_probe.make_trainer`, which builds the callback list itself [code: g_probe.py:92-95]. Its file is not in A6's
   list. Also, `experiments/phase2/profile.py` already computes the GPU busy share and the host-wait counts from a
   torch.profiler trace (`GPU_CATS`, `SYNCS`) [code: profile.py:41-43]. A6 rebuilds this from scratch.
   Fix: add `g_probe.py` (`make_trainer`: append `callbacks.NvtxSteps(a.profile_window)` when
   `getattr(a, "profile_window", None)`, and the same for `TraceStep` from `a.trace`) to the file list. Add `--trace DIR`
   to `train_rl` next to `--profile-window`. For the fallback, reuse `profile.py`'s busy-share and sync code: move the
   function to `src/objdel/train/trace_stats.py` and test it, or call it. Do not write a third copy.

3. **major: A6 step 4 (`NvtxSteps`): the `step` range leaves out the data loader and the log, so it is not the 23.4 s step.**
   Problem: in transformers 5.20.1, `get_batch_samples` (all 16 micro-batch fetches, OpenOnce image opens) runs *before*
   `on_step_begin`. `_maybe_log_save_evaluate` (FastGRPO.log: the `is_stats.tolist()` read, TRL's completions table) runs
   *after* `on_step_end` [code: .venv transformers/trainer.py:1856,1877,1929-1931]. The CPU phase the BRIEF names
   (3.0-3.2 s) is partly outside the range. Also, backward has no range: `_compute_loss` covers the forward only, and
   `accelerator.backward` runs in `training_step`. The 90 % coverage gate then fails by design.
   Fix: define `step` from `on_step_end(s-1)` to `on_step_end(s)`: push in `on_step_end` and pop at the next
   `on_step_end`, or emit `mark("step_end")` and measure between marks. Start the profiler at `on_step_end` when
   `state.global_step == s - 1`. Add a `micro` range by overriding the public `training_step` in FastGRPO
   (`with nvtx.range("micro"): return super().training_step(...)`). Add the names `micro`, `loader` (wrap
   `get_batch_samples`; it is public) and `log` to the fixed list.

4. **major: A7 step 2 (check order): the signature checks never fire.**
   Problem: `check()` raises on the version first. When the version equals 1.17.1 (PyPI wheel), the signatures equal
   the ones recorded. When it differs, the version check has already raised. So the method, global and attribute checks
   only run under monkeypatched tests. The "pin + signature" ADR then protects nothing beyond `trl.__version__ == PIN`.
   Fix: make the version mismatch a soft failure with an explicit opt-in. Order: methods, globals, attrs, then version.
   On a version mismatch with all contract checks passing, raise unless `OBJDEL_TRL_UNPINNED=1`. With the variable set,
   print a warning naming the version. Then a TRL upgrade can be tested by the contract alone (a planned job,
   Architecture.md §6.2). Add a test: version patched to "1.18.0", env var set, contract intact -> passes with a warning.

5. **major: A7 "Verified signatures" list: it leaves out the names that break silently.**
   Problem: FastGRPO also reads private or undocumented trainer state and input keys that a signature check cannot see:
   `self.use_vllm`, `self.loss_type`, `self.num_iterations`, `self.top_entropy_quantile`,
   `self.vllm_importance_sampling_{correction,mode,clip_min,clip_max}`, `self.vllm_generation.{llm,model,sync_weights}`,
   `llm.{generate,wake_up,sleep}`, and the input keys `sampling_per_token_logps`, `completion_mask`, `tool_mask`,
   `importance_sampling_ratio` [code: src/objdel/train/fast.py:524-603,652-690]. C10 (`skip_old_logps`) depends on a
   behaviour: TRL skips the old-logprob pass when `vllm_importance_sampling_correction` is False [code: fast.py:540-547].
   `is_ratio` claims to copy grpo_trainer.py 2776-2824 [code: fast.py:351-363]. All exist in 1.17.1
   [measured: substring counts in `inspect.getsource(GRPOTrainer)`].
   Fix: add these to `ATTRS` (substring in `inspect.getsource(GRPOTrainer)` or `VLLMGeneration`) and to a new `KEYS`
   tuple (substring `"<key>"` in the source of `_generate_and_score_completions` / `_compute_loss`). Add one behavioural
   test in `tests/test_trl_adapter.py`: TRL's own ratio code path and `fast.is_ratio` give equal tensors on a toy
   (logps, sampling, mask) for `sequence_mask` and `token_truncate`. Update the "10 items" count.

6. **major: A6 step 1 (`nvtx.enabled()` evaluated at import): it can break fork safety for the scoring pool.**
   Problem: `torch.cuda.is_available()` calls `cudaGetDeviceCount`, which runs `cuInit` (it poisons fork) unless
   `PYTORCH_NVML_BASED_CUDA_CHECK=1` [code: .venv torch/cuda/__init__.py:213-227]. `train_rl` imports `g_probe` (so
   `fast`, so `nvtx`) *before* `serve.fork_pool` [code: train_rl.py:40-50]. The pool workers run the CUDA checker
   (`metrics` -> `score_many`). `fork_pool`'s guard checks only `torch.cuda.is_initialized()` [code: serve.py:167-169],
   and that stays False after `is_available`. So the guard does not catch this case.
   Fix: at import, read only `os.environ.get("OBJDEL_NVTX") == "1"`. Run the CUDA check lazily at the first `range()`
   call and cache the result (after the fork). Test: importing `objdel.train.nvtx` with the variable set leaves
   `torch.cuda` uncalled (spy on `torch.cuda.is_available`).

7. **major: A1 step 5 and A-G1 (R5 fixture): the exact-CI tests have no committed input.**
   Problem: `/runs/` is git-ignored [code: .gitignore:23]. `runs/r5/eval/eval_test/items.jsonl` does not exist in a clean
   checkout or on CI, so the "exactly [0.6334..., 0.6561...]" tests cannot run.
   Fix: commit `tests/fixtures/r5_test_items.jsonl.gz` with only `item_id` and `iou_tc` (2,200 rows, about 60 KB
   gzipped [estimate]). Tests read that file. State in A1 that `group_id` for the patient test comes from a small
   synthetic map, not the data folder.

8. **major: A3 `final_eval.main` (write order): the answers can be lost.**
   Problem: A3 lists `answers.jsonl` among the outputs but gives no order. Today all files are written after
   `score_items` [code: final_eval.py:75-84]. If scoring or the pool crashes (the R5 log ends with a `cumem_allocator`
   error), the 17.5 min of sampling is lost. The reason for A3 is that answers are never resampled.
   Fix: write `<out>/answers.jsonl` (and a partial `final_eval.json` with `view_manifest_sha256`, seed, K, `code_sha`)
   immediately after `sample.sample` returns, before `score_items`. Then `--rescore` can recover the run.

9. **major: A3 plan file / `nsys_capture.sh`: an nsys hang has no time limit, and the rental bills until the cap.**
   Problem: the step is `"~"` (a failure does not stop the chain), but a hang is not a failure. Capture-range runs with
   colocated vLLM are untested (A.md §7 says so). The cap is 1.5 h.
   Fix: wrap the profile command in `timeout --signal=INT 900 nsys profile ...` (15 min: model load + 6 steps at <= 31 s
   + nsys flush [estimate]). After it, run `teardown`-safe cleanup (`pkill -f nsys`). Put the nsys step *last*, as now,
   and make `pull.sh` copy `eval_*` before it starts if the runner allows it.

10. **major: A6 step 8 (nsys install): "pick the URL later" is vague, and Blackwell support is not checked.**
    Problem: the install step is left open. The RTX 5090 (sm_120) needs a recent Nsight Systems. An old build may run
    but not trace the GPU.
    Fix: install from NVIDIA's devtools apt repo, CLI only, at a pinned version:
    `apt-get install -y nsight-systems-cli-<YYYY.N.N>` after adding the repo key
    [web: developer.download.nvidia.com/devtools/repo-deb/amd64/Packages]. Choose a version whose release notes list
    Blackwell. Run `nsys status -e` and `nsys --version` first, and abort the nsys step (fallback) if the status
    reports a problem. Write the exact commands in `nsys_capture.sh` and run `curl -sI` against the package URL on the
    Mac before the rental.

11. **major: A6 step 7 (`zero_share`): the counter mixes two different failures.**
    Problem: in `sequence_mask` mode, `r` is 0 either when it is masked above the cap (`r > hi`) or when
    `exp(sum diff)` underflows in fp32 (sum < about -103, trainer << vLLM) [code: fast.py:355-363]. B2 needs to know
    which side loses the answers.
    Fix: have `is_ratio` also return the raw ratio, or count before masking. Add two device counts, `masked_hi`
    = `(r_raw > hi).sum()` and `masked_lo` = `(r_raw < lo).sum() + (r_raw == 0).sum()`, and log
    `sampling/importance_sampling_ratio/{zero_share,hi_share,lo_share}`. The single host read is unchanged.

12. **major: A-G3 (sanity rule): the gate is too weak to catch a regression.**
    Problem: "new test mean inside the R5 nodule CI" accepts any shift up to about ±0.017. Both runs use the same items
    and the same model, so the population-sampling CI is the wrong yardstick. A box image with a different vLLM or
    torch could lower IoU_tc by 0.01-0.015 and still pass.
    Fix: compare paired by `item_id` against `runs/r5/eval/eval_test/items.jsonl`. Pass when the nodule-bootstrap CI
    of the per-item difference contains 0 and |mean diff| <= 0.01. Also report the per-item correlation. Keep the
    current rule as a secondary line.

13. **major: A8 (scope): it counts doctor outlines only, not model answers.**
    Problem: the limit question (BRIEF §2: 38 % of R5 steps hit 1,536, runaways masked) is about the *model's*
    lengths. A3 now saves `tokens[K]` per answer, but A8 does not use them.
    Fix: add `--answers runs/a3/eval_{val,test}/answers.jsonl` to `outline_lengths.py`. For each L, report the share of
    model answers > L, their legal share, and their mean IoU_tc from `metrics.jsonl`. Then the decision table shows
    "limit L cuts N model answers, of which M were legal with IoU_tc >= 0.5" beside the doctor count.

14. **major: A-G8 and A8 step 4: the gate cannot close inside Phase A, which conflicts with the "dataset first" rule.**
    Problem: the shifted and rescaled train rows wait for B1, but building them is free on the Mac
    (`build_eval_views --split train --out data/views --procs 6`) [code: experiments/phase2/configs/rl_shifted.yaml:10].
    The owner's top rule is dataset quality, and B1 needs these views anyway.
    Fix: move the train-view build into A8 (Mac, no money). Run A8 on all three roots now and close A-G8 in Phase A.
    B1 then reuses the same views and only audits them.

15. **minor: A-G4 last sentence: a tautology used as a gate.** `sel_medoid_legal` equals the share of items with
    >= 1 legal answer by construction. Fix: replace it with a real check. The test medoid's IoU_tc on test (with a
    nodule CI) moves in the same direction as on val. Report it, and never choose the method on test.

16. **minor: A-G2 (legacy check): the escape clause makes it impossible to fail.** "Within ±0.005 ... or the note
    explains" always passes. The pooling order is also unstated (mean per item, then over items, versus pooled over
    raters). Fix: compute `legacy_loo_dmax` both ways. Pass if either is within ±0.005 of 0.734 / 0.732. Otherwise
    the gate fails and A2 stops until the owner accepts the difference.

17. **minor: A4 `medoid` (fallback to parsed samples): dead logic.** `sample_iou` is built on `drawn` masks, where
    illegal samples are empty [code: metrics.py:306-307]. Two empty masks have IoU 1, so the parsed fallback always
    returns the first parsed index, and its IoU_tc is 0 anyway (gated). Fix: "no legal sample -> return 0". Drop the
    `parsed` parameter and its test. Note in the docstring that legal empty answers have IoU 1 with each other, so a
    majority of empty answers wins the medoid (that is intended).

18. **minor: A7 step 4 (preflight check name).** `preflight` groups its stages by `name.split()[0]` [code:
    experiments/phase2/preflight.py:266-280]. `("trl_names", ...)` would form its own stage, and the stamp key would
    not be `CODE ...`. The function must return a list of problems with signature `fn(a, cfg)`. Fix:
    `("CODE TRL private names", lambda a, cfg: _trl_problems())`, where `_trl_problems()` returns `[str(e)]` on
    `TRLContractError`. A-G7 reads `checks["CODE TRL private names"]["pass"]`.

19. **minor: A8 step 1: `configs.processor` does not exist.** Only `MODELS` and `processor_kwargs` exist
    [code: src/objdel/train/configs.py:70,134]. Fix:
    `AutoProcessor.from_pretrained(configs.MODELS[model], **configs.processor_kwargs(model, dataset)).tokenizer`, as
    `sample.export_prompts` does [code: src/objdel/train/sample.py:115-116]. Load it from the local HF cache with
    `HF_HUB_OFFLINE=1`.

20. **minor: A3 step 3 / plan header: the preflight command is incomplete.** `make preflight` needs `CONFIG=`.
    Fix: write `make preflight CONFIG=experiments/phase2/configs/r1_qwen3vl2b.yaml` (the config `launch.sh --config`
    checks, as in `rl_shifted.yaml:19-21`), and the NA-only `offers.sh RTX_5090 --pick ID`.

21. **minor: A3 `code_sha`: always null on the box.** The box gets code by `sync_up.sh`, with no `.git`, so
    `git rev-parse` fails. Fix: have the bundle or `sync_up.sh` write `CODE_SHA` (the Mac's `git rev-parse HEAD` plus
    a `git diff --stat` hash). `final_eval` reads that file first.

22. **minor: A9 (docs list) and A1 schema change: existing docs quote item CIs.** `docs/guides/evaluation.md:53` and
    `docs/research/notes/final_eval_2026-11-08.md:14` quote `0.645 [0.633, 0.656]` [measured: grep]. After A1,
    `ci95` means the nodule CI. Fix: add both files to A9 and relabel the old numbers "item CI". Better, keep `ci95`
    unambiguous: name the keys `ci95_nodule`, `ci95_patient`, `ci95_item`, and drop the bare `ci95` in schema 2, so
    an old reader fails loudly rather than mixing meanings.

23. **minor: §1 evidence item 1 versus the code: "CT series" against "patient".** The CI labelled "CT series"
    ([0.6248, 0.6634]) is not what A1 computes. A1 clusters on `group_id` = patient. On val, series = patients (78)
    [measured], but on test they may differ. Fix: re-measure with `key = group_id` and quote that number.

24. **minor: A5 step 1 (boundary IoU d): `d = NSD_TAU` makes boundary IoU nearly equal mask IoU on small nodules.**
    With d = 2 px, a 2 px band covers most of a nodule about 8 px wide. Boundary IoU then adds little over IoU there.
    Fix: keep d = 2 for comparison with NSD. Also report `boundary_iou` by size tercile in the result note, so the
    reader sees where it differs from IoU. Clip the d + 1 margin to the frame, so zero padding happens only at the real
    image edge (as the official API does).

25. **minor: A2 (fairness, checked): doctor masks versus polygon text.** I checked whether doctors gain from raw
    masks while the model is rasterized from 0-1000 grid text. On 200 val items (628 targets), each target text
    rasterized gives IoU 1.000 with its rater mask [measured: `parse` + `to_masks`]. No action needed. Add this as a
    test (`targets[r]` rasterizes to `cache.raters[r]`), so the A2 doctor side is backed by a check.

26. **minor: A1 step 5, test 2: it can be made exact.** With one item per cluster, `signal.draws` reproduces
    `rng.integers(0, n, (B, n))` row for row (n = 2,200 and n = 7) [measured]. Fix: assert the nodule CI equals the
    old item CI *exactly*, not only "computed from the same draws".

27. **minor: A6 step 3 (host-read marks): redundant with the trace.** D2H copies and syncs already appear in
    `cuda_api_trace` (item 1). Fix: drop `mark("host_read")`, or keep one mark per Python call site only if the
    sqlite query cannot attribute syncs to a range. That leaves fewer edits in `fast.py`.

28. **minor: §3.2 GPU vs CPU: set metrics read host masks from a second rasterizer.** `metrics()` uses host
    `to_masks` [code: metrics.py:255], while the CUDA checker can return the same masks (`check_polygons`)
    [code: gpu_scoring.md:16]. The CPU-oracle reasoning is fine for A. Fix: state that `to_masks` is bit-equal to the
    device raster and name the test that proves it, or switch `_samples` to `check_polygons` masks now. Then Phase D
    swaps only the arithmetic.

29. **minor: A6 step 8: the profile command runs Python without `-I`, and graph kernels are hidden.** Fix: use
    `python -I -m experiments.phase2.train_rl` (the repo rule for data reads). Keep `--cuda-graph-trace=graph`
    (default, low overhead) for the idle share. Note in `profiling.md` that kernels inside vLLM's CUDA graphs show up
    as single graph launches. Use `node` only for a separate kernel breakdown [web: Nsight Systems user guide].

30. **minor: §8 cost / A-G6: no check that nsys slows the step.** Fix: compare `time/batch_s` and the step wall time of
    steps 3 and 6 (outside the capture window) with R5's 23.4 s. If the window steps are > 10 % slower than the
    steps around them, the note marks the baseline as overhead-affected.

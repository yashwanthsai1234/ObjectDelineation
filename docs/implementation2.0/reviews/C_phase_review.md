# Review of Phase C (GPU speed): `docs/implementation2.0/C.md`

Reviewer: a fresh critical reviewer, 2026-11-08. I opened every code location C.md cites in REPO (fast.py, configs.py,
TRL 1.17.1 `grpo_trainer.py`, `vllm_generation.py`, `utils.py`, `data_utils.py`, transformers 5.20.1 Qwen3-VL modeling
and processing) and recomputed the R5 numbers from `REPO/runs/r5/shared_rl.json`.

## Skills used

| Skill | What it changed in this review |
|---|---|
| nvidia-skills:nvidia-skill-finder | Checked the NVIDIA catalog route (C.md's saved catalog check). There is no new catalog skill for vLLM weight sync or for NVTX gates. Nsight Systems stays the tool. Items 9 and 10 cover how nsys is used, and item 21 covers the nvImageCodec PNG claim. Nothing was installed. |
| engineering:code-review | Line-by-line check of each cited `file:line`. It found the reference-pass gap in C3 (item 1), the CUDA-grid sync in C1 (item 3) and the pass-2 text work (item 2). |
| engineering:debug | Reproduce-first check of the gates. Two gates fail no matter what C1-C5 do, because vLLM's own syncs and H2D copies fall inside the step (items 6 and 7). The baseline arm is also distorted by profiling (item 9). |
| engineering:system-design | Re-ran the GPU/CPU placement table against the real data flow. vLLM's CPU-side image processing and its duplicate vision encode are missing (item 5). |
| engineering:architecture | ADR check of "delta push vs one-call full push" (item 4: delta push adds risk and saves milliseconds) and of "awake vs asleep" (item 8: the R4 evidence is not reconciled). |
| engineering:documentation | Found mislabelled gate values, numbering errors and internal contradictions (items 12, 17, 18, 19). |
| huggingface-skills:huggingface-best | No leaderboard run: Phase C does not choose a model. One risk follows: C2, C3 and C5 are specific to Qwen3-VL (DeepStack, M-RoPE, 196 tensors), so a model change from M1-M3 voids them (item 20). |
| engineering:testing-strategy | Coverage of the equality gates. With the planned 15-step box arms and checks every 50 steps, G-EQ2 and G-EQ3 get only one real batch on the box (item 11). The sync counter undercounts because of Python's warning dedup (item 6). |
| engineering:tech-debt | Checked the C10 deletions against the owner's archive rule. Checked the C4 step-6 copy of 265 lines of TRL (item 16) and the leftover CPU fallbacks (item 15). |

Web checks done: vLLM PR 10980 (V1 parallel sampling: `n>1` spawns `n` child requests, merged 2025-02-24) [web:
https://github.com/vllm-project/vllm/pull/10980] (C.md is correct). The vLLM 0.32.0 weight-transfer API (an `ipc` backend
for colocated training, a 4-phase protocol, `WeightTransferConfig` on `LLM`) [web:
https://docs.vllm.ai/en/v0.32.0/training/weight_transfer/]. The `torch.cuda.set_sync_debug_mode` docs ("experimental,
not all synchronizing operations will trigger") [web: https://docs.pytorch.org/docs/stable/generated/torch.cuda.set_sync_debug_mode.html].
nvImageCodec: GPU codecs listed are nvJPEG, nvJPEG2000 and nvTIFF, with no GPU PNG codec [web: https://developer.nvidia.com/nvimagecodec].
vLLM Qwen3-VL precomputed `image_embeds` with DeepStack: a community thread only, so the input format is unverified
[web: https://discuss.vllm.ai/t/how-do-i-precompute-multimodal-embeddings/2294].

Facts I re-checked and found correct (no action): R5 step 23.43 s; sync 1.267 s, wake 0.106 s, sleep 0.292 s; generate
11.11 s; batch 15.76 s; max-length mean 1,271; 5 of steps 3-15 hit 1,536; 6 of 15 steps have a naturally ended answer
above 1,200 [measured: runs/r5/shared_rl.json, recomputed now]. The C7 saving of 0.765 s on R5 steps is correct
(2,313 / 13 tokens x 4.3 ms). dr_grpo divides by `max_completion_length` [code: grpo_trainer.py:3286]. With `beta == 0` or
PEFT there is no reference model [code: grpo_trainer.py:1031-1037]. TRL sleeps at level 2 and pushes per parameter
[code: vllm_generation.py:377, 520-522, 772]. Both `mm_encoder_outputs` and "raise if both are given" exist
[code: modeling_qwen3_vl.py:1155-1163]. `get_rope_index` starts from zeros [code: modeling_qwen3_vl.py:966].
`split_tensor_dict` slices lists [code: trl/trainer/utils.py:1076-1078]. The G-SP1 midpoint arithmetic is internally
consistent (see item 12 for the scope error).

---

## Findings (most severe first)

1. **blocking - C3 Steps 1-3 and §3.2 step 7: the reference pass runs inside TRL before C3 swaps `pixel_values` for
   `image_ids`.**
   Problem: C3 says "after TRL builds the batch, replace `pixel_values` by `image_ids`" in
   `FastGRPO._generate_and_score_completions`. But the reference pass runs inside
   `super()._generate_and_score_completions` [code: grpo_trainer.py:2827-2846]. Its `forward_kwargs` still hold the
   256-row `pixel_values` from processor pass 2 (grpo_trainer.py:2645). So the reference pass never sees `image_ids`.
   Either it keeps the old pixel path (no C3 saving there, and G-EQ4 "reference pass" tests nothing new) or it breaks.
   C3 also never says where `EncoderCache.encode` is called "before the first forward", because that forward is inside
   TRL.
   Fix: State the mechanism. In `FastGRPO._generate` (fast.py:618), after vLLM returns, set `self.meta` (C4 RowMeta
   with `image_id`) and call `EncoderCache.encode` on the 32 distinct images. In
   `FastForward._get_per_token_logps_and_entropies`, map row index to image id from `self.meta` (or from a `host_meta`
   row field) whenever `pixel_values` arrives. Drop `pixel_values` before calling the model in both the reference and
   update passes. Make C3 depend on C4 step 1 in §8.3. Add one G-EQ4 test that monkeypatches the vision tower's
   `forward` to count calls: exactly 1 encode per generation batch.

2. **major - C2 and §3.3 row "Chat template, tokenization": processor pass 2 still renders and tokenizes 256 prompts on
   the CPU.**
   Problem: C2 memoizes only pass 1 (`_tokenize_prompts`). Pass 2 runs `apply_chat_template` once per row (256 jinja
   renders) [code: grpo_trainer.py:2640], then `self.processing_class(images=..., text=prompts_text, ...)`
   [code: grpo_trainer.py:2645]. That call expands each image pad to t·h·w/4 pads, tokenizes 256 rows of about 325
   tokens and builds `mm_token_type_ids`. TRL then throws away `input_ids`/`attention_mask` (it keeps only
   `forward_kwargs`). This is the CPU work the owner rule forbids, and it sits inside the "prepare" range that G-SP3
   caps at 0.5 s. §3.3 claims "once per prompt layout", which is false for pass 2.
   Fix: Add C2 step 5, a `PromptMemo` wrapper on `processing_class.__call__` (same install pattern as `ImageMemo`).
   Key: (prompt text, image (H, W)). Value: `mm_token_type_ids` row and length. It assembles the 256-row output by one
   pad-and-stack and delegates `pixel_values`/`image_grid_thw` to `ImageMemo`. Key the memo by the message text parts
   plus image size, not by "rendered prompt text": rendering is itself the cost (C2 step 1 has the same flaw). Mac test:
   wrapper output `torch.equal` to the plain processor on 256 rows. Count it in G-SP3.

3. **major - C1 Step 2 and Interface: a CUDA `image_grid_thw` adds 256 host syncs inside the processor.**
   Problem: `patchify` returns the grid on `device`, and `ImageMemo` stacks the per-image grid rows. Qwen3-VL's
   processor then computes `image_inputs["image_grid_thw"][image_idx].prod() // merge_length` per image to build the pad
   string [code: transformers/models/qwen3_vl/processing_qwen3_vl.py:78]. With a CUDA grid, that is one D2H sync per
   image: 256 per step in pass 2. This undoes part of G-SP4.
   Fix: State in the C1 interface that `image_grid_thw` is always a CPU int64 tensor, computed from `(H, W)` with no
   device work. Only `pixel_values` goes to the device. Add a Mac test that asserts `grid.device.type == "cpu"`. Add a
   box assertion: `count_syncs()` around the pass-2 processor call reads 0.

4. **major - C5 Interface and Step 2: the delta push adds staleness risk but saves only milliseconds over a one-call
   full push.**
   Problem: The 196 touched tensors are about 1.41 B of the model's parameters (per layer: q 2048x2048, k and v
   1024x2048, o 2048x2048, gate/up/down 3 x 6144x2048, which is about 50.3 M, times 28) [estimate: Qwen3-VL-2B config
   shapes]. A full push adds embeddings and vision, about 0.5-0.7 B more. On the device these are D2D copies of about
   2.8 GB vs about 4-4.5 GB at more than 1 TB/s, so the gap is a few ms [estimate]. R4 shows the 1.01 s awake push comes
   from about 700 Python `load_weights` calls (each rebuilds loader maps), not from bytes [measured: r4 note]. In
   return, the delta push brings: `lora_touched_names` mapping code, a "delta while asleep raises" branch,
   `--full-push-every 50`, resume handling, an untouched-hash check, and the G-EQ5 "3 delta pushes" protocol.
   Fix: Make the C5 design "ONE `load_weights` call, full push, every step" in both modes. Delete `lora_touched_names`,
   `full=`, `--full-push-every`, the sleeping-delta `ValueError` and the G-EQ5 delta clause. Keep G-EQ5 as "after a
   push, every vLLM parameter `torch.equal` to the trainer's `lora_merged_params`". Add a box micro-benchmark in
   `c_speed_check --sync` that times one-call full vs one-call 196-tensor push. Re-add delta only if the gap is
   >= 50 ms. Also record the vLLM 0.32 `weight_transfer` IPC backend as considered and rejected: it targets separate
   processes, while TRL colocate runs vLLM in-process, so one direct call is less code.

5. **major - §3.1, §3.3: vLLM's own CPU image processing, H2D copy and duplicate vision encode are missing.**
   Problem: In colocate mode, vLLM receives 256 PIL images per step (`multi_modal_data`, vllm_generation.py:723-730).
   Its input processor hashes each one, runs the HF image processor on the CPU for the distinct ones, copies their
   fp32 patches H2D, and runs its own vision tower [SP/B_gpu_audit.md row 1d: "vLLM multimodal processing of 256 PIL
   images"]. After C3, the trainer encodes the same 32 images again with the same frozen weights. The placement table
   has no row for this, so it violates "every CPU item justified". It also breaks G-SP5 (item 7).
   Fix: Add §3.3 rows "vLLM image hashing + preprocessing (CPU)" and "vLLM vision encode (GPU, duplicate)", with the
   C0 measured time. Add an investigation task C6b: feed vLLM precomputed `image_embeds` + `image_grid_thw` from the
   C3 `EncoderCache`. This removes vLLM's CPU preprocessing, its H2D copy and its encode, and it makes vLLM and the
   trainer use identical image features, which may also help B2's mismatch. Gate it on (a) vLLM 0.32 Qwen3-VL
   accepting embeds with DeepStack levels (unverified: only a community thread [web: discuss.vllm.ai/t/2294]) and
   (b) a statistical gate like G-ST1. Until then, C6 (`n=8`) is the cheap partial fix: it cuts hashing and processing
   from 256 to 32 requests.

6. **major - C0 Step 1, G-SP4: the sync counter undercounts, and its scope includes vLLM.**
   Problem: (a) `warnings.catch_warnings(record=True)` alone keeps Python's "default" action. A warning from the same
   code location is recorded only once, so a loop of 1,000 `.item()` calls at one line counts as 1. (b)
   `set_sync_debug_mode` is process-global, and vLLM runs in-process in colocate mode. Counting across `generate`
   includes vLLM's per-decode-step D2H reads of sampled tokens (one per decode step, about 1,300-1,500 per step at the
   R5 longest answer) [estimate]. Then "< 200 per step" is unreachable whatever C4 does. (c) The baseline "900-1,250"
   is an estimate, not a C0 number.
   Fix: In `count_syncs()` call `warnings.simplefilter("always")` inside the `catch_warnings` block. Restrict G-SP4 to
   the non-generate ranges ("batch minus generate" plus the update). Restate the target as "< 200 and <= 25 % of the
   C0-measured trainer-side count". Validate the counter once on the box against nsys `cudaStreamSynchronize` +
   `cudaMemcpy*` (sync) API counts for the same range. Report the vLLM-range count separately, for information only.

7. **major - G-SP5 (H2D < 50 MB per step): unreachable as scoped, because vLLM copies its own pixels.**
   Problem: vLLM's preprocessing output for 32 distinct images is about 32 x 1,024 x 1,536 x 4 B = 201 MB fp32 (100 MB if
   bf16), copied H2D each step inside `generate` [estimate; B_gpu_audit row 1d]. A whole-step H2D count therefore stays
   above 50 MB even with C1 perfect.
   Fix: Scope G-SP5 to trainer ranges ("prepare" + forward/backward): "trainer-side H2D < 50 MB". Report vLLM H2D
   separately. If C6b (item 5) is adopted, add "vLLM H2D < 10 MB".

8. **major - C5 Step 2, §8.1 row C5 mode A (1.45-1.6 s): the R4 evidence against it is not reconciled.**
   Problem: In R4, removing wake and sleep and shortening sync cut handovers by about 0.64 s (asleep 1.28 + 0.12 + 0.25
   vs awake 1.01) [measured: r4 note]. Yet the step fell only 0.05 s. Generate minus handovers rose from 8.83 s to
   9.22 s (+0.39 s) [estimate: from r4 table]. C.md blames all of R4's small gain on the 700-call push. But part of the
   handover saving vanished elsewhere (longest-answer noise, or slower awake decode under memory pressure); the cause is
   unknown. Also, `time/wake_s` 0.106 is the mean of two calls, and one is nested inside `sync_s` [SP/B_gpu_audit.md
   row 1c]. So "sync + wake + sleep" double-counts about 0.1 s.
   Fix: In §8.1, set mode A to "sync 1.27 + kv-wake 0.11 + sleep 0.29 - new sync; R4 shows an unexplained +0.39 s
   residual-decode rise when awake, so the low end is about 1.0 s". Re-derive the G-SP1 bound from that. In the C5 A/B,
   log generate minus handovers and the longest answer per step. Compare residual decode using the 4.3 ms/token fit, so
   the awake decode cost is measured and not assumed.

9. **major - §6 protocol and C0 Step 3: the speed baseline is the same arm that runs nsys on steps 4-5.**
   Problem: §8.2 runs "C0 baseline arm with nsys capture on steps 4-5". §6 takes means over steps 3-15 of that baseline
   arm. `nsys -t cuda,nvtx,osrt --cuda-memory-usage=true` slows the captured steps, so the baseline is inflated and G-SP1
   passes too easily.
   Fix: Capture after the timed window: run 17 steps with `--profile-steps 16-17` in both the baseline and the fast arm.
   Means stay steps 3-15. Profile both arms the same way, because G-SP3 and G-SP6 compare traces. Add 2 x about 1 min
   to §8.2.

10. **major - G-SP6 (GPU-active >= 70 %) and G-SP3: no measurement method, and the fallback path cannot produce
    them.**
    Problem: "GPU-active share" is not defined (union of kernel intervals? per stream?). The nsys-absent fallback
    (`g_probe --trace`, a torch.profiler trace of one step) gives no `nsys stats` reports. nsys presence on the
    `vllm_native` image is [unverified].
    Fix: Define GPU-active share = (union of all kernel and memcpy intervals on the device over the step) / step wall
    time, computed by a small script in `experiments/phase2/` from the `cuda_gpu_trace` CSV (nsys) or from the Chrome
    trace (torch.profiler). The same script serves both paths, with a Mac unit test on a synthetic trace. Add to box
    setup: `apt-get install -y nsight-systems-cli` from NVIDIA's apt repo, or bake nsys into the Dockerfile before the
    rental. That makes the fallback rare.

11. **major - G-EQ2, G-EQ3: "50 recorded batches" is not reachable in the planned box session.**
    Problem: The box arms are 15 steps, and `--check-rope-every` / `--check-meta-every` default to 50. So the box
    checks only step 1: one real batch. The Mac recordings come from a tiny model whose answer lengths and runaways do
    not look like real ones.
    Fix: (a) Run the fast arm with `--check-rope-every 1 --check-meta-every 1` (cost: one HF `get_rope_index` pass per
    step, which is CPU work but only in the check arm), giving 15 real batches. (b) Build 50+ Mac batches with the real
    Qwen3-VL processor (CPU is enough for ids, rope and meta) from Phase A's saved answers (final_eval answers on
    shifted and rescaled views, including length-1,536 truncated rows). (c) Add `--record-batches DIR` to the C0
    baseline arm so real ids, masks and grids come home. Restate both gates as ">= 50 batches, of which >= 15 come from
    the box".

12. **major - G-SP1 vs §8.1: the gate counts C6's 0-1 s, but the gated arm is C1-C5 only.**
    Problem: The C1-C6 range in §8.1 (mode A 2.7-5.6 s) includes C6's 0-1 s. G-SP1 applies to the "fast arm with C1-C5"
    (C6 is a separate optional arm, default off). For C1-C5, mode A is 2.7-4.6 s (midpoint 3.65 -> 19.78 -> 20.0 s) and
    mode B is 1.55-4.0 s (midpoint 2.78 -> 20.65 -> 21.0 s).
    Fix: Recompute G-SP1 from C1-C5 only: mode A <= 20.0 s, mode B <= 21.0 s. The fail thresholds stay at the low ends
    (20.7 / 21.9 s; after item 8 the mode A low end moves). Keep a separate C6 increment in G-ST1.

13. **major - C1 Step 2 + C3: after C3, the 256-row `pixel_values` (1.61 GB) is still built on the GPU and then
    discarded.**
    Problem: C1 assembles the 256-row output by `index_select` into a 1.61 GB CUDA tensor. TRL then holds it in
    `forward_kwargs` and in the batch until C3 replaces it. So C3's memory claim ("frees 1.6-3.2 GB") holds only after
    the reference pass. The GPU write of 1.61 GB plus `split_pixel_values_by_grid` work happen every step for nothing.
    This also eats the memory headroom mode A needs (G-MEM).
    Fix: With `--encoder-once` on, make `ImageMemo` return the 32-image patch tensor (201 MB) and a per-row image-index
    list. TRL's InternVL check reads only `pixel_values.ndim == 4` [code: grpo_trainer.py:2664], and a 2-D tensor passes
    it. Have the patched `split/unsplit_pixel_values_by_grid` pass it through. Measure
    `torch.cuda.max_memory_allocated()` both ways in G-MEM.

14. **minor - C5 Tests / G-SP2 mode B: the gate does not separate a good result from a poor one.**
    Problem: "Mode B: `time/sync_s` < 1.01 s". Mode B's baseline `sync_s` is 1.27 s (it includes the nested weights
    wake). A one-call push that saves only 0.27 s passes. Mode B's saving in §8.1 (0.3-1.0 s) is [unverified].
    Fix: Time one `load_weights` call on the box first (`c_speed_check --sync`), then set G-SP2 B to "<= measured
    one-call push + weights-wake + 0.1 s", expected about 0.3-0.4 s [estimate]. Write the measured value into §8.1
    before the A/B.

15. **minor - C1 Step 4: a silent CPU fallback hides a failed GPU path.**
    Problem: If neither backend passes the start-up self-check, training continues on the CPU processor and only logs
    `gpu_pixels/backend=cpu`. That conflicts with "GPU everywhere" and makes the speed arm measure the old path without
    saying so. The self-check also uses one image only.
    Fix: Under `--gpu-pixels auto` in a speed or B arm, raise at start-up unless `--allow-cpu-pixels` is given. Run the
    self-check on 8 images (centred, shifted, rescaled, and one digit image built in memory).

16. **minor - C4 Step 6: copying TRL's `_compute_loss` (265 lines) conflicts with "prefer a library; simple code".**
    Problem: A private copy plus a source-hash guard is code we own. It goes stale on every TRL bump.
    Fix: Keep the ">= 0.1 s" adoption rule. Also add: before copying, open an upstream TRL issue/PR proposing that the
    metric `.item()` calls be deferred to `log`. Prefer a 5-line monkeypatch of `self._metrics[mode][k].append` that
    accepts device tensors (read once in `log`), if TRL's code allows it without the copy. Record the choice in the C10
    ADR.

17. **minor - §1 Goal vs C6: "No task in C1-C6 may change one sampled token" contradicts C6 step 4 ("Not bit-equal").**
    Fix: Change to "No task in C1-C5 ...; C6 is a statistical arm (G-ST1)".

18. **minor - G-EQ6: the values are mislabelled and the count is wrong.**
    Problem: 5 metrics are named and 4 values are given. In R5 step 1, `reward` is about -3e-8 (advantage-centred, so
    equality on it means nothing). 0.5467 is `rewards/iou_tc/mean`, not reward. Grad norm is 0.003264
    [measured: runs/r5/shared_rl.json step 1].
    Fix: "IoU_tc 0.5467, legal 0.8984, mean answer length 117.89, grad norm 0.003264 (4 significant digits), plus
    `sampling/sampling_logp_difference/mean` equal to the baseline arm on the same host."

19. **minor - C5 Steps numbering and cross-references.**
    Problem: Steps run 1-6, then "5. Rejected". Tests say "the step-3 parameter equality", but it is step 5.
    Fix: Renumber "Rejected" as step 7. Change the Tests reference to "step 5".

20. **minor - §7 Risks: no row for a model change.**
    Problem: The owner is open to any model up to 4B. C2 (M-RoPE template), C3 (DeepStack outputs), C5 (196 tensors)
    and the K3 kernels are specific to Qwen3-VL. M1-M3 may change the model.
    Fix: Add a risk row: "Model changes (M1-M3)". Signal: the start-up self-checks fail or `RopeCache.install`/
    `EncoderCache.install` return None. Mitigation: each wrapper returns None for other families, so TRL's path is used.
    Rollback: none needed. Gate: re-run G-EQ1..5 on the new model.

21. **minor - §3.3 PNG row and C1 Step 1: smaller transfer and allocation details.**
    Problem: (a) The 25.2 MB figure assumes 3 channels. CT images are grey, so a 1-channel plane is 8.4 MB, with the RGB
    expand done on the device (the T1 view-bank row already says this). (b) `torch.empty(..., pin_memory=True)` every
    step calls cudaHostAlloc each time (ms-scale, serializing). (c) The nvImageCodec claim is right (GPU codecs: nvJPEG,
    nvJPEG2000, nvTIFF only) [web: developer.nvidia.com/nvimagecodec].
    Fix: Use one persistent pinned buffer (sized for 32 images), reused after a CUDA event confirms the previous copy
    finished. Send grey planes when the source is grey (`pil.mode == "L"`) and expand on the device. Keep the PNG-on-CPU
    justification with the web tag.

22. **minor - C6 Interface: `expand_outputs` is not needed; TRL already flattens `n>1` outputs.**
    Problem: TRL's colocate path builds completions as `[o.token_ids for outputs in all_outputs for o in
    outputs.outputs]` and discards `all_prompt_ids` [code: vllm_generation.py:745-746; grpo_trainer.py:1844]. TRL's own
    server mode already dedupes with `all_prompts[::num_generations]` [code: vllm_generation.py:618].
    Fix: The wrapper sends `prompts[::G]` with `sampling_params.clone()` and `n=G`, and returns the 32
    `RequestOutput`s as they are. Keep only `group_requests` with its run-length check. Adapt the `one_pad` check to
    compare `prompts[::G]`. C6 step 2: `max_num_seqs` is known from code
    (`per_device_train_batch_size x TP x steps_per_generation` = 16 x 1 x 16 = 256 [code: grpo_trainer.py:1186-1188;
    configs.py:185]). Note also TRL's fixed `max_num_batched_tokens=4096` [code: vllm_generation.py:368].

23. **minor - §4 order / missing task: the largest idle item (decode tail 4-5 s) has no Phase C task except C7.**
    Problem: B_gpu_audit names the decode tail (1-2 sequences decoding alone) as the largest idle source. C.md defers
    XGrammar and speculative decoding and leaves only the owner-gated limit (C7).
    Fix: Add a C0 output, "tail seconds per step" (time with fewer than 8 active sequences, taken from vLLM's per-step
    scheduler stats or the NVTX trace). Add one line to §8.3 handing the tail to Phase B/D with the measured number.
    Truncated answers are masked (`mask_truncated_completions` True, configs.py:188), so every token past the limit is
    pure waste. That makes the owner's C7 decision a data-backed one.

24. **minor - C10: the deletion list misses pieces C4 replaces, and the rollback text goes stale.**
    Problem: C10 lists "device reads replaced in C4" without naming them. After C10 removes `--gpu-pixels` etc., the §7
    rollback column ("`--gpu-pixels off`") is no longer valid.
    Fix: Name the deletions: `window`'s device `argmax/sum().tolist()` (fast.py:41-42), `by_length` /
    `by_group_length` device sorts (fast.py:379-381, 394-397), `prompt_groups`' `same.tolist()` (fast.py:74),
    `lengths.tolist()` (fast.py:125), `grid.prod(-1).tolist()` (fast.py:111), and the `packed_chunk` per-row
    `torch.full`/`torch.tensor` loops (fast.py:175-204). Add a §7 note: "after C10, rollback = owner runs
    `git revert <sha>`".

# Packed shared-prefix attention, async and graphs: research 2026-11-08

Follows `shared_prefix_2026-11-08.md` (R5: 31.0 -> 23.6 s/step). Not repeated here: beta=0, sleep mode, weight
sync, ImageMemo, n=8 requests (see `grpo_utilization_research`, `rl_code_optimality_review`, `r4_speed`).
**[V]** = read in source or a primary page. **[I]** = my estimate, not measured. No repo code was changed.

Step budget after R5 [I]: 11 s generation + 1.3 s sync leaves about 11 s for processor, ref and update. The "16 s
training" figure fits the base arm, not the packed one. From the 0.361 token ratio: update about 6-7 s, ref about 1.5 s.

## 1. Attention path for the packed row on sm_120

Sizes: a packed row is W ≈ 324 + 8 x 131 = 1,372 tokens, 2 rows per micro-batch, 32 rows per step. Pairs needed:
prompt causal 52k + 8 x (131 x 324 + 131²/2) = 461k, against W² = 1.88M that a dense-mask kernel computes (**24.5 %**).

| Path | sm_120 today | Qwen3-VL (transformers 5.19/5.20) | Ckpt + LoRA | Tiles computed | Est. saving vs (a) |
|---|---|---|---|---|---|
| (a) SDPA + 4-D bool mask (now) | yes, but flash refuses masks, so it runs mem-efficient; transformers then also runs `repeat_kv` (GQA only without a mask) [V `integrations/sdpa_attention.py` `use_gqa_in_sdpa`] | yes | yes | 100 % | — |
| (b) FlexAttention + BlockMask | yes. Inductor has sm12x defaults; bwd for head_dim 128 is 64x64 tiles with **1 pipeline stage** (99 KB SMEM) [V torch `template_heuristics/triton.py`, present at v2.13.0]. Flex reaches 90 % fwd / 85 % bwd of FA2 on A100 [V PyTorch blog] | **`flex_attention` is refused**: `_supports_flex_attn` is False, on main too [V `modeling_utils.py` `_flex_attn_can_dispatch`]. Works only as a custom `AttentionInterface` name (`_supports_attention_backend=True`). A 4-D BlockMask passes `_preprocess_mask_arguments` unchanged [V] | yes (HOP autograd) | ≈45-50 % (128-blocks, answers ≈ 1 block, misaligned) [I] | 0.3-0.6 s; recompiles on new widths; eager fallback is very slow |
| (c) Flash varlen, trunk+branch | **yes, no new dependency**: `torch.nn.attention.varlen.varlen_attn` (torch ≥ 2.12: separate `cu_seq_q`/`cu_seq_k`, `window_size=(-1,0)` causal, `enable_gqa`, autograd) [V v2.13.0 source]. cuDNN is chosen only on sm90/sm100, so sm_120 gets the flash kernel [V `_should_use_cudnn`]. F.sdpa-flash on a 5090: 187 TFLOPS fwd [V gau-nernst] | custom `AttentionInterface` function | yes, if boundaries travel as forward kwargs (non-reentrant recompute replays them) | ≈ exact (25-30 %) | **0.7-1.1 s** [I] |

(d) A prefix pass, then an answer pass with KV reuse (Prefix Grouper: two calls, prefix KV joined by
`batch_repeat_cat` [V arXiv 2506.05433]) is (c) with one extra call per layer. Cascade attention (FlashInfer) and
Hydragen are serving/decode kernels with no backward pass, so they are not usable for training.

The arithmetic behind (c) [I]. Dense attention costs 32 rows x 28 layers x 4·W²·128·16 (15.4 GFLOP) x 5.5 fwd-units
(ref 1 + fwd 1 + recompute 1 + bwd 2.5) ≈ 76 TFLOP per step. At 60-90 TFLOPS (mem-efficient kernel with a bias) that
is 0.85-1.3 s. Exact varlen at ~30 % of the tiles and ~120 TFLOPS is about 0.2 s. The linear layers (≈ 490 TFLOP)
still dominate. So attention is about 9 % of the training time: a real gain, but not the largest one.

**Recommended design** (same layout as ours, verified by verl RFC #8299 / PR #8300 on H100): flatten the packed rows.
Make one varlen call per layer. Query segments: the trunk (prompt minus its last token, causal over itself), then each
branch `[p_last, a_k]`. Key segments: trunk for the trunk; `[trunk ; branch_k]` for each branch, built by a K/V gather
after RoPE. Flash's causal rule is bottom-right aligned when keys outnumber queries, so branch token i sees the whole
trunk plus branch tokens ≤ i. In the backward pass, the gather sums each branch's gradient into the trunk [V RFC
8299 §3.2]. Build it with `expand`+`cat`: `index_select`'s backward uses atomic bf16 adds [I]. With
`attention_mask=None` and a custom name, `create_causal_mask` builds no mask [V `masking_utils.py` 800-806]. Count
the calls per forward and raise if a layer skipped the function (verl does this). The K/V gather is about 16 MB per
group per layer. The n²/PAIR term in `packing_wins` goes away, along with the W×W mask and the `repeat_kv` copies.

## 2. What production frameworks do

- **verl**: PrefixGrouper merged (#4368, 2026-01-05). It forked the actor, needs `use_remove_padding=False`, and has
  had no caller since the engine refactor. Open: #6689 (Megatron, MAGI/Flex tree masks, 5.5k lines); #8300 (FSDP,
  **our exact layout with the duplicated last prompt token**, varlen flash + K/V gather; Qwen3-4B with a 2k prompt: step
  2.7x, actor 3.1x, ref 3.9x; with 83-token prompts 1.18x). #8300 rejects VLMs. [V gh]
- **Tree Training** (arXiv 2511.00413): a custom FA3-based kernel (Hopper) on Megatron; 6.2-6.3x on agentic trees;
  not evaluated on single-prompt GRPO [V]. MiniMax Forge "prefix tree merging" and AReaL-DTA are cited in verl
  #6401/#7561.
- **TRL**: AsyncGRPO trie packing, text-only (earlier note). OpenRLHF, slime, ROLL, NeMo-RL: I found no shared-prefix
  path in an issue/PR search; they pack with varlen and repeat the prompt [V gh search, absence only].
- **Verdict**: no better layout is known for one-level sharing. Ours equals verl #8300's layout. The only upgrade is
  the attention kernel (dense mask -> varlen with K/V gather).

## 3. CUDA graphs vs async overlap

- **Graphs**: possible only with packed-token buckets (for example multiples of 256), per-layer
  `make_graphed_callables`, non-reentrant checkpoint inside the capture, and `max_q` fixed per bucket. Blockers: the
  graph conflict with vLLM's cumem allocator already seen (P1 reduce-overhead crash), one graph pool per bucket, and
  TRL's `.item()` syncs. The gain is only launch overhead: about 0.5-1 s [I]. **Not now.**
- **Async (one-step off-policy)**: the ceiling is the GPU's idle time. That is about 30 % of 23.6 s ≈ 7 s; realistic
  3-5 s [I]. TRL has no colocate async, and AsyncGRPO is text-only. verl's +23-40 % uses separate GPU pools. On one
  GPU it needs vLLM awake (R4: costs nothing), a generation thread or second process, and stale samples (TIS corrects
  them; C10's old=current shortcut then clips against θ_t, not the sampler). Larger ceiling, high effort and risk:
  **after** the items below.

## 4. Ranked next optimizations

| # | Change | Est. s/step | Effort | Risk |
|---|---|---|---|---|
| 1 | **Checkpointing off for packed micro-batches.** Packing cut tokens per micro-batch to ≈ 2.7k (P1's no-ckpt OOM at 30.1 GB was 16 padded rows, ≈ 13-30k tokens). Activations ≈ 2.2 MB/token [I] ≈ 6 GB. Toggle the per-layer `gradient_checkpointing` flag per micro-batch on a packed-token budget (≈ 6k), so rare long groups keep checkpointing. Removes the recompute, about 30 % of update time. | **1.5-2** | ~15 lines + one memory arm | OOM on long groups (the budget guards it) |
| 2 | **Varlen trunk+branch attention** (§1 c) as a custom `AttentionInterface` function | **0.7-1.1** | ~80 lines; current parity tests apply | `varlen_attn` is new API; check flash backward on sm_120 first with a 30-second micro-benchmark against (a) |
| 3 | **Wider chunks**: 4 groups per micro-batch (32 x 8, the same update under dr_grpo) and 64-row no-grad ref chunks. Halves launches and the ~29 host syncs per micro-batch; bigger LoRA GEMMs | **0.5-1.5** | config + ref batch override | competes with #1 for memory: test {ckpt off, 2 groups} against {ckpt on, 4 groups} |

Expected after 1+2: ≈ 23.6 -> ~20-21 s [I]. Measure in one rental: arms base / #1 / #1+#3 / #1+#2.

## Sources

- transformers 5.20.1 installed: `masking_utils.py` 746-800, `modeling_utils.py` 1695, 1843-1880,
  `integrations/{sdpa,flex}_attention.py`, `models/qwen3_vl/modeling_qwen3_vl.py` 497-535, 610-625; main:
  https://github.com/huggingface/transformers/blob/main/src/transformers/models/qwen3_vl/modeling_qwen3_vl.py
- torch varlen: https://github.com/pytorch/pytorch/blob/v2.13.0/torch/nn/attention/varlen.py ; lower-right causal
  `torch/nn/attention/bias.py`; flex sm12x configs `torch/_inductor/template_heuristics/triton.py` (v2.13.0)
- FlexAttention blog: https://pytorch.org/blog/flexattention/ ; 5090 SDPA numbers: https://gau-nernst.github.io/fa-5090/
- FA4 on sm_120 (community, not upstream): https://github.com/sgl-project/sglang/discussions/10564 ,
  https://huggingface.co/SecondNatureComputing/flash-attn-4-sm120 ; FA2 sm_120 status:
  https://github.com/Dao-AILab/flash-attention/issues/1987
- Prefix Grouper: https://arxiv.org/abs/2506.05433 , https://github.com/johncaged/PrefixGrouper
- Tree Training: https://arxiv.org/abs/2511.00413 ; Hydragen: https://arxiv.org/abs/2402.05099 ;
  FlashInfer: https://github.com/flashinfer-ai/flashinfer
- verl: https://github.com/verl-project/verl/issues/8299 , /pull/8300 , /pull/4368 , /issues/6401 , /pull/6689 ,
  /issues/7561 ; one-step off-policy: https://verl.readthedocs.io/en/latest/advance/one_step_off.html
- Repo: `experiments/phase2/perf/p1_results.md` (no_ckpt OOM 30,146 MiB; graphs crash), `src/objdel/train/fast.py`

"""B5 / H-A5-1 + H-A5-2(a), step 2 (vLLM env: the 0.32.0 macOS CPU wheel, transformers 5.19.0, torch 2.13.0).

A standalone sampler as Phase 2 would run it: untokenized chat text (rendered by the training env, export_prompts.py)
plus the PIL image, mm_processor_kwargs = objdel.train.configs.processor_kwargs minus padding_side.
(1) prompt parity: RequestOutput.prompt_token_ids vs the training env's GRPO input_ids, every exported item;
(2) samples for the logprob check: the first S items of each dataset x n = 8, T = 1, top_p = 1, top_k = 0 (off),
    min_p = 0, max_tokens 96, logprobs = 0 (the sampled token's processed logprob), seed 1.
The model's generation_config is left at vLLM's default ("auto"), as a real run would; explicit SamplingParams must win.
macOS CPU wheel workaround (Qwen3.5 only, CPU-only packaging gap, irrelevant on CUDA): vllm._custom_ops defines
causal_conv1d_weight_pack but the macOS _C library lacks the op, so the ARM-bf16 prepack in
layers/utils.dispatch_cpu_unquantized_gemm dies at load; with VLLM_ENABLE_V1_MULTIPROCESSING=0 (engine in this
process) the wrapper is removed so vLLM keeps the unpacked conv weight (its own fallback path).
    VLLM_CPU_KVCACHE_SPACE=1 HF_HUB_OFFLINE=1 <vllm venv>/bin/python -I vllm_parity.py prompts.json <model> <dtype>
        <out.json> [S]
"""
import json
import os
import sys
import time

from PIL import Image
from vllm import LLM, SamplingParams


def _patch_conv_pack() -> bool:
    import torch
    import vllm._custom_ops as ops

    if hasattr(ops, "causal_conv1d_weight_pack") and not hasattr(torch.ops._C, "causal_conv1d_weight_pack"):
        del ops.causal_conv1d_weight_pack
        return True
    return False


def main(prompts: str, model: str, dtype: str, out: str, s: int) -> None:
    P = json.load(open(prompts))["models"][model]
    repo = P["repo"]
    res = {"model": model, "repo": repo, "dtype": dtype, "datasets": {}, "conv_pack_patched": _patch_conv_pack()}
    t0 = time.time()
    # one engine per model; SmolVLM's per-dataset size goes in each request's mm_processor_kwargs
    # profiling dummy image capped at 512 x 512 (our inputs are <= 512 px); Qwen's processor default (longest_edge
    # 16.8 M px) otherwise profiles a 16k-token encoder pass, which takes > 30 min on 2 CPU threads. Profiling only.
    llm = LLM(model=repo, dtype=dtype, max_model_len=2048, max_num_seqs=8,
              limit_mm_per_prompt={"image": {"count": 1, "width": 512, "height": 512}},
              logprobs_mode="processed_logprobs", seed=0)
    res["load_s"] = round(time.time() - t0, 1)
    for dataset, D in P["datasets"].items():
        mm = {k: v for k, v in D["processor_kwargs"].items() if k != "padding_side"}
        reqs = [{"prompt": x["text"], "multi_modal_data": {"image": Image.open(x["image"]).convert("RGB")},
                 **({"mm_processor_kwargs": mm} if mm else {})} for x in D["items"]]
        outs = llm.generate(reqs, SamplingParams(n=1, temperature=0.0, max_tokens=1), use_tqdm=False)
        parity = [{"item_id": x["item_id"], "hf_len": len(x["gen_ids"]), "vllm_len": len(o.prompt_token_ids),
                   "equal": list(o.prompt_token_ids) == x["gen_ids"]} for x, o in zip(D["items"], outs, strict=True)]
        sp = SamplingParams(n=8, temperature=1.0, top_p=1.0, top_k=0, min_p=0.0, repetition_penalty=1.0,
                            max_tokens=96, logprobs=0, seed=1)
        t0 = time.time()
        gens = llm.generate(reqs[:s], sp, use_tqdm=False)
        samples = [{"item_id": x["item_id"], "samples": [
            {"ids": list(c.token_ids), "lp": [next(iter(d.values())).logprob for d in c.logprobs],
             "finish": c.finish_reason, "text": c.text} for c in o.outputs]}
            for x, o in zip(D["items"][:s], gens, strict=True)]
        res["datasets"][dataset] = {"mm_processor_kwargs": mm, "sample_s": round(time.time() - t0, 1),
                                    "sampling_params": repr(sp), "parity": parity, "samples": samples}
        print(model, dataset, "prompt ids equal", sum(p["equal"] for p in parity), "/", len(parity),
              "sampled", sum(len(c["ids"]) for x in samples for c in x["samples"]), "tokens", flush=True)
    json.dump(res, open(out, "w"))


if __name__ == "__main__":
    os.environ.setdefault("VLLM_CPU_KVCACHE_SPACE", "1")
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], int(sys.argv[5]) if len(sys.argv) > 5 else 4)

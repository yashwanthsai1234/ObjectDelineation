"""A5 sandbox (vLLM 0.32.0 CPU wheel, macOS): SmolVLM-256M through vLLM directly, outside TRL.

(a) prompt-length parity: vLLM's processed prompt vs the HF processor's input_ids, longest_edge 512 and 1024
(b) TRL's current colocate hand-off (processor-expanded ids + image): what vLLM does with it (#6294)
(c) group of G = 8 at temperature 1 from one request (n=8) vs 8 duplicate requests: num_cached_tokens per child
(d) structured outputs (xgrammar, D10 regex): parse rate, and the logprobs vLLM returns under the constraint
Writes JSON to argv[2] for the HF teacher-forcing comparison (hf_check.py).
"""
import json
import os
import sys
import time

os.environ.setdefault("VLLM_CPU_KVCACHE_SPACE", "1")
from PIL import Image, ImageDraw  # noqa: E402
from transformers import AutoProcessor  # noqa: E402
from vllm import LLM, SamplingParams  # noqa: E402
from vllm.sampling_params import StructuredOutputsParams  # noqa: E402

REPO = "HuggingFaceTB/SmolVLM-256M-Instruct"
EDGE = int(os.environ.get("EDGE", "512"))
QUESTION = "Outline the nodule as a polygon. Answer with JSON only."


def image():
    img = Image.new("RGB", (512, 512), (40, 40, 40))
    ImageDraw.Draw(img).ellipse((200, 180, 320, 300), fill=(200, 200, 200))
    return img


def main():
    D = json.load(open(sys.argv[1]))
    OUT = sys.argv[2]
    res = {"edge": EDGE}
    proc = AutoProcessor.from_pretrained(REPO, size={"longest_edge": EDGE})
    msgs = [{"role": "user", "content": [{"type": "image"}, {"type": "text", "text": QUESTION}]}]
    text = proc.apply_chat_template(msgs, add_generation_prompt=True)
    img = image()
    hf_ids = proc(text=text, images=[img], return_tensors="np")["input_ids"][0].tolist()
    res["prompt_text"] = text
    res["hf_prompt_len"] = len(hf_ids)
    t0 = time.time()
    llm = LLM(model=REPO, dtype="float32", max_model_len=4096, max_num_seqs=16,
              mm_processor_kwargs={"size": {"longest_edge": EDGE}}, limit_mm_per_prompt={"image": 1},
              logprobs_mode="processed_logprobs", seed=0)
    res["load_s"] = round(time.time() - t0, 1)
    # (a) unexpanded text prompt + image (what PR #6718 sends, and what a standalone sampler sends)
    sp = SamplingParams(n=1, temperature=0.0, max_tokens=8)
    o = llm.generate([{"prompt": text, "multi_modal_data": {"image": img}}], sp, use_tqdm=False)[0]
    res["vllm_prompt_len_text"] = len(o.prompt_token_ids)
    res["prompt_ids_equal_hf"] = list(o.prompt_token_ids) == hf_ids
    # (b) TRL's current hand-off: expanded ids + image
    try:
        o = llm.generate([{"prompt_token_ids": hf_ids, "multi_modal_data": {"image": img}}], sp, use_tqdm=False)[0]
        res["expanded_handoff"] = {"prompt_len": len(o.prompt_token_ids), "text": o.outputs[0].text}
    except Exception as e:  # noqa: BLE001
        res["expanded_handoff"] = {"error": f"{type(e).__name__}: {str(e)[:300]}"}
    # (c) G = 8: one n=8 request vs 8 duplicates (fresh prompt each so the first is a cache miss)
    sp8 = SamplingParams(n=8, temperature=1.0, top_p=1.0, top_k=0, max_tokens=48, seed=1)
    llm.reset_prefix_cache()
    t0 = time.time()
    o = llm.generate([{"prompt": text, "multi_modal_data": {"image": img}}], sp8, use_tqdm=False)[0]
    res["n8"] = {"num_cached_tokens": o.num_cached_tokens, "s": round(time.time() - t0, 2),
                 "distinct": len({c.text for c in o.outputs})}
    llm.reset_prefix_cache()
    sp1 = SamplingParams(n=1, temperature=1.0, top_p=1.0, top_k=0, max_tokens=48)
    t0 = time.time()
    outs = llm.generate([{"prompt": text, "multi_modal_data": {"image": img}}] * 8, sp1, use_tqdm=False)
    res["dup8"] = {"num_cached_tokens": [x.num_cached_tokens for x in outs], "s": round(time.time() - t0, 2)}
    o = llm.generate([{"prompt": text, "multi_modal_data": {"image": img}}], sp1, use_tqdm=False)[0]
    res["repeat_after"] = {"num_cached_tokens": o.num_cached_tokens, "prompt_len": len(o.prompt_token_ids)}
    try:
        res["block_size"] = llm.llm_engine.vllm_config.cache_config.block_size
    except Exception as e:  # noqa: BLE001
        res["block_size"] = repr(e)[:100]
    # (d) structured outputs with the D10 regex (lidc schema); unconstrained twin for comparison
    pat = D["patterns"]["lidc"]
    for name, so in (("free", None), ("regex", StructuredOutputsParams(regex=pat))):
        spc = SamplingParams(n=8, temperature=1.0, top_p=1.0, top_k=0, max_tokens=160, logprobs=0, seed=2,
                             structured_outputs=so)
        t0 = time.time()
        o = llm.generate([{"prompt": text, "multi_modal_data": {"image": img}}], spc, use_tqdm=False)[0]
        samples = []
        for c in o.outputs:
            lp = [next(iter(d.values())).logprob for d in c.logprobs]
            samples.append({"text": c.text, "ids": list(c.token_ids), "lp": lp, "finish": c.finish_reason})
        res[f"so_{name}"] = {"s": round(time.time() - t0, 2), "samples": samples,
                             "prompt_ids": list(o.prompt_token_ids)}
    json.dump(res, open(OUT, "w"))
    print(json.dumps({k: v for k, v in res.items() if not k.startswith("so_")}, indent=1))
    for name in ("free", "regex"):
        s = res[f"so_{name}"]
        print(name, s["s"], [(x["finish"], len(x["ids"]), x["text"][:70]) for x in s["samples"][:4]])


if __name__ == "__main__":
    main()

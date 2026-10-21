"""Sample K completions per prompt with vLLM offline (plan D13, C-BRIDGE); runs in /opt/vllm, never in /opt/venv.

    /opt/vllm/bin/python -I experiments/phase2/sample_vllm.py --prompts P.jsonl --model DIR --out C.jsonl

Input (objdel.train.sample.export_prompts in the training env): one JSON object per prompt, {"id", "prompt" (chat
text rendered with tokenize=False), "image" (path), "processor_kwargs" (configs.processor_kwargs minus
padding_side, sent as mm_processor_kwargs), "sampling" (every SamplingParams field explicit: configs.sampling plus
seed and stop_token_ids), "n"}; each request samples with request_seed(seed, id). Output: one line per completion, {"id", "k", "text", "token_ids", "finish_reason"}
(+ "logprobs" with --logprobs: the sampled token's processed logprob, for the H-A5-2 teacher-forcing check), and
<out>.meta.json with every distinct SamplingParams as vLLM built it (the echo). vLLM's token_ids end with the stop
id when finish_reason is "stop" (B5 CPU runs: LFM id 7, SmolVLM id 49279), so the training env applies TRL's
truncation rule (objdel.train.rewards.truncated) to them in objdel.train.sample.import_completions.

Imports nothing from objdel (its numpy / numba / torch pins differ from vLLM 0.32's, B5 Q1); vLLM, torch and PIL
are imported inside main, so --help works anywhere. GPU sequencing (D13): the launcher (select_ckpt) has not
initialised CUDA; this process refuses to start unless torch.cuda.mem_get_info() reports >= 90 % free, then builds
LLM(gpu_memory_utilization=0.85). The profiling image is capped at 512 x 512 (our inputs are <= 512 px; Qwen's
default profile is a 16k-token encoder pass, B5).

References:
- vLLM offline inference (LLM.generate, multi_modal_data, mm_processor_kwargs, SamplingParams):
  https://docs.vllm.ai/en/latest/serving/offline_inference.html , https://docs.vllm.ai/en/latest/api/vllm/sampling_params.html
- experiments/phase2/lab/B5/vllm_parity.py and REPORT.md (prompt-id parity 64 / 64, logprob parity, CPU notes)
"""
import argparse
import hashlib
import json
import sys
from pathlib import Path

FREE = 0.90  # of total device memory, before the engine starts (D13)


def check_memory(mem_get_info, need: float = FREE) -> None:
    free, total = mem_get_info()
    if free < need * total:
        raise SystemExit(f"sample_vllm: {100 * free / total:.0f} % of GPU memory free, need {100 * need:.0f} % "
                         "(another process holds the GPU)")


def request_seed(seed: int, item_id: str) -> int:
    """One seed per request from (the exported base seed, the item id): one shared seed would correlate the sample
    streams across items; a hash keeps reruns identical (31 bits, a valid SamplingParams seed)."""
    return int.from_bytes(hashlib.sha256(f"{seed}/{item_id}".encode()).digest()[:4], "big") >> 1


def requests(lines: list[dict], SamplingParams, open_image, logprobs: bool = False) -> tuple[list, list]:
    """vLLM requests and one SamplingParams each, built from the exported fields only (the seed per item)."""
    reqs, params = [], []
    for x in lines:
        req = {"prompt": x["prompt"], "multi_modal_data": {"image": open_image(x["image"])}}
        if x["processor_kwargs"]:
            req["mm_processor_kwargs"] = x["processor_kwargs"]
        reqs.append(req)
        sp = x["sampling"] | {"seed": request_seed(x["sampling"]["seed"], x["id"])}
        params.append(SamplingParams(n=x["n"], **sp, **({"logprobs": 0} if logprobs else {})))
    return reqs, params


def outputs(lines: list[dict], results):
    for x, res in zip(lines, results, strict=True):
        for k, c in enumerate(res.outputs):
            rec = {"id": x["id"], "k": k, "text": c.text, "token_ids": list(c.token_ids),
                   "finish_reason": c.finish_reason}
            if c.logprobs:
                rec["logprobs"] = [next(iter(d.values())).logprob for d in c.logprobs]
            yield rec


def write(recs, path) -> None:
    with open(path, "w") as f:
        for r in recs:
            f.write(json.dumps(r) + "\n")


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--prompts", required=True, help="JSONL from objdel.train.sample.export_prompts")
    ap.add_argument("--model", required=True, help="merged checkpoint directory (with its processor)")
    ap.add_argument("--out", required=True, help="completions JSONL")
    ap.add_argument("--dtype", default="bfloat16")
    ap.add_argument("--gpu-memory-utilization", type=float, default=0.85)
    ap.add_argument("--max-model-len", type=int, default=2048, help="prompt (<= ~420 tokens) + max_tokens")
    ap.add_argument("--seed", type=int, default=0, help="engine seed (each request carries its own)")
    ap.add_argument("--logprobs", action="store_true", help="also write the sampled tokens' processed logprobs")
    a = ap.parse_args(argv)
    import torch
    check_memory(torch.cuda.mem_get_info)
    import vllm
    from PIL import Image

    lines = [json.loads(line) for line in Path(a.prompts).read_text().splitlines()]
    reqs, params = requests(lines, vllm.SamplingParams, lambda p: Image.open(p).convert("RGB"), a.logprobs)
    llm = vllm.LLM(model=a.model, dtype=a.dtype, gpu_memory_utilization=a.gpu_memory_utilization,
                   max_model_len=a.max_model_len, enable_prefix_caching=True, seed=a.seed,
                   limit_mm_per_prompt={"image": {"count": 1, "width": 512, "height": 512}},
                   **({"logprobs_mode": "processed_logprobs"} if a.logprobs else {}))
    write(outputs(lines, llm.generate(reqs, params, use_tqdm=False)), a.out)
    meta = {"vllm": vllm.__version__, "model": a.model, "prompts": len(lines),
            "sampling_params": sorted({repr(p) for p in params})}
    Path(a.out + ".meta.json").write_text(json.dumps(meta, indent=1) + "\n")
    print(json.dumps(meta), file=sys.stderr)


if __name__ == "__main__":
    main()

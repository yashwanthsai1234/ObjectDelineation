"""Sampling: one sampler interface for evaluation, the trainer's own generation and the JSONL bridge to the sampler
child. Completions go to the evaluation harness as rows {"item_id", "samples": [text, ...], "truncated": [bool, ...]}.

- `Sampler.sample(rows, params) -> [Sample]`; params_for(model, n, seed, stop_token_ids) is the one source of the
  sampling fields (objdel.policy.models.sampling), and each row's request seed is request_seed(seed, key).
- VllmSampler: vllm.LLM offline in this process (the vLLM and objdel venv is one). The caller has not initialised
  CUDA, and it refuses to start unless torch.cuda.mem_get_info() reports at least FREE of the device free. ENV is set
  before vLLM is imported: FlashInfer's top-k sampler compiles with nvcc, and vLLM's PyTorch sampler draws from the
  same distribution (top_k 0, top_p 1). With `logprobs`, records also hold the sampled tokens' processed logprobs and
  vLLM's prompt ids (parity.py reads both). With `adapter` (a LoRA folder over the ckpt base) the engine serves the
  adapter unmerged, as in training.
- The HF path (load, generate, batches, todo, main --kind hf): transformers generate on left-padded batches of
  `batch_items` prompts, every sampling field passed explicitly over a fresh GenerationConfig. Batch b is seeded
  seed + b * batch_items, so the same lines, model and library versions give the same samples.
- build_prompt: TRL GRPO's generation prompt as text and the row's image path.
- Truncation is the reward's rule: models.truncated on the stop ids (`stop_ids(path, proc)` for a checkpoint).
- export_prompts / import_completions: the JSONL bridge to the sampler child (prompt text, image path, processor
  kwargs minus padding_side, every sampling field and n; completions come back as token ids). export reads TRL's
  completions tables back into rows; write(rows, path) writes JSONL.
- main: the child-process entry, `python -I -m objdel.policy.sampler --kind vllm|hf --prompts P --model DIR --out C`
  (CUDA untouched in the parent; the caller sets a timeout, kills the process group and resumes).

References: TRL GRPOTrainer (https://github.com/huggingface/trl/blob/v1.17.1/trl/trainer/grpo_trainer.py),
vLLM offline inference (https://docs.vllm.ai/en/latest/serving/offline_inference.html), HF generate
(https://huggingface.co/docs/transformers/main_classes/text_generation).
"""

import argparse
import hashlib
import json
import os
import sys
import time
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from functools import cached_property
from pathlib import Path
from typing import Any, Protocol

import pandas as pd
from trl.data_utils import prepare_multimodal_messages

from objdel.core.io import read_jsonl
from objdel.policy import models
from objdel.policy.errors import Refusal
from objdel.verify.reward.spec import REWARD_V4

G = REWARD_V4.group


def export(paths: Iterable, keys: Iterable[str] | None = None) -> list[dict]:
    """Rows of TRL completions tables, in file and table order; ValueError unless the table is whole groups."""
    keys = None if keys is None else set(keys)
    rows = []
    for path in paths:
        df = pd.read_parquet(path, columns=["step", "completion", "key", "truncated", "group"])
        last = None
        for i in range(0, len(df), G):
            b = df.iloc[i:i + G]
            if len(b) != G or b["key"].nunique() != 1 or b["group"].nunique() != 1:
                raise ValueError(f"{path}: rows {i}..{i + len(b) - 1} are not one group of {G} completions")
            key = str(b["key"].iloc[0])
            if keys is not None and key not in keys:
                last = None
                continue
            texts = [models._text(c if isinstance(c, str) else list(c)) for c in b["completion"]]
            if last is None or last["item_id"] != key:
                last = {"item_id": key, "step": int(b["step"].iloc[0]), "samples": [], "truncated": []}
                rows.append(last)
            last["samples"] += texts
            last["truncated"] += [bool(t) for t in b["truncated"]]
    return rows


def write(rows: Iterable[dict], path) -> None:
    """Write rows as JSON lines to `path`, creating its folder."""
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        for r in rows:
            f.write(json.dumps(r) + "\n")


def _path(image) -> str:
    return image if isinstance(image, str) else image["path"]


def export_prompts(rows, path, model: str, n: int, processor=None, stop_ids=None, seed: int = 0,
                   max_tokens: int | None = None) -> None:
    """rows: prompt, image (a path; grpo_rows' decoded Image column is read as paths), key (unique) -> the prompts
    JSONL. processor: AutoProcessor of the sampled checkpoint (default: the entry's hub id with processor_kwargs).
    max_tokens: a completion cap below the model's budget (the CPU smoke), else the entry's max_completion."""
    kw = models.processor_kwargs(model)
    if processor is None:
        processor = models.processor(model)
    if hasattr(rows, "cast_column"):  # a datasets.Dataset: paths, not decoded images
        from datasets import Image
        rows = rows.cast_column("image", Image(decode=False))
    sampling = models.sampling(model) | {"seed": seed} | ({"max_tokens": max_tokens} if max_tokens else {})
    if stop_ids is not None:
        sampling["stop_token_ids"] = sorted(stop_ids)
    write(prompt_lines(rows, processor, kw, sampling, n), path)


def prompt_lines(rows, processor, processor_kwargs: dict, sampling: dict, n: int) -> list[dict]:
    """Prompt rows (prompt, image, key) -> the prompts JSONL lines; ValueError on a repeated key."""
    mm = {k: v for k, v in processor_kwargs.items() if k != "padding_side"}
    seen, out = set(), []
    for r in rows:
        if r["key"] in seen:
            raise ValueError(f"{r['key']}: exported twice (ids are item keys)")
        seen.add(r["key"])
        text, (image,) = build_prompt(r, processor)
        out.append({"id": r["key"], "prompt": text, "image": image, "processor_kwargs": mm, "sampling": sampling,
                    "n": n})
    return out


def build_prompt(row, processor) -> tuple[str, list]:
    """TRL GRPO's generation prompt as text and the row's image paths (module docstring)."""
    image = _path(row["image"])
    text = processor.apply_chat_template([prepare_multimodal_messages(row["prompt"], images=[image])],
                                         add_generation_prompt=True, tokenize=False)
    return (text[0] if isinstance(text, list) else text), [image]


def import_completions(path, stop_ids, n: int = G, ids: Iterable[str] | None = None) -> list[dict]:
    """The sampler child's completions -> eval_run rows {"item_id", "samples", "truncated", "tokens"} (tokens:
    completion length in token ids, for verify.metrics' token stats; samples in k order, items in first-seen order).
    ValueError unless every id has exactly the exported n samples k = 0 .. n-1, each once, and (with ids, the
    exported keys) the ids are exactly those (a cut or doubled file would change K silently)."""
    by: dict[str, dict[int, dict]] = {}
    for c in read_jsonl(path):
        if by.setdefault(c["id"], {}).setdefault(c["k"], c) is not c:
            raise ValueError(f"{c['id']}: sample k = {c['k']} twice")
    if ids is not None and set(by) != set(ids):
        raise ValueError(f"{path}: missing: {sorted(set(ids) - set(by))[:5]}, not exported: "
                         f"{sorted(set(by) - set(ids))[:5]}")
    rows = []
    for key, cs in by.items():
        if sorted(cs) != list(range(n)):
            raise ValueError(f"{key}: {len(cs)} samples k = {sorted(cs)[:3]}.., want {n} samples k = 0 .. {n - 1}")
        cs = [cs[k] for k in range(n)]
        rows.append({"item_id": key, "samples": [c["text"] for c in cs],
                     "truncated": [models.truncated(c["token_ids"], stop_ids) for c in cs],
                     "tokens": [len(c["token_ids"]) for c in cs]})
    return rows


def stop_ids(path: Path, proc) -> set[int]:
    """models.stop_ids of a checkpoint directory. Base weights may have no generation_config.json (eos only in
    text_config): then the model config's, as from_pretrained falls back for a loaded model."""
    from types import SimpleNamespace

    from transformers import AutoConfig, GenerationConfig

    try:
        gen = GenerationConfig.from_pretrained(path)
    except OSError:
        gen = GenerationConfig.from_model_config(AutoConfig.from_pretrained(path))  # reads text_config's eos too
    return models.stop_ids(SimpleNamespace(generation_config=gen), proc)


def params_for(model: str, *, n: int, seed: int, stop_token_ids: Sequence[int]) -> dict:
    """Every vLLM SamplingParams field of one sampling (module docstring). No vLLM import; HfSampler maps it with
    models.as_names(..., "hf")."""
    return models.sampling(model, "vllm") | {"n": n, "seed": seed, "stop_token_ids": list(stop_token_ids)}


@dataclass(frozen=True)
class Sample:
    """One prompt row's samples in k order; `truncated` is models.truncated on the stop ids (the reward's rule)."""

    item_id: str
    texts: tuple[str, ...]
    token_ids: tuple[tuple[int, ...], ...]
    truncated: tuple[bool, ...]


class Sampler(Protocol):
    """Rows (prompt, image path, key) and params_for's dict -> one Sample per row, in row order."""

    def sample(self, rows: Sequence[Mapping], params: Mapping[str, Any]) -> list[Sample]:
        """One Sample per row, in row order."""


def _rows_lines(rows: Sequence[Mapping], processor, model: str, params: Mapping[str, Any]) -> list[dict]:
    """Prompt rows of one dataset (the first key's first part) and params_for's dict -> prompts JSONL lines."""
    kw = models.processor_kwargs(model) if rows else {}
    return prompt_lines(rows, processor, kw, {k: v for k, v in params.items() if k != "n"}, params["n"])


def samples(lines: list[dict], recs: Iterable[dict], stop) -> list[Sample]:
    """Output records ({"id", "k", "text", "token_ids"}) -> one Sample per prompt line, in line order."""
    by: dict[str, dict[int, dict]] = {}
    for c in recs:
        by.setdefault(c["id"], {})[c["k"]] = c
    out = []
    for x in lines:
        cs = [by[x["id"]][k] for k in range(x["n"])]
        out.append(Sample(x["id"], tuple(c["text"] for c in cs), tuple(tuple(c["token_ids"]) for c in cs),
                          tuple(models.truncated(c["token_ids"], stop) for c in cs)))
    return out


FREE = 0.90  # of total device memory, before the vLLM engine starts

ENV = {"VLLM_USE_FLASHINFER_SAMPLER": "0"}  # the PyTorch sampler: no JIT build (module docstring)


def check_memory(mem_get_info, need: float = FREE) -> None:
    """Refusal unless mem_get_info() reports at least `need` of the device free (D13: no other process on the GPU)."""
    free, total = mem_get_info()
    if free < need * total:
        raise Refusal(f"sampler: {100 * free / total:.0f} % of GPU memory free, need {100 * need:.0f} % "
                      "(another process holds the GPU)")


def request_seed(seed: int, item_id: str) -> int:
    """One seed per request from (the base seed, the item id): one shared seed would correlate the sample streams
    across items; a hash keeps reruns identical (31 bits, a valid SamplingParams seed)."""
    return int.from_bytes(hashlib.sha256(f"{seed}/{item_id}".encode()).digest()[:4], "big") >> 1


def requests(lines: list[dict], SamplingParams, open_image, logprobs: bool = False) -> tuple[list, list]:
    """vLLM requests and one SamplingParams each, built from the prompt lines' fields only (the seed per item)."""
    extra = {"logprobs": 0} if logprobs else {}
    reqs = [{"prompt": x["prompt"], "multi_modal_data": {"image": open_image(x["image"])}}
            | ({"mm_processor_kwargs": x["processor_kwargs"]} if x["processor_kwargs"] else {}) for x in lines]
    params = [SamplingParams(n=x["n"], **x["sampling"] | {"seed": request_seed(x["sampling"]["seed"], x["id"])},
                             **extra) for x in lines]  # fmt: skip
    return reqs, params


def outputs(lines: list[dict], results):
    """One output record per completion of vLLM's results, in prompt-line and k order (module docstring)."""
    for x, res in zip(lines, results, strict=True):
        for k, c in enumerate(res.outputs):
            rec = {"id": x["id"], "k": k, "text": c.text, "token_ids": list(c.token_ids),
                   "finish_reason": c.finish_reason}
            if c.logprobs:  # parity runs (--logprobs): the sampled token's entry, and the prompt ids vLLM built
                rec["logprobs"] = [d[t].logprob for t, d in zip(c.token_ids, c.logprobs, strict=True)]
                rec["prompt_token_ids"] = list(res.prompt_token_ids)
            yield rec


def _rgb(path):
    from PIL import Image

    return Image.open(path).convert("RGB")


class VllmSampler(Sampler):
    """vllm.LLM offline, in this process (module docstring). `model` (a models.MODELS key) is needed by sample() and
    with an adapter (its LoRA rank); records() samples prompt lines as they are. `sampling_params` collects repr of
    every SamplingParams built. The engine's memory share and length default to the model's vllm_eval (the CLI
    passes both and no model)."""

    def __init__(self, ckpt: Path, model: str | None = None, gpu_memory_utilization: float | None = None,
                 max_model_len: int | None = None, *, adapter: Path | None = None, dtype: str = "bfloat16",
                 seed: int = 0, logprobs: bool = False):
        import torch

        check_memory(torch.cuda.mem_get_info)
        os.environ.update(ENV)  # before vLLM reads its environment
        import vllm

        self.vllm, self.ckpt, self.model, self.logprobs, self.sampling_params = vllm, ckpt, model, logprobs, set()
        mem, length = models.get(model).vllm_eval if model else (None, None)  # explicit values win
        self.lora = None if adapter is None else vllm.lora.request.LoRARequest("adapter", 1, str(adapter))
        self.llm = vllm.LLM(model=str(ckpt), dtype=dtype, gpu_memory_utilization=gpu_memory_utilization or mem,
                            max_model_len=max_model_len or length, enable_prefix_caching=True, seed=seed,
                            limit_mm_per_prompt={"image": {"count": 1, "width": 512, "height": 512}},
                            **({"logprobs_mode": "processed_logprobs"} if logprobs else {}),
                            **({"enable_lora": True, "max_lora_rank": models.get(model).lora_r}
                               if adapter is not None else {}))

    @cached_property
    def processor(self):
        """The checkpoint's processor (build_prompt's chat template)."""
        from transformers import AutoProcessor

        return AutoProcessor.from_pretrained(self.ckpt)

    def records(self, lines: list[dict]) -> list[dict]:
        """vLLM's output records of prompt lines (module docstring)."""
        reqs, params = requests(lines, self.vllm.SamplingParams, _rgb, self.logprobs)
        self.sampling_params |= {repr(p) for p in params}
        return list(outputs(lines, self.llm.generate(reqs, params, use_tqdm=False, lora_request=self.lora)))

    def sample(self, rows: Sequence[Mapping], params: Mapping[str, Any]) -> list[Sample]:
        """One Sample per row, in row order (all rows in one generate call)."""
        lines = _rows_lines(rows, self.processor, self.model, params)
        return samples(lines, self.records(lines), set(params["stop_token_ids"]))


def load(path, device: str):
    """(model, processor): library generation defaults, never the card's; left padding for batches. CUDA samples in
    bf16 whatever the checkpoint holds (full FT saves fp32), the dtype of D13's vLLM path; CPU runs use fp32."""
    import torch
    from transformers import GenerationConfig

    dtype = torch.bfloat16 if device == "cuda" else torch.float32
    model, proc = models.load_pair(path, dtype=dtype, device=device, proc_kw={"padding_side": "left"})
    model.generation_config = GenerationConfig()
    if proc.tokenizer.padding_side != "left":
        raise Refusal(f"sampler: {path}: the tokenizer pads on the {proc.tokenizer.padding_side}")
    return model, proc


def records(rows: list[list[int]], prompts: list[dict], stop, decode) -> list[dict]:
    """Generated ids (n rows per prompt, in order) -> output records, each row cut after its own first stop id."""
    out, n = [], prompts[0]["n"]
    for j, row in enumerate(rows):
        cut = next((i + 1 for i, t in enumerate(row) if t in stop), None)
        toks = row[:cut] if cut else row
        out.append({"id": prompts[j // n]["id"], "k": j % n, "token_ids": toks,
                    "finish_reason": "stop" if cut else "length", "text": decode(toks)})  # fmt: skip
    return out


def generate(model, proc, dev, batch: list[tuple[int, dict]]) -> list[dict]:
    """One generate call on a batch of (line index, prompt line), seeded by its first index."""
    import torch

    ps = [p for _, p in batch]
    s, kw = ps[0]["sampling"], {k: v for k, v in ps[0]["processor_kwargs"].items() if k != "padding_side"}
    gen = models.as_names(s, "hf")
    x = proc(text=[p["prompt"] for p in ps], images=[[_rgb(p["image"])] for p in ps], padding=True,
             return_tensors="pt", **kw)
    torch.manual_seed(s["seed"] + batch[0][0])
    with torch.no_grad():
        ids = model.generate(**x.to(dev), **gen, num_return_sequences=ps[0]["n"],
                             pad_token_id=proc.tokenizer.pad_token_id)  # fmt: skip
    rows = ids[:, x["input_ids"].shape[1] :].tolist()
    return records(rows, ps, set(gen["eos_token_id"]), lambda t: proc.tokenizer.decode(t, skip_special_tokens=True))


def batches(lines: list[dict], batch_items: int) -> list[list[tuple[int, dict]]]:
    """(line index, line) batches of batch_items lines; Refusal unless all lines share sampling, n and kwargs."""
    if len({json.dumps([p["sampling"], p["n"], p["processor_kwargs"]], sort_keys=True) for p in lines}) > 1:
        raise Refusal("sampler: prompts differ in sampling, n or processor_kwargs; a batch needs one setting")
    indexed = list(enumerate(lines))
    return [indexed[i : i + batch_items] for i in range(0, len(lines), batch_items)]


class HfSampler(Sampler):
    """transformers generate (module docstring); the Mac reference and the CPU test path."""

    def __init__(self, ckpt: Path, model: str, device: str = "cpu", batch_items: int = 8):
        self.model, self.device, self.batch_items = model, device, batch_items
        self.net, self.processor = load(ckpt, device)

    def sample(self, rows: Sequence[Mapping], params: Mapping[str, Any]) -> list[Sample]:
        """One Sample per row, in row order (batch_items rows per generate call)."""
        lines = _rows_lines(rows, self.processor, self.model, params)
        recs = [c for b in batches(lines, self.batch_items) for c in generate(self.net, self.processor, self.device, b)]
        return samples(lines, recs, set(params["stop_token_ids"]))


def todo(out: Path, bs: list[list[tuple[int, dict]]]) -> list:
    """Resume: rewrite `out` with the lines of every whole finished batch and return the other batches."""
    have = {}
    for x in out.read_text().splitlines() if out.exists() else []:
        try:
            c = json.loads(x)
            have.setdefault((c["id"], c["k"]), x)
        except (ValueError, KeyError, TypeError):
            continue  # a line cut by a kill
    keys = [[(p["id"], k) for _, p in b for k in range(p["n"])] for b in bs]
    done = [all(key in have for key in ks) for ks in keys]
    out.write_text("".join(have[key] + "\n" for ks, d in zip(keys, done, strict=True) if d for key in ks))
    return [b for b, d in zip(bs, done, strict=True) if not d]


def _main_vllm(a, lines: list[dict]) -> None:
    s = VllmSampler(a.model, None, a.gpu_memory_utilization, a.max_model_len, dtype=a.dtype, seed=a.seed,
                    logprobs=a.logprobs)
    write(s.records(lines), a.out)
    meta = {"vllm": s.vllm.__version__, "model": str(a.model), "prompts": len(lines),
            "sampling_params": sorted(s.sampling_params)}
    Path(str(a.out) + ".meta.json").write_text(json.dumps(meta, indent=1) + "\n")
    print(json.dumps(meta), file=sys.stderr)


def _main_hf(a, lines: list[dict]) -> None:
    import torch

    left = todo(a.out, batches(lines, a.batch_items))
    t0 = time.perf_counter()
    dev = "cuda" if torch.cuda.is_available() else "cpu"
    model, proc = load(a.model, dev) if left else (None, None)
    t1 = time.perf_counter()
    recs = []
    with open(a.out, "a") as f:
        for b in left:
            got = generate(model, proc, dev, b)
            f.write("".join(json.dumps(c) + "\n" for c in got))  # whole batches: a kill cuts at most one line
            f.flush()
            recs += got
    if a.timing:
        gen_s, n = time.perf_counter() - t1, len(recs)
        a.timing.write_text(json.dumps({
            "items": sum(map(len, left)), "samples": n, "batch_items": a.batch_items, "load_s": t1 - t0,
            "generate_s": gen_s, "s_per_sample": gen_s / n if n else None,
            "tokens_mean": sum(len(c["token_ids"]) for c in recs) / n if n else None,
            "length_share": sum(c["finish_reason"] == "length" for c in recs) / n if n else None}) + "\n")  # fmt: skip


def main(argv=None) -> None:
    """Child-process entry (module docstring). Input: export_prompts' lines. Output: one record per completion (vLLM:
    written once at the end, with <out>.meta.json holding every distinct SamplingParams as vLLM built it; HF: appended
    per whole batch, and --out is resumed: lines of whole finished batches stay, a half-written batch or a cut last
    line is dropped and resampled with its own seed, so a resumed file equals an uninterrupted one). HF --timing
    writes generate seconds (model load apart) per item-sample."""
    ap = argparse.ArgumentParser(description="Sample K completions per prompt line with vLLM or HF generate.")
    ap.add_argument("--kind", choices=("vllm", "hf"), required=True)
    ap.add_argument("--prompts", type=Path, required=True, help="JSONL from objdel.policy.sampler.export_prompts")
    ap.add_argument("--model", type=Path, required=True, help="checkpoint directory (with its processor)")
    ap.add_argument("--out", type=Path, required=True, help="completions JSONL")
    ap.add_argument("--dtype", default="bfloat16", help="vllm")
    ap.add_argument("--gpu-memory-utilization", type=float, default=0.85, help="vllm")
    ap.add_argument("--max-model-len", type=int, default=2048, help="vllm: prompt (<= ~420 tokens) + max_tokens")
    ap.add_argument("--seed", type=int, default=0, help="vllm engine seed (each request carries its own)")
    ap.add_argument("--logprobs", action="store_true", help="vllm: also write the sampled tokens' processed logprobs")
    ap.add_argument("--batch_items", type=int, default=8, help="hf: prompts per generate call (x n samples each)")
    ap.add_argument("--timing", type=Path, help="hf: write generate seconds per item-sample here (R0 rate)")
    a = ap.parse_args(argv)
    try:
        (_main_vllm if a.kind == "vllm" else _main_hf)(a, read_jsonl(a.prompts))
    except Refusal as e:
        raise e.system_exit() from None


if __name__ == "__main__":
    main()

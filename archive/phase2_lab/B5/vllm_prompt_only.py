"""B5 / H-A5-1 without running the model (vLLM env: 0.32.0 macOS CPU wheel, transformers 5.19.0).

The vLLM macOS CPU wheel cannot run Qwen3.5 (its _C library lacks causal_conv1d_fwd_cpu / _weight_pack, used by the
Gated-DeltaNet layers; vllm_parity.py log). Prompt processing does not need the model: LLM.generate hands every
prompt to the engine's renderer (vllm.renderers.renderer_from_config(vllm_config).render_cmpl -> process_for_engine:
tokenize the text, parse the image, apply the HF processor through vLLM's multimodal processor, expand the
placeholders). This script builds the same VllmConfig (EngineArgs.create_engine_config, no weights loaded) and renderer
and compares its prompt_token_ids with the training env's GRPO input_ids (export_prompts.py), for every exported item.
    HF_HUB_OFFLINE=1 <vllm venv>/bin/python -I vllm_prompt_only.py prompts.json <out.json> [model ...]
"""
import json
import sys

from PIL import Image
from vllm.engine.arg_utils import EngineArgs
from vllm.renderers import renderer_from_config


def main(prompts: str, out: str, models: list[str]) -> None:
    P = json.load(open(prompts))["models"]
    res = {}
    for model in models or list(P):
        repo = P[model]["repo"]
        cfg = EngineArgs(model=repo, max_model_len=2048, limit_mm_per_prompt={"image": 1}).create_engine_config()
        renderer = renderer_from_config(cfg)
        res[model] = {}
        for dataset, D in P[model]["datasets"].items():
            mm = {k: v for k, v in D["processor_kwargs"].items() if k != "padding_side"}
            rows = []
            for x in D["items"]:
                req = {"prompt": x["text"], "multi_modal_data": {"image": Image.open(x["image"]).convert("RGB")},
                       **({"mm_processor_kwargs": mm} if mm else {})}
                ids = list(renderer.render_cmpl([req])[0]["prompt_token_ids"])
                rows.append({"item_id": x["item_id"], "hf_len": len(x["gen_ids"]), "vllm_len": len(ids),
                             "equal": ids == x["gen_ids"]})
            res[model][dataset] = rows
            print(model, dataset, "renderer prompt ids equal", sum(r["equal"] for r in rows), "/", len(rows),
                  "lens", sorted({r["vllm_len"] for r in rows}), flush=True)
    json.dump(res, open(out, "w"))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3:])

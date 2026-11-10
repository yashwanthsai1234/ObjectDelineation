"""Test 5 + offline checks (no model weights needed).

(a) TRL #6294 diagnostic: processor-expanded prompt -> number of separate image-token runs per model and setting.
    vLLM re-expands placeholders; a model whose expanded prompt has >1 run (Idefics3/SmolVLM tiling) is the one
    that breaks. Without vLLM TRL feeds the expanded ids + pixel_values straight to model.generate (no re-expansion).
(b) GDPO (normalize_then_sum) arithmetic on a constant and a near-constant term, using the exact TRL 1.16 formula.

    HF_HOME=./hf HF_HUB_OFFLINE=1 .venv/bin/python -I t5_checks.py
"""
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np  # noqa: E402
import torch  # noqa: E402
from PIL import Image  # noqa: E402

from common import HERE, model_path  # noqa: E402


def runs(ids, tid):
    n, prev = 0, False
    for t in ids:
        cur = t == tid
        n += cur and not prev
        prev = cur
    return n


def diag_6294():
    from transformers import AutoProcessor
    img = Image.new("RGB", (448, 448), (120, 60, 30))
    out = {}
    cases = [("smolvlm-256m", {}), ("smolvlm-256m", {"size": {"longest_edge": 512}}),
             ("smolvlm-256m", {"do_image_splitting": False}), ("smolvlm-500m", {}),
             ("lfm2.5-vl-450m", {}), ("qwen3-vl-2b", {}), ("qwen3.5-0.8b", {})]
    for name, kw in cases:
        proc = AutoProcessor.from_pretrained(model_path(name), **kw)
        tok = proc.tokenizer
        msgs = [{"role": "user", "content": [{"type": "image"}, {"type": "text", "text": "Outline the disc."}]}]
        text = proc.apply_chat_template(msgs, add_generation_prompt=True)
        raw = tok(text, add_special_tokens=False)["input_ids"]
        exp = proc(text=[text], images=[[img]], return_tensors="pt")["input_ids"][0].tolist()
        tid = getattr(proc, "image_token_id", None)
        if tid is None:
            tid = tok.convert_tokens_to_ids(proc.image_token)
        out[f"{name} {json.dumps(kw)}"] = {
            "tokenizer_only_len": len(raw), "image_tokens_tokenizer_only": raw.count(tid),
            "expanded_len": len(exp), "image_tokens_expanded": exp.count(tid), "image_token_runs_expanded": runs(exp, tid)}
        print(name, kw, out[f"{name} {json.dumps(kw)}"], flush=True)
    return out


def gdpo_checks():
    """TRL normalize_then_sum: z_k = (r_k - mean_g) / (std_g + 1e-4); sum w_k z_k; then batch z-score."""
    from trl.trainer.utils import nanstd

    def trl_ns(r, w, G):
        r = torch.tensor(r, dtype=torch.float32)
        g = r.view(-1, G, r.shape[1])
        z = (g - torch.nanmean(g, 1, keepdim=True)) / (nanstd(g, dim=1, keepdim=True) + 1e-4)
        s = (z.view(-1, r.shape[1]) * torch.tensor(w)).nansum(1)
        return ((s - torch.nanmean(s)) / (nanstd(s) + 1e-4)).tolist()

    G = 4
    iou = [0.80, 0.82, 0.85, 0.90]
    res = {}
    for label, legal in [("constant 1.0", [1.0, 1.0, 1.0, 1.0]), ("near-constant (one 0.999)", [1.0, 1.0, 1.0, 0.999]),
                         ("near-constant (one 0.99999)", [1.0, 1.0, 1.0, 0.99999])]:
        r = [[a, b] for a, b in zip(iou, legal)]
        adv = trl_ns(r, [1.0, 1.0], G)
        adv_iou_only = trl_ns([[a] for a in iou], [1.0], G)
        res[label] = {"advantages": [round(x, 3) for x in adv], "iou_only": [round(x, 3) for x in adv_iou_only],
                      "any_nan": any(np.isnan(adv))}
        print("GDPO", label, res[label], flush=True)
    return res


if __name__ == "__main__":
    out = {"diag_6294": diag_6294(), "gdpo": gdpo_checks()}
    (HERE / "results").mkdir(exist_ok=True)
    (HERE / "results" / "t5_checks.json").write_text(json.dumps(out, indent=1))

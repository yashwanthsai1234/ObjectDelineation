"""Test 1: tokens per polygon answer by tokenizer and format; image tokens per 448 px crop.

    HF_HOME=./hf .venv/bin/python -I t1_tokens.py
"""
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np  # noqa: E402
from PIL import Image  # noqa: E402

from common import HERE, MODELS, model_path  # noqa: E402


def polygon(n, seed=0):
    """Realistic disc outline: noisy ellipse, centre ~ (480, 520), radii ~ (170, 190) in a 0-1000 frame."""
    rng = np.random.default_rng(seed)
    t = np.linspace(0, 2 * np.pi, n, endpoint=False)
    r = 1 + 0.03 * rng.standard_normal(n)
    x = 480 + 170 * r * np.cos(t)
    y = 520 + 190 * r * np.sin(t)
    return np.stack([x, y], 1)


def formats(n):
    p = polygon(n)
    hole = polygon(max(8, n // 2), seed=1)
    hole = (hole - [480, 520]) * 0.45 + [470, 500]
    pi, hi = np.rint(p).astype(int).tolist(), np.rint(hole).astype(int).tolist()
    pf = (p / 1000).round(3).tolist()
    p100 = np.rint(p / 10).astype(int).tolist()
    flat = lambda pts: " ".join(f"{x} {y}" for x, y in pts)
    deltas = [pi[0]] + [[b[0] - a[0], b[1] - a[1]] for a, b in zip(pi[:-1], pi[1:])]
    return {
        "a_json_min": json.dumps({"disc": [{"outer": pi, "holes": []}]}, separators=(",", ":")),
        "a_json_spaced": json.dumps({"disc": [{"outer": pi, "holes": []}]}),
        "a_json_min_1hole": json.dumps({"disc": [{"outer": pi, "holes": [hi]}]}, separators=(",", ":")),
        "b_compact": f"disc: {flat(pi)}",
        "b_compact_1hole": f"disc: {flat(pi)} ; hole: {flat(hi)}",
        "c_qwen_native": json.dumps([{"label": "optic disc", "polygon": pi}]),  # Qwen3-VL style, 0-1000 ints
        "c_lfm_native_float3": json.dumps([{"label": "optic disc", "polygon": pf}]),  # LFM style, 0-1 floats
        "d_compact_grid100": f"disc: {flat(p100)}",  # coarser 0-100 grid
        "d_delta_ints": f"disc: {flat(deltas)}",
    }


def main():
    from transformers import AutoProcessor, AutoTokenizer
    names = list(MODELS)
    rows = []
    img = Image.new("RGB", (448, 448), (120, 60, 30))
    image_tokens = {}
    for name in names:
        path = model_path(name)
        try:
            tok = AutoTokenizer.from_pretrained(path)
        except Exception as e:  # noqa: BLE001
            print(f"{name}: tokenizer failed: {e!r}")
            continue
        print(f"{name}: {type(tok).__name__} vocab={len(tok)}")
        for n in (16, 32, 64):
            for fmt, s in formats(n).items():
                rows.append({"model": name, "n": n, "format": fmt, "chars": len(s),
                             "tokens": len(tok(s, add_special_tokens=False)["input_ids"])})
        # digit splitting probe
        probe = {k: tok.tokenize(k) for k in ("1000", "523", "0.523", " 523")}
        print("   tokenize probe:", probe)
        # image tokens for a 448x448 crop, via the processor + chat template
        try:
            proc = AutoProcessor.from_pretrained(path)
            msgs = [{"role": "user", "content": [{"type": "image"}, {"type": "text", "text": "hi"}]}]
            text = proc.apply_chat_template(msgs, add_generation_prompt=True)
            enc = proc(text=[text], images=[[img]], return_tensors="pt")
            ids = enc["input_ids"][0].tolist()
            tid = getattr(proc, "image_token_id", None) or tok.convert_tokens_to_ids(getattr(proc, "image_token", "<image>"))
            image_tokens[name] = {"image_token_count": ids.count(tid), "prompt_len": len(ids)}
        except Exception as e:  # noqa: BLE001
            image_tokens[name] = {"error": repr(e)[:300]}
        print("   image tokens (448x448):", image_tokens[name])
    out = HERE / "results"
    out.mkdir(exist_ok=True)
    (out / "t1_tokens.json").write_text(json.dumps({"rows": rows, "image_tokens": image_tokens}, indent=1))
    # pivot table: format x model for each n
    fmts = list(formats(16))
    for n in (16, 32, 64):
        print(f"\n### n={n} vertices (tokens)")
        models = [m for m in names if any(r["model"] == m for r in rows)]
        print("| format | " + " | ".join(models) + " |")
        print("|---|" + "---|" * len(models))
        for f in fmts:
            vals = [next(r["tokens"] for r in rows if r["model"] == m and r["n"] == n and r["format"] == f) for m in models]
            print(f"| {f} | " + " | ".join(map(str, vals)) + " |")


if __name__ == "__main__":
    main()

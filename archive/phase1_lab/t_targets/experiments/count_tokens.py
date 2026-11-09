"""Add token counts (Qwen3-VL-2B, LFM2.5-VL-450M, SmolVLM-256M) to a results jsonl.
Run with the repo venv (transformers):  python experiments/count_tokens.py in.jsonl out.jsonl"""
import json
import sys

from transformers import AutoTokenizer

TOKS = {"qwen3vl": "Qwen/Qwen3-VL-2B-Instruct", "lfm25": "LiquidAI/LFM2.5-VL-450M",
        "smolvlm": "HuggingFaceTB/SmolVLM-256M-Instruct"}

recs = [json.loads(l) for l in open(sys.argv[1])]
texts = []
for r in recs:
    texts.append(r["text"])
    texts += [o["text"] for o in r["objs"].values()]
uniq = sorted(set(texts))
counts = {}
for name, path in TOKS.items():
    tok = AutoTokenizer.from_pretrained(path)
    ids = tok(uniq, add_special_tokens=False)["input_ids"]
    counts[name] = {t: len(i) for t, i in zip(uniq, ids)}
with open(sys.argv[2], "w") as f:
    for r in recs:
        r["tokens"] = {n: counts[n][r["text"]] for n in TOKS}
        for o in r["objs"].values():
            o["tokens"] = {n: counts[n][o["text"]] for n in TOKS}
        f.write(json.dumps(r) + "\n")

"""Token length of every certified target (Qwen3-VL-2B tokenizer, local cache only). Run from REPO with -I."""
import json, numpy as np, os
os.environ["HF_HUB_OFFLINE"]="1"
from transformers import AutoTokenizer
tok = AutoTokenizer.from_pretrained("Qwen/Qwen3-VL-2B-Instruct")
srcs = {"prepared/train (centred)": ("data/prepared/lidc/train.jsonl", "rec"), "train-shifted": ("data/views/lidc-train-shifted/lidc/train.jsonl", "rec"),
        "train-rescaled": ("data/views/lidc-train-rescaled/lidc/train.jsonl", "rec"), "val-shifted": ("data/eval_views/lidc-val-shifted/lidc/val.jsonl", "rec"),
        "val-rescaled": ("data/eval_views/lidc-val-rescaled/lidc/val.jsonl", "rec"), "sft bundle v3": ("data/views/lidc-bundle-v3/rows.jsonl", "row")}
for n, (p, kind) in srcs.items():
    texts = []
    for l in open(p):
        r = json.loads(l)
        texts += [t for t in r["targets"] if t] if kind == "rec" else [r["target"]]
    L = np.array([len(x) for x in tok(texts, add_special_tokens=False)["input_ids"]])
    print(f"{n:26s} n={len(L):6d} mean={L.mean():6.1f} p50={np.median(L):5.0f} p90={np.percentile(L,90):5.0f} p99={np.percentile(L,99):5.0f} max={L.max():5d} >600={int((L>600).sum())} >1200={int((L>1200).sum())} >1536={int((L>1536).sum())}")

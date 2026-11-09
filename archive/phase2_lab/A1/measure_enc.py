"""Tokens per vertex and round-trip legality of candidate polygon encodings on the real SFT targets.
Run: .venv/bin/python -I measure_enc.py <repo>/data/targets  (data is untrusted: -I)."""

import glob
import json
import os
import sys

import numpy as np
from huggingface_hub import hf_hub_download
from tokenizers import Tokenizer

from objdel.geometry.text import SCHEMAS, parse

TOKS = {
    "qwen3-vl-2b": "Qwen/Qwen3-VL-2B-Instruct",
    "qwen3.5-0.8b": "Qwen/Qwen3.5-0.8B",
    "smolvlm": "HuggingFaceTB/SmolVLM-256M-Instruct",
    "lfm2.5-vl": "LiquidAI/LFM2.5-VL-450M",
    "gemma3n": "unsloth/gemma-3n-E2B-it",
    "florence2": "microsoft/Florence-2-base",
    "phi4-mm": "microsoft/Phi-4-multimodal-instruct",
    "minicpm-v4": "openbmb/MiniCPM-V-4",
    "internvl3.5": "OpenGVLab/InternVL3_5-1B-HF",
}
root = sys.argv[1]
recs = []
for f in sorted(glob.glob(os.path.join(root, "*", "*.jsonl"))):
    ds, split = f.split(os.sep)[-2], os.path.basename(f)[:-6]
    for line in open(f):
        r = json.loads(line)
        r["ds"], r["split"] = ds, split
        recs.append(r)
print("records", len(recs))


def rings(objs, names):
    for n in names:
        for o, hs in objs[n]:
            yield n, o, hs


def enc(r, mode):
    names = SCHEMAS[r["ds"]]
    p = parse(r["text"], names)
    assert p.ok
    k = 1000 // r["W"]

    def pt(v):
        return v // k if mode.startswith("px") or "_px" in mode else v

    def ring(a):
        a = np.asarray(a)
        if "delta" in mode:
            a = np.vstack([a[:1], np.diff(a, axis=0)])
        a = [[int(pt(x)) if "delta" not in mode else int(x // k if "_px" in mode else x) for x in v] for v in a]
        if "flat" in mode:
            return "[" + ",".join(f"{x},{y}" for x, y in a) + "]"
        if "loc" in mode:
            return "[" + "".join(f"<x{x}><y{y}>" for x, y in a) + "]"
        if "space" in mode:
            return "[" + " ".join(f"{x} {y}" for x, y in a) + "]"
        return "[" + ",".join(f"[{x},{y}]" for x, y in a) + "]"

    parts = []
    for n in names:
        ps = ",".join(
            '{"outer":' + ring(o) + ',"holes":[' + ",".join(ring(h) for h in hs) + "]}" for o, hs in p.objs[n]
        )
        parts.append(f'"{n}":[{ps}]')
    nv = sum(len(o) + sum(len(h) for h in hs) for _, o, hs in rings(p.objs, names))
    return "{" + ",".join(parts) + "}", nv


MODES = ["json1k", "json_px", "flat1k", "flat_px", "space1k", "delta1k", "delta_px", "loc1k"]
texts = {m: [] for m in MODES}
nverts = []
for r in recs:
    for m in MODES:
        t, nv = enc(r, m)
        texts[m].append(t)
    nverts.append(nv)
assert texts["json1k"] == [r["text"] for r in recs], "serializer mismatch"
nverts = np.array(nverts)
ds = np.array([r["ds"] for r in recs])
# coordinate stats
vals = []
for r in recs:
    p = parse(r["text"], SCHEMAS[r["ds"]])
    for _, o, hs in rings(p.objs, SCHEMAS[r["ds"]]):
        for a in [o, *hs]:
            vals.append((r["ds"], r["W"], a.ravel()))
for d in ("riga", "lidc"):
    v = np.concatenate([a for dd, W, a in vals if dd == d])
    Ws = [W for dd, W, a in vals if dd == d]
    print(
        d,
        "coords",
        v.size,
        "min",
        v.min(),
        "max",
        v.max(),
        "n==0",
        (v == 0).sum(),
        "n==1000",
        (v == 1000).sum(),
        "frame sides",
        dict(zip(*np.unique(Ws, return_counts=True))),
    )
out = {}
for name, repo in TOKS.items():
    tk = Tokenizer.from_file(hf_hub_download(repo, "tokenizer.json"))
    if "loc" in " ".join(MODES):
        tk.add_special_tokens([f"<x{i}>" for i in range(1001)] + [f"<y{i}>" for i in range(1001)])
    row = {}
    for m in MODES:
        encs = tk.encode_batch(texts[m], add_special_tokens=False)
        n = np.array([len(e.ids) for e in encs])
        rt = np.mean(
            [
                tk.decode(e.ids, skip_special_tokens=False).replace(" ", "") == t.replace(" ", "")
                if m == "loc1k"
                else tk.decode(e.ids, skip_special_tokens=False) == t
                for e, t in zip(encs, texts[m])
            ]
        )
        for d in ("riga", "lidc"):
            s = ds == d
            row[f"{m}/{d}"] = dict(
                med=float(np.median(n[s])), p95=float(np.percentile(n[s], 95)), tpv=float(n[s].sum() / nverts[s].sum())
            )
        row[f"{m}/rt"] = float(rt)
    out[name] = row
    print(
        name,
        " ".join(
            f"{m}:{row[m + '/riga']['tpv']:.2f}/{row[m + '/lidc']['tpv']:.2f} rt={row[m + '/rt']:.4f}" for m in MODES
        ),
    )
json.dump(out, open(sys.argv[2] if len(sys.argv) > 2 else "enc_tokens.json", "w"), indent=1)

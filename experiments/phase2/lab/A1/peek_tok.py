import re

from transformers import AutoTokenizer

REPOS = {
    "smolvlm-256m": "HuggingFaceTB/SmolVLM-256M-Instruct",
    "lfm2.5-vl-450m": "LiquidAI/LFM2.5-VL-450M",
    "lfm2.5-vl-1.6b": "LiquidAI/LFM2.5-VL-1.6B",
    "qwen3-vl-2b": "Qwen/Qwen3-VL-2B-Instruct",
    "qwen3.5-0.8b": "Qwen/Qwen3.5-0.8B",
    "internvl3.5-1b": "OpenGVLab/InternVL3_5-1B-HF",
    "gemma3n-e2b": "unsloth/gemma-3n-E2B-it",
    "medgemma-4b": "unsloth/medgemma-4b-it",
    "florence2": "microsoft/Florence-2-base",
    "moondream2": "vikhyatk/moondream2",
    "phi4-mm": "microsoft/Phi-4-multimodal-instruct",
    "minicpm-v4": "openbmb/MiniCPM-V-4",
    "fastvlm-0.5b": "apple/FastVLM-0.5B",
}
s = '{"nodule":[{"outer":[[500,440],[560,440],[1000,7],[42,9]],"holes":[]}]}'
for k, r in REPOS.items():
    try:
        t = AutoTokenizer.from_pretrained(r)
    except Exception as e:
        print(k, "LOAD FAIL", type(e).__name__, str(e)[:150])
        continue
    ids = t(s, add_special_tokens=False)["input_ids"]
    toks = [t.decode([i]) for i in ids]
    v = t.get_vocab()
    num = [w for w in v if re.fullmatch(r"(Ġ|▁| )?[0-9]+", w)]
    maxlen = max((len(re.sub(r"\D", "", w)) for w in num), default=0)
    print(
        f"{k:16s} n={len(ids):3d} vocab={len(v)} numtoks={len(num)} maxdigits={maxlen} rt={t.decode(ids) == s} |"
        + "|".join(toks)
    )

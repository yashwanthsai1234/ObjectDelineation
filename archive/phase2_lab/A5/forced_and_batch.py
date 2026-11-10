"""A5: forced-token share (jump-forward potential) and batched mask cost (trainer-side grammar masking estimate).

Usage (scratch venv with xgrammar 0.2.8, transformers, torch; HF cache with the three tokenizers):
    HF_HUB_OFFLINE=1 venv/bin/python -I forced_and_batch.py targets_texts.json
targets_texts.json: {"patterns": {schema: regex}, "texts": {schema: [target text, ...]}} exported from data/targets
with objdel.geometry.text.pattern (see docs/research/notes/Phase 2 SFT/A5_decoding_serving.md §11).
"""
import json
import random
import sys
import time

import numpy as np
import torch
import xgrammar as xgr
from transformers import AutoTokenizer

MODELS = {
    "smolvlm-256m": ("HuggingFaceTB/SmolVLM-256M-Instruct", [49279]),
    "lfm2.5-vl-450m": ("LiquidAI/LFM2.5-VL-450M", [7]),
    "qwen3-vl-2b": ("Qwen/Qwen3-VL-2B-Instruct", [151645, 151643]),
}


def n_allowed(bitmask, vocab: int) -> int:
    return int(np.unpackbits(bitmask[0].numpy().view(np.uint8), bitorder="little")[:vocab].sum())


def main():
    data = json.load(open(sys.argv[1]))
    rng = random.Random(0)
    for name, (repo, stop) in MODELS.items():
        tok = AutoTokenizer.from_pretrained(repo)
        vocab = max(len(tok), max(tok.get_vocab().values()) + 1)
        info = xgr.TokenizerInfo.from_huggingface(tok, vocab_size=vocab, stop_token_ids=stop)
        compiler = xgr.GrammarCompiler(info, max_threads=1)
        for ds in ("riga", "lidc"):
            grammar = compiler.compile_regex(data["patterns"][ds])
            texts = rng.sample(data["texts"][ds], 256)
            seqs = [tok.encode(t, add_special_tokens=False) + [stop[0]] for t in texts]
            bm = xgr.allocate_token_bitmask(1, vocab)
            forced = total = 0
            jump = []
            for ids in seqs[:100]:
                m = xgr.GrammarMatcher(grammar, override_stop_tokens=stop)
                for tid in ids:
                    m.fill_next_token_bitmask(bm)
                    forced += n_allowed(bm, vocab) == 1
                    total += 1
                    jump.append(len(m.find_jump_forward_string()))
                    m.accept_token(tid)
            # trainer-side masking for one 256-rollout batch: a bitmask per completion position
            length = max(map(len, seqs))
            t0 = time.perf_counter()
            masks = torch.zeros((len(seqs), length, (vocab + 31) // 32), dtype=torch.int32)
            for b, ids in enumerate(seqs):
                m = xgr.GrammarMatcher(grammar, override_stop_tokens=stop)
                for k, tid in enumerate(ids):
                    m.fill_next_token_bitmask(masks[b], k)
                    m.accept_token(tid)
            dt = time.perf_counter() - t0
            ntok = sum(map(len, seqs))
            print(json.dumps({"model": name, "ds": ds, "forced_share": round(forced / total, 4),
                              "jump_forward_chars_mean": round(sum(jump) / len(jump), 3), "batch256_tokens": ntok,
                              "batch256_mask_s_1thread": round(dt, 3), "us_per_token": round(1e6 * dt / ntok, 2),
                              "mask_MB": round(masks.numel() * 4 / 2**20, 1)}), flush=True)
            del masks


if __name__ == "__main__":
    main()

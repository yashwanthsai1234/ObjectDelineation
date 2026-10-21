"""A6 (a): SFT data path on CPU: objdel.train.dataset.sft_rows (row building), then per batch the dataset fetch
(lazy PNG decode by datasets.Image) and TRL's DataCollatorForVisionLanguageModeling with each screened model's
processor (cached in ~/.cache/huggingface; HF_HUB_OFFLINE=1). cProfile splits the collator into image processing
and chat template / tokenization.

    HF_HUB_OFFLINE=1 .venv/bin/python -I -m experiments.phase2.lab.A6.prof_collate --out <dir>
"""
import argparse
import time

import numpy as np
from transformers import AutoProcessor
from trl.trainer.sft_trainer import DataCollatorForVisionLanguageModeling

from objdel.train import configs
from objdel.train.dataset import sft_rows

from .common import PREPARED, by_module, dump, env, profile, stats_ms


def collate_bench(ds_rows, collator, bs, n_batches, rng):
    fetch, coll, tokens, pix = [], [], [], []
    for _ in range(n_batches):
        idx = rng.choice(len(ds_rows), bs, replace=False).tolist()
        t0 = time.perf_counter()
        ex = [dict(ds_rows[i]) for i in idx]  # datasets decodes the PNG to PIL here
        t1 = time.perf_counter()
        out = collator(ex)
        t2 = time.perf_counter()
        fetch.append(1e3 * (t1 - t0))
        coll.append(1e3 * (t2 - t1))
        tokens.append(int(out["attention_mask"].sum()))
        pix.append(tuple(out["pixel_values"].shape))
    return {"fetch_ms": stats_ms(fetch), "collate_ms": stats_ms(coll), "tokens_per_batch": float(np.mean(tokens)),
            "pixel_values_shape": str(pix[0]), "keys": None}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--batches", type=int, default=6)
    ap.add_argument("--models", default="smolvlm-256m,lfm2.5-vl-450m,qwen3-vl-2b")
    a = ap.parse_args()
    res = {"env": env(), "rows": {}, "collate": {}}
    rng = np.random.default_rng(0)
    rows = {}
    tok_model = "smolvlm-256m"  # rows are built once, with SmolVLM's tokenizer and budget (I5 filter)
    tok = AutoProcessor.from_pretrained(configs.MODELS[tok_model]).tokenizer
    for ds in ("riga", "lidc"):
        t0 = time.perf_counter()
        rows[ds] = sft_rows(PREPARED, ds, "train", tok, configs.MAX_COMPLETION[tok_model])
        dt = time.perf_counter() - t0
        res["rows"][f"{ds}/train"] = {"rows": len(rows[ds]), "s": dt, "ms_per_row": 1e3 * dt / len(rows[ds])}
        print("rows", ds, len(rows[ds]), round(dt, 1), "s", flush=True)
    _, prow, _ = profile(sft_rows, PREPARED, "riga", "train", tok, configs.MAX_COMPLETION[tok_model], top=400)
    res["rows"]["riga/train"]["cprofile_by_file_ms"] = {k: round(v, 1) for k, v in list(by_module(prow).items())[:10]}
    res["rows"]["riga/train"]["cprofile_top"] = [(round(tt, 1), round(ct, 1), nc, nm) for tt, ct, nc, nm in prow[:10]]
    for model in a.models.split(","):
        for ds in ("riga", "lidc"):
            proc = AutoProcessor.from_pretrained(configs.MODELS[model], **configs.processor_kwargs(model, ds))
            d = rows[ds]
            ip = getattr(proc, "image_processor", None)
            for aol in (True, False):
                try:
                    collator = DataCollatorForVisionLanguageModeling(proc, assistant_only_loss=aol)
                    collator([dict(d[0]), dict(d[1])])
                    break
                except Exception as e:  # noqa: BLE001 - recorded
                    res["collate"][f"{model}/{ds}/assistant_only_loss={aol}"] = f"error: {type(e).__name__}: {e}"[:300]
            for bs in (8, 16):
                b = collate_bench(d, collator, bs, a.batches, rng)
                b["image_processor"] = type(ip).__name__ if ip is not None else None
                b["assistant_only_loss"] = aol
                res["collate"][f"{model}/{ds}/bs{bs}"] = b
                print(model, ds, bs, "fetch", round(b["fetch_ms"]["median"], 1), "collate",
                      round(b["collate_ms"]["median"], 1), "tokens", b["tokens_per_batch"], b["pixel_values_shape"],
                      b["image_processor"], flush=True)
            idx = rng.choice(len(d), 16, replace=False).tolist()
            ex = [dict(d[i]) for i in idx]
            _, crow, ctot = profile(collator, ex, top=400)
            res["collate"][f"{model}/{ds}/bs16"]["cprofile_by_file_ms"] = {
                k: round(v, 2) for k, v in list(by_module(crow).items())[:12]}
            res["collate"][f"{model}/{ds}/bs16"]["cprofile_top"] = [
                (round(tt, 2), round(ct, 2), nc, nm) for tt, ct, nc, nm in crow[:12]]
    res["env_end"] = env()
    dump(res, f"{a.out}/collate.json")


if __name__ == "__main__":
    main()

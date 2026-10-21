"""G-probe (plan §5.1, §3.2 ④): can GRPO improve this checkpoint? 20 GRPO steps at the pre-registered probe LR.

    python -I -m experiments.phase2.g_probe --model lfm2.5-vl-450m (--ckpt runs/r1/branches/e1 | --select runs/r1/select)
        --train data/prepared --probe data/eval_views/lidc-val-shifted --out runs/r1/g_probe

configs.grpo_config with prereg.G_PROBE: learning rate 5e-6 (full FT) or 5e-5 (LoRA), constant (grpo_config's
default), 20 steps, G = 8, 32 LIDC train prompts per step with disagreement items x4 (dataset.grpo_rows oversample),
HF generation (use_vllm=False), sampling fields from configs.sampling(model, "grpo") (D13). Reward: TRLReward on the
preloaded caches (serve.loader) and a fork pool made before the model loads. Fixed G-probe set: a seeded subset of
`items` (128) shifted-val records (select_ckpt.view_rows), K = 8 at T = 1, sampled through the trainer's own
generation path (objdel.train.sample.sample): twice at step 0 with different seeds (the null, half-width h0) and once
at step 20. Per item: the raw composite R (objdel.eval.signal). Verdict and CIs: prereg.g_probe, with the step-0
legal rate = the mean of both step-0 draws. Also reported: per-item Spearman of step-0 sd_geom, r_std and mixed with
the gain (D12 veto check, §6.2), the training-log form (mean logged raw composite R = signal.WEIGHTS on TRL's
per-term means, steps 1-5 vs 11-20) and the step timing.

--select (select_ckpt's output dir): probe the selected branch; when it is not PASS, the D12(c) repair probes the
next-earlier eligible branch once. <out>/g_probe.json is the last probe's result plus "probed" (branch names in
order) and "earlier" (the results before it); its "ckpt" is the checkpoint that received the final verdict (C2 depth).
"""

import argparse
import gc
import json
import sys
import time
from pathlib import Path

import numpy as np

from experiments.phase2.prereg import G_PROBE, g_probe
from experiments.phase2.select_ckpt import _score, view_rows  # imported before the pool forks (no import in a child)
from objdel.eval import signal

SEEDS = (11, 12, 13)  # step-0 draw a, step-0 draw b, step 20


def score(trainer, rows, probe_root, pool, seed: int) -> dict:
    """K = 8 samples per probe row -> per-item panel metrics and the group records."""
    import torch

    from objdel.train import sample

    torch.manual_seed(seed)
    done = sample.sample(trainer, rows, 8)
    meta = {r["key"]: r for r in rows}
    jobs = [(d["item_id"], d["samples"], d["truncated"], *(meta[d["item_id"]][k] for k in ("W", "H", "names")),
             probe_root) for d in done]  # fmt: skip
    res = list(pool.map(_score, jobs))
    groups = [signal.group(r["item_id"], d["samples"], r) for r, d in zip(res, done, strict=True)]
    m = signal.per_item(groups)
    ids = [g["item_id"] for g in groups]
    return {"ids": ids, "per_item": m, "legal": float(np.mean(m["legal"])), "sigma_share": float(np.mean(m["sigma"]))}


def verdict(s0a: dict, s0b: dict, s20: dict) -> dict:
    """prereg.g_probe: gain = R20 - mean of the two step-0 draws; step-0 legal = mean of both draws."""
    r0a, r0b, r20 = (s["per_item"]["r_mean"] for s in (s0a, s0b, s20))
    return g_probe(r20, r0a, r0b, s20["ids"], (s0a["legal"] + s0b["legal"]) / 2, s20["legal"], s20["sigma_share"])


def train_log(log: list[dict]) -> dict:
    """§5.1 training-log form: mean logged raw composite R, steps 1-5 vs 11-20 (TRL logs rewards/<term>/mean)."""
    R = [sum(w * h.get(f"rewards/{t}/mean", np.nan) for t, w in signal.WEIGHTS.items()) for h in log
         if "rewards/iou_tc/mean" in h]  # fmt: skip
    return {"steps_1_5": float(np.mean(R[:5])), "steps_11_20": float(np.mean(R[10:20]))}


def candidates(sel: dict) -> list[str]:
    """D12(c): the selected branch, then the next-earlier (shallower) eligible branch."""
    if not sel["selected"]:
        raise SystemExit("g_probe: select_ckpt found no eligible branch (a §5.1 gate failure): nothing to probe")
    b, depth = sel["branches"], sel["branches"][sel["selected"]]["depth"]
    earlier = sorted((x["depth"], n) for n, x in b.items() if x["eligible"] and x["depth"] < depth)
    return [sel["selected"]] + [n for _, n in earlier[-1:]]


def probe(a, ckpt: Path, train, rows, pool) -> dict:
    """One G-probe run on one checkpoint."""
    from transformers import AutoModelForImageTextToText, AutoProcessor
    from trl import GRPOTrainer

    from objdel.train import configs, rewards, serve

    model = AutoModelForImageTextToText.from_pretrained(ckpt)
    proc = AutoProcessor.from_pretrained(ckpt, **configs.processor_kwargs(a.model, "lidc"))
    peft = configs.lora(a.model)
    lr = G_PROBE["lr_lora"] if peft else G_PROBE["lr_full"]
    over = {"learning_rate": lr, "max_steps": G_PROBE["steps"], "use_vllm": False, "save_strategy": "no"}
    over |= {"logging_steps": 1, "log_completions": False} | configs.sampling(a.model, "grpo") | json.loads(a.grpo)
    cfg = configs.grpo_config(a.model, str(a.out / "trainer"), **over)
    funcs, _ = rewards.TRLReward(serve.loader(a.train), rewards.stop_ids(model, proc), pool=pool).funcs()
    trainer = GRPOTrainer(
        model=model, reward_funcs=funcs, args=cfg, train_dataset=train, processing_class=proc, peft_config=peft
    )
    root = str(a.probe.resolve())
    t0 = time.time()
    s0a, s0b = (score(trainer, rows, root, pool, s) for s in SEEDS[:2])
    t1 = time.time()
    trainer.train()
    t2 = time.time()
    s20 = score(trainer, rows, root, pool, SEEDS[2])
    r0a, r0b, r20 = (s["per_item"]["r_mean"] for s in (s0a, s0b, s20))
    gain = r20 - (r0a + r0b) / 2
    step0 = {k: (s0a["per_item"][k] + s0b["per_item"][k]) / 2 for k in ("sd_geom", "r_std", "mixed")}
    log = trainer.state.log_history
    return verdict(s0a, s0b, s20) | {
        "lr": lr,
        "steps": G_PROBE["steps"],
        "items": len(s20["ids"]),
        "ckpt": str(ckpt),
        "seeds": SEEDS,
        "legal": {"step0": [s0a["legal"], s0b["legal"]], "step20": s20["legal"]},
        "spearman_gain": {k: signal.spearman_ci(v, gain, s20["ids"]) for k, v in step0.items()},
        "per_item": {"ids": s20["ids"], "r0a": r0a.tolist(), "r0b": r0b.tolist(), "r20": r20.tolist()},
        "timing_s": {"step0_sampling": t1 - t0, "train": t2 - t1, "per_step": (t2 - t1) / G_PROBE["steps"]},
        "train_log_R": train_log(log),
        "log": log,
    }


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--model", required=True, help="configs.MODELS key")
    ck = ap.add_mutually_exclusive_group(required=True)
    ck.add_argument("--ckpt", type=Path, help="one model-only checkpoint (no repair)")
    ck.add_argument("--select", type=Path, help="select_ckpt output dir: selected branch, then the D12(c) repair")
    ap.add_argument("--train", type=Path, required=True, help="prepared root with lidc/train.jsonl and caches")
    ap.add_argument("--probe", type=Path, required=True, help="shifted-val eval view")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--items", type=int, help="probe items (default prereg.G_PROBE['items'])")
    ap.add_argument("--procs", type=int, default=4)
    ap.add_argument("--grpo", default="{}", help="JSON overrides of configs.grpo_config (smoke tests)")
    a = ap.parse_args(argv)
    if a.select:
        sel = json.loads((a.select / "selected.json").read_text())
        ckpts = [Path(sel["run"]) / "branches" / n for n in candidates(sel)]
    else:
        ckpts = [a.ckpt]
    from objdel.train import serve

    for root in (a.train, a.probe):
        serve.preload(root, ["lidc"])
    pool = serve.fork_pool(a.procs)  # before the model and CUDA (D14)
    from PIL import Image

    from objdel.train import dataset

    train = dataset.grpo_rows(a.train, "lidc", "train", oversample=dataset.DISAGREEMENT_OVERSAMPLE)
    rows = [
        r | {"image": Image.open(r["image"]).convert("RGB")} for r in view_rows(a.probe, a.items or G_PROBE["items"])
    ]
    docs = []
    for ckpt in ckpts:  # at most the selected branch and one repair branch
        docs.append(probe(a, ckpt, train, rows, pool))
        gc.collect()  # the trainer's reference cycles hold the model
        if docs[-1]["verdict"] == "PASS":
            break
    pool.shutdown()
    doc = docs[-1] | {"probed": [c.name for c in ckpts[: len(docs)]], "earlier": docs[:-1]}
    a.out.mkdir(parents=True, exist_ok=True)
    (a.out / "g_probe.json").write_text(json.dumps(doc, indent=1, default=float) + "\n")
    print(json.dumps({k: doc[k] for k in ("verdict", "h0", "gain", "probed")}))
    return 0


if __name__ == "__main__":
    sys.exit(main())

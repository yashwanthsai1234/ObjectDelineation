"""G-probe (PLAN §5.1, §3.2 ④): can GRPO improve this checkpoint? 20 GRPO steps at the pre-registered probe LR.

    python -I -m experiments.phase2.g_probe --model qwen3-vl-2b (--ckpt runs/r1/branches/e1 | --select runs/r1/select)
        --train data/prepared --probe data/eval_views/lidc-val-shifted --out runs/r1/g_probe
    python -I -m experiments.phase2.g_probe --spec specs/r1_qwen3vl2b.yaml --stage probe [--rd RD] [--data DATA]

The flags are an objdel.core.spec.ProbeStage (cli_stage), with the trainer flags train_rl has (TrainerFlags: --lr,
--grpo JSON, --shared-prefix, --compile, --plain, --trace, the C flags); objdel.rl.build.build_trainer builds the
trainer from them.

configs.grpo_config with prereg.G_PROBE: LR 5e-5 (LoRA; every model, ADR-0008), constant, 20 steps, G = 8, 32 LIDC train
prompts per step (disagreement items x4), vLLM colocate generation (owner 2026-11-06, perf F1: HF generate of 256
sequences needs 20-51 GiB of KV cache; needs the Phase 3 single-venv image, D15), sampling from
configs.sampling(model, "grpo") (D13);
fast.FastGRPO (per-chunk trimmed forwards, a per-step vision cache: perf C1, C6), TRLReward on the
preloaded caches; the probe metrics run in this process (D7). Probe set: a seeded subset of `items`
(128) shifted-val records (policy.evaluate.eval_rows), K = 8 through the trainer's own generation (sample.sample), twice
at step 0 (the null, h0) and once at step 20; per item the raw composite R. Verdict: prereg.g_probe, step-0 legal =
the mean of both draws. Also: Spearman of step-0 sd_geom / r_std / mixed with the gain (D12 veto, §6.2), the
training-log form (mean logged R, steps 1-5 vs 11-20) and the step timing.

--select (select_ckpt's dir): probe the selected branch, then (not PASS) the D12(c) repair: the next-earlier
eligible branch, once. <out>/g_probe.json = the last probe's result + "probed" (names) + "earlier" (results before);
its "ckpt" received the final verdict (C2 depth); NaN (e.g. Spearman of constant gains) is written as null. When
R1's §5.1 gates failed (selected.json "gates_pass" not true) nothing is probed: {"verdict": null, "skipped": "§5.1
gates failed: <gate names>", "probed": [], "gates_branch"}, exit 0, so the chain (train_sft && select_ckpt &&
g_probe) ends cleanly with the reason on disk.
"""
import gc
import sys
import time
from pathlib import Path

import numpy as np

from experiments.phase2.prereg import G_PROBE, g_probe
from objdel.eval import signal
from objdel.io import read_json
from objdel.policy.evaluate import eval_rows, score_items
from objdel.rl.build import build_trainer
from objdel.rl.run import save_json

SEEDS = (11, 12, 13)  # step-0 draw a, step-0 draw b, step 20


def score(trainer, rows, probe_root: Path, seed: int) -> dict:
    """K = 8 samples per probe row -> per-item panel metrics."""
    import torch

    from objdel.train import sample

    torch.manual_seed(seed)
    _, groups = score_items(rows, {d["item_id"]: d for d in sample.sample(trainer, rows, 8)}, probe_root)
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
    b, depth = sel["branches"], sel["branches"][sel["selected"]]["depth"]
    earlier = sorted((x["depth"], n) for n, x in b.items() if x["eligible"] and x["depth"] < depth)
    return [sel["selected"]] + [n for _, n in earlier[-1:]]


def probe(a, ckpt: Path, train, rows) -> dict:
    """One G-probe run on one checkpoint."""
    trainer, _, lr = build_trainer(a.model, ckpt, train, Path(a.out), [Path(a.train)], G_PROBE["steps"], a)
    cfg = trainer.args
    t0 = time.time()
    s0a, s0b = (score(trainer, rows, Path(a.probe).resolve(), s) for s in SEEDS[:2])
    t1 = time.time()
    trainer.train()
    t2 = time.time()
    s20 = score(trainer, rows, Path(a.probe).resolve(), SEEDS[2])
    r0a, r0b, r20 = (s["per_item"]["r_mean"] for s in (s0a, s0b, s20))
    step0 = {k: (s0a["per_item"][k] + s0b["per_item"][k]) / 2 for k in ("sd_geom", "r_std", "mixed")}
    log = trainer.state.log_history
    return verdict(s0a, s0b, s20) | {
        "lr": lr, "steps": G_PROBE["steps"], "items": len(s20["ids"]), "ckpt": str(ckpt), "seeds": SEEDS,
        "legal": {"step0": [s0a["legal"], s0b["legal"]], "step20": s20["legal"]},
        "spearman_gain": {k: signal.spearman_ci(v, r20 - (r0a + r0b) / 2, s20["ids"]) for k, v in step0.items()},
        "per_item": {"ids": s20["ids"], "r0a": r0a.tolist(), "r0b": r0b.tolist(), "r20": r20.tolist()},
        "timing_s": {"step0_sampling": t1 - t0, "train": t2 - t1, "per_step": (t2 - t1) / cfg.max_steps},
        "train_log_R": train_log(log), "log": log,
    }  # fmt: skip


def main(argv=None) -> int:
    """CLI: the ProbeStage -> one probe per candidate checkpoint (at most two) -> <out>/g_probe.json."""
    from objdel.core.spec import ProbeStage, cli_stage
    from objdel.train import serve

    a = cli_stage(ProbeStage, sys.argv[1:] if argv is None else argv)
    out = Path(a.out)
    ckpts = [Path(a.ckpt)] if a.ckpt else []
    if a.select:
        sel = read_json(Path(a.select) / "selected.json")
        if not (sel["gates_pass"] and sel["selected"]):  # R1 failed §5.1: a probe would not change the outcome
            why = ", ".join(sel.get("gates_failed") or ["no eligible branch (D12)"])
            return save(out, {"verdict": None, "skipped": f"§5.1 gates failed: {why}", "probed": [],
                            "gates_branch": sel["gates_branch"]})
        ckpts = [Path(sel["run"]) / "branches" / n for n in candidates(sel)]
    for root in (a.train, a.probe):
        serve.preload(root, ["lidc"])
    from PIL import Image

    from objdel.train import dataset

    train = dataset.grpo_rows(a.train, "lidc", "train", oversample=dataset.DISAGREEMENT_OVERSAMPLE)
    rows = eval_rows(a.probe, a.items or G_PROBE["items"])
    rows = [r | {"image": Image.open(r["image"]).convert("RGB")} for r in rows]
    docs = []
    for ckpt in ckpts:  # at most the selected branch and one repair branch
        docs.append(probe(a, ckpt, train, rows))
        gc.collect()  # the trainer's reference cycles hold the model
        if docs[-1]["verdict"] == "PASS":
            break
    return save(out, docs[-1] | {"probed": [c.name for c in ckpts[: len(docs)]], "earlier": docs[:-1]})


def save(out: Path, doc: dict) -> int:
    """<out>/g_probe.json and its one-line summary."""
    return save_json(out / "g_probe.json", doc, ("verdict", "skipped", "h0", "gain", "probed"))


if __name__ == "__main__":
    sys.exit(main())

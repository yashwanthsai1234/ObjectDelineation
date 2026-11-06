"""R0 stop rule (PLAN §4.7 step 9, §4.9): project R1 GPU-hours from the 30-step timing arm; stop above the PLAN's stop value.

    python -I -m experiments.phase2.project_r1 --model qwen3-vl-2b --hf_rate runs/<run>/hf-rate-qwen3vl/hf_rate.json \
        runs/<run>/r0-qwen3vl-arm/manifest.json [...]

Basis: of the given train_sft manifests of --model, the one timing arm with deterministic off and 6 loader workers
(R1's loader), trained as R1 is (freeze_vision, gradient checkpointing equal); other arms are listed. R1 steps: train_sft.plan on R1's config (its branches) at
R1's steps per epoch = the manifest's rows_per_epoch (the same bundle view, uncapped) // its effective batch, which
must equal R1's (s/step is per batch); the stable run plus every decay.

Evaluation: select_ckpt's item-samples in policy.chain's R1 call (eval_samples; with hf, --items SUBSET for
selection) x the sampler's s per sample, plus one model load per sampler child (sampler_calls). hf: the R0 hf_rate
step's s_per_sample and load_s. vllm: VLLM_S_PER_SAMPLE and VLLM_LOAD_S (measured in R0 v3), plus select_ckpt's HF
check (HF_CHECK_ITEMS x K at the hf rate, one load). Both read --hf_rate. Not counted: scoring, G-vision, G-probe.
Projected R1 h = steps x s/step / 3600 + eval h. Prints a verdict line and JSON; exit 0 = go, 1 = stop, 2 = usage or
data error (no usable basis, an arm trained unlike R1, no or another model's HF rate for an HF R1).
"""

import argparse
import json
import sys
from pathlib import Path

from experiments.phase2 import train_sft as TS

CONFIGS = Path(__file__).parent / "configs"
# vLLM 0.32 on the R0 v3 box (RTX 4090, 2026-11-03, runs/v3/r0.log): 256 base-weight answers of Qwen3-VL-2B
# (148,760 tokens, 581 per answer) generated in 23 s, engine start 72 s. Trained answers are shorter: conservative.
VLLM_S_PER_SAMPLE, VLLM_LOAD_S = 0.09, 95.0
HF_CHECK_ITEMS = 64  # select_ckpt.HF_CHECK_ITEMS (a test keeps them equal)
# R1 stop values in hours (PLAN §4.9): LFM 3.4 (2 x its 1.70 h high estimate); Qwen3-VL-2B 5.0 (expected 3.1 h;
# 5.0 h is about $2.1 of the $3.18 credit after R0)
R1 = {"lfm2.5-vl-450m": ("r1_lfm.yaml", 3.4), "qwen3-vl-2b": ("r1_qwen3vl2b.yaml", 6.0)}  # 6.0: owner 2026-11-05, 3 epochs
# select_ckpt's sampling (lidc/val records of the shifted, rescaled and rotated views; probe-precision 1,000 and
# probe-topology 680 records): K = 8 gates, K = 16 on shifted; policy.chain passes --items 256 with sampler hf
VAL, PROBES, K, K16, SUBSET = 1084, 1000 + 680, 8, 16, 256


def eval_samples(branches: int, items: int | None) -> int:
    """Item-samples of select_ckpt: shifted at K on every branch (the --items subset), then the selected branch's
    full shifted at K (only after a subset), rescaled and rotated at K, shifted at K16, both probes at K."""
    full = VAL * K if items else 0
    return branches * (items or VAL) * K + full + 2 * VAL * K + VAL * K16 + PROBES * K


def sampler_calls(branches: int, items: int | None) -> int:
    """Sampler children (one model load each): one per branch, then judge's views and probes."""
    return branches + (1 if items else 0) + 5


def _ckpt(c: dict) -> bool:
    from objdel.train import configs

    return c["sft"].get("gradient_checkpointing", c["model"] in configs._GRAD_CKPT)  # sft_config's default


def project(a) -> dict:
    import yaml

    docs = {str(m): json.loads(m.read_text()) for m in a.manifests}
    docs = {m: d for m, d in docs.items() if d["config"]["model"] == a.model}
    basis = [m for m, d in docs.items() if not d["config"].get("deterministic")
             and d["config"]["sft"].get("dataloader_num_workers", 6) == 6]  # fmt: skip
    cfg, stop = R1[a.model]
    r1 = TS.RunArgs(**yaml.safe_load((CONFIGS / cfg).read_text()))
    batch = r1.sft.get("per_device_train_batch_size", 8) * r1.sft.get("gradient_accumulation_steps", 1)
    if len(basis) != 1 or (t := docs[basis[0]].get("timing", {})).get("s_per_step") is None or t["batch"] != batch:
        raise ValueError(f"need one finished deterministic-off, workers-6 timing arm of {a.model} at batch {batch}")
    arm, want = docs[basis[0]]["config"], {"freeze_vision": r1.freeze_vision}
    got = {k: arm.get(k, False) for k in want} | {"gradient_checkpointing": _ckpt(arm)}
    if got != (want := want | {"gradient_checkpointing": _ckpt({"model": a.model, "sft": r1.sft})}):
        raise ValueError(f"{basis[0]}: timing arm trained with {got}, R1 ({cfg}) with {want}")
    p = TS.plan(r1, t["rows_per_epoch"] // batch)
    steps = p["stable"] + sum(b["decay"] for b in p["branches"].values())
    sampler = r1.sampler or "vllm"
    if a.hf_rate is None:
        raise ValueError(f"{cfg}: pass --hf_rate (R0's hf-rate-*/hf_rate.json; the HF sampler or the HF check)")
    rate = json.loads(a.hf_rate.read_text())
    if rate["model"] != a.model or not rate["s_per_sample"] > 0:
        raise ValueError(f"{a.hf_rate}: rate of {rate['model']} ({rate['s_per_sample']} s), want {a.model}")
    items = SUBSET if sampler == "hf" else None
    n, calls = eval_samples(len(r1.branches), items), sampler_calls(len(r1.branches), items)
    s, load = (rate["s_per_sample"], rate["load_s"]) if sampler == "hf" else (VLLM_S_PER_SAMPLE, VLLM_LOAD_S)
    check = 0.0 if sampler == "hf" else HF_CHECK_ITEMS * K * rate["s_per_sample"] + rate["load_s"]
    ev = {"sampler": sampler, "eval_samples": n, "s_per_sample": s, "sampler_calls": calls, "hf_check_s": check,
          "eval_h": (n * s + calls * load + check) / 3600}
    hours = steps * t["s_per_step"] / 3600 + ev["eval_h"]
    doc = {"model": a.model, "basis": basis[0], "s_per_step": t["s_per_step"], "r1_steps": steps, **ev,
           "projected_h": hours, "stop_h": stop, "go": hours <= stop,
           "arms": {m: d.get("timing", {}).get("s_per_step") for m, d in docs.items()}}  # fmt: skip
    return doc


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--model", required=True, choices=sorted(R1))
    ap.add_argument("--hf_rate", type=Path, help="R0's hf_rate.json for --model (the HF sampler or the HF check)")
    ap.add_argument("manifests", nargs="+", type=Path, help="train_sft manifest.json of R0 timing arms")
    a = ap.parse_args(argv)
    try:
        doc = project(a)
    except (OSError, ValueError, KeyError, TypeError) as e:
        print(f"project_r1: {e}", file=sys.stderr)
        return 2
    go, hours = doc["go"], doc["projected_h"]
    print(f"{'GO' if go else 'STOP'}: {a.model} R1 projected {hours:.2f} h "
          f"{'<=' if go else '>'} stop {doc['stop_h']:.2f} h (PLAN §4.9)")  # fmt: skip
    print(json.dumps(doc))
    return 0 if go else 1


if __name__ == "__main__":
    sys.exit(main())

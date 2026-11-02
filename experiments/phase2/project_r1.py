"""R0 stop rule (PLAN §4.7 step 9, §4.9): project R1 GPU-hours from the timing prelude; stop above 2x the high end.

    python -I -m experiments.phase2.project_r1 --model lfm2.5-vl-450m --hf_rate runs/v1/hf-rate-lfm/hf_rate.json \
        runs/v1/r0-timing-w6-det-off/manifest.json [...]

Basis: of the given train_sft manifests of --model, the one timing arm with deterministic off and 6 loader workers
(R1's loader; the Qwen3.5 arm is one too), trained as R1 is (lora, freeze_vision, gradient checkpointing equal); the
other arms are listed for the determinism / workers cost. R1 steps: train_sft.plan on R1's config (its branches) at
R1's steps per epoch = the manifest's rows_per_epoch (the same bundle view, uncapped) // its effective batch, which
must equal R1's (s/step is per batch); the stable run plus every decay.

Evaluation, by R1's `sampler`: vllm -> EVAL_H_VLLM (§4.9, UNVERIFIED); hf -> select_ckpt's item-samples in
policy.chain's R1 call (eval_samples: --items SUBSET for selection) x the R0 hf_rate step's s_per_sample, plus one
model load (load_s) per sampler child (sampler_calls). Not counted: CPU scoring, G-vision, G-probe.
Projected R1 h = steps x s/step / 3600 + eval h. Prints a verdict line and JSON; exit 0 = go, 1 = stop, 2 = usage or
data error (no usable basis, an arm trained unlike R1, no or another model's HF rate for an HF R1).
"""

import argparse
import json
import sys
from pathlib import Path

from experiments.phase2 import train_sft as TS

CONFIGS = Path(__file__).parent / "configs"
EVAL_H_VLLM = 0.3  # §4.9: evaluation allowance of R1 with vLLM (UNVERIFIED there)
# §4.9 R1 high estimates (LFM 0.76-1.70 h; Qwen3.5 taken as 2x LFM); the stop rule is 2x the high end
R1 = {"lfm2.5-vl-450m": ("r1_lfm.yaml", 1.70), "qwen3.5-0.8b": ("r1_qwen35.yaml", 3.40)}
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


def project(a) -> tuple[dict, float]:
    import yaml

    docs = {str(m): json.loads(m.read_text()) for m in a.manifests}
    docs = {m: d for m, d in docs.items() if d["config"]["model"] == a.model}
    basis = [m for m, d in docs.items() if not d["config"].get("deterministic")
             and d["config"]["sft"].get("dataloader_num_workers", 6) == 6]  # fmt: skip
    cfg, high = R1[a.model]
    r1 = TS.RunArgs(**yaml.safe_load((CONFIGS / cfg).read_text()))
    batch = r1.sft.get("per_device_train_batch_size", 8) * r1.sft.get("gradient_accumulation_steps", 1)
    if len(basis) != 1 or (t := docs[basis[0]].get("timing", {})).get("s_per_step") is None or t["batch"] != batch:
        raise ValueError(f"need one finished deterministic-off, workers-6 timing arm of {a.model} at batch {batch}")
    arm, want = docs[basis[0]]["config"], {"lora": r1.lora, "freeze_vision": r1.freeze_vision}
    got = {k: arm.get(k, False) for k in want} | {"gradient_checkpointing": _ckpt(arm)}
    if got != (want := want | {"gradient_checkpointing": _ckpt({"model": a.model, "sft": r1.sft})}):
        raise ValueError(f"{basis[0]}: timing arm trained with {got}, R1 ({cfg}) with {want}")
    p = TS.plan(r1, t["rows_per_epoch"] // batch)
    steps = p["stable"] + sum(b["decay"] for b in p["branches"].values())
    sampler = r1.sampler or "vllm"
    ev = {"sampler": sampler, "eval_h": EVAL_H_VLLM}
    if sampler == "hf":
        if a.hf_rate is None:
            raise ValueError(f"{cfg} samples with hf: pass --hf_rate (R0's hf-rate-*/hf_rate.json)")
        rate = json.loads(a.hf_rate.read_text())
        if rate["model"] != a.model or not rate["s_per_sample"] > 0:
            raise ValueError(f"{a.hf_rate}: rate of {rate['model']} ({rate['s_per_sample']} s), want {a.model}")
        n, calls = eval_samples(len(r1.branches), SUBSET), sampler_calls(len(r1.branches), SUBSET)
        ev = {"sampler": "hf", "eval_samples": n, "s_per_sample": rate["s_per_sample"], "sampler_calls": calls,
              "eval_h": (n * rate["s_per_sample"] + calls * rate["load_s"]) / 3600}  # fmt: skip
    hours = steps * t["s_per_step"] / 3600 + ev["eval_h"]
    doc = {"model": a.model, "basis": basis[0], "s_per_step": t["s_per_step"], "r1_steps": steps, **ev,
           "projected_h": hours, "high_h": high, "stop_h": 2 * high, "go": hours <= 2 * high,
           "arms": {m: d.get("timing", {}).get("s_per_step") for m, d in docs.items()}}  # fmt: skip
    return doc, high


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--model", required=True, choices=sorted(R1))
    ap.add_argument("--hf_rate", type=Path, help="R0's hf_rate.json for --model (required when R1 samples with hf)")
    ap.add_argument("manifests", nargs="+", type=Path, help="train_sft manifest.json of R0 timing arms")
    a = ap.parse_args(argv)
    try:
        doc, high = project(a)
    except (OSError, ValueError, KeyError, TypeError) as e:
        print(f"project_r1: {e}", file=sys.stderr)
        return 2
    go, hours = doc["go"], doc["projected_h"]
    print(f"{'GO' if go else 'STOP'}: {a.model} R1 projected {hours:.2f} h "
          f"{'<=' if go else '>'} stop {2 * high:.2f} h (2 x the PLAN §4.9 high estimate {high:.2f} h)")  # fmt: skip
    print(json.dumps(doc))
    return 0 if go else 1


if __name__ == "__main__":
    sys.exit(main())

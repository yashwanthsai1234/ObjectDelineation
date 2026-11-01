"""R0 stop rule (PLAN §4.7 step 9, §4.9): project R1 GPU-hours from the timing prelude; stop above 2x the high end.

    python -I -m experiments.phase2.project_r1 --model lfm2.5-vl-450m runs/v1/r0-timing-w6-det-off/manifest.json [...]

Basis: of the given train_sft manifests of --model, the one timing arm with deterministic off and 6 loader workers
(R1's loader; the Qwen3.5 arm is one too); the other arms are listed for the determinism / workers cost. R1 steps:
train_sft.plan on R1's config (its branches) at R1's steps per epoch = the manifest's rows_per_epoch (the same
bundle view, uncapped) // its effective batch, which must equal R1's (s/step is per batch); the stable run plus
every decay. Projected R1 h = steps x s/step / 3600 + EVAL_H. Prints a verdict line and JSON; exit 0 = go, 1 = stop.
"""

import argparse
import json
import sys
from pathlib import Path

from experiments.phase2 import train_sft as TS

CONFIGS = Path(__file__).parent / "configs"
EVAL_H = 0.3  # §4.9: evaluation allowance of R1 (UNVERIFIED there)
# §4.9 R1 high estimates (LFM 0.76-1.70 h; Qwen3.5 taken as 2x LFM); the stop rule is 2x the high end
R1 = {"lfm2.5-vl-450m": ("r1_lfm.yaml", 1.70), "qwen3.5-0.8b": ("r1_qwen35.yaml", 3.40)}


def main(argv=None) -> int:
    import yaml

    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--model", required=True, choices=sorted(R1))
    ap.add_argument("manifests", nargs="+", type=Path, help="train_sft manifest.json of R0 timing arms")
    a = ap.parse_args(argv)
    docs = {str(m): json.loads(m.read_text()) for m in a.manifests}
    docs = {m: d for m, d in docs.items() if d["config"]["model"] == a.model}
    basis = [m for m, d in docs.items() if not d["config"].get("deterministic")
             and d["config"]["sft"].get("dataloader_num_workers", 6) == 6]  # fmt: skip
    cfg, high = R1[a.model]
    r1 = TS.RunArgs(**yaml.safe_load((CONFIGS / cfg).read_text()))
    batch = r1.sft.get("per_device_train_batch_size", 8) * r1.sft.get("gradient_accumulation_steps", 1)
    if len(basis) != 1 or (t := docs[basis[0]].get("timing", {})).get("s_per_step") is None or t["batch"] != batch:
        raise SystemExit(f"need one finished deterministic-off, workers-6 timing arm of {a.model} at batch {batch}")
    p = TS.plan(r1, t["rows_per_epoch"] // batch)
    steps = p["stable"] + sum(b["decay"] for b in p["branches"].values())
    hours = steps * t["s_per_step"] / 3600 + EVAL_H
    go = hours <= 2 * high
    print(f"{'GO' if go else 'STOP'}: {a.model} R1 projected {hours:.2f} h "
          f"{'<=' if go else '>'} stop {2 * high:.2f} h (2 x the PLAN §4.9 high estimate {high:.2f} h)")  # fmt: skip
    print(json.dumps({"model": a.model, "basis": basis[0], "s_per_step": t["s_per_step"], "r1_steps": steps,
                      "eval_h": EVAL_H, "projected_h": hours, "high_h": high, "stop_h": 2 * high, "go": go,
                      "arms": {m: d.get("timing", {}).get("s_per_step") for m, d in docs.items()}}))  # fmt: skip
    return 0 if go else 1


if __name__ == "__main__":
    sys.exit(main())

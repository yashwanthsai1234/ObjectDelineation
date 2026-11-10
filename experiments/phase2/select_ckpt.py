"""D12 checkpoint selection CLI (PLAN D12, D13, §3.2 ②, §5.1): objdel.sft.select scores the branches, this script
applies rule R-B and the gates; never touches CUDA.

    python -I -m experiments.phase2.select_ckpt --run runs/r1 --shifted data/eval_views/lidc-val-shifted \
        --vision '{"prepared": root, "reference": R0 base value}' | --no-vision \
        [--sampler vllm|hf] [--items 256] [--rescaled V] [--rotated V] [--prepared P] [--probes DIR]
    python -I -m experiments.phase2.select_ckpt --run runs/c2 --probes DIR --only topology [--sampler hf]
    python -I -m experiments.phase2.select_ckpt --spec specs/r1_qwen3vl2b.yaml --stage select [--rd RD] [--data DATA]

The flags are an objdel.core.spec.SelectStage (cli_stage): K is -k and the sampler timeout --timeout-s; --max-tokens
caps every sampled completion below the model's budget (the CPU smoke on a random model, which never stops). The order of
work, the sampler child and the outputs: objdel.sft.select.
"""

import argparse
import json
import sys
from pathlib import Path

from objdel.geom import engine as _gpucheck
from objdel.sft import select
from objdel.verify import gates
from objdel.verify import store as serve


def judge(a, sel: dict, out: Path) -> dict:
    """The §5.1 report of the selected (or gated) branch, on the full views (--items never applies here, D13)."""
    from experiments.phase2 import baselines
    from objdel.verify import stats as signal

    best, whole = a.run / "branches" / sel["name"], argparse.Namespace(**vars(a) | {"items": None})

    def full(view, sub, k=a.k):
        return select.score_view(whole, best, view, out / sub, k)

    scored = {"shifted": full(a.shifted, "shifted-full") if a.items else sel}
    scored |= {"rescaled": full(a.rescaled, "rescaled")} if a.rescaled else {}
    report = {"views": {v: {"summary": r["summary"]} for v, r in scored.items()}, "vision": sel["vision"]}
    report |= {"signal": scored["shifted"]["panel"], "k16": full(a.shifted, "shifted-k16", 16)["summary"]}
    if a.rotated:  # free attribution check: rotated ~ identity means the dihedral training transferred
        rot = full(a.rotated, "rotated")
        diff = signal.paired(rot["panel"], scored["shifted"]["panel"], "iou_median")
        report["attribution"] = {"summary": rot["summary"], "iou_median_rotated_minus_identity": diff}
    if a.prepared:
        fit = baselines.fit_prior(a.prepared, "lidc", "train")
        views = {"shifted": a.shifted, "rescaled": a.rescaled}
        report["prior"] = {v: select.prior_ci(r["iou"], baselines.score_prior(views[v], "lidc", "val", fit))
                           for v, r in scored.items()}  # fmt: skip
    if a.probes:
        report["probes"] = {p: select.score_probe(whole, best, p, out) for p in ("precision", "topology")}
    return report


def dump(path: Path, doc: dict) -> None:
    """Indented JSON to path (parents created)."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(doc, indent=1, default=str) + "\n")


def topology_only(a, out: Path) -> int:
    """--only topology (C2): the topology probe on the run's single branch -> <out>/report.json."""
    names = sorted(p.name for p in (a.run / "branches").iterdir())
    if a.probes is None or len(names) != 1:
        raise SystemExit(f"select_ckpt --only topology: needs --probes and one branch, found {names}")
    topo = select.score_probe(a, a.run / "branches" / names[0], "topology", out / names[0])
    dump(out / "report.json", {"only": "topology", "branch": names[0], "sampler": a.sampler,
        "checker": _gpucheck.build(), "k": a.k,
                               "probes": {"topology": topo}})
    print(json.dumps({"branch": names[0], "pass": topo.get("pass")}))
    return 0


def vision_config(a) -> dict | None:
    """The G-vision config of the stage (None with --no-vision); SystemExit when prepared or reference is bad."""
    if a.no_vision:
        return None
    cfg = a.vision
    ref = cfg.get("reference")
    ok = Path(ref).is_file() if isinstance(ref, str) else isinstance(ref, int | float) and ref > 0
    if not isinstance(cfg.get("prepared"), str) or not ok:
        raise SystemExit(f"select_ckpt: --vision needs prepared (a path) and reference (> 0, or R0's JSON): {cfg}")
    return cfg


PATHS = ("run", "shifted", "rescaled", "rotated", "prepared", "probes", "out")


def main(argv=None) -> int:
    """CLI: the SelectStage -> score every branch, judge, write <out>/selected.json."""
    from objdel.core.spec import SelectStage, cli_stage

    st = cli_stage(SelectStage, sys.argv[1:] if argv is None else argv)
    a = argparse.Namespace(**{k: Path(v) if k in PATHS and v is not None else v
                              for k, v in st.model_dump(exclude={"kind", "name", "soft"}).items()})
    a.sampler = st.sampler or "vllm"
    vision_cfg = None if a.only else vision_config(a)
    a.model, out = a.model or select.exited(a.run)["config"]["model"], a.out or a.run / "select"
    if a.only:
        return topology_only(a, out)
    vision_cfg = vision_cfg and {"model": a.model} | vision_cfg  # sft.train.reference checks R0's JSON model
    for v in filter(None, (a.shifted, a.rescaled, a.rotated)):
        serve.preload(v, ["lidc"])
    branches = []
    for b in sorted((a.run / "branches").iterdir(), key=lambda p: float(p.name[1:])):
        res = select.score_view(a, b, a.shifted, out / b.name / "shifted", a.k)
        res["vision"] = select.vision_check(b, vision_cfg)
        res["eligible"] = gates.eligible(res["summary"], res["vision"])  # D12: G-format, G-signal, G-vision
        branches.append({"name": b.name, "depth": float(b.name[1:]), **res})
    from experiments.phase2 import prereg

    name, why = prereg.select_rb(branches)
    gate, judged = None, name or (branches[-1]["name"] if branches else None)  # none eligible: the deepest (§6.3)
    if judged:
        report = judge(a, next(x for x in branches if x["name"] == judged), out / judged)
        gate = prereg.entry(report, ["lidc"])  # phase2_entry with the pre-registered G-prior (a) slot
        if a.sampler == "vllm":  # D13: the gates stand only if vLLM's answers score like HF's
            report["hf_check"] = check = select.hf_check(a, judged, out)
            gate["gates"]["G-hf-check"] = {"requirement": "every vLLM - HF metric difference: 0 inside its 95 % CI",
                                           "values": check["vllm_minus_hf"], "pass": check["pass"]}
            gate["passed"]["G-hf-check"], gate["pass"] = check["pass"], gate["pass"] and check["pass"]
        dump(out / "report.json", report)
        dump(out / "gates.json", gate)
    torch = sys.modules.get("torch")
    dump(out / "selected.json", {
        "cuda_initialized": bool(torch and torch.cuda.is_initialized()),
        "run": str(a.run),  # g_probe --select reads the branches from here
        "selected": name, "reasons": why, "gates_branch": judged, "gates_pass": gate and gate["pass"],
        "gates_failed": [n for n, g in gate["gates"].items() if not g["pass"]] if gate else [],
        "sampler": a.sampler, "checker": _gpucheck.build(), "k": a.k, "items": a.items,
        "branches": {x["name"]: {k: x[k] for k in ("depth", "eligible", "vision")}
                     | {"summary": x["summary"]["datasets"]["lidc"], "panel": x["panel"]["summary"]} for x in branches},
    })  # fmt: skip
    print(json.dumps({"selected": name, "eligible": why["eligible"]}))
    return 0

if __name__ == "__main__":
    sys.exit(main())

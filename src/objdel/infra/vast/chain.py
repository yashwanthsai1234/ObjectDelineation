"""A run spec's (or an old plan's) stages on the box as argv lists (`chain`), the G-vision references it reads
(`vision_refs`) and the folder `push_hf` uploads (`push_target`)."""

import json
from pathlib import Path

from objdel.core.provenance import repo_root
from objdel.core.spec import RunSpec, fill, load_spec, to_argv
from objdel.infra.vast.remote import die

DATA = "/workspace/data/lidc"  # {data} on the box


def chain(cfg: Path, run_id: str) -> list[list[str]]:
    """A config's stages on the box, as argv lists (`-m module ...` runs in /opt/venv, `bash script ...` as is).
    `steps:` (R0) lists them literally with {rd} / {data} filled in; a leading "~" lets a step fail without stopping
    the chain. Otherwise `stages:` (default train, select, probe; C2: train, topology) picks from train_sft, then on
    the saved branches select_ckpt (G-vision block, rotated view, `sampler:`) or its topology-only form, and g_probe;
    a config without branches (timing arms) only trains. A chain that scores (select, topology or probe) first builds
    the CUDA checker (run_box.sh --build-only) and checks the scorers on it against the host scorer
    (tests/test_gpu_scoring.py),
    so a failed build or a difference stops the chain before training: scoring runs on the GPU only. Box
    paths: data /workspace/data/lidc, runs /workspace/runs/<id>."""
    import yaml

    c, data, rd = yaml.safe_load(Path(cfg).read_text()), DATA, f"/workspace/runs/{run_id}"
    if "schema_version" in c:  # a run spec (specs/*.yaml)
        return [to_argv(s, rd, data) for s in load_spec(Path(cfg)).stages]
    if "steps" in c:
        return [[str(x).format(rd=rd, data=data) for x in s] for s in c["steps"]]
    run, ev, prep = f"{rd}/{c['output_dir']}", f"{data}/eval_views", f"{data}/prepared"
    v, samp = c.get("vision"), ["--sampler", c["sampler"]] if c.get("sampler") else []
    items = ["--items", "256"] if c.get("sampler") == "hf" else []  # D13: HF path selects on 256 shifted-val items
    ref = v and (f"{rd}/{v['reference']}" if isinstance(v["reference"], str) else v["reference"])  # a number as is
    vis = ["--vision", json.dumps({"prepared": f"{data}/{v['prepared']}", "reference": ref})] if v else ["--no-vision"]
    select = ["-m", "experiments.phase2.select_ckpt", "--run", run, "--prepared", prep, "--probes", ev, *samp]
    table = {
        "train": ["-m", "experiments.phase2.train_sft", "--config",
                  Path(cfg).resolve().relative_to(repo_root()).as_posix()],
        "select": [*select, "--shifted", f"{ev}/lidc-val-shifted", "--rescaled", f"{ev}/lidc-val-rescaled",
                   "--rotated", f"{ev}/lidc-val-rotated", *items, *vis],
        "topology": [*select, "--only", "topology"],
        "probe": ["-m", "experiments.phase2.g_probe", "--model", c["model"], "--select", f"{run}/select", "--train",
                  prep, "--probe", f"{ev}/lidc-val-shifted", "--out", f"{run}/g_probe"],
    }  # fmt: skip
    stages = c.get("stages", ["train", "select", "probe"]) if c.get("branches") or c.get("branches_from") else ["train"]
    gpu = [["bash", "-c", "bash experiments/cuda_checker/run_box.sh --build-only && /opt/venv/bin/python -m pytest -q "
                          "-p no:cacheprovider tests/test_gpu_scoring.py"]] if set(stages) - {"train"} else []
    return gpu + [table[s] for s in stages]


def vision_refs(cfg: Path, run_id: str) -> list[str]:
    """The G-vision reference JSON a config reads on the box ([] for a number or `vision: off`)."""
    import yaml

    c = yaml.safe_load(Path(cfg).read_text())
    if "schema_version" in c:  # a run spec: the reference of each sft and select stage, in stage order, once
        rd, refs = f"/workspace/runs/{run_id}", []
        for s in load_spec(Path(cfg)).stages:
            if s.kind == "sft" and s.vision and isinstance(s.vision["reference"], str):
                refs.append(f"{rd}/{s.vision['reference']}")
            elif s.kind == "select" and s.vision and isinstance(s.vision.reference, str):
                refs.append(fill(s.vision.reference, rd, DATA))
        return list(dict.fromkeys(refs))
    ref = (c.get("vision") or {}).get("reference")
    return [f"/workspace/runs/{run_id}/{ref}"] if isinstance(ref, str) else []


def push_target(spec: RunSpec, run_id: str) -> tuple[str, str, bool]:
    """(folder on the box, repo revision, is an rl output) of a spec's last rl or sft stage. For rl the revision is the
    last part of its out folder; for sft the output_dir."""
    last = next((s for s in reversed(spec.stages) if s.kind in ("rl", "sft")), None)
    if last is None:
        die("nothing to push: the spec has no rl or sft stage")
    rd = f"/workspace/runs/{run_id}"
    if last.kind == "rl":
        path = fill(last.out, rd, DATA)
        return path, path.rsplit("/", 1)[1], True
    return f"{rd}/{last.output_dir}", last.output_dir, False

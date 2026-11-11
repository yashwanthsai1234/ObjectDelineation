"""A run spec's stages on the box as argv lists (`chain`), the G-vision references it reads (`vision_refs`) and the
folder `push_hf` uploads (`push_target`)."""

from pathlib import Path

from objdel.core.provenance import repo_root
from objdel.core.spec import RunSpec, fill, load_spec, to_argv
from objdel.infra.errors import InfraError

DATA = "/workspace/data/lidc"  # {data} on the box


def chain(cfg: Path, run_id: str) -> list[list[str]]:
    """A run spec's stages on the box, as argv lists (`-m module ...` runs in /opt/venv, `bash script ...` as is; a
    leading "~" lets a stage fail without stopping the chain). Box paths: data /workspace/data/lidc, runs
    /workspace/runs/<id>. Only a spec in specs/: an old `steps:` plan raises ValidationError (load_spec), a spec
    elsewhere raises ValueError, both before any remote call."""
    spec = load_spec(Path(cfg))
    if Path(cfg).resolve().parent != (repo_root() / "specs").resolve():
        raise ValueError(f"{cfg}: the box chain runs only a run spec in specs/")
    return [to_argv(s, f"/workspace/runs/{run_id}", DATA) for s in spec.stages]


def vision_refs(cfg: Path, run_id: str) -> list[str]:
    """The G-vision reference JSON a run spec reads on the box: the reference of each sft and select stage, in stage
    order, once ([] for a number or no vision)."""
    rd, refs = f"/workspace/runs/{run_id}", []
    for s in load_spec(Path(cfg)).stages:
        if s.kind == "sft" and s.vision and isinstance(s.vision["reference"], str):
            refs.append(f"{rd}/{s.vision['reference']}")
        elif s.kind == "select" and s.vision and isinstance(s.vision.reference, str):
            refs.append(fill(s.vision.reference, rd, DATA))
    return list(dict.fromkeys(refs))


def push_target(spec: RunSpec, run_id: str) -> tuple[str, str, bool]:
    """(folder on the box, repo revision, is an rl output) of a spec's last rl or sft stage. For rl the revision is the
    last part of its out folder; for sft the output_dir."""
    last = next((s for s in reversed(spec.stages) if s.kind in ("rl", "sft")), None)
    if last is None:
        raise InfraError("nothing to push: the spec has no rl or sft stage")
    rd = f"/workspace/runs/{run_id}"
    if last.kind == "rl":
        path = fill(last.out, rd, DATA)
        return path, path.rsplit("/", 1)[1], True
    return f"{rd}/{last.output_dir}", last.output_dir, False

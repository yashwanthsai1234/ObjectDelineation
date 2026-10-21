"""Phase 2 SFT with WSD decay branches (plan D5, D7, §3.2 ①): one stable run, then one decay continuation per branch.

    python -I -m experiments.phase2.train_sft --config experiments/phase2/configs/r1_lfm.yaml [--resume true]
    python -I -m experiments.phase2.train_sft --config <same> --vision_reference_only true   # R0: G-vision base ref

Schedule (D5, configs.sft_config): warmup (3 % of the stable steps), constant LR, and for each branch point b
(epochs, default 0.5 / 1 / 2) a linear decay to 0 over 15 % of b. The stable run uses the deepest branch's schedule
and stops at its branch point, so every saved branch point sits at LR factor 1; it saves a full checkpoint (model,
optimizer, scheduler, RNG, data position) at every b, and each decay branch is a new trainer resumed from it with
the same warmup_stable_decay schedule ending at b + d (sft_config(steps, branch_at=b)): on a tiny model every
branch is bit-identical to an uninterrupted WSD run that decays at b (impl/M8.md). Branch points that round to one
step, or a decay of 0 steps, are refused. `branches: []` is a stable-only run of one epoch that saves nothing (the
R0 timing prelude and Qwen3.5 arm). Kept: `stable/e<b>` and `branches/e<b>`, model-only with the processor (the
sampler loads them); LoRA runs (D6: Qwen3-VL-2B, or `lora: true`) are saved merged. Rolling checkpoints (optimizer
included, limit 3) every `rolling_minutes`; at each one the G-vision check (objdel.eval.vision, §5.1) on the
training device and dtype (recorded; the reference is CPU float32, B1): ratio < 0.5 **and** TV < 0.002 aborts with
the reason in manifest.json (exit code 3; D7 repair: rerun with tower and projector frozen), one alone warns.

G-vision is required (D7, D12): `vision: {prepared, reference}` with reference = the base model's projector
distance, a number or the JSON that `--vision_reference_only true` writes in R0 (§4.7 step 8: base model, CPU,
float32, on vision.fixed_set). `vision: off` only for the discarded timing runs; `vision: {}` is refused.

Paths: view and vision.prepared resolve under data_root (default $OBJDEL_DATA_ROOT, else /workspace/data/lidc on a
Vast box, else data); output_dir and a vision.reference JSON under run_dir (default $OBJDEL_RUN_DIR, else runs).
On a box (/workspace/runs exists) the output must land under /workspace/runs, where pull.sh reads.

Data: objdel.train.dataset.view_rows over the view directory; the rows of view epoch e are visited in epoch order
(length-grouped within each epoch, seeded by data_seed + e), so "b epochs" means the first b passes over fresh
views (D10). A branch's decay continues the same stream; the view must hold enough epochs for max(b) + its decay
(R1: build_views --epochs 3), checked before anything loads. The process never samples: select_ckpt starts after
it exits (D13). manifest.json records the config, the view hash, SHA-256 of the code, the step plan, the G-vision
log, the pid and the exit status. `--resume true` continues a crashed run (status "running", same training config
and view): the newest complete rolling checkpoint, then each branch not yet saved (from its own decay checkpoint
when there is one).

References:
- transformers 5.20 get_wsd_schedule, Trainer resume (optimizer, scheduler, RNG, data skip), LengthGroupedSampler:
  https://github.com/huggingface/transformers/blob/v5.20.1/src/transformers/optimization.py ,
  https://github.com/huggingface/transformers/blob/v5.20.1/src/transformers/trainer.py
- TRL 1.16 SFTTrainer, TrlParser (--config YAML): https://huggingface.co/docs/trl/main/en/sft_trainer
- PEFT merge_and_unload: https://huggingface.co/docs/peft/main/en/package_reference/lora
- WSD: Hu et al., MiniCPM (2404.06395); Hägele et al. (2405.18392)
"""

import hashlib
import json
import math
import os
import shutil
import socket
import subprocess
import sys
import time
from dataclasses import asdict, dataclass, field
from pathlib import Path

ABORTED = 3  # exit code of a G-vision abort
ROOT = Path(__file__).resolve().parents[2]
BOX_DATA, BOX_RUNS = Path("/workspace/data/lidc"), Path("/workspace/runs")  # sync_up / pull.sh on a Vast box (M7)
SAME_RUN = ("model", "branches", "a0", "freeze_vision", "lora", "seed", "steps_per_epoch", "sft")  # resume check
TRAINER_FILES = ("optimizer", "scheduler", "rng_state", "trainer_state", "adapter_")


def _data_root() -> str:
    return os.environ.get("OBJDEL_DATA_ROOT") or (str(BOX_DATA) if BOX_DATA.is_dir() else "data")


@dataclass
class RunArgs:
    model: str = field(metadata={"help": "objdel.train.configs.MODELS key"})
    view: str = field(metadata={"help": "training view directory, relative to data_root (views/<view>)"})
    output_dir: str = field(metadata={"help": "run directory, relative to run_dir"})
    data_root: str = field(default_factory=_data_root, metadata={"help": "$OBJDEL_DATA_ROOT, Vast bundle dir or data"})
    run_dir: str = field(
        default_factory=lambda: os.environ.get("OBJDEL_RUN_DIR", "runs"),
        metadata={"help": "$OBJDEL_RUN_DIR (Vast: /workspace/runs/<RUN_ID>, set by run.sh) or runs"},
    )
    model_path: str | None = field(default=None, metadata={"help": "local weights (default MODELS[model])"})
    branches: list[float] = field(default_factory=lambda: [0.5, 1.0, 2.0], metadata={"help": "epochs ([] = timing)"})
    a0: bool = field(default=False, metadata={"help": "empty every holes list (A0)"})
    freeze_vision: bool = field(default=False, metadata={"help": "D7 repair: tower and projector frozen"})
    lora: bool = field(default=False, metadata={"help": "D6: LM-only LoRA r16 (Qwen3.5 if R0's peak memory says so)"})
    seed: int = 1
    rolling_minutes: float = 10.0
    vision: dict = field(default_factory=dict, metadata={"help": "{prepared, reference} for G-vision, or off"})
    vision_reference_only: bool = field(default=False, metadata={"help": "R0: write vision.reference and exit"})
    resume: bool = field(default=False, metadata={"help": "continue a crashed run in output_dir"})
    steps_per_epoch: int | None = field(default=None, metadata={"help": "override (smoke tests, timing)"})
    sft: dict = field(default_factory=dict, metadata={"help": "overrides of configs.sft_config"})

    def __post_init__(self):  # a bad G-vision config fails here, not at the first rolling checkpoint (BK3-5)
        if self.vision in (False, "off"):
            self.vision = False
        elif self.vision:
            ref = self.vision.get("reference")
            num = isinstance(ref, (int, float)) and not isinstance(ref, bool) and math.isfinite(ref) and ref > 0
            if set(self.vision) != {"prepared", "reference"} or not (num or str(ref).endswith(".json")):
                raise ValueError(
                    f"vision {self.vision!r}: need {{prepared, reference}}, reference = R0's base projector distance "
                    "(a positive number, or the .json that --vision_reference_only writes)"
                )


def plan(args: RunArgs, steps_per_epoch: int) -> dict:
    """The step plan in optimizer steps: the stable run's length, and per branch its point and decay (as
    configs.sft_config computes them: warmup = 3 % of the stable steps, decay = 15 % of the branch point)."""
    from objdel.train.configs import DECAY, WARMUP

    at = {f"e{b:g}": max(1, round(b * steps_per_epoch)) for b in args.branches}
    br = {n: {"at": s, "decay": round(DECAY * s)} for n, s in at.items()}
    if len(set(at.values())) < len(at) or any(b["decay"] < 1 for b in br.values()):
        raise ValueError(f"branches {br} at {steps_per_epoch} steps/epoch: points collide or a decay is 0 steps")
    stable = max(at.values(), default=steps_per_epoch)
    return {"steps_per_epoch": steps_per_epoch, "stable": stable, "warmup": round(WARMUP * stable), "branches": br}


def resolve(args: RunArgs) -> RunArgs:
    """view and vision.prepared under data_root; output_dir and a vision.reference JSON under run_dir."""
    args.view, args.output_dir = str(Path(args.data_root, args.view)), str(Path(args.run_dir, args.output_dir))
    if args.vision:
        ref = args.vision["reference"]
        args.vision = {
            "prepared": str(Path(args.data_root, args.vision["prepared"])),
            "reference": str(Path(args.run_dir, ref)) if isinstance(ref, str) else ref,
        }
    return args


def vision_set(prepared) -> tuple[list[dict], str]:
    from objdel.eval import vision

    return vision.fixed_set(prepared)


def reference(ref, model: str | None, set_sha256: str) -> float:
    """The G-vision base reference: a number, or R0's JSON for this model (when given) and fixed set."""
    if not isinstance(ref, str):
        return float(ref)
    doc = json.loads(Path(ref).read_text())
    if (doc["model"], doc["set_sha256"]) != (model or doc["model"], set_sha256):
        raise ValueError(f"{ref}: reference of {doc['model']} on set {doc['set_sha256'][:12]}, not this run's")
    return doc["reference"]


def vision_reference(args: RunArgs) -> int:
    """R0 (§4.7 step 8): projector_distance of the base model on the fixed set, CPU float32 (B1) -> vision.reference."""
    import torch
    from transformers import AutoModelForImageTextToText, AutoProcessor

    from objdel.eval import vision
    from objdel.train import configs

    if not (args.vision and isinstance(args.vision["reference"], str)):
        raise ValueError("--vision_reference_only writes vision.reference: set it to a .json path")
    src = args.model_path or configs.MODELS[args.model]
    proc = AutoProcessor.from_pretrained(src, **configs.processor_kwargs(args.model, "lidc"))
    model = AutoModelForImageTextToText.from_pretrained(src, dtype=torch.float32).eval()
    items, sha = vision_set(args.vision["prepared"])
    with torch.no_grad():
        d = vision.projector_distance(model, proc, [it["image"] for it in items[: vision.N_IMAGES]])
    dest = Path(args.vision["reference"])
    dest.parent.mkdir(parents=True, exist_ok=True)
    doc = {"model": args.model, "source": src, "reference": d, "set_sha256": sha, "device": "cpu", "dtype": "float32"}
    write_json(dest, doc | {"code": code_hashes()})
    return 0


def vision_check(model, processor, cfg: dict) -> dict:
    """objdel.eval.vision.check on the fixed set: cfg = train's preloaded {"items", "set_sha256", "reference"}, or
    {"prepared", "reference"} (select_ckpt)."""
    import torch

    from objdel.eval import vision

    if "items" not in cfg:
        items, sha = vision_set(cfg["prepared"])
        cfg = {"items": items, "set_sha256": sha, "reference": reference(cfg["reference"], cfg.get("model"), sha)}
    images = [it["image"] for it in cfg["items"][: vision.N_IMAGES]]
    with torch.no_grad():
        res = vision.check(model, processor, images, cfg["items"], reference=cfg["reference"])
    return res | {"set_sha256": cfg["set_sha256"]}


def record_vision(man: dict, step: int, res: dict) -> bool:
    """Log one G-vision result in the manifest; True (stop training) on "abort": ratio and TV both fail (D7)."""
    keys = ("action", "ratio", "tv", "tv_black", "tv_other", "set_sha256", "device", "dtype")
    man["vision"].append({"step": step} | {k: res.get(k) for k in keys})
    if res["action"] == "abort":
        man["abort"] = {
            "step": step,
            "reason": f"G-vision (a) ratio {res['ratio']:.4g} < 0.5 and (b) TV "
            f"{res['tv']:.4g} < 0.002 (D7 repair: rerun with freeze_vision)",
        }
    return res["action"] == "abort"


def sha256(path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def code_hashes() -> dict:
    files = sorted((ROOT / "src" / "objdel").rglob("*.py")) + [Path(__file__)]
    h = hashlib.sha256(b"".join(p.relative_to(ROOT).as_posix().encode() + bytes.fromhex(sha256(p)) for p in files))
    git = subprocess.run(["git", "-C", str(ROOT), "rev-parse", "HEAD"], capture_output=True, text=True)
    return {"objdel_and_script_sha256": h.hexdigest(), "git_head": git.stdout.strip() or None, "files": len(files)}


def write_json(path, doc) -> None:
    tmp = Path(f"{path}.tmp")
    tmp.write_text(json.dumps(doc, indent=1, sort_keys=True, default=str) + "\n")
    os.replace(tmp, path)


def latest(d) -> str | None:
    """The newest complete Trainer checkpoint in d (trainer_state.json is written last)."""
    done = [c for c in Path(d).glob("checkpoint-*") if (c / "trainer_state.json").exists()]
    return str(max(done, key=lambda c: int(c.name.rsplit("-", 1)[1]))) if done else None


def replace_dir(write, dest: Path) -> None:
    """write(tmp) then rename: a crash never leaves a half-written kept directory (resume trusts what exists)."""
    tmp = dest.with_name(f".{dest.name}.tmp")
    shutil.rmtree(tmp, ignore_errors=True)
    write(tmp)
    shutil.rmtree(dest, ignore_errors=True)
    os.replace(tmp, dest)


def main(argv=None) -> int:
    from trl import TrlParser

    (args,) = TrlParser([RunArgs]).parse_args_and_config(argv)
    return train(args)


def train(args: RunArgs) -> int:
    import torch
    from transformers import AutoModelForImageTextToText, AutoProcessor, TrainerCallback
    from transformers.trainer_pt_utils import get_length_grouped_indices
    from trl import SFTTrainer

    from objdel.train import configs, dataset

    resolve(args)
    if args.vision_reference_only:
        return vision_reference(args)
    if args.vision == {}:
        raise ValueError("vision: {} skips G-vision (D7, D12): give {prepared, reference}; `vision: off` for timing")
    out = Path(args.output_dir)
    if BOX_RUNS.is_dir() and BOX_RUNS.resolve() not in out.resolve().parents:
        raise ValueError(f"{out}: on a Vast box outputs go under {BOX_RUNS}/<RUN_ID> (pull.sh): set OBJDEL_RUN_DIR")
    view = json.loads(Path(args.view, "manifest.json").read_text())
    deep = max((b * (1 + configs.DECAY) for b in args.branches), default=1.0)
    if view["config"]["epochs"] < math.ceil(deep):  # S1-8, before anything loads
        raise ValueError(
            f"{args.view}: {view['config']['epochs']} view epoch(s), the deepest branch decays to {deep:g} epochs: "
            f"rebuild it with build_views --epochs {math.ceil(deep)}"
        )
    vis = None
    if args.vision:
        items, sha = vision_set(args.vision["prepared"])
        vis = {"items": items, "set_sha256": sha, "reference": reference(args.vision["reference"], args.model, sha)}
    out.mkdir(parents=True, exist_ok=True)
    old = json.loads((out / "manifest.json").read_text()) if (out / "manifest.json").exists() else None
    if old and not args.resume:
        raise ValueError(f"{out} holds a run: --resume true continues it, or pick another output_dir")
    man = {
        "config": asdict(args),
        "view": {
            "dir": args.view,
            "view_hash": view.get("view_hash"),
            "manifest_sha256": sha256(Path(args.view, "manifest.json")),
        },
        "code": code_hashes(),
        "vision": [],
        "log": {},
    }
    if args.resume:
        same = old and all(old["config"][k] == man["config"][k] for k in SAME_RUN)
        if not (same and old["status"] == "running" and old["view"]["view_hash"] == man["view"]["view_hash"]):
            raise ValueError(f"--resume: {out} is not a crashed run of this config and view")
        man |= {k: old[k] for k in ("vision", "log", "started")} | {"resumed": old.get("resumed", []) + [time.time()]}
    man |= {
        "pid": os.getpid(),
        "host": socket.gethostname(),
        "status": "running",
        "started": man.get("started", time.time()),
    }
    write_json(out / "manifest.json", man)
    src = args.model_path or configs.MODELS[args.model]
    proc = AutoProcessor.from_pretrained(src, **configs.processor_kwargs(args.model, "lidc"))
    rows = dataset.view_rows(args.view, proc.tokenizer, configs.MAX_COMPLETION[args.model], a0=args.a0)
    if "epoch" in rows.column_names:
        epochs = list(rows["epoch"])
    else:  # row_id = "<epoch>/<item_id>/<rater>/<rep>"
        epochs = [int(i.split("/", 1)[0]) for i in rows["row_id"]]
    lengths = list(rows["length"])

    configs.training_template(args.model, proc)  # ValueError if Qwen3.5's template is not TRL's no-think one
    peft = configs.lora(args.model, force=args.lora)
    dtype = getattr(torch, configs.dtype(args.model, lora=peft is not None))  # fp32 for full FT (bf16 autocast)

    def sft_args(out_dir, stable, branch_at=None):
        o = {"seed": args.seed, "data_seed": args.seed, "save_strategy": "no", "save_total_limit": 3}
        o |= {"logging_steps": 10}
        return configs.sft_config(args.model, str(out_dir), stable, branch_at, lora=peft is not None, **(o | args.sft))

    class Trainer(SFTTrainer):
        def _get_train_sampler(self, train_dataset=None):  # view epochs in order, length-grouped inside each
            bs = self.args.train_batch_size * self.args.gradient_accumulation_steps
            order = []
            for e in sorted(set(epochs)):
                idx = [i for i, x in enumerate(epochs) if x == e]
                g = torch.Generator().manual_seed(self.args.data_seed + e)
                order += [idx[j] for j in get_length_grouped_indices([lengths[i] for i in idx], bs, generator=g)]
            return order

    class Rolling(TrainerCallback):
        def __init__(self, keep_at=(), stop_at=None):
            self.keep_at, self.stop_at, self.last = dict(keep_at), stop_at, time.monotonic()

        def on_step_end(self, a, state, control, **kw):
            due = time.monotonic() - self.last >= 60 * args.rolling_minutes
            control.should_save = due or state.global_step in self.keep_at
            control.should_training_stop |= state.global_step == self.stop_at

        def on_save(self, a, state, control, model=None, processing_class=None, **kw):
            self.last, step = time.monotonic(), state.global_step
            ckpt = Path(a.output_dir, f"checkpoint-{step}")
            if step in self.keep_at:
                replace_dir(lambda tmp: shutil.copytree(ckpt, tmp), out / "stable" / self.keep_at[step])
            if vis:
                was, prm = model.training, next(model.parameters())
                res = vision_check(model, proc, vis) | {"device": str(prm.device), "dtype": str(prm.dtype)}
                control.should_training_stop |= record_vision(man, step, res)
                model.train(was)
                write_json(out / "manifest.json", man)

    def trainer(out_dir, stable, branch_at=None, keep_at=(), stop_at=None):
        model = AutoModelForImageTextToText.from_pretrained(src, dtype=dtype)
        cfg = sft_args(out_dir, stable, branch_at)
        opt = None
        if not peft:  # full FT: vision tower + projector at 0.1x (D7), or frozen (the D7 repair)
            lm, vision = configs.param_groups(model, cfg.learning_rate)
            for prm in vision["params"] if args.freeze_vision else ():
                prm.requires_grad_(False)
            groups = [lm] if args.freeze_vision else [lm, vision]
            opt = torch.optim.AdamW(groups, lr=cfg.learning_rate, weight_decay=cfg.weight_decay)
        kw = {"train_dataset": rows, "processing_class": proc, "peft_config": peft, "optimizers": (opt, None)}
        return Trainer(model=model, args=cfg, callbacks=[Rolling(keep_at, stop_at)], **kw)

    def save(model, dest):  # D6: adapters merged, so every kept directory is a full model with its processor
        model = model.merge_and_unload() if hasattr(model, "merge_and_unload") else model
        replace_dir(lambda tmp: (model.save_pretrained(tmp), proc.save_pretrained(tmp)), dest)

    def finish(d):  # a kept stable checkpoint: model-only (LoRA merged) with the processor, for the depth curve
        if (d / "adapter_config.json").exists():
            from peft import PeftModel

            model = AutoModelForImageTextToText.from_pretrained(src, dtype=dtype)
            PeftModel.from_pretrained(model, d).merge_and_unload().save_pretrained(d)
        for f in d.glob("*"):
            if f.name.startswith(TRAINER_FILES) or f.suffix == ".pt":
                f.unlink()
        proc.save_pretrained(d)

    gbs = args.sft.get("per_device_train_batch_size", 8) * args.sft.get("gradient_accumulation_steps", 1)  # 1 GPU
    p = plan(args, args.steps_per_epoch or sum(e == min(epochs) for e in epochs) // gbs)
    man["plan"] = p
    need = max([b["at"] + b["decay"] for b in p["branches"].values()] or [p["stable"]]) * gbs
    if need > len(rows):
        raise ValueError(f"{args.view}: {len(rows)} rows, the deepest branch and its decay need {need}")
    keep = {b["at"]: name for name, b in p["branches"].items()}
    if not (keep and (out / "stable" / keep[p["stable"]]).is_dir()):  # stable phase not finished yet
        t = trainer(out / "rolling", p["stable"], p["stable"] if keep else None, keep, p["stable"])
        t.train(resume_from_checkpoint=latest(out / "rolling") if args.resume else None)
        man["log"]["stable"] = t.state.log_history
        write_json(out / "manifest.json", man)
    for name, b in p["branches"].items() if "abort" not in man else ():
        if not (out / "branches" / name).is_dir():
            t = trainer(out / "decay" / name, p["stable"], branch_at=b["at"])
            t.train(resume_from_checkpoint=(args.resume and latest(out / "decay" / name)) or str(out / "stable" / name))
            save(t.model, out / "branches" / name)
            man["log"][name] = t.state.log_history
        finish(out / "stable" / name)  # frees the optimizer state once its branch exists
        write_json(out / "manifest.json", man)
        if "abort" in man:
            break
    for d in (out / "stable").glob("e*"):  # after an abort too
        finish(d)
    shutil.rmtree(out / "rolling", ignore_errors=True)
    shutil.rmtree(out / "decay", ignore_errors=True)
    man |= {"status": "aborted" if "abort" in man else "done", "finished": time.time()}
    write_json(out / "manifest.json", man)
    return ABORTED if "abort" in man else 0


if __name__ == "__main__":
    sys.exit(main())

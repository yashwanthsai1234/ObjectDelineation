"""SFT with WSD decay branches: one stable run, then one decay continuation per branch.

`train(stage)` runs one objdel.core.spec.SftStage; src/objdel/cli/train/sft.py is its CLI. The stage is frozen, so
resolve() returns a resolved copy.

Schedule: warmup 3 % of the stable steps, constant LR, and per branch point b (epochs, default 0.5 / 1 / 2) a linear
decay to 0 over 15 % of b. The stable run uses the deepest branch's schedule and stops at its point, so every branch
point sits at LR factor 1; it saves a full checkpoint (model, optimizer, scheduler, RNG, data position) at every b,
and each branch is a new trainer resumed from it whose schedule ends at b + d: bit-identical to an uninterrupted WSD
run decaying at b. Colliding branch points and 0-step decays are refused. `branches: []` is the stable run only, one
epoch, nothing saved. Kept: `stable/e<b>` and `branches/e<b>`, model-only with the processor, LoRA saved merged.

Rolling checkpoints (optimizer included, limit 3) every `rolling_minutes`, each with the G-vision check
(objdel.policy.vision) on the training device and dtype: ratio < 0.5 and TV < 0.002 aborts with the reason in
manifest.json (exit 3), one alone warns. G-vision is required: `vision: {prepared, reference}`, the reference a number
or the JSON `--vision-reference-only` writes from the base model; `vision: off` only for discarded timing runs.

Paths: view and vision.prepared under data_root ($OBJDEL_DATA_ROOT, else the Vast bundle dir, else data); output_dir
and a vision.reference JSON under run_dir ($OBJDEL_RUN_DIR, else runs).

Data: view_rows in view-epoch order, length-grouped within each epoch (seed data_seed + e), so "b epochs" is the first
b passes over fresh views; the view must hold max(b) + its decay epochs (checked before anything loads). manifest.json
holds the config, view hash, code SHA-256, step plan, G-vision log, pid, status, `timing` and `peak_memory_gib` per
phase. `--resume` continues a crashed run (status "running", same config and view) from the newest complete rolling
checkpoint, then each unsaved branch. Every model is LoRA over a frozen bf16 base.

`deterministic: true` sets CUBLAS_WORKSPACE_CONFIG=:4096:8 unless set and torch.use_deterministic_algorithms(True,
warn_only=True) before anything loads. `branches_from: <run_dir-relative g_probe.json>` takes branches = [the depth of
the checkpoint with the final G-probe verdict] (Refusal when there is none). `profile: <dir>` (timing runs only)
writes torch.profiler traces of optimizer step WARM + 2 to <dir>/trace.json.gz and <dir>/ops.txt.

References: transformers get_wsd_schedule and Trainer resume
(https://github.com/huggingface/transformers/blob/v5.20.1/src/transformers/trainer.py); WSD: Hu et al., MiniCPM
(2404.06395); Hagele et al. (2405.18392).
"""

import gc
import math
import os
import re
import shutil
import socket
import time
from dataclasses import dataclass
from pathlib import Path

from objdel.core.hashing import file_sha256
from objdel.core.io import read_json, replace_dir, write_json
from objdel.core.provenance import code_hashes
from objdel.policy.errors import Refusal

ABORTED = 3  # exit code of a G-vision abort


WARM = 3  # timing: steps left out of the mean s/step (compilation, allocator and loader warm-up)


BOX_DATA, BOX_RUNS = Path("/workspace/data/lidc"), Path("/workspace/runs")  # sync_up / pull.sh on a Vast box (M7)


SAME_RUN = ("model", "branches", "a0", "seed", "steps_per_epoch", "sft")  # resume check


TRAINER_FILES = ("optimizer", "scheduler", "rng_state", "trainer_state", "adapter_")


def plan(args, steps_per_epoch: int) -> dict:
    """The step plan in optimizer steps: the stable run's length, and per branch its point and decay (as
    configs.sft_config computes them: warmup = 3 % of the stable steps, decay = 15 % of the branch point)."""
    from objdel.sft.config import DECAY, WARMUP

    at = {f"e{b:g}": max(1, round(b * steps_per_epoch)) for b in args.branches}
    br = {n: {"at": s, "decay": round(DECAY * s)} for n, s in at.items()}
    if len(set(at.values())) < len(at) or any(b["decay"] < 1 for b in br.values()):
        raise ValueError(f"branches {br} at {steps_per_epoch} steps/epoch: points collide or a decay is 0 steps")
    stable = max(at.values(), default=steps_per_epoch)
    return {"steps_per_epoch": steps_per_epoch, "stable": stable, "warmup": round(WARMUP * stable), "branches": br}


def peak_gib() -> float | None:
    """Peak CUDA memory allocated since the last reset_peak_memory_stats (each trainer resets it), None on CPU."""
    import torch

    return torch.cuda.max_memory_allocated() / 2**30 if torch.cuda.is_available() else None


def timing(stamps: list[float], batch: int) -> dict:
    """Step-end perf_counter stamps of one trainer -> mean s/step after the first WARM steps and peak CUDA memory."""
    s = (stamps[-1] - stamps[WARM - 1]) / (len(stamps) - WARM) if len(stamps) > WARM else None
    return {"steps": len(stamps), "warmup_excluded": WARM, "s_per_step": s, "samples_per_s": s and batch / s,
            "batch": batch, "max_memory_gib": peak_gib()}  # fmt: skip


def probed_depth(path: Path) -> float:
    """The depth (epochs) of the checkpoint that received g_probe's final verdict: its "ckpt" is <run>/branches/e<b>."""
    doc = read_json(path) if path.is_file() else {}
    if doc.get("verdict") and (m := re.fullmatch(r"e(\d+(\.\d+)?)", Path(doc.get("ckpt") or "").name)):
        return float(m[1])
    why = doc.get("skipped") or "missing or unfinished R2"
    raise Refusal(f"branches_from {path}: no G-probe verdict on a branch checkpoint ({why}): C2 must not start")


def resolve(args):
    """A copy of the SftStage with every path resolved: data_root ($OBJDEL_DATA_ROOT, else the Vast bundle dir, else
    data) and run_dir ($OBJDEL_RUN_DIR, else runs) when unset; view and vision.prepared under data_root; output_dir and
    a vision.reference JSON under run_dir. The input stage is unchanged."""
    data = args.data_root or os.environ.get("OBJDEL_DATA_ROOT") or (str(BOX_DATA) if BOX_DATA.is_dir() else "data")
    runs = args.run_dir if args.run_dir is not None else os.environ.get("OBJDEL_RUN_DIR", "runs")
    new = {"data_root": data, "run_dir": runs, "view": str(Path(data, args.view)),
           "output_dir": str(Path(runs, args.output_dir))}
    if args.vision:
        ref = args.vision["reference"]
        new["vision"] = {"prepared": str(Path(data, args.vision["prepared"])),
                         "reference": str(Path(runs, ref)) if isinstance(ref, str) else ref}
    return args.model_copy(update=new)


def reference(ref, model: str | None, set_sha256: str) -> float:
    """The G-vision base reference: a number, or R0's JSON for this model (when given) and fixed set."""
    if not isinstance(ref, str):
        return float(ref)
    doc = read_json(ref)
    if (doc["model"], doc["set_sha256"]) != (model or doc["model"], set_sha256):
        raise ValueError(f"{ref}: reference of {doc['model']} on set {doc['set_sha256'][:12]}, not this run's")
    return doc["reference"]


def load_items(cfg: dict, model: str | None, dataset: str) -> dict:
    """{"prepared", "reference"} and the dataset -> the fixed set, its checked reference and the dataset's prompt
    {"items", "set_sha256", "reference", "prompt"}."""
    from objdel.core import datasets
    from objdel.policy import vision

    items, sha = vision.fixed_set(cfg["prepared"], dataset)
    return {"items": items, "set_sha256": sha, "reference": reference(cfg["reference"], model, sha),
            "prompt": datasets.get(dataset).prompt}


def vision_reference(args) -> int:
    """projector_distance of the base model on the fixed set, CPU float32 -> vision.reference."""
    import torch

    from objdel.core import datasets
    from objdel.policy import models as policy_models
    from objdel.policy import vision

    if not (args.vision and isinstance(args.vision["reference"], str)):
        raise ValueError("--vision-reference-only writes vision.reference: set it to a .json path")
    src = args.model_path or policy_models.get(args.model).hf_id
    model, proc = policy_models.load_pair(src, dtype=torch.float32, proc_kw=policy_models.processor_kwargs(args.model))
    dataset = read_json(Path(args.view, "manifest.json"))["dataset"]  # the training view's
    items, sha = vision.fixed_set(args.vision["prepared"], dataset)
    d = vision.projector_distance(model, proc, [it["image"] for it in items[: vision.N_IMAGES]],
                                  datasets.get(dataset).prompt)  # inference_mode
    dest = Path(args.vision["reference"])
    dest.parent.mkdir(parents=True, exist_ok=True)
    write_json(dest, {"model": args.model, "source": src, "reference": d, "set_sha256": sha, "device": "cpu",
                      "dtype": "float32", "code": code_hashes()})
    return 0


def vision_check(model, processor, cfg: dict) -> dict:
    """objdel.policy.vision.check (inference mode, module modes restored) on the fixed set: cfg = load_items' output, or
    {"prepared", "reference", "model"?, "dataset"} (select_ckpt)."""
    from objdel.policy import vision

    cfg = cfg if "items" in cfg else load_items(cfg, cfg.get("model"), cfg["dataset"])
    images = [it["image"] for it in cfg["items"][: vision.N_IMAGES]]
    res = vision.check(model, processor, images, cfg["items"], reference=cfg["reference"], prompt=cfg["prompt"])
    return res | {"set_sha256": cfg["set_sha256"]}


def record_vision(man: dict, step: int, res: dict) -> bool:
    """Log one G-vision result in the manifest; True (stop training) on "abort": ratio and TV both fail."""
    keys = ("action", "ratio", "tv", "tv_black", "tv_other", "set_sha256", "device", "dtype")
    man["vision"].append({"step": step} | {k: res.get(k) for k in keys})
    if res["action"] == "abort":
        man["abort"] = {"step": step, "reason": f"G-vision (a) ratio {res['ratio']:.4g} < 0.5 and (b) TV "
                                                f"{res['tv']:.4g} < 0.002"}
    return res["action"] == "abort"


def latest(d) -> str | None:
    """The newest complete Trainer checkpoint in d (trainer_state.json is written last)."""
    done = [c for c in Path(d).glob("checkpoint-*") if (c / "trainer_state.json").exists()]
    return str(max(done, key=lambda c: int(c.name.rsplit("-", 1)[1]))) if done else None


@dataclass
class Base:
    """What every trainer of a run loads: the weights (src), their processor and the load dtype (the bf16 LoRA base)."""
    src: str
    proc: object
    dtype: object


def validate(args) -> dict:
    """The refusals before anything loads (profile outside a timing run, vision {}, outputs outside /workspace/runs on
    a Vast box, a view shorter than the deepest decay: S1-8) -> the view's manifest."""
    from objdel.sft import config as sft_config

    if args.profile and args.branches:
        raise ValueError("profile is for timing runs: branches: []")
    if args.vision == {}:
        raise ValueError("vision: {} skips G-vision: give {prepared, reference}; `vision: off` for timing")
    out = Path(args.output_dir)
    if BOX_RUNS.is_dir() and BOX_RUNS.resolve() not in out.resolve().parents:
        raise ValueError(f"{out}: on a Vast box outputs go under {BOX_RUNS}/<RUN_ID> (pull.sh): set OBJDEL_RUN_DIR")
    view = read_json(Path(args.view, "manifest.json"))
    deep = max((b * (1 + sft_config.DECAY) for b in args.branches), default=1.0)
    if view["config"]["epochs"] < math.ceil(deep):
        raise ValueError(f"{args.view}: {view['config']['epochs']} view epoch(s), the deepest branch decays to "
                         f"{deep:g} epochs: rebuild it with build_views --epochs {math.ceil(deep)}")
    return view


def resume_manifest(args, out: Path, view: dict) -> dict:
    """The run manifest (config, code, view identity, pid, status "running"). A run already in out is refused unless
    --resume, which needs a crashed run ("running") of the same config (SAME_RUN) and view, and carries its G-vision
    log, logs, start, timing and peak memory over."""
    old = read_json(out / "manifest.json") if (out / "manifest.json").exists() else None
    if old and not args.resume:
        raise ValueError(f"{out} holds a run: --resume continues it, or pick another output_dir")
    man = {"config": args.model_dump(mode="json", exclude={"kind", "name", "soft"}), "code": code_hashes(),
           "vision": [], "log": {},
           "view": {"dir": args.view, "view_hash": view.get("view_hash"),
                    "manifest_sha256": file_sha256(Path(args.view, "manifest.json"))}}
    if args.resume:
        same = old and all(old["config"][k] == man["config"][k] for k in SAME_RUN)
        if not (same and old["status"] == "running" and old["view"]["view_hash"] == man["view"]["view_hash"]):
            raise ValueError(f"--resume: {out} is not a crashed run of this config and view")
        man |= {k: old[k] for k in ("vision", "log", "started", "timing", "peak_memory_gib") if k in old}
        man["resumed"] = old.get("resumed", []) + [time.time()]
    return man | {"pid": os.getpid(), "host": socket.gethostname(), "status": "running",
                  "started": man.get("started", time.time()), "peak_memory_gib": man.get("peak_memory_gib", {})}


class Trainer:
    """train's overrides of trl's SFTTrainer, mixed in by trainer() with epochs and lengths (one per row) as class
    attributes (trl and transformers load late: ~4 s of imports, this module ~0.03 s)."""

    def floating_point_ops(self, inputs):
        """0: perf F3, no per-micro-batch walk over ~2.5k PEFT modules."""
        return 0

    def _get_train_sampler(self, train_dataset=None):  # view epochs in order, length-grouped inside each
        import torch
        from transformers.trainer_pt_utils import get_length_grouped_indices

        bs = self.args.train_batch_size * self.args.gradient_accumulation_steps
        order = []
        for e in sorted(set(self.epochs)):
            idx = [i for i, x in enumerate(self.epochs) if x == e]
            g = torch.Generator().manual_seed(self.args.data_seed + e)
            order += [idx[j] for j in get_length_grouped_indices([self.lengths[i] for i in idx], bs, generator=g)]
        return order


class Rolling:
    """Step-end stamps, a save every `minutes` and at each step of keep_at ({step: branch name}, hard-linked to
    out/stable/<name>), G-vision after every save (vis: load_items' output, or False) logged in man; stops at stop_at.
    train() mixes in transformers' TrainerCallback."""

    def __init__(self, man: dict, out: Path, proc, vis, minutes: float, keep_at=(), stop_at=None):
        self.man, self.out, self.proc, self.vis, self.minutes = man, out, proc, vis, minutes
        self.keep_at, self.stop_at, self.last, self.stamps = dict(keep_at), stop_at, time.monotonic(), []

    def on_step_end(self, a, state, control, **kw):
        """Stamp the step; save when due or at a keep_at step; stop at stop_at."""
        self.stamps.append(time.perf_counter())
        due = time.monotonic() - self.last >= 60 * self.minutes
        control.should_save = due or state.global_step in self.keep_at
        control.should_training_stop |= state.global_step == self.stop_at

    def on_save(self, a, state, control, model=None, **kw):
        """Hard-link a keep_at checkpoint to stable/<name>, then run G-vision and log it in the manifest."""
        self.last, step = time.monotonic(), state.global_step
        if step in self.keep_at:
            ckpt = Path(a.output_dir, f"checkpoint-{step}")
            replace_dir(lambda tmp: shutil.copytree(ckpt, tmp, copy_function=os.link),  # hard links: no copy
                        self.out / "stable" / self.keep_at[step])
        if self.vis:  # vision.check restores the module modes; training_step calls model.train() anyway
            prm = next(model.parameters())
            res = vision_check(model, self.proc, self.vis) | {"device": str(prm.device), "dtype": str(prm.dtype)}
            control.should_training_stop |= record_vision(self.man, step, res)
            write_json(self.out / "manifest.json", self.man)


def trainer(args, base: Base, rows, stable: int, out_dir, rolling: Rolling, branch_at=None):
    """One trainer from base: SFTConfig (stable steps; branch_at = a branch's decay), the LM-only LoRA, callbacks
    rolling and the profiler."""
    import torch
    from transformers import AutoModelForImageTextToText
    from trl import SFTTrainer

    from objdel.policy import callbacks
    from objdel.policy import models as policy_models
    from objdel.sft import config as sft_config

    if torch.cuda.is_available():  # this phase's peak (Qwen3-VL LoRA R1: 6.1 GiB)
        torch.cuda.reset_peak_memory_stats()
    model = AutoModelForImageTextToText.from_pretrained(base.src, dtype=base.dtype)
    o = {"seed": args.seed, "data_seed": args.seed, "save_strategy": "no", "save_total_limit": 3}
    o |= {"logging_steps": 10}
    cfg = sft_config.sft_config(args.model, str(out_dir), stable, branch_at, **(o | args.sft))
    kw = {"train_dataset": rows, "processing_class": base.proc, "peft_config": policy_models.lora(args.model)}
    calls = [rolling] + ([callbacks.TraceStep(args.profile, WARM)] if args.profile else [])
    cls = type("Trainer", (Trainer, SFTTrainer), {"epochs": list(rows["epoch"]), "lengths": list(rows["length"])})
    return cls(model=model, args=cfg, callbacks=calls, **kw)


def finish(base: Base, d: Path) -> None:
    """A kept stable checkpoint, model-only (LoRA merged) with the processor, for the depth curve."""
    if (d / "adapter_config.json").exists():
        from objdel.policy import models as policy_models

        policy_models.merged(base.src, d, base.dtype).save_pretrained(d)
    for f in d.glob("*"):
        if f.name.startswith(TRAINER_FILES) or f.suffix == ".pt":
            f.unlink()
    base.proc.save_pretrained(d)


def train(args) -> int:
    """One SFT run of an objdel.core.spec.SftStage."""
    os.environ.setdefault("PYTORCH_ALLOC_CONF", "expandable_segments:True")  # read at the first CUDA allocation
    import torch
    from transformers import TrainerCallback

    from objdel.policy import models as policy_models
    from objdel.sft import rows as dataset

    args = resolve(args)
    if args.vision_reference_only:
        return vision_reference(args)
    if args.branches_from:
        args = args.model_copy(update={"branches": (probed_depth(Path(args.run_dir, args.branches_from)),)})
    if args.deterministic:  # before CUDA and cuBLAS initialise
        os.environ.setdefault("CUBLAS_WORKSPACE_CONFIG", ":4096:8")
        torch.use_deterministic_algorithms(True, warn_only=True)
    view, out = validate(args), Path(args.output_dir)
    vis = args.vision and load_items(args.vision, args.model, view["dataset"])
    out.mkdir(parents=True, exist_ok=True)
    man = resume_manifest(args, out, view)
    write_json(out / "manifest.json", man)
    src = args.model_path or policy_models.get(args.model).hf_id
    proc = policy_models.processor(args.model, src)
    base = Base(src, proc, torch.bfloat16)  # the frozen LoRA base
    rows = dataset.view_rows(args.view, proc.tokenizer, policy_models.get(args.model).max_completion, a0=args.a0)
    epochs = list(rows["epoch"])
    gbs = args.sft.get("per_device_train_batch_size", 8) * args.sft.get("gradient_accumulation_steps", 1)  # 1 GPU
    rows0 = epochs.count(min(epochs, default=0))  # rows of one view epoch
    p = man["plan"] = plan(args, args.steps_per_epoch or rows0 // gbs)
    need = max([b["at"] + b["decay"] for b in p["branches"].values()] or [p["stable"]]) * gbs
    if need > len(rows):
        raise ValueError(f"{args.view}: {len(rows)} rows, the deepest branch and its decay need {need}")
    last = latest(out / "rolling")
    for name, b in p["branches"].items():  # a crash between a branch point's rolling save and its stable copy
        if last and int(last.rsplit("-", 1)[1]) >= b["at"] and not (out / "stable" / name).is_dir():
            raise ValueError(f"--resume: {last} is past branch {name} (step {b['at']}) but {out}/stable/{name} is "
                             f"missing: copy {out}/rolling/checkpoint-{b['at']} there if it exists, else start over")
    callback = type("Rolling", (Rolling, TrainerCallback), {})

    keep = {b["at"]: name for name, b in p["branches"].items()}
    if not (keep and (out / "stable" / keep[p["stable"]]).is_dir()):  # stable phase not finished yet
        rolling = callback(man, out, proc, vis, args.rolling_minutes, keep, p["stable"])
        t = trainer(args, base, rows, p["stable"], out / "rolling", rolling, p["stable"] if keep else None)
        t.train(resume_from_checkpoint=latest(out / "rolling") if args.resume else None)
        man["log"]["stable"], man["peak_memory_gib"]["stable"] = t.state.log_history, peak_gib()
        man.setdefault("timing", timing(rolling.stamps, gbs) | {"rows_per_epoch": rows0})
        write_json(out / "manifest.json", man)
    if all((out / "stable" / n).is_dir() for n in keep.values()):  # rolling/ holds full AdamW checkpoints
        shutil.rmtree(out / "rolling", ignore_errors=True)
    for name, b in p["branches"].items():
        if "abort" in man:
            break
        if not (out / "branches" / name).is_dir():
            t = None  # V1 R1: the previous trainer's model, grads and AdamW state (~7 GB) made e1 OOM
            gc.collect()
            torch.cuda.empty_cache()
            rolling = callback(man, out, proc, vis, args.rolling_minutes)
            t = trainer(args, base, rows, p["stable"], out / "decay" / name, rolling, branch_at=b["at"])
            t.train(resume_from_checkpoint=(args.resume and latest(out / "decay" / name)) or str(out / "stable" / name))
            man["log"][name], man["peak_memory_gib"][name] = t.state.log_history, peak_gib()
            if "abort" in man:  # G-vision stopped the decay: the branch is unfinished
                break
            policy_models.save_merged(t.model, proc, out / "branches" / name)  # D6: merged
            shutil.rmtree(out / "decay" / name, ignore_errors=True)
        finish(base, out / "stable" / name)  # frees the optimizer state once its branch exists
        write_json(out / "manifest.json", man)
    for d in (out / "stable").glob("e*"):  # after an abort too
        finish(base, d)
    for d in ("rolling", "decay"):
        shutil.rmtree(out / d, ignore_errors=True)
    man |= {"status": "aborted" if "abort" in man else "done", "finished": time.time()}
    write_json(out / "manifest.json", man)
    return ABORTED if "abort" in man else 0

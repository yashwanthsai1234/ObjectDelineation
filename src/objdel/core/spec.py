"""One typed run spec: the stages of a run (what runs, in order, with which flags) as a strict pydantic model.

`load_spec(path)` reads a YAML file (duplicate keys and unknown keys are errors; scalars are not coerced, so write
`5.0e-5`, not `5e-5`), `to_argv(stage, rd, data)` is the one definition of every CLI flag, and `cli_parse(Stage, argv)`
reads the same flags back into a stage. Rules that need other objdel packages (a model name, the reward group size)
do not live here.

{rd} and {data} in a string field are the run folder and the data folder on the box. A stage starts with `~` in the
argv list (`soft`) when its failure must not stop the chain.
"""

import argparse
import json
import math
import re
import string
import sys
from collections.abc import Sequence
from pathlib import Path
from typing import Annotated, Literal

import yaml
from pydantic import BaseModel, BeforeValidator, Field, JsonValue, ValidationError, model_validator
from pydantic_settings import CliApp, CliSettingsSource, CliSuppress

HEX40, HEX64 = re.compile(r"[0-9a-f]{40}"), re.compile(r"[0-9a-f]{64}")
PLACEHOLDERS = {"rd", "data"}


def _tuple(v):
    return tuple(v) if isinstance(v, list) else v


TUPLE_FIELDS = {"view", "split", "out", "argv", "branches", "profile_window"}  # YAML lists; strict mode wants tuples


class Base(BaseModel, frozen=True, extra="forbid"):
    """Strict, frozen model: unknown keys are errors."""


class StageBase(Base):
    """`name` and `soft` belong to the spec only: a CLI gives name="" and soft=False."""

    name: CliSuppress[str] = ""
    soft: CliSuppress[bool] = False

    @model_validator(mode="before")
    @classmethod
    def _lists_are_tuples(cls, d):
        return {k: _tuple(v) if k in TUPLE_FIELDS else v for k, v in d.items()} if isinstance(d, dict) else d


class TrainerFlags(BaseModel, frozen=True, extra="forbid"):
    """Flags of the trainer that train_rl and g_probe build alike (rl.build.build_trainer reads each one)."""

    lr: float | None = None  # None: configs.rl_lr(model)
    shared_prefix: bool = False
    compile: Literal["default", "reduce-overhead"] | None = None
    plain: bool = False  # TRL's GRPOTrainer as shipped, not FastGRPO (profile arm)
    trace: str | None = None
    grpo: dict[str, JsonValue] = {}  # GRPOConfig overrides
    log_by_prompt: bool = False
    fast_sync: bool = False
    group_requests: bool = False
    no_ckpt_layers: int = Field(0, ge=0)  # C9 arm (a): <= the model's text decoder layers (rl.run.checked); no clamp
    runaway_cap: float = Field(0, ge=0)
    profile_window: tuple[int, int] | None = None  # S,E: NvtxSteps; 1,steps traces every step; S >= 2 with trace
    count_syncs: bool = False
    host_meta: bool = False
    check_meta_every: int = Field(50, ge=0)
    check_rope_every: int = Field(50, ge=0)
    checks_from_step: int = Field(1, ge=1)  # C0: the check flags (and train_rl's record_batches) act from step S on
    prompt_memo: bool = False
    encoder_once: bool = False  # C3, with host_meta
    gpu_pixels: Literal["off", "auto"] = "off"
    allow_cpu_pixels: bool = False

    @model_validator(mode="after")
    def _window(self):
        w = self.profile_window
        if w and not (1 <= w[0] <= w[1] and (w[0] >= 2 or not self.trace)):
            raise ValueError(f"profile_window {w}: needs 1 <= S <= E, and S >= 2 with --trace (torch.profiler warms "
                             "up one step before step S)")
        return self


Money = Annotated[float, Field(allow_inf_nan=False)]


SAVE_STEPS = 50  # train_rl CLI default; E2.7 refuses any other save_steps above steps


class RlStage(StageBase, TrainerFlags, frozen=True, extra="forbid"):  # train_rl
    """A train_rl stage."""

    kind: CliSuppress[Literal["rl"]] = "rl"
    merge_last: bool = False
    model: str
    ckpt: str
    train: str | None = None
    out: str
    steps: int | None = Field(None, ge=1)  # TRL reads max_steps <= 0 as "train num_train_epochs"
    save_steps: int = Field(SAVE_STEPS, ge=1)  # a short smoke run on the default writes none
    resume: bool = False
    rows_revision: str | None = None
    rows_manifest_sha256: str | None = None
    budget_usd: Money | None = Field(None, gt=0)
    session_start: Money | None = Field(None, ge=0)  # None: OBJDEL_SESSION_START at run time (train_rl)
    rate_usd_per_h: Money | None = Field(None, gt=0)
    fixed_usd: Money | None = Field(None, ge=0)
    heldout: str | None = None
    heldout_items: int = Field(64, ge=1)
    heldout_every: int = Field(50, ge=1)
    reward: Literal["v4", "v5"] = "v4"
    numerics: str | None = None  # a configs.NUMERICS key; None: the B1 control (train_rl checks the name)
    f_check: bool = False  # B2.2 F step 5; neither flag: on under f_check.PRESETS (train_rl resolves it)
    no_f_check: bool = False
    b1_adapter: str | None = None  # the B1 final adapter of the F comparator (--f-check under lora-seqmask)
    opsm_start_delta: float | None = None
    record_batches: str | None = None  # C0: each generation batch from checks_from_step on, to DIR/step<N>.pt
    expect_step1: str | None = None  # C0 (E-EQ6): another arm's rl.json; needs encode_check
    encode_check: str | None = None  # C0: c_speed_check --encoder's encoder.json
    dataset: str = "lidc"  # a key of objdel.core.datasets; its data folders must hold it (policy.rows.check_root)

    @model_validator(mode="after")
    def _rules(self):
        if not self.merge_last and (self.train is None or self.steps is None):
            raise ValueError("train and steps are required (except with merge_last)")
        if self.steps is not None and self.save_steps > self.steps and self.save_steps != SAVE_STEPS:
            raise ValueError(f"save_steps {self.save_steps} > steps {self.steps} writes no checkpoint (E2.7)")
        if self.f_check and self.no_f_check:
            raise ValueError("f_check and no_f_check exclude each other")
        if self.expect_step1 and not self.encode_check:
            raise ValueError("expect_step1 needs encode_check (its bits_equal picks the grad_norm rule)")
        if (self.rows_revision is None) != (self.rows_manifest_sha256 is None):
            raise ValueError("rows_revision and rows_manifest_sha256 are both set or both unset")
        if self.rows_revision is not None and not (HEX40.fullmatch(self.rows_revision)
                                                   and HEX64.fullmatch(self.rows_manifest_sha256)):
            raise ValueError("rows_revision needs 40 and rows_manifest_sha256 64 lowercase hex digits")
        return self


class SftStage(StageBase, frozen=True, extra="forbid"):  # train_sft, minus the runner-only stages, sampler, rows_plan
    """A train_sft stage."""

    kind: CliSuppress[Literal["sft"]] = "sft"
    model: str
    view: str
    output_dir: str
    data_root: str | None = None  # None: the CLI's environment-dependent default
    run_dir: str | None = None
    model_path: str | None = None
    branches: tuple[float, ...] = (0.5, 1.0, 2.0)
    a0: bool = False
    seed: int = 1
    rolling_minutes: float = 10.0
    vision: dict[str, JsonValue] | Literal[False] = {}
    vision_reference_only: bool = False
    resume: bool = False
    steps_per_epoch: int | None = None
    sft: dict[str, JsonValue] = {}
    deterministic: bool = False
    branches_from: str | None = None
    profile: str | None = None

    @model_validator(mode="before")
    @classmethod
    def _vision(cls, d):
        if isinstance(d, dict) and d.get("vision") == "off":
            d = d | {"vision": False}
        return d

    @model_validator(mode="after")
    def _rules(self):
        if self.vision:  # the check RunArgs makes in __post_init__
            ref = self.vision.get("reference")
            num = isinstance(ref, int | float) and not isinstance(ref, bool) and math.isfinite(ref) and ref > 0
            if set(self.vision) != {"prepared", "reference"} or not (num or str(ref).endswith(".json")):
                raise ValueError(f"vision {self.vision!r}: need {{prepared, reference}}, reference = R0's base "
                                 "projector distance (a positive number, or a .json file)")
        return self


class VisionRef(Base):
    """A vision check: the prepared file and R0's reference projector distance (number or .json file)."""

    prepared: str
    reference: str | float


class SelectStage(StageBase, frozen=True, extra="forbid"):  # select_ckpt
    """A select_ckpt stage."""

    kind: CliSuppress[Literal["select"]] = "select"
    run: str
    prepared: str | None = None  # the image-blind prior's training rows (not read under only)
    probes: str | None = None
    sampler: Literal["vllm", "hf"] | None = None
    shifted: str | None = None
    rescaled: str | None = None
    rotated: str | None = None
    items: int | None = None
    vision: VisionRef | None = None
    no_vision: bool = False
    only: Literal["topology"] | None = None
    model: str | None = None
    k: int = 8
    seed: int = 0
    timeout_s: float | None = None
    max_tokens: int | None = Field(None, ge=1)  # sampled completion cap (CPU smoke); None: the model's max_completion
    out: str | None = None
    dataset: str = "lidc"  # a key of objdel.core.datasets; its data folders must hold it (policy.rows.check_root)

    @model_validator(mode="after")
    def _rules(self):
        if self.only is None and (self.vision is None) == (not self.no_vision):
            raise ValueError("exactly one of vision and no_vision (unless only)")
        if self.only is not None and (self.vision is not None or self.no_vision):
            raise ValueError("with only set, neither vision nor no_vision")
        if self.only is None and self.shifted is None:
            raise ValueError("shifted is required unless only")
        return self


class ProbeStage(StageBase, TrainerFlags, frozen=True, extra="forbid"):  # g_probe
    """A g_probe stage."""

    kind: CliSuppress[Literal["probe"]] = "probe"
    model: str
    ckpt: str | None = None
    select: str | None = None
    train: str
    probe: str
    out: str
    items: int | None = None
    dataset: str = "lidc"  # a key of objdel.core.datasets; its data folders must hold it (policy.rows.check_root)

    @model_validator(mode="after")
    def _rules(self):
        if (self.ckpt is None) == (self.select is None):
            raise ValueError("exactly one of ckpt and select")
        return self


class EvalStage(StageBase, frozen=True, extra="forbid"):  # final_eval
    """A final_eval stage."""

    kind: CliSuppress[Literal["eval"]] = "eval"
    ckpt: str | None = None
    view: tuple[str, ...] = Field(min_length=1)
    split: tuple[str, ...] = ()  # () -> every view uses "val"; else one split per view
    out: tuple[str, ...] = Field(min_length=1)  # one per view, or one
    model: str = "qwen3-vl-2b"
    adapter: str | None = None  # ckpt's base plus this LoRA folder (a lora-seqmask run), served unmerged
    k: int = Field(8, ge=2)
    seed: int = 11
    rescore: str | None = None
    limit: int | None = Field(None, ge=1)  # smoke runs: at most N items per view
    dataset: str = "lidc"  # a key of objdel.core.datasets; its data folders must hold it (policy.rows.check_root)

    @model_validator(mode="after")
    def _rules(self):
        if (self.ckpt is None) == (self.rescore is None):
            raise ValueError("exactly one of ckpt and rescore")
        if self.adapter is not None and self.ckpt is None:
            raise ValueError("adapter needs ckpt (its base)")
        if len(self.split) not in (0, len(self.view)) or len(self.out) not in (1, len(self.view)):
            raise ValueError("split: none or one per view; out: one, or one per view")
        if any("," in x for x in (*self.view, *self.out)):
            raise ValueError("a comma splits a repeated flag value: no comma in view or out")
        if not set(self.split) <= {"val", "test"}:
            raise ValueError("split is val or test")
        return self


class CmdStage(StageBase, frozen=True, extra="forbid"):  # escape hatch: a literal argv (env_check.sh, shell steps)
    """A literal argv stage."""

    kind: CliSuppress[Literal["cmd"]] = "cmd"
    argv: tuple[str, ...] = Field(min_length=1)


Stage = Annotated[RlStage | SftStage | SelectStage | ProbeStage | EvalStage | CmdStage, Field(discriminator="kind")]
STAGE_CLASSES = {c.model_fields["kind"].default: c for c in (RlStage, SftStage, SelectStage, ProbeStage, EvalStage,
                                                              CmdStage)}
# The package of the command-line modules (`python -I -m objdel.cli.<group>.<command>`, or `objdel <group> <command>`)
# and the module each stage kind runs; the box chain, preflight and the Vast CLI read both.
CLI = "objdel.cli"
MODULE = {k: f"{CLI}.{m}" for k, m in {"rl": "train.rl", "sft": "train.sft", "select": "train.select",
                                       "probe": "train.probe", "eval": "eval.run"}.items()}

# The order of the flags in today's plans; `views` stands for the interleaved view, split and out of an eval stage.
TRAINER_ORDER = ("fast_sync", "group_requests", "no_ckpt_layers", "runaway_cap", "profile_window", "trace",
                 "count_syncs", "host_meta", "check_meta_every", "check_rope_every", "checks_from_step",
                 "prompt_memo", "encoder_once", "gpu_pixels", "allow_cpu_pixels", "compile", "plain")
FLAG_ORDER = {
    "rl": ("merge_last", "model", "ckpt", "train", "out", "steps", "save_steps", "shared_prefix", "log_by_prompt",
           "lr", "resume", "reward", "numerics", "f_check", "no_f_check", "b1_adapter", "opsm_start_delta", "grpo",
           "rows_revision", "rows_manifest_sha256", "budget_usd", "session_start", "rate_usd_per_h", "fixed_usd",
           "heldout", "heldout_items", "heldout_every", *TRAINER_ORDER, "record_batches", "expect_step1",
           "encode_check", "dataset"),
    "sft": ("model", "view", "output_dir", "data_root", "run_dir", "model_path", "branches", "a0", "seed",
            "rolling_minutes", "vision", "vision_reference_only", "resume", "steps_per_epoch", "sft",
            "deterministic", "branches_from", "profile"),
    "select": ("run", "prepared", "probes", "sampler", "shifted", "rescaled", "rotated", "items", "vision",
               "no_vision", "only", "model", "k", "seed", "timeout_s", "max_tokens", "out", "dataset"),
    "probe": ("model", "ckpt", "select", "train", "probe", "out", "items", "grpo", "compile", "shared_prefix", "plain",
              "trace", "lr", "log_by_prompt", *(f for f in TRAINER_ORDER if f not in ("compile", "plain",
                                                                                               "trace")), "dataset"),
    "eval": ("ckpt", "views", "model", "adapter", "k", "seed", "rescore", "limit", "dataset"),
}


def spec_dataset(spec: "RunSpec") -> str:
    """The one dataset the spec's stages name (the stages with a dataset field; default lidc when none has one)."""
    found = {s.dataset for s in spec.stages if "dataset" in type(s).model_fields}
    if len(found) > 1:
        raise ValueError(f"{spec.name}: stages name {sorted(found)}; a run reads one dataset")
    return found.pop() if found else RlStage.model_fields["dataset"].default


class ComputeSpec(Base):
    """The rented GPU, the hour cap and the run budget (inherited by the rl stages that train) of a run."""

    gpu: str
    cap_hours: float = Field(gt=0)
    budget_usd: Money | None = Field(None, gt=0)


class RowsRef(Base):
    """A rows root on the Hub: its folder name under {data}/rows, its revision and its manifest's SHA-256."""

    name: str
    revision: str
    manifest_sha256: str


def _inherit(data: dict) -> dict:
    """The run-level model, dataset, rows and compute.budget_usd copied into each stage (a dict) whose class has the
    field and that does not set it, when the value is not the field's default (to_argv emits only set fields, so a
    spec written either way gives the same argv). The budget goes to the rl stages that train only; rows (train
    and both pins, as one unit) to those that also leave train unset."""
    rows = data.get("rows")
    top = {"model": data.get("model"), "dataset": data.get("dataset"),
           "budget_usd": (data.get("compute") or {}).get("budget_usd")}
    stages = []
    for st in data.get("stages") or ():
        cls = STAGE_CLASSES.get(st.get("kind")) if isinstance(st, dict) else None
        trains = cls is RlStage and not st.get("merge_last")
        add = {} if cls is None else {
            k: v for k, v in top.items() if k in cls.model_fields and k not in st and v is not None
            and v != cls.model_fields[k].default and (k != "budget_usd" or trains)}
        if trains and isinstance(rows, dict) and "train" not in st:
            add |= {k: v for k, v in {"train": f"{{data}}/rows/{rows.get('name')}", "rows_revision": rows.get(
                "revision"), "rows_manifest_sha256": rows.get("manifest_sha256")}.items() if k not in st}
        stages.append(st | add if add else st)
    return data | {"stages": stages}


class RunSpec(Base):
    """One run: its compute, the settings its stages inherit (model, dataset, rows; `_inherit`) and its ordered
    stages."""

    schema_version: Literal[1]
    name: str
    model: str | None = None
    dataset: str | None = None
    rows: RowsRef | None = None
    compute: ComputeSpec
    stages: Annotated[tuple[Stage, ...], BeforeValidator(_tuple), Field(min_length=1)]
    bundle_view: str | None = None  # the SFT training view under the data root, read by the preflight DATA checks

    @model_validator(mode="before")
    @classmethod
    def _defaults(cls, data):
        return _inherit(data) if isinstance(data, dict) else data

    @model_validator(mode="after")
    def _rules(self):
        names = [s.name for s in self.stages]
        if "" in names or len(set(names)) != len(names):
            raise ValueError(f"stage names are non-empty and unique: {names}")
        kinds = {s.kind for s in self.stages}
        if "rl" in kinds and "sft" not in kinds and self.bundle_view is None:
            raise ValueError("a spec with an rl stage and no sft stage needs bundle_view")
        if "sft" in kinds and self.bundle_view is not None:
            raise ValueError("bundle_view beside an sft stage: two training views")
        return self


class _UniqueKeyLoader(yaml.SafeLoader):
    pass


def _construct_mapping(loader, node, deep=False):
    seen = set()
    for key_node, _ in node.value:
        if key_node.tag == "tag:yaml.org,2002:merge":
            continue
        key = loader.construct_object(key_node, deep=True)
        if key in seen:
            raise yaml.constructor.ConstructorError(None, None, f"duplicate key {key!r}", key_node.start_mark)
        seen.add(key)
    return yaml.SafeLoader.construct_mapping(loader, node, deep)


_UniqueKeyLoader.add_constructor(yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG, _construct_mapping)


def load_spec(path: Path) -> RunSpec:
    """The spec in a YAML file. Raises on a duplicate or unknown key, a wrong type, a broken rule or an unknown
    placeholder."""
    spec = RunSpec.model_validate(yaml.load(Path(path).read_text(), Loader=_UniqueKeyLoader), strict=True)
    for s in spec.stages:
        to_argv(s, "rd", "data")
    return spec


def num(v: float) -> str:
    """A float as the text plans write: 5e-5, not 5e-05."""
    return re.sub(r"e([+-])0*(\d)", r"e\1\2", repr(v)).replace("e+", "e")


def fill(s: str, rd: str, data: str) -> str:
    """str.format with {rd} and {data}; {{ and }} are literal braces, any other placeholder raises."""
    unknown = {f for _, f, _, _ in string.Formatter().parse(s) if f is not None} - PLACEHOLDERS
    if unknown:
        raise ValueError(f"{s!r}: unknown placeholder {sorted(unknown)}; only {sorted(PLACEHOLDERS)}")
    return s.format(rd=rd, data=data)


def _flag(name: str) -> str:
    return f"-{name}" if len(name) == 1 else "--" + name.replace("_", "-")


def _value(v, rd: str, data: str) -> list[str]:
    if isinstance(v, str):
        return [fill(v, rd, data)]
    if isinstance(v, float):
        return [num(v)]
    if isinstance(v, VisionRef):  # its paths are the only filled part of a JSON value
        ref = fill(v.reference, rd, data) if isinstance(v.reference, str) else v.reference
        return [json.dumps({"prepared": fill(v.prepared, rd, data), "reference": ref})]  # as chain() wrote it
    if isinstance(v, tuple | dict):
        return [json.dumps(list(v) if isinstance(v, tuple) else v, separators=(",", ":"))]
    return [str(v)]


def to_argv(stage: Stage, rd: str, data: str) -> list[str]:
    """The argv of a stage on the box: `~` first when soft, then `-m <MODULE of its kind>` and the flags that
    the spec sets (in FLAG_ORDER), or a cmd stage's argv. {rd} and {data} are filled in string fields."""
    if stage.kind == "cmd":
        body = [fill(x, rd, data) for x in stage.argv]
    else:
        body = ["-m", MODULE[stage.kind]]
        for name in FLAG_ORDER[stage.kind]:
            v = None if name == "views" else getattr(stage, name)
            if name == "views":
                body += _views(stage, rd, data)
            elif name not in stage.model_fields_set or v is None or v == ():
                continue
            elif name == "vision" and v is False:
                body += ["--vision", "off"]
            elif v is True:
                body.append(_flag(name))
            elif v is not False:
                body += [_flag(name), *_value(v, rd, data)]
    return (["~"] if stage.soft else []) + body


def _views(s: EvalStage, rd: str, data: str) -> list[str]:
    """--view/--split/--out in groups, then one shared --out after the last group."""
    out = []
    for i, v in enumerate(s.view):
        out += ["--view", fill(v, rd, data)]
        if s.split:
            out += ["--split", s.split[i]]
        if len(s.out) == len(s.view):
            out += ["--out", fill(s.out[i], rd, data)]
    if len(s.out) != len(s.view):
        out += ["--out", fill(s.out[0], rd, data)]
    return out


def cli_parse(cls: type[Stage], argv: Sequence[str]) -> Stage:
    """The stage a CLI reads from its flags, in the one parser form (kebab flags, bare boolean flags). "toggle", not
    True: True gives argparse's BooleanOptionalAction, which reads every --no-* option as False (--no-vision)."""
    return CliApp.run(cls, cli_args=list(argv),
                      cli_settings_source=CliSettingsSource(cls, cli_kebab_case=True, cli_implicit_flags="toggle"))


def same_stage(a: Stage, b: Stage) -> bool:
    """Equal but for `name` and `soft`, which belong to the spec only."""
    return a.model_dump(exclude={"name", "soft"}) == b.model_dump(exclude={"name", "soft"})


def cli_stage(cls: type[Stage], argv: Sequence[str]) -> Stage:
    """The stage a CLI runs. `--spec PATH --stage NAME [--rd RD] [--data DATA]` takes the stage NAME of a spec file
    with {rd} and {data} filled (default runs/<spec name> and data, the Mac folders); any other argv is the CLI's own
    flags (cli_parse). A broken rule exits 2 with the reason, as an argparse error does."""
    argv = list(argv)
    try:
        if argv[:1] == ["--spec"]:
            ap = argparse.ArgumentParser(prog=cls.__name__)
            ap.add_argument("--spec", type=Path, required=True)
            ap.add_argument("--stage", required=True)
            ap.add_argument("--rd")
            ap.add_argument("--data", default="data")
            a = ap.parse_args(argv)
            spec = load_spec(a.spec)
            stage = next((s for s in spec.stages if s.name == a.stage), None)
            if not isinstance(stage, cls):
                ap.error(f"{a.spec} has no {cls.model_fields['kind'].default} stage named {a.stage!r}")
            argv = to_argv(stage.model_copy(update={"soft": False}), a.rd or f"runs/{spec.name}", a.data)[2:]
        return cli_parse(cls, argv)
    except ValidationError as e:
        why = "; ".join(f"{_flag(str(x['loc'][0])) + ': ' if x['loc'] else ''}{x['msg']}" for x in e.errors())
        print(f"{cls.__name__}: error: {why}", file=sys.stderr)
        raise SystemExit(2) from e

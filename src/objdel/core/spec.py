"""One typed run spec: the stages of a run (what runs, in order, with which flags) as a strict pydantic model.

`load_spec(path)` reads a YAML file (duplicate keys and unknown keys are errors; scalars are not coerced, so write
`5.0e-5`, not `5e-5`), `to_argv(stage, rd, data)` is the one definition of every CLI flag, and `cli_parse(Stage, argv)`
reads the same flags back into a stage. Rules that need other objdel packages (a model name, the reward group size)
do not live here.

{rd} and {data} in a string field are the run folder and the data folder on the box. A stage starts with `~` in the
argv list (`soft`) when its failure must not stop the chain.
"""

import json
import re
import string
from collections.abc import Sequence
from pathlib import Path
from typing import Annotated, Literal

import yaml
from pydantic import BaseModel, BeforeValidator, Field, JsonValue, model_validator
from pydantic_settings import CliApp, CliSettingsSource, CliSuppress

HEX40, HEX64 = re.compile(r"[0-9a-f]{40}"), re.compile(r"[0-9a-f]{64}")
PLACEHOLDERS = {"rd", "data"}


def _tuple(v):
    return tuple(v) if isinstance(v, list) else v


TUPLE_FIELDS = {"view", "split", "out", "argv", "branches"}  # YAML lists; strict mode wants tuples


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


class TrainerFlags(BaseModel):
    """Flags that train_rl and g_probe share."""

    shared_prefix: bool = False
    grpo: dict[str, JsonValue] | None = None


class RlStage(StageBase, TrainerFlags, frozen=True, extra="forbid"):  # train_rl
    """A train_rl stage."""

    kind: CliSuppress[Literal["rl"]] = "rl"
    merge_last: bool = False
    model: str
    ckpt: str
    train: str | None = None
    out: str
    steps: int | None = None
    save_steps: int = 50
    resume: bool = False
    lr: float | None = None
    rows_revision: str | None = None
    rows_manifest_sha256: str | None = None
    budget_usd: float | None = None
    session_start: float | None = None
    rate_usd_per_h: float | None = None
    fixed_usd: float | None = None
    heldout: str | None = None
    heldout_items: int = 64
    heldout_every: int = 50

    @model_validator(mode="after")
    def _rules(self):
        if not self.merge_last and (self.train is None or self.steps is None):
            raise ValueError("train and steps are required (except with merge_last)")
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
    freeze_vision: bool = False
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
            num = isinstance(ref, int | float) and not isinstance(ref, bool) and ref > 0
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
    prepared: str
    probes: str
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
    out: str | None = None

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
    compile: Literal["default", "reduce-overhead"] | None = None
    plain: bool = False
    trace: str | None = None

    @model_validator(mode="after")
    def _rules(self):
        if (self.ckpt is None) == (self.select is None):
            raise ValueError("exactly one of ckpt and select")
        return self


class EvalStage(StageBase, frozen=True, extra="forbid"):  # final_eval
    """A final_eval stage."""

    kind: CliSuppress[Literal["eval"]] = "eval"
    ckpt: str | None = None
    train: str | None = None  # E4 deletes final_eval's --train
    view: tuple[str, ...] = Field(min_length=1)
    split: tuple[str, ...] = ()  # () -> every view uses "val"; else one split per view
    out: tuple[str, ...] = Field(min_length=1)  # one per view, or one
    model: str = "qwen3-vl-2b"
    k: int = 8
    seed: int = 11
    per_call: int = Field(64, gt=0)
    rescore: str | None = None

    @model_validator(mode="after")
    def _rules(self):
        if (self.ckpt is None) == (self.rescore is None):
            raise ValueError("exactly one of ckpt and rescore")
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
MODULE = {"rl": "train_rl", "sft": "train_sft", "select": "select_ckpt", "probe": "g_probe", "eval": "final_eval"}

# The order of the flags in today's plans; `views` stands for the interleaved view, split and out of an eval stage.
FLAG_ORDER = {
    "rl": ("merge_last", "model", "ckpt", "train", "out", "steps", "save_steps", "resume", "shared_prefix", "lr",
           "rows_revision", "rows_manifest_sha256", "budget_usd", "session_start", "rate_usd_per_h", "fixed_usd",
           "heldout", "heldout_items", "heldout_every", "grpo"),
    "sft": ("model", "view", "output_dir", "data_root", "run_dir", "model_path", "branches", "a0", "freeze_vision",
            "seed", "rolling_minutes", "vision", "vision_reference_only", "resume", "steps_per_epoch", "sft",
            "deterministic", "branches_from", "profile"),
    "select": ("run", "prepared", "probes", "sampler", "shifted", "rescaled", "rotated", "items", "vision",
               "no_vision", "only", "model", "k", "seed", "timeout_s", "out"),
    "probe": ("model", "ckpt", "select", "train", "probe", "out", "items", "grpo", "compile", "shared_prefix", "plain",
              "trace"),
    "eval": ("ckpt", "train", "views", "model", "k", "seed", "per_call", "rescore"),
}


class ComputeSpec(Base):
    """The rented GPU and the hour cap of a run."""

    gpu: str
    cap_hours: float = Field(gt=0)


class RunSpec(Base):
    """One run: its compute and its ordered stages."""

    schema_version: Literal[1]
    name: str
    compute: ComputeSpec
    stages: Annotated[tuple[Stage, ...], BeforeValidator(_tuple), Field(min_length=1)]
    bundle_view: str | None = None  # the SFT training view under the data root, read by the preflight DATA checks

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
    """The argv of a stage on the box: `~` first when soft, then `-m experiments.phase2.<module>` and the flags that
    the spec sets (in FLAG_ORDER), or a cmd stage's argv. {rd} and {data} are filled in string fields."""
    if stage.kind == "cmd":
        body = [fill(x, rd, data) for x in stage.argv]
    else:
        body = ["-m", f"experiments.phase2.{MODULE[stage.kind]}"]
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
    """The stage a CLI reads from its flags, in the one parser form (kebab flags, bare boolean flags)."""
    return CliApp.run(cls, cli_args=list(argv),
                      cli_settings_source=CliSettingsSource(cls, cli_kebab_case=True, cli_implicit_flags=True))


def same_stage(a: Stage, b: Stage) -> bool:
    """Equal but for `name` and `soft`, which belong to the spec only."""
    return a.model_dump(exclude={"name", "soft"}) == b.model_dump(exclude={"name", "soft"})

# Phase E spec: code structure (for test writers)

Revision 1 (2026-11-08): kept consistent with docs/implementation2.0/E.md after docs/implementation2.0/reviews/E_phase_review.md.

## Goal

`src/objdel` is one layered library with the layer order (high to low)
`infra > {rl | sft} > policy > data > verify > geom > core`. A layer imports only layers below it. `rl` and `sft`
never import each other. All touches of private TRL names are in `objdel.rl.trl_adapter`. Run plans are typed specs.
Evaluation samples without a trainer. The verifier has a versioned batch API. Rewards, targets, trainer updates and
run argv are identical to the pre-Phase-E code (goldens of B10).

## Required behaviours

B1. **Layers.**
- Packages that exist after Phase E: `objdel.core`, `objdel.geom`, `objdel.verify`, `objdel.data`, `objdel.policy`,
  `objdel.rl`, `objdel.sft`, `objdel.infra`, and the compiled extension `objdel._gpucheck`.
- Packages that no longer exist: `objdel.checker`, `objdel.geometry`, `objdel.reward`, `objdel.eval`,
  `objdel.targets`, `objdel.train`, and the module `objdel.train.fast`. There is no compatibility re-export module.
- `lint-imports` passes with a `layers` contract over container `objdel` in the order above. The contract is
  `exhaustive = true`, its only exhaustive ignore is `_gpucheck`, and `ignore_imports` is empty. Every contract has a
  `broken_contract_guidance` text.
- Before the end of Phase E, `lint-imports` passes after every commit with contracts that name only packages that
  exist at that commit, and the config has at least one contract.
- Only `objdel.geom.engine` imports `objdel._gpucheck`. `objdel.geom.engine.__all__` holds every `_gpucheck` name
  used anywhere in `src/`, `tests/` and `experiments/`.
- CI runs `lint-imports` (not only `ruff check`).

B2. **No ML stack below policy.** Import of `objdel.core`, `objdel.geom`, `objdel.verify` and `objdel.data` (all
submodules) loads none of: `torch`, `trl`, `transformers`, `vllm`, `peft`, `datasets`, `accelerate`.

B3. **Dependency direction.**
- No module under `src/objdel` imports `experiments`, `tests`, `oracles` or a module whose name starts with `test_`.
  No module under `src/objdel` changes `sys.path`.
- No test file outside `tests/experiments/` imports `experiments`.
- The box import-check string (`IMPORT_CHECK` in `objdel.infra.vast`) imports no `tests.*` and no `experiments.*`
  module, and it imports `pydantic`, `vllm` and `objdel.policy.sampler`.
- Running any live CLI under `experiments/` with `--help` loads no `tests.` or `oracles.` module.
- The CUDA parity harness lives under `tests/oracles/` (not in `src/`). The K3 reward harness lives under
  `tests/harness/`.

B4. **TRL adapter.**
- Only `src/objdel/rl/trl_adapter.py` imports `trl.trainer.grpo_trainer` or `trl.generation.vllm_generation`
  (test files of the adapter are exempt).
- Only that file accesses these names on TRL objects (as `self.<name>` in a subclass of `GRPOTrainer` or
  `VLLMGeneration`, as `trainer.<name>`, as `super().<name>`, or on a TRL module): `_get_per_token_logps_and_entropies`,
  `_generate_and_score_completions`, `_prepare_inputs`, `_compute_loss`, `_generate`, `_calculate_rewards`,
  `shuffle_sequence_dict`, `_iter_named_params`, `_metrics`, `_step`, `_buffered_inputs`. Other classes may use the
  same names.
- `check_private_names()` returns `None` on trl 1.17.1. It raises `TrlDriftError` when any of these is true: the TRL
  version is not `1.17.1`; one of the six overridden methods is missing; a method's parameter names differ from the
  recorded ones; the source SHA-256 of an overridden TRL method differs from the recorded one; a built trainer lacks
  one of the attributes `_metrics`, `_step`, `_buffered_inputs`, `num_iterations`. The error message contains the
  name and the pin `1.17.1`.
- `FastGRPO.__init__` calls `check_private_names()` before TRL's `__init__`.
- `env_check.sh` records the result of `check_private_names()` under the key `trl_private_names` in `env_check.json`.
- No code assigns to the module global `trl.trainer.grpo_trainer.shuffle_sequence_dict`.
- With the FastGRPO length ordering on (dr_grpo, shared prefix on or off), the micro-batches and the loss equal the
  recorded pre-change values on the tiny offline test model. In eval mode, and when the fallback condition holds
  (not training, loss not `dr_grpo`, entropy cut on, or `gradient_accumulation_steps % steps_per_generation != 0`),
  the output equals TRL's own `_prepare_inputs`.

B5. **Single sources.**
- Code reads the RL learning-rate default only from `objdel.rl.config.rl_lr(model)`: `5e-5` for LoRA models and
  `5e-6` for full fine-tuning. `grpo_config` and `build_trainer` use it when no learning rate is given.
- `experiments.phase2.prereg.G_PROBE` keeps `lr_full == 5e-6` and `lr_lora == 5e-5` (pre-registered record).
  For every model key, `rl_lr(m) == G_PROBE["lr_lora" if LoRA else "lr_full"]`. No code outside `prereg.py` reads
  `G_PROBE["lr_full"]` or `G_PROBE["lr_lora"]`.
- `objdel.verify.stats.WEIGHTS == dict(zip(objdel.verify.reward.TERMS, objdel.verify.reward.REWARD_V4.weights))`.
  No second weight literal exists.
- The host region rule is the constant `objdel.infra.vast.offers.REGION == ("US", "CA")`; no run-spec field sets it.
- `experiments` RL training creates no worker pool.

B6. **Trainer builder.**
- `objdel.rl.build.build_trainer` takes explicit typed arguments (model, ckpt, rows, out, roots, steps, flags). It
  accepts no argparse Namespace.
- `flags.plain = True` gives TRL's `GRPOTrainer`; `False` gives `FastGRPO`.
- For the same inputs, the resolved `GRPOConfig.to_dict()` (path and time fields removed) equals the pre-Phase-E value
  for the `train_rl` stage of `rl_shifted` and for a default G-probe (golden 6).
- No function named `make_trainer` exists in `src`, `experiments` or `tests`.

B7. **Sampler and evaluation.**
- `objdel.policy.sampler.build_prompt` gives the same prompt token ids as TRL's own prompt path
  (`prepare_multimodal_messages` + processor `apply_chat_template(add_generation_prompt=True)`) for `qwen3-vl-2b` on
  64 validation rows.
- `params_for(model, n=n, seed=seed, stop_token_ids=s)` returns a dict equal to
  `objdel.policy.models.sampling(model, "vllm") | {"n": n, "seed": seed, "stop_token_ids": list(s)}`, including
  `min_p`, `repetition_penalty`, `presence_penalty` and `frequency_penalty`.
- `requests(...)` gives each request the seed `request_seed(seed, item_id)` and `mm_processor_kwargs` when the line
  has processor kwargs. `ENV["VLLM_USE_FLASHINFER_SAMPLER"] == "0"` is set before vLLM is imported.
- `HfSampler.sample` returns one `Sample` per row, with `len(texts) == len(token_ids) == len(truncated) == n`.
  `truncated[i] == truncated(token_ids[i], stop_token_ids)`.
- `python -I -m objdel.policy.sampler --kind vllm|hf ...` is the sampler child process; `select_ckpt.sampler_cmd`
  builds that command.
- `experiments/phase2/sample_vllm.py` and `experiments/phase2/hf_sample.py` exist only under `archive/`.
- The evaluation path (`objdel.policy.evaluate.evaluate`) builds no trainer and calls no `_generate`.
- Importing `objdel.policy.sampler` and `objdel.policy.evaluate` imports neither `trl.trainer` nor `vllm`.

B8. **Run spec.**
- `load_spec(path)` returns a frozen `RunSpec`. Unknown keys at any level raise `pydantic.ValidationError`
  (including `compute.regions`). An unknown `kind` raises. `schema_version` other than `1` raises. Duplicate stage
  names raise.
- `RlStage` without `save_steps` has `save_steps == 50`; `save_steps > steps` raises. `ProbeStage` with both or
  neither of `ckpt` and `select` raises.
- Every file in `specs/*.yaml` loads.
- `to_argv(stage, rd, data)` on the ported `specs/rl_shifted.yaml`, `specs/rl_full.yaml` and `specs/r5_shared.yaml`
  gives, string for string, the argv lists that the pre-Phase-E `policy.chain()` gave for the old YAML files (golden
  5, after the recorded module rename map). For `specs/r1_qwen3vl2b.yaml`, the select, probe and cmd stages match
  golden 5 string for string, and the `SftStage` field values equal `train_sft.RunArgs` built from the old YAML file.
- Each live CLI parses its old flags and `--spec PATH --stage NAME` into the same stage object.
- `CmdStage.argv` is passed through literally, after `{rd}` / `{data}` substitution.
- Preflight takes a run spec. It passes on `specs/rl_shifted.yaml` with fixture data. After the spec file changes,
  the launch refusal fires (stamp keyed on the spec path and its SHA-256).
- Each run writes `runs/<id>/run_spec.json`. It holds the spec, the git sha, the trl, torch and transformers versions,
  and the reward fingerprint.

B9. **Verifier API.**
- `Verifier("reward_v4", view_root).score(items)` gives the same floats (bit-equal via `float.hex`) as
  `objdel.verify.reward.score_batch` on the same inputs.
- `Verifier.stamp` holds `api_version == API_VERSION`, `reward_name`, `reward_fingerprint == fingerprint()`,
  `checker_build in {"cuda", "cpu"}` and `code_version`.
- A key that is not in the store raises `VerifyError` with `code == "unknown_key"` for the whole batch. An outline
  that fails to parse gives a normal `ScoreResult` with `error` set (not an exception).
- HTTP (only when the `serve` extra is installed): `POST /v1/score` returns the same floats as the in-process call.
  Every response body carries `api_version`. A batch larger than `MAX_BATCH` returns HTTP 413 with
  `error.code == "batch_too_large"`. `GET /v1/spec`, `GET /livez` and `GET /readyz` return 200 when ready. Error body:
  `{"error": {"code", "message", "retryable"}}`, and `retryable` is true only for `device_error`.
- The CI coverage job installs the `serve` extra, so the HTTP tests run there.

B10. **No result change (move invariance).** These goldens, recorded on the pre-Phase-E commit, stay equal:
1. `sha256(repr((astuple(REWARD_V4), [min_area(a) for a in (0, 999, 1001, 10**4, 10**6)])))`. Also the full
   `fingerprint()` hex string, checked only when numpy, OpenCV and scipy have the recorded versions.
2. `target_set_hash` of the targets of the fixture items in the target tests.
3. `score_batch` output on a fixed corpus of 64 polygon texts (CPU build), each float as `float.hex`. Per platform
   (`invariance.<sys.platform>-<machine>.json`).
4. Loss, gradient L2 norm and parameter-change L2 norm (as `float.hex`) of one FastGRPO update on the tiny offline
   test model, CPU, seed 0. Per platform.
5. The chain argv lists for `rl_shifted`, `rl_full`, `r5_shared` and `r1_qwen3vl2b`, compared after the recorded
   module rename map of the experiments move.
6. The resolved `GRPOConfig.to_dict()` (path and time fields removed) for the `train_rl` stage of `rl_shifted` and
   for a default G-probe on `qwen3-vl-2b`.

A platform with no recorded golden file skips that file with a message. darwin-arm64 is the binding platform.

B11. **No new host syncs.** For each module under `src/objdel` that imports torch, the count of regex matches
`\.(item|cpu|tolist|numpy)\(\)` is <= its recorded count in `tests/goldens/phase_e/host_syncs.json` (a moved file is
matched through the recorded old -> new path map).

B12. **No stale paths or strings.**
- Every string literal matching `objdel(\.\w+)+` in `src/`, live `experiments/` and `tests/` resolves statically to
  a module file (or `objdel._gpucheck`) plus, optionally, a top-level name defined in it.
- Every string matching `tests/[\w/]+\.py(::\w+)?` in `pyproject.toml`, `.github/`, live `experiments/`, `specs/`,
  `src/` and `Makefile` names an existing file (and test function, when given).
- `objdel.infra.paths.repo_root()` returns the directory that holds `pyproject.toml` (env `OBJDEL_REPO` first, then
  git). With `OBJDEL_SECRETS` unset, `objdel.infra.secrets.path() == repo_root() / ".secrets"`. No module under
  `src/objdel/infra/` derives a root from `__file__`.
- The box code snapshot (`code_files()`) lists no path under `.claude/` or `archive/`, and it contains
  `tests/test_gpu_scoring.py`.

B13. **Verifier store and seams.**
- `objdel.verify.store` holds `preload`, `loader`, `VIEW_SEP`; `objdel.verify.pool` holds `fork_pool`, `pool_size`.
  No module under `objdel.verify` names a module of a higher layer, in an import or in a string.
- `fork_pool(n, modules)` imports exactly the modules passed (no built-in list).
- `objdel.geom.polygonize(mask, offset=(0, 0))` is the one polygonizer entry point that production code calls.

B14. **Every name has a home.** For each pre-Phase-E module that is split or moved (`train/dataset.py`,
`train/configs.py`, `train/fast.py`, `train/rewards.py`, `train/callbacks.py`, `train/serve.py`), its set of top-level
names (recorded on the base commit) equals the union of the top-level names of its new modules, minus the names
listed as archived or renamed.

## Interfaces

```python
# objdel.core.spec  (pydantic v2; all models frozen=True, extra="forbid"; stages are never hashed)
class TrainerFlags(BaseModel):
    lr: float | None = None; shared_prefix: bool = False
    compile: Literal["default", "reduce-overhead"] | None = None; plain: bool = False; trace: str | None = None
    grpo: Mapping[str, JsonValue] = {}
class RlStage(TrainerFlags):
    kind: Literal["rl"] = "rl"; name: str; model: str; ckpt: str; train: str; train_views: tuple[str, ...] = ()
    out: str; steps: int  # > 0
    save_steps: int = 50  # > 0, <= steps
    resume: bool = False
class SftStage(BaseModel): kind: Literal["sft"]; name: str; ...      # the fields of the former train_sft.RunArgs
class SelectStage(BaseModel):
    kind: Literal["select"]; name: str; run: str; prepared: str; probes: str
    sampler: Literal["vllm", "hf"] | None = None; items: int | None = None
    shifted: str | None = None; rescaled: str | None = None; rotated: str | None = None
    vision: VisionRef | None = None; only: Literal["topology"] | None = None
class ProbeStage(TrainerFlags):
    kind: Literal["probe"]; name: str; model: str; train: str; probe: str; out: str
    ckpt: str | None = None; select: str | None = None   # exactly one
    items: int | None = None
class EvalStage(BaseModel):
    kind: Literal["eval"]; name: str; ckpt: str; view: str; split: Literal["val", "test"]; train: str; out: str
    k: int = 8
class CmdStage(BaseModel): kind: Literal["cmd"]; name: str; argv: tuple[str, ...]
Stage = Annotated[RlStage | SftStage | SelectStage | ProbeStage | EvalStage | CmdStage, Field(discriminator="kind")]
class ComputeSpec(BaseModel): gpu: str; cap_hours: float  # > 0 ; no regions field
class RunSpec(BaseModel): schema_version: Literal[1]; name: str; compute: ComputeSpec; stages: tuple[Stage, ...]
def load_spec(path: Path) -> RunSpec
def to_argv(stage: Stage, rd: str, data: str) -> list[str]

# objdel.rl.build
class Built(NamedTuple): trainer: Any; processor: Any; lr: float
def build_trainer(model: str, ckpt: Path, rows: "datasets.Dataset", out: Path, roots: Sequence[Path], steps: int,
                  flags: TrainerFlags, **over) -> Built

# objdel.rl.run
def run_rl(stage: RlStage, rd: Path) -> dict
def save_json(obj: Any, path: Path) -> None      # NaN/inf rejected; numpy scalars and arrays serialised

# objdel.rl.trl_adapter
PIN: str = "1.17.1"
class TrlDriftError(RuntimeError): ...
def check_private_names() -> None
class FastGRPO(GRPOTrainer): ...

# objdel.rl.config
def rl_lr(model: str) -> float
def grpo_config(model: str, output_dir: str, **overrides) -> "trl.GRPOConfig"

# objdel.policy.rows
def prompts(root, dataset: str, split: str, rl_only: bool = True, oversample: int = 1, subset=None) -> list[dict]
def records(root, dataset: str, split: str) -> list[dict]
def cache_loader(root, maxsize: int | None = 64) -> Callable[[str], dict]
DISAGREEMENT_OVERSAMPLE: int

# objdel.policy.sampler
@dataclass(frozen=True)
class Sample:
    item_id: str; texts: tuple[str, ...]; token_ids: tuple[tuple[int, ...], ...]; truncated: tuple[bool, ...]
class Sampler(Protocol):
    def sample(self, rows: Sequence[Mapping], params: Mapping[str, Any]) -> list[Sample]: ...
class VllmSampler:
    def __init__(self, ckpt: Path, model: str, gpu_memory_utilization: float, max_model_len: int): ...
    def sample(self, rows, params) -> list[Sample]: ...
class HfSampler:
    def __init__(self, ckpt: Path, model: str, device: str = "cpu"): ...
    def sample(self, rows, params) -> list[Sample]: ...
ENV: dict[str, str]                               # {"VLLM_USE_FLASHINFER_SAMPLER": "0"}
def build_prompt(row: Mapping, processor: Any, model: str) -> tuple[str, list]
def params_for(model: str, *, n: int, seed: int, stop_token_ids: Sequence[int]) -> dict
def request_seed(seed: int, item_id: str) -> int
def check_memory(mem_get_info, need: float = ...) -> None
def truncated(ids: Sequence[int], stop_token_ids: Sequence[int]) -> bool
def main(argv: Sequence[str] | None = None) -> None   # child-process entry

# objdel.policy.evaluate
def evaluate(sampler: Sampler, rows: Sequence[Mapping], params: Mapping[str, Any], view_root: Path) -> dict

# objdel.geom
def polygonize(mask, offset: tuple[int, int] = (0, 0)) -> list   # the one polygonizer entry point
# objdel.geom.engine: explicit re-exports of objdel._gpucheck with __all__

# objdel.verify.store / objdel.verify.pool
def preload(root, datasets: Iterable[str], threads: int = 4) -> int
def loader(*roots) -> Callable[[str], Mapping]
VIEW_SEP: str = "::"
def fork_pool(n: int, modules: Iterable[str] = ()) -> "ProcessPoolExecutor"
def pool_size(cgroup: Path = ...) -> int

# objdel.verify.api  (no pydantic, no torch)
API_VERSION: str = "1.0"
@dataclass(frozen=True)
class ScoreItem: id: str; key: str; text: str; truncated: bool
@dataclass(frozen=True)
class ScoreResult:
    id: str; reward: float | None; terms: Mapping[str, float]; defects: Mapping[str, int]; error: str | None
@dataclass(frozen=True)
class Stamp: api_version: str; reward_name: str; reward_fingerprint: str; checker_build: str; code_version: str
class VerifyError(Exception):
    code: Literal["unknown_key", "spec_mismatch", "batch_too_large", "device_error"]
class Verifier:
    stamp: Stamp
    def __init__(self, reward: str, view_root: Path, device: int = 0): ...
    def score(self, items: Sequence[ScoreItem]) -> list[ScoreResult]: ...
    def spec(self) -> dict: ...

# objdel.verify.http  (optional extra "serve": fastapi==0.136.3, uvicorn==0.54.0)
app: fastapi.FastAPI       # POST /v1/score {"reward": str, "items": [ScoreItem fields]} ; GET /v1/spec /livez /readyz
MAX_BATCH: int

# objdel.infra.paths
def repo_root() -> Path
```

CLIs (thin; run as modules with `python -I -m`): each live script under `experiments/<group>/` takes
`--spec PATH --stage NAME` or its old flags (both build the same stage object) and answers `--help` with exit 0.

## Constraints

- Python 3.12. Pins: trl==1.17.1, import-linter==2.16 (dev), ruff==0.17.1. Base dependencies as ranges:
  `pydantic>=2.14,<3`, `pydantic-settings>=2.16,<3` (exact versions only in the image locks). `serve` extra:
  fastapi==0.136.3, uvicorn==0.54.0; dev: httpx==0.28.1.
- `make test` (= `pytest -q -m "not slow"`) and `make lint` (= ruff check + `lint-imports`) pass after every commit.
- Nothing in Phase E needs a GPU. All tests above run on a Mac CPU with the CPU build of `objdel._gpucheck`. GPU-only
  checks are marked and skipped without a GPU.
- The vLLM sampler runs only in the vllm_native image (vLLM and objdel in one venv).
- Python that reads `data/` runs with `-I`.
- No shims: the old module paths raise `ModuleNotFoundError`.
- Retired code lives under `archive/<name>/` with a README. `archive/` is not imported, tested, linted or shipped to
  rented boxes.
- No Claude/AI attribution in commits.

## Acceptance gates

| # | Gate | Pass condition |
|---|---|---|
| G1 | Fast suite | 0 failed; passed >= 1447 minus the removed tests that the commit lists (baseline: 1447 passed, 6 skipped, 31 deselected) |
| G2 | Lint | ruff 0 errors; lint-imports: all contracts kept, locally and in CI; final contract exhaustive with 0 `ignore_imports` |
| G3 | Move invariance | goldens 1-6 of B10 equal on darwin-arm64 (other platforms: where a file is recorded) |
| G4 | No ML stack below policy | B2 holds |
| G5 | TRL in one place | B4 holds; `check_private_names()` passes and fails in 3 negative cases (missing method, changed signature, changed source hash); it finishes in < 10 s on the Mac, including the TRL import |
| G6 | No global patch | 0 assignments to `shuffle_sequence_dict`; micro-batches and loss equal to the recorded values; eval-mode and fallback output equal TRL's |
| G7 | Single sources | B5 holds |
| G8 | Direction | B3 holds |
| G9 | Typed plans | B8 holds |
| G10 | Evaluation without a trainer | B7 holds |
| G11 | No new host syncs | B11 holds |
| G12 | Verifier API | B9 holds |
| G13 | Box smoke (at the next approved rental) | `env_check.json` has `trl_private_names` passed, and `IMPORT_CHECK` exits 0, before training. Box sampler check: prompt ids and SamplingParams fields of the new sampler equal TRL's on 64 val-shifted rows (binding); mean IoU_tc difference trainer path vs sampler path on 1,084 val-shifted items, K 8: 90 % nodule-level bootstrap CI inside ±0.01 (else "not shown") |
| G14 | Old code gone | the old packages of B1 do not exist; no `make_trainer`; no `fast.py`; `sample_vllm.py` and `hf_sample.py` only under `archive/` |
| G15 | No stale strings | B12 holds |
| G16 | Every name has a home | B14 holds |

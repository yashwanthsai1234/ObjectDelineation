"""Start-up check of the private TRL names and the silent couplings that FastGRPO relies on (TRL pinned to 1.17.1).

`check()` runs first in `FastGRPO.__init__` and in the preflight ("CODE TRL private names"). It reads signatures and
source text only, so it runs on the Mac without vLLM or a GPU. A TRL change that breaks one item fails in seconds,
naming the item and the place in `fast.py` that overrides or reads it. The version is checked last: with
OBJDEL_TRL_UNPINNED=1 a different TRL version only warns, so an upgrade can be tried on the contract alone.
"""

import inspect
import os
import warnings

PIN = "1.17.1"

# name -> (parameter names in order, the place in fast.py that overrides or calls it)
METHODS = {  # GRPOTrainer
    "_generate_and_score_completions": (("self", "inputs"), "FastGRPO._generate_and_score_completions"),
    "_prepare_inputs": (("self", "generation_batch"), "FastGRPO._prepare_inputs"),
    "_compute_loss": (("self", "model", "inputs"), "FastGRPO._compute_loss"),
    "_get_per_token_logps_and_entropies": (
        ("self", "model", "input_ids", "attention_mask", "logits_to_keep", "batch_size", "compute_entropy",
         "compute_aux_loss", "pixel_values", "image_grid_thw", "num_images", "pixel_attention_mask", "spatial_shapes",
         "num_tiles", "image_sizes", "token_type_ids", "mm_token_type_ids", "image_position_ids", "sampling_mask"),
        "FastForward._get_per_token_logps_and_entropies"),  # fmt: skip
    "_generate": (("self", "prompts"), "FastGRPO._generate"),
    "_calculate_rewards": (("self", "inputs", "prompts", "completions", "completion_ids_list"),
                           "FastGRPO._calculate_rewards"),  # fmt: skip
}
GLOBALS = {"shuffle_sequence_dict": (("seq_dict",), "FastGRPO._prepare_inputs")}  # trl.trainer.grpo_trainer
VLLM_METHODS = {"_iter_named_params": (("self",), "FastGRPO.sync_without_merge")}  # VLLMGeneration

# name -> (sources that must contain the needle, needle, the place in fast.py that relies on it)
PRIVATE_ATTRS = {
    "_metrics": (("GRPOTrainer.__init__",), "self._metrics", "FastGRPO.timed, FastGRPO.log"),
    "_dist": (("VLLMGeneration",), "self._dist", "FastGRPO.sync_without_merge"),
}
ATTRS = {
    **{n: (("GRPOTrainer",), f"self.{n} =", p) for n, p in (
        ("use_vllm", "FastGRPO.skip_old_logps"), ("loss_type", "FastGRPO._prepare_inputs"),
        ("num_iterations", "FastGRPO.skip_old_logps"), ("num_generations", "FastGRPO._prepare_inputs"),
        ("top_entropy_quantile", "FastGRPO._prepare_inputs"),
        ("vllm_importance_sampling_correction", "FastGRPO._generate_and_score_completions"),
        ("vllm_importance_sampling_mode", "FastGRPO.set_is_ratio"),
        ("vllm_importance_sampling_clip_min", "FastGRPO.set_is_ratio"),
        ("vllm_importance_sampling_clip_max", "FastGRPO.set_is_ratio"),
        ("vllm_generation", "FastGRPO.sync_without_merge"))},
    "llm": (("VLLMGeneration",), "self.llm =", "FastGRPO.time_vllm_phases, FastGRPO.patch_vllm_prompts"),
    "model": (("VLLMGeneration",), "self.model =", "FastGRPO.sync_without_merge"),
    "sync_weights": (("VLLMGeneration",), "def sync_weights", "FastGRPO.time_vllm_phases"),
}  # fmt: skip
_GASC = "GRPOTrainer._generate_and_score_completions"
KEYS = {
    "sampling_per_token_logps": ((_GASC,), '"sampling_per_token_logps"', "FastGRPO.set_is_ratio"),
    "completion_mask": ((_GASC,), '"completion_mask"', "FastGRPO.set_is_ratio"),
    "tool_mask": ((_GASC,), '"tool_mask"', "FastGRPO.set_is_ratio"),
    "importance_sampling_ratio": ((_GASC, "GRPOTrainer._compute_loss"), '"importance_sampling_ratio"',
                                  "FastGRPO.set_is_ratio, FastGRPO._compute_loss"),
}  # fmt: skip
BEHAVIOURS = {  # C10: TRL skips the old-logprob pass when the IS correction is off; fast.is_ratio copies this block
    "old-logprob pass gated on the IS correction": (
        (_GASC,), "self.use_vllm and self.vllm_importance_sampling_correction", "FastGRPO.skip_old_logps"),
    "IS ratio block": ((_GASC,), "per_token_logps_diff = (old_per_token_logps - sampling_per_token_logps) * mask",
                       "fast.is_ratio"),
    "log_metric to reward functions": (("GRPOTrainer._calculate_rewards",), 'reward_kwargs["log_metric"]',
                                       "TRLReward.funcs"),
    "log_extra to reward functions": (("GRPOTrainer._calculate_rewards",), 'reward_kwargs["log_extra"]',
                                      "TRLReward.funcs"),
    "completion_ids to reward functions": (("GRPOTrainer._calculate_rewards",), "completion_ids=completion_ids_list",
                                           "TRLReward.funcs"),
}  # fmt: skip
_SIGNED = {"method": METHODS, "global": GLOBALS, "VLLMGeneration method": VLLM_METHODS}
_TEXT = {"private attribute": PRIVATE_ATTRS, "attribute": ATTRS, "key": KEYS, "behaviour": BEHAVIOURS}


class TRLContractError(RuntimeError):
    """A private TRL name or coupling that fast.py relies on changed."""


def _fail(kind: str, name: str, expected, found, place: str):
    raise TRLContractError(f"TRL {kind} {name}: expected {expected!r}, found {found!r}; relied on by fast.py {place}. "
                           f"Pin is {PIN}.")  # fmt: skip


def check() -> None:
    """Raise TRLContractError on the first broken item (methods, globals, private attributes, attributes, keys,
    behaviours), then on a version other than PIN unless OBJDEL_TRL_UNPINNED=1 (then one warning)."""
    import trl
    import trl.trainer.grpo_trainer as grpo
    from trl.generation.vllm_generation import VLLMGeneration  # imports no vLLM module

    owners = {"method": grpo.GRPOTrainer, "global": grpo, "VLLMGeneration method": VLLMGeneration}
    for kind, items in _SIGNED.items():
        for name, (params, place) in items.items():
            fn = getattr(owners[kind], name, None)
            found = "missing" if fn is None else tuple(inspect.signature(fn).parameters)
            if found != params:
                _fail(kind, name, params, found, place)
    sources = {}
    for kind, items in _TEXT.items():
        for name, (where, needle, place) in items.items():
            for w in where:
                if w not in sources:
                    obj = {"GRPOTrainer": grpo.GRPOTrainer, "VLLMGeneration": VLLMGeneration}[w.split(".")[0]]
                    for attr in w.split(".")[1:]:
                        obj = getattr(obj, attr)
                    sources[w] = inspect.getsource(obj)
                if needle not in sources[w]:
                    _fail(kind, name, f"{needle} in the source of {w}", "absent", place)
    if trl.__version__ != PIN:
        msg = f"TRL {trl.__version__} is installed; every contract item passed. Pin is {PIN}."
        if os.environ.get("OBJDEL_TRL_UNPINNED") != "1":
            raise TRLContractError(msg + " Set OBJDEL_TRL_UNPINNED=1 to try this version.")
        warnings.warn(msg, UserWarning, stacklevel=2)

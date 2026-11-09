"""Reward v4 as TRL (version pinned in pyproject.toml) GRPOTrainer reward functions: weight-0 per-term loggers
plus one aggregator (plan D7, §4.6).

TRL calls every reward function once per generation batch with the same arguments (prompts, completions,
completion_ids, every dataset column, trainer_state, log_metric, log_extra) and gathers their outputs into
rewards_per_func. With multi_objective_aggregation="sum_then_normalize" and scale_rewards="none" its advantage is
    A = sum_k w_k r_k - group_mean(sum_k w_k r_k).
The aggregator (weight 1) returns objdel.reward.advantages_v4, whose group means are already 0, and the loggers
(weight 0) return one term each, so A = advantages_v4 to float32 precision (R2, tests/test_train_trl.py) and every
term is logged as rewards/<term>/mean. The geometric loggers (GATED: iou_tc, f5, g_area) return the term times
objdel.reward.legal_gate, as advantages_v4 sees it: an illegal rollout logs 0, so neither the curves nor the
held-out eval_rewards/iou_tc/mean credit illegal output (R2 follow-up). A weight-0 logger adds 0 * value = 0 to
TRL's sum, so the advantages do not depend on what the loggers return. A logger given weight w (GRPOConfig
reward_weights) adds w (term - group mean), with the gated term. TRL's own
`reward` metric is then the mean advantage, about 0: dashboards use the per-term logs.

All functions share one computation per batch: the first call scores the batch (objdel.reward.score_batch: one CUDA
checker call checks and scores all its answers against the item caches resident on the device), the others reuse it
(the same `completions` list object is passed to each function).

Groups: advantages_v4 needs each prompt's G completions contiguous within one call. TRL computes rewards on each
process's local slice before gathering; TRL's RepeatSampler yields each prompt G times in a row and accelerate
gives process p the p-th block of per_device_train_batch_size * steps_per_generation indices, so a local slice
holds whole groups exactly when that size is a multiple of G (objdel.train.configs guarantees it). The call
checks it anyway (one key per block of G) and raises otherwise.

Truncation: TRL's HF and vLLM generation cut completion_ids after the first EOS; a completion is truncated when
its last id is not an EOS or pad id (GRPOTrainer: `ids[-1] not in eos_and_pad`). The reward uses the same rule
(`truncated` on `stop_ids`, also objdel.train.sample's), so a truncated completion is a format failure (plan D10)
even when its prefix parses.

Typed logs (log_metric, averaged by TRL per logging step, eval-prefixed in evaluation): error/<parser code> is
the share of completions with that error (codes enter the set when first seen and are logged every batch after),
defect/<class> the share with that checker defect, allmiss_share the share of groups whose answers are all
missed answers (gated IoU_tc < MISS_IOU); rows with the kind and cell columns (a B1 rows root) add kind/<k>/iou_tc,
cell/<ij>/iou_tc and cell/<ij>/allmiss_groups. log_extra adds the columns key (item id), truncated and group
(index of the completion's group within its reward call) to the completions table (log_completions=True), so
objdel.train.sample.export can turn it into eval_run JSONL (break-it B4-3), and error (the parser code).

Group size: TRL's num_generations must be REWARD_V4.group, or v4 would silently score each TRL group as several
independent groups (break-it B2-5). The keys cannot show it (an oversampled item drawn twice in a row gives 16
equal keys), so objdel.train.configs.grpo_config refuses any other value.

References:
- TRL GRPOTrainer docs (reward function kwargs, log_metric / log_extra, multi_objective_aggregation,
  scale_rewards, mask_truncated_completions): https://huggingface.co/docs/trl/main/en/grpo_trainer
- TRL (version pinned in pyproject.toml) source, trl/trainer/grpo_trainer.py: `_calculate_rewards` (local call,
  then gather),
  `_generate_and_score_completions` (sum_then_normalize, process_slice), `_generate_single_turn` (mask after the
  first EOS), eos_token_ids (generation_config.eos_token_id + tokenizer.eos_token_id); trl/trainer/utils.py
  RepeatSampler: https://github.com/huggingface/trl/blob/v1.17.1/trl/trainer/grpo_trainer.py
- accelerate DataLoaderShard / BatchSamplerShard (process p takes batch p of every num_processes batches):
  https://huggingface.co/docs/accelerate/main/en/package_reference/torch_wrappers
- Dr. GRPO (no std scaling): https://arxiv.org/abs/2503.20783 ; GDPO (per-term normalization, the arm TRL's
  normalize_then_sum implements and v4 replaces): https://arxiv.org/abs/2601.05242
- Lab: archive/phase1_lab/t_stack/prod.py (TRLReward), r2_trl_equiv.py, results/INTEGRATION.md.
"""
from collections.abc import Callable, Iterable, Mapping
from concurrent.futures import ThreadPoolExecutor

import numpy as np

from objdel._gpucheck import DEFECTS
from objdel.reward import REWARD_V4, ExpertCache, advantages_v4, legal_gate, score_batch, topology_reward
from objdel.train import nvtx

LOGGED = ("format_ok", "legal", "iou_tc", "f5", "g_area", "topo_any", "r_topo")
GATED = ("iou_tc", "f5", "g_area")  # logged times legal_gate: 0 for an illegal rollout, as advantages_v4 sees them
AGGREGATOR = "advantage_v4"
MISS_IOU = 0.1  # a missed answer: gated IoU_tc below this (glossary)


def stop_ids(model, processing_class) -> set[int]:
    """The ids that end a completion for TRL: generation_config.eos_token_id, the tokenizer's EOS and pad."""
    tok = getattr(processing_class, "tokenizer", processing_class)
    eos = model.generation_config.eos_token_id
    return {*(eos if isinstance(eos, list) else [eos]), tok.eos_token_id, tok.pad_token_id} - {None}


def truncated(ids, stop) -> bool:
    """TRL's rule: a completion is truncated when its last id is not a stop id (an empty one has no end either)."""
    return not ids or ids[-1] not in stop


def _score(jobs) -> dict:
    """score_batch inside the NVTX range `checker`, on whichever thread scores (the prefetch thread or the reward
    call). The range lives here, not in objdel.reward, which never imports objdel.train (S12)."""
    with nvtx.range("checker"):
        return score_batch(jobs)


def _text(completion) -> str:
    """Conversational completions are [{"role": "assistant", "content": str | [parts]}]."""
    if isinstance(completion, list):
        content = completion[0]["content"]
        return content if isinstance(content, str) else "".join(p.get("text", "") for p in content)
    return completion


class TRLReward:
    """caches(key) -> {object name: ExpertCache} for one prompt (objdel.train.dataset.cache_loader, or
    serve.loader(root) on preloaded caches); stop_ids: the ids TRL treats as the end of a completion
    (`stop_ids(model, processor)`). Each reward call (256 answers) is one CUDA checker call that checks and scores
    them all; each cache's tables are uploaded to the device the first time it is scored (score_batch)."""

    def __init__(self, caches: Callable[[str], Mapping[str, ExpertCache]], stop_ids: Iterable[int]):
        self.caches, self.stop = caches, frozenset(stop_ids)
        self.codes: set[str] = set()
        self._last = (None, None)
        self.by_prompt: list[float] | None = None  # the last call's gated IoU_tc, mean per prompt (B2.3 bootstrap)
        self._ahead, self._thread = None, None

    def funcs(self) -> tuple[list[Callable], list[float]]:
        """(reward_funcs, reward_weights) for GRPOTrainer / GRPOConfig."""
        def make(name):
            def f(completions, completion_ids, key, W, H, names, log_metric=None, log_extra=None, kind=None,
                  cell=None, **_):
                return self._terms(completions, completion_ids, key, W, H, names, log_metric, log_extra, kind,
                                   cell)[name]
            f.__name__ = name
            f.trl_reward = self  # lets FastGRPO.prefetch_rewards find the reward among trainer.reward_funcs
            return f
        return [make(n) for n in (*LOGGED, AGGREGATOR)], [0.0] * len(LOGGED) + [1.0]

    def prefetch(self, completions, completion_ids, key, W, H, names) -> None:
        """Start scoring these answers on a background thread (the checker releases the GIL), so scoring overlaps
        TRL's old-logprob and reference forwards. The next reward call uses the
        result only for exactly these answers; otherwise it scores as usual."""
        keys = list(key)
        jobs = self._jobs(completions, completion_ids, keys, W, H, names)
        if self._thread is None:
            self._thread = ThreadPoolExecutor(1, thread_name_prefix="reward-prefetch")
        self._ahead = (self._signature(keys, jobs), self._thread.submit(_score, jobs))

    def _jobs(self, completions, completion_ids, keys, W, H, names) -> list[tuple]:
        return [(_text(c), truncated(ids, self.stop), tuple(n), self.caches(k), w, h)
                for c, ids, k, w, h, n in zip(completions, completion_ids, keys, W, H, names, strict=True)]

    @staticmethod
    def _signature(keys, jobs) -> list[tuple]:
        return [(k, t, tr, n, w, h) for k, (t, tr, n, _, w, h) in zip(keys, jobs, strict=True)]

    def _terms(self, completions, completion_ids, key, W, H, names, log_metric, log_extra, kind=None,
               cell=None) -> dict:
        if self._last[0] is completions:
            return self._last[1]
        g = REWARD_V4.group
        keys = list(key)
        if len(keys) % g or any(keys[i:i + g].count(keys[i]) != g for i in range(0, len(keys), g)):
            raise ValueError(f"reward call must receive whole groups of {g} completions per prompt, in order; got "
                             f"{len(keys)} completions (per_device_train_batch_size * steps_per_generation and "
                             f"per_device_eval_batch_size must be multiples of {g})")
        jobs = self._jobs(completions, completion_ids, keys, W, H, names)
        ahead, self._ahead = self._ahead, None
        same = ahead is not None and ahead[0] == self._signature(keys, jobs)
        T = ahead[1].result() if same else _score(jobs)
        legal = legal_gate(T)
        self.by_prompt = (T["iou_tc"] * legal).reshape(-1, g).mean(1).tolist()
        out = {n: (T[n] * legal if n in GATED else T[n]).tolist() for n in LOGGED[:-1]}
        out["r_topo"] = topology_reward(T).tolist()
        out[AGGREGATOR] = advantages_v4(T).tolist()
        if log_metric is not None:
            self.codes.update(e for e in T["error"] if e)
            for code in sorted(self.codes):
                log_metric(f"error/{code}", T["error"].count(code) / len(jobs))
            for d in DEFECTS:
                log_metric(f"defect/{d}", sum(d in x for x in T["defects"]) / len(jobs))
            self._log_views(log_metric, T["iou_tc"] * legal, kind, cell)
        if log_extra is not None:
            log_extra("key", keys)
            log_extra("truncated", [j[1] for j in jobs])
            log_extra("group", [i // g for i in range(len(jobs))])
            log_extra("error", T["error"])
        self._last = (completions, out)
        return out

    @staticmethod
    def _log_views(log_metric, iou, kind, cell) -> None:
        """allmiss_share: the share of groups whose G answers are all missed answers (gated IoU_tc < MISS_IOU); with
        the rows root's kind and cell columns also kind/<k>/iou_tc, cell/<ij>/iou_tc and cell/<ij>/allmiss_groups."""
        g = REWARD_V4.group
        allmiss = (iou < MISS_IOU).reshape(-1, g).all(1)
        log_metric("allmiss_share", float(allmiss.mean()))
        for col, vals in (("kind", kind), ("cell", cell)):
            vals = None if vals is None else np.asarray(vals)
            for v in [] if vals is None else sorted(set(vals)):
                log_metric(f"{col}/{v}/iou_tc", float(iou[vals == v].mean()))
                if col == "cell":
                    log_metric(f"cell/{v}/allmiss_groups", float(allmiss[vals[::g] == v].mean()))

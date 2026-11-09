"""Reward v4 (plan §3.2, §4.3): per-prompt expert caches, per-rollout terms, group advantages.

    cache = ExpertCache.build(rater_stack)                       # offline, once per item and object
    terms = score_batch([(text, truncated, names, caches, W, H), ...])   # on the CUDA checker
    adv = advantages_v4(terms)                                     # (B*G,), prompt by prompt
"""
from .advantages import advantages_v4, legal_gate, loo_bias, topology_reward, topology_valid
from .cache import ExpertCache, StaleCacheError, rl_eligible, voters
from .score import (
    NUMERIC,
    ObjectTerms,
    rater_answers,
    score_batch,
    score_completion,
    score_many,
    score_mask,
    sig_id,
)
from .spec import CACHE_SPEC, REWARD_V4, TERMS, fingerprint

__all__ = [
    "CACHE_SPEC", "NUMERIC", "REWARD_V4", "TERMS", "ExpertCache", "ObjectTerms", "StaleCacheError", "advantages_v4",
    "fingerprint", "legal_gate", "loo_bias", "rater_answers", "rl_eligible", "score_batch", "score_completion",
    "score_many", "score_mask", "sig_id", "topology_reward", "topology_valid", "voters",
]

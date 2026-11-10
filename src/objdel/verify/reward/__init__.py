"""Reward v4 (plan §3.2, §4.3) and v5 (Implementation.md B3.2-B3.4: location term, RLOO; advantages_v5):
per-prompt expert caches, per-rollout terms, group advantages.

    cache = ExpertCache.build(rater_stack)                       # offline, once per item and object
    terms = score_batch([(text, truncated, names, caches, W, H), ...])   # on the CUDA checker
    adv = advantages_v4(terms)                                     # (B*G,), prompt by prompt
"""

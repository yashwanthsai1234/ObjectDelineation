"""Synthetic GRPO batches built from the fixtures (realistic crop sizes and rollout errors)."""
from __future__ import annotations

import numpy as np

from . import fixtures, perturb


def rollout(stack: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    """A plausible policy output: a rater's mask with a random boundary error, shift or defect."""
    m = stack[int(rng.integers(len(stack)))]
    kind = rng.integers(0, 5)
    if kind == 0:
        return perturb.jitter(m, float(rng.uniform(0.5, 3.0)), rng)
    if kind == 1:
        return perturb.shift(m, int(rng.integers(-3, 4)), int(rng.integers(-3, 4)))
    if kind == 2:
        return perturb.grow(m, float(rng.uniform(-2.0, 2.0)))
    if kind == 3:
        out = perturb.add_speck(m, int(rng.integers(1, 10)), float(rng.uniform(3, 10)), rng)
        return m if out is None else out
    return m.copy()


def make_batch(seed: int, n_prompts: int = 32, group: int = 8, dataset: str | None = None):
    """(stacks, prompt_idx, masks): n_prompts fixture items x group rollouts.

    dataset=None samples prompts uniformly over all fixture items (about 40 % RIGA, 60 % LIDC);
    'riga' or 'lidc' restricts to one dataset.
    """
    rng = np.random.default_rng(seed)
    items = [st for ds, _, st in fixtures.load_all() if dataset in (None, ds)]
    pick = rng.choice(len(items), n_prompts, replace=False)
    stacks = [items[j] for j in pick]
    prompt_idx, masks = [], []
    for i, st in enumerate(stacks):
        for _ in range(group):
            prompt_idx.append(i)
            masks.append(rollout(st, rng))
    return stacks, prompt_idx, masks

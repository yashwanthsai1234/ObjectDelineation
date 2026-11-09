"""Task 6, KPI K10: identical scores bit for bit across runs, threads, processes, interpreters."""
import subprocess
import sys
from pathlib import Path

import cv2
import numpy as np
import pytest

from rewardcore.batch import BatchScorer, advantages_v4, digest
from rewardcore.reward import ExpertCache, RewardConfig, score
from rewardcore.workload import make_batch

ROOT = Path(__file__).resolve().parents[1]
CFG = RewardConfig(band_formula="F5")


@pytest.fixture(scope="module")
def batch():
    stacks, idx, masks = make_batch(seed=7, n_prompts=12, group=8)
    return stacks, [ExpertCache.build(s) for s in stacks], idx, masks


def run(caches, idx, masks, workers, backend):
    """Per-rollout terms followed by the reward-v2 group advantages."""
    with BatchScorer(caches, CFG, workers, backend) as sc:
        t = sc(idx, masks)[0]
    return np.concatenate([t.ravel(), advantages_v4(t)])


def test_repeat_threads_processes(batch):
    _, caches, idx, masks = batch
    ref = digest(run(caches, idx, masks, 1, "thread"))
    assert digest(run(caches, idx, masks, 1, "thread")) == ref
    assert digest(run(caches, idx, masks, 8, "thread")) == ref
    assert digest(run(caches, idx, masks, 1, "process")) == ref
    assert digest(run(caches, idx, masks, 8, "process")) == ref


def test_opencv_thread_count_does_not_matter(batch):
    _, caches, idx, masks = batch
    out = {}
    for n in (1, 8):
        cv2.setNumThreads(n)
        out[n] = digest(np.array([[score(m, caches[i]).band["F5"], score(m, caches[i]).depth] for i, m in zip(idx, masks)]))
    cv2.setNumThreads(1)
    assert out[1] == out[8]


def test_cache_build_is_deterministic(batch):
    stacks, caches, _, _ = batch
    for st, c in zip(stacks, caches):
        c2 = ExpertCache.build(st)
        assert (c2.tau, c2.tol4, c2.tol5, c2.perimeter) == (c.tau, c.tol4, c.tol5, c.perimeter)
        assert np.array_equal(c2.dist_out_union, c.dist_out_union) and np.array_equal(c2.dist_in_core, c.dist_in_core)


def test_fresh_interpreter_gives_same_digest(batch):
    _, caches, idx, masks = batch
    ref = digest(run(caches, idx, masks, 1, "thread"))
    code = (
        "import sys; sys.path.insert(0, %r)\n"
        "import numpy as np\n"
        "from rewardcore.batch import BatchScorer, advantages_v4, digest\n"
        "from rewardcore.reward import ExpertCache, RewardConfig\n"
        "from rewardcore.workload import make_batch\n"
        "stacks, idx, masks = make_batch(seed=7, n_prompts=12, group=8)\n"
        "with BatchScorer([ExpertCache.build(s) for s in stacks], RewardConfig(band_formula='F5'), 4, 'thread') as sc:\n"
        "    t = sc(idx, masks)[0]\n"
        "print(digest(np.concatenate([t.ravel(), advantages_v4(t)])))\n" % str(ROOT)
    )
    out = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, check=True, env={"PYTHONHASHSEED": "123"})
    assert out.stdout.strip() == ref

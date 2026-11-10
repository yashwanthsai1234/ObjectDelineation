"""E.9.4 PF2 and PF3: the verifier API and its HTTP app add almost nothing to the shared scoring function (CPU build).
The time bound of PF2 is recorded per platform in tests/goldens/perf.json; a platform without a record skips that
part. Overheads are compared on the minimum of repeated runs, which load on the machine does not lower."""
import json
import platform
import resource
import statistics
import sys
import time
from pathlib import Path

import pytest
from conftest import KEY, SIDE

from objdel.verify.api import ScoreItem, Verifier
from objdel.verify.reward.cache import ExpertCache
from objdel.verify.reward.score import score_batch

PERF = Path(__file__).resolve().parents[1] / "goldens/perf.json"
PLATFORM = f"{sys.platform}-{platform.machine()}"
GRPO_BATCH = 128  # answers per RL step on one card


def seconds(f, runs: int) -> list[float]:
    out = []
    for _ in range(runs):
        t = time.perf_counter()
        f()
        out.append(time.perf_counter() - t)
    return out


def batch(answers, n: int) -> list[ScoreItem]:
    return [ScoreItem(f"i{i}", KEY, answers[i % len(answers)], False) for i in range(n)]


def raw_jobs(view_root: Path, items: list[ScoreItem]) -> list[tuple]:
    """The score_batch jobs of `items` (as tests/verify/test_api.py builds them), outside the timed part."""
    cache = {"nodule": ExpertCache.load(view_root / "caches" / KEY / "nodule.npz")}
    return [(i.text, i.truncated, ("nodule",), cache, SIDE, SIDE) for i in items]


def max_rss_mb() -> float:
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / (2**20 if sys.platform == "darwin" else 2**10)


def test_pf2_grpo_batch_under_the_recorded_bound(view_root, answers):
    bound = json.loads(PERF.read_text()).get(PLATFORM, {}).get("verifier_128_s")
    if bound is None:
        pytest.skip(f"no verifier_128_s recorded for {PLATFORM} in {PERF.name}")
    v, items = Verifier("reward_v4", view_root), batch(answers, GRPO_BATCH)
    v.score(items)  # warm
    assert statistics.median(seconds(lambda: v.score(items), 5)) < bound * 1.25


# slow, like PF1: a ratio of wall times is not stable beside 8 xdist workers (2026-11-10, 1024 items: 1.058 alone,
# 1.11-1.16 inside make test); make test-all runs it
@pytest.mark.slow
@pytest.mark.parametrize(("n", "rel", "abs_s"), [(GRPO_BATCH, 1.05, 0.002), (1024, 1.10, 0.0)])
def test_pf2_api_overhead_over_score_batch(view_root, answers, n, rel, abs_s):
    v, items = Verifier("reward_v4", view_root), batch(answers, n)
    jobs = raw_jobs(view_root, items)
    api, raw = [], []
    for _ in range(10):  # interleaved, so both see the same machine load
        api += seconds(lambda: v.score(items), 1)
        raw += seconds(lambda: score_batch(jobs), 1)
    assert min(api) < min(raw) * rel + abs_s, (min(api), min(raw))


def test_pf2_memory_does_not_grow_over_50_calls(view_root, answers):
    v, items = Verifier("reward_v4", view_root), batch(answers, GRPO_BATCH)
    v.score(items)
    before = max_rss_mb()
    for _ in range(50):
        v.score(items)
    assert max_rss_mb() - before < 50


def test_pf3_http_overhead(view_root, answers):
    pytest.importorskip("fastapi")
    pytest.importorskip("httpx")
    from fastapi.testclient import TestClient

    from objdel.verify import http

    v = Verifier("reward_v4", view_root)
    client, items = TestClient(http.make_app(v)), batch(answers, GRPO_BATCH)
    body = {"reward": "reward_v4", "items": [{"id": i.id, "key": i.key, "text": i.text, "truncated": False}
                                             for i in items]}
    client.post("/v1/score", json=body)  # warm
    served = min(seconds(lambda: client.post("/v1/score", json=body), 5))
    local = min(seconds(lambda: v.score(items), 5))
    assert served - local < 0.050, (served, local)

    big = {"reward": "reward_v4", "items": [{"id": f"i{i}", "key": KEY, "text": "x" * 1024, "truncated": False}
                                            for i in range(http.MAX_BATCH + 1)]}
    t = time.perf_counter()
    r = client.post("/v1/score", json=big)
    assert r.status_code == 413 and time.perf_counter() - t < 0.100

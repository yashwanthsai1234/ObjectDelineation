"""The verifier API: score answers against the expert caches of a prepared view root, stamped with the reward spec.

API_VERSION: the minor number grows with additive fields only; the major number grows with a breaking change.
Scoring goes through objdel.train.rewards.score_jobs and logged_terms, the functions TRLReward calls, so the
service and the trainer have one scoring path (one CUDA checker call per batch). A model outline that does not parse
is a normal result (ScoreResult.error = the parser code); a batch the service cannot score is a VerifyError.

ScoreResult.reward is the weighted sum of the legality-gated terms, sum_k w_k term_k over (iou_tc, f5, legal, g_area)
with REWARD_V4.weights (they sum to 1): a per-answer number in [0, 1]. The group-relative advantage of the
trainer needs whole groups and is not a per-answer value, so it is not served.
"""
import dataclasses
from collections.abc import Mapping, Sequence
from importlib.metadata import version
from pathlib import Path

from objdel.geom import engine as _gpucheck
from objdel.reward import REWARD_V4, TERMS, fingerprint
from objdel.train import rewards, serve

API_VERSION = "1.0"
REWARDS = {"reward_v4": REWARD_V4}
CODES = ("unknown_key", "spec_mismatch", "batch_too_large", "device_error")


@dataclasses.dataclass(frozen=True)
class ScoreItem:
    """One answer to score: the model's outline text for the item `key` of the view root."""

    id: str
    key: str  # item id in the view root
    text: str  # the model's outline text
    truncated: bool  # generation hit the token limit (a format failure, plan D10)


@dataclasses.dataclass(frozen=True)
class ScoreResult:
    """The score of one ScoreItem: reward, the logged terms, the checker defects and the parser code."""

    id: str
    reward: float | None  # None when the outline did not parse
    terms: Mapping[str, float]
    defects: Mapping[str, int]
    error: str | None  # the parser code


@dataclasses.dataclass(frozen=True)
class Stamp:
    """Which API, reward, checker build and code version produced a result."""

    api_version: str
    reward_name: str
    reward_fingerprint: str
    checker_build: str  # "cuda" or "cpu"
    code_version: str


class VerifyError(Exception):
    """A batch the service cannot score; `code` is one of CODES."""

    def __init__(self, code: str, message: str = ""):
        assert code in CODES, code
        super().__init__(message or code)
        self.code = code


class Verifier:
    """Scores batches against the caches of one prepared view root (every dataset under it, preloaded at once). The
    CUDA engine runs on the process's current device: `device` is the index spec() reports, and the process is
    pinned to it with CUDA_VISIBLE_DEVICES."""

    def __init__(self, reward: str, view_root: Path, device: int = 0):
        if reward not in REWARDS:
            raise VerifyError("spec_mismatch", f"unknown reward {reward!r}; this service has {sorted(REWARDS)}")
        self.reward, self.device = REWARDS[reward], device
        root = Path(view_root)
        serve.preload(root, sorted({p.parent.name for p in root.glob("*/*.jsonl")}))
        self._caches = serve.loader(root)
        self.stamp = Stamp(API_VERSION, reward, fingerprint(), _gpucheck.build(), version("objdel"))

    def score(self, items: Sequence[ScoreItem]) -> list[ScoreResult]:
        """Score the batch in one checker call; results keep the order of `items`."""
        try:
            caches = [self._caches(i.key) for i in items]
        except KeyError as e:
            raise VerifyError("unknown_key", f"not in the view root: {e.args[0]}") from e
        jobs = [(i.text, i.truncated, tuple(c), c, *c[next(iter(c))].union.shape[::-1])
                for i, c in zip(items, caches, strict=True)]
        try:
            T = rewards.score_jobs(jobs)
        except RuntimeError as e:
            raise VerifyError("device_error", str(e)) from e
        terms = rewards.logged_terms(T)
        w = dict(zip(TERMS, self.reward.weights, strict=True))
        return [ScoreResult(i.id, None if T["error"][j] else float(sum(w[k] * terms[k][j] for k in w)),
                            {k: float(v[j]) for k, v in terms.items()}, T["defects"][j], T["error"][j] or None)
                for j, i in enumerate(items)]

    def spec(self) -> dict:
        """The stamp, the reward spec and the device index."""
        return {"stamp": dataclasses.asdict(self.stamp), "reward": dataclasses.asdict(self.reward),
                "device": self.device}

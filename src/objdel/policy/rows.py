"""Training data: TRL rows built from the prepared set (objdel.data.prepare).

Rows (Hugging Face datasets, conversational, one {"type": "image"} placeholder that TRL fills from the "image"
column, decoded lazily from the PNG path):
- sft_rows: messages = [user: image + PROMPTS[dataset], assistant: target text], one row per voting rater with a
  target. A rater whose whole mask is below the reward's size filter abstains from the topology vote (break-it B2-4,
  objdel.verify.reward.score.rater_answers): its own answer can never be topo-any, so SFT does not teach it (10 / 23,469
  LIDC train targets; RL keeps the items). All voting raters are rows (not one random rater per epoch): an epoch then
  shows each item's whole rater distribution, which is the topology histogram the RL reward matches; it needs no custom
  sampler. Issue I5: with tokenizer and
  max_tokens (the model's max_completion_length), targets longer than the budget (+1 end-of-turn token) are
  dropped, so SFT never teaches an answer that GRPO generation would truncate.
- prompt_rows / grpo_rows: prompt (the same user turn), image, key (item_id, the cache key), W, H, names; only
  RL-eligible items (prompt_rows(rl_only=False): every item, for probe views).
  oversample > 1 repeats disagreement items (plan §10). TRL's RepeatSampler permutes rows uniformly and takes no
  weights, so a repeated row is the simplest sampler-free weighting; the factor is decided in
  experiments/phase1/impl/M7_train.md (DISAGREEMENT_OVERSAMPLE).
Both row builders decode every image they reference and compare its pixels with the record's image_sha256
(break-it B3-5): rows are the only way images reach a trainer, so a swapped or edited PNG cannot be trained on
(ValueError). Cost: about 1 ms per image with 3 threads.
cache_loader checks every cache against the record's digest when loading it (ExpertCache.load(expected_digest));
a cache built under another reward fingerprint raises StaleCacheError ("re-run prepare_rl"), distinct from a
corrupted or swapped file (break-it B2-8).

Flip augmentation (issue I7) is not implemented: Phase 1 runs no training. If Phase 2 adds it, flip the image and
the rater masks and rebuild targets (build_target) and caches (ExpertCache.build) from the flipped masks; never
mirror polygon vertices, because the top-left pixel-centre tie rule is not mirror-symmetric.

References:
- TRL dataset formats (conversational, vision: "image" column + {"type": "image"}):
https://huggingface.co/docs/trl/main/en/dataset_formats
- TRL GRPO (image column, extra columns reach reward functions): https://huggingface.co/docs/trl/main/en/grpo_trainer
- datasets Image feature (paths decoded on access): https://huggingface.co/docs/datasets/main/en/image_load
- OpenCV Python bindings release the GIL around every call (PyAllowThreads, modules/python/src2/cv2_util.hpp):
  https://github.com/opencv/opencv/blob/4.x/modules/python/src2/cv2_util.hpp
- Plan §4.3 (caches stored with a hash), §4.4 (targets), §10 (oversample disagreement items).
"""

from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import cv2
import numpy as np

from objdel.core.hashing import pixels_sha256
from objdel.core.io import read_jsonl
from objdel.geom.text import SCHEMAS
from objdel.policy.prompt import _user

DISAGREEMENT_OVERSAMPLE = 4  # decided in experiments/phase1/impl/M7_train.md


def records(root, dataset: str, split: str) -> list[dict]:
    """The rows of `<root>/<dataset>/<split>.jsonl` as dicts."""
    return read_jsonl(Path(root, dataset, f"{split}.jsonl"))


def _image(root, item_id: str) -> str:
    return str(Path(root, "images", item_id + ".png").resolve())


def _check(root, rec: dict) -> None:
    im = cv2.imread(_image(root, rec["item_id"]), cv2.IMREAD_COLOR)
    if im is None or pixels_sha256(cv2.cvtColor(im, cv2.COLOR_BGR2RGB)) != rec["image_sha256"]:
        raise ValueError(f"{rec['item_id']}: prepared image missing or not the one recorded (image_sha256)")


def _verified(root, recs: list[dict]) -> list[dict]:
    """recs, after checking every image's pixels against its record (cv2 decodes without the GIL: 3 threads)."""
    with ThreadPoolExecutor(3) as pool:
        list(pool.map(lambda r: _check(root, r), recs))
    return recs


def target_tokens(tokenizer, text: str) -> int:
    """Token count of the text under the tokenizer, without special tokens."""
    return len(tokenizer(text, add_special_tokens=False)["input_ids"])


def prompts(root, dataset: str, split: str, rl_only: bool = True, oversample: int = 1, subset=None) -> list[dict]:
    """Prompt rows (prompt, image path, key, W, H, names) of the records (images checked against them): RL-eligible
    items only, or with rl_only=False every item (probe views, whose records are rl: false); disagreement items
    repeated `oversample` times; subset = (n, seed): n records drawn once with that seed, in record order. Window
    records of a rows root (build_rl_rows: `cell` and view.kind) add the columns kind and cell (TRLReward logs per
    kind and cell); other records have no such columns."""
    recs = [r for r in records(root, dataset, split) if r["rl"] or not rl_only]
    if subset and subset[0] < len(recs):
        n, seed = subset
        recs = [recs[i] for i in sorted(np.random.default_rng(seed).choice(len(recs), n, replace=False))]
    return [{"prompt": [_user(dataset)], "image": _image(root, r["item_id"]), "key": r["item_id"], "W": r["W"],
             "H": r["H"], "names": list(SCHEMAS[dataset]),
             **({"kind": r["view"]["kind"], "cell": r["cell"]} if "cell" in r else {})}
            for r in _verified(root, recs) for _ in range(oversample if r["disagreement"] else 1)]

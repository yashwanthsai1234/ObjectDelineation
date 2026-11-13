"""Training data: TRL rows built from the prepared set (objdel.data.prepare).

Rows are Hugging Face datasets in conversational form with one {"type": "image"} placeholder that TRL fills from the
"image" column (decoded lazily from the PNG path):
- sft_rows: messages = [user: image + the dataset's prompt, assistant: target text], one row per voting rater with a
  target. A rater whose whole mask is below the reward's size filter abstains from the topology vote, so SFT does not
  teach its answer (RL keeps the item). Every voting rater is a row: an epoch shows each item's whole rater
  distribution, the topology histogram the RL reward matches, and needs no custom sampler. With tokenizer and
  max_tokens, targets longer than the budget (+1 end-of-turn token) are dropped, so SFT never teaches an answer that
  GRPO generation would truncate.
- prompt_rows: prompt (the same user turn), image, key (item_id, the cache key), W, H, names; only RL-eligible items
  (rl_only=False: every item, for probe views). oversample > 1 repeats disagreement items: TRL's RepeatSampler takes
  no weights, so a repeated row is the simplest weighting (DISAGREEMENT_OVERSAMPLE).
Both builders decode every image they reference and compare its pixels with the record's image_sha256 (ValueError):
rows are the only way images reach a trainer, so a swapped or edited PNG cannot be trained on. cache_loader checks
every cache against the record's digest; a cache built under another reward fingerprint raises StaleCacheError
("re-run prepare_rl"), distinct from a corrupted file.

Flip augmentation is not implemented. If added, flip the image and the rater masks and rebuild targets and caches from
the flipped masks; never mirror polygon vertices, because the top-left pixel-centre tie rule is not mirror-symmetric.

References: TRL dataset formats (https://huggingface.co/docs/trl/main/en/dataset_formats), TRL GRPO (extra columns
reach reward functions: https://huggingface.co/docs/trl/main/en/grpo_trainer).
"""

from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

from objdel.core import datasets
from objdel.core.hashing import pixels_sha256  # noqa: F401 (tests read policy.rows.pixels_sha256)
from objdel.core.io import read_json, read_jsonl
from objdel.data.prepare import image_matches
from objdel.policy.prompt import _user

DISAGREEMENT_OVERSAMPLE = 4  # decided in experiments/phase1/impl/M7_train.md


def check_root(root, dataset: str, split: str) -> None:
    """ValueError unless the folder's manifest lists the set `<dataset>/<split>` (a stage's dataset against its data
    folder)."""
    sets = read_json(Path(root, "manifest.json"))["sets"]
    if f"{dataset}/{split}" not in sets:
        raise ValueError(f"{root}: holds {sorted(sets)}, not {dataset}/{split} (the stage's dataset and split)")


def records(root, dataset: str, split: str) -> list[dict]:
    """The rows of `<root>/<dataset>/<split>.jsonl` as dicts."""
    return read_jsonl(Path(root, dataset, f"{split}.jsonl"))


def image_path(root, item_id: str) -> str:
    """The resolved path of the prepared image of `item_id`: <root>/images/<item_id>.png."""
    return str(Path(root, "images", item_id + ".png").resolve())


_image = image_path  # the pinned old name (tests/goldens/phase_e/old_names.json)


def rgb(path):
    """The image at `path` as a PIL RGB image (Pillow decode; the processor's input)."""
    from PIL import Image

    return Image.open(path).convert("RGB")


def _check(root, rec: dict) -> None:
    if not image_matches(image_path(root, rec["item_id"]), rec["image_sha256"]):
        raise ValueError(f"{rec['item_id']}: prepared image missing or not the one recorded (image_sha256)")


def _verified(root, recs: list[dict]) -> list[dict]:
    """recs, after checking every image's pixels against its record (cv2 decodes without the GIL: 3 threads)."""
    with ThreadPoolExecutor(3) as pool:
        list(pool.map(lambda r: _check(root, r), recs))
    return recs


def draw_subset(seq, n: int, seed: int) -> list:
    """n elements of seq drawn once without replacement with numpy's default_rng(seed), kept in seq order."""
    return [seq[i] for i in sorted(np.random.default_rng(seed).choice(len(seq), n, replace=False))]


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
        recs = draw_subset(recs, *subset)
    return [{"prompt": [_user(dataset)], "image": _image(root, r["item_id"]), "key": r["item_id"], "W": r["W"],
             "H": r["H"], "names": list(datasets.get(dataset).names),
             **({"kind": r["view"]["kind"], "cell": r["cell"]} if "cell" in r else {})}
            for r in _verified(root, recs) for _ in range(oversample if r["disagreement"] else 1)]

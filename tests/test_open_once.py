"""Open-once GRPO rows (dataset.OpenOnce): each step opens each of its 32 images once, not 16 x 256 times.

TRL's RepeatSampler hands the DataLoader the same 256 indices (32 images x 8 answers) once per micro-batch, 16 times a
step, and transformers fetches all 16 before the step starts; a datasets.Dataset opens a PNG for every row of every
fetch. The loader here uses the trainer's settings (configs.grpo_config: 8 answers per image, 16 x 16, so TRL's sampler
gets 32 prompts per generation batch and a repeat count of 16, and the DataLoader 256 rows per fetch)."""

import numpy as np
import PIL.Image
import pytest
import torch
from datasets import Dataset, Image
from trl.trainer.utils import RepeatSampler

from objdel.train.dataset import OpenOnce

G, PROMPTS, GA = 8, 32, 16


@pytest.fixture
def rows(tmp_path):
    rng = np.random.default_rng(0)
    out = []
    for i in range(2 * PROMPTS):  # two steps' worth
        p = tmp_path / f"{i}.png"
        PIL.Image.fromarray(rng.integers(0, 256, (24, 20, 3), dtype=np.uint8)).save(p)
        out.append({"prompt": [{"role": "user", "content": f"q{i}"}], "image": str(p), "key": f"k{i}", "W": 20,
                    "H": 24})
    return Dataset.from_list(out).cast_column("image", Image())


def loader(data):
    sampler = RepeatSampler(data, mini_repeat_count=G, batch_size=PROMPTS, repeat_count=GA, seed=0)
    return torch.utils.data.DataLoader(data, batch_size=G * PROMPTS, sampler=sampler, collate_fn=lambda b: b)


@pytest.fixture
def opened(monkeypatch):
    calls, real = [], PIL.Image.open
    monkeypatch.setattr(PIL.Image, "open", lambda f, *a, **k: (calls.append(str(f)), real(f, *a, **k))[1])
    return calls


def test_open_once_rows_equal_the_decoded_rows_fetch_for_fetch(rows):
    for a, b in zip(loader(rows), loader(OpenOnce(rows)), strict=True):
        assert len(a) == len(b) == G * PROMPTS
        for x, y in zip(a, b, strict=True):
            assert x.keys() == y.keys()
            assert (x["prompt"], x["key"], x["W"], x["H"]) == (y["prompt"], y["key"], y["W"], y["H"])
            assert (x["image"].mode, x["image"].size) == (y["image"].mode, y["image"].size)
            assert np.array_equal(np.asarray(x["image"]), np.asarray(y["image"]))


def test_a_step_opens_each_of_its_32_images_once_where_the_plain_rows_open_4096(rows, opened):
    plain = iter(loader(rows))
    for _ in range(GA):
        next(plain)
    assert len(opened) == GA * G * PROMPTS == 4096  # the waste: every row of every fetch
    opened.clear()

    it = iter(loader(OpenOnce(rows)))
    step1 = [next(it) for _ in range(GA)]
    assert len(opened) == len(set(opened)) == PROMPTS
    assert all(b is step1[0] for b in step1)  # one stack, handed out for all 16 micro-batches
    step2 = [next(it) for _ in range(GA)]
    assert len(opened) == len(set(opened)) == 2 * PROMPTS  # the next step opens its own 32
    assert {r["key"] for r in step2[0]}.isdisjoint({r["key"] for r in step1[0]})


def test_the_8_rows_of_an_image_share_one_opened_image(rows):
    batch = next(iter(loader(OpenOnce(rows))))
    by_key = {}
    for r in batch:
        by_key.setdefault(r["key"], set()).add(id(r["image"]))
    assert len(by_key) == PROMPTS and all(len(ids) == 1 for ids in by_key.values())


def test_single_rows_still_read(rows):
    data = OpenOnce(rows)
    assert len(data) == len(rows)
    assert np.array_equal(np.asarray(data[3]["image"]), np.asarray(rows[3]["image"]))

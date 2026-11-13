"""D11 single sources (MT-12, MT-13, MT-19): registry values instead of literals, one RL stage refusal, one miss
line for train and eval."""
from pathlib import Path

import pytest

from objdel.core import datasets


def test_baseline_margins_cover_every_library_dataset():
    from objdel.cli.eval import baselines

    assert set(baselines.MARGIN) >= set(datasets.DATASETS)  # the tests' riga included


def test_box_data_is_one_template():
    from objdel.infra.vast import chain
    from objdel.sft import train as sft_train

    assert chain.BOX_DATA == datasets.BOX_DATA == "/workspace/data/{dataset}"
    assert sft_train.BOX_DATA == Path("/workspace/data/lidc")


def test_ceiling_reads_the_dataset_and_its_structure_from_the_registry(tmp_path):
    from objdel.cli.eval import ceiling

    for d in ("lidc", "riga"):  # riga (disc + cup) is not a ceiling dataset
        (tmp_path / d).mkdir()
        (tmp_path / d / "val.jsonl").write_text("")
    assert ceiling.view_dataset(tmp_path, "val", None) == ("lidc", "nodule")
    with pytest.raises(ValueError, match="pass --dataset"):
        ceiling.view_dataset(tmp_path, "test", None)
    with pytest.raises(ValueError, match="one structure"):
        ceiling.view_dataset(tmp_path, "val", "riga")  # disc + cup


def test_one_stage_refusal_for_train_rl_and_the_probe():
    from objdel.rl import run as rl_run

    assert rl_run.stage_refusal("qwen3-vl-2b", "lidc", None) is None
    assert rl_run.stage_refusal("qwen3-vl-2b", "lidc", 28) is None
    assert "28 text decoder layers" in rl_run.stage_refusal("qwen3-vl-2b", "lidc", 29)
    assert rl_run.stage_refusal("qwen3-vl-2b", "nope", None) == rl_run.rl_refusal("nope")


def test_the_train_miss_line_is_the_eval_miss_line():
    from objdel.rl import reward
    from objdel.verify import stats

    assert reward.MISS_IOU == stats.MISS == 0.1

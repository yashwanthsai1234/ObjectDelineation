"""The stage field `dataset`: default lidc (today's argv unchanged), read by RL, the probe, select and final eval;
policy.rows.check_root holds a stage's dataset against its data folder; RL refuses a dataset whose items have fewer
raters than the reward needs (no single-rater reward mode)."""
import dataclasses
import json

import pytest

from objdel.core import datasets, spec
from objdel.policy import rows as policy_rows

BASE = ["--model", "qwen3-vl-2b", "--ckpt", "c", "--train", "t", "--out", "o", "--steps", "20"]


def test_the_four_stages_default_to_lidc_and_emit_dataset_only_when_set():
    stages = {spec.RlStage: dict(model="m", ckpt="c", train="t", out="o", steps=2),
              spec.ProbeStage: dict(model="m", ckpt="c", train="t", probe="p", out="o"),
              spec.SelectStage: dict(run="r", shifted="s", no_vision=True),
              spec.EvalStage: dict(ckpt="c", view=("v",), out=("o",))}
    for cls, kw in stages.items():
        st = cls(**kw)
        assert st.dataset == "lidc" and "--dataset" not in spec.to_argv(st, "rd", "data"), cls.__name__
        assert spec.to_argv(cls(**kw, dataset="other"), "rd", "data")[-2:] == ["--dataset", "other"], cls.__name__
    assert "dataset" not in spec.SftStage.model_fields  # SFT reads its view's dataset (the view manifest)


def test_check_root_holds_the_dataset_against_the_folder(tmp_path):
    (tmp_path / "manifest.json").write_text(json.dumps({"sets": {"lidc/train": {}, "lidc/val": {}}}))
    policy_rows.check_root(tmp_path, "lidc", "val")
    with pytest.raises(ValueError, match=r"not riga/val"):
        policy_rows.check_root(tmp_path, "riga", "val")
    with pytest.raises(ValueError, match=r"not lidc/test"):
        policy_rows.check_root(tmp_path, "lidc", "test")


def test_rl_refuses_unknown_and_single_rater_datasets(monkeypatch, capsys):
    from objdel.rl import run as rl_run

    def parse(argv):  # train_rl.parse
        return rl_run.checked(spec.cli_stage(spec.RlStage, argv))

    assert rl_run.rl_refusal("lidc") is None and parse(BASE).dataset == "lidc"
    one = dataclasses.replace(datasets.LIDC, name="polyp", min_raters=1)
    monkeypatch.setitem(datasets.DATASETS, one.name, one)
    for name, why in (("polyp", "1 rater(s) per item; the RL reward needs 2"), ("nope", "not one of")):
        assert why in rl_run.rl_refusal(name)
        with pytest.raises(SystemExit) as e:
            parse([*BASE, "--dataset", name])
        assert e.value.code == 2 and why in capsys.readouterr().err


def test_a_spec_reads_one_dataset_and_its_box_folder():
    from objdel.core.provenance import repo_root
    from objdel.infra.vast import chain

    run = spec.load_spec(repo_root() / "specs/b1_rows.yaml")
    assert spec.spec_dataset(run) == "lidc" and chain.box_data(run) == "/workspace/data/lidc"
    rl = next(s for s in run.stages if s.kind == "rl")
    two = run.model_copy(update={"stages": (rl, rl.model_copy(update={"dataset": "x", "name": "b"}))})
    with pytest.raises(ValueError, match="a run reads one dataset"):
        spec.spec_dataset(two)

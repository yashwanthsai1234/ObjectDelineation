"""C0 flags --record-batches DIR and --checks-from-step S (U-CF1, U-CF2): the check and recording steps of a
FastGRPO run, and the file format that G-EQ2 and G-EQ3 read."""

import collections

import pytest

torch = pytest.importorskip("torch")
pytest.importorskip("trl")
from test_row_meta import EOS, device_masks, generation_trainer, rows  # noqa: E402

from objdel.policy import caches as policy_caches  # noqa: E402
from objdel.rl import record  # noqa: E402
from objdel.rl import trl_adapter as rl_trl_adapter  # noqa: E402
from tests.experiments.test_trl_adapter_equivalence import group_batch, tiny_model  # noqa: E402

BASE = ["--model", "m", "--ckpt", "c", "--train", "t", "--out", "o", "--steps", "20"]


def test_flags_default_to_the_control_and_validate():
    from experiments.phase2 import train_rl

    a = train_rl.parse(BASE)
    assert (a.record_batches, a.checks_from_step) == (None, 1)
    a = train_rl.parse([*BASE, "--record-batches", "d", "--checks-from-step", "18"])
    assert (str(a.record_batches), a.checks_from_step) == ("d", 18)
    with pytest.raises(SystemExit):
        train_rl.parse([*BASE, "--checks-from-step", "0"])


@pytest.mark.parametrize(("every", "S", "step", "due"), [
    (1, 1, 1, True), (50, 1, 51, True), (50, 1, 2, False), (0, 1, 1, False),  # S = 1: steps 1, N + 1, ... as before
    (1, 18, 17, False), (1, 18, 18, True), (1, 18, 32, True),  # the check arm: every step from 18
    (50, 18, 18, True), (50, 18, 19, False), (50, 18, 68, True), (0, 18, 18, False)])
def test_checks_act_from_step_s_on(monkeypatch, every, S, step, due):
    """U-CF1: with S the check flags count their steps from S (the step is global_step + 1); spy on both checks."""
    t, inputs, _ = generation_trainer(monkeypatch, every, step - 1)
    t.check_rope_every, t.checks_from_step = every, S
    calls = []
    monkeypatch.setattr(rl_trl_adapter, "check_meta", lambda batch: calls.append("meta"))
    t._generate_and_score_completions(inputs)
    assert (calls == ["meta"], t.rope_check) == (due, due)


def recording_trainer(monkeypatch, tmp_path, step, S=18, training=True):
    """generation_trainer whose batch holds what TRL's builder returns (ids, masks, grids)."""
    t, inputs, G = generation_trainer(monkeypatch, 0, step - 1)
    prompt, comp, pm, P, T, _ = rows(5)
    cm = device_masks(prompt, comp, pm, T, True)
    n = len(prompt)

    def build(self, inputs):
        self._generate(prompt)
        return {"prompt_ids": torch.randint(10, 99, (n, P)), "prompt_mask": pm,
                "completion_ids": torch.zeros(n, T).long(), "completion_mask": cm,
                "image_grid_thw": torch.tensor([[1, 4, 4]] * n), "num_images": [1] * n,
                "advantages": torch.zeros(n), "pixel_values": torch.zeros(n, 3)}

    monkeypatch.setattr(rl_trl_adapter.GRPOTrainer, "_generate_and_score_completions", build)
    t.record_dir, t.checks_from_step = record.open_dir(tmp_path / "rec"), S
    t.model.train(training)
    return t, inputs


def test_a_batch_is_recorded_from_step_s_on_and_loads_with_weights_only(monkeypatch, tmp_path):
    """U-CF2: one file per recorded step; ids, masks, grids and RowMeta, no pixels; weights_only loads it."""
    t, inputs = recording_trainer(monkeypatch, tmp_path, step=17)
    t._generate_and_score_completions(inputs)
    assert list((tmp_path / "rec").iterdir()) == []  # step 17 < 18
    t.state.global_step = 17  # step 18
    batch = t._generate_and_score_completions(inputs)
    t.state.global_step = 18
    t._generate_and_score_completions(inputs)
    assert sorted(f.name for f in (tmp_path / "rec").iterdir()) == ["step0018.pt", "step0019.pt"]
    rec = torch.load(tmp_path / "rec/step0018.pt", weights_only=True)
    assert set(rec) == {"prompt_ids", "prompt_mask", "completion_ids", "completion_mask", "image_grid_thw",
                        "num_images", "meta"}
    assert torch.equal(rec["prompt_ids"], batch["prompt_ids"]) and rec["num_images"].tolist() == [1] * 16
    assert rec["meta"]["prompt_len"].tolist() == [r.prompt_len for r in batch["host_meta"]]
    assert rec["meta"]["completion_len"].tolist() == [r.completion_len for r in batch["host_meta"]]
    assert set(rec["meta"]) == {"prompt_len", "completion_len", "image_id", "group", "layout", "rope_bad", "pixel_rows"}


def test_an_eval_batch_is_not_recorded(monkeypatch, tmp_path):
    t, inputs = recording_trainer(monkeypatch, tmp_path, step=20, training=False)
    t._metrics["eval"] = collections.defaultdict(list)
    t._generate_and_score_completions(inputs)
    assert list((tmp_path / "rec").iterdir()) == []


def test_record_works_without_host_meta_and_leaves_the_batch_alone(monkeypatch, tmp_path):
    t, inputs = recording_trainer(monkeypatch, tmp_path, step=20)
    t.host_meta = False
    batch = t._generate_and_score_completions(inputs)
    assert "host_meta" not in batch and (tmp_path / "rec/step0020.pt").is_file()
    assert "meta" in torch.load(tmp_path / "rec/step0020.pt", weights_only=True)


def test_a_non_empty_directory_raises_and_nothing_is_written(tmp_path):
    d = tmp_path / "rec"
    assert record.open_dir(d) == d and d.is_dir()
    assert record.open_dir(d) == d  # an empty folder is fine
    (d / "old.pt").write_bytes(b"x")
    with pytest.raises(FileExistsError, match="not empty"):
        record.open_dir(d)
    assert [f.name for f in d.iterdir()] == ["old.pt"]


def test_positions_are_hf_get_rope_index_on_the_full_rows():
    """The recorded positions (3, rows, width) are HF's own on prompt + answer, masked holes included."""
    m = tiny_model()
    ids, mask, T, kw = group_batch(groups=2, G=4, lengths=[3, 7, 1, 5, 2, 6, 4, 10], pads=[2, 0])
    P = ids.size(1) - T
    batch = {"prompt_ids": ids[:, :P], "prompt_mask": mask[:, :P], "completion_ids": ids[:, P:],
             "completion_mask": mask[:, P:], "mm_token_type_ids": kw["mm_token_type_ids"],
             "image_grid_thw": kw["image_grid_thw"], "num_images": kw["num_images"]}
    rec = record.batch_record(batch, None, policy_caches.RopeCache.install(m).base)
    want, _ = m.model.get_rope_index(ids, mm_token_type_ids=kw["mm_token_type_ids"],
                                     image_grid_thw=kw["image_grid_thw"], attention_mask=mask)
    assert torch.equal(rec["positions"], want) and "meta" not in rec
    assert EOS  # the shared row helpers are imported for the fixtures above

"""objdel.policy.models.save_merged (MT-02, IO-08): the one writer of a merged model folder (RL's final, SFT's
branches), with fakes in place of the model and processor."""
from objdel.policy import models


class _Saves:
    def __init__(self, name, log):
        self.name, self.log = name, log

    def save_pretrained(self, d):
        d.mkdir(parents=True, exist_ok=True)
        (d / f"{self.name}.bin").write_text(self.name)
        self.log.append((self.name, d.name))


class _Lora(_Saves):
    def merge_and_unload(self):
        self.log.append(("merge", None))
        return _Saves("merged", self.log)


def test_save_merged_writes_adapter_then_merged_model_through_the_named_temp(tmp_path):
    log, out = [], tmp_path / "out"
    (out / "final").mkdir(parents=True)
    (out / "final" / "half.bin").write_text("an older final")
    models.save_merged(_Lora("adapter", log), _Saves("proc", log), out / "final", adapter_to=out / "final_adapter",
                       tmp=out / "final.tmp")
    assert log == [("adapter", "final_adapter"), ("merge", None), ("merged", "final.tmp"), ("proc", "final.tmp")]
    assert sorted(p.name for p in (out / "final").iterdir()) == ["merged.bin", "proc.bin"]
    assert sorted(p.name for p in out.iterdir()) == ["final", "final_adapter"]


def test_save_merged_of_a_plain_model_and_default_temp(tmp_path):
    log = []
    models.save_merged(_Saves("plain", log), _Saves("proc", log), tmp_path / "branches" / "e1")
    assert log == [("plain", ".e1.tmp"), ("proc", ".e1.tmp")]
    assert sorted(p.name for p in (tmp_path / "branches" / "e1").iterdir()) == ["plain.bin", "proc.bin"]


def test_rl_save_final_keeps_its_names(tmp_path, monkeypatch):
    from objdel.rl import run as rl_run

    seen = {}
    monkeypatch.setattr(models, "save_merged", lambda m, p, dest, **k: seen.update(dest=dest, **k))
    rl_run.save_final(None, None, tmp_path)
    assert seen == {"dest": tmp_path / "final", "adapter_to": tmp_path / "final_adapter", "tmp": tmp_path / "final.tmp"}


def test_rl_last_checkpoint_takes_the_highest_step_without_a_completeness_check(tmp_path):
    from objdel.rl import run as rl_run

    assert rl_run.last_checkpoint(tmp_path) == (None, 0)
    for n in (5, 40, 250):
        (tmp_path / "trainer" / f"checkpoint-{n}").mkdir(parents=True)
    (tmp_path / "trainer" / "checkpoint-999.txt").write_text("not a folder")
    assert rl_run.last_checkpoint(tmp_path) == (str(tmp_path / "trainer" / "checkpoint-250"), 250)

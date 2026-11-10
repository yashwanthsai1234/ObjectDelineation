"""D8's one GPU skip rule (tests/conftest.py): a `gpu` test skips unless the extension is the CUDA build, whatever
nvidia-smi says; no test file keeps its own nvidia-smi skipif."""
from pathlib import Path

import pytest

from objdel.geom import engine
from tests import conftest


class Item:
    def __init__(self, *marks):
        self.marks = list(marks)

    def get_closest_marker(self, name):
        return name if name in self.marks else None

    def add_marker(self, marker):
        self.marks.append(marker)


@pytest.mark.parametrize(("build", "skipped"), [("cpu", True), ("cuda", False)])
def test_gpu_tests_skip_on_the_cpu_build_only(monkeypatch, build, skipped):
    monkeypatch.setattr(engine, "build", lambda: build)
    gpu, other = Item("gpu"), Item("slow")
    conftest.pytest_collection_modifyitems(None, [gpu, other])
    assert len(gpu.marks) == 1 + skipped and other.marks == ["slow"]


def test_no_test_file_skips_on_nvidia_smi():
    offenders = [p.name for p in Path(__file__).parents[1].rglob("test_*.py")
                 if p != Path(__file__) and "nvidia-smi -L)" in p.read_text()]
    assert offenders == []

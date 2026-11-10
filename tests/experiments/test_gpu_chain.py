"""The GPU check chain on the Mac at $0 (Implementation.md D8): `make gpu-check DRY=1` prints every box step in order
and collects D_GPU_TESTS; each entry of the list collects a test, and the sanitizer selection is the corpus. Reads no
data/. The JUnit skip check fails the box on any skip; d_gpu_test.yaml calls the chain and A's two box tests."""
import re
import subprocess
from pathlib import Path

from experiments.cuda_checker.junit_skips import skipped
from objdel.infra.vast import chain as C

ROOT = Path(__file__).resolve().parents[2]
BOX = ROOT / "experiments/cuda_checker/run_box.sh"
PY = ROOT / ".venv/bin/python"


def _var(name: str) -> str:
    return re.search(rf'^{name}="([^"]*)"', BOX.read_text().replace("\\\n", " "), re.M).group(1)


def _collect(ids: list[str], sel: str) -> list[str]:
    out = subprocess.run([PY, "-I", "-m", "pytest", "--collect-only", "-q", "-p", "no:cacheprovider", *ids, "-m", sel],
                         cwd=ROOT, capture_output=True, text=True)
    assert out.returncode == 0, out.stdout[-2000:] + out.stderr[-2000:]
    return [line for line in out.stdout.splitlines() if "::" in line]


def test_dry_chain_prints_the_box_steps_in_order_and_collects_the_list():
    out = subprocess.run(["make", "gpu-check", "DRY=1"], cwd=ROOT, capture_output=True, text=True)
    assert out.returncode == 0, out.stdout[-2000:] + out.stderr[-2000:]
    steps = [line[2:] for line in out.stdout.splitlines() if line.startswith("+ ")]
    want = ["uv pip install", "CUDA_HOME=", '[ -e "$CUDA_HOME/lib64" ]', "ln -sf libcudart", "export CUDACXX=",
            "make cuda", '"$PY" -I -c "import objdel._gpucheck', '[ -f "$OBJDEL_DATA_ROOT/checker_eval/masks.npz" ]',
            '"$PY" -I -m pytest -q -p no:cacheprovider $D_GPU_TESTS -m "$SEL" --junitxml "$OUT/gpu.xml"',
            '"$PY" -I -m experiments.cuda_checker.junit_skips "$OUT/gpu.xml"',
            *(f"--tool {t} --error-exitcode 99" for t in ("memcheck", "initcheck", "synccheck", "racecheck"))]
    assert len(steps) == len(want) and all(w in s for w, s in zip(want, steps, strict=True)), steps
    assert not any("--leak-check" in s for s in steps)
    assert sum("::" in line for line in out.stdout.splitlines()) > 100  # the collect-only listing


def test_every_listed_entry_collects_a_test_and_the_sanitizer_corpus_is_selected():
    ids, sel = _var("D_GPU_TESTS").split(), _var("SEL")
    got = _collect(ids, sel)
    assert [i for i in ids if not any(g.startswith(i) for g in got)] == []
    for node in ("tests/geom/test_cuda_checker.py::test_validity_adversarial_families",
                 "tests/experiments/test_scorer_sums.py::test_g01_golden", "tests/geom/test_engine_counters.py::"):
        assert any(g.startswith(node) for g in got), node
    archive = ("tests/geom/test_cuda_checker.py::test_slot_polygonizer_is_archived",
               "tests/geom/test_cuda_checker.py::test_shape_totals_is_archived")
    assert not set(archive) & set(got)  # git-index tests: the box's source has no .git or archive/ (they would skip)
    assert set(archive) <= set(_collect(["tests/geom/test_cuda_checker.py"], "not slow"))  # make test still runs them
    san = _collect(ids, "gpu and sanitizer")
    assert {g.split("[")[0] for g in san} == {"tests/geom/test_gpu_sanitize.py::test_san",
                                             "tests/geom/test_gpu_sanitize.py::test_san_no_leak",
                                             "tests/geom/test_gpu_sanitize.py::test_san_pool",
                                             "tests/experiments/test_scorer_sums.py::test_g07_sanitize_sums"}
    assert len(san) == 11


def test_junit_skips_lists_each_skip(tmp_path):
    x = tmp_path / "gpu.xml"
    x.write_text('<testsuites><testsuite skipped="1"><testcase classname="tests.a" name="t1"/>'
                 '<testcase classname="tests.a" name="t2"><skipped message="no input"/></testcase>'
                 '</testsuite></testsuites>')
    assert skipped(str(x)) == ["tests.a::t2: no input"]


def test_d_gpu_test_runs_the_chain_then_a_box_tests():
    s = C.chain(ROOT / "specs/d_gpu_test.yaml", "r")
    assert s[0] == ["bash", "-c", "OBJDEL_DATA_ROOT=/workspace/data/lidc PY=/opt/venv/bin/python "
                                  "bash experiments/cuda_checker/run_box.sh --gpu-check"]
    ids = [t for t in s[1][2].split() if "::" in t]
    assert ids == ["tests/experiments/test_a_gpu.py::test_g02b_builds",
                   "tests/experiments/test_a_gpu.py::test_g03_device_bytes"]
    assert len(_collect(ids, "gpu")) == 2


def test_cuda_host_code_does_not_contract_fmas():
    """docs/research/notes cuda-kernels sections 2 and 10: nvcc's host compiler contracts FMAs by default."""
    cmake = (ROOT / "CMakeLists.txt").read_text()
    assert "$<$<COMPILE_LANGUAGE:CUDA>:--fmad=false -Xcompiler=-ffp-contract=off>" in cmake

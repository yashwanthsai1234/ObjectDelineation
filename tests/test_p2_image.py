"""Phase 2 image size and causal-conv1d (V0, 2026-11-01): the first push held two copies of the CUDA wheels both venvs
pin (3.94 + 3.40 GB compressed), and causal-conv1d fell back silently because the image has no nvcc."""

import importlib.util
import os
import re
import shutil
import tomllib
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
DOCKER = ROOT / "experiments" / "phase2" / "docker"
LOCKS = DOCKER / "locks"


def load(name):
    spec = importlib.util.spec_from_file_location(name, DOCKER / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_dedupe_hardlinks_identical_files_across_venvs(tmp_path):
    dedupe = load("dedupe")
    a, b = tmp_path / "vllm", tmp_path / "venv"
    for root in (a, b):
        (root / "nvidia/cublas/lib").mkdir(parents=True)
        (root / "nvidia/cublas/lib/libcublas.so.13").write_bytes(b"same" * 300_000)
    (a / "torch.so").write_bytes(b"2.13" * 300_000)  # same size, different bytes
    (b / "torch.so").write_bytes(b"2.14" * 300_000)
    (a / "small.py").write_text("x")
    (b / "small.py").write_text("x")
    (b / "link.so").symlink_to(b / "torch.so")
    saved = dedupe.dedupe([a, b], min_size=1000)
    assert saved == 1_200_000
    assert os.path.samefile(a / "nvidia/cublas/lib/libcublas.so.13", b / "nvidia/cublas/lib/libcublas.so.13")
    assert not os.path.samefile(a / "torch.so", b / "torch.so") and (b / "torch.so").read_bytes()[:4] == b"2.14"
    assert not os.path.samefile(a / "small.py", b / "small.py")  # below min_size
    assert (b / "link.so").is_symlink()
    assert dedupe.dedupe([a, b], min_size=1000) == 0  # idempotent


def test_dedupe_keeps_different_modes_apart(tmp_path):
    dedupe = load("dedupe")
    (tmp_path / "a").mkdir(), (tmp_path / "b").mkdir()
    (tmp_path / "a/x").write_bytes(b"1" * 5000), (tmp_path / "b/x").write_bytes(b"1" * 5000)
    os.chmod(tmp_path / "b/x", 0o755)
    assert dedupe.dedupe([tmp_path / "a", tmp_path / "b"], min_size=1000) == 0


def test_patch_arch_keeps_only_our_cards(tmp_path):
    patch = load("patch_arch")
    src = (
        "def append_nvcc_threads(nvcc_extra_args):\n"
        "    return nvcc_extra_args + [\"--threads\", \"4\"]\n"
        "    else:\n"
        "        cc_flag.append(\"-gencode\")\n"
        "        cc_flag.append(\"arch=compute_75,code=sm_75\")\n"
        "        cc_flag.append(\"-gencode\")\n"
        "        cc_flag.append(\"arch=compute_80,code=sm_80\")\n"
        "        if bare_metal_version >= Version(\"13.0\"):\n"
        "            cc_flag.append(\"arch=compute_121,code=sm_121\")\n"
        "\n"
        "    # HACK: The compiler flag -D_GLIBCXX_USE_CXX11_ABI is set\n"
    )
    out = patch.patch(src)
    assert "compute_75" not in out and "compute_121" not in out and "bare_metal_version >=" not in out
    assert re.findall(r"arch=compute_\d+,code=\w+", out) == [
        "arch=compute_80,code=sm_80", "arch=compute_120,code=sm_120", "arch=compute_120,code=compute_120"]
    assert '["--threads", "1"]' in out and '"--threads", "4"' not in out  # V0: --threads 4 OOM-killed bwd at 4 GB
    with pytest.raises(ValueError):
        patch.patch("no gencode block here")


def run_steps(text):
    """RUN instructions with their line continuations joined."""
    return [re.sub(r"\\\n\s*", " ", r) for r in re.findall(r"^RUN .*?(?<!\\)$", text, re.M | re.S)]


def test_dockerfile_builds_both_venvs_in_one_layer_and_dedupes():
    text = (DOCKER / "Dockerfile.phase2").read_text()
    phase2 = text[text.index(" AS phase2"):]
    steps = [s for s in run_steps(phase2) if "/opt/venv" in s or "/opt/vllm" in s]
    one = [s for s in steps if "uv venv" in s and "/opt/vllm" in s and "/opt/venv" in s and "dedupe.py" in s]
    assert len(one) == 1, "both venvs and the dedupe must be one RUN: a hardlink across layers is a copy"


def test_dockerfile_compiles_causal_conv1d_in_a_builder_and_fails_loudly():
    text = (DOCKER / "Dockerfile.phase2").read_text()
    assert "MISSING" not in text and "||" not in text  # the silent fallback hid the missing nvcc
    builder = text[text.index(" AS ccd"):text.index(" AS sampler")]
    assert 'cuda-toolkit[nvcc,cccl]==13.0.3' in builder  # torch 2.14.1 pins cuda-toolkit==13.0.3
    assert re.search(r"sha256sum -c", builder) and "patch_arch.py" in builder
    assert "CAUSAL_CONV1D_FORCE_BUILD=TRUE" in builder
    assert "ln -sf libcudart.so.13 $CUDA_HOME/lib/libcudart.so" in builder  # the PyPI runtime ships only .so.13
    phase2 = text[text.index(" AS phase2"):]
    assert "COPY --from=ccd" in phase2 and "causal_conv1d" in phase2
    assert "cuda-toolkit" not in phase2  # the compiler stays in the builder


@pytest.mark.skipif(not shutil.which("bash"), reason="bash")
def test_build_images_copies_the_helpers_into_the_context():
    text = (DOCKER / "build_images.sh").read_text()
    for name in ("dedupe.py", "patch_arch.py"):
        assert name in text


# The box never runs these: Mac dev tools, raw-DICOM fetchers, Phase 0's Qwen2.5-VL helper and its video decoder
JUNK_TRAIN = {"ruff", "pytest", "pytest-cov", "coverage", "hypothesis", "pylidc", "pydicom", "sqlalchemy", "matplotlib",
              "scikit-image", "qwen-vl-utils", "av"}
JUNK_VLLM = {"nvidia-cutlass-dsl-libs-cu12", "pynvvideocodec"}  # CUDA 12 libs; video decoding (we sample images)


def pins(path):
    out = {}
    for line in path.read_text().splitlines():
        line = line.split("#")[0].strip()
        if "==" in line:
            name, ver = line.split("==", 1)
            out[name.strip().lower()] = ver.split(";")[0].strip()
    return out


def test_train_lock_is_pyproject_train_for_linux_pinned_to_the_mac():
    """V0 audit: the lock was the Mac venv's freeze, so the image carried dev and download tools (31 of 109 packages)
    behind a hand-kept exclude list. It is now pyproject's base + [train] for Linux x86-64, constrained to the Mac's
    versions (lock mode of build_images.sh), plus the build backend the box needs for `pip install -e .`."""
    lock = pins(LOCKS / "train.txt")
    assert not JUNK_TRAIN & lock.keys()
    project = tomllib.loads((ROOT / "pyproject.toml").read_text())
    reqs = (project["project"]["dependencies"] + project["project"]["optional-dependencies"]["train"]
            + project["build-system"]["requires"])
    for req in reqs:
        name = re.split(r"[<>=!~\[ ;]", req, maxsplit=1)[0].lower()
        assert name in lock, name
        if "==" in req:
            assert lock[name] == req.split("==")[1].strip(), name
    assert "torchvision" in project["project"]["optional-dependencies"]["train"][-1]  # transformers/TRL vision path
    assert any(k.startswith("nvidia-") for k in lock)  # the Linux CUDA wheels are pinned too, not resolved at build
    assert lock.keys() >= {"pandas"}  # objdel.train.sample


def test_dockerfile_installs_the_locks_and_imports_what_it_keeps():
    """`uv pip install --excludes` drops the vLLM junk even though the lock pins it; the build imports each venv's
    entry points, so a bad lock or exclusion fails on the Mac or in CI, not on a paid box."""
    names = lambda f: {ln.strip() for ln in (DOCKER / f).read_text().splitlines() if ln.strip() and ln[0] != "#"}  # noqa: E731
    assert names("exclude-vllm.txt") == JUNK_VLLM and not (DOCKER / "exclude-venv.txt").exists()
    assert not JUNK_VLLM - pins(LOCKS / "vllm.txt").keys()  # still pinned by vllm, so the exclusion matters
    text = (DOCKER / "Dockerfile.phase2").read_text()
    assert "exclude-venv" not in text and "qwen_vl_utils" not in text
    phase2 = " ".join(run_steps(text[text.index(" AS phase2"):]))
    assert "-r /opt/requirements.txt &&" in phase2.replace("  ", " ")
    assert "-r /opt/requirements-vllm.txt --excludes /opt/exclude-vllm.txt" in phase2
    assert re.search(r"/opt/vllm/bin/python -c \"import vllm", phase2)
    keep = r"/opt/venv/bin/python -c \"import [^\"]*torch[^\"]*torchvision[^\"]*transformers[^\"]*trl"
    assert re.search(keep, phase2)


def test_build_images_builds_from_the_committed_locks():
    """CI has no Mac venv: `--lock` (Mac only) writes locks/, a build only copies them into the context."""
    sh = (DOCKER / "build_images.sh").read_text()
    start = sh.index('if [ "$LOCK" = 1 ]')
    lock_branch = sh[start:sh.index("\nfi\n", start)]
    assert "uv pip freeze" in lock_branch and "pyproject.toml --extra train" in lock_branch
    assert sh.count("uv pip freeze") == 1  # only in lock mode
    assert 'cp "$LOCKS"/' in sh and "exclude-venv" not in sh


def test_image_workflow_is_manual_pinned_and_pushes_with_the_repo_token():
    wf = (ROOT / ".github" / "workflows" / "image.yml").read_text()
    on = wf[wf.index("\non:"):wf.index("\npermissions:")]
    assert "workflow_dispatch" in on and "push" not in on and "pull_request" not in on  # each run is a choice
    assert re.search(r"packages: write", wf) and re.search(r"contents: read", wf)
    uses = re.findall(r"uses: (\S+)", wf)
    assert uses and all(re.search(r"@[0-9a-f]{40}$", u) for u in uses)
    assert "secrets.GITHUB_TOKEN" in wf and "secrets." not in wf.replace("secrets.GITHUB_TOKEN", "")
    assert "build_images.sh --push phase2" in wf

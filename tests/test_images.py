"""Phase 2 image: locks, layers and the CI build. 2026-11-03: one 5.61 GB layer (both venvs, hardlinked) kept
restarting on weak hosts; each venv now installs in four layers."""

import re
import tomllib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DOCKER = ROOT / "experiments" / "phase2" / "docker"
LOCKS = DOCKER / "locks"


def run_steps(text):
    """RUN instructions with their line continuations joined."""
    return [re.sub(r"\\\n\s*", " ", r) for r in re.findall(r"^RUN .*?(?<!\\)$", text, re.M | re.S)]


def test_dockerfile_installs_each_venv_in_small_layers():
    """No RUN builds both venvs (a 5.61 GB layer of both kept restarting on weak hosts); each venv takes four RUNs:
    cuBLAS + cuDNN, the other CUDA wheels, torch, the rest."""
    text = (DOCKER / "Dockerfile.phase2").read_text()
    steps = run_steps(text[text.index(" AS phase2"):])
    assert not [s for s in steps if "/opt/venv/" in s and "/opt/vllm/" in s]
    for venv, lock in (("/opt/venv/", "requirements.txt"), ("/opt/vllm/", "requirements-vllm.txt")):
        mine = [s for s in steps if venv in s]
        assert len(mine) == 4, venv
        assert "nvidia-(cublas|cudnn)" in mine[0] and "(nvidia|cuda)-" in mine[1]
        assert "(torch|torchvision|triton)" in mine[2]
        assert all("--no-deps" in s for s in mine[:3]) and f"-r /opt/{lock}" in mine[3]
    installs = " ".join(run_steps(text))
    assert "causal" not in installs and "flash-linear-attention" not in installs and "dedupe" not in installs


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
    assert 'build_images.sh --push "$WHICH"' in wf and "options: [phase2, vllm_native]" in on
    runs = re.findall(r"run: (.*)", wf)
    assert not any("${{" in r for r in runs)  # inputs reach the shell through env, never pasted into the script

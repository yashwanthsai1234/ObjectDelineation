"""Box inputs ship (Implementation.md front 6.1 step 4): every path a box step reads is in the code sync
(policy.code_files) or in the LIDC bundle (bundle.members). S09: the K3 runbook (D1 step 0; the R-D1
plan it also checked is finished: archived in E6.6, local only).
Reads no data/: the bundle is listed for a temp tree."""
import os
import re
import shutil
import subprocess
from pathlib import Path

import pytest

from objdel.infra import bundle as BU
from objdel.infra import preflight as PF
from objdel.infra.vast import remote as P

ROOT = Path(__file__).resolve().parents[2]
KERNELS = ROOT / "docs/implementation2.0/kernels"
SRC = "/workspace/src/"


def unshipped(text: str, code: list[str], bundle: set[str]) -> list[str]:
    """The paths `text` reads that neither ship nor come from the bundle. `NAME=/workspace/src/...` variables are
    expanded; paths under a build/ folder are made by an earlier step; a relative data/ path never exists on the box
    (the code sync drops data/; the bundle sits at {data})."""
    for name, value in re.findall(r"\b(\w+)=(/workspace/src/[\w./-]+)", text):
        text = text.replace(f"${name}/", f"{value}/").replace(f"${name} ", f"{value} ")
    bad = [f"relative data path: {m}" for m in re.findall(r"(?<![\w{}/.-])(data/[\w./-]+)", text)]
    for kind, rel in re.findall(r"(/workspace/src/|\{data\}/)([\w./-]*\w)", text):
        if "/build/" in f"/{rel}/":
            continue
        have = code if kind == SRC else bundle
        if rel not in have and not any(f.startswith(rel + "/") for f in have):
            bad.append(f"{kind}{rel}")
    return bad


@pytest.fixture(scope="module")
def code() -> list[str]:
    return P.code_files()


@pytest.fixture(scope="module")
def bundle(tmp_path_factory) -> set[str]:
    """The member names of a LIDC bundle built from a minimal temp data tree."""
    d = tmp_path_factory.mktemp("data")
    files = ["prepared/lidc/val.jsonl", "prepared/images/lidc/a.png", "prepared/caches/lidc/1.2/n000/z001/nodule.npz",
             "views/lidc-val-a/rows.jsonl", "checker_eval/manifest.json", "checker_eval/masks.npz",
             *(f"eval_views/{v}/rows.jsonl" for v in BU.LIDC_EVAL)]
    for f in files:
        (d / f).parent.mkdir(parents=True, exist_ok=True)
        (d / f).write_text("x")
    (d / "prepared/manifest.json").write_text('{"files": {}, "sets": {}}')
    return {name for name, _ in BU.members("lidc", d)}


def test_k3_reads_only_shipped_paths(code, bundle):
    assert unshipped((KERNELS / "K3/GPU_RUN.md").read_text(), code, bundle) == []


def test_a_planted_unshipped_read_fails(code, bundle):
    assert unshipped('- [bash, -c, "cat data/k1_bundle/lidc/in.bits"]', code, bundle) == [
        "relative data path: data/k1_bundle/lidc/in.bits"]
    k9 = "/workspace/src/docs/implementation2.0/kernels/K9"
    assert unshipped(f"K9={k9}; cat $K9/x.cpp {{data}}/k1_bundle/in", code, bundle) == [
        k9, f"{k9}/x.cpp", "{data}/k1_bundle/in"]


def test_kernel_drafts_ship_without_build_or_results(code):
    k1 = "docs/implementation2.0/kernels/K1/"
    assert f"{k1}k1poly.cpp" in code and f"{k1}GPU_RUN.md" in code
    assert not [f for f in code if f.startswith((f"{k1}build/", f"{k1}results/"))]
    assert {f.split("/")[2] for f in code if f.startswith("docs/")} <= {"kernels", "frozen_views.json",
                                                                      "val_halves.json"}


def test_runbooks_have_no_scratch_path_hand_copy_or_stale_torch(tmp_path):
    assert PF.runbooks(None) == []
    k = tmp_path / "docs/implementation2.0/kernels/K9"
    k.mkdir(parents=True)
    (k / "GPU_RUN.md").write_text("torch 2.14.1 in /opt/venv\nscp root@h:/x .\nrsync -a a b\n")
    (k / "a.cpp").write_text("// /" + "private/tmp/x\n")
    assert len(PF.runbooks(None, tmp_path)) == 4


# ------------------------------------------------------------------ S14: the old-engine tar ships (D-G7, D8)

BOX = ROOT / "experiments/cuda_checker/run_box.sh"
TAR = ("old_engine.tar", "old_engine.tar.sha256")
TAR_PATHS = ("cuda", "src", "CMakeLists.txt", "pyproject.toml", "Makefile", "README.md")


def _git(cwd: Path, *args: str) -> subprocess.CompletedProcess:
    return subprocess.run(["git", "-c", "user.name=t", "-c", "user.email=t@t", "-c", "commit.gpgsign=false",
                           "-c", "tag.gpgsign=false", *args], cwd=cwd, capture_output=True, text=True)


def _tagged_repo(tmp_path: Path, tag: bool = True) -> Path:
    """A tmp git repo with a copy of run_box.sh, the paths old-engine-tar archives and (tag) a pre-d2 tag."""
    r = tmp_path / "repo"
    (r / "experiments/cuda_checker").mkdir(parents=True)
    shutil.copy(BOX, r / "experiments/cuda_checker/run_box.sh")
    for p in TAR_PATHS:
        f = r / p / "x.txt" if p in ("cuda", "src") else r / p
        f.parent.mkdir(parents=True, exist_ok=True)
        f.write_text(p)
    _git(r, "init", "-q")
    _git(r, "add", ".")
    assert _git(r, "commit", "-q", "-m", "old tree").returncode == 0
    if tag:
        _git(r, "tag", "pre-d2")
    return r


def test_code_files_ship_the_untracked_tar(tmp_path, monkeypatch):
    for f in (*TAR, "tests/verify/test_gpu_scoring.py"):
        (tmp_path / f).parent.mkdir(parents=True, exist_ok=True)
        (tmp_path / f).write_text("x")
    listed = "\0".join([*TAR, "tests/verify/test_gpu_scoring.py", "data/lidc/x.npz"]) + "\0"
    monkeypatch.setattr(P, "repo_root", lambda: tmp_path)
    monkeypatch.setattr(P.subprocess, "run", lambda *a, **k: subprocess.CompletedProcess(a, 0, stdout=listed))
    assert P.code_files() == [*TAR, "tests/verify/test_gpu_scoring.py"]


def test_tar_is_neither_ignored_nor_tracked():
    ignored = subprocess.run(["git", "check-ignore", *TAR], cwd=ROOT, capture_output=True, text=True)  # -q: one path
    assert ignored.returncode == 1 and ignored.stdout == ""
    assert subprocess.run(["git", "ls-files", "--error-unmatch", TAR[0]], cwd=ROOT, capture_output=True).returncode
    target = re.search(r"^old-engine-tar:.*?(?=^\S)", (ROOT / "Makefile").read_text(), re.M | re.S).group(0)
    assert "pre-d2" in target and not re.search(r"\b[0-9a-f]{7,40}\b", target)


def test_box_step_5_builds_the_old_engine_from_the_root_tar():
    text = BOX.read_text()
    assert "sha256sum -c old_engine.tar.sha256" in text and "tar -xf old_engine.tar -C /opt/old_src" in text
    assert "make -C /opt/old_src cuda PY=/opt/old/bin/python" in text
    assert not re.search(r"(runs|build)/old_engine\.tar", text) and "<dir>" not in text and "<old dir>" not in text


def test_make_old_engine_tar_needs_the_tag_and_writes_both_root_files(tmp_path):
    make = ["make", "-s", "-f", str(ROOT / "Makefile"), "old-engine-tar"]
    r = _tagged_repo(tmp_path / "untagged", tag=False)
    p = subprocess.run(make, cwd=r, capture_output=True, text=True)
    assert p.returncode and "owner: git tag pre-d2" in p.stderr and not (r / TAR[0]).exists()
    r = _tagged_repo(tmp_path)
    assert subprocess.run(make, cwd=r, capture_output=True, text=True).returncode == 0
    names = subprocess.run(["tar", "-tf", TAR[0]], cwd=r, capture_output=True, text=True).stdout.split()
    assert sorted(names) == sorted(["cuda/", "cuda/x.txt", "src/", "src/x.txt", *TAR_PATHS[2:]])
    assert subprocess.run(["sha256sum", "-c", TAR[1]], cwd=r, capture_output=True).returncode == 0


def test_gpu_check_without_the_tar_stops_before_the_toolkit(tmp_path):
    r = _tagged_repo(tmp_path)
    env = {**os.environ, "OBJDEL_RUN_DIR": str(tmp_path / "run"), "OBJDEL_DATA_ROOT": str(tmp_path)}
    box = ["bash", "experiments/cuda_checker/run_box.sh", "--gpu-check"]
    for args in ([], ["--dry"]):  # the box mode always needs the tar; the dry run needs it once pre-d2 exists
        p = subprocess.run([*box, *args], cwd=r, env=env, capture_output=True, text=True)
        assert p.returncode == 1 and "old_engine.tar missing" in p.stderr and "uv pip install" not in p.stdout
    assert not (tmp_path / "run").exists()

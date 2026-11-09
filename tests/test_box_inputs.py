"""Box inputs ship (Implementation.md front 6.1 step 4): every path a box step reads is in the code sync
(policy.code_files) or in the LIDC bundle (bundle.members). S09: the R-D1 run config and the K3 runbook (D1 step 0).
Reads no data/: the bundle is listed for a temp tree."""
import re
from pathlib import Path

import pytest

from experiments.phase2 import preflight as PF
from experiments.phase2.vast import bundle as BU
from experiments.phase2.vast import policy as P

ROOT = Path(__file__).resolve().parents[1]
R_D1 = ROOT / "experiments/phase2/configs/r_d1.yaml"
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


def test_r_d1_and_k3_read_only_shipped_paths(code, bundle):
    text = R_D1.read_text()
    assert unshipped(text, code, bundle) == []
    assert unshipped((KERNELS / "K3/GPU_RUN.md").read_text(), code, bundle) == []
    assert "/workspace/src/docs/implementation2.0/kernels/K1" in text and "{data}/checker_eval/masks.npz" in text
    assert "OBJDEL_DATA_ROOT={data}" in next(s for s in text.split("\n  - ") if "bench_checker" in s)


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
    assert PF.runbooks(None, None) == []
    k = tmp_path / "docs/implementation2.0/kernels/K9"
    k.mkdir(parents=True)
    (k / "GPU_RUN.md").write_text("torch 2.14.1 in /opt/venv\nscp root@h:/x .\nrsync -a a b\n")
    (k / "a.cpp").write_text("// /" + "private/tmp/x\n")
    assert len(PF.runbooks(None, None, tmp_path)) == 4

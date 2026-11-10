"""Content tests for the documentation tasks DOC1, DOC2, DOC8 and A9 (plan section DOC.9: I15, I18, X1, U19, I02, S08).

Bad strings are built at run time so that this file never trips the scans it implements.
"""

from __future__ import annotations

import re
import subprocess
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
LOCAL_ONLY = pytest.mark.skipif(
    not (REPO / "docs/plan").exists(), reason="planning records are local-only (git-ignored)"
)

# Frozen paths (DOC.spec "Definitions"); experiments/phase1 is live after DOC4, except its lab folder.
FROZEN_PARTS = ("/lab/", "/impl/")
FROZEN_PREFIXES = ("docs/research/", "experiments/cuda_checker/", "docs/implementation2.0/")
BANNER = "> Frozen record (written "


def _read(rel: str) -> str:
    return (REPO / rel).read_text(encoding="utf-8")


def _git(*args: str) -> str:
    return subprocess.run(["git", "-C", str(REPO), *args], capture_output=True, text=True, check=True).stdout


def _tracked_or_new(*globs: str) -> list[str]:
    out = _git("ls-files", "-z", "--cached", "--others", "--exclude-standard", "--", *globs)
    return [p for p in out.split("\0") if p and (REPO / p).is_file()]


def _is_frozen(rel: str) -> bool:
    return rel.startswith(FROZEN_PREFIXES) or any(part in "/" + rel for part in FROZEN_PARTS)


# ---------------------------------------------------------------- DOC1


@LOCAL_ONLY
def test_plan_moved_out_of_frozen_folder():
    names = sorted(p.name for p in (REPO / "docs/plan").iterdir())
    assert names == ["Architecture.md", "Direction.md", "Research2.0.md", "previous-vs-new.md"]
    assert not (REPO / "docs/research" / ("research" + "2.0")).exists()


def test_no_reference_to_the_old_plan_path():
    old = ("research" + "2.0", "previous " + "vs new", "previous%20" + "vs%20new")
    targets = ["README.md", "ROADMAP.md", "Makefile", "pyproject.toml", "docs", "src", "tests", "experiments"]
    targets.append(".github")
    hits = []
    for rel in _tracked_or_new(*targets):
        if rel.startswith(("docs/research/", "docs/implementation2.0/")):
            continue  # frozen, and the implementation plan names the old path as history
        try:
            text = _read(rel)
        except UnicodeDecodeError:
            continue
        hits += [f"{rel}: {o}" for o in old if o in text]
    assert hits == []


@LOCAL_ONLY
def test_research_page_corrected_ceiling():
    row = next(ln for ln in _read("docs/plan/Research2.0.md").splitlines() if ln.startswith("| All images"))
    assert "0.723" in row and "0.734" not in row


def test_gitignore_has_docs_build_rules():
    rules = set(_read(".gitignore").splitlines())
    wanted = ("/site/", "/.cache/", ".lycheecache", "docs/research/*.pdf", "docs/implementation2.0/kernels/*/build/")
    for rule in wanted:
        assert rule in rules, rule
    ignored = subprocess.run(
        ["git", "-C", str(REPO), "check-ignore", "-q", "docs/research/2212.07020v1.pdf"], capture_output=True
    )
    assert ignored.returncode == 0


def test_no_scratch_or_home_paths_in_implementation_docs():
    patterns = [
        "/" + "private/tmp",
        "/" + "Users/",
        "S" + "P/impl",
        "S" + "P/[A-Z]",
    ]
    rx = re.compile("|".join(patterns))
    hits = []
    for rel in _tracked_or_new("docs/implementation2.0"):
        if rel.startswith("docs/implementation2.0/notes/research/"):
            continue  # written by the research wave, scanned by its own owner
        try:
            lines = _read(rel).splitlines()
        except UnicodeDecodeError:
            continue
        hits += [f"{rel}:{i}" for i, ln in enumerate(lines, 1) if rx.search(ln)]
    assert hits == []


def test_kernel_build_output_is_not_tracked_or_cited():
    assert [p for p in _git("ls-files", "docs/implementation2.0/kernels").splitlines() if "/build/" in p] == []
    for k in ("K2", "K3"):
        assert not re.search(r"\bbuild/", _read(f"docs/implementation2.0/kernels/{k}/DESIGN.md"))
    results = REPO / "docs/implementation2.0/kernels"
    assert not list(results.rglob("*.log"))  # *.log is git-ignored: evidence is kept as .log.txt


# ---------------------------------------------------------------- DOC2


def test_policy_docstring_has_no_dollar_value():
    import importlib

    policy = importlib.import_module("experiments.phase2.vast.policy")
    assert re.search(r"\$\d", policy.__doc__) is None
    assert 0 < policy.ALARM < policy.HARD_STOP


def test_gpu_run_runbook_has_no_stale_policy_sentences():
    text = _read("docs/runbooks/gpu-run.md")
    assert "docstring still says" not in text
    sweeper = next(ln for ln in text.splitlines() if ln.startswith("| Sweeper | Destroys"))
    assert "`HARD_STOP`" in sweeper and "`ALARM`" in sweeper and "$" not in sweeper.split("|")[2]


def test_adr_0009_action_item_3_ticked_and_history_kept():
    text = _read("docs/architecture/decisions/0009-vast-ai-north-america-and-preflight.md")
    assert "3. [x] Update the docstring of `policy.py`" in text
    assert "1. [x]" in text and "2. [x]" in text
    assert "The module docstring still says $15." in text  # history line is left as it was


def test_no_stale_trl_version_text_outside_frozen_paths():
    bad_name = "TRL " + "1.16"
    bad_url = "blob/v" + "1.16.0/"
    hits = []
    for rel in _tracked_or_new("src", "tests", "experiments"):
        name = Path(rel).name
        live_type = (
            rel.endswith(".py") or name.startswith("Dockerfile") or (rel.endswith(".yaml") and "/configs/" in rel)
        )
        if not live_type or _is_frozen(rel):
            continue
        text = _read(rel)
        hits += [rel for s in (bad_name, bad_url, "trl==" + "1.16") if s in text]
    assert hits == []


def test_r2_config_comment_is_not_history_wording():
    text = _read("experiments/phase2/configs/r2_qwen3vl2b.yaml")
    assert "first run on a GPU here" not in text
    assert "first GPU run: r2 (see docs/results/runs.md)" in text


FROZEN_RECORDS = (
    ["experiments/phase2/PLAN.md", "experiments/phase2/V0_V1.md", "experiments/phase1/PLAN.md"]
    + [f"experiments/phase2/JUDGE_R{i}.md" for i in range(1, 5)]
    + ["docs/research/records/phase0/PLAN.md", "docs/research/records/phase0/PLAN_0b.md"]
)


@LOCAL_ONLY
def test_frozen_records_carry_the_banner_on_line_one_only():
    for rel in FROZEN_RECORDS:
        lines = _read(rel).splitlines()
        assert re.fullmatch(
            r"> Frozen record \(written \d{4}-\d{2}-\d{2}\)\. Versions, budgets and paths were true then\. "
            r"Current facts: docs/README\.md\.",
            lines[0],
        ), rel
        assert sum(ln.startswith(BANNER) for ln in lines) == 1, rel


def test_no_banner_in_python_or_dockerfiles():
    for rel in _tracked_or_new("src", "tests", "experiments"):
        if rel != "tests/test_docs_content.py" and (rel.endswith(".py") or Path(rel).name.startswith("Dockerfile")):
            assert BANNER not in _read(rel), rel


# ---------------------------------------------------------------- A9 pages (S08) and DOC8


def _ticks(rel: str) -> set[str]:
    """The code spans of a page."""
    return set(re.findall(r"`([^`\n]+)`", _read(rel)))


def _a9_items():
    """Key sets of the schema-2 files, taken from the code: final_eval.summarize on fake items, one metrics() call."""
    import numpy as np

    from experiments.phase2 import final_eval as FE
    from objdel.eval import metrics

    ids = [f"lidc/p{p}/n001/z{s}" for p in range(2) for s in range(2)]
    items = [{"item_id": i, **{m: 0.5 + k / 10 for m in FE.METRICS}, **{m: 0.5 for m in FE.SELECTION}}
             for k, i in enumerate(ids)]
    items[0]["iou_tc"] = 0.05  # one miss, so that iou_found has a subset
    summary = FE.summarize(items, {i: i.split("/")[1] for i in ids}, draws=20)
    yy, xx = np.ogrid[:50, :50]
    raters = {"nodule": np.stack([(yy - 25) ** 2 + (xx - 25) ** 2 <= r * r for r in (12, 14)])}
    return summary, metrics(['{"nodule":[]}'] * 2, raters, 50, 50)


def test_metrics_page_names_every_metric_and_the_ci_rules():
    page = _read("docs/reference/metrics.md")
    for word in ("fair ceiling", "`ci95_nodule`", "`ci95_patient`", "`ci95_item`", "`gate_ci`", "D²_GED", "HM-IoU8",
                 "soft Dice", "NSD", "boundary IoU", "miss rate", "`iou_found`", "`sel_medoid`", "`sel_vote`",
                 "`noise_floor.item_sd`", "`noise_floor.mean_sd`", "`cluster_t`", "CR1", "no BCa",
                 "CIs exclude training-seed variance; a gain is a single-run result"):
        assert word in page, word
    for path in re.findall(r"`((?:src|experiments)/[\w/]+\.py)`", page):
        assert (REPO / path).is_file(), path
    summary, _ = _a9_items()
    assert {"ged", "hm_iou", "soft_dice", "nsd", "boundary_iou", "miss_rate"} <= set(summary) and "`ged`" in page


def test_output_format_page_lists_every_schema_2_key():
    summary, m = _a9_items()
    page, ticks = _read("docs/reference/output-format.md"), _ticks("docs/reference/output-format.md")
    keys = {*summary, *summary["iou_tc"], *summary["miss_rate"], *summary["iou_found"], *summary["paired"],
            *summary["paired"]["medoid_minus_mean"], *m, *m["samples"], *m["objects"]["nodule"], *m["groups"],
            "item_id", "dataset", "group_id", "W", "H", "samples", "truncated", "tokens", "partial", "seconds",
            "noise_floor", "item_sd", "mean_sd", "limit", "rescored_from", "view_manifest_sha256",
            "view_content_sha256", "code_sha", "checker_build", "ckpt", "view", "split", "K", "seed", "metrics_absent",
            "resummarized_from", "ci95_item", "ci95_nodule", "ci95_patient", "gate_ci"}
    assert sorted(k for k in keys if k not in ticks) == []
    assert "no bare `ci95`" in page and "CODE_SHA" in page


def test_profiling_page_lists_the_range_names_and_the_gate_numbers():
    from objdel.train import nvtx

    page = _read("docs/guides/profiling.md")
    names = re.search(r"Range names are fixed .*?: (.*?)\.\n", nvtx.__doc__, re.S)[1].replace("\n", " ").split(", ")
    from test_nvtx import NAMES

    assert set(names) == NAMES and len(names) == len(NAMES)
    assert sorted(n for n in names if f"`{n}`" not in page) == []
    for word in ("OBJDEL_NVTX", "OBJDEL_NVTX_SYNC", "nsys_capture.sh", "nsys_baseline", "GRAPH_TABLES", "graph_rows",
                 "upper bound", "--profile-window", "--trace", "--count-syncs", "`coverage`", "0.90", "fallback"):
        assert word in page, word


def test_glossary_has_the_a9_terms():
    ticks = {t.strip("*") for t in re.findall(r"^\| \*\*(.+?)\*\* \|", _read("docs/glossary.md"), re.M)}
    wanted = {"fair ceiling", "medoid", "miss", "missed answer", "found image", "gate CI", "cluster", "doctor",
              "sample"}
    assert wanted - ticks == set()


def test_the_item_ci_is_labelled_where_the_headline_is_quoted():
    text = _read("docs/guides/evaluation.md")
    assert "0.645 [0.633, 0.656]" in text and "item CI" in text
    assert "[0.625, 0.664]" in text and "[0.628, 0.662]" in text


def test_contributing_states_the_archive_rules():
    text = _read("docs/contributing.md")
    for part in ("What it was", "What replaced it", "Retiring commit", "How to restore"):  # test_docs.ARCHIVE_PARTS
        assert part in text, part
    assert "not built, not imported, not tested, not linted" in text
    assert re.search(r'extend-exclude = \[.*"archive"', _read("pyproject.toml"))

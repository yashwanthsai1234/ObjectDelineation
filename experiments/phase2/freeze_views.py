"""Freeze eval views into docs/reference/frozen_views.json (B1.6, DG15): one entry per view, keyed by its name.

    .venv/bin/python -I -m experiments.phase2.freeze_views VIEW_DIR [VIEW_DIR ...] [--out FILE] [--known FILE]
    .venv/bin/python -I -m experiments.phase2.freeze_views --check [--out FILE]

An entry holds the content hash (objdel.data.quality.view_content_sha256: the records' item_id, view.offset,
view.side and image_sha256, so another Python or library build that renders the same windows and pixels gives the
same hash; the manifest SHA does not), the target-set hash, the build arguments (seed, kinds, store manifest SHA-256),
the fallback and padded counts, the commit the build ran on, each render file with its last commit (the content
hash holds image_sha256 and the target-set hash, which come from objdel.train.dataset._prepare), the R5 code commit
(runs/r5/CODE_SHA once A3's sync writes it, else null) and the view's known label issues (B1.7 report). Render files
with uncommitted edits stop the freeze. --check lists the render files whose last commit moved since each entry was
frozen (a report: the measured check is the content hash). final_eval refuses a frozen view that differs.
"""
import argparse
import json
import sys
from pathlib import Path

from objdel.data import quality as Q
from objdel.io import code_sha, git, git_sha

ROOT = Path(__file__).resolve().parents[2]
RENDER_FILES = ("src/objdel/targets/augment.py", "src/objdel/targets/frame.py", "src/objdel/targets/build.py",
                "src/objdel/targets/simplify.py", "src/objdel/train/dataset.py", "experiments/phase2/build_eval_views.py")
R5_IDENTITY = {"lidc-test-shifted": (
    "equal to the R5 box view by code identity, not measured: targets/frame.py unchanged since 734b1ec; augment.py and "
    "build_eval_views.py changed in B1.6 without moving a window (their Mac rebuild of the R5 val shifted, rescaled and "
    "rotated views gives the stored content and target-set hashes); other render files as listed")}


def _git(*args, root=ROOT) -> str:
    return git(*args, root=root)  # the box has no .git: a missing one raises


def render_files(root=ROOT) -> dict[str, str]:
    return {f: _git("log", "-1", "--format=%H", "--", f, root=root) for f in RENDER_FILES}


def dirty(root=ROOT) -> list[str]:
    """Render files with uncommitted changes."""
    return _git("diff", "--name-only", "HEAD", "--", *RENDER_FILES, root=root).splitlines()


def drift(entry: dict, root=ROOT) -> list[str]:
    """Render files whose last commit differs from the entry's."""
    now = render_files(root)
    return [f for f, c in entry["render_files"].items() if now.get(f) != c]


def entry(view_dir, known: list[str]) -> dict:
    d = Path(view_dir)
    m = json.loads((d / "manifest.json").read_text())
    (s,) = m["sets"].values()
    ids = {json.loads(x)["item_id"] for x in next((d / "lidc").glob("*.jsonl")).read_text().splitlines()}
    return {"view": m["view"], "items": s["items"], "content_sha256": Q.view_content_sha256(d),
            "target_set_sha256": s["target_set_sha256"], "seed": m["seed"], "kinds": m["kinds"],
            "store_manifest_sha256": m["store_manifest_sha256"], "fallback": s["fallback"], "padded": s["padded"],
            "generator_commit": git_sha(ROOT), "render_files": render_files(),
            "r5_code_commit": code_sha(ROOT / "runs/r5", required=False),
            "known_label_issues": [i for i in known if i in ids]} | (
        {"r5_identity": R5_IDENTITY[m["view"]]} if m["view"] in R5_IDENTITY else {})


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("views", nargs="*", type=Path)
    ap.add_argument("--out", type=Path, default=ROOT / "docs/reference/frozen_views.json")
    ap.add_argument("--known", type=Path, default=ROOT / "runs/b1/series_z.json",
                    help="check_series_z report; its known_label_issues (val and test ids)")
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args(argv)
    old = json.loads(a.out.read_text()) if a.out.is_file() else {}
    if a.check:
        for v, e in old.items():
            print(v, "render files moved:", drift(e) or "none")
        return
    if bad := dirty():
        sys.exit(f"freeze_views: commit the render files first: {bad}")
    report = json.loads(a.known.read_text()) if a.known.is_file() else {}
    known = [i for ids in report.get("known_label_issues", {}).values() for i in ids]
    new = {(e := entry(d, known))["view"]: e for d in a.views}
    a.out.write_text(json.dumps(old | new, indent=1, sort_keys=True) + "\n")
    print(json.dumps({v: e["content_sha256"] for v, e in new.items()}, indent=1))


if __name__ == "__main__":
    main()

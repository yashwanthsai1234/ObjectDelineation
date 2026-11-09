"""Pure helpers of the docs tests. DOC6 uses headline_strings; the DOC4 drift checks (tests/test_docs.py) bring
their own helpers back here when DOC4 is merged into this branch (unused helpers fail the dead-code gate)."""

from __future__ import annotations


def headline_strings(h: dict) -> list[str]:
    """The test numbers as README.md writes them (rounded only here), plus '<unit>-level' when not item-level."""
    m = {k: v["mean"] for k, v in h["test"]["metrics"].items()}
    out = [f"IoU {m['iou_tc']:.3f}", f"boundary F-score {m['f5']:.3f}", f"{m['legal'] * 100:.1f} %"]
    return out + ([f"{h['ci_unit']}-level"] if h["ci_unit"] != "item" else [])

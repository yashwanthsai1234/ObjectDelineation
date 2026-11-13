"""Outline text of the raters' masks, built the way the library builds an SFT target."""
from objdel.data.targets.build import build_target, min_areas


def rater_texts(stacks, names=None) -> list[str]:
    """One outline text per rater: build_target of rater r's masks of `stacks` {name: (k, H, W)}, at the size rule
    min_areas gives for the item."""
    names = tuple(stacks if names is None else names)
    ma = min_areas({n: stacks[n] for n in names})
    k = len(stacks[names[0]])
    return [build_target({n: stacks[n][r] for n in names}, names, min_area=ma).text for r in range(k)]


def mask_text(mask, name="nodule") -> str:
    """The outline text of one mask, sized by its own rule (exact for a one-rater nodule)."""
    return build_target({name: mask}, (name,), None, min_areas({name: mask[None]})).text

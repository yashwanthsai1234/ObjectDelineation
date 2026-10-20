"""The single geometry spec for model text output, shared by SFT targets, the reward and evaluation.

text:  GRID, SCHEMAS, caps, the compact grammar, serialize / parse (strict, never raises), pattern (regex).
shape: canonicalize (clockwise outers, counter-clockwise holes, top-left starts, sorted parts) and to_mask
       (even-odd mask + legality through objdel.checker.polycore only).
"""
from .shape import Geometry, canonicalize, to_mask
from .text import GRID, MAX_CHARS, MAX_RINGS, MAX_VERTICES, SCHEMAS, Objects, Parsed, Part, parse, pattern, serialize

__all__ = [
    "GRID", "MAX_CHARS", "MAX_RINGS", "MAX_VERTICES", "SCHEMAS", "Geometry", "Objects", "Parsed", "Part",
    "canonicalize", "parse", "pattern", "serialize", "to_mask",
]

"""Model text <-> polygon objects: the compact grammar, its serializer and its strict parser.

Grammar (no whitespace inside; surrounding JSON whitespace is stripped; names in schema order, all present once; a
dataset's schema is its structure names, objdel.core.datasets):

    doc   := '{' '"' name '":' parts (',' '"' name '":' parts)* '}'
    parts := '[' ']' | '[' part (',' part)* ']'
    part  := '{"outer":' ring ',"holes":' rings '}'          "holes" is always present
    rings := '[' ']' | '[' ring (',' ring)* ']'
    ring  := '[' pt ',' pt ',' pt (',' pt)* ']'               >= 3 points
    pt    := '[' int ',' int ']'                              int := 0 | [1-9][0-9]{0,2} | 1000

Every text in the grammar is also compact JSON (json.dumps with separators=(",", ":")). Caps: <= MAX_VERTICES
points per ring and <= MAX_RINGS rings per object. A completion cut at the token limit is a format failure even
if its prefix parses.

Surrounding whitespace is exactly JSON's insignificant whitespace, space, tab, LF and CR (RFC 8259 §2); str.strip()
with no argument would also strip Unicode spaces and the separators U+001C-U+001F.
Texts longer than MAX_CHARS = 16,384 characters are refused first ("too_long"), so parse time is bounded by the
cap and not by the input. The densest grammar text measured costs 2.14 characters per token (a retired 450M
tokenizer; Qwen3-VL 1.7; real targets and adversarial texts), so a 1,536-token completion is at most ~3,300
characters and the longest of the 30,843 targets in data/targets (summary.json) has 2,864: the cap leaves 5x
headroom over the largest completion budget and still rejects a 10 MB input in microseconds (implementation note
fixes/F1.md).

Parsing is a regular-language check followed by json.loads: `pattern(names)` is a regex for the grammar,
compiled once per schema; `re.fullmatch` decides acceptance and json.loads (C) builds the lists. The grammar is
deterministic (every choice is decided by the next one or two characters), so a failed match backtracks O(1)
per quantifier step and the check is linear; tests feed 200 k-character near-misses. Measured on 50-1,000-vertex
completions it is 5-6x faster than a hand-written recursive-descent parser.
Only rejected texts are walked again, by `_diagnose`, to give a short stable error code for typed logging.
The same regex string has no Python-only syntax (no lookarounds, possessive or atomic groups), so it can be passed
to vLLM / xgrammar constrained decoding (TRL `GRPOConfig.vllm_structured_outputs_regex`).

References:
    https://docs.python.org/3/library/re.html (fullmatch; [0-9] rather than \\d, which matches Unicode digits)
    https://docs.python.org/3/library/json.html (object_pairs_hook keeps duplicate keys and their order;
        parse_int receives the digit string, so huge integers never reach int() and its 4,300-digit limit)
    https://docs.python.org/3/library/stdtypes.html#int-max-str-digits (int-string conversion length limit)
    https://datatracker.ietf.org/doc/html/rfc8259#section-2 (insignificant whitespace: %x20 / %x09 / %x0A / %x0D)
    https://docs.python.org/3/library/stdtypes.html#str.strip (default strips every str.isspace() character)
    https://www.regular-expressions.info/catastrophic.html (backtracking blow-up needs ambiguous nested repeats)
    https://docs.python.org/3/whatsnew/3.11.html#re (possessive quantifiers: faster, but not portable, so unused)
    https://docs.vllm.ai/en/latest/features/structured_outputs/ (regex structured outputs; Rust-style regex)
    https://xgrammar.mlc.ai/docs/latest/api/python/grammar.html (Grammar.from_regex; compiles this pattern)
    https://github.com/huggingface/trl/blob/main/trl/trainer/grpo_config.py (vllm_structured_outputs_regex)
"""
from __future__ import annotations

import json
import re
from functools import cache
from itertools import chain
from typing import NamedTuple

import numpy as np

GRID = 1000
MAX_VERTICES, MAX_RINGS = 1000, 10
MAX_CHARS = 16_384  # >= 5x a 1,536-token completion at 2.14 characters per token (module docstring)
_JSON_WS = " \t\n\r"  # RFC 8259 §2

Part = tuple[np.ndarray, list[np.ndarray]]  # (outer (k, 2) int64, holes)
Objects = dict[str, list[Part]]


class Parsed(NamedTuple):
    """Result of parsing completion text: success flag, the objects (None on failure) and an error code (empty on
    success).
    """
    ok: bool
    objs: Objects | None
    error: str  # "" when ok, else a short stable code (see parse)


_INT = "(?:1000|[1-9][0-9]{0,2}|0)"
_PT = rf"\[{_INT},{_INT}\]"
_RING = rf"\[{_PT}(?:,{_PT}){{2,}}\]"
_RINGS = rf"\[(?:{_RING}(?:,{_RING})*)?\]"
_PART = rf'\{{"outer":{_RING},"holes":{_RINGS}\}}'
_PARTS = rf"\[(?:{_PART}(?:,{_PART})*)?\]"


def pattern(names) -> str:
    """Regex of the grammar for one schema (caps excluded); portable to Rust-style regex engines."""
    return r"\{" + ",".join(f'"{re.escape(n)}":{_PARTS}' for n in names) + r"\}"


@cache
def _regex(names: tuple[str, ...]) -> re.Pattern:
    return re.compile(pattern(names))


def serialize(objs: Objects, names) -> str:
    """Objects -> compact text, names in the given order (canonicalize first for SFT targets)."""
    doc = {n: [{"outer": np.asarray(o).tolist(), "holes": [np.asarray(h).tolist() for h in hs]} for o, hs in objs[n]]
           for n in names}
    return json.dumps(doc, separators=(",", ":"))


def _ring(points: list) -> np.ndarray:
    return np.fromiter(chain.from_iterable(points), np.int64, 2 * len(points)).reshape(-1, 2)


def parse(text, names, truncated: bool = False) -> Parsed:
    """Strict parser; never raises. Error codes: truncated, not_text, too_long, empty, not_ascii, whitespace,
    leading_zero, incomplete, trailing_text, syntax, keys, part_keys, schema, not_int, range, few_points, cap_rings,
    cap_vertices."""
    if truncated:
        return Parsed(False, None, "truncated")
    if not isinstance(text, str):
        return Parsed(False, None, "not_text")
    if len(text) > MAX_CHARS:
        return Parsed(False, None, "too_long")
    names = tuple(names)
    s = text.strip(_JSON_WS)
    if not _regex(names).fullmatch(s):
        return Parsed(False, None, _diagnose(s, names))
    doc = json.loads(s)
    for parts in doc.values():
        if sum(1 + len(p["holes"]) for p in parts) > MAX_RINGS:
            return Parsed(False, None, "cap_rings")
        if any(len(r) > MAX_VERTICES for p in parts for r in (p["outer"], *p["holes"])):
            return Parsed(False, None, "cap_vertices")
    objs = {n: [(_ring(p["outer"]), [_ring(h) for h in p["holes"]]) for p in doc[n]] for n in names}
    return Parsed(True, objs, "")


class _Obj(list):
    """A JSON object as its (key, value) pairs, so duplicates and key order stay visible."""


def _int(digits: str) -> int:
    """JSON integer literal -> int; 6+ characters (sign included) is out of [0, GRID] anyway, so it becomes GRID + 1
    ("range") without int()'s 4,300-digit ValueError."""
    return int(digits) if len(digits) < 6 else GRID + 1


def _diagnose(s: str, names: tuple[str, ...]) -> str:
    """Error code for a text the grammar rejected (the first problem found, coarse to fine)."""
    if not s:
        return "empty"
    if not s.isascii():
        return "not_ascii"
    if re.search(r"\s", s):
        return "whitespace"
    try:
        doc = json.loads(s, object_pairs_hook=_Obj, parse_int=_int)
    except json.JSONDecodeError as e:
        if re.search(r"(?<![0-9.eE])0[0-9]", s):  # JSON forbids leading zeros too
            return "leading_zero"
        return "incomplete" if e.pos >= len(s) else "trailing_text" if e.msg == "Extra data" else "syntax"
    except RecursionError:  # nesting deeper than the decoder's recursion limit
        return "syntax"
    if not isinstance(doc, _Obj) or [k for k, _ in doc] != list(names):
        return "keys"
    for _, parts in doc:
        if type(parts) is not list:
            return "schema"
        for part in parts:
            if not isinstance(part, _Obj) or [k for k, _ in part] != ["outer", "holes"]:
                return "part_keys"
            outer, holes = part[0][1], part[1][1]
            if type(holes) is not list:
                return "schema"
            for ring in (outer, *holes):
                if type(ring) is not list or any(type(p) is not list or len(p) != 2 for p in ring):
                    return "schema"
                values = [v for p in ring for v in p]
                if any(type(v) is not int for v in values):
                    return "not_int"
                if any(not 0 <= v <= GRID for v in values):
                    return "range"
                if len(ring) < 3:
                    return "few_points"
    return "syntax"  # valid JSON of the right shape but not in the grammar's spelling (e.g. escapes, -0)

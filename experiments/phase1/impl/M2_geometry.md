# M2: geometry spec module (`objdel.geometry`)

Gate: R1 parser part (plan §7). Promotes `lab/t_stack/prod.py` (serializer, parser, geometry) and the canonical
form of `lab/t_targets/targets/builder.py` (`canon_ring`, `canonical_parts`).

## What was built

`src/objdel/geometry/` (about 240 lines, half of them docstrings):

| File | Contents |
|---|---|
| `text.py` | `GRID`, `SCHEMAS`, `MAX_VERTICES = 1000`, `MAX_RINGS = 10`, `Part`, `Objects`, `Parsed`; `pattern(names)` (regex of the grammar); `serialize`; `parse` (strict, never raises, error codes) |
| `shape.py` | `Geometry`; `canonicalize`; `to_mask` (polycore `pack` / `validate(roles)` / `rasterize` only) |
| `__init__.py` | re-exports |

Deleted: `src/objdel/checker/polygon_io.py` and `tests/test_polygon_io.py` (Phase 0 PIL path: no pixel-centre
sampling, inflates small nodules 1.56×). `src/objdel/checker/__init__.py` never imported it.
Still importing it: `experiments/phase0/baseline_vlm.py` (frozen record, left to the coordinator) and lab files
under `experiments/phase1/lab/` (frozen evidence). *(Decided: Phase 0 scripts run against commit `faaacc4`,
`experiments/phase0/README.md`.)*

## Grammar (unchanged from the lab)

```
doc   := '{' '"' name '":' parts (',' '"' name '":' parts)* '}'     names in schema order, each once
parts := '[' ']' | '[' part (',' part)* ']'
part  := '{"outer":' ring ',"holes":' rings '}'                      "holes" always present
rings := '[' ']' | '[' ring (',' ring)* ']'
ring  := '[' pt ',' pt ',' pt (',' pt)* ']'                          >= 3 points
pt    := '[' int ',' int ']'          int := 0 | [1-9][0-9]{0,2} | 1000
```

No whitespace inside; surrounding whitespace is stripped (as in the lab; since F1 / B1-1 only JSON whitespace:
space, tab, LF, CR). Caps: ≤ 1,000 points per ring,
≤ 10 rings per object (outer + holes, over all parts). `truncated=True` is always a format failure.
`serialize` is `json.dumps(..., separators=(",", ":"))`: every grammar text is compact JSON and the output equals
the lab's f-string serializer on all 6,338 fixture targets.

Error codes (stable, short; the first problem found, coarse to fine): `truncated`, `not_text`, `empty`,
`not_ascii`, `whitespace`, `leading_zero`, `incomplete` (JSON ends early), `trailing_text`, `syntax`, `keys`
(object names: order, missing, extra, duplicate), `part_keys` (`outer`/`holes` missing, extra, reordered,
duplicated), `schema` (wrong nesting/arity), `not_int` (floats, booleans, strings), `range` (< 0 or > 1000),
`few_points`, `cap_rings`, `cap_vertices`. *Since F1 (I4): `too_long` first, for texts over `MAX_CHARS` = 16,384.*

## Parser: research, alternatives, numbers

Sources:
- Python `re` (`fullmatch`; `\d` matches Unicode digits, so the grammar uses `[0-9]`): https://docs.python.org/3/library/re.html
- Possessive quantifiers / atomic groups exist in `re` since 3.11: https://docs.python.org/3/whatsnew/3.11.html#re ,
  https://github.com/python/cpython/pull/31982 , https://www.regular-expressions.info/possessive.html
- Catastrophic backtracking needs ambiguous nested repetition: https://www.regular-expressions.info/catastrophic.html
- `json` `object_pairs_hook` (keeps duplicate keys and order, used only to name errors): https://docs.python.org/3/library/json.html
- vLLM structured outputs (regex / JSON schema / grammar; xgrammar, guidance and outlines use Rust-style regex;
  `guided_regex` became `structured_outputs={"regex": ...}` in v0.12): https://docs.vllm.ai/en/latest/features/structured_outputs/
- TRL 1.16 `GRPOConfig.vllm_structured_outputs_regex` (checked in `.venv`): https://github.com/huggingface/trl/blob/main/trl/trainer/grpo_config.py
- xgrammar `Grammar.from_regex` (ECMA-262 regex syntax, converted to EBNF): https://xgrammar.mlc.ai/docs/api/python/grammar.html

Sandbox scripts: `scratchpad/impl/M2/bench_parse.py`, `bench_conv.py`, `bench_final.py` (Apple M3, other agents
running; best of 5).

Building blocks, 2-object RIGA-like completion (disc n + cup n + one hole n/4 vertices), µs:

| n | chars | lab parser | regex gate (plain) | regex gate (possessive) | gate + json.loads + np.array | json walk (no regex) | json.loads alone |
|---|---|---|---|---|---|---|---|
| 50 | 1,185 | 120 | 8.4 | 5.4 | 29 | 34 | 10 |
| 100 | 2,315 | 238 | 16 | 10.5 | 57 | 66 | 21 |
| 200 | 4,565 | 471 | 29 | 21 | 113 | 133 | 41 |
| 400 | 9,065 | 931 | 48 | 40 | 218 | 251 | 79 |
| 1,000 | 22,565 | 2,317 | 108 | 103 | 570 | 639 | 200 |

List → array conversion (400-point ring): `np.array(list, int64)` 42 µs, `np.fromiter(chain(...), int64)` 15 µs
(chosen); regex `findall` of all integers 180 µs (rejected). orjson is not installed and json.loads is not the
bottleneck after the regex gate, so no new dependency.

Production `parse` vs the lab parser, µs per completion:

| vertices | one ring (LIDC-like): lab / production | disc + cup (RIGA-like): lab / production |
|---|---|---|
| 50 | 77 / 16 | 165 / 33 |
| 100 | 161 / 32 | 312 / 66 |
| 200 | 295 / 63 | 696 / 113 |
| 400 | 678 / 116 | 1,255 / 220 |
| 1,000 | 1,555 / 287 | 3,487 / 537 |

On the 6,338 fixture targets (lab builder, both modes, every rater): median 30 µs, p99 214 µs per parse.

Adversarial rejections, µs (lab / production, code):

| input | chars | lab | production | code |
|---|---|---|---|---|
| 1,000-vertex doc + `x` | 10,043 | 1,617 | 420 | `trailing_text` |
| 1,000-vertex doc cut 3 chars early | 10,039 | 1,772 | 462 | `incomplete` |
| 100 k points, never closed | 600,019 | 1,125 | 40,086 | `incomplete` |
| `[` × 100,000 | 100,000 | 1 | 659 | `syntax` |
| a 100,000-digit integer | 100,020 | 4,944 | 543 | `syntax` |
| 5,000 parts, junk at the end | 205,020 | 31,318 | 12,468 | `trailing_text` |

Scaling of the plain regex on failing input (1k/2k/4k/8k parts with holes): 0.53 / 1.05 / 2.5 / 5.1 ms, linear.
The grammar is deterministic (each choice is decided by the next one or two characters), so giving up a
repetition fails in O(1) and the backtracking engine stays linear; possessive quantifiers were 1.0–1.5× faster
but are not valid in Rust-style regex, so one plain pattern serves both the parser and constrained decoding.
The 100 k-point case is the only one where the lab is faster (it stops at the vertex cap while reading); 40 ms
on 600 kB is linear and far beyond any `max_new_tokens`.

**Decision:** `re.fullmatch(pattern(names))` decides acceptance, then `json.loads` and `np.fromiter` build the
arrays and the caps are checked. Rejected texts are walked again (`_diagnose`) only to name the error. Exactness
rests on one definition (the regex); the diagnosis can only change the code, and falls back to `syntax`.

**Constrained decoding (optional, done because it was free):** `pattern(names)` has no lookarounds, possessive
or atomic groups and no `\d`. xgrammar 0.2.8 (`Grammar.from_regex`, scratch venv) compiles it in 2 ms and agrees
with the parser on 9 accept/reject cases (`scratchpad/impl/M2/xg_check.py`). It can be passed as
`GRPOConfig(vllm_structured_outputs_regex=pattern(names))`; caps are not in the regex. Whether to constrain
decoding at all is an M7 / Phase 2 decision (it removes format errors from RL, which changes the reward mix).

**Lab bug found:** the lab parser tested digits with `str.isdigit` and converted with `int()`, so it accepted
Unicode digits (`[1١,1]` parsed as 11). The production parser rejects them (`not_ascii`); the oracle copy in
`tests/test_geometry_oracle.py` uses ASCII digits.

## Canonical form: research and decision

Sources:
- PolyFormer (CVPR 2023) §3: vertices clockwise, start "closest to the top left corner of the image", polygons of
  one object "ordered based on the distance between their starting points and the image origin";
  canonical order +12.4 mIoU over random: https://arxiv.org/abs/2302.07387 (text read at https://ar5iv.labs.arxiv.org/html/2302.07387)
- VectorLLM: start at the top-left corner point, clockwise; random start helps pretraining, a fixed start rule
  "brings slight performance improvements and stability" in SFT: https://arxiv.org/abs/2507.04664
- GeoJSON RFC 7946 §3.1.6, right-hand rule (exterior counter-clockwise, holes clockwise, y up):
  https://datatracker.ietf.org/doc/html/rfc7946#section-3.1.6
- SimpleSeg (clockwise ordering required for valid polygons) and the survey in
  `docs/research/notes/Phase 1 stack and algorithms/output_format.md` §1, §2, recommendation 3.

Decision (plan D10, as the lab built the targets):
- screen coordinates, y down; `polycore.signed_area2 > 0` = clockwise on screen. Outer rings clockwise (> 0),
  holes counter-clockwise (< 0). Numerically this is the GeoJSON rule (exterior has positive shoelace area);
  drawn with y down it looks clockwise, which is the PolyFormer / VectorLLM / Florence-2 convention.
- start at the top-most, then left-most vertex (minimum (y, x), row-major). PolyFormer and VectorLLM use the
  vertex nearest the image corner (Euclidean); row-major is the repo's `polygonize` and lab convention, is exact
  on the integer grid and has no distance ties ((3, 4) vs (4, 3)).
- parts and holes sorted by the (y, x) sequence of their canonical ring (parts: outer, then holes). For legal
  shapes the starts differ, so this is the lab's sort by start vertex; the full sequence only breaks ties.
- made total so canonicalize is a true canonical form on any input (needed for the invariance property): a
  zero-area ring has no orientation, so both directions are candidates; a start vertex that occurs more than once
  gives several candidates; the smallest (y, x) sequence wins.
- signed area comes from `polycore.signed_area2` (one implementation per operation); tests check it against an
  independent Python-int shoelace.

## Empty object

`to_mask([], W, H)` is an all-false mask with `legal=True`, `defects={}`: "no object" (`"nodule":[]`) is a
well-formed, geometrically legal answer. Whether it is right is the reward's topology and IoU terms' job;
zeroing it through the legality gate would push the policy to always draw something. polycore handles zero rings
itself, so `to_mask` has no special case.

## Tests (written first, run red, then green)

- `tests/test_geometry.py` (text + canonical form): Hypothesis round trips (`parse∘serialize` exact on any
  objects; `serialize∘parse` identity on canonical objects); 39 parametrized rejections, each with its code
  (whitespace, leading zeros, > 1000, negative, floats, booleans, strings, Unicode digits, missing / reordered /
  extra / duplicate `holes`, wrong / missing / extra / duplicate names, < 3 points, wrong arity, trailing text,
  incomplete, escapes, 100 k nested brackets, code fences); non-str input; truncated flag; empty string; caps at
  exactly the limit and one beyond; 1,000-vertex ring parse < 2 ms (measured ~0.3 ms); four ~200 k-character
  near-miss inputs < 100 ms each (linear); the regex accepts every serialization and has no non-portable syntax;
  canonical-form properties (orientation vs an independent int shoelace, top-left start, sorted parts and holes);
  idempotence; invariance to start rotation, ring reversal, part and hole order (500 Hypothesis cases); tie
  breaks (zero-area, figure-eight, repeated top vertex, parts sharing a start); input not modified; four lab
  builder hand cases (one blob, three parts, two holes, RIGA disc + cup) are canonical and come back from 20
  random perturbations each.
- `tests/test_geometry_oracle.py`: differential test against the lab recursive-descent parser on 3,000
  near-miss texts (serializations with up to 3 random edits, biased to land after digits) — same accept/reject
  and same values; `st.text()` fuzz, never raises. It caught the Unicode-digit bug above (the original lab
  oracle fails on its first such example).
- `tests/test_geometry_mask.py`: hand case (square with a hole on a 4×4 crop, pixel centres); mask shape (H, W);
  empty part list; roles (two parts with swapped holes give `hole_outside_outer` = 2); winding is a legality
  defect but not a mask change; bow-tie is illegal but rasterized; Hypothesis: `to_mask` equals
  `polycore.rasterize(...)["evenodd"]` and `polycore.validate(roles)` totals, and canonicalization keeps legal
  masks.

Run: 76 tests green in ~14 s (against the lab polycore build first, then M1's module).

## Gate: R1 parser part

`scratchpad/impl/M2/r1_parser.py`, all fixture items (RIGA + LIDC experts), lab builder in both modes, every
rater, `pad_div` crop, 3 processes, 36 s:

| check | result (6,338 targets) |
|---|---|
| `serialize` text equals the lab's `prod.serialize` text | 100 % |
| `parse` ok | 100 % |
| `serialize(parse(text)) == text` | 100 % |
| lab target is already canonical (`canonicalize` returns it) | 100 % |
| reversed rings + rotated starts + reversed part/hole order canonicalize back to the same text | 100 % |
| `to_mask` legal on every object | 100 % |

Vertices per target: median 17, p99 188, max 239. **R1 parser part: PASS.**

## Status of to_mask (M1 dependency)

`to_mask` and `canonicalize` call `objdel.checker.polycore` at the contract names (`pack`,
`validate(shape, roles=, G=)["total"]`, `rasterize(shape, W, H, G=)["evenodd"]`, `DEFECTS`, `signed_area2`).
They were first developed against the lab polycore build (`lab/t_polycore`, built in the scratch dir and mapped to
`objdel.checker.polycore` by a dev-only pytest plugin, never in the repo). After M1 landed:

- the 76 geometry tests pass unchanged against `src/objdel/checker/polycore.py` (14 s);
- the R1 parser run above was repeated with M1's polycore: 6,338 / 6,338 on every check (parse median 19 µs,
  p99 106 µs; 24 s with 3 processes);
- `make test`: 238 passed, 7 skipped, 6 deselected (slow); `make lint`: clean.

## Open issues

- A 100,000-digit integer reports `syntax`, not `range` (Python's int-string limit fires inside json.loads).
- Inputs far beyond any completion length (600 kB) take ~40 ms to reject; a length cap was not added.
- Constrained decoding with `pattern(names)` is available but not adopted (M7 / Phase 2 decision).

# Output format: polygon JSON on the 0-1000 grid

The model writes one JSON text per image. SFT targets, the reward and evaluation all use this one format.
Source: `src/objdel/geometry/text.py` (grammar, parser), `shape.py` (canonical form), `cuda/gpucheck.h` (legality).
Rationale: [ADR-0010](../architecture/decisions/0010-compact-json-grammar-and-strict-parser.md) and [ADR-0003](../architecture/decisions/0003-exact-integer-geometry.md).

## 1. The prompt

The user turn has one image and one of these texts (`objdel.train.dataset.PROMPTS`).

LIDC:

```text
Outline the lung nodule in this CT crop. Answer with JSON only: {"nodule":[{"outer":[[x,y],...],"holes":[[[x,y],...]]}]}, integers 0-1000 relative to the crop, outer rings clockwise, holes counter-clockwise.
```

RIGA:

```text
Outline the optic disc and the optic cup in this fundus crop. Answer with JSON only: {"disc":[{"outer":[[x,y],...],"holes":[]}],"cup":[{"outer":[[x,y],...],"holes":[]}]}, integers 0-1000 relative to the crop, outer rings clockwise.
```

The schemas are `SCHEMAS = {"riga": ("disc", "cup"), "lidc": ("nodule",)}`.
Names must appear in schema order, each once.

## 2. Grammar

```text
doc   := '{' '"' name '":' parts (',' '"' name '":' parts)* '}'
parts := '[' ']' | '[' part (',' part)* ']'
part  := '{"outer":' ring ',"holes":' rings '}'        "holes" is always present
rings := '[' ']' | '[' ring (',' ring)* ']'
ring  := '[' pt ',' pt ',' pt (',' pt)* ']'            at least 3 points
pt    := '[' int ',' int ']'
int   := 0 | [1-9][0-9]{0,2} | 1000
```

Rules of the text:

- No whitespace inside. Whitespace around the whole text is stripped. Only JSON whitespace counts: space, tab, LF, CR.
- Every text of the grammar is also compact JSON: `json.dumps(obj, separators=(",", ":"))`.
- Keys inside a part are in the order `outer`, then `holes`.
- Integers have no leading zero, no sign, no fraction and no exponent. The range is 0 to 1000.
- Caps: at most 1,000 points per ring (`MAX_VERTICES`), at most 10 rings per object (`MAX_RINGS`, counting each outer ring and each hole), at most 16,384 characters (`MAX_CHARS`).

`pattern(names)` returns the grammar as a regex with no Python-only syntax. It could be passed to a regex-constrained decoder. No config uses constrained decoding today.

## 3. Parser rules (`parse(text, names, truncated=False) -> Parsed`)

`Parsed` has three fields: `ok`, `objs` (`{name: [(outer int64 (k, 2), [hole, ...]), ...]}` or `None`) and `error`.
`parse` never raises. The order of the checks is:

1. `truncated` is true: `truncated`.
2. `text` is not a `str`: `not_text`.
3. More than `MAX_CHARS` characters: `too_long`. This check is first so that parse time is bounded by the cap.
4. Strip JSON whitespace. `re.fullmatch` against the grammar. On a match, build lists with `json.loads`. Then check the ring cap and the vertex cap.
5. On a mismatch, `_diagnose` walks the text once and returns the first problem found.

A completion that was cut at the token limit is a format failure even if its prefix parses.

### Error codes

| Code | Meaning | Example |
|---|---|---|
| `truncated` | The caller says the answer was cut | any text with `truncated=True` |
| `not_text` | Not a string | `5` |
| `too_long` | Over 16,384 characters | 20,000 characters |
| `empty` | Empty after stripping | `""` |
| `not_ascii` | Non-ASCII character | |
| `whitespace` | Whitespace inside the text | `{"nodule": []}` (space after the colon) |
| `leading_zero` | A number starts with 0 and has more digits | `[0,00]` |
| `incomplete` | Valid so far, but the text ends too early | a missing final `]` |
| `trailing_text` | Valid JSON followed by more text | |
| `syntax` | Not valid JSON, or valid JSON outside the grammar's spelling | escapes, `-0` |
| `keys` | Top-level keys are not exactly the schema names in order, or the top level is not an object | `[1]`, or `{"nodule":[],"cup":[]}` for the schema `("nodule",)` |
| `part_keys` | A part does not have exactly `outer`, `holes` in that order | missing `holes`, or `holes` first |
| `schema` | Wrong container type | |
| `not_int` | A coordinate is not an integer | `10.5` |
| `range` | A coordinate is outside 0 to 1000 | `1001` |
| `few_points` | A ring has fewer than 3 points | |
| `cap_rings` | More than 10 rings in one object | |
| `cap_vertices` | More than 1,000 points in a ring | |

Real parser results (`parse(..., ("nodule",))`):

| Text | `ok` | `error` |
|---|---|---|
| `{"nodule":[]}` | true | |
| ` {"nodule":[]} ` (with outer spaces) | true | |
| `{"nodule": []}` | false | `whitespace` |
| `{"nodule":[{"outer":[[0,0],[10,0],[10,10]]}]}` | false | `part_keys` |
| `{"nodule":[{"outer":[[0,0],[10.5,0],[10,10]],"holes":[]}]}` | false | `not_int` |
| `{"nodule":[{"outer":[[0,0],[10,0],[10,10]],"holes":[]}` | false | `incomplete` |

`{"nodule":[]}` is a well-formed answer: "no object". It is a legal empty mask. The reward's topology terms judge it.

## 4. Coordinates

- A frame is W x H pixels. LIDC and RIGA frames are square, with W = H = S (`S` divides 1000, at most 500).
- Grid point (x, y) is pixel point (x * W / 1000, y * H / 1000). One pixel is 1000 / S grid units. For S = 100, a pixel is 10 grid units.
- Screen coordinates: x to the right, y down.
- A pixel (c, r) is inside the outline when its centre (c + 1/2, r + 1/2) is inside under the even-odd rule with a half-open crossing rule. A centre on a left or top edge is inside. A centre on a right or bottom edge is outside.

## 5. Legality rules

An answer is legal when it parses and the checker finds none of the 13 defects in any ring of any object.
The defects are listed with their meaning in [CUDA engine, section 2](cuda-engine.md#the-13-defect-classes).
In short:

| Rule | Defect when broken |
|---|---|
| Each ring has at least 3 distinct points and non-zero area | `fewer_than_3_points`, `zero_area` |
| A ring does not cross, touch or overlap itself, and has no spike | `self_crossing`, `vertex_pinch`, `collinear_overlap`, `spike` |
| Rings of one object do not cross or touch | `ring_ring_crossing`, `ring_ring_touching` |
| A hole lies inside its outer ring | `hole_outside_outer` |
| Sibling holes and sibling parts do not overlap or nest | `nested_holes`, `nested_outers` |
| Outer rings run clockwise on screen (positive shoelace sum), holes counter-clockwise | `wrong_winding` |
| Every coordinate is within the grid, 0 to 1000 | `out_of_range` |

The shoelace sum is `sum(x_i * y_(i+1) - x_(i+1) * y_i)` with y down. A positive sum means clockwise on screen.
Hole roles come from nesting in the text: the holes of a part belong to that part's outer ring.

## 6. Canonical form (`canonicalize`)

SFT targets are written in canonical form so that the model learns one spelling.

1. An outer ring runs clockwise. A hole runs counter-clockwise. A ring with zero area has no orientation, so both directions are tried.
2. A ring starts at its top-most, then left-most vertex (minimum `(y, x)`). If that vertex occurs more than once, every occurrence is a candidate. The ring with the smallest `(y, x)` vertex sequence wins.
3. Holes of a part, and parts of an object, are sorted by the `(y, x)` sequence of their canonical ring.
4. The function copies. It never modifies its input.

Targets are built with a tolerance ladder and certified by the checker: legal, at most 1,000 points per ring, the same signature as the rater's mask under the size rule, and IoU of at least 0.99 with the rater's mask. For LIDC nodules the target must also match the mask pixel for pixel.

## 7. Examples

All examples are real output of `serialize(canonicalize(...))` on a 100 x 100 frame.

**One part without a hole** (a square):

```json
{"nodule":[{"outer":[[300,300],[700,300],[700,700],[300,700]],"holes":[]}]}
```

**One part with one hole.** The ring order of the hole is counter-clockwise on screen. The checker finds it legal. Its mask has 3,200 pixels: 60 x 60 minus 20 x 20.

```json
{"nodule":[{"outer":[[200,200],[800,200],[800,800],[200,800]],"holes":[[[400,400],[400,600],[600,600],[600,400]]]}]}
```

**Two parts:**

```json
{"nodule":[{"outer":[[100,100],[300,100],[300,300],[100,300]],"holes":[]},{"outer":[[600,600],[900,600],[900,900],[600,900]],"holes":[]}]}
```

**RIGA, disc and cup** (schema `("disc", "cup")`):

```json
{"disc":[{"outer":[[500,200],[800,500],[500,800],[200,500]],"holes":[]}],"cup":[{"outer":[[500,400],[600,500],[500,600],[400,500]],"holes":[]}]}
```

**Nothing there:**

```json
{"nodule":[]}
```

**Illegal examples** (scored with `score_completion`; the frame is 100 x 100):

| Text | `format_ok` | `legal` | `defects` |
|---|---|---|---|
| Square written counter-clockwise: `{"nodule":[{"outer":[[300,300],[300,700],[700,700],[700,300]],"holes":[]}]}` | 1 | 0 | `wrong_winding: 1` |
| Bow-tie: `{"nodule":[{"outer":[[300,300],[700,700],[700,300],[300,700]],"holes":[]}]}` | 1 | 0 | `zero_area: 1`, `self_crossing: 1` |
| Coordinate 1001: `{"nodule":[{"outer":[[0,0],[1001,0],[5,5]],"holes":[]}]}` | 0 | 0 | none. Parse error `range` |

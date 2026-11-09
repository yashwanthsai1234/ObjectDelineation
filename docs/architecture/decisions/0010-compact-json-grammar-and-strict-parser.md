# ADR-0010: The model writes compact JSON on a 0-1000 grid, and a strict parser reads it

**Status:** Accepted
**Date:** 2026-10-19 (commit `1c54bc4`, `src/objdel/geometry/text.py`; format rationale `experiments/phase2/PLAN.md` D2)
**Deciders:** Repository owner

## Context

The model must write polygons as text. The format affects three things: how many tokens an answer costs, how hard it is to write legally, and how easily the reward can be fooled.
The same format is used by SFT targets, by the reward and by evaluation. Two readers of one format would disagree on edge cases.

## Decision

- The output is one JSON object. Keys are the structure names in schema order (`nodule`, or `disc` then `cup`). Each value is a list of parts. Each part is `{"outer": ring, "holes": [rings]}`. `holes` is always present.
- Coordinates are integers from 0 to 1000 (`GRID`). Outer rings run clockwise on screen. Holes run counter-clockwise. A ring starts at its top-most, then left-most vertex.
- No whitespace inside the text. Surrounding JSON whitespace is stripped.
- Caps: at most 1,000 points per ring, at most 10 rings per object, at most 16,384 characters in all.
- `objdel.geometry.text.parse` is strict and never raises. It is one regular expression followed by `json.loads`. A rejected text gets a short stable error code from `_diagnose`.
- A completion cut at the token limit is a format failure even if its prefix parses.
- `objdel.geometry.canonicalize` is the single canonical form. SFT targets are serialised from it.

## Options considered

### Option A: Compact JSON, absolute integers on a 1000 grid (chosen)

| Dimension | Assessment |
|---|---|
| Tokens | Qwen3-VL reads about 1.7 characters per token on this grammar. The longest of 30,843 targets has 2,864 characters |
| Strictness | One regex decides acceptance. It is deterministic, so a failed match backtracks in linear time |
| Portability | The pattern has no lookarounds or atomic groups, so it can go to vLLM or xgrammar regex decoding |

### Option B: Delta coordinates, pixel units, coarse grids, coordinate tokens, random start

`experiments/phase2/PLAN.md` D2 lists these as rejected, with the reasons measured in lab A1.
Delta cut 25 % of tokens, but a one-pixel slip cost 13 IoU points on LIDC. Pixel units cut 23 % of LIDC tokens, but cannot express some bevel vertices (0.34 to 0.47 % of targets).
A coarse grid cost IoU and legality (G = 100 on RIGA: -2.7 IoU points and -2.2 points legality). A random start cost 40 IoU points.

### Option C: A hand-written recursive-descent parser

The lab parser. The regex parser is 5 to 6 times faster on 50- to 1,000-vertex answers (`experiments/phase1/impl/M2_geometry.md`).
A differential test keeps the two equal (`tests/test_geometry_oracle.py`, slow).

## Trade-off analysis

Strictness gives the policy no room to hide format errors. A strict parser also makes error codes a useful training signal: the trainer logs `error/<code>` shares per batch.
The cost is token length. The model also needs to learn winding and start-point rules, which the legality checker enforces.

## Consequences

- Easier: one spec for targets, reward and evaluation. `parse` bounds its time by the character cap, not by the input size.
- Easier: the typed error codes drive logs and the hack audit.
- Harder: the model must write rings clockwise from the top-left. A wrong direction is a defect (`wrong_winding`), not a format error.
- Harder: a change to the grammar changes every SFT target, every target-set hash and the prompt text. It needs a new prepared set.
- Not used: constrained decoding. The regex can feed TRL's `vllm_structured_outputs_regex`, but no config sets it.

## Action items

1. [x] `parse`, `serialize`, `canonicalize`, `pattern`.
2. [x] Differential test against the lab parser.
3. [ ] Decide if constrained decoding should remove format failures. Legal answers are already 99.4 % on the test split (`final_eval_2026-11-08.md`).

"""objdel.geometry text spec: serializer, strict parser, canonical form (plan D10, R1 parser part)."""
import json
import re
import time

import numpy as np
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

from objdel.geometry import GRID, MAX_CHARS, MAX_RINGS, MAX_VERTICES, SCHEMAS, canonicalize, parse, pattern, serialize

RIGA, LIDC = SCHEMAS["riga"], SCHEMAS["lidc"]

# Canonical texts produced by the lab target builder (experiments/phase1/lab/t_targets, BuildConfig() and the
# lossless mode give the same text here) for hand-made masks of 25 x 20 px, serialized in the lab grammar.
LAB_TARGETS = [
    (LIDC, '{"nodule":[{"outer":[[240,200],[600,200],[600,400],[800,400],[800,800],[400,800],[400,600],[240,600]],'
           '"holes":[]}]}'),
    (LIDC, '{"nodule":[{"outer":[[120,100],[320,100],[320,350],[120,350]],"holes":[]},{"outer":[[640,150],'
           '[760,150],[760,300],[640,300]],"holes":[]},{"outer":[[560,500],[880,500],[880,850],[560,850]],"holes":[]}]}'),
    (LIDC, '{"nodule":[{"outer":[[160,150],[840,150],[840,850],[160,850]],"holes":[[[320,350],[320,550],[480,550],'
           '[480,350]],[[560,600],[560,750],[720,750],[720,600]]]}]}'),
    (RIGA, '{"disc":[{"outer":[[480,100],[560,150],[640,150],[800,250],[880,450],[880,600],[800,800],[640,900],'
           '[560,900],[480,950],[480,900],[360,900],[200,800],[120,600],[120,450],[200,250],[360,150],[480,150]],'
           '"holes":[]}],"cup":[{"outer":[[480,250],[480,300],[520,300],[640,400],[640,500],[560,650],[480,650],'
           '[480,700],[400,650],[280,550],[280,400],[400,300]],"holes":[]}]}'),
]


def area2(ring):
    """Independent oracle: twice the signed area with Python ints (> 0 = clockwise on screen, y down)."""
    pts = [(int(x), int(y)) for x, y in ring]
    return sum(x0 * y1 - x1 * y0 for (x0, y0), (x1, y1) in zip(pts, pts[1:] + pts[:1], strict=True))


def same(a, b):
    """Objects equal (serialize is injective on integer rings) and every ring is int64."""
    rings = [r for parts in a.values() for o, hs in parts for r in (o, *hs)]
    return serialize(a, list(a)) == serialize(b, list(b)) and all(r.dtype == np.int64 for r in rings)


# ------------------------------------------------------------------ strategies

coord = st.integers(0, GRID)
ring = st.lists(st.tuples(coord, coord), min_size=3, max_size=30).map(lambda p: np.array(p, np.int64))
part = st.tuples(ring, st.lists(ring, max_size=3))
parts = st.lists(part, max_size=3).filter(lambda ps: sum(1 + len(hs) for _, hs in ps) <= MAX_RINGS)
schema = st.sampled_from([RIGA, LIDC])


@st.composite
def objects(draw):
    names = draw(schema)
    return names, {n: draw(parts) for n in names}


@st.composite
def perturbed(draw):
    """Objects and a copy with rotated starts, reversed rings, shuffled parts and holes."""
    names, objs = draw(objects())
    rnd = draw(st.randoms(use_true_random=False))

    def move(r):
        r = np.roll(r, -rnd.randrange(len(r)), axis=0)
        return r[::-1].copy() if rnd.random() < 0.5 else r

    out = {}
    for n in names:
        ps = [(move(o), rnd.sample([move(h) for h in hs], len(hs))) for o, hs in objs[n]]
        out[n] = rnd.sample(ps, len(ps))
    return names, objs, out


# ------------------------------------------------------------------ serializer + parser

@settings(max_examples=300, deadline=None)
@given(objects())
def test_parse_inverts_serialize(case):
    names, objs = case
    text = serialize(objs, names)
    p = parse(text, names)
    assert p.ok and p.error == "" and same(p.objs, objs)
    assert text.count('"outer"') == text.count('"holes"')  # "holes" always present
    assert " " not in text


@settings(max_examples=300, deadline=None)
@given(objects())
def test_serialize_inverts_parse_on_canonical_objects(case):
    names, objs = case
    text = serialize(canonicalize(objs), names)
    p = parse(text, names)
    assert p.ok and serialize(p.objs, names) == text


def test_empty_object_lists():
    assert serialize({"nodule": []}, LIDC) == '{"nodule":[]}'
    assert parse('{"disc":[],"cup":[]}', RIGA) == (True, {"disc": [], "cup": []}, "")


def test_surrounding_whitespace_is_stripped():
    assert parse('\n {"nodule":[]}\n', LIDC).ok
    assert parse(' \t\r\n{"nodule":[]}\r\n\t ', LIDC).ok  # exactly JSON whitespace (RFC 8259 section 2)


@pytest.mark.parametrize("ws,code", [("\u3000", "not_ascii"), ("\u2028", "not_ascii"), ("\xa0", "not_ascii"),
                                     ("\x85", "not_ascii"), ("\x0b", "whitespace"), ("\x0c", "whitespace"),
                                     ("\x1c", "whitespace"), ("\x1f", "whitespace"), ("\x00", "")])
def test_only_json_whitespace_is_stripped(ws, code):
    """B1-1: str.strip() also strips Unicode spaces and the separators U+001C-U+001F."""
    assert parse(ws + '{"nodule":[]}', LIDC) == (False, None, code or "syntax")
    assert parse('{"nodule":[]}' + ws, LIDC) == (False, None, code or "trailing_text")


def test_length_cap():
    """I4: texts longer than MAX_CHARS are refused before any parsing; the cap is far above any completion budget
    (1,536 tokens at <= 2.14 characters per token for the three tokenizers is <= 3,300 characters)."""
    assert MAX_CHARS >= 4 * 1536 * 2.14
    doc = '{"nodule":[]}'
    assert parse(doc + " " * (MAX_CHARS - len(doc)), LIDC).ok
    assert parse(doc + " " * (MAX_CHARS - len(doc) + 1), LIDC) == (False, None, "too_long")
    assert parse("[" * 10 ** 6, LIDC) == (False, None, "too_long")
    assert _best_ms(lambda: parse("x" * 10 ** 7, LIDC)) < 1


@pytest.mark.parametrize("digits", [5, 4300, 4301, 10_000])
def test_huge_integers_are_out_of_range(digits):
    """I4: Python's int-string limit (4,300 digits) used to surface as 'syntax'; any integer of 5+ digits is out of
    range, as is a negative one."""
    for v in ("9" * digits, "-" + "9" * digits, "1" + "0" * (digits - 1)):
        assert parse(OK.replace("[5,5]", f"[{v},5]"), LIDC) == (False, None, "range")


RING = "[[1,1],[5,1],[5,5]]"
PART = '{"outer":' + RING + ',"holes":[]}'
OK = '{"nodule":[' + PART + "]}"


@pytest.mark.parametrize("text,names,code", [
    ("", LIDC, "empty"),
    ("   \n", LIDC, "empty"),
    (OK.replace(",[5,1]", ", [5,1]"), LIDC, "whitespace"),
    (OK.replace("[1,1]", "[1,\n1]"), LIDC, "whitespace"),
    (OK.replace("[5,5]", "[05,5]"), LIDC, "leading_zero"),
    (OK.replace("[5,5]", "[5,00]"), LIDC, "leading_zero"),
    (OK.replace("[5,5]", "[1001,5]"), LIDC, "range"),
    (OK.replace("[5,5]", "[5,99999]"), LIDC, "range"),
    (OK.replace("[5,5]", "[-5,5]"), LIDC, "range"),
    (OK.replace("[5,5]", "[5.0,5]"), LIDC, "not_int"),
    (OK.replace("[5,5]", "[5e2,5]"), LIDC, "not_int"),
    (OK.replace("[5,5]", "[true,5]"), LIDC, "not_int"),
    (OK.replace("[5,5]", '["5",5]'), LIDC, "not_int"),
    (OK.replace("[5,5]", "[５,5]"), LIDC, "not_ascii"),
    (OK.replace("[5,5]", "[٥,5]"), LIDC, "not_ascii"),
    (OK.replace(',"holes":[]', ""), LIDC, "part_keys"),
    (OK.replace('{"outer":' + RING + ',"holes":[]}', '{"holes":[],"outer":' + RING + "}"), LIDC, "part_keys"),
    (OK.replace(',"holes":[]', ',"holes":[],"extra":[]'), LIDC, "part_keys"),
    (OK.replace(',"holes":[]', ',"holes":[],"holes":[]'), LIDC, "part_keys"),
    ('{"cup":[],"disc":[]}', RIGA, "keys"),
    ('{"disc":[]}', RIGA, "keys"),
    ('{"disc":[],"cup":[],"rim":[]}', RIGA, "keys"),
    ('{"nodule":[],"nodule":[]}', LIDC, "keys"),
    ('{"Nodule":[]}', LIDC, "keys"),
    ('[]', LIDC, "keys"),
    (OK.replace(RING, "[[1,1],[5,1]]"), LIDC, "few_points"),
    (OK.replace(',"holes":[]', ',"holes":[[[1,1],[2,2]]]'), LIDC, "few_points"),
    (OK.replace(RING, "[[1,1],[5,1],[5]]"), LIDC, "schema"),
    (OK.replace(RING, "[[1,1],[5,1],[5,5,5]]"), LIDC, "schema"),
    (OK.replace(',"holes":[]', ',"holes":{}'), LIDC, "schema"),
    ('{"nodule":{}}', LIDC, "schema"),
    (OK + "x", LIDC, "trailing_text"),
    (OK + OK, LIDC, "trailing_text"),
    (OK[:-2], LIDC, "incomplete"),
    (OK[:20], LIDC, "incomplete"),
    ('{"nodule":[' + PART + ",]}", LIDC, "syntax"),
    ('{"nod\\u0075le":[]}', LIDC, "syntax"),
    ("[" * 10_000, LIDC, "syntax"),
    ("```json" + OK + "```", LIDC, "syntax"),
])
def test_rejections(text, names, code):
    p = parse(text, names)
    assert p == (False, None, code)


@pytest.mark.parametrize("value", [None, b'{"nodule":[]}', 3, ["x"]])
def test_non_text_is_rejected(value):
    assert parse(value, LIDC) == (False, None, "not_text")


def test_truncated_completion_fails_even_if_it_parses():
    assert parse(OK, LIDC, truncated=True) == (False, None, "truncated")


def test_caps():
    pts = ",".join(f"[{i % 1000},{i // 1000}]" for i in range(MAX_VERTICES))
    assert parse('{"nodule":[{"outer":[' + pts + '],"holes":[]}]}', LIDC).ok
    assert parse('{"nodule":[{"outer":[' + pts + ',[7,7]],"holes":[]}]}', LIDC).error == "cap_vertices"
    nine = '{"outer":' + RING + ',"holes":[' + ",".join([RING] * (MAX_RINGS - 2)) + "]}"
    assert parse('{"nodule":[' + nine + "," + PART + "]}", LIDC).ok
    assert parse('{"nodule":[' + nine + "," + PART + "," + PART + "]}", LIDC).error == "cap_rings"
    # the cap is per object: 10 rings in each of disc and cup is fine
    assert parse('{"disc":[' + nine + "," + PART + '],"cup":[' + nine + "," + PART + "]}", RIGA).ok


def _best_ms(f, repeat=5):
    best = float("inf")
    for _ in range(repeat):
        t = time.perf_counter()
        f()
        best = min(best, time.perf_counter() - t)
    return best * 1e3


def test_parse_speed_on_a_1000_vertex_ring():
    t = np.linspace(0, 2 * np.pi, MAX_VERTICES, endpoint=False)
    ring = np.rint(np.stack([500 + 400 * np.cos(t), 500 + 400 * np.sin(t)], 1)).astype(np.int64)
    text = serialize({"nodule": [(ring, [])]}, LIDC)
    assert parse(text, LIDC).ok
    assert _best_ms(lambda: parse(text, LIDC)) < 2.0  # measured ~0.25 ms on an Apple M3


@pytest.mark.parametrize("make", [
    lambda m: '{"disc":[' + ",".join([PART] * m) + '],"cup":[]}x',              # fails at the very end
    lambda m: '{"disc":[{"outer":[' + "[1,1]," * (20 * m),                        # unterminated, no caps hit
    lambda m: '{"disc":[{"outer":' + RING + ',"holes":[' + ",".join([RING] * m) + ']}],"cup":[]}x',
    lambda m: '{"disc":[' + ",".join(['{"outer":[[1,1],[5,1],[5,05]],"holes":[]}'] * m) + '],"cup":[]}',
])
def test_long_adversarial_inputs_are_linear(make):
    """Near-miss texts at the length cap: a quadratic parser or a backtracking regex would take far longer."""
    m = 1
    while len(make(2 * m)) <= MAX_CHARS:
        m *= 2
    text = make(m)
    assert len(text) > MAX_CHARS // 4 and not parse(text, RIGA).ok and parse(text, RIGA).error != "too_long"
    assert _best_ms(lambda: parse(text, RIGA), repeat=3) < 20


def test_pattern_is_the_grammar():
    """The exported regex (for constrained decoding) accepts exactly what the parser accepts, caps aside."""
    rx = re.compile(pattern(RIGA))
    assert rx.fullmatch(serialize({"disc": [], "cup": []}, RIGA))
    assert not rx.fullmatch(OK)
    assert pattern(LIDC) == pattern(list(LIDC))
    for bad in ("(?>", "*+", "++", "\\d", "(?=", "(?!", "(?<"):  # portable: no Python-only or lookaround syntax
        assert bad not in pattern(RIGA)


@settings(max_examples=300, deadline=None)
@given(objects())
def test_pattern_accepts_every_serialization(case):
    names, objs = case
    assert re.fullmatch(pattern(names), serialize(objs, names))


# ------------------------------------------------------------------ canonical form

@settings(max_examples=300, deadline=None)
@given(objects())
def test_canonical_form_properties(case):
    names, objs = case
    can = canonicalize(objs)
    assert list(can) == list(names)
    for n in names:
        for o, hs in can[n]:
            for r, hole in [(o, False), *[(h, True) for h in hs]]:
                a = area2(r)
                assert a == 0 or (a > 0) != hole  # outer clockwise on screen, holes counter-clockwise
                assert tuple(r[0][::-1]) == min(tuple(p[::-1]) for p in r.tolist())  # top-most, then left-most
            starts = [tuple(h[0][::-1]) for h in hs]
            assert starts == sorted(starts)
        starts = [tuple(o[0][::-1]) for o, _ in can[n]]
        assert starts == sorted(starts)


@settings(max_examples=300, deadline=None)
@given(objects())
def test_canonicalize_is_idempotent(case):
    names, objs = case
    once = canonicalize(objs)
    assert same(canonicalize(once), once)


@settings(max_examples=500, deadline=None)
@given(perturbed())
def test_canonicalize_ignores_start_direction_and_order(case):
    names, objs, moved = case
    assert serialize(canonicalize(moved), names) == serialize(canonicalize(objs), names)


def test_canonical_tie_breaks():
    # degenerate ring (zero area): either direction; ring visiting its top-left vertex twice
    flat = np.array([[0, 0], [5, 0], [9, 0]])
    eight = np.array([[0, 0], [4, 4], [8, 0], [8, 8], [4, 4], [0, 8]])  # zero net area, pinched
    pinch = np.array([[2, 0], [4, 2], [2, 0], [0, 2], [2, 4]])          # top vertex visited twice
    for r in (flat, eight, pinch):
        texts = {serialize(canonicalize({"nodule": [(np.roll(v, k, 0), [])]}), LIDC)
                 for v in (r, r[::-1]) for k in range(len(r))}
        assert len(texts) == 1
    # parts with the same start vertex are ordered by the rest of the ring, then by their holes
    a = (np.array([[0, 0], [9, 0], [9, 9]]), [])
    b = (np.array([[0, 0], [9, 0], [9, 5]]), [])
    c = (np.array([[0, 0], [9, 0], [9, 5]]), [np.array([[5, 1], [6, 2], [6, 1]])])
    assert serialize(canonicalize({"nodule": [c, a, b]}), LIDC) == serialize({"nodule": [b, c, a]}, LIDC)


def test_canonicalize_does_not_modify_its_input():
    r = np.array([[5, 5], [0, 0], [9, 0]])
    objs = {"nodule": [(r, [])]}
    canonicalize(objs)
    assert r.tolist() == [[5, 5], [0, 0], [9, 0]] and objs["nodule"][0][0] is r


@pytest.mark.parametrize("names,text", LAB_TARGETS)
def test_lab_targets_are_canonical_and_round_trip(names, text):
    """R1 parser part on lab targets: parse ok, canonical, and canonicalization of any perturbation gives it back."""
    p = parse(text, names)
    assert p.ok and serialize(p.objs, names) == text
    assert serialize(canonicalize(p.objs), names) == text
    rng = np.random.default_rng(0)
    for _ in range(20):
        moved = {}
        for n in names:
            ps = []
            for o, hs in p.objs[n]:
                rings = [np.roll(r[::-1] if rng.random() < 0.5 else r, -int(rng.integers(len(r))), 0)
                         for r in [o, *hs]]
                ps.append((rings[0], rings[1:][::-1]))
            moved[n] = ps[::-1]
        assert serialize(canonicalize(moved), names) == text


def test_lab_grammar_compatibility():
    """The lab texts are also plain compact JSON (json.dumps with separators)."""
    for _, text in LAB_TARGETS:
        assert json.dumps(json.loads(text), separators=(",", ":")) == text

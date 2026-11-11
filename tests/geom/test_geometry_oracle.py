"""Differential test: objdel.geom.text.parse vs the lab's hand-written recursive-descent parser (slow oracle,
archive/phase1_lab/t_stack/prod.py, trimmed to accept/reject and the parsed values).

One fix to the lab code: it tested digits with str.isdigit and converted with int(), so it accepted Unicode digits
("1\u0661" parsed as 11); the oracle uses ASCII digits, as the spec says."""
import numpy as np
from hypothesis import given, settings
from hypothesis import strategies as st

from objdel.geom.text import GRID, MAX_RINGS, MAX_VERTICES, parse, serialize

SCHEMAS = (("disc", "cup"), ("nodule",))  # a two-structure schema and LIDC's


class _Fail(Exception):
    pass


class _Lab:
    def __init__(self, s):
        self.s, self.i = s, 0

    def lit(self, t):
        if not self.s.startswith(t, self.i):
            raise _Fail
        self.i += len(t)

    def peek(self):
        return self.s[self.i] if self.i < len(self.s) else ""

    def digit(self, k):
        return k < len(self.s) and "0" <= self.s[k] <= "9"

    def int_(self):
        s, j = self.s, self.i
        if j < len(s) and s[j] == "0":
            self.i = j + 1
            if self.digit(self.i):
                raise _Fail
            return 0
        if not (j < len(s) and "1" <= s[j] <= "9"):
            raise _Fail
        k = j
        while self.digit(k):
            k += 1
        if k - j > 4 or int(s[j:k]) > GRID:
            raise _Fail
        self.i = k
        return int(s[j:k])

    def ring(self):
        self.lit("[")
        pts = []
        while True:
            self.lit("[")
            x = self.int_()
            self.lit(",")
            y = self.int_()
            self.lit("]")
            pts.append((x, y))
            if len(pts) > MAX_VERTICES:
                raise _Fail
            if self.peek() != ",":
                break
            self.i += 1
        self.lit("]")
        if len(pts) < 3:
            raise _Fail
        return np.array(pts, np.int64)

    def seq(self, item):
        self.lit("[")
        out = []
        if self.peek() == "]":
            self.i += 1
            return out
        while True:
            out.append(item())
            if self.peek() != ",":
                break
            self.i += 1
        self.lit("]")
        return out

    def part(self):
        self.lit('{"outer":')
        o = self.ring()
        self.lit(',"holes":')
        hs = self.seq(self.ring)
        self.lit("}")
        return o, hs


def lab_parse(text, names):
    p = _Lab(text.strip())
    try:
        p.lit("{")
        objs = {}
        for j, n in enumerate(names):
            if j:
                p.lit(",")
            p.lit(f'"{n}":')
            objs[n] = p.seq(p.part)
            if sum(1 + len(hs) for _, hs in objs[n]) > MAX_RINGS:
                raise _Fail
        p.lit("}")
        if p.i != len(p.s):
            raise _Fail
    except _Fail:
        return None
    return objs


coord = st.integers(0, GRID)
ring = st.lists(st.tuples(coord, coord), min_size=1, max_size=8).map(lambda p: np.array(p, np.int64))
parts = st.lists(st.tuples(ring, st.lists(ring, max_size=4)), max_size=4)
EXTRA = st.sampled_from([" ", "\n", "-", ".", "e", "0", "00", "1000", "1001", "9999", "١", "５", '"', ":", ",", "[",
                         "]", "{", "}", '"holes":[]', '"outer":', "null", "true", "\\u0030"])


@st.composite
def near_miss(draw):
    """A serialized (possibly cap-violating, < 3-point) object set with a few random edits."""
    names = draw(st.sampled_from(SCHEMAS))
    text = serialize({n: draw(parts) for n in names}, names)
    digits = [k + 1 for k, c in enumerate(text) if c.isdigit()]  # edits right after a digit are the subtle ones
    for _ in range(draw(st.integers(0, 3))):
        i = draw(st.integers(0, len(text)) | st.sampled_from(digits)) if digits else draw(st.integers(0, len(text)))
        i = min(i, len(text))
        j = draw(st.integers(i, min(len(text), i + 3)))
        text = text[:i] + draw(st.one_of(st.just(""), EXTRA)) + text[j:]
    return names, text


@settings(max_examples=2000, deadline=None)
@given(near_miss())
def test_parse_equals_lab_parser(case):
    names, text = case
    got, want = parse(text, names), lab_parse(text, names)
    assert got.ok == (want is not None), (text, got.error)
    assert (got.objs is None) == (not got.ok) and (got.error == "") == got.ok
    if got.ok:
        assert serialize(got.objs, names) == serialize(want, names)


def test_unicode_digits_are_rejected():
    text = '{"nodule":[{"outer":[[1\u0661,1],[5,1],[5,5]],"holes":[]}]}'
    assert lab_parse(text, ("nodule",)) is None and parse(text, ("nodule",)).error == "not_ascii"


@settings(max_examples=1000, deadline=None)
@given(st.text(), st.sampled_from(SCHEMAS))
def test_parse_never_raises(text, names):
    p = parse(text, names)
    assert p.ok == (lab_parse(text, names) is not None) and isinstance(p.error, str)

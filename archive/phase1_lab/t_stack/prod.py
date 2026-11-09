"""Production single-spec path: serializer, strict parser, polycore geometry, reward-v2 terms, advantages, TRL rewards.

One spec for SFT targets and RL rewards:
- Text grammar (compact, no whitespace; surrounding whitespace is stripped, nothing else is tolerated):
      doc   := '{' '"' name '":' parts (',' '"' name '":' parts)* '}'   names = schema, in schema order, all present
      parts := '[' ']' | '[' part (',' part)* ']'
      part  := '{"outer":' ring ',"holes":' rings '}'                    "holes" is ALWAYS present
      rings := '[' ']' | '[' ring (',' ring)* ']'
      ring  := '[' pt ',' pt ',' pt (',' pt)* ']'                         >= 3 points
      pt    := '[' int ',' int ']'      int := '0' | [1-9][0-9]{0,3}, value <= 1000
  Caps: <= MAX_VERTS vertices per ring, <= MAX_RINGS rings per object. A completion cut at max_new_tokens
  (no EOS) is a format failure even if the prefix happens to parse.
- Geometry: polycore only (rasterize even-odd at the crop size, G = 1000; validate with roles, per object).
- Terms: rewardcore v2 (ExpertCache, score -> F5 band, signature, IoU_tc). IoU_tc fallback for a topology no
  rater drew is a parameter: "vote" (default, IoU vs strict majority) or "max_all" (iteration-2 behaviour).
- Advantages: reward v2 (legality gate, GDPO with sigma floors, batch z, DM-loo topology mean-centred, c_topo = 2)
  plus "drop low-spread groups": a group whose geometric terms all span < drop_spread and whose rollouts all give
  the same topology answer gets advantage 0 and is excluded from the batch-level z-score.
"""
from __future__ import annotations

import json
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE / "vendor"))  # vendored copies of rewardcore / targets (unedited)

import polycore  # noqa: E402
from rewardcore import advantages as RA  # noqa: E402
from rewardcore.reward import ExpertCache, score as rc_score  # noqa: E402

G_GRID = 1000
MAX_VERTS, MAX_RINGS = 1000, 10
RIGA_NAMES = ("disc", "cup")
LIDC_NAMES = ("nodule",)


# ------------------------------------------------------------------ serializer (SFT targets)

def _ring_txt(r):
    return "[" + ",".join(f"[{int(x)},{int(y)}]" for x, y in np.asarray(r).reshape(-1, 2)) + "]"


def serialize(objs: dict, names) -> str:
    """objs: {name: [(outer (k,2), [holes])]} -> compact text in the production grammar."""
    out = []
    for n in names:
        parts = ",".join('{"outer":' + _ring_txt(o) + ',"holes":[' + ",".join(_ring_txt(h) for h in hs) + "]}"
                         for o, hs in objs[n])
        out.append(f'"{n}":[{parts}]')
    return "{" + ",".join(out) + "}"


# ------------------------------------------------------------------ strict parser

class ParseError(Exception):
    pass


@dataclass
class Parsed:
    ok: bool
    objs: dict | None  # {name: [(outer int64 (k,2), [holes])]}
    error: str = ""


class _P:
    __slots__ = ("s", "i")

    def __init__(self, s):
        self.s, self.i = s, 0

    def lit(self, t):
        if not self.s.startswith(t, self.i):
            raise ParseError(f"expected {t!r} at {self.i}" if self.i < len(self.s) else "truncated")
        self.i += len(t)

    def peek(self):
        return self.s[self.i] if self.i < len(self.s) else ""

    def int_(self):
        j = self.i
        s = self.s
        if j >= len(s):
            raise ParseError("truncated")
        if s[j] == "0":
            self.i = j + 1
            if self.peek().isdigit():
                raise ParseError("leading zero")
            return 0
        if not ("1" <= s[j] <= "9"):
            raise ParseError(f"expected int at {j}")
        k = j
        while k < len(s) and s[k].isdigit():
            k += 1
        if k - j > 4:
            raise ParseError("int too long")
        v = int(s[j:k])
        if v > G_GRID:
            raise ParseError("int > 1000")
        self.i = k
        return v

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
            if len(pts) > MAX_VERTS:
                raise ParseError("cap: vertices")
            if self.peek() == ",":
                self.i += 1
                continue
            self.lit("]")
            break
        if len(pts) < 3:
            raise ParseError("ring < 3 points")
        return np.array(pts, np.int64)

    def rings(self):
        self.lit("[")
        out = []
        if self.peek() == "]":
            self.i += 1
            return out
        while True:
            out.append(self.ring())
            if self.peek() == ",":
                self.i += 1
                continue
            self.lit("]")
            return out

    def parts(self):
        self.lit("[")
        out = []
        if self.peek() == "]":
            self.i += 1
            return out
        while True:
            self.lit('{"outer":')
            o = self.ring()
            self.lit(',"holes":')
            hs = self.rings()
            self.lit("}")
            out.append((o, hs))
            if self.peek() == ",":
                self.i += 1
                continue
            self.lit("]")
            return out


def parse(text: str, names, truncated: bool = False) -> Parsed:
    """Strict production parser. Never raises; failures carry an error code."""
    if truncated:
        return Parsed(False, None, "truncated")
    if not isinstance(text, str):
        return Parsed(False, None, "not_text")
    p = _P(text.strip())
    try:
        p.lit("{")
        objs = {}
        for j, n in enumerate(names):
            if j:
                p.lit(",")
            p.lit(f'"{n}":')
            objs[n] = p.parts()
            if sum(1 + len(hs) for _, hs in objs[n]) > MAX_RINGS:
                raise ParseError("cap: rings")
        p.lit("}")
        if p.i != len(p.s):
            raise ParseError("trailing text")
    except ParseError as e:
        return Parsed(False, None, str(e).split(" at ")[0])
    return Parsed(True, objs)


# ------------------------------------------------------------------ geometry (polycore only)

def shape_of(parts):
    rings, roles = [], []
    for o, hs in parts:
        k = len(rings)
        rings.append(o)
        roles.append(-1)
        for h in hs:
            rings.append(h)
            roles.append(k)
    return rings, np.asarray(roles, np.int32)


def geometry(parts, W, H):
    """-> (mask bool (H, W) even-odd, legal bool, defect totals dict)."""
    if not parts:
        return np.zeros((H, W), bool), True, {}
    rings, roles = shape_of(parts)
    xy, off = polycore.pack(rings)
    v = polycore.validate((xy, off), roles=roles, G=G_GRID)
    tot = {d: int(c) for d, c in zip(polycore.DEFECTS, v["total"]) if c}
    m = polycore.rasterize((xy, off), W, H, G=G_GRID)["evenodd"]
    return np.asarray(m, bool), not tot, tot


# ------------------------------------------------------------------ reward-v2 terms

def build_cache(stack) -> ExpertCache:
    return ExpertCache.build(np.asarray(stack, bool))


def terms(mask, cache: ExpertCache, fallback: str = "vote") -> dict:
    """Per-object reward-v2 terms. fallback: IoU_tc for a topology no rater drew ("vote" | "max_all")."""
    t = rc_score(mask, cache)
    drawn = cache.signatures.get(t.sig, 0) > 0
    iou_tc = t.iou_tc if (drawn or fallback == "max_all") else t.iou
    return {"iou_tc": float(iou_tc), "f5": float(t.band["F5"]), "iou_vote": float(t.iou), "sig": t.sig,
            "topo_any": float(drawn), "topo_freq": float(t.topo_freq)}


def sig_id(sig) -> float:
    return float(hash(tuple(sig)) % 1_000_003)  # tuple-of-int hashes are not salted (as rewardcore.batch)


# ------------------------------------------------------------------ advantages (reference implementation)

W_GEO = np.array([0.3, 0.3, 0.1])
FLOORS = np.array([0.05, 0.05, 0.25])


def advantages(iou_tc, f5, legal, sig_ids, topo_freq, topo_any, group=8, c_topo=2.0, lam=0.5,
               drop_spread: float | None = None):
    """Reward-v2 advantages (identical to rewardcore.batch.advantages_v2) + optional low-spread group drop.
    All inputs are flat (B*G,) arrays ordered prompt by prompt."""
    a = [np.asarray(x, float).reshape(-1, group) for x in (iou_tc, f5, legal, sig_ids, topo_freq, topo_any)]
    iou_tc, f5, legal, sids, freq, tany = a
    B = iou_tc.shape[0]
    keep = np.ones(B, bool)
    if drop_spread is not None:
        geo = np.stack([iou_tc, f5, legal], -1)
        spread = (geo.max(1) - geo.min(1)).max(-1)
        same_topo = (sids == sids[:, :1]).all(1)
        keep = ~((spread < drop_spread) & same_topo)
    out = np.zeros((B, group))
    if keep.any():
        k = keep
        valid = (tany[k] > 0) & (iou_tc[k] >= 0.5)
        r_topo = RA.dm_loo_topology(sids[k], freq[k], valid, lam=lam)
        R = np.stack([iou_tc[k], f5[k], legal[k]], -1)
        adv, _ = RA.reward_v2(R, W_GEO, FLOORS, ("iou", "band", "legal"), r_topo, c_topo=c_topo)
        out[k] = adv
    return out.reshape(-1)


# ------------------------------------------------------------------ TRL reward functions (aggregator + loggers)

LOG_TERMS = ("format_ok", "legal", "iou_tc", "f5", "topo_any", "r_topo")


class TRLReward:
    """Builds TRL reward callables that share one per-batch computation.

    provider(prompts, completions, completion_ids, **columns) -> dict of flat arrays with keys
    iou_tc, f5, legal, sig_id, topo_freq, topo_any, format_ok (one value per completion).
    Loggers (weight 0) return one term each; `aggregator` (weight 1) returns the v2 advantage, so with
    multi_objective_aggregation="sum_then_normalize" and scale_rewards="none" TRL's advantage is
    aggregator - group_mean(aggregator) = aggregator (group means are already 0).
    """

    def __init__(self, provider, group=8, drop_spread=1e-2, c_topo=2.0, lam=0.5):
        self.provider, self.G, self.drop, self.c_topo, self.lam = provider, group, drop_spread, c_topo, lam
        self._key, self._val = None, None
        self.history = []  # every batch: terms + advantages (for tests / audit)

    def _compute(self, prompts, completions, completion_ids, **kw):
        key = (len(completions), tuple(map(str, completions)), tuple(map(str, kw.get("key", []))),
               getattr(kw.get("trainer_state"), "global_step", None))
        if key == self._key:
            return self._val
        T = self.provider(prompts=prompts, completions=completions, completion_ids=completion_ids, **kw)
        T = {k: np.asarray(v, float) for k, v in T.items()}
        g = self.G
        sids = T["sig_id"].reshape(-1, g)
        valid = (T["topo_any"].reshape(-1, g) > 0) & (T["iou_tc"].reshape(-1, g) >= 0.5)
        r_topo = RA.dm_loo_topology(sids, T["topo_freq"].reshape(-1, g), valid, lam=self.lam)
        T["r_topo"] = (r_topo * (T["legal"].reshape(-1, g) >= 0.5)).reshape(-1)
        T["advantage"] = advantages(T["iou_tc"], T["f5"], T["legal"], T["sig_id"], T["topo_freq"], T["topo_any"],
                                    group=g, c_topo=self.c_topo, lam=self.lam, drop_spread=self.drop)
        self._key, self._val = key, T
        self.history.append(T)
        return T

    def funcs(self):
        def make(name):
            def f(prompts, completions, completion_ids=None, **kw):
                return self._compute(prompts, completions, completion_ids, **kw)[name].tolist()
            f.__name__ = name
            return f
        loggers = [make(n) for n in LOG_TERMS]
        agg = make("advantage")
        agg.__name__ = "v2_aggregate"
        return loggers + [agg], [0.0] * len(loggers) + [1.0]


# ------------------------------------------------------------------ real provider (text -> terms)

def score_text(text, truncated, names, caches: dict, W, H, fallback="vote"):
    """One completion -> per-rollout terms (objects combined by mean; legality = all objects legal).
    Unparseable / truncated -> format_ok 0, legal 0, all terms 0 (the legality gate zeroes it anyway)."""
    p = parse(text, names, truncated)
    if not p.ok:
        return {"format_ok": 0.0, "legal": 0.0, "iou_tc": 0.0, "f5": 0.0, "topo_any": 0.0, "topo_freq": 0.0,
                "sig_id": -1.0, "error": p.error}
    legal_all, ts, sigs = True, [], []
    for n in names:
        m, legal, _ = geometry(p.objs[n], W, H)
        legal_all &= legal
        t = terms(m, caches[n], fallback)
        ts.append(t)
        sigs.append(t["sig"])
    # joint (item-level) topology answer: share of raters whose per-object signatures all equal the sample's
    k = caches[names[0]].k
    joint = float(np.mean([all(caches[n].rater_sigs[r] == s for n, s in zip(names, sigs)) for r in range(k)]))
    return {"format_ok": 1.0, "legal": float(legal_all),
            "iou_tc": float(np.mean([t["iou_tc"] for t in ts])), "f5": float(np.mean([t["f5"] for t in ts])),
            "topo_any": float(joint > 0), "topo_freq": joint,
            "sig_id": sig_id(tuple(len(s) for s in sigs) + tuple(x for s in sigs for x in s)), "error": ""}


def json_dumps_compact(o):
    return json.dumps(o, separators=(",", ":"))

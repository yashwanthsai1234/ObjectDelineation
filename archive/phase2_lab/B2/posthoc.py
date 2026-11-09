"""B2 post-hoc (EXPLORATORY, added after seeing analyze.py's results; no verdict attached): where does the group
reward variance come from?

    .venv/bin/python -I -m experiments.phase2.lab.B2.posthoc

Illegal samples score R = 0 (legality gate), so per group of G = 8 with legal share p, legal-sample mean mu_L and
legal-sample variance var_L (population variances), the law of total variance gives exactly
    Var(R) = p (1 - p) mu_L^2   (legality part: which samples are legal)
           + p var_L            (geometric part: how legal samples differ).
GRPO with a legality gate learns legality from the first part and geometry from the second; G-format (legal >= 90 %)
caps the first at 0.09 mu_L^2. Reported per run and view, plus paired depth (150 - 36) and interpolation differences
of the geometric part's standard deviation sqrt(p var_L), with nodule-cluster CIs.
"""

import json

import numpy as np

from . import prereg as P
from .analyze import INTERP, ci, draws
from .score import B2


def parts(v: dict) -> dict[str, np.ndarray]:
    R, L = np.asarray(v["R"], float), np.asarray(v["legal"], float) > 0
    p = L.mean(1)
    n = np.maximum(L.sum(1), 1)
    mu = (R * L).sum(1) / n
    var_l = ((R - mu[:, None]) ** 2 * L).sum(1) / n
    out = {"var": R.var(1), "var_legality": p * (1 - p) * mu**2, "var_geom": p * var_l}
    out["sd_geom"] = np.sqrt(out["var_geom"])
    out["geom_share"] = np.where(out["var"] > 0, out["var_geom"] / np.where(out["var"] > 0, out["var"], 1), 0.0)
    return out


def main():
    runs = sorted({r for p in P.DEPTH_RUNS.values() for r in p})
    runs += [r for r in INTERP.values() if (B2 / "scored" / f"{r}.json").exists()]
    data = {r: json.loads((B2 / "scored" / f"{r}.json").read_text()) for r in runs}
    res = {"tables": {}, "paired": {}}
    for r, d in data.items():
        res["tables"][r] = {vw: {k: ci(x, draws(v["item_id"])) for k, x in parts(v).items()} for vw, v in d.items()}
    pairs = {f"{arm} 150-36": (deep, shallow) for arm, (shallow, deep) in P.DEPTH_RUNS.items()}
    pairs |= {f"lfm a{a:g} - a1": (r, P.INTERP_RUN) for a, r in INTERP.items() if r in data}
    for name, (a, b) in pairs.items():
        res["paired"][name] = {}
        for vw in ("centred", "shifted"):
            ids = [i for i in data[a][vw]["item_id"] if i in set(data[b][vw]["item_id"])]
            ia = [data[a][vw]["item_id"].index(i) for i in ids]
            ib = [data[b][vw]["item_id"].index(i) for i in ids]
            pa, pb = parts(data[a][vw]), parts(data[b][vw])
            res["paired"][name][vw] = {
                k: ci(pa[k][ia] - pb[k][ib], draws(ids)) for k in ("sd_geom", "var_geom", "var_legality")
            }
    (B2 / "posthoc.json").write_text(json.dumps(res, indent=1))
    for r, t in res["tables"].items():
        print(
            r.ljust(24),
            {vw: {k: t[vw][k][0] for k in ("var", "var_legality", "var_geom", "sd_geom", "geom_share")} for vw in t},
        )
    for name, t in res["paired"].items():
        print(name.ljust(28), t)


if __name__ == "__main__":
    main()

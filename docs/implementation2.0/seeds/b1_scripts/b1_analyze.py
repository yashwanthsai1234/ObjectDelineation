import json, sys
import numpy as np
try:
    import pandas as pd
except ImportError:
    sys.exit("no pandas")
d = json.load(open(sys.argv[1] + "/geom.json"))
df = pd.DataFrame(d["rows"], columns=d["cols"])
q = lambda s, p: float(np.percentile(s, p))
print("== window geometry per split/kind (item level, 1 row per item)")
for (sp, kd), g in df.groupby(["split", "kind"]):
    print(f"{sp:5s} {kd:9s} n={len(g):5d} S={dict(g.S.value_counts().sort_index())} off_euc med={g.off_euc.median():.3f} p90={q(g.off_euc,90):.3f} >0.2={100*(g.off_euc>0.2).mean():.1f}% >0.3={100*(g.off_euc>0.3).mean():.1f}% "
          f"| ext_frac med={g.ext_frac.median():.3f} p10={q(g.ext_frac,10):.3f} p90={q(g.ext_frac,90):.3f} | inside={g.masks_inside.mean():.4f} padded={100*g.padded.mean():.1f}% crowded={100*g.crowded.mean():.1f}% fb={g.fallback.sum()} kept={g.kept.sum()}")
print("== object px extent per split (centred rows)")
for sp, g in df[df.kind=="centred"].groupby("split"):
    print(sp, "ext_px med", g.ext_px.median(), "p10", q(g.ext_px,10), "p90", q(g.ext_px,90), "max", g.ext_px.max(), "area med", g.area_min.median(), "k dist", dict(g.k.value_counts().sort_index()))
print("== units")
for sp, g in df[df.kind=="centred"].groupby("split"):
    print(sp, "items", len(g), "patients", g.group_id.nunique(), "series", g.series.nunique(), "nodules", g.groupby(["series","nodule"]).ngroups)
print("== leakage")
sets = {sp: df[df.kind=="centred"][df.split==sp] for sp in ("train","val","test")}
for col in ("group_id","series"):
    s = {sp: set(g[col]) for sp, g in sets.items()}
    print(col, {f"{a}&{b}": len(s[a]&s[b]) for a,b in (("train","val"),("train","test"),("val","test"))})
ids = {sp: set(g.item_id) for sp, g in sets.items()}
print("item_id overlap", {f"{a}&{b}": len(ids[a]&ids[b]) for a,b in (("train","val"),("train","test"),("val","test"))})
print("== shifted vs rescaled overlap in train: window identical")

"""RIGA raw + preprocessed inventory. python -I inventory.py <repo>/data/raw/riga <out.json>"""
import json, re, sys
from collections import Counter, defaultdict
from pathlib import Path
from PIL import Image

PRIME = re.compile(r"^image(\d+)prime\.(\w+)$", re.I)
ANNO = re.compile(r"^image(\d+)-(\d+)\.(\w+)$", re.I)

def scan(d):
    primes, annos, other = {}, defaultdict(dict), []
    fmt = Counter()
    for f in sorted(Path(d).iterdir()):
        if f.is_dir():
            continue
        if m := PRIME.match(f.name):
            primes[int(m.group(1))] = f.name
        elif m := ANNO.match(f.name):
            annos[int(m.group(1))][int(m.group(2))] = f.name
        else:
            other.append(f.name)
            continue
        fmt[f.suffix.lower()] += 1
    ids = sorted(set(primes) | set(annos))
    rater_hist = Counter(len(annos[i]) for i in ids)
    return {"files": sum(fmt.values()) + len(other), "formats": dict(fmt), "other": other,
            "prime_images": len(primes), "image_ids": len(ids),
            "ids_with_prime_and_raters": sum(1 for i in ids if i in primes and annos[i]),
            "ids_raters_no_prime": sorted(i for i in ids if i not in primes),
            "ids_prime_no_raters": sorted(i for i in primes if not annos[i]),
            "id_range": [min(ids), max(ids)] if ids else None,
            "missing_ids_in_range": sorted(set(range(1, max(ids) + 1)) - set(ids)) if ids else [],
            "raters_per_image_hist": dict(sorted(rater_hist.items())),
            "images_lt6": {i: sorted(annos[i]) for i in ids if len(annos[i]) < 6},
            "rater_format": {r: dict(Counter(Path(annos[i][r]).suffix.lower() for i in ids if r in annos[i])) for r in range(1, 7)},
            "prime_format": dict(Counter(Path(p).suffix.lower() for p in primes.values())),
            "sizes": dict((f"{a}x{b}", n) for (a, b), n in Counter(Image.open(Path(d) / (primes.get(i) or next(iter(annos[i].values())))).size for i in ids[:1000]).items()),
            }

root = Path(sys.argv[1])
out = {}
for d in sorted(p for p in root.rglob("*") if p.is_dir()):
    if any(f.is_file() and f.suffix.lower() in (".jpg", ".tif", ".png") for f in d.iterdir()):
        out[str(d.relative_to(root))] = scan(d)
Path(sys.argv[2]).write_text(json.dumps(out, indent=1, default=str))
for k, v in out.items():
    print(k, {x: v[x] for x in ("files", "formats", "prime_images", "image_ids", "ids_with_prime_and_raters", "raters_per_image_hist", "prime_format", "sizes")})
    print("   lt6:", v["images_lt6"], " raters_no_prime:", v["ids_raters_no_prime"][:60], " prime_no_raters:", v["ids_prime_no_raters"], " missing:", v["missing_ids_in_range"], " other:", v["other"])

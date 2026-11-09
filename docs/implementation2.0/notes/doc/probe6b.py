import re, sys
from fnmatch import fnmatch
from pathlib import Path
root = Path(sys.argv[1])
FROZEN = ["docs/research/**", "experiments/phase0/**", "experiments/phase1/**", "experiments/**/lab/**",
          "experiments/**/impl/**", "experiments/cuda_checker/**"]
pat = re.compile(r"(?i)\b(trl|transformers|torch|vllm)(?:\[\w+\])?(?:==|\s+)v?(\d+\.\d+(?:\.\d+)?)\b(?!\s*(?:s|ms|GB|MB|%|x)\b)")
url = re.compile(r"github\.com/huggingface/(trl|transformers)/blob/v(\d+\.\d+\.\d+)/")
files = [p for p in (root/"src").rglob("*.py")] + list((root/"tests").rglob("*.py"))
files += list((root/"experiments").rglob("*.py")) + list((root/"experiments").rglob("Dockerfile*")) + list((root/"experiments").rglob("configs/*.yaml"))
files += [root/"README.md", root/"ROADMAP.md"] + [p for p in (root/"docs").rglob("*.md")]
import tomllib
pins = {n: set() for n in ("trl","transformers","torch","vllm")}
for lock in (root/"experiments/phase2/docker/locks").glob("*.txt"):
    for l in lock.read_text().splitlines():
        m = re.match(r"(trl|transformers|torch|vllm)==(\S+)", l)
        if m: pins[m.group(1)].add(m.group(2))
pp = tomllib.loads((root/"pyproject.toml").read_text())
for dep in sum(pp["project"]["optional-dependencies"].values(), []):
    m = re.match(r"(trl|transformers|torch|vllm)==(\S+)", dep)
    if m: pins[m.group(1)].add(m.group(2))
print(pins)
for f in files:
    rel = f.relative_to(root).as_posix()
    if "__pycache__" in rel or any(fnmatch(rel, g) for g in FROZEN):
        continue
    for i, line in enumerate(f.read_text(errors="ignore").splitlines(), 1):
        for m in list(pat.finditer(line)) + list(url.finditer(line)):
            n, v = m.group(1).lower(), m.group(2)
            if not any(p == v or p.startswith(v + ".") for p in pins[n]):
                print(f"{rel}:{i}: {n} {v}  | {line.strip()[:110]}")

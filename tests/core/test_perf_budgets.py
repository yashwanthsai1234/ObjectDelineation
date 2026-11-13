"""E.9.4 PF1 and PF4: the layer contract check stays fast, and the pure layers import fast without torch.
PF1's recorded times (tests/goldens/perf.json, per platform: the minimum of 3 warm runs on the Mac) bound a run at
2x; a platform without a record checks only the plan's fixed limits."""
import json
import platform
import subprocess
import sys
import time

import pytest

from tests.support.paths import ROOT

PERF = ROOT / "tests/goldens/perf.json"
PLATFORM = f"{sys.platform}-{platform.machine()}"
FAST_DEPS = ("numpy", "cv2", "scipy")  # only their modules may take more than 0.5 s of self time


def best_of_3(argv: list[str]) -> float:
    times = []
    for _ in range(3):
        t = time.perf_counter()
        subprocess.run(argv, cwd=ROOT, check=True, capture_output=True)
        times.append(time.perf_counter() - t)
    return min(times)


@pytest.mark.slow
@pytest.mark.parametrize(("key", "argv", "limit"), [
    ("lint_imports_s", [str(ROOT / ".venv/bin/lint-imports"), "--no-cache", "--no-logo"], 60),
    ("make_lint_s", ["make", "-s", "lint"], 120),
])
def test_pf1_lint_time(key, argv, limit):
    took = best_of_3(argv)
    recorded = json.loads(PERF.read_text()).get(PLATFORM, {}).get(key)
    assert took < limit
    if recorded is not None:
        assert took <= 2 * recorded, f"{key}: {took:.2f} s, recorded {recorded} s"


def test_pf4_pure_layers_import_time():
    r = subprocess.run([sys.executable, "-X", "importtime", "-I", "-c",
                        "import sys, objdel.core, objdel.geom, objdel.verify; print('torch' in sys.modules)"],
                       cwd=ROOT, check=True, capture_output=True, text=True)
    assert r.stdout.strip() == "False"
    rows = [line.split("|") for line in r.stderr.splitlines() if line.startswith("import time:") and "self" not in line]
    self_us = {name.strip(): int(us.split(":")[1]) for us, _, name in rows}
    assert sum(self_us.values()) < 2e6, f"{sum(self_us.values()) / 1e6:.2f} s"
    slow = {m: us for m, us in self_us.items() if us > 5e5 and m.split(".")[0] not in FAST_DEPS}
    assert not slow, slow

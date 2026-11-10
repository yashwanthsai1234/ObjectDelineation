"""E8 step 3: a live CLI loads no test code: `python -I -X importtime -m <module> --help` imports no `tests.*` and no
`oracles.*` module."""

import re
import subprocess
import sys

import pytest

from objdel.infra.vast.remote import LIVE_CLIS


@pytest.mark.parametrize("module", LIVE_CLIS)
def test_cli_loads_no_test_module(module):
    p = subprocess.run([sys.executable, "-I", "-X", "importtime", "-m", module, "--help"], capture_output=True,
                       text=True, timeout=300)  # fmt: skip
    assert p.returncode == 0, p.stderr[-2000:]
    loaded = re.findall(r"^import time:\s+\d+ \|\s+\d+ \|\s*(\S+)$", p.stderr, re.M)
    assert loaded and [m for m in loaded if m.split(".")[0] in ("tests", "oracles")] == []

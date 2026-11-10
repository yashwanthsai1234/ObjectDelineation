"""E8 step 3: every live CLI (objdel.infra.vast.remote.LIVE_CLIS) answers --help as `python -I -m <module>`."""

import subprocess
import sys

import pytest

from objdel.infra.vast.remote import LIVE_CLIS


@pytest.mark.parametrize("module", LIVE_CLIS)
def test_cli_help(module):
    p = subprocess.run([sys.executable, "-I", "-m", module, "--help"], capture_output=True, text=True, timeout=300)
    assert p.returncode == 0 and "usage" in p.stdout.lower(), p.stderr[-2000:]

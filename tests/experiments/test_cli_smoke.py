"""E8 step 3: every live CLI (objdel.infra.vast.remote.LIVE_CLIS) answers --help as `python -I -m <module>`."""


import pytest

from objdel.infra.vast.remote import LIVE_CLIS
from tests.support.proc import run_isolated


@pytest.mark.parametrize("module", LIVE_CLIS)
def test_cli_help(module):
    p = run_isolated("-m", module, "--help", timeout=300)
    assert p.returncode == 0 and "usage" in p.stdout.lower(), p.stderr[-2000:]

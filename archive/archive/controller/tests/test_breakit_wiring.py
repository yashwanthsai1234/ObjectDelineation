"""Controller tests moved from tests/test_breakit_wiring.py (CLEAN6, from commit f531f86). Archived: never collected, linted or
run. To run them, restore the Controller from the tag pre-clean-2026-11-08 (README.md)."""
import re
# chain: the fixture of tests/test_breakit_wiring.py with the Controller attached (alarm_inputs collected
# by `Controller(trainer, alarm=lambda m: alarm_inputs.append(dict(m)) and None)`), as at f531f86.


def test_b4_4_alarm_receives_hack_audit_statistics(chain):
    m = chain["alarm_inputs"][0]
    assert any(re.search(r"a_amax|above_tail|components|dup", k) for k in m)

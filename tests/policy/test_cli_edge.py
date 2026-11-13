"""The CLI edge helpers: Refusal -> exit (policy.errors.exits_on_refusal, Refusal.usage), InfraError -> exit
(infra.errors.exits), and the stage CLIs that use them."""
import pytest

from objdel.infra.errors import InfraError, exits
from objdel.policy.errors import Refusal, exits_on_refusal


def test_refusal_usage_is_argparse_worded_with_exit_code_2():
    e = Refusal.usage("train_rl", "--steps: bad")
    assert str(e) == "train_rl: error: --steps: bad" and e.code == 2


def test_exits_on_refusal_turns_a_refusal_into_the_exit_the_cli_always_gave(capsys):
    @exits_on_refusal
    def main(argv=None, fail=None):
        if fail:
            raise fail
        return ("ok", argv)

    assert main(["x"]) == ("ok", ["x"]) and main.__name__ == "main"
    with pytest.raises(SystemExit) as e:  # code 1: Python prints the message itself
        main(fail=Refusal("refused: a"))
    assert e.value.code == "refused: a"
    with pytest.raises(SystemExit) as e:  # another code: the message goes to stderr, the code is the exit status
        main(fail=Refusal.usage("p", "m"))
    assert e.value.code == 2 and capsys.readouterr().err == "p: error: m\n"
    with pytest.raises(ValueError):  # anything else passes through
        main(fail=ValueError("v"))


def test_exits_turns_an_infra_error_into_a_prefixed_exit():
    @exits("tool")
    def main(fail=None):
        if fail:
            raise fail
        return 7

    assert main() == 7 and main.__name__ == "main"
    with pytest.raises(SystemExit) as e:
        main(InfraError("no box"))
    assert e.value.code == "tool: no box"
    with pytest.raises(KeyError):
        main(KeyError("k"))


def test_the_stage_clis_exit_on_a_refusal_instead_of_a_traceback(monkeypatch, capsys):
    from objdel.cli.train import probe

    def refuse(cfg, prereg):
        raise Refusal("probe: refused")

    monkeypatch.setattr("objdel.rl.probe.run", refuse)
    monkeypatch.setattr("objdel.core.spec.cli_stage", lambda cls, argv=None: object())
    with pytest.raises(SystemExit) as e:
        probe.main([])
    assert e.value.code == "probe: refused"

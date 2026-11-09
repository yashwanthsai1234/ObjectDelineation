"""E1.2: one source for each reward constant."""
from objdel.eval import signal
from objdel.reward.spec import REWARD_V4, TERMS


def test_signal_weights_equal_reward_spec():
    assert tuple(signal.WEIGHTS) == TERMS
    assert tuple(signal.WEIGHTS.values()) == REWARD_V4.weights

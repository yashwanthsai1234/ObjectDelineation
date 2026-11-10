"""Archived by CLEAN8 from tests/test_p2_breakit_scoring.py (retired scripts and plans; nothing runs these)."""


def test_topology_probe_c2_input_keys_match_c2_analyze():
    seed = {"bins": [[0.5] * 6] * 3, "false_nonsimple": 0.0, "legal": 1.0}  # score_topology's keys
    res = C2A.analyze([seed] * 3, boot=50)
    assert res["pooled"] == [0.5, 0.5, 0.5] and res["outcome"] == "pass"


def _bins(hits_per_bin, n=120, K=8):
    """Per-item K-sample means whose bin totals are exactly the given hit counts."""
    out = []
    for h in hits_per_bin:
        x = np.zeros(n * K)
        x[:h] = 1
        out.append(x.reshape(n, K).mean(1).tolist())
    return out


def test_c2_pooled_exactly_on_the_floor_is_not_a_fail():
    seed = {"bins": _bins((349, 349, 44)), "false_nonsimple": 0.0, "legal": 1.0}
    # exact pooled = 0.8 * 349 / 960 + 0.2 * 44 / 960 = 288 / 960 = 0.30; float: 0.29999999999999993
    res = C2A.analyze([seed] * 3, boot=50)
    assert res["outcome"] != "fail"


def test_c2_pooled_exactly_on_the_pass_line_counts_as_050():
    seed = {"bins": _bins((597, 597, 12)), "false_nonsimple": 0.0, "legal": 1.0}
    res = C2A.analyze([seed] * 3, boot=50)
    assert res["pooled"][0] >= 0.5

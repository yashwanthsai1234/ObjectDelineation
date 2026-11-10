"""Archived by CLEAN8 from tests/test_p2_probes.py (retired scripts and plans; nothing runs these)."""


def test_topology_output_feeds_c2_analyze(topology_dir):
    c2 = pytest.importorskip("experiments.phase2.c2_analyze")
    seeds = [P.score_topology(_rows(topology_dir, M.model_noise(0.5), k=8, seed=s), topology_dir) for s in (1, 2, 3)]
    for s in seeds:
        assert [len(b) for b in s["bins"]] == [8, 8, 8] and all(0 <= v <= 1 for b in s["bins"] for v in b)
    out = c2.analyze([{k: s[k] for k in ("bins", "false_nonsimple", "legal")} for s in seeds])
    assert out["seeds"] == 3 and 0 < out["pooled"][0] < 1 and out["outcome"]
    exact = P.score_topology(_rows(topology_dir, M.model_exact), topology_dir)
    assert c2.analyze([exact] * 3)["pooled"][0] == pytest.approx(1.0)

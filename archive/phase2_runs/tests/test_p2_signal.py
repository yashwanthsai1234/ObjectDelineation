"""Archived by CLEAN8 from tests/test_p2_signal.py (retired scripts and plans; nothing runs these)."""


def test_c2_stratified_bootstrap_equals_power_c2():
    rng = np.random.default_rng(11)
    bins = [rng.random(120) for _ in range(3)]
    assert C2.pooled_ci(bins, seed=0) == pytest.approx(LAB["c2_pooled"], abs=1e-12)  # rounded to 12 decimals (BK2-4)
    assert C2.WEIGHTS == tuple(LAB["W"])


def fake_seed(rng, level, legal=0.95, fns=0.01):
    return {
        "bins": [np.clip(level + 0.1 * rng.standard_normal(120), 0, 1).tolist() for _ in range(3)],
        "false_nonsimple": fns,
        "legal": legal,
    }


def test_c2_analyze_cli_and_outcomes(tmp_path, monkeypatch, capsys):
    rng = np.random.default_rng(2)
    paths = []
    for i, level in enumerate((0.7, 0.68, 0.72)):
        p = tmp_path / f"s{i}.json"
        p.write_text(json.dumps({"probes": {"topology": fake_seed(rng, level) | {"ids": list(range(120))}}}))
        paths.append(str(p))
    out = tmp_path / "c2.json"
    monkeypatch.setattr(sys, "argv", ["c2_analyze", "--seeds", *paths, "--out", str(out)])
    C2.main()
    monkeypatch.setattr(sys, "argv", ["prereg"])
    P.main()
    assert '"pass"' in capsys.readouterr().out
    r = json.loads(out.read_text())
    assert r["outcome"] == "pass" and r["conditional"] is None and r["clauses"]["2"] and not r["seed_discordant"]
    r = C2.analyze([fake_seed(rng, v) for v in (0.7, 0.45, 0.4)])  # 1 of 3 seeds >= 0.50
    assert r["outcome"] == "partial" and r["conditional"] == "R-x8" and r["fourth_seed_eligible"]
    r = C2.analyze([fake_seed(rng, v) for v in (0.2, 0.25, 0.2)])
    assert r["outcome"] == "fail" and r["conditional"] == "R-ctrl"
    r4 = C2.analyze([fake_seed(rng, v) for v in (0.7, 0.6, 0.45, 0.65)])  # 3 of 4 >= 0.50
    assert r4["clauses"]["2"] and r4["outcome"] == "pass"

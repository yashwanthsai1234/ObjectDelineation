"""Archived by CLEAN8 from tests/test_p2_hf_eval.py (retired scripts and plans; nothing runs these)."""


def test_hf_rate_step_samples_base_weights_and_writes_the_rate(tmp_path, monkeypatch):
    pytest.importorskip("trl")
    from experiments.phase2 import hf_rate as HR
    from experiments.phase2 import r0_parity

    seen = []

    def sample_view(a, branch, rows, out, k, extra=()):
        seen.append({"branch": branch, "rows": len(rows), "k": k, "sampler": a.sampler, "extra": extra})
        Path(extra[1]).write_text(json.dumps({"items": len(rows), "samples": len(rows) * k, "generate_s": 64.0,
                                              "s_per_sample": 0.25, "load_s": 9.0, "batch_items": 16,
                                              "tokens_mean": 90.0, "length_share": 0.0}))
        return {}

    monkeypatch.setattr(r0_parity, "weights", lambda model: tmp_path / "base")
    monkeypatch.setattr(evaluate, "view_rows",
                        lambda view, items=None, rl_only=True: [{"key": str(i)} for i in range(items)])
    monkeypatch.setattr(sft_select, "sample_view", sample_view)
    out = tmp_path / "hf-rate-lfm"
    assert HR.main(["--model", "lfm2.5-vl-450m", "--view", str(tmp_path / "v"), "--out", str(out)]) == 0
    (s,) = seen
    assert s["branch"] == tmp_path / "base" and s["rows"] == 32 and s["k"] == 8 and s["sampler"] == "hf"
    assert s["extra"][0] == "--timing"
    doc = json.loads((out / "hf_rate.json").read_text())
    assert doc["model"] == "lfm2.5-vl-450m" and doc["s_per_sample"] == 0.25 and doc["load_s"] == 9.0


def test_r0_measures_the_hf_rate_of_the_rung():
    steps = _yaml(CONFIGS / "r0.yaml")["steps"]
    rate = [s for s in steps if "experiments.phase2.hf_rate" in s]
    assert [s[s.index("--model") + 1] for s in rate] == ["qwen3-vl-2b"]
    assert all(s[0] == "~" for s in rate)  # a failure leaves project_r1 without a rate (exit 2), not R0 stopped
    assert all(s[s.index("--view") + 1] == "{data}/eval_views/lidc-val-shifted" for s in rate)


def _report(tmp_path, name, level, ids, rng):
    topo = {"bins": [np.clip(level + 0.1 * rng.standard_normal(120), 0, 1).tolist() for _ in range(3)],
            "false_nonsimple": 0.01, "legal": 0.95, "ids": ids}
    p = tmp_path / f"{name}.json"
    p.write_text(json.dumps({"only": "topology", "branch": "e1", "sampler": "hf", "k": 8,
                             "probes": {"topology": topo}}))
    return str(p)


def test_c2_analyze_reads_select_ckpts_topology_report(tmp_path, monkeypatch):
    pytest.importorskip("scipy")
    from experiments.phase2 import c2_analyze as C2

    rng, ids = np.random.default_rng(2), [f"lidc/p{i}" for i in range(680)]
    paths = [_report(tmp_path, f"s{i}", v, ids, rng) for i, v in enumerate((0.7, 0.68, 0.72))]
    out = tmp_path / "c2.json"
    monkeypatch.setattr(sys, "argv", ["c2_analyze", "--seeds", *paths, "--out", str(out)])
    C2.main()
    assert json.loads(out.read_text())["outcome"] == "pass"
    bad = _report(tmp_path, "other", 0.7, ids[::-1], rng)  # other items: per-item means could not be paired
    monkeypatch.setattr(sys, "argv", ["c2_analyze", "--seeds", *paths[:2], bad, "--out", str(out)])
    with pytest.raises(SystemExit, match="items"):
        C2.main()

"""Archived by CLEAN8 from tests/test_p2_v1_select.py (retired scripts and plans; nothing runs these)."""


class _Proc:
    """A processor whose prompt ids are the prompt's characters (vLLM echoes them unless told otherwise)."""

    def __call__(self, text, images, add_special_tokens, return_tensors):
        assert add_special_tokens is False and len(images[0]) == 1
        return {"input_ids": np.array([[ord(ch) for ch in text[0]]])}


def _r0(tmp_path, monkeypatch, *, bad_ids=(), shift=0.0):
    from PIL import Image

    from experiments.phase2 import r0_parity as R

    Image.new("RGB", (8, 8)).save(tmp_path / "im.png")
    seen = {"teacher": [], "sample": []}

    def sample_view(a, branch, rows, out, k, extra=()):
        seen["sample"].append({"branch": branch, "rows": [r["key"] for r in rows], "k": k, "extra": tuple(extra),
                               "sampler": a.sampler, "model": a.model})
        out.mkdir(parents=True, exist_ok=True)
        ps = [{"id": r["key"], "prompt": f"<{r['key']}>", "image": str(tmp_path / "im.png"),
               "processor_kwargs": {}, "n": k} for r in rows]
        cs = [{"id": p["id"], "k": j, "token_ids": [9, 9, 9], "logprobs": [-1.0 + shift] * 3,
               "prompt_token_ids": [ord(ch) for ch in p["prompt"]][: -1 if p["id"] in bad_ids else None]}
              for p in ps for j in range(k)]
        (out / "prompts.jsonl").write_text("".join(json.dumps(p) + "\n" for p in ps))
        (out / "completions.jsonl").write_text("".join(json.dumps(c) + "\n" for c in cs))
        return {}

    def teacher(model, processor, x, ids):
        seen["teacher"].append(x["prompt"])
        return np.full(len(ids), -1.0)

    monkeypatch.setattr(R, "weights", lambda model: tmp_path / "base")
    monkeypatch.setattr(evaluate, "view_rows", lambda view, items=None, rl_only=True: [{"key": f"lidc/n{i:02d}/z0"}
                                                                                 for i in range(items)])
    monkeypatch.setattr(sft_select, "sample_view", sample_view)
    monkeypatch.setattr(parity, "load", lambda model_dir, prompts, dtype: ("model", _Proc()))
    monkeypatch.setattr(parity, "teacher", teacher)
    out = tmp_path / "parity-lfm"
    argv = ["--model", "lfm2.5-vl-450m", "--view", str(tmp_path / "lidc-val-shifted"), "--out", str(out)]
    return R, argv, out, seen


def test_r0_parity_passes_on_matching_ids_and_logprobs(tmp_path, monkeypatch):
    R, argv, out, seen = _r0(tmp_path, monkeypatch)
    assert R.main(argv) == 0
    (s,) = seen["sample"]
    assert s["branch"] == tmp_path / "base" and len(s["rows"]) == 32 and s["k"] == 8
    assert s["extra"] == ("--logprobs",) and s["sampler"] == "vllm" and s["model"] == "lfm2.5-vl-450m"
    assert seen["teacher"] == [f"<{k}>" for k in s["rows"][:8] for _ in range(8)]  # 8 items x n = 8
    doc = json.loads((out / "parity.json").read_text())
    assert doc["pass"] is True and doc["dtype"] == "bfloat16" and doc["sequences"] == 64
    assert doc["prompt_ids"] == {"items": 32, "equal": 32, "mismatched": [], "pass": True}
    assert doc["thresholds"] == {"median_abs": 0.02} and doc["report_only"] == {"seq_p99_abs": 0.5}  # D13, 2026-11-05
    assert doc["model"] == "lfm2.5-vl-450m" and doc["weights"] == str(tmp_path / "base")


def test_r0_parity_fails_on_one_prompt_id_mismatch(tmp_path, monkeypatch):
    R, argv, out, _ = _r0(tmp_path, monkeypatch, bad_ids={"lidc/n20/z0"})  # outside the logprob items
    assert R.main(argv + ["--items", "24"]) == 1
    doc = json.loads((out / "parity.json").read_text())
    assert doc["prompt_ids"] == {"items": 24, "equal": 23, "mismatched": ["lidc/n20/z0"], "pass": False}
    assert doc["pass"] is False and doc["median_abs"] == 0.0


def test_r0_parity_fails_on_a_logprob_shift(tmp_path, monkeypatch):
    R, argv, out, _ = _r0(tmp_path, monkeypatch, shift=0.05)
    assert R.main(argv) == 1
    doc = json.loads((out / "parity.json").read_text())
    assert doc["prompt_ids"]["pass"] is True and doc["median_abs"] == pytest.approx(0.05) and doc["pass"] is False


def test_r0_parity_help_and_no_prepared(tmp_path, monkeypatch):
    R, argv, _, _ = _r0(tmp_path, monkeypatch)
    for extra, code in ((["--help"], 0), (["--prepared", "p"], 2)):  # --prepared was never used: not accepted
        with pytest.raises(SystemExit) as e:
            R.main(argv + extra)
        assert e.value.code == code

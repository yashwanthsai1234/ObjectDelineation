"""objdel.cli.eval.gates: the Gate B1 / B3 CLI on planted final_eval folders gives the library's verdict."""
import json

from objdel.cli.eval import gates as run_gates
from objdel.verify import rungates as G
from tests.verify.test_rungates import FROZEN, converted, world


def _view(root, name, meta):
    (root / name / "lidc").mkdir(parents=True)
    rows = [{"item_id": i, "group_id": m["group"],
             "view": {"offset_frac": m["offset"], "padded": m["padded"], "fallback": m["fallback"]}}
            for i, m in meta.items()]
    (root / name / "lidc" / "val.jsonl").write_text("".join(json.dumps(r) + "\n" for r in rows))
    return root / name


def _eval(root, name, iou, view, sha=None):
    out = root / name
    out.mkdir()
    (out / "items.jsonl").write_text("".join(json.dumps({"item_id": i, "iou_tc": v, "legal": 1.0}) + "\n"
                                             for i, v in iou.items()))
    sha = FROZEN[view]["content_sha256"] if sha is None else sha
    (out / "final_eval.json").write_text(json.dumps({"view": f"data/eval_views/{view}", "view_content_sha256": sha}))
    return out


def _world(tmp_path):
    meta, base = world()
    frozen = tmp_path / "frozen_views.json"
    frozen.write_text(json.dumps(FROZEN))
    views = ["--view", str(_view(tmp_path, "lidc-val-shifted", meta)),
             "--edge-view", str(_view(tmp_path, "lidc-val-edge", meta)), "--frozen", str(frozen)]
    return meta, base, views


def test_b1_cli_writes_the_library_verdict(tmp_path, capsys):
    meta, base, views = _world(tmp_path)
    b1 = _eval(tmp_path, "b1", converted(base), "lidc-val-shifted")
    b1e, a3 = _eval(tmp_path, "b1e", base, "lidc-val-edge"), _eval(tmp_path, "a3", base, "lidc-val-shifted")
    a3e, r5 = _eval(tmp_path, "a3e", base, "lidc-val-edge"), _eval(tmp_path, "r5", base, "lidc-val-shifted", "r5")
    out = tmp_path / "gate_b1.json"
    assert run_gates.main(["b1", *views, "--b1", str(b1), "--b1-edge", str(b1e), "--r5", str(r5), "--a3", str(a3),
                           "--a3-edge", str(a3e), "--identity", "holds", "--out", str(out)]) == 0
    got = json.loads(out.read_text())
    want = G.evaluate_b1(G.load_run(b1), G.load_run(b1e), a3=G.load_run(a3), a3_edge=G.load_run(a3e),
                         r5=G.load_run(r5), meta=meta, edge_meta=meta, frozen=FROZEN, chain_ok=False,
                         identity_ok=True)
    assert got["result"] == want["result"] == "pass" and got["rules"] == want["rules"]
    assert json.loads(capsys.readouterr().out)["result"] == "pass"


def test_b1_without_a3_and_without_pre_waits_for_v2(tmp_path):
    _, base, views = _world(tmp_path)
    b1 = _eval(tmp_path, "b1", converted(base), "lidc-val-shifted")
    b1e, r5 = _eval(tmp_path, "b1e", base, "lidc-val-edge"), _eval(tmp_path, "r5", base, "lidc-val-shifted", "r5")
    out = tmp_path / "g.json"
    run_gates.main(["b1", *views, "--b1", str(b1), "--b1-edge", str(b1e), "--r5", str(r5), "--out", str(out)])
    assert json.loads(out.read_text())["result"] == "wait for V2 item 1b"


def test_b3_cli_and_a_missing_run_exits(tmp_path):
    meta, base, views = _world(tmp_path)
    b3 = _eval(tmp_path, "b3", converted(base, 3), "lidc-val-shifted")
    b3e, b1 = _eval(tmp_path, "b3e", base, "lidc-val-edge"), _eval(tmp_path, "b1", base, "lidc-val-shifted")
    b1e = _eval(tmp_path, "b1e", base, "lidc-val-edge")
    sel = tmp_path / "select.json"
    sel.write_text(json.dumps(sorted({m["group"] for m in meta.values()})))
    out = tmp_path / "g3.json"
    run_gates.main(["b3", *views, "--b3", str(b3), "--b3-edge", str(b3e), "--b1", str(b1), "--b1-edge", str(b1e),
                    "--select", str(sel), "--out", str(out)])
    assert json.loads(out.read_text())["result"] in ("pass", "fail")
    (tmp_path / "empty").mkdir()
    try:
        run_gates.main(["b3", *views, "--b3", str(tmp_path / "empty"), "--b3-edge", str(b3e), "--b1", str(b1),
                        "--b1-edge", str(b1e), "--select", str(sel), "--out", str(out)])
    except SystemExit as e:
        assert "no items.jsonl" in str(e)
    else:
        raise AssertionError("an empty run folder must exit")

"""B3.7: the lab report on synthetic scored groups (the real round needs A3's saved val answers)."""
import shutil

import numpy as np
import pytest

from objdel.verify.reward.spec import REWARD_V5
from tests.support.paths import ROOT
from tools.studies.reward_lab_v5 import (
    allmiss_grad_share,
    b3_grpo,
    lab_gate,
    report,
    rloo,
    spearman,
    variance_share,
    write_grpo,
)


def groups():
    """Two all-miss groups (location falls with distance), one mixed group, one all-hit group."""
    rows = []
    d = np.r_[1.0, 2.0, 3.0, 3.5, 4.0, 5.0, 6.0, 8.0]
    for _ in range(2):
        loc = np.exp(-d**2 / 2)
        rows.append(dict(iou_tc=np.zeros(8), loc=loc, d=d))
    rows.append(dict(iou_tc=np.r_[0.6, 0.5, 0.05, 0.0, 0.7, 0.4, 0.0, 0.3], loc=np.zeros(8), d=np.full(8, 0.2)))
    rows.append(dict(iou_tc=np.linspace(0.5, 0.9, 8), loc=np.zeros(8), d=np.full(8, 0.1)))
    n = 8 * len(rows)
    T = dict(format_ok=np.ones(n), legal=np.ones(n), iou_tc=np.concatenate([r["iou_tc"] for r in rows]),
             f5=np.full(n, 0.5), g_area=np.ones(n), topo_any=np.zeros(n), topo_freq=np.ones(n), sig_id=np.full(n, 7.0),
             loc=np.concatenate([r["loc"] for r in rows]))
    T["iou_min"] = T["iou_tc"]
    return T, np.concatenate([r["d"] for r in rows])


def test_spearman_of_constant_input_is_zero():
    assert spearman([1, 1, 1], [1, 2, 3]) == 0.0 and spearman([1, 2, 3], [3, 4, 5]) == 1.0


def test_report_sees_the_location_signal_and_its_gates():
    T, d = groups()
    rep = report(T, d, ["p1", "p2", "p3", "p4"], np.zeros(32, bool))
    assert rep["allmiss_groups"] == 2 and rep["loc_signal"]["groups"] == 2
    assert rep["loc_signal"]["v5"]["mean"] > 0.9 and rep["loc_signal"]["v5_minus_v4"]["mean"] > 0
    assert rep["drop_share"] == {"v4": 1.0, "v5": 0.0}
    assert rep["miss_beats_hit"] == 0.0 and rep["adv_size"]["r"] > 0
    assert rep["allmiss_grad_share"]["v5"] > 0 and rep["runaway"]["share"] == 0.0


def test_report_splits_the_advantage_variance_by_term():
    """Rule 3: in all-miss groups only loc varies, so it carries the whole variance; rule 8: the better half of the
    all-hit group has a valid topology, so r_topo carries part of that group's variance and the jump is c_topo x 1."""
    T, d = groups()
    T["topo_any"][28:32] = 1.0  # the four best IoUs, so topology and IoU move together
    rep = report(T, d, ["p1", "p2", "p3", "p4"], np.zeros(32, bool))
    share = rep["variance_share"]["allmiss"]
    assert np.isclose(share["loc"], 1.0) and np.isclose(share["area"], 0.0)
    assert np.isclose(sum(share.values()), 1.0)
    topo = rep["topology"]
    assert 0 < topo["hit_share"] < 1 and np.isclose(topo["jump"], 1.0)


def test_report_bounds_c_topo_by_the_fallback_margin():
    """Rule 8: an answer no doctor drew (topo_any 0) beats the best valid answer of its group by iou_tc - best valid
    iou_tc; c_topo must be at least q95 of that margin over eps, so the fallback never beats a valid topology."""
    n = 8
    T = dict(format_ok=np.ones(n), legal=np.ones(n), f5=np.full(n, 0.5), g_area=np.ones(n), loc=np.zeros(n),
             iou_tc=np.array([0.6, 0.55, 0.5, 0.52, 0.7, 0.65, 0.4, 0.62]),
             topo_any=np.r_[np.ones(4), np.zeros(4)], topo_freq=np.r_[np.full(4, 0.5), np.zeros(4)],
             sig_id=np.arange(n, dtype=float))
    T["iou_min"] = T["iou_tc"]
    topo = report(T, np.ones(n), ["p1"], np.zeros(n, bool))["topology"]
    q95 = np.percentile([0.1, 0.05, -0.2, 0.02], 95)
    assert np.isclose(topo["fallback_margin_q95"], q95)
    assert np.isclose(topo["c_topo_min"], q95 / REWARD_V5.eps)


def test_report_has_no_fallback_bound_without_fallback_answers():
    T, d = groups()
    topo = report(T, d, ["p1", "p2", "p3", "p4"], np.zeros(32, bool))["topology"]
    assert np.isnan(topo["c_topo_min"])


def test_u65_gradient_share_hand_value_and_flag():
    """U65 (rule 10): an all-miss group with A = +-0.05 and a mixed group with A = +-0.5: S = 0.4 / 4.4; a group
    with A = 0 adds 0 to the denominator; S_v5 < 0.5 S_v4 fails the pass flag."""
    pm = np.r_[1, -1, 1, -1, 1, -1, 1, -1]
    a = np.stack([0.05 * pm, 0.5 * pm, np.zeros(8)])
    allmiss = np.array([True, False, True])
    s = allmiss_grad_share(a, a, allmiss)
    assert s["v4"] == pytest.approx(0.4 / 4.4) and s["v4"] == pytest.approx(0.0909, abs=1e-4)
    assert s["ratio"] == pytest.approx(1.0) and s["pass"]
    low = np.stack([0.02 * pm, 0.5 * pm, np.zeros(8)])  # S_v5 = 0.16 / 4.16 < 0.5 x 0.0909
    assert not allmiss_grad_share(a, low, allmiss)["pass"]
    assert allmiss_grad_share(np.zeros((1, 8)), np.zeros((1, 8)), np.array([True]))["v4"] == 0.0


def test_u66_topology_share_is_one_or_zero():
    """U66 (rule 8): c_topo r_topo the only varying term in hit groups carries the whole variance; constant, none."""
    rows = np.ones(1, bool)
    topo = np.linspace(0, 1, 8)[None]
    parts = {"iou": np.full((1, 8), 0.7), "topo": topo}
    assert np.isclose(variance_share(parts, rloo(topo), rows)["topo"], 1.0)
    iou = np.linspace(0.5, 0.9, 8)[None]
    parts = {"iou": iou, "topo": np.ones((1, 8))}
    assert np.isclose(variance_share(parts, rloo(iou), rows)["topo"], 0.0)


def lab_report(**change):
    """A report that passes every pre-run Gate B3 lab rule."""
    rep = {"loc_signal": {"v5": {"mean": 0.31, "ci": [0.01, 0.6]}, "v5_minus_v4": {"mean": 0.01, "ci": [-0.1, 0.2]}},
           "drop_share": {"v4": 0.5, "v5": 0.5}, "miss_beats_hit": 0.01,
           "allmiss_grad_share": {"v4": 0.2, "v5": 0.1, "ratio": 0.5, "pass": True},
           "topology": {"hit_share": 0.5, "jump": 1.0}}
    for k, v in change.items():
        top, _, sub = k.partition("__")
        rep[top] = (rep[top] | {sub: v}) if sub else v
    return rep


def test_lab_gate_rules_and_boundaries():
    """Gate B3 before the run (B.6.2): loc mean > 0.3, CI lower bound > 0, (v5 - v4) mean > 0; drop share v5 <= v4;
    miss beats hit <= 1 %; S_v5 >= 0.5 S_v4 (U65); topology share <= 0.5 passes, 0.5001 fails (U66)."""
    assert all(lab_gate(lab_report()).values())
    for k, bad in [("loc_signal__v5", {"mean": 0.3, "ci": [0.01, 0.6]}),
                   ("loc_signal__v5", {"mean": 0.31, "ci": [0.0, 0.6]}),
                   ("loc_signal__v5_minus_v4", {"mean": 0.0, "ci": [-0.1, 0.2]}),
                   ("drop_share__v5", 0.5001), ("miss_beats_hit", 0.0101),
                   ("allmiss_grad_share__pass", False), ("topology__hit_share", 0.5001)]:
        assert not all(lab_gate(lab_report(**{k: bad})).values()), k


def test_u52_area_gate_and_loc_width_rules():
    """U52: area_gate turns True exactly when the miss-beats-hit share > 1 %; loc_width keeps the kernel width when
    the median missed d / r <= 3, else the width that gives the median missed answer loc 0.1."""
    T, d = groups()
    rep = report(T, d, ["p1", "p2", "p3", "p4"], np.zeros(32, bool))
    assert rep["miss_beats_hit"] == 0.0 and rep["area_gate"] is False
    assert rep["distance"]["median"] > 3  # the all-miss groups' d reach 8
    w = rep["distance"]["loc_width"]
    assert np.isclose(np.exp(-rep["distance"]["median"] ** 2 / (2 * w**2)), 0.1)
    near = report(T, d / 10, ["p1", "p2", "p3", "p4"], np.zeros(32, bool))["distance"]
    assert near["median"] <= 3 and near["loc_width"] == REWARD_V5.loc_width
    T["iou_tc"][16:24] = np.r_[0.1, 0.5, 0.05, 0.0, 0.7, 0.4, 0.0, 0.3]  # the mixed group: a hit at 0.1
    T["f5"][16], T["f5"][18] = 0.0, 1.0  # a miss at 0.05 with F5 1 beats it: 1 of 1 mixed group
    T["iou_min"] = T["iou_tc"]
    beats = report(T, d, ["p1", "p2", "p3", "p4"], np.zeros(32, bool))
    assert beats["miss_beats_hit"] == 1.0 and beats["area_gate"] is True


def test_u71_adam_eps_scales_with_the_gradient():
    """U71: AdamW with (g / r, eps / r) takes the steps of (g, eps), r = 10."""
    torch = pytest.importorskip("torch")
    r = 10.0
    out = []
    for scale, eps in ((1.0, 1e-8), (1 / r, 1e-8 / r)):
        torch.manual_seed(0)
        model = torch.nn.Linear(4, 3).double()
        opt = torch.optim.AdamW(model.parameters(), lr=1e-2, eps=eps)
        x = torch.randn(16, 4, dtype=torch.float64)
        for _ in range(5):
            opt.zero_grad()
            (model(x) ** 2).mean().mul(1e-6 * scale).backward()  # small gradients, so eps matters
            opt.step()
        out.append(torch.cat([p.detach().flatten() for p in model.parameters()]))
    assert torch.allclose(out[0], out[1], rtol=1e-6, atol=0)


def test_u52_u71_b3_config_gets_beta_and_adam_eps_from_one_r(tmp_path):
    """U52 and U71: the builder writes beta = 1e-3 / r (0.0 with the owner flag) and adam_epsilon = 1e-8 / r into
    the B3 spec's train stage; adam_epsilon reaches the GRPOConfig train_rl builds."""
    from objdel.core.spec import RlStage, load_spec, to_argv

    assert b3_grpo(4.0) == {"beta": 1e-3 / 4, "adam_epsilon": 1e-8 / 4}
    assert b3_grpo(4.0, beta_zero=True) == {"beta": 0.0, "adam_epsilon": 1e-8 / 4}
    spec = tmp_path / "b3_v5.yaml"
    shutil.copy(ROOT / "specs/b3_v5.yaml", spec)
    write_grpo(spec, b3_grpo(4.0))
    train = next(s for s in load_spec(spec).stages if isinstance(s, RlStage) and s.name == "train")
    assert train.grpo == {"mask_truncated_completions": False, "beta": 2.5e-4, "adam_epsilon": 2.5e-9}
    write_grpo(spec, b3_grpo(4.0, beta_zero=True))
    train = next(s for s in load_spec(spec).stages if isinstance(s, RlStage) and s.name == "train")
    assert train.grpo["beta"] == 0.0 and train.grpo["adam_epsilon"] == 2.5e-9
    pytest.importorskip("trl")
    from tests.experiments.test_move_invariance import grpo_config_dict

    argv = to_argv(train, "rd", "data")
    cfg = grpo_config_dict(argv[argv.index("objdel.cli.train.rl") + 1:], save_strategy="steps",
                           save_steps=50, shuffle_dataset=False)
    assert cfg["adam_epsilon"] == 2.5e-9 and cfg["beta"] == 0.0

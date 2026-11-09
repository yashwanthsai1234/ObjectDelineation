"""B3.7: the lab report on synthetic scored groups (the real round needs A3's saved val answers)."""
import numpy as np

from experiments.phase2.reward_lab_v5 import report, spearman


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

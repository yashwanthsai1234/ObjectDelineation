"""experiments.phase2.gate_power (B1.9): power of the Gate B1 rules (1)-(4) under a simulated far-miss effect, with
K = 8 sampling noise on both arms; Gate B3 power as a report line (E10, E10b, E10c)."""

import numpy as np
import pytest

from experiments.phase2 import gate_power as GP
from objdel.eval import signal

RUNS = 300


def synthetic(patients=60, seed=0):
    """Per patient: 4 near items (offset 0.05), 4 mid (0.25), 6 far (0.4: 3 misses at 0.02, 3 found); size 15 px.
    The first half of the patients is val-select."""
    rng = np.random.default_rng(seed)
    ids, mu, groups, off = [], [], [], []
    for p in range(patients):
        for j, (o, v) in enumerate([(0.05, rng.uniform(0.5, 0.9)) for _ in range(4)]
                                   + [(0.25, rng.uniform(0.4, 0.9)) for _ in range(4)]
                                   + [(0.4, 0.02)] * 3 + [(0.4, rng.uniform(0.4, 0.9)) for _ in range(3)]):
            ids.append(f"lidc/p{p}/n000/z{j}")
            mu.append(v)
            groups.append(p)
            off.append(o)
    off, groups = np.array(off), np.array(groups)
    return GP.Items(ids=ids, mu=np.array(mu), groups=groups, cell=GP.cell_of(off, np.full(len(off), 15)),
                    far=off > GP.FAR, near=off < GP.NEAR, select=groups < patients // 2)


@pytest.fixture(scope="module")
def data():
    items = synthetic()
    return items, GP.weights(items.ids, items.groups)


def test_cells_follow_the_b1_edges_right_closed():
    """B1.1 cells: searchsorted side="right" on SIZE_EDGES (10, 20, 40) px and OFFSET_EDGES (0.15, 0.2, 0.3)."""
    assert GP.cell_of(np.array([0.0, 0.15, 0.31, 0.2]), np.array([9, 10, 40, 25])).tolist() == [0, 5, 15, 10]


def test_weights_reproduce_the_gate_patient_bootstrap(data):
    """The vectorized CI equals signal.ci on signal.draws' patient draws (same B and seed as the gate)."""
    items, W = data
    x = np.random.default_rng(1).normal(size=len(items.ids))
    D = signal.draws(items.ids, key=dict(zip(items.ids, items.groups.tolist(), strict=True)).__getitem__)
    reps = W @ x / W.sum(1)
    lo, hi = np.quantile(reps, [0.025, 0.975])
    assert W.shape == (signal.BOOT, len(x)) and np.allclose([x.mean(), lo, hi], signal.ci(x, D), rtol=0, atol=1e-12)


def test_noise_sd_of_an_item_at_one_half_is_sqrt_c_quarter_over_8():
    z = np.random.default_rng(2).standard_normal(5)
    for c in (0.077, 1.0):
        assert np.allclose((GP.noisy(np.full(5, 0.5), c, z) - 0.5) / z, np.sqrt(c * 0.25 / 8))
    assert np.array_equal(GP.noisy(np.full(5, 0.5), 0.0, z), np.full(5, 0.5))


def test_the_arms_draw_independent_noise():
    """E10: the correlation of the two arms' noise vectors over 10,000 items is within 0.03 of 0."""
    _, z = GP.streams(0, 1, 10_000)
    assert abs(np.corrcoef(z["base"][0], z["b1"][0])[0, 1]) < 0.03


def test_baseline_without_noise_is_the_one_arm_run(data):
    """E10b: with c = 0 on the baseline arm and c = 1 on the B1 arm, the baseline is the R5 K-means and the B1 arm
    is the one it is when the baseline is noisy too (same seed): the one-arm run of revision 5."""
    items, _ = data
    t0, b1_0, base_0 = GP.arms(items, GP.EFFECT, 0.0, 1.0, 50, 7)
    t1, b1_1, base_1 = GP.arms(items, GP.EFFECT, 1.0, 1.0, 50, 7)
    assert np.array_equal(base_0, np.broadcast_to(items.mu, base_0.shape))
    assert np.array_equal(t0, t1) and np.array_equal(b1_0, b1_1) and not np.array_equal(base_1, base_0)


def test_convert_takes_the_share_of_far_misses_patient_by_patient_from_found_items_of_the_cell(data):
    items, _ = data
    far_miss = items.far & (items.mu < GP.MISS)
    true, _, _ = GP.arms(items, GP.EFFECT, 0.0, 0.0, 40, 3)
    found_far = set(items.mu[items.far & (items.mu >= GP.MISS)])
    for row in true:
        moved = row != items.mu
        assert moved.sum() == round(GP.EFFECT * far_miss.sum()) and not (moved & ~far_miss).any()
        assert set(row[moved]) <= found_far
        pats = np.unique(items.groups[moved])
        assert len(pats) <= int(np.ceil(moved.sum() / 3)) + 1  # whole patients (3 far misses each), one partial


def test_huge_effect_passes_every_rule(data):
    items, W = data
    s = GP.summary(GP.simulate(items, W, 1.0, GP.C_MEASURED, GP.C_MEASURED, RUNS, 0)["rules"])
    assert s["joint"] == 1.0 and set(s["p_pass"]) == {"iou", "far_miss", "near_iou", "did"}


def test_no_effect_keeps_the_iou_rule_at_its_level(data):
    items, W = data
    noisy = GP.summary(GP.simulate(items, W, 0.0, 1.0, 1.0, 1000, 0)["rules"])["p_pass"]["iou"]
    assert noisy <= 0.05 + 0.03
    assert GP.summary(GP.simulate(items, W, 0.0, 0.0, 0.0, 20, 0)["rules"])["p_pass"]["iou"] == 0.0


def test_the_difference_in_differences_rule_is_applied():
    """A far gain smaller than the near gain fails rule (4) while the IoU rule passes."""
    items = synthetic()
    W = GP.weights(items.ids, items.groups)
    base = np.broadcast_to(items.mu, (1, len(items.mu)))
    b1 = base + np.where(items.far, 0.05, np.where(items.near, 0.2, 0.0))
    flags = GP.b1_rules(b1, base, items.far, items.near, W)
    assert flags["iou"].tolist() == [True] and flags["did"].tolist() == [False]
    flags = GP.b1_rules(base + np.where(items.far, 0.2, 0.0), base, items.far, items.near, W)
    assert flags["did"].tolist() == [True] and flags["near_iou"].tolist() == [True]


def test_far_miss_point_rule_reads_20_5_percent():
    """Rule (2): the B1 far-miss share <= 20.5 % and the paired CI upper end < 0."""
    items = synthetic()
    W = GP.weights(items.ids, items.groups)
    base = np.broadcast_to(items.mu, (1, len(items.mu))).copy()
    b1 = base.copy()
    fm = np.flatnonzero(items.far & (items.mu < GP.MISS))  # 180 of 360 far items
    b1[0, fm[:106]] = 0.8  # 74 / 360 = 20.6 % left
    assert GP.b1_rules(b1, base, items.far, items.near, W)["far_miss"].tolist() == [False]
    b1[0, fm[106]] = 0.8  # 73 / 360 = 20.3 %
    assert GP.b1_rules(b1, base, items.far, items.near, W)["far_miss"].tolist() == [True]


def test_report_joint_is_the_four_gate_b1_rules_only(data):
    """E10c: joint P(pass) is recomputed from the per-run flags; the Gate B3 report never moves it; the other Gate
    B1 rules are listed as not simulated; Gate B3 is printed per rule at each B3 share; the smallest effect line
    exists for every noise line."""
    items, W = data
    lines = {"c=0.077": (0.077, 0.077), "noise 0": (0.0, 0.0)}
    doc = GP.report(items, W, 100, 0, lines, b3_shares=(0.25, 0.5))
    for line in doc["lines"].values():
        flags = np.array([line["flags"][r] for r in GP.B1_RULES], bool)
        assert line["joint"] == flags.all(0).mean() <= min(line["p_pass"].values())
        assert set(line["b3"]) == {"0.25", "0.5"}
        assert set(line["b3"]["0.5"]) == {"p_pass", "joint"} and set(line["b3"]["0.5"]["p_pass"]) == set(GP.B3_RULES)
        assert "smallest_effect" in line
    hopeless = GP.report(items, W, 100, 0, lines, b3_shares=(0.0,))
    assert hopeless["lines"]["c=0.077"]["joint"] == doc["lines"]["c=0.077"]["joint"]
    assert hopeless["lines"]["c=0.077"]["b3"]["0.0"]["joint"] == 0.0
    assert doc["not_simulated"] == ["edge view", "legal share", "frozen-view hash"]


def test_c_from_a3_item_sd_halves_the_retest_variance_per_arm():
    """math-stats §11 item 1: A3's item_sd is the SD of a retest difference, so one arm has item_sd^2 / 2."""
    mu = np.array([0.5, 0.5])
    assert GP.c_from_item_sd(0.1, mu) == pytest.approx(0.005 / (0.25 / 8))

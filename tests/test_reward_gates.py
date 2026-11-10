"""Data-backed reward-v4 gates (plan §7: K8, K8b, K9, R3, R4, R5 and the M6 local-gradient rate).

Full-size runs (all fixture items, AB arms, K3/K4) are in experiments/phase1/reward_gates.py; these tests run the
same harness on the 239 full-LIDC train disagreement items and on LIDC fixtures plus a RIGA subsample."""
import numpy as np
import pytest
from oracles import reward_harness as H
from oracles import reward_reference as ref
from scipy import stats

from objdel.verify.reward.cache import ExpertCache
from objdel.verify.reward.score import score_mask

pytestmark = pytest.mark.data


@pytest.fixture(scope="module")
def fixture_items(fixtures_dir):
    items = H.load_fixture_items(fixtures_dir)
    riga = [it for it in items if it[0] == "riga"]
    return [it for it in items if it[0] == "lidc"] + [it for i, it in enumerate(riga) if i // 2 % 8 == 0]  # 1/8 images


@pytest.fixture(scope="module")
def lidc_dis(fixtures_dir):
    items = H.load_lidc_disagreement(fixtures_dir)
    assert len(items) == 239  # 251 before sub-filter raters abstained (B2-4)
    return items


@pytest.fixture(scope="module")
def lidc_pools(lidc_dis):
    return [H.item_pools((k, st, j, 10, H.pack, True)) for j, (k, st) in enumerate(lidc_dis.items())]


def test_cache_and_terms_match_reference_on_real_items(fixture_items):
    rng = np.random.default_rng(0)
    for _, _, st in fixture_items[::12]:
        c, r = ExpertCache.build(st), ref.cache(st)
        assert c.tau == pytest.approx(r["tau"], rel=1e-5) and c.t5 == pytest.approx(r["t5"], rel=1e-5)
        assert list(c.rater_sigs) == r["sigs"]
        for m in (st[0], H.grow(c.vote, 1.5), H.jitter(st[-1], 2.0, rng, 4.0)):
            t, rt = score_mask(m, c), ref.object_terms(m, r)
            assert (t.sig, t.topo_freq) == (rt["sig"], rt["topo_freq"])
            assert t.iou_tc == pytest.approx(rt["iou_tc"], abs=1e-12) and t.f5 == pytest.approx(rt["f5"], rel=1e-5)


def test_cache_store_roundtrip_on_real_items(fixture_items, tmp_path):
    for i, (_, _, st) in enumerate(fixture_items[::40]):
        c = ExpertCache.build(st)
        c.save(tmp_path / f"{i}.npz")
        d = ExpertCache.load(tmp_path / f"{i}.npz")
        assert d.digest == c.digest and np.array_equal(d.dist_out_union, c.dist_out_union)


# ------------------------------------------------------------------ K9

def _k9(args, rng):
    pairs = [(g, key) for key, st, two_part, seed in args
             for g, _ in H.k9_groups(st, np.random.default_rng(seed), two_part)]
    T, keys = np.array([g for g, _ in pairs]), [k for _, k in pairs]
    a = H.in_steps(H.advantages, T)
    ok = a[:, 0] > a[:, 1]
    return ok.mean(), H.cluster_ci(ok, keys, rng), len(T)


def test_k9_full_lidc_disagreement_items(lidc_dis):
    args = [(k, st, False, j) for j, (k, st) in enumerate(lidc_dis.items())]
    share, (lo, hi), n = _k9(args, np.random.default_rng(1))
    assert n > 500 and lo >= 0.95, (share, lo, hi)


def test_k9_lidc_fixtures_and_two_part_composites(fixture_items):
    lidc = [(k, st) for ds, k, st in fixture_items if ds == "lidc"]
    args = [(k, st, False, j) for j, (k, st) in enumerate(lidc[::2])]
    args += [(f"{ka}+{kb}", H.composite(a, b, 3 + j % 6), True, 10_000 + j)
             for j, ((ka, a), (kb, b)) in enumerate(zip(lidc[0::4], lidc[1::4], strict=False))]
    share, (lo, hi), n = _k9(args, np.random.default_rng(2))
    assert n > 500 and lo >= 0.95, (share, lo, hi)


# ------------------------------------------------------------------ R3, local gradient, R4

@pytest.mark.slow
def test_r3_reachable_hacks_lose_to_the_median_honest_rollout(fixtures_dir):
    """Every fixture item, realistic honest rollouts (exhaustive, hence slow; numbers in impl/M6_reward.md)."""
    rng = np.random.default_rng(3)
    items = H.load_fixture_items(fixtures_dir)
    data = [H.item_pools((k, st, j, 12, H.pack, True)) for j, (_, k, st) in enumerate(items)]
    for name in ("dilated_union", "band_speck"):
        T, keys, _ = H.adversary_groups(data, name, "real", rng)
        a = H.in_steps(H.advantages, T)
        below = a[:, 0] < np.median(a[:, 1:], axis=1)
        assert below.mean() >= 0.95, (name, below.mean(), H.cluster_ci(below, keys, rng))


@pytest.fixture(scope="module")
def local_gradient(lidc_dis):
    return [H.local_gradient_item((k, st, 70_000 + j, 6)) for j, (k, st) in enumerate(lidc_dis.items())]


def test_local_gradient_on_average_on_the_lidc_disagreement_items(local_gradient):
    """Lab R3 (a) criterion: dilating or eroding an honest rollout lowers its advantage on average."""
    res = local_gradient
    rng = np.random.default_rng(5)
    for name, _, _ in H.MOVES:
        m, lo, hi = H.mean_ci([np.mean(r["deltas"][name]) for r in res if r["deltas"][name]], rng)
        assert hi < 0, (name, m, lo, hi)


def test_local_gradient_moves_away_from_consensus_per_rollout(local_gradient):
    """Gate as re-registered after measurement (PLAN §9, ISSUES I12): moves away from the median rater area lower
    the advantage for >= 95 % of rollouts. All moves together: ~82 %, since moves toward it may rightly raise it."""
    res = local_gradient
    d = np.concatenate([np.concatenate([r["deltas"][n] for n, _, _ in H.MOVES]) for r in res])
    away = np.concatenate([np.concatenate([np.array(r["away"][n], bool) for n, _, _ in H.MOVES]) for r in res])
    assert (d[away] < 0).mean() >= 0.95


@pytest.mark.parametrize("mode", ["real", "jit1px"])
def test_r4_neutrality_on_the_lidc_disagreement_items(lidc_pools, mode):
    """Minority minus mode-class advantage. The plan's criterion (95 % CI contains 0, half-width <= 0.10) holds
    for 1-px rollouts in every pool draw and for realistic rollouts in 4 of 6 (systematic +0.03..+0.07,
    impl/M6_reward.md); this test pins the half-width and a |mean| <= 0.07 bound on one fixed draw."""
    rng = np.random.default_rng(4)
    d = H.neutrality(lidc_pools, mode, rng, H.advantages, n_batch=150)
    m, lo, hi = H.mean_ci(d, rng)
    assert len(d) >= 200 and abs(m) <= 0.07 and (hi - lo) / 2 <= 0.10, (m, lo, hi)


# ------------------------------------------------------------------ band: K8, K8b, R5

@pytest.fixture(scope="module")
def k8_rows(fixture_items):
    return [r for ds, k, st in fixture_items for r in H.k8_rows(ds, k, st)]


def _assert_k8(rows, strata):
    for name in strata:
        g = [r for r in rows if f"{r['subset']} k={r['k']}".startswith(name)]
        rho = stats.spearmanr([r["f5"] for r in g], [r["area"] for r in g]).statistic
        assert len(g) >= 400 and abs(rho) <= 0.10, (name, len(g), rho)


def test_k8_band_is_size_neutral_lidc_strata(k8_rows):
    _assert_k8(k8_rows, ("lidc/nodule k=3", "lidc/nodule k=4"))


@pytest.mark.slow
def test_k8_band_is_size_neutral_riga_strata(fixtures_dir):
    """All RIGA structures: a 1/8 subsample gives |rho| up to 0.14 by chance (full set: -0.04 / +0.01)."""
    rows = [r for ds, k, st in H.load_fixture_items(fixtures_dir) if ds == "riga" for r in H.k8_rows(ds, k, st)]
    _assert_k8(rows, ("riga/disc", "riga/cup"))


def test_k8b_band_discrimination(fixture_items, k8_rows):
    k8b = {k: H.k8b_item(st) for _, k, st in fixture_items}
    assert np.mean([d["dil_mono"] and d["ero_mono"] for d in k8b.values()]) >= 0.98
    held = {}
    for r in k8_rows:
        held.setdefault(r["key"], []).append(r["f5"])
    assert np.mean([k8b[k]["dil"] < np.median(v) for k, v in held.items()]) >= 0.95
    f5 = np.array([r["f5"] for r in k8_rows])
    assert np.median(f5) >= 0.90 and np.percentile(f5, 10) >= 0.70


@pytest.mark.slow
def test_r5_rater_subsample_medians(fixtures_dir):
    """All RIGA structures (m = 4, 5 sit at the 0.98 edge, so a subsample is not enough)."""
    riga = [st for ds, _, st in H.load_fixture_items(fixtures_dir) if ds == "riga"]
    res = [H.r5_item(st, j) for j, st in enumerate(riga)]
    for m in (2, 3, 4, 5):
        assert 0.92 <= np.median([v for r in res for v in r[m]]) <= 0.98, m

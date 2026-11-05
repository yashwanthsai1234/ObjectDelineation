"""Scoring on the CUDA checker: objdel.geometry.to_masks (objdel._gpucheck) and every scorer above it give the same
results as the polycore reference (to_mask). On a CUDA box the chain runs this file on the GPU build before training
(policy.chain); elsewhere it runs on the same source compiled for the CPU."""
import json

import numpy as np
import pytest

from experiments.cuda_checker import compare as C
from experiments.phase2.vast import policy as P
from objdel import _gpucheck
from objdel.eval import metrics
from objdel.geometry import parse, serialize, to_mask, to_masks
from objdel.reward import ExpertCache, score_batch, score_many
from objdel.reward import score as score_mod

ROOT = P.ROOT


def disk(n, r, cy, cx):
    y, x = np.mgrid[:n, :n]
    return np.hypot(y - cy, x - cx) <= r


def parts_of(rings, roles):
    """A fuzz case's rings as to_mask parts, or None when a hole's role does not name an outer ring."""
    roles = [-1] * len(rings) if roles is None else list(roles)
    if any(p != -1 and not (0 <= p < len(rings) and roles[p] == -1) for p in roles):
        return None
    return [(np.asarray(rings[i]).reshape(-1, 2), [np.asarray(rings[j]).reshape(-1, 2)
                                                    for j, p in enumerate(roles) if p == i])
            for i, p in enumerate(roles) if p == -1]  # fmt: skip


def fuzz_objects():
    out = []
    for _, rings, roles, W, H, G in C.fuzz_polygons(20, seed=3):
        parts = parts_of(rings, roles)
        if G == 1000 and parts is not None and sum(1 + len(h) for _, h in parts) <= 10:
            out.append((parts, W, H))
    return out


def same_geometry(a, b) -> bool:
    return a.legal == b.legal and a.defects == b.defects and a.mask.shape == b.mask.shape and \
        np.array_equal(a.mask, b.mask)


def reference(objs):
    return [to_mask(parts, W, H) for parts, W, H in objs]


def test_to_masks_equals_polycore():
    objs = fuzz_objects() + [([], 5, 4)]  # + an object with no parts: a legal empty mask
    assert len(objs) > 300
    cpu, gpu = reference(objs), to_masks(objs)
    assert [i for i, (a, b) in enumerate(zip(cpu, gpu, strict=True)) if not same_geometry(a, b)] == []
    assert any(not g.legal for g in cpu) and any(g.legal and g.mask.any() for g in cpu)


def test_to_masks_of_nothing():
    assert to_masks([]) == [] and _gpucheck.build() in ("cuda", "cpu")


@pytest.fixture(scope="module")
def jobs():
    """LIDC (one object) and RIGA-like (two objects) jobs: fuzz outlines on the 0-1000 grid, an empty answer, an
    illegal bow-tie, a hole, garbage, a truncated and a wrong-schema completion."""
    nod = ExpertCache.build(np.stack([disk(40, r, 20, 21) for r in (8, 9, 10)]))
    disc = ExpertCache.build(np.stack([disk(60, r, 30, 30) for r in (20, 21, 22)]))
    cup = ExpertCache.build(np.stack([disk(60, r, 31, 30) for r in (8, 9, 9)]))
    lidc, riga = {"nodule": nod}, {"disc": disc, "cup": cup}
    texts = []
    for parts, _, _ in fuzz_objects():
        if all(np.asarray(r).min() >= 0 and np.asarray(r).max() <= 1000 and len(r) >= 3
               for o, hs in parts for r in (o, *hs)):
            texts.append(serialize({"nodule": parts}, ("nodule",)))
    texts = [t for t in texts if parse(t, ("nodule",)).ok][:150]
    texts += ['{"nodule":[]}', '{"nodule":[{"outer":[[100,100],[900,900],[900,100],[100,900]],"holes":[]}]}',
              '{"nodule":[{"outer":[[200,200],[800,200],[800,800],[200,800]],"holes":[[[400,400],[400,600],'
              '[600,600],[600,400]]]}]}', "garbage", '{"disc":[]}']
    out = [(t, False, ("nodule",), lidc, 40, 40) for t in texts]
    out.append((texts[0], True, ("nodule",), lidc, 40, 40))
    def square(a, b):
        return [(np.array([[a, a], [b, a], [b, b], [a, b]]), [])]

    out += [(serialize({"disc": square(150, 850), "cup": square(a, b)}, ("disc", "cup")), False, ("disc", "cup"), riga,
             60, 60) for a, b in ((350, 650), (400, 600), (100, 900))]
    return out


def test_score_batch_equals_polycore(jobs, monkeypatch):
    gpu = score_batch(jobs)
    monkeypatch.setattr(score_mod, "to_masks", reference)
    cpu = score_batch(jobs)
    assert cpu.keys() == gpu.keys()
    for k, v in cpu.items():
        assert (np.array_equal(v, gpu[k]) if isinstance(v, np.ndarray) else v == gpu[k]), k
    assert 0 < cpu["legal"].sum() < len(jobs) and cpu["format_ok"].sum() < len(jobs)


def test_score_many_keeps_job_order_and_masks(jobs):
    one = [score_many([j])[0] for j in jobs]
    many = score_many(jobs)
    for (t1, m1, _), (t2, m2, _) in zip(one, many, strict=True):
        assert json.dumps(t1, sort_keys=True) == json.dumps(t2, sort_keys=True)
        assert (m1 is None) == (m2 is None)
        if m1 is not None:
            assert all(np.array_equal(m1[n], m2[n]) for n in m1)


def test_exam_metrics_equal_polycore(jobs, monkeypatch):
    lidc = [j for j in jobs if j[2] == ("nodule",)][:16]
    texts, flags, caches = [j[0] for j in lidc], [j[1] for j in lidc], lidc[0][3]
    b = metrics(texts, {}, 40, 40, ("nodule",), flags, caches)
    monkeypatch.setattr(score_mod, "to_masks", reference)
    a = metrics(texts, {}, 40, 40, ("nodule",), flags, caches)
    assert json.dumps(a, sort_keys=True) == json.dumps(b, sort_keys=True)


def test_every_scoring_chain_builds_and_checks_the_cuda_checker_first():
    gpu, train, select, probe = P.chain(ROOT / "experiments/phase2/configs/r1_qwen3vl2b.yaml", "r1")
    assert gpu[:2] == ["bash", "-c"] and "run_box.sh --build-only && " in gpu[2]
    assert gpu[2].endswith("tests/test_gpu_scoring.py")  # this file, on the GPU build, before training
    assert [train[1], select[1], probe[1]] == ["experiments.phase2.train_sft", "experiments.phase2.select_ckpt",
                                               "experiments.phase2.g_probe"]
    (timing,) = P.chain(ROOT / "experiments/phase2/configs/r0_qwen3vl_arm.yaml", "r1")  # train only: nothing to score
    assert timing[1] == "experiments.phase2.train_sft"

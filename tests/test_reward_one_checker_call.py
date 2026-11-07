"""C3 (experiments/phase2/perf/implementation.md; issues RL-5, CK-5): one CUDA checker call per reward call, in the
main process; pool workers only score masks on the CPU and never touch objdel._gpucheck."""
import json
import pickle
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor

import numpy as np

from objdel import _gpucheck
from objdel.reward import ExpertCache
from objdel.reward import score as S
from objdel.train import rewards as R
from objdel.train import serve

G, W = 8, 50
EOS, PAD = 7, 0


def _box(x0, y0, x1, y1):
    return [[x0, y0], [x1, y0], [x1, y1], [x0, y1]]


def _compact(objs) -> str:
    return json.dumps(objs, separators=(",", ":"))


def _nodule(*ring, hole=None):
    return _compact({"nodule": [{"outer": _box(*ring), "holes": [] if hole is None else [_box(*hole)]}]})


def _riga(disc, cup):
    return _compact({"disc": [{"outer": _box(*disc), "holes": []}], "cup": [{"outer": _box(*cup), "holes": []}]})


def _stack(boxes):
    m = np.zeros((len(boxes), W, W), bool)
    for i, (a, b, c, d) in enumerate(boxes):
        m[i, b:d, a:c] = True
    return m


CACHES = {
    "a": {"nodule": ExpertCache.build(_stack([(10, 10, 30, 30), (11, 11, 31, 31), (10, 10, 32, 32)]))},
    "b": {"nodule": ExpertCache.build(_stack([(15, 12, 35, 30), (16, 12, 36, 31), (15, 13, 34, 33)]))},
    "r": {"disc": ExpertCache.build(_stack([(5, 5, 45, 45), (6, 5, 45, 44), (5, 6, 44, 45)])),
          "cup": ExpertCache.build(_stack([(18, 18, 32, 32), (19, 18, 31, 32), (18, 19, 32, 31)]))},
}
NAMES = {"a": ("nodule",), "b": ("nodule",), "r": ("disc", "cup")}
BOWTIE = '{"nodule":[{"outer":[[200,200],[600,600],[600,200],[200,600]],"holes":[]}]}'
TEXTS = {
    "a": [_nodule(200, 200, 600, 600), _nodule(220, 220, 620, 620), _nodule(200, 200, 640, 640), BOWTIE, "nonsense",
          _nodule(100, 100, 900, 900), _nodule(200, 200, 610, 600), '{"nodule":[]}'],
    "b": [_nodule(300, 240, 700, 600), _nodule(320, 240, 720, 620), _nodule(300, 260, 680, 660),
          _nodule(300, 240, 700, 600, hole=(400, 400, 500, 500)), _nodule(0, 0, 1000, 1000), '{"disc":[]}',
          _nodule(310, 250, 690, 610), _nodule(300, 240, 700, 600)],
    "r": [_riga((100, 100, 900, 900), (360, 360, 640, 640)), _riga((120, 100, 900, 880), (380, 360, 620, 640)),
          _riga((100, 120, 880, 900), (360, 380, 640, 620)), _riga((100, 100, 900, 900), (0, 0, 100, 100)),
          _riga((100, 100, 900, 900), (300, 300, 700, 700)), "garbage",
          '{"disc":[{"outer":[[100,100],[900,900],[900,100],[100,900]],"holes":[]}],"cup":[]}',
          _riga((90, 90, 910, 910), (350, 350, 650, 650))],
}


def batch(keys=("a", "b", "r"), truncate=(1, 17)):
    """TRL's reward kwargs for one group of G answers per key, in order; answer i is truncated when i in truncate."""
    texts = [t for k in keys for t in TEXTS[k]]
    ks = [k for k in keys for _ in range(G)]
    ids = [[5, 6, 6] if i in truncate else [5, 6, EOS] for i in range(len(texts))]
    return dict(prompts=[[{"role": "user", "content": "p"}]] * len(texts),
                completions=[[{"role": "assistant", "content": t}] for t in texts], completion_ids=ids, key=ks,
                W=[W] * len(texts), H=[W] * len(texts), names=[list(NAMES[k]) for k in ks])


def call_all(pool=None, kw=None):
    rw = R.TRLReward(CACHES.__getitem__, stop_ids={EOS, PAD}, pool=pool)
    funcs, _ = rw.funcs()
    kw = batch() if kw is None else kw
    return {f.__name__: f(**kw) for f in funcs}


# Values recorded from the code before C3 (commit e9f23d7: one checker call per pool chunk, in the worker).
EXPECTED = {
    "format_ok": [1.0, 0.0, 1.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 1.0, 1.0, 1.0, 0.0, 1.0, 1.0,
                  1.0, 0.0, 1.0, 1.0],
    "legal": [1.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1.0, 1.0, 0.0, 1.0, 1.0, 1.0,
              0.0, 0.0, 1.0],
    "iou_tc": [1.0, 0.0, 1.0, 0.0, 0.0, 0.3025, 1.0, 0.0, 1.0, 1.0, 1.0, 0.0, 0.152, 0.0, 0.95, 1.0, 1.0, 0.0, 1.0, 0.5,
               0.745, 0.0, 0.0, 0.9114627536519267],
    "f5": [1.0, 0.0, 1.0, 0.0, 0.0, 0.032258537578351974, 1.0, 0.16951174443254863, 1.0, 1.0, 1.0, 0.0,
           0.016914010942196488, 0.0, 1.0, 1.0, 1.0, 0.0, 1.0, 0.5505019983106747, 0.6092762746618532, 0.0, 0.0,
           0.869142931320022],
    "g_area": [1.0, 0.0, 1.0, 0.0, 0.0, 0.3025, 1.0, 0.0, 1.0, 1.0, 1.0, 0.0, 0.152, 0.0, 0.95, 1.0, 1.0, 0.0, 1.0,
               0.5744047619047619, 0.745, 0.0, 0.0, 0.9114627536519267],
    "topo_any": [1.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 1.0, 1.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1.0, 1.0, 0.0, 1.0, 1.0,
                 1.0, 0.0, 0.0, 1.0],
    "r_topo": [1.4236489301936017, 0.0, 1.4236489301936017, 0.0, 0.0, 0.0, 1.4236489301936017, 0.0, 1.1682361183106065,
               1.1682361183106065, 1.1682361183106065, 0.0, 0.0, 0.0, 1.1682361183106065, 1.1682361183106065,
               1.2798078939677113, 0.0, 1.2798078939677113, 0.0, 0.0, 0.0, 0.0, 1.2798078939677113],
    "advantage_v4": [2.9999291474211436, -1.986578267415168, 2.9999291474211436, -1.986578267415168, -1.986578267415168,
                     -1.3688139292807404, 2.9999291474211436, -1.6712387107371847, 1.6417370151746964,
                     1.6417370151746964, 1.6417370151746964, -2.8613472098603996, -2.420365003858335,
                     -2.8613472098603996, 1.576111362880349, 1.6417370151746964, 2.69745210240513, -2.18197016196783,
                     2.69745210240513, -0.8306511090607538, -0.4983515753096066, -2.18197016196783, -2.18197016196783,
                     2.4800089654635906],
}


def test_rewards_equal_the_previous_path():
    assert call_all() == EXPECTED
    with ThreadPoolExecutor(3) as pool:
        assert call_all(pool) == EXPECTED
    pool = serve.fork_pool(2)
    try:
        assert call_all(pool) == EXPECTED
    finally:
        pool.shutdown()


def test_reward_calls_the_checker_once_per_reward_call(monkeypatch):
    """256 answers (32 groups), serially and on a thread pool (where a worker's call would be counted too)."""
    calls = []
    real = _gpucheck.check_polygons
    monkeypatch.setattr(_gpucheck, "check_polygons", lambda *a: calls.append(len(a[4])) or real(*a))
    kw = batch(keys=("a", "b", "r", "a") * 8)
    assert len(kw["key"]) == 256
    serial = call_all(kw=kw)
    assert len(calls) == 1 and calls[0] > 256  # one call for every object of every parsed answer
    with ThreadPoolExecutor(4) as pool:
        assert call_all(pool, kw) == serial
    assert len(calls) == 2


CHILD = """
import json, pickle, sys
from objdel.reward import score
rows = score._score_checked(pickle.load(sys.stdin.buffer))
print(json.dumps({"rows": rows, "gpucheck": "objdel._gpucheck" in sys.modules}))
"""


def test_pool_jobs_never_touch_cuda():
    """A pool job (the main process's checked masks) is scored in a fresh interpreter that never imports
    objdel._gpucheck, with the same terms as in this process."""
    kw = batch()
    jobs = [(c[0]["content"], R.truncated(i, {EOS, PAD}), tuple(n), CACHES[k], w, h)
            for c, i, k, w, h, n in zip(kw["completions"], kw["completion_ids"], kw["key"], kw["W"], kw["H"],
                                        kw["names"], strict=True)]
    work = S._checked(jobs)
    out = subprocess.run([sys.executable, "-I", "-c", CHILD], input=pickle.dumps(work), capture_output=True,
                         check=True)
    child = json.loads(out.stdout)
    assert child["gpucheck"] is False
    assert child["rows"] == json.loads(json.dumps(S._score_checked(work)))
    assert [r["format_ok"] for r in child["rows"]] == EXPECTED["format_ok"]

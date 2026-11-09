"""B4 probe: forkserver (and spawn) workers re-import the parent's __main__ as __mp_main__, so a trainer script with
top-level torch / datasets imports re-runs them in every reward worker. HEAVY=1 imports both at top level.

    HEAVY=0|1 .venv/bin/python -I experiments/phase2/lab/B4/fs_main_probe.py
(run as a file, like a training script; stderr lists every process that imported the main module)
"""
import os
import sys
import time

HEAVY = os.environ.get("HEAVY") == "1"
if HEAVY:
    import datasets  # noqa: F401
    import torch  # noqa: F401
print("main imported in pid", os.getpid(), "as", __name__, file=sys.stderr, flush=True)

if __name__ == "__main__":
    from experiments.phase2.lab.B4 import serve
    from experiments.phase2.lab.B4.common import ROOT, first_keys
    t = time.time()
    try:
        fs = serve.forkserver_pool(2, ROOT, first_keys(32, 32))
        print("ok", round(time.time() - t, 1), flush=True)
        fs.shutdown()
    except Exception as e:
        print("FAIL", type(e).__name__, flush=True)

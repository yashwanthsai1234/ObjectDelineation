"""Imported by the fork server (serve.forkserver_pool): fills serve.STORE from B4_ROOT / B4_KEYS before any worker
is forked, so the workers share the caches copy-on-write with the single-threaded server, not with the trainer.
Errors propagate (multiprocessing's forkserver swallows only ImportError, so a missing key or a bad digest kills
the server loudly). Single-threaded on purpose: with a 2-thread preload the macOS fork server's children crashed
in the Objective-C runtime (+[NSCharacterSet initialize] in progress in another thread at fork; B4 trl_pattern
run 1)."""
import gc
import json
import os

from experiments.phase2.lab.B4 import serve

if os.environ.get("B4_ROOT") and os.environ.get("B4_KEYS"):
    _root = os.environ["B4_ROOT"]
    _keys = set(json.loads(open(os.environ["B4_KEYS"]).read()))
    serve.preload(_root, {k: v for k, v in serve.rl_digests(_root).items() if k in _keys}, threads=1)
    gc.freeze()

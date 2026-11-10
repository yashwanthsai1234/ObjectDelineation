"""B4 pre-registered pass criteria (this file was criteria.py until the coordinator asked for prereg.py; content unchanged) (fixed 2026-10-19 before any timed run; calibration before this: riga_stats only).

Every B4 script evaluates its results against these thresholds (`check`), so the verdict cannot drift after the
numbers are seen. Workloads: experiments.phase1.reward_gates.make_jobs, 32 prompts x G = 8 = 256 completions,
RIGA-only / mixed (40 % RIGA) / LIDC-only on data/prepared train caches; wall = advantages_v4(score_batch(...)),
warm (second of two runs of a batch); hash = reward_gates.digest.

H-4 preload (preload.py): every RL cache of riga/train + lidc/train into a dict, 2 threads; then 10 TRL-like steps
(seeded permutation of the RIGA train rows, 32 prompts per step) through TRLReward with the production
cache_loader(maxsize=64) vs the preloaded dict.
H-2 pool (serve.py, pool_bench.py, trl_pattern.py): CacheRef(key) jobs, fork pool created after the preload and
before torch work; forkserver + preloading module as the post-CUDA-safe alternative; 10 batches per workload,
serial and 2 / 4 / 8 workers interleaved per batch.
H-K3SYN (synth.py, k3syn.py, k3_harness.py): 64 real RIGA K3 prompts vs 64 synthetic ones, 20 repetitions, order
alternating per repetition.
H-3 decision: C++ batch scorer needed if P2-1 or P2-2 fail on the Mac, or if the x86 RIGA-only wall projected at 4
workers (Mac speed-up applied to the x86 serial estimate) exceeds 125 ms; otherwise not for Phase 2, with a design
for Phase 3 if larger groups would exceed the gate.
"""
CRITERIA = {
    # H-4
    "P4-1 preload wall s (Mac, 2 threads)": ("<=", 60.0),
    "P4-2 RSS after preload GB": ("<=", 6.0),
    "P4-3a preloaded loads during 10 steps": ("==", 0),
    "P4-3b LRU misses per step (mean)": (">=", 28.0),
    "P4-4 LRU minus preloaded step reward call, median s": (">=", 1.5),
    "P4-5 preloaded TRLReward / plain score_batch, median": ("<=", 1.10),
    "P4-6 step hashes identical LRU vs preloaded": ("==", True),
    # H-2
    "P2-1 RIGA-only serial / 4-worker median (Mac)": (">=", 3.0),
    "P2-2 RIGA-only best-pool median wall ms (Mac)": ("<=", 125.0),
    "P2-3 hashes identical to serial (all)": ("==", True),
    "P2-4 USS growth per worker MB (max)": ("<", 200.0),
    "P2-5 no-op pool overhead / 4-worker wall": ("<=", 0.10),
    "P2-6 TRL pattern: one dispatch per batch, outputs ==, 10 steps": ("==", True),
    "P2-7a Docker RIGA-only speed-up 2 workers": (">=", 1.6),
    "P2-7b Docker RIGA-only speed-up 4 workers": (">=", 2.2),
    # H-K3SYN
    "PS-1a RIGA-only synthetic/real median ratio (Mac)": ("in", (0.90, 1.10)),
    "PS-1b its 95 % bootstrap CI": ("within", (0.85, 1.15)),
    "PS-2 mixed synthetic/real median ratio (Mac)": ("in", (0.90, 1.10)),
    "PS-3 |p99 synthetic - p99 real| per rollout ms": ("<=", 1.0),
    "PS-4 RIGA-only ratio in Docker": ("in", (0.90, 1.10)),
    "PS-5 synthetic corpus digest Mac == x86": ("==", True),
}


def check(name: str, value) -> bool:
    op, ref = CRITERIA[name]
    if op == "<=":
        return value <= ref
    if op == "<":
        return value < ref
    if op == ">=":
        return value >= ref
    if op == "==":
        return value == ref
    if op == "in":
        return ref[0] <= value <= ref[1]
    if op == "within":
        return ref[0] <= value[0] and value[1] <= ref[1]
    raise ValueError(op)


def verdicts(values: dict) -> dict:
    """{criterion: {"value": v, "pass": bool}} for the criteria present in values."""
    return {k: {"value": v, "pass": bool(check(k, v))} for k, v in values.items()}

# ADR-0009: Rent from North American Vast AI hosts only, and prove everything on the Mac first

**Status:** Accepted
**Date:** 2026-11-01 (region rule, commit `787c096`); 2026-11-02 (preflight gate, commit `734b1ec`)
**Deciders:** Repository owner

## Context

GPUs are rented per hour on Vast AI. A box costs money from the moment it starts.
Early runs lost money and time on avoidable problems:

- A Thailand host could not pull the image (`policy.py` comment, V1).
- A busy, power-capped card was found only after a 12-minute bundle upload (commit `999969b`, 2026-11-03, "GPU check before the bundle upload"). Machine 28360 is now blocked (`BAD_MACHINES`).

The owner also wants a hard limit on total spend and no forgotten instance.

## Decision

Two rules, both enforced in code and not by habit.

**1. Region and offer filter (`policy.terms`, `policy.reasons`).**
Only hosts in the US or Canada (`REGION`). One GPU. Verified host. Reliability above 0.98. At least 12 effective cores. CPU RAM at least 32 GB (48 GB for RIGA and for the RTX PRO 6000 cards).
CUDA at least 13.0. `inet_down` of 500 or more. Two direct ports. A price ceiling per card (RTX 5090: $1.50 per hour).
The filter is both sent as a query and re-checked on the returned JSON.

**2. Preflight before every rental (`experiments/phase2/preflight.py`).**
`make preflight CONFIG=...` runs three stages, and a failed stage stops the later ones:

| Stage | Checks |
|---|---|
| DATA | bundle hashes, no patient in both training and an exam view, coverage of every exam view by the training view, same prompt text, consistent CT window, every training target parses and is legal |
| CHECKER | the checker and scorer tests, the exact oracles, the archived CPU checker's recorded results |
| CODE | a quick test subset, a CPU dry run of train to select on a tiny model, the projected GPU hours and dollars |

It writes a stamp, `runs/preflight/<config>.json`, with hashes of code, data and config.
`launch.sh` calls `preflight.refusal` and refuses unless the stamp passes and all three hashes still match.

Around these rules the tooling adds guards:

- Cost: a launch is refused above the hard stop (`HARD_STOP = 17.5` USD) and needs `--owner-approved` at or above the alarm (`ALARM = 10.0`) (`cmd_launch`).
- Time: the on-box watchdog destroys the instance at the cap, or when the GPU idles after the pull, or when no job ever started.
- Cleanup: `teardown.sh` destroys the box even if the final pull fails. `sweep_orphans.sh` destroys unknown or overdue instances.
- Governance: RIGA only with recorded permission and a datacenter offer (`governance`). RIGA never enters `run_box.sh`.

## Options considered

### Option A: North America only, preflight gate (chosen)

| Dimension | Assessment |
|---|---|
| Complexity | Medium. A filter, a stamp and a watchdog |
| Cost | Fewer offers, so some hours cost more |
| Safety | Failures happen on the Mac, where hours are free |

### Option B: Any region, cheapest offer

Lower hourly price. The V1 run showed a host that could not pull the image. Data governance also becomes harder.

### Option C: Test on the box

Simplest. The box bills while you debug. Not chosen: the preflight docstring states the goal as "everything that can be proven on the Mac, before any GPU is rented".

## Trade-off analysis

The filter narrows the market but removes whole classes of failure.
Preflight costs minutes on the Mac and saves hours on the box.
Recorded costs: the full 770-step run billed about $4.56. Phase 2 committed spend was $15.59 of $17.50 after R4 (`docs/research/notes/r4_speed_2026-11-08.md`).

## Consequences

- Easier: a launch cannot happen on stale code, stale data or a changed config.
- Easier: one command, `make preflight`, tells you if a rental is safe.
- Harder: every code change after a preflight makes the stamp stale, so preflight runs last.
- Harder: the hard stop needs a code change to raise. It rose from $15 to $16.60, $17.00 and $17.50 across commits `cf2f507`, `74a0135` and `edf3129`. The module docstring still says $15.
- Revisit: the `inet_up_cost` and `inet_down_cost` ceilings carry an unverified unit note in `policy.py`.

## Action items

1. [x] `terms`, `reasons`, `governance` and `refusal`.
2. [x] Watchdog, teardown and orphan sweep.
3. [ ] Update the docstring of `policy.py` to the current hard stop.

"""The offer filter as data: `query` sends it to `vastai search offers`, `reasons` re-checks it on the returned JSON;
`offer_cols` are the offer fields a launch records."""

import operator
from pathlib import Path

from objdel.infra.vast.ledger import INET_CEIL

# §4.2 / §4.8: gpu_ram bounds in GB (stock + 2, or >= 79 for the 80 GB A100), price ceilings = 1.5 x snapshot median.
CARDS = {
    "RTX_3090": {"ram": (None, 26), "ceiling": 0.50},
    "RTX_4090": {"ram": (None, 26), "ceiling": 1.00},
    "RTX_5090": {"ram": (None, 34), "ceiling": 1.50},
    "RTX_PRO_6000_WS": {"ram": (None, 98), "ceiling": 2.40, "cpu_ram": 48},
    "RTX_PRO_6000_S": {"ram": (None, 98), "ceiling": 2.40, "cpu_ram": 48},
    "A100_SXM4": {"ram": (79, 82), "ceiling": 1.80},
}
CPU_RAM = 32  # GB, the least host RAM of any card (a card's own "cpu_ram" raises it)
OPS = {"=": operator.eq, ">": operator.gt, ">=": operator.ge, "<=": operator.le}
REGION = ("US", "CA")  # North American hosts only: a Thailand host could not pull the image
# Machines that failed the GPU check on a paid box; never offered again
BAD_MACHINES = {28360: "GPU busy before our job, power limit 150 W, 7.4 of 165 bf16 TFLOPS"}


def is_raw(p: Path) -> bool:
    """True when `p` resolves under a data/raw/ folder (never bundled or synced)."""
    return "/data/raw/" in f"{Path(p).resolve().as_posix()}/"


def terms(gpu: str, disk: int = 80) -> list[tuple]:
    """(reason, query term, offer JSON key, op, limit, scale) per term. Vast reports gpu_ram (MiB) and cpu_ram (MB) and
    the query takes GB; cpu_ram>=32 returns 32160, so cpu_ram counts 1 GB = 1000 MB. `verification` is the
    string "verified" in the JSON (the query term is `verified=true`). `duration` is query-only: its unit is
    UNVERIFIED."""
    if gpu not in CARDS:
        raise ValueError(f"unknown card {gpu}; one of {sorted(CARDS)}")
    c, mb = CARDS[gpu], 1 / 1024
    (lo, hi), cpu = c["ram"], max(CPU_RAM, c.get("cpu_ram", 0))
    rows = [
        ("num_gpus", "num_gpus=1", "num_gpus", "=", 1, 1),
        ("gpu_name", f"gpu_name={gpu}", "gpu_name", "=", gpu, 1),
        ("gpu_ram", f"gpu_ram>={lo}", "gpu_ram", ">=", lo, mb) if lo else None,
        ("gpu_ram", f"gpu_ram<={hi}", "gpu_ram", "<=", hi, mb),
        ("verified", "verified=true", "verification", "=", "verified", 1),
        ("rentable", "rentable=true", "rentable", "=", True, 1),
        ("reliability", "reliability>0.98", "reliability2", ">", 0.98, 1),
        ("cores", "cpu_cores_effective>=12", "cpu_cores_effective", ">=", 12, 1),
        ("cpu_ram", f"cpu_ram>={cpu}", "cpu_ram", ">=", cpu, 1 / 1000),
        ("cuda", "cuda_vers>=13.0", "cuda_max_good", ">=", 13.0, 1),
        ("disk", f"disk_space>={disk}", "disk_space", ">=", disk, 1),
        ("inet_down", "inet_down>=500", "inet_down", ">=", 500, 1),
        ("ports", "direct_port_count>=2", "direct_port_count", ">=", 2, 1),
        ("duration", "duration>=1", None, None, None, None),
        ("price", f"dph<={c['ceiling']}", "dph_total", "<=", c["ceiling"], 1),
        ("bandwidth", f"inet_up_cost<={INET_CEIL}", "inet_up_cost", "<=", INET_CEIL, 1),
        ("bandwidth", f"inet_down_cost<={INET_CEIL}", "inet_down_cost", "<=", INET_CEIL, 1),
        ("region", f"geolocation in [{','.join(REGION)}]", None, None, None, None),  # checked in reasons()
    ]
    return [r for r in rows if r]


def query(gpu: str, disk: int = 80) -> str:
    """The `vastai search offers` query string of the filter."""
    return " ".join(t[1] for t in terms(gpu, disk))


def reasons(o: dict, gpu: str, disk: int = 80) -> list[str]:
    """Why an offer JSON fails the filter ([] = admissible); a missing or mistyped field fails."""
    bad = []
    for why, _, key, op, limit, scale in terms(gpu, disk):
        v = o.get(key)
        if key and (v is None or isinstance(v, str) != isinstance(limit, str) or not OPS[op](
                v.replace(" ", "_") if isinstance(v, str) else v * scale, limit)):  # fmt: skip
            bad.append(why)
    if (o.get("geolocation") or "").rsplit(", ", 1)[-1] not in REGION:  # "Texas, US"; none or unknown fails
        bad.append("region")
    if o.get("machine_id") in BAD_MACHINES:
        bad.append(f"machine {o['machine_id']} ({BAD_MACHINES[o['machine_id']]})")
    return [f"{r} fails the §4.2 filter" for r in dict.fromkeys(bad)]


def offer_cols(o: dict) -> dict:
    """The offer fields the ledger records with a launch."""
    return {k: o.get(k) for k in ("dph_total", "inet_up_cost", "inet_down_cost", "hosting_type", "host_id",
                                  "machine_id")}

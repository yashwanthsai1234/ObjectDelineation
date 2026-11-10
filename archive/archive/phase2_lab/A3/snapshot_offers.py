"""Vast.ai price/availability snapshot for the Phase 2 GPU short-list (lab evidence for A3_vast_gpus.md).

    python3 -I experiments/phase2/lab/A3/snapshot_offers.py --out experiments/phase2/lab/A3/vast_snapshot.json

Uses the public offer search that `vastai search offers` calls (POST https://console.vast.ai/api/v0/bundles/, the
CLI's default query `verified=true external=false rentable=true`; vastai 1.8.3 runs it without an API key). No
account, no authenticated call, ~90 requests, 0.5 s apart. Filters = what our runs need (ROADMAP "When renting",
RESULTS section 8): 1 GPU, >= 8 effective vCPUs (the reward scores 256 rollouts per step on the CPU), reliability >
0.98, disk >= 50 GB, CUDA >= 12.8 (Blackwell needs it; torch 2.14.1's PyPI wheel needs driver >= 580 = CUDA 13.0),
inet_down >= 200 Mb/s. Prices are dph_total with 100 GB of storage included (--storage 100).

References: https://docs.vast.ai/cli/reference/search-offers , https://docs.vast.ai/guides/instances/pricing
"""
import argparse
import json
import statistics as st
import time
import urllib.request

GPUS = ["RTX 3090", "RTX 4090", "RTX 5090", "L40S", "RTX 6000Ada", "RTX A6000", "A100 SXM4", "A100 PCIE",
        "H100 SXM", "H100 NVL", "H200", "H200 NVL", "B200", "RTX PRO 6000 WS", "RTX PRO 6000 S"]
BASE = {"verified": {"eq": True}, "external": {"eq": False}, "rentable": {"eq": True}, "num_gpus": {"eq": 1},
        "cpu_cores_effective": {"gte": 8}, "reliability": {"gt": 0.98}, "disk_space": {"gte": 50},
        "cuda_max_good": {"gte": 12.8}, "inet_down": {"gte": 200}, "allocated_storage": 100.0, "limit": 500,
        "order": [["dph_total", "asc"]]}


def search(q):
    req = urllib.request.Request("https://console.vast.ai/api/v0/bundles/", data=json.dumps(q).encode(),
                                 headers={"Content-Type": "application/json", "User-Agent": "objdel-a3-snapshot"})
    with urllib.request.urlopen(req, timeout=60) as r:
        return json.load(r).get("offers", [])


def summarize(offers):
    if not offers:
        return {"n": 0}
    p = sorted(o["dph_total"] for o in offers)
    q = lambda v, f: v[min(len(v) - 1, int(f * len(v)))]  # noqa: E731
    return {
        "n": len(p), "min": round(p[0], 3), "p25": round(q(p, .25), 3), "median": round(st.median(p), 3),
        "max": round(p[-1], 3),
        "vram_gb": sorted({round(o["gpu_ram"] / 1024) for o in offers}),
        "cores_eff_median": st.median(o["cpu_cores_effective"] for o in offers),
        "ram_gb_per_gpu_median": round(st.median(o["cpu_ram"] * o.get("gpu_frac", 1) for o in offers) / 1024, 1),
        "reliability_median": round(st.median(o["reliability2"] for o in offers), 4),
        "dlperf_median": round(st.median(o["dlperf"] for o in offers), 1),
        "gpu_mem_bw_median": round(st.median(o["gpu_mem_bw"] for o in offers)),
        "cuda_max_good_min": min(o["cuda_max_good"] for o in offers),
        "storage_usd_gb_month_median": round(st.median(o["storage_cost"] for o in offers), 3),
        "inet_up_usd_tb_median": round(st.median(o["internet_up_cost_per_tb"] for o in offers), 2),
        "countries": sorted({(o.get("geolocation") or "?").split(",")[-1].strip() for o in offers}),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    out = {"taken_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), "filters": BASE, "offers": {}}
    for g in GPUS:
        for typ in ("on-demand", "bid"):
            for dc in (False, True):
                q = dict(BASE, gpu_name={"eq": g}, type=typ)
                if dc:
                    q["datacenter"] = {"eq": True}
                out["offers"][f"{g}|{typ}|{'secure_cloud' if dc else 'any'}"] = summarize(search(q))
                time.sleep(0.5)
    with open(a.out, "w") as f:
        json.dump(out, f, indent=1, sort_keys=True)
    for k, v in out["offers"].items():
        if v["n"]:
            print(f"{k:40s} n={v['n']:3d} ${v['min']:.2f} / {v['median']:.2f} / {v['max']:.2f}  cores {v['cores_eff_median']}"
                  f"  dlperf {v['dlperf_median']}")


if __name__ == "__main__":
    main()

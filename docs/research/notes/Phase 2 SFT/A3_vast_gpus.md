# A3: Vast.ai operations, data governance and GPUs for Phase 2

Written 2026-10-19 by research agent A3 (Phase 2 planning, `experiments/phase2/RESEARCH.md`).

**Question.** How do we use Vast.ai for Phase 2 SFT screening (SmolVLM-256M → LFM2.5-VL-450M → Qwen3-VL-2B, ~6
GPU-hours) and for the GRPO phases? May RIGA go onto a rented host? Which five GPU types do we use, and at what cost?
What image, scripts and guardrails do we deploy, and what goes into the first-session checklist?

**Method.**
- Read the Vast docs from their index ([llms.txt](https://docs.vast.ai/llms.txt), 55 pages fetched as markdown on
  2026-10-19), the Terms of Service and the DPA, and the `vastai` CLI 1.8.3 help for every command we would use.
- Took a price and availability snapshot through the public offer search. The `vastai search offers` CLI runs it
  without an API key; no account, no authenticated call, no spend.
- Read the dataset licences (RIGA DOI metadata, MESSIDOR terms, TCIA / LIDC) and the NVIDIA architecture
  whitepapers.
- Built the candidate image for linux/amd64 on the Mac and checked it (§6.3).
- Wrote a bottom-up cost estimator.
- Lab evidence is in `experiments/phase2/lab/A3/`:
  - `snapshot_offers.py` and `vast_snapshot.json` (the price snapshot);
  - `estimate_cost.py`;
  - `Dockerfile.vast` and `image_check.sh`.

**Limits.**
- The session's WebSearch budget ran out early in this task. Community experience (Reddit, Discord, blog post-mortems)
  was therefore **not** surveyed.
- Pitfalls below come from Vast's own docs, the CLI help and the official `vastai` agent guide
  ([SKILL.md](https://raw.githubusercontent.com/vast-ai/vast-cli/master/vastai/SKILL.md)).
- Measured cross-GPU fine-tuning benchmarks for small VLMs were not found. Lambda's benchmark page renders no data
  without JavaScript ([lambda.ai/gpu-benchmarks](https://lambda.ai/gpu-benchmarks)). Throughput below rests on:
  - datasheet peaks;
  - Vast's per-machine measured DLPerf and memory bandwidth;
  - one published small-VLM datapoint (nanoVLM).
- Nothing here is legal advice. §5 is a risk analysis for the coordinator.

---

## 0. Recommendations (ranked)

1. **Phase 2 runs on on-demand, verified, single-GPU Vast instances.** Each run is launched from one pinned image.
   - The image holds dependencies only. A code tarball and the data are rsynced in per session.
   - The Mac **pulls** every artifact back. No cloud credentials ever sit on a host.
   - Every instance is **destroyed**, never stopped:
     - storage keeps billing while an instance is stopped;
     - a stopped instance's GPU can be taken by someone else, so restart may hang in "scheduling"
       ([pricing](https://docs.vast.ai/guides/instances/pricing), [troubleshooting](https://docs.vast.ai/guides/reference/troubleshooting)).
2. **The five GPUs:**

   | # | GPU | Role |
   |---|---|---|
   | 1 | RTX 4090 24 GB | default screening card for SmolVLM-256M and LFM2.5-VL-450M; named by the RESULTS §8 cost gate |
   | 2 | RTX 5090 32 GB | best DLPerf per $ among consumer cards, and the largest Secure Cloud supply of any type |
   | 3 | RTX PRO 6000 Blackwell 96 GB | Qwen3-VL-2B full SFT, and Phase 3 GRPO with vLLM colocate on one card |
   | 4 | A100 SXM4 80 GB | mature-stack datacenter fallback for Qwen; Secure Cloud supply exists |
   | 5 | RTX 3090 24 GB | cheapest box for the first session's x86 gates, smoke tests and Phase 8 CUDA work |

   H100 and H200 are kept only for Phases 3–4, and only if wall-clock time binds (§3).
3. **Money is not the constraint in Phase 2.** One screening run costs **$0.3–5** on any of the five cards. The whole
   Phase 2 package costs **$2–17** (§4).
   - What matters is iteration speed, stack risk (Blackwell needs CUDA ≥ 12.8; our torch wheel needs driver ≥ 580)
     and data governance.
   - Guardrails therefore target *forgotten instances* and *runaway runs*, not the hourly price (§7).
4. **RIGA governance:**
   - **LIDC may go anywhere.** Its licence is CC BY 3.0.
   - **RIGA never goes to community hosts.**
   - Its MESSIDOR-derived images carry MESSIDOR's terms: "Copy, redistribution … prohibited". Use is "restricted to
     individuals or organizations that obtained the database directly".
   - On Vast, hosts "can technically access files on their machines".
   - Therefore:
     - Phase 2 screening uses **LIDC only**;
     - RIGA SFT runs only on **Secure Cloud** (`datacenter=True`, ISO 27001) hosts, with application-level encryption
       at rest, no Vast Cloud Sync, and destroy on completion;
     - **in parallel, ask the RIGA depositor and the MESSIDOR partners for written permission** (§5).
5. **Stack finding to resolve before Phase 3.** Each item is backed by the PyPI metadata linked in §6.1.
   - The pinned `torch==2.14.1` Linux wheel is the **CUDA 13.0** build, which needs host driver ≥ 580.
   - **TRL 1.16's `[vllm]` extra caps vLLM at ≤ 0.31.0, and vLLM 0.31/0.32 pin `torch==2.13.0`.**
   - So the vLLM colocate path (LFM, Qwen; RESULTS §8 I22) cannot live in the training venv as pinned.
   - Phase 2 SFT does not need vLLM. A separate vLLM venv or image covers offline K-sampling and Phase 3 (§6.1, H-VLLM).
6. **The first session** costs about 1.5 h and about $1 (§8):
   - image smoke test;
   - K3/K4 x86 timing on LIDC plus synthetic RIGA-like items (the existing harness reads `data/fixtures/riga_experts.npz`,
     which is RIGA data);
   - `xarch_hashes` compare;
   - SFT and GRPO step-time measurement.

---

## 1. Vast facts we depend on

| # | Fact | Consequence for us | Source |
|---|---|---|---|
| V1 | Marketplace: hosts set prices. Billing is per second while running, plus storage per GB for every second the instance exists (stopped included), plus bandwidth per byte. No charge while "Loading". | Destroy, don't stop. Pull artifacts first. | [pricing](https://docs.vast.ai/guides/instances/pricing), [billing](https://docs.vast.ai/guides/reference/billing) |
| V2 | Three rental types:<br>• **on-demand**: fixed price, never interrupted;<br>• **reserved**: on-demand prepaid for a discount, by default 20 % for 1 month and 30 % for 3 months, with clawback plus a 5 % fee on early cancel;<br>• **interruptible**: bid. It is stopped when outbid, then resumes. Its data is kept while stopped. | Phase 2: on-demand. Reserved is pointless for hour-long runs. Interruptible only for idempotent jobs (eval sampling, seed replicates with resume). | [instance types](https://docs.vast.ai/guides/instances/choosing/instance-types), [reserved](https://docs.vast.ai/guides/instances/choosing/reserved-instances), [rental-types FAQ](https://docs.vast.ai/guides/reference/faq/rental-types) |
| V3 | `create instance` defaults to **on-demand at dph_total** unless `--bid_price` is passed, even after a `--type bid` search. When outbid, the instance goes to `stopped` and storage keeps billing. | Scripts must pass `--bid_price` explicitly for spot. | [SKILL.md](https://raw.githubusercontent.com/vast-ai/vast-cli/master/vastai/SKILL.md) |
| V4 | Stopping releases the GPU. Restart waits in "Scheduling", possibly forever, if someone else rented it. | Never stop to "pause". | [manage instances](https://docs.vast.ai/guides/instances/manage-instances), [troubleshooting](https://docs.vast.ai/guides/reference/troubleshooting) |
| V5 | Every offer has a max duration. At expiry the instance stops. Expired instances "may be deleted 48 hours after expiration". | Filter `duration ≥ 1 day` (snapshot median 42 days). | [instance types](https://docs.vast.ai/guides/instances/choosing/instance-types), [instances FAQ](https://docs.vast.ai/guides/reference/faq/instances) |
| V6 | Container disk size is fixed at creation and cannot be resized. Volumes are **local to one machine**. | Size the disk per model (§6.4). Volumes are useless across hosts; ours are short-lived. | [storage types](https://docs.vast.ai/guides/instances/storage/types), [volumes](https://docs.vast.ai/guides/instances/storage/volumes) |
| V7 | CPU and RAM baselines are proportional to the GPU fraction of the machine. Bursting is allowed. RAM above the allocation risks the OOM killer. `/dev/shm` is proportional too. | Filter on `cpu_cores_effective ≥ 8` and RAM per GPU. Keep dataloader workers within the baseline. | [Docker environment](https://docs.vast.ai/guides/instances/docker-environment) |
| V8 | Three launch modes:<br>• **ssh** and **jupyter** replace the image entrypoint with Vast's (use `onstart`);<br>• **entrypoint/args** runs the image as is;<br>• `--direct` gives a direct SSH port. The proxy SSH is slow for large transfers (use it for < 1 GB only). | ssh mode with `--direct`; data over direct SSH. | [connect overview](https://docs.vast.ai/guides/instances/connect/overview), [data movement](https://docs.vast.ai/guides/instances/storage/data-movement) |
| V9 | Env vars:<br>• `-e` in the docker options;<br>• account-level env vars are injected into every container;<br>• in ssh mode they are not visible in SSH sessions unless written to `/etc/environment`;<br>• `CONTAINER_API_KEY` and `CONTAINER_ID` let the container stop or destroy **itself**. | On-box watchdog can self-destroy (§7). Never put secrets in templates. | [Docker environment](https://docs.vast.ai/guides/instances/docker-environment), [template settings](https://docs.vast.ai/guides/templates/template-settings) |
| V10 | Private images: the template's registry login (`--login`) is used by the **host** to pull. | Use a read-only, single-repo registry token (GHCR `read:packages`). | [template settings](https://docs.vast.ai/guides/templates/template-settings), CLI `create instance --help` |
| V11 | `vastai copy` (rsync over SSH; local ↔ instance, instance ↔ instance, volume ↔ instance) works even on stopped instances. Cloud Sync (S3, GDrive, Backblaze, Dropbox, HF buckets) **copies the cloud credentials onto the host**, and Vast advises Secure Cloud plus scoped keys. | Use rsync from the Mac. No Cloud Sync. | [data movement](https://docs.vast.ai/guides/instances/storage/data-movement), [cloud sync](https://docs.vast.ai/guides/instances/storage/cloud-sync), [CLI copy](https://docs.vast.ai/cli/reference/copy) |
| V12 | **Vast takes its own backups** of instances and volumes, for example before a machine is decommissioned. "A backup outlives the instance it came from" and can be downloaded with tokened URLs. | Copies of whatever was on the disk can exist after we destroy. This matters for RIGA (§5). | CLI 1.8.3 `vastai show instance-backups --help` ([vast-cli](https://github.com/vast-ai/vast-cli)) |
| V13 | Verification is automated: reliability, infrastructure, DLPerf, supply and demand. The **Datacenter / Secure Cloud** label needs:<br>• ISO/IEC 27001;<br>• a registered business;<br>• a signed Datacenter Hosting Agreement;<br>• owner identity verification;<br>• at least 5 servers.<br>It is filterable as `datacenter=True`. | The governance tier for RIGA (§5). | [understanding verification](https://docs.vast.ai/host/understanding-verification), [datacenter status](https://docs.vast.ai/host/datacenter-status), [CLI search](https://docs.vast.ai/cli/reference/search-offers) |
| V14 | Hosts' data access, per Vast:<br>• "Hosts can technically access files on their machines";<br>• containers are isolated from *other renters*, but protection from *providers* "varies significantly";<br>• advice: Secure Cloud, encrypt at rest, no credentials on instances.<br>The ToS makes providers independent contractors, has the **renter** warrant all rights to the data, and contains no host-confidentiality clause. | §5 | [manage instances](https://docs.vast.ai/guides/instances/manage-instances), [security FAQ](https://docs.vast.ai/guides/reference/faq/security), [ToS (2026-09-01)](https://vast.ai/terms) |
| V15 | Balance at $0: instances are stopped. With a saved card, the card is auto-charged. With no card, instances and data are **destroyed** after a short grace period. New accounts have a `spend_rate_limit` that stays very small until the email is verified. | Guardrail design (§7). Verify the email before session 1. | [billing](https://docs.vast.ai/guides/reference/billing), [troubleshooting](https://docs.vast.ai/guides/reference/troubleshooting) |
| V16 | DLPerf is Vast's per-machine estimate of deep-learning throughput (for example V100 ≈ 21, 2080 Ti ≈ 14). Reliability starts at 60 % for new machines. | Use DLPerf only as a relative proxy (§3). Filter reliability > 0.98. | [rental-types FAQ](https://docs.vast.ai/guides/reference/faq/rental-types), [concepts](https://docs.vast.ai/guides/concepts) |
| V17 | Docker-in-Docker is unsupported. Instances are unprivileged containers (a minority are VMs). | No nested builds and no LUKS or FUSE assumptions on the box. | [instances FAQ](https://docs.vast.ai/guides/reference/faq/instances), [security FAQ](https://docs.vast.ai/guides/reference/faq/security) |
| V18 | API keys can be **scoped**: `instance_read`, `instance_write`, `misc` (search), `billing_*`, and so on, with per-ID constraints. Webhooks and email fire on `client:low_credit` and `client:outbid`. | The Mac's automation key gets no billing scope. Low-credit alert. | [CLI permissions](https://docs.vast.ai/cli/permissions), [notification webhooks](https://docs.vast.ai/guides/reference/notification-webhooks) |
| V19 | Status values: `exited`, `unknown` and `offline` never reach `running`, so poll loops need a timeout plus destroy. Storage charges start at creation. | `launch.sh` polls with a timeout. | [SKILL.md](https://raw.githubusercontent.com/vast-ai/vast-cli/master/vastai/SKILL.md) |
| V20 | `vastai run benchmarks --template_hash … --gpus RTX_4090,RTX_5090,…` rents one instance per GPU, measures, and tears down (it costs money). `search benchmarks` needs login. | An option for H-STEP, but its "perf" metric is serverless-oriented. A plain scripted run is more transparent. | CLI 1.8.3 `run benchmarks --help` |
| V21 | Serverless (endpoints, workergroups, PyWorker) autoscales inference workers and is billed as plain instances. | Not useful for Phase 2. Possibly for Phase 9 deployment demos. | [concepts](https://docs.vast.ai/guides/concepts), [pricing](https://docs.vast.ai/guides/instances/pricing) |

---

## 2. How we use Vast: step by step

**One-time account setup** (coordinator, before Phase 2; RESULTS §8 says the account must be funded first):
1. Create the account and **verify the email**. This lifts the new-account `spend_rate_limit` (V15).
2. Add the Mac's ed25519 public key at the [keys page](https://cloud.vast.ai/manage-keys/)
   ([SSH guide](https://docs.vast.ai/guides/instances/connect/ssh)).
3. Create a **scoped API key** for the scripts, with `misc`, `instance_read`, `instance_write` and `user_read` and no
   billing or team scopes, and store it in the macOS keychain (V18). The full-power key stays in the browser only.
4. Set up money:
   - prepaid credit in small tranches ($25 to start; §7);
   - a low-credit email or webhook at $10;
   - decide card versus no card (§7, G1).
5. Container registry: a private GHCR package for the image, plus a read-only PAT (`read:packages`) used only as the
   Vast `--login` (V10).
6. Write the `.vast-allowlist`: the offer filters (§3.3) and the governance tier per dataset (§5).

**Per session** (scripts proposed in §6.5):
1. `offers.sh <gpu> [--secure]`: runs `vastai search offers` with the §3.3 filter and keeps the cheapest 5. A human
   picks one; the script records the offer JSON (host id, machine id, datacenter flag, driver, price) in the run
   ledger.
2. `launch.sh <offer> <run_id>`:
   - runs `vastai create instance <offer> --image ghcr.io/<owner>/objdel-vast:deps-<pinhash> --login '…' --ssh --direct --disk <GB> --label <run_id> --onstart-cmd '<watchdog>' --cancel-unavail`;
   - polls `show instance --raw` until `running`, with a 15-minute timeout; on `exited`, `unknown` or `offline` it
     destroys (V19).
3. `sync_up.sh`:
   - rsyncs the code tarball (`git ls-files`, as `x86.sh` makes it; no `data/`, `.venv`, `build/`) and the dataset
     bundle the governance tier allows, over direct SSH;
   - the box verifies SHA-256 against `data/prepared/manifest.json`.
4. On the box, in tmux:
   - install the code editable (about 20 s on x86; M9);
   - run `image_check.sh`;
   - download the model at a pinned HF revision;
   - start `run.sh <config>`, which wraps training in `timeout <cap>` and writes a heartbeat file.
5. `pull.sh <run_id>`:
   - runs from the Mac every N minutes (launchd) and at the end;
   - rsyncs `runs/<run_id>/` (logs, trainer_state, eval jsonl, last plus best checkpoint) to the Mac;
   - verifies checksums.
6. `teardown.sh <run_id>`:
   - final pull and checksum;
   - `vastai destroy instance <id> -y`;
   - `show instances` must not list it;
   - append `show invoices-v1` charges to the ledger.

**Run types:**
- **Screening SFT:** on-demand.
- **Eval K-sampling of finished checkpoints:** interruptible is fine. The job is idempotent and resumes per item.
- **Phase 3+ GRPO seeds:** on-demand for decisive runs; interruptible with `save_steps` plus `resume_from_checkpoint`
  only for replicates.

---

## 3. GPUs

### 3.1 Specifications (primary sources)

Dense BF16 is tensor-core throughput with FP32 accumulate, without sparsity. GeForce cards run FP32-accumulate at half
the FP16-accumulate rate (whitepaper tables).

| GPU | Arch (cc) | Memory | Mem BW (datasheet) | Dense BF16 TFLOPS | FP8 | FlashAttention | Source |
|---|---|---|---|---|---|---|---|
| RTX 3090 | Ampere (8.6) | 24 GB GDDR6X | 936 GB/s | 71 | no | FA2 yes | [GA102 whitepaper](https://www.nvidia.com/content/PDF/nvidia-ampere-ga-102-gpu-architecture-whitepaper-v2.pdf) (Table 9) |
| RTX 4090 | Ada (8.9) | 24 GB GDDR6X | 1,008 GB/s | 165.2 | yes (330 FP32-acc) | FA2 yes | [Ada whitepaper](https://images.nvidia.com/aem-dam/Solutions/geforce/ada/nvidia-ada-gpu-architecture.pdf) (App. A) |
| RTX 5090 | Blackwell (12.0) | 32 GB GDDR7 | 1,792 GB/s | 209.5 | yes, plus FP4 | FA2 not listed for sm_120; SDPA works | [RTX Blackwell whitepaper](https://images.nvidia.com/aem-dam/Solutions/geforce/blackwell/nvidia-rtx-blackwell-gpu-architecture.pdf), [5090 page](https://www.nvidia.com/en-us/geforce/graphics-cards/50-series/rtx-5090/) |
| RTX PRO 6000 Blackwell | Blackwell (12.0) | 96 GB GDDR7 ECC | 1,792 GB/s | ≈ 500, **derived** (4,000 sparse FP4 TOPS ÷ 8, no GeForce halving; UNVERIFIED) | yes, plus FP4 | as 5090 | [RTX PRO 6000 page](https://www.nvidia.com/en-us/products/workstations/professional-desktop-gpus/rtx-pro-6000/) |
| L40S | Ada (8.9) | 48 GB GDDR6 ECC | 864 GB/s | 362 | yes (733) | FA2 yes | [L40S page](https://www.nvidia.com/en-us/data-center/l40s/) |
| A100 SXM4 80 GB | Ampere (8.0) | 80 GB HBM2e | 2,039 GB/s | 312 | no | FA2 yes | [A100 page](https://www.nvidia.com/en-us/data-center/a100/) |
| H100 SXM | Hopper (9.0) | 80 GB HBM3 | 3.35 TB/s | 989 (1,979 sparse) | yes | FA2, FA3, FA4 | [H100 page](https://www.nvidia.com/en-us/data-center/h100/) |

FlashAttention support:
- FA2 supports "Ampere, Ada, or Hopper GPUs", with bf16 on those.
- FA3 targets Hopper.
- FA4 (CuTeDSL) targets "Hopper and Blackwell (e.g. H100, B200)".
- Source: [flash-attention README](https://github.com/Dao-AILab/flash-attention).
- A2 already recommends SDPA by default and Hub FA2 kernels only if profiling shows attention matters (A2 §stack).
  FlashAttention is therefore **not** a GPU selection criterion.

Software compatibility:
- PyPI `torch==2.14.1` for Linux x86_64 depends on `cuda-toolkit==13.0.3` and cuDNN cu13
  ([PyPI JSON](https://pypi.org/pypi/torch/2.14.1/json)).
- CUDA 13.0 needs Linux driver **≥ 580**, and CUDA 13 dropped pre-Turing GPUs
  ([CUDA release notes](https://docs.nvidia.com/cuda/cuda-toolkit-release-notes/index.html)).
- torch 2.14.1 also ships cu126 and cu132 builds ([download index](https://download.pytorch.org/whl/torch/)). cu126
  would accept drivers ≥ 560 but predates Blackwell support (sm_100 and sm_120 need CUDA ≥ 12.8; architecture lists
  in §6.3).

### 3.2 Vast snapshot, 2026-10-19 about 24:00 UTC

Filters:
- verified, single GPU, `cpu_cores_effective ≥ 8`, reliability > 0.98;
- disk ≥ 50 GB, CUDA ≥ 12.8, inet down ≥ 200 Mb/s;
- prices are `dph_total` including 100 GB of storage.

The snapshot is in `experiments/phase2/lab/A3/vast_snapshot.json`; queries are listed in its `filters` field. It is
one moment of a moving market. Re-run `snapshot_offers.py` before each phase.

| GPU | On-demand n · min / median / max $/h | Interruptible min / median | Secure Cloud on-demand n · min–max | Median eff. vCPUs | Median DLPerf | Measured mem BW (median) |
|---|---|---|---|---|---|---|
| RTX 3090 | 21 · 0.20 / 0.32 / 0.82 | 0.15 / 0.19 | 2 · 0.39–0.40 | 16 | 44 | 805 GB/s |
| RTX 4090 | 23 · 0.51 / 0.66 / 1.57 | 0.30 / 0.49 | 1 · 1.57 | 16 | 97 | 879 GB/s |
| RTX 5090 | 37 · 0.63 / 1.02 / 2.29 | 0.40 / 0.75 | **9 · 0.80–2.29** | 24 | 198 | 1,455 GB/s |
| RTX PRO 6000 WS / Server | 17 · 1.44 / 1.60 / 2.05 and 10 · 1.36 / 1.56 / 1.82 | 0.56–0.69 / 1.03–1.36 | 1 + 1 · 1.61–2.05 | 32–34 | 280–283 | 1,400 GB/s |
| A100 SXM4 (40 + 80 GB) | 14 · 0.48 / 0.92 / 1.32; **80 GB only: 7 · 1.09 / 1.19 / 1.30** | 0.34 / 0.53 | 4 · 0.69–1.19 | 28 | 107–129 | 1,317 GB/s |
| L40S | 6 · 0.56 / 0.82 / 0.86 | 0.43 / 0.60 | 0 | 24 | 94 | 625 GB/s |
| RTX 6000 Ada / A6000 | 2 · 0.67–0.68 / 1 · 0.70 | — | 0 | 15–24 | 113 / 54 | — |
| H100 SXM | 2 · 2.17 / 3.75 / 5.32 | 1.60 | 0 | 25 | 320 | 2,880 GB/s |
| H200 NVL (141 GB) | 5 · 3.74 / 4.10 / 4.34 | 2.79 / 3.50 | 0 | 64 | 452 | 3,883 GB/s |
| B200 | 3 · 6.28 / 8.73 / 10.76 | 5.96 / 7.18 | 1 · 10.76 | 24 | 681 | 6,080 GB/s |

Other observations from the same pull (raw per-offer JSON, 1,500+ offers, `rentable` and not `rented`):
- **16 of 104 single "RTX 4090" offers report 48 GB of VRAM**: modified cards, some in CN. Filter `gpu_ram ≤ 26` for
  a stock 4090 so that memory results transfer.
- 85–94 % of 3090/4090/5090 offers, 83 % of A100 SXM4 and 100 % of RTX PRO 6000 report `cuda_max_good ≥ 13.0`, which
  the driver ≥ 580 our wheel needs.
- Only 5 % of offers in the pull are datacenter-hosted.
- Bandwidth in eligible offers:
  - median download 951 Mb/s, upload 925 Mb/s;
  - upload cost median $8/TB, p90 $40/TB;
  - storage median $0.27 per GB-month, p90 $0.67.
- 51 % of offers have a static IP. 100 % have direct ports.

### 3.3 The search filter we standardise on

```
vastai search offers 'num_gpus=1 gpu_name=RTX_4090 gpu_ram<=26 verified=true rentable=true reliability>0.98
  cpu_cores_effective>=8 cpu_ram>=32 cuda_vers>=13.0 disk_space>=80 inet_down>=500 direct_port_count>=2
  duration>=1 dph<=0.90' -o 'dph_total' --storage 80 --raw            # add datacenter=True for RIGA (§5)
```

Fields and syntax are from the [search-offers reference](https://docs.vast.ai/cli/reference/search-offers). Notes:
- `cpu_cores_effective ≥ 8` because the reward scores 256 rollouts per step on the CPU (ROADMAP "When renting").
- `cpu_ram ≥ 32` because the RIGA cache LRU holds about 9 MB per item and the prepare step peaks at 1.8–2.4 GB (M7).

### 3.4 Selection: five GPUs

The estimates in the table come from `estimate_cost.py`. Assumptions:
- 2 epochs over every voting-rater target: LIDC 23,469 plus RIGA 3,282 rows;
- MFU 0.10–0.30;
- plus 0.25 h of eval sampling per run;
- on-demand min–median prices from §3.2.

| Pick | GPU | Phase 2 role | Phase 3–4 role | Per-run hours and $ | Why | Risks |
|---|---|---|---|---|---|---|
| 1 | **RTX 4090 24 GB** | Default screening card for SmolVLM-256M and LFM2.5-VL-450M (full FT fits: about 6 / 10 GB). The RESULTS §8 cost gate is defined on it. | SmolVLM / LFM GRPO (about 10–20 GB, ROADMAP table) | Smol 0.6–1.2 h, $0.3–0.8; LFM 0.7–1.5 h, $0.3–1.0; Qwen LoRA 2.4–6.6 h, $1.1–4.4 | Mature Ada stack (FA2, FP8), 23 eligible offers, 16 median vCPUs | Modified 48 GB cards; consumer reliability; Qwen full FT does not fit |
| 2 | **RTX 5090 32 GB** | Same models, 2× the 4090's DLPerf (198 vs 97) at about 1.5× the price, so the best DLPerf per $. **The RIGA card**: 9 Secure Cloud offers, the most of any type. | SmolVLM / LFM GRPO with headroom; Qwen LoRA GRPO | Smol 0.5–1.0 h, $0.3–1.0 | Blackwell sm_120 tests the CUDA 13 stack early; 1.8 TB/s helps small, bandwidth-bound models | FA2 not listed for sm_120 (use SDPA); youngest driver stack |
| 3 | **RTX PRO 6000 Blackwell 96 GB** | Qwen3-VL-2B full SFT (about 38 GB) with room for larger batches | **Qwen 2B GRPO with vLLM colocate on one card** (50–60 GB plus the KV cache); LFM / Qwen best-of-8 sampling | Qwen 0.95–2.35 h, $1.3–3.7 | 27 offers at $1.36–2.05, cheaper than an A100 80 GB per DLPerf (283 vs 129), 32 vCPUs | Derived BF16 peak; Blackwell stack as for the 5090 |
| 4 | **A100 SXM4 80 GB** | Fallback for Qwen full SFT if the Blackwell stack misbehaves. The ROADMAP's reference card. | Qwen GRPO (the ROADMAP plan); Secure Cloud RIGA runs on 80 GB | Qwen 1.4–3.6 h, $1.5–4.3 | Most mature stack (Ampere, FA2); 4 Secure Cloud offers; 28–32 vCPUs | No FP8; lower DLPerf per $ than the PRO 6000; 40 GB variants mixed into the listings (filter `gpu_ram ≥ 79`) |
| 5 | **RTX 3090 24 GB** | First session: x86 K3/K4 and `xarch_hashes`, image smoke test, data-path dry runs. SmolVLM dev and debug. | Phase 8 CUDA checker development ("cheapest NVIDIA card", ROADMAP) | Smol 1.0–2.5 h, $0.2–0.8 | $0.20–0.33 on-demand, $0.15 spot; Ampere bf16 (`is_torch_bf16_gpu_available`, RESULTS §8) | Slowest; some hosts older |

Rejected:
- **L40S / RTX 6000 Ada / A6000 (48 GB):** thin supply (1–6 offers) and no Secure Cloud. The L40S matches the A100's
  DLPerf, but its 48 GB is tight for Qwen GRPO. The PRO 6000 dominates at 96 GB.
- **H100 / H200 / B200:**
  - 1.5–8× the throughput but $2.2–10.8/h;
  - only 1–5 single-GPU offers;
  - no Secure Cloud single-GPU offers in the snapshot.
  
  Phase 2 jobs are short enough that wall-clock does not bind. Reconsider for Phases 3–4 only if a decisive-arm batch
  (about 50 GPU-h) must finish within a day: 50 h on an H100 at $2.2–3.8 costs $110–190, against 50–80 h on a
  PRO 6000 at $1.4–1.6, which costs $70–130.
- **Multi-GPU:** our models fit on one card. Single-process TRL also avoids I22's group-split check.

**Phase 3–4 implications:**
- GRPO step time is dominated by generation (K12 on the Mac: 13–17 s generation plus scoring, 28–33 s update per step).
- Generation is memory-bandwidth-bound. The PRO 6000 (1.8 TB/s, 96 GB for vLLM's KV cache) and the H100 (3.35 TB/s)
  gain most.
- The reward needs at least 8 vCPUs. All five picks have a median of 16 or more.
- At an assumed $1.5/h on a PRO 6000:
  - Track A (about 70 GPU-h, ROADMAP) is about $105;
  - the full ~150 GPU-h program is about $225 on-demand, roughly half that with interruptible replicates.

---

## 4. Cost estimate for Phase 2

| Scenario | 4090 | 5090 | PRO 6000 | A100 80 GB | 3090 |
|---|---|---|---|---|---|
| Phase 2 package if SmolVLM wins: 3 screening runs plus A0 and A0′ on the winner (5 runs) | 4.8–11.8 GPU-h, $2.3–7.8 | 4.0–9.6 h, $2.3–9.2 | 2.4–4.7 h, $3.3–7.4 | 3.1–6.8 h, $3.4–8.1 | 9.4–26 h, $2.3–8.5 |
| Same if Qwen3-VL-2B is needed (Qwen on 24–32 GB means LoRA) | 8.3–22.5 h, $4.0–14.9 | 6.8–18 h, $3.8–17.3 | 3.6–8.3 h, $4.9–13.0 | 5.0–12.5 h, $5.5–14.9 | — |

Overheads:
- About 0.2–0.4 h per session: image pull (3.6 GB at the median 950 Mb/s is ≈ 30 s plus extraction), data upload, setup.
- Bandwidth: pulling 20 GB of Qwen checkpoints at the p90 $40/TB costs $0.80.
- Storage: 100 GB at $0.27 per GB-month is $0.04/h.

These are negligible. The ROADMAP's ~6 GPU-hours is plausible for the SmolVLM path. A2's independent estimate is in the
same range: SmolVLM 7 min per epoch, Qwen full FT 85 min per epoch on a 4090 (A2 §1, which counts LM tokens only). The
first session's measured s/step replaces all of this (RESULTS §8 "Cost").

---

## 5. Data governance: RIGA and LIDC on rented machines

### 5.1 What the licences say

- **LIDC-IDRI:**
  - licence CC BY 3.0; citation DOI 10.7937/K9/TCIA.2015.LO9QL9SX; an NCI and FNIH acknowledgement is required;
  - the TCIA policy forbids re-identification and requires attribution;
  - it does not restrict cloud or third-party storage;
  - sources: [TCIA LIDC-IDRI](https://www.cancerimagingarchive.net/collection/lidc-idri/),
    [TCIA data usage policy](https://www.cancerimagingarchive.net/data-usage-policies-and-restrictions/).
- **RIGA (Almazroa 2018, Deep Blue, DOI 10.7302/Z23R0R29):**
  - The **registered rights in the DataCite metadata are CC BY-NC 4.0**
    ([DataCite API record](https://api.datacite.org/dois/10.7302/Z23R0R29)).
  - The Deep Blue landing page (403 to scripts; quoted in Phase 1 via a search snippet:
    `docs/research/notes/Phase 1 stack and algorithms/data_engineering.md`) says "free for research and educational
    purposes; copy, redistribution, and any unauthorized commercial use are prohibited".
  - That sentence is **verbatim the MESSIDOR terms**: "Data included in this database can be used, free of charge,
    only for research and educational purposes. Copy, redistribution, and any unauthorized commercial use are
    prohibited."
  - MESSIDOR adds that "the use of this database is restricted to individuals or organizations that obtained the
    database directly from this website" ([ADCIS Messidor](https://www.adcis.net/en/third-party/messidor/)).
  - RIGA's MESSIDOR file holds 460 of the 750 originals ([DataCite](https://api.datacite.org/dois/10.7302/Z23R0R29)).
  - Our train split is MESSIDOR plus BinRushed (ROADMAP Phase 1; PLAN D13), so **most RIGA training images are
    MESSIDOR-derived**.
- **CC BY-NC 4.0** defines "Share" as providing material "to the public". Reproduction for NonCommercial purposes is
  granted ([legal code](https://creativecommons.org/licenses/by-nc/4.0/legalcode.en)).

### 5.2 Is a rented host "redistribution"?

- **Under CC BY-NC (BinRushed, Magrabia, the annotations):** probably not.
  - A private, rented, single-tenant container is not "the public".
  - Making working copies for our NonCommercial research is granted.
- **Under the MESSIDOR terms (460 of the originals):** this is the risk.
  - The terms forbid "copy" and restrict use to direct obtainers.
  - A **community host is an independent contractor** (ToS). It has root on the machine and "can technically access
    files" (V14).
  - It has signed no confidentiality obligation towards us. Vast's DPA lists no hosts as subprocessors
    ([DPA](https://vast.ai/data-processing-agreement)).
  - Vast itself may keep **backups that outlive the instance** (V12).
  - Placing the images there arguably gives an unvetted third party access to a copy. That is the
    "redistribution / non-direct obtainer" situation the terms target, even though nothing is published.
- **Secure Cloud hosts change the risk, not the principle.** They are ISO 27001-certified businesses with verified
  owners and a signed Datacenter Hosting Agreement (V13). That is comparable to an ordinary cloud processor, which
  research groups routinely use for licence-restricted data. Root access still exists.
- **Encryption:**
  - We can store an encrypted tarball (age or GPG) on the container disk and decrypt only into `/dev/shm` or process
    memory.
  - This protects against disk remnants after destroy, against Vast backups and against casual browsing.
  - It does **not** protect against a malicious root on a live machine: the key and the plaintext are in memory.
  - LUKS and FUSE (gocryptfs) usually need privileges that unprivileged containers lack (V17).
- **Confidential computing** (H100 or Blackwell TEE): not offered on Vast as far as the docs show. UNVERIFIED.

### 5.3 Options

| Option | Licence risk | Cost / effort | Verdict |
|---|---|---|---|
| A. LIDC anywhere; RIGA never leaves the Mac. RIGA SFT on MPS: 3,282 rows × 2 epochs at about 1.2 s/row (K12: 4.8 s per 4-row LoRA step) ≈ 2.2 h for SmolVLM. | none | Mac-bound. Qwen full FT on RIGA is infeasible on a 16 GB Mac. MPS fp32 ≠ CUDA bf16, so arms are not comparable. | Fallback only |
| B. LIDC anywhere; **RIGA only on Secure Cloud**. Encrypted tarball → `/dev/shm`; no Cloud Sync; destroy; ledger entry per session (host id, datacenter flag, file hashes, destroy time). | low (residual: datacenter root access, Vast backups of the encrypted blob) | Secure Cloud premium: 5090 $0.80–0.96 vs $0.63; A100 $0.69–1.19. Thin supply: 1–9 offers per type. | **Recommended default** |
| C. Ask permission: email the RIGA depositor (Almazroa, via Deep Blue) and the Messidor partners (Decencière, Laÿ) asking whether processing on rented cloud GPUs counts as permitted use. | eliminates it if granted | One email. Unknown latency. | **Do now, in parallel with B** |
| D. Split RIGA: BinRushed and Magrabia images (CC BY-NC 4.0) on any verified host; MESSIDOR images on Secure Cloud or the Mac. | low | More plumbing, and only 290 of 750 originals are non-MESSIDOR | Only if B's supply fails |
| E. RIGA on community hosts with encryption | **high** (V14, V12) | cheap | Reject |

**Recommendation.**
- Phase 2 model screening and the A0 vs A0′ arms use **LIDC only**. LIDC carries the topology claim; RIGA tests the
  band (ROADMAP "Phase 0 changed").
- RIGA SFT and eval runs go to **Secure Cloud** hosts under option B, preferably 5090 or A100 80 GB.
- Send option C's emails before Phase 2 starts.
- Keep two never-mix data bundles:
  - `bundle-lidc.tar` for any host;
  - `bundle-riga.tar.age` for Secure Cloud only.
- `sync_up.sh` refuses to send the RIGA bundle to an instance whose offer JSON lacks `hosting_type = 1` /
  `datacenter = True`.
- **The x86 K3/K4 harness reads `data/fixtures/riga_experts.npz`.** On community boxes, run it LIDC-only plus with
  synthetic RIGA-like items (H-K3SYN).

---

## 6. Deployment architecture

### 6.1 Image and environments

```
 Mac (arm64)                                   registry (private GHCR)            Vast host (x86-64, NVIDIA driver ≥ 580)
 ───────────                                   ───────────────────────            ───────────────────────────────────────
 uv pip freeze (.venv) ─► requirements.txt ─►  objdel-vast:deps-<pinhash>  ──pull─► container (ssh mode, --direct)
 docker buildx --platform linux/amd64           (3.6 GB pull, 6.4 GB venv; rebuilt only when       /opt/venv  torch 2.14.1 (CUDA 13.0 wheels)
                                                pins change)                        /workspace/src  ◄─ rsync code tarball, `uv pip install -e .` (20 s)
 git ls-files tarball ──────────── rsync over direct SSH ─────────────────────────► /workspace/data ◄─ bundle-lidc (any) / bundle-riga.age (Secure Cloud)
 runs/<run_id>/ ◄──────────────── rsync pull (Mac-initiated, every N min) ───────── /workspace/runs/<run_id> (logs, ckpt last+best, eval jsonl)
 ledger.csv ◄── vastai show invoices-v1 / instance JSON
```

Image decisions:
- **Base: `python:3.12-slim-trixie` at the same multi-arch digest as `experiments/phase1/bench/Dockerfile`.** No CUDA
  base image is needed:
  - the PyPI torch wheel brings the CUDA 13.0 user-space libraries as pip packages
    ([torch 2.14.1 metadata](https://pypi.org/pypi/torch/2.14.1/json));
  - the host's NVIDIA runtime injects the driver.
  
  `NVIDIA_VISIBLE_DEVICES` and `NVIDIA_DRIVER_CAPABILITIES` are set in the image in case the runtime does not set them.
- **Pins:** exactly the Mac's `uv pip freeze`, the same file `x86.sh` uses. CPU reward hashes then match the Mac / x86
  comparison of M9. `reward.spec.fingerprint` includes numpy, scipy and OpenCV versions.
- **Code is not baked in.** It arrives per session as the same tarball `x86.sh` builds, then installs editable. A code
  change never rebuilds or re-pushes the large layers.
  - The run manifest records the git commit, the dirty-diff hash, the image digest, the host id and the driver.
  - The cibuildwheel manylinux wheel stays the CI release artifact. `Dockerfile.vast` can optionally install it.
- **Two venvs** once vLLM is needed: `/opt/venv` (train, torch 2.14.1) and `/opt/vllm` (vLLM ≤ 0.31.0, torch 2.13.0).
  - TRL 1.16 requires `vllm<=0.31.0,>=0.22.0` for its `[vllm]` extra
    ([TRL 1.16.0 metadata](https://pypi.org/pypi/trl/1.16.0/json)).
  - vLLM 0.31.0 and 0.32.0 pin `torch==2.13.0` ([vLLM 0.32.0](https://pypi.org/pypi/vllm/0.32.0/json),
    [0.31.0](https://pypi.org/pypi/vllm/0.31.0/json)).
  - Offline K-sampling of saved checkpoints works across venvs.
  - GRPO colocate needs one venv, so either Phase 3 pins torch 2.13.0 or TRL / vLLM move. **Decision for the
    coordinator** (H-VLLM).
- **Tools in the image:** openssh-server (ssh launch mode expects a "typical ssh daemon setup"; V8), rsync, tmux,
  procps, build-essential, cmake and ninja (objdel builds from source as in M9).

### 6.2 Data, models, checkpoints

- **Data bundles** come from `data/prepared/` (D13 / M7): jsonl, PNG images and caches.
  - LIDC images are 1.3 GB, jsonl 8 MB, caches 42 MB.
  - RIGA images are 148 MB.
  - The box verifies each file against `manifest.json` before training. A data hash mismatch aborts the run.
- **Base models** come from the HF Hub at a pinned commit (revision SHA in the run manifest). No HF token is needed for
  public models. If one is ever needed, it is a read-only fine-grained token sent as an env var at create time, never
  in a template (V9).
- **Checkpoints:**
  - TRL `save_steps` with `save_total_limit=2` (last two) plus best (the Controller's best-checkpoint rollback, D9).
  - The Mac pulls them every 15–30 min (launchd), so a host failure loses at most one interval.
  - For SmolVLM-256M full FT, a checkpoint with optimizer is about 3 GB. For Qwen 2B it is about 25 GB, so pull
    weights only, plus the final optimizer state.
  - Resume with `trainer.train(resume_from_checkpoint=...)` on a new box after rsyncing the checkpoint back.
- **Outputs never go to third-party clouds from the host.** If off-Mac durability is wanted, the **Mac** pushes to a
  private HF repo or B2.

### 6.3 Docker image test (H-IMG, this sandbox)

Build:
- `docker buildx build --platform linux/amd64` on a dedicated builder capped at 2 CPUs and 4 GB;
- the legacy builder failed on `COPY --from` with a platform image (`failed to get destination image`), so buildx it is.

Run: `docker run --platform linux/amd64 --cpus 2 --memory 4g --network none` with `image_check.sh` on a fresh source
snapshot (`git ls-files`, no `data/`). Results:

| Check | Result |
|---|---|
| Build | **OK** in 1,066 s, of which about 200 s was `uv pip install` of the pins and about 740 s was exporting and loading the image into Docker Desktop. Image content **3.61 GB** compressed; `/opt/venv` is **6.4 GB** on disk. A Vast pull at the median 950 Mb/s is about 30 s plus extraction. |
| torch | `2.14.1+cu130`, CUDA build 13.0, cuDNN 9.24.0 |
| Kernel archs (`torch._C._cuda_getArchFlags()`) | **sm_75 sm_80 sm_86 sm_90 sm_100 sm_120**. All five picks are covered: 3090 = 8.6; 4090 = 8.9, which runs the 8.6 cubins because cubins are binary-compatible within a major version; A100 = 8.0; 5090 and PRO 6000 = 12.0; plus H100 = 9.0 and B200 = 10.0. |
| Stack imports | transformers 5.20.1, trl 1.16.0, peft 0.21.2, accelerate 1.15.0, datasets 5.1.0, numpy 2.5.4, scipy 1.18.1, OpenCV 5.0.0, numba 0.68.0, `objdel._polycore`. Without network or GPU, `cuda.is_available()` is False, as expected. |
| objdel editable build from the snapshot | rc 0. It took 197 s under Rosetta; M9 measured about 20 s on x86. |
| `xarch_hashes` against the bench x86 CPU image (`objdel-x86:deps`, CPU torch), `--strict` | **26 / 26 sections identical: PASS.** The CUDA torch build does not touch the reward or eval path. |
| `xarch_hashes` against the Mac | 22 / 26 identical. The 4 accepted differences are the same `splits.*` ones as M9: BLAS sdot, FFT, ORB. **PASS.** |
| Timing | Emulation only: imports 120 s, hashes 109 s. **Not x86 timing**; K3/K4 stay a first-box gate. |

So H-IMG passes on everything checkable without a GPU. The GPU half (driver ≥ 580 handshake, bf16 matmul) is
first-session item 3.

### 6.4 Disk sizing

| Model | Container disk | Why |
|---|---|---|
| SmolVLM-256M, LFM2.5-VL-450M | 60–80 GB | image about 7 GB unpacked; data 1.5 GB; 2–3 checkpoints at about 3–5 GB; HF cache; eval outputs |
| Qwen3-VL-2B full FT | 150 GB | 2 checkpoints with optimizer at about 25 GB; merged weights; vLLM venv |

Disk cannot grow later (V6). Storage at about $0.27 per GB-month is negligible (§4).

### 6.5 Scripts (to write in `experiments/phase2/vast/`; proposal)

| Script | Does | Guardrail built in |
|---|---|---|
| `offers.sh` | §3.3 search → top 5 → offer JSON into the ledger | price ceiling per GPU; `gpu_ram` bounds; `cuda_vers ≥ 13.0`; `datacenter=True` when the RIGA bundle is requested |
| `launch.sh` | `create instance … --ssh --direct --cancel-unavail --label run_id`, then poll with a timeout | destroys on timeout or `exited` / `offline`; refuses a 2nd concurrent instance unless `--allow-parallel` |
| `sync_up.sh` | rsync code and the allowed bundle; remote SHA-256 verify | refuses RIGA to non-datacenter offers; refuses any path under `data/raw` |
| `onstart` (in `--onstart-cmd`) | env → `/etc/environment` (V9); starts the watchdog | watchdog self-destroys via `CONTAINER_API_KEY` when uptime exceeds the cap, or when the GPU has been idle for 30 min and `PULLED_OK` exists |
| `run.sh` | `timeout <cap> python -I -m experiments.phase2.<train> …` in tmux; heartbeat | wall-clock cap per config |
| `pull.sh` | rsync `runs/<run_id>/` to the Mac, verify, then `touch PULLED_OK` on the box | — |
| `teardown.sh` | final pull → destroy → verify absent → invoices into the ledger | — |
| `sweep_orphans.sh` (launchd, hourly) | `show instances --raw`: destroy anything not in `active_runs.txt` or older than its cap; alert | catches forgotten instances |

---

## 7. Cost guardrails

| # | Guardrail | Mechanism | Source |
|---|---|---|---|
| G1 | **Prepaid credit as the hard cap.** Top up in $25–50 tranches. No autobilling during Phase 2. | A saved card is auto-charged when the balance goes negative (V15), which defeats a cap. Without a card, instances are destroyed at $0 after a grace period. That is acceptable because everything of value is pulled to the Mac every 15–30 min. The coordinator may prefer the card for deletion protection; then set the autobill threshold to about one day's spend and the alert at 75 % of it ([billing](https://docs.vast.ai/guides/reference/billing)). | [billing](https://docs.vast.ai/guides/reference/billing) |
| G2 | Low-credit email or webhook at $10 | `client:low_credit` | [webhooks](https://docs.vast.ai/guides/reference/notification-webhooks) |
| G3 | Scoped automation key: no billing or team scopes | Permission JSON | [CLI permissions](https://docs.vast.ai/cli/permissions) |
| G4 | Destroy, never stop | `teardown.sh`; storage bills while stopped (V1) | [pricing](https://docs.vast.ai/guides/instances/pricing) |
| G5 | Per-run wall-clock cap: on-box `timeout` plus a self-destroying watchdog | `CONTAINER_API_KEY` (V9) | [Docker environment](https://docs.vast.ai/guides/instances/docker-environment) |
| G6 | Orphan sweep from the Mac, hourly | `sweep_orphans.sh` | — |
| G7 | Price ceiling in every search (`dph<=`); interruptible only with an explicit `--bid_price` (V3) | `offers.sh` | [SKILL.md](https://raw.githubusercontent.com/vast-ai/vast-cli/master/vastai/SKILL.md) |
| G8 | Poll loops with a timeout; destroy on `exited`, `unknown` or `offline` (V19) | `launch.sh` | [SKILL.md](https://raw.githubusercontent.com/vast-ai/vast-cli/master/vastai/SKILL.md) |
| G9 | Ledger: one row per session with offer JSON, start and stop times, `show invoices-v1` charges, datasets sent and the destroy confirmation | `experiments/phase2/vast/ledger.csv` (no secrets) | [CLI invoices](https://docs.vast.ai/cli/reference/show-invoices-v1) |
| G10 | Budget gate per phase. Phase 2 alarm at $40, about 3× the high estimate in §4. | The coordinator approves top-ups | — |

---

## 8. First-session checklist (one RTX 3090 or 4090, LIDC only; about 1.5 h, about $1)

1. **Account preconditions:** verified email, SSH key, scoped key, $25 credit, low-credit alert (§2).
2. **Choose and launch.**
   - `offers.sh RTX_4090`: stock 24 GB, ≥ 8 vCPUs, `cuda_vers ≥ 13.0`. Use a 3090 if the 4090 price is above $0.70.
   - `launch.sh` with a 60–80 GB disk.
   - Record host id, CPU model (`lscpu`), `nproc`, cgroup `cpu.max`, driver and `nvidia-smi -q` in the ledger.
3. **Image smoke test.** Run `image_check.sh`. On a GPU box it also prints the GPU, the cc, bf16 support and an 8192³
   bf16 matmul TFLOPS figure.
   - Pass: `torch.cuda.is_available()`;
   - the matmul reaches at least 60 % of the dense peak in §3.1 (sanity check on host throttling or a wrong card);
   - `objdel._polycore` imports.
4. **x86 determinism.**
   - Run `python -I -m experiments.phase1.bench.xarch_hashes --out vast.json`.
   - Pull it and run `--compare build/x86/hashes-mac.json vast.json` (RESULTS §8 step 2).
   - Also compare against the Docker x86 run. Any non-ACCEPTED difference blocks training.
   - If `reward.advantages` differs (AVX-512 SVML), rerun with `NPY_DISABLE_CPU_FEATURES="X86_V4"` and decide
     (M9 residual risk).
5. **K3/K4 x86 timing.**
   - Time 256 completions, serially and with `ThreadPoolExecutor(8)` (M6 open issue), on LIDC fixtures plus synthetic
     RIGA-like items (H-K3SYN).
   - Gate: batch ≤ 250 ms and p99 ≤ 5 ms (RESULTS §8).
   - Also log the batch time with the GPU busy, because training contends for the CPU.
6. **SFT step time (H-STEP).**
   - SmolVLM-256M, production SFT config, LIDC, 100 steps.
   - Log s/step (median after warm-up), rows/s, `torch.cuda.max_memory_allocated`, GPU utilisation (`nvidia-smi dmon`),
     CPU utilisation and dataloader wait.
   - Repeat with LFM2.5-VL-450M if time allows.
7. **GRPO step time** (the RESULTS §8 "Cost" gate): SmolVLM, G = 8, 32 prompts, `max_completion_length` 1,536, HF
   generation, production `TRLReward`; 5 steps. Log the generation, scoring and update split.
8. **Pull, destroy, verify, ledger.** Write the measured numbers into `estimate_cost.py` (replace the MFU range) and
   into ROADMAP §Compute (coordinator).

---

## 9. Hypotheses to test

| ID | Hypothesis | Test design | Pass criterion |
|---|---|---|---|
| **H-IMG** | The candidate image (`Dockerfile.vast`) builds for linux/amd64, carries CUDA 13.0 torch with sm_80 / 86 / 89 / 90 / 100 / 120 kernels, and its CPU reward path hashes identically to the validated x86 CPU image. | Mac sandbox (done, §6.3): build, `image_check.sh`, `xarch_hashes --compare` against Mac hashes and the bench x86 image. On the first box: the same script on a real GPU. | Build OK; arch flags include sm_80 and sm_120; compare shows only M9-ACCEPTED differences against the Mac and **none** against x86-CPU; on the GPU, `cuda.is_available()` and bf16 matmul ≥ 60 % of peak |
| **H-STEP** | Measured SFT throughput lies inside `estimate_cost.py`'s MFU 0.10–0.30 band. | First-session item 6 on the 4090. Then the same 100-step probe on the 5090 and PRO 6000 when they are first rented. | Measured s/step × steps per epoch within the band. Otherwise replace the band, and the GPU ranking follows measured $/epoch. |
| **H-5090** | On the RTX 5090 (sm_120, CUDA 13), our stack (torch 2.14.1 cu130, SDPA, transformers 5.20, TRL 1.16) trains SmolVLM and LFM without errors and with ≥ 1.4× the 4090's rows/s. | The same 100-step probe on both cards, same config and seed. | No kernel or arch errors; loss curves within run-to-run noise of the 4090; rows/s ratio ≥ 1.4 |
| **H-PRO** | One RTX PRO 6000 runs Qwen3-VL-2B full-FT SFT at batch ≥ 8 without OOM, and a Qwen GRPO step with vLLM colocate (G = 8, 32 prompts) fits in 96 GB. | 50 SFT steps; then 3 GRPO steps in the vLLM venv once H-VLLM is resolved. | Peak memory ≤ 85 GB; no OOM; s/step recorded |
| **H-VLLM** | vLLM (≤ 0.31.0, torch 2.13.0) can K-sample our SFT checkpoints from a separate venv, giving the same parse and legal rates as HF generate at matched sampling settings. | Generate K = 8 on 200 LIDC val items with both engines, fixed seeds, temperature 1. Compare format and legal rates with a two-proportion test. | Rates within 2 pp, p > 0.05. vLLM ≥ 3× faster in tokens/s. Otherwise keep HF generate for Phase 2 eval. |
| **H-K3SYN** | Synthetic RIGA-like items (as `xarch_hashes` generates: 6 raters, disc and cup) time like real RIGA fixtures in the K3 harness, so the x86 gate can run without RIGA data on community hosts. | On the Mac: K3 batch wall time on real `riga_experts.npz` items versus size-matched synthetic items, serial, 20 repetitions, interleaved. | Median ratio within 0.8–1.25 and p99 per rollout within 1 ms |
| **H-GOV** | `sync_up.sh` cannot send the RIGA bundle to a non-datacenter offer. | Unit test with recorded offer JSONs (`hosting_type` 0 and 1) and a dry-run rsync. | Refuses 100 % of non-datacenter offers; allows datacenter offers |
| **H-PULL** | A killed host loses at most one pull interval of work. Resume on a new box reproduces the loss trajectory. | On the 3090: SmolVLM 300 steps; destroy at step 200; relaunch; rsync the last checkpoint; resume. | Resumed loss at steps 201–300 within ±2 % of an uninterrupted control run (same seed, data order from the trainer state) |
| **H-WATCH** | The on-box watchdog destroys the instance at the cap. | Launch with cap = 10 min and no job. Confirm the instance is gone and billing stops. | `show instances` empty at 11 min; ledger charge ≤ 12 min |

---

## 10. Rejected alternatives

- **Vast Cloud Sync, or S3/B2 pushed from the host:** puts storage credentials on the host (V11). Mac pull is enough at
  our data sizes.
- **Volumes:** local to one machine (V6) and billed separately. Our sessions are short and do not return to the same
  machine.
- **Reserved instances:** a month-scale prepay with clawback (V2). Phase 2 needs hours.
- **Baking data into the image:** the registry and Vast snapshots (`take snapshot` pushes a container to a registry)
  are further third parties, and RIGA governance forbids it.
- **CUDA base images (`nvidia/cuda`, `vastai/pytorch`):** the pip wheels already carry CUDA. `vastai/pytorch` bundles
  its own torch, which conflicts with our pins
  ([advanced setup](https://docs.vast.ai/guides/templates/advanced-setup)), and its instance portal and Caddy are
  unneeded.
- **Jupyter launch mode:** for interactive notebooks. We run scripted jobs in tmux.
- **Serverless:** for autoscaled inference endpoints (V21).
- **Multi-GPU and H100 for Phase 2:** §3.4.

---

## 11. Sources

**Vast docs and CLI**
- docs index: https://docs.vast.ai/llms.txt
- concepts: https://docs.vast.ai/guides/concepts
- pricing: https://docs.vast.ai/guides/instances/pricing
- instance types: https://docs.vast.ai/guides/instances/choosing/instance-types
- reserved instances: https://docs.vast.ai/guides/instances/choosing/reserved-instances
- rental-types FAQ: https://docs.vast.ai/guides/reference/faq/rental-types
- security FAQ: https://docs.vast.ai/guides/reference/faq/security
- billing: https://docs.vast.ai/guides/reference/billing
- storage types: https://docs.vast.ai/guides/instances/storage/types
- volumes: https://docs.vast.ai/guides/instances/storage/volumes
- data movement: https://docs.vast.ai/guides/instances/storage/data-movement
- cloud sync: https://docs.vast.ai/guides/instances/storage/cloud-sync
- cloud backups: https://docs.vast.ai/guides/instances/storage/cloud-backups
- managing instances: https://docs.vast.ai/guides/instances/manage-instances
- Docker environment: https://docs.vast.ai/guides/instances/docker-environment
- connect overview: https://docs.vast.ai/guides/instances/connect/overview
- SSH: https://docs.vast.ai/guides/instances/connect/ssh
- template settings: https://docs.vast.ai/guides/templates/template-settings
- advanced setup: https://docs.vast.ai/guides/templates/advanced-setup
- instances FAQ: https://docs.vast.ai/guides/reference/faq/instances
- technical FAQ: https://docs.vast.ai/guides/reference/faq/technical
- troubleshooting: https://docs.vast.ai/guides/reference/troubleshooting
- notification webhooks: https://docs.vast.ai/guides/reference/notification-webhooks
- CLI permissions: https://docs.vast.ai/cli/permissions
- CLI search-offers: https://docs.vast.ai/cli/reference/search-offers
- CLI copy: https://docs.vast.ai/cli/reference/copy
- CLI show-invoices-v1: https://docs.vast.ai/cli/reference/show-invoices-v1
- understanding verification: https://docs.vast.ai/host/understanding-verification
- datacenter status: https://docs.vast.ai/host/datacenter-status
- vastai agent guide: https://raw.githubusercontent.com/vast-ai/vast-cli/master/vastai/SKILL.md
- vast-cli repository (CLI 1.8.3 `--help` for create instance, show instance-backups, take snapshot, run benchmarks):
  https://github.com/vast-ai/vast-cli
- offers API: `POST https://console.vast.ai/api/v0/bundles/`

**Vast legal**
- Terms of Service (dated 2026-09-01): https://vast.ai/terms
- DPA: https://vast.ai/data-processing-agreement

**Licences**
- RIGA (DataCite record): https://api.datacite.org/dois/10.7302/Z23R0R29
- RIGA (Deep Blue record, 403 to scripts): https://deepblue.lib.umich.edu/data/concern/data_sets/3b591905z
- MESSIDOR: https://www.adcis.net/en/third-party/messidor/
- TCIA LIDC-IDRI: https://www.cancerimagingarchive.net/collection/lidc-idri/
- TCIA data usage policy: https://www.cancerimagingarchive.net/data-usage-policies-and-restrictions/
- CC BY-NC 4.0 legal code: https://creativecommons.org/licenses/by-nc/4.0/legalcode.en

**GPUs**
- GA102 (Ampere) whitepaper: https://www.nvidia.com/content/PDF/nvidia-ampere-ga-102-gpu-architecture-whitepaper-v2.pdf
- Ada whitepaper: https://images.nvidia.com/aem-dam/Solutions/geforce/ada/nvidia-ada-gpu-architecture.pdf
- RTX Blackwell whitepaper: https://images.nvidia.com/aem-dam/Solutions/geforce/blackwell/nvidia-rtx-blackwell-gpu-architecture.pdf
- RTX 5090: https://www.nvidia.com/en-us/geforce/graphics-cards/50-series/rtx-5090/
- RTX PRO 6000: https://www.nvidia.com/en-us/products/workstations/professional-desktop-gpus/rtx-pro-6000/
- L40S: https://www.nvidia.com/en-us/data-center/l40s/
- A100: https://www.nvidia.com/en-us/data-center/a100/
- H100: https://www.nvidia.com/en-us/data-center/h100/
- FlashAttention README: https://github.com/Dao-AILab/flash-attention
- CUDA release notes: https://docs.nvidia.com/cuda/cuda-toolkit-release-notes/index.html
- Small-VLM throughput datapoint (222M nanoVLM, about 1.7M samples in about 6 h on one H100):
  https://huggingface.co/blog/nanovlm
- SmolVLM fine-tuning memory (QLoRA on an L4, about 16 GB): https://huggingface.co/blog/smolvlm

**Stack**
- torch 2.14.1 metadata: https://pypi.org/pypi/torch/2.14.1/json
- torch wheel index: https://download.pytorch.org/whl/torch/
- TRL 1.16.0 metadata: https://pypi.org/pypi/trl/1.16.0/json
- vLLM 0.32.0 metadata: https://pypi.org/pypi/vllm/0.32.0/json
- vLLM 0.31.0 metadata: https://pypi.org/pypi/vllm/0.31.0/json

**Repo**
- `experiments/phase1/RESULTS.md` §8
- `experiments/phase1/impl/M6_reward.md` (K3/K4)
- `experiments/phase1/impl/M7_train.md` (K12, token counts)
- `experiments/phase1/bench/`
- `docs/research/notes/Phase 2 SFT/A2_sft_stack.md` §1

"""B5 pre-registered hypotheses and pass criteria (fixed 2026-10-19 18:12 PDT, before the vLLM parity runs and the
container tests; the harness refused REPORT.md, so the text went to this file and is unchanged from the 18:12 draft).
The report (experiments/phase2/lab/B5/REPORT.md) quotes them and adds the results.

    python experiments/phase2/lab/B5/prereg.py    # prints the criteria
"""
PREREG = {
    "Q1/H-A5-7": {
        "candidates": {
            "a": "one venv: TRL + vLLM + training stack (Phase 3 colocate)",
            "b": "training venv as pinned + separate vLLM env (Phase 2 offline sampler; Phase 3 TRL vllm_mode=server)",
            "c": "SGLang next to the training pins",
        },
        "resolution_test": "uv pip compile, x86_64-manylinux_2_40, Python 3.12 (resolve_matrix.py)",
        "a_pass": [
            "a lock resolves with TRL 1.16, vLLM >= 0.30, transformers >= 5.18 (assistant-only loss on vision data) "
            "keeping numpy 2.5.4 and numba 0.68.0 (hash pins) unless hashes are shown unchanged at other versions; "
            "every uv override named and justified",
            "the linux/amd64 image of that lock imports vllm, vllm.LLM, trl.trainer.grpo_trainer, "
            "trl.generation.vllm_generation and objdel.train.* without error",
            "pytest -q -m 'not slow' (as x86_inside.sh): same pass/skip/fail counts as the pinned training venv, "
            "0 failures",
            "xarch_hashes --compare --strict against the pinned venv: identical in every section",
        ],
        "b_pass_phase2": "the two-venv image builds and both venvs import their stacks",
        "b_pass_phase3": [
            "(i) the trainer process needs no vLLM import",
            "(ii) weight sync works with trainer and server on one GPU (single-GPU Vast box, A3)",
            "(iii) TRL 1.16's client is wire-compatible with the vLLM version the server env resolves",
            "(static check: TRL 1.16 / vLLM 0.32 source)",
        ],
        "c_rule": "considered only if it resolves next to the training pins without overrides; even then a sampler "
                  "only (no TRL integration)",
    },
    "Q2/H-A5-1": {
        "design": "LFM2.5-VL-450M, Qwen3.5-0.8B, SmolVLM-256M x 32 val items per dataset (LIDC, RIGA), evenly "
                  "spaced over sorted ids; HF reference = training env (.venv, transformers 5.20.1, TRL 1.16 GRPO "
                  "prompt path, configs.processor_kwargs); vLLM 0.32.0 macOS CPU wheel (transformers 5.19.0) with "
                  "untokenized chat text + PIL image + mm_processor_kwargs = processor_kwargs minus padding_side",
        "pass": "prompt_token_ids == HF input_ids for 100 % of the 192 (model, item) pairs, and the SFT collator's "
                "token sequence starts with the GRPO generation prompt for 100 % of items (loss-span differences "
                "are reported, not failed)",
    },
    "Q2/H-A5-2a (CPU proxy)": {
        "design": "8 items (4 LIDC + 4 RIGA) x n = 8 per model, T = 1, top_p = 1, top_k = 0, max_tokens 96, "
                  "processed logprobs, base weights; fp32 for LFM and SmolVLM, bf16 for Qwen3.5 in both engines "
                  "(fp32 weights ~3.5 GB exceed this agent's 3 GB budget)",
        "pass_fp32": "median per-token |delta| <= 1e-3 nats and p99 per-sequence |sum delta| <= 0.5 nats",
        "pass_bf16": "median per-token |delta| <= 0.02 nats and p99 per-sequence |sum delta| <= 0.5 nats",
        "leak_check": "no Qwen generation_config default (T 0.7, top-k / top-p) leaks: it would show as a "
                      "systematic delta far above these bounds",
        "not_tested": "H-A5-2(b) distribution parity on an SFT checkpoint at K = 64 (needs a checkpoint and a GPU)",
    },
    "Q3/images": {
        "pass": "each image (Phase 2: train venv + vLLM sampler venv; Phase 3 candidate: one venv) builds for "
                "linux/amd64 on a builder capped at 2 CPUs / 4 GB; per venv torch.version.cuda == '13.0', the arch "
                "list includes sm_80 and sm_120, the stack imports without a GPU; sizes reported; images and builder "
                "deleted afterwards",
    },
}

if __name__ == "__main__":
    import json

    print(json.dumps(PREREG, indent=1))

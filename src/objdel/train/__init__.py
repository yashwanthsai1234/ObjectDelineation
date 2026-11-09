"""TRL integration (plan §4.6): import the submodules directly, so that using one does not import the others' stack.

dataset:   prepare frames / SFT targets / expert caches once (data/prepared); SFT and GRPO rows; cache_loader.
rewards:   TRLReward: reward v4 as weight-0 per-term loggers plus one aggregator; stop_ids.
callbacks: TraceStep (torch.profiler step), NvtxSteps (NVTX ranges, nsys window), BudgetStop (session spend).
serve:     preloaded expert caches, CacheRef (pickles as its key) and the warmed fork pool for reward scoring.
sample:    checkpoint -> K samples per item, completions table -> eval_run JSONL.
configs:   sft_config / grpo_config / lora / processor_kwargs for the three screened models (plan D8, D12).
"""

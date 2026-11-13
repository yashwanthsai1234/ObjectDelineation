"""The run specs and the stage builders the SFT and RL runner tests share."""
from tests.support.paths import ROOT

# the sft stages of the run specs (the old train_sft configs left git in E6.6)
SPECS = sorted((ROOT / "specs").glob("*.yaml"))
CPU_SFT = {"per_device_train_batch_size": 1, "gradient_accumulation_steps": 1, "use_cpu": True, "bf16": False,
           "logging_steps": 1, "dataloader_num_workers": 0, "dataloader_persistent_workers": False,
           "dataloader_prefetch_factor": None}  # tiny CPU runs of train_sft
GRPO_SMOKE = {"per_device_train_batch_size": 8, "gradient_accumulation_steps": 1, "per_device_eval_batch_size": 8,
              "max_completion_length": 8, "use_cpu": True, "bf16": False, "beta": 0.0, "use_vllm": False}


def stage(model, view, output_dir, **kw):
    """train_sft's SftStage in the old RunArgs order (model, view, output_dir)."""
    from objdel.core.spec import SftStage

    return SftStage(model=model, view=view, output_dir=output_dir, **kw)


def sft_stages():
    """(spec name, SftStage) for each sft stage of the run specs."""
    from objdel.core.spec import load_spec

    return [(p.stem, st) for p in SPECS for st in load_spec(p).stages if st.kind == "sft"]

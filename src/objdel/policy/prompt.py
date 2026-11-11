"""The user turn of a dataset (its entry's verbatim prompt, objdel.core.datasets), shared by the RL and SFT rows,
G-vision and the prepared-set manifest (PROMPTS: the table the builders store)."""
from objdel.core import datasets

PROMPTS = {name: d.prompt for name, d in datasets.DATASETS.items()}  # the library's datasets, derived at import


def _user(dataset: str) -> dict:
    return {"role": "user", "content": [{"type": "image"}, {"type": "text", "text": datasets.get(dataset).prompt}]}

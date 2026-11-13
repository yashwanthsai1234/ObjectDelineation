"""The user turn of a dataset (its entry's verbatim prompt, objdel.core.datasets), shared by the RL and SFT rows,
G-vision and the prepared-set manifest (PROMPTS: the table the builders store)."""
from objdel.core import datasets

PROMPTS = {name: d.prompt for name, d in datasets.DATASETS.items()}  # the library's datasets, derived at import


def user_turn(text: str) -> dict:
    """The chat user turn: the image, then the prompt text (RL, SFT and eval rows, and G-vision)."""
    return {"role": "user", "content": [{"type": "image"}, {"type": "text", "text": text}]}


def assistant_turn(text: str) -> dict:
    """The chat assistant turn of an SFT row: the target text."""
    return {"role": "assistant", "content": [{"type": "text", "text": text}]}


def _user(dataset: str) -> dict:
    return user_turn(datasets.get(dataset).prompt)

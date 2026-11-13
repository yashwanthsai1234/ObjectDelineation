"""GRPO rows from the prepared set: prompt_rows, grpo_rows and OpenOnce (objdel.policy.rows holds the records and
the image check)."""

from datasets import Dataset, Image

from objdel.policy.rows import as_dataset, prompts


def prompt_rows(root, dataset: str, split: str, rl_only: bool = True, oversample: int = 1) -> Dataset:
    """prompts() as a datasets.Dataset with decoded images (TRL GRPO)."""
    return as_dataset(prompts(root, dataset, split, rl_only, oversample))


def grpo_rows(root, dataset: str, split: str, oversample: int = 1) -> Dataset:
    """RL prompts: RL-eligible items, disagreement items repeated `oversample` times."""
    return prompt_rows(root, dataset, split, True, oversample)


class OpenOnce:
    """GRPO rows whose images are opened once per generation batch (perf: the 16 x 256 opens per step).

    TRL's RepeatSampler hands the DataLoader the same 256 indices (32 images x 8 answers) once per micro-batch, 16 times
    a step, and transformers fetches all 16 before the step (get_batch_samples); TRL generates from the first and reuses
    its answers for the other 15. A datasets.Dataset opens a PNG for every row of every fetch: 4,096 opens for 32 images
    (P1: 12.96 s a step with the GPU idle). OpenOnce reads the rows without their images, opens each distinct image once
    with the column's own decoder (the same pixels), shares it among the rows that name it, and hands the same list back
    while the indices repeat. TRL copies the messages it edits (prepare_multimodal_messages), so a shared row is never
    changed in place. Main process only (dataloader_num_workers=0): each worker would keep its own last stack."""

    def __init__(self, rows: Dataset):
        self.feature = rows.features["image"]
        self.rows = rows.cast_column("image", Image(decode=False))
        self.last: tuple[tuple, list] = ((), [])

    def __len__(self) -> int:
        return len(self.rows)

    def __getitem__(self, i: int) -> dict:
        return self.__getitems__([i])[0]

    def __getitems__(self, indices: list[int]) -> list[dict]:
        key = tuple(indices)
        if key and key == self.last[0]:
            return self.last[1]
        opened, out = {}, []
        for row in self.rows.__getitems__(list(indices)):
            v = row["image"]
            k = v["path"] or v["bytes"]
            if k not in opened:
                opened[k] = self.feature.decode_example(v)
            out.append(row | {"image": opened[k]})
        self.last = (key, out)
        return out

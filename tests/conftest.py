"""Shared pytest setup: data-backed fixtures are optional, so CI (which has no data) skips those tests; synthetic
items and a tiny offline TRL model shared by the training tests (test_train_*, test_breakit_wiring)."""
import os
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
FIXTURES = Path(os.environ.get("OBJDEL_FIXTURES", ROOT / "data" / "fixtures"))
STORE = Path(os.environ.get("OBJDEL_STORE", ROOT / "data" / "store"))


@pytest.fixture(autouse=True)
def _no_real_secrets(tmp_path_factory, monkeypatch):
    """No test reads the owner's <repo>/.secrets (objdel.secrets): a path that does not exist, and a fresh cache."""
    from objdel import secrets

    monkeypatch.setenv("OBJDEL_SECRETS", str(tmp_path_factory.mktemp("secrets") / ".secrets"))
    secrets.get.cache_clear()
    yield
    secrets.get.cache_clear()


@pytest.fixture(scope="session")
def fixtures_dir() -> Path:
    if not FIXTURES.is_dir():
        pytest.skip(f"fixtures not found at {FIXTURES}")
    return FIXTURES


@pytest.fixture(scope="session")
def store_dir() -> Path:
    if not STORE.is_dir():
        pytest.skip(f"data store not found at {STORE}")
    return STORE


def _disc(h, w, cy, cx, r, hole=0):
    y, x = np.ogrid[:h, :w]
    d2 = (y - cy) ** 2 + (x - cx) ** 2
    return (d2 <= r * r) & (d2 > hole * hole)


def _synthetic_items(splits=("train",), unanimous=False):
    """Per split (numbered s = 1, 2, ...): riga/MESSIDOR/image00<s> (3 raters, nested disc and cup) and LIDC items
    lidc/1.2.<s>/n000/z001 (raters disagree: a solid disc vs a ring with a hole), z002 (no pixels), z003 (one rater)
    and, with `unanimous`, z004 (two raters, one topology). Sorted as objdel.data.store.write needs."""
    from objdel.data import Item
    rng = np.random.default_rng(0)
    hu = lambda: rng.integers(-1000, 400, (100, 100), dtype=np.int16)
    out = []
    for s, split in enumerate(splits, 1):
        out.append(Item(f"riga/MESSIDOR/image00{s}", "riga", f"riga:MESSIDOR/image00{s}", split,
                        rng.integers(0, 256, (120, 160, 3), dtype=np.uint8),
                        {"disc": np.stack([_disc(120, 160, 60, 80, 30 + d) for d in (0, 1, 2)]),
                         "cup": np.stack([_disc(120, 160, 60, 80, 14 + d) for d in (0, 1, 2)])}, ["1", "2", "3"], {}))
        g, z = f"LIDC-IDRI-000{s}", f"lidc/1.2.{s}/n000/z00"
        nod = lambda *discs: {"nodule": np.stack([_disc(100, 100, *d) for d in discs])}
        out += [Item(z + "1", "lidc", g, split, hu(), nod((50, 50, 12), (50, 50, 12, 5)), ["1", "2"], {}),
                Item(z + "2", "lidc", g, split, None, nod((40, 40, 8), (40, 40, 8)), ["1", "2"], {}),
                Item(z + "3", "lidc", g, split, hu(), nod((30, 60, 6)), ["1"], {})]
        if unanimous:
            out.append(Item(z + "4", "lidc", g, split, hu(), nod((60, 40, 9), (60, 40, 10)), ["1", "2"], {}))
    return sorted(out, key=lambda i: (i.dataset, i.split, i.item_id))


@pytest.fixture(scope="session")
def disc_mask():
    """disc_mask(h, w, cy, cx, r, hole=0) -> bool (h, w): a disc, a ring when hole > 0."""
    return _disc


@pytest.fixture(scope="session")
def synthetic_items():
    """synthetic_items(splits=("train",), unanimous=False) -> [Item] (see _synthetic_items)."""
    return _synthetic_items


@pytest.fixture(scope="session")
def char_model(tmp_path_factory):
    """(checkpoint path, tokenizer): a 1-layer Llama (seed 0) with a character-level tokenizer over the answer
    grammar's characters, <pad> = 0 and <eos> = 1, saved offline for TRL tests."""
    torch = pytest.importorskip("torch")
    pytest.importorskip("transformers")
    from tokenizers import Regex, Tokenizer, decoders, models, pre_tokenizers
    from transformers import LlamaConfig, LlamaForCausalLM, PreTrainedTokenizerFast
    chars = sorted(set('{}[]":,0123456789abcdefghijklmnopqrstuvwxyz'))
    vocab = {w: i for i, w in enumerate(["<pad>", "<eos>", *chars])}
    tk = Tokenizer(models.WordLevel(vocab, unk_token="<pad>"))
    tk.pre_tokenizer = pre_tokenizers.Split(Regex("."), behavior="isolated")
    tk.decoder = decoders.Fuse()
    tok = PreTrainedTokenizerFast(tokenizer_object=tk, pad_token="<pad>", eos_token="<eos>", padding_side="left")
    torch.manual_seed(0)
    path = tmp_path_factory.mktemp("char") / "model"
    LlamaForCausalLM(LlamaConfig(vocab_size=len(tok), hidden_size=16, intermediate_size=32, num_hidden_layers=1,
                                 num_attention_heads=2, num_key_value_heads=2, pad_token_id=0, eos_token_id=1,
                                 bos_token_id=1, max_position_embeddings=4096)).save_pretrained(path)
    return path, tok


@pytest.fixture
def clean_store():
    """An empty objdel.train.serve STORE before and after the test, and gc unfrozen (fork_pool freezes it)."""
    import gc

    from objdel.train import serve

    serve.STORE.clear()
    yield
    serve.STORE.clear()
    gc.unfreeze()

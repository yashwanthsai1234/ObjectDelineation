"""E-G15 (E8.1): the Apache-2.0 text ships inside objdel.infra; the card fallback returns exactly those bytes."""

import hashlib

import pytest

from objdel.infra import hf

APACHE_SHA256 = "cfc7749b96f63bd31c3c42b5c471bf756814053e847c10f3eb003417bc523d30"


def test_base_card_falls_back_to_the_packaged_apache_licence(monkeypatch):
    hub = pytest.importorskip("huggingface_hub")
    from huggingface_hub.errors import EntryNotFoundError

    card = type("Card", (), {"data": type("D", (), {"to_dict": lambda self: {"license": "apache-2.0"}})()})()
    monkeypatch.setattr(hub.ModelCard, "load", classmethod(lambda cls, base, token=None: card))

    def missing(*a, **k):
        raise EntryNotFoundError("no LICENSE")

    monkeypatch.setattr(hub, "hf_hub_download", missing)
    lic, meta = hf.base_card("Qwen/Qwen3-VL-2B-Instruct")
    assert hashlib.sha256(lic).hexdigest() == APACHE_SHA256 and meta == {"license": "apache-2.0"}

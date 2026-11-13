"""Fakes for the Hugging Face push and pull tests."""

REPO = "svenk029/objdel-sft"


def token(role="fineGrained", scoped=None, global_=()):
    """A whoami() answer: a fine-grained token scoped to REPO by default."""
    scoped = [{"entity": {"type": "model", "name": REPO}, "permissions": ["repo.content.read", "repo.write"]}] \
        if scoped is None else scoped  # fmt: skip
    return {"auth": {"type": "access_token", "accessToken": {
        "role": role, "fineGrained": {"global": list(global_), "scoped": scoped}}}}  # fmt: skip


def fake_mac(monkeypatch, who) -> list:
    """Patch the Mac-side box commands (get_run, secret, remote) and hf.whoami to answer `who`; return the list of
    (cmd, stdin, what) the fake `remote` was sent."""
    from objdel.infra import hf as H
    from objdel.infra.vast import __main__ as P

    sent = []
    monkeypatch.setattr(P, "get_run", lambda rid: {"instance": 7, "ip": "1.2.3.4", "port": 22})
    monkeypatch.setattr(P, "secret", lambda s: {"objdel-hf-write": "hf_secret"}[s])
    monkeypatch.setattr(H, "whoami", lambda tok: who)
    monkeypatch.setattr(P, "remote", lambda run, cmd, stdin=None, what="": sent.append((cmd, stdin, what)))
    return sent

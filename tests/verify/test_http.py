"""The verifier HTTP app (E7): the same floats as in process, typed errors, an api_version in every response."""
import pytest
from conftest import KEY

pytest.importorskip("fastapi")
pytest.importorskip("httpx")

from fastapi.testclient import TestClient  # noqa: E402

from objdel.verify import http  # noqa: E402
from objdel.verify.api import API_VERSION, ScoreItem, Verifier, VerifyError  # noqa: E402


@pytest.fixture
def client(view_root):
    return TestClient(http.make_app(Verifier("reward_v4", view_root)))


def body(answers, n=None):
    n = n or len(answers)
    return {"reward": "reward_v4", "items": [{"id": f"i{i}", "key": KEY, "text": answers[i % len(answers)],
                                              "truncated": False} for i in range(n)]}


def test_score_returns_the_floats_of_the_process(client, view_root, answers):
    r = client.post("/v1/score", json=body(answers))
    want = Verifier("reward_v4", view_root).score([ScoreItem(f"i{i}", KEY, a, False) for i, a in enumerate(answers)])
    assert r.status_code == 200 and r.json()["api_version"] == API_VERSION
    assert r.json()["stamp"]["reward_name"] == "reward_v4"
    for got, w in zip(r.json()["results"], want, strict=True):
        assert got["id"] == w.id and got["error"] is None
        assert float(got["reward"]).hex() == w.reward.hex()
        assert {k: float(x).hex() for k, x in got["terms"].items()} == {k: x.hex() for k, x in w.terms.items()}


def test_batch_over_max_is_413(client, answers):
    r = client.post("/v1/score", json=body(answers, http.MAX_BATCH + 1))
    assert r.status_code == 413 and r.json()["api_version"] == API_VERSION
    assert r.json()["error"]["code"] == "batch_too_large" and r.json()["error"]["retryable"] is False


def test_unknown_key_is_404_and_wrong_reward_409(client, answers):
    b = body(answers)
    b["items"][0]["key"] = "lidc/nope"
    r = client.post("/v1/score", json=b)
    assert r.status_code == 404 and r.json()["error"]["code"] == "unknown_key"
    r = client.post("/v1/score", json=body(answers) | {"reward": "reward_v3"})
    assert r.status_code == 409 and r.json()["error"]["code"] == "spec_mismatch"


def test_device_error_is_the_only_retryable(client, monkeypatch, answers):
    def boom(self, items):
        raise VerifyError("device_error", "out of memory")
    monkeypatch.setattr(Verifier, "score", boom)
    r = client.post("/v1/score", json=body(answers))
    assert r.status_code == 503 and r.json()["error"]["retryable"] is True


def test_spec_and_probes(client):
    s = client.get("/v1/spec")
    assert s.status_code == 200 and s.json()["api_version"] == API_VERSION and s.json()["stamp"]["reward_fingerprint"]
    for path in ("/livez", "/readyz"):
        r = client.get(path)
        assert r.status_code == 200 and r.json()["api_version"] == API_VERSION

"""HTTP face of the verifier (extra `serve`; training does not use it): FastAPI over objdel.verify.api.

    python -I -m objdel.verify.http --view data/eval_views/<view> --port 8000

POST /v1/score {"reward", "items": [{"id", "key", "text", "truncated"}]} (at most MAX_BATCH items), GET /v1/spec,
GET /livez, GET /readyz. Every JSON body carries api_version. An error is {"api_version", "error": {"code",
"message", "retryable"}}; only device_error is retryable. pydantic models live only here.

MAX_BATCH: one RL step scores 32 prompts x 8 answers = 256 answers (ADR-0022); 512 leaves 2x headroom.

References: https://fastapi.tiangolo.com/tutorial/handling-errors/#install-custom-exception-handlers
"""
import argparse
import dataclasses
from pathlib import Path

import uvicorn
from fastapi import FastAPI, Request
from fastapi.responses import JSONResponse
from pydantic import BaseModel

from objdel.verify.api import API_VERSION, ScoreItem, Verifier, VerifyError

MAX_BATCH = 512
STATUS = {"unknown_key": 404, "spec_mismatch": 409, "batch_too_large": 413, "device_error": 503}


class Item(BaseModel):
    id: str
    key: str
    text: str
    truncated: bool


class ScoreRequest(BaseModel):
    reward: str
    items: list[Item]


def make_app(verifier: Verifier) -> FastAPI:
    app = FastAPI(title="objdel verifier", version=API_VERSION)

    @app.exception_handler(VerifyError)
    def refuse(_: Request, e: VerifyError) -> JSONResponse:
        error = {"code": e.code, "message": str(e), "retryable": e.code == "device_error"}
        return JSONResponse({"api_version": API_VERSION, "error": error}, STATUS[e.code])

    @app.post("/v1/score")
    def score(req: ScoreRequest) -> dict:
        if req.reward != verifier.stamp.reward_name:
            raise VerifyError("spec_mismatch", f"this service scores {verifier.stamp.reward_name}, not {req.reward}")
        if len(req.items) > MAX_BATCH:
            raise VerifyError("batch_too_large", f"{len(req.items)} items; the limit is {MAX_BATCH}")
        results = verifier.score([ScoreItem(**i.model_dump()) for i in req.items])
        return {"api_version": API_VERSION, "stamp": dataclasses.asdict(verifier.stamp),
                "results": [dataclasses.asdict(r) for r in results]}

    @app.get("/v1/spec")
    def spec() -> dict:
        return {"api_version": API_VERSION, "max_batch": MAX_BATCH} | verifier.spec()

    @app.get("/livez")
    def livez() -> dict:
        return {"api_version": API_VERSION, "status": "ok"}

    @app.get("/readyz")  # the app exists only after the Verifier loaded its caches
    def readyz() -> dict:
        return {"api_version": API_VERSION, "status": "ready"}

    return app


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--view", type=Path, required=True, help="prepared view root (data/prepared or an eval view)")
    p.add_argument("--port", type=int, default=8000)
    a = p.parse_args()
    uvicorn.run(make_app(Verifier("reward_v4", a.view)), host="127.0.0.1", port=a.port)


if __name__ == "__main__":
    main()

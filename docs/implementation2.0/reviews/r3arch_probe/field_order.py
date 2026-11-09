from typing import Literal, Mapping
from pydantic import BaseModel, Field, JsonValue
class StageBase(BaseModel, frozen=True, extra="forbid"):
    kind: str
    name: str
    soft: bool = False
class TrainerFlags(BaseModel, frozen=True, extra="forbid"):
    lr: float | None = None
    shared_prefix: bool = False
    grpo: Mapping[str, JsonValue] = {}
class RlStage(StageBase, TrainerFlags, frozen=True, extra="forbid"):
    kind: Literal["rl"] = "rl"
    model: str
    ckpt: str
    steps: int = Field(gt=0)
print(list(RlStage.model_fields))
s = RlStage(name="a", model="m", ckpt="c", steps=1)
try:
    hash(s); print("hash ok")
except TypeError as e:
    print("hash TypeError", e)

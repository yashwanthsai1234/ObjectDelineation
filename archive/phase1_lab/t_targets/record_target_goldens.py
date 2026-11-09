"""Record tests/goldens/targets_lab.npz: what the lab target builder makes of the fixture crops (CLEAN4).

    .venv/bin/python -I experiments/phase1/lab/t_targets/record_target_goldens.py data/fixtures tests/goldens/targets_lab.npz

Only 1-D or 0-D values are stored, never a mask (RIGA terms forbid copying it). Per LIDC case `lidc/<item>/<rater>`:
the lab vertex count, the SHA-256 of the lab text and the SHA-256 of the lab-cropped masks. For RIGA only three
totals: `riga_vertices`, `riga_cases`, `riga_sha` (one hash over all lab-cropped masks in case order).
"""

import hashlib
import sys
from pathlib import Path

import numpy as np

sys.path[:0] = [str(Path(__file__).parent), str(Path(__file__).parents[4])]
from targets.builder import RIGA_DERIVED, BuildConfig, build_target  # noqa: E402
from targets.crops import apply_crop  # noqa: E402

from objdel.geometry import canonicalize, serialize  # noqa: E402
from tests.oracles.crops_oracle import mask_sha  # noqa: E402

CASE = np.dtype([("vertices", "i4"), ("text_sha", "U64"), ("crop_sha", "U64")])


def record(fixtures: Path) -> dict:
    out, riga_v, riga_n, riga_h = {}, 0, 0, hashlib.sha256()
    for ds, step in (("riga", 20), ("lidc", 6)):
        z = np.load(fixtures / f"{ds}_experts.npz")
        names = ("disc", "cup") if ds == "riga" else ("nodule",)
        for it in sorted({k.rsplit("__", 1)[0] for k in z.files})[::step]:
            st = {n: z[f"{it}__{n}"].astype(bool) for n in names}
            for r in range(st[names[0]].shape[0]):
                m = apply_crop({n: st[n][r] for n in names}, "pad_div")[0]
                _, lo, _ = build_target(m, {n: "parts" for n in names}, BuildConfig(), RIGA_DERIVED if ds == "riga" else ())
                if ds == "lidc":
                    text = serialize(canonicalize({n: lo[n].parts for n in names}), names)
                    sha = hashlib.sha256(text.encode()).hexdigest()
                    out[f"lidc/{it}/{r}"] = np.array((lo["nodule"].vertices, sha, mask_sha(m)), CASE)
                else:
                    riga_v += sum(lo[n].vertices for n in names)
                    riga_n += 1
                    riga_h.update(mask_sha(m).encode())
    return out | {"riga_vertices": np.array(riga_v), "riga_cases": np.array(riga_n), "riga_sha": np.array(riga_h.hexdigest())}


if __name__ == "__main__":
    np.savez_compressed(sys.argv[2], **record(Path(sys.argv[1])))

"""The dataset registry: every fact the shared code needs about a dataset, in one entry per dataset.

Adding a dataset = one `Dataset` entry in DATASETS plus its ingest code (an `objdel.data.store.Item` per image,
written by an offline build script). The shared layers read the entry at call time through `get` (structure names,
pixel contract, window, frame policy, ops, prompt, metric and gate thresholds), so nothing else names a dataset.
Plain data, below every layer that reads it (no numpy here).
"""
from dataclasses import dataclass


@dataclass(frozen=True)
class Structure:
    """One outlined structure: its answer key and store mask name, and its metric and gate thresholds."""
    name: str
    nsd_tau: float  # NSD@tau of verify.metrics, px at frame resolution
    biou_d: int  # boundary-IoU band width, px at frame resolution
    gate_median_iou: float | None = None  # G-geometry median IoU_tc; None: the gate cannot pass


@dataclass(frozen=True)
class FrameConfig:
    """Frame policy: largest and smallest side and margins (pixels, percent) and the model image side."""
    max_side: int = 500  # a divisor of 1000
    min_side: int = 50
    margin_px: int = 8
    margin_pct: int = 25
    image_side: int | None = None  # upscale the model image to this side (LIDC), else the frame side


@dataclass(frozen=True)
class Dataset:
    """A dataset's facts. `prompt` is verbatim: changing it changes the model's input."""
    name: str  # item id prefix and folder name
    structures: tuple[Structure, ...]  # in answer order
    prompt: str
    pixels: tuple[str, int]  # the store's image dtype and ndim: ("int16", 2) | ("uint8", 3)
    window: tuple[int, int] | None = None  # intensity window to 0..255; None: the images are uint8 already
    window_jitter: tuple[tuple[int, int], tuple[int, int]] | None = None  # (level, width) centre and half-range
    noise_max: float | None = None  # Gaussian noise sigma upper bound of the training views, image units
    frame: FrameConfig = FrameConfig()
    ops: tuple[str, ...] = ("id",)  # the training views' dihedral ops (objdel.data.targets.augment.OPS)
    rotated_op: str | None = None  # the op of the rotated eval view; None: no rotated view
    min_raters: int = 2  # raters per RL item; RL refuses an entry below the reward's 2 (no single-rater reward)
    exams: tuple[str, ...] = ()  # the eval views (data/eval_views/<name>) a bundle carries and preflight checks
    checker_eval: bool = False  # data/checker_eval (the CUDA checker's archived references) rides in its bundle
    restricted: bool = False  # its licence forbids uploads (no rows, no model trained on it, on the Hub)

    @property
    def names(self) -> tuple[str, ...]:
        """The structure names in answer order (the parse schema)."""
        return tuple(s.name for s in self.structures)


LIDC = Dataset(
    name="lidc",
    structures=(Structure("nodule", nsd_tau=2.0, biou_d=2, gate_median_iou=0.60),),
    prompt='Outline the lung nodule in this CT crop. Answer with JSON only: '
           '{"nodule":[{"outer":[[x,y],...],"holes":[[[x,y],...]]}]}, '
           'integers 0-1000 relative to the crop, outer rings clockwise, holes counter-clockwise.',
    pixels=("int16", 2),  # Hounsfield units
    window=(-1000, 400),  # the lung window
    window_jitter=((-300, 60), (1400, 80)),  # CT level and width, HU
    noise_max=10.0,  # HU
    frame=FrameConfig(image_side=512),
    ops=("id", "r90", "r180", "r270", "fh", "fv", "tr", "at"),
    rotated_op="r90",
    exams=("lidc-val-shifted", "lidc-val-rescaled", "lidc-val-rotated", "probe-precision", "probe-topology"),
    checker_eval=True,
)

DATASETS: dict[str, Dataset] = {d.name: d for d in (LIDC,)}


def get(name: str) -> Dataset:
    """The entry of `name`; KeyError naming the registered datasets."""
    if name not in DATASETS:
        raise KeyError(f"unknown dataset {name!r}; registered: {sorted(DATASETS)} (objdel.core.datasets)")
    return DATASETS[name]


def structure(name: str) -> Structure:
    """The structure `name` of any registered dataset; KeyError when none or two different entries have it."""
    found = {s for d in DATASETS.values() for s in d.structures if s.name == name}
    if len(found) != 1:
        raise KeyError(f"structure {name!r}: {len(found)} different entries in the registered datasets")
    return found.pop()

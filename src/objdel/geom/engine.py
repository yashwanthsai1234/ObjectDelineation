"""The CUDA geometry engine: the only importer of the `objdel._gpucheck` extension (E3a).

Callers import it as `from objdel.geom import engine`. A test checks that `__all__` covers every
`_gpucheck.<name>` used in the repo.
"""

from objdel._gpucheck import (  # noqa: TID251
    DEFECTS,
    add_tables,
    add_tables_from_raters,
    add_tables_many,
    build,
    certify_polygons,
    check_masks,
    check_polygons,
    device_count,
    drop_tables,
    edt,
    eval_sets,
    reset_stats,
    resident,
    score_masks,
    score_polygons,
    stats,
    test_hooks,
    topology,
    trim,
)

__all__ = [
    "DEFECTS",
    "add_tables",
    "add_tables_from_raters",
    "add_tables_many",
    "build",
    "certify_polygons",
    "check_masks",
    "check_polygons",
    "device_count",
    "drop_tables",
    "edt",
    "eval_sets",
    "reset_stats",
    "resident",
    "score_masks",
    "score_polygons",
    "stats",
    "test_hooks",
    "topology",
    "trim",
]

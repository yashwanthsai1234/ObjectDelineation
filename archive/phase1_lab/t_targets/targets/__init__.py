"""SFT target builder: expert mask -> compact, legal, canonical polygon JSON (see REPORT.md).

    from targets import build_target, BuildConfig, RIGA_SCHEMA, RIGA_DERIVED, LIDC_SCHEMA
    text, objs, info = build_target({"disc": d, "cup": c}, RIGA_SCHEMA, BuildConfig(), RIGA_DERIVED)
"""
from .builder import (LIDC_SCHEMA, RIGA_DERIVED, RIGA_SCHEMA, BuildConfig, build_item, build_target, canonical_parts,
                      is_canonical, parse, rasterize_parts, serialize)
from .crops import apply_crop

__all__ = ["BuildConfig", "LIDC_SCHEMA", "RIGA_DERIVED", "RIGA_SCHEMA", "apply_crop", "build_item", "build_target",
           "canonical_parts", "is_canonical", "parse", "rasterize_parts", "serialize"]

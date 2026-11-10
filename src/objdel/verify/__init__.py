"""Verify layer: reward, metrics, statistics, gates and the doctor-store loader. No torch, no TRL, no data layer.

The versioned batch service (plan E7) sits on top: api (in-process, no web or torch stack) and http (FastAPI, extra
`serve`). Training scores through verify.reward directly and does not use the service."""

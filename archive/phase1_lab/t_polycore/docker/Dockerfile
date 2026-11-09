# Linux portability check for polycore (build from source + full test suite).
FROM python:3.12-slim
RUN apt-get update && apt-get install -y --no-install-recommends build-essential cmake \
    && rm -rf /var/lib/apt/lists/*
RUN pip install --no-cache-dir numpy pytest hypothesis shapely scikit-image numba nanobind scikit-build-core \
    && (pip install --no-cache-dir rasterio || echo "rasterio unavailable: its test will be skipped")
WORKDIR /w

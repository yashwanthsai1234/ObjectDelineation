# x86-64 (linux/amd64) test image for objdel (Phase 1, M9). Built and run by tools/xarch/x86.sh.
#
# Only the dependencies live in the image; the repo source is streamed in at run time (x86.sh), so a re-run on new
# code reuses this layer. The build context is a generated directory holding only requirements.txt, the exact
# package pins of the Mac's .venv (uv pip freeze), so the container runs the same numpy / scipy / OpenCV / numba /
# torch versions as the Mac: reward-cache fingerprints and hashes are only comparable at equal versions (plan §6).
#
# On Apple silicon, Docker Desktop runs linux/amd64 images under Rosetta 2 (or QEMU) emulation; correctness and
# determinism are meaningful there, timing is not.
#
# References:
#   Docker multi-platform builds: https://docs.docker.com/build/building/multi-platform/
#   Docker Desktop Rosetta emulation: https://docs.docker.com/desktop/settings-and-maintenance/settings/#general
#   Dockerfile --platform / pinning by digest: https://docs.docker.com/reference/dockerfile/#from
#   uv in Docker: https://docs.astral.sh/uv/guides/integration/docker/
#   uv + CPU-only PyTorch (--torch-backend): https://docs.astral.sh/uv/guides/integration/pytorch/
#   scikit-build-core (no build isolation needs the backend in the env):
#     https://scikit-build-core.readthedocs.io/en/latest/configuration/index.html
#   auditwheel: https://github.com/pypa/auditwheel ; abi3audit: https://github.com/pypa/abi3audit

# python:3.12-slim-trixie, multi-arch index digest (resolved 2026-10-19), x86.sh builds and runs it with --platform linux/amd64
FROM python:3.12-slim-trixie@sha256:a6e34c598f2467ed0e9a8d349809fcd8b5c603269512df273a0bb1784edc11b1

# g++ and cmake build objdel._gpucheck from source (CMakeLists.txt: Release, -O3 -DNDEBUG, nanobind NOMINSIZE);
# patchelf for auditwheel; git is not needed (the source arrives as a tarball).
RUN apt-get update && apt-get install -y --no-install-recommends build-essential cmake ninja-build patchelf \
    && rm -rf /var/lib/apt/lists/*

COPY --from=ghcr.io/astral-sh/uv:0.12.23 /uv /usr/local/bin/uv
ENV UV_LINK_MODE=copy UV_NO_CACHE=1 PYTHONDONTWRITEBYTECODE=1

# The Mac's pins (core + dev + train + data + models extras; scikit-build-core and nanobind included, so the
# extension builds with --no-build-isolation and the run needs no network). torch / torchvision from the CPU index.
COPY requirements.txt /opt/requirements.txt
RUN uv venv --python /usr/local/bin/python3.12 /opt/venv \
    && uv pip install --python /opt/venv/bin/python --torch-backend cpu -r /opt/requirements.txt

# Wheel tools in their own env, so they cannot change a pinned package of /opt/venv.
RUN uv venv --python /usr/local/bin/python3.12 /opt/tools \
    && uv pip install --python /opt/tools/bin/python auditwheel==6.8.2 abi3audit==0.0.26

ENV PATH=/opt/venv/bin:$PATH
WORKDIR /work

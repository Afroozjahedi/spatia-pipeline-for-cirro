# ============================================================
# Dockerfile — SPATIA pipeline (CRC TMA / Nature Protocols / Cirro)
#
# Python 3.9, updated 2026-07-21: matches the one confirmed real environment
# spacec has actually run in for Afrouz — BTC's JupyterHub. Everything else
# in this file (llvmlite/numba/conda-forge ordering) is still unverified
# against that environment (see trailing comment) — Python version is the
# one piece we now have real confidence in, not a guess.
#
# Written fresh on 2026-07-21 (Q12). NOTE: there is an older Dockerfile in
# OneDrive-InsideMDAnderson/BTC/SPATIA_pipeline/ — Afrouz explicitly asked
# NOT to reuse it ("ignore docker file. create one heaere"), so this file
# does not reference or copy from it. Only this file is the source of truth.
#
# Packaging target: single black-box container (Q10 — one Nextflow process
# calls run_pipeline.py end-to-end inside this image; no per-step containers).
#
# Dependency ordering below mirrors setup_local_env.sh exactly, because
# Day 2's research (SPATIA_PIPELINE_LOG.md) found that `pip install spacec`
# alone fails on x86_64: llvmlite (pulled in transitively via numba) tries
# to build from source and fails. The fix is installing llvmlite/numba from
# conda-forge (prebuilt binaries) BEFORE pip ever touches spacec. This is
# the single most important ordering constraint in this file — do not
# collapse the RUN steps below into one big `pip install -r requirements`
# without preserving conda-forge-first.
# ============================================================

FROM continuumio/miniconda3:latest

WORKDIR /app

# ---- system-level deps that setup_local_env.sh installs via conda ----
# (graphviz/libvips/openslide are native libraries, not pip-installable)
RUN conda create -y -n spatia python=3.9 graphviz libvips openslide -c conda-forge \
    && conda clean -afy

SHELL ["conda", "run", "-n", "spatia", "/bin/bash", "-c"]

# ---- llvmlite/numba from conda-forge FIRST (prebuilt binaries) ----
# This must happen before spacec is installed via pip, or the pip install
# will attempt to build llvmlite from source and fail on x86_64.
RUN conda install -y -c conda-forge llvmlite numba \
    && conda clean -afy

# ---- spacec + pipeline dependencies ----
# spacec pinned to 0.0.10 — confirmed via `pip show spacec` on Afrouz's
# real working environment (BTC JupyterHub, Python 3.9). That output shows
# deepcell as a normal transitive dependency of spacec 0.0.10 (under
# "Requires"), so plain `pip install spacec==0.0.10` should pull it in
# automatically — no separate/manual deepcell install needed. (Earlier
# version of this file excluded deepcell on purpose, based on
# setup_local_env.sh; Afrouz has since said that script was never actually
# validated, so that guidance is corrected here now that a real source
# exists.)
RUN pip install --no-cache-dir spacec==0.0.10 lifelines kneed seaborn pyyaml

# ---- scanpy (added 2026-07-29, Q24) ----
# spatia.analysis.cell_typing.run_cell_typing() raises ImportError immediately
# if scanpy isn't importable -- this is true for BOTH "automatic" and
# "semi_automatic" mode (see cell_typing.py line ~371: `if not HAS_SCANPY:
# raise ImportError(...)`), not just semi_automatic as the try/except import
# comment there might suggest at a glance. Confirmed by reading the code
# directly, not assumed. spacec's own dependency tree was not confirmed to
# include scanpy (never actually installed end-to-end in this sandbox), so
# this is pinned explicitly rather than relied on as a transitive dependency.
RUN pip install --no-cache-dir scanpy

# ---- SPATIA package itself ----
COPY pyproject.toml ./
COPY spatia/ ./spatia/
COPY run_pipeline.py ./
RUN pip install --no-cache-dir -e . --no-deps

# ---- experiment configs + cell-type definitions ----
# (added cell_type_definitions/ 2026-07-29, Q24 -- it was missing from this
# COPY list even though experiments/*.yaml already reference files under it,
# e.g. crc_tma_celltyping.yaml -> cell_type_definitions/crc_tma.yaml. Any
# config using cell_type_definitions_file would have failed inside the
# container with a FileNotFoundError until this was added.)
COPY experiments/ ./experiments/
COPY cell_type_definitions/ ./cell_type_definitions/

# Make the conda env the default python on PATH for the container's
# entrypoint, so `docker run <image> --config ...` works without needing
# callers to know about `conda run`.
ENV PATH=/opt/conda/envs/spatia/bin:$PATH

ENTRYPOINT ["python", "run_pipeline.py"]

# ============================================================
# Build/verification status (2026-07-21, updated 2026-07-29): written and
# reviewed against setup_local_env.sh's documented ordering, but NOT built
# or run in this sandbox — `pip install spacec` alone pulls ~900MB of torch
# as a transitive dependency, and the sandbox used for this work ran out of
# disk space attempting exactly this install on 2026-07-20 (see
# SPATIA_PIPELINE_LOG.md Day 4, Part C). A real build/test needs to happen
# on a machine with sufficient disk (Afrouz's machine, Seadragon, or an
# actual Cirro/AWS Batch build step) before this is trusted as correct —
# flagging honestly rather than claiming it's verified.
#
# 2026-07-29 (Q24): Afrouz confirmed she already has a real spacec
# container (ghcr.io/break-through-cancer/btc-spatial-proteomics/spacec)
# rather than needing this from-scratch conda/pip build. See
# Dockerfile.spacec-base in this same directory for an alternate version
# that builds FROM that image instead of continuumio/miniconda3 — lower
# risk if that container's spacec/deepcell install already works, but
# unverified from here: this sandbox cannot pull ghcr.io images (network
# allowlist blocks it, confirmed directly) to inspect the base image's
# actual Python path / conda env name. Treat Dockerfile.spacec-base as a
# draft to adjust once you've inspected the real container, not as tested.
# ============================================================

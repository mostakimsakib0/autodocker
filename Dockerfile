# ============================================================
# Stage 0 — Build Image Base
# ============================================================
FROM docker.io/debian:bookworm-slim AS base_builder

COPY tools/builder/apt.sh /scripts/
RUN /scripts/apt.sh --init builder
# copy other tools
COPY --exclude=apt.sh tools/builder/* /scripts/
WORKDIR /src

# ============================================================
# Stage 1 — Build Rootless Mode ACL Manager
# ============================================================

FROM base_builder AS aclmgr_builder
COPY tools/aclmgr ./aclmgr
RUN /scripts/buildpkg.sh aclmgr

# ============================================================
# Stage 1 — Build AutoDock Vina
# ============================================================
FROM base_builder AS vina_builder
COPY tools/autodock-vina/deps ./autodock-vina/deps
RUN /scripts/apt.sh dev autodock-vina

COPY --exclude=deps tools/autodock-vina ./autodock-vina
RUN /scripts/buildpkg.sh autodock-vina

# ============================================================
# Stage 1 — Build SMINA
# ============================================================
FROM base_builder AS smina_builder
COPY tools/smina/deps ./smina/deps
RUN /scripts/apt.sh dev smina

COPY --exclude=deps tools/smina ./smina
RUN /scripts/buildpkg.sh smina
RUN /scripts/shlib.sh boost smina/build/fromsmina smina/build/tosmina smina/build/smina > deps

# ============================================================
# Stage 1 — Build QuickVina2
# ============================================================
FROM base_builder AS qvina_builder
COPY tools/qvina/deps ./qvina/deps
RUN /scripts/apt.sh dev qvina

COPY --exclude=deps tools/qvina ./qvina
RUN /scripts/buildpkg.sh qvina

# ============================================================
# Stage 1 — Build NGL
# ============================================================
# Pinned to a fixed pnpm release (not :latest) so builds are reproducible.
FROM ghcr.io/pnpm/pnpm:12.8.0 AS ngl_builder
COPY --from=base_builder /scripts /scripts
WORKDIR /src
COPY tools/ngl ./ngl
RUN /scripts/buildpkg.sh ngl

# ============================================================
# Stage 1 — Build FPocket
# ============================================================

FROM base_builder AS fpocket_builder
COPY tools/fpocket/deps ./fpocket/deps
RUN /scripts/apt.sh dev fpocket

COPY --exclude=deps tools/fpocket ./fpocket
RUN /scripts/patch.sh fpocket
RUN make -j"$(nproc)" -C fpocket/repo CXX=g++
RUN strip fpocket/repo/bin/*

# ============================================================
# Stage 1.5 — Accumulate deps
# ============================================================
FROM docker.io/debian:bookworm-slim as deps

WORKDIR /acc
COPY deps base
COPY --from=smina_builder /src/deps smina
RUN cat * > /deps
# ============================================================
# Stage 2 — Minimal runtime image
# ============================================================
FROM docker.io/python:3.12-slim-bookworm

LABEL org.opencontainers.image.title="AutoDocker" \
      org.opencontainers.image.description="Containerized virtual-screening pipeline (protein/ligand prep, Vina docking, ranked reports)." \
      org.opencontainers.image.url="https://github.com/mostakimsakib0/autodocker" \
      org.opencontainers.image.source="https://github.com/mostakimsakib0/autodocker" \
      org.opencontainers.image.documentation="https://github.com/mostakimsakib0/autodocker#readme" \
      org.opencontainers.image.vendor="Mostakim Sakib" \
      org.opencontainers.image.licenses="MIT" \
      org.opencontainers.image.version="1.1.4"

COPY --from=deps /deps /deps
COPY tools/builder/apt.sh /apt.sh
RUN /apt.sh --init run .

COPY requirements.txt /requirements.txt
RUN pip install --no-cache-dir --root-user-action=ignore -r /requirements.txt

COPY autodocker /autodocker

COPY --from=vina_builder  /pkg/ /
COPY --from=smina_builder /pkg/ /
COPY --from=qvina_builder /pkg/ /

COPY --from=fpocket_builder \
    /src/fpocket/repo/bin/fpocket \
    /src/fpocket/repo/bin/dpocket \
    /src/fpocket/repo/bin/tpocket \
    /usr/local/bin/
COPY --from=fpocket_builder \
    /src/fpocket/repo/plugins/LINUXAMD64/molfile/ \
    /usr/local/lib/

COPY --from=ngl_builder /pkg/ /

COPY --from=aclmgr_builder /pkg/ /

RUN useradd --create-home --uid 1000 dockuser
RUN mkdir -p /workspace && chown -R dockuser:dockuser /workspace /autodocker

WORKDIR /workspace
VOLUME ["/workspace"]

ENTRYPOINT ["/autodocker/aclmgr"]

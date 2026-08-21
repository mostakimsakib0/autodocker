# ============================================================
# Stage 0 — Build Image Base
# ============================================================
FROM docker.io/debian:bookworm-slim AS builder

ENV DEBIAN_FRONTEND=noninteractive
COPY scripts/apt.sh /scripts/apt.sh
COPY scripts/patch.sh /scripts/patch.sh

RUN <<EOF
set -eu
apt-get update
apt-get install \
	-y --no-install-recommends \
build-essential git
EOF

# ============================================================
# Stage 1 — Build AutoDock Vina
# ============================================================
FROM builder AS vina_builder

WORKDIR /src
COPY tools/autodock-vina ./autodock-vina
RUN /scripts/apt.sh dev autodock-vina/deps
RUN make -j"$(nproc)" -C autodock-vina/repo/build/linux/release
RUN strip autodock-vina/repo/build/linux/release/vina

# ============================================================
# Stage 1 — Build SMINA
# ============================================================
FROM builder AS smina_builder

WORKDIR /src
COPY tools/smina ./smina
RUN /scripts/apt.sh dev smina/deps
RUN cmake -S smina/repo -B smina/build
RUN make -j"$(nproc)" -C smina/build
RUN strip smina/build/fromsmina smina/build/tosmina smina/build/smina

# ============================================================
# Stage 1 — Build QuickVina2
# ============================================================
FROM builder AS qvina_builder
WORKDIR /src
COPY tools/qvina ./qvina
RUN /scripts/apt.sh dev qvina/deps
RUN make -j"$(nproc)" -C qvina/repo/build/linux/release
RUN <<EOF
strip qvina/repo/build/linux/release/vina
mv qvina/repo/build/linux/release/vina qvina/repo/build/linux/release/qvina
EOF

# ============================================================
# Stage 1 — Build NGL
# ============================================================
FROM ghcr.io/pnpm/pnpm:latest AS ngl_builder
WORKDIR /src
COPY tools/ngl .
RUN pnpm install
RUN mv node_modules/ngl/dist ngl

# ============================================================
# Stage 1 — Build FPocket
# ============================================================

FROM builder AS fpocket_builder
WORKDIR /src
COPY tools/fpocket ./fpocket
RUN /scripts/patch.sh fpocket
RUN /scripts/apt.sh dev fpocket/deps
RUN make -j"$(nproc)" -C fpocket/repo CXX=g++
RUN strip fpocket/repo/bin/*

# ============================================================
# Stage 1 — Build Rootless Mode ACL Manager
# ============================================================

FROM builder AS aclmgr_builder
WORKDIR /src
COPY tools/aclmgr.c .
RUN gcc aclmgr.c -o aclmgr
RUN strip aclmgr

# ============================================================
# Stage 2 — Minimal runtime image
# ============================================================
FROM docker.io/python:3.12-slim-bookworm

COPY deps /deps
COPY scripts/apt.sh /apt.sh
RUN <<EOF
export DEBIAN_FRONTEND=noninteractive
apt update
/apt.sh run /deps
rm -f /apt.sh
EOF

COPY requirements.txt /requirements.txt
RUN pip install --no-cache-dir --root-user-action=ignore -r /requirements.txt

COPY autodocker /autodocker

COPY --from=vina_builder \
    /src/autodock-vina/repo/build/linux/release/vina \
    /usr/local/bin/vina


COPY --from=smina_builder \
    /src/smina/build/smina \
    /src/smina/build/tosmina \
    /src/smina/build/fromsmina \
    /usr/local/bin/

COPY --from=qvina_builder \
    /src/qvina/repo/build/linux/release/qvina \
    /usr/local/bin/qvina

COPY --from=fpocket_builder \
    /src/fpocket/repo/bin/fpocket \
    /src/fpocket/repo/bin/dpocket \
    /src/fpocket/repo/bin/tpocket \
    /usr/local/bin/
COPY --from=fpocket_builder \
    /src/fpocket/repo/plugins/LINUXAMD64/molfile/ \
    /usr/local/lib/

COPY --from=ngl_builder \
    src/ngl /ngl

COPY --from=aclmgr_builder \
	src/aclmgr /autodocker/aclmgr

RUN useradd --create-home --uid 1000 dockuser
RUN mkdir -p /workspace && chown -R dockuser:dockuser /workspace /autodocker

WORKDIR /workspace
VOLUME ["/workspace"]

ENTRYPOINT ["/autodocker/aclmgr"]

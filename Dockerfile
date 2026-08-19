# ============================================================
# Stage 0 — Build Image Base
# ============================================================
FROM docker.io/debian:bookworm-slim AS builder

ENV DEBIAN_FRONTEND=noninteractive
COPY scripts/apt.sh /scripts/apt.sh
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
# Stage 1 — Build OpenBabel
# ============================================================
FROM builder AS obabel_builder
WORKDIR /src
COPY tools/openbabel ./openbabel
RUN /scripts/apt.sh dev openbabel/deps
RUN cmake -S openbabel/repo -B openbabel/build \
    -DCMAKE_BUILD_TYPE=Release \
    -DBUILD_GUI=OFF \
    -DBUILD_SHARED=OFF \
    -DWITH_COORDGEN=OFF \
    -DWITH_MAEPARSER=OFF \
    -DWITH_AVALON=OFF \
    -DWITH_SWIG=OFF \
    -DWITH_PYTHON=OFF \
    -DWITH_OPENMP=ON \
    -DCMAKE_INSTALL_PREFIX=/usr/local
RUN make -j"$(nproc)" -C openbabel/build
RUN strip openbabel/build/bin/obabel
RUN <<EOF
make -C openbabel/build DESTDIR=/src/openbabel/pfx install
mkdir /src/openbabel/pfx/share
mv /src/openbabel/pfx/usr/local/bin /src/openbabel/pfx/bin
mv /src/openbabel/pfx/usr/local/share/openbabel /src/openbabel/pfx/share/openbabel
rm -rf /src/openbabel/pfx/usr
EOF
# ============================================================
# Stage 1 — Build FPocket
# ============================================================

FROM builder AS fpocket_builder
WORKDIR /src
COPY tools/fpocket ./fpocket
RUN /scripts/apt.sh dev fpocket/deps
RUN make -j"$(nproc)" -C fpocket/repo
RUN strip fpocket/repo/bin/*

# ============================================================
# Stage 2 — Minimal runtime image
# ============================================================
FROM docker.io/python:3.12-slim-bookworm

COPY requirements.txt /requirements.txt

RUN pip install --no-cache-dir --root-user-action=ignore -r /requirements.txt

COPY autodocker /autodocker

COPY --from=vina_builder \
    /src/autodock-vina/repo/build/linux/release/vina \
    /usr/local/bin/vina

COPY --from=qvina_builder \
    /src/qvina/repo/build/linux/release/qvina \
    /usr/local/bin/qvina

COPY --from=obabel_builder \
    /src/openbabel/pfx/ \
    /usr/local/

COPY --from=fpocket_builder \
    /src/fpocket/repo/bin/fpocket \
    /src/fpocket/repo/bin/dpocket \
    /src/fpocket/repo/bin/tpocket \
    /usr/local/bin/

# Run as a non-root user for safety
RUN useradd --create-home --uid 1000 dockuser
RUN mkdir -p /workspace && chown -R dockuser:dockuser /workspace /autodocker

USER dockuser
WORKDIR /workspace
VOLUME ["/workspace"]

# Smallest practical default
ENTRYPOINT ["/autodocker/entry.sh"]

# Reproducibility pins (override with --build-arg only on purpose):
#  - Base images are pinned by multi-arch index digest (each index covers
#    linux/amd64 and linux/arm64); the tag is kept for readability.
#
# Architectures: linux/amd64 and linux/arm64, each built natively (see
# .github/workflows/publish.yml). Floating-point code is compiled with
# -ffp-contract=off (AD_FPFLAGS) on both, so the compiler never fuses
# multiply-adds; GCC does so by default on arm64 but cannot on baseline
# x86-64, and the fusion changes rounding and hence scores.
#  - Debian packages (toolchain, Boost, Open Babel, runtime libraries) come
#    from snapshot.debian.org at a fixed timestamp, so every build resolves the
#    exact same package versions. Old snapshots are never removed, unlike
#    superseded versions on deb.debian.org, so version pins cannot rot.
#    Set DEBIAN_SNAPSHOT="" to use the live mirrors instead.
#  - NGL is installed from a committed pnpm lockfile (--frozen-lockfile).
#  - Python packages are pinned exactly in requirements.txt.
#  - No "# syntax=" directive on purpose: the builder's built-in Dockerfile
#    frontend (heredoc support) is used instead of pulling a floating one.
#
# Version metadata (single source of truth: autodocker/VERSION):
#   docker build \
#     --build-arg AUTODOCKER_VERSION="$(cat autodocker/VERSION)" \
#     --build-arg VCS_REF="$(git rev-parse HEAD)" -t autodocker .
ARG DEBIAN_IMAGE=docker.io/library/debian:bookworm-slim@sha256:abd67ffcfa541b485a3dff59865ab629aa048a6c613e639d36e7456b0b229241
ARG PYTHON_IMAGE=docker.io/library/python:3.12-slim-bookworm@sha256:a116514e19457bcb7af7efe9c3dd0b9b71e85b317694e7882a1c52aa15a78134
ARG PNPM_IMAGE=ghcr.io/pnpm/pnpm:12.8.0@sha256:e6247edbaab71cccb0a81ca2b74d3a576efd0200cd61d2f1927de3e62d98bdde
ARG DEBIAN_SNAPSHOT=20261010T000000Z

# ============================================================
# Stage 0 — Build Image Base
# ============================================================
FROM ${DEBIAN_IMAGE} AS builder

ARG DEBIAN_SNAPSHOT
ENV DEBIAN_FRONTEND=noninteractive
# Passed explicitly to every native build below (see header).
ENV AD_FPFLAGS="-ffp-contract=off"
COPY scripts/apt.sh /scripts/apt.sh
COPY scripts/shlib.sh /scripts/shlib.sh
COPY scripts/patch.sh /scripts/patch.sh

RUN <<EOF
set -eu
if [ -n "${DEBIAN_SNAPSHOT}" ]; then
	sed -i -E "s#^URIs: http://deb\.debian\.org/(debian|debian-security)\$#URIs: http://snapshot.debian.org/archive/\1/${DEBIAN_SNAPSHOT}#" \
		/etc/apt/sources.list.d/debian.sources
	test "$(grep -c "^URIs: http://snapshot.debian.org/archive/.*/${DEBIAN_SNAPSHOT}\$" /etc/apt/sources.list.d/debian.sources)" -eq 2
	# Snapshot Release files keep their original Valid-Until.
	echo 'Acquire::Check-Valid-Until "false";' > /etc/apt/apt.conf.d/99snapshot
fi
apt-get update -o Acquire::Retries=5
apt-get install \
	-y --no-install-recommends \
	-o Acquire::Retries=5 \
	build-essential git
EOF

# ============================================================
# Stage 1 — Build AutoDock Vina
# ============================================================
FROM builder AS vina_builder

WORKDIR /src
COPY tools/autodock-vina ./autodock-vina
RUN /scripts/apt.sh dev autodock-vina/deps
RUN make -j"$(nproc)" -C autodock-vina/repo/build/linux/release CFLAGS="$AD_FPFLAGS"
RUN strip autodock-vina/repo/build/linux/release/vina

# ============================================================
# Stage 1 — Build SMINA
# ============================================================
FROM builder AS smina_builder

WORKDIR /src
COPY tools/smina ./smina
RUN /scripts/apt.sh dev smina/deps
RUN cmake -S smina/repo -B smina/build \
	-DCMAKE_C_FLAGS="$AD_FPFLAGS" -DCMAKE_CXX_FLAGS="$AD_FPFLAGS"
RUN make -j"$(nproc)" -C smina/build
RUN strip smina/build/fromsmina smina/build/tosmina smina/build/smina
RUN /scripts/shlib.sh boost smina/build/fromsmina smina/build/tosmina smina/build/smina > deps

# ============================================================
# Stage 1 — Build QuickVina2
# ============================================================
FROM builder AS qvina_builder
WORKDIR /src
COPY tools/qvina ./qvina
RUN /scripts/apt.sh dev qvina/deps
RUN make -j"$(nproc)" -C qvina/repo/build/linux/release CFLAGS="$AD_FPFLAGS"
RUN <<EOF
set -eu
strip qvina/repo/build/linux/release/vina
mv qvina/repo/build/linux/release/vina qvina/repo/build/linux/release/qvina
EOF

# ============================================================
# Stage 1 — Build NGL
# ============================================================
# pnpm image pinned by version and digest; the exact NGL release is fixed by
# tools/ngl/pnpm-lock.yaml, so a lockfile drift fails the build.
FROM ${PNPM_IMAGE} AS ngl_builder
WORKDIR /src
COPY tools/ngl .
RUN pnpm install --frozen-lockfile
RUN mv node_modules/ngl/dist ngl

# ============================================================
# Stage 1 — Build FPocket
# ============================================================

FROM builder AS fpocket_builder
WORKDIR /src
COPY tools/fpocket ./fpocket
RUN /scripts/patch.sh fpocket
RUN /scripts/apt.sh dev fpocket/deps
# fpocket links VMD's prebuilt libmolfile_plugin.a, shipped only for
# LINUXAMD64. On arm64 a stub (tools/fpocket/molfile_stub.c) stands in: the
# plugins serve only mmCIF/topology/trajectory input, never the PDB in/out
# AutoDocker uses, and any plugin use exits with a clear error.
RUN <<EOF
set -eu
case "$(uname -m)" in
	x86_64)
		arch=LINUXAMD64
		;;
	aarch64)
		arch=LINUXARM64
		dir="fpocket/repo/plugins/$arch/molfile"
		mkdir -p "$dir"
		cp fpocket/repo/plugins/LINUXAMD64/molfile/libmolfile_plugin.h "$dir/"
		# shellcheck disable=SC2086
		gcc -O2 $AD_FPFLAGS -I fpocket/repo/plugins/include \
			-c fpocket/molfile_stub.c -o "$dir/molfile_stub.o"
		ar rcs "$dir/libmolfile_plugin.a" "$dir/molfile_stub.o"
		rm "$dir/molfile_stub.o"
		;;
	*)
		echo "unsupported build architecture: $(uname -m)" >&2
		exit 1
		;;
esac
make -j"$(nproc)" -C fpocket/repo CXX=g++ ARCH="$arch" \
	CC="gcc $AD_FPFLAGS" CCQHULL="gcc $AD_FPFLAGS"
strip fpocket/repo/bin/*
mkdir -p /out/molfile
cp -a "fpocket/repo/plugins/$arch/molfile/." /out/molfile/
EOF

# ============================================================
# Stage 1 — Build Rootless Mode ACL Manager
# ============================================================

FROM builder AS aclmgr_builder
WORKDIR /src
COPY tools/aclmgr.c .
RUN gcc aclmgr.c -o aclmgr
RUN strip aclmgr

# ============================================================
# Stage 1.5 — Accumulate deps
# ============================================================
FROM ${DEBIAN_IMAGE} AS deps

WORKDIR /acc
COPY deps base
COPY --from=smina_builder /src/deps smina
RUN cat ./* > /deps
# ============================================================
# Stage 2 — Minimal runtime image
# ============================================================
FROM ${PYTHON_IMAGE}

ARG PYTHON_IMAGE
ARG DEBIAN_SNAPSHOT
# The Open Babel build the published results were produced with. The build
# fails if the package set above resolves to anything else.
ARG OPENBABEL_VERSION=3.1.1+dfsg-9+b3

LABEL org.opencontainers.image.title="AutoDocker" \
      org.opencontainers.image.description="Containerized virtual-screening pipeline (protein/ligand prep, Vina docking, ranked reports)." \
      org.opencontainers.image.url="https://github.com/mostakimsakib0/autodocker" \
      org.opencontainers.image.source="https://github.com/mostakimsakib0/autodocker" \
      org.opencontainers.image.documentation="https://github.com/mostakimsakib0/autodocker#readme" \
      org.opencontainers.image.vendor="Mostakim Sakib" \
      org.opencontainers.image.licenses="MIT" \
      org.opencontainers.image.base.name="${PYTHON_IMAGE}" \
      io.github.mostakimsakib0.autodocker.debian-snapshot="${DEBIAN_SNAPSHOT}"

COPY --from=deps /deps /deps
COPY scripts/apt.sh /apt.sh
RUN <<EOF
set -eu
export DEBIAN_FRONTEND=noninteractive
if [ -n "${DEBIAN_SNAPSHOT}" ]; then
	sed -i -E "s#^URIs: http://deb\.debian\.org/(debian|debian-security)\$#URIs: http://snapshot.debian.org/archive/\1/${DEBIAN_SNAPSHOT}#" \
		/etc/apt/sources.list.d/debian.sources
	test "$(grep -c "^URIs: http://snapshot.debian.org/archive/.*/${DEBIAN_SNAPSHOT}\$" /etc/apt/sources.list.d/debian.sources)" -eq 2
	echo 'Acquire::Check-Valid-Until "false";' > /etc/apt/apt.conf.d/99snapshot
fi
apt-get update -o Acquire::Retries=5
/apt.sh run /deps
rm -f /apt.sh
# dpkg-query -W prints "openbabel<TAB>VERSION"; keep the last field.
installed="$(dpkg-query -W openbabel)"
installed="${installed##*[[:space:]]}"
if [ "$installed" != "${OPENBABEL_VERSION}" ]; then
	echo "openbabel ${installed} installed, expected ${OPENBABEL_VERSION}" >&2
	exit 1
fi
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
    /out/molfile/ \
    /usr/local/lib/

COPY --from=ngl_builder \
    /src/ngl /ngl

COPY --from=aclmgr_builder \
	/src/aclmgr /autodocker/aclmgr

RUN useradd --create-home --uid 1000 dockuser
# /autodocker stays root-owned and read-only for dockuser: aclmgr runs as root
# at every start, so a dockuser-writable copy could be replaced. Bytecode is
# compiled now because dockuser cannot write __pycache__ there at run time.
RUN python3 -m compileall -q /autodocker \
	&& chmod -R a+rX,go-w /autodocker \
	&& mkdir -p /workspace && chown -R dockuser:dockuser /workspace

# Release metadata, declared last so a new commit (VCS_REF) does not
# invalidate the cached layers above. autodocker/VERSION is the single source
# of truth: a build-arg that disagrees with it fails the build instead of
# labelling the image with the wrong release.
ARG AUTODOCKER_VERSION=""
ARG VCS_REF=""
RUN <<EOF
set -eu
test -s /autodocker/VERSION
file_version="$(cat /autodocker/VERSION)"
if [ -n "${AUTODOCKER_VERSION}" ] && [ "${AUTODOCKER_VERSION}" != "${file_version}" ]; then
	echo "AUTODOCKER_VERSION=${AUTODOCKER_VERSION} does not match autodocker/VERSION (${file_version})" >&2
	exit 1
fi
EOF

LABEL org.opencontainers.image.version="${AUTODOCKER_VERSION}" \
      org.opencontainers.image.revision="${VCS_REF}"
ENV AUTODOCKER_VERSION="${AUTODOCKER_VERSION}"

WORKDIR /workspace
VOLUME ["/workspace"]

ENTRYPOINT ["/autodocker/aclmgr"]

# Reproducibility pins (override with --build-arg only on purpose):
#  - Base images are pinned by multi-arch index digest (each index covers
#    linux/amd64 and linux/arm64); the tag is kept for readability.
#
# Architectures: linux/amd64 and linux/arm64, each built natively (see
# .github/workflows/publish.yml). Floating-point code is compiled with
# -ffp-contract=off (AD_FPFLAGS) on both, so the compiler never fuses
# multiply-adds; GCC does so by default on arm64 but cannot on baseline
# x86-64, and the fusion changes rounding and hence scores.
#  - Debian packages (toolchain, Boost, Eigen, runtime libraries) come from
#    snapshot.debian.org at a fixed timestamp, so every build resolves the
#    exact same package versions. Old snapshots are never removed, unlike
#    superseded versions on deb.debian.org, so version pins cannot rot.
#    Set DEBIAN_SNAPSHOT="" to use the live mirrors instead.
#  - Open Babel is built from tools/openbabel/repo (submodule pinned to the
#    openbabel-3-2-1 tag), not taken from Debian (3.1.1); SMINA links it too.
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
# Stage 1 — Build Open Babel
# ============================================================
# Open Babel 3.2.1 (tools/openbabel/repo, tag openbabel-3-2-1 = commit
# 0e94434fa75c9f61095023e3c12e0d5f2ac035ff). 3.2 honours OB_RANDOM_SEED,
# which the pipeline sets on every call, and no longer adds a hydrogen to
# pre-protonated Lys/Arg (see docs/openbabel.md). Only the command-line tools
# and the library are built: no GUI, no language bindings, no tests, and no
# optional component that CMake would download (MaeParser, CoordGen,
# RapidJSON) or AutoDocker does not use (Cairo PNG depiction, libxml2 CML).
# OpenMP stays off (single-threaded force fields). USE_RANDOMMT=ON selects
# the generator OB_RANDOM_SEED seeds. SOURCE_DATE_EPOCH (the commit time of
# the tag, 2026-07-11T19:21:27Z) fixes the __DATE__/__TIME__ that "obabel -V"
# prints, so the binary does not record when the image was built.
FROM builder AS openbabel_builder
ARG OPENBABEL_VERSION=3.2.1
ARG OPENBABEL_SOURCE_DATE_EPOCH=1783797687
WORKDIR /src
COPY tools/openbabel ./openbabel
RUN /scripts/apt.sh dev openbabel/deps
RUN <<EOF
set -eu
test -f openbabel/repo/CMakeLists.txt || {
	echo "tools/openbabel/repo is empty: run git submodule update --init" >&2
	exit 1
}
export SOURCE_DATE_EPOCH="${OPENBABEL_SOURCE_DATE_EPOCH}"
cmake -S openbabel/repo -B openbabel/build \
	-DCMAKE_BUILD_TYPE=Release \
	-DCMAKE_INSTALL_PREFIX=/usr/local \
	-DCMAKE_C_FLAGS="$AD_FPFLAGS" -DCMAKE_CXX_FLAGS="$AD_FPFLAGS" \
	-DBUILD_GUI=OFF -DPYTHON_BINDINGS=OFF -DRUN_SWIG=OFF -DENABLE_TESTS=OFF \
	-DENABLE_OPENMP=OFF -DOPTIMIZE_NATIVE=OFF -DUSE_RANDOMMT=ON \
	-DWITH_MAEPARSER=OFF -DWITH_COORDGEN=OFF -DWITH_JSON=OFF \
	-DCMAKE_DISABLE_FIND_PACKAGE_Boost=ON \
	-DCMAKE_DISABLE_FIND_PACKAGE_Cairo=ON \
	-DCMAKE_DISABLE_FIND_PACKAGE_LibXml2=ON
make -j"$(nproc)" -C openbabel/build
make -C openbabel/build install DESTDIR=/out/full
prefix=/out/full/usr/local
strip "$prefix"/bin/* "$prefix"/lib/lib*.so.*.* "$prefix"/lib/openbabel/"${OPENBABEL_VERSION}"/*.so
# Runtime subset: tools, shared libraries, plugins and data (no headers,
# CMake or pkg-config files).
mkdir -p /out/runtime/lib /out/runtime/share
cp -a "$prefix"/bin /out/runtime/
cp -a "$prefix"/lib/lib*.so.* "$prefix"/lib/openbabel /out/runtime/lib/
cp -a "$prefix"/share/openbabel /out/runtime/share/
# The build must be the pinned release, load its plugins and write PDBQT.
export LD_LIBRARY_PATH="$prefix/lib"
export BABEL_LIBDIR="$prefix/lib/openbabel/${OPENBABEL_VERSION}"
export BABEL_DATADIR="$prefix/share/openbabel/${OPENBABEL_VERSION}"
version="$("$prefix"/bin/obabel -V)"
formats="$("$prefix"/bin/obabel -L formats)"
pdbqt="$("$prefix"/bin/obabel -:CCO -opdbqt -h)"
case "$version" in "Open Babel ${OPENBABEL_VERSION} "*) ;; *) echo "built: $version" >&2; exit 1 ;; esac
case "$formats" in *"pdbqt -- "*) ;; *) echo "no PDBQT format plugin" >&2; exit 1 ;; esac
case "$pdbqt" in *ROOT*ENDROOT*) ;; *) echo "PDBQT conversion failed" >&2; exit 1 ;; esac
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
# SMINA links the same Open Babel 3.2.1 (headers and library from the
# openbabel_builder stage; the runtime image ships the library).
COPY --from=openbabel_builder /out/full/usr/local/ /usr/local/
RUN ldconfig
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
COPY tools/aclmgr.c tools/aclmgr.*.h .
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
# The Open Babel release built in openbabel_builder (tools/openbabel/repo);
# the build fails if obabel reports anything else.
ARG OPENBABEL_VERSION=3.2.1

LABEL org.opencontainers.image.title="AutoDocker" \
      org.opencontainers.image.description="Containerized virtual-screening pipeline (protein/ligand prep, Vina docking, ranked reports)." \
      org.opencontainers.image.url="https://github.com/mostakimsakib0/autodocker" \
      org.opencontainers.image.source="https://github.com/mostakimsakib0/autodocker" \
      org.opencontainers.image.documentation="https://github.com/mostakimsakib0/autodocker#readme" \
      org.opencontainers.image.vendor="Mostakim Sakib" \
      org.opencontainers.image.licenses="MIT" \
      org.opencontainers.image.base.name="${PYTHON_IMAGE}" \
      io.github.mostakimsakib0.autodocker.debian-snapshot="${DEBIAN_SNAPSHOT}" \
      io.github.mostakimsakib0.autodocker.openbabel-version="${OPENBABEL_VERSION}"

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
# Debian's Open Babel (3.1.1) must not come in through a dependency.
for pkg in openbabel libopenbabel7; do
	if dpkg -s "$pkg" > /dev/null 2>&1; then
		echo "Debian package $pkg is installed; the image ships Open Babel ${OPENBABEL_VERSION}" >&2
		exit 1
	fi
done
EOF

# Open Babel (tools, library, plugins, data) under /usr/local. These are
# also the compiled-in plugin and data directories; the variables make them
# explicit.
COPY --from=openbabel_builder /out/runtime/ /usr/local/
ENV BABEL_LIBDIR=/usr/local/lib/openbabel/${OPENBABEL_VERSION} \
    BABEL_DATADIR=/usr/local/share/openbabel/${OPENBABEL_VERSION}
RUN <<EOF
set -eu
ldconfig
test "$(command -v obabel)" = /usr/local/bin/obabel
version="$(obabel -V)"
case "$version" in
	"Open Babel ${OPENBABEL_VERSION} "*) ;;
	*)
		echo "obabel -V reports '${version}', expected Open Babel ${OPENBABEL_VERSION}" >&2
		exit 1
		;;
esac
formats="$(obabel -L formats)"
case "$formats" in *"pdbqt -- "*) ;; *) echo "obabel finds no PDBQT format plugin" >&2; exit 1 ;; esac
# The compiled-in directories must work too (an environment without the
# variables, e.g. setpriv --reset-env).
formats="$(env -u BABEL_LIBDIR -u BABEL_DATADIR obabel -L formats)"
case "$formats" in *"pdbqt -- "*) ;; *) echo "obabel needs BABEL_LIBDIR" >&2; exit 1 ;; esac
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

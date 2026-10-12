#!/usr/bin/env bash
# Run one command inside the benchmark image with the fixed mount layout of
# .github/workflows/benchmarks.yml:
#
#   /repo   the repository checkout (read-only; scripts/ and autodocker/)
#   /data   BENCH_DATA, the prepared benchmark inputs (writable: DUD-E
#           docking writes <target>/dock next to its inputs)
#   /out    BENCH_OUT, this job's results (also the working directory)
#   /tools  BENCH_TOOLS, optional: baseline virtual environments, raw Vina
#
# Every job uses the same container paths, so the absolute paths recorded in
# logs, run manifests and summary.json files agree across shards. The
# command runs as the calling user (files in the bind mounts stay owned by
# the runner) with HOME=/tmp, bypassing the image's entrypoint.
#
# Usage: scripts/ci/in_image.sh PROGRAM [ARGS...]
#   BENCH_IMAGE  image (default autodocker:bench)
#   BENCH_DATA, BENCH_OUT (required), BENCH_TOOLS (optional) host directories
set -euo pipefail

if [ "$#" -lt 1 ]; then
	echo "usage: $0 PROGRAM [ARGS...]" >&2
	exit 2
fi
: "${BENCH_DATA:?set BENCH_DATA to the host input directory}"
: "${BENCH_OUT:?set BENCH_OUT to the host output directory}"
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
IMAGE="${BENCH_IMAGE:-autodocker:bench}"

mkdir -p "$BENCH_DATA" "$BENCH_OUT"
mounts=(-v "$REPO:/repo:ro" -v "$BENCH_DATA:/data" -v "$BENCH_OUT:/out")
if [ -n "${BENCH_TOOLS:-}" ]; then
	mkdir -p "$BENCH_TOOLS"
	mounts+=(-v "$BENCH_TOOLS:/tools")
fi

program="$1"
shift
exec docker run --rm \
	--user "$(id -u):$(id -g)" \
	--entrypoint "$program" \
	"${mounts[@]}" \
	-w /out \
	-e HOME=/tmp \
	-e PYTHONDONTWRITEBYTECODE=1 \
	-e PYTHONUNBUFFERED=1 \
	"$IMAGE" "$@"

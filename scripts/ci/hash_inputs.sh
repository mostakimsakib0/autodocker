#!/usr/bin/env bash
# SHA-256 manifest of prepared benchmark inputs, in the format the workflow's
# prepare job uploads (inputs/SHA256SUMS): `sha256sum` lines
# "<hex>  <path>", path relative to ROOT, every regular file below the given
# sub-directories, sorted bytewise (LC_ALL=C). Run it on any machine to get a
# manifest comparable with the CI one (scripts/ci/compare_inputs.sh):
#
#   scripts/ci/hash_inputs.sh ~/Downloads/Research/benchmark_data/v2 \
#       casf2016_full dude dude_prepare > laptop.sha256
#
# Usage: scripts/ci/hash_inputs.sh ROOT SUBDIR...
set -euo pipefail
export LC_ALL=C

if [ "$#" -lt 2 ]; then
	echo "usage: $0 ROOT SUBDIR..." >&2
	exit 2
fi
root="$1"
shift
cd "$root"
for d in "$@"; do
	if [ ! -d "$d" ]; then
		echo "$root/$d: not a directory" >&2
		exit 1
	fi
done
find "$@" -type f -print0 | sort -z | xargs -0 -r sha256sum --

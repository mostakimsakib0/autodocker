#!/usr/bin/env bash
# Load the image built by the workflow's build job and check that it is that
# image: its ID must equal build-info.json's image_id.
#
# Usage: scripts/ci/load_image.sh DIR   (DIR holds image.tar.zst and build-info.json)
#   BENCH_IMAGE  tag to give the loaded image (default autodocker:bench)
set -euo pipefail

dir="${1:?usage: $0 DIR}"
tag="${BENCH_IMAGE:-autodocker:bench}"
want="$(jq -r .image_id "$dir/build-info.json")"
sum="$(jq -r .archive_sha256 "$dir/build-info.json")"

echo "${sum}  $dir/image.tar.zst" | sha256sum -c -
zstd -d --stdout "$dir/image.tar.zst" | docker load
got="$(docker image inspect --format '{{.Id}}' "$want")"
if [ "$got" != "$want" ]; then
	echo "loaded image $got, build job recorded $want" >&2
	exit 1
fi
docker tag "$want" "$tag"
echo "image $want loaded as $tag"

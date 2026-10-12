#!/usr/bin/env bash
# Save the benchmark image for the other jobs of benchmarks.yml and record
# what it is: DIR/image.tar.zst (`docker save | zstd -19`) and
# DIR/build-info.json (image ID, archive SHA-256, OCI labels, and the
# versions of AutoDocker and of every tool and Python package in the image).
#
# Usage: scripts/ci/save_image.sh IMAGE DIR
set -euo pipefail

image="${1:?usage: $0 IMAGE DIR}"
dir="${2:?usage: $0 IMAGE DIR}"
mkdir -p "$dir"

in_image() {
	docker run --rm --entrypoint sh "$image" -c "$1"
}

image_id="$(docker image inspect --format '{{.Id}}' "$image")"
labels="$(docker image inspect --format '{{json .Config.Labels}}' "$image")"
arch="$(docker image inspect --format '{{.Architecture}}' "$image")"
runner_version="$(docker run --rm --entrypoint python3 "$image" /autodocker/runner.py --version)"
vina="$(in_image 'vina --version 2>&1 | head -n 1')"
obabel="$(in_image 'obabel -V')"
smina="$(in_image 'smina --version 2>&1 | head -n 1')"
qvina="$(in_image 'qvina --version 2>&1 | head -n 1')"
python="$(in_image 'python3 --version')"
rdkit="$(docker run --rm --entrypoint python3 "$image" -c 'from rdkit import rdBase; print(rdBase.rdkitVersion)')"
freeze="$(in_image 'pip freeze --all --disable-pip-version-check' | jq -R . | jq -s .)"
binaries="$(in_image 'sha256sum /usr/local/bin/vina /usr/local/bin/obabel /usr/local/bin/smina /usr/local/bin/qvina /usr/local/bin/fpocket' |
	jq -R 'split("  ") | {(.[1]): .[0]}' | jq -s add)"

docker save "$image" | zstd -T0 -19 -q -o "$dir/image.tar.zst"
archive_sha256="$(sha256sum "$dir/image.tar.zst" | cut -d ' ' -f 1)"

jq -n \
	--arg image_id "$image_id" \
	--arg tag "$image" \
	--arg arch "$arch" \
	--arg archive_sha256 "$archive_sha256" \
	--argjson archive_bytes "$(stat -c %s "$dir/image.tar.zst")" \
	--argjson labels "$labels" \
	--arg runner_version "$runner_version" \
	--arg vina "$vina" \
	--arg obabel "$obabel" \
	--arg smina "$smina" \
	--arg qvina "$qvina" \
	--arg python "$python" \
	--arg rdkit "$rdkit" \
	--argjson pip_freeze "$freeze" \
	--argjson binaries_sha256 "$binaries" \
	--arg git_sha "${GITHUB_SHA:-}" \
	--arg built_at "$(date -u +%Y-%m-%dT%H:%M:%SZ)" \
	'{image_id: $image_id, tag: $tag, architecture: $arch, archive: "image.tar.zst",
	  archive_sha256: $archive_sha256, archive_bytes: $archive_bytes, labels: $labels,
	  git_sha: $git_sha, built_at: $built_at,
	  versions: {autodocker: $runner_version, vina: $vina, openbabel: $obabel, smina: $smina,
	             qvina: $qvina, python: $python, rdkit: $rdkit},
	  binaries_sha256: $binaries_sha256, pip_freeze: $pip_freeze}' > "$dir/build-info.json"
cat "$dir/build-info.json"

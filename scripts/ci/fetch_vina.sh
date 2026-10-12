#!/usr/bin/env bash
# The upstream AutoDock Vina release binary (statically linked) for the J5
# baseline, verified by SHA-256: the same release and checksum as
# .github/workflows/ci.yml.
#
# Usage: scripts/ci/fetch_vina.sh DEST
#   VINA_VERSION (1.2.7), VINA_SHA256 (checksum of vina_<version>_linux_x86_64)
set -euo pipefail

dest="${1:?usage: $0 DEST}"
version="${VINA_VERSION:-1.2.7}"
sha256="${VINA_SHA256:-f31f774f723bba7bbe6e9d1c47577020eea9a8da16424284c043d22593570644}"
url="https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v${version}/vina_${version}_linux_x86_64"

mkdir -p "$(dirname "$dest")"
tmp="$dest.part"
curl -fsSL --retry 5 -o "$tmp" "$url"
echo "${sha256}  $tmp" | sha256sum -c -
chmod 0755 "$tmp"
mv "$tmp" "$dest"
"$dest" --version

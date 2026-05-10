#!/bin/bash
set -euo pipefail
env="$1"
shift
for tool in "$@"; do
	sed -nE '/^'"$env"'/s/^'"$env"'[[:space:]]+//p' "$tool"
done | xargs apt-get install -y --no-install-recommends
rm -rf /var/lib/apt/lists/*

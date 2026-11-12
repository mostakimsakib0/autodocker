#!/bin/bash
set -euo pipefail
export DEBIAN_FRONTEND=noninteractive
add()
{
	apt-get install -y \
		--no-install-recommends \
		-o Acquire::Retries=5 \
		-o Acquire::http::Timeout=60 "$@"
}

if [ "$1" = "--init" ]; then
	apt-get update
	shift
fi
env="$1"
shift
case "$env" in
	builder) add build-essential git quilt ;;
	*)
		for tool in "$@"; do
			sed -nE '/^'"$env"'/s/^'"$env"'[[:space:]]+//p' "$tool/deps"
		done | (
			add $(xargs)
			rm "$0"
		)
		;;
esac

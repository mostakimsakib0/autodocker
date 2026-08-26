#!/bin/bash

declare -i eid=0 len=0
if [[ -f /proc/self/uid_map ]]; then
    read -r _ eid len < <(
        grep -E '^\s+0\s+' /proc/self/uid_map
    )
fi
args=()

for arg in "$@"; do
	case "$arg" in
		-s | --shell | --debug) args=(--pty) ;;
        -r | --root) exec bash -i;;
	esac
done

if ((eid >= 1000 && eid <= 60000 && len == 1)); then
    chown -R dockuser:dockuser /workspace
fi

runuser "${args[@]}" -u dockuser -- "${BASH_SOURCE[0]%/*}/entry.sh" "$@"

if ((eid >= 1000 && eid <= 60000 && len == 1)); then
    chown -R root:root /workspace
fi

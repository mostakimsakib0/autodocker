#!/usr/bin/env bash
# Compare two input manifests written by scripts/ci/hash_inputs.sh (e.g. the
# CI prepare job's inputs/SHA256SUMS and one made on the laptop).
#
# Files whose content is expected to differ between machines are listed but
# do not fail the comparison:
#   *_ligand_modelserver.sdf   RCSB ModelServer response metadata (request
#                              id, timestamp, timings)
#   dude_prepare/summary.json  records the host's worker count and the
#                              runner/Vina versions of that machine
# COMPARE_IGNORE (extended regex on the path) replaces this list.
#
# Usage: scripts/ci/compare_inputs.sh A.sha256 B.sha256
# Exit status: 0 = every other file is present in both with equal SHA-256.
set -euo pipefail
export LC_ALL=C

if [ "$#" -ne 2 ]; then
	echo "usage: $0 A.sha256 B.sha256" >&2
	exit 2
fi
export IGNORE_RE="${COMPARE_IGNORE:-(_ligand_modelserver\.sdf|^dude_prepare/summary\.json)\$}"

table() {
	# "<path>\t<sha256>", sorted by path
	awk '{ print substr($0, 67) "\t" $1 }' "$1" | sort -t "$(printf '\t')" -k1,1
}

tab="$(printf '\t')"
join -t "$tab" -a 1 -a 2 -e MISSING -o 0,1.2,2.2 <(table "$1") <(table "$2") |
	awk -F '\t' '
		$2 == $3 { same++; next }
		{
			kind = ($2 == "MISSING") ? "only-in-B" : (($3 == "MISSING") ? "only-in-A" : "differs")
			if ($1 ~ ENVIRON["IGNORE_RE"]) { ignored++; printf "%-10s %s (expected)\n", kind, $1; next }
			bad++
			printf "%-10s %s\n", kind, $1
		}
		END {
			printf "identical: %d, expected differences: %d, unexpected differences: %d\n", same, ignored, bad
			exit (bad > 0)
		}
	'

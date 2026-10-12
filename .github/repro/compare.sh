#!/usr/bin/env bash
# Compare a case-study workspace written by run_case.sh with the laptop
# reference in 1hsg/expected/ (ranking.csv and the SHA-256 of every pose
# file) and print a Markdown report. Exit 0 if ranking.csv and every pose
# file are byte-identical, 1 otherwise (measured, not guaranteed).
#
# Usage: .github/repro/compare.sh LABEL WORKSPACE
set -uo pipefail

label="${1:?usage: $0 LABEL WORKSPACE}"
out="${2:?usage: $0 LABEL WORKSPACE}/output"
ref="$(cd "$(dirname "$0")" && pwd)/1hsg/expected"

echo "### $label"
echo
if [ ! -f "$out/ranking.csv" ]; then
	echo "Result: **no output** (the run produced no ranking.csv)"
	exit 1
fi
status=0
if cmp -s "$ref/ranking.csv" "$out/ranking.csv"; then
	echo "- ranking.csv: byte-identical to the laptop run"
else
	status=1
	echo "- ranking.csv: **differs** from the laptop run"
fi
n="$(wc -l < "$ref/poses.sha256")"
same="$(cd "$out" && sha256sum -c --quiet "$ref/poses.sha256" 2>/dev/null; true)"
bad="$(printf '%s\n' "$same" | grep -c 'FAILED' || true)"
echo "- pose files: $((n - bad)) of $n byte-identical"
[ "$bad" -eq 0 ] || status=1
extra="$(find "$out/docked" -maxdepth 1 -name "*_out.pdbqt" 2>/dev/null | wc -l)"
[ "$extra" -eq "$n" ] || { status=1; echo "- pose files written: $extra (laptop: $n)"; }
if ! cmp -s "$ref/ranking.csv" "$out/ranking.csv"; then
	echo
	echo '```'
	diff "$ref/ranking.csv" "$out/ranking.csv" | head -n 20
	echo '```'
fi
echo
[ "$status" -eq 0 ] && echo "Result: **identical**" || echo "Result: **differs**"
exit "$status"

#!/usr/bin/env bash
# Check AutoDocker output directories with scripts/verify_run.py (consistency
# of ranking.csv, poses, logs and run_manifest.json hashes; no re-docking).
# Writes OUT/<n>_<name>.json per run and OUT/verify.tsv ("<run dir>\t<exit
# status>"). Without scripts/verify_run.py (checkouts before 2.0.0's
# verify_run.py) it writes OUT/SKIPPED and succeeds.
#
# Usage: scripts/ci/verify_runs.sh OUT RUN_DIR...
#   PYTHON (python3)
# Exit status: 0 when every run passes (warnings allowed), 1 otherwise.
set -uo pipefail

out="${1:?usage: $0 OUT RUN_DIR...}"
shift
python="${PYTHON:-python3}"
verify="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/verify_run.py"
mkdir -p "$out"
if [ ! -f "$verify" ]; then
	echo "scripts/verify_run.py not in this checkout" > "$out/SKIPPED"
	echo "[!] scripts/verify_run.py not in this checkout: runs not verified"
	exit 0
fi

status=0
n=0
: > "$out/verify.tsv"
for run in "$@"; do
	n=$((n + 1))
	name="$(printf '%04d_%s' "$n" "$(basename "$run")")"
	"$python" "$verify" "$run" --json > "$out/$name.json"
	rc=$?
	printf '%s\t%s\n' "$run" "$rc" >> "$out/verify.tsv"
	if [ "$rc" -ne 0 ]; then
		echo "[x] verify_run.py $run: exit $rc"
		status=1
	fi
done
echo "verified $n run(s), exit $status"
exit "$status"

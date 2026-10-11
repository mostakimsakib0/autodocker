#!/usr/bin/env bash
# Smoke screen + determinism check for a built AutoDocker image.
#
# Screens the bundled fixtures (crambin, PDB 1CRN, protein atoms only; phenol
# and benzoic acid) twice through the image's entrypoint (aclmgr -> entry.sh
# -> runner.py) and requires ranking.csv and every docked pose
# (docked/*_out.pdbqt) to be byte-identical between the two runs. Both runs
# mount their workspace at the same container path, so the absolute paths
# recorded in ranking.csv are identical as well.
#
# Usage: .github/smoke/run_smoke.sh IMAGE [WORKDIR]
#   CONTAINER_ENGINE    docker (default) or podman
#   SMOKE_ENGINE_ARGS   extra "run" options, e.g. "--cpus 2"
#
# Exit status: 0 = both runs succeeded and are identical, 1 = otherwise.
set -euo pipefail

IMAGE="${1:?usage: $0 IMAGE [WORKDIR]}"
ROOT="${2:-$(mktemp -d)}"
ENGINE="${CONTAINER_ENGINE:-docker}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
read -r -a ENGINE_ARGS <<< "${SMOKE_ENGINE_ARGS:-}"

screen() {
	local ws="$1"
	mkdir -p "$ws/ligs"
	cp "$HERE/protein.pdb" "$ws/protein.pdb"
	cp "$HERE"/ligs/*.sdf "$ws/ligs/"
	# Under rootful Docker the pipeline runs as the image's dockuser
	# (uid 1000), which must be able to create output/ in the bind mount.
	chmod -R a+rwX "$ws"
	"$ENGINE" run --rm "${ENGINE_ARGS[@]}" \
		-v "$ws:/workspace" \
		-e EXHAUSTIVENESS=2 \
		-e BINDING_MODES=3 \
		-e PROCESSES=1 \
		-e RESUME=0 \
		"$IMAGE"
}

for run in 1 2; do
	echo "[*] smoke run #$run -> $ROOT/run$run"
	screen "$ROOT/run$run"
done

out1="$ROOT/run1/output"
out2="$ROOT/run2/output"
status=0

for f in ranking.csv Top_hits.txt metrics.txt; do
	if [ ! -s "$out1/$f" ]; then
		echo "[x] missing or empty output: $f" >&2
		status=1
	fi
done

shopt -s nullglob
poses=("$out1"/docked/*_out.pdbqt)
if [ "${#poses[@]}" -ne 2 ]; then
	echo "[x] expected 2 docked poses, found ${#poses[@]}" >&2
	status=1
fi

for f in "$out1/ranking.csv" "${poses[@]}"; do
	rel="${f#"$out1"/}"
	if ! cmp -s "$f" "$out2/$rel"; then
		echo "[x] NON-DETERMINISTIC: $rel differs between runs" >&2
		diff "$f" "$out2/$rel" | head -20 >&2 || true
		status=1
	fi
done

if [ "$status" -eq 0 ]; then
	echo "[ok] DETERMINISTIC: ranking.csv and ${#poses[@]} docked poses are byte-identical"
fi
exit "$status"

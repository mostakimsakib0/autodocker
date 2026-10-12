#!/usr/bin/env bash
# Run the paper's 1HSG case-study screen (--chain all --pockets 2) in IMAGE on
# this machine, as AutoDocker_Paper/rerun/run_v200.sh does on the laptop: the
# workspace gets protein.pdb, ligs/, output/ and a log whose first line is the
# exact container command, plus host.txt describing the machine.
#
# Usage: .github/repro/run_case.sh IMAGE WORKSPACE
set -euo pipefail

image="${1:?usage: $0 IMAGE WORKSPACE}"
ws="${2:?usage: $0 IMAGE WORKSPACE}"
here="$(cd "$(dirname "$0")" && pwd)"

rm -rf "$ws"
mkdir -p "$ws"
cp "$here/1hsg/protein.pdb" "$ws/"
cp -r "$here/1hsg/ligs" "$ws/"
# The image runs as dockuser (UID 1000), not as the runner's user: let it
# write output/ (as .github/smoke/run_smoke.sh does).
chmod -R a+rwX "$ws"

{
	echo "uname: $(uname -srm)"
	echo "runner: ${RUNNER_NAME:-?} (${RUNNER_OS:-?} ${RUNNER_ARCH:-?}, ${ImageOS:-?} ${ImageVersion:-?})"
	echo "image: $image $(docker image inspect --format '{{.Id}} {{.Architecture}}' "$image")"
	echo
	lscpu
} > "$ws/host.txt"

timer=()
[ -x /usr/bin/time ] && timer=(/usr/bin/time -v)
(
	set +e
	cd "$ws" || exit 1
	echo "\$ docker run --rm -v \"\$PWD:/workspace\" $image --chain all --pockets 2"
	date
	"${timer[@]}" docker run --rm -v "$PWD:/workspace" "$image" --chain all --pockets 2
	echo "container exit $?"
	date
) > "$ws/log" 2>&1
tail -n 40 "$ws/log"
grep -q '^container exit 0$' "$ws/log"

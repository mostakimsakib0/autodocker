#!/usr/bin/env bash
# Virtual environments of the J5 baseline (scripts/benchmark_baseline.py),
# every package pinned exactly, including the transitive ones. The pins are
# the environments benchmark_baseline.py was validated with (BENCHMARKS.md,
# J5): Meeko 0.8.0 with RDKit 2026.09.1 and reduce 4.15 (reduce-binary
# 4.15.0 wheel), and EasyDock 1.3.1 with the Vina 1.2.7 Python bindings.
# Installed from wheels only, without dependency resolution (--no-deps), then
# checked with `pip check`; `pip freeze --all` is written to PREFIX.freeze.txt.
#
# Run it with the Python that will run the tools: in the workflow, inside the
# benchmark image (Python 3.12), e.g.
#   scripts/ci/in_image.sh bash /repo/scripts/ci/baseline_env.sh baseline /tools/baseline
#
# Usage: scripts/ci/baseline_env.sh {baseline|easydock} PREFIX
set -euo pipefail

kind="${1:?usage: $0 {baseline|easydock} PREFIX}"
prefix="${2:?usage: $0 {baseline|easydock} PREFIX}"
python="${PYTHON:-python3}"

# Meeko + reduce (receptor and ligand preparation of the standard workflow).
baseline_pins='
biopython==1.88
gemmi==0.7.5
markdown-it-py==4.2.0
mdurl==0.1.2
meeko==0.8.0
molscrub==0.3.0
numpy==2.5.4
pillow==12.3.0
prody==2.6.1
pygments==2.21.0
pyparsing==3.1.1
rdkit==2026.9.1
reduce-binary==4.15.0
rich==15.0.0
scipy==1.18.1
'

# EasyDock 1.3.1 and its runtime (Dask, the Vina Python bindings, Meeko).
easydock_pins='
biopython==1.88
click==8.5.0
cloudpickle==3.1.2
dask==2026.8.0
distributed==2026.8.0
easydock==1.3.1
fsspec==2026.9.0
gemmi==0.7.5
jinja2==3.1.6
joblib==1.6.0
locket==1.0.0
markdown-it-py==4.2.0
markupsafe==3.0.4
mdurl==0.1.2
meeko==0.8.0
molscrub==0.3.0
msgpack==1.2.3
narwhals==2.27.1
numpy==2.5.4
packaging==26.3
partd==1.4.2
pillow==12.3.0
prody==2.6.1
psutil==7.2.2
pygments==2.21.0
pyparsing==3.1.1
pyyaml==6.0.3
rdkit==2026.9.1
rich==15.0.0
scikit-learn==1.9.1
scipy==1.18.1
setuptools==84.0.0
sortedcontainers==2.4.0
tblib==3.2.2
threadpoolctl==3.7.0
toolz==1.2.0
tornado==6.5.10
vina==1.2.7
wheel==0.48.0
zict==3.0.0
'

case "$kind" in
	baseline) pins="$baseline_pins" ;;
	easydock) pins="$easydock_pins" ;;
	*)
		echo "usage: $0 {baseline|easydock} PREFIX" >&2
		exit 2
		;;
esac

"$python" -m venv "$prefix"
req="$(mktemp)"
trap 'rm -f "$req"' EXIT
printf '%s\n' "$pins" | sed '/^$/d' > "$req"
"$prefix/bin/pip" install --no-cache-dir --disable-pip-version-check \
	--only-binary=:all: --no-deps -r "$req"
"$prefix/bin/pip" check
"$prefix/bin/pip" freeze --all > "$prefix.freeze.txt"
cat "$prefix.freeze.txt"

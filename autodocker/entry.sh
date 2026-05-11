#!/bin/bash
for arg in "$@"; do
	case "$arg" in
		-s | --shell | --debug) exec bash -i ;;
	esac
done

set -euo pipefail

DIR="$(dirname "$(realpath "${BASH_SOURCE[0]}")")"
WORK="/workspace"
OUT="$WORK/output"
SCRIPT="$DIR/runner.py"

mkdir -p "$OUT"

if [ -n "${PDB:-}" ]; then
	: # Use PDB from the environment.
elif [ -f "$WORK/protein.pdb" ]; then
	PDB="$WORK/protein.pdb"
elif [ -f "$WORK/protein_clean.pdb" ]; then
	PDB="$WORK/protein_clean.pdb"
else
	PDB="$WORK/protein.pdb"
fi

LIGS="$WORK/ligs"

echo "======================================"
echo " PenDrive Virtual Screening Engine"
echo "======================================"

# Check if required files exist
if [ ! -f "$PDB" ]; then
	echo "[!] Protein file not found: $PDB"
	echo "[!] Please place protein.pdb in workspace/ or set PDB=/path/to/file.pdb"
	exit 1
fi

if [ ! -d "$LIGS" ]; then
	echo "[!] Ligands directory not found: $LIGS"
	echo "[!] Please create workspace/ligs/ and add .sdf files"
	exit 1
fi

# Check for ligand files
LIG_COUNT=$(find "$LIGS" -name "*.sdf" | wc -l)
if [ "$LIG_COUNT" -eq 0 ]; then
	echo "[!] No .sdf files found in $LIGS"
	exit 1
fi

echo "[*] Found $LIG_COUNT ligand(s)"
echo "[*] Protein: $PDB"

# Run the pipeline
CHAIN="${CHAIN:-A}"
PROCESSES="${PROCESSES:-${MAX_WORKERS:--1}}"
EXHAUSTIVENESS="${EXHAUSTIVENESS:-${EXHAUST:-8}}"
BINDING_MODES="${BINDING_MODES:-9}"
ENERGY_RANGE="${ENERGY_RANGE:-3.0}"
PADDING="${PADDING:-6}"
TOP_N="${TOP_N:-20}"
TIMEOUT="${TIMEOUT:-900}"

CMD=(
	"$SCRIPT"
	-r "$PDB"
	-l "$LIGS"
	-o "$OUT"
	--chain "$CHAIN"
	--no-admet
	--padding "$PADDING"
	--exhaustiveness "$EXHAUSTIVENESS"
	--binding-modes "$BINDING_MODES"
	--energy-range "$ENERGY_RANGE"
	--timeout "$TIMEOUT"
	--top-n "$TOP_N"
	--no-resume
	-p "$PROCESSES"
)

if [ -n "${POCKETS:-}" ]; then
	CMD+=(--pockets "$POCKETS")
fi

if [ "${NO_MINIMIZE:-0}" = "1" ]; then
	CMD+=(--no-minimize)
fi

python3 "${CMD[@]}"

echo "======================================"
echo " Pipeline Completed Successfully"
echo "======================================"
echo ""
echo "Results saved to: $OUT"
echo "  - ranking.csv: Docking scores"
echo "  - Top_hits.txt: Top-ranked hits"
echo "  - Results_full.txt: Full text report"
echo "  - metrics.txt: Pose metrics"
echo "  - pocket_summary.csv: fpocket pocket scores"
echo "  - ligand_metadata.json: ligand descriptors and source files"
echo "  - grid_box.py: PyMOL grid visualization"
echo "  - docked/: Docked poses"

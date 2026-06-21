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

if [ -n "$INPUT" ] && [ -f "$INPUT" ]; then
    echo "[INFO] Single ligand mode detected: $INPUT"
    LIGAND_MODE="single"
else
    echo "[INFO] Batch mode: scanning $LIGS"
    LIGAND_MODE="batch"
fi


if [ -n "${PDB:-}" ]; then
	: # Use PDB from the environment.
elif [ -f "$WORK/protein.pdb" ]; then
	PDB="$WORK/protein.pdb"
elif [ -f "$WORK/protein_clean.pdb" ]; then
	PDB="$WORK/protein_clean.pdb"
else
	PDB="$WORK/protein.pdb"
fi

INPUT="${1:-}"
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
SDF_COUNT=$(find "$LIGS" -type f -name "*.sdf" | wc -l)
PDBQT_COUNT=$(find "$LIGS" -type f -name "*.pdbqt" | wc -l)
MOL2_COUNT=$(find "$LIGS" -type f -name "*.mol2" | wc -l)
PDB_COUNT=$(find "$LIGS" -type f -name "*.pdb" | wc -l)

TOTAL_LIGS=$((SDF_COUNT + PDBQT_COUNT + MOL2_COUNT + PDB_COUNT))

if [ "$TOTAL_LIGS" -eq 0 ]; then
    echo "[!] No ligand files (.sdf/.pdbqt) found in $LIGS"
    exit 1
fi

echo "[*] Found $TOTAL_LIGS ligands"
echo "    - SDF: $SDF_COUNT"
echo "    - PDBQT: $PDBQT_COUNT"
echo "[*] Found $TOTAL_LIGS ligand(s)"
echo "[*] Protein loaded: $(basename "$PDB")"

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

"${CMD[@]}"

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

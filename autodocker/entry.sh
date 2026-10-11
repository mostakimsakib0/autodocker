#!/bin/bash
# Container entrypoint (started by /autodocker/aclmgr as the unprivileged
# user): resolves receptor/ligand inputs from the workspace and the
# environment, then runs runner.py. All other arguments are passed through to
# runner.py unchanged.

# An interactive shell (as the unprivileged user) is opened only when the
# FIRST argument is -s, --shell or --debug, so these words can still appear as
# values later on the command line. (A root shell is aclmgr's --root-shell.)
case "${1:-}" in
	-s | --shell | --debug) exec bash -i ;;
esac

set -euo pipefail

SELF="$(realpath "${BASH_SOURCE[0]}")"
DIR="$(dirname "$SELF")"
WORK="/workspace"
OUT="$WORK/output"
SCRIPT="$DIR/runner.py"

LIGS="${LIGS:-/workspace/ligs}"
PDB="${PDB:-}"

# -r/-l/-o given on the command line override the workspace defaults
# (argparse keeps the last value); validate and report what will be used.
USER_PDB="" USER_LIGS="" USER_OUT=""
for ((i = 1; i <= $#; i++)); do
	next=$((i + 1))
	case "${!i}" in
		-r|--receptor) USER_PDB="${!next:-}" ;;
		--receptor=*) USER_PDB="${!i#*=}" ;;
		-l|--ligands) USER_LIGS="${!next:-}" ;;
		--ligands=*) USER_LIGS="${!i#*=}" ;;
		-o|--output) USER_OUT="${!next:-}" ;;
		--output=*) USER_OUT="${!i#*=}" ;;
	esac
done
if [ -n "$USER_OUT" ]; then OUT="$USER_OUT"; fi

mkdir -p "$OUT"

INPUT="${INPUT:-}"

# -l on the command line replaces INPUT/LIGS from the environment.
if [ -n "$USER_LIGS" ]; then
	LIGS="$USER_LIGS"
	if [ -f "$LIGS" ]; then INPUT="$LIGS"; else INPUT=""; fi
fi

if [ -n "${INPUT}" ] && [ ! -f "${INPUT}" ]; then
	echo "[!] INPUT is set but the file does not exist: $INPUT"
	exit 1
fi

if [ -n "${INPUT}" ] && [ -f "${INPUT}" ]; then
	echo  "[INFO] Single ligand mode detected: $INPUT"
	LIGAND_MODE="single"
	LIGS="$INPUT"
else
	echo  "[INFO] Batch mode: scanning $LIGS"
	LIGAND_MODE="batch"
fi

if [ -n "$USER_PDB" ]; then
	PDB="$USER_PDB"
elif [ -n "${PDB:-}" ]; then
	: # Use PDB from the environment.
elif [ -f "$WORK/protein.pdb" ]; then
	PDB="$WORK/protein.pdb"
elif [ -f "$WORK/protein_clean.pdb" ]; then
	PDB="$WORK/protein_clean.pdb"
else
	PDB="$WORK/protein.pdb"
fi

echo "======================================"
echo " AutoDocker Virtual Screening Pipeline"
echo "======================================"

# Check if required files exist
if [ ! -f "$PDB" ]; then
	echo "[!] Protein file not found: $PDB"
	echo "[!] Please place protein.pdb in workspace/ or set PDB=/path/to/file.pdb"
	exit 1
fi

# PubChem library modes (--library fda|custom) download their own ligands,
# so no local ligand files are required.
PUBCHEM_LIBRARY=0
for ((i = 1; i <= $#; i++)); do
	case "${!i}" in
		--library=fda|--library=custom) PUBCHEM_LIBRARY=1 ;;
		--library)
			next=$((i + 1))
			case "${!next:-}" in fda|custom) PUBCHEM_LIBRARY=1 ;; esac
			;;
	esac
done

if [ "$PUBCHEM_LIBRARY" = "1" ]; then
	mkdir -p "$LIGS"
	echo "[*] PubChem library mode: ligands will be downloaded"
elif [ "$LIGAND_MODE" = "single" ]; then
	if [ ! -f "$LIGS" ]; then
		echo "[!] Ligand file not found: $LIGS"
		exit 1
	fi

	ext="${LIGS##*.}"
	case "${ext,,}" in
		sdf)
			SDF_COUNT=1
			PDBQT_COUNT=0
			MOL2_COUNT=0
			PDB_COUNT=0
			;;
		pdbqt)
			SDF_COUNT=0
			PDBQT_COUNT=1
			MOL2_COUNT=0
			PDB_COUNT=0
			;;
		mol2)
			SDF_COUNT=0
			PDBQT_COUNT=0
			MOL2_COUNT=1
			PDB_COUNT=0
			;;
		pdb)
			SDF_COUNT=0
			PDBQT_COUNT=0
			MOL2_COUNT=0
			PDB_COUNT=1
			;;
		*)
			echo "[!] Unsupported ligand file type: $LIGS"
			exit 1
			;;
	esac
	TOTAL_LIGS=1
else
	if [ ! -d "$LIGS" ]; then
		echo "[!] Ligands directory not found: $LIGS"
		echo "[!] Please create workspace/ligs/ and add .sdf files"
		exit 1
	fi

	# Check for ligand files
	SDF_COUNT=$(find "$LIGS" -type f -iname "*.sdf" | wc -l)
	PDBQT_COUNT=$(find "$LIGS" -type f -iname "*.pdbqt" | wc -l)
	MOL2_COUNT=$(find "$LIGS" -type f -iname "*.mol2" | wc -l)
	PDB_COUNT=$(find "$LIGS" -type f -iname "*.pdb" | wc -l)

	TOTAL_LIGS=$((SDF_COUNT + PDBQT_COUNT + MOL2_COUNT + PDB_COUNT))
fi

if [ "$PUBCHEM_LIBRARY" != "1" ]; then
	if [ "$TOTAL_LIGS" -eq 0 ]; then
		echo "[!] No ligand files (.sdf/.pdbqt/.mol2/.pdb) found in $LIGS"
		exit  1
	fi

	echo "[*] Found $TOTAL_LIGS ligand(s)"
	echo "    - SDF: $SDF_COUNT"
	echo "    - PDBQT: $PDBQT_COUNT"
	echo "    - MOL2: $MOL2_COUNT"
	echo "    - PDB: $PDB_COUNT"
fi
echo "[*] Protein loaded: $(basename "$PDB")"

# Run the pipeline
# The container docks on all available CPUs unless PROCESSES (or MAX_WORKERS)
# says otherwise; runner.py itself defaults to one process.
PROCESSES="${PROCESSES:-${MAX_WORKERS:--1}}"
EXHAUSTIVENESS="${EXHAUSTIVENESS:-${EXHAUST:-}}"

CMD=(
	"$SCRIPT"
	-r "$PDB"
	-l "$LIGS"
	-o "$OUT"
	-p "$PROCESSES"
)

# Other settings are passed only when their variable is set, so the defaults
# live in one place (autodocker/vspipeline/defaults.py, see runner.py --help).
add_opt() {  # add_opt VALUE OPTION
	if [ -n "$1" ]; then
		CMD+=("$2" "$1")
	fi
}
add_opt "$EXHAUSTIVENESS" --exhaustiveness
add_opt "${BINDING_MODES:-}" --binding-modes
add_opt "${ENERGY_RANGE:-}" --energy-range
add_opt "${PADDING:-}" --padding
add_opt "${TOP_N:-}" --top-n

# Per-command timeout: only pass --timeout when TIMEOUT is set, so that
# runner.py's own default (VS_COMMAND_TIMEOUT, else 900 s) applies otherwise.
if [ -n "${TIMEOUT:-}" ]; then
	CMD+=(--timeout "$TIMEOUT")
fi

# Chain selection: only pass --chain when explicitly set, so runner.py's
# auto-detection (first chain) is the default. Set CHAIN=A (or A,B/all) to
# override.
if [ -n "${CHAIN:-}" ]; then
	CMD+=(--chain "$CHAIN")
fi

# The Lipinski rule-of-five filter is ON by default (matches runner.py); set
# LIPINSKI=0 to disable it. ADMET=0 is the deprecated name of the same switch.
if [ "${ADMET:-1}" = "0" ]; then
	echo "[!] ADMET=0 is deprecated; use LIPINSKI=0"
fi
if [ "${LIPINSKI:-1}" = "0" ] || [ "${ADMET:-1}" = "0" ]; then
	CMD+=(--no-lipinski)
fi

# Checkpoint resume is ON by default; set RESUME=0 to force a fresh run.
if [ "${RESUME:-1}" = "0" ]; then
	CMD+=(--no-resume)
fi

if [ -n "${POCKETS:-}" ]; then
	CMD+=(--pockets "$POCKETS")
fi

if [ "${NO_MINIMIZE:-0}" = "1" ]; then
	CMD+=(--no-minimize)
fi

# Exit code 3: outputs were written, but more ligands failed to dock than
# --max-failure-fraction allows; report the outputs, then exit with 3.
rc=0
"${CMD[@]}" "$@" || rc=$?
if [ "$rc" -ne 0 ] && [ "$rc" -ne 3 ]; then
	exit "$rc"
fi

echo "======================================"
if [ "$rc" -eq 3 ]; then
	echo " Pipeline completed, but too many ligands failed to dock (exit 3)"
else
	echo " Pipeline Completed Successfully"
fi
echo "======================================"
echo ""
echo "Results saved to: $OUT"
# List only what this run actually produced (e.g. pocket_summary.csv needs
# fpocket pockets, results_report.html is skipped by --no-html-report).
OUTPUTS=(
	"ranking.csv|Docking scores"
	"Top_hits.txt|Top-ranked hits"
	"Results_full.txt|Full text report"
	"metrics.txt|Pose metrics"
	"pocket_summary.csv|fpocket pocket scores"
	"ligand_metadata.json|ligand descriptors and source files"
	"run_manifest.json|versions, settings and checksums of this run"
	"consensus_ranking.csv|Vina and SMINA (Vinardo) scores"
	"consensus_summary.json|Vina/SMINA Spearman rho"
	"grid_box.pml|PyMOL grid visualization"
	"results_report.html|Self-contained HTML report"
	"docked/|Docked poses"
)
for entry in "${OUTPUTS[@]}"; do
	name="${entry%%|*}"
	if [ -e "$OUT/$name" ]; then
		echo "  - $name: ${entry#*|}"
	fi
done
exit "$rc"

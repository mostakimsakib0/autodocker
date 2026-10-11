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

# Invalid input exits with 2, like an invalid option of runner.py.
EXIT_USAGE=2
usage_error() {
	local line
	for line in "$@"; do echo "[!] $line"; done
	exit "$EXIT_USAGE"
}

# --help and --version need no workspace: hand them straight to the runner.
for arg in "$@"; do
	case "$arg" in
		--) break ;;
		-h | --help | --version) exec "$SCRIPT" "$@" ;;
	esac
done

LIGS="${LIGS:-/workspace/ligs}"
PDB="${PDB:-}"

# -r/-l/-o and --library given on the command line override the workspace
# defaults (argparse keeps the last value); validate and report what will be
# used. Recognised in every form runner.py accepts: "-r FILE", "-rFILE",
# "-r=FILE", "-vrFILE", "--receptor FILE" and "--receptor=FILE" (runner.py
# rejects abbreviated long options). A token starting with "-" is always an
# option to argparse, so option values cannot be mistaken for options here.
USER_PDB="" USER_LIGS="" USER_OUT="" LIBRARY=""
ARGS=("$@")
for ((i = 0; i < ${#ARGS[@]}; i++)); do
	arg="${ARGS[i]}"
	value="${ARGS[i + 1]:-}"
	case "$arg" in
		--) break ;;
		--receptor) USER_PDB="$value"; i=$((i + 1)) ;;
		--receptor=*) USER_PDB="${arg#*=}" ;;
		--ligands) USER_LIGS="$value"; i=$((i + 1)) ;;
		--ligands=*) USER_LIGS="${arg#*=}" ;;
		--output) USER_OUT="$value"; i=$((i + 1)) ;;
		--output=*) USER_OUT="${arg#*=}" ;;
		--library) LIBRARY="$value"; i=$((i + 1)) ;;
		--library=*) LIBRARY="${arg#*=}" ;;
		--*) ;;
		-?*)
			# Short options; -v takes no value and may precede another one.
			short="${arg#-}"
			while [ "${short:0:1}" = "v" ]; do short="${short:1}"; done
			case "$short" in
				r) USER_PDB="$value"; i=$((i + 1)) ;;
				r=*) USER_PDB="${short#r=}" ;;
				r?*) USER_PDB="${short#r}" ;;
				l) USER_LIGS="$value"; i=$((i + 1)) ;;
				l=*) USER_LIGS="${short#l=}" ;;
				l?*) USER_LIGS="${short#l}" ;;
				o) USER_OUT="$value"; i=$((i + 1)) ;;
				o=*) USER_OUT="${short#o=}" ;;
				o?*) USER_OUT="${short#o}" ;;
			esac
			;;
	esac
done
if [ -n "$USER_OUT" ]; then OUT="$USER_OUT"; fi

# PubChem library modes (--library fda|custom) download their own ligands,
# so no local ligand files are required and INPUT/LIGS/-l are not used.
PUBCHEM_LIBRARY=0
case "$LIBRARY" in fda | custom) PUBCHEM_LIBRARY=1 ;; esac

INPUT="${INPUT:-}"

if [ "$PUBCHEM_LIBRARY" = "1" ]; then
	if [ -n "$INPUT" ]; then
		echo "[!] INPUT=$INPUT is ignored with --library $LIBRARY (ligands come from PubChem)"
		INPUT=""
	fi
elif [ -n "$USER_LIGS" ]; then
	# -l on the command line replaces INPUT/LIGS from the environment.
	LIGS="$USER_LIGS"
	if [ -f "$LIGS" ]; then INPUT="$LIGS"; else INPUT=""; fi
fi

if [ -n "${INPUT}" ] && [ ! -f "${INPUT}" ]; then
	usage_error "INPUT is set but the file does not exist: $INPUT"
fi

if [ -n "${INPUT}" ]; then
	echo  "[INFO] Single ligand mode detected: $INPUT"
	LIGAND_MODE="single"
	LIGS="$INPUT"
elif [ "$PUBCHEM_LIBRARY" != "1" ]; then
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
	usage_error "Protein file not found: $PDB" \
		"Please place protein.pdb in workspace/ or set PDB=/path/to/file.pdb"
fi

if [ "$PUBCHEM_LIBRARY" = "1" ]; then
	echo "[*] PubChem library mode: ligands will be downloaded"
elif [ "$LIGAND_MODE" = "single" ]; then
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
			usage_error "Unsupported ligand file type: $LIGS"
			;;
	esac
	TOTAL_LIGS=1
else
	if [ ! -d "$LIGS" ]; then
		usage_error "Ligands directory not found: $LIGS" \
			"Please create workspace/ligs/ and add .sdf files"
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
		usage_error "No ligand files (.sdf/.pdbqt/.mol2/.pdb) found in $LIGS"
	fi

	echo "[*] Found $TOTAL_LIGS ligand(s)"
	echo "    - SDF: $SDF_COUNT"
	echo "    - PDBQT: $PDBQT_COUNT"
	echo "    - MOL2: $MOL2_COUNT"
	echo "    - PDB: $PDB_COUNT"
fi

if ! mkdir -p -- "$OUT"; then
	usage_error "Cannot create the output directory: $OUT"
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

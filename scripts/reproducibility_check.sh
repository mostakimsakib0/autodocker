#!/usr/bin/env bash
# Reproducibility check for the autodocker pipeline.
#
# Runs the pipeline twice on the SAME input with the SAME seed and verifies
# that the ranked results (ranking.csv) are byte-identical. This is the
# end-to-end determinism guarantee advertised in README.md and the JOSS
# paper.
#
# Usage:
#   scripts/reproducibility_check.sh <protein.pdb> <ligand_dir> [output_dir]
#
# Exit code 0 = deterministic, 1 = results differ, 2 = usage/missing tools.
set -euo pipefail

if [ "$#" -lt 2 ] || [ "$#" -gt 3 ]; then
    echo "Usage: $0 <protein.pdb> <ligand_dir> [output_dir]" >&2
    exit 2
fi

PDB="$1"
LIGS="$2"
ROOT="${3:-$(mktemp -d)}"
SEED="${VS_SEED:-42}"
PROCS="${PROCESSES:-4}"

for tool in vina obabel; do
    command -v "$tool" >/dev/null 2>&1 || {
        echo "[!] required binary not found on PATH: $tool" >&2
        exit 2
    }
done

RUN1="$ROOT/run1"
RUN2="$ROOT/run2"
mkdir -p "$RUN1" "$RUN2"

echo "[*] run #1 -> $RUN1"
python3 autodocker/runner.py \
    -r "$PDB" -l "$LIGS" -o "$RUN1" \
    --seed "$SEED" --exhaustiveness "${EXHAUSTIVENESS:-4}" \
    --binding-modes 1 --energy-range 1.0 --processes "$PROCS" \
    --no-resume >/dev/null

echo "[*] run #2 -> $RUN2"
python3 autodocker/runner.py \
    -r "$PDB" -l "$LIGS" -o "$RUN2" \
    --seed "$SEED" --exhaustiveness "${EXHAUSTIVENESS:-4}" \
    --binding-modes 1 --energy-range 1.0 --processes "$PROCS" \
    --no-resume >/dev/null

echo "[*] comparing ranking.csv ..."
if diff -q "$RUN1/ranking.csv" "$RUN2/ranking.csv" >/dev/null; then
    echo "[✔] DETERMINISTIC: both runs produced identical ranking.csv"
    exit 0
else
    echo "[✘] NON-DETERMINISTIC: ranking.csv differs between runs" >&2
    diff "$RUN1/ranking.csv" "$RUN2/ranking.csv" | head -40 >&2 || true
    exit 1
fi
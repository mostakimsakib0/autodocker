#!/usr/bin/env bash
# End-to-end reproducibility check for the AutoDocker pipeline.
#
# Runs the pipeline three times on the same input with the same seed:
#   run1  -p 1
#   run2  -p 1          (repeatability)
#   runN  -p PROCESSES  (worker count / Vina threads per process differ)
# and verifies that run2 and runN match run1:
#   * ranking.csv is byte-identical once the run directory itself is
#     normalised away (versions that write absolute Docked_File paths)
#   * the same set of docked pose files exists and each one is
#     byte-identical (docked/*_out.pdbqt, docked/*_smina.pdbqt)
#   * results_report.html is byte-identical (it holds no timestamp; the NGL
#     viewer it embeds must come from the same source in every run, e.g. the
#     image's /ngl/ngl.js)
#
# Usage:
#   scripts/reproducibility_check.sh <protein.pdb> <ligand_dir> [output_dir] [-- runner args...]
# Environment: VS_SEED (read by runner.py; its default otherwise), PROCESSES (4),
#              EXHAUSTIVENESS (4), PYTHON (python3)
#
# Runner arguments after "--" are passed to every run, e.g. a fixed box
# (--grid-center X,Y,Z --grid-size S,S,S) to keep pocket detection out of
# the comparison.
#
# Exit codes: 0 = all runs identical, 1 = results differ,
#             2 = usage error / missing tools, 3 = a pipeline run failed.
set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RUNNER="$SCRIPT_DIR/../autodocker/runner.py"
PYTHON="${PYTHON:-python3}"

usage() {
    echo "Usage: $0 <protein.pdb> <ligand_dir> [output_dir] [-- runner args...]" >&2
    exit 2
}

[ "$#" -ge 2 ] || usage
PDB="$1"; LIGS="$2"; shift 2
ROOT=""
if [ "$#" -gt 0 ] && [ "$1" != "--" ]; then ROOT="$1"; shift; fi
EXTRA=()
if [ "$#" -gt 0 ]; then
    [ "$1" = "--" ] || usage
    shift
    EXTRA=("$@")
fi
[ -n "$ROOT" ] || ROOT="$(mktemp -d)"

PROCS="${PROCESSES:-4}"
EXH="${EXHAUSTIVENESS:-4}"

[ -f "$PDB" ] || { echo "[!] receptor not found: $PDB" >&2; exit 2; }
[ -e "$LIGS" ] || { echo "[!] ligands not found: $LIGS" >&2; exit 2; }
[ -f "$RUNNER" ] || { echo "[!] runner not found: $RUNNER" >&2; exit 2; }
for tool in vina obabel "$PYTHON"; do
    command -v "$tool" >/dev/null 2>&1 || {
        echo "[!] required binary not found on PATH: $tool" >&2
        exit 2
    }
done

run_pipeline() {  # run_pipeline <outdir> <processes>
    local out="$1" procs="$2"
    rm -rf "$out"
    mkdir -p "$out"
    echo "[*] $(basename "$out"): -p $procs -> $out"
    if ! "$PYTHON" "$RUNNER" -r "$PDB" -l "$LIGS" -o "$out" \
            --exhaustiveness "$EXH" \
            --binding-modes 1 --energy-range 1.0 --processes "$procs" \
            --no-resume ${EXTRA[@]+"${EXTRA[@]}"} \
            >"$out.log" 2>&1; then
        echo "[!] pipeline run failed: $out (log: $out.log)" >&2
        tail -20 "$out.log" >&2
        exit 3
    fi
    [ -f "$out/ranking.csv" ] || { echo "[!] no ranking.csv in $out" >&2; exit 3; }
}

compare_runs() {  # compare_runs <reference> <other>; prints differences
    "$PYTHON" - "$1" "$2" <<'EOF'
import glob, os, sys
ref, other = sys.argv[1], sys.argv[2]
problems = []

def ranking(run):
    with open(os.path.join(run, "ranking.csv"), "rb") as f:
        data = f.read()
    for prefix in (os.path.abspath(run), run):  # absolute Docked_File paths
        data = data.replace(prefix.encode() + b"/", b"<RUN>/")
    return data

if ranking(ref) != ranking(other):
    problems.append("ranking.csv differs")

def report(run):
    path = os.path.join(run, "results_report.html")
    if not os.path.exists(path):
        return None
    with open(path, "rb") as f:
        return f.read()

if report(ref) is None and report(other) is None:
    pass  # --no-html-report given after "--"
elif report(ref) != report(other):
    problems.append("results_report.html differs")

def poses(run):
    files = []
    for pattern in ("*_out.pdbqt", "*_smina.pdbqt"):
        files += glob.glob(os.path.join(run, "docked", pattern))
    return {os.path.basename(p): p for p in files}

a, b = poses(ref), poses(other)
if not a:
    problems.append("no docked poses in reference run")
for name in sorted(set(a) ^ set(b)):
    problems.append(f"pose file only in one run: {name}")
for name in sorted(set(a) & set(b)):
    with open(a[name], "rb") as fa, open(b[name], "rb") as fb:
        if fa.read() != fb.read():
            problems.append(f"pose file differs: {name}")
for p in problems[:40]:
    print("    " + p)
sys.exit(1 if problems else 0)
EOF
}

run_pipeline "$ROOT/run1" 1
run_pipeline "$ROOT/run2" 1
run_pipeline "$ROOT/runN" "$PROCS"

status=0
for other in run2 runN; do
    echo "[*] comparing run1 vs $other (ranking.csv, docked poses, results_report.html) ..."
    if compare_runs "$ROOT/run1" "$ROOT/$other"; then
        echo "[ok] run1 and $other are identical"
    else
        echo "[x] run1 and $other differ" >&2
        status=1
    fi
done

if [ "$status" -eq 0 ]; then
    echo "DETERMINISTIC: ranking.csv, all docked poses and the HTML report identical (-p 1, -p 1, -p $PROCS)"
else
    echo "NON-DETERMINISTIC: see differences above" >&2
fi
exit "$status"

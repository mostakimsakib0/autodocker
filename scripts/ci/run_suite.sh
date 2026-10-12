#!/usr/bin/env bash
# One job of .github/workflows/benchmarks.yml: run a benchmark shard inside
# the benchmark image with the fixed protocol of BENCHMARKS.md. Each step's
# exit status is appended to $BENCH_OUT/exit_codes.txt ("step=status"); a
# failing step does not stop the following ones, so a shard always keeps the
# results it has. The script exits 1 if any step failed.
#
# Usage: scripts/ci/run_suite.sh SUITE [ARG]
#   prepare               fetch CASF (CASF_SET), prepare the DUD-E targets
#                         (DUDE_TARGETS, empty = skip), split the CASF IDs
#                         into SHARDS shards, write SHA256SUMS
#   casf SHARD            J1 re-docking of shard SHARD (e.g. 007); on the
#                         shard's IDs of the seeded sample also the
#                         crystal-start comparison and, with RUN_PARITY=1, J3
#   dude TARGET           J2 enrichment of one DUD-E target
#   baseline-redock       J5 baseline (Meeko + raw Vina) on the CASF sample
#   baseline-enrichment TARGET
#                         J5 baseline on one DUD-E target
#   easydock              J5 EasyDock on the CASF sample
#
# Environment: BENCH_DATA, BENCH_OUT, BENCH_TOOLS (host directories, see
# scripts/ci/in_image.sh), CASF_SET (full|sample), BENCH_IMAGE.
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"
: "${BENCH_DATA:?}" "${BENCH_OUT:?}"
export BENCH_DATA BENCH_OUT
CASF_SET="${CASF_SET:-full}"
case "$CASF_SET" in
	full) CASF_DIR=casf2016_full ;;
	sample) CASF_DIR=casf2016 ;;
	*) echo "CASF_SET must be full or sample" >&2; exit 2 ;;
esac
NPROC="$(nproc)"

# Benchmark protocol (BENCHMARKS.md J1-J5); fixed, not workflow inputs.
SEED=42
SAMPLE_SIZE=40
SAMPLE_SEED=42
REDOCK_EXHAUSTIVENESS=8
REDOCK_MODES=9
REDOCK_ENERGY_RANGE=3
ENRICH_EXHAUSTIVENESS=4
MAX_ACTIVES=50
MAX_DECOYS=500
# CASF-2016 core-set IDs left out of J1: the core-set SMILES puts an sp3
# stereocentre where the crystal carbon is planar (2xys, 4ih7: no RDKit start
# conformer with the crystal InChIKey; 3gr2, 3n86: same tautomer artefact),
# see BENCHMARKS.md. SKIPPED: peptide ligands fetch_casf_rcsb.py cannot build.
CASF_EXCLUDED="${CASF_EXCLUDED:-2xys,4ih7,3gr2,3n86}"
CASF_SKIPPED="${CASF_SKIPPED:-1a30,3bv9,3uri}"
# Wall-time limit of one DUD-E target's runner call (s); not a docking
# parameter. 4.5 h stays below the job's 6 h limit.
ENRICH_TIMEOUT=16200

mkdir -p "$BENCH_OUT"
codes="$BENCH_OUT/exit_codes.txt"
status=0

step() {
	local name="$1"
	shift
	echo "::group::$name"
	"$@"
	local rc=$?
	echo "::endgroup::"
	echo "$name=$rc" >> "$codes"
	if [ "$rc" -ne 0 ]; then
		echo "::error::$name failed (exit $rc)"
		status=1
	fi
	return 0
}

# shellcheck disable=SC2317  # called through step()
img() {
	"$HERE/in_image.sh" "$@"
}

join_ids() {
	paste -sd, "$1"
}

# The IDs of one file that are also in another (both sorted, one per line).
common_ids() {
	LC_ALL=C comm -12 "$1" "$2"
}

run_dirs() {
	# Host run directories -> container paths for verify_run.py.
	local host_root="$1" container_root="$2" d
	for d in "$host_root"/*/; do
		[ -d "$d" ] && printf '%s\n' "$container_root/$(basename "$d")"
	done
}

suite_prepare() {
	local manifest
	if [ "$CASF_SET" = full ]; then
		step fetch_casf img python3 /repo/scripts/fetch_casf_rcsb.py \
			--outdir "/data/$CASF_DIR" --sample 0
	else
		step fetch_casf img python3 /repo/scripts/fetch_casf_rcsb.py \
			--outdir "/data/$CASF_DIR" --sample "$SAMPLE_SIZE" --sample-seed "$SAMPLE_SEED"
	fi
	manifest="$BENCH_DATA/$CASF_DIR/manifest.csv"
	step split_casf python3 "$REPO/scripts/merge_shards.py" split \
		--set "$CASF_SET" --exclude "$CASF_EXCLUDED" --skipped "$CASF_SKIPPED" \
		--manifest "$manifest" --shards "${SHARDS:-20}" --outdir "$BENCH_DATA/casf_ids"
	local subdirs=("$CASF_DIR")
	if [ -n "${DUDE_TARGETS:-}" ]; then
		step prepare_dude img python3 /repo/scripts/benchmark_enrichment.py --prepare-only \
			--targets "$DUDE_TARGETS" --workdir /data/dude --outdir /data/dude_prepare \
			--box-source crystal --max-actives "$MAX_ACTIVES" --max-decoys "$MAX_DECOYS" \
			--seed "$SEED" --runner /autodocker/runner.py
		subdirs+=(dude dude_prepare)
	fi
	"$HERE/hash_inputs.sh" "$BENCH_DATA" "${subdirs[@]}" > "$BENCH_DATA/SHA256SUMS"
	echo "hash_inputs=$?" >> "$codes"
	wc -l < "$BENCH_DATA/SHA256SUMS"
}

suite_casf() {
	local shard="$1" ids sample_ids
	local ids_file="$BENCH_DATA/casf_ids/shard_$shard.txt"
	if [ ! -s "$ids_file" ]; then
		echo "::error::no IDs for shard $shard ($ids_file)"
		echo "ids=1" >> "$codes"
		return 1
	fi
	ids="$(join_ids "$ids_file")"
	common_ids "$ids_file" "$BENCH_DATA/casf_ids/sample.txt" > "$BENCH_OUT/sample_ids.txt"
	cp "$ids_file" "$BENCH_OUT/ids.txt"
	sample_ids="$(join_ids "$BENCH_OUT/sample_ids.txt")"
	local docking=(--runner /autodocker/runner.py --processes "$NPROC" --threads 1
		--seed "$SEED" --exhaustiveness "$REDOCK_EXHAUSTIVENESS"
		--binding-modes "$REDOCK_MODES" --energy-range "$REDOCK_ENERGY_RANGE"
		--receptor-hetatm metals)
	step redock img python3 /repo/scripts/merge_shards.py capture redock -- \
		--data-dir "/data/$CASF_DIR" --pdbids "$ids" --outdir /out/redock \
		--start-conformer rdkit --conformer-seed "$SEED" "${docking[@]}"
	if [ -n "$sample_ids" ]; then
		step redock_crystal img python3 /repo/scripts/merge_shards.py capture redock -- \
			--data-dir "/data/$CASF_DIR" --pdbids "$sample_ids" --outdir /out/redock_crystal \
			--start-conformer crystal "${docking[@]}"
		if [ "${RUN_PARITY:-0}" = 1 ]; then
			step parity img python3 /repo/scripts/benchmark_vina_parity.py \
				--workdir /out/redock/work --pdbids "$sample_ids" --outdir /out/parity
		fi
	fi
	local runs
	mapfile -t runs < <(run_dirs "$BENCH_OUT/redock/work" /out/redock/work
		run_dirs "$BENCH_OUT/redock_crystal/work" /out/redock_crystal/work)
	step verify img bash /repo/scripts/ci/verify_runs.sh /out/verify "${runs[@]}"
}

suite_dude() {
	local target="$1"
	step enrichment img python3 /repo/scripts/benchmark_enrichment.py \
		--targets "$target" --workdir /data/dude --outdir /out/enrichment \
		--box-source crystal --max-actives "$MAX_ACTIVES" --max-decoys "$MAX_DECOYS" \
		--exhaustiveness "$ENRICH_EXHAUSTIVENESS" --seed "$SEED" \
		--processes "$NPROC" --runner /autodocker/runner.py --timeout "$ENRICH_TIMEOUT"
	step verify img bash /repo/scripts/ci/verify_runs.sh /out/verify "/data/dude/$target/dock"
	# The runner's output directory lives next to the inputs; keep it.
	mkdir -p "$BENCH_OUT/dude_work/$target"
	cp -a "$BENCH_DATA/dude/$target/dock" "$BENCH_DATA/dude/$target/runner.log" \
		"$BENCH_OUT/dude_work/$target/" 2>/dev/null || true
}

baseline_tools() {
	: "${BENCH_TOOLS:?}"
	export BENCH_TOOLS
	step fetch_vina "$HERE/fetch_vina.sh" "$BENCH_TOOLS/vina"
	step env_baseline img bash /repo/scripts/ci/baseline_env.sh baseline /tools/baseline
	mkdir -p "$BENCH_OUT/env"
	cp "$BENCH_TOOLS"/*.freeze.txt "$BENCH_OUT/env/" 2>/dev/null || true
}

suite_baseline_redock() {
	baseline_tools
	step baseline_redock img python3 /repo/scripts/merge_shards.py capture baseline -- \
		redock --data-dir "/data/$CASF_DIR" --sample "$SAMPLE_SIZE" --sample-seed "$SAMPLE_SEED" \
		--start-conformer rdkit --conformer-seed "$SEED" \
		--meeko-bin /tools/baseline/bin --vina-bin /tools/vina \
		--outdir /out/baseline/redock --processes "$NPROC" --threads 1 --seed "$SEED" \
		--exhaustiveness "$REDOCK_EXHAUSTIVENESS" --binding-modes "$REDOCK_MODES" \
		--energy-range "$REDOCK_ENERGY_RANGE"
}

suite_baseline_enrichment() {
	local target="$1"
	baseline_tools
	step baseline_enrichment img python3 /repo/scripts/benchmark_baseline.py enrichment \
		--workdir /data/dude --targets "$target" --ligands-from /data/dude_prepare \
		--max-actives "$MAX_ACTIVES" --max-decoys "$MAX_DECOYS" \
		--meeko-bin /tools/baseline/bin --vina-bin /tools/vina \
		--outdir /out/baseline/enrichment --processes "$NPROC" --threads 1 --seed "$SEED" \
		--exhaustiveness "$ENRICH_EXHAUSTIVENESS"
}

suite_easydock() {
	baseline_tools
	step env_easydock img bash /repo/scripts/ci/baseline_env.sh easydock /tools/easydock
	cp "$BENCH_TOOLS"/*.freeze.txt "$BENCH_OUT/env/" 2>/dev/null || true
	step easydock_redock img python3 /repo/scripts/merge_shards.py capture baseline -- \
		redock --engine easydock --easydock-bin /tools/easydock/bin/easydock \
		--data-dir "/data/$CASF_DIR" --sample "$SAMPLE_SIZE" --sample-seed "$SAMPLE_SEED" \
		--start-conformer rdkit --conformer-seed "$SEED" --meeko-bin /tools/baseline/bin \
		--outdir /out/baseline/easydock_redock --processes "$NPROC" --threads 1 --seed "$SEED" \
		--exhaustiveness "$REDOCK_EXHAUSTIVENESS" --binding-modes "$REDOCK_MODES" \
		--energy-range "$REDOCK_ENERGY_RANGE"
}

suite="${1:?usage: $0 SUITE [ARG]}"
case "$suite" in
	prepare) suite_prepare ;;
	casf) suite_casf "${2:?shard}" || status=1 ;;
	dude) suite_dude "${2:?target}" ;;
	baseline-redock) suite_baseline_redock ;;
	baseline-enrichment) suite_baseline_enrichment "${2:?target}" ;;
	easydock) suite_easydock ;;
	*) echo "unknown suite: $suite" >&2; exit 2 ;;
esac
exit "$status"

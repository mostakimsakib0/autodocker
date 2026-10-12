#!/usr/bin/env bash
# Aggregate job of .github/workflows/benchmarks.yml: merge the downloaded
# shard artifacts with scripts/merge_shards.py into the layout of a single
# local run (the layout jcim_2.0/numbers/config.json reads as results_2.0/):
#
#   RESULTS/redock  redock_crystal  parity  enrichment
#   RESULTS/baseline/{redock,enrichment,easydock_redock}
#   RESULTS/merge/<name>.json      merge report per result (missing IDs, ...)
#   RESULTS/inputs/                SHA256SUMS and casf_ids/ of the prepare job
#   RESULTS/build-info.json        the image record of the build job
#   RESULTS/benchmarks_summary.json
#
# Usage: scripts/ci/aggregate.sh SHARDS RESULTS
#   SHARDS holds one directory per artifact (casf-shard-000, dude-akt1, ...,
#   build-info, inputs-manifest). SUITES (comma list), CASF_SET and the
#   GITHUB_* variables describe the run.
# Exit status: 0 when every requested result merged completely, else 1.
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MERGE="$HERE/../merge_shards.py"
shards="${1:?usage: $0 SHARDS RESULTS}"
results="${2:?usage: $0 SHARDS RESULTS}"
suites=",${SUITES:-casf,dude,baseline,easydock,parity},"
mkdir -p "$results/merge"
status=0

manifest="$shards/inputs-manifest"
ids="$manifest/casf_ids/ids.txt"
sample="$manifest/casf_ids/sample.txt"
plan="$manifest/casf_ids/plan.json"
selection=""
[ -f "$plan" ] && selection="$(jq -r .selection "$plan")"
targets="${DUDE_TARGETS:-akt1,braf,cdk2,egfr,fa10,hivpr,andr,aldr,hmdh,aces}"

wanted() {
	case "$suites" in *",$1,"*) return 0 ;; *) return 1 ;; esac
}

# merge NAME SUITE OUTDIR SUBDIR GLOB EXPECTED [EXTRA ARGS...]
merge() {
	local name="$1" suite="$2" outdir="$3" subdir="$4" glob="$5" expected="$6"
	shift 6
	local dirs=() d
	for d in "$shards"/$glob/"$subdir"; do
		[ -f "$d/summary.json" ] && dirs+=("$d")
	done
	local args=(--outdir "$results/$outdir" --report "$results/merge/$name.json"
		--name "$name" --allow-missing)
	# EXPECTED: an ID file (skipped when the prepare job left none) or a list.
	if [[ "$expected" == */* ]]; then
		if [ -f "$expected" ]; then
			args+=(--expected-ids "$expected")
		else
			echo "::warning::$name: $expected missing; completeness not checked"
		fi
	elif [ -n "$expected" ]; then
		args+=(--expected-ids "$expected")
	fi
	if [ "${#dirs[@]}" -eq 0 ]; then
		echo "::error::$name: no shard produced results"
		jq -n --arg name "$name" --arg suite "$suite" \
			'{name: $name, suite: $suite, complete: false, error: "no shard produced results"}' \
			> "$results/merge/$name.json"
		status=1
		return
	fi
	echo "::group::merge $name (${#dirs[@]} shard(s))"
	python3 "$MERGE" "$suite" --shards "${dirs[@]}" "${args[@]}" "$@"
	local rc=$?
	echo "::endgroup::"
	if [ "$rc" -ne 0 ]; then
		echo "::error::$name: merge failed (exit $rc)"
		status=1
	elif [ "$(jq -r '.complete' "$results/merge/$name.json")" = false ]; then
		echo "::warning::$name: incomplete, missing $(jq -r '.missing | join(",")' "$results/merge/$name.json")"
		status=1
	fi
}

if wanted casf; then
	sel=()
	[ "${CASF_SET:-full}" = sample ] && [ -n "$selection" ] && sel=(--selection "$selection")
	merge redock redock redock redock 'casf-shard-*' "$ids" "${sel[@]}"
	crystal_sel=()
	[ -n "$selection" ] && crystal_sel=(--selection "$selection")
	merge redock_crystal redock redock_crystal redock_crystal 'casf-shard-*' "$sample" \
		"${crystal_sel[@]}"
fi
if wanted parity; then
	merge parity parity parity parity 'casf-shard-*' "$sample"
fi
if wanted dude; then
	merge enrichment enrichment enrichment enrichment 'dude-*' "$targets" --targets "$targets"
fi
if wanted baseline; then
	merge baseline_redock baseline-redock baseline/redock baseline/redock \
		'baseline-redock' "$sample"
	merge baseline_enrichment baseline-enrichment baseline/enrichment baseline/enrichment \
		'baseline-enrichment-*' "$targets" --targets "$targets"
fi
if wanted easydock; then
	merge easydock_redock baseline-redock baseline/easydock_redock \
		baseline/easydock_redock 'easydock-redock' "$sample"
fi

if [ -d "$manifest" ]; then
	mkdir -p "$results/inputs"
	cp -a "$manifest/." "$results/inputs/"
fi
[ -f "$shards/build-info/build-info.json" ] && cp "$shards/build-info/build-info.json" "$results/"

reports=()
for r in "$results"/merge/*.json; do
	[ -f "$r" ] && reports+=("$r")
done
python3 "$MERGE" report --out "$results/benchmarks_summary.json" --shards-root "$shards" \
	--merge-report "${reports[@]}" \
	--build-info "$shards/build-info/build-info.json" \
	--inputs-sums "$manifest/SHA256SUMS" --plan "$plan" \
	--meta "git_sha=${GITHUB_SHA:-}" --meta "git_ref=${GITHUB_REF:-}" \
	--meta "repository=${GITHUB_REPOSITORY:-}" --meta "run_id=${GITHUB_RUN_ID:-}" \
	--meta "run_attempt=${GITHUB_RUN_ATTEMPT:-}" --meta "workflow=${GITHUB_WORKFLOW:-}" \
	--meta "suites=${SUITES:-}" --meta "casf_set=${CASF_SET:-}" \
	--meta "shards=${SHARDS:-}" || status=1
exit "$status"

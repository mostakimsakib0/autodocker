#!/usr/bin/env bash
# Record of the machine a benchmark job ran on (shard-info.json), collected
# by `merge_shards.py report` into benchmarks_summary.json.
#
# Usage: scripts/ci/shard_info.sh OUT.json NAME [KEY=VALUE...]
#   STARTED_AT   job start (UTC, ISO 8601), recorded as started_at
#   KEY=VALUE    extra string fields, e.g. exit codes (redock_exit=0)
set -euo pipefail

out="${1:?usage: $0 OUT.json NAME [KEY=VALUE...]}"
name="${2:?usage: $0 OUT.json NAME [KEY=VALUE...]}"
shift 2

cpu_model="$(awk -F': *' '/^model name/ { print $2; exit }' /proc/cpuinfo 2>/dev/null || true)"
mem_kb="$(awk '/^MemTotal:/ { print $2 }' /proc/meminfo 2>/dev/null || echo 0)"
extra='{}'
for kv in "$@"; do
	extra="$(jq -c --arg k "${kv%%=*}" --arg v "${kv#*=}" '. + {($k): $v}' <<< "$extra")"
done

mkdir -p "$(dirname "$out")"
jq -n \
	--arg name "$name" \
	--arg cpu_model "$cpu_model" \
	--argjson nproc "$(nproc)" \
	--argjson mem_kb "${mem_kb:-0}" \
	--arg kernel "$(uname -r)" \
	--arg arch "$(uname -m)" \
	--arg runner_name "${RUNNER_NAME:-}" \
	--arg runner_os "${RUNNER_OS:-}" \
	--arg image_os "${ImageOS:-}" \
	--arg image_version "${ImageVersion:-}" \
	--arg job "${GITHUB_JOB:-}" \
	--arg started_at "${STARTED_AT:-}" \
	--arg finished_at "$(date -u +%Y-%m-%dT%H:%M:%SZ)" \
	--argjson extra "$extra" \
	'{name: $name, job: $job, cpu_model: $cpu_model, nproc: $nproc, mem_kb: $mem_kb,
	  kernel: $kernel, arch: $arch, runner_name: $runner_name, runner_os: $runner_os,
	  runner_image: {os: $image_os, version: $image_version},
	  started_at: $started_at, finished_at: $finished_at} + $extra' > "$out"
cat "$out"

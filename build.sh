#!/bin/sh
cd "$(realpath "$(dirname "$0")")"
sed -nE 's/^FROM .+ AS\s+([^_]+)_builder$/Tool: \1/p' Dockerfile
while read -rp 'Tool? ' tool; do
	podman build --target "${tool}_builder" -t "${PWD##*/}:$tool-builder" .
done

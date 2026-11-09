#!/bin/sh
cd "$(realpath "$(dirname "$0")")"
sed -nE 's/^FROM .+ AS\s+([^_]+)(_builder)$/\1\2 \1-builder/p' Dockerfile |
	while read -r stage tag; do
		podman build --target "$stage" -t "${PWD##*/}:$tag" .
	done

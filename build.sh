#!/bin/sh
cd "$(realpath "$(dirname "$0")")" || exit
list()
{
	sed -nE 's/^FROM .+ AS\s+([^_]+)_builder$/Tool: \1/p' Dockerfile
}

build()
{
	podman build --target "${1}_builder" -t "${PWD##*/}:$1-builder" .
}

list
while read -rp 'Tool? ' tool; do
	case "$tool" in
		all)
			list | while read -r _ all; do build "$all"; done
			podman build -t "${PWD##*/}:latest" .
			break
			;;
		*) build "$tool" ;;
	esac
	list
done

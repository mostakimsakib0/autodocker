#!/bin/sh
die()
{
	echo >&2 "ERROR: $1"
	exit 1
}

try_patch()
{
	tree="$1"
	patch="$2"
	shift 2

	if patch --dry-run -d "$tree" "$@" < "$patch"; then
		patch -d "$tree" "$@" < "$patch"
	fi
}

test -d "$1/repo" || die "Target repo for patching does not exist"
test -d "$1/patches" || die "Patch tree does not exist"

for patch in "$1/patches"/*.patch; do
	! test "$1/patches/*.patch" = "$patch" || {
		echo >&2 "Empty patch tree"
		break
	}

	try_patch "$1/repo" "$patch" -p1
done

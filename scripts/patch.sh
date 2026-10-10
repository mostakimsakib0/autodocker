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

	if patch --dry-run -d "$tree" "$@" < "$patch" >/dev/null; then
		patch -d "$tree" "$@" < "$patch"
	elif patch --dry-run -R -d "$tree" "$@" < "$patch" >/dev/null; then
		echo >&2 "Already applied: $patch"
	else
		# A patch that no longer applies (e.g. after a submodule bump) must
		# stop the build instead of silently shipping unpatched code.
		die "Patch does not apply: $patch"
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

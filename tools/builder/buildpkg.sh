#!/bin/sh
die()
{
	echo >&2 "ERROR: $1"
	exit 1
}

addbin()
{
	install -TDm755 "$1" "$PKG/$2"
}

addfile()
{
	install -TDm644 "$1" "$PKG/$2"
}

addtree()
{
	cp -av --no-preserve=ownership "$1/" "$PKG/$2/"
}

set -eu

cd "$1" || die "Source directory does not exist"
test -f pkg.sh || die "Package script does not exist"

readonly PKG=/pkg
export DESTDIR="$PKG"

if test -f "$PWD/patches/series"; then
	test -d "$PWD/repo" || die "Patch destination does not exist"
	(
		export QUILT_PATCHES="$PWD/patches"
		cd "$PWD/repo"
		set -x
		if test -n "$(quilt unapplied)"; then
			quilt push -a -q
		fi
	) || die "Failed to apply patch series"
elif test -d "$PWD/patches"; then
	die "Patch tree exists, but not a quilt series"
fi

# shellcheck disable=SC1091
. "$PWD/pkg.sh"
(
	set -x
	build
) || die "Build failed..."

(
	set -x
	package
) || die "Package failed..."

/scripts/shlibs.sh "$1"

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

PKGVER='#default'

set -eu
test 1 -eq "$#" || die "Incorrect number of arguments\nUsage: $0 <buildtree>"
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

pkgver=$(
	{
		if ! test -d repo/.git; then
			exit
		fi

		cd repo || exit
		git describe --tags --long 2> /dev/null ||
			printf "r%s.g%s" "$(git rev-list --count HEAD)" "$(git rev-parse --short HEAD)"
	} | sed -E -e "$PKGVER" -e 's/([^-]*-g)/r\1/;s/-/./g' -e 's/^[^r0-9]*//g' -e '/^r\.$/d'
)

if test -n "$pkgver"; then
	test -d "/deps/$1" || {
		rm -f "/deps/$1"
		mkdir -p "/deps/$1"
	}
	echo "$1: $pkgver" >&2
	echo "$pkgver" > "/deps/$1/version"
fi

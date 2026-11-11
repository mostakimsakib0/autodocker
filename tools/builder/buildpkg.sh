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
cd "$1" || die "source dir not present"
test -f pkg.sh || die "pkg script does not exist"

# shellcheck disable=SC2034
readonly PKG=/pkg
export DESTDIR="$PKG"
# shellcheck disable=SC1091
. "$PWD/pkg.sh"
(
	set -x
	build
)
(
	set -x
	package
)

/scripts/shlibs.sh "$1"

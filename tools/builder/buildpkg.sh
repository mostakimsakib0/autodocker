#!/bin/sh
die()
{
	echo >&2 "ERROR: $1"
	exit 1
}

set -eu
cd "$1" || die "source dir not present"
test -f pkg.sh || die "pkg script does not exist"

# shellcheck disable=SC2034
readonly PKG=/pkg
# shellcheck disable=SC1091
(. ./pkg.sh)

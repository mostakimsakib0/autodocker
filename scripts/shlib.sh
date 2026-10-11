#!/bin/bash
# Usage: shlib.sh FILTER BINARY...
#
# Print one "run <package>" line (apt.sh deps format) for every Debian package
# that ships a shared library whose name matches FILTER (an ERE fragment, e.g.
# "boost") and that any BINARY links against. Used at build time to give the
# runtime image exactly the libraries the compiled tools need.
#
# Any failure (ldd, realpath, a library no package owns, or no match at all)
# exits non-zero, so the image build stops instead of shipping a binary
# without its libraries.
set -euo pipefail

if [ "$#" -lt 2 ]; then
	echo "usage: $0 FILTER BINARY..." >&2
	exit 2
fi

filter="$1"
shift

libs="$(ldd "$@" |
	sed -nE '/\s*'"$filter"'.+ => /s/^.+ => ([^[:blank:]]+) \(0x[0-9a-f]+\)/\1/p' |
	sort -u)"

if [ -z "$libs" ]; then
	echo "$0: no shared library matching '$filter' found in: $*" >&2
	exit 1
fi

# Resolve symlinks (soname -> real file) so dpkg-query finds the owner.
resolved="$(printf '%s\n' "$libs" | xargs realpath -e)"
mapfile -t paths <<< "$resolved"

owners="$(dpkg-query -S "${paths[@]}")"
printf '%s\n' "$owners" | sed -E 's/: .+$//' | sort -u | sed 's/^/run /'

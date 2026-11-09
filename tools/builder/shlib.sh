#!/bin/sh
filter="$1"
shift
ldd "$@" |
	sed -nE '/\s*'"$filter"'.+ => /s/^.+ => ([^[:blank:]]+) \(0x[0-9a-f]+\)/\1/p' |
	sort -u |
	xargs realpath -q |
	xargs dpkg-query -S 2> /dev/null |
	sed -E 's/: .+$//' |
	sed 's/^/run /'

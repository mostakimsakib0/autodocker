#!/bin/bash
set -euo pipefail

filterelf()
{
	for f; do
		magic=$(od -An -N4 -tx1 "$f" 2> /dev/null | tr -d " \n")
		[ "$magic" = "7f454c46" ] && printf "%s\n" "$f"
	done
}
export -f filterelf

lsdeps()
{
	find /pkg -type f -exec bash -c 'filterelf "$@"' -- {} + | while read -r file; do
		printf 'shlib-deps: %s\n' "${file#/pkg/}" >&2
		ldd "$file" 2> /dev/null
	done \
		| sed -nE '/\s*.+ => /s/^.+ => ([^[:blank:]]+) \(0x[0-9a-f]+\)/\1/p' \
		| sort -u \
		| xargs -r realpath -q \
		| xargs -r dpkg-query -S 2> /dev/null \
		| sed -E 's/: .+$//' \
		| sed -E 's/:[a-z0-9]+$//' \
		| sort -u
}

lsbase()
{
	dpkg-query -W -f='${binary:Package} ${Priority} ${Essential}\n' | while read -r pkg prio req; do
		case "$prio:$req" in
			required:*) ;;
			important:*) ;;
			*:yes) ;;
			*:*) continue ;;
		esac
		echo "$pkg"
	done | sort -u
}

pkgdeps()
{
	apt-cache depends --installed "$1" \
		| sed -nE 's/^[[:space:]]*(PreDepends|Depends): ([^ ]+)$/\2/p'
}

declare -A candidates=()
declare -A referenced=()

# Filter out system deps
while IFS= read -r pkg; do
	candidates["$pkg"]=1
done < <(
	comm -23 \
		<(lsdeps) \
		<(lsbase)
)

# Mark transitive dependencies
for pkg in "${!candidates[@]}"; do
	while IFS= read -r dep; do
		((${candidates[$dep]:-0})) || continue
		referenced["$dep"]=1
	done < <(pkgdeps "$pkg")
done

if (($# == 0)); then
	set -- /dev/stdout
else
	mkdir -p "/deps/$1"
	set -- "/deps/$1/deps" "/src/$1/deps"
fi

# Collect dependencies referenced by graph root
{
	(($# < 2)) || [[ ! -f $2 ]] || sed -nE '/^run\s+/p' "$2"
	for pkg in "${!candidates[@]}"; do
		((${referenced[$pkg]:-0})) || echo "run $pkg"
	done
} | sort -u > "$1"

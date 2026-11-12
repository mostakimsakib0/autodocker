#!/bin/sh
cd -P "$(dirname "$0")" || exit 1
cd "$(git rev-parse --show-toplevel)" || exit 1

#set -x
git submodule status --recursive |
	sed -nE 's/^\s*[[:punct:]]?[0-9a-f]+\s+([^[:blank:]]+)(\s+.+)?$/\1/p' |
	tac | while read -r path; do
	test -f "$path/.git" || continue
	gitdir=$(git -C "${path}" rev-parse --absolute-git-dir)
	test -d "$gitdir" || continue

	rm -f "$path/.git"
	mv -v "$gitdir" "$path/.git"
	git --work-tree="$path" --git-dir="$path/.git" config --local --path --unset core.worktree
done

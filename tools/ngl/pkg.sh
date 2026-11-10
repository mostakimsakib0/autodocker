build()
{
	pnpm install
}

package()
{
	for ext in js js.map; do
		for file in ngl ngl.umd ngl.esm; do
			addfile "node_modules/ngl/dist/$file.$ext" "ngl/$file.$ext"
		done
	done
}

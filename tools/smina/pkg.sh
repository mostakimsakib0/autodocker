build()
{
	cmake -S repo -B build
	cmake --build build -j"$NPROC"

	strip build/fromsmina build/tosmina build/smina
}

package()
{
	for bin in fromsmina tosmina smina; do
		addbin "build/$bin" "usr/local/bin/$bin"
	done
}

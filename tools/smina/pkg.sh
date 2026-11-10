build()
{
	cmake -S repo -B build
	make -j"$(nproc)" -C build
	strip build/fromsmina build/tosmina build/smina
}

package()
{
	for bin in fromsmina tosmina smina; do
		addbin "build/$bin" "usr/local/bin/$bin"
	done
}

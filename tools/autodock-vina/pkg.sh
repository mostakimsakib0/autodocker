build()
{
	make -j"$NPROC" -C repo/build/linux/release
	strip repo/build/linux/release/vina
}

package()
{
	addbin repo/build/linux/release/vina usr/local/bin/vina
}

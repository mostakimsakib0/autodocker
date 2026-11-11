build()
{
	cmake -B build -S repo \
		-DCMAKE_INSTALL_PREFIX=/usr/local \
		-DCMAKE_SKIP_RPATH=TRUE

	cmake --build build -j $(nproc)
}

package()
{
	cmake --install build
}

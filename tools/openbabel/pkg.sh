PKGVER='s/openbabel-//;s/-/./;s/\.([0-9]+)\.(g[0-9a-f]+)$/-\1-\2/'

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

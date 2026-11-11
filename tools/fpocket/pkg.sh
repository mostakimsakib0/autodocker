build()
{
	make -j"$(nproc)" -C repo
	strip repo/bin/*
}

package()
{
	addbin repo/bin/fpocket usr/local/bin/fpocket
	addbin repo/bin/dpocket usr/local/bin/dpocket
	addbin repo/bin/tpocket usr/local/bin/tpocket
	addbin repo/bin/mdpocket usr/local/bin/mdpocket
	addtree repo/plugins/LINUXAMD64/molfile usr/local/lib
}

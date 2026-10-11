/*
 * Stand-in for VMD's libmolfile_plugin.a on architectures for which fpocket
 * ships no prebuilt molfile plugins (fpocket only bundles LINUXAMD64 among
 * the Linux 64-bit targets; the VMD plugin sources need a registered download
 * and cannot be fetched reproducibly at build time).
 *
 * fpocket uses the plugins only for mmCIF input (.cif), topology files
 * (-y/--topology_file) and trajectories (mdpocket). PDB input and PDB output,
 * which is all AutoDocker uses, never call them. Any use of a plugin here
 * stops the program with a clear message instead of misbehaving.
 *
 * Built by the Dockerfile's fpocket stage into
 * plugins/<ARCH>/molfile/libmolfile_plugin.a when ARCH != LINUXAMD64.
 */
#include <stdio.h>
#include <stdlib.h>

#include "vmdplugin.h"

static void unavailable(const char *plugin)
{
	fprintf(stderr,
		"fpocket: the VMD molfile plugin '%s' is not available on this "
		"architecture; only PDB input is supported (no mmCIF, topology or "
		"trajectory files).\n", plugin);
	exit(EXIT_FAILURE);
}

#define MOLFILE_STUB(name)                                                   \
	int molfile_##name##_init(void) { return VMDPLUGIN_SUCCESS; }        \
	int molfile_##name##_register(void *v, vmdplugin_register_cb cb)     \
	{                                                                    \
		(void) v;                                                    \
		(void) cb;                                                   \
		unavailable(#name);                                          \
		return VMDPLUGIN_ERROR;                                      \
	}                                                                    \
	int molfile_##name##_fini(void) { return VMDPLUGIN_SUCCESS; }

MOLFILE_STUB(crdplugin)
MOLFILE_STUB(dcdplugin)
MOLFILE_STUB(dtrplugin)
MOLFILE_STUB(gromacsplugin)
MOLFILE_STUB(lammpsplugin)
MOLFILE_STUB(netcdfplugin)
MOLFILE_STUB(parm7plugin)
MOLFILE_STUB(pdbxplugin)
MOLFILE_STUB(pdbplugin)
MOLFILE_STUB(psfplugin)
MOLFILE_STUB(xyzplugin)

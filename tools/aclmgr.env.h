#ifndef __ACLMGR_ENV_H__
#define __ACLMGR_ENV_H__

static const char* const forwarded_env[] = {
	"ADMET", "AUTODOCKER_VERSION", "BINDING_MODES", "CHAIN", "ENERGY_RANGE",
	"EXHAUST", "EXHAUSTIVENESS", "INPUT", "LIGAND_MODE", "LIGS", "LIPINSKI",
	"MAX_WORKERS", "NO_MINIMIZE", "OMP_NUM_THREADS", "PADDING", "PDB",
	"POCKETS", "PROCESSES", "PUBCHEM_REST", "RESUME", "TIMEOUT", "TOP_N",
	/* Open Babel's plugin and data directories (Dockerfile ENV; the same
	 * as its compiled-in defaults) */
	"BABEL_DATADIR", "BABEL_LIBDIR",
	/* proxies, so PubChem library modes work behind one */
	"HTTP_PROXY", "HTTPS_PROXY", "NO_PROXY", "http_proxy", "https_proxy",
	"no_proxy", NULL
};

__attribute__((always_inline))
static inline bool is_forwarded(const char* entry)
{
	const char* eq = strchr(entry, '=');
	if (!eq) return false;
	size_t name_len = (size_t) (eq - entry);

	if (name_len > 3 && strncmp(entry, "VS_", 3) == 0) return true;
	for (size_t i = 0; forwarded_env[i] != NULL; i++)
		if (strlen(forwarded_env[i]) == name_len
		    && strncmp(entry, forwarded_env[i], name_len) == 0)
			return true;
	return false;
}


#endif /* __ACLMGR_ENV_H__ */

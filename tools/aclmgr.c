#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <pwd.h>
#include <sys/wait.h>
#include <unistd.h>
#include <stdbool.h>

#define  DROP_USER  "dockuser"
#define  DROP_GROUP "dockuser"

extern char **environ;

/* setpriv --reset-env clears the environment for the unprivileged child, so
 * the pipeline's own settings (README "Environment variables", every variable
 * entry.sh reads, the variables runner.py records in run_manifest.json, plus
 * any VS_* knob) are forwarded explicitly through /usr/bin/env.
 * tests/test_container_env.py checks this list against entry.sh and README. */
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

static bool is_forwarded(const char* entry)
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

bool is_rootless = false;
/* Run as capability-less container root instead of dockuser: set when the
 * container is rootless and /workspace is owned by container root, i.e. a
 * bind mount of the invoking user's own directory. */
bool run_as_root = false;
volatile sig_atomic_t caught = 0;
/* Owner of /workspace (stat in main). Under rootful Docker or Podman the
 * pipeline runs as dockuser, which must be able to write that bind mount. */
bool have_workspace = false;
uid_t workspace_uid = 0;
gid_t workspace_gid = 0;

static const char* const drop_caps_runner[]  = {
	"/usr/bin/setpriv",
	"--reuid", DROP_USER,
	"--regid", DROP_GROUP,
	"--init-groups", "--nnp",
	"--inh-caps=-all",
	"--ambient-caps=-all",
	"--bounding-set=-all",
	"--reset-env", "--"
};
static const size_t drop_caps_len = sizeof(drop_caps_runner) / sizeof(drop_caps_runner[0]);

/* Rootless podman: container root is the invoking host user, without any
 * privilege on the host. Staying uid 0 (with every capability dropped) makes
 * the pipeline write files the host user owns, so /workspace never has to be
 * chown'ed to dockuser and back -- a step that left the user's files owned by
 * a sub-uid whenever the container was killed. */
static const char* const drop_caps_rootless[]  = {
	"/usr/bin/setpriv",
	"--nnp",
	"--inh-caps=-all",
	"--ambient-caps=-all",
	"--bounding-set=-all",
	"--reset-env", "--"
};
static const size_t drop_caps_rootless_len =
	sizeof(drop_caps_rootless) / sizeof(drop_caps_rootless[0]);

int direct_run(const char* args[])
{
	pid_t pid = fork();  // Create a new child process
	if (pid < 0) {
		perror("fork failed");
		return -1;
	}

	if (pid == 0) {
		execv(args[0], (char* const*) args);
		_exit(127);
	}
	return pid;
}

int dropcap_run(const char* args[], size_t length)
{
	const char* const* prefix = run_as_root ? drop_caps_rootless : drop_caps_runner;
	size_t prefix_len = run_as_root ? drop_caps_rootless_len : drop_caps_len;
	const char* dropcap_args[prefix_len + length + 1];

	for(size_t i = 0; i < prefix_len; i++)
		dropcap_args[i] = prefix[i];
	for(size_t i = 0; i < length; i++)
		dropcap_args[prefix_len + i] = args[i];
	dropcap_args[prefix_len + length] = NULL;

	return direct_run(dropcap_args);
}

int run(bool drop_caps, const char* args[])
{

	size_t i = 0;
	printf("init$ ");

	if (drop_caps) printf("dropcaps ");

	do {
		printf("%s%c", args[i], (args[i + 1] == NULL)?'\n':' ');
	} while(args[++ i] != NULL);

	pid_t pid = -1;

	if (drop_caps) pid = dropcap_run(args, i);
	else pid = direct_run(args);

	if (pid == -1) return 1;

	int status;
	if (waitpid(pid, &status, 0) < 0)
		return 1;

	if (WIFEXITED(status))
		return WEXITSTATUS(status);

	if (WIFSIGNALED(status))
		return 128 + WTERMSIG(status);

	return 0;
}


static void handle_sigint(int sig)
{
	if (caught) return;
	caught = 1;

	printf("Received signal %d, forwarding...\n", sig);

	kill(0, sig);
	while (waitpid(-1, NULL, 0) > 0) {}
}


/* Rootful engine: give dockuser the UID (and, unless it is root's group, the
 * GID) of the owner of /workspace, in this container's /etc/passwd only. The
 * host user's directory is then writable without a chmod, and every output
 * belongs to that user. dockuser keeps its name, its home and the same
 * privilege drop (setpriv, no capabilities); only the numbers change. A
 * root-owned /workspace (e.g. a directory Docker created) or a failed
 * usermod leaves dockuser as it is. */
static void adopt_workspace_owner(void)
{
	if (!have_workspace || workspace_uid == 0) return;
	struct passwd* pw = getpwnam(DROP_USER);
	if (!pw) return;

	char uid[24], gid[24];
	snprintf(uid, sizeof uid, "%lu", (unsigned long) workspace_uid);
	snprintf(gid, sizeof gid, "%lu", (unsigned long) workspace_gid);
	if (workspace_gid != 0 && pw->pw_gid != workspace_gid) {
		const char* groupmod[] = {"/usr/sbin/groupmod", "-o", "-g", gid, DROP_GROUP, NULL};
		if (run(false, groupmod) != 0)
			fprintf(stderr, "warning: could not give %s the group of /workspace\n", DROP_GROUP);
	}
	if (pw->pw_uid != workspace_uid) {
		const char* usermod[] = {"/usr/sbin/usermod", "-o", "-u", uid, DROP_USER, NULL};
		if (run(false, usermod) != 0)
			fprintf(stderr, "warning: could not give %s the owner of /workspace; "
			        "outputs may not be writable\n", DROP_USER);
	}
}


int run_child_with_args(int argc, char** argv)
{

	size_t env_count = 0;
	for (char** e = environ; *e != NULL; e++)
		if (is_forwarded(*e)) env_count++;

	/* /usr/bin/env NAME=value... /autodocker/entry.sh args... NULL */
	const char* cmdline[1 + env_count + argc + 1];
	size_t n = 0;
	cmdline[n++] = "/usr/bin/env";
	for (char** e = environ; *e != NULL; e++)
		if (is_forwarded(*e)) cmdline[n++] = *e;
	cmdline[n++] = "/autodocker/entry.sh";
	for(int i = 1; i < argc; i++)
		cmdline[n++] = argv[i];
	cmdline[n] = NULL;

	if (run_as_root) {
		/* Hand back to the user any files an older image left owned by
		 * dockuser's sub-uid (it chown'ed /workspace and could be killed
		 * before restoring it). Container root is the invoking user, so
		 * this only ever gives files to their owner. */
		static const char* repair[] = {"/usr/bin/chown", "-R", "0:0", "/workspace", NULL};
		run(false, repair);
	} else if (!is_rootless) {
		adopt_workspace_owner();
	}
	return run(true, cmdline);
}

int main(int argc, char** argv)
{
	signal(SIGINT, handle_sigint);
	signal(SIGTERM, handle_sigint);

	setvbuf(stdout, NULL, _IONBF, 0);  // Disable buffering for stdout
	setvbuf(stderr, NULL, _IONBF, 0);  // Disable buffering for stderr
	setvbuf(stdin, NULL, _IONBF, 0);   // Disable buffering for stdin

	if (access("/proc/self/uid_map", F_OK) == 0) {
		FILE *umf = fopen("/proc/self/uid_map", "r");
		long long int iuid = -1, euid = -1, length = -1;

		if(umf) {
			while (fscanf(umf, "%lld %lld %lld", &iuid, &euid, &length) == 3)
				if (iuid == 0 && euid >= 1000 && euid <= 60000 && length == 1) {
					is_rootless = true;
					break;
				}

			fclose(umf);
		}
	}

	struct stat ws;
	if (stat("/workspace", &ws) == 0) {
		have_workspace = true;
		workspace_uid = ws.st_uid;
		workspace_gid = ws.st_gid;
	}
	run_as_root = is_rootless && have_workspace && workspace_uid == 0;

	/* Debug shell as root, only on an explicit first argument: "-r" is
	 * the pipeline's --receptor flag and must reach the runner. */
	static char* const debug_shell[] = {"/usr/bin/bash", "-i", NULL};
	if (argc > 1 && strcmp(argv[1], "--root-shell") == 0)
		execv("/usr/bin/bash", debug_shell);

	int ret = run_child_with_args(argc, argv);
	printf(ret == 0 ? "OK.\n" : "ERROR.\n");
	/* Propagate the pipeline's exit status so callers can detect failures. */
	return ret;
}

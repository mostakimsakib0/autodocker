#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/wait.h>
#include <unistd.h>
#include <stdbool.h>

#define  DROP_USER  "dockuser"
#define  DROP_GROUP "dockuser"

extern char **environ;

/* setpriv --reset-env clears the environment for the unprivileged child, so
 * the pipeline's own settings (README "Environment variables", plus any VS_*
 * knob) are forwarded explicitly through /usr/bin/env. */
static const char* const forwarded_env[] = {
	"ADMET", "BINDING_MODES", "CHAIN", "ENERGY_RANGE", "EXHAUST",
	"EXHAUSTIVENESS", "INPUT", "LIGAND_MODE", "LIGS", "MAX_WORKERS",
	"NO_MINIMIZE", "PADDING", "PDB", "POCKETS", "PROCESSES", "PUBCHEM_REST",
	"RESUME", "TIMEOUT", "TOP_N", NULL
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
volatile sig_atomic_t caught = 0;

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

int direct_run(const char* args[])
{
	pid_t pid = fork();  // Create a new child process
	if (pid < 0) {
		perror("fork failed");
		return -1;
	}

	if (pid == 0) {
		if (execv(args[0], (char* const*) args))
			_exit(127);
	} else {
		return pid;
	}
}

int dropcap_run(const char* args[], size_t length)
{
	const char* dropcap_args[drop_caps_len + length + 1];

	dropcap_args[drop_caps_len + length] = NULL;
	for(size_t i = 0; i < drop_caps_len; i++) {
		dropcap_args[i] = drop_caps_runner[i];

		if(i < length)
			dropcap_args[drop_caps_len + i] = args[i];
	}

	for(size_t i = drop_caps_len; i < length; i++)
		dropcap_args[drop_caps_len + i] = args[i];

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


int run_child_with_args(int argc, char** argv)
{

	static const char* pre_chown[]  = {"/usr/bin/chown", "-R", DROP_USER":"DROP_GROUP, "/workspace", NULL};
	static const char* post_chown[] = {"/usr/bin/chown", "-R", "root:root", "/workspace", NULL};

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
	for(size_t i = 1; i < argc; i++)
		cmdline[n++] = argv[i];
	cmdline[n] = NULL;

	if (is_rootless)
		run(false, pre_chown);

	int ret = run(true, cmdline);

	if (is_rootless)
		run(false, post_chown);
	return ret;
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

	static char* const debug_shell[] = {"/usr/bin/bash", "-i", NULL};
	for(size_t i = 1; i < argc; i++) {
		if(strcmp(argv[i], "--root") == 0
		   ||strcmp(argv[i], "-r") == 0)
			execv("/usr/bin/bash", debug_shell);
		;
	}

	if(run_child_with_args(argc, argv) == 0) printf("OK.\n");
	else printf("ERROR.\n");
	return 0;
}

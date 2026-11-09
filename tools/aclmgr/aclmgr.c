#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/wait.h>
#include <unistd.h>
#include <stdbool.h>
#include "aclmgr.env.h"

#define  DROP_USER  "dockuser"
#define  DROP_GROUP "dockuser"

bool is_rootless = false;
volatile sig_atomic_t caught = 0;

#define PREARGS(name, ...) \
	static const char* const name##_args[]  = { \
		__VA_ARGS__ \
	}; \
	static const size_t name##_len = sizeof(name##_args) / sizeof(name##_args[0])

PREARGS(drop_caps,
        "/usr/bin/setpriv",
        "--reuid", DROP_USER,
        "--regid", DROP_GROUP,
        "--init-groups", "--nnp",
        "--inh-caps=-all",
        "--ambient-caps=-all",
        "--bounding-set=-all",
        "--reset-env", "--"
       );

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
		dropcap_args[i] = drop_caps_args[i];

		if(i < length)
			dropcap_args[drop_caps_len + i] = args[i];
	}

	for(size_t i = drop_caps_len; i < length; i++)
		dropcap_args[drop_caps_len + i] = args[i];

	return direct_run(dropcap_args);
}

int run_with_env(bool drop_caps, size_t length, const char* args[], const char* envs[])
{

	if (! drop_caps) return direct_run(args);

	if(envs == NULL) return -1; // dropcap must have envp

	size_t env_count = 1; // /usr/bin/env takes one slot
	for (size_t i = 0; envs[i] != NULL; i++)
		if (is_forwarded(envs[i])) env_count++;

	const char* cmdline[env_count + length + 1];

	cmdline[0] = "/usr/bin/env";
	cmdline[env_count + length] = NULL;

	for (size_t i = 0, j = 1; envs[i] != NULL; i++)
		if (is_forwarded(envs[i])) cmdline[j++] = envs[i];

	for(size_t i = 0; i < length; i++)
		cmdline[env_count + i] = args[i];

	return dropcap_run(cmdline, env_count + length);
}

int run(bool drop_caps, const char* args[], const char* envs[])
{

	size_t i = 0;
	printf("init$ ");

	if (drop_caps) printf("dropcaps ");

	do {
		printf("%s%c", args[i], (args[i + 1] == NULL)?'\n':' ');
	} while(args[++ i] != NULL);

	pid_t pid = -1;

	pid = run_with_env(drop_caps, i, args, envs);
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


int run_child_with_args(int argc, char** argv, char** envp)
{

	static const char* pre_chown[]  = {"/usr/bin/chown", "-R", DROP_USER":"DROP_GROUP, "/workspace", NULL};
	static const char* post_chown[] = {"/usr/bin/chown", "-R", "root:root", "/workspace", NULL};

	const char* cmdline[argc + 1];
	cmdline[0] = "/autodocker/entry.sh";
	cmdline[argc] = NULL;

	for(size_t i = 1; i < argc; i++)
		cmdline[i] = argv[i];

	if (is_rootless)
		run(false, pre_chown, NULL);

	int ret = run(true, cmdline, (const char**) envp);

	if (is_rootless)
		run(false, post_chown, NULL);
	return ret;
}

int main(int argc, char** argv, char** envp)
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
	for(size_t i = 1; i < argc; i++)
		if(strcmp(argv[i], "--root") == 0
		   ||i == 1 && strcmp(argv[i], "-r") == 0)
			execv("/usr/bin/bash", debug_shell);

	size_t ret = run_child_with_args(argc, argv, envp);
	if(ret == 0) printf("OK.\n");
	else printf("ERROR.\n");
	return ret;
}

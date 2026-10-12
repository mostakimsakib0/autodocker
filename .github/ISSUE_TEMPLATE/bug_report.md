---
name: Bug report
about: Something failed, gave a wrong result or was confusing
labels: bug
---

**What happened**
What you expected, and what you got instead (paste error messages exactly).

**Command**
```
docker run ... autodocker ...      (or: python3 autodocker/runner.py ...)
```

**Environment**
- AutoDocker version (`--version`, or `autodocker_version` in `run_manifest.json`):
- How you run it: Docker / Podman, rootless or not / native
- Operating system and CPU architecture (x86-64 or arm64):
- Your user ID (`id -u`, Linux only):

**Files**
Attach `output/run_manifest.json` and, if it ran, the output of
`python3 scripts/verify_run.py output`. If you can share them, attach the
receptor and one ligand that shows the problem.

# What AutoDocker 2.x needs from the build

The build (Dockerfile, `tools/`, `deps` files, `scripts/{apt,shlib,patch}.sh`,
`requirements.txt`, `.github/workflows/`) is back to its v1.1.4 state, plus
Mubashshir's aclmgr from `0cb767b`. This page lists what the build has to
provide for the 2.x pipeline to meet its goals. It also says which
pipeline behaviour depends on each item. The build maintainer decides how to
implement each one. The commits listed show one way it was done before; they
stay in the history and can be cherry-picked or used as reference.

The AutoDocker 2.0.0 results (benchmarks, case study and paper) were produced with the
image published from `e5cc36e`:

    ghcr.io/mostakimsakib0/autodocker:sha-e5cc36e
    ghcr.io/mostakimsakib0/autodocker@sha256:5f75af76d644bd7515b1f6f808ca126cf0c0d30a485f59f47d517d917e8aaf78

To reproduce 2.0.0 results, pull it by that digest. Don't use `:latest` or `:main`;
those now come from the v1.1.4-style build.

## Goals

1. **Same seed, same results.** A run with a given `--seed` gives byte-identical
   poses, scores and reports, on any machine.
2. **Known tool versions.** Every tool in the image is a fixed, recorded release.
3. **Builds that fail loudly** instead of shipping a broken or unpatched tool.
4. **amd64 and arm64** images.
5. **CI that checks all of the above**, so they can't regress unnoticed.

## Required changes

### R1. Open Babel 3.2.1 instead of Debian's 3.1.1 (goal 1, 2)

The pipeline needs Open Babel 3.2.0 or later, and 3.2.1 is what the 2.0.0 results were made with
(`autodocker/vspipeline/utils.py`, `OBABEL_VALIDATED_VERSION`). 3.1.1 ignores
`OB_RANDOM_SEED`, so 3D generation and conformers are not seeded. It also adds a
second hydrogen to pre-protonated Lys/Arg. The pipeline accepts 3.1.1 with a
loud warning, and the results then differ from the 2.0.0 results.

What has to change:
- Build Open Babel from a pinned source (the `tools/openbabel/repo` submodule at the
  `openbabel-3-2-1` tag) with `USE_RANDOMMT=ON`.
- Stop installing Debian's `openbabel`/`libopenbabel-dev` (root `deps`,
  `tools/smina/deps`), and fail the build if a dependency pulls one in.
- Link SMINA against the same 3.2.1 library.
- Set `BABEL_LIBDIR`/`BABEL_DATADIR` for the runtime image.

Reference: `97a3f16`.

### R2. RDKit and fully pinned Python packages (goal 1, 2)

Without RDKit, 2D ligands get 3D coordinates from `obabel --gen3d`, which is not
reproducible. The pipeline logs a warning and records the fallback in
`run_manifest.json`. Without RDKit, a 2D PDB ligand whose bond orders cannot be
read is also docked, when it should be rejected. `requirements.txt` needs `rdkit==2023.9.3` (and `pillow`) plus
exact pins for all transitive packages.

Until then, 4 tests fail in CI because RDKit is missing:
`tests/test_charges.py::test_fix_charges_real_benzene_keeps_coordinates` and
three in `tests/test_ligand_fixes.py` (seeded 3D coordinates for 2D
MOL2 and PDB ligands, and rejection of a 2D PDB without bond orders).

Reference: `f021b13`, `557ff1f` (`requirements.txt`).

### R3. Pinned inputs (goal 2)

- Base images pinned by digest, not by tag.
- Debian packages from snapshot.debian.org at a fixed timestamp.
- AutoDock Vina 1.2.7: the `tools/autodock-vina/repo` submodule at `8eb4040`, not `3c65c0b`.
- NGL 2.5.0 exactly, from a committed lockfile installed with `--frozen-lockfile`.
  The HTML report loads `ngl@2.5.0` (`autodocker/vspipeline/report.py`).
  `tests/test_code_conventions.py::test_report_ngl_version_matches_the_image_build`
  fails until the build pins 2.5.0 again.

Reference: `557ff1f`.

### R4. Floating-point settings that are the same on every architecture (goal 1, 4)

Compile every native tool with `-ffp-contract=off`. GCC fuses multiply-adds by
default on arm64 but not on baseline x86-64, which changes rounding and therefore
scores.

Reference: `df57f46`.

### R5. arm64 (goal 4)

- Build each architecture natively. Under QEMU the SMINA build takes hours.
- fpocket ships its molfile plugin only for amd64. On arm64 a stub is needed; the
  pipeline uses only PDB input and output, which never call the plugins.
- `shlib.sh` must strip the `:arch` qualifier, so the runtime package list
  is the same on both architectures.

Reference: `df57f46` (Dockerfile, `publish.yml`, `tools/fpocket/molfile_stub.c`).

### R6. Fail-loud build scripts (goal 3)

- `patch.sh`: stop the build when a patch no longer applies (after a submodule
  bump), instead of skipping it silently.
- `shlib.sh`: fail when `ldd`, `realpath` or `dpkg-query` fails, or when nothing
  matches. Otherwise the image ships a binary without its libraries.
- Each tool stage checks its result: the version string, the PDBQT plugin
  and a test conversion.

Reference: `557ff1f`, `5b7506d`.

### R7. Version metadata (goal 2)

`autodocker/VERSION` is the single source of the version. The build should
fail when `AUTODOCKER_VERSION` disagrees with it. Set the OCI labels and
`AUTODOCKER_VERSION` in the image; the run manifest records them.

Reference: `557ff1f`.

### R8. Container entry (aclmgr)

Mubashshir's aclmgr (`0cb767b`) is kept. Three things need his decision:
- `-r` as the first argument opens a root shell, but `-r` is the pipeline's
  `--receptor` flag (`autodocker/runner.py`). `docker run IMAGE -r protein.pdb …`
  must start a docking run. `README.md` and `autodocker/entry.sh` document
  `--root-shell`.
- Under a rootful engine, `dockuser` (UID 1000) can write a bind-mounted
  `/workspace` only when the host user is also UID 1000.
- Under a rootless engine, if the container is killed between the two `chown -R`
  calls, the user's files stay owned by a sub-UID.

`gcc -Wall -Wextra` warnings: a missing return in `direct_run`, `&&`/`||`
without parentheses, signed/unsigned comparisons, and an exit status held in a `size_t`.

### R9. CI (goal 5)

The v1.1.4 workflows run pytest only. 2.x relied on these checks:
- pytest on Python 3.9–3.12 with Open Babel 3.2.1 and Vina 1.2.7 installed
  (`a3b04d5`)
- ruff, shellcheck and hadolint (`557ff1f`, `a5c3a73`)
- a Docker smoke test plus a determinism check on amd64 and arm64
  (`.github/smoke/`, `df57f46`, `2fb71c2`)
- the 1HSG case study on a second machine and on arm64, run by hand
  (`.github/repro/`, `154b2be`)
- sharded benchmarks (`scripts/ci/`, `a3d50fe`)

The helper scripts for these (`.github/smoke`, `.github/repro`, `scripts/ci`) are
still in the repository; no workflow calls them now.

### R10. Bit-for-bit image reproducibility (optional)

The goals above need reproducible inputs and results, not an identical image.
The ways an image differs between builds are `.pyc` source mtimes, layer
timestamps, apt/dpkg logs, `useradd` dates and the `created` label. Each has a
small fix: `SOURCE_DATE_EPOCH` with `rewrite-timestamp`, `checked-hash` bytecode,
and log cleanup. None of them affects results.

## Documentation that no longer matches the image

These describe the 2.0.0 build: `README.md` (install, architectures,
environment), `docs/reproducibility.md`, `docs/openbabel.md`, `BENCHMARKS.md`,
`CHANGELOG.md` (2.0.0 section) and `paper/paper.md`. Update them once the
build decisions are made, or point them at the 2.0.0 image digest above.

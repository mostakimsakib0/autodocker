# External user testing

This template records feedback from users who are **not** the developers,
run AutoDocker on a fresh machine, and report results. It is part of the
JOSS submission evidence: independent installation + first-run success.

## The 15-minute test

Please do this on your own machine, using only the README and this page; do
not ask the developers until you have filled in the report. Note every place
where you got stuck or had to guess. That is the main thing we want to learn.

1. **Install.** Install Docker (or use Podman), then follow the README's
   *Quick start* to get the image (`docker pull ...` and `docker tag ...`).
2. **Prepare a small screen.** In an empty directory:

   ```bash
   mkdir -p works/ligs
   curl -fsSL https://files.rcsb.org/download/1HSG.pdb -o works/protein.pdb
   git clone --depth 1 https://github.com/mostakimsakib0/autodocker.git src
   cp src/.github/smoke/ligs/*.sdf works/ligs/
   ```

   (1HSG is HIV-1 protease with the inhibitor indinavir, residue MK1; the two
   ligands are small test molecules.)
3. **Run it**, with the box on the crystal ligand:

   ```bash
   docker run --rm -v "$PWD/works:/workspace" autodocker --chain all --box-from-het MK1
   ```

   The outputs should belong to you (`ls -ln works/output` shows your
   `id -u`). If the run stops with `Permission denied`, run
   `chmod -R a+rwX works`, retry, and **record that this happened**, with
   your `id -u` and how you run containers (Docker, Podman, rootless or not).
4. **Look at the results.** Open `works/output/ranking.csv` and
   `works/output/results_report.html`. Does the box drawn in the report sit
   on the inhibitor?
5. **Check that it is reproducible.** Copy the inputs to a second directory,
   run the same command there, and compare:

   ```bash
   mkdir -p works2 && cp -r works/protein.pdb works/ligs works2/
   docker run --rm -v "$PWD/works2:/workspace" autodocker --chain all --box-from-het MK1
   sha256sum works/output/ranking.csv works2/output/ranking.csv
   python3 src/scripts/verify_run.py works/output
   python3 src/scripts/verify_run.py works/output works2/output
   ```

   The two `ranking.csv` hashes should be equal and both `verify_run.py`
   calls should end with exit status 0 (`echo $?`).

## Report template

```
Machine:
  OS / kernel: 
  Python version:           (native install only)
  Docker version:           (container install only)
  RAM / CPU cores:
  Installation method:      Docker | native

Steps performed (from docs/installation.md):
  1) 
  2) 
  3) 

Time taken (wall clock):
  build / install:
  first screening run:

Receptor used:        PDB id / source:
Ligands:             count / format:
Command used:
  docker run ...   or   python3 autodocker/runner.py ...

Outputs produced (check that apply):
  [ ] ranking.csv     [ ] Top_hits.txt   [ ] metrics.txt
  [ ] pocket_summary.csv    [ ] docked/*_out.pdbqt
  [ ] grid_box.pml    [ ] report.html (when requested)

15-minute test (above):
  user ID (`id -u`, Linux only):
  needed `chmod -R a+rwX works`:        yes | no
  box on the inhibitor in the report:   yes | no | could not tell
  ranking.csv hashes equal:             yes | no
  verify_run.py exit status (one run / two runs):

Optional, longer determinism check (docs/reproducibility.md):
  scripts/reproducibility_check.sh result: DETERMINISTIC | NON-DETERMINISTIC

Anything that failed or was confusing:
  
Any suggestions:
  
Rating (1 easy - 5 hard):
  install:
  first run:
```

## Completed reports

| Date | Tester | Method | Install | First run | Determinism | Issue/link |
|---|---|---|---|---|---|---|
| (add rows as they are reported) | | | | | | |
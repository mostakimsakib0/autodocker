# External user testing

This template records feedback from users who are **not** the developers,
run AutoDocker on a fresh machine, and report results. It is part of the
JOSS submission evidence: independent installation + first-run success.

## How to participate

1. Pick a machine you have not used AutoDocker on before (or a fresh VM).
2. Follow `docs/installation.md` (Docker or native).
3. Run your own protein (`protein.pdb`) and ligand files (`ligs/*.sdf`) or
   the small bundled example.
4. Fill in the report below and open an issue or link it here.

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
  [ ] grid_box.py     [ ] report.html (when requested)

Determinism check (docs/reproducibility.md):
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
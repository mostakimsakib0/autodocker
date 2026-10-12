#!/usr/bin/env python3
"""Baseline: the standard AutoDock Vina workflow on the AutoDocker benchmarks.

Runs the preparation recommended by the AutoDock Vina documentation (Meeko)
and raw Vina on the SAME complexes/targets, boxes, seeds and docking
settings as ``benchmark_redock.py`` (J1) and ``benchmark_enrichment.py``
(J2), and writes results in their CSV/summary.json schema, so AutoDocker
and the baseline can be put side by side.

Preparation chain
-----------------
  * Receptor : PDB element columns filled from the atom names where blank
               (DUD-E receptors have none) and DUD-E's own names of
               standard amino acids restored (LEV, GLZ, ASQ, MEU -> LEU,
               GLY, ASP, MET; Meeko has no template for them and cannot
               delete them from inside a chain); no other change. Then
               ``reduce -Trim`` (remove hydrogens), ``reduce -BUILD`` (add
               hydrogens, flip Asn/Gln/His) and ``mk_prepare_receptor.py -i
               receptorH.pdb -o receptor -p`` (Meeko, ProDy reader, Gasteiger
               charges from templates). If Meeko refuses residues that match
               no template, it is re-run with ``-x`` (delete them) and the
               deleted residues are recorded in ``prep.json``.
  * Ligand   : an RDKit molecule with explicit hydrogens written as SDF, then
               ``mk_prepare_ligand.py -i ligand.sdf -o ligand.pdbqt``.
               ``redock``: the start conformer of benchmark_redock.py
               (``--start-conformer rdkit``: the same seeded ETKDGv3 + MMFF
               conformer from ``_benchlib.rdkit_start_conformer``, InChIKey
               checked against the crystal ligand; ``crystal``: the fetch
               script's ``{pdbid}_ligand.sdf``). ``enrichment``: the DUD-E
               entry docked by benchmark_enrichment.py (same compounds, same
               entry; read by RDKit with its hydrogens, DUD-E coordinates).
  * Docking  : raw Vina (``--vina-bin``, 1.2.7 for the paper) with the box
               from the crystal ligand (``_benchlib.ligand_box``: centroid,
               edge max(22.5, extent + 12) A), ``--seed``, exhaustiveness,
               ``--num_modes`` and ``--energy_range`` of the AutoDocker
               protocol (redock: pipeline defaults; enrichment: 4, 1, 1.0),
               ``--cpu --threads``. The exact command is logged at the top
               of ``<name>_vina.log`` like AutoDocker's logs.
  * Metrics  : redock: obrms RMSD of every mode vs. the crystal MOL2
               (``benchmark_redock.mode_rmsds``/``summarize``); enrichment:
               AUC, EF, BEDROC, ROCe, bootstrap CIs and Mann-Whitney
               (``benchmark_enrichment.compute_stats``/``save_tables``).
               Poses stay in the input frame, so no frame shift is applied.

``--engine easydock`` docks the same prepared receptor and start ligands
with EasyDock (``easydock --program vina``; EasyDock prepares the ligand
with Meeko itself, keeps 3D input coordinates and docks with the Vina
Python API, searching 50 poses and returning ``--binding-modes``). If the
``easydock`` executable is missing or does not start, the run stops with
exit status 3 and a summary.json that says why.

The Meeko/reduce tools may live in their own virtual environment
(``--meeko-bin DIR``); their package versions are read with that
environment's interpreter and recorded in summary.json.

Results: ``outdir/results.csv`` + ``outdir/summary.json`` (redock),
``outdir/enrichment.csv`` + ``outdir/summary.json`` + ``outdir/<target>/``
(enrichment), exactly as the AutoDocker benchmarks write them.
"""

import argparse
import csv
import hashlib
import json
import multiprocessing
import os
import random
import shlex
import shutil
import sqlite3
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _benchlib  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))


def _load_script(name):
    """A sibling benchmark script as a module (their metric functions are
    reused so both sides are scored by the same code)."""
    import importlib.util
    spec = importlib.util.spec_from_file_location(f"_baseline_{name}",
                                                  os.path.join(HERE, f"{name}.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


br = _load_script("benchmark_redock")
be = _load_script("benchmark_enrichment")
vr = _load_script("verify_run")

# results.csv columns of benchmark_redock.py (checked in tests).
REDOCK_COLUMNS = ["pdbid", "status", "affinity", "affinity_mode1",
                  "rmsd_top1", "rmsd_best", "rmsd_modes",
                  "n_modes_rmsd_failed", "reused", "error"]
ENGINES = ("vina", "easydock")
PREP_TIMEOUT = 600            # one reduce / Meeko call
EXIT_ENGINE_UNAVAILABLE = 3
# Python packages whose versions define the preparation (read with the
# interpreter of the Meeko scripts).
PREP_PACKAGES = ("meeko", "rdkit", "numpy", "scipy", "gemmi", "prody",
                 "reduce-binary", "molscrub")
EASYDOCK_PACKAGES = ("easydock", "vina", "meeko", "rdkit", "numpy", "dask")
ENRICHMENT_MODES = 1           # benchmark_enrichment.py: --binding-modes 1
ENRICHMENT_ENERGY_RANGE = 1.0  # benchmark_enrichment.py: --energy-range 1.0
ENRICHMENT_EXHAUSTIVENESS = 4  # benchmark_enrichment.py default
TWO_LETTER_ELEMENTS = {"FE", "ZN", "MG", "MN", "CA", "CL", "BR", "NA", "CU", "CO",
                       "NI", "CD", "HG", "SE", "LI", "AL", "SR", "BA", "CS", "RB"}


# ---------------------------------------------------------------------------
# Tools and versions
# ---------------------------------------------------------------------------

def find_tool(name, bin_dir=None):
    """``bin_dir/name`` when given and present, else ``name`` on PATH."""
    if bin_dir:
        path = os.path.join(bin_dir, name)
        return path if os.path.exists(path) else None
    return shutil.which(name)


def script_interpreter(script):
    """Python interpreter that runs a console script: the ``python3`` next to
    it (virtual environment), else its shebang (incl. pip's ``exec`` form)."""
    if not script:
        return None
    for cand in ("python3", "python"):
        path = os.path.join(os.path.dirname(script), cand)
        if os.path.exists(path):
            return path
    try:
        with open(script, encoding="utf-8", errors="replace") as f:
            head = [f.readline() for _ in range(3)]
    except OSError:
        return None
    for line in head:
        if line.startswith("#!") and "python" in line:
            return line[2:].split()[0]
        if line.startswith("'''exec' "):
            return shlex.split(line[len("'''exec' "):])[0]
    return None


def package_versions(interpreter, packages):
    """``{package: version or None}`` as seen by ``interpreter``."""
    if not interpreter:
        return {p: None for p in packages}
    code = ("import json, sys, importlib.metadata as m\n"
            "out = {}\n"
            "for p in sys.argv[1:]:\n"
            "    try:\n        out[p] = m.version(p)\n"
            "    except Exception:\n        out[p] = None\n"
            "print(json.dumps(out))\n")
    try:
        proc = subprocess.run([interpreter, "-I", "-c", code, *packages], capture_output=True,
                              text=True, timeout=_benchlib.VERSION_TIMEOUT)
        return json.loads(proc.stdout)
    except (OSError, ValueError, subprocess.SubprocessError):
        return {p: None for p in packages}


def reduce_het_dict(reduce_bin, explicit=None):
    """reduce's het dictionary: ``--reduce-het-dict``, else the one the
    reduce-binary wheel ships next to its binary (``../reduce_wwPDB_het_dict.txt``)."""
    if explicit:
        return explicit
    if not reduce_bin:
        return None
    real = os.path.realpath(reduce_bin)
    cand = os.path.join(os.path.dirname(os.path.dirname(real)), "reduce_wwPDB_het_dict.txt")
    return cand if os.path.exists(cand) else None


def wheel_reduce(interpreter):
    """The reduce binary of the reduce-binary wheel installed for
    ``interpreter`` (None when it has none)."""
    if not interpreter:
        return None
    try:
        proc = subprocess.run(
            [interpreter, "-I", "-c",
             "import reduce_binary; print(reduce_binary.REDUCE_BIN_PATH)"],
            capture_output=True, text=True, timeout=_benchlib.VERSION_TIMEOUT)
    except (OSError, subprocess.SubprocessError):
        return None
    path = proc.stdout.strip()
    return path if proc.returncode == 0 and path and os.path.exists(path) else None


def resolve_tools(args):
    """Paths of every external tool of the chosen chain (None = missing).
    reduce: --reduce-bin, else ``reduce`` in --meeko-bin, else the
    reduce-binary wheel of the Meeko interpreter, else PATH."""
    meeko_bin = args.meeko_bin
    receptor_script = find_tool("mk_prepare_receptor.py", meeko_bin)
    reduce_bin = (args.reduce_bin or (find_tool("reduce", meeko_bin) if meeko_bin else None)
                  or wheel_reduce(script_interpreter(receptor_script))
                  or shutil.which("reduce"))
    tools = {
        "reduce": reduce_bin,
        "reduce_het_dict": reduce_het_dict(reduce_bin, args.reduce_het_dict),
        "mk_prepare_receptor": receptor_script,
        "mk_prepare_ligand": find_tool("mk_prepare_ligand.py", meeko_bin),
        "vina": args.vina_bin or shutil.which("vina"),
        "obabel": shutil.which("obabel"),
        "easydock": (args.easydock_bin or shutil.which("easydock")
                     if args.engine == "easydock" else None),
    }
    return tools


def tool_record(tools):
    """Versions and SHA-256 of the tools, recorded in summary.json."""
    rec = {}
    for name in ("reduce", "vina", "obabel"):
        path = tools.get(name)
        args = ("-version",) if name == "reduce" else (("-V",) if name == "obabel"
                                                      else ("--version",))
        rec[name] = {"path": path, "version": _benchlib.tool_version(path, args) if path else None,
                     "sha256": _benchlib.file_sha256(os.path.realpath(path)) if path else None}
    rec["reduce_het_dict"] = {"path": tools.get("reduce_het_dict"),
                              "sha256": _benchlib.file_sha256(tools.get("reduce_het_dict"))}
    interp = script_interpreter(tools.get("mk_prepare_receptor"))
    rec["meeko_python"] = {"interpreter": interp,
                           "packages": package_versions(interp, PREP_PACKAGES)}
    if tools.get("easydock"):
        interp = script_interpreter(tools["easydock"])
        rec["easydock"] = {"path": tools["easydock"], "interpreter": interp,
                           "packages": package_versions(interp, EASYDOCK_PACKAGES)}
    rec["rdkit_start_conformers"] = _benchlib.rdkit_version()
    return rec


def missing_tools(tools, engine):
    needed = ["reduce", "mk_prepare_receptor", "mk_prepare_ligand", "obabel"]
    needed.append("vina" if engine == "vina" else "easydock")
    return [n for n in needed if not tools.get(n)]


def easydock_available(tools):
    """``(ok, reason)``: the easydock executable exists and starts."""
    exe = tools.get("easydock")
    if not exe or not os.path.exists(exe):
        return False, ("easydock executable not found (pip install easydock plus its "
                       "dependencies rdkit, meeko, gemmi, vina, dask, distributed, "
                       "scikit-learn, prody; or pass --easydock-bin)")
    try:
        proc = subprocess.run([exe, "--help"], capture_output=True, text=True,
                              timeout=_benchlib.RUNNER_VERSION_TIMEOUT)
    except (OSError, subprocess.SubprocessError) as e:
        return False, f"easydock does not start: {e}"
    if proc.returncode != 0:
        tail = (proc.stderr or proc.stdout).strip().splitlines()
        return False, f"easydock does not start: {tail[-1] if tail else proc.returncode}"
    return True, ""


# ---------------------------------------------------------------------------
# Receptor preparation
# ---------------------------------------------------------------------------

def element_from_name(name):
    """PDB element of an atom name (columns 13-16) whose element column is
    blank: right-justified two-letter elements (``FE``, ``ZN``, ...) start
    in column 13; otherwise the first letter."""
    name = name.ljust(4)
    two = name[:2].strip().upper()
    if name[0].isalpha() and two in TWO_LETTER_ELEMENTS:
        return two
    letters = [c for c in name if c.isalpha()]
    return letters[0].upper() if letters else ""


# DUD-E's names for standard amino acids in some receptors (cdk2 LEV, hivpr
# GLZ/ASQ, egfr MEU; the atom names are those of the standard residue).
DUDE_RESIDUE_NAMES = {"LEV": "LEU", "GLZ": "GLY", "ASQ": "ASP", "MEU": "MET"}


def fill_pdb_elements(src, dst, renamed=None):
    """Copy a PDB file, writing the element (columns 77-78) of ATOM/HETATM
    records that have none and restoring the standard name of the residues
    in DUDE_RESIDUE_NAMES. Returns the number of records whose element was
    filled; ``renamed`` (a dict) receives {old name: records renamed}."""
    changed = 0
    with open(src, encoding="utf-8", errors="replace") as fin, \
            open(dst, "w", encoding="utf-8") as fout:
        for line in fin:
            if line.startswith(("ATOM", "HETATM")):
                body = line.rstrip("\r\n").ljust(80)
                if not body[76:78].strip():
                    body = body[:76] + element_from_name(body[12:16]).rjust(2) + body[78:]
                    changed += 1
                standard = DUDE_RESIDUE_NAMES.get(body[17:20])
                if standard:
                    if renamed is not None:
                        renamed[body[17:20]] = renamed.get(body[17:20], 0) + 1
                    body = body[:17] + standard + body[20:]
                line = body.rstrip() + "\n"
            fout.write(line)
    return changed


def _run(cmd, log, stdout_path=None, timeout=PREP_TIMEOUT, ok_codes=(0,)):
    """Run one preparation command, appending it and its output to ``log``.
    ``stdout_path`` receives stdout (reduce writes the PDB there)."""
    log.write(f"# Command: {shlex.join(cmd)}\n")
    log.flush()
    env = {**os.environ, "OMP_NUM_THREADS": "1"}
    if stdout_path:
        with open(stdout_path, "w", encoding="utf-8") as out:
            proc = subprocess.run(cmd, stdout=out, stderr=subprocess.PIPE, text=True,
                                  timeout=timeout, env=env)
        stdout = ""
    else:
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout, env=env)
        stdout = proc.stdout or ""
    log.write(f"# Exit status: {proc.returncode}\n{stdout}{proc.stderr or ''}\n")
    if ok_codes is not None and proc.returncode not in ok_codes:
        raise RuntimeError(f"{os.path.basename(cmd[0])} failed (exit {proc.returncode}): "
                           f"{(proc.stderr or stdout).strip()[-300:]}")
    return proc.returncode, stdout + (proc.stderr or "")


def _has_atoms(pdb):
    if not os.path.exists(pdb):
        return False
    with open(pdb, encoding="utf-8", errors="replace") as f:
        return any(line.startswith(("ATOM", "HETATM")) for line in f)


def deleted_residues(meeko_output):
    """Residue lines Meeko reports as deleted/unmatched (best effort)."""
    keep = []
    for line in meeko_output.splitlines():
        low = line.lower()
        if ("delet" in low or "no template" in low or "unmatched" in low
                or "did not match" in low):
            keep.append(line.strip())
    return keep[:50]


def prepare_receptor(pdb, workdir, tools):
    """reduce -Trim, reduce -BUILD, mk_prepare_receptor.py. Returns
    ``(receptor.pdbqt, record)``; the record holds the commands and whether
    unmatched residues were deleted."""
    os.makedirs(workdir, exist_ok=True)
    rec = {"input": os.path.abspath(pdb), "input_sha256": _benchlib.file_sha256(pdb)}
    with open(os.path.join(workdir, "receptor_prep.log"), "w", encoding="utf-8") as log:
        elem = os.path.join(workdir, "receptor_elements.pdb")
        renamed = {}
        rec["elements_filled"] = fill_pdb_elements(pdb, elem, renamed)
        rec["residues_renamed"] = renamed
        trim = os.path.join(workdir, "receptor_trim.pdb")
        built = os.path.join(workdir, "receptorH.pdb")
        # reduce's exit status is not an error indicator (it is non-zero when
        # it changed nothing, e.g. -Trim on a structure without hydrogens);
        # the output is checked instead.
        rec["reduce_trim_exit"], _ = _run([tools["reduce"], "-Trim", elem], log, trim,
                                          ok_codes=None)
        if not _has_atoms(trim):
            raise RuntimeError("reduce -Trim wrote no atoms")
        cmd = [tools["reduce"], "-BUILD"]
        if tools.get("reduce_het_dict"):
            cmd += ["-DB", tools["reduce_het_dict"]]
        rec["reduce_build_exit"], _ = _run(cmd + [trim], log, built, ok_codes=None)
        if not _has_atoms(built):
            raise RuntimeError("reduce -BUILD wrote no atoms")
        base = os.path.join(workdir, "receptor")
        out = base + ".pdbqt"
        cmd = [tools["mk_prepare_receptor"], "-i", built, "-o", base, "-p"]
        rec["mk_prepare_receptor"] = shlex.join(cmd)
        rec["deleted_unmatched_residues"] = False
        code, text = _run(cmd, log, ok_codes=None)
        if code != 0 or not os.path.exists(out):
            code, text = _run(cmd + ["-x"], log)
            rec["mk_prepare_receptor"] = shlex.join(cmd + ["-x"])
            rec["deleted_unmatched_residues"] = True
            rec["meeko_messages"] = deleted_residues(text)
    if not os.path.exists(out):
        raise RuntimeError("mk_prepare_receptor.py wrote no PDBQT")
    rec["receptor_pdbqt_sha256"] = _benchlib.file_sha256(out)
    return out, rec


# ---------------------------------------------------------------------------
# Ligand preparation
# ---------------------------------------------------------------------------

def mol2_to_sdf(mol2, sdf, seed=None):
    """Open Babel MOL2 -> SDF (coordinates, bond orders, charges and
    hydrogens kept), with the benchmarks' Open Babel environment."""
    proc = subprocess.run(["obabel", "-imol2", mol2, "-osdf", "-O", sdf], capture_output=True,
                          text=True, timeout=_benchlib.TOOL_TIMEOUT,
                          env=_benchlib.obabel_env(seed))
    if proc.returncode != 0 or not os.path.exists(sdf) or not os.path.getsize(sdf):
        raise RuntimeError(f"obabel mol2->sdf failed: {proc.stderr[-300:]}")
    return sdf


def mol2_block_to_sdf(block, sdf, name):
    """One DUD-E MOL2 entry -> SDF with all its hydrogens, via RDKit (no
    coordinate change). Raises RuntimeError when RDKit cannot read it."""
    from rdkit import Chem
    mol = Chem.MolFromMol2Block(block, removeHs=False, sanitize=True)
    if mol is None:
        raise RuntimeError("RDKit could not read the MOL2 entry")
    if not any(a.GetAtomicNum() == 1 for a in mol.GetAtoms()):
        mol = Chem.AddHs(mol, addCoords=True)
    mol.SetProp("_Name", name)
    with open(sdf, "w", encoding="utf-8") as f:
        f.write(Chem.MolToMolBlock(mol))
        f.write("$$$$\n")
    return sdf


def prepare_ligand_pdbqt(sdf, pdbqt, tools, log_path):
    with open(log_path, "w", encoding="utf-8") as log:
        _run([tools["mk_prepare_ligand"], "-i", sdf, "-o", pdbqt], log)
    if not os.path.exists(pdbqt) or not os.path.getsize(pdbqt):
        raise RuntimeError("mk_prepare_ligand.py wrote no PDBQT")
    return pdbqt


def redock_start_sdf(pdbid, complex_dir, work, opts):
    """The start ligand of benchmark_redock.py as SDF (with hydrogens).

    rdkit: benchmark_redock.prepare_start_ligand writes the same MOL2 the
    AutoDocker benchmark docks (seeded ETKDGv3 + MMFF, InChIKey checked);
    it is converted to SDF without changing coordinates and checked again.
    crystal: the fetch script's ``{pdbid}_ligand.sdf``."""
    if opts["start_conformer"] == "crystal":
        src = os.path.join(complex_dir, f"{pdbid}_ligand.sdf")
        sdf = os.path.join(work, "input", f"{pdbid}_ligand.sdf")
        os.makedirs(os.path.dirname(sdf), exist_ok=True)
        shutil.copy(src, sdf)
        meta = {"method": "crystal"}
    else:
        mol2, meta = br.prepare_start_ligand(pdbid, complex_dir, work, opts)
        sdf = mol2[:-len(".mol2")] + ".sdf"
        mol2_to_sdf(mol2, sdf, opts["conformer_seed"])
    br.check_same_molecule(os.path.join(complex_dir, f"{pdbid}_ligand.mol2"), sdf)
    return sdf, meta


# ---------------------------------------------------------------------------
# Docking
# ---------------------------------------------------------------------------

def vina_command(tools, receptor, ligand, out, center, size, opts):
    cmd = [tools["vina"], "--receptor", receptor, "--ligand", ligand, "--out", out]
    for axis, value in zip("xyz", center):
        cmd += [f"--center_{axis}", f"{value:.3f}"]
    for axis in "xyz":
        cmd += [f"--size_{axis}", f"{size:.3f}"]
    cmd += ["--exhaustiveness", str(opts["exhaustiveness"]),
            "--num_modes", str(opts["binding_modes"]),
            "--energy_range", str(opts["energy_range"]),
            "--seed", str(opts["seed"])]
    if opts.get("threads"):
        cmd += ["--cpu", str(opts["threads"])]
    return cmd


def dock_vina(tools, receptor, ligand, out, log_path, center, size, opts):
    """Raw Vina; the command is the first line of the log (like AutoDocker's
    ``# Command:`` header). Returns the mode-1 score token."""
    cmd = vina_command(tools, receptor, ligand, out, center, size, opts)
    with open(log_path, "w", encoding="utf-8") as log:
        log.write(f"# Command: {shlex.join(cmd)}\n")
        log.flush()
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=opts["timeout"],
                              env={**os.environ, "OMP_NUM_THREADS": "1"})
        log.write(proc.stdout or "")
        if proc.stderr:
            log.write("# STDERR:\n" + proc.stderr)
    if proc.returncode != 0 or not os.path.exists(out):
        raise RuntimeError(f"vina failed (exit {proc.returncode}): {(proc.stderr or '')[-300:]}")
    token, _ = _benchlib.first_vina_result(out)
    if token is None:
        raise RuntimeError("no REMARK VINA RESULT in the Vina output")
    return token


def write_easydock_config(workdir, receptor, center, size, opts):
    setup = os.path.join(workdir, "easydock_box.txt")
    with open(setup, "w", encoding="utf-8") as f:
        for axis, value in zip("xyz", center):
            f.write(f"center_{axis} = {value:.3f}\n")
        for axis in "xyz":
            f.write(f"size_{axis} = {size:.3f}\n")
    config = os.path.join(workdir, "easydock_config.yml")
    with open(config, "w", encoding="utf-8") as f:
        f.write(f"protein: {os.path.abspath(receptor)}\n"
                f"protein_setup: {os.path.abspath(setup)}\n"
                f"exhaustiveness: {opts['exhaustiveness']}\n"
                f"seed: {opts['seed']}\n"
                f"n_poses: {opts['binding_modes']}\n"
                f"ncpu: {opts.get('threads') or 1}\n")
    return config


def dock_easydock(tools, receptor, sdf, workdir, center, size, opts, ncpu=1):
    """Dock every molecule of ``sdf`` with EasyDock (Vina). Returns
    ``{name: pose PDBQT path}`` for the docked molecules."""
    os.makedirs(workdir, exist_ok=True)
    config = write_easydock_config(workdir, receptor, center, size, opts)
    db = os.path.join(workdir, "easydock.db")
    if os.path.exists(db):
        os.remove(db)
    cmd = [tools["easydock"], "-i", sdf, "-o", db, "--program", "vina",
           "--config", config, "-c", str(ncpu)]
    with open(os.path.join(workdir, "easydock.log"), "w", encoding="utf-8") as log:
        _run(cmd, log, timeout=opts["timeout"])
    poses = {}
    with sqlite3.connect(db) as conn:
        rows = conn.execute("SELECT id, raw_block FROM mols WHERE docking_score IS NOT NULL "
                            "AND raw_block IS NOT NULL ORDER BY id, stereo_id").fetchall()
    for name, raw in rows:
        if name in poses:
            continue  # first stereoisomer (EasyDock -s 1 enumerates one)
        path = os.path.join(workdir, f"{name}_out.pdbqt")
        with open(path, "w", encoding="utf-8") as f:
            f.write(raw if raw.endswith("\n") else raw + "\n")
        poses[name] = path
    return poses


# ---------------------------------------------------------------------------
# Redock mode
# ---------------------------------------------------------------------------

def baseline_params(opts):
    keys = ("engine", "seed", "exhaustiveness", "binding_modes", "energy_range", "padding",
            "min_box", "threads", "start_conformer", "conformer_seed", "rdkit_version",
            "tool_versions")
    return {k: opts.get(k) for k in keys}


def redock_complex(job):
    """One complex: same record layout as benchmark_redock.process_complex."""
    pdbid, data_dir, workroot, opts = job
    rec = {"pdbid": pdbid, "status": "error", "affinity": None, "affinity_mode1": None,
           "rmsd_top1": None, "rmsd_best": None, "rmsd_modes": [], "reused": False,
           "error": None}
    complex_dir = os.path.join(data_dir, pdbid)
    protein = os.path.join(complex_dir, f"{pdbid}_protein.pdb")
    crystal = os.path.join(complex_dir, f"{pdbid}_ligand.mol2")
    work = os.path.join(workroot, pdbid)
    out = os.path.join(work, f"{pdbid}_ligand" + _benchlib.POSE_SUFFIX)
    tools = opts["tools"]
    try:
        if not (os.path.exists(protein) and os.path.exists(crystal)):
            rec["status"] = "missing"
            rec["error"] = opts["manifest"].get(pdbid) or "not in --data-dir"
            return rec
        params = {**baseline_params(opts), "inputs": br.complex_inputs(pdbid, complex_dir)}
        if os.path.exists(out) and br.params_match(work, params):
            rec["reused"] = True
        else:
            shutil.rmtree(work, ignore_errors=True)
            os.makedirs(work)
            cx, cy, cz, size = br.ligand_grid(crystal, opts["padding"], opts["min_box"])
            receptor, rec_prep = prepare_receptor(protein, os.path.join(work, "receptor"), tools)
            sdf, conf_meta = redock_start_sdf(pdbid, complex_dir, work, opts)
            if opts["engine"] == "easydock":
                poses = dock_easydock(tools, receptor, sdf, os.path.join(work, "easydock"),
                                      (cx, cy, cz), size, opts)
                if not poses:
                    raise RuntimeError("EasyDock returned no pose")
                shutil.copy(next(iter(poses.values())), out)
            else:
                pdbqt = os.path.join(work, "input", f"{pdbid}_ligand.pdbqt")
                prepare_ligand_pdbqt(sdf, pdbqt, tools,
                                     os.path.join(work, "input", "ligand_prep.log"))
                dock_vina(tools, receptor, pdbqt, out,
                          os.path.join(work, f"{pdbid}_ligand" + _benchlib.VINA_LOG_SUFFIX),
                          (cx, cy, cz), size, opts)
            with open(os.path.join(work, "params.json"), "w", encoding="utf-8") as f:
                json.dump({"docking": params, "start_conformer": conf_meta,
                           "receptor_preparation": rec_prep,
                           "box": {"center_input_frame": [cx, cy, cz], "edge": size}},
                          f, indent=2)
        token, value = _benchlib.first_vina_result(out)
        rec["affinity"] = rec["affinity_mode1"] = value
        # Poses are in the input frame: the crystal MOL2 is the reference as is.
        br.apply_rmsds(rec, *br.mode_rmsds(crystal, out))
    except Exception as e:
        rec["status"] = "error"
        rec["error"] = str(e)[:500]
    return rec


def write_redock_results(results, outdir):
    path = os.path.join(outdir, "results.csv")
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(REDOCK_COLUMNS)
        for r in results:
            writer.writerow([
                r["pdbid"], r["status"],
                "" if r["affinity"] is None else r["affinity"],
                "" if r["affinity_mode1"] is None else r["affinity_mode1"],
                "" if r["rmsd_top1"] is None else f"{r['rmsd_top1']:.3f}",
                "" if r["rmsd_best"] is None else f"{r['rmsd_best']:.3f}",
                " ".join("NA" if x is None else f"{x:.3f}" for x in r["rmsd_modes"]),
                r.get("n_modes_rmsd_failed", ""), r["reused"], r["error"] or ""])
    return path


def _map(func, jobs, processes):
    if processes <= 1 or len(jobs) <= 1:
        for job in jobs:
            yield func(job)
        return
    with multiprocessing.Pool(processes=processes) as pool:
        yield from pool.imap_unordered(func, jobs)


def run_redock(args, tools, tool_versions):
    pdbids = br.select_ids(args)
    opts = {"engine": args.engine, "padding": args.padding, "min_box": args.min_box,
            "seed": args.seed, "exhaustiveness": args.exhaustiveness,
            "binding_modes": args.binding_modes, "energy_range": args.energy_range,
            "threads": args.threads, "start_conformer": args.start_conformer,
            "conformer_seed": args.conformer_seed if args.start_conformer == "rdkit" else None,
            "rdkit_version": (_benchlib.rdkit_version()
                              if args.start_conformer == "rdkit" else None),
            "tool_versions": _versions_only(tool_versions),
            "timeout": args.timeout, "manifest": br.read_manifest(args.data_dir),
            "tools": tools}
    workroot = os.path.join(args.outdir, "work")
    jobs = [(p, args.data_dir, workroot, opts) for p in pdbids]
    results = []
    for rec in _map(redock_complex, jobs, args.processes):
        results.append(rec)
        print(f"[{rec['pdbid']}] {rec['status']}{' (reused)' if rec['reused'] else ''} "
              f"aff={rec['affinity']} rmsd_top1={rec['rmsd_top1']} "
              f"rmsd_best={rec['rmsd_best']} err={rec['error']}", flush=True)
    results.sort(key=lambda r: r["pdbid"])
    write_redock_results(results, args.outdir)
    summary = br.summarize(results, args.rmsd_cutoff)
    summary["protocol"] = {
        "method": f"baseline: standard Vina workflow ({args.engine})",
        "selection": ("pdbids" if args.pdbids else
                      f"first {args.first}" if args.first else
                      f"sample {args.sample} seed {args.sample_seed} of "
                      f"{os.path.basename(args.ids_csv)}"),
        "pdbids": pdbids,
        **{k: v for k, v in baseline_params(opts).items() if k != "tool_versions"},
        "receptor_hetatm": "kept only if Meeko has a template (else deleted with -x)",
        "preparation": PREPARATION_CHAIN,
        "rmsd": "obrms heavy-atom, symmetry-corrected, no superposition",
        "top1": "mode 1 (first MODEL)",
        "tools": tool_versions,
        "vina_version": tool_versions["vina"]["version"],
    }
    with open(os.path.join(args.outdir, "summary.json"), "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)
    print("\n=== BASELINE REDOCK SUMMARY ===")
    for k, v in summary.items():
        if k not in ("excluded", "protocol"):
            print(f"  {k}: {v}")
    return 0


# ---------------------------------------------------------------------------
# Enrichment mode
# ---------------------------------------------------------------------------

def target_seed(target):
    """benchmark_enrichment.py's per-target seed (actives; decoys use
    ``seed ^ DECOY_SEED_SALT``)."""
    return int(hashlib.md5(target.encode()).hexdigest()[:8], 16)


def read_ligand_table(path):
    with open(path, newline="", encoding="utf-8") as f:
        return [{"ligand": r["ligand"], "compound": r["compound"],
                 "entry_index": int(r["entry_index"]), "n_entries": int(r["n_entries"])}
                for r in csv.DictReader(f)]


def select_entries(mol2_gz, prefix, max_n, seed, convert, table=None):
    """The compounds benchmark_enrichment.py docks, each as an SDF.

    With ``table`` (rows of the AutoDocker run's ligands.csv with this
    prefix) exactly those entries are used, and a compound RDKit cannot
    read stays in the set as a failure (ranked last). Otherwise the
    selection is repeated: compounds grouped by ID, IDs sorted and shuffled
    with ``random.Random(seed)``, first entry per compound, the next
    compound replacing one that fails to convert.

    ``convert(block, name)`` writes the SDF and returns its path. Returns
    ``(selected, failures, info)``; ``selected`` rows carry ``sdf`` (None
    when the conversion failed), ``info`` the entry/compound counts of the
    file."""
    blocks = be.read_molecules(mol2_gz)
    groups = be.group_compounds(blocks)
    selected, failures = [], []
    info = {"file_entries": len(blocks), "file_compounds": len(groups)}
    if table is not None:
        for row in table:
            block = blocks[row["entry_index"]]
            if be.molecule_name(block) != row["compound"]:
                raise RuntimeError(f"{row['ligand']}: entry {row['entry_index']} of "
                                   f"{os.path.basename(mol2_gz)} is not {row['compound']}")
            try:
                sdf = convert(block, row["ligand"])
            except Exception as e:
                failures.append({"ligand": row["ligand"], "error": str(e)[-300:]})
                sdf = None
            selected.append({**row, "sdf": sdf})
        return selected, failures, info
    order = sorted(groups)
    random.Random(seed).shuffle(order)
    for compound in order:
        if len(selected) >= max_n:
            break
        idx = groups[compound][0]
        stem = be._file_stem(prefix, compound)
        try:
            sdf = convert(blocks[idx], stem)
        except Exception as e:
            failures.append({"ligand": f"{stem} (entry {idx})", "error": str(e)[-300:]})
            continue
        selected.append({"ligand": stem, "compound": compound, "entry_index": idx,
                         "n_entries": len(groups[compound]), "sdf": sdf})
    return selected, failures, info


def _dock_one_ligand(job):
    name, sdf, lig_dir, dock_dir, receptor, center, size, opts = job
    pdbqt = os.path.join(lig_dir, name + ".pdbqt")
    out = os.path.join(dock_dir, name + _benchlib.POSE_SUFFIX)
    try:
        prepare_ligand_pdbqt(sdf, pdbqt, opts["tools"], os.path.join(lig_dir, name + ".log"))
        dock_vina(opts["tools"], receptor, pdbqt, out,
                  os.path.join(dock_dir, name + _benchlib.VINA_LOG_SUFFIX), center, size, opts)
        return name, out, None
    except Exception as e:
        return name, None, str(e)[-300:]


def write_ranking(path, poses, errors, engine):
    """A ranking.csv in AutoDocker's format (status rules of the pipeline:
    OK/WEAK in the window, OUT_OF_WINDOW outside, FAILED without a pose)."""
    rows = []
    for name in sorted(set(poses) | set(errors)):
        pose = poses.get(name)
        token, value = _benchlib.first_vina_result(pose) if pose else (None, None)
        if value is None:
            rows.append((2, float("inf"), name, "FAILED", "FAILED", ""))
            continue
        status = vr.expected_status(value, vr.DEFAULT_MIN_VALID_AFFINITY)
        rows.append((vr.STATUS_GROUP[status], value, name, token, status,
                     os.path.relpath(pose, os.path.dirname(path))))
    rows.sort()
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=vr.RANKING_COLUMNS)
        writer.writeheader()
        for _, _, name, token, status, docked in rows:
            writer.writerow({"Ligand": name, "Binding_Affinity": token, "Docked_File": docked,
                             "Status": status, "Engine": engine})
    return path


def enrichment_target(target, opts):
    """One DUD-E target; same record layout as benchmark_enrichment."""
    rec = {"target": target, "family": be.TARGET_FAMILY.get(target, "?"),
           "status": "error", "n_actives": 0, "n_decoys": 0, "box_source": "crystal",
           "error": None}
    tdir = os.path.join(opts["workdir"], target)
    work = os.path.join(opts["outdir"], "work", target)
    tools = opts["tools"]
    try:
        receptor_pdb = os.path.join(tdir, be.DUD_E_RECEPTOR)
        crystal = os.path.join(tdir, be.DUD_E_CRYSTAL_LIGAND)
        for path in (receptor_pdb, crystal, os.path.join(tdir, be.DUD_E_ACTIVES),
                     os.path.join(tdir, be.DUD_E_DECOYS)):
            if not os.path.exists(path):
                raise RuntimeError(f"missing {path} (run benchmark_enrichment.py "
                                   "--prepare-only first)")
        shutil.rmtree(work, ignore_errors=True)
        sdf_dir = os.path.join(work, "ligands_sdf")
        lig_dir = os.path.join(work, "ligands")
        dock_dir = os.path.join(work, "docked")
        for d in (sdf_dir, lig_dir, dock_dir):
            os.makedirs(d)

        def convert(block, name):
            return mol2_block_to_sdf(block, os.path.join(sdf_dir, name + ".sdf"), name)

        table = None
        if opts.get("ligands_from"):
            path = os.path.join(opts["ligands_from"], target, "ligands.csv")
            if not os.path.exists(path):
                raise RuntimeError(f"--ligands-from: {path} not found")
            table = read_ligand_table(path)
            rec["ligand_source"] = os.path.abspath(path)
        seed = target_seed(target)
        actives, fail_a, info_a = select_entries(
            os.path.join(tdir, be.DUD_E_ACTIVES), "A_", opts["max_actives"], seed, convert,
            None if table is None else [r for r in table if r["ligand"].startswith("A_")])
        decoys, fail_d, info_d = select_entries(
            os.path.join(tdir, be.DUD_E_DECOYS), "D_", opts["max_decoys"],
            seed ^ be.DECOY_SEED_SALT, convert,
            None if table is None else [r for r in table if r["ligand"].startswith("D_")])
        rec["n_actives"], rec["n_decoys"] = len(actives), len(decoys)
        rec["n_convert_failed_actives"], rec["n_convert_failed_decoys"] = len(fail_a), len(fail_d)
        rec["convert_failures"] = fail_a + fail_d
        for label, info in (("actives", info_a), ("decoys", info_d)):
            rec[f"dude_{label}_entries"] = info["file_entries"]
            rec[f"dude_{label}_compounds"] = info["file_compounds"]
        selected = actives + decoys
        be.write_ligand_table(target, opts["outdir"],
                              [{k: r[k] for k in ("ligand", "compound", "entry_index",
                                                  "n_entries")} for r in selected])
        if not actives or not decoys:
            raise RuntimeError("no usable actives/decoys")
        center, size = be.crystal_box(crystal, opts["padding"], opts["min_box"])
        rec["box_center"] = [round(c, 3) for c in center]
        rec["box_size"] = round(size, 3)
        receptor, rec["receptor_preparation"] = prepare_receptor(
            receptor_pdb, os.path.join(work, "receptor"), tools)

        dockable = [r for r in selected if r["sdf"]]
        errors = {r["ligand"]: "ligand conversion failed" for r in selected if not r["sdf"]}
        poses = {}
        if opts["engine"] == "easydock":
            combined = os.path.join(work, "ligands_all.sdf")
            with open(combined, "w", encoding="utf-8") as out:
                for r in dockable:
                    with open(r["sdf"], encoding="utf-8") as f:
                        out.write(f.read())
            got = dock_easydock(tools, receptor, combined, os.path.join(work, "easydock"),
                                center, size, opts, ncpu=opts["processes"])
            for r in dockable:
                if r["ligand"] in got:
                    poses[r["ligand"]] = got[r["ligand"]]
                else:
                    errors[r["ligand"]] = "EasyDock returned no pose"
        else:
            jobs = [(r["ligand"], r["sdf"], lig_dir, dock_dir, receptor, center, size, opts)
                    for r in dockable]
            for name, pose, err in _map(_dock_one_ligand, jobs, opts["processes"]):
                if pose:
                    poses[name] = pose
                else:
                    errors[name] = err
        ranking = write_ranking(os.path.join(work, _benchlib.RANKING_CSV), poses, errors,
                                "Vina" if opts["engine"] == "vina" else "EasyDock")
        expected = [r["ligand"] for r in selected]
        names, scores, info = be.collect_scores(ranking, expected, opts["out_of_window"])
        rec["n_out_of_window"] = sum(1 for s in info["status"].values() if s == "OUT_OF_WINDOW")
        rec["n_missing_from_ranking"] = sum(1 for s in info["status"].values()
                                            if s == "MISSING")
        rec["docking_errors"] = [{"ligand": k, "error": v} for k, v in sorted(errors.items())]
        rec["metrics"] = be.write_metrics(target, opts["outdir"], names, scores,
                                          info["status"], ranking)
        stats = be.compute_stats(names, scores, seed, opts["bootstrap"], opts["ci"])
        stats["warnings"] = rec.get("warnings", []) + stats["warnings"]
        rec.update(stats)
        rec["status"] = "ok"
    except Exception as e:
        rec["error"] = str(e)
    return rec


def run_enrichment(args, tools, tool_versions):
    targets = [t.strip() for t in args.targets.split(",") if t.strip()]
    opts = {"engine": args.engine, "workdir": args.workdir, "outdir": args.outdir,
            "max_actives": args.max_actives, "max_decoys": args.max_decoys,
            "ligands_from": args.ligands_from, "padding": args.padding,
            "min_box": args.min_box, "seed": args.seed,
            "exhaustiveness": args.exhaustiveness, "binding_modes": ENRICHMENT_MODES,
            "energy_range": ENRICHMENT_ENERGY_RANGE, "threads": args.threads,
            "processes": args.processes, "timeout": args.timeout,
            "bootstrap": args.bootstrap, "ci": args.ci,
            "out_of_window": args.out_of_window, "tools": tools}
    protocol = {
        "method": f"baseline: standard Vina workflow ({args.engine})",
        "targets": targets, "box_source": "crystal", "padding": args.padding,
        "min_box": args.min_box, "max_actives": args.max_actives,
        "max_decoys": args.max_decoys,
        "ligand_selection": ("the AutoDocker run's <target>/ligands.csv from "
                             f"{os.path.abspath(args.ligands_from)}" if args.ligands_from else
                             "benchmark_enrichment.py's seeded compound sample, repeated"),
        "sampling_unit": "compound (DUD-E compound ID), not MOL2 entry",
        "entry_per_compound": be.ENTRY_RULE,
        "primary_metrics": list(be.PRIMARY_METRICS),
        "secondary_metrics_caveats": be.SECONDARY_METRIC_CAVEATS,
        "seed": args.seed, "exhaustiveness": args.exhaustiveness,
        "binding_modes": ENRICHMENT_MODES, "energy_range": ENRICHMENT_ENERGY_RANGE,
        "threads": args.threads, "processes": args.processes,
        "bootstrap": args.bootstrap, "ci": args.ci, "bedroc_alpha": be.BEDROC_ALPHA,
        "out_of_window": args.out_of_window,
        "failed_ligands": "kept in N, ranked last (tied)",
        "preparation": PREPARATION_CHAIN, "tools": tool_versions,
        "vina_version": tool_versions["vina"]["version"],
    }
    results = []
    for t in targets:
        print(f"\n=== [{t}] baseline ({args.engine}) ===", flush=True)
        rec = enrichment_target(t, opts)
        results.append(rec)
        print(f"[{t}] {rec['status']} actives={rec['n_actives']} decoys={rec['n_decoys']} "
              f"auc={rec.get('auc')} ef5={rec.get('ef5')} err={rec['error']}", flush=True)
        be.save_tables(results, args.outdir, protocol)
    be.save_tables(results, args.outdir, protocol)
    return 0 if all(r["status"] == "ok" for r in results) else 1


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

PREPARATION_CHAIN = {
    "receptor": ["fill blank PDB element columns from atom names", "reduce -Trim",
                 "reduce -BUILD [-DB het dict]",
                 "mk_prepare_receptor.py -i receptorH.pdb -o receptor -p "
                 "(re-run with -x only if residues match no template)"],
    "ligand": ["RDKit molecule with explicit hydrogens -> SDF",
               "mk_prepare_ligand.py -i ligand.sdf -o ligand.pdbqt"],
    "docking": "vina --receptor --ligand --out --center_* --size_* --exhaustiveness "
               "--num_modes --energy_range --seed --cpu",
}


def _versions_only(tool_versions):
    """The tool facts that change results (for the resume fingerprint)."""
    out = {k: [v.get("version"), v.get("sha256")] for k, v in tool_versions.items()
           if isinstance(v, dict) and "sha256" in v}
    out["meeko_python"] = tool_versions.get("meeko_python", {}).get("packages")
    if "easydock" in tool_versions:
        out["easydock"] = tool_versions["easydock"].get("packages")
    return out


def write_unavailable(outdir, engine, reason, tool_versions):
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "summary.json"), "w", encoding="utf-8") as f:
        json.dump({"engine": engine, "available": False, "reason": reason,
                   "tools": tool_versions}, f, indent=2)


def _add_common(ap, defaults):
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--engine", choices=ENGINES, default="vina",
                    help="raw Vina (default) or EasyDock (Vina engine)")
    ap.add_argument("--vina-bin", default=None, help="Vina binary (default: vina on PATH)")
    ap.add_argument("--meeko-bin", default=None,
                    help="Directory with mk_prepare_receptor.py/mk_prepare_ligand.py (and "
                         "reduce), e.g. a virtual environment's bin/ (default: PATH)")
    ap.add_argument("--reduce-bin", default=None, help="reduce binary (default: --meeko-bin "
                    "or PATH)")
    ap.add_argument("--reduce-het-dict", default=None,
                    help="reduce het dictionary (default: the reduce-binary wheel's)")
    ap.add_argument("--easydock-bin", default=None,
                    help="easydock executable for --engine easydock (default: PATH)")
    ap.add_argument("--seed", type=int, default=_benchlib.PIPELINE_DEFAULTS.SEED)
    ap.add_argument("--exhaustiveness", type=int, default=defaults["exhaustiveness"])
    ap.add_argument("--padding", type=float, default=_benchlib.BENCH_BOX_PADDING)
    ap.add_argument("--min-box", type=float, default=_benchlib.BENCH_MIN_BOX_EDGE)
    ap.add_argument("--threads", type=int, default=1,
                    help="Vina --cpu per docking (default 1)")
    ap.add_argument("--processes", type=int, default=_benchlib.available_cpus(),
                    help="Parallel dockings (default: available CPUs)")
    ap.add_argument("--timeout", type=int, default=defaults["timeout"])


def build_parser():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="mode", required=True)
    pipeline = _benchlib.PIPELINE_DEFAULTS

    rd = sub.add_parser("redock", help="CASF re-docking (benchmark_redock.py's set)")
    rd.add_argument("--data-dir", required=True,
                    help="fetch_casf_rcsb.py output (one sub-dir per PDB ID)")
    rd.add_argument("--pdbids", help="Comma-separated PDB IDs (overrides --sample)")
    rd.add_argument("--sample", type=int, default=_benchlib.BENCH_SAMPLE_SIZE)
    rd.add_argument("--sample-seed", type=int, default=_benchlib.BENCH_SEED)
    rd.add_argument("--ids-csv", default=br.DEFAULT_IDS)
    rd.add_argument("--first", type=int, default=None)
    rd.add_argument("--binding-modes", type=int, default=pipeline.BINDING_MODES)
    rd.add_argument("--energy-range", type=float, default=pipeline.ENERGY_RANGE)
    rd.add_argument("--start-conformer", choices=["rdkit", "crystal"], default="rdkit")
    rd.add_argument("--conformer-seed", type=int, default=_benchlib.BENCH_SEED)
    rd.add_argument("--rmsd-cutoff", type=float, default=2.0)
    _add_common(rd, {"exhaustiveness": pipeline.EXHAUSTIVENESS,
                     "timeout": _benchlib.DOCKING_RUN_TIMEOUT})

    en = sub.add_parser("enrichment", help="DUD-E enrichment (benchmark_enrichment.py's set)")
    en.add_argument("--workdir", required=True,
                    help="benchmark_enrichment.py --workdir (DUD-E files per target)")
    en.add_argument("--targets", default=",".join(be.DEFAULT_TARGETS))
    en.add_argument("--ligands-from", default=None, metavar="DIR",
                    help="AutoDocker enrichment --outdir: dock exactly the compounds of "
                         "its <target>/ligands.csv (recommended)")
    en.add_argument("--max-actives", type=int, default=be.DEFAULT_MAX_ACTIVES)
    en.add_argument("--max-decoys", type=int, default=be.DEFAULT_MAX_DECOYS)
    en.add_argument("--bootstrap", type=int, default=1000)
    en.add_argument("--ci", type=float, default=be.BOOTSTRAP_CI)
    en.add_argument("--out-of-window", choices=["score", "fail"], default="score")
    _add_common(en, {"exhaustiveness": ENRICHMENT_EXHAUSTIVENESS,
                     "timeout": _benchlib.DOCKING_RUN_TIMEOUT})
    return ap


def main(argv=None):
    args = build_parser().parse_args(argv)
    os.makedirs(args.outdir, exist_ok=True)
    tools = resolve_tools(args)
    tool_versions = tool_record(tools)
    if args.engine == "easydock":
        ok, reason = easydock_available(tools)
        if not ok:
            write_unavailable(args.outdir, "easydock", reason, tool_versions)
            print(f"[!] --engine easydock not available: {reason}", file=sys.stderr)
            return EXIT_ENGINE_UNAVAILABLE
    missing = missing_tools(tools, args.engine)
    if missing:
        reason = "missing tools: " + ", ".join(missing)
        write_unavailable(args.outdir, args.engine, reason, tool_versions)
        print(f"[!] {reason}", file=sys.stderr)
        return EXIT_ENGINE_UNAVAILABLE
    if args.mode == "redock":
        return run_redock(args, tools, tool_versions)
    return run_enrichment(args, tools, tool_versions)


if __name__ == "__main__":
    sys.exit(main())

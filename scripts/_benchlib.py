"""Helpers shared by the benchmark scripts and their tests.

The benchmark scripts import this module as ``_benchlib`` (each script puts
its own directory on ``sys.path`` first), so it works both when a script is
run directly and when a test loads a script with ``importlib``.

Only standard-library code lives here; optional dependencies (RDKit) are
imported inside the functions that need them.

The scripts run standalone (they need not import the pipeline package), so
the Vina/SMINA output parser below is a copy of the pipeline's rules
(``vspipeline.docking.parse_pdbqt_results``): the score of a model is its
first ``REMARK minimizedAffinity`` (SMINA) if present, else its first
``REMARK VINA RESULT``. tests/test_benchmarks.py checks that both parsers
give the same scores.
"""

import csv
import math
import os
import random
import re
import shutil
import subprocess
import sys

# Names of AutoDocker's output files read by the benchmarks. They mirror
# autodocker/vspipeline/outputs.py (the scripts also run outside the package,
# e.g. from /repo/scripts in the container); tests/test_benchmarks.py checks
# that the two agree.
RANKING_CSV = "ranking.csv"
DOCKED_DIR = "docked"
FRAME_OFFSET_TXT = "frame_offset.txt"
COL_LIGAND = "Ligand"
COL_AFFINITY = "Binding_Affinity"
COL_STATUS = "Status"
COL_DOCKED_FILE = "Docked_File"
COL_ENGINE = "Engine"

# Timeouts (seconds) for the external helpers and downloads used by the
# benchmark scripts.
TOOL_TIMEOUT = 120          # obrms, obabel conversions
VERSION_TIMEOUT = 30        # `<tool> --version`
RUNNER_VERSION_TIMEOUT = 60  # `runner.py --version` (imports the pipeline)
DOWNLOAD_TIMEOUT = 120      # one HTTP download (RCSB)


def available_cpus():
    """CPUs this process may use (respects affinity/cgroup pinning); the
    default number of parallel workers of the benchmark scripts."""
    try:
        return max(1, len(os.sched_getaffinity(0)))
    except (AttributeError, OSError):
        return max(1, os.cpu_count() or 1)


# ---------------------------------------------------------------------------
# Docking frame
# ---------------------------------------------------------------------------


def read_frame_offset(path):
    """Return the (x, y, z) offset the runner wrote to ``frame_offset.txt``.

    ``path`` is the file itself or the pipeline output directory holding it.
    Poses are in the docking frame; input frame = docking frame + offset.
    """
    if os.path.isdir(path):
        path = os.path.join(path, FRAME_OFFSET_TXT)
    values = {}
    with open(path, encoding="utf-8") as f:
        for line in f:
            if "=" in line and not line.lstrip().startswith("#"):
                key, val = line.split("=", 1)
                values[key.strip()] = float(val)
    return values["offset_x"], values["offset_y"], values["offset_z"]


# ---------------------------------------------------------------------------
# MOL2
# ---------------------------------------------------------------------------


def parse_mol2_atoms(mol2_path):
    """Atoms of the first molecule in a MOL2 file.

    Returns a list of dicts with ``name``, ``x``, ``y``, ``z``, ``type``
    (SYBYL type) and ``element`` (the part of the type before the dot).
    """
    atoms, in_atoms, n_molecules = [], False, 0
    with open(mol2_path, encoding="utf-8") as f:
        for line in f:
            stripped = line.strip()
            if stripped.startswith("@<TRIPOS>MOLECULE"):
                n_molecules += 1
                if n_molecules > 1:
                    break
                continue
            if stripped.startswith("@<TRIPOS>ATOM"):
                in_atoms = True
                continue
            if stripped.startswith("@<TRIPOS>"):
                in_atoms = False
                continue
            if in_atoms and stripped:
                parts = stripped.split()
                if len(parts) < 6:
                    continue
                atype = parts[5]
                atoms.append({
                    "name": parts[1],
                    "x": float(parts[2]), "y": float(parts[3]),
                    "z": float(parts[4]),
                    "type": atype,
                    "element": atype.split(".")[0],
                })
    return atoms


def mol2_heavy_coords(mol2_path):
    """Coordinates of the non-hydrogen atoms of a MOL2 ligand."""
    coords = [(a["x"], a["y"], a["z"]) for a in parse_mol2_atoms(mol2_path)
              if a["element"].upper() not in ("H", "D")]
    if not coords:
        raise ValueError(f"no heavy atoms in {mol2_path}")
    return coords


def shift_mol2(mol2_path, shift, out_path):
    """Copy a MOL2 file with every ATOM coordinate translated by ``-shift``.

    Used to bring a crystal ligand (input frame) into the docking frame.
    """
    with open(mol2_path, encoding="utf-8") as src, open(out_path, "w", encoding="utf-8") as dst:
        in_atoms = False
        for line in src:
            stripped = line.strip()
            if stripped.startswith("@<TRIPOS>"):
                in_atoms = stripped.startswith("@<TRIPOS>ATOM")
                dst.write(line)
                continue
            if in_atoms and stripped:
                parts = line.split()
                x = float(parts[2]) - shift[0]
                y = float(parts[3]) - shift[1]
                z = float(parts[4]) - shift[2]
                parts[2], parts[3], parts[4] = f"{x:.4f}", f"{y:.4f}", f"{z:.4f}"
                dst.write(" ".join(parts) + "\n")
            else:
                dst.write(line)


# ---------------------------------------------------------------------------
# Docking box
# ---------------------------------------------------------------------------

# Benchmark protocol for boxes around a reference (crystal) ligand: edge =
# max(BENCH_MIN_BOX_EDGE, ligand extent + 2 x BENCH_BOX_PADDING) per axis,
# i.e. max(22.5, extent + 12) Angstrom. This is the benchmarks' own rule; it
# is unrelated to the pipeline's --padding around fpocket pockets (also 6 A).
BENCH_BOX_PADDING = 6.0
BENCH_MIN_BOX_EDGE = 22.5


def ligand_box(coords, padding=BENCH_BOX_PADDING, min_size=BENCH_MIN_BOX_EDGE):
    """Cubic docking box around a reference ligand.

    Center = heavy-atom centroid; edge = max(min_size, largest extent along
    x/y/z + 2*padding). With the defaults this is max(22.5, extent + 12) A.
    Returns ``((cx, cy, cz), edge)``.
    """
    coords = list(coords)
    if not coords:
        raise ValueError("no coordinates for the docking box")
    n = float(len(coords))
    center = tuple(sum(c[i] for c in coords) / n for i in range(3))
    extent = max(max(c[i] for c in coords) - min(c[i] for c in coords)
                 for i in range(3))
    return center, max(min_size, extent + 2.0 * padding)


def box_cli_args(center, size):
    """Runner arguments for an explicit box given in input-PDB coordinates."""
    return [
        "--no-fpocket",
        f"--grid-center={center[0]:.3f},{center[1]:.3f},{center[2]:.3f}",
        "--grid-center-frame=input",
        f"--grid-size={size:.3f},{size:.3f},{size:.3f}",
    ]


# ---------------------------------------------------------------------------
# Vina / PDBQT output
# ---------------------------------------------------------------------------

VINA_RESULT_RE = re.compile(
    r"^REMARK VINA RESULT:\s+(\S+)\s+(\S+)\s+(\S+)")
MINIMIZED_AFFINITY_RE = re.compile(r"^REMARK minimizedAffinity\s+(\S+)")
_HYDROGEN_TYPES = ("H", "HD", "HS")


def _float_or_none(text):
    try:
        value = float(text)
    except (TypeError, ValueError):
        return None
    return value if math.isfinite(value) else None


def parse_vina_models(pdbqt_path):
    """Every MODEL of a Vina/SMINA output file, in file order.

    Returns a list of dicts: ``model`` (1-based position in the file),
    ``affinity`` (float or None), ``affinity_str`` (the token as printed:
    the model's first ``REMARK minimizedAffinity`` when present, which is the
    score SMINA computed, else its first ``REMARK VINA RESULT``; None when
    the model has neither),
    ``lines`` (the ligand lines of the model, without MODEL/ENDMDL and
    without flexible-residue blocks) and ``heavy_coords``. A file without
    MODEL records is treated as a single model.
    """
    models, current, in_flex = [], None, False
    seen_model = False

    def new_model():
        return {"model": len(models) + 1, "affinity": None,
                "affinity_str": None, "lines": [], "heavy_coords": [],
                "_minimized": False}

    def close(model):
        model.pop("_minimized", None)
        models.append(model)

    with open(pdbqt_path, encoding="utf-8") as f:
        for raw in f:
            line = raw.rstrip("\n")
            if line.startswith("MODEL"):
                current, in_flex = new_model(), False
                seen_model = True
                continue
            if line.startswith("ENDMDL"):
                if current is not None:
                    close(current)
                current = None
                continue
            if current is None:
                if seen_model or not line.strip():
                    continue  # text between/after MODEL blocks
                current = new_model()  # file without MODEL records
            if line.startswith("BEGIN_RES"):
                in_flex = True
                continue
            if line.startswith("END_RES"):
                in_flex = False
                continue
            if in_flex:
                continue
            m = MINIMIZED_AFFINITY_RE.match(line)
            if m and not current["_minimized"]:
                current["_minimized"] = True
                current["affinity_str"] = m.group(1)
                current["affinity"] = _float_or_none(m.group(1))
            m = VINA_RESULT_RE.match(line)
            if m and current["affinity_str"] is None:
                current["affinity_str"] = m.group(1)
                current["affinity"] = _float_or_none(m.group(1))
            current["lines"].append(line)
            if line.startswith(("ATOM", "HETATM")):
                atype = line.split()[-1] if line.split() else ""
                if atype not in _HYDROGEN_TYPES:
                    try:
                        current["heavy_coords"].append(
                            (float(line[30:38]), float(line[38:46]),
                             float(line[46:54])))
                    except ValueError:
                        pass
    if current is not None and not seen_model:
        close(current)
    return models


def first_vina_result(pdbqt_path):
    """``(token, value)`` of the first scored model (see parse_vina_models),
    or ``(None, None)``. For Vina output this is the mode-1 score exactly as
    Vina printed it, the value the pipeline writes to ranking.csv."""
    for model in parse_vina_models(pdbqt_path):
        if model["affinity_str"] is not None:
            return model["affinity_str"], model["affinity"]
    return None, None


def write_model(model, path):
    """Write one parsed model as a standalone PDBQT file."""
    with open(path, "w", encoding="utf-8") as out:
        out.write("\n".join(model["lines"]) + "\n")
    return path


# ---------------------------------------------------------------------------
# ranking.csv
# ---------------------------------------------------------------------------

AFFINITY_COLUMNS = (COL_AFFINITY, "Affinity")  # "Affinity": other tools' CSVs
# Statuses that mean "no usable score", whatever the affinity column holds.
# Status values without a usable score: FAILED (pipeline Status.FAILED), the
# pre-2.0 REJECTED, and ERROR from other tools' CSVs. Checked against the
# pipeline in tests/test_benchmarks.py.
FAILED_STATUS_PREFIXES = ("FAILED", "REJECTED", "ERROR")


def ligand_key(name):
    """Ligand name as used for matching (file extension dropped)."""
    name = (name or "").strip()
    for ext in (".pdbqt", ".mol2", ".sdf", ".pdb"):
        if name.lower().endswith(ext):
            return name[:-len(ext)]
    return name


def read_ranking(ranking_csv):
    """Rows of a pipeline (or external) ranking.csv, read by column name.

    Works with both the 1.1.x format (``FAILED`` in Binding_Affinity,
    absolute Docked_File) and newer ones (numeric out-of-window scores with
    Status ``OUT_OF_WINDOW``, relative Docked_File, an Engine column).

    Each row is a dict with ``name``, ``affinity_str`` (as written),
    ``affinity`` (float, or None when the value is not a finite number or
    the status marks a failure), ``status`` (upper case, may be ""),
    ``docked_file`` and ``engine``.
    """
    rows = []
    with open(ranking_csv, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        fields = reader.fieldnames or []
        key = next((k for k in AFFINITY_COLUMNS if k in fields), None)
        if key is None or COL_LIGAND not in fields:
            raise ValueError(
                f"{ranking_csv}: need a Ligand column and one of "
                f"{', '.join(AFFINITY_COLUMNS)} (found {fields})")
        for row in reader:
            text = (row.get(key) or "").strip()
            status = (row.get(COL_STATUS) or "").strip().upper()
            affinity = _float_or_none(text)
            if status.startswith(FAILED_STATUS_PREFIXES):
                affinity = None
            rows.append({
                "name": ligand_key(row.get(COL_LIGAND)),
                "affinity_str": text,
                "affinity": affinity,
                "status": status,
                "docked_file": (row.get(COL_DOCKED_FILE) or "").strip(),
                "engine": (row.get(COL_ENGINE) or "").strip(),
            })
    return rows


def resolve_docked_file(ranking_csv, docked_file):
    """Absolute path of a Docked_File entry (absolute, or relative to the
    directory holding ranking.csv). Falls back to ``docked/<basename>`` next
    to ranking.csv when a recorded absolute path no longer exists (e.g. a
    run made inside a container)."""
    if not docked_file:
        return None
    base = os.path.dirname(os.path.abspath(ranking_csv))
    path = docked_file if os.path.isabs(docked_file) else os.path.join(base, docked_file)
    if os.path.exists(path):
        return path
    alt = os.path.join(base, DOCKED_DIR, os.path.basename(docked_file))
    return alt if os.path.exists(alt) else path


# ---------------------------------------------------------------------------
# External tools
# ---------------------------------------------------------------------------


def obrms(ref, query, timeout=TOOL_TIMEOUT):
    """Heavy-atom, symmetry-aware RMSD (Open Babel ``obrms``) without
    superposition. Raises RuntimeError when obrms fails or prints no value."""
    proc = subprocess.run(["obrms", ref, query], capture_output=True,
                          text=True, timeout=timeout)
    if proc.returncode != 0:
        raise RuntimeError(f"obrms failed: {proc.stderr.strip()[-300:]}")
    for line in proc.stdout.splitlines():
        if line.startswith("RMSD"):
            value = _float_or_none(line.split()[-1])
            if value is not None:
                return value
    raise RuntimeError(f"no RMSD in obrms output: {proc.stdout[-300:]!r}")


def tool_version(binary, args=("--version",), timeout=VERSION_TIMEOUT):
    """First non-empty output line of ``binary --version`` (or None)."""
    path = shutil.which(binary) if binary and not os.path.isabs(binary) else binary
    if not path or not os.path.exists(path):
        return None
    try:
        proc = subprocess.run([path, *args], capture_output=True, text=True,
                              timeout=timeout)
    except (OSError, subprocess.TimeoutExpired):
        return None
    for line in (proc.stdout + "\n" + proc.stderr).splitlines():
        if line.strip():
            return line.strip()
    return None


def runner_version(runner_script, timeout=RUNNER_VERSION_TIMEOUT):
    """``AutoDocker x.y.z`` as printed by ``runner.py --version``."""
    import sys
    try:
        proc = subprocess.run([sys.executable, runner_script, "--version"],
                              capture_output=True, text=True, timeout=timeout)
    except (OSError, subprocess.TimeoutExpired):
        return None
    out = (proc.stdout or proc.stderr).strip()
    return out.splitlines()[-1] if out else None


# ---------------------------------------------------------------------------
# Subset selection
# ---------------------------------------------------------------------------


def select_sample(ids, n, seed):
    """Seeded random sample of ``n`` IDs, returned sorted.

    The input order does not matter (IDs are sorted and de-duplicated first),
    so the same (n, seed) always selects the same IDs. ``n`` <= 0 or >= the
    number of IDs selects all of them.
    """
    pool = sorted(set(ids))
    if n is None or n <= 0 or n >= len(pool):
        return pool
    return sorted(random.Random(seed).sample(pool, n))


# ---------------------------------------------------------------------------
# Start conformers (RDKit)
# ---------------------------------------------------------------------------

def _pipeline_module(name):
    """A module of the pipeline package (autodocker/vspipeline/<name>.py), so
    the benchmarks use the pipeline's own defaults and conformer recipe. The
    scripts live next to autodocker/ in the repository (also in the image's
    /repo mount)."""
    root = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "autodocker")
    if root not in sys.path:
        sys.path.insert(0, root)
    import importlib
    return importlib.import_module(f"vspipeline.{name}")


def _pipeline_conformers():
    return _pipeline_module("conformers")


# The pipeline's defaults (vspipeline/defaults.py): docking settings that the
# benchmarks use unless told otherwise.
PIPELINE_DEFAULTS = _pipeline_module("defaults")

# Seed of the benchmarks' own random choices (which complexes, ligands and
# start conformers are used); recorded in every summary.json.
BENCH_SEED = 42
# Timeouts (seconds) of whole benchmark steps.
DOCKING_RUN_TIMEOUT = 1800        # one re-docked complex / one parity re-run
ENRICHMENT_RUN_TIMEOUT = 7200     # docking one DUD-E target's ligand set


def rdkit_start_conformer(src, out_mol2, seed, mmff_iters=None):
    """Write a crystal-independent start conformer for re-docking.

    The molecule (topology, bond orders, formal charges) is read from
    ``src`` (SDF preferred; MOL2 is converted to SDF with Open Babel first),
    its coordinates are discarded and one ETKDGv3 conformer is embedded with
    a fixed ``randomSeed`` on a single thread, then MMFF94-optimised with a
    fixed iteration cap (UFF when MMFF lacks parameters). The result is
    deterministic for a given RDKit version. Written as MOL2 via Open Babel.

    Returns a metadata dict (seed, RDKit version, force field, converged).
    Raises ImportError when RDKit is missing, RuntimeError on failure.
    """
    from rdkit import Chem, rdBase

    conformers = _pipeline_conformers()
    if mmff_iters is None:
        mmff_iters = conformers.FORCE_FIELD_MAX_ITERS

    tmp_dir = os.path.dirname(os.path.abspath(out_mol2))
    sdf = src
    if not src.lower().endswith(".sdf"):
        sdf = os.path.join(tmp_dir, os.path.basename(src) + ".tmp.sdf")
        subprocess.run(["obabel", src, "-osdf", "-O", sdf], check=True,
                       capture_output=True, timeout=TOOL_TIMEOUT)
    mol = Chem.MolFromMolFile(sdf, removeHs=False, sanitize=True)
    if mol is None:
        raise RuntimeError(f"RDKit could not read {sdf}")
    mol = Chem.RemoveHs(mol)
    mol = Chem.AddHs(mol)
    mol.RemoveAllConformers()
    try:
        ff, converged = conformers.embed_and_relax(mol, seed=int(seed), max_iters=mmff_iters)
    except ValueError as e:
        raise RuntimeError(str(e)) from None
    tmp_sdf = out_mol2 + ".tmp.sdf"
    with open(tmp_sdf, "w", encoding="utf-8") as f:
        f.write(Chem.MolToMolBlock(mol))
        f.write("$$$$\n")
    try:
        subprocess.run(["obabel", "-isdf", tmp_sdf, "-omol2", "-O", out_mol2],
                       check=True, capture_output=True, timeout=TOOL_TIMEOUT)
    finally:
        os.remove(tmp_sdf)
        if sdf != src and os.path.exists(sdf):
            os.remove(sdf)
    return {"method": "ETKDGv3+" + ff, "seed": int(seed),
            "max_iters": mmff_iters, "converged": converged,
            "rdkit": rdBase.rdkitVersion}

#!/usr/bin/env python3
"""Shared helpers: subprocess execution, Open Babel calls, chain parsing.

The tool paths and the command timeout come from ``tools.settings()``, read
at call time. Every external command of the pipeline goes through ``run``
(the Open Babel calls through ``run_obabel``), so a test can replace
``utils.run`` to fake them.
"""
import hashlib
import logging
import os
import re
import shlex
import subprocess
from functools import lru_cache
from multiprocessing import cpu_count as mp_cpu_count
from typing import Dict, Iterable, List, Optional, Tuple, Union

from vspipeline import defaults, tools

logger = logging.getLogger(__name__)

# Chain ID assumed for a PDB whose ATOM records carry no chain identifier
# (blank column 22, common in DUD-E-style receptors).
DEFAULT_CHAIN = "A"

# Open Babel's force fields and conformer search use OpenMP; a single thread
# keeps their floating-point summation order, and so their output, repeatable.
OBABEL_ENV = {"OMP_NUM_THREADS": "1"}

# Open Babel >= 3.2.0 seeds its random generator from this environment
# variable (upstream PR #2241); unset, it draws a random seed. It decides the
# hydrogens whose direction the geometry leaves free (waters, isolated atoms,
# some -h/-p cases) and --gen3d coordinates. Every Open Babel call of the
# pipeline sets it from the run's --seed (run_obabel). Open Babel 3.1.x
# ignores it and seeds from the clock; hydrogens.py makes the hydrogen
# positions deterministic with either version.
OB_RANDOM_SEED_ENV = "OB_RANDOM_SEED"
# First Open Babel release that honours OB_RANDOM_SEED and no longer adds a
# second hydrogen to pre-protonated Lys/Arg (upstream PR #2380).
OBABEL_SEEDABLE_VERSION = (3, 2, 0)
# The release the container image builds (tools/openbabel/repo) and the 2.0
# results were produced with.
OBABEL_VALIDATED_VERSION = (3, 2, 1)

_OBABEL_VERSION = re.compile(r"Open Babel (\d+)\.(\d+)\.(\d+)")

# Lines of Open Babel stderr quoted in an error message.
_STDERR_TAIL_LINES = 5

# Elements Open Babel 3.1 has Gasteiger charge parameters for. Open Babel
# skips a molecule containing any other element when --partialcharge
# gasteiger is requested ("0 molecules converted", exit code 0).
GASTEIGER_ELEMENTS = frozenset({"H", "D", "C", "N", "O", "F", "CL", "BR", "I", "S", "P"})

_CONVERTED = re.compile(r"(\d+) molecules? converted")


def file_sha256(path: str) -> Optional[str]:
    """SHA-256 (hex) of a file's bytes; None when it cannot be read."""
    try:
        digest = hashlib.sha256()
        with open(path, "rb") as f:
            for block in iter(lambda: f.read(1 << 20), b""):
                digest.update(block)
        return digest.hexdigest()
    except OSError:
        return None


def run(cmd_list: List[str], capture: bool = False,
        extra_env: Optional[Dict[str, str]] = None) -> Optional[Tuple[str, str, int]]:
    """Run a command (no shell) with the run's command timeout
    (``tools.settings().command_timeout``).

    ``extra_env`` adds or overrides environment variables for this command only.
    """
    env = {**os.environ, **extra_env} if extra_env else None
    timeout = tools.settings().command_timeout
    try:
        if capture:
            result = subprocess.run(
                cmd_list, capture_output=True, text=True, timeout=timeout,
                env=env)
            if result.returncode != 0:
                raise subprocess.CalledProcessError(
                    result.returncode, cmd_list, output=result.stdout, stderr=result.stderr
                )
            return result.stdout, result.stderr, result.returncode
        subprocess.run(cmd_list, check=True, timeout=timeout, env=env)
    except subprocess.TimeoutExpired:
        logger.error(f"Command timed out after {timeout}s: {shlex.join(cmd_list)}")
        raise
    except subprocess.CalledProcessError as e:
        logger.error(
            f"Command failed with return code {e.returncode}: {shlex.join(cmd_list)}")
        raise
    except FileNotFoundError:
        logger.error(f"Command not found: {shlex.join(cmd_list)}")
        raise
    return None


def _tail(text: Optional[str]) -> str:
    lines = [line for line in (text or "").strip().splitlines() if line.strip()]
    return " | ".join(lines[-_STDERR_TAIL_LINES:])


def log_tool_stderr(tool: str, stderr: str) -> None:
    """Log a tool's stderr at WARNING, collapsing repeated messages.

    Open Babel prints one block per atom for some PDB warnings (e.g. a blank
    element column), separated by '=====' lines; identical blocks (ignoring
    the atom number) are logged once with a count.
    """
    if not stderr or not stderr.strip():
        return
    blocks = [b.strip() for b in re.split(r"^=+\s*$", stderr, flags=re.M) if b.strip()]
    counts = {}
    for block in blocks:
        key = re.sub(r"\(atom \d+\)", "(atom N)", " ".join(block.split()))
        counts[key] = counts.get(key, 0) + 1
    for message, count in counts.items():
        if re.fullmatch(r"\d+ molecules? converted", message):
            logger.debug(f"{tool}: {message}")
            continue
        suffix = f" [x{count}]" if count > 1 else ""
        logger.warning(f"[!] {tool}: {message}{suffix}")


class ObabelNoOutput(RuntimeError):
    """Open Babel exited 0 but converted nothing or wrote no/empty output."""


def obabel_env(seed: Optional[int] = None) -> Dict[str, str]:
    """Environment overrides of every Open Babel call: one OpenMP thread and
    ``OB_RANDOM_SEED`` = ``seed`` (default: ``defaults.SEED``)."""
    value = defaults.SEED if seed is None else int(seed)
    return {**OBABEL_ENV, OB_RANDOM_SEED_ENV: str(value)}


def parse_obabel_version(text: Optional[str]) -> Optional[Tuple[int, int, int]]:
    """(major, minor, patch) from ``obabel -V`` output, None if absent."""
    match = _OBABEL_VERSION.search(text or "")
    return tuple(int(g) for g in match.groups()) if match else None


@lru_cache(maxsize=None)
def _obabel_version_of(path: str) -> Optional[Tuple[int, int, int]]:
    try:
        proc = subprocess.run([path, "-V"], capture_output=True, text=True, timeout=60)
    except (OSError, subprocess.SubprocessError):
        return None
    return parse_obabel_version((proc.stdout or "") + (proc.stderr or ""))


def obabel_version(path: Optional[str] = None) -> Optional[Tuple[int, int, int]]:
    """Version of the Open Babel executable ``path`` (default: the obabel
    found on PATH, ``tools.settings().obabel``), None when it cannot be run
    or read."""
    path = path or tools.settings().obabel
    return _obabel_version_of(path) if path else None


def run_obabel(args: List[str], outputs: Iterable[str] = (),
               return_stderr: bool = False,
               seed: Optional[int] = None) -> Union[str, Tuple[str, str]]:
    """Run Open Babel single-threaded and return its stdout (with
    ``return_stderr``: ``(stdout, stderr)``).

    Every Open Babel call of the pipeline (ligands, charges, receptor) goes
    through here, so each one gets ``OMP_NUM_THREADS=1``,
    ``OB_RANDOM_SEED=<seed>`` (the run's --seed; ``defaults.SEED`` when
    ``seed`` is None), the global command timeout and a DEBUG log line with
    the exact command. On failure the files
    in ``outputs`` are deleted (a partial output must never be mistaken for a
    prepared ligand or receptor) and a RuntimeError with Open Babel's last
    stderr lines is raised. Open Babel exits 0 when it skips every molecule
    (e.g. an element without Gasteiger parameters: "0 molecules converted"),
    so a run that reports 0 molecules converted or leaves an output of
    ``outputs`` missing or empty is a failure too (ObabelNoOutput, a
    RuntimeError).
    """
    obabel = tools.settings().obabel
    if not obabel:
        raise RuntimeError("Open Babel (obabel) not found in PATH")
    cmd = [obabel, *args]
    env = obabel_env(seed)
    logger.debug("Running: " + " ".join(f"{k}={v}" for k, v in env.items())
                 + f" {shlex.join(cmd)}")
    try:
        result = run(cmd, capture=True, extra_env=env)
    except BaseException as e:
        for path in outputs:
            try:
                os.remove(path)
            except FileNotFoundError:
                pass
        if isinstance(e, subprocess.TimeoutExpired):
            raise RuntimeError(
                f"obabel timed out after {tools.settings().command_timeout}s") from e
        if isinstance(e, subprocess.CalledProcessError):
            raise RuntimeError(
                f"obabel exited with code {e.returncode}: {_tail(e.stderr)}") from e
        raise
    stdout, stderr = (result[0], result[1]) if result else ("", "")
    problem = _conversion_problem(stderr, outputs)
    if problem:
        # Open Babel exits 0 when it skips every molecule (e.g. an element
        # without Gasteiger parameters); never let an empty file through.
        for path in outputs:
            try:
                os.remove(path)
            except FileNotFoundError:
                pass
        detail = _tail(_CONVERTED.sub("", stderr or ""))
        raise ObabelNoOutput(f"obabel {problem}" + (f": {detail}" if detail else ""))
    for path in outputs:
        if path.lower().endswith((".sdf", ".sd", ".mol")):
            _strip_molfile_timestamps(path)
        elif path.lower().endswith(".pdbqt"):
            _strip_pdbqt_header_records(path)
    if return_stderr:
        return stdout or "", stderr or ""
    if stderr and stderr.strip():
        logger.debug(f"obabel: {_tail(stderr)}")
    return stdout or ""


def _conversion_problem(stderr: Optional[str], outputs: Iterable[str]) -> Optional[str]:
    """Why a successful (exit code 0) Open Babel run produced nothing usable:
    it reported converting 0 molecules, or an expected output file is missing
    or empty. None if the run looks fine."""
    match = None
    for match in _CONVERTED.finditer(stderr or ""):
        pass
    if match is not None and int(match.group(1)) == 0:
        return "converted 0 molecules (exit code 0)"
    for path in outputs:
        if not os.path.exists(path):
            return f"wrote no output file {os.path.basename(path)}"
        if os.path.getsize(path) == 0:
            return f"wrote an empty {os.path.basename(path)}"
    return None


def elements_without_gasteiger(elements: Iterable[str]) -> List[str]:
    """Sorted element symbols (as given) that Open Babel cannot assign
    Gasteiger charges to."""
    return sorted({e for e in elements if e and e.upper() not in GASTEIGER_ELEMENTS})


def molecule_elements(path: str, fmt: Optional[str] = None,
                      seed: Optional[int] = None) -> List[str]:
    """Element symbols of the molecule(s) in ``path`` from Open Babel's
    formula (best effort for error messages; [] when it cannot be read)."""
    args = ([f"-i{fmt}"] if fmt else []) + [path, "-osmi", "--append", "formula"]
    try:
        stdout, _ = run_obabel(args, return_stderr=True, seed=seed)
    except (RuntimeError, OSError):
        return []
    elements = set()
    for line in (stdout or "").splitlines():
        if line.split():
            # The formula is the last token (after the SMILES and the title).
            elements.update(re.findall(r"[A-Z][a-z]?", line.split()[-1]))
    return sorted(elements)


# PDB header records the Open Babel >= 3.2 PDBQT writer copies from the input
# molecule (src/formats/pdbqtformat.cpp) right after "REMARK  Name = ...".
# AutoDock Vina rejects every one of them except REMARK ("Unknown or
# inappropriate tag"); COMPND carries the input's path (the title of a PDB
# that Open Babel wrote), and a PDBQT input gets its REMARK lines (torsion
# list included) written twice. Open Babel 3.1 never wrote them.
PDBQT_PASSTHROUGH_RECORDS = frozenset({
    "HEADER", "OBSLTE", "TITLE", "SPLIT", "CAVEAT", "COMPND", "SOURCE", "KEYWDS",
    "EXPDTA", "NUMMDL", "MDLTYP", "AUTHOR", "REVDAT", "SPRSDE", "JRNL", "REMARK",
    "DBREF", "DBREF1", "DBREF2", "SEQADV", "SEQRES", "MODRES", "HET", "HETNAM",
    "HETSYN", "FORMUL", "HELIX", "SHEET", "SSBOND", "LINK", "CISPEP", "SITE",
    "ORIGX1", "ORIGX2", "ORIGX3", "SCALE1", "SCALE2", "SCALE3", "MATRIX1",
    "MATRIX2", "MATRIX3", "MODEL",
})
_PDBQT_TORSIONS_REMARK = re.compile(r"^REMARK  \d+ active torsions:")


def _strip_pdbqt_header_records(path: str) -> None:
    """Remove the PDB header records Open Babel >= 3.2 copies into a PDBQT
    (see PDBQT_PASSTHROUGH_RECORDS).

    After each ``REMARK  Name = `` line the writer puts the copied records,
    then (torsion-tree output) its own block starting with the *last*
    ``REMARK  N active torsions:`` line before the first atom or ROOT
    record; a copy of an input PDBQT's REMARK lines can contain such a line
    too. Everything between the name line and that block (rigid output:
    the first non-header record) is removed, so the file is laid out as
    Open Babel 3.1 wrote it."""
    try:
        with open(path, encoding="utf-8") as fh:
            lines = fh.readlines()
    except (FileNotFoundError, UnicodeDecodeError):
        return
    kept = []
    i = 0
    while i < len(lines):
        line = lines[i]
        kept.append(line)
        i += 1
        if not line.startswith("REMARK  Name = "):
            continue
        j = i
        while j < len(lines) and lines[j][:6].rstrip() in PDBQT_PASSTHROUGH_RECORDS:
            j += 1
        own = [k for k in range(i, j) if _PDBQT_TORSIONS_REMARK.match(lines[k])]
        i = own[-1] if own else j
    if len(kept) != len(lines):
        with open(path, "w", encoding="utf-8") as fh:
            fh.writelines(kept)


# MDL molfile header line 2: program name followed by a MMDDYYHHmm date.
_MOLFILE_STAMP = re.compile(r"^(\s*OpenBabel)\d{10}")


def _strip_molfile_timestamps(path: str) -> None:
    """Zero the date Open Babel writes into line 2 of every molfile record,
    so identical structures written in different minutes give identical
    files (the date carries no chemical information)."""
    try:
        with open(path, encoding="utf-8") as fh:
            lines = fh.read().split("\n")
    except (FileNotFoundError, UnicodeDecodeError):
        return
    changed = False
    header = 1  # line 2 of the first record
    for i, line in enumerate(lines):
        if i == header:
            new = _MOLFILE_STAMP.sub(lambda m: m.group(1) + "0" * 10, line)
            changed |= new != line
            lines[i] = new
        if line.strip() == "$$$$":
            header = i + 2
    if changed:
        with open(path, "w", encoding="utf-8") as fh:
            fh.write("\n".join(lines))


def ligand_name(path: str) -> str:
    """Ligand ID used for output names: the file name without a .pdbqt
    extension (matched case-insensitively, so X.PDBQT -> X)."""
    base = os.path.basename(path)
    return base[:-6] if base.lower().endswith(".pdbqt") else base


def available_cpus() -> int:
    """CPUs this process may actually use.

    ``multiprocessing.cpu_count()`` reports every host core, ignoring the
    process's CPU affinity and a container's cgroup quota (docker --cpus), so
    defaults based on it oversubscribe a limited container.
    """
    try:
        n = len(os.sched_getaffinity(0))
    except (AttributeError, OSError):
        n = mp_cpu_count()
    try:
        with open("/sys/fs/cgroup/cpu.max", encoding="utf-8") as f:  # cgroup v2: "<quota> <period>"
            quota, period = f.read().split()[:2]
        if quota != "max":
            n = min(n, max(1, int(int(quota) // int(period))))
    except (OSError, ValueError):
        pass
    return max(1, n)


def parse_pdb_coords(pdb_path: str) -> Tuple[List[float], List[float], List[float]]:
    """Extract XYZ coordinates from PDB file."""
    xs, ys, zs = [], [], []
    try:
        with open(pdb_path, encoding="utf-8") as f:
            for line in f:
                if line.startswith(("ATOM", "HETATM")):
                    try:
                        xs.append(float(line[30:38]))
                        ys.append(float(line[38:46]))
                        zs.append(float(line[46:54]))
                    except (ValueError, IndexError):
                        logger.warning(
                            f"Skipping malformed PDB line: {line.strip()}")
                        continue
    except IOError as e:
        logger.error(f"Cannot read PDB file {pdb_path}: {e}")
        raise

    if not xs:
        raise ValueError(f"No ATOM records found in {pdb_path}")

    return xs, ys, zs


def _safe_ligand_id(name: str) -> str:
    """Sanitize a compound name into a safe file basename token."""
    token = re.sub(r"[^A-Za-z0-9_.-]", "_", name.strip().lower())
    return token or "compound"


class ChainAmbiguityError(ValueError):
    """ATOM records with a blank chain ID next to records of chain
    ``DEFAULT_CHAIN``: the blank chain is read as that chain, so the two
    cannot be told apart."""


def get_chains(pdb_file: str) -> List[str]:
    """Sorted unique chain IDs of the ATOM records of a PDB file.

    ATOM records with a blank chain ID (column 22) are reported as
    ``DEFAULT_CHAIN`` ("A"), the ID receptor preparation reads them as
    (``pdb_records.record_chain``), so ``--chain A`` and ``--chain all`` select
    them and a blank-chain receptor (DUD-E style) is reported as chain A. A
    file in which blank-chain ATOM records coexist with records of a real
    chain A raises ChainAmbiguityError (the two could not be told apart);
    give the blank-chain records a chain ID. An unreadable file raises
    (OSError) instead of silently defaulting. Chain IDs are case-sensitive
    and returned as written.
    """
    chains = set()
    blank = named_default = 0
    with open(pdb_file, errors="replace", encoding="utf-8") as f:
        for line in f:
            if line.startswith("ATOM") and len(line) > 21:
                chain = line[21]
                if chain.strip():
                    chains.add(chain)
                    named_default += chain == DEFAULT_CHAIN
                else:
                    blank += 1
    if blank and named_default:
        raise ChainAmbiguityError(
            f"{pdb_file}: {blank} ATOM record(s) have a blank chain ID and "
            f"{named_default} are in chain {DEFAULT_CHAIN}; a blank chain is read as "
            f"chain {DEFAULT_CHAIN}, so they cannot be told apart. Give the "
            "blank-chain records a chain ID (PDB column 22) other than "
            f"{DEFAULT_CHAIN}.")
    if blank:
        chains.add(DEFAULT_CHAIN)
    return sorted(chains) if chains else [DEFAULT_CHAIN]


def _parse_chain_selection(chain: str, available: List[str]) -> List[str]:
    """Resolve a --chain value ("A", "A,B", "all") against the available IDs.

    PDB chain IDs are case-sensitive (large complexes use both "A" and "a"),
    so IDs are compared exactly; a case-insensitive match is suggested in the
    error message.
    """
    choice = (chain or "").strip()
    if not choice or choice.lower() in {"all", "*"}:
        return available

    selected = [c.strip() for c in re.split(r"[,\s]+", choice) if c.strip()]
    invalid = [c for c in selected if c not in available]
    if invalid:
        hints = []
        for c in invalid:
            matches = [a for a in available if a.lower() == c.lower()]
            if matches:
                hints.append(f"'{c}' -> '{matches[0]}'")
        message = (f"Invalid chain(s): {', '.join(invalid)}. "
                   f"Available chains: {', '.join(available)}")
        if hints:
            message += (" (chain IDs are case-sensitive; did you mean "
                        + ", ".join(hints) + "?)")
            logger.warning(f"[!] {message}")
        raise ValueError(message)
    return selected

#!/usr/bin/env python3
"""Shared helpers: subprocess execution, Open Babel calls, chain parsing.

Cross-module references go through the ``runner`` namespace at call time so
monkeypatching ``runner.<name>`` (e.g. ``runner.run`` in tests) keeps working.
"""
import logging
import os
import re
import shlex
import subprocess
from multiprocessing import cpu_count as mp_cpu_count
from typing import Dict, Iterable, List, Optional, Tuple, Union

# The runner module owns the tool-path globals (OBABEL, VINA*, SMINA,
# COMMAND_TIMEOUT); during the initial import of runner this binds the
# partially-initialized module, which is fine because attributes are only
# accessed at call time.
import runner  # noqa: E402

logger = logging.getLogger(__name__)

# Chain ID assumed for a PDB whose ATOM records carry no chain identifier
# (blank column 22, common in DUD-E-style receptors).
DEFAULT_CHAIN = "A"

# Open Babel's force fields and conformer search use OpenMP; a single thread
# keeps their floating-point summation order, and so their output, repeatable.
OBABEL_ENV = {"OMP_NUM_THREADS": "1"}

# Lines of Open Babel stderr quoted in an error message.
_STDERR_TAIL_LINES = 5

# Elements Open Babel 3.1 has Gasteiger charge parameters for. Open Babel
# skips a molecule containing any other element when --partialcharge
# gasteiger is requested ("0 molecules converted", exit code 0).
GASTEIGER_ELEMENTS = frozenset({"H", "D", "C", "N", "O", "F", "CL", "BR", "I", "S", "P"})

_CONVERTED = re.compile(r"(\d+) molecules? converted")


def run(cmd_list: List[str], capture: bool = False,
        extra_env: Optional[Dict[str, str]] = None) -> Optional[Tuple[str, str, int]]:
    """Run a command (no shell) with the global ``runner.COMMAND_TIMEOUT``.

    ``extra_env`` adds or overrides environment variables for this command only.
    """
    env = {**os.environ, **extra_env} if extra_env else None
    try:
        if capture:
            result = subprocess.run(
                cmd_list, capture_output=True, text=True, timeout=runner.COMMAND_TIMEOUT,
                env=env)
            if result.returncode != 0:
                raise subprocess.CalledProcessError(
                    result.returncode, cmd_list, output=result.stdout, stderr=result.stderr
                )
            return result.stdout, result.stderr, result.returncode
        subprocess.run(cmd_list, check=True, timeout=runner.COMMAND_TIMEOUT, env=env)
    except subprocess.TimeoutExpired:
        logger.error(f"Command timed out after {runner.COMMAND_TIMEOUT}s: {shlex.join(cmd_list)}")
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


class ObabelNoOutput(RuntimeError):
    """Open Babel exited 0 but converted nothing or wrote no/empty output."""


def run_obabel(args: List[str], outputs: Iterable[str] = (),
               return_stderr: bool = False) -> Union[str, Tuple[str, str]]:
    """Run Open Babel single-threaded and return its stdout (with
    ``return_stderr``: ``(stdout, stderr)``).

    Every Open Babel call of the pipeline (ligands, charges, receptor) goes
    through here, so each one gets ``OMP_NUM_THREADS=1``, the global command
    timeout and a DEBUG log line with the exact command. On failure the files
    in ``outputs`` are deleted (a partial output must never be mistaken for a
    prepared ligand or receptor) and a RuntimeError with Open Babel's last
    stderr lines is raised. Open Babel exits 0 when it skips every molecule
    (e.g. an element without Gasteiger parameters: "0 molecules converted"),
    so a run that reports 0 molecules converted or leaves an output of
    ``outputs`` missing or empty is a failure too (ObabelNoOutput, a
    RuntimeError).
    """
    if not runner.OBABEL:
        raise RuntimeError("Open Babel (obabel) not found in PATH")
    cmd = [runner.OBABEL, *args]
    logger.debug(f"Running: OMP_NUM_THREADS=1 {shlex.join(cmd)}")
    try:
        result = runner.run(cmd, capture=True, extra_env=dict(OBABEL_ENV))
    except BaseException as e:
        for path in outputs:
            try:
                os.remove(path)
            except FileNotFoundError:
                pass
        if isinstance(e, subprocess.TimeoutExpired):
            raise RuntimeError(
                f"obabel timed out after {runner.COMMAND_TIMEOUT}s") from e
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


def molecule_elements(path: str, fmt: Optional[str] = None) -> List[str]:
    """Element symbols of the molecule(s) in ``path`` from Open Babel's
    formula (best effort for error messages; [] when it cannot be read)."""
    args = ([f"-i{fmt}"] if fmt else []) + [path, "-osmi", "--append", "formula"]
    try:
        stdout, _ = run_obabel(args, return_stderr=True)
    except (RuntimeError, OSError):
        return []
    elements = set()
    for line in (stdout or "").splitlines():
        if line.split():
            # The formula is the last token (after the SMILES and the title).
            elements.update(re.findall(r"[A-Z][a-z]?", line.split()[-1]))
    return sorted(elements)


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
    (``protein.record_chain``), so ``--chain A`` and ``--chain all`` select
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

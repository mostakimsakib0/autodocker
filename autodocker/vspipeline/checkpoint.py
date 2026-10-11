#!/usr/bin/env python3
"""Checkpoint/resume store for docking results.

Results are appended to a JSON-lines journal as they arrive and consolidated
into ``.docking_checkpoint.json`` by :meth:`DockingCheckpoint.flush`, which
writes atomically (temporary file + ``os.replace``) and then truncates the
journal. An advisory lock (``fcntl``) on ``.docking.lock`` keeps two runs
from sharing one output directory.

The stored timestamp is bookkeeping only: it never reaches ranking.csv or any
other result file.
"""
import fcntl
import json
import logging
import os
from datetime import datetime
from typing import Dict, List, Optional, Tuple

logger = logging.getLogger(__name__)

LOCK_FILENAME = ".docking.lock"

# Lock file descriptors held by this process, keyed by real path of the
# output directory. POSIX record locks belong to the process and are dropped
# when *any* descriptor of the file is closed, so each directory is opened
# exactly once per process and kept open until release_lock().
_HELD_LOCKS: Dict[str, int] = {}


class OutputDirLockedError(RuntimeError):
    """Another process holds the output directory's docking lock."""


def acquire_lock(outdir: str) -> None:
    """Take the advisory lock of ``outdir`` (re-entrant within a process).

    Raises OutputDirLockedError when another process holds it.
    """
    key = os.path.realpath(outdir)
    if key in _HELD_LOCKS:
        return
    os.makedirs(outdir, exist_ok=True)
    path = os.path.join(outdir, LOCK_FILENAME)
    fd = os.open(path, os.O_RDWR | os.O_CREAT, 0o644)
    try:
        fcntl.lockf(fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except OSError:
        os.close(fd)
        raise OutputDirLockedError(
            f"Output directory {outdir} is in use by another AutoDocker run "
            f"(lock held on {path}). Wait for it to finish or use a different "
            "--output directory.")
    os.ftruncate(fd, 0)
    os.write(fd, f"{os.getpid()}\n".encode())
    _HELD_LOCKS[key] = fd


def release_lock(outdir: str) -> None:
    """Release the lock taken by acquire_lock (no-op when not held)."""
    fd = _HELD_LOCKS.pop(os.path.realpath(outdir), None)
    if fd is not None:
        try:
            fcntl.lockf(fd, fcntl.LOCK_UN)
        finally:
            os.close(fd)


def _atomic_write_json(path: str, data) -> None:
    """Write JSON to ``path`` so readers see either the old or the new file
    (insertion order is kept)."""
    tmp = f"{path}.tmp.{os.getpid()}"
    try:
        with open(tmp, "w", encoding="utf-8") as f:
            json.dump(data, f, indent=2)
            f.flush()
            os.fsync(f.fileno())
        os.replace(tmp, path)
    finally:
        if os.path.exists(tmp):
            os.remove(tmp)


class DockingCheckpoint:
    """Manages docking progress and resume capability."""

    def __init__(self, outdir: str, fingerprint: Optional[str] = None,
                 resume: bool = True):
        self.outdir = outdir
        self.fingerprint = fingerprint
        self.checkpoint_file = os.path.join(outdir, ".docking_checkpoint.json")
        self.journal_file = os.path.join(outdir, ".docking_checkpoint.jsonl")
        self.meta_file = os.path.join(outdir, ".docking_checkpoint.meta.json")
        acquire_lock(outdir)
        if not resume:
            self.completed = {}
            self.clear()
        else:
            self.completed = self._load_checkpoint()
        if fingerprint is not None:
            self._check_fingerprint(fingerprint)

    def _check_fingerprint(self, fingerprint: str):
        """Discard cached results produced with a different receptor, box,
        docking parameters or engine (the fingerprint covers all of them)."""
        stored = None
        if os.path.exists(self.meta_file):
            try:
                with open(self.meta_file, encoding="utf-8") as f:
                    stored = json.load(f).get("fingerprint")
            except Exception:
                stored = None
        if self.completed and stored != fingerprint:
            reason = ("was written by an older version without a run fingerprint"
                      if stored is None else
                      "was made with a different receptor, box, docking parameters or engine")
            logger.warning(
                f"[!] Checkpoint {reason}; discarding {len(self.completed)} cached "
                "result(s) and docking every ligand again")
            self.clear()
        try:
            _atomic_write_json(self.meta_file, {"fingerprint": fingerprint})
        except Exception as e:
            logger.warning(f"Could not write checkpoint fingerprint: {e}")

    def clear(self):
        """Forget all cached results, in memory and on disk."""
        self.completed = {}
        for path in (self.checkpoint_file, self.journal_file):
            if os.path.exists(path):
                os.remove(path)

    def _load_checkpoint(self) -> Dict:
        """Load previous docking progress.

        Reads the consolidated ``.json`` (if any) and replays the append-only
        ``.jsonl`` journal so results saved during the current/last run are
        recovered without a full-file rewrite per ligand.
        """
        data: Dict = {}
        if os.path.exists(self.checkpoint_file):
            try:
                with open(self.checkpoint_file, "r", encoding="utf-8") as f:
                    data = json.load(f)
            except Exception as e:
                logger.warning(
                    f"Could not load checkpoint JSON: {e}. Starting fresh.")
                data = {}
        if os.path.exists(self.journal_file):
            try:
                with open(self.journal_file, "r", encoding="utf-8") as f:
                    for line in f:
                        line = line.strip()
                        if not line:
                            continue
                        try:
                            rec = json.loads(line)
                            data[rec["name"]] = {
                                "score": rec["score"],
                                "timestamp": rec.get("timestamp"),
                                "metrics": rec.get("metrics", {}),
                                "ligand_hash": rec.get("ligand_hash"),
                            }
                        except Exception:
                            continue
            except Exception as e:
                logger.warning(f"Could not replay checkpoint journal: {e}")
        if data:
            logger.info(
                f"[*] Loaded checkpoint: {len(data)} ligands already docked")
        return data

    def is_completed(self, ligand_name: str, ligand_hash: Optional[str] = None) -> bool:
        """Check if ligand already docked (from the same file, when a hash is given)."""
        entry = self.completed.get(ligand_name)
        if entry is None:
            return False
        return ligand_hash is None or entry.get("ligand_hash") == ligand_hash

    def save_result(self, ligand_name: str, score: Optional[float],
                    metrics: Optional[Dict] = None,
                    ligand_hash: Optional[str] = None):
        """Save one docking result (``score`` is None for an OUT_OF_WINDOW
        result, whose affinity is kept in ``metrics["raw_score"]``).

        Updates the in-memory record and appends one line to the journal;
        :meth:`flush` consolidates the journal into the JSON file at run end.
        """
        self.completed[ligand_name] = {
            "score": score,
            "timestamp": datetime.now().isoformat(),
            "metrics": metrics or {},
            "ligand_hash": ligand_hash,
        }
        try:
            with open(self.journal_file, "a", encoding="utf-8") as f:
                f.write(json.dumps({
                    "name": ligand_name,
                    "score": score,
                    "timestamp": self.completed[ligand_name]["timestamp"],
                    "metrics": metrics or {},
                    "ligand_hash": ligand_hash,
                }) + "\n")
                f.flush()
        except Exception as e:
            logger.warning(f"Failed to journal checkpoint: {e}")
            self._persist()

        # Keep the consolidated .json present immediately (so external checks
        # see the file right away) but only write it on the first save;
        # later saves rely on the journal and the final flush().
        if not os.path.exists(self.checkpoint_file):
            self._persist()

    def _persist(self) -> bool:
        """Atomically write the in-memory state to the JSON checkpoint."""
        try:
            _atomic_write_json(self.checkpoint_file, self.completed)
            return True
        except Exception as e:
            logger.warning(f"Failed to save checkpoint: {e}")
            return False

    def flush(self):
        """Write the consolidated JSON checkpoint and, once it is safely on
        disk, truncate the journal it now contains (call at run end)."""
        if self._persist() and os.path.exists(self.journal_file):
            try:
                with open(self.journal_file, "w", encoding="utf-8"):
                    pass
            except Exception as e:
                logger.warning(f"Could not truncate checkpoint journal: {e}")

    def close(self):
        """Flush and release the output-directory lock."""
        self.flush()
        release_lock(self.outdir)

    def get_results(self, names: Optional[List[str]] = None) -> List[Tuple[str, float]]:
        """Return completed results as (ligand_name, score) tuples, limited to
        ``names`` (the current library) when given."""
        keep = None if names is None else set(names)
        return [(name, data["score"]) for name, data in self.completed.items()
                if keep is None or name in keep]

    def get_metrics(self, names: Optional[List[str]] = None) -> Dict:
        """Return cached metrics keyed by ligand name (limited to ``names``)."""
        keep = None if names is None else set(names)
        return {name: data.get("metrics", {}) for name, data in self.completed.items()
                if keep is None or name in keep}

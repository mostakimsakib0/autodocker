#!/usr/bin/env python3
"""CHECKPOINT/RESUME SYSTEM.

Split out of runner.py (behavior-preserving refactor). All cross-module
references go through the ``runner`` namespace at call time so the public
entry point and its monkeypatch surface stay unchanged.
"""
import os
import re
import csv
import json
import math
import time
import shlex
import uuid
import shutil
import fnmatch
import glob
import html
import statistics
import subprocess
import logging
from collections import defaultdict
from datetime import datetime
from pathlib import Path
from typing import List, Tuple, Optional, Dict
from urllib.parse import quote
from multiprocessing import Pool, cpu_count as mp_cpu_count

# The runner module owns the tool-path globals (OBABEL, VINA*, SMINA,
# COMMAND_TIMEOUT) and every public name; import it lazily-safe: during the
# initial import of runner this binds the partially-initialized module
# object, which is fine because attributes are only accessed at call time.
import runner  # noqa: E402

logger = logging.getLogger(__name__)

class DockingCheckpoint:
    """Manages docking progress and resume capability."""

    def __init__(self, outdir: str):
        self.outdir = outdir
        self.checkpoint_file = os.path.join(outdir, ".docking_checkpoint.json")
        self.journal_file = os.path.join(outdir, ".docking_checkpoint.jsonl")
        self.completed = self._load_checkpoint()

    def _load_checkpoint(self) -> Dict:
        """Load previous docking progress.

        Reads the consolidated ``.json`` (if any) and replays the append-only
        ``.jsonl`` journal so results saved during the current/last run are
        recovered without a full-file rewrite per ligand.
        """
        data: Dict = {}
        if os.path.exists(self.checkpoint_file):
            try:
                with open(self.checkpoint_file, "r") as f:
                    data = json.load(f)
            except Exception as e:
                logger.warning(
                    f"Could not load checkpoint JSON: {e}. Starting fresh.")
                data = {}
        if os.path.exists(self.journal_file):
            try:
                with open(self.journal_file, "r") as f:
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
                            }
                        except Exception:
                            continue
            except Exception as e:
                logger.warning(f"Could not replay checkpoint journal: {e}")
        if data:
            logger.info(
                f"[✔] Loaded checkpoint: {len(data)} ligands already docked")
        return data

    def is_completed(self, ligand_name: str) -> bool:
        """Check if ligand already docked."""
        return ligand_name in self.completed

    def save_result(self, ligand_name: str, score: float, metrics: Optional[Dict] = None):
        """Save docking result to checkpoint.

        Updates the in-memory record and appends one O(1) line to the journal
        (instead of rewriting the whole JSON each call), eliminating the
        O(N**2) write amplification for large libraries. Call :meth:`flush`
        at the end of a run to consolidate the journal into the JSON file.
        """
        self.completed[ligand_name] = {
            "score": score,
            "timestamp": datetime.now().isoformat(),
            "metrics": metrics or {},
        }
        try:
            with open(self.journal_file, "a") as f:
                f.write(json.dumps({
                    "name": ligand_name,
                    "score": score,
                    "timestamp": self.completed[ligand_name]["timestamp"],
                    "metrics": metrics or {},
                }) + "\n")
        except Exception as e:
            logger.warning(f"Failed to journal checkpoint: {e}")
            self._persist()

        # Keep the consolidated .json present immediately (so external checks
        # see the file right away) but only rewrite it on the first save;
        # subsequent saves rely on the O(1) journal append and the final
        # flush(), avoiding a full-dict rewrite per ligand (was O(N**2) I/O).
        if not os.path.exists(self.checkpoint_file):
            try:
                self._persist()
            except Exception as e:
                logger.warning(f"Failed to seed checkpoint file: {e}")

    def _persist(self):
        """Consolidate the in-memory state into the JSON checkpoint file."""
        try:
            with open(self.checkpoint_file, "w") as f:
                json.dump(self.completed, f, indent=2)
        except Exception as e:
            logger.warning(f"Failed to save checkpoint: {e}")

    def flush(self):
        """Write the consolidated JSON checkpoint (call once at run end)."""
        self._persist()

    def get_results(self) -> List[Tuple[str, float]]:
        """Return all completed results as (ligand_name, score) tuples."""
        return [(name, data["score"]) for name, data in self.completed.items()]

    def get_metrics(self) -> Dict:
        """Return cached metrics keyed by ligand name."""
        return {name: data.get("metrics", {}) for name, data in self.completed.items()}

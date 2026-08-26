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
        self.checkpoint_file = os.path.join(outdir, ".docking_checkpoint.json")
        self.completed = self._load_checkpoint()

    def _load_checkpoint(self) -> Dict:
        """Load previous docking progress."""
        if os.path.exists(self.checkpoint_file):
            try:
                with open(self.checkpoint_file, "r") as f:
                    data = json.load(f)
                logger.info(
                    f"[✔] Loaded checkpoint: {len(data)} ligands already docked")
                return data
            except Exception as e:
                logger.warning(
                    f"Could not load checkpoint: {e}. Starting fresh.")
                return {}
        return {}

    def is_completed(self, ligand_name: str) -> bool:
        """Check if ligand already docked."""
        return ligand_name in self.completed

    def save_result(self, ligand_name: str, score: float, metrics: Optional[Dict] = None):
        """Save docking result to checkpoint."""
        self.completed[ligand_name] = {
            "score": score,
            "timestamp": datetime.now().isoformat(),
            "metrics": metrics or {},
        }
        self._persist()

    def _persist(self):
        """Write checkpoint to disk."""
        try:
            with open(self.checkpoint_file, "w") as f:
                json.dump(self.completed, f, indent=2)
        except Exception as e:
            logger.warning(f"Failed to save checkpoint: {e}")

    def get_results(self) -> List[Tuple[str, float]]:
        """Return all completed results as (ligand_name, score) tuples."""
        return [(name, data["score"]) for name, data in self.completed.items()]

    def get_metrics(self) -> Dict:
        """Return cached metrics keyed by ligand name."""
        return {name: data.get("metrics", {}) for name, data in self.completed.items()}

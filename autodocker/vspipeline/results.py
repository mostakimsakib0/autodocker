#!/usr/bin/env python3
"""RESULTS ANALYSIS (ADVANCED).

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

def _score_sort_key(result: Tuple[str, Optional[float]]) -> Tuple[bool, float]:
    score = result[1]
    return score is None, score if score is not None else float("inf")


def _score_csv(score: Optional[float]):
    return "FAILED" if score is None else score


def _score_text(score: Optional[float]) -> str:
    return "FAILED" if score is None else f"{score:.2f}"

# =============================
# RESULTS ANALYSIS (ADVANCED)
# =============================


class ResultsAnalyzer:
    """Advanced results analysis with SimScore metrics."""

    def __init__(self, outdir: str, top_n: int = 20):
        self.outdir = outdir
        self.top_n = top_n
        self.results_file = os.path.join(outdir, "Results_full.txt")
        self.top_hits_file = os.path.join(outdir, "Top_hits.txt")
        self.csv_file = os.path.join(outdir, "ranking.csv")
        self.metrics_file = os.path.join(outdir, "metrics.txt")
        self.metadata_file = os.path.join(outdir, "ligand_metadata.json")
        self.metadata = self._load_metadata()

    def _load_metadata(self) -> Dict:
        if not os.path.exists(self.metadata_file):
            return {}
        try:
            with open(self.metadata_file) as f:
                return json.load(f)
        except Exception as e:
            logger.warning(f"Could not load ligand metadata: {e}")
            return {}

    def save_ranking(self, results: List[Tuple[str, Optional[float]]], mode: str = "simple",
                     metrics_dict: Optional[Dict] = None):
        """Save docking results with metrics."""
        results_sorted = sorted(results, key=_score_sort_key)

        try:
            with open(self.csv_file, "w", newline="") as f:
                writer = csv.writer(f)
                if metrics_dict:
                    writer.writerow([
                        "Ligand", "Binding_Affinity", "SimScore", "Binding_Modes",
                        "Mode_Coverage", "MW", "LogP", "HBA", "HBD", "TPSA",
                        "Rotors", "ZINC_URL", "Docked_File", "Status"
                    ])
                    for lig, score in results_sorted:
                        metrics = metrics_dict.get(lig, {})
                        meta = self.metadata.get(lig, {})
                        status = metrics.get(
                            "status", "FAILED" if score is None else "OK")
                        writer.writerow([lig, _score_csv(score), metrics.get('simscore', 0.0),
                                         metrics.get('binding_modes', 1),
                                         metrics.get('mode_coverage', ""),
                                         meta.get("mw", ""), meta.get(
                                             "logp", ""),
                                         meta.get("hba", ""), meta.get(
                                             "hbd", ""),
                                         meta.get("tpsa", ""), meta.get(
                                             "rotors", ""),
                                         meta.get("source_url", ""),
                                         metrics.get("dock_file", os.path.join(
                                             self.outdir, "docked", f"{lig}_out.pdbqt")),
                                         status])
                else:
                    writer.writerow(["Ligand", "Binding_Affinity", "Status"])
                    for lig, score in results_sorted:
                        writer.writerow(
                            [lig, _score_csv(score), "FAILED" if score is None else "OK"])
            logger.info(f"[✔] Ranking saved: {self.csv_file}")
        except IOError as e:
            logger.error(f"Failed to write ranking file: {e}")
            raise

        self._save_full_results(results_sorted)
        self._save_top_hits(results_sorted[:self.top_n], metrics_dict)

    def _save_full_results(self, results: List[Tuple[str, Optional[float]]]):
        """Save full results with header."""
        try:
            with open(self.results_file, "w") as f:
                f.write("=" * 70 + "\n")
                f.write("VIRTUAL SCREENING RESULTS\n")
                f.write("=" * 70 + "\n")
                f.write(f"Date: {datetime.now().isoformat()}\n")
                f.write(f"Total compounds screened: {len(results)}\n")
                f.write(f"Top-N report size: {self.top_n}\n")
                f.write("=" * 70 + "\n\n")

                f.write("RANKING (by binding affinity):\n")
                f.write("-" * 70 + "\n")
                f.write(
                    f"{'Rank':<6} {'Ligand':<35} {'Affinity (kcal/mol)':<15}\n")
                f.write("-" * 70 + "\n")

                for rank, (name, score) in enumerate(results, 1):
                    f.write(f"{rank:<6} {name:<35} {_score_text(score):>10}\n")

                f.write("=" * 70 + "\n")
            logger.info(f"[✔] Full results saved: {self.results_file}")
        except IOError as e:
            logger.error(f"Failed to save results: {e}")
            raise

    def _save_top_hits(self, top_results: List[Tuple[str, Optional[float]]],
                       metrics_dict: Optional[Dict] = None):
        """Save top 10 hits with metrics."""
        try:
            with open(self.top_hits_file, "w") as f:
                f.write("=" * 70 + "\n")
                f.write(f"TOP {self.top_n} HITS FOR FURTHER ANALYSIS\n")
                f.write("=" * 70 + "\n")
                f.write(f"Date: {datetime.now().isoformat()}\n")
                f.write("=" * 70 + "\n\n")

                for rank, (name, score) in enumerate(top_results, 1):
                    f.write(f"\n{'─' * 70}\n")
                    f.write(f"{rank}. {name}\n")
                    if score is None:
                        f.write("   Binding Affinity: FAILED\n")
                    else:
                        f.write(f"   Binding Affinity: {score:.2f} kcal/mol\n")

                    if metrics_dict and name in metrics_dict:
                        metrics = metrics_dict[name]
                        meta = self.metadata.get(name, {})
                        f.write(
                            f"   Binding Modes: {metrics.get('binding_modes', 1)}\n")
                        simscore = metrics.get('simscore', 0.0)
                        f.write(
                            f"   SimScore (pose consistency): {simscore:.2f}\n")
                        if metrics.get("mode_coverage") != "":
                            f.write(
                                f"   Mode Coverage: {metrics.get('mode_coverage', 0):.2f}\n")
                        if meta:
                            f.write(
                                f"   MW: {meta.get('mw', 0):.2f}; LogP: {meta.get('logp', 0):.2f}\n")
                        if meta.get("source_url"):
                            f.write(f"   Source: {meta['source_url']}\n")
                        f.write(
                            f"   Docked file: {metrics.get('dock_file', os.path.join(self.outdir, 'docked', name + '_out.pdbqt'))}\n")

                    if score is None:
                        interp = "Docking failed - check ligand/receptor logs"
                    elif score < -7:
                        interp = "Strong binder - HIGH PRIORITY"
                    elif score < -5:
                        interp = "Moderate binder - MEDIUM PRIORITY"
                    else:
                        interp = "Weak binder - LOW PRIORITY"
                    f.write(f"   Interpretation: {interp}\n")
                    if score is None:
                        f.write(
                            f"   Next steps: inspect docked/{name}.log and rerun if needed\n")
                    else:
                        f.write(f"   Next steps: MD/ADMET screening\n")

                f.write("\n" + "=" * 70 + "\n")
            logger.info(f"[✔] Top hits saved: {self.top_hits_file}")
        except IOError as e:
            logger.error(f"Failed to save top hits: {e}")
            raise

    def save_metrics_report(self, results: List[Tuple[str, Optional[float]]], metrics_dict: Dict):
        """Save detailed metrics report."""
        try:
            with open(self.metrics_file, "w") as f:
                f.write("=" * 70 + "\n")
                f.write("ADVANCED METRICS ANALYSIS\n")
                f.write("=" * 70 + "\n")
                f.write(f"Generated: {datetime.now().isoformat()}\n\n")

                f.write("SUMMARY STATISTICS\n")
                f.write("-" * 70 + "\n")

                affinities = [score for _,
                              score in results if score is not None]
                failed = sum(1 for _, score in results if score is None)
                if affinities:
                    f.write(f"Total compounds: {len(results)}\n")
                    f.write(f"Successful dockings: {len(affinities)}\n")
                    f.write(f"Failed dockings: {failed}\n")
                    f.write(f"Best affinity: {min(affinities):.2f} kcal/mol\n")
                    f.write(
                        f"Worst affinity: {max(affinities):.2f} kcal/mol\n")
                    f.write(
                        f"Mean affinity: {sum(affinities)/len(affinities):.2f} kcal/mol\n\n")
                else:
                    f.write(f"Total compounds: {len(results)}\n")
                    f.write(f"Successful dockings: 0\n")
                    f.write(f"Failed dockings: {failed}\n\n")

                f.write("DETAILED METRICS\n")
                f.write("-" * 70 + "\n")
                f.write(
                    f"{'Ligand':<30} {'Affinity':<12} {'SimScore':<12} {'Modes':<8} {'MW':<10} {'LogP':<8}\n")
                f.write("-" * 70 + "\n")

                for name, score in sorted(results, key=_score_sort_key):
                    metrics = metrics_dict.get(name, {})
                    meta = self.metadata.get(name, {})
                    simscore = metrics.get('simscore', 0.0)
                    modes = metrics.get('binding_modes', 1)
                    mw = meta.get("mw", "")
                    logp = meta.get("logp", "")
                    mw_text = f"{mw:.2f}" if isinstance(
                        mw, (int, float)) and mw else ""
                    logp_text = f"{logp:.2f}" if isinstance(
                        logp, (int, float)) or isinstance(logp, float) else ""
                    f.write(
                        f"{name:<30} {_score_text(score):>10}  {simscore:>10.2f}  {modes:>6}  {mw_text:>8}  {logp_text:>6}\n")

                f.write("\n" + "=" * 70 + "\n")
                f.write("INTERPRETATION GUIDE:\n")
                f.write("  Affinity: Lower (more negative) = stronger binding\n")
                f.write(
                    "  SimScore: protocol-style RMSD convergence; higher = more consistent poses\n")
                f.write("  Mode coverage: generated modes / requested modes\n")
                f.write("  Modes: Number of distinct binding poses found\n")
                f.write("=" * 70 + "\n")

            logger.info(f"[✔] Metrics report saved: {self.metrics_file}")
        except IOError as e:
            logger.error(f"Failed to save metrics: {e}")
            raise

#!/usr/bin/env python3
"""Result files: ranking.csv, Results_full.txt, Top_hits.txt, metrics.txt.

All files are a pure function of the docking results (no timestamps, paths
relative to the output directory), so identical inputs give byte-identical
files. Run-level information (date, versions, settings) belongs in the run
manifest written by runner.py.

Ranking order: in-window results by affinity (best first), then
OUT_OF_WINDOW results by affinity, then FAILED ligands; the ligand name
breaks ties so the order does not depend on worker completion order.
"""
import csv
import json
import logging
import os
from typing import Dict, List, Optional, Tuple

from vspipeline import defaults  # noqa: E402
from vspipeline import outputs  # noqa: E402
from vspipeline.docking import (SIMSCORE_RMSD_LB_CUTOFF, SIMSCORE_RMSD_UB_CUTOFF,
                                 Status)

logger = logging.getLogger(__name__)

# Bins used for the plain-language note in Top_hits.txt (kcal/mol).
STRONG_AFFINITY_THRESHOLD = -7.0
MODERATE_AFFINITY_THRESHOLD = -5.0

_GROUP = {Status.OK: 0, Status.WEAK: 0, Status.OUT_OF_WINDOW: 1, Status.FAILED: 2}


def _score_sort_key(result: Tuple[str, Optional[float]]) -> Tuple[bool, float, str]:
    """Sort key for (name, score) pairs: scored first by affinity, then
    unscored; the name breaks ties."""
    score = result[1]
    return score is None, score if score is not None else float("inf"), str(result[0])


def _score_csv(score: Optional[float]):
    return Status.FAILED.value if score is None else score


def _score_text(score: Optional[float]) -> str:
    return Status.FAILED.value if score is None else f"{score:.2f}"


def _blank_if_none(value):
    return "" if value is None else value


def _minimized_text(meta: Dict) -> str:
    """ranking.csv Minimized column: yes/no from the ligand metadata, blank
    when unknown (no metadata for the ligand)."""
    value = meta.get("minimized")
    if value is None:
        return ""
    return "yes" if value else "no"


class ResultsAnalyzer:
    """Writes the ranking and the text reports of one screen."""

    def __init__(self, outdir: str, top_n: int = defaults.TOP_N):
        self.outdir = outdir
        self.top_n = top_n
        self.results_file = os.path.join(outdir, outputs.RESULTS_FULL_TXT)
        self.top_hits_file = os.path.join(outdir, outputs.TOP_HITS_TXT)
        self.csv_file = os.path.join(outdir, outputs.RANKING_CSV)
        self.metrics_file = os.path.join(outdir, outputs.METRICS_TXT)
        self.metadata_file = os.path.join(outdir, outputs.LIGAND_METADATA_JSON)
        self.metadata = self._load_metadata()

    def _load_metadata(self) -> Dict:
        if not os.path.exists(self.metadata_file):
            return {}
        try:
            with open(self.metadata_file, encoding="utf-8") as f:
                return json.load(f)
        except Exception as e:
            logger.warning(f"Could not load ligand metadata: {e}")
            return {}

    # -- helpers -----------------------------------------------------------

    @staticmethod
    def _classify(score: Optional[float], metrics: Optional[Dict]
                  ) -> Tuple[Status, Optional[float]]:
        """(status, affinity to report) of one ligand. OUT_OF_WINDOW rows
        report their raw engine score."""
        metrics = metrics or {}
        status = Status.of(metrics.get("status"), score)
        if status is Status.OUT_OF_WINDOW:
            raw = metrics.get("raw_score")
            return status, float(raw) if isinstance(raw, (int, float)) else None
        if score is None:
            return Status.FAILED, None
        if status is Status.FAILED:  # inconsistent record: trust the score
            status = Status.OK
        return status, score

    def _rows(self, results: List[Tuple[str, Optional[float]]],
              metrics_dict: Optional[Dict]) -> List[Tuple[str, Status, Optional[float]]]:
        """[(name, status, affinity)] in ranking order."""
        metrics_dict = metrics_dict or {}
        rows = []
        for name, score in results:
            status, affinity = self._classify(score, metrics_dict.get(name))
            rows.append((name, status, affinity))
        rows.sort(key=lambda r: (_GROUP[r[1]],
                                 r[2] if r[2] is not None else float("inf"),
                                 str(r[0])))
        return rows

    def _relative_path(self, path: str) -> str:
        """``path`` relative to the output directory (stable across runs)."""
        if not path:
            return ""
        return os.path.relpath(os.path.abspath(path), os.path.abspath(self.outdir))

    def _dock_file(self, name: str, metrics: Dict) -> str:
        return self._relative_path(metrics.get(
            "dock_file", os.path.join(self.outdir, outputs.DOCKED_DIR, name + outputs.POSE_SUFFIX)))

    @staticmethod
    def _affinity_csv(status: Status, affinity: Optional[float]):
        return Status.FAILED.value if status is Status.FAILED or affinity is None else affinity

    @staticmethod
    def _affinity_text(status: Status, affinity: Optional[float]) -> str:
        return Status.FAILED.value if status is Status.FAILED or affinity is None else f"{affinity:.2f}"

    # -- writers -----------------------------------------------------------

    def save_ranking(self, results: List[Tuple[str, Optional[float]]], mode: str = "simple",
                     metrics_dict: Optional[Dict] = None):
        """Write ranking.csv, Results_full.txt and Top_hits.txt.

        ``mode`` is accepted for backwards compatibility; the extended column
        set is written whenever ``metrics_dict`` is given.
        """
        rows = self._rows(results, metrics_dict)

        try:
            with open(self.csv_file, "w", newline="", encoding="utf-8") as f:
                writer = csv.writer(f)
                if metrics_dict:
                    writer.writerow(outputs.RANKING_COLUMNS)
                    for lig, status, affinity in rows:
                        metrics = metrics_dict.get(lig, {})
                        meta = self.metadata.get(lig, {})
                        writer.writerow([
                            lig, self._affinity_csv(status, affinity),
                            _blank_if_none(metrics.get("simscore")),
                            _blank_if_none(metrics.get("binding_modes")),
                            _blank_if_none(metrics.get("mode_coverage")),
                            meta.get("mw", ""), meta.get("logp", ""),
                            meta.get("hba", ""), meta.get("hbd", ""),
                            meta.get("tpsa", ""), meta.get("rotors", ""),
                            meta.get("source_url", ""),
                            self._dock_file(lig, metrics),
                            status.value,
                            metrics.get("engine") or metrics.get("vina_tool", ""),
                            _minimized_text(meta),
                        ])
                else:
                    writer.writerow([outputs.COL_LIGAND, outputs.COL_AFFINITY, outputs.COL_STATUS])
                    for lig, status, affinity in rows:
                        writer.writerow([lig, self._affinity_csv(status, affinity),
                                         status.value])
            logger.info(f"[*] Ranking saved: {self.csv_file}")
        except IOError as e:
            logger.error(f"Failed to write ranking file: {e}")
            raise

        self._save_full_results(rows)
        self._save_top_hits(rows[:self.top_n], metrics_dict)

    def _save_full_results(self, rows: List[Tuple[str, Status, Optional[float]]]):
        """Write Results_full.txt (every ligand, in ranking order)."""
        try:
            with open(self.results_file, "w", encoding="utf-8") as f:
                f.write("=" * 70 + "\n")
                f.write("VIRTUAL SCREENING RESULTS\n")
                f.write("=" * 70 + "\n")
                f.write(f"Total compounds screened: {len(rows)}\n")
                f.write(f"Top-N report size: {self.top_n}\n")
                f.write("=" * 70 + "\n\n")

                f.write("RANKING (by binding affinity):\n")
                f.write("-" * 70 + "\n")
                f.write(
                    f"{'Rank':<6} {'Ligand':<35} {'Affinity (kcal/mol)':<20} Status\n")
                f.write("-" * 70 + "\n")

                for rank, (name, status, affinity) in enumerate(rows, 1):
                    f.write(f"{rank:<6} {name:<35} "
                            f"{self._affinity_text(status, affinity):>10}{'':10} "
                            f"{status.value}\n")

                f.write("=" * 70 + "\n")
            logger.info(f"[*] Full results saved: {self.results_file}")
        except IOError as e:
            logger.error(f"Failed to save results: {e}")
            raise

    @staticmethod
    def _interpretation(status: Status, affinity: Optional[float]) -> str:
        if status is Status.FAILED or affinity is None:
            return "docking failed - see the per-ligand log"
        if status is Status.OUT_OF_WINDOW:
            return "score outside the accepted affinity window - not ranked with the hits"
        if affinity < STRONG_AFFINITY_THRESHOLD:
            return f"strong predicted affinity (< {STRONG_AFFINITY_THRESHOLD:g} kcal/mol)"
        if affinity < MODERATE_AFFINITY_THRESHOLD:
            return (f"moderate predicted affinity ({STRONG_AFFINITY_THRESHOLD:g} to "
                    f"{MODERATE_AFFINITY_THRESHOLD:g} kcal/mol)")
        return f"weak predicted affinity (>= {MODERATE_AFFINITY_THRESHOLD:g} kcal/mol)"

    def _save_top_hits(self, top_rows: List[Tuple[str, Status, Optional[float]]],
                       metrics_dict: Optional[Dict] = None):
        """Write Top_hits.txt: the first ``top_n`` rows of the ranking."""
        try:
            with open(self.top_hits_file, "w", encoding="utf-8") as f:
                f.write("=" * 70 + "\n")
                f.write(f"TOP {self.top_n} HITS\n")
                f.write("=" * 70 + "\n\n")

                for rank, (name, status, affinity) in enumerate(top_rows, 1):
                    f.write(f"\n{'-' * 70}\n")
                    f.write(f"{rank}. {name}\n")
                    if status is Status.FAILED or affinity is None:
                        f.write("   Binding Affinity: FAILED\n")
                    else:
                        f.write(f"   Binding Affinity: {affinity:.2f} kcal/mol\n")
                    f.write(f"   Status: {status.value}\n")

                    if metrics_dict and name in metrics_dict:
                        metrics = metrics_dict[name]
                        meta = self.metadata.get(name, {})
                        if metrics.get("binding_modes") is not None:
                            f.write(f"   Binding Modes: {metrics['binding_modes']}\n")
                        simscore = metrics.get("simscore")
                        if isinstance(simscore, (int, float)):
                            f.write(f"   SimScore (pose consistency): {simscore:.2f}\n")
                        if isinstance(metrics.get("mode_coverage"), (int, float)):
                            f.write(f"   Mode Coverage: {metrics['mode_coverage']:.2f}\n")
                        if (isinstance(meta.get("mw"), (int, float))
                                and isinstance(meta.get("logp"), (int, float))):
                            f.write(f"   MW: {meta['mw']:.2f}; LogP: {meta['logp']:.2f}\n")
                        if meta.get("source_url"):
                            f.write(f"   Source: {meta['source_url']}\n")
                        if metrics.get("engine"):
                            f.write(f"   Engine: {metrics['engine']}\n")
                        f.write(f"   Docked file: {self._dock_file(name, metrics)}\n")

                    f.write(f"   Interpretation: {self._interpretation(status, affinity)}\n")
                    if status is Status.FAILED:
                        f.write(f"   Log: {outputs.DOCKED_DIR}/{name}{outputs.VINA_LOG_SUFFIX}\n")

                f.write("\n" + "=" * 70 + "\n")
            logger.info(f"[*] Top hits saved: {self.top_hits_file}")
        except IOError as e:
            logger.error(f"Failed to save top hits: {e}")
            raise

    def save_metrics_report(self, results: List[Tuple[str, Optional[float]]], metrics_dict: Dict):
        """Write metrics.txt: counts, affinity summary and per-ligand metrics.

        "Successful" counts in-window results; OUT_OF_WINDOW results are
        counted separately and are not failures. The affinity summary covers
        in-window results only.
        """
        metrics_dict = metrics_dict or {}
        rows = self._rows(results, metrics_dict)
        accepted = [a for _, s, a in rows if s in (Status.OK, Status.WEAK)]
        out_of_window = sum(1 for _, s, _ in rows if s is Status.OUT_OF_WINDOW)
        failed = sum(1 for _, s, _ in rows if s is Status.FAILED)
        try:
            with open(self.metrics_file, "w", encoding="utf-8") as f:
                f.write("=" * 70 + "\n")
                f.write("METRICS SUMMARY\n")
                f.write("=" * 70 + "\n\n")

                f.write("SUMMARY STATISTICS\n")
                f.write("-" * 70 + "\n")
                f.write(f"Total compounds: {len(rows)}\n")
                f.write(f"Successful dockings: {len(accepted)}\n")
                f.write(f"Out-of-window results: {out_of_window}\n")
                f.write(f"Failed dockings: {failed}\n")
                if accepted:
                    f.write(f"Best affinity: {min(accepted):.2f} kcal/mol\n")
                    f.write(f"Worst affinity: {max(accepted):.2f} kcal/mol\n")
                    f.write(
                        f"Mean affinity: {sum(accepted)/len(accepted):.2f} kcal/mol\n")
                f.write("\n")

                f.write("DETAILED METRICS\n")
                f.write("-" * 70 + "\n")
                f.write(
                    f"{'Ligand':<30} {'Affinity':<12} {'SimScore':<12} {'Modes':<8} "
                    f"{'MW':<10} {'LogP':<8} Status\n")
                f.write("-" * 70 + "\n")

                for name, status, affinity in rows:
                    metrics = metrics_dict.get(name, {})
                    meta = self.metadata.get(name, {})
                    simscore = metrics.get("simscore")
                    simscore_text = (f"{simscore:.2f}" if isinstance(simscore, (int, float))
                                     else "-")
                    modes = metrics.get("binding_modes")
                    modes_text = "-" if modes is None else str(modes)
                    mw, logp = meta.get("mw"), meta.get("logp")
                    mw_text = f"{mw:.2f}" if isinstance(mw, (int, float)) and mw else ""
                    logp_text = f"{logp:.2f}" if isinstance(logp, (int, float)) else ""
                    f.write(
                        f"{name:<30} {self._affinity_text(status, affinity):>10}  "
                        f"{simscore_text:>10}  {modes_text:>6}  {mw_text:>8}  "
                        f"{logp_text:>6}  {status.value}\n")

                f.write("\n" + "=" * 70 + "\n")
                f.write("INTERPRETATION GUIDE:\n")
                f.write("  Affinity: Lower (more negative) = stronger predicted binding\n")
                f.write(
                    "  SimScore: mean of the fractions of non-best modes with Vina "
                    f"rmsd_lb < {SIMSCORE_RMSD_LB_CUTOFF:g} A and rmsd_ub < "
                    f"{SIMSCORE_RMSD_UB_CUTOFF:g} A; '-' = undefined (single mode)\n")
                f.write("  Mode coverage: generated modes / requested modes\n")
                f.write("  Modes: Number of distinct binding poses found\n")
                f.write("  OUT_OF_WINDOW: valid engine score outside the accepted window\n")
                f.write("=" * 70 + "\n")

            logger.info(f"[*] Metrics report saved: {self.metrics_file}")
        except IOError as e:
            logger.error(f"Failed to save metrics: {e}")
            raise

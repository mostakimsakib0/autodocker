#!/usr/bin/env python3
"""Verify an AutoDocker output directory, byte for byte and against raw Vina.

Modes
-----
``verify_run.py OUTDIR``
    Consistency: the run is complete and its outputs agree with each other.

    * ``run_manifest.json`` exists, is strict JSON (no NaN/Infinity), has the
      required keys and a finished status (``completed``; ``too_many_failures``
      is reported as a warning; ``failed``/``interrupted`` fail);
    * ``ranking.csv`` has the pipeline's columns in order, one row per
      ligand, rows in the pipeline's order (OK/WEAK by score, then
      OUT_OF_WINDOW, then FAILED; ties by name);
    * every scored row has its pose file and its engine log in ``docked/``;
      its Binding_Affinity equals the pose file's first ``REMARK VINA
      RESULT`` (or SMINA ``minimizedAffinity``) at full precision; its Status
      follows from the score (window and WEAK threshold of the run);
    * the log's last ``# Command:`` writes that pose file and uses the run's
      seed, exhaustiveness and number of modes; its ``# Result:`` line agrees
      with ranking.csv; the pose file has no more models than requested;
    * the manifest's docking counts equal the ranking's;
    * the docking binary named in the logs, when present on this machine,
      has the SHA-256 recorded in the manifest, and the log's ``# Engine:``
      line (AutoDocker >= 2.0.0) records that SHA-256;
    * every file in the manifest's ``outputs`` ({relative path: sha256})
      exists with that SHA-256 (so an edit of a pose file AND ranking.csv
      that keeps them consistent is detected); hashed outputs the manifest
      does not list are warnings;
    * the manifest's ``inputs`` (receptor, ligand files, --box-ligand) that
      exist at their recorded path have the recorded SHA-256 (inputs not
      present on this machine are skipped);
    * with --replicates N > 1: the replicate columns of ranking.csv, each
      replicate pose's first score equals replicate_scores.csv, and each
      replicate's log uses seed --seed + k;
    * pose files in ``docked/`` that no ranking row refers to are warnings.

``verify_run.py OUTDIR --redock [--sample K] [--seed S]``
    The consistency checks, then re-execute the exact engine command logged
    at the top of ``docked/<name>_vina.log`` for K seeded-random ligands
    (all with ``--sample 0``). Paths under the original output directory are
    mapped to OUTDIR (so a copied/moved run can be checked) and ``--out`` is
    redirected to a temporary directory. The new pose file must be
    byte-identical to the stored one, every mode's score must match, and the
    ranking.csv affinity must equal raw mode 1.

``verify_run.py RUN_A RUN_B``
    Compare two output directories: identical / different / one-sided files
    (ignoring run_manifest.json timestamps and the ``outputs`` hashes of
    ``*.log`` files, ``*.log``, the Monte-Carlo
    ``Volume`` column of pocket_summary.csv and the ``.docking*`` checkpoint
    files; the run directory path is normalised), the Spearman rank
    correlation of the two rankings and per-ligand score differences.

``--json`` prints one machine-readable JSON document instead of text.
Exit status: 0 = every check passed (warnings allowed), 1 = a check failed,
2 = usage error.

The script needs only the standard library (plus the Vina binary for
``--redock``); it runs outside the package, like the benchmark scripts.
"""

import argparse
import fnmatch
import json
import math
import os
import shlex
import subprocess
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _benchlib  # noqa: E402

# Names mirrored from autodocker/vspipeline/outputs.py, docking.py and
# results.py (tests/test_verify_run.py checks that they agree).
RUN_MANIFEST_JSON = "run_manifest.json"
RECEPTOR_PDBQT = "receptor.pdbqt"
POCKET_SUMMARY_CSV = "pocket_summary.csv"
SMINA_POSE_SUFFIX = "_smina.pdbqt"
SMINA_LOG_SUFFIX = "_smina.log"
RANKING_COLUMNS = (
    "Ligand", "Binding_Affinity", "SimScore", "Binding_Modes", "Mode_Coverage",
    "MW", "LogP", "HBA", "HBD", "TPSA", "Rotors", "ZINC_URL",
    "Docked_File", "Status", "Engine", "Minimized",
)
RANKING_REPLICATE_COLUMNS = ("n_replicates", "score_mean", "score_sd", "score_min",
                             "score_max")
REPLICATE_SCORES_CSV = "replicate_scores.csv"
REPLICATE_TAG = "_rep"
# Outputs whose SHA-256 the manifest records ("outputs"), as in outputs.py.
HASHED_OUTPUTS = (
    "ranking.csv", "Top_hits.txt", "metrics.txt", "Results_full.txt",
    "results_report.html", "pose_clusters.json", "consensus_ranking.csv",
    "consensus_summary.json", "replicate_scores.csv", "replicate_summary.json",
    "receptor.pdbqt", "protein_clean.pdb", "rigid.pdbqt", "flex.pdbqt",
    "receptor_no_reference.pdb", "reference_ligand.pdb",
    "grid.conf", "grid_box.pml", "frame_offset.txt", "ligand_metadata.json",
)
HASHED_OUTPUT_DIRS = ("docked", "ligands", "ligands_minimized")
STATUSES = ("OK", "WEAK", "OUT_OF_WINDOW", "FAILED")
STATUS_GROUP = {"OK": 0, "WEAK": 0, "OUT_OF_WINDOW": 1, "FAILED": 2}
AFFINITY_FLOOR = -20.0
WEAK_AFFINITY_THRESHOLD = -4.0
DEFAULT_MIN_VALID_AFFINITY = -1.0
MANIFEST_REQUIRED_KEYS = ("autodocker_version", "argv", "args", "seed", "status",
                          "exit_code", "tools", "started_at", "finished_at")
MANIFEST_FINISHED = ("completed",)
MANIFEST_FINISHED_WITH_WARNING = ("too_many_failures",)
# Compare mode: files never compared, and manifest keys that differ run to run.
COMPARE_IGNORE = ("*.log", ".docking*", "*/.docking*")
MANIFEST_VOLATILE_KEYS = ("started_at", "finished_at")
POCKET_VOLATILE_COLUMNS = ("Volume",)

PASS, FAIL, WARN, SKIP = "pass", "fail", "warn", "skip"


class Report:
    """Ordered list of check results."""

    def __init__(self):
        self.checks = []

    def add(self, check, status, detail="", **extra):
        self.checks.append({"check": check, "status": status, "detail": detail, **extra})
        return status

    @property
    def failed(self):
        return [c for c in self.checks if c["status"] == FAIL]

    def counts(self):
        out = {s: 0 for s in (PASS, FAIL, WARN, SKIP)}
        for c in self.checks:
            out[c["status"]] += 1
        return out


# ---------------------------------------------------------------------------
# Manifest
# ---------------------------------------------------------------------------

def _reject_constant(token):
    raise ValueError(f"non-standard JSON token {token}")


def load_strict_json(path):
    """Parse ``path`` as strict JSON (NaN/Infinity are errors)."""
    with open(path, encoding="utf-8") as f:
        return json.load(f, parse_constant=_reject_constant)


def check_manifest(outdir, report):
    """Manifest checks; returns the manifest dict (or None)."""
    path = os.path.join(outdir, RUN_MANIFEST_JSON)
    if not os.path.exists(path):
        report.add("manifest.present", FAIL, f"{RUN_MANIFEST_JSON} missing")
        return None
    try:
        manifest = load_strict_json(path)
    except ValueError as e:
        report.add("manifest.strict_json", FAIL, str(e))
        return None
    report.add("manifest.strict_json", PASS)
    missing = [k for k in MANIFEST_REQUIRED_KEYS if k not in manifest]
    report.add("manifest.keys", FAIL if missing else PASS,
               f"missing: {', '.join(missing)}" if missing else "")
    status = manifest.get("status")
    if status in MANIFEST_FINISHED and manifest.get("exit_code") == 0:
        report.add("manifest.status", PASS, status)
    elif status in MANIFEST_FINISHED_WITH_WARNING:
        report.add("manifest.status", WARN,
                   f"{status} (exit {manifest.get('exit_code')}): the run finished but "
                   "more ligands failed than --max-failure-fraction allows")
    else:
        report.add("manifest.status", FAIL,
                   f"status {status!r}, exit code {manifest.get('exit_code')!r}, "
                   f"error {manifest.get('error')!r}")
    return manifest


def manifest_args(manifest):
    return (manifest or {}).get("args") or {}


def expected_engine_sha(manifest, binary):
    """SHA-256 the manifest records for ``binary`` (by path or resolved
    path), or None."""
    if not manifest or not binary:
        return None
    candidates = {binary, os.path.realpath(binary)}
    for entry in manifest.get("docking_engines") or []:
        if entry.get("path") in candidates and entry.get("sha256"):
            return entry["sha256"]
    for info in (manifest.get("tools") or {}).values():
        if not isinstance(info, dict):
            continue
        if (info.get("path") in candidates or info.get("resolved_path") in candidates) \
                and info.get("sha256"):
            return info["sha256"]
    return None


# ---------------------------------------------------------------------------
# Engine logs
# ---------------------------------------------------------------------------

def parse_log(log_path):
    """The last ``# Command:`` (as a list) and the last ``# Result:`` line
    of a per-ligand engine log, i.e. the attempt that wrote the output."""
    command, result = None, None
    with open(log_path, encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("# Command:"):
                command = shlex.split(line.split(":", 1)[1].strip())
                result = None  # a result belongs to the command before it
            elif line.startswith("# Result:"):
                result = line.split(":", 1)[1].strip()
    return command, result


def parse_engine_line(log_path):
    """``(binary, sha256)`` of the ``# Engine:`` line that follows the last
    ``# Command:`` of a per-ligand log (AutoDocker >= 2.0.0 writes it), or
    ``(None, None)``."""
    engine = (None, None)
    with open(log_path, encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("# Command:"):
                engine = (None, None)
            elif line.startswith("# Engine:"):
                binary, _, sha = line.split(":", 1)[1].strip().rpartition(" sha256=")
                engine = (binary or None, sha or None)
    return engine


def command_option(cmd, option):
    """Value following ``option`` in an argument list (last occurrence)."""
    value = None
    for i, arg in enumerate(cmd or []):
        if arg == option and i + 1 < len(cmd):
            value = cmd[i + 1]
        elif arg.startswith(option + "="):
            value = arg.split("=", 1)[1]
    return value


def parse_result_line(text):
    """``(status, affinity)`` of a ``# Result: OK - Affinity: -6.147 kcal/mol``
    line (affinity None when absent)."""
    if not text:
        return None, None
    status = text.split()[0]
    affinity = None
    if "Affinity:" in text:
        token = text.split("Affinity:", 1)[1].split()[0]
        affinity = _benchlib._float_or_none(token)
    return status, affinity


def log_path_for(outdir, name, docked_file):
    """The engine log of a ranking row (Vina or SMINA naming). The log is
    named after the pose file (``<stem>_out.pdbqt`` -> ``<stem>_vina.log``;
    with --replicates the stem is ``<name>_rep<k>``)."""
    docked = os.path.join(outdir, _benchlib.DOCKED_DIR)
    if docked_file.endswith(SMINA_POSE_SUFFIX):
        return os.path.join(docked, name + SMINA_LOG_SUFFIX)
    base = os.path.basename(docked_file)
    stem = base[:-len(_benchlib.POSE_SUFFIX)] if base.endswith(_benchlib.POSE_SUFFIX) else name
    return os.path.join(docked, stem + _benchlib.VINA_LOG_SUFFIX)


def replicate_index(name, docked_file):
    """k of a --replicates pose file ``<name>_rep<k>_out.pdbqt``, else 0."""
    base = os.path.basename(docked_file or "")
    prefix = name + REPLICATE_TAG
    if base.startswith(prefix) and base.endswith(_benchlib.POSE_SUFFIX):
        k = base[len(prefix):-len(_benchlib.POSE_SUFFIX)]
        if k.isdigit():
            return int(k)
    return 0


# ---------------------------------------------------------------------------
# (a) consistency
# ---------------------------------------------------------------------------

def _inside(path, root):
    path, root = os.path.realpath(path), os.path.realpath(root)
    return path == root or path.startswith(root + os.sep)


def expected_status(score, min_valid):
    if not (AFFINITY_FLOOR <= score < min_valid):
        return "OUT_OF_WINDOW"
    return "OK" if score < WEAK_AFFINITY_THRESHOLD else "WEAK"


def read_header(csv_path):
    import csv
    with open(csv_path, newline="", encoding="utf-8") as f:
        return next(csv.reader(f), [])


def check_row(outdir, row, manifest, report, binaries):
    """Checks of one ranking.csv row. Returns the per-ligand record."""
    name, status = row["name"], row["status"]
    args = manifest_args(manifest)
    rec = {"ligand": name, "status": status, "affinity": row["affinity_str"]}
    tag = f"ligand[{name}]"
    if status not in STATUSES:
        report.add(f"{tag}.status", FAIL, f"unknown status {status!r}")
        return rec
    docked_file = row["docked_file"]
    pose = (_benchlib.resolve_docked_file(os.path.join(outdir, _benchlib.RANKING_CSV),
                                          docked_file) if docked_file else None)
    log = log_path_for(outdir, name, docked_file or "")
    rec["pose"], rec["log"] = pose, log
    if status == "FAILED":
        ok = row["affinity_str"] == "FAILED"
        report.add(f"{tag}.failed_row", PASS if ok else FAIL,
                   "" if ok else f"FAILED row with affinity {row['affinity_str']!r}")
        if not os.path.exists(log):
            report.add(f"{tag}.log", WARN, "no engine log (preparation failure?)")
        return rec

    if not pose or not os.path.exists(pose):
        report.add(f"{tag}.pose", FAIL, f"pose file missing: {docked_file!r}")
        return rec
    if not _inside(pose, outdir):
        report.add(f"{tag}.pose", FAIL, f"pose file outside the output dir: {pose}")
        return rec
    if not os.path.exists(log):
        report.add(f"{tag}.log", FAIL, f"engine log missing: {os.path.relpath(log, outdir)}")
        return rec

    token, value = _benchlib.first_vina_result(pose)
    rec["pose_score"] = token
    ranking_value = _benchlib._float_or_none(row["affinity_str"])
    if value is None or ranking_value is None or value != ranking_value:
        report.add(f"{tag}.score", FAIL,
                   f"ranking.csv {row['affinity_str']!r} != pose file {token!r}")
    else:
        report.add(f"{tag}.score", PASS, token)
        min_valid = args.get("min_valid_affinity", DEFAULT_MIN_VALID_AFFINITY)
        want = expected_status(value, float(min_valid))
        report.add(f"{tag}.status", PASS if want == status else FAIL,
                   "" if want == status else f"status {status} but score {token} "
                   f"implies {want}")

    models = _benchlib.parse_vina_models(pose)
    n_modes = args.get("binding_modes")
    if n_modes and len(models) > int(n_modes):
        report.add(f"{tag}.models", FAIL,
                   f"{len(models)} models > --binding-modes {n_modes}")

    cmd, result = parse_log(log)
    if not cmd:
        report.add(f"{tag}.log_command", FAIL, "no '# Command:' line in the log")
        return rec
    rec["command"] = cmd
    binaries.add(cmd[0])
    problems = []
    out = command_option(cmd, "--out")
    if not out or os.path.basename(out) != os.path.basename(pose):
        problems.append(f"--out {out!r} is not {os.path.basename(pose)}")
    k = replicate_index(name, docked_file) if (args.get("replicates") or 1) > 1 else 0
    for option, key in (("--seed", "seed"), ("--exhaustiveness", "exhaustiveness"),
                        ("--num_modes", "binding_modes")):
        want = manifest.get("seed") if key == "seed" and manifest else args.get(key)
        if key == "seed" and want is not None and k:
            want = int(want) + k  # replicate k is docked with --seed + k
        got = command_option(cmd, option)
        if want is not None and got is not None and str(got) != str(want):
            problems.append(f"{option} {got} != manifest {want}")
    res_status, res_aff = parse_result_line(result)
    if res_status is None:
        problems.append("no '# Result:' line after the last command")
    else:
        if res_status != status:
            problems.append(f"log result {res_status} != ranking status {status}")
        if value is not None and res_aff is not None and round(value, 3) != res_aff:
            problems.append(f"log affinity {res_aff} != pose score {token}")
    report.add(f"{tag}.log_command", FAIL if problems else PASS, "; ".join(problems))
    logged_binary, logged_sha = parse_engine_line(log)
    if logged_sha is not None:
        want = expected_engine_sha(manifest, cmd[0])
        if logged_binary != cmd[0]:
            report.add(f"{tag}.log_engine", FAIL,
                       f"'# Engine:' names {logged_binary}, the command runs {cmd[0]}")
        elif want is not None and logged_sha != want:
            report.add(f"{tag}.log_engine", FAIL,
                       f"logged engine sha256 {logged_sha} != manifest {want}")
        else:
            report.add(f"{tag}.log_engine", PASS, logged_sha)
    return rec


def check_ranking_order(rows, report):
    def key(r):
        aff = r["affinity"] if r["status"] != "FAILED" else None
        if aff is None:
            aff = _benchlib._float_or_none(r["affinity_str"])
        return (STATUS_GROUP.get(r["status"], 3),
                aff if aff is not None else math.inf, r["name"])
    ok = [key(r) for r in rows] == sorted(key(r) for r in rows)
    report.add("ranking.order", PASS if ok else FAIL,
               "" if ok else "rows are not in the pipeline's ranking order")


def check_consistency(outdir, report):
    """Mode (a). Returns ``(manifest, records)``."""
    manifest = check_manifest(outdir, report)
    ranking = os.path.join(outdir, _benchlib.RANKING_CSV)
    if not os.path.exists(ranking):
        report.add("ranking.present", FAIL, f"{_benchlib.RANKING_CSV} missing")
        return manifest, []
    header = tuple(read_header(ranking))
    n_replicates = (manifest_args(manifest).get("replicates") or 1) if manifest else 1
    columns = RANKING_COLUMNS + (RANKING_REPLICATE_COLUMNS if n_replicates > 1 else ())
    report.add("ranking.columns", PASS if header == columns else FAIL,
               "" if header == columns else f"columns {list(header)}")
    try:
        rows = _benchlib.read_ranking(ranking)
    except ValueError as e:
        report.add("ranking.read", FAIL, str(e))
        return manifest, []
    names = [r["name"] for r in rows]
    dups = sorted({n for n in names if names.count(n) > 1})
    report.add("ranking.unique", FAIL if dups else PASS,
               f"duplicate ligands: {dups}" if dups else f"{len(rows)} rows")
    check_ranking_order(rows, report)

    binaries = set()
    records = [check_row(outdir, row, manifest or {}, report, binaries) for row in rows]

    if manifest:
        failures = manifest.get("docking_failures")
        if failures:
            n_failed = sum(1 for r in rows if r["status"] == "FAILED")
            ok = failures.get("docked") == len(rows) and failures.get("failed") == n_failed
            report.add("manifest.counts", PASS if ok else FAIL,
                       "" if ok else f"manifest docking_failures {failures} vs ranking "
                       f"{len(rows)} rows / {n_failed} FAILED")
        else:
            report.add("manifest.counts", WARN, "no docking_failures in the manifest")
        for binary in sorted(binaries):
            check_binary(binary, manifest, report)

    referenced = {os.path.realpath(r["pose"]) for r in records if r.get("pose")}
    if n_replicates > 1:
        referenced |= check_replicate_scores(outdir, report)
    check_outputs(outdir, manifest, report)
    check_inputs(manifest, report)
    docked = os.path.join(outdir, _benchlib.DOCKED_DIR)
    if os.path.isdir(docked):
        orphans = sorted(
            f for f in os.listdir(docked)
            if f.endswith((_benchlib.POSE_SUFFIX, SMINA_POSE_SUFFIX))
            and os.path.realpath(os.path.join(docked, f)) not in referenced)
        if orphans:
            report.add("docked.orphans", WARN,
                       f"{len(orphans)} pose file(s) not in ranking.csv: {orphans[:5]}")
    return manifest, records


def check_replicate_scores(outdir, report):
    """--replicates N > 1: every replicate pose of replicate_scores.csv
    exists and its first score equals the Score column. Returns the
    realpaths of those pose files (they are not orphans)."""
    import csv
    path = os.path.join(outdir, REPLICATE_SCORES_CSV)
    if not os.path.exists(path):
        report.add("replicates.scores", FAIL, f"{REPLICATE_SCORES_CSV} missing")
        return set()
    with open(path, newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    poses, problems = set(), []
    for row in rows:
        rel = row.get("Docked_File") or ""
        if not rel:
            continue
        pose = os.path.join(outdir, rel)
        poses.add(os.path.realpath(pose))
        if not os.path.exists(pose):
            problems.append(f"{rel} missing")
            continue
        token, value = _benchlib.first_vina_result(pose)
        score = _benchlib._float_or_none(row.get("Score") or "")
        if score is not None and value != score:
            problems.append(f"{rel}: replicate_scores.csv {row.get('Score')} != pose {token}")
    report.add("replicates.scores", FAIL if problems else PASS,
               "; ".join(problems[:5]) if problems else f"{len(rows)} replicate rows")
    return poses


def _hashed_files(outdir):
    """Relative paths of the outputs a manifest's ``outputs`` covers."""
    found = [n for n in HASHED_OUTPUTS if os.path.isfile(os.path.join(outdir, n))]
    for top in HASHED_OUTPUT_DIRS:
        for root, dirs, files in os.walk(os.path.join(outdir, top)):
            dirs[:] = [d for d in dirs if not d.startswith(".")]
            found += [os.path.relpath(os.path.join(root, f), outdir).replace(os.sep, "/")
                      for f in files if not f.startswith(".")]
    return sorted(found)


def check_outputs(outdir, manifest, report):
    """The manifest's ``outputs`` hashes (AutoDocker >= 2.0.0) against the
    files in ``outdir``."""
    recorded = (manifest or {}).get("outputs")
    if not isinstance(recorded, dict):
        report.add("outputs.sha256", SKIP, "the manifest records no output hashes")
        return
    missing, changed = [], []
    for rel, want in sorted(recorded.items()):
        path = os.path.join(outdir, rel)
        if not os.path.isfile(path):
            missing.append(rel)
        elif _benchlib.file_sha256(path) != want:
            changed.append(rel)
    problems = ([f"changed: {', '.join(changed[:10])}"] if changed else []) + \
               ([f"missing: {', '.join(missing[:10])}"] if missing else [])
    report.add("outputs.sha256", FAIL if problems else PASS,
               "; ".join(problems) if problems else f"{len(recorded)} files",
               changed=changed, missing=missing)
    extra = [rel for rel in _hashed_files(outdir) if rel not in recorded]
    if extra:
        report.add("outputs.unrecorded", WARN,
                   f"{len(extra)} output file(s) not in the manifest (added after the "
                   f"run?): {extra[:5]}")


def check_inputs(manifest, report):
    """The manifest's ``inputs`` that exist at their recorded path have the
    recorded SHA-256; inputs not on this machine are skipped."""
    inputs = (manifest or {}).get("inputs")
    if not isinstance(inputs, dict):
        report.add("inputs.sha256", SKIP, "the manifest records no input hashes")
        return
    entries = [("receptor", inputs.get("receptor")), ("box_ligand", inputs.get("box_ligand"))]
    entries += [(f"ligand[{name}]", entry)
                for name, entry in sorted((inputs.get("ligands") or {}).items())]
    checked, absent, changed = 0, 0, []
    for label, entry in entries:
        if not entry or not entry.get("path"):
            continue
        if not os.path.isfile(entry["path"]):
            absent += 1
            continue
        checked += 1
        if _benchlib.file_sha256(entry["path"]) != entry.get("sha256"):
            changed.append(f"{label} {entry['path']}")
    if changed:
        report.add("inputs.sha256", FAIL, "changed since the run: " + ", ".join(changed[:5]))
    elif checked:
        report.add("inputs.sha256", PASS,
                   f"{checked} input file(s) unchanged, {absent} not on this machine")
    else:
        report.add("inputs.sha256", SKIP, f"none of the {absent} input files is on this "
                   "machine")


def check_binary(binary, manifest, report, required=False):
    """The logged engine binary has the manifest's SHA-256."""
    want = expected_engine_sha(manifest, binary)
    tag = f"engine[{binary}]"
    if want is None:
        return report.add(f"{tag}.sha256", FAIL if required else WARN,
                          "the manifest records no SHA-256 for this binary")
    if not os.path.exists(binary):
        return report.add(f"{tag}.sha256", FAIL if required else SKIP,
                          "binary not present on this machine")
    got = _benchlib.file_sha256(binary)
    return report.add(f"{tag}.sha256", PASS if got == want else FAIL,
                      "" if got == want else f"sha256 {got} != manifest {want}")


# ---------------------------------------------------------------------------
# (b) re-docking
# ---------------------------------------------------------------------------

def original_outdir(cmd):
    """Output directory of the original run, from the logged ``--out``
    (``<outdir>/docked/<name>_out.pdbqt``)."""
    out = command_option(cmd, "--out")
    return os.path.dirname(os.path.dirname(out)) if out else None


def relocate_command(cmd, outdir, new_out):
    """The logged command with paths under the original output directory
    mapped to ``outdir`` and ``--out`` set to ``new_out``.

    Returns ``(command, missing_inputs)``; ``missing_inputs`` lists input
    files (receptor, ligand, flex, config) that do not exist."""
    old = original_outdir(cmd)
    new = []
    for arg in cmd:
        if old and (arg == old or arg.startswith(old.rstrip(os.sep) + os.sep)):
            arg = os.path.join(outdir, os.path.relpath(arg, old))
        new.append(arg)
    for i, arg in enumerate(new):
        if arg == "--out" and i + 1 < len(new):
            new[i + 1] = new_out
    missing = []
    for option in ("--receptor", "--ligand", "--flex", "--config"):
        value = command_option(new, option)
        if value is not None and not os.path.exists(value):
            missing.append(f"{option} {value}")
    return new, missing


def run_engine(cmd, timeout):
    env = {**os.environ, "OMP_NUM_THREADS": os.environ.get("OMP_NUM_THREADS", "1")}
    proc = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout, env=env)
    if proc.returncode != 0:
        raise RuntimeError(f"engine exit {proc.returncode}: {(proc.stderr or '')[-300:]}")


def redock(outdir, manifest, records, report, sample, seed, vina_bin=None,
           timeout=_benchlib.DOCKING_RUN_TIMEOUT, engine=run_engine):
    """Mode (b). Returns the per-ligand re-docking records."""
    candidates = {r["ligand"]: r for r in records if r.get("command") and r.get("pose")}
    if not candidates:
        report.add("redock.candidates", FAIL, "no docked ligand with a logged command")
        return []
    chosen = _benchlib.select_sample(list(candidates), sample, seed)
    report.add("redock.sample", PASS, f"{len(chosen)} of {len(candidates)} ligands "
               f"(sample {sample}, seed {seed})", ligands=chosen)
    checked_binaries = set()
    results = []
    with tempfile.TemporaryDirectory(prefix="verify_redock_") as tmp:
        for name in chosen:
            rec = candidates[name]
            tag = f"redock[{name}]"
            out = os.path.join(tmp, os.path.basename(rec["pose"]))
            cmd, missing = relocate_command(rec["command"], outdir, out)
            if vina_bin:
                cmd = [vina_bin] + cmd[1:]
            elif cmd[0] not in checked_binaries:
                checked_binaries.add(cmd[0])
                check_binary(cmd[0], manifest, report, required=True)
            res = {"ligand": name, "command": shlex.join(cmd), "identical": None,
                   "modes": None, "modes_matched": None, "ranking_match": None}
            results.append(res)
            if missing:
                report.add(tag, FAIL, "input files not found: " + ", ".join(missing))
                continue
            try:
                engine(cmd, timeout)
            except Exception as e:  # recorded as a failed check
                report.add(tag, FAIL, f"re-run failed: {e}")
                continue
            if not os.path.exists(out):
                report.add(tag, FAIL, "re-run wrote no pose file")
                continue
            with open(out, "rb") as a, open(rec["pose"], "rb") as b:
                res["identical"] = a.read() == b.read()
            new_models = _benchlib.parse_vina_models(out)
            old_models = _benchlib.parse_vina_models(rec["pose"])
            res["modes"] = len(new_models)
            res["modes_matched"] = sum(
                1 for x, y in zip(new_models, old_models)
                if x["affinity"] is not None and x["affinity"] == y["affinity"])
            raw1 = new_models[0]["affinity"] if new_models else None
            ranking = _benchlib._float_or_none(rec["affinity"])
            res["ranking_match"] = raw1 is not None and ranking == raw1
            problems = []
            if not res["identical"]:
                problems.append("pose file differs from the re-run")
            if len(new_models) != len(old_models) or res["modes_matched"] != len(new_models):
                problems.append(f"{res['modes_matched']}/{len(old_models)} mode scores match")
            if not res["ranking_match"]:
                problems.append(f"ranking.csv {rec['affinity']} != raw mode 1 "
                                f"{new_models[0]['affinity_str'] if new_models else None}")
            report.add(tag, FAIL if problems else PASS,
                       "; ".join(problems) if problems else "byte-identical, scores equal")
    if vina_bin:
        report.add("redock.binary", WARN,
                   f"--vina-bin {vina_bin} replaces the logged binary (sha256 "
                   f"{_benchlib.file_sha256(vina_bin)})")
    return results


# ---------------------------------------------------------------------------
# (c) comparing two runs
# ---------------------------------------------------------------------------

def _ignored(rel, patterns):
    return any(fnmatch.fnmatch(rel, p) or fnmatch.fnmatch(os.path.basename(rel), p)
               for p in patterns)


def list_files(root, patterns):
    files = []
    for dirpath, dirnames, filenames in os.walk(root):
        dirnames.sort()
        for fname in sorted(filenames):
            rel = os.path.relpath(os.path.join(dirpath, fname), root)
            if not _ignored(rel, patterns):
                files.append(rel)
    return files


def _normalise_manifest(data, run_dir):
    manifest = json.loads(data.decode("utf-8"), parse_constant=_reject_constant)
    for key in MANIFEST_VOLATILE_KEYS:
        manifest.pop(key, None)
    if isinstance(manifest.get("outputs"), dict):
        # Logs name the run directory (and are not compared either).
        manifest["outputs"] = {rel: sha for rel, sha in manifest["outputs"].items()
                               if not _ignored(rel, COMPARE_IGNORE)}
    text = json.dumps(manifest, sort_keys=True, indent=1)
    return _normalise_path(text.encode(), run_dir)


def _normalise_pocket_summary(data):
    import csv
    import io
    rows = list(csv.reader(io.StringIO(data.decode("utf-8"))))
    if not rows:
        return data
    drop = {i for i, col in enumerate(rows[0]) if col in POCKET_VOLATILE_COLUMNS}
    return "\n".join(",".join(c for i, c in enumerate(r) if i not in drop)
                     for r in rows).encode()


def _normalise_path(data, run_dir):
    for prefix in (os.path.realpath(run_dir), os.path.abspath(run_dir)):
        data = data.replace(prefix.encode() + b"/", b"<RUN>/").replace(prefix.encode(), b"<RUN>")
    return data


def normalised_bytes(path, rel, run_dir):
    with open(path, "rb") as f:
        data = f.read()
    if rel == RUN_MANIFEST_JSON:
        try:
            return _normalise_manifest(data, run_dir)
        except ValueError:
            return data
    if rel == POCKET_SUMMARY_CSV:
        data = _normalise_pocket_summary(data)
    return _normalise_path(data, run_dir)


def spearman(xs, ys):
    """Spearman rank correlation (average ranks for ties); None for < 2
    pairs or a constant input."""
    if len(xs) < 2:
        return None
    rx, ry = _average_ranks(xs), _average_ranks(ys)
    mx, my = sum(rx) / len(rx), sum(ry) / len(ry)
    sxy = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    sxx = sum((a - mx) ** 2 for a in rx)
    syy = sum((b - my) ** 2 for b in ry)
    if sxx == 0 or syy == 0:
        return None
    return sxy / math.sqrt(sxx * syy)


def _average_ranks(values):
    order = sorted(range(len(values)), key=lambda i: values[i])
    ranks = [0.0] * len(values)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and values[order[j + 1]] == values[order[i]]:
            j += 1
        for k in range(i, j + 1):
            ranks[order[k]] = (i + j) / 2.0 + 1.0
        i = j + 1
    return ranks


def compare_runs(run_a, run_b, report, ignore=()):
    patterns = COMPARE_IGNORE + tuple(ignore)
    files_a, files_b = set(list_files(run_a, patterns)), set(list_files(run_b, patterns))
    identical, different = [], []
    for rel in sorted(files_a & files_b):
        pa, pb = os.path.join(run_a, rel), os.path.join(run_b, rel)
        if normalised_bytes(pa, rel, run_a) == normalised_bytes(pb, rel, run_b):
            identical.append(rel)
        else:
            different.append(rel)
    only_a, only_b = sorted(files_a - files_b), sorted(files_b - files_a)
    report.add("compare.files", FAIL if (different or only_a or only_b) else PASS,
               f"{len(identical)} identical, {len(different)} different, "
               f"{len(only_a)} only in A, {len(only_b)} only in B")
    result = {"identical": identical, "different": different,
              "only_in_a": only_a, "only_in_b": only_b,
              "ignored_patterns": list(patterns)}

    ra, rb = (os.path.join(r, _benchlib.RANKING_CSV) for r in (run_a, run_b))
    if not (os.path.exists(ra) and os.path.exists(rb)):
        report.add("compare.ranking", FAIL, "ranking.csv missing in one run")
        return result
    rows_a = {r["name"]: r for r in _benchlib.read_ranking(ra)}
    rows_b = {r["name"]: r for r in _benchlib.read_ranking(rb)}
    common = sorted(set(rows_a) & set(rows_b))
    scored = [n for n in common if rows_a[n]["affinity"] is not None
              and rows_b[n]["affinity"] is not None]
    diffs = []
    for name in common:
        a, b = rows_a[name]["affinity"], rows_b[name]["affinity"]
        diffs.append({"ligand": name, "a": a, "b": b,
                      "diff": None if a is None or b is None else b - a,
                      "status_a": rows_a[name]["status"], "status_b": rows_b[name]["status"]})
    rho = spearman([rows_a[n]["affinity"] for n in scored],
                   [rows_b[n]["affinity"] for n in scored])
    changed = [d for d in diffs if d["diff"] != 0 or d["status_a"] != d["status_b"]]
    max_abs = max((abs(d["diff"]) for d in diffs if d["diff"] is not None), default=None)
    result.update({
        "ranking": {"common": len(common), "scored_in_both": len(scored),
                    "only_in_a": sorted(set(rows_a) - set(rows_b)),
                    "only_in_b": sorted(set(rows_b) - set(rows_a)),
                    "spearman_rho": rho, "max_abs_diff": max_abs,
                    "n_changed": len(changed), "per_ligand": diffs}})
    ok = not changed and set(rows_a) == set(rows_b)
    report.add("compare.ranking", PASS if ok else FAIL,
               f"rho={rho if rho is None else round(rho, 4)}, {len(changed)} of "
               f"{len(common)} ligands changed, max |diff| {max_abs}")
    return result


# ---------------------------------------------------------------------------
# Manifest gaps (what a verifier would need that the manifest lacks; runs
# made before the 2.0.0 manifest recorded output/input hashes)
# ---------------------------------------------------------------------------

def manifest_gaps(manifest, records, outdir):
    """Information missing for a complete independent check of the run."""
    gaps = []
    if manifest is None:
        return ["no run_manifest.json"]
    if "outputs" not in manifest:
        gaps.append("no SHA-256 of the output files (e.g. manifest['outputs'] = "
                    "{relative path: sha256} for ranking.csv, docked/*, receptor.pdbqt, "
                    "ligands*/): a consistent edit of a pose file and ranking.csv is "
                    "undetectable without re-docking")
    if "inputs" not in manifest:
        gaps.append("no SHA-256 of the input receptor and ligand files")
    outside = sorted({ligand for ligand, old in (
        (command_option(r["command"], "--ligand"), original_outdir(r["command"]))
        for r in records if r.get("command"))
        if ligand and old and not ligand.startswith(old.rstrip(os.sep) + os.sep)})
    if outside:
        gaps.append(f"{len(outside)} docked ligand file(s) live outside the output "
                    f"directory (e.g. {outside[0]}): the run is not self-contained "
                    "(runs before 2.0.0 docked PDBQT inputs in place)")
    if not manifest.get("docking_engines"):
        gaps.append("no docking_engines (path, sha256, version) entry")
    return gaps


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _print_text(doc):
    for c in doc["checks"]:
        if c["status"] == PASS and not doc.get("verbose"):
            continue
        detail = f": {c['detail']}" if c["detail"] else ""
        print(f"[{c['status'].upper():4}] {c['check']}{detail}")
    if doc.get("compare"):
        cmp_ = doc["compare"]
        for rel in cmp_["different"]:
            print(f"  differs   : {rel}")
        for rel in cmp_["only_in_a"]:
            print(f"  only in A : {rel}")
        for rel in cmp_["only_in_b"]:
            print(f"  only in B : {rel}")
        rk = cmp_.get("ranking")
        if rk:
            print(f"  ranking: {rk['common']} common ligands, Spearman rho "
                  f"{rk['spearman_rho']}, max |diff| {rk['max_abs_diff']}")
            for d in rk["per_ligand"]:
                if d["diff"] != 0:
                    print(f"    {d['ligand']}: {d['a']} -> {d['b']}")
    for gap in doc.get("manifest_gaps", []):
        print(f"[NOTE] manifest gap: {gap}")
    counts = doc["counts"]
    print(f"\n{doc['result']}: {counts['pass']} passed, {counts['fail']} failed, "
          f"{counts['warn']} warnings, {counts['skip']} skipped")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("runs", nargs="+", metavar="OUTDIR",
                    help="One output directory to verify, or two to compare")
    ap.add_argument("--redock", action="store_true",
                    help="Re-execute the logged engine commands and compare poses/scores")
    ap.add_argument("--sample", type=int, default=5,
                    help="Ligands re-docked with --redock (default 5; 0 = all)")
    ap.add_argument("--seed", type=int, default=_benchlib.BENCH_SEED,
                    help="Seed of the --redock sample (default 42)")
    ap.add_argument("--vina-bin", default=None,
                    help="Use this engine binary for --redock instead of the logged one "
                         "(its SHA-256 is reported; the manifest check is skipped)")
    ap.add_argument("--timeout", type=int, default=_benchlib.DOCKING_RUN_TIMEOUT,
                    help="Seconds per re-docked ligand")
    ap.add_argument("--ignore", action="append", default=[], metavar="GLOB",
                    help="Compare mode: extra file pattern to skip (repeatable)")
    ap.add_argument("--json", action="store_true", help="Machine-readable output")
    ap.add_argument("-v", "--verbose", action="store_true", help="Also list passed checks")
    args = ap.parse_args(argv)

    if len(args.runs) > 2:
        ap.error("give one output directory (verify) or two (compare)")
    for run in args.runs:
        if not os.path.isdir(run):
            ap.error(f"not a directory: {run}")
    if len(args.runs) == 2 and args.redock:
        ap.error("--redock takes a single output directory")

    report = Report()
    doc = {"tool": "verify_run", "runs": [os.path.abspath(r) for r in args.runs]}
    if len(args.runs) == 2:
        doc["mode"] = "compare"
        doc["compare"] = compare_runs(args.runs[0], args.runs[1], report, args.ignore)
    else:
        outdir = os.path.abspath(args.runs[0])
        doc["mode"] = "redock" if args.redock else "verify"
        manifest, records = check_consistency(outdir, report)
        doc["autodocker_version"] = (manifest or {}).get("autodocker_version")
        if args.redock:
            doc["redock"] = redock(outdir, manifest, records, report, args.sample,
                                   args.seed, args.vina_bin, args.timeout)
        doc["manifest_gaps"] = manifest_gaps(manifest, records, outdir)
    doc["checks"] = report.checks
    doc["counts"] = report.counts()
    doc["result"] = "FAIL" if report.failed else "PASS"
    doc["verbose"] = args.verbose
    if args.json:
        doc.pop("verbose")
        json.dump(doc, sys.stdout, indent=2, allow_nan=False, default=str)
        sys.stdout.write("\n")
    else:
        _print_text(doc)
    return 1 if report.failed else 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except BrokenPipeError:
        sys.exit(1)
    except KeyboardInterrupt:
        sys.exit(130)

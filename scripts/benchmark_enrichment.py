#!/usr/bin/env python3
"""J2: DUD-E enrichment test (screening power).

Protocol
--------
  * Targets : curated set of DUD-E targets spanning several protein
              families: protein kinases (akt1, braf, cdk2, egfr), serine
              protease (fa10), HIV-1 aspartic protease (hivpr), nuclear
              receptor (andr), reductases (aldr, hmdh) and the esterase
              superfamily (aces)
  * Ligands : per target a seeded sample of 50 unique active and 500
              unique decoy compounds (1:10; ``--max-actives``,
              ``--max-decoys``) from DUD-E (``actives_final.mol2.gz``,
              ``decoys_final.mol2.gz``). DUD-E lists a compound once per
              protonation state/tautomer under the same ID, so entries are
              grouped by compound ID, compounds are sampled, and each one
              is docked as its first entry in the file (ENTRY_RULE; counts
              of entries and compounds go to summary.json, the entry used
              per compound to ``<target>/ligands.csv``). Converted to PDBQT
              with Open Babel; MOL2 intermediates are written to a scratch
              directory, so the ligand directory handed to the runner holds
              only the final PDBQT files (``A_<CHEMBL id>.pdbqt``,
              ``D_<ZINC id>.pdbqt``). A compound whose conversion fails is
              recorded and replaced by the next one in the seeded order.
  * Receptor: DUD-E ``receptor.pdb``, all chains (``--chain all``).
  * Box     : ``--box-source crystal`` (default): cubic box centred on the
              heavy-atom centroid of DUD-E's ``crystal_ligand.mol2``, edge
              max(22.5, extent + 2*6) A (the same rule as the re-docking
              benchmark), passed to the runner as --grid-center/--grid-size
              in input-PDB coordinates. ``--box-source fpocket`` reproduces
              the 1.1.x protocol (runner's fpocket pocket choice).
  * Docking : autodocker ``runner.py`` (Vina), fixed seed, exhaustiveness 4,
              one binding mode.
  * Ranking : lower score = better. Every sampled ligand stays in N: a
              ligand without a usable score (docking failure, not in
              ranking.csv) is ranked last, tied with the other failures.
  * Metrics : primary: ROC AUC (tie-aware) with its bootstrap CI, and EF
              at 5 % with the maximum achievable EF. Secondary (computed,
              reported with caveats, see SECONDARY_METRIC_CAVEATS): EF at
              1 % (n_top = 6 of 550) and BEDROC (alpha = 80.5, Truchon &
              Bayly 2007; alpha * Ra = 7.3 > 1 at 50/550). Also
              ROC enrichment at 1 % / 5 % FPR, stratified bootstrap CIs
              (``--ci``, default 95 %), one-sided Mann-Whitney U test
              (normal approximation, tie-corrected). Ties are resolved by
              their expected value under a random tie-break, so no metric
              depends on input order.

Definitions (N ligands, n actives, n_top = ceil(x N))::

    EF_x     = (TP_x / n_top) / (n / N)
    EF_x,max = (min(n_top, n) / n_top) / (n / N)
    ROCe_x   = TPR(FPR = x) / x

Actives and decoys are labelled by filename prefix ``A_`` / ``D_``.

Modes
-----
* default: full pipeline (download + prep + dock with ``runner.py``)
* ``--prepare-only``: stop after ligand prep; useful for handing the exact
  same prepared ligands to another tool
* ``--ingest DIR``: skip docking; read an external ranking.csv (columns
  Ligand, Binding_Affinity|Affinity) found under ``DIR/<target>/`` and score
  it against the prepared ligand set
"""

import argparse
import csv
import gzip
import hashlib
import http.client
import json
import math
import os
import random
import shutil
import subprocess
import sys
import time
import urllib.error
import urllib.request

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _benchlib  # noqa: E402

DUD_E_BASE = "https://dude.docking.org/targets/{t}/{f}"
# File names on the DUD-E server (also used for the local copies).
DUD_E_RECEPTOR = "receptor.pdb"
DUD_E_CRYSTAL_LIGAND = "crystal_ligand.mol2"
DUD_E_ACTIVES = "actives_final.mol2.gz"
DUD_E_DECOYS = "decoys_final.mol2.gz"

# Functional diversity across protein families so screening power is not
# claimed on a single receptor class.
DEFAULT_TARGETS = [
    "akt1", "braf", "cdk2", "egfr",  # kinases
    "fa10",                          # serine protease
    "hivpr",                         # aspartic protease
    "andr",                          # nuclear hormone receptor
    "aldr", "hmdh",                  # reductases
    "aces",                          # esterase superfamily
]

TARGET_FAMILY = {
    "akt1": "kinase", "braf": "kinase", "cdk2": "kinase", "egfr": "kinase",
    "fa10": "serine protease", "hivpr": "aspartic protease",
    "andr": "nuclear receptor", "aldr": "reductase", "hmdh": "reductase",
    "aces": "esterase",
}

# Confidence level of the bootstrap intervals.
BOOTSTRAP_CI = 0.95
# Conventional significance marks for the Mann-Whitney p-values.
SIGNIFICANCE_STARS = ((0.001, "***"), (0.01, "**"), (0.05, "*"))
# Mixed into the target's seed so decoys are drawn from a different random
# stream than actives (fixed: changing it changes the sampled decoys).
DECOY_SEED_SALT = 0xDEAD
# Seconds for converting one DUD-E MOL2 file to PDBQT with Open Babel.
MOL2_CONVERSION_TIMEOUT = 300
BEDROC_ALPHA = 80.5
EF_FRACTIONS = (0.01, 0.05)
ROCE_FRACTIONS = (0.01, 0.05)
MIN_STABLE_N_TOP = 10
# Default library per target (protocol "option B"): 50 unique active
# compounds + 500 unique decoys, i.e. DUD-E's 1:50 active:decoy ratio is
# reduced to 1:10 to keep the docking cost manageable.
DEFAULT_MAX_ACTIVES = 50
DEFAULT_MAX_DECOYS = 500
# Metrics reported as the result of the benchmark. The others are computed
# and written as well but are secondary, for the reasons below.
PRIMARY_METRICS = ("auc", "ef5")
SECONDARY_METRIC_CAVEATS = {
    "ef1": (f"EF1% rests on n_top = ceil(0.01 N) ligands (6 of 550 with the "
            f"defaults, < {MIN_STABLE_N_TOP}): one active more or less in the top "
            "changes it by a large step, so it is unstable at this library size"),
    "bedroc": (f"BEDROC (alpha = {BEDROC_ALPHA}) assumes alpha * Ra << 1; with 50 "
               "actives in 550 ligands alpha * Ra = 7.3, so BEDROC saturates and is "
               "not comparable with values from full DUD-E sets (Ra about 0.02)"),
}


# ---------------------------------------------------------------------------
# Download / ligand preparation
# ---------------------------------------------------------------------------

def _gzip_ok(path):
    """True for non-.gz files, and for .gz files that are a complete,
    readable gzip stream. The DUD-E server sometimes returns a truncated
    body that only surfaces as 'end-of-stream marker' errors at read time."""
    if not path.endswith(".gz"):
        return True
    try:
        with gzip.open(path, "rb") as f:
            while f.read(1 << 20):
                pass
        return True
    except (OSError, EOFError, gzip.BadGzipFile):
        return False


# DUD-E downloads are large (decoy files of tens of MB) and the server is slow.
DOWNLOAD_TIMEOUT = 180  # seconds per attempt


def download(url, dest, retries=5, backoff=10.0):
    if os.path.exists(dest) and os.path.getsize(dest) > 0 and _gzip_ok(dest):
        return
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    req = urllib.request.Request(
        url, headers={"User-Agent": "autodocker-bench",
                      # The legacy DUD-E server drops/truncates plain GETs on
                      # large files; an open-ended Range header makes it
                      # stream the complete body (verified empirically).
                      "Range": "bytes=0-"},
    )
    last_err = None
    for attempt in range(1, retries + 1):
        try:
            with urllib.request.urlopen(req, timeout=DOWNLOAD_TIMEOUT) as r, open(dest, "wb") as out:
                shutil.copyfileobj(r, out)
            if not _gzip_ok(dest):
                raise RuntimeError("gzip integrity check failed (truncated body)")
            return
        except (urllib.error.URLError, ConnectionError,
                TimeoutError, http.client.IncompleteRead,
                OSError, RuntimeError) as e:
            last_err = e
            if os.path.exists(dest):
                os.remove(dest)  # never leave a truncated file behind
            print(f"    download attempt {attempt}/{retries} failed for "
                  f"{os.path.basename(dest)}: {e}; retrying in "
                  f"{backoff * attempt:.0f}s", flush=True)
            if attempt == retries:
                break
            time.sleep(backoff * attempt)
    raise RuntimeError(f"download failed after {retries} attempts: {last_err}")


def mol2_to_pdbqt(mol2_path, pdbqt_path):
    proc = subprocess.run(
        ["obabel", "-imol2", mol2_path, "-opdbqt", "-O", pdbqt_path],
        capture_output=True, text=True, timeout=MOL2_CONVERSION_TIMEOUT,
        env=_benchlib.obabel_env())
    if (proc.returncode != 0 or not os.path.exists(pdbqt_path)
            or os.path.getsize(pdbqt_path) == 0):
        if os.path.exists(pdbqt_path):
            os.remove(pdbqt_path)
        raise RuntimeError(f"obabel mol2->pdbqt failed for {mol2_path}: "
                           f"{proc.stderr[-300:]}")


def read_molecules(mol2_gz):
    """Return the list of @<TRIPOS>MOLECULE blocks from a compound file."""
    with gzip.open(mol2_gz, "rt", errors="replace") as f:
        text = f.read()
    blocks = text.split("@<TRIPOS>MOLECULE")
    return [("@<TRIPOS>MOLECULE" + b) for b in blocks[1:] if b.strip()]


def molecule_name(block):
    """Compound ID of a MOL2 block (first token of the name line), or None."""
    lines = block.splitlines()
    for i, line in enumerate(lines):
        if line.strip().startswith("@<TRIPOS>MOLECULE"):
            name = lines[i + 1].split() if i + 1 < len(lines) else []
            return name[0] if name else None
    return None


# DUD-E lists each compound under one ID (CHEMBL... for actives, ZINC... for
# decoys) once per protonation state/tautomer, so a compound can have
# several MOL2 entries. The benchmark samples compounds, not entries, and
# docks one entry per compound: the first one in the file (see
# group_compounds). Recorded in summary.json.
ENTRY_RULE = ("one entry per compound ID: the first of its entries in DUD-E's "
              "*_final.mol2.gz (DUD-E documents no preferred state; the first "
              "entry is a fixed, reproducible choice)")


def group_compounds(blocks):
    """``{compound_id: [entry indices in file order]}`` (insertion-ordered).

    Entries without a name get a unique placeholder ID (``entry<i>``), so
    they are never merged."""
    groups = {}
    for idx, block in enumerate(blocks):
        name = molecule_name(block) or f"entry{idx}"
        groups.setdefault(name, []).append(idx)
    return groups


def _file_stem(prefix, compound):
    safe = "".join(ch if ch.isalnum() or ch in "-." else "_" for ch in compound)
    return f"{prefix}{safe}"


def sample_ligands(mol2_gz, out_dir, prefix, max_n, seed, scratch_dir,
                   skip_existing=True):
    """Seeded sample of up to ``max_n`` unique compounds, converted to PDBQT.

    Entries are grouped by compound ID (group_compounds); the IDs, sorted,
    are shuffled with ``random.Random(seed)`` and taken in that order, so
    the sample depends only on the set of IDs and the seed. Each compound is
    docked as one entry (its first in the file, ENTRY_RULE) in the file
    ``<prefix><compound ID>.pdbqt``. MOL2 intermediates go to
    ``scratch_dir``; ``out_dir`` receives only the PDBQT files (the runner
    docks every ligand file it finds there). A compound whose conversion
    fails is recorded and replaced by the next compound in the seeded order.

    Returns ``(pdbqt_paths, failures, info)``; ``failures`` is a list of
    ``{"ligand", "error"}`` dicts; ``info`` holds the entry/compound counts
    of the file and one record per docked compound.
    """
    os.makedirs(out_dir, exist_ok=True)
    os.makedirs(scratch_dir, exist_ok=True)
    paths, failures, selected = [], [], []
    blocks = read_molecules(mol2_gz)
    groups = group_compounds(blocks)
    stems = {}
    for compound in groups:
        stem = _file_stem(prefix, compound)
        if stem in stems:
            raise RuntimeError(f"compound IDs {stems[stem]!r} and {compound!r} map to "
                               f"the same file name {stem}")
        stems[stem] = compound
    order = sorted(groups)
    random.Random(seed).shuffle(order)
    for compound in order:
        if len(paths) >= max_n:
            break
        entries = groups[compound]
        idx = entries[0]
        stem = _file_stem(prefix, compound)
        pdbqt_path = os.path.join(out_dir, stem + ".pdbqt")
        record = {"ligand": stem, "compound": compound, "entry_index": idx,
                  "n_entries": len(entries)}
        if not (skip_existing and os.path.exists(pdbqt_path) and os.path.getsize(pdbqt_path)):
            mol2_path = os.path.join(scratch_dir, stem + ".mol2")
            try:
                with open(mol2_path, "w", encoding="utf-8") as out:
                    out.write(blocks[idx])
                mol2_to_pdbqt(mol2_path, pdbqt_path)
            except Exception as e:  # recorded, never fatal
                failures.append({"ligand": f"{stem} (entry {idx})",
                                 "error": str(e)[-300:]})
                sys.stderr.write(f"  [!] convert failed [{stem}]: {e}\n")
                continue
        paths.append(pdbqt_path)
        selected.append(record)
    info = {"file_entries": len(blocks), "file_compounds": len(groups),
            "file_compounds_with_several_entries": sum(1 for v in groups.values()
                                                       if len(v) > 1),
            "selected_with_several_entries": sum(1 for r in selected if r["n_entries"] > 1),
            "selected": selected}
    return paths, failures, info


def clean_ligand_dir(lig_dir, keep):
    """Remove every file in ``lig_dir`` that is not one of ``keep``.

    The runner docks every ligand file in the directory, so leftovers from
    an earlier run (MOL2 intermediates written by versions <= 1.1.5, or
    PDBQTs from a larger sample) would silently join the screen.
    """
    keep = {os.path.abspath(p) for p in keep}
    removed = []
    for fname in sorted(os.listdir(lig_dir)):
        path = os.path.abspath(os.path.join(lig_dir, fname))
        if os.path.isfile(path) and path not in keep:
            os.remove(path)
            removed.append(fname)
    if removed:
        print(f"    removed {len(removed)} stale file(s) from {lig_dir} "
              f"(e.g. {removed[0]})", flush=True)
    return removed


# ---------------------------------------------------------------------------
# Ranking input
# ---------------------------------------------------------------------------

def is_active(name):
    return name.startswith("A_")


def collect_scores(ranking_csv, expected=None, out_of_window="score"):
    """Scores for the screened ligand set from a ranking.csv.

    ``expected`` is the list of ligand names that were screened; ligands in
    it but absent from ranking.csv, or present without a usable score, get
    score None (ranked last). When ``expected`` is None the ranking's own
    ligands define the set. Rows whose Status is ``OUT_OF_WINDOW`` keep
    their numeric score with ``out_of_window="score"`` and count as failed
    with ``"fail"``.

    Returns ``(names, scores, info)``; ``info`` holds per-ligand status
    strings and counts of unexpected/duplicate rows.
    """
    rows = _benchlib.read_ranking(ranking_csv)
    by_name, duplicates = {}, 0
    for row in rows:
        if row["name"] in by_name:
            duplicates += 1
            continue
        by_name[row["name"]] = row
    names = sorted(expected) if expected is not None else sorted(by_name)
    scores, status = [], {}
    for name in names:
        row = by_name.get(name)
        if row is None:
            scores.append(None)
            status[name] = "MISSING"
            continue
        value = row["affinity"]
        if row["status"] == "OUT_OF_WINDOW" and out_of_window == "fail":
            value = None
        scores.append(value)
        status[name] = row["status"] or ("OK" if value is not None else "FAILED")
    unexpected = sorted(set(by_name) - set(names))
    return names, scores, {"status": status, "unexpected": unexpected,
                           "duplicate_rows": duplicates}


# ---------------------------------------------------------------------------
# Metrics (lower score = better; None = failed, ranked last)
# ---------------------------------------------------------------------------

def _sort_key(score):
    return math.inf if score is None else score


def tie_groups(scores, labels):
    """Ranked tie groups as ``[(n_actives, n_decoys), ...]``, best first.

    All failed ligands (score None) form the last group."""
    groups = {}
    for s, lab in zip(scores, labels):
        g = groups.setdefault(_sort_key(s), [0, 0])
        g[0 if lab else 1] += 1
    return [tuple(groups[k]) for k in sorted(groups)]


def _counts(groups):
    n_pos = sum(a for a, _ in groups)
    n_neg = sum(d for _, d in groups)
    return n_pos, n_neg


def _auc_groups(groups):
    n_pos, n_neg = _counts(groups)
    if not n_pos or not n_neg:
        return None
    pairs, dec_before = 0.0, 0
    for a, d in groups:
        dec_after = n_neg - dec_before - d
        pairs += a * dec_after + 0.5 * a * d
        dec_before += d
    return pairs / (n_pos * n_neg)


def roc_auc(scores, labels):
    """P(random active scores better than random decoy); ties count 1/2."""
    return _auc_groups(tie_groups(scores, labels))


def _n_top(fraction, n_total):
    # round() guards against 0.05 * 100 = 5.000000000000001 -> 6
    return max(1, math.ceil(round(fraction * n_total, 9)))


def _ef_groups(groups, fraction):
    n_pos, n_neg = _counts(groups)
    n_total = n_pos + n_neg
    if not n_pos or not n_total:
        return None
    n_top = _n_top(fraction, n_total)
    remaining, tp = n_top, 0.0
    for a, d in groups:
        if remaining <= 0:
            break
        size = a + d
        take = min(size, remaining)
        tp += a * take / size  # expected actives under random tie-break
        remaining -= take
    base = n_pos / n_total
    return {"ef": (tp / n_top) / base,
            "ef_max": (min(n_top, n_pos) / n_top) / base,
            "n_top": n_top, "tp": tp}


def enrichment_factor(scores, labels, fraction):
    """EF at ``fraction`` of the ranked list (see module docstring).

    Returns a dict with ``ef``, ``ef_max``, ``n_top`` and ``tp`` (expected
    number of actives in the top n_top; fractional only with ties at the
    cut-off), or None without actives."""
    return _ef_groups(tie_groups(scores, labels), fraction)


def _bedroc_groups(groups, alpha=BEDROC_ALPHA):
    n_pos, n_neg = _counts(groups)
    n_total = n_pos + n_neg
    if not n_pos or not n_neg:
        return None
    q = math.exp(-alpha / n_total)
    s, start = 0.0, 1
    for a, d in groups:
        m = a + d
        if a:
            # mean of exp(-alpha r / N) over ranks start..start+m-1
            s += a * (q ** start) * (1.0 - q ** m) / ((1.0 - q) * m)
        start += m
    ra = n_pos / n_total
    rie = s / (ra * (1.0 - math.exp(-alpha)) / (math.exp(alpha / n_total) - 1.0))
    return (rie * ra * math.sinh(alpha / 2.0)
            / (math.cosh(alpha / 2.0) - math.cosh(alpha / 2.0 - alpha * ra))
            + 1.0 / (1.0 - math.exp(alpha * (1.0 - ra))))


def bedroc(scores, labels, alpha=BEDROC_ALPHA):
    """BEDROC (Truchon & Bayly, J. Chem. Inf. Model. 2007, 47, 488).
    Tied actives contribute the mean of their possible exponential terms."""
    return _bedroc_groups(tie_groups(scores, labels), alpha)


def _roce_groups(groups, fpr):
    n_pos, n_neg = _counts(groups)
    if not n_pos or not n_neg or fpr <= 0:
        return None
    fp = tp = 0.0
    for a, d in groups:
        if d == 0:
            tp += a / n_pos
            continue
        if fp >= fpr - 1e-12:
            break
        nfp, ntp = fp + d / n_neg, tp + a / n_pos
        if nfp > fpr + 1e-12:
            tp += (ntp - tp) * (fpr - fp) / (nfp - fp)  # tie = diagonal ROC step
            break
        fp, tp = nfp, ntp
    return tp / fpr


def roc_enrichment(scores, labels, fpr):
    """ROC enrichment: TPR at the given false-positive rate divided by it
    (Jain & Nicholls 2008); ties are linear segments of the ROC curve."""
    return _roce_groups(tie_groups(scores, labels), fpr)


def _ranked_with_ties(values):
    """Average ranks for a list of values (ties share the mean rank)."""
    order = sorted(range(len(values)), key=lambda i: values[i])
    ranks = [0.0] * len(values)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and values[order[j + 1]] == values[order[i]]:
            j += 1
        avg = (i + 1 + j + 1) / 2.0
        for k in range(i, j + 1):
            ranks[order[k]] = avg
        i = j + 1
    return ranks


def mannwhitney(active_scores, decoy_scores):
    """One- and two-sided Mann-Whitney U, normal approximation with tie
    correction. Lower score = better; None (failed) ranks last.

    H1 (one-sided): actives score better (lower) than decoys.
    Returns ``(U_actives, p_one_sided, p_two_sided)``; p-values are kept at
    full double precision (erfc keeps small tails representable).
    """
    n1, n2 = len(active_scores), len(decoy_scores)
    if n1 == 0 or n2 == 0:
        return None, None, None
    allv = [_sort_key(s) for s in list(active_scores) + list(decoy_scores)]
    ranks = _ranked_with_ties(allv)
    u1 = sum(ranks[:n1]) - n1 * (n1 + 1) / 2.0
    mu = n1 * n2 / 2.0
    n = n1 + n2
    ties = {}
    for v in allv:
        ties[v] = ties.get(v, 0) + 1
    tie_adj = sum(t ** 3 - t for t in ties.values() if t > 1)
    var_u = n1 * n2 / 12.0 * ((n + 1) - tie_adj / (n * (n - 1.0)))
    if var_u <= 0:
        return u1, 1.0, 1.0
    z = (u1 - mu) / math.sqrt(var_u)
    p_one = 0.5 * math.erfc(-z / math.sqrt(2.0))       # P(Z <= z)
    p_two = min(1.0, math.erfc(abs(z) / math.sqrt(2.0)))
    return u1, p_one, p_two


def _percentile(values, p):
    values = sorted(values)
    k = (len(values) - 1) * p
    f, c = math.floor(k), math.ceil(k)
    if f == c:
        return values[int(k)]
    return values[int(f)] * (c - k) + values[int(c)] * (k - f)


def metric_set(scores, labels):
    """All point metrics for one ranked list (dict of floats or None)."""
    groups = tie_groups(scores, labels)
    out = {"auc": _auc_groups(groups),
           "bedroc": _bedroc_groups(groups)}
    for frac in EF_FRACTIONS:
        tag = f"{frac * 100:g}"
        ef = _ef_groups(groups, frac)
        out[f"ef{tag}"] = ef["ef"] if ef else None
        out[f"ef{tag}_max"] = ef["ef_max"] if ef else None
        out[f"ef{tag}_n_top"] = ef["n_top"] if ef else None
    for frac in ROCE_FRACTIONS:
        out[f"roce{frac * 100:g}"] = _roce_groups(groups, frac)
    return out


BOOTSTRAP_METRICS = ("auc", "bedroc", "ef1", "ef5", "roce1", "roce5")


def bootstrap_ci(scores, labels, n_boot, seed, ci=BOOTSTRAP_CI):
    """Stratified percentile-bootstrap confidence intervals.

    Each replicate draws ``n_actives`` (score, label) pairs with
    replacement from the actives and ``n_decoys`` from the decoys, by
    index, so duplicated draws stay duplicated, and recomputes every metric.
    Returns ``{metric: (lo, hi)}`` for the central ``ci`` band, or None.
    """
    if n_boot <= 0 or not 0.0 < ci < 1.0:
        return None
    pos = [s for s, lab in zip(scores, labels) if lab]
    neg = [s for s, lab in zip(scores, labels) if not lab]
    if not pos or not neg:
        return None
    rng = random.Random(seed)
    labels_b = [True] * len(pos) + [False] * len(neg)
    draws = {m: [] for m in BOOTSTRAP_METRICS}
    for _ in range(n_boot):
        sample = ([pos[rng.randrange(len(pos))] for _ in range(len(pos))]
                  + [neg[rng.randrange(len(neg))] for _ in range(len(neg))])
        values = metric_set(sample, labels_b)
        for m in BOOTSTRAP_METRICS:
            if values.get(m) is not None:
                draws[m].append(values[m])
    alpha = (1.0 - ci) / 2.0
    return {m: (_percentile(v, alpha), _percentile(v, 1.0 - alpha))
            for m, v in draws.items() if v}


def _sig(p):
    if p is None:
        return ""
    for cutoff, stars in SIGNIFICANCE_STARS:
        if p < cutoff:
            return stars
    return ""


def compute_stats(names, scores, seed, n_boot, ci=BOOTSTRAP_CI):
    """Point metrics, bootstrap CIs and the Mann-Whitney test for a target.

    ``names`` carry the A_/D_ label; ``scores`` may contain None (failed)."""
    labels = [is_active(n) for n in names]
    point = metric_set(scores, labels)
    act = [s for s, lab in zip(scores, labels) if lab]
    dec = [s for s, lab in zip(scores, labels) if not lab]
    _u, p1, p2 = mannwhitney(act, dec)
    bands = bootstrap_ci(scores, labels, n_boot, seed, ci) or {}
    row = dict(point)
    for m in BOOTSTRAP_METRICS:
        lo_hi = bands.get(m)
        row[f"{m}_lo"] = lo_hi[0] if lo_hi else None
        row[f"{m}_hi"] = lo_hi[1] if lo_hi else None
    warnings = []
    for frac in EF_FRACTIONS:
        n_top = point.get(f"ef{frac * 100:g}_n_top")
        if n_top is not None and n_top < MIN_STABLE_N_TOP:
            warnings.append(
                f"EF{frac * 100:g}% uses only n_top={n_top} ligands "
                f"(< {MIN_STABLE_N_TOP}); it is unstable at this library size")
    if names and act and BEDROC_ALPHA * len(act) / len(names) > 1.0:
        warnings.append(
            f"BEDROC: alpha*Ra = {BEDROC_ALPHA * len(act) / len(names):.1f} > 1 "
            "(active fraction too high for alpha=80.5); BEDROC saturates and "
            "is not comparable with published DUD-E values")
    for frac in ROCE_FRACTIONS:
        if dec and frac * len(dec) < 1:
            warnings.append(f"ROCe{frac * 100:g}% is below one decoy "
                            f"({len(dec)} decoys)")
    row.update({
        "p_mw_one": p1, "p_mw_two": p2, "sig": _sig(p1),
        "total": len(names), "n_actives": len(act), "n_decoys": len(dec),
        "n_failed_actives": sum(1 for s in act if s is None),
        "n_failed_decoys": sum(1 for s in dec if s is None),
        "ci_level": ci, "bootstrap": n_boot, "bedroc_alpha": BEDROC_ALPHA,
        "warnings": warnings,
    })
    return row


# ---------------------------------------------------------------------------
# Per-target pipeline
# ---------------------------------------------------------------------------

def write_metrics(target, outdir, names, scores, status, ranking):
    """Per-ligand ranked table (failed ligands last) plus the raw ranking."""
    os.makedirs(os.path.join(outdir, target), exist_ok=True)
    metrics = os.path.join(outdir, target, "metrics.csv")
    order = sorted(range(len(names)), key=lambda i: (_sort_key(scores[i]), names[i]))
    with open(metrics, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["rank", "name", "label", "affinity", "status"])
        for rank, i in enumerate(order, 1):
            writer.writerow([rank, names[i],
                             "active" if is_active(names[i]) else "decoy",
                             "" if scores[i] is None else repr(scores[i]),
                             status.get(names[i], "")])
    dest = os.path.join(outdir, target, _benchlib.RANKING_CSV)
    if os.path.abspath(ranking) != os.path.abspath(dest):
        shutil.copy(ranking, dest)
    return metrics


def write_ligand_table(target, outdir, selected):
    """``<outdir>/<target>/ligands.csv``: the docked compounds, with the DUD-E
    entry used (0-based position in the .mol2.gz) and its number of entries."""
    os.makedirs(os.path.join(outdir, target), exist_ok=True)
    path = os.path.join(outdir, target, "ligands.csv")
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=["ligand", "compound", "entry_index",
                                               "n_entries"])
        writer.writeheader()
        writer.writerows(selected)
    return path


def crystal_box(crystal_mol2, padding, min_size):
    center, size = _benchlib.ligand_box(
        _benchlib.mol2_heavy_coords(crystal_mol2), padding=padding,
        min_size=min_size)
    return center, size


def process_target(target, workdir, opts):
    rec = {"target": target, "family": TARGET_FAMILY.get(target, "?"),
           "status": "error", "n_actives": 0, "n_decoys": 0,
           "box_source": opts["box_source"], "error": None}
    tdir = os.path.join(workdir, target)
    target_seed = int(hashlib.md5(target.encode()).hexdigest()[:8], 16)
    try:
        receptor = os.path.join(tdir, DUD_E_RECEPTOR)
        actives_gz = os.path.join(tdir, DUD_E_ACTIVES)
        decoys_gz = os.path.join(tdir, DUD_E_DECOYS)
        download(DUD_E_BASE.format(t=target, f=DUD_E_RECEPTOR), receptor)
        download(DUD_E_BASE.format(t=target, f=DUD_E_ACTIVES), actives_gz)
        download(DUD_E_BASE.format(t=target, f=DUD_E_DECOYS), decoys_gz)

        lig_dir = os.path.join(tdir, "ligands")
        scratch = os.path.join(tdir, "ligands_mol2")
        actives, fail_a, info_a = sample_ligands(actives_gz, lig_dir, "A_",
                                                 opts["max_actives"], target_seed, scratch)
        decoys, fail_d, info_d = sample_ligands(decoys_gz, lig_dir, "D_",
                                                opts["max_decoys"],
                                                target_seed ^ DECOY_SEED_SALT, scratch)
        clean_ligand_dir(lig_dir, actives + decoys)
        rec["n_actives"], rec["n_decoys"] = len(actives), len(decoys)
        rec["n_convert_failed_actives"] = len(fail_a)
        rec["n_convert_failed_decoys"] = len(fail_d)
        rec["convert_failures"] = fail_a + fail_d
        for label, info in (("actives", info_a), ("decoys", info_d)):
            rec[f"dude_{label}_entries"] = info["file_entries"]
            rec[f"dude_{label}_compounds"] = info["file_compounds"]
            rec[f"dude_{label}_compounds_with_several_entries"] = (
                info["file_compounds_with_several_entries"])
            rec[f"sampled_{label}_with_several_entries"] = info["selected_with_several_entries"]
        write_ligand_table(target, opts["outdir"], info_a["selected"] + info_d["selected"])
        for label, n, wanted in (("actives", len(actives), opts["max_actives"]),
                                 ("decoys", len(decoys), opts["max_decoys"])):
            if n < wanted:
                rec.setdefault("warnings", []).append(
                    f"only {n} unique {label} usable (requested {wanted})")
        if not actives or not decoys:
            raise RuntimeError("no usable actives/decoys after conversion")
        expected = [os.path.basename(p)[:-len(".pdbqt")] for p in actives + decoys]

        if opts["box_source"] == "crystal":
            crystal = os.path.join(tdir, DUD_E_CRYSTAL_LIGAND)
            download(DUD_E_BASE.format(t=target, f=DUD_E_CRYSTAL_LIGAND), crystal)
            center, size = crystal_box(crystal, opts["padding"], opts["min_box"])
            rec["box_center"] = [round(c, 3) for c in center]
            rec["box_size"] = round(size, 3)

        if opts["prepare_only"]:
            rec["status"] = "prepared"
            return rec

        outdir = os.path.join(tdir, "dock")
        if not opts["ingest"]:
            cmd = [
                sys.executable, opts["runner"],
                "-r", receptor, "-l", lig_dir, "-o", outdir,
                f"--chain={opts['chain']}",
                "--seed", str(opts["seed"]),
                "--exhaustiveness", str(opts["exhaustiveness"]),
                "--binding-modes", "1",
                "--energy-range", "1.0",
                "--processes", str(opts["processes"]),
            ]
            if opts["box_source"] == "crystal":
                cmd += _benchlib.box_cli_args(center, size)
            rec["runner_cmd"] = " ".join(cmd)
            proc = subprocess.run(cmd, capture_output=True, text=True,
                                  timeout=opts["timeout"])
            with open(os.path.join(tdir, "runner.log"), "w", encoding="utf-8") as f:
                f.write(proc.stdout or "")
                f.write(proc.stderr or "")
            ranking = os.path.join(outdir, _benchlib.RANKING_CSV)
            if proc.returncode != 0 or not os.path.exists(ranking):
                raise RuntimeError(
                    f"docking failed (rc={proc.returncode}): "
                    f"{(proc.stdout or '')[-400:]}{(proc.stderr or '')[-400:]}")
        else:
            ranking = os.path.join(opts["ingest"], target, _benchlib.RANKING_CSV)
            if not os.path.exists(ranking):
                raise RuntimeError(f"ingest ranking.csv missing: {ranking}")

        names, scores, info = collect_scores(ranking, expected,
                                             opts["out_of_window"])
        if info["unexpected"]:
            rec["n_unexpected_ranking_rows"] = len(info["unexpected"])
            print(f"    [!] {len(info['unexpected'])} ranking row(s) are not "
                  f"in the sampled set and were ignored "
                  f"(e.g. {info['unexpected'][0]})", flush=True)
        rec["n_out_of_window"] = sum(1 for s in info["status"].values()
                                     if s == "OUT_OF_WINDOW")
        rec["n_missing_from_ranking"] = sum(1 for s in info["status"].values()
                                            if s == "MISSING")
        rec["metrics"] = write_metrics(target, opts["outdir"], names, scores,
                                       info["status"], ranking)
        stats = compute_stats(names, scores, target_seed, opts["bootstrap"], opts["ci"])
        stats["warnings"] = rec.get("warnings", []) + stats["warnings"]
        rec.update(stats)
        rec["status"] = "ok"
    except Exception as e:
        rec["error"] = str(e)
    return rec


SUMMARY_COLS = [
    "target", "family", "status", "box_source", "total", "n_actives",
    "n_decoys", "n_failed_actives", "n_failed_decoys",
    "n_convert_failed_actives", "n_convert_failed_decoys",
    "dude_actives_entries", "dude_actives_compounds",
    "dude_decoys_entries", "dude_decoys_compounds",
    # primary metrics
    "auc", "auc_lo", "auc_hi",
    "ef5", "ef5_max", "ef5_n_top", "ef5_lo", "ef5_hi",
    # secondary metrics (see SECONDARY_METRIC_CAVEATS)
    "ef1", "ef1_max", "ef1_n_top", "ef1_lo", "ef1_hi",
    "bedroc", "bedroc_lo", "bedroc_hi",
    "roce1", "roce1_lo", "roce1_hi", "roce5", "roce5_lo", "roce5_hi",
    "p_mw_one", "p_mw_two", "sig", "ci_level", "bootstrap", "error",
]


def save_tables(results, outdir, protocol):
    rows_path = os.path.join(outdir, "enrichment.csv")
    with open(rows_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=SUMMARY_COLS, extrasaction="ignore")
        writer.writeheader()
        for r in results:
            writer.writerow(r)
    with open(os.path.join(outdir, "summary.json"), "w", encoding="utf-8") as f:
        json.dump({"protocol": protocol,
                   "primary_metrics": list(PRIMARY_METRICS),
                   "secondary_metrics_caveats": SECONDARY_METRIC_CAVEATS,
                   "targets": results}, f, indent=2)
    print(f"\nResults written to {rows_path}")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--targets", default=",".join(DEFAULT_TARGETS),
                    help="Comma-separated DUD-E target codes "
                         f"(default: {','.join(DEFAULT_TARGETS)})")
    ap.add_argument("--max-actives", type=int, default=DEFAULT_MAX_ACTIVES,
                    help="Unique active compounds sampled per target "
                         f"(default {DEFAULT_MAX_ACTIVES})")
    ap.add_argument("--max-decoys", type=int, default=DEFAULT_MAX_DECOYS,
                    help="Unique decoy compounds sampled per target "
                         f"(default {DEFAULT_MAX_DECOYS})")
    ap.add_argument("--bootstrap", type=int, default=1000,
                    help="Bootstrap replicates for CIs (default 1000)")
    ap.add_argument("--ci", type=float, default=BOOTSTRAP_CI,
                    help="Confidence level of the bootstrap intervals (default 0.95)")
    ap.add_argument("--box-source", choices=["crystal", "fpocket"], default="crystal",
                    help="Docking box: DUD-E crystal ligand (default) or the "
                         "runner's fpocket pocket (1.1.x protocol)")
    ap.add_argument("--padding", type=float, default=_benchlib.BENCH_BOX_PADDING,
                    help="Crystal box padding per side, A (default 6.0)")
    ap.add_argument("--min-box", type=float, default=_benchlib.BENCH_MIN_BOX_EDGE,
                    help="Minimum crystal box edge, A (default 22.5)")
    ap.add_argument("--chain", default="all",
                    help="Receptor chain(s) passed to the runner (default all)")
    ap.add_argument("--out-of-window", choices=["score", "fail"], default="score",
                    help="Rows with Status OUT_OF_WINDOW: rank by their score "
                         "(default) or count them as failed")
    ap.add_argument("--prepare-only", action="store_true",
                    help="Download+prep ligands only (no docking)")
    ap.add_argument("--ingest", default=None, metavar="DIR",
                    help="Do not dock; read DIR/<target>/ranking.csv instead")
    ap.add_argument("--workdir", default="/tmp/autodocker-bench/dude")
    ap.add_argument("--outdir", default="benchmarks/results/enrichment")
    ap.add_argument("--runner", default="autodocker/runner.py")
    ap.add_argument("--seed", type=int, default=_benchlib.PIPELINE_DEFAULTS.SEED)
    ap.add_argument("--exhaustiveness", type=int, default=4)
    ap.add_argument("--processes", type=int, default=_benchlib.available_cpus(),
                    help="Ligands docked in parallel (default: available CPUs)")
    ap.add_argument("--timeout", type=int, default=_benchlib.ENRICHMENT_RUN_TIMEOUT)
    args = ap.parse_args(argv)
    if not 0.0 < args.ci < 1.0:
        ap.error("--ci must be between 0 and 1")

    targets = [t.strip() for t in args.targets.split(",") if t.strip()]
    os.makedirs(args.outdir, exist_ok=True)
    runner = os.path.abspath(args.runner)
    protocol = {
        "targets": targets, "box_source": args.box_source,
        "padding": args.padding, "min_box": args.min_box, "chain": args.chain,
        "max_actives": args.max_actives, "max_decoys": args.max_decoys,
        "sampling_unit": "compound (DUD-E compound ID), not MOL2 entry",
        "entry_per_compound": ENTRY_RULE,
        "primary_metrics": list(PRIMARY_METRICS),
        "secondary_metrics_caveats": SECONDARY_METRIC_CAVEATS,
        "seed": args.seed, "exhaustiveness": args.exhaustiveness,
        "binding_modes": 1, "processes": args.processes,
        "bootstrap": args.bootstrap, "ci": args.ci,
        "bedroc_alpha": BEDROC_ALPHA, "out_of_window": args.out_of_window,
        "failed_ligands": "kept in N, ranked last (tied)",
        "ingest": args.ingest,
        "runner_version": (None if args.ingest
                           else _benchlib.runner_version(runner)),
        "vina_version": _benchlib.tool_version("vina"),
    }
    results = []
    for t in targets:
        print(f"\n=== [{t}] ({TARGET_FAMILY.get(t, '?')}) ===", flush=True)
        rec = process_target(t, args.workdir, {
            "max_actives": args.max_actives, "max_decoys": args.max_decoys,
            "seed": args.seed, "exhaustiveness": args.exhaustiveness,
            "processes": args.processes, "timeout": args.timeout,
            "runner": runner, "bootstrap": args.bootstrap, "ci": args.ci,
            "box_source": args.box_source, "padding": args.padding,
            "min_box": args.min_box, "chain": args.chain,
            "out_of_window": args.out_of_window,
            "prepare_only": args.prepare_only,
            "ingest": args.ingest, "outdir": args.outdir})
        results.append(rec)
        print(f"[{t}] {rec['status']} actives={rec['n_actives']} "
              f"decoys={rec['n_decoys']} failed="
              f"{rec.get('n_failed_actives')}/{rec.get('n_failed_decoys')} "
              f"auc={rec.get('auc')} ({rec.get('auc_lo')}-{rec.get('auc_hi')}) "
              f"ef5={rec.get('ef5')} (max {rec.get('ef5_max')}) "
              f"p_mw={rec.get('p_mw_one')} {rec.get('sig')} | secondary: "
              f"ef1={rec.get('ef1')} bedroc={rec.get('bedroc')} "
              f"err={rec['error']}", flush=True)
        for w in rec.get("warnings", []):
            print(f"    [!] {w}", flush=True)
        save_tables(results, args.outdir, protocol)  # partial runs stay usable

    save_tables(results, args.outdir, protocol)
    aucs = [r["auc"] for r in results if r["status"] == "ok" and r.get("auc") is not None]
    if aucs:
        print(f"\nMean AUC over {len(aucs)} targets: {sum(aucs) / len(aucs):.3f}")
    return 0 if all(r["status"] in ("ok", "prepared") for r in results) else 1


if __name__ == "__main__":
    sys.exit(main())

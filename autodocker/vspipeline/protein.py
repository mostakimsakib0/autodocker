#!/usr/bin/env python3
"""Receptor preparation: the ``ProteinPreparation`` class used by the
pipeline, and receptor cleaning.

ProteinPreparation runs the receptor stage on top of these modules:

* pdb_records: PDB record reading (first model, one alternate location per
  residue, element perception),
* protonation: Open Babel protonation, deterministic hydrogens and residue
  label restoration,
* fpocket: fpocket pockets (running, parsing, ranking, selection),
* gridbox: box geometry and the grid/PyMOL/frame-offset writers.

Cleaning (chain selection, hydrogens, HETATM records, selenomethionine) is
defined here.
"""
import logging
import os
import tempfile
from collections import defaultdict
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Set, Tuple

from vspipeline import (  # noqa: E402
    charges,
    defaults,
    fpocket,
    gridbox,
    outputs,
    pocket,
    protonation,
    residue_templates,
    utils,
)
from vspipeline.pdb_records import (  # noqa: E402
    BLANK_CHAIN_ID,
    _has_blank_chain,
    _heavy_atom_coords,
    _residue_position,
    is_hydrogen_record,
    read_pdb_records,
    record_chain,
)

logger = logging.getLogger(__name__)

# Modified amino acids deposited as HETATM; they are part of the protein and
# dropping them would leave holes in the backbone or side chains.
MODIFIED_RESIDUES = {
    "MSE", "SEP", "TPO", "PTR", "CSO", "CSD", "CME", "MLY", "HYP", "KCX",
    "PCA", "CSS", "OCS", "SCH", "NEP", "ALY", "M3L", "TYS",
}

# Selenomethionine is written as methionine (Se -> S): AutoDock Vina has no
# selenium atom type and Open Babel no Gasteiger parameters for Se, so a SeMet
# structure could not be converted at all. The side chains are isosteric
# (C-Se 1.95 A vs C-S 1.81 A); coordinates are kept as deposited.
SELENOMETHIONINE = "MSE"
_SEMET_TO_MET_ATOM = {"SE": " SD "}

# A receptor with fewer ATOM records than this is rejected as truncated.
MIN_PROTEIN_ATOM_LINES = 100
# Default padding (A) added on each side of the pocket atoms' bounding box.
DEFAULT_PADDING = defaults.PADDING
# Minimum edge (A) of the whole-receptor box used with --no-fpocket.
MIN_RECEPTOR_BOX_SIZE = 24.0


# ---------------------------------------------------------------------------
# Cleaning
# ---------------------------------------------------------------------------

def _semet_to_met(line: str) -> Tuple[str, Tuple[str, str]]:
    """A selenomethionine (MSE) record as methionine: residue name MET, the
    selenium as ' SD ' with element S. Returns (line, residue position)."""
    name = line[12:16]
    element = line[76:78] if len(line) > 76 else ""
    if name.strip() in _SEMET_TO_MET_ATOM:
        name = _SEMET_TO_MET_ATOM[name.strip()]
        element = " S"
    line = line[:12] + name + line[16:17] + "MET" + line[20:]
    if len(line) > 76:
        line = line[:76] + element + line[78:]
    return line, _residue_position(line)


@dataclass
class _CleaningCounts:
    """What receptor cleaning kept and dropped (for the log)."""
    kept_atoms: int = 0
    skipped_hetero: int = 0
    removed_h: int = 0
    semet: Set[Tuple[str, str]] = field(default_factory=set)  # residue positions
    # chain -> [ATOM, HETATM] records of the chains not selected
    other_chains: Dict[str, List[int]] = field(
        default_factory=lambda: defaultdict(lambda: [0, 0]))


def _clean_records(lines: Sequence[str], selected_chains: Sequence[str],
                   keep_hetero: bool) -> Tuple[List[str], _CleaningCounts]:
    """The records of ``pdb_clean`` (see ProteinPreparation._write_clean_pdb)
    from the first-model, one-altloc ``lines``, and the counts."""
    counts = _CleaningCounts()
    last_was_atom = False
    out = []
    for line in lines:
        is_atom = line.startswith("ATOM")
        is_het = line.startswith("HETATM")
        if is_het and line[17:20].strip() in MODIFIED_RESIDUES:
            line = "ATOM  " + line[6:]
            is_atom, is_het = True, False
        if is_atom or (keep_hetero and is_het):
            if len(line) <= 21 or record_chain(line) not in selected_chains:
                if len(line) > 21:
                    counts.other_chains[record_chain(line)][0 if is_atom else 1] += 1
                continue
            # Input hydrogens are dropped: Open Babel's pH model would
            # protonate the already protonated groups once more and put
            # the extra hydrogen in a random direction.
            if is_hydrogen_record(line):
                counts.removed_h += 1
                continue
            if is_atom and line[17:20] == SELENOMETHIONINE:
                line, residue = _semet_to_met(line)
                counts.semet.add(residue)
            out.append(line)
            counts.kept_atoms += 1
            last_was_atom = True
        elif is_het:
            counts.skipped_hetero += 1
        elif line.startswith("TER") and last_was_atom:
            out.append(line)
            last_was_atom = False
    out.append("END\n")
    return out, counts


def _log_cleaning(counts: _CleaningCounts, skipped_altloc: int, keep_hetero: bool) -> None:
    """Report what cleaning dropped or changed."""
    if skipped_altloc:
        logger.info(f"[*] Kept one alternate location per residue (highest "
                    f"occupancy); dropped {skipped_altloc} alternate-location "
                    "atom records")
    if counts.removed_h:
        logger.warning(
            f"[!] Removed {counts.removed_h} hydrogen records from the input; the receptor "
            "is re-protonated by Open Babel at the requested pH (input "
            "protonation states are not kept)")
    if counts.semet:
        logger.warning(
            f"[!] Wrote {len(counts.semet)} selenomethionine (MSE) residue(s) as "
            "methionine (MET, Se -> S): AutoDock Vina has no selenium atom type "
            "and Open Babel no Gasteiger parameters for Se; coordinates are kept")
    if counts.skipped_hetero and not keep_hetero:
        logger.info(
            f"[*] Removed {counts.skipped_hetero} HETATM records during receptor cleaning")
    het_only = sorted(c for c, (n_atom, _) in counts.other_chains.items() if n_atom == 0)
    if keep_hetero and het_only:
        logger.warning(
            "[!] --keep-hetero: not keeping HETATM records of chain(s) without "
            "protein atoms, which --chain never selects: " + ", ".join(
                f"{c} ({counts.other_chains[c][1]} records)" for c in het_only))
    dropped_het = sum(n_het for c, (n_atom, n_het) in counts.other_chains.items() if n_atom)
    if keep_hetero and dropped_het:
        logger.info(f"[*] --keep-hetero: {dropped_het} HETATM records of chains "
                    "outside the selection are not kept")


# ---------------------------------------------------------------------------
# Receptor files
# ---------------------------------------------------------------------------

def _translate_pdbqt(path: str, offset: Tuple[float, float, float]) -> None:
    """Subtract ``offset`` from every ATOM/HETATM coordinate of a PDBQT file."""
    dx, dy, dz = offset
    out = []
    with open(path, encoding="utf-8") as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")) and len(line) >= 54:
                x = float(line[30:38]) - dx
                y = float(line[38:46]) - dy
                z = float(line[46:54]) - dz
                line = f"{line[:30]}{x:8.3f}{y:8.3f}{z:8.3f}{line[54:]}"
            out.append(line)
    with open(path, "w", encoding="utf-8") as f:
        f.writelines(out)


def _write_without_residues(src: str, dst: str, residues: set) -> None:
    """Copy a PDB without the HETATM residues whose ``pocket._residue_key``
    is in ``residues``."""
    with open(src, encoding="utf-8") as f:
        lines = f.readlines()
    kept = []
    for line in lines:
        parsed = pocket._parse_hetatm(line)
        if parsed is not None and pocket._residue_key(parsed) in residues:
            continue
        kept.append(line)
    with open(dst, "w", encoding="utf-8") as f:
        f.writelines(kept)


class ProteinPreparation:
    """Receptor preparation: chain selection, cleaning, protonation, pocket
    detection and the docking box."""

    def __init__(self, pdb_file: str, outdir: str, seed: int = defaults.SEED):
        self.pdb_file = pdb_file
        self.outdir = outdir
        # OB_RANDOM_SEED of the receptor conversion (the run's --seed).
        self.seed = seed
        self.pdb_clean = os.path.join(outdir, outputs.PROTEIN_CLEAN_PDB)
        self.receptor_pdbqt = os.path.join(outdir, outputs.RECEPTOR_PDBQT)
        self.grid_conf = os.path.join(outdir, outputs.GRID_CONF)
        self.grid_box_script = os.path.join(outdir, outputs.GRID_BOX_PML)
        self.pocket_summary_file = os.path.join(outdir, outputs.POCKET_SUMMARY_CSV)
        # Trp/His ring-hydrogen decisions (residue_templates.apply_residue_templates);
        # None until prepare_receptor has run.
        self.residue_templates: Optional[Dict] = None

    def validate(self):
        """Validate protein PDB file."""
        if not os.path.exists(self.pdb_file):
            raise FileNotFoundError(f"Protein file not found: {self.pdb_file}")

        if not os.access(self.pdb_file, os.R_OK):
            raise PermissionError(f"Cannot read protein file: {self.pdb_file}")

        try:
            with open(self.pdb_file, encoding="utf-8") as f:
                count = sum(1 for line in f if line.startswith("ATOM"))
        except IOError as e:
            raise IOError(f"Error reading protein file: {e}")

        if count < MIN_PROTEIN_ATOM_LINES:
            raise ValueError(
                f"Invalid protein: only {count} ATOM lines "
                f"(need >= {MIN_PROTEIN_ATOM_LINES})")

        logger.info(f"[*] Protein validated: {count} ATOM lines")

    def select_chain(self) -> str:
        """Choose a chain for docking when --chain is not given.

        Non-interactive: with several protein chains the first one (in
        sorted order) is used and a warning lists the alternatives.
        """
        chains = utils.get_chains(self.pdb_file)

        if len(chains) == 1:
            note = ""
            if chains[0] == BLANK_CHAIN_ID and _has_blank_chain(self.pdb_file):
                note = (f" (blank chain ID in the PDB, read as {BLANK_CHAIN_ID}; "
                        f"--chain {BLANK_CHAIN_ID} and residue keys such as "
                        f"{BLANK_CHAIN_ID}45 or 45 address it)")
            logger.info(f"[*] Single chain detected: {chains[0]}{note}")
            return chains[0]

        logger.warning(
            f"[!] WARNING: {len(chains)} protein chains present "
            f"({', '.join(chains)}); docking against chain '{chains[0]}' only. "
            "Pass --chain (e.g. --chain all or --chain A,B) to choose explicitly.")
        return chains[0]

    def _write_clean_pdb(self, selected_chains: Sequence[str], keep_hetero: bool) -> None:
        """Write ``pdb_clean``: first model, one altloc per residue, the
        selected chains, no hydrogens, HETATM records only with
        ``keep_hetero`` (modified amino acids are always kept as ATOM;
        selenomethionine is written as methionine, see SELENOMETHIONINE)."""
        try:
            lines, skipped_altloc = read_pdb_records(self.pdb_file)
        except IOError as e:
            raise IOError(f"Failed to read {self.pdb_file}: {e}")

        out, counts = _clean_records(lines, selected_chains, keep_hetero)
        if counts.kept_atoms == 0:
            raise ValueError(
                f"Selected chain(s) '{', '.join(selected_chains)}' have no atoms "
                f"in {self.pdb_file}")
        try:
            with open(self.pdb_clean, "w", encoding="utf-8") as dst:
                dst.writelines(out)
        except IOError as e:
            raise IOError(
                f"Failed to filter chain(s) {', '.join(selected_chains)}: {e}")
        _log_cleaning(counts, skipped_altloc, keep_hetero)

    def _metal_ions(self) -> List[Tuple[str, str, str]]:
        """Residue keys (see ``pocket._residue_key``) of the single-atom metal
        ions in ``pdb_clean`` that AutoDock Vina can type."""
        ions = []
        for key, atoms in pocket._hetatm_residues(self.pdb_clean).items():
            if len(atoms) == 1 and atoms[0]["element"] in pocket.VINA_METAL_TYPE_MAP:
                ions.append(key)
        return ions

    def _protonate(self, ph: float, source: Optional[str] = None) -> None:
        """Convert ``source`` (default ``pdb_clean``) to ``receptor_pdbqt``
        with Open Babel, adding hydrogens for ``ph``, with deterministic
        hydrogen positions (protonation.protonate)."""
        protonation.protonate(source or self.pdb_clean, self.receptor_pdbqt, ph, self.seed)

    def prepare_receptor(self, chain: str = BLANK_CHAIN_ID, keep_hetero: bool = False,
                         ph: float = defaults.PH) -> List[str]:
        """Clean protein and convert to PDBQT."""
        available_chains = utils.get_chains(self.pdb_file)
        selected_chains = utils._parse_chain_selection(chain, available_chains)
        logger.info(
            f"[*] Preparing receptor (chain(s): {', '.join(selected_chains)})...")

        if os.path.realpath(self.pdb_file) == os.path.realpath(self.pdb_clean):
            raise ValueError(
                f"Receptor input {self.pdb_file} is the file the cleaned receptor "
                "is written to; move it out of the output directory")

        self._write_clean_pdb(selected_chains, keep_hetero)
        self.selected_chains = selected_chains

        try:
            # Metal ions kept by --keep-hetero stay in pdb_clean (and in the
            # frame offset) but are not given to Open Babel, which has no
            # Gasteiger parameters for them and would skip the whole receptor;
            # they are appended with their Vina metal types afterwards.
            ions = self._metal_ions() if keep_hetero else []
            with tempfile.TemporaryDirectory(prefix="autodocker_receptor_") as tmp:
                source = self.pdb_clean
                if ions:
                    # Same file name, so the PDBQT's REMARK Name is unchanged.
                    source = os.path.join(tmp, os.path.basename(self.pdb_clean))
                    _write_without_residues(self.pdb_clean, source, set(ions))
                self._protonate(ph, source)
                # obabel -p re-perceives residues; put the input's labels back.
                relabelled = protonation._restore_residue_labels(self.receptor_pdbqt, source)
            logger.info(f"[*] Restored input residue labels on {relabelled} receptor atoms")
            # Trp HE1 and His tautomers (residue_templates): in the input frame,
            # after label restoration, before the translation to the docking frame.
            self.residue_templates = residue_templates.apply_residue_templates(
                self.receptor_pdbqt, self.pdb_file)
            # Centre on the heavy atoms of pdb_clean (the frame offset written
            # to frame_offset.txt).
            _translate_pdbqt(self.receptor_pdbqt, self._receptor_centroid())
            if ions:
                appended = pocket._append_hetatm_to_receptor(
                    self.receptor_pdbqt, self.pdb_clean, ions, selected_chains,
                    offset=self.frame_offset())
                logger.info(
                    f"[*] --keep-hetero: added {appended} metal ion(s) with their Vina "
                    "metal types (not charged by Open Babel): "
                    + ", ".join(f"{n} {c}{r}" for c, r, n in ions))

            # Receptor PDBQT must not contain ligand torsion tags (ROOT/BRANCH/TORSDOF).
            charges._sanitize_receptor_pdbqt(self.receptor_pdbqt)

            if not charges._ensure_pdbqt_has_charges(self.receptor_pdbqt):
                raise RuntimeError(
                    f"Open Babel assigned no partial charges to {self.receptor_pdbqt}")
            n_hd, n_polar = protonation._count_polar_atoms(self.receptor_pdbqt)
            if n_hd == 0 and n_polar > 0:
                raise RuntimeError(
                    f"Receptor {self.receptor_pdbqt} has no polar hydrogens after "
                    f"protonation at pH {ph:g}; Vina would see no H-bond donors")

            logger.info(f"[*] Receptor prepared: {self.receptor_pdbqt} "
                        f"(protonated at pH {ph:g}, {n_hd} polar hydrogens)")
        except Exception as e:
            logger.error(f"Failed to prepare receptor: {e}")
            raise

        return selected_chains

    def rebuild_receptor_keep_hetatm(self, keep_residues) -> str:
        """Append the selected HETATM residues to the receptor PDBQT.

        ``keep_residues`` is a sequence of ``(chain, resnum, resname)``
        triples identifying the HETATM residues to retain (waters, metal
        ions, cofactors). Appending (rather than re-converting the whole
        receptor) preserves the protein charges already assigned by
        OpenBabel. Returns the receptor PDBQT path.
        """
        if not os.path.exists(self.receptor_pdbqt):
            raise FileNotFoundError(
                f"Receptor PDBQT not found: {self.receptor_pdbqt}")

        selected_chains = getattr(self, "selected_chains", None) or [BLANK_CHAIN_ID]
        appended = pocket._append_hetatm_to_receptor(
            self.receptor_pdbqt, self.pdb_file, keep_residues, selected_chains,
            offset=self.frame_offset())

        if appended == 0:
            logger.warning(
                "[!] No HETATM atoms were added to the receptor (none found, or "
                "all already present)")
            return self.receptor_pdbqt

        charges._sanitize_receptor_pdbqt(self.receptor_pdbqt)
        logger.info(
            f"[*] Receptor updated keeping {appended} HETATM atoms: {self.receptor_pdbqt}")
        return self.receptor_pdbqt

    def _box_source(self) -> str:
        if self.receptor_pdbqt and os.path.exists(self.receptor_pdbqt):
            return self.receptor_pdbqt
        return self.pdb_clean if os.path.exists(self.pdb_clean) else self.pdb_file

    def centroid_grid(self, padding: float = DEFAULT_PADDING
                      ) -> Tuple[float, float, float, float, float, float]:
        """Box around the whole receptor (used with --no-fpocket).

        Centred on the midpoint of the receptor's heavy-atom bounding box,
        ``padding`` on each side and every edge at least
        ``MIN_RECEPTOR_BOX_SIZE``. Raises if the receptor cannot be read.

        Public API (formerly the private ``_protein_centroid_grid``).
        """
        source = self._box_source()
        try:
            xs, ys, zs = _heavy_atom_coords(source)
        except (ValueError, IOError) as e:
            raise RuntimeError(f"Cannot build the receptor box from {source}: {e}") from e
        cx, cy, cz, sx, sy, sz = gridbox.bounding_box(xs, ys, zs, padding,
                                                      MIN_RECEPTOR_BOX_SIZE)
        logger.info(
            f"[*] Whole-receptor grid centred on the receptor bounding box "
            f"({cx:.2f}, {cy:.2f}, {cz:.2f}), size ({sx:.1f}, {sy:.1f}, {sz:.1f})")
        return cx, cy, cz, sx, sy, sz

    def detect_pocket(self, pocket_spec: Optional[str] = None,
                      padding: float = DEFAULT_PADDING
                      ) -> Tuple[float, float, float, float, float, float]:
        """Detect binding pockets with fpocket and build a grid from one or
        more of them. Raises when fpocket is unavailable, fails or finds no
        pocket (use ``centroid_grid`` for a deliberate whole-receptor box)."""
        logger.info("[*] Running fpocket...")
        fpocket_target = self.pdb_clean if os.path.exists(
            self.pdb_clean) else self.pdb_file
        pocket_root = fpocket.run_fpocket(fpocket_target)
        pockets = fpocket.rank_pockets(fpocket.collect_pockets(pocket_root, fpocket_target))
        fpocket.write_pocket_summary(self.pocket_summary_file, pockets)
        fpocket.log_top_pockets(pockets)

        selected = fpocket.select_pockets(pocket_spec, pockets)
        selected_names = ", ".join(f"pocket{p['number']}" for p in selected)
        logger.info(
            f"[*] Pocket selection: {selected_names} (padding={padding:g} A)")
        cx, cy, cz, sx, sy, sz = gridbox.pocket_box(
            [p["path"] for p in selected], padding=padding)
        dx, dy, dz = self._receptor_centroid()
        return cx - dx, cy - dy, cz - dz, sx, sy, sz

    def frame_offset(self) -> Tuple[float, float, float]:
        """Translation from the receptor PDBQT/docking frame to the input PDB frame.

        Add it to docked-pose coordinates to overlay them on the input PDB;
        subtract it to move input-PDB coordinates into the docking frame.
        """
        return self._receptor_centroid()

    def write_frame_offset(self) -> Optional[str]:
        """Record the frame offset in ``frame_offset.txt`` next to grid.conf."""
        path = os.path.join(self.outdir, outputs.FRAME_OFFSET_TXT)
        try:
            gridbox.write_frame_offset(path, *self.frame_offset())
        except (ValueError, IOError) as e:
            logger.warning(f"[!] Could not write frame offset: {e}")
            return None
        return path

    def _receptor_centroid(self) -> Tuple[float, float, float]:
        """Unweighted mean of the cleaned receptor's heavy atoms.

        ``prepare_receptor`` translates the receptor PDBQT by this vector so
        it is centred at the origin, while ``pdb_clean`` and fpocket's pocket
        files keep the original coordinates; the fpocket-derived grid centre
        is shifted by it to land in the receptor PDBQT frame.
        """
        source = self.pdb_clean if os.path.exists(
            self.pdb_clean) else self.pdb_file
        xs, ys, zs = _heavy_atom_coords(source)
        return sum(xs) / len(xs), sum(ys) / len(ys), sum(zs) / len(zs)

    def write_grid(self, cx: float, cy: float, cz: float,
                   sx: float, sy: float, sz: float):
        """Write AutoDock Vina grid configuration, the PyMOL box script and
        the frame offset."""
        try:
            gridbox.write_grid_conf(self.grid_conf, cx, cy, cz, sx, sy, sz)
            logger.info(f"[*] Grid ready: {self.grid_conf}")
            self.write_grid_box_script(cx, cy, cz, sx, sy, sz)
            self.write_frame_offset()
        except IOError as e:
            logger.error(f"Failed to write grid file: {e}")
            raise

    def write_grid_box_script(self, cx: float, cy: float, cz: float,
                              sx: float, sy: float, sz: float):
        """Write a PyMOL script that visualizes the final docking box."""
        gridbox.write_grid_box_script(self.grid_box_script, self.receptor_pdbqt,
                                      cx, cy, cz, sx, sy, sz)

    # Backward-compatible alias (deprecated): the former private name.
    _protein_centroid_grid = centroid_grid

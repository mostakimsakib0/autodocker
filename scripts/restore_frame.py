#!/usr/bin/env python3
"""Move docked poses back into the input PDB's coordinate frame.

The pipeline centers the receptor at the origin (``obabel -c``), so
receptor.pdbqt, grid.conf and every docked pose share that centered frame.
This is fine for the bundled PyMOL script and HTML report, but poses will not
overlay on the original protein.pdb (e.g. for MD set-up or comparison with a
crystal ligand). The pipeline records the translation in
``<output_dir>/frame_offset.txt``; this script adds it to every ATOM/HETATM
record and writes ``<name>_inputframe.<ext>`` next to each input file.

Usage:
    python3 scripts/restore_frame.py works/output works/output/docked/*_out.pdbqt
"""

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _benchlib import FRAME_OFFSET_TXT, read_frame_offset  # noqa: E402


def read_offset(output_dir):
    """(x, y, z) offset recorded in ``<output_dir>/frame_offset.txt``."""
    return read_frame_offset(os.path.join(output_dir, FRAME_OFFSET_TXT))


def shift_file(src, dst, offset):
    dx, dy, dz = offset
    with open(src, encoding="utf-8") as fin, open(dst, "w", encoding="utf-8") as fout:
        for line in fin:
            if line.startswith(("ATOM", "HETATM")) and len(line) >= 54:
                x = float(line[30:38]) + dx
                y = float(line[38:46]) + dy
                z = float(line[46:54]) + dz
                line = f"{line[:30]}{x:8.3f}{y:8.3f}{z:8.3f}{line[54:]}"
            fout.write(line)


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("output_dir", help="pipeline output directory (contains frame_offset.txt)")
    parser.add_argument("files", nargs="+", help="PDBQT/PDB files in the docking frame")
    args = parser.parse_args()

    try:
        offset = read_offset(args.output_dir)
    except (OSError, KeyError, ValueError) as e:
        sys.exit(f"[!] Cannot read frame offset from {args.output_dir}: {e}")

    for src in args.files:
        root, ext = os.path.splitext(src)
        dst = f"{root}_inputframe{ext}"
        shift_file(src, dst, offset)
        print(dst)


if __name__ == "__main__":
    main()

import os
import shutil
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "autodocker"))

import runner


def _write(path, content):
    with open(path, "w") as f:
        f.write(content)
    return path


@pytest.fixture
def outdir(tmp_path):
    d = tmp_path / "out"
    d.mkdir(exist_ok=True)
    return str(d)


def _atom_line(serial, x, y, z, chain="A", resname="ALA", atom="CA",
               element="C", charge="0.100"):
    """Build a PDBQT-compatible ATOM line with a columnar charge."""
    return (
        f"ATOM  {serial:5d} {atom:>4s} {resname:3s} {chain:1s}{serial:4d}    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00 20.00      {element:>2s} {charge:>6s}\n"
    )


@pytest.fixture
def protein_pdb(tmp_path):
    """A synthetic PDB with >100 ATOM lines across two chains."""
    lines = ["REMARK synthetic test protein\n"]
    serial = 0
    for chain in ("A", "B"):
        for i in range(70):
            serial += 1
            x = float(i % 10)
            y = float(i // 10)
            z = 5.0
            lines.append(_atom_line(serial, x, y, z, chain=chain))
    lines.append("END\n")
    return _write(str(tmp_path / "protein.pdb"), "".join(lines))


@pytest.fixture
def tiny_protein_pdb(tmp_path):
    """A PDB with fewer than 100 ATOM lines (must be rejected)."""
    lines = [
        _atom_line(i, float(i), 0.0, 0.0) for i in range(1, 10)
    ]
    return _write(str(tmp_path / "tiny.pdb"), "".join(lines))


SDF_METHANE = """methane
made-by-tests

  5  4  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0870    0.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.5435   -0.9414    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.5435    0.4707    0.8153 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.5435    0.4707   -0.8153 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  1  3  1  0  0  0  0
  1  4  1  0  0  0  0
  1  5  1  0  0  0  0
M  END
$$$$
"""


@pytest.fixture
def small_sdf(tmp_path):
    return _write(str(tmp_path / "methane.sdf"), SDF_METHANE)


@pytest.fixture
def charged_pdbqt(tmp_path):
    """A PDBQT with nonzero charges (B3 'good' case)."""
    content = "\n".join([
        _atom_line(1, 0.0, 0.0, 0.0, element="C", charge="0.120"),
        _atom_line(2, 1.0, 0.0, 0.0, element="H", charge="0.000"),
        "END",
    ])
    return _write(str(tmp_path / "charged.pdbqt"), content)


@pytest.fixture
def zero_charge_pdbqt(tmp_path):
    """A PDBQT whose atoms all carry 0.000 charge (B3 regression)."""
    content = "\n".join([
        _atom_line(1, 0.0, 0.0, 0.0, element="C", charge="0.000"),
        _atom_line(2, 1.0, 0.0, 0.0, element="H", charge="0.000"),
        "END",
    ])
    return _write(str(tmp_path / "zero.pdbqt"), content)


@pytest.fixture
def empty_pdbqt(tmp_path):
    return _write(str(tmp_path / "empty.pdbqt"), "")


@pytest.fixture
def mol2_file(tmp_path):
    return _write(str(tmp_path / "methane.mol2"), """@<TRIPOS>MOLECULE
methane
 5 4 1 0 0
SMALL
NO_CHARGES
****
@<TRIPOS>ATOM
      1 C1          0.0000    0.0000    0.0000 C.3     1  MOL1     0.0000
      2 H1         -1.0870    0.0000    0.0000 H       1  MOL1     0.0000
      3 H2          0.5435   -0.9414    0.0000 H       1  MOL1     0.0000
      4 H3          0.5435    0.4707    0.8153 H       1  MOL1     0.0000
      5 H4          0.5435    0.4707   -0.8153 H       1  MOL1     0.0000
@<TRIPOS>BOND
     1     1     2    1
     2     1     3    1
     3     1     4    1
     4     1     5    1
""")
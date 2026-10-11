#!/usr/bin/env python3
"""Deterministic 3D coordinates for 2D ligand input.

Open Babel's ``--gen3d`` draws from a random generator in its
conformer/rotor search: Open Babel 3.1 seeds it from the clock, so
heavy-atom coordinates differ between runs of the same 2D input (hydrogen
re-placement cannot fix that); 3.2 seeds it from OB_RANDOM_SEED, which the
pipeline sets, but the result still depends on the seed and on the Open
Babel build. When RDKit is importable, ``generate_3d`` embeds with ETKDGv3
using a fixed random seed and one thread, then relaxes with MMFF94 (UFF
when MMFF has no parameters) for a fixed number of iterations; the output
is byte-identical for identical input. Without RDKit it falls back to
``obabel --gen3d`` and logs a warning.
"""
import logging
from typing import Callable, Optional
from vspipeline import defaults  # noqa: E402

logger = logging.getLogger(__name__)

# Seed for ETKDG embedding; fixed so identical input gives identical output.
EMBED_SEED = defaults.SEED
# Force-field iterations after embedding (a fixed count, not a convergence
# criterion that could depend on timing).
FORCE_FIELD_MAX_ITERS = 2000

METHOD_RDKIT = "rdkit-etkdgv3-mmff94"
METHOD_OBABEL = "obabel-gen3d"


def rdkit_available() -> bool:
    try:
        import rdkit  # noqa: F401
        from rdkit.Chem import AllChem  # noqa: F401
    except ImportError:
        return False
    return True


def _rdkit_version() -> str:
    import rdkit
    return rdkit.__version__


def embed_and_relax(mol, seed: int = EMBED_SEED,
                    max_iters: int = FORCE_FIELD_MAX_ITERS):
    """The deterministic conformer recipe, shared by the pipeline and the
    re-docking benchmark (scripts/_benchlib.py): ETKDGv3 with a fixed seed on
    one thread (retried from random coordinates if the default embedding
    fails), then MMFF94, or UFF when MMFF lacks parameters, for a fixed number
    of iterations. ``mol`` must already carry the hydrogens wanted; it gets
    one conformer in place. Returns (force_field, converged); raises
    ValueError if embedding fails."""
    from rdkit.Chem import AllChem

    params = AllChem.ETKDGv3()
    params.randomSeed = int(seed)
    params.numThreads = 1
    if AllChem.EmbedMolecule(mol, params) != 0:
        params.useRandomCoords = True
        if AllChem.EmbedMolecule(mol, params) != 0:
            raise ValueError("RDKit ETKDGv3 could not embed the molecule")
    if AllChem.MMFFHasAllMoleculeParams(mol):
        status = AllChem.MMFFOptimizeMolecule(mol, maxIters=max_iters)
        return "MMFF94", status == 0
    status = AllChem.UFFOptimizeMolecule(mol, maxIters=max_iters)
    return "UFF", status == 0


def embed_3d_rdkit(src_sdf: str, dst_sdf: str, seed: int = EMBED_SEED,
                   max_iters: int = FORCE_FIELD_MAX_ITERS) -> int:
    """Write 3D coordinates for every molecule of ``src_sdf`` to ``dst_sdf``
    (all hydrogens explicit). Raises ValueError if a molecule cannot be read
    or embedded. Returns the number of molecules written."""
    from rdkit import Chem

    supplier = Chem.SDMolSupplier(src_sdf, removeHs=False)
    molecules = []
    for index, mol in enumerate(supplier):
        if mol is None:
            raise ValueError(f"RDKit cannot read molecule {index + 1} of {src_sdf}")
        mol = Chem.AddHs(mol)
        try:
            embed_and_relax(mol, seed=seed, max_iters=max_iters)
        except ValueError:
            raise ValueError(
                f"RDKit ETKDG could not embed molecule {index + 1} of {src_sdf}") from None
        molecules.append(mol)
    if not molecules:
        raise ValueError(f"No molecules in {src_sdf}")
    writer = Chem.SDWriter(dst_sdf)
    try:
        for mol in molecules:
            writer.write(mol)
    finally:
        writer.close()
    return len(molecules)


def generate_3d(src_sdf: str, dst_sdf: str, seed: int = EMBED_SEED,
                run_obabel: Optional[Callable] = None) -> str:
    """3D coordinates for a 2D SDF; returns the method used (``METHOD_RDKIT``
    or ``METHOD_OBABEL``) so the caller can record it.

    ``run_obabel(args, outputs=[...], seed=...)`` runs the Open Babel
    fallback (default: ``utils.run_obabel``, which runs obabel with
    OMP_NUM_THREADS=1), with ``OB_RANDOM_SEED`` = ``seed``.
    """
    if rdkit_available():
        embed_3d_rdkit(src_sdf, dst_sdf, seed=seed)
        logger.debug(f"3D coordinates for {src_sdf}: RDKit {_rdkit_version()} ETKDGv3 "
                     f"(seed {seed}) + force field ({FORCE_FIELD_MAX_ITERS} iterations)")
        return METHOD_RDKIT
    logger.warning(
        f"[!] RDKit not available: 3D coordinates for {src_sdf} come from obabel "
        "--gen3d, which is NOT reproducible with Open Babel < 3.2 (heavy atoms "
        "differ between runs; 3.2+ seeds it from OB_RANDOM_SEED) and differs "
        "from the RDKit coordinates of the container image")
    from vspipeline.hydrogens import _default_obabel
    run = run_obabel or _default_obabel
    run(["-isdf", src_sdf, "-osdf", "-O", dst_sdf, "--gen3d"], outputs=[dst_sdf], seed=seed)
    return METHOD_OBABEL

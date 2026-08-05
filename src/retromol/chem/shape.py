"""Module for principal-moments-of-inertia (PMI) molecular shape analysis.

See https://pubs.acs.org/doi/full/10.1021/ci025599w for the normalized-PMI-ratio
convention used here (the "shape triangle": rod/disc/sphere-like character).
"""

from dataclasses import dataclass

from rdkit import Chem

from retromol.chem.geom import calculate_moments_of_inertia_and_axes, get_conformer_energy, get_conformers
from retromol.chem.mol import standardize_from_smiles

# Above this heavy-atom count, conformer embedding/optimization is too slow to run
# synchronously inside a shared request thread -- reject rather than hang a request.
MAX_HEAVY_ATOMS = 200


@dataclass(frozen=True)
class PmiShape:
    """Normalized principal moments of inertia for one molecule's lowest-energy conformer."""

    npr1: float
    npr2: float
    conformers_generated: int


def compute_pmi_shape(smiles: str, num_confs: int = 25, max_iters: int = 500) -> PmiShape:
    """
    Compute normalized PMI ratios (NPR1 = I1/I3, NPR2 = I2/I3) for a molecule's
    lowest-energy conformer.

    Generates `num_confs` conformers (MMFF-optimized), picks the lowest-energy one,
    and computes its mass-weighted principal moments of inertia.

    :param smiles: the input SMILES string
    :param num_confs: how many conformers to embed and optimize
    :param max_iters: max MMFF optimization iterations per conformer
    :return: the computed PmiShape
    :raises ValueError: if the SMILES is invalid, the molecule is too large, or no conformers could be embedded
    """
    mol = standardize_from_smiles(smiles, keep_stereo=True, neutralize=True, tautomer_canon=False)
    if mol is None:
        raise ValueError(f"could not parse or standardize SMILES: {smiles}")

    if mol.GetNumAtoms() > MAX_HEAVY_ATOMS:
        raise ValueError(f"molecule has {mol.GetNumAtoms()} heavy atoms, exceeding the {MAX_HEAVY_ATOMS} limit")

    mol = Chem.AddHs(mol)

    confs = get_conformers(
        mol,
        num_confs=num_confs,
        use_random_coords=True,
        optimize_confs=True,
        max_iters_optimization=max_iters,
    )
    if not confs:
        raise ValueError("could not embed any conformers")

    conf = min(confs, key=get_conformer_energy)

    coords = conf.GetPositions()
    weights = [atom.GetMass() for atom in mol.GetAtoms()]

    moments, _ = calculate_moments_of_inertia_and_axes(coords, weights)

    return PmiShape(
        npr1=float(moments[0] / moments[2]),
        npr2=float(moments[1] / moments[2]),
        conformers_generated=len(confs),
    )

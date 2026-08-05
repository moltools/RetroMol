"""Module for shape cheminformatics: conformer generation and principal moments of inertia."""

import typing as ty

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from scipy.linalg import eig


def get_conformers(
    mol: Chem.Mol,
    num_confs: int = 10,
    max_attempts: int = 0,
    random_seed: int = -1,
    clear_confs: bool = True,
    use_random_coords: bool = False,
    box_size_mult: float = 2.0,
    random_neg_eig: bool = True,
    num_zero_fail: int = 1,
    prune_rms_thresh: float = -1.0,
    coord_map: ty.Optional[ty.Dict[int, ty.Tuple[float, float, float]]] = {},
    force_tol: float = 1e-3,
    ignore_smoothing_failures: bool = False,
    enforce_chirality: bool = True,
    num_threads: int = 1,
    use_exp_torsion_angle_prefs: bool = True,
    use_basic_knowledge: bool = True,
    print_exp_torsion_angles: bool = False,
    use_small_ring_torsions: bool = False,
    use_macrocycle_torsions: bool = False,
    et_version: int = 1,
    optimize_confs: bool = True,
    max_iters_optimization: int = 200,
    mmff_variant: str = "MMFF94",
    non_bonded_thresh: float = 100.0,
    ignore_interfrag_interactions: bool = True,
) -> ty.List[Chem.Conformer]:
    """
    Embed and optionally optimize conformers for a molecule.

    :param mol: The molecule to embed conformers for.
    :param num_confs: The number of conformers to embed.
    :param max_attempts: The maximum number of embedding attempts. If 0, defaults to num_confs.
    :param random_seed: The random seed to use for embedding. If -1, a random seed is generated.
    :param clear_confs: Whether to clear the conformers before embedding.
    :param use_random_coords: Whether to use random coordinates for embedding (vs. eigenvalues of the distance matrix).
    :param box_size_mult: The box size multiplier for embedding.
    :param random_neg_eig: Whether to use random eigenvalues for embedding.
    :param num_zero_fail: The number of zero eigenvalue failures allowed.
    :param prune_rms_thresh: The RMS threshold for pruning near-duplicate conformers.
    :param coord_map: The coordinate map for embedding.
    :param force_tol: The force tolerance for embedding.
    :param ignore_smoothing_failures: Whether to ignore smoothing failures.
    :param enforce_chirality: Whether to enforce chirality.
    :param num_threads: The number of threads to use for embedding and optimization.
    :param use_exp_torsion_angle_prefs: Whether to use experimental torsion angle preferences.
    :param use_basic_knowledge: Whether to use basic knowledge.
    :param print_exp_torsion_angles: Whether to print experimental torsion angles.
    :param use_small_ring_torsions: Whether to use small ring torsions.
    :param use_macrocycle_torsions: Whether to use macrocycle torsions.
    :param et_version: The experimental torsion angle version.
    :param optimize_confs: Whether to optimize the conformers.
    :param max_iters_optimization: The maximum number of iterations for optimization.
    :param mmff_variant: The MMFF variant to use for optimization.
    :param non_bonded_thresh: The non-bonded threshold for optimization.
    :param ignore_interfrag_interactions: Whether to ignore inter-fragment interactions.
    :return: The conformers of the molecule.
    """
    conf_ids = AllChem.EmbedMultipleConfs(
        mol,
        numConfs=num_confs,
        maxAttempts=max_attempts,
        randomSeed=random_seed,
        clearConfs=clear_confs,
        useRandomCoords=use_random_coords,
        boxSizeMult=box_size_mult,
        randNegEig=random_neg_eig,
        numZeroFail=num_zero_fail,
        pruneRmsThresh=prune_rms_thresh,
        coordMap=coord_map,
        forceTol=force_tol,
        ignoreSmoothingFailures=ignore_smoothing_failures,
        enforceChirality=enforce_chirality,
        numThreads=num_threads,
        useExpTorsionAnglePrefs=use_exp_torsion_angle_prefs,
        useBasicKnowledge=use_basic_knowledge,
        printExpTorsionAngles=print_exp_torsion_angles,
        useSmallRingTorsions=use_small_ring_torsions,
        useMacrocycleTorsions=use_macrocycle_torsions,
        ETversion=et_version,
    )

    if optimize_confs:
        AllChem.MMFFOptimizeMoleculeConfs(
            mol,
            numThreads=num_threads,
            maxIters=max_iters_optimization,
            mmffVariant=mmff_variant,
            nonBondedThresh=non_bonded_thresh,
            ignoreInterfragInteractions=ignore_interfrag_interactions,
        )

    conformers = [mol.GetConformer(conf_id) for conf_id in conf_ids]

    return conformers


def get_conformer_energy(conf: Chem.Conformer) -> float:
    """
    Get the MMFF94 energy of a conformer.

    :param conf: The conformer to get the energy of.
    :return: The energy of the conformer.
    """
    energy = AllChem.MMFFGetMoleculeForceField(
        conf.GetOwningMol(),
        AllChem.MMFFGetMoleculeProperties(
            conf.GetOwningMol(),
            mmffVariant="MMFF94"
        ),
        confId=conf.GetId(),
    ).CalcEnergy()

    return energy


def calculate_center_of_mass(coords: np.ndarray, weights: ty.List[float]) -> np.ndarray:
    """
    Calculate the center of mass of a set of coordinates.

    :param coords: The coordinates to calculate the center of mass of.
    :param weights: The weights of the coordinates.
    :return: The center of mass of the coordinates.
    """
    center_of_mass = np.average(coords, axis=0, weights=weights)

    return center_of_mass


def calculate_moments_of_inertia_and_axes(
    coords: np.ndarray,
    weights: ty.List[float]
) -> ty.Tuple[np.ndarray, np.ndarray]:
    """
    Calculate the moment of inertia tensor and principal axes.

    :param coords: The coordinates to calculate the moment of inertia tensor and principal axes of, shape (n, 3).
    :param weights: The weights of the coordinates, length n.
    :return: A tuple of (principal moments of inertia, sorted ascending; principal axes).
    """
    if not isinstance(coords, np.ndarray):
        raise TypeError("coords must be a numpy array.")
    if coords.shape[1] != 3:
        raise ValueError("coords must be of shape (n, 3).")

    if not isinstance(weights, list):
        raise TypeError("weights must be a list.")
    if len(weights) != coords.shape[0]:
        raise ValueError("weights must be of length n.")

    center_of_mass = calculate_center_of_mass(coords, weights)

    coords = coords - center_of_mass

    tensor = np.zeros((3, 3))
    for i, (x, y, z) in enumerate(coords):
        tensor[0, 0] += weights[i] * (y**2 + z**2)
        tensor[1, 1] += weights[i] * (x**2 + z**2)
        tensor[2, 2] += weights[i] * (x**2 + y**2)
        tensor[0, 1] -= weights[i] * x * y
        tensor[0, 2] -= weights[i] * x * z
        tensor[1, 2] -= weights[i] * y * z
    tensor[1, 0] = tensor[0, 1]
    tensor[2, 0] = tensor[0, 2]
    tensor[2, 1] = tensor[1, 2]

    eigvals, eigvecs = eig(tensor)
    idx = eigvals.argsort()
    eigvals = eigvals[idx]
    eigvecs = eigvecs[:, idx]

    return eigvals, eigvecs

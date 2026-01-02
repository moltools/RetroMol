"""Module defining the Submission data class."""

from dataclasses import dataclass, field
from typing import Any

from rdkit.Chem.rdchem import Mol

from retromol.chem.mol import standardize_from_smiles, mol_to_inchikey
from retromol.chem.tagging import tag_mol


@dataclass(frozen=True)
class Submission:
    """
    Represents a submission of a molecule for retrosynthetic analysis.

    :var smiles: str: SMILES representation of the submitted molecule
    :var name: str | None: optional name of the submitted molecule
    :var props: dict[str, Any] | None: optional additional properties associated with the submission
    """

    smiles: str
    name: str | None = None
    props: dict[str, Any] | None = None

    mol: Mol = field(init=False, repr=False)
    inchikey: str = field(init=False, repr=False)

    def __post_init__(self) -> None:
        """
        """
        # Sanitize SMILES
        smiles = self.smiles.replace("[N]", "N")  # avoid parsing issues with RDKit

        # Generate standardized molecule
        mol = standardize_from_smiles(
            smiles,
            keep_stereo=True,
            neutralize=True,
            tautomer_canon=True,
        )

        # Generate InChIKey
        inchikey = mol_to_inchikey(self.mol)

        # Tag molecule
        tag_mol(mol)

        object.__setattr__(self, "smiles", smiles)
        object.__setattr__(self, "mol", mol)
        object.__setattr__(self, "inchikey", inchikey)
